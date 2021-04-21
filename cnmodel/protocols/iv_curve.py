from __future__ import print_function
from neuron import h
import numpy as np
import scipy
import scipy.integrate
import scipy.stats
import scipy.optimize

try:
    import pyqtgraph as pg

    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from ..util.stim import make_pulse
from ..util import fitting
from ..util import custom_init
from .protocol import Protocol


class IVCurve(Protocol):
    def __init__(self):
        super(IVCurve, self).__init__()
        self.reset()

    def reset(self):
        super(IVCurve, self).reset()
        self.voltage_traces = []
        self.durs = None  # durations of current steps
        self.current_cmd = None  # Current command levels
        self.current_traces = []
        self.axon_voltage_traces = None
        self.time_values = None
        self.dt = None
        self.initdelay = 0.0

    def run(
        self,
        ivrange,
        cell,
        durs=None,
        sites=None,
        reppulse=None,
        temp=22,
        dt=0.025,
        initdelay=0.0,
        stimdict=None,
    ):
        """
        Run a current-clamp I/V curve on *cell*.
        
        Parameters
        ----------
        ivrange : dict of list of tuples
            Each item in the list is (min, max, step) describing a range of 
            levels to test. Range values are inclusive, so the max value may
            appear in the test values. Using multiple ranges allows finer 
            measurements in some ranges.
            For example::
    
                {'pulse': [(-1., 0., 1.), (-0.1, 0., 0.02)], 'prepulse': [(-0.5, 0, 0.1)]}
            
            will generate a set of negative pulses between -1 and 0 nA, AND
            a set of negative pulses from -0.1 to 0 nA in steps of 0.1 nA, AND
            each of these will be prefaced by a set of prepulses, in the range between
            -0.5 and 0 nA, in steps of 0.1 nA
    
            Optional keys include:
    
            * 'pulsedur' : the duration of the pulse, in ms
            * 'prepulsecur: the duration of the prepulse, in ms
            
            The prepulse or the pulse can have a single value if the other is ranged.
    
        cell : Cell
            The Cell instance to test.
        durs : tuple
            durations of (pre, pulse, post) regions of the command (in msec)
        sites : list
            Sections to add recording electrodes
        reppulse : 
            stimulate with pulse train
        temp : 
            temperature of simulation (32) (Celsius)
        dt : 
            timestep of simulation (0.025) (msec)
        stim: dict
            a dict with::
            
                stimdict = {
                    'NP': 10,  # number of pulses
                    'Sfreq': 400.0,  # stimulus frequency (pulse rate, Hz)
                    'delay': 10.0,  # delay to first pulses, msec
                    'dur': 1.0,  # duration of each pulse, msec
                    'amp': 1.0,  # amplitude of each pulse, nA
                    'PT': 0.0,  # post test pulse (after train), delay in msec
                }


        """

        # print('iv_curve:run')
        self.reset()
        self.cell = cell
        self.initdelay = initdelay
        self.dt = dt
        self.temp = temp
        # Calculate current pulse levels
        icur = []
        precur = [0.0]
        self.pre_current_cmd = []
        npresteps = 0
        if isinstance(ivrange["pulse"], tuple):
            icmd = [ivrange["pulse"]]  # convert to a list with tuple(s) embedded
        else:
            icmd = ivrange["pulse"]  # was already a list with multiple tuples
        for c in icmd:  # unpack current levels for the main pulse
            try:
                (imin, imax, istep) = c  # unpack a tuple... or list
            except:
                raise TypeError(
                    "run_iv arguments must be a dict with tuples {'pulse': (imin, imax, istep), 'prepulse': ...}"
                )
            nstep = np.floor((imax - imin) / istep) + 1
            icur.extend(imin + istep * np.arange(nstep))  # build the list
        self.current_cmd = np.array(sorted(icur))
        nsteps = self.current_cmd.shape[0]
        # Configure IClamp
        if durs is None:
            durs = [10.0, 100.0, 50.0]  # set default durs

        if "prepulse" in ivrange.keys():
            if isinstance(ivrange["prepulse"], tuple):
                icmd = [ivrange["prepulse"]]  # convert to a list with tuple(s) embedded
            else:
                icmd = ivrange["prepulse"]  # was already a list with multiple tuples
            precur = []
            for c in icmd:
                try:
                    (imin, imax, istep) = c  # unpack a tuple... or list
                except:
                    raise TypeError(
                        "run_iv arguments must be a dict with tuples {'pulse': (imin, imax, istep), 'prepulse': ...}"
                    )
                nstep = np.floor((imax - imin) / istep) + 1
                precur.extend(imin + istep * np.arange(nstep))  # build the list
            self.pre_current_cmd = np.array(sorted(precur))
            npresteps = self.pre_current_cmd.shape[0]
            durs.insert(1, 50.0)

        self.durs = durs
        # set up stimulation with a pulse train
        if reppulse is not None and stimdict is None:
            stim = {
                "NP": 10,
                "Sfreq": 400.0,
                "delay": 10.0,
                "dur": 1.0,
                "amp": 1.0,
                "PT": 0.0,
                "dt": self.dt,
            }
            self.p_start = stim['delay']
            self.p_dur = stim['dur']
            self.p_end = self.p_start + self.p_dur + 10.
        elif stimdict is not None:
            # print('stimdict: ', stimdict)
            stim = stimdict
            stim["dt"] = self.dt
        elif "prepulse" in ivrange.keys():
            stim = {
                "NP": 2,
                "delay": durs[0],
                "predur": durs[1],
                "dur": durs[2],
                "amp": 1.0,
                "preamp": 0.0,
                "dt": self.dt,
            }
            self.p_start = durs[0] + durs[1]
            self.p_end = self.p_start + durs[2]
            self.p_dur = durs[2]
        else:
            stim = {
                "NP": 1,
                "delay": durs[0],
                "dur": durs[1],
                "amp": 1.0,
                "dt": self.dt,
            }
            self.p_start = durs[0]
            self.p_end = self.p_start + durs[1]
            self.p_dur = durs[1]

        # print stim
        # print('p_: ', self.p_start, self.p_end, self.p_dur)
        istim = h.iStim(0.5, sec=cell.soma)
        istim.delay = 5.0
        istim.dur = 1e9  # these actually do not matter...
        istim.iMax = 0.0
        self.tend = np.sum(durs)  # maxt + len(iextend)*stim['dt']
        self.axon_voltage_traces = []

        self.cell = cell
        for i in range(nsteps):
            # Generate current command for this level
            stim["amp"] = self.current_cmd[i]
            if npresteps > 0:
                for j in range(npresteps):
                    stim["preamp"] = self.pre_current_cmd[j]
                    self.run_one(istim, stim, initflag=(i == 0 and j == 0), sites=sites)
            else:
                self.run_one(istim, stim, initflag=(i == 0), sites=sites)

    def run_one(self, istim, stim, initflag=True, sites=None):
        """
        Perform one run in current-clamp for the selected cell
        and add the data to the traces
        
        Parameters
        ----------
        istim : Stimulus electrode instance
        stim : waveform information
        initflag : boolean (default: True)
            If true, force initialziation of the cell and computation of 
            point Rin, tau and Vm
        sites : list of sections to monitor in addition to the soma
        """
        # print('iv_curve:run_one')
        (secmd, maxt, tstims) = make_pulse(stim)
        # print('maxt, dt*lencmd: ', maxt, len(secmd)*self.dt)# secmd = np.append(secmd, [0.])
        # print('stim: ', stim, self.tend)

        # connect current command vector
        playvector = h.Vector(secmd)
        playvector.play(istim._ref_i, h.dt, 0, sec=self.cell.soma)

        # Connect recording vectors
        self["v_soma"] = self.cell.soma(0.5)._ref_v
        # self['q10'] = self.cell.soma(0.5).ihpyr_adj._ref_q10
        # self['ih_ntau'] = self.cell.soma(0.5).ihpyr_adj._ref_kh_n_tau
        self["i_inj"] = istim._ref_i
        self["time"] = h._ref_t
        if sites is not None:
            recvec = [h.Vector() for x in sites]
            for i, r in enumerate(recvec):
                r.record(sites[i](0.5)._ref_v, h.dt, 0, sec=sites[i])
        # h('secondorder=0')  # direct call fails; let hoc do the work
        h.celsius = self.cell.status["temperature"]
        # print("iv_curve:run_one:calling cell_initialize")
        self.cell.cell_initialize()
        h.dt = self.dt
        # print("iv_curve:run_one:calling custom_init")
        custom_init(v_init=self.cell.vm0)

        h.t = 0.0
        h.tstop = self.tend
        while h.t < h.tstop:
            h.fadvance()

        self.voltage_traces.append(self["v_soma"])
        self.current_traces.append(self["i_inj"])
        self.time_values = np.array(self["time"] - self.initdelay)
        if sites is not None:
            # for i, r in enumerate(recvec):
            #     print('r2: ', np.array(recvec[i]))
            self.axon_voltage_traces.append(np.array(recvec.copy()))
            # print('# axon site traces: ', len(self.axon_voltage_traces))
        # self.mon_q10 = np.array(self['q10'])
        # self.mon_ih_ntau = np.array(self['ih_ntau'])

    def peak_vm(self, window=0.5):
        """
        Parameters
        ----------
        window : float (default: 0.5)
            fraction of trace to look at to find peak value
        
        Returns
        -------
        peak membrane voltage for each trace.

        """
        Vm = self.voltage_traces
        Icmd = self.current_cmd
        steps = len(Icmd)
        peakStart = int(self.p_start / self.dt)
        peakStop = int(
            peakStart + (self.p_dur * window) / self.dt
        )  # peak can be in first half
        Vpeak = []
        for i in range(steps):
            if Icmd[i] > 0:
                Vpeak.append(Vm[i][peakStart:peakStop].max())
            else:
                Vpeak.append(Vm[i][peakStart:peakStop].min())
        return np.array(Vpeak)

    def steady_vm(self, window=0.1):
        """
        Parameters
        ----------
        window: (float) default: 0.1
        fraction of window to use for steady-state measurement, taken
        immediately before the end of the step
        
        Returns
        -------
        steady-state membrane voltage for each trace.
        
        """
        Vm = self.voltage_traces
        steps = len(Vm)
        steadyStop = int((self.p_end) / self.dt)
        steadyStart = int(
            steadyStop - (self.p_end * window) / self.dt
        )  # measure last 10% of trace
        Vsteady = [Vm[i][steadyStart:steadyStop].mean() for i in range(steps)]
        return np.array(Vsteady)

    def spike_times(self, threshold=None):
        """
        Return an array of spike times for each trace.
        
        Parameters
        ----------
        threshold: float (default: None)
            Optional threshold at which to detect spikes. By 
            default, this queries cell.spike_threshold.
        
        Returns
        -------
        list of spike times.
        
        """
        if threshold is None:
            threshold = self.cell.spike_threshold

        Vm = self.voltage_traces
        steps = len(Vm)
        spikes = []
        for i in range(steps):
            # dvdt = np.diff(Vm[i]) / self.dt
            # mask = (dvdt > 40).astype(int)
            mask = (Vm[i] > threshold).astype(int)
            indexes = np.argwhere(np.diff(mask) == 1)[:, 0] + 1
            times = indexes.astype(float) * self.dt
            spikes.append(times)
        return spikes

    def spike_filter(self, spikes, window=(0.0, np.inf)):
        """Filter the spikes to only those occurring in a defined window.
        
        Required to compute input resistance in traces with no spikes during
        the stimulus, because some traces will have anodal break spikes.
        
        Parameters
        ----------
        spikes : list
            the list of spike trains returned from the spike_times method
        window : (start, stop)
            the window over which to look for spikes (in msec: default is
            the entire trace).

        Returns
        -------
        the spikes in a list
        
        """
        filteredspikes = []
        for i in range(len(spikes)):
            winspikes = []  # spikes is arranged by current; so this is for one level
            for j in range(len(spikes[i])):
                if spikes[i][j] >= window[0] and spikes[i][j] <= window[1]:
                    winspikes.append(spikes[i][j])
            filteredspikes.append(winspikes)  # now build filtered spike list
        return filteredspikes

    def rest_vm(self):
        """
        Parameters
        ----------
        None
        
        Returns
        -------
        The mean resting membrane potential.
        """
        d = int(self.durs[0] / self.dt)
        rvm = np.array(
            [
                np.array(self.voltage_traces[i][d // 2 : d]).mean()
                for i in range(len(self.voltage_traces))
            ]
        ).mean()
        return rvm

    def input_resistance_tau(self, vmin=-10.0, imax=0):
        """
        Estimate resting input resistance and time constant.
        
        Parameters
        ----------
        vmin : float
            minimum voltage to use in computation relative to resting
        imax : float
            maximum current to use in computation.
        return_eval : bool
            If True, return lmfit.ModelFit instances for the subthreshold trace
            fits as well.
            
        Returns
        -------
        dict :
            Dict containing:
            
            * 'slope' and 'intercept' keys giving linear 
              regression for subthreshold traces near rest
            * 'tau' giving the average first-exponential fit time constant
            * 'fits' giving a record array of exponential fit data to subthreshold
              traces.
        
        Analyzes only traces hyperpolarizing pulse traces near rest, with no 
        spikes.
        
        """
        self.fits = None
        Vss = self.steady_vm()
        vmin += self.rest_vm()
        Icmd = self.current_cmd
        rawspikes = self.spike_times()
        spikes = self.spike_filter(rawspikes, window=[self.p_start, self.p_end])
        steps = len(Icmd)

        nSpikes = np.array([len(s) for s in spikes])
        # find traces with Icmd < 0, Vm > -70, and no spikes.
        vmask = Vss >= vmin
        imask = Icmd <= imax
        smask = nSpikes > 0
        mask = vmask & imask & ~smask
        if mask.sum() < 2:
            print(
                "WARNING: Not enough traces to do linear regression in "
                "IVCurve.input_resistance_tau()."
            )
            print(
                "{0:<15s}: {1:s}".format(
                    "vss", ", ".join(["{:.2f}".format(v) for v in Vss])
                )
            )
            print(
                "{0:<15s}: {1:s}".format(
                    "Icmd", ", ".join(["{:.2f}".format(i) for i in Icmd])
                )
            )
            print("{0:<15s}: {1:s}".format("vmask", repr(vmask.astype(int))))
            print("{0:<15s}: {1:s} ".format("imask", repr(imask.astype(int))))
            print("{0:<15s}: {1:s}".format("spikemask", repr(smask.astype(int))))
            slope = np.nan
            intercept = np.nan
            tau = np.nan
            fit_data = []
            ret = {"slope": slope, "intercept": intercept, "tau": tau, "fits": fit_data}
            return ret

        # Use these to measure input resistance by linear regression.
        reg = scipy.stats.linregress(Icmd[mask], Vss[mask])
        (slope, intercept, r, p, stderr) = reg

        # also measure the tau in the same traces:
        pulse_start = int(self.p_start / self.dt)
        pulse_stop = int((self.p_end) / self.dt)

        fits = []
        fit_inds = []
        tx = self.time_values[pulse_start:pulse_stop].copy()
        for i, m in enumerate(mask):
            if not m or (self.rest_vm() - Vss[i]) <= 1:
                continue
            trace = self.voltage_traces[i][pulse_start:pulse_stop]

            # find first minimum in the trace
            min_ind = np.argmin(trace)
            min_val = trace[min_ind]
            min_diff = trace[0] - min_val
            tau_est = min_ind * self.dt * (1.0 - 1.0 / np.e)
            # print ('minind: ', min_ind, tau_est)
            fit = fitting.Exp1().fit(
                trace[:min_ind],
                method="nelder",
                x=tx[:min_ind],
                xoffset=(tx[0], "fixed"),
                yoffset=(min_val, -120.0, -10.0),
                amp=(min_diff, 0.0, 50.0),
                tau=(tau_est, 0.5, 50.0),
            )

            # find first maximum in the trace (following with first minimum)
            max_ind = np.argmax(trace[min_ind:]) + min_ind
            max_val = trace[max_ind]
            max_diff = min_val - max_val
            tau2_est = max_ind * self.dt * (1.0 - 1.0 / np.e)
            amp1_est = fit.params["amp"].value
            tau1_est = fit.params["tau"].value
            amp2_est = fit.params["yoffset"] - max_val
            # print('tau1, tau2est: ', tau1_est, tau2_est)
            # fit up to first maximum with double exponential, using prior
            # fit as seed.
            fit = fitting.Exp2().fit(
                trace[:max_ind],
                method="nelder",
                x=tx[:max_ind],
                xoffset=(tx[0], "fixed"),
                yoffset=(max_val, -120.0, -10.0),
                amp1=(amp1_est, 0.0, 200.0),
                tau1=(tau1_est, 0.5, 50.0),
                amp2=(amp2_est, -200.0, -0.5),
                tau_ratio=(tau2_est / tau1_est, 2.0, 50.0),
                tau2="tau_ratio * tau1",
            )

            fits.append(fit)
            fit_inds.append(i)
        # convert fits to record array
        #        print len(fits) # fits[0].params
        if len(fits) > 0:
            key_order = sorted(
                fits[0].params
            )  # to ensure that unit tests remain stable
            dtype = [(k, float) for k in key_order] + [("index", int)]
            fit_data = np.empty(len(fits), dtype=dtype)
            for i, fit in enumerate(fits):
                for k, v in fit.params.items():
                    fit_data[i][k] = v.value
                fit_data[i]["index"] = fit_inds[i]

            if "tau" in fit_data.dtype.fields:
                tau = fit_data["tau"].mean()
            else:
                tau = fit_data["tau1"].mean()
        else:
            slope = 0.0
            intercept = 0.0
            tau = 0.0
            fit_data = []
        ret = {"slope": slope, "intercept": intercept, "tau": tau, "fits": fit_data}
        self.fits = fits
        return ret

    def show(self, cell=None, rmponly=False):
        """
        Plot results from run_iv()
        
        Parameters
        ----------
        cell : cell object (default: None)
        
        """
        if not HAVE_PG:
            raise Exception("Requires pyqtgraph")

        #
        # Generate figure with subplots
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow(
            "%s  %s (%s)"
            % (cell.status["name"], cell.status["modelType"], cell.status["species"])
        )
        self.win = win
        win.resize(1000, 800)
        Vplot = win.addPlot(labels={"left": "Vm (mV)", "bottom": "Time (ms)"})
        rightGrid = win.addLayout(rowspan=2)
        win.nextRow()
        Iplot = win.addPlot(labels={"left": "Iinj (nA)", "bottom": "Time (ms)"})

        IVplot = rightGrid.addPlot(labels={"left": "Vm (mV)", "bottom": "Icmd (nA)"})
        IVplot.showGrid(x=True, y=True)
        rightGrid.nextRow()
        spikePlot = rightGrid.addPlot(
            labels={"left": "Iinj (nA)", "bottom": "Spike times (ms)"}
        )
        rightGrid.nextRow()
        FIplot = rightGrid.addPlot(
            labels={"left": "Spike count", "bottom": "Iinj (nA)"}
        )

        win.ci.layout.setRowStretchFactor(0, 10)
        win.ci.layout.setRowStretchFactor(1, 5)

        #
        # Plot simulation and analysis results
        #
        Vm = self.voltage_traces
        Iinj = self.current_traces
        Icmd = self.current_cmd
        DVm = self.axon_voltage_traces
        t = self.time_values
        steps = len(Icmd)
        # plot I, V traces
        colors = [(i, steps * 3.0 / 2.0) for i in range(steps)]
        for i in range(steps):
            Vplot.plot(t, Vm[i], pen=colors[i])
            Iplot.plot(t, Iinj[i], pen=colors[i])
            if len(DVm) == 0:
                continue
            nnodes = len(DVm[i])
            axcolors = [pg.intColor(k, nnodes) for k in range(nnodes)]
            for j in range(len(DVm[i])):  # for each site
                if j in range(nnodes):  # [1, nnodes-1]:
                    ilen = len(DVm[i][j])
                    Vplot.plot(t[:ilen], DVm[i][j], pen=axcolors[j])

        if rmponly:
            return
        # I/V relationships
        IVplot.plot(
            Icmd,
            self.peak_vm(),
            symbol="o",
            symbolBrush=(50, 150, 50, 255),
            symbolSize=4.0,
        )
        IVplot.plot(Icmd, self.steady_vm(), symbol="s", symbolSize=4.0)

        # F/I relationship and raster plot
        spikes = self.spike_times()
        for i, times in enumerate(spikes):
            spikePlot.plot(
                x=times,
                y=[Icmd[i]] * len(times),
                pen=None,
                symbol="d",
                symbolBrush=colors[i],
                symbolSize=4.0,
            )
        FIplot.plot(x=Icmd, y=[len(s) for s in spikes], symbol="o", symbolSize=4.0)

        # Print Rm, Vrest
        rmtau = self.input_resistance_tau()
        if self.fits is None:
            return
        s = rmtau["slope"]
        i = rmtau["intercept"]
        # print(s)
        # tau1 = rmtau['fits']['tau1'].mean()
        # tau2 = rmtau['fits']['tau2'].mean()
        # print ("\nMembrane resistance (chord): {0:0.1f} MOhm  Taum1: {1:0.2f}  Taum2: {2:0.2f}".format(s, tau1, tau2))

        # Plot linear subthreshold I/V relationship
        ivals = np.array([Icmd.min(), Icmd.max()])
        vvals = s * ivals + i
        line = pg.QtGui.QGraphicsLineItem(ivals[0], vvals[0], ivals[1], vvals[1])
        line.setPen(pg.mkPen(255, 0, 0, 70))
        line.setZValue(-10)
        IVplot.addItem(line, ignoreBounds=True)

        # plot exponential fits
        for fit in fits:
            t = np.linspace(self.p_start, self.p_end, 1000)
            y = fit.eval(x=t)
            Vplot.plot(
                t, y, pen={"color": (100, 100, 0), "style": pg.QtCore.Qt.DashLine}
            )

            # plot initial guess
            # y = fit.eval(x=t, **fit.init_params.valuesdict())
            # Vplot.plot(t, y, pen={'color': 'b', 'style': pg.QtCore.Qt.DashLine})

        print("Resting membrane potential: %0.1f mV\n" % self.rest_vm())
