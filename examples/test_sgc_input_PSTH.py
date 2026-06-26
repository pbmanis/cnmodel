import sys
import argparse
import numpy as np
import timeit
import pyqtgraph as pg
import pyqtgraph.multiprocess as mp

# import dask
# import dask.delayed

from neuron import h
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
from cnmodel.util.pynrnutilities import custom_init
import cnmodel.util.pynrnutilities as PU
from cnmodel import data

species = "mouse"  # tables for other species do not yet exist



def buildCell(cell):
    """
    Make the postsynaptic cell
    
    Parameters
    ----------
    cell : str
        Name of a cell type to instantiate
    
    Return
    ------
    cnmodel postsynaptic cell object
    """
    if cell == "bushy":
        post_cell = cells.Bushy.create(species=species)
    elif cell == "tstellate":
        post_cell = cells.TStellate.create(species=species)
    elif cell == "octopus":
        post_cell = cells.Octopus.create(species=species)
    elif cell == "dstellate":
        post_cell = cells.DStellate.create(species=species)
    elif cell == "tuberculoventral":
        post_cell = cells.DStellate.create(species=species)
    elif cell == "pyramidal":
        post_cell = cells.DStellate.create(species=species)
    else:
        raise ValueError("cell %s is not yet implemented for PSTH testing" % self.cell)
    # print('made postcell')
    return post_cell

def makeSynapses(post_cell, info):
    """
    Attach synapses to the postsynaptic cell by the rules in info
    
    Parameters
    post_cell : cnmodel cell object
    
    info: dict 
    
    Returns
    -------
    dict of:
        pre_cells (list of cnmodel dummy SGC cells)
        synapses (list of snmodel synapses) 
        and the xmtr variable (hoc)
    """
    pre_cells = []
    synapses = []
    j = 0
    xmtr = {}

    for _ in range(info["n_sgc"]):
        pre_cells.append(cells.DummySGC(cf=info["cf"], sr=info["sr"]))
        if info["synapse_type"] == "simple":
            synapses.append(pre_cells[-1].connect(post_cell, type=info["synapse_type"]))
            synapses[-1].terminal.netcon.weight[0] = info["gmax"]
        elif info["synapse_type"] == "multisite":
            synapses.append(
                pre_cells[-1].connect(
                    post_cell,
                    post_opts={"AMPAScale": 1.0, "NMDAScale": 1.0},
                    type=info["synapse_type"],
                )
            )
            for i in range(synapses[-1].terminal.n_rzones):
                xmtr["xmtr%04d" % j] = h.Vector()
                xmtr["xmtr%04d" % j].record(synapses[-1].terminal.relsite._ref_XMTR[i])
                j = j + 1
            synapses[
                -1
            ].terminal.relsite.Dep_Flag = False  # no depression in these simulations
        # Claude fixed 2026-06-26: set_sound_stim moved to setStimuli(); called per rep
        # so the same pre_cells can be reused across repetitions without recreating cells.

    return {'pre_cells':pre_cells, 'synapses': synapses, 'xmtr':xmtr}


def setStimuli(psx, info):
    """Update the AN spike train seed on existing pre_cells for one repetition.

    Call this once per repetition (after updating info["seed"]) instead of
    recreating pre_cells and synapses from scratch each time.
    """
    for nsgc, pc in enumerate(psx['pre_cells']):
        pc.set_sound_stim(
            info["stim"], seed=info["seed"] + nsgc, simulator=info["simulator"]
        )


def runNeuron(post_cell, info, psx):

    Vm = h.Vector()
    Vm.record(post_cell.soma(0.5)._ref_v)
    rtime = h.Vector()
    rtime.record(h._ref_t)
    tstop = 1e3 * info["run_duration"]  # duration of a run
    h.celsius = info["temp"]
    h.dt = info["dt"]
    post_cell.cell_initialize()
    info["init"]()
    # Claude changed 2026-06-26: batch_run (no filename, no prior batch_save) leaves stale
    # NetCon events in the queue in NEURON 9; the next rep's cell_initialize→finitialize hits
    # "events out of order". fadvance loop avoids batch-mode side effects entirely.
    # custom_init already sets h.t=0, so no explicit reset needed.
    # TO REVERT: uncomment h.t/h.batch_run below and remove the while loop.
    # h.t = 0.0
    # h.batch_run(tstop, h.dt)
    while h.t < tstop:
        h.fadvance()
    pre_trains = [psx['pre_cells'][k]._spiketrain for k in range(len(psx['pre_cells']))]
    # print(psx['xmtr'].keys())
    # print([psx['xmtr'][k].to_python() for k in list(psx['xmtr'].keys())])
    # exit()
    return {
        "time": np.array(rtime),
        "vm": Vm.to_python(),
        "xmtr": [psx['xmtr'][k].to_python() for k in list(psx['xmtr'].keys())],
        "pre_cells": pre_trains,
        # "post_cell": post_cell,
        # "synapses": psx['synapses'],
    }


def runTrial(cell, info):
    # Used by parallel workers — each worker handles one rep so no accumulation occurs.
    post_cell = buildCell(cell)
    psx = makeSynapses(post_cell, info)
    setStimuli(psx, info)  # Claude fixed 2026-06-26: set stim separately from synapse setup
    res = runNeuron(post_cell, info, psx)
    return res
    
    
class SGCInputTestPSTH(Protocol):
    def __init__(self):
        self.nrep = 10
        self.stimulus = "tone"
        self.run_duration = 0.20  # in seconds
        self.pip_duration = 0.05  # in seconds
        self.pip_start = [0.1]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = 4000.0  # stimulus in Hz
        self.cf = 4000.0  # SGCs in Hz
        self.fMod = 100.0  # mod freq, Hz
        self.dMod = 0.0  # % mod depth, Hz
        self.dbspl = 50.0
        self.simulator = "py3" # "cochlea"
        self.sr = 2  # set SR group
        self.synapse_type = "simple"  # Claude fixed 2026-06-26: was module-level global

    def set_cell(self, cell="bushy"):
        self.cell = cell

    def set_freq(self, fr=4000.0):
        self.f0 = fr
    
    def set_cf(self, cf=4000.0):
        self.cf = fr
    
    def set_db(self, db=50.):
        self.dbspl = db
    
    def set_sr(self, srgroup=1):
        self.sr = srgroup

    def set_synapse_type(self, syntype="simple"):
        self.synapse_type = syntype
        
    def run(
        self,
        temp=34.0,
        dt=0.025,
        seed=575982035,
        reps=10,
        stimulus="tone",
        # Claude fixed 2026-06-26: was "cochlea"; default matches __init__ and CLI default
        simulator="py3",
        # Claude fixed 2026-06-26: default changed from 'serial' to 'mp'
        parallelmode='mp',
    ):
        self.simulator = simulator  # Claude fixed 2026-06-26: was never assigned, so info dict always used __init__ value
        self.stimulus = stimulus
        assert self.stimulus in ["tone", "SAM", "noise", "clicks"]  # cases available
        assert self.cell in [
            "bushy",
            "tstellate",
            "octopus",
            "dstellate",
            "tuberculoventral",
            "pyramidal",
        ]

        self.nrep = reps
        if self.stimulus == "SAM":
            self.dMod = 100.0
            self.stim = sound.SAMTone(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                fmod=self.fMod,
                dmod=self.dMod,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
            )
        if self.stimulus == "tone":
            self.f0 = 4000.0
            self.cf = 4000.0
            self.stim = sound.TonePip(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
            )
        if self.stimulus == "noise": # will be frozen, as only generated once.
            self.f0 = 4000.0
            self.cf = 4000.0
            self.stim = sound.NoisePip(
                rate=self.Fs,
                duration=self.run_duration,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
                seed=9,
            )

        if self.stimulus == "clicks":
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                dbspl=self.dbspl,
                click_start=0.010,
                click_duration=100.0e-6,
                click_interval=self.click_rate,
                nclicks=int((self.run_duration - 0.01) / self.click_rate),
                ramp_duration=2.5e-3,
            )

        n_sgc = data.get(
            "convergence", species=species, post_type=self.cell, pre_type="sgc"
        )[0]
        self.n_sgc = int(np.round(n_sgc))
        # for simple synapses, need this value:
        self.AMPA_gmax = (
            data.get(
                "sgc_synapse", species=species, post_type=self.cell, field="AMPA_gmax"
            )[0]
            / 1e3
        )  # convert nS to uS for NEURON
        self.vms = [None for n in range(self.nrep)]
        self.synapses = [None for n in range(self.nrep)]
        self.xmtrs = [None for n in range(self.nrep)]
        self.pre_cells = [None for n in range(self.nrep)]
        self.time = [None for n in range(self.nrep)]
        info = {
            "n_sgc": self.n_sgc,
            "gmax": self.AMPA_gmax,
            "stim": self.stim,
            "simulator": self.simulator,
            "cf": self.cf,
            "sr": self.sr,
            "seed": seed,
            "run_duration": self.run_duration,
            "temp": temp,
            "dt": dt,
            "init": custom_init,
            "synapse_type": self.synapse_type,
        }
        if parallelmode == 'serial':
            start_time = timeit.default_timer()
            # Claude fixed 2026-06-26: build cell and synapses ONCE outside the rep loop.
            # Creating new NEURON sections per rep accumulated them in NEURON's global
            # section list, causing a Bus error (~rep 12) when finitialize iterated
            # over partially-freed sections from previous reps.
            post_cell = buildCell(self.cell)
            psx = makeSynapses(post_cell, info)
            for nr in range(self.nrep):
                info["seed"] = seed + 3 * self.n_sgc * nr
                info["rep_no"] = nr
                print(f"\r  rep {nr+1}/{self.nrep}", end="", flush=True)
                setStimuli(psx, info)
                res = runNeuron(post_cell, info, psx)
                self.pre_cells[nr] = res["pre_cells"]
                self.time[nr] = res["time"]
                self.xmtrs[nr] = res['xmtr']
                self.vms[nr] = res["vm"]
            print()  # close the \r progress line
            elapsed = timeit.default_timer() - start_time
            print(f"Not Parallel Elapsed time for {self.nrep:d} stimuli: {elapsed:f} secs")

        if parallelmode in ['mp', 'multiprocessing']:
            results = {}
            workers = mp.Parallelize.suggestedWorkerCount()  # use suggested #
            tasks = list(range(self.nrep))
            start_time = timeit.default_timer()
            # Claude fixed 2026-06-26: create cell once per worker (same root cause as
            # serial-mode fix). Each worker handles several tasks sequentially; creating
            # new NEURON sections per task accumulates them and causes Bus error.
            # Mutable sentinel: forked worker creates its cell on the first task and
            # reuses it for all subsequent tasks assigned to that worker.
            _cell_ref: list = [None]
            _psx_ref: list = [None]
            with mp.Parallelize(enumerate(tasks), results=results, workers=workers, progressDialog='Running repetitions..') as tasker:
                for i, task in tasker:
                    if _cell_ref[0] is None:
                        _cell_ref[0] = buildCell(self.cell)
                        _psx_ref[0] = makeSynapses(_cell_ref[0], info)
                    info["seed"] = seed + 3 * self.n_sgc * i
                    info["rep_no"] = i
                    setStimuli(_psx_ref[0], info)
                    repres = runNeuron(_cell_ref[0], info, _psx_ref[0])
                    tasker.results[i] = repres
            # get time of run before display
            elapsed = timeit.default_timer() - start_time
            print('Elapsed time for %d stimuli: %f  (%f sec per stim)' % (len(tasks), elapsed, elapsed/len(tasks)))
            res = results
            for nr in range(self.nrep):
                self.pre_cells[nr] = res[nr]["pre_cells"]
                self.time[nr] = res[nr]["time"]
                self.xmtrs[nr] = res[nr]["xmtr"] # {k: v.to_python() for k, v in list(res[nr]["xmtr"].items())}
                self.vms[nr] = res[nr]["vm"]
                # self.synapses[nr] = res[nr]["synapses"]
                
        # if parallelmode == 'dask':
        #     # this does not work...
        #     # dask does not allow neuron section objects cannot be pickled to be passed between
        #     # elements of the task.
        #     # import dask.multiprocessing
        #     # dask.config.set(scheduler='processes')
        #     from dask.distributed import Client, LocalCluster
        #     cluster = LocalCluster()
        #     client = Client(cluster) #client = Client() # Client(cluster) # client = Client('127.0.0.1:8786')
        #     start_time = timeit.default_timer()
        #     def _tasks(nrep, nsgc, info, cell):
        #         res = []
        #         for nr in range(nrep):
        #             info["seed"] = seed + 3 * nsgc * nr
        #             # repres = runTrial_dask(self.cell, info)
        #             post_cell = buildCell_dask(cell)
        #             psx = makeSynapses_dask(post_cell, info)
        #             # nvars = setupNeuron_dask(post_cell, info)
        #             repres = runNeuron_dask(post_cell, info, psx) # , nvars['rtime'], nvars['Vm'])
        #             res.append(repres)
        #         return(res)
        #
        #     res = dask.compute(_tasks(self.nrep, self.n_sgc, info, self.cell))# , scheduler='single-threaded')
        #     elapsed = timeit.default_timer() - start_time
        #     print(f"Parallel Elapsed time for {self.nrep:d} stimuli: {elapsed:f} secs")
        #     # res = results
        #     res = res[0]
        #     for nr in range(self.nrep):
        #         # print('nr: ', nr)
        #         # print(res[nr])
        #         self.pre_cells[nr] = res[nr]["pre_cells"]
        #         self.time[nr] = res[nr]["time"]
        #         self.xmtrs[nr] = res[nr]["xmtr"] # {k: v.to_python() for k, v in list(res[nr]["xmtr"].items())}
        #         self.vms[nr] = res[nr]["vm"]
        #         # self.synapses[nr] = res[nr]["synapses"]
        #


    def show(self):
        self.win = pg.GraphicsLayoutWidget()
        self.win.setBackground("w")
        Fs = self.Fs
        p1 = self.win.addPlot(
            title="Stimulus", row=0, col=0, labels={"bottom": "T (ms)", "left": "V"}
        )
        p1.plot(self.stim.time * 1000, self.stim.sound, pen=pg.mkPen("k", width=0.75))
        p1.setXLink(p1)

        p2 = self.win.addPlot(
            title="AN spikes",
            row=1,
            col=0,
            labels={"bottom": "T (ms)", "left": "AN spikes (first trial)"},
        )
        for nr in range(self.nrep):
            xan = []
            yan = []
            for k in range(len(self.pre_cells[nr])):
                r = self.pre_cells[nr][k] # [k]._spiketrain
                xan.extend(r)
                yr = k + np.zeros_like(r) + 0.2
                yan.extend(yr)
            c = pg.PlotCurveItem()
            xp = np.repeat(np.array(xan), 2)
            yp = np.repeat(np.array(yan), 2)
            yp[1::2] = yp[::2] + 0.6
            c.setData(
                xp.flatten(),
                yp.flatten(),
                connect="pairs",
                pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0),
            )
            p2.addItem(c)
        p2.setXLink(p1)

        p3 = self.win.addPlot(
            title="%s Spikes" % self.cell,
            row=2,
            col=0,
            labels={"bottom": "T (ms)", "left": "Trial #"},
        )
        xcn = []
        ycn = []
        xspks = []
        for k in range(self.nrep):
            bspk = PU.findspikes(self.time[k], self.vms[k], -35.0)
            xcn.extend(bspk)
            yr = k + np.zeros_like(bspk) + 0.2
            ycn.extend(yr)
        d = pg.PlotCurveItem()
        xp = np.repeat(np.array(xcn), 2)
        yp = np.repeat(np.array(ycn), 2)
        yp[1::2] = yp[::2] + 0.6
        d.setData(
            xp.flatten(), yp.flatten(), connect="pairs", pen=pg.mkPen("k", width=1.5)
        )
        p3.addItem(d)
        p3.setXLink(p1)

        p4 = self.win.addPlot(
            title="%s Vm" % self.cell,
            row=3,
            col=0,
            labels={"bottom": "T (ms)", "left": "Vm (mV)"},
        )
        for nr in range(self.nrep):
            p4.plot(
                self.time[nr],
                self.vms[nr],
                pen=pg.mkPen(pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0),
            )
        p4.setXLink(p1)

        p5 = self.win.addPlot(
            title="xmtr", row=0, col=1, labels={"bottom": "T (ms)", "left": "gSyn"}
        )
        if self.synapse_type == "multisite":
            # for nr in [0]:
            #     syn = self.synapses[nr]
            #     if syn is None:
            #         continue
            for nr in range(self.nrep):
                # for k in range(self.n_sgc):
                    # print(len(self.xmtrs[k]))
                    # synapse = syn[k]
                    for j in range(len(self.xmtrs[nr])):
                        if j > 10:  # just show a few events... 
                            continue
                        p5.plot(
                            self.time[nr],
                            self.xmtrs[nr][j], # ["xmtr%04d" % j],
                            pen=pg.mkPen(
                                pg.intColor(nr, self.nrep), hues=self.nrep, width=1.0
                            ),
                        )
        p5.setXLink(p1)

        p6 = self.win.addPlot(
            title="AN PSTH",
            row=1,
            col=1,
            labels={"bottom": "T (ms)", "left": "Sp/ms/trial"},
        )
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xan, bins)
        curve6 = p6.plot(
            binedges,
            hist,
            stepMode=True,
            fillBrush=(0, 0, 0, 255),
            brush=pg.mkBrush("k"),
            fillLevel=0,
        )

        p7 = self.win.addPlot(
            title="%s PSTH" % self.cell,
            row=2,
            col=1,
            labels={"bottom": "T (ms)", "left": "Sp/ms/trial"},
        )
        bins = np.arange(0, 200, 1)
        (hist, binedges) = np.histogram(xcn, bins)
        curve7 = p7.plot(
            binedges,
            hist,
            stepMode=True,
            fillBrush=(0, 0, 0, 255),
            brush=pg.mkBrush("k"),
            fillLevel=0,
        )

        self.win.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute AN only PSTH in postsynaptic cell"
    )
    parser.add_argument(
        type=str,
        dest="cell",
        default="bushy",
        choices=[
            "bushy",
            "tstellate",
            "dstellate",
            "octopus",
            "tuberculoventral",
            "pyramidal",
        ],
        help="Select target cell",
    )
    parser.add_argument(
        "-p",
        "--parallel",
        type=str,
        dest="parallelmode",
        # Claude fixed 2026-06-26: default changed from 'serial' to 'mp'
        default="mp",
        choices=["serial", "mp", "multiprocessing"],
        help="Select stimulus from ['serial', 'mp' or 'multiprocessing']",
    )
    
    parser.add_argument(
        "-s",
        "--stimulus",
        type=str,
        dest="stimulus",
        default="tone",
        choices=["tone", "SAM", "noise", "clicks"],
        help="Select stimulus from ['tone', 'noise', 'SAM', 'clicks']",
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        dest="syntype",
        default="simple",
        choices=["simple", "multisite"],
        help="Set synapse type (simple, multisite)",
    )
    parser.add_argument(
        "-n",
        "--nrep",
        type=int,
        dest="nrep",
        default=10,
        help="Set number of repetitions",
    )

    parser.add_argument(
        "-d",
        "--dB",
        type=float,
        dest="dbspl",
        default=50.,
        help="Set sound intensity, dbSPL",
    )
    parser.add_argument(
        "--sr",
        type=int,
        dest="srgroup",
        default=2,
        help="SR group (1-high, 2-medium, 3-low)",
    )
    # Claude fixed 2026-06-26: added --simulator and --temp; script was missing both
    parser.add_argument(
        "--simulator",
        type=str,
        dest="simulator",
        default="py3",
        choices=["py3", "matlab", "cochlea"],
        help="AN model simulator",
    )
    parser.add_argument(
        "--temp",
        type=float,
        dest="temp",
        default=34.0,
        help="Simulation temperature (°C)",
    )

    args = parser.parse_args()

    cell = args.cell

    print("cell type: ", cell)
    print('nrep: ', args.nrep)
    print('parallelmode: ', args.parallelmode)
    prot = SGCInputTestPSTH()
    prot.set_cell(cell)
    prot.set_synapse_type(args.syntype)  # Claude fixed 2026-06-26: was module-level global
    prot.set_db(args.dbspl)
    prot.set_sr(args.srgroup)
    pg.mkQApp()  # Claude fixed 2026-06-26: must exist before mp.Parallelize creates ProgressDialog
    # Claude fixed 2026-06-26: use args.simulator and args.temp instead of hardcoded values
    prot.run(stimulus=args.stimulus, reps=args.nrep, simulator=args.simulator, temp=args.temp, parallelmode=args.parallelmode)
    prot.show()

    import sys

    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()
