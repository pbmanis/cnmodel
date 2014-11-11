__author__='pbmanis'

"""
Sync_Freq.py does a measurement of phase locking as a function of frequency for the selected model
Attempts to uses multiprocessing to take advantage of the processors in a given system

"""


import faulthandler
import sys
import numpy as np
from neuron import h
#import pyqtgraph.multiprocess as mproc
import multiprocessing
import pyqtgraph as pg
import pylibrary.pyqtgraphPlotHelpers as pgh
import nrnlibrary.util.pynrnutilities as PU
from nrnlibrary.protocols import Protocol
from nrnlibrary import cells
from nrnlibrary.util import sound
import cPickle
from PyQt4 import QtGui
faulthandler.enable()

def plotResults(results):

    win = pgh.figure(title='RIFunc')
    layout = pgh.LayoutMaker(rows=len(results), cols=2, win=win, labelEdges=True, ticks='talbot')

    for i, result in enumerate(results):
        if result is None:
            continue
        layout.plot((i,0), result['r']['t'], result['r']['vm'])
        layout.title((i,0), '%d dBSPL'%result['pars']['SPL'])
        layout.title((i,1), title='AN spikes at %d dBSPL' % result['pars']['SPL'])
        nconverge = result['pars']['nConverge']
        dist = 1./nconverge
        ytick = np.linspace(0, dist, nconverge+1)
        for k in range(nconverge):
            spiketrain = result['r']['preCell%03d'%k]
            nsp = len(spiketrain)
            vt = pg.ScatterPlotItem(spiketrain, [ytick[k+1]]*nsp, symbol='+', pen=(k,nconverge))
            layout.getPlot((i,1)).addItem(vt)
    pgh.show()

def saveResults(results):
    """
    Store the results to a disk file.
    """
    pass

def retrieveResults():
    pass


class SGCInputTestPL(Protocol):
    def __init__(self):
        self.run_duration = 1.
        self.pip_duration = 0.8
        self.pip_start = [0.02]
        self.Fs = 100e3
        self.f0 = 4000.  # stimulus frequency
        self.cf = self.f0  # fiber cf location (frequency)
        self.dBSPL = 60.
        self.sr = ['l']  # list of sr group per convergent input
        self.nconverge = len(self.sr)
        self.depFlag = False  # depressing or not depressing synapses

    def setSR(self, sr):
        self.sr = sr

    def setCF(self, cf):
        self.cf = cf
        self.f0 = cf

    def setdB(self, dB):
        self.dBSPL = dB

    def info(self):
        return{'Fs': self.Fs, 'F0': self.f0, 'CF': self.cf, 'SPL': self.dBSPL, 'sr': self.sr, 'nconverge': self.nconverge}

    def getPars(self):
        return({'CF': self.cf, 'dB': self.dBSPL, 'temp': 34.0,
                            'dt': 0.025, 'seed': 5759820, 'SRs': self.sr})

    def run_one(self, pars={'CF': 4000., 'dB': 60, 'temp': 34.0,
                            'dt': 0.025, 'seed': 5759820, 'SRs': {'l': 1, 'm': 2, 'h': 1}}):
        self.setCF(pars['CF'])
        self.setdB(pars['dB'])
        self.setSR(pars['SRs'])
        self.nconverge = len(pars['SRs'])
        srmap = {'h': 1, 'm': 2, 'l': 3}
        temp = pars['temp']
        dt = pars['dt']
        seed = pars['seed']
        postCell = cells.Bushy.create()
        preCell = [None]*self.nconverge
        synapse = [None]*self.nconverge
        for i in range(self.nconverge):
            print 'making sgc of type: ', self.sr[i], srmap[self.sr[i]]
            preCell[i] = cells.DummySGC(cf=self.cf, sr=srmap[self.sr[i]])
            synapse[i] = preCell[i].connect(postCell)
        self.pre_cell = preCell
        self.post_cell = postCell
        self.pre_spk = [None]*self.nconverge
        self.synapse = synapse
        self.stim = sound.TonePip(rate=self.Fs, duration=self.run_duration, f0=self.f0, dbspl=self.dBSPL,
                                  ramp_duration=2.5e-3, pip_duration=self.pip_duration,
                                  pip_start=self.pip_start)

        for i in range(self.nconverge):
            preCell[i].set_sound_stim(self.stim, seed=seed*i)

        self.res = {'vm': h.Vector()}
        self.res['vm'].record(postCell.soma(0.5)._ref_v)
        #self['prevm'] = preCell.soma(0.5)._ref_v
        for i in range(30):
            self.res['xmtr%d'%i] = h.Vector()
            self.res['xmtr%d'%i].record(synapse[0].terminal.relsite._ref_XMTR[i])
            synapse[0].terminal.relsite.Dep_Flag = self.depFlag
        self.res['t'] = h.Vector()
        self.res['t'].record(h._ref_t)

        h.tstop = 1e3*self.run_duration # duration of a run
        h.celsius = temp
        h.dt = dt

        self.custom_init()
        h.run()
        for i in range(self.nconverge):
            self.pre_spk[i] = self.pre_cell[i]._spiketrain
#        print 'vm: ', self.res['vm']
        self.post_spk = PU.findspikes(np.array(self.res['t']), np.array(self.res['vm']), -30.)
#        print 'bu spk: ', self.post_s


    def get_VS(self, text_ident):
        if text_ident == 'AN':
            return self.computeVS(self.pre_cell[0]._spiketrain, text_ident=text_ident)
        elif text_ident == 'post':
            return self.computeVS(self.post_spk, text_ident=text_ident)

    def computeVS(self, spikes, text_ident='??'):
        if len(spikes) == 0:
            return {'r': 0., 'ph': 0., 'd': 0., 'n': 0, 'p': 1.0, 'R': 0.}
        phasewin = [self.pip_start[0] + 0.25*self.pip_duration, self.pip_start[0] + self.pip_duration]
        spkin = spikes[np.where(spikes > phasewin[0]*1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
        vs = PU.vector_strength(spikesinwin, self.f0)
        print ('{:4s} [cf={:.1f}] Hz Vector Strength: {:7.3f}, d={:.2f} (us) Rayleigh: {:7.3f}  p = {:.3e}  n = {:d}'
              .format(text_ident, self.cf, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n']))
        return vs

    def show(self):
        self.win = pg.GraphicsWindow()
        Fs = self.Fs
        p1 = self.win.addPlot(title='stim', row=0, col=0)
        p1.plot(self.stim.time * 1000, self.stim.sound)
        p1.setXLink(p1)

        p2 = self.win.addPlot(title='AN spikes', row=1, col=0)
        vt = pg.VTickGroup(self.pre_spk[0])
        p2.addItem(vt)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title='Bushy Spikes', row=2, col=0)
        bspk = PU.findspikes(self.res['t'], self.res['vm'], -30.)
        bspktick = pg.VTickGroup(bspk)
        p3.addItem(bspktick)
        p3.setXLink(p1)

        p4 = self.win.addPlot(title='Bushy Vm', row=3, col=0)
        p4.plot(self.res['t'], self.res['vm'])
        p4.setXLink(p1)

        p5 = self.win.addPlot(title='xmtr', row=0, col=1)
        for i in range(30):
            p5.plot(self.res['t'], self.res['xmtr%d'%i], pen=(i, 15))
        p5.setXLink(p1)

        p6 = self.win.addPlot(title='AN phase', row=1, col=1)
        vs = self.computeVS(self.pre_spk[0], 'AN')
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p6.setXRange(0., 2*np.pi)

        p7 = self.win.addPlot(title='Bushy phase', row=2, col=1)
        vs = self.computeVS(bspk, 'BU')
        (hist, binedges) = np.histogram(vs['ph'])
        curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
        p7.setXRange(0., 2*np.pi)
        p7.setXLink(p6)

        self.win.show()


def run_oneparset(pars):
    s = SGCInputTestPL()
    s.setCF(pars['CF'])
    s.setdB(pars['dB'])
    s.setSR(pars['SRs'])
    s.run_one(pars)
    result = {'AN_VS': s.get_VS('AN'),
                                 'post_VS': s.get_VS('post'),
                                 'pars': s.info(),
                                 't': np.array(s.res['t']),
                                 'vm': np.array(s.res['vm']),
                                 'post_spk': s.post_spk,
                                 'pre_spk': s.pre_spk,
                                 'stim': s.stim,
                                  }
    return result

def computeVS(self, spikes, text_ident='??'):
    if len(spikes) == 0:
        return {'r': 0., 'ph': 0., 'd': 0., 'n': 0, 'p': 1.0, 'R': 0.}
    phasewin = [self.pip_start[0] + 0.25*self.pip_duration, self.pip_start[0] + self.pip_duration]
    spkin = spikes[np.where(spikes > phasewin[0]*1e3)]
    spikesinwin = spkin[np.where(spkin <= phasewin[1]*1e3)]
    vs = PU.vector_strength(spikesinwin, self.f0)
    print ('{:4s} [cf={:.1f}] Hz Vector Strength: {:7.3f}, d={:.2f} (us) Rayleigh: {:7.3f}  p = {:.3e}  n = {:d}'
          .format(text_ident, self.cf, vs['r'], vs['d']*1e6, vs['R'], vs['p'], vs['n']))
    return vs

def run_multifreq(srlist):
    #cflist = np.logspace(np.log10(200), np.log10(6000), 5)
    cflist = [200, 333., 500., 1000., 1200., 1500., 2000.,
              2500., 3000.,4000., 5000.]
    nCFs = len(cflist)

    nWorkers = 6
    pars = [None]*nCFs
    for i in range(nCFs):
        pars[i] = {'CF': cflist[i], 'dB': 60., 'i': i, 'temp': 34.0,
                    'dt': 0.025, 'seed': 57598201*i, 'SRs': srlist}
    pool = multiprocessing.Pool(processes = nWorkers)
    results = pool.map(run_oneparset, (pars[j] for j in range(nCFs)))

    #results = [p.get() for p in res]  # from async...
    pool.close()
    pool.join()

    # TASKS = [s for s in range(len(cflist))]
    # results = [None]*len(TASKS)
    # sgc_post = [None]*len(TASKS)
    # with mproc.Parallelize(enumerate(TASKS), results=results, workers=nWorkers) as tasker:
    #     for i, x in tasker:
    #         sgc_post[i] = SGCInputTestPL()
    #         print i
    #         sgc_post[i].setCF(cflist[i])
    #         sgc_post[i].run_one(seed=57598*i)
    #         tasker.results[i] = {'AN_VS': sgc_post[i].get_VS('AN'),
    #                              'post_VS': sgc_post[i].get_VS('post'),
    #                              'pars': sgc_post[i].info(),
    #                              't': np.array(sgc_post[i].res['t']),
    #                              'vm': np.array(sgc_post[i].res['vm']),
    #                              'post_spk': sgc_post[i].post_spk,
    #                              'pre_spk': sgc_post[i].pre_spk,
    #                              }
    return(results)

def show_vs(results, title=''):
    print 'show_vs'
    win = pgh.figure(title=title)
    layout = pgh.LayoutMaker(cols=1,rows=2, win=win, labelEdges=True, ticks='talbot')
    p1 = layout.getPlot(0)
    p2 = layout.getPlot(1)
    layout.title(0, 'VS')
    layout.title(1, 'd')
    #talbotTicks(layout.getPlot(1))
    #layout.columnAutoScale(col=3, axis='left')
    # reassemble results
    fl = [x['pars']['F0']*1e-3 for x in results]

    vs_pre = [x['AN_VS']['r'] for x in results]
    vs_post = [x['post_VS']['r'] for x in results]
    pp = pg.PlotDataItem(fl, vs_pre, pen=pg.mkPen('k'), symbol='o', symbolPen=pg.mkPen('k'), symbolBrush=pg.mkBrush('k'))
    p1.addItem(pp)
    ps = pg.PlotDataItem(fl, vs_post, pen = pg.mkPen('r'), symbol='s', symbolPen=pg.mkPen('r'), symbolBrush=pg.mkBrush('r'))
    p1.addItem(ps)
    p1.setLogMode(x=True)
    pgh.labelAxes(p1, 'F (Hz)', 'Vector Strength')
    d_pre = [x['AN_VS']['d']*1e6 for x in results]
    d_post = [x['post_VS']['d']*1e6 for x in results]
    pp = pg.PlotDataItem(fl, d_pre, pen=pg.mkPen('g'), symbol='o', symbolPen=pg.mkPen('g'), symbolBrush=pg.mkBrush('g'))
    p2.addItem(pp)
    ps = pg.PlotDataItem(fl, d_post, pen = pg.mkPen('r'), symbol='s', symbolPen=pg.mkPen('r'), symbolBrush=pg.mkBrush('r'))
    p2.addItem(ps)
    pgh.labelAxes(p2, 'F (Hz)', 'Dispersion (uS)')
    p2.setLogMode(x=True)
    pgh.show()
    return win

def show_traces(result):

    win = pg.GraphicsWindow()
    Fs = result['pars']['Fs']
    p1 = win.addPlot(title='stim', row=0, col=0)

    p1.plot(result['stim'].time * 1000, result['stim'].sound)
    p1.setXLink(p1)

    p2 = win.addPlot(title='AN spikes', row=1, col=0)
    vt = pg.VTickGroup(result['pre_spk'][0])
    p2.addItem(vt)
    p2.setXLink(p1)

    p3 = win.addPlot(title='Bushy Spikes', row=2, col=0)
    bspk = PU.findspikes(result['t'], result['vm'], -30.)
    bspktick = pg.VTickGroup(bspk)
    p3.addItem(bspktick)
    p3.setXLink(p1)

    p4 = win.addPlot(title='Bushy Vm', row=3, col=0)
    p4.plot(result['t'], result['vm'])
    p4.setXLink(p1)

    # p5 = self.win.addPlot(title='xmtr', row=0, col=1)
    # for i in range(30):
    #     p5.plot(result['t'], self.res['xmtr%d'%i], pen=(i, 15))
    # p5.setXLink(p1)

    p6 = win.addPlot(title='AN phase', row=1, col=1)
    vs = computeVS(result['pre_spk'][0], 'AN')
    (hist, binedges) = np.histogram(vs['ph'])
    curve = p6.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
    p6.setXRange(0., 2*np.pi)

    p7 = win.addPlot(title='Bushy phase', row=2, col=1)
    vs = computeVS(bspk, 'BU')
    (hist, binedges) = np.histogram(vs['ph'])
    curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
    p7.setXRange(0., 2*np.pi)
    p7.setXLink(p6)

    win.show()


if __name__ == '__main__':
    single = False
    showdata = True
    dorun = False

    species = 'guineapig'
    celltype = 'bushy'
    srlist = ['l', 'm', 'm', 'h']

    convergence_pattern='LS%03dMS%03dHS%03d' % (srlist.count('l'), srlist.count('m'), srlist.count('h'))
    filename = 'Phase-Runs_%s_%s_%s.p' % (species, celltype, convergence_pattern )

    if single:
        u = SGCInputTestPL()
        u.setSR(srlist)
        u.run_one(pars=u.getPars())
        u.show()

    if dorun:
        results = run_multifreq(srlist)
        rf = open(filename, 'w')
        cPickle.dump(results, rf)
        rf.close()

    if showdata:
        rf = open(filename, 'r')
        results = cPickle.load(rf)
        rf.close()
        w=show_vs(results, filename)
        show_traces(results[0])






