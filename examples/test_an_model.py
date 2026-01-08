
"""
Test the embedded auditory nerve model with a set of tone pips.
(Zilany et al. 2014; requires MATLAB)

This demonstrates the lowest-level access to the auditory nerve model that is
available in the cnmodel API. For higher-level tools, see test_sound_stim.py
(which uses an_model.get_spiketrain) and test_sgc_input.py (which uses 
cells.DummySGC).

Usage:
python test_an_model.py 
(no arguments)

(Adapted from makeANF_CF_RI.m)
"""
import time
import numpy as np
import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound
from pyzbc2014 import pyzbc2014

zbc = pyzbc2014.pyzbc2014()

def test_an_model():
    # model fiber parameters
    CF    = 16e3   # CF in Hz   
    cohc  = 1.0    # normal ohc function
    cihc  = 1.0    # normal ihc function
    species = 1    # 1 for cat (2 for human with Shera et al. tuning 3 for human with Glasberg & Moore tuning)
    noiseType = 0  # 1 for variable fGn (0 for fixed fGn)
    fiberType = 3  # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects "1" = Low "2" = Medium "3" = High
    implnt = 0     # "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
    if fiberType == 1:
        sr = 'lsr'
    elif fiberType == 2:
        sr = 'msr'
    elif fiberType == 3:
        sr = 'hsr'
    else:
        raise ValueError("Invalid fiberType: %d" % fiberType)
    if noiseType == 0:
        noisetype = 'none'
    elif noiseType == 1:
        noisetype = 'fresh'
    else:
        raise ValueError("Invalid noiseType: %d" % noiseType)
    if species == 1:
        species_str = 'cat'
    elif species == 2:
        species_str = 'human'
    elif species == 3:
        species_str = 'human_glasberg'
    else:
        raise ValueError("Invalid species: %d" % species)
    # stimulus parameters
    F0 = CF     # stimulus frequency in Hz
    Fs = 100e3  # sampling rate in Hz (must be 100, 200 or 500 kHz)
    # T  = 150e-3  # stimulus duration in seconds
    pdur = 0.3  # pip duration
    pstart = [0.1]  # pip start times
    T = pstart[-1] + pdur + 0.1  # total stimulus duration
    rt = 2.5e-3 # rise/fall time in seconds
    stimdb = 30 # stimulus intensity in dB SPL

    # PSTH parameters
    nrep = 1           # number of stimulus repetitions (e.g., 50)
    psthbinwidth = 0.5e-3 # binwidth in seconds

    stim = sound.TonePip(rate=Fs, duration=T, f0=F0, dbspl=stimdb, 
                         pip_duration=pdur, pip_start=pstart, ramp_duration=rt)
    t = stim.time
    pin = stim.sound
    db = stim.measure_dbspl(rt, T-rt)


    an_model.seed_rng(34978)
    start = time.time()
    # vihc = an_model.model_ihc(pin, CF, nrep, 1/Fs, T+1e-3, cohc, cihc, species) # , _transfer=False) 
    vihc = zbc.sim_ihc_zbc2014(pin, CF, nrep, Fs, cohc, cihc, species=species_str)
    print("IHC time:", time.time() - start)
    start = time.time()
#    m, v, psth = an_model.model_synapse(vihc, CF, nrep, 1/Fs, fiberType, noiseType, implnt)
    fiberType = 'hsr'
    powerlaw = 'true'
    noiseType = 'fresh'

    an_drive = zbc.sim_anrate_zbc2014(vihc, CF, nrep, Fs, fiberType, powerlaw, noiseType)
    print("Syn time:", time.time() - start)

    win = pg.GraphicsLayoutWidget()
    win.setWindowTitle('AN Model Test')
    win.resize(800, 800)
    p1 = win.addPlot(title='Input Stimulus (%0.1f dBSPL)' % db)
    p1.plot(t, pin)

    p2 = win.addPlot(col=0, row=1, title='IHC voltage')
    p2.setXLink(p1)
    #vihc = vihc.get()[0]
    vihc = vihc[:len(vihc) // nrep]
    t = np.arange(len(vihc)) * 1./Fs
    p2.plot(t, vihc)

    p3 = win.addPlot(col=0, row=2, title='ANF drive')
    p3.setXLink(p2)
    ds = 1
    size = an_drive.size // ds
    an_drive = an_drive[:size*ds].reshape(size, ds).sum(axis=1)
    t = np.arange(len(an_drive)) * 1./Fs * ds
    p3.plot(t, an_drive[:-1], stepMode=True, fillLevel=0, fillBrush='w')
    
    p4 = win.addPlot(col=0, row=3, title='PSTH')
    p4.setXLink(p3)
    nreps = 500
    spikes = []
    isi_times = []
    total_stim = t.max() * nreps

    for nr in range(nreps):
        if nr == 0:
            ihcout = zbc.sim_ihc_zbc2014(pin, cf=CF, nrep=1, species=species_str)
            anout = zbc.sim_anrate_zbc2014(
                ihcout, cf=CF, fibertype=fiberType, nrep=1, noisetype=noisetype
            )
        spiketimes = zbc.sim_spike_generator_zbc2014(anout, fs=Fs, totalstim=total_stim, nrep=1)
        spiketimes = spiketimes[spiketimes > 0.0]
        spikes.append(spiketimes)
        stim_spikes = spiketimes[(spiketimes > pstart[0]) & (spiketimes < pstart[0] + pdur)]
        isi_times.append(np.diff(stim_spikes))
    spks = np.concatenate(spikes)
    isi_times = np.concatenate(isi_times)

    binwidth = 0.0005  # 0.5 ms
    numbins = int(np.ceil(t.max() / binwidth))
    rounded_max = np.ceil(t.max())
    xbins = np.arange(0, rounded_max + binwidth, binwidth)
    (hist, binedges) = np.histogram(spks, xbins)
    p4.plot(binedges, hist, stepMode=True, fillLevel=0, fillBrush='w', brush=pg.mkBrush('w'))
    
    p5 = win.addPlot(col=0, row=4, title='ISI Histogram')
    binwidth = 0.0001  # 0.1 ms
    max_isi = 0.05  # 50 ms
    xbins = np.arange(0, max_isi + binwidth, binwidth)
    (hist, binedges) = np.histogram(isi_times, xbins)
    p5.plot(binedges, hist, stepMode=True, fillLevel=0, fillBrush='w', brush=pg.mkBrush('w'))
    
    
    win.show()
    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()
    
if __name__ == "__main__":
    import sys
    test_an_model()
