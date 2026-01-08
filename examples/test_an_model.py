
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
import argparse
import sys
import time
import numpy as np
import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound
from pyzbc2014 import pyzbc2014

zbc = pyzbc2014.pyzbc2014()

def test_an_model(fibertype:str='hsr', cf:float=16e3, dbspl:float=30., noiseType=0, species=1   ):
    assert fibertype in ['hsr', 'msr', 'lsr'], "Invalid fiber type: %s" % fibertype
    assert species in ['cat', 'human', 'human_glasberg'], "Invalid species: %s" % species


    # model fiber parameters
    CF    = cf # 16e3   # CF in Hz   
    cohc  = 1.0    # normal ohc function
    cihc  = 1.0    # normal ihc function
    # species = species    # 1 for cat (2 for human with Shera et al. tuning 3 for human with Glasberg & Moore tuning)
    noiseType = 0  # 1 for variable fGn (0 for fixed fGn)
    fiberType = fibertype  # spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects "1" = Low "2" = Medium "3" = High
    implnt = 0     # "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

    if noiseType == 0:
        noisetype = 'none'
    elif noiseType == 1:
        noisetype = 'fresh'
    else:
        raise ValueError("Invalid noiseType: %d" % noiseType)

    F0 = CF     # stimulus frequency in Hz
    Fs = 100e3  # sampling rate in Hz (must be 100, 200 or 500 kHz)
    pdur = 0.3  # pip duration
    pstart = [0.1]  # pip start times
    T = pstart[-1] + pdur + 0.1  # total stimulus duration
    rt = 2.5e-3 # rise/fall time in seconds
    stimdb = dbspl # stimulus intensity in dB SPL

    # PSTH parameters
    nrep = 1           # number of stimulus repetitions (e.g., 50)
    psthbinwidth = 0.5e-3 # binwidth in seconds

    stim = sound.TonePip(rate=Fs, duration=T, f0=F0, dbspl=stimdb, 
                         pip_duration=pdur, pip_start=pstart, ramp_duration=rt)
    t = stim.time
    pin = stim.sound
    db = stim.measure_dbspl(pstart[-1]+rt, pdur-rt) # only measure during steady state portion of pip


    an_model.seed_rng(34978)
    start = time.time()
    # vihc = an_model.model_ihc(pin, CF, nrep, 1/Fs, T+1e-3, cohc, cihc, species) # , _transfer=False) 
    vihc = zbc.sim_ihc_zbc2014(pin, CF, nrep, Fs, cohc, cihc, species=species)
    print("IHC time:", time.time() - start)
    start = time.time()
#    m, v, psth = an_model.model_synapse(vihc, CF, nrep, 1/Fs, fiberType, noiseType, implnt)
    powerlaw = 'true'
    noiseType = 'fresh'

    an_drive = zbc.sim_anrate_zbc2014(vihc, CF, nrep, Fs, fiberType, powerlaw, noiseType)
    print("Syn time:", time.time() - start)

    win = pg.GraphicsLayoutWidget()
    win.setWindowTitle('AN Model Test')
    win.resize(800, 800)
    p1 = win.addPlot(title=f'Species: {species} Freq: {CF} Hz Level {db:0.1f} dBSPL, Fiber: {fibertype}')
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
            ihcout = zbc.sim_ihc_zbc2014(pin, cf=CF, nrep=1, species=species)
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
    parser = argparse.ArgumentParser(description="test an model basics")

    parser.add_argument(
        "--species",
        type=str,
        choices=["cat", "human", "human_glasberg"],
        default="cat",
        help="Species",
    )
    parser.add_argument(
        "-S",
        "--stimulus",
        type=str,
        choices=["tone", "noise", "SAM", "clicks"],
        default="tone",
        help="stimulus type",
    )
    parser.add_argument(
        "--dB",
        "--dBSPL",
        type=float,
        default=30.,
        help="Sound pressure level, SPL",
    )
    parser.add_argument(
        "-f",
        "--fibertype",
        type=str,
        choices=["hsr", "msr", "lsr"],
        default="hsr",
        help="Fiber type (spontaneous rate) hsr, msr, lsr",
    )
    parser.add_argument(
        "--dmod",
        type=float,
        default=100.,
        help="Modulation depth for SAM (percent)",
    )
    parser.add_argument(
        "--fmod",
        type=float,
        default=200.0,
        help="Modulation Frequency for SAM (Hz)",
    )
    parser.add_argument(
        "--CF",
        type=float,
        default=16000.,
        help="Carrier Frequency for SAM (Hz)",
    )
    
    parser.add_argument(
        "--RI",
        action="store_true",
        default=False,
        dest="RI",
        help="Run Rate-intensity with these parameters",
    )

    args = parser.parse_args()
    fibertype = args.fibertype

    test_an_model(fibertype=args.fibertype, noiseType=1, species=args.species, cf=args.CF, dbspl=args.dB)
