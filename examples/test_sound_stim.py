"""
Test using sound stimulation to generate SGC spike trains.

This script uses an_model.get_spiketrain(), which internally calls MATLAB to 
generate spike trains and caches the output. A higher-level approach is to use
DummySGC, which will automatically feed the spike train into a synapse for 
input to downstream neurons (see test_sgc_input.py). Lower-level access to the
auditory nerve model is demonstrated in test_an_model.py.

This script also can measure the time required to generate the spike trains, or
to retrieve them from the cache. In addition, it can access a second interface
to the Zilany et al. mode that is in pure python ("cochlea") for comparison.
In general the speed to compute the spike trains for the same sets of stimuli
is faster with the pure Python interface than with the Matlab interface, and
fastest for retrieval of pre-computed trains. Note that changing the value
of "seed" will force a recomputation of the spike trains.
"""

import argparse
import sys
import numpy as np
import time

import pyqtgraph as pg
from cnmodel import an_model
from cnmodel.util import sound
try:
    import cochlea
    HAVE_COCHLEA = True
except ImportError:
    HAVE_COCHLEA = False

seed = 34986


def time_usage(func):
    def wrapper(*args, **kwargs):
        beg_ts = time.time()
        res = func(*args, **kwargs)
        end_ts = time.time()
        print(("** Elapsed time: %f" % (end_ts - beg_ts)))
        return res

    return wrapper


def set_dbspl(signal, dbspl):
    """Scale the level of `signal` to the given dB_SPL."""
    p0 = 20e-6
    rms = np.sqrt(np.sum(signal ** 2) / signal.size)
    scaled = signal * 10 ** (dbspl / 20.0) * p0 / rms
    return scaled


@time_usage
def sound_stim(seed, reps=10, pip_duration:float = 0.5, useMatlab=False):
    cf = 1.5e3
    levels = list(range(-10, 101, 10))

    result = {}
    spikes: dict = {}
    if useMatlab:
        simulator = "matlab"
    else:
        simulator = "py3"
    for sr in 0, 1, 2:
        spikes = {}
        for rep in range(reps):
            spikes[rep] = []
            for level in levels:
                stim = sound.TonePip(
                        rate=100e3,
                        duration=pip_duration+0.005,
                        f0=cf,
                        dbspl=level,
                        pip_duration=pip_duration,
                        pip_start=[0.005],
                        ramp_duration=2.5e-3,
                    )
                if simulator == "py3":
                        stim._sound = set_dbspl(stim.sound, level)  # adjust scaling here
                spikes[rep].append(
                        an_model.get_spiketrain(
                            cf=cf, sr=sr, seed=seed, stim=stim, simulator=simulator
                            )
                        )
                seed += 1
        result[sr] = {"levels": levels, "spikes": spikes}
    return result


def runtest():
    parser = argparse.ArgumentParser(description="Sound stimulation spike-train test")
    parser.add_argument(
        "--simulator",
        type=str,
        choices=["py3", "matlab"],
        default="py3",
        help="AN simulator backend (default: py3)",
    )
    args = parser.parse_args()
    usematlab = args.simulator == "matlab"
    reps = 10
    pip_duration = 0.5 # seconds
    result = sound_stim(seed, reps=reps, pip_duration=pip_duration, useMatlab=usematlab)

    win = pg.GraphicsLayoutWidget()
    p1 = win.addPlot(title="Rate-level function")
    p1.setLabel("bottom", "Tone level (dB SPL)")
    p1.setLabel("left", "Firing Rate (sp/s)")
    group = {0: "LSR", 1: "MSR", 2: "HSR"}
    p1.addLegend()
    for i, x in enumerate(result.keys()):  # across srs
        for r in range(reps):
            if r == 0:
                name = group[x]
            else:
                name = None
            p1.plot(result[x]['levels'], [
                s.size/pip_duration for s in result[x]["spikes"][r]], 
                symbol = 'o', symbolBrush=0.25, symbolPen=(x,6), symbolSize=3,
                pen=(x, 6),
                name=name)

    return win


if __name__ == "__main__":
    win = runtest()
    win.show()
    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()
