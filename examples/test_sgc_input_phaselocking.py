"""
test_sgc_input_phaselocking.py

Test phase locking from an input sgc to a target cell type. Runs simulations
with AN input, and plots the results, including PSTH and phase histogram.

usage: test_sgc_input_phaselocking.py [-h]
                                      [-c {bushy,tstellate,octopus,dstellate}]
                                      [-S {tone,SAM,clicks}]
                                      [-s {guineapig,rat,mouse}]

test sgc input phaselocking

optional arguments:
  -h, --help            show this help message and exit
  -c {bushy,tstellate,octopus,dstellate}, --celltype {bushy,tstellate,octopus,dstellate}
                        cell type
  -S {tone,SAM,clicks}, --stimulus {tone,SAM,clicks}
                        stimulus type
  -s {guineapig,rat,mouse}, --species {guineapig,rat,mouse}
                        species

Note: Not all combinations of inputs are valid (not all cell types are
known for each species)

"""

import sys
import numpy as np
import argparse

from neuron import h
import pyqtgraph as pg
from cnmodel.protocols import Protocol
from cnmodel import cells
from cnmodel.util import sound
import cnmodel.util.pynrnutilities as PU
from cnmodel import data


class SGCInputTestPL(Protocol):
    def set_cell(self, cell="bushy"):
        self.cell = cell

    def run(
        self, args, temp=34.0, dt=0.025, seed=575982035, dB=None,
    ):
        if self.cell == "bushy":
            postCell = cells.Bushy.create(species=args.species)
        elif self.cell == "tstellate":
            postCell = cells.TStellate.create(species=args.species)
        elif self.cell == "octopus":
            postCell = cells.Octopus.create(species=args.species)
        elif self.cell == "dstellate":
            postCell = cells.DStellate.create(species=args.species)
        else:
            raise ValueError(
                "cell %s is not yet implemented for phaselocking" % self.cell
            )
        self.post_cell = postCell
        self.species = args.species
        self.stimulus = args.stimulus
        self.run_duration = 1.0  # in seconds
        self.pip_duration = 0.8  # in seconds
        self.pip_start = [0.02]  # in seconds
        self.Fs = 100e3  # in Hz
        self.f0 = args.CF  # stimulus in Hz
        self.cf = args.CF  # SGCs in Hz
        self.fMod = args.fmod  # mod freq, Hz
        self.dMod = args.dmod  # % mod depth, percentage
        if dB is None:
            self.dbspl = args.dB
        else:
            self.dbspl = dB
        if self.stimulus == "SAM":
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
            self.vs_freq = self.fMod
        if self.stimulus == "tone":
            self.f0 = 1000.0
            self.cf = 1000.0
            self.stim = sound.TonePip(
                rate=self.Fs,
                duration=self.run_duration,
                f0=self.f0,
                dbspl=self.dbspl,
                ramp_duration=2.5e-3,
                pip_duration=self.pip_duration,
                pip_start=self.pip_start,
            )
            self.vs_freq = self.f0
        if self.stimulus == "clicks":
            self.click_rate = 0.020  # msec
            self.stim = sound.ClickTrain(
                rate=self.Fs,
                duration=self.run_duration,
                dbspl=self.dbspl,
                click_starts=np.linspace(
                    0.01,
                    self.run_duration - 0.01,
                    int((self.run_duration) / self.click_rate),
                ),
                click_duration=100.0e-6,
                # click_interval=self.click_rate, nclicks=int((self.run_duration-0.01)/self.click_rate),
                ramp_duration=2.5e-3,
            )
        n_sgc = data.get(
            "convergence", species=self.species, post_type=postCell.celltype, pre_type="sgc"
        )[0]
        self.n_sgc = int(np.round(n_sgc))

        self.pre_cells = []
        self.synapses = []
        j = 0
        for k in range(self.n_sgc):
            seed = seed + k
            preCell = cells.DummySGC(cf=self.cf, sr=2)
            synapse = preCell.connect(postCell)
            for i in range(synapse.terminal.n_rzones):
                self["xmtr%03d" % j] = synapse.terminal.relsite._ref_XMTR[i]
                j = j + 1
            synapse.terminal.relsite.Dep_Flag = False
            preCell.set_sound_stim(self.stim, seed=seed)
            self.pre_cells.append(preCell)
            self.synapses.append(synapse)

        self["vm"] = postCell.soma(0.5)._ref_v
        # self['prevm'] = preCell.soma(0.5)._ref_v
        self["t"] = h._ref_t
        postCell.cell_initialize()
        h.tstop = 1e3 * self.run_duration  # duration of a run
        h.celsius = temp
        h.dt = dt

        self.custom_init()
        h.run()

    def window_spikes(self, spiketrain):
        phasewin = [
            self.pip_start[0] + 0.25 * self.pip_duration,
            self.pip_start[0] + self.pip_duration,
        ]
        spkin = spiketrain[np.where(spiketrain > phasewin[0] * 1e3)]
        spikesinwin = spkin[np.where(spkin <= phasewin[1] * 1e3)]
        return spikesinwin
    
    def compute_vs(self):
        self.post_spikes = PU.findspikes(self["t"], self["vm"], -30.0)
        self.post_spikes_win = self.window_spikes(self.post_spikes)
        self.an_spikes_win = self.window_spikes(self.pre_cells[0]._spiketrain)  # just sample one...

        # set freq for VS calculation
        if self.stimulus == "tone":
            f0 = self.f0
            print(
                "Tone: f0=%.3f at %3.1f dbSPL, cell CF=%.3f"
                % (self.f0, self.dbspl, self.cf)
            )
        if self.stimulus == "SAM":
            f0 = self.fMod
            print(
                (
                    "SAM Tone: f0=%.3f at %3.1f dbSPL, fMod=%3.1f  dMod=%5.2f, cell CF=%.3f"
                    % (self.f0, self.dbspl, self.fMod, self.dMod, self.cf)
                )
            )
        if self.stimulus == "clicks":
            f0 = 1.0 / self.click_rate
            print(
                "Clicks: interval %.3f at %3.1f dbSPL, cell CF=%.3f "
                % (self.click_rate, self.dbspl, self.cf)
            )
        self.an_vs = PU.vector_strength(self.an_spikes_win*1e-3, f0)
        andiff = self.an_spikes_win*1e-3
        print(
            "AN Vector Strength at %.1f: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
            % (f0, self.an_vs["r"], self.an_vs["d"] * 1e6, self.an_vs["R"], self.an_vs["p"], self.an_vs["n"])
        )
        self.post_cell_vs = PU.vector_strength(self.post_spikes_win*1e-3, f0)
        print(
            "%s Vector Strength: %7.3f, d=%.2f (us) Rayleigh: %7.3f  p = %.3e  n = %d"
            % (self.cell, self.post_cell_vs["r"], self.post_cell_vs["d"] * 1e6,
             self.post_cell_vs["R"], self.post_cell_vs["p"], self.post_cell_vs["n"])
        )
        

    def show(self):
        self.compute_vs()
        self.win = pg.GraphicsWindow()
        p1 = self.win.addPlot(title="stim", row=0, col=0)
        p1.plot(self.stim.time * 1000, self.stim.sound)
        p1.setXLink(p1)

        p2 = self.win.addPlot(title="AN spikes", row=1, col=0)
        vt = pg.VTickGroup(self.pre_cells[0]._spiketrain)
        p2.addItem(vt)
        p2.setXLink(p1)

        p3 = self.win.addPlot(title="%s Spikes" % self.cell, row=2, col=0)
        bspktick = pg.VTickGroup(self.post_spikes)
        p3.addItem(bspktick)
        p3.setXLink(p1)

        p4 = self.win.addPlot(title="%s Vm" % self.cell, row=3, col=0)
        p4.plot(self["t"], self["vm"])
        p4.setXLink(p1)

        p5 = self.win.addPlot(title="xmtr", row=0, col=1)
        j = 0
        for k in range(self.n_sgc):
            synapse = self.synapses[k]
            for i in range(synapse.terminal.n_rzones):
                p5.plot(self["t"], self["xmtr%03d" % j], pen=(i, 15))
                j = j + 1
        p5.setXLink(p1)

        p6 = self.win.addPlot(title="AN phase", row=1, col=1)
        # phasewin = [
        #     self.pip_start[0] + 0.25 * self.pip_duration,
        #     self.pip_start[0] + self.pip_duration,
        # ]
        print("\nCell type: %s" % self.cell)
        print("Stimulus: ")


        (hist, binedges) = np.histogram(self.an_vs["ph"])
        p6.plot(
            binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0
        )
        p6.setXRange(0.0, 2 * np.pi)

        p7 = self.win.addPlot(title="%s phase" % self.cell, row=2, col=1)
        
        (hist, binedges) = np.histogram(self.post_cell_vs["ph"])
        p7.plot(
            binedges, hist, stepMode=True, fillBrush=(255, 100, 100, 150), fillLevel=0
        )
        p7.setXRange(0.0, 2 * np.pi)
        p7.setXLink(p6)

        self.win.show()


def main():
    parser = argparse.ArgumentParser(description="test sgc input phaselocking")
    parser.add_argument(
        "-c",
        "--celltype",
        type=str,
        choices=["bushy", "tstellate", "octopus", "dstellate"],
        default="bushy",
        help="cell type",
    )
    parser.add_argument(
        "--species",
        type=str,
        choices=["guineapig", "mouse", "rat"],
        default="mouse",
        help="Species",
    )
    parser.add_argument(
        "-S",
        "--stimulus",
        type=str,
        choices=["tone", "SAM", "clicks",],
        default="tone",
        help="stimulus type",
    )
    parser.add_argument(
        "--dB",
        "--dBSPL",
        type=float,
        default=60.,
        help="Sound pressure level, SPL",
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

    prot = SGCInputTestPL()
    prot.set_cell(args.celltype)
    
    if not args.RI:
        prot.run(args) # stimulus=args.stimulus, species=args.species)
        prot.show()
    else:
        an_vs = []
        post_vs = []
        dbrange = np.linspace(0, 70, 15)
        for db in dbrange:
            prot.run(args, dB=db)
            prot.compute_vs()
            an_vs.append(prot.an_vs["r"])
            post_vs.append(prot.post_cell_vs["r"])
        print(f" {'dB':3s}  {'vsan':6s}  {'vsbu':6s}")
        for i, db in enumerate(dbrange):
            print(f" {int(db):3d}  {an_vs[i]:5.2f}  {post_vs[i]:5.2f}")

    import sys

    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()


if __name__ == "__main__":
    main()
