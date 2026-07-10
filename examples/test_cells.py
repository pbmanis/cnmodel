#!/usr/bin/python

"""
Test the basic membrane physiology of cell types.

Basic Usage:   python test_cells.py celltype species [--cc | --vc]

This script generates a cell of the specified type and species, then tests the
cell with a series of current/voltage pulses to produce I/V, F/I, and spike
latency analyses.

"""
import argparse
import os
import sys

import numpy as np
import pyqtgraph as pg
from neuron import h

import cnmodel
from cnmodel import cells as cells
from cnmodel import data as cnmodel_data
from cnmodel.protocols import IVCurve, VCCurve

debugFlag = True
ax = None
h.celsius = 22
cclamp = False

cellinfo = {
    "types": [
        "bushy",
        "bushycoop",
        "tstellate",
        "tstellatenav11",
        "dstellate",
        "dstellateeager",
        "sgc",
        "cartwheel",
        "pyramidal",
        "pyramidalceballos",
        "octopus",
        "tuberculoventral",
        "mso",
        "granule",
    ],
    "morphology": ["point", "waxon", "stick"],
    "nav": ["std", "jsrna", "nav11", "nacncoop", "grc_na"],
    "species": ["guineapig", "cat", "rat", "mouse"],
    "pulse": ["step", "pulse"],
}

# Format for ivranges is list of tuples. This allows finer increments in selected ranges, such as close to rest
ccivrange = {
    "mouse": {
        "bushy": {"pulse": [(-1, 1.2, 0.05)]},
        "bushycoop": {"pulse": [(-0.5, 0.7, 0.02)]},
        "tstellate": {"pulse": [(-1.0, 1.01, 0.05), (-0.015, 0, 0.005)]},
        "tstellatenav11": {"pulse": [(-1, 1.0, 0.1)]},
        "dstellate": {"pulse": [(-0.3, 0.301, 0.015)]},
        "octopus": {"pulse": [(-1.0, 1.0, 0.05)]},
        "sgc": {"pulse": [(-0.3, 0.6, 0.02)]},
        "cartwheel": {"pulse": [(-0.5, 0.5, 0.05)]},
        "pyramidal": {"pulse": [(-0.3, 0.3, 0.025), (-0.040, 0.025, 0.005)]},
        "pyramidalceballos": {
            "pulse": [(-0.09, 0.00, 0.09), (0, 0.008, 0.008)]
        },  # , 'prepulse': [(-0.25, -0.25, 0.25)]},
        "tuberculoventral": {"pulse": [(-0.35, 1.0, 0.05), (-0.040, 0.01, 0.005)]},
        "granule": {"pulse": [(-0.02, 0.02, 0.005)]},
    },
    "guineapig": {
        "bushy": {"pulse": [(-1, 1.2, 0.05)]},
        "tstellate": {"pulse": [(-0.15, 0.15, 0.01)]},
        "dstellate": {"pulse": [(-0.25, 0.25, 0.025)]},
        "dstellateeager": {"pulse": [(-0.6, 1.0, 0.025)]},
        "octopus": {"pulse": [(-2.0, 6.0, 0.2)]},
        "sgc": {"pulse": [(-0.3, 0.6, 0.02)]},
        "mso": {"pulse": [(-1, 1.2, 0.05)]},
    },
    "rat": {
        "pyramidal": {
            "pulse": [(-0.3, 0.3, 0.025), (-0.040, 0.025, 0.005)]
        },  # 'prepulse': [(-0.25, -0.25, 0.25)]},
    },
}

# scales holds some default scaling to use in the cciv plots
# argument is {cellname: (xmin, xmax, IVymin, IVymax, FIspikemax,
# offset(for spikes), crossing (for IV) )}
## the "offset" refers to setting the axes back a bit
scale = {
    "bushy": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "bushycoop": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "tstellate": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "tstellatenav11": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "tstellatedend": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "dstellate": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "dstellateeager": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "sgc:": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "cartwheel": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "pyramidal": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "pyramidalceballos": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "tuberculoventral": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "granule": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "octopus": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
    "mso": (-1.0, -160.0, 1.0, -40, 0, 40, "offset", 5, "crossing", [0, -60]),
}


def _count_iv_steps(ivrange: dict) -> int:
    """Count total current-command steps from an IVCurve range dict."""
    total = 0
    for c in ivrange.get('pulse', []):
        imin, imax, istep = c
        total += int(np.floor((imax - imin) / istep) + 1)
    return total


class Tests:
    """
    Class to select cells for tests
    """

    def __init__(self):
        pass

    def selectCell(self, args):
        """
        Parameters
        ----------
        args : argparse args from command line

        Returns
        -------
        cell
            Instantiated cell of the selected celltype
        """
        h.celsius = float(args.temp)
        #
        # Spiral Ganglion cell tests
        #
        if args.celltype == "sgc":  # morphology is always "point" for SGCs
            cell = cells.SGC.create(
                debug=debugFlag,
                species=args.species,
                nach=args.nav,
                ttx=args.ttx,
                modelType=args.type,
            )
        #
        # Bushy tests
        #
        elif args.celltype == "bushy" and args.morphology == "point":
            cell = cells.Bushy.create(
                model="RM03",
                species=args.species,
                modelType=args.type,
                ttx=args.ttx,
                nach=args.nav,
                debug=debugFlag,
            )
        #            cell.soma().klt.gbar = 0.0003

        elif args.celltype == "bushy" and args.morphology == "waxon":
            cell = cells.Bushy.create(
                model="RM03",
                species=args.species,
                modelType=args.type,
                nach=args.nav,
                ttx=args.ttx,
                debug=debugFlag,
            )
            cell.add_axon()

        elif args.celltype == "bushy" and args.morphology == "stick":
            # XM13/mouse is the only combo with a compartments table;
            # RM03/guineapig has no compartments table so the decorator cannot work.
            cell = cells.Bushy.create(
                model="XM13",
                species="mouse",
                modelType=args.type,
                morphology="cnmodel/morphology/bushy_stick.hoc",
                decorator=True,
                nach=args.nav,
                ttx=args.ttx,
                debug=debugFlag,
            )
            h.topology()

        elif args.celltype == "bushycoop" and args.morphology == "point":
            cell = cells.Bushy.create(
                model="RM03",
                species=args.species,
                modelType=args.type,
                ttx=args.ttx,
                nach=args.nav,
                debug=debugFlag,
            )

        #
        # T-stellate tests
        #
        elif args.celltype == "tstellate" and args.morphology == "point":
            cell = cells.TStellate.create(
                model="RM03",
                species=args.species,
                modelType=args.type,
                nach=args.nav,
                ttx=args.ttx,
                debug=debugFlag,
            )

        elif args.celltype == "tstellate" and args.morphology == "stick":
            # no tstellate compartments table; uniform fallback used.
            cell = cells.TStellate.create(
                model="RM03",
                species=args.species,
                modelType=args.type,
                nach=args.nav,
                ttx=args.ttx,
                debug=debugFlag,
                morphology="cnmodel/morphology/tstellate_stick.hoc",
                decorator=True,
            )

        elif (
            args.celltype == "tstellatenav11" and args.morphology == "point"
        ):  # note this uses a different model...
            print("test_cells: Stellate NAV11")
            cell = cells.TStellateNav11.create(
                model="Nav11", species=args.species, modelType=None, ttx=args.ttx, debug=debugFlag
            )

        elif (
            args.celltype == "tstellatenav11" and args.morphology == "stick"
        ):  # note this uses a different model...
            #  no tstellatenav11 compartments table; uniform fallback used.
            cell = cells.TStellateNav11.create(
                model="Nav11",
                species=args.species,
                modelType=None,
                morphology="cnmodel/morphology/tstellate_stick.hoc",
                decorator=True,
                ttx=args.ttx,
                debug=debugFlag,
            )
            h.topology()

        #
        # Octopus cell tests
        #
        elif args.celltype == "octopus" and args.morphology == "point":
            cell = cells.Octopus.create(
                species=args.species,
                modelType="RM03",  # args.type,
                nach=args.nav,
                ttx=args.ttx,
                debug=debugFlag,
            )

        elif (
            args.celltype == "octopus" and args.morphology == "stick"
        ):  # Go to spencer et al. model
            # no octopus/Spencer compartments table; uniform fallback used.
            # Spencer model only supports mouse; nach defaults to 'nacn' if not specified.
            cell = cells.Octopus.create(
                modelType="Spencer",
                species="mouse",
                morphology="cnmodel/morphology/octopus_spencer_stick.hoc",
                decorator=True,
                nach=args.nav if args.nav is not None else 'nacn',
                ttx=args.ttx,
                debug=debugFlag,
            )
            h.topology()

        #
        # D-stellate tests
        #
        elif args.celltype == "dstellate":
            cell = cells.DStellate.create(
                debug=debugFlag, species=args.species, ttx=args.ttx, modelType=args.type
            )

        elif args.celltype == "dstellateeager":
            cell = cells.DStellateEager.create(debug=debugFlag, ttx=args.ttx, modelType=args.type)

        #
        # DCN pyramidal cell tests
        #
        elif args.celltype == "pyramidal":
            cell = cells.Pyramidal.create(
                model="POK", modelType=args.type, ttx=args.ttx, debug=debugFlag
            )

        elif args.celltype == "pyramidalceballos":
            cell = cells.Pyramidal.create(
                modelType=args.type, model="Ceballos", ttx=args.ttx, debug=debugFlag
            )
            cell.vm0 = -65.0  # cell may not have a stable state, so force rmp

        #
        # DCN tuberculoventral cell tests
        #
        elif args.celltype == "tuberculoventral" and args.morphology == "point":
            cell = cells.Tuberculoventral.create(
                species="mouse", modelType="TVmouse", ttx=args.ttx, nach=args.nav, debug=debugFlag
            )

        elif args.celltype == "tuberculoventral" and args.morphology == "stick":
            # no TVmouse compartments table; uniform fallback used.
            cell = cells.Tuberculoventral.create(
                species="mouse",
                modelType="TVmouse",
                morphology="cnmodel/morphology/tv_stick.hoc",
                decorator=True,
                ttx=args.ttx,
                debug=debugFlag,
            )
            h.topology()

        #
        # DCN granule cell tests
        #
        elif args.celltype == "granule" and args.morphology == "point":
            print("Dir cells: ", dir(cells))
            cell = cells.Granule.create(
                species="mouse", modelType="GRC", ttx=args.ttx, nach=args.nav, debug=debugFlag
            )

        # added decorator=True, because that channel_manager is wired for GRC
        elif args.celltype == 'granule' and args.morphology == 'stick':
            cell = cells.Granule.create(species='mouse', modelType='GRC',
                    morphology='cnmodel/morphology/granule_stick_diwakar.hoc',
                    decorator=True,
                    ttx=args.ttx, nach=args.nav, debug=debugFlag)
            # h.topology()

        #
        # DCN cartwheel cell tests
        #
        elif args.celltype == "cartwheel":
            cell = cells.Cartwheel.create(modelType=args.type, ttx=args.ttx, debug=debugFlag)

        #
        # MSO principal neuron tests
        #
        elif args.celltype == "mso" and args.morphology == "point":
            print("MSO Creation")
            cell = cells.MSO.create(
                model="MSO-principal",
                species=args.species,
                modelType=args.type,
                ttx=args.ttx,
                nach=args.nav,
                debug=debugFlag,
            )

        else:
            avail = cnmodel_data.report_available_configurations(
                celltype=args.celltype, species=args.species
            )
            raise ValueError(
                f"Cell type '{args.celltype}' with nav='{args.nav}', "
                f"morphology='{args.morphology}', species='{args.species}' "
                f"is not a defined configuration.\n"
                f"Available cell types: "
                f"{sorted(k for k in cells.__dict__ if not k.startswith('_'))}\n"
                + avail
            )

        print(cell.__doc__)
        self.cell = cell

    def run_test(self, sites, ptype, args):
        """
        Run either vc or cc test, and plot the result

        Parameters
        ----------
        args : argparse args from command line

        """
        self.cell.set_temperature(float(args.temp))
        print("Cell status dictionary:\n", self.cell.status)
        print("Resting potential: ", self.cell.vm0)
        durations = eval(args.durations)
        if self.cell.vm0 is None:
            V0 = self.cell.find_i0(showinfo=True)
        else:
            V0 = self.cell.vm0
        #        self.cell.cell_initialize()
        print("Currents at nominal Vrest= %.2f I = 0: I = %g " % (V0, self.cell.i_currents(V=V0)))
        self.cell.print_mechs(self.cell.soma)
        instant = self.cell.compute_rmrintau(auto_initialize=False, vrange=None)
        print(
            "    From Inst: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}".format(
                instant["Rin"], instant["tau"], instant["v"]
            )
        )
        if args.cc is True:
            # define the current clamp electrode and default settings
            self.iv = IVCurve()
            self._run_cc_with_progress(
                ivrange=ccivrange[args.species][args.celltype],
                cell=self.cell,
                durs=durations,
                sites=sites,
                reppulse=ptype,
                temp=float(args.temp),
            )
            ret = self.iv.input_resistance_tau()
            if not np.isnan(ret["slope"]):
                print(
                    "    From IV: Rin = {:7.1f}  Tau = {:7.1f}  Vm = {:7.1f}".format(
                        ret["slope"], ret["tau"], ret["intercept"]
                    )
                )
            self.iv.show(cell=self.cell)

        elif args.rmp is True:
            print("temperature: ", self.cell.status["temperature"])
            self.iv = IVCurve()
            self.iv.run(
                {"pulse": (0, 0, 1)},
                self.cell,
                durs=durations,
                sites=sites,
                reppulse=ptype,
                temp=float(args.temp),
            )
            self.iv.show(cell=self.cell, rmponly=True)

        elif args.vc is True:
            # define the voltage clamp electrode and default settings
            self.vc = VCCurve()
            self.vc.run((-120, 40, 5), self.cell)
            self.vc.show(cell=self.cell)

        else:
            print("No test mode specified (--cc/--vc/--rmp); showing Current Plots button only.")

        self._add_current_plots_button(args, durations)


    def _run_cc_with_progress(self, ivrange, cell, durs, sites, reppulse, temp):
        """
        Call self.iv.run() while showing a progress bar + Cancel button.

        custom_init() calls h.finitialize() exactly twice per IV step, so
        h.FInitializeHandler fires twice per step.  We count every 2nd call
        as one completed step.

        Cancel sets h.stoprun=1, which causes h.batch_run() (called inside
        run_one) to return immediately for the current and all subsequent steps.
        The remaining steps still run their warmup but complete near-instantly.
        """
        app = pg.QtWidgets.QApplication.instance()
        nsteps = _count_iv_steps(ivrange)

        # Build the progress panel.
        panel = pg.QtWidgets.QWidget()
        panel.setWindowTitle("Running IV curve")
        vlayout = pg.QtWidgets.QVBoxLayout(panel)
        vlayout.setSpacing(8)
        vlayout.setContentsMargins(16, 12, 16, 12)

        hdr = pg.QtWidgets.QLabel(f"Running {nsteps} current-clamp steps…")
        vlayout.addWidget(hdr)

        pbar = pg.QtWidgets.QProgressBar()
        pbar.setRange(0, nsteps)
        pbar.setValue(0)
        vlayout.addWidget(pbar)

        step_lbl = pg.QtWidgets.QLabel(f"Step 0 of {nsteps}")
        vlayout.addWidget(step_lbl)

        cancel_btn = pg.QtWidgets.QPushButton("Cancel")
        vlayout.addWidget(cancel_btn)

        panel.resize(440, 155)
        panel.show()
        if app:
            app.processEvents()

        state = {'finit_count': 0, 'step': 0, 'cancelled': False}

        def _on_cancel():
            state['cancelled'] = True
            h.stoprun = 1
            cancel_btn.setEnabled(False)
            cancel_btn.setText("Cancelling…")

        cancel_btn.clicked.connect(_on_cancel)

        # FInitializeHandler fires at each h.finitialize() call.
        # custom_init() makes exactly 2 per IV step; every 2nd call = 1 step done.
        def _on_finit():
            if state['step'] >= nsteps:  # Claude fixed 2026-07-10: ignore extra finitialize() calls after run completes
                return
            state['finit_count'] += 1
            if state['finit_count'] % 2 == 0:
                state['step'] += 1
                pbar.setValue(state['step'])
                step_lbl.setText(f"Step {state['step']} of {nsteps}")
            if app:
                app.processEvents()
            if state['cancelled']:
                h.stoprun = 1

        fih = h.FInitializeHandler(1, _on_finit)

        try:
            self.iv.run(
                ivrange, cell,
                durs=durs, sites=sites, reppulse=reppulse, temp=temp,
            )
        finally:
            del fih
            panel.close()

    def _add_current_plots_button(self, args, durations):
        from cnmodel.util.current_plots import CurrentPlotsWindow
        iamp_nA = args.iamp / 1000.0

        self._cur_plot_win = None
        cell = self.cell
        tests_ref = self

        # Build a horizontal bar: [Level: <combo>] [Current Plots]
        # Claude fixed 2026-07-09: combo lets the user pick any simulated current level.
        bar = pg.QtWidgets.QWidget()
        bar_layout = pg.QtWidgets.QHBoxLayout(bar)
        bar_layout.setContentsMargins(4, 2, 4, 2)

        combo = None
        cmd_levels = getattr(getattr(self, 'iv', None), 'current_cmd', None)
        if cmd_levels is not None:
            bar_layout.addWidget(pg.QtWidgets.QLabel('Level:'))
            combo = pg.QtWidgets.QComboBox()
            for level in cmd_levels:
                combo.addItem(f'{level * 1000:.1f} pA')
            # Pre-select the entry closest to the command-line --iamp value.
            best = int(np.argmin(np.abs(cmd_levels - iamp_nA)))
            combo.setCurrentIndex(best)
            bar_layout.addWidget(combo)

        btn = pg.QtWidgets.QPushButton('Current Plots')
        bar_layout.addWidget(btn)

        def _on_click():
            if combo is not None and cmd_levels is not None:
                level = float(cmd_levels[combo.currentIndex()])
            else:
                level = iamp_nA
            tests_ref._cur_plot_win = CurrentPlotsWindow(
                cell, iamp_nA=level, durations=durations)

        btn.clicked.connect(_on_click)

        # Embed bar in the existing summary window (iv or vc) when available;
        # otherwise show it as a standalone toolbar.
        summary_win = None
        if hasattr(self, 'iv') and hasattr(self.iv, 'win'):
            summary_win = self.iv.win
        elif hasattr(self, 'vc') and hasattr(self.vc, 'win'):
            summary_win = self.vc.win

        if summary_win is not None:
            win_w = summary_win.width()
            win_h = summary_win.height()
            outer = pg.QtWidgets.QWidget()
            outer.setWindowTitle(summary_win.windowTitle())
            vbox = pg.QtWidgets.QVBoxLayout(outer)
            vbox.setContentsMargins(0, 0, 0, 0)
            vbox.setSpacing(2)
            vbox.addWidget(summary_win)  # reparents the plot into this container
            vbox.addWidget(bar)
            outer.resize(win_w, win_h + 40)
            outer.show()
            self._tools_widget = outer
        else:
            bar.setWindowTitle('Tools')
            bar.resize(320, 50)
            bar.show()
            self._tools_widget = bar


def main():
    parser = argparse.ArgumentParser(
        description=(
            "test_cells.py:",
            " Biophysical representations of neurons (mostly auditory), test file",
        )
    )
    parser.add_argument("celltype", action="store")
    parser.add_argument("species", action="store", default="guineapig")
    parser.add_argument("--type", action="store", default=None)
    parser.add_argument("--temp", action="store", default=22.0, help=("Temp DegC (22 default)"))
    parser.add_argument(
        "-m",
        action="store",
        dest="morphology",
        default="point",
        help=("Set morphology: %s " % [morph for morph in cellinfo["morphology"]]),
    )
    parser.add_argument(
        "--nav",
        action="store",
        dest="nav",
        default=None,
        help=("Choose sodium channel: %s " % [ch for ch in cellinfo["nav"]]),
    )
    parser.add_argument(
        "--ttx", action="store_true", dest="ttx", default=False, help=("Use TTX (no sodium current")
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="pulsetype",
        default="step",
        help=("Set CCIV pulse to step or repeated pulse"),
    )
    parser.add_argument(
        "-d",
        "--durations",
        action="store",
        dest="durations",
        default="[10, 100, 50]",
        help=("Set pulse durations in format '[10, 100, 20]' (as a string)"),
    )
    parser.add_argument(
        "--iamp",
        action="store",
        type=float,
        dest="iamp",
        default=20.0,
        help="Current injection amplitude for Current Plots (pA, default 20)",
    )
    clampgroup = parser.add_mutually_exclusive_group()
    clampgroup.add_argument("--vc", action="store_true", help="Run in voltage clamp mode")
    clampgroup.add_argument("--cc", action="store_true", help="Run in current clamp mode")
    clampgroup.add_argument(
        "--rmp", action="store_true", help="Run to get RMP in current clamp mode"
    )

    args = parser.parse_args()
    try:
        eval(args.durations)
    except:
        raise ValueError(
            "Durations values could not be parsed\nFor example, use: -d '[10,100,10]' in quotes"
        )

    if args.celltype not in cellinfo["types"]:
        print("cell: %s is not in our list of cell types" % (args.celltype))
        print("celltypes: ", cellinfo["types"])
        sys.exit(1)

    path = os.path.dirname(cnmodel.__file__)
    #    h.load_file("stdrun.hoc")
    #    h.load_file(os.path.join(path, "custom_init.hoc")) # replace init with one that gets closer to steady state

    # print 'Species: ', args.species
    # print 'Morphology configuration: ', args.morphology
    sites = None
    if args.pulsetype == "step":
        ptype = None
    else:
        ptype = "pulses"
    if args.morphology in cellinfo["morphology"]:
        print("Morphological configuration %s is ok" % args.morphology)

    t = Tests()
    t.selectCell(args)
    app = pg.mkQApp()
    t.run_test(sites, ptype, args)

    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()


if __name__ == "__main__":
    main()
