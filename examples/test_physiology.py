"""
Test principal cell responses to tone pips of varying frequency and intensity.

This is an example of model construction from a very high level--we specify
only which populations of cells are present and which ones should be connected.
The population and cell classes take care of all the details of generating the
network needed to support a small number of output cells.

Note: run time for this example can be very long. To speed things up, reduce
n_frequencies or n_levels, or reduce the number of selected output cells (see
cells_per_band).
 
"""
import argparse
import os
import pickle
import sys
import timeit
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Any, Callable, Dict

import numpy as np
import pyqtgraph as pg
from neuron import h
from pyqtgraph import multiprocess as mp

# from pyqtgraph.Qt import QtGui, QtCore
from cnmodel import populations
from cnmodel.protocols import Protocol
from cnmodel.util import random_seed, sound
from cnmodel.util.network_widgets import NetworkControlPanel


CELLTYPES = ["bushy", "tstellate", "dstellate", "pyramidal", "tuberculoventral"]


# ---------------------------------------------------------------------------
# Stimulus type registry
# ---------------------------------------------------------------------------

@dataclass
class StimDef:
    """Describes one stimulus type: how to build it and label it."""
    name: str                           # identifier used in cache keys
    label: str                          # display string in GUI
    factory: Callable                   # factory(stimpar, f0, dbspl, seed, **extra) → Sound
    extra: Dict[str, Any] = field(default_factory=dict)  # type-specific params


def _tone_pip_factory(stimpar, f0, dbspl, seed, **kw):
    return sound.TonePip(
        rate=100e3,
        duration=stimpar["dur"],
        f0=f0,
        dbspl=dbspl,
        ramp_duration=2.5e-3,
        pip_duration=stimpar["pip"],
        pip_start=stimpar["start"],
    )


def _bl_noise_factory(stimpar, f0, dbspl, seed, octave_width=1.0, **kw):
    return sound.BandlimitedNoisePip(
        rate=100e3,
        duration=stimpar["dur"],
        center_freq=f0,
        dbspl=dbspl,
        ramp_duration=2.5e-3,
        pip_duration=stimpar["pip"],
        pip_start=stimpar["start"],
        seed=seed,
        octave_width=octave_width,
    )


def _wb_noise_factory(stimpar, f0, dbspl, seed, **kw):
    stim = sound.NoisePip(
        rate=100e3,
        duration=stimpar["dur"],
        dbspl=dbspl,
        ramp_duration=2.5e-3,
        pip_duration=stimpar["pip"],
        pip_start=stimpar["start"],
        seed=seed,
    )
    stim.opts['f0'] = f0  # store centre frequency so _stim_freq() returns a valid value
    return stim


STIM_DEFS: Dict[str, StimDef] = {
    "tone_pip":          StimDef("tone_pip",          "Tone Pip",           _tone_pip_factory),
    "bandlimited_noise": StimDef("bandlimited_noise", "Band-Limited Noise", _bl_noise_factory),
    "wideband_noise":    StimDef("wideband_noise",    "Wide-Band Noise",    _wb_noise_factory),
}


class CNSoundStim(Protocol):
    def __init__(self, seed, temp=34.0, dt=0.025, synapsetype="simple",
                 celltype="bushy", cf=16e3, species='mouse'):
        Protocol.__init__(self)

        self.seed = seed
        self.temp = temp
        self.dt = dt
        self.celltype = celltype
        self.species = species
        self.record_syn_conductance = False  # set from main() via control panel checkbox

        # Seed now to ensure network generation is stable
        random_seed.set_seed(seed)

        # SGC population is always present
        self.sgc = populations.SGC(model="dummy", species=species)
        # set synapse type to use in the sgc population - simple is fast, multisite is slower
        self.sgc._synapsetype = synapsetype

        frequencies = [cf]
        cells_per_band = 1
        self._build_topology(celltype, frequencies, cells_per_band, species=species)

    def _build_topology(self, celltype, frequencies, cells_per_band, species='mouse'):
        """Create populations and symbolic connections for the given target cell type.

        Sets self.target, self.recording_pops, self.pre_pops, and self.populations.
        Does NOT call resolve_inputs(); call _resolve_network() after configuring
        per-connection synapse types via _synapsetype_per_pre on each post-population.
        """
        if celltype == "bushy":
            self.dstellate = populations.DStellate(species=species)
            self.tstellate = populations.TStellate(species=species)
            self.tuberculoventral = populations.Tuberculoventral(species=species)
            self.bushy = populations.Bushy(species=species)
            self.target = self.bushy

            self.sgc.connect(self.bushy, self.dstellate, self.tuberculoventral, self.tstellate)
            self.dstellate.connect(self.bushy, self.tstellate)
            self.tuberculoventral.connect(self.bushy, self.tstellate)
            self.tstellate.connect(self.bushy)

            pops = [self.sgc, self.dstellate, self.tuberculoventral, self.tstellate, self.bushy]
            self.pre_pops = [self.sgc, self.dstellate, self.tuberculoventral]
            self.recording_pops = [self.bushy, self.dstellate, self.tstellate, self.tuberculoventral]

        elif celltype == "tstellate":
            self.dstellate = populations.DStellate(species=species)
            self.tuberculoventral = populations.Tuberculoventral(species=species)
            self.tstellate = populations.TStellate(species=species)
            self.target = self.tstellate

            self.sgc.connect(self.tstellate, self.dstellate, self.tuberculoventral)
            self.dstellate.connect(self.tstellate)
            self.tuberculoventral.connect(self.tstellate)

            pops = [self.sgc, self.dstellate, self.tuberculoventral, self.tstellate]
            self.pre_pops = [self.sgc, self.dstellate, self.tuberculoventral]
            self.recording_pops = [self.tstellate, self.dstellate, self.tuberculoventral]

        elif celltype == "dstellate":
            self.dstellate = populations.DStellate(species=species)
            self.target = self.dstellate

            self.sgc.connect(self.dstellate)

            pops = [self.sgc, self.dstellate]
            self.pre_pops = [self.sgc]
            self.recording_pops = [self.dstellate]

        elif celltype == "pyramidal":
            self.dstellate = populations.DStellate(species=species)
            self.tuberculoventral = populations.Tuberculoventral(species=species)
            self.pyramidal = populations.Pyramidal(species=species)
            self.target = self.pyramidal

            self.sgc.connect(self.pyramidal, self.dstellate, self.tuberculoventral)
            self.dstellate.connect(self.pyramidal, self.tuberculoventral)  # Claude fixed 2026-06-25: DS→TV was missing
            self.tuberculoventral.connect(self.pyramidal)

            pops = [self.sgc, self.dstellate, self.tuberculoventral, self.pyramidal]
            self.pre_pops = [self.sgc, self.dstellate, self.tuberculoventral]
            self.recording_pops = [self.pyramidal, self.dstellate, self.tuberculoventral]

        elif celltype == "tuberculoventral":
            self.dstellate = populations.DStellate(species=species)
            self.tuberculoventral = populations.Tuberculoventral(species=species)
            self.target = self.tuberculoventral

            self.sgc.connect(self.tuberculoventral, self.dstellate)
            self.dstellate.connect(self.tuberculoventral)

            pops = [self.sgc, self.dstellate, self.tuberculoventral]
            self.pre_pops = [self.sgc, self.dstellate]
            self.recording_pops = [self.tuberculoventral, self.dstellate]

        else:
            raise ValueError(
                f"Unknown celltype {celltype!r}. Choose from: {CELLTYPES}"
            )

        self.populations = OrderedDict([(pop.type, pop) for pop in pops])

        # Instantiate target cells now so the dialog can inspect real cell objects.
        # resolve_inputs() is deferred to _resolve_network() so that synapse types
        # can be configured before any NEURON synapse objects are created.
        for f in frequencies:
            self.target.select(cells_per_band, cf=f, create=True)

    def _resolve_network(self):
        """Create NEURON synapses for all real target cells.

        Call after _build_topology() and after any _synapsetype_per_pre overrides
        have been set on post-populations.  depth=2 resolves inputs for the target
        cells (level 1) and their pre-synaptic cells (level 2).
        """
        self.target.resolve_inputs(depth=2)

    def custom_init(self, vinit=-60.):
        from cnmodel.util.pynrnutilities import custom_init as _custom_init
        _custom_init(v_init=vinit, dur=10.)  # 10 ms warmup; 50 ms default is overkill

    def run(self, stim, seed):
        """Run the network simulation with *stim* as the sound source and a unique
        *seed* used to configure the random number generators.
        """
        self.reset()

        # Generate 2 new seeds for the SGC spike generator and for the NEURON simulation
        rs = np.random.RandomState()
        rs.seed(self.seed ^ seed)
        seed1, seed2 = rs.randint(0, 2 ** 32, 2)
        random_seed.set_seed(seed1)
        self.sgc.set_seed(seed2)

        self.sgc.set_sound_stim(stim, parallel=False)

        # set up recording vectors
        for pop in self.recording_pops:
            for ind in pop.real_cells():
                cell = pop.get_cell(ind)
                self[cell] = cell.soma(0.5)._ref_v
        self["t"] = h._ref_t

        _tstop = stim.duration * 1000
        h.celsius = self.temp
        h.dt = self.dt

        # Claude fixed 2026-06-26: recording setup moved before custom_init so
        # h.finitialize() (inside custom_init) resets these vectors at the same
        # moment as t/Vm vectors, making all array lengths identical for plotting.
        # Duck-typing on psd attributes avoids importing PSD subclasses here.
        _syn_vecs: dict = {}  # (pre_type, pre_ind) -> list[h.Vector]
        if self.record_syn_conductance:
            for _tind in self.target.real_cells():
                _tcell = self.target.get_cell(_tind)
                for _inp in _tcell.inputs:
                    _syn = _inp[0]
                    _pre_cell = _syn.terminal.cell
                    _pre_type = _pre_cell.celltype
                    if _pre_type == "sgc":
                        continue
                    _pre_pop = self.populations.get(_pre_type)
                    if _pre_pop is None:
                        continue
                    _pre_ind = int(_pre_pop.get_cell_index(_pre_cell))
                    _key = (_pre_type, _pre_ind)
                    _psd = _syn.psd
                    _vecs: list = []
                    if hasattr(_psd, "all_psd"):          # GlyPSD
                        for _m in _psd.all_psd:
                            _v = h.Vector()
                            _v.record(_m._ref_g)
                            _vecs.append(_v)
                    elif hasattr(_psd, "ampa_psd"):       # GluPSD
                        for _m in _psd.ampa_psd:
                            _v = h.Vector()
                            _v.record(_m._ref_g)
                            _vecs.append(_v)
                    elif hasattr(_psd, "syn"):             # Exp2PSD
                        _v = h.Vector()
                        _v.record(_psd.syn._ref_g)
                        _vecs.append(_v)
                    if _vecs:
                        _syn_vecs[_key] = _vecs

        # set_maxstep runs after finitialize (correct NEURON order); warmup fadvance is also inside.
        self.custom_init()

        h.batch_run(_tstop, h.dt)

        # collect conductance traces (sum release zones → one array per pre-cell)
        vec_syn_g: dict = {}
        if self.record_syn_conductance:
            for _key, _vecs in _syn_vecs.items():
                # Claude fixed 2026-06-26: h.Vector must be converted via .to_python() before np.array()
                vec_syn_g[_key] = sum(np.array(_v.to_python()) for _v in _vecs)

        # record vsoma and spike times for all cells
        vec = {}
        for k in self._vectors:
            v = self[k].copy()
            if k == "t":
                vec[k] = v
                continue
            spike_inds = np.argwhere((v[1:] > -20) & (v[:-1] <= -20))[:, 0]
            spikes = self["t"][spike_inds]
            pop = k.celltype
            # print('pop: ', pop)
            assert isinstance(pop, str)
            cell_ind = getattr(self, pop).get_cell_index(k)
            vec[(pop, cell_ind)] = [v, spikes]

        # record SGC spike trains
        for ind in self.sgc.real_cells():
            cell = self.sgc.get_cell(ind)
            vec[("sgc", ind)] = [None, cell._spiketrain]

        vec["syn_g"] = vec_syn_g   # empty dict when record_syn_conductance is False
        return vec


def _stim_freq(stim) -> float:
    """Return the x-axis frequency for any stim type (tone f0 or noise center_freq)."""
    opts = stim.opts
    return float(opts.get("f0") or opts.get("center_freq") or 0.0)


class NetworkSimDisplay(pg.QtWidgets.QSplitter):
    def __init__(self, prot, results, baseline, response):
        pg.QtWidgets.QSplitter.__init__(self, pg.QtCore.Qt.Orientation.Horizontal)
        self.selected_cell = None

        self.prot = prot
        self.baseline = baseline  # (start, stop)
        self.response = response  # (start, stop)

        self.ctrl = pg.QtWidgets.QWidget()
        self.layout = pg.QtWidgets.QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.ctrl.setLayout(self.layout)
        self.addWidget(self.ctrl)

        self.nv = NetworkVisualizer(prot.populations)
        self.layout.addWidget(self.nv)
        self.nv.cell_selected.connect(self.nv_cell_selected)

        self.stim_type_combo = pg.QtWidgets.QComboBox()
        self.layout.addWidget(self.stim_type_combo)
        self.stim_combo = pg.QtWidgets.QComboBox()
        self.layout.addWidget(self.stim_combo)
        self.trial_combo = pg.QtWidgets.QComboBox()
        self.layout.addWidget(self.trial_combo)

        # Group incoming results by stimulus type.
        # Raw key structure: (stim_name, f0, dbspl, iteration)
        self.all_results: Dict[str, OrderedDict] = {}
        for k, v in list(results.items()):
            sname, f0, dbspl, iteration = k
            stim, result = v
            stype_res = self.all_results.setdefault(sname, OrderedDict())
            label = "f0: %0.0f  dBspl: %0.0f" % (f0, dbspl)
            stype_res.setdefault(label, [stim, {}])
            stype_res[label][1][iteration] = result

        for sname in self.all_results:
            display = STIM_DEFS[sname].label if sname in STIM_DEFS else sname
            self.stim_type_combo.addItem(display, sname)

        # Placeholder attributes — filled by stim_type_selected() below.
        self.results = OrderedDict()
        self.stim_order: list = []
        self.freqs: list = []
        self.levels: list = []
        self.iterations = 1
        self._df = 0.25
        self._dl = 10.0

        self.stim_type_combo.currentIndexChanged.connect(self.stim_type_selected)
        self.stim_combo.currentIndexChanged.connect(self.stim_selected)
        self.trial_combo.currentIndexChanged.connect(self.trial_selected)

        self.tuning_plot = pg.PlotWidget()
        self.tuning_plot.setLogMode(x=True, y=False)
        self.tuning_plot.setLabel('bottom', 'Frequency (Hz)')
        self.tuning_plot.setLabel('left', 'Level (dB SPL)')
        self.tuning_plot.scene().sigMouseClicked.connect(self.tuning_plot_clicked)
        self.layout.addWidget(self.tuning_plot)

        self.tuning_img = pg.ImageItem()
        self.tuning_plot.addItem(self.tuning_img)
        # Claude fixed 2026-06-25: diverging colormap — cool→black→hot centered at 0 (spontaneous)
        self.tuning_cmap = pg.ColorMap(
            pos=np.array([0.0, 0.3, 0.5, 0.7, 1.0]),
            color=np.array([
                [  0,   0, 160, 255],   # deep blue  (max decrease below spont)
                [  0, 200, 220, 255],   # cyan        (moderate decrease)
                [  0,   0,   0, 255],   # black       (no change = spontaneous)
                [220,  80,   0, 255],   # orange      (moderate increase)
                [255, 240,   0, 255],   # yellow      (max increase above spont)
            ], dtype=np.uint8),
        )
        self.tuning_img.setColorMap(self.tuning_cmap)
        self.tuning_cbar = pg.ColorBarItem(
            interactive=False,
            colorMap=self.tuning_cmap,
            label="Net rate (sp/s)  [black = spont]",
        )
        self.tuning_cbar.setImageItem(self.tuning_img, insert_in=self.tuning_plot.plotItem)

        # Rect sized to one pixel; size updated by stim_type_selected().
        self.stim_rect = pg.QtWidgets.QGraphicsRectItem(
            pg.QtCore.QRectF(0, 0, self._df, self._dl)
        )
        self.stim_rect.setPen(pg.mkPen("c"))
        self.stim_rect.setZValue(20)
        self.tuning_plot.addItem(self.stim_rect)
        self._grid_items: list = []  # faint frequency/level grid; rebuilt by update_tuning()

        # self.network_tree = NetworkTree(self.prot)
        # self.layout.addWidget(self.network_tree)

        self.pw = pg.GraphicsLayoutWidget()
        self.addWidget(self.pw)

        self.stim_plot = self.pw.addPlot(row=0, col=0)
        self.pw.ci.layout.setRowFixedHeight(0, 100)
        self.cell_plot = self.pw.addPlot(row=1, col=0, labels={"left": "Vm"})
        self.input_plot = self.pw.addPlot(
            row=2, col=0,
            labels={"left": "input #", "bottom": "time"}, title="Input spike times",
        )
        self.input_plot.setXLink(self.cell_plot)
        self.stim_plot.setXLink(self.cell_plot)
        # Claude changed 2026-06-26: rowspan=2 so ri_plot and g_plot share vertical space evenly
        self.ri_plot = self.pw.addPlot(
            row=0, col=1, rowspan=2,
            labels={"left": "Rate (sp/s)", "bottom": "Level (dB SPL)"},
            title="Rate-Intensity",
        )
        self.g_plot = self.pw.addPlot(
            row=2, col=1, rowspan=1,
            labels={"left": "Est. g (norm.)", "bottom": "Time (ms)"},
            title="Estimated input conductance",
        )
        self.g_plot.setXLink(self.cell_plot)
        self.pw.ci.layout.setColumnFixedWidth(1, 250)

        self.stim_type_selected()

    def update_stim_plot(self):
        stim = self.selected_stim
        self.stim_plot.plot(stim.time * 1000, stim.sound, clear=True, antialias=True)

    def update_raster_plot(self):
        self.input_plot.clear()
        if self.selected_cell is None:
            return
        pop, ind = self.selected_cell

        rec = pop._cells[ind]
        i = 0
        plots = []
        # plot spike times for all presynaptic cells
        labels = []
        if rec["connections"] == 0:
            return

        pop_colors = {
            "sgc": "g",
            "dstellate": "y",
            "tstellate": "b",
            "tuberculoventral": "r",
            "bushy": "w",
            "pyramidal": "m",
        }
        pop_order = self.prot.pre_pops
        trials = self.selected_trials()
        for pop in pop_order:
            pre_inds = rec["connections"].get(pop, [])
            for preind in pre_inds:
                # iterate over all trials
                for j in trials:
                    result = self.selected_results[j]
                    rkey = (pop.type, preind)
                    if rkey not in result:
                        continue
                    spikes = result[rkey][1]
                    y = np.ones(len(spikes)) * i + j / (
                        2.0 * len(self.selected_results)
                    )
                    self.input_plot.plot(
                        spikes,
                        y,
                        pen=None,
                        symbolBrush=pop_colors[pop.type],
                        symbol="+",
                        symbolPen=None,
                    )
                i += 1
                labels.append(pop.type + " " + str(preind))
        self.input_plot.getAxis("left").setTicks([list(enumerate(labels))])

    def update_cell_plot(self):
        self.cell_plot.clear()
        if self.selected_cell is None:
            return
        pop, cell_ind = self.selected_cell

        self.cell_plot.setTitle(
            "%s %d   %s" % (pop.type, cell_ind, str(self.stim_combo.currentText()))
        )
        key = (pop.type, cell_ind)
        trials = self.selected_trials()
        for i in trials:
            result = self.selected_results[i]
            if key not in result:
                continue
            y = result[key][0]
            if y is not None:
                p = self.cell_plot.plot(
                    self.selected_results[0]["t"],
                    y,
                    name="%s-%d" % self.selected_cell,
                    antialias=True,
                    pen=(i, len(self.selected_results) * 1.5),
                )
        # p.curve.setClickable(True)
        # p.sigClicked.connect(self.cell_curve_clicked)
        # p.cell_ind = ind

    def tuning_plot_clicked(self, event):
        spos = event.scenePos()
        stimpos = self.tuning_plot.plotItem.vb.mapSceneToView(spos)
        # stimpos is already in log10(Hz) / dB space (tuning_plot has logMode x=True)
        log_x = stimpos.x()
        y = stimpos.y()

        # Find the stim point nearest the click in (log-freq, dB) space
        best = None
        best_dist = float("inf")
        for stim, result in list(self.results.values()):
            f0 = _stim_freq(stim)
            dbspl = stim.opts["dbspl"]
            dist = (log_x - np.log10(f0)) ** 2 + (y - dbspl) ** 2
            if dist < best_dist:
                best = stim
                best_dist = dist

        if best is None:
            return
        self.select_stim(_stim_freq(best), best.opts["dbspl"])

    def nv_cell_selected(self, nv, cell):
        self.select_cell(*cell)

    def stim_selected(self):
        key = str(self.stim_combo.currentText())
        if not key or key not in self.results:
            return
        results = self.results[key]
        self.selected_results = results[1]
        self.selected_stim = results[0]
        self.update_stim_plot()
        self.update_raster_plot()
        self.update_cell_plot()
        self.update_ri_plot()
        self.update_g_plot()

        # Centre the rect on the stim point (image pixels are also centred on stim points)
        self.stim_rect.setPos(
            np.log10(_stim_freq(results[0])) - self._df / 2,
            results[0].opts["dbspl"] - self._dl / 2,
        )

    def trial_selected(self):
        self.update_raster_plot()
        self.update_cell_plot()
        self.update_tuning()
        self.update_g_plot()

    def stim_type_selected(self):
        """Switch to the currently selected stimulus type and rebuild dependent widgets."""
        sname = self.stim_type_combo.currentData()
        if sname is None or sname not in self.all_results:
            return

        self.results = self.all_results[sname]

        freqs: set = set()
        levels: set = set()
        max_iter = 0
        for _, (stim, iters) in self.results.items():
            freqs.add(_stim_freq(stim))
            levels.add(stim.opts["dbspl"])
            for it in iters:
                max_iter = max(max_iter, it)

        self.freqs = sorted(freqs)
        self.levels = sorted(levels)
        self.iterations = max_iter + 1
        self._df = (
            (np.log10(self.freqs[1]) - np.log10(self.freqs[0]))
            if len(self.freqs) > 1 else 0.15  # matches single-freq pixel width in update_tuning
        )
        self._dl = (self.levels[1] - self.levels[0]) if len(self.levels) > 1 else 10.0
        self.stim_order = [
            (_stim_freq(stim), stim.opts["dbspl"])
            for _, (stim, _) in self.results.items()
        ]

        # Rebuild stim_combo without triggering stim_selected mid-update.
        self.stim_combo.blockSignals(True)
        self.stim_combo.clear()
        for label in self.results:
            self.stim_combo.addItem(label)
        self.stim_combo.blockSignals(False)

        # Rebuild trial_combo.
        self.trial_combo.blockSignals(True)
        self.trial_combo.clear()
        self.trial_combo.addItem("all trials")
        for i in range(self.iterations):
            self.trial_combo.addItem(str(i))
        self.trial_combo.blockSignals(False)

        self.stim_rect.setRect(pg.QtCore.QRectF(0, 0, self._df, self._dl))
        self.stim_selected()

    def selected_trials(self):
        if self.trial_combo.currentIndex() == 0:
            return list(range(self.iterations))
        else:
            return [self.trial_combo.currentIndex() - 1]

    def select_stim(self, f0, dbspl):
        i = self.stim_order.index((f0, dbspl))
        self.stim_combo.setCurrentIndex(i)

    def select_cell(self, pop, cell_id):
        self.selected_cell = pop, cell_id
        self.update_tuning()
        self.update_cell_plot()
        self.update_raster_plot()
        self.update_g_plot()

    # def cell_curve_clicked(self, c):
    # if self.selected is not None:
    # pen = self.selected.curve.opts['pen']
    # pen.setWidth(1)
    # self.selected.setPen(pen)

    # pen = c.curve.opts['pen']
    # pen.setWidth(3)
    # c.setPen(pen)
    # self.selected = c

    # self.show_cell(c.cell_ind)

    def update_tuning(self):
        # update matrix image
        if self.selected_cell is None:
            return

        pop, ind = self.selected_cell
        fvals = set()
        lvals = set()

        # first get lists of all frequencies and levels in the matrix
        for stim, vec in list(self.results.values()):
            fvals.add(_stim_freq(stim))
            lvals.add(stim.key()["dbspl"])
        fvals = sorted(list(fvals))
        lvals = sorted(list(lvals))

        # Claude added 2026-06-26: faint grey grid at each (freq, level) point
        for item in self._grid_items:
            self.tuning_plot.removeItem(item)
        self._grid_items = []
        _grey_pen = pg.mkPen(color=(140, 140, 140, 160), width=1)
        _lv0 = lvals[0] - 5
        _lv1 = lvals[-1] + 5
        _lf0 = np.log10(fvals[0]) - 0.1
        _lf1 = np.log10(fvals[-1]) + 0.1
        for _f in fvals:
            _ln = pg.QtWidgets.QGraphicsLineItem(np.log10(_f), _lv0, np.log10(_f), _lv1)
            _ln.setPen(_grey_pen)
            _ln.setZValue(1)
            self.tuning_plot.addItem(_ln)
            self._grid_items.append(_ln)
        for _lv in lvals:
            _ln = pg.QtWidgets.QGraphicsLineItem(_lf0, _lv, _lf1, _lv)
            _ln.setPen(_grey_pen)
            _ln.setZValue(1)
            self.tuning_plot.addItem(_ln)
            self._grid_items.append(_ln)

        # Get spontaneous rate statistics
        key = (pop.type, ind)
        spont_spikes = 0
        spont_time = 0
        for stim, iterations in list(self.results.values()):
            for vec in list(iterations.values()):
                if key not in vec:
                    self.tuning_img.setImage(np.zeros((1, 1)))
                    return
                spikes = vec[key][1]
                spont_spikes += (
                    (spikes >= self.baseline[0]) & (spikes < self.baseline[1])
                ).sum()
                spont_time += self.baseline[1] - self.baseline[0]
        spont_rate = spont_spikes / spont_time if spont_time > 0 else 0.0

        # next count the number of spikes for the selected cell at each point in the matrix
        matrix = np.zeros((len(fvals), len(lvals)))
        trials = self.selected_trials()
        for stim, iteration in list(self.results.values()):
            for i in trials:
                vec = iteration[i]
                if key not in vec:
                    continue
                spikes = vec[key][1]
                n_spikes = (
                    (spikes >= self.response[0]) & (spikes < self.response[1])
                ).sum()
                i = fvals.index(_stim_freq(stim))
                j = lvals.index(stim.key()["dbspl"])
                matrix[i, j] += n_spikes - spont_rate * (
                    self.response[1] - self.response[0]
                )
        matrix /= self.iterations

        # Convert net spikes/trial → net rate (sp/s) for display
        response_dur_s = (self.response[1] - self.response[0]) / 1000.0
        rate_matrix = matrix / response_dur_s

        # Scale and position the image so each pixel occupies its correct
        # (log10-freq, dB) region. sx / sy = size of one pixel in data units.
        # For single-frequency (noise) stimuli, use a fixed pixel width of 0.15 log10
        # decades (~15 % of a ±0.5-decade visible window) so the column is visible.
        sx = (np.log10(fvals[-1]) - np.log10(fvals[0])) / (len(fvals) - 1) if len(fvals) > 1 else 0.15
        sy = (lvals[-1] - lvals[0]) / (len(lvals) - 1) if len(lvals) > 1 else 1.0
        tr = pg.QtGui.QTransform()
        tr.scale(sx, sy)
        self.tuning_img.setImage(rate_matrix)
        self.tuning_img.setTransform(tr)
        # Shift by -sx/2, -sy/2 so pixel CENTRES align with stim data points,
        # matching the stim_rect centred positioning in stim_selected().
        self.tuning_img.setPos(np.log10(fvals[0]) - sx / 2, lvals[0] - sy / 2)
        # For single-frequency stimuli, centre the x-axis view on the CF.
        if len(fvals) == 1:
            self.tuning_plot.setXRange(
                np.log10(fvals[0]) - 0.5,
                np.log10(fvals[0]) + 0.5,
                padding=0,
            )
        # Claude fixed 2026-06-25: symmetric range so 0 (spontaneous) is always the midpoint
        # of the diverging colormap (maps to black).
        sym_max = max(float(abs(rate_matrix.min())), float(rate_matrix.max()), 1.0)
        self.tuning_cbar.setLevels(values=(-sym_max, sym_max))
        try:
            self.tuning_cbar.axis.setTicks([
                [(-sym_max, f'{-sym_max:.0f}'), (0, '0\n(spont)'), (sym_max, f'+{sym_max:.0f}')],
            ])
        except Exception:
            pass
        self.update_ri_plot()

    def update_ri_plot(self):
        """Rate-intensity function for the selected cell.

        For tone pip: uses the currently selected stimulus frequency.
        For noise types: uses the nearest test frequency to the cell's actual CF.
        Spontaneous rate is estimated from the pre-stimulus baseline across ALL trials
        and ALL stimulus conditions (maximum statistical power).
        """
        self.ri_plot.clear()
        if self.selected_cell is None:
            return

        pop, ind = self.selected_cell
        key = (pop.type, ind)
        trials = self.selected_trials()
        all_trials = list(range(self.iterations))

        cell_cf = float(pop._cells[ind]['cf'])

        fvals = sorted(set(_stim_freq(s) for s, _ in self.results.values()))
        lvals = sorted(set(s.key()["dbspl"] for s, _ in self.results.values()))
        if not fvals or not lvals:
            return

        # Determine which test frequency to use for the RI plot.
        sname = self.stim_type_combo.currentData()
        if (sname == "tone_pip"
                and getattr(self, 'selected_stim', None) is not None):
            cf_test = _stim_freq(self.selected_stim)
        else:
            cf_arr = np.array(fvals)
            cf_idx = int(np.argmin(np.abs(cf_arr - cell_cf)))
            cf_test = fvals[cf_idx]

        response_dur_s = (self.response[1] - self.response[0]) / 1000.0
        baseline_dur_s = (self.baseline[1] - self.baseline[0]) / 1000.0

        # Spontaneous rate: pool baseline spikes across ALL results and ALL trials
        # for the most stable estimate (pre-stimulus rate is stimulus-independent).
        spont_n = spont_tr = 0
        for _, (stim, iteration) in self.results.items():
            for i in all_trials:
                if i not in iteration:
                    continue
                vec = iteration[i]
                if key not in vec:
                    continue
                spikes = vec[key][1]
                spont_n += ((spikes >= self.baseline[0]) & (spikes < self.baseline[1])).sum()
                spont_tr += 1
        mean_spont = (spont_n / spont_tr / baseline_dur_s) if spont_tr > 0 else 0.0

        # Response rate vs level at cf_test.
        ri_rates = []
        valid_levels = []
        for lv in lvals:
            res_key = "f0: %0.0f  dBspl: %0.0f" % (cf_test, lv)
            if res_key not in self.results:
                continue
            stim, iteration = self.results[res_key]
            n_resp = n_tr = 0
            for i in trials:
                if i not in iteration:
                    continue
                vec = iteration[i]
                if key not in vec:
                    continue
                spikes = vec[key][1]
                n_resp += ((spikes >= self.response[0]) & (spikes < self.response[1])).sum()
                n_tr += 1
            if n_tr > 0:
                ri_rates.append(n_resp / n_tr / response_dur_s)
                valid_levels.append(lv)

        if not valid_levels:
            return

        self.ri_plot.plot(
            valid_levels, ri_rates,
            pen=pg.mkPen('w', width=2),
            symbol='o', symbolBrush='w', symbolSize=6,
        )
        if mean_spont > 0:
            self.ri_plot.addLine(
                y=mean_spont,
                pen=pg.mkPen('y', style=pg.QtCore.Qt.PenStyle.DashLine, width=1),
            )

        ymax = max(max(ri_rates) if ri_rates else 0.0, mean_spont) * 1.15
        self.ri_plot.setYRange(0, max(ymax, 1.0), padding=0)

        if sname == "tone_pip":
            self.ri_plot.setTitle(
                "%s %d   f0: %.0f Hz" % (pop.type, ind, cf_test),
                size="9pt",
            )
        else:
            self.ri_plot.setTitle(
                "%s %d   CF: %.0f Hz" % (pop.type, ind, cell_cf),
                size="9pt",
            )

    # Claude added 2026-06-26: actual synaptic conductance recorded during simulation
    def update_g_plot(self):
        """Plot recorded synaptic conductances for non-SGC inputs to the selected cell.

        Conductance data is only present when 'Record synaptic conductances' was
        checked in the control panel before the run.  When absent, a prompt is
        shown instead of traces.
        """
        self.g_plot.clear()
        if not self.selected_results:
            return

        # Check whether conductance data was collected this run
        _sample = self.selected_results.get(0) or next(iter(self.selected_results.values()))
        _syn_g_sample = _sample.get("syn_g", {})
        if not _syn_g_sample:
            self.g_plot.setTitle(
                "Synaptic conductance  [enable 'Record synaptic conductances' in control panel]",
                size="8pt",
            )
            return

        if self.selected_cell is None:
            self.g_plot.setTitle(
                "Synaptic conductance  [click the target cell in the network diagram]",
                size="8pt",
            )
            return

        pop, ind = self.selected_cell
        rec = pop._cells[ind]
        if rec["connections"] == 0:
            return

        pop_colors = {
            "sgc": "g", "dstellate": "y", "tstellate": "b",
            "tuberculoventral": "r", "bushy": "w", "pyramidal": "m",
        }

        t = np.array(self.selected_results[0]["t"])
        trials = self.selected_trials()

        for pre_pop in self.prot.pre_pops:
            if pre_pop.type == "sgc":
                continue
            pre_inds = rec["connections"].get(pre_pop, [])
            # Claude fixed 2026-06-26: numpy array([0]) is falsy — must check length explicitly
            if len(pre_inds) == 0:
                continue
            color = pop_colors.get(pre_pop.type, "w")

            g_total = np.zeros(len(t))
            n_contrib = 0
            for preind in pre_inds:
                key = (pre_pop.type, int(preind))
                for j in trials:
                    g_data = self.selected_results[j].get("syn_g", {})
                    if key not in g_data:
                        continue
                    g_arr = g_data[key]
                    if len(g_arr) == len(t):
                        g_total += g_arr
                        n_contrib += 1

            if n_contrib == 0:
                continue
            g_avg = g_total / n_contrib
            peak = float(np.abs(g_avg).max())
            if peak > 0:
                g_avg = g_avg / peak
            self.g_plot.plot(
                t, g_avg,
                pen=pg.mkPen(color, width=1.5),
                name=pre_pop.type,
            )
        self.g_plot.setTitle("Synaptic conductance (normalised)", size="9pt")


class NetworkTree(pg.QtWidgets.QTreeWidget):
    def __init__(self, prot):
        self.prot = prot
        QtGui.QTreeWidget.__init__(self)
        self.setColumnCount(2)

        self.update_tree()

    def update_tree(self):
        for pop_name in ["bushy", "tstellate", "dstellate", "tuberculoventral", "sgc"]:
            if not hasattr(self.prot, pop_name):
                continue
            pop = getattr(self.prot, pop_name)
            grp = QtGui.QTreeWidgetItem([pop_name])
            self.addTopLevelItem(grp)
            for cell in pop.real_cells():
                self.add_cell(grp, pop, cell)

    def add_cell(self, grp_item, pop, cell):
        item = QtGui.QTreeWidgetItem([str(cell)])
        grp_item.addChild(item)
        all_conns = pop.cell_connections(cell)
        if all_conns == 0:
            return
        for cpop, conns in list(all_conns.items()):
            pop_grp = QtGui.QTreeWidgetItem([cpop.celltype, str(conns)])
            item.addChild(pop_grp)


class NetworkVisualizer(pg.PlotWidget):

    cell_selected = pg.QtCore.Signal(object, object)

    def __init__(self, populations):
        self.pops = populations
        pg.PlotWidget.__init__(self)
        self.setLogMode(x=True, y=False)

        self.cells = pg.ScatterPlotItem(clickable=True)
        self.cells.setZValue(10)
        self.addItem(self.cells)
        self.cells.sigClicked.connect(self.cells_clicked)

        self.selected = pg.ScatterPlotItem()
        self.selected.setZValue(20)
        # pyqtgraph's GraphicsScene.sendClickEvent() calls mouseClickEvent() on
        # items by z-order regardless of Qt's acceptedMouseButtons. self.selected
        # sits on top of self.cells and its spots cover every highlighted cell
        # (including presynaptic cells shown in red). Without this override,
        # clicking any highlighted cell hits self.selected first, which accepts
        # the event, and self.cells never receives it.
        # Claude fixed 2026-06-24: override mouseClickEvent to always ignore.
        self.selected.mouseClickEvent = lambda ev: ev.ignore()
        self.addItem(self.selected)

        self.connections = pg.PlotCurveItem()
        self.connections.setAcceptedMouseButtons(pg.QtCore.Qt.MouseButton.NoButton)
        self.addItem(self.connections)

        # first assign positions of all cells
        pop_brushes = {
            "sgc": pg.mkBrush(80, 200, 80),
            "dstellate": pg.mkBrush(200, 200, 0),
            "tstellate": pg.mkBrush(80, 80, 220),
            "tuberculoventral": pg.mkBrush(200, 60, 60),
            "bushy": pg.mkBrush(0, 200, 220),
            "pyramidal": pg.mkBrush(180, 80, 220),
        }
        cells = []
        for y, pop in enumerate(self.pops.values()):
            pop.cell_spots = []
            pop.fwd_connections = {}
            brush = pop_brushes.get(pop.type, pg.mkBrush("w"))
            for i, cell in enumerate(pop._cells):
                pos = (np.log10(cell["cf"]), y)
                real = cell["cell"] != 0
                if not real:
                    pop.cell_spots.append(None)
                    continue
                spot = {
                    "x": pos[0],
                    "y": pos[1],
                    "size": 12,
                    "symbol": "o",
                    "brush": brush,
                    "pen": pg.mkPen("w", width=1),
                    "data": (pop, i),
                }
                cells.append(spot)
                pop.cell_spots.append(spot)

        self.cells.setData(cells)

        self.getAxis("left").setTicks([list(enumerate(self.pops.keys()))])

        # now assign connection lines and record forward connectivity
        con_x = []
        con_y = []
        for pop in list(self.pops.values()):
            for i, cell in enumerate(pop._cells):
                conns = cell["connections"]
                if conns == 0:
                    continue
                for prepop, precells in list(conns.items()):
                    spot = pop.cell_spots[i]
                    if spot is None:
                        continue
                    p1 = spot["x"], spot["y"]
                    for j in precells:
                        prepop.fwd_connections.setdefault(j, [])
                        prepop.fwd_connections[j].append((pop, i))
                        spot2 = prepop.cell_spots[j]
                        if spot2 is None:
                            continue  # Claude fixed 2026-06-24: was `return`, which exited __init__ early
                        p2 = spot2["x"], spot2["y"]
                        con_x.extend([p1[0], p2[0]])
                        con_y.extend([p1[1], p2[1]])
        self.connections.setData(
            x=con_x, y=con_y, connect="pairs", pen=(255, 255, 255, 60)
        )

    def cells_clicked(self, *args):
        selected = None
        for spot in args[1]:
            # find the first real cell
            pop, i = spot.data()
            if pop._cells[i]["cell"] != 0:
                selected = spot
                break
        if selected is None:
            self.selected.hide()
            return

        rec = pop._cells[i]
        pos = selected.pos()
        spots = [
            {
                "x": pos.x(),
                "y": pos.y(),
                "size": 15,
                "symbol": "o",
                "pen": "y",
                "brush": "b",
            }
        ]

        # display presynaptic cells
        if rec["connections"] != 0:
            for prepop, preinds in list(rec["connections"].items()):
                for preind in preinds:
                    # Claude fixed 2026-06-24: guard against None (virtual/uninstantiated cell)
                    s = prepop.cell_spots[preind]
                    if s is None:
                        continue
                    s = s.copy()
                    s["size"] = 15
                    s["brush"] = pg.mkBrush((255, 0, 0, 75))
                    spots.append(s)

        # display postsynaptic cells
        for postpop, postind in pop.fwd_connections.get(i, []):
            # Claude fixed 2026-06-24: guard against None (virtual/uninstantiated cell)
            s = postpop.cell_spots[postind]
            if s is None:
                continue
            s = s.copy()
            s["size"] = 15
            s["brush"] = pg.mkBrush((0, 255, 0, 75))
            spots.append(s)

        self.selected.setData(spots)
        self.selected.show()

        self.cell_selected.emit(self, selected.data())


def main():
    parser = argparse.ArgumentParser(description="CN network physiology simulation")
    parser.add_argument(
        "--celltype",
        choices=CELLTYPES,
        default="bushy",
        help="Target cell type (default: bushy)",
    )
    parser.add_argument(
        "--species",
        choices=["mouse", "guineapig", "rat"],
        default="mouse",
        help="Species for cell parameters and connectivity (default: mouse)",
    )
    parser.add_argument(
        "--simulator",
        choices=["py3", "matlab"],
        default=None,
        help="AN spike-train simulator (default: auto-detect)",
    )
    parser.add_argument(
        "--ignore-cache",
        action="store_true",
        help="Ignore cached results and rerun simulations",
    )
    parser.add_argument("--fmin", type=float, default=4000.,
                        help="Minimum frequency (Hz, default 4000)")
    parser.add_argument("--fmax", type=float, default=32000.,
                        help="Maximum frequency (Hz, default 32000)")
    parser.add_argument("--octave-spacing", type=float, default=0.25,
                        help="Octave spacing between test frequencies (default 0.25)")
    parser.add_argument("--db-min", type=float, default=20.,
                        help="Minimum sound level (dB SPL, default 20)")
    parser.add_argument("--db-max", type=float, default=100.,
                        help="Maximum sound level (dB SPL, default 100)")
    parser.add_argument("--db-step", type=float, default=10.,
                        help="Sound level step size (dB, default 10)")
    parser.add_argument("--execution", choices=["parallel", "serial"], default="parallel",
                        help="Run simulations in parallel or serial (default: parallel)")
    parser.add_argument("--cf", type=float, default=16000.,
                        help="Characteristic frequency of the target cell (Hz, default 16000)")
    parser.add_argument(
        "--stim-types",
        nargs="+",
        choices=list(STIM_DEFS.keys()),
        default=["tone_pip"],
        metavar="TYPE",
        help="One or more stimulus types to run: %s (default: tone_pip)"
             % ", ".join(STIM_DEFS.keys()),
    )
    parser.add_argument(
        "--noise-octave-width",
        type=float,
        default=1.0,
        help="Bandwidth for band-limited noise in octaves (default: 1.0)",
    )
    args = parser.parse_args()

    pg.mkQApp()

    parallel = (args.execution == "parallel")

    nreps = 1
    n_frequencies = int(np.log2(args.fmax / args.fmin) / args.octave_spacing) + 1
    fvals = (
        np.logspace(
            np.log2(args.fmin / 1000.0),
            np.log2(args.fmax / 1000.0),
            num=n_frequencies,
            endpoint=True,
            base=2,
        )
        * 1000.0
    )

    n_levels = int(round((args.db_max - args.db_min) / args.db_step)) + 1
    levels = np.linspace(args.db_min, args.db_max, n_levels)

    print(("Frequencies (kHz):", fvals / 1000.0))
    print(("Levels (dB SPL):", levels))

    syntype = "multisite"
    path = os.path.dirname(__file__)
    cachepath = os.path.join(path, "cache")
    if not os.path.isdir(cachepath):
        os.mkdir(cachepath)

    seed = 34657845
    prot = CNSoundStim(seed=seed, synapsetype=syntype, celltype=args.celltype,
                       cf=args.cf, species=args.species)

    # Show connection configuration dialog before running any simulations.
    # The user can inspect the convergence table and scale/silence individual
    # connection types.  Changes are applied live to the network's synapse
    # gmax values so they take effect when prot.run() is called below.
    net_dlg = NetworkControlPanel()
    net_dlg.update_table(args.species)
    net_dlg.set_synapse_types(prot)  # populate Synapse Types tab; Build Network fires _resolve_network internally
    if net_dlg.exec() != pg.QtWidgets.QDialog.DialogCode.Accepted:
        print("Simulation cancelled.")
        return
    _conn_modified = net_dlg.any_modified()

    # Claude added 2026-06-26: propagate conductance-recording flag to protocol
    prot.record_syn_conductance = net_dlg.get_record_syn_conductance()

    # Claude added 2026-06-26: zero Na channels on target cells when requested
    if net_dlg.get_na_zero():
        _NA_MECHS = ["nav11", "na", "nabu", "jsrna", "nacncoop", "nacn"]
        _target_pop = prot.populations.get(args.celltype)
        if _target_pop is not None:
            for _cid in _target_pop.real_cells():
                _cobj = _target_pop._cells[_cid]["cell"]
                if _cobj == 0:
                    continue
                for _sec_list in _cobj.all_sections.values():
                    for _sec in _sec_list:
                        for _seg in _sec.allseg():
                            for _mech in _NA_MECHS:
                                if hasattr(_seg, _mech):
                                    getattr(_seg, _mech).gbar = 0.0

    start_time = timeit.default_timer()

    stimpar = {
        # "dur": 0.35,        # Claude fixed 2026-06-25: was 350 ms (100 pre + 100 pip + 150 post)
        "dur": 0.25,        # 250 ms total: 100 ms pre + 100 ms pip + 50 ms post
        "pip": 0.1,         # 100 ms pip
        "start": [0.1],     # pip starts at 100 ms
        "baseline": [50, 100],   # 50 ms pre-stimulus window (ms)
        "response": [100, 200],  # 100 ms response window matching pip duration (ms)
    }

    # Noise stimuli run only at the target CF (no frequency sweep).
    _stim_freqs: Dict[str, np.ndarray] = {
        "tone_pip":          fvals,
        "bandlimited_noise": np.array([args.cf]),
        "wideband_noise":    np.array([args.cf]),
    }
    # Type-specific extra parameters forwarded to each factory.
    _stim_extra: Dict[str, Dict[str, Any]] = {
        "tone_pip":          {},
        "bandlimited_noise": {"octave_width": args.noise_octave_width},
        "wideband_noise":    {},
    }

    tasks = []
    for sname in args.stim_types:
        for f in _stim_freqs[sname]:
            for db in levels:
                for i in range(nreps):
                    tasks.append((sname, float(f), db, i))

    results = {}

    workers = 1 if not parallel else mp.Parallelize.suggestedWorkerCount()
    tot_runs = len(tasks)
    with mp.Parallelize(
        enumerate(tasks),
        results=results,
        progressDialog="Running parallel simulation.. (may take a while)",
        workers=workers,
    ) as tasker:
        for i, task in tasker:
            sname, f, db, iteration = task
            sdef = STIM_DEFS[sname]
            extra = _stim_extra[sname]
            stim = sdef.factory(stimpar, f, db, iteration, **extra)

            # Stable extra-param suffix so different configurations get separate caches.
            extra_str = "".join("_%s=%s" % (k, v) for k, v in sorted(extra.items()))
            print(f"=== Start run {i+1:5d}/{tot_runs:5d} ({sname}) === ", end="")
            cachefile = os.path.join(
                cachepath,
                "seed=%d_stim=%s_f0=%f_dbspl=%f_syntype=%s_celltype=%s%s_iter=%d.pk"
                % (seed, sname, f, db, syntype, args.celltype, extra_str, iteration),
            )
            if args.ignore_cache or _conn_modified or not os.path.isfile(cachefile):
                result = prot.run(stim, seed=i)
                if not _conn_modified:  # don't cache non-default configurations
                    pickle.dump(result, open(cachefile, "wb"))
            else:
                print("  (Loading cached results)", end="")
                result = pickle.load(open(cachefile, "rb"))
            tasker.results[(sname, f, db, iteration)] = (stim, result)
            print(f"  --- finished run {i+1:5d}/{tot_runs:5d} ---")

    elapsed = timeit.default_timer() - start_time
    print(
        "Elapsed time for %d stimuli: %f  (%f sec per stim), celltype: %s, synapses: %s"
        % (len(tasks), elapsed, elapsed / len(tasks), args.celltype, prot.sgc._synapsetype)
    )

    nd = NetworkSimDisplay(
        prot, results, baseline=stimpar["baseline"], response=stimpar["response"]
    )
    nd.show()

    if sys.flags.interactive == 0:
        pg.QtWidgets.QApplication.exec()


if __name__ == "__main__":
    main()