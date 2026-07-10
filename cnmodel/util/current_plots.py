"""
Per-mechanism current/conductance traces for cnmodel cells.

CurrentPlotsWindow(cell, iamp_nA, durations) opens a pyqtgraph window
with Vm (row 0, one column per compartment) and per-mechanism ionic current
or conductance traces (rows 1..N), recorded from the first (most-proximal)
section of each compartment type present in cell.all_sections.

A toolbar button toggles between current (mA/cm²) and conductance (mho/cm²)
views without re-running the simulation.  Mechanisms with zero conductance
show a 'g = 0' label instead of a trace.  All x-axes are linked.
"""
from neuron import h
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets

# ── NEURON RANGE variable for each mechanism's ionic current ────────────────
MECH_CURRENT_VAR = {
    'GRC_NA':   'ina',
    'GRC_KV':   'ik',
    'GRC_KA':   'ik',
    'GRC_KM':   'ik',
    'GRC_KIR':  'ik',
    'GRC_KCA':  'ik',
    'GRC_CA':   'ica',
    'GRC_LKG1': 'i',
    'GRC_LKG2': 'i',
    'nacn':     'ina',
    'na':       'ina',
    'nav11':    'ina',
    'nacncoop': 'ina',
    'nabu':     'ina',
    'kht':      'ik',
    'klt':      'ik',
    'ka':       'ik',
    'kd':       'ik',
    'ihvcn':    'ih',
    'hcno':     'ih',
    'hcnobo':   'ih',
    'leak':     'i',
    'pas':      'i',
}

# ── NEURON RANGE variable for each mechanism's time-varying conductance ─────
MECH_COND_RANGE_VAR = {
    'GRC_NA':   'gna',
    'GRC_KV':   'gk',
    'GRC_KA':   'gk',
    'GRC_KM':   'gk',
    'GRC_KIR':  'gk',
    'GRC_KCA':  'gk',
    'GRC_CA':   'gca',
    'GRC_LKG1': 'gl',     # parameter — constant flat line
    'GRC_LKG2': 'ggaba',  # parameter — constant flat line
    'nacn':     'gna',
    'na':       'gna',
    'nav11':    'gna',
    'nacncoop': 'gna',
    'nabu':     'gna',
    'kht':      'gk',
    'klt':      'gk',
    'ka':       'gk',
    'ihvcn':    'g',
    'hcno':     'g',
    'hcnobo':   'g',
    'leak':     'g',
    'pas':      'g',
}

# ── gbar variable name used to test whether a mechanism is truly active ─────
MECH_COND_VAR = {
    'GRC_LKG1': 'gl',
    'GRC_LKG2': 'ggaba',
}

# Mechanisms with no plottable ionic current (skip entirely)
SKIP_MECHS = {'GRC_CALC', 'cadiff', 'cap', 'extracellular'}

# Cycling color palette for compartment columns (RGB tuples, 0-255).
_COLUMN_COLORS = [
    (86,  180, 233),   # sky blue
    (230, 159,   0),   # orange
    (  0, 158, 115),   # teal/green
    (213,  94,   0),   # vermillion
    (204, 121, 167),   # pink
    (240, 228,  66),   # yellow
    (  0, 114, 178),   # blue
    (150, 150, 150),   # grey
]


# ── helpers ─────────────────────────────────────────────────────────────────

def _cond_var(mech):
    return MECH_COND_VAR.get(mech, 'gbar')


def _gbar_nonzero(sec, mech):
    """True if mech is inserted in sec with a non-zero conductance parameter."""
    try:
        g = getattr(getattr(sec(0.5), mech), _cond_var(mech))
        return float(g) > 0.0
    except (AttributeError, RuntimeError):
        return False


def _try_ref(seg, var, mech):
    """Return `_ref_<var>_<mech>` on seg, or None if unavailable."""
    try:
        return getattr(seg, f'_ref_{var}_{mech}')
    except (AttributeError, RuntimeError):
        return None


def _find_var(seg, mech, preferred, fallbacks):
    """
    Try `preferred` first, then each fallback, returning the first name that
    produces a valid `_ref_` on seg, or None.
    """
    for v in ([preferred] if preferred else []) + fallbacks:
        if _try_ref(seg, v, mech) is not None:
            return v
    return None


def _mechs_in_sec(sec):
    return list(sec.psection().get('density_mechs', {}).keys())


# ── main window ─────────────────────────────────────────────────────────────

class CurrentPlotsWindow(QtWidgets.QMainWindow):
    """
    Grid window: Vm row (one per compartment) + per-mechanism rows.
    Columns = compartment types (first proximal section each).
    A toolbar button toggles between current (mA/cm²) and conductance (mho/cm²).
    """

    def __init__(self, cell, iamp_nA=0.020, durations=(10, 100, 50)):
        super().__init__()
        self.setWindowTitle('Current Plots')
        self.cell = cell
        self.iamp_nA = iamp_nA
        self.durations = list(durations)
        self._vecs = []           # keep all h.Vector refs alive
        self._mode = 'current'    # 'current' | 'conductance'
        self._t_arr = None
        self._vm_data = {}        # comp_name -> np.ndarray (Vm per compartment)
        self._cur_data = {}       # (comp_name, mech) -> np.ndarray | None
        self._cond_data = {}      # (comp_name, mech) -> np.ndarray | None
        self._curves = {}         # (comp_name, mech) -> PlotDataItem | None
        self._plots = {}          # (comp_name, mech) -> PlotItem
        self._row_y_roots = {}    # mech -> PlotItem (y-axis root for each row)
        self._row_plots = {}      # mech -> [PlotItem, ...] (all live plots per row)
        self._vm_plots = {}       # comp_name -> PlotItem for Vm row (plots with data)
        self._vm_y_root = None    # y-axis root for the Vm row
        self._ref_p = None        # x-axis root for the whole grid
        self._all_mechs = []      # ordered mechanism list (set in _build_layout)
        self._run_and_plot()

    # ── compartment discovery ────────────────────────────────────────────────

    def _comp_sections(self):
        """
        Return {comp_name: section} — one (first/proximal) section per
        compartment type, or {'soma': cell.soma} for point models.
        """
        all_sec = getattr(self.cell, 'all_sections', None)
        if all_sec:
            return {name: secs[0] for name, secs in all_sec.items() if secs}
        return {'soma': self.cell.soma}

    # ── simulation & recording ───────────────────────────────────────────────

    def _run_and_plot(self):
        cell = self.cell
        tdel, tdur, tpost = self.durations
        tstop = tdel + tdur + tpost
        soma = cell.soma

        comp_secs = self._comp_sections()

        # Collect all unique mechanisms across compartment sections (order preserved).
        all_mechs = []
        seen = set()
        for sec in comp_secs.values():
            for m in _mechs_in_sec(sec):
                if m not in SKIP_MECHS and m not in seen:
                    all_mechs.append(m)
                    seen.add(m)

        # Discover recordable RANGE variables per (comp, mech).
        cur_var = {}    # (comp, mech) -> str | None
        cond_var = {}   # (comp, mech) -> str | None
        has_mech = {}   # (comp, mech) -> bool

        for comp_name, sec in comp_secs.items():
            sec_mechs = set(_mechs_in_sec(sec))
            seg = sec(0.5)
            for mech in all_mechs:
                key = (comp_name, mech)
                if mech not in sec_mechs:
                    has_mech[key] = False
                    continue
                has_mech[key] = True
                cur_var[key] = _find_var(
                    seg, mech,
                    MECH_CURRENT_VAR.get(mech),
                    ['ina', 'ik', 'ica', 'ih', 'i'],
                )
                cond_var[key] = _find_var(
                    seg, mech,
                    MECH_COND_RANGE_VAR.get(mech),
                    ['gna', 'gk', 'gca', 'g', 'gl', 'ggaba'],
                )

        # Time and Vm recordings — one Vm vector per compartment section.
        t_vec = h.Vector().record(h._ref_t)
        self._vecs.append(t_vec)

        vm_vecs = {}  # comp_name -> h.Vector
        for comp_name, sec in comp_secs.items():
            vec = h.Vector()
            vec.record(sec(0.5)._ref_v)
            vm_vecs[comp_name] = vec
            self._vecs.append(vec)

        # Ionic current and conductance recordings.
        cur_vecs = {}   # (comp, mech) -> h.Vector
        cond_vecs = {}  # (comp, mech) -> h.Vector

        for comp_name, sec in comp_secs.items():
            seg = sec(0.5)
            for mech in all_mechs:
                key = (comp_name, mech)
                if not has_mech.get(key) or not _gbar_nonzero(sec, mech):
                    continue
                for var, store in ((cur_var.get(key), cur_vecs),
                                   (cond_var.get(key), cond_vecs)):
                    if var is None:
                        continue
                    ref = _try_ref(seg, var, mech)
                    if ref is not None:
                        vec = h.Vector()
                        vec.record(ref)
                        store[key] = vec
                        self._vecs.append(vec)

        # Silence pre-existing IClamps, inject our own, run.
        iclamp_snap = h.List('IClamp')
        n_prev = int(iclamp_snap.count())
        saved_amps = [float(iclamp_snap.o(i).amp) for i in range(n_prev)]
        for i in range(n_prev):
            iclamp_snap.o(i).amp = 0.0

        our_clamp = h.IClamp(soma(0.5))
        our_clamp.delay = tdel
        our_clamp.dur   = tdur
        our_clamp.amp   = self.iamp_nA

        v0 = cell.vm0 if cell.vm0 is not None else -60.0
        h.finitialize(v0)
        h.continuerun(tstop)

        our_clamp.amp = 0.0
        for i in range(n_prev):
            iclamp_snap.o(i).amp = saved_amps[i]

        # Cache numpy arrays.
        self._t_arr = np.array(t_vec)

        for comp_name in comp_secs:
            self._vm_data[comp_name] = np.array(vm_vecs[comp_name])

        def _arr(store, key):
            v = store.get(key)
            return np.array(v) if (v is not None and len(v) > 0) else None

        for comp_name in comp_secs:
            for mech in all_mechs:
                key = (comp_name, mech)
                self._cur_data[key]  = _arr(cur_vecs, key)
                self._cond_data[key] = _arr(cond_vecs, key)

        self._build_layout(comp_secs, all_mechs, has_mech, tstop)

    # ── layout ───────────────────────────────────────────────────────────────

    def _build_layout(self, comp_secs, all_mechs, has_mech, tstop):
        comp_names = list(comp_secs.keys())
        n_comps    = len(comp_names)
        n_mechs    = len(all_mechs)

        # Toolbar with toggle + reset buttons.
        toolbar = self.addToolBar('View')
        toolbar.setMovable(False)
        self._toggle_btn = QtWidgets.QPushButton('Show Conductance')
        self._toggle_btn.clicked.connect(self._toggle_mode)
        toolbar.addWidget(self._toggle_btn)
        reset_btn = QtWidgets.QPushButton('Reset Scale')
        reset_btn.clicked.connect(self._reset_scale)
        toolbar.addWidget(reset_btn)

        # Graphics layout grid:
        # Col 0:       row-label column (mechanism name, rotated)
        # Cols 1..N:   one column per compartment
        # Row 0:       compartment name headers (coloured labels)
        # Row 1:       Vm plots (one per compartment column, coloured)
        # Rows 2..M+1: one mechanism row per mechanism
        win_h = max(800, 100 + 90 * (n_mechs + 1))
        win_w = max(600, 80  + 210 * n_comps)

        # Assign one colour per compartment column.
        colors = [_COLUMN_COLORS[ci % len(_COLUMN_COLORS)] for ci in range(n_comps)]

        gw = pg.GraphicsLayoutWidget()
        self.setCentralWidget(gw)
        self.resize(win_w, win_h)

        # Row 0: compartment name column headers, coloured to match traces.
        gw.addLabel('', row=0, col=0)
        for ci, name in enumerate(comp_names):
            gw.addLabel(name, row=0, col=ci + 1,
                        bold=True, size='11pt',
                        color=colors[ci])

        # Row 1: one Vm plot per compartment; x and y axes linked across columns.
        gw.addLabel('Vm', row=1, col=0, angle=-90, bold=True, size='11pt')
        ref_p = None      # x-axis root for the whole grid
        ref_vm_p = None   # y-axis root for the Vm row
        vmin = 200.
        vmax = -200.
        self.vm_plots = {}
        for ci, comp_name in enumerate(comp_names):
            col = ci + 1
            pen = pg.mkPen(color=colors[ci], width=1)
            vm_p = gw.addPlot(row=1, col=col)
            vm_p.showGrid(x=False, y=True, alpha=0.3)
            vm_data = self._vm_data.get(comp_name)
            print(f"compartment: {comp_name}, vm_data: {vm_data[0]} to {vm_data[-1]}")
            if vm_data is not None:
                if np.min(vm_data) < vmin: vmin = np.min(vm_data)
                if np.max(vm_data) > vmax: vmax = np.max(vm_data)
                vm_p.plot(self._t_arr, vm_data, pen=pen)
                self._vm_plots[comp_name] = vm_p
            if ci == 0:
                vm_p.setLabel('left', 'Vm (mV)')
                vm_p.setTitle(
                    f'I = {self.iamp_nA * 1000:.0f} pA  '
                    f'(delay={self.durations[0]} ms, dur={self.durations[1]} ms)')
                ref_p    = vm_p
                ref_vm_p = vm_p
            else:
                vm_p.hideAxis('left')
                if ref_p    is not None: vm_p.setXLink(ref_p)
                # if ref_vm_p is not None: vm_p.setYLink(ref_vm_p)

        # ref_vm_p.setYRange(vmin, vmax, padding=0.05)
        print("vmin, vmax: ", vmin, vmax)
        for ci, comp_name in enumerate(comp_names):
            vm_p = self._vm_plots.get(comp_name)
            if vm_p is not None:
                vm_p.setYRange(vmin, vmax, padding=0.05)

        self._ref_p    = ref_p
        self._vm_y_root = ref_vm_p

        # Rows 2+: one per mechanism.
        # Within each row, plots with live curves are y-linked; the shared range
        # is computed explicitly in _reset_scale() after layout is complete.
        self._all_mechs = list(all_mechs)
        for ri, mech in enumerate(all_mechs):
            row = ri + 2
            gw.addLabel(mech, row=row, col=0, angle=-90,
                        size='10pt', color=(180, 180, 180))

            row_ref_y = None   # y-axis root for this mechanism row

            for ci, comp_name in enumerate(comp_names):
                key = (comp_name, mech)
                col = ci + 1
                pen = pg.mkPen(color=colors[ci], width=1)
                p = gw.addPlot(row=row, col=col)
                p.showGrid(x=False, y=True, alpha=0.3)
                if ref_p is not None:
                    p.setXLink(ref_p)
                self._plots[key] = p

                if ci > 0:
                    p.hideAxis('left')

                if not has_mech.get(key, False):
                    # Mechanism absent from this compartment — grey dash.
                    txt = pg.TextItem('—', anchor=(0.5, 0.5),
                                      color=(90, 90, 90))
                    p.addItem(txt)
                    txt.setPos(tstop / 2.0, 0.0)
                    p.setYRange(-0.01, 0.01)
                    self._curves[key] = None
                else:
                    arr = self._cur_data.get(key)
                    if arr is not None:
                        curve = p.plot(self._t_arr, arr, pen=pen)
                        self._curves[key] = curve
                        if ci == 0:
                            p.setLabel('left', 'mA/cm²')
                        # Link y-axis to first live plot in this row.
                        if row_ref_y is None:
                            row_ref_y = p
                        else:
                            p.setYLink(row_ref_y)
                        self._row_plots.setdefault(mech, []).append(p)
                    else:
                        txt = pg.TextItem('g = 0', anchor=(0.5, 0.5),
                                          color=(220, 200, 0))
                        p.addItem(txt)
                        txt.setPos(tstop / 2.0, 0.0)
                        p.setYRange(-0.01, 0.01)
                        self._curves[key] = None

            if row_ref_y is not None:
                self._row_y_roots[mech] = row_ref_y

        self._tstop = tstop
        # Apply correct x range and global per-row y ranges.
        if ref_p is not None:
            ref_p.setXRange(0, tstop, padding=0)
        self._reset_scale()
        self.show()

    # ── toggle ───────────────────────────────────────────────────────────────

    def _toggle_mode(self):
        self._mode = 'conductance' if self._mode == 'current' else 'current'
        use_cond  = (self._mode == 'conductance')
        data_dict = self._cond_data if use_cond else self._cur_data
        y_label   = 'mho/cm²'     if use_cond else 'mA/cm²'

        for key, curve in self._curves.items():
            if curve is None:
                continue
            arr = data_dict.get(key)
            p   = self._plots.get(key)
            if arr is not None:
                curve.setData(self._t_arr, arr)
                if p is not None:
                    p.setLabel('left', y_label)

        self._toggle_btn.setText(
            'Show Current' if use_cond else 'Show Conductance')
        # Re-apply per-row y ranges (conductance mode clamps min to 0).
        self._reset_scale()

    # ── y-range helpers ──────────────────────────────────────────────────────

    def _compute_row_y_range(self, mech: str) -> tuple[float, float]:
        """
        Return (ylo, yhi) for one mechanism row, spanning all compartment columns.

        In conductance mode the lower bound is always 0 (conductance ≥ 0).
        A 5 % padding is added above (and below for current mode).
        Falls back to (-0.01, 0.01) when no data exist for this row.
        """
        use_cond  = (self._mode == 'conductance')
        data_dict = self._cond_data if use_cond else self._cur_data

        ymin =  float('inf')
        ymax = -float('inf')
        for (_, m), arr in data_dict.items():
            if m != mech or arr is None:
                continue
            ymin = min(ymin, float(np.nanmin(arr)))
            ymax = max(ymax, float(np.nanmax(arr)))

        if not np.isfinite(ymin) or not np.isfinite(ymax) or ymin >= ymax:
            return (-0.01, 0.01)

        if use_cond:
            ymin = 0.0
            return (0.0, ymax * 1.05)

        span    = ymax - ymin
        padding = 0.05 * span
        return (ymin - padding, ymax + padding)

    def _reset_scale(self):
        """
        Apply the canonical y range for every mechanism row and the Vm row,
        then restore the x range.  Called on initial build, mode toggle, and
        by the Reset Scale toolbar button.
        """
        # Mechanism rows: apply range to every live plot in the row directly.
        for mech, plots in self._row_plots.items():
            lo, hi = self._compute_row_y_range(mech)
            for p in plots:
                p.setYRange(lo, hi, padding=0)

        # Vm row: compute union range and apply to every Vm plot directly.
        if self._vm_plots:
            vmin =  float('inf')
            vmax = -float('inf')
            for arr in self._vm_data.values():
                if arr is not None:
                    vmin = min(vmin, float(np.nanmin(arr)))
                    vmax = max(vmax, float(np.nanmax(arr)))
            if np.isfinite(vmin) and np.isfinite(vmax) and vmin < vmax:
                pad = 0.05 * (vmax - vmin)
                lo, hi = vmin - pad, vmax + pad
                # for p in self._vm_plots.values():
                #     p.setYRange(lo, hi, padding=0)

        # X range.
        if self._ref_p is not None:
            self._ref_p.setXRange(0, self._tstop, padding=0)
