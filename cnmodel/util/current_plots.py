"""
Per-mechanism current traces for cnmodel cells.

CurrentPlotsWindow(cell, iamp_nA, durations) opens a pyqtgraph window
with Vm (row 0) and per-mechanism ionic current traces (rows 1..N),
all recorded from the soma.  Mechanisms with zero conductance show a
'g = 0' label instead of a trace.  All x-axes are linked.
"""
from neuron import h
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets

# NEURON current-variable RANGE name for each mechanism
MECH_CURRENT_VAR = {
    'GRC_NA':    'ina',
    'GRC_KV':    'ik',
    'GRC_KA':    'ik',
    'GRC_KM':    'ik',
    'GRC_KIR':   'ik',
    'GRC_KCA':   'ik',
    'GRC_CA':    'ica',
    'GRC_LKG1':  'i',
    'GRC_LKG2':  'i',
    'nacn':      'ina',
    'na':        'ina',
    'nav11':     'ina',
    'nacncoop':  'ina',
    'kht':       'ik',
    'klt':       'ik',
    'ka':        'ik',
    'kd':        'ik',
    'ihvcn':     'ih',
    'hcno':      'ih',
    'leak':      'i',
    'pas':       'i',
}

# Conductance variable name per mechanism (default 'gbar')
MECH_COND_VAR = {
    'GRC_LKG1': 'gl',
    'GRC_LKG2': 'ggaba',
}

# Mechanisms that have no plottable ionic current
SKIP_MECHS = {'GRC_CALC', 'cadiff', 'cap', 'extracellular'}


def _cond_var(mech):
    return MECH_COND_VAR.get(mech, 'gbar')


def _gbar_nonzero(sec, mech):
    """True if mech is inserted in sec with conductance > 0."""
    try:
        g = getattr(getattr(sec(0.5), mech), _cond_var(mech))
        return float(g) > 0.0
    except (AttributeError, RuntimeError):
        return False


def _current_var(mech, seg):
    """Return NEURON RANGE variable for mech's current at seg, or None."""
    known = MECH_CURRENT_VAR.get(mech, '__probe__')
    if known is None:
        return None
    candidates = [] if known == '__probe__' else [known]
    for c in ('ina', 'ik', 'ica', 'ih', 'i'):
        if c not in candidates:
            candidates.append(c)
    for cvar in candidates:
        try:
            getattr(seg, f'_ref_{cvar}_{mech}')
            return cvar
        except (AttributeError, RuntimeError):
            pass
    return None


def _mechs_in_sec(sec):
    """Return list of density mechanism names inserted in sec."""
    return list(sec.psection().get('density_mechs', {}).keys())


class CurrentPlotsWindow(QtWidgets.QMainWindow):
    """
    Pyqtgraph window: Vm (row 0) + per-mechanism currents (rows 1..N),
    all recorded from cell.soma during a single current-clamp run.
    """

    def __init__(self, cell, iamp_nA=0.020, durations=(10, 100, 50)):
        super().__init__()
        self.setWindowTitle('Current Plots')
        self.cell = cell
        self.iamp_nA = iamp_nA
        self.durations = list(durations)
        self._vecs = []  # keep h.Vector refs alive for the lifetime of the window
        self._run_and_plot()

    def _run_and_plot(self):
        cell = self.cell
        tdel, tdur, tpost = self.durations
        tstop = tdel + tdur + tpost
        soma = cell.soma

        # Discover mechanisms inserted in soma
        all_mechs = _mechs_in_sec(soma)
        seg = soma(0.5)

        plot_mechs = []
        cvar_map = {}
        for mech in all_mechs:
            if mech in SKIP_MECHS:
                continue
            cvar = _current_var(mech, seg)
            if cvar is not None:
                plot_mechs.append(mech)
                cvar_map[mech] = cvar

        # Set up time and Vm recordings
        t_vec = h.Vector().record(h._ref_t)
        vm_vec = h.Vector().record(soma(0.5)._ref_v)
        self._vecs += [t_vec, vm_vec]

        zero_g = {}
        cur_vecs = {}
        for mech in plot_mechs:
            if not _gbar_nonzero(soma, mech):
                zero_g[mech] = True
                continue
            zero_g[mech] = False
            cvar = cvar_map[mech]
            vec = h.Vector()
            try:
                vec.record(getattr(seg, f'_ref_{cvar}_{mech}'))
                cur_vecs[mech] = vec
                self._vecs.append(vec)
            except (AttributeError, RuntimeError):
                zero_g[mech] = True

        # Null pre-existing IClamps so they don't interfere with our run
        iclamp_snap = h.List('IClamp')
        n_prev = int(iclamp_snap.count())
        saved_amps = []
        for i in range(n_prev):
            ic = iclamp_snap.o(i)
            saved_amps.append(float(ic.amp))
            ic.amp = 0.0

        our_clamp = h.IClamp(soma(0.5))
        our_clamp.delay = tdel
        our_clamp.dur = tdur
        our_clamp.amp = self.iamp_nA

        v0 = cell.vm0 if cell.vm0 is not None else -60.0
        h.finitialize(v0)
        h.continuerun(tstop)

        # Disable our clamp and restore the originals
        our_clamp.amp = 0.0
        for i in range(n_prev):
            iclamp_snap.o(i).amp = saved_amps[i]

        # Build the plot layout
        t_arr = np.array(t_vec)
        n_rows = 1 + len(plot_mechs)
        win_h = max(1440, 140 * n_rows)

        gw = pg.GraphicsLayoutWidget()
        self.setCentralWidget(gw)
        self.resize(600, win_h)

        ref_p = None

        # Row 0: Vm
        p = gw.addPlot(row=0, col=0)
        p.plot(t_arr, np.array(vm_vec), pen='w')
        p.setLabel('left', 'Vm (mV)')
        p.setTitle(f'Soma  —  I = {self.iamp_nA * 1000:.0f} pA  '
                   f'(delay={tdel} ms, dur={tdur} ms)')
        p.showGrid(x=False, y=True, alpha=0.3)
        ref_p = p

        # Rows 1+: one per mechanism
        for ri, mech in enumerate(plot_mechs):
            p = gw.addPlot(row=1 + ri, col=0)
            if zero_g.get(mech, True):
                txt = pg.TextItem('g = 0', anchor=(0.5, 0.5), color=(220, 200, 0))
                p.addItem(txt)
                txt.setPos(tstop / 2.0, 0.0)
                p.setXRange(0, tstop)
                p.setYRange(-0.01, 0.01)
                p.hideAxis('left')
            else:
                vec = cur_vecs.get(mech)
                if vec is not None and len(vec) > 0:
                    p.plot(t_arr, np.array(vec), pen='w')
                p.setLabel('left', 'mA/cm²')
            cvar_label = cvar_map.get(mech, '?')
            p.setTitle(f'{mech}  ({cvar_label})')
            p.showGrid(x=False, y=True, alpha=0.3)
            if ref_p is not None:
                p.setXLink(ref_p)

        self.show()
