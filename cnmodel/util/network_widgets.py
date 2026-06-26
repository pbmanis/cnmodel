"""Reusable Qt widgets for inspecting and controlling a cnmodel network.

Designed to be imported from test_physiology.py (or other scripts) as a
floating dialog that operates on a live CNSoundStim / NetworkSimProtocol
object without touching any on-disk data files.

Public API
----------
NetworkControlPanel  — QDialog with two tabs:
    "Convergence Table"   — read-only grid from cnmodel.data
    "Connection Isolation" — per-connection checkbox + scale spinbox
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore, QtGui

# Cell types treated as inhibitory pre-synaptic sources for colour coding.
INHIBITORY_PRE = {'dstellate', 'tuberculoventral', 'cartwheel'}


def _collect_synapse_handles(prot):
    """Enumerate all conductance handles in the live network.

    Walks every real cell's ``outputs`` list, inspects the PSD type, and
    records a ``(setter_fn, original_value)`` tuple for every individual
    release-site gmax (GluPSD / GlyPSD) or NetCon weight (Exp2PSD).

    Parameters
    ----------
    prot : CNSoundStim (or any object with a ``.populations`` OrderedDict)

    Returns
    -------
    dict
        ``{(pre_type, post_type): [(setter_fn, original_value), ...]}``
        Calling ``setter_fn(new_value)`` immediately changes the live NEURON
        mechanism; no data file is touched.
    """
    from cnmodel.synapses import GluPSD, GlyPSD, Exp2PSD

    handles: dict = {}

    for pre_pop in prot.populations.values():
        for i in pre_pop.real_cells():
            pre_cell = pre_pop._cells[i]['cell']
            if pre_cell == 0:
                continue
            for syn in getattr(pre_cell, 'outputs', []):
                psd = syn.psd
                term = syn.terminal

                try:
                    post_type = psd.cell.celltype
                except Exception:
                    continue

                key = (pre_pop.type, post_type)
                bucket = handles.setdefault(key, [])

                if isinstance(psd, GluPSD):
                    for mech in psd.ampa_psd:
                        g0 = float(mech.gmax)
                        bucket.append(
                            (lambda v, m=mech: setattr(m, 'gmax', v), g0)
                        )
                    for mech in psd.nmda_psd:
                        g0 = float(mech.gmax)
                        bucket.append(
                            (lambda v, m=mech: setattr(m, 'gmax', v), g0)
                        )
                elif isinstance(psd, GlyPSD):
                    for mech in psd.all_psd:
                        g0 = float(mech.gmax)
                        bucket.append(
                            (lambda v, m=mech: setattr(m, 'gmax', v), g0)
                        )
                elif isinstance(psd, Exp2PSD) and hasattr(term, 'netcon'):
                    nc = term.netcon
                    g0 = float(nc.weight[0])
                    bucket.append(
                        (lambda v, n=nc: n.weight.__setitem__(0, v), g0)
                    )

    return handles


class NetworkControlPanel(QtWidgets.QDialog):
    """Pre-simulation modal dialog for inspecting and configuring connections.

    Show with ``exec()`` before the simulation loop.  The dialog blocks until
    the user clicks "Run Simulation" (accept) or "Cancel" (reject).  Call
    ``get_scale_config()`` after ``exec()`` returns to find out whether any
    connections were modified (if so, skip the result cache).

    Tab 1 — Convergence Table
        Read-only grid from cnmodel.data for the current species.
        Green = excitatory, red = inhibitory pre-types.

    Tab 2 — Connection Isolation
        One row per active (pre_type → post_type) pair.  Checkbox unchecked
        → scale 0; spinbox sets scale 0.00–2.00.  Changes are applied
        immediately to the live NEURON synapse gmax values so they are in
        effect when prot.run() is called.

    Usage
    -----
    ::

        dlg = NetworkControlPanel()
        dlg.update_table(species)
        dlg.set_network(prot)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            any_modified = dlg.any_modified()
            # run simulations ...
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Configure Network — Build, then Run Simulation")
        self.resize(760, 560)

        self._prot = None                  # set by set_synapse_types()
        self._syn_handles: dict = {}       # {(pre, post): [(setter, orig), ...]}
        self._row_widgets: dict = {}       # {(pre, post): (QCheckBox, QDoubleSpinBox)}
        self._syn_type_widgets: dict = {}  # {(pre_type, post_type): QComboBox}

        root = QtWidgets.QVBoxLayout(self)
        self.tabs = QtWidgets.QTabWidget()
        root.addWidget(self.tabs)

        self._build_conv_tab()
        self._build_syntype_tab()   # tab 1 — must be configured before Build Network
        self._build_iso_tab()       # tab 2 — populated after Build Network fires

        # Button bar: Build Network → Run Simulation | Cancel
        btn_row = QtWidgets.QHBoxLayout()

        self._build_btn = QtWidgets.QPushButton("Build Network")
        self._build_btn.setToolTip(
            "Apply the synapse type choices above and create NEURON synapses.\n"
            "This must be clicked before running simulations."
        )
        self._build_btn.clicked.connect(self._on_build_network)
        btn_row.addWidget(self._build_btn)

        btn_row.addStretch()

        self._run_btn = QtWidgets.QPushButton("Run Simulation")
        self._run_btn.setEnabled(False)   # enabled after Build Network
        self._run_btn.clicked.connect(self.accept)
        btn_row.addWidget(self._run_btn)

        cancel_btn = QtWidgets.QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)

        root.addLayout(btn_row)

    # ── Tab builders ──────────────────────────────────────────────────────────

    def _build_conv_tab(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        self._species_label = QtWidgets.QLabel("Species: —")
        lay.addWidget(self._species_label)
        self.conv_table = QtWidgets.QTableWidget()
        self.conv_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        lay.addWidget(self.conv_table)
        self.tabs.addTab(w, "Convergence Table")

    def _build_syntype_tab(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)

        note = QtWidgets.QLabel(
            "Choose the synapse implementation for each pre→post connection.\n"
            "•  multisite — stochastic release zones with short-term plasticity (default).\n"
            "•  simple    — single Exp2Syn, no depression or facilitation.\n"
            "Click ‘Build Network’ to apply these choices and create NEURON synapses."
        )
        note.setWordWrap(True)
        lay.addWidget(note)

        self.syntype_table = QtWidgets.QTableWidget()
        self.syntype_table.setColumnCount(2)
        self.syntype_table.setHorizontalHeaderLabels(["Connection", "Synapse type"])
        hdr = self.syntype_table.horizontalHeader()
        hdr.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.syntype_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        lay.addWidget(self.syntype_table)

        self.tabs.addTab(w, "Synapse Types")

    def _build_iso_tab(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)

        note = QtWidgets.QLabel(
            "Changes take effect on the next simulation run.  "
            "Scale 1.0 = full strength · 0.0 = silenced · "
            "unchecking a row forces scale to 0."
        )
        note.setWordWrap(True)
        lay.addWidget(note)

        self.iso_table = QtWidgets.QTableWidget()
        self.iso_table.setColumnCount(3)
        self.iso_table.setHorizontalHeaderLabels(["Active", "Connection", "Scale"])
        hdr = self.iso_table.horizontalHeader()
        hdr.setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        hdr.setSectionResizeMode(
            1, QtWidgets.QHeaderView.ResizeMode.Stretch
        )
        hdr.setSectionResizeMode(
            2, QtWidgets.QHeaderView.ResizeMode.ResizeToContents
        )
        lay.addWidget(self.iso_table)

        btn_row = QtWidgets.QHBoxLayout()
        reset_btn = QtWidgets.QPushButton("Reset all to 1.0")
        reset_btn.clicked.connect(self._reset_all)
        btn_row.addWidget(reset_btn)
        btn_row.addStretch()
        lay.addLayout(btn_row)

        # Claude added 2026-06-26: Na-zero checkbox for conductance diagnostics
        na_row = QtWidgets.QHBoxLayout()
        self._na_zero_chk = QtWidgets.QCheckBox(
            "Zero Na channels on target cell (reveals subthreshold conductances)"
        )
        self._na_zero_chk.setToolTip(
            "Sets all Na-channel gbar to 0 on target population cells before running.\n"
            "Removes action potentials so pre-synaptic conductance waveforms are visible."
        )
        na_row.addWidget(self._na_zero_chk)
        na_row.addStretch()
        lay.addLayout(na_row)

        # Claude added 2026-06-26: opt-in conductance recording (bypasses cache)
        syn_g_row = QtWidgets.QHBoxLayout()
        self._rec_syn_g_chk = QtWidgets.QCheckBox(
            "Record synaptic conductances on target cell (debug; slower)"
        )
        self._rec_syn_g_chk.setToolTip(
            "Before each run, attaches h.Vector recorders to every non-SGC\n"
            "PSD conductance on target-population cells.  Results are shown\n"
            "in the conductance panel.  Forces a fresh simulation (no cache)."
        )
        syn_g_row.addWidget(self._rec_syn_g_chk)
        syn_g_row.addStretch()
        lay.addLayout(syn_g_row)

        self.tabs.addTab(w, "Connection Isolation")

    # ── Public API ────────────────────────────────────────────────────────────

    def set_synapse_types(self, prot):
        """Populate the Synapse Types tab from the live network topology in *prot*.

        Stores a reference to *prot* so that the 'Build Network' button can call
        prot._resolve_network() when clicked.  Call before exec().
        """
        self._prot = prot
        self._syn_type_widgets = {}
        self.syntype_table.setRowCount(0)

        for post_pop in prot.populations.values():
            for pre_pop in post_pop._pre_connections:
                key = (pre_pop.type, post_pop.type)
                if key in self._syn_type_widgets:
                    continue  # skip duplicates (same pair via multiple paths)

                row = self.syntype_table.rowCount()
                self.syntype_table.insertRow(row)

                label = QtWidgets.QTableWidgetItem(
                    f"{pre_pop.type.title()} → {post_pop.type.title()}"
                )
                label.setFlags(label.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.syntype_table.setItem(row, 0, label)

                combo = QtWidgets.QComboBox()
                combo.addItems(["multisite", "simple"])
                combo.setCurrentText("multisite")
                self.syntype_table.setCellWidget(row, 1, combo)
                self._syn_type_widgets[key] = combo

        self.syntype_table.resizeColumnsToContents()

    def _on_build_network(self):
        """Apply synapse type choices, build NEURON synapses, populate gmax tab."""
        if self._prot is None:
            return

        # Push per-pre synapse type overrides onto each post-population.
        for (pre_type, post_type), combo in self._syn_type_widgets.items():
            post_pop = self._prot.populations.get(post_type)
            if post_pop is not None:
                post_pop._synapsetype_per_pre[pre_type] = combo.currentText()

        # Create all NEURON synapses with those types.
        self._prot._resolve_network()

        # Now that synapses exist, populate the Connection Isolation tab.
        self.set_network(self._prot)

        # Lock the synapse type widgets — can't change type after build.
        for combo in self._syn_type_widgets.values():
            combo.setEnabled(False)
        self._build_btn.setEnabled(False)
        self._build_btn.setText("Network Built ✓")

        # Enable Run and switch view to Connection Isolation tab.
        self._run_btn.setEnabled(True)
        self._run_btn.setDefault(True)
        self.tabs.setCurrentIndex(2)

    def get_synapse_type_config(self) -> dict:
        """Return {(pre_type, post_type): synapse_type_string} for all connections."""
        return {key: combo.currentText() for key, combo in self._syn_type_widgets.items()}

    def update_table(self, species: str):
        """Populate the convergence table for *species*."""
        from cnmodel import data as cdata
        from cnmodel.data._db import get_table_info

        self._species_label.setText(f"Species: {species}")

        info = get_table_info('convergence')
        all_pre = info.get('pre_type', [])
        all_post = info.get('post_type', [])

        # Only show rows/cols that have at least one value for this species.
        valid_pre: list = []
        valid_post: list = []
        value_cache: dict = {}
        for pre in all_pre:
            for post in all_post:
                try:
                    v = cdata.get('convergence', species=species,
                                  pre_type=pre, post_type=post)
                    value_cache[(pre, post)] = v
                    if pre not in valid_pre:
                        valid_pre.append(pre)
                    if post not in valid_post:
                        valid_post.append(post)
                except KeyError:
                    pass

        self.conv_table.setRowCount(len(valid_pre))
        self.conv_table.setColumnCount(len(valid_post))
        self.conv_table.setHorizontalHeaderLabels(valid_post)
        self.conv_table.setVerticalHeaderLabels(valid_pre)

        for r, pre in enumerate(valid_pre):
            for c, post in enumerate(valid_post):
                v = value_cache.get((pre, post))
                if v is None:
                    self.conv_table.setItem(
                        r, c, QtWidgets.QTableWidgetItem("")
                    )
                    continue
                numeric = v[0] if isinstance(v, tuple) else v
                if isinstance(v, tuple):
                    text = f"{v[0]:.1f}±{v[1]:.1f}"
                else:
                    text = f"{v:.1f}" if v != 0 else "0"
                item = QtWidgets.QTableWidgetItem(text)
                item.setTextAlignment(
                    QtCore.Qt.AlignmentFlag.AlignRight
                    | QtCore.Qt.AlignmentFlag.AlignVCenter
                )
                if numeric != 0:
                    if pre in INHIBITORY_PRE:
                        item.setBackground(QtGui.QColor(90, 30, 30))
                    else:
                        item.setBackground(QtGui.QColor(30, 80, 30))
                self.conv_table.setItem(r, c, item)

        self.conv_table.resizeColumnsToContents()

    def set_network(self, prot):
        """Rebuild the isolation tab from the live network in *prot*.

        Safe to call multiple times (e.g., after a network rebuild).
        """
        self._syn_handles = _collect_synapse_handles(prot)
        self._row_widgets = {}
        self.iso_table.setRowCount(0)

        for key in sorted(self._syn_handles):
            pre_type, post_type = key
            handles = self._syn_handles[key]

            row = self.iso_table.rowCount()
            self.iso_table.insertRow(row)

            # col 0 — checkbox (centred in a wrapper widget)
            chk = QtWidgets.QCheckBox()
            chk.setChecked(True)
            chk_cell = QtWidgets.QWidget()
            chk_lay = QtWidgets.QHBoxLayout(chk_cell)
            chk_lay.addWidget(chk)
            chk_lay.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            chk_lay.setContentsMargins(0, 0, 0, 0)
            self.iso_table.setCellWidget(row, 0, chk_cell)

            # col 1 — label
            label_str = (
                f"{pre_type.title()} → {post_type.title()}"
                f"  ({len(handles)} handles)"
            )
            label_item = QtWidgets.QTableWidgetItem(label_str)
            label_item.setFlags(
                label_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable
            )
            self.iso_table.setItem(row, 1, label_item)

            # col 2 — spinbox
            spin = QtWidgets.QDoubleSpinBox()
            spin.setRange(0.0, 2.0)
            spin.setSingleStep(0.05)
            spin.setDecimals(2)
            spin.setValue(1.0)
            self.iso_table.setCellWidget(row, 2, spin)

            self._row_widgets[key] = (chk, spin)
            chk.toggled.connect(lambda checked, k=key: self._on_change(k))
            spin.valueChanged.connect(lambda val, k=key: self._on_change(k))

        self.iso_table.resizeColumnsToContents()

    # ── Private helpers ───────────────────────────────────────────────────────

    def _on_change(self, key):
        chk, spin = self._row_widgets[key]
        scale = spin.value() if chk.isChecked() else 0.0
        for setter, orig in self._syn_handles[key]:
            setter(orig * scale)

    def _reset_all(self):
        """Restore all connections to scale 1.0 and re-enable checkboxes."""
        for key, (chk, spin) in self._row_widgets.items():
            chk.blockSignals(True)
            spin.blockSignals(True)
            chk.setChecked(True)
            spin.setValue(1.0)
            chk.blockSignals(False)
            spin.blockSignals(False)
        for key, handles in self._syn_handles.items():
            for setter, orig in handles:
                setter(orig)

    def get_scale_config(self) -> dict:
        """Return the effective scale for every connection type.

        Returns
        -------
        dict  {(pre_type, post_type): float}
            1.0  = full strength (default)
            0.0  = silenced (checkbox unchecked or spinbox at 0)
            other = scaled up/down
        """
        return {
            key: (spin.value() if chk.isChecked() else 0.0)
            for key, (chk, spin) in self._row_widgets.items()
        }

    def get_na_zero(self) -> bool:
        """Return True if Na channels should be zeroed on the target cell."""
        return self._na_zero_chk.isChecked()

    def get_record_syn_conductance(self) -> bool:
        """Return True if synaptic conductances should be recorded during the run."""
        return self._rec_syn_g_chk.isChecked()

    def any_modified(self) -> bool:
        """Return True if anything differs from defaults (forces a fresh simulation)."""
        return (
            self._na_zero_chk.isChecked()
            or self._rec_syn_g_chk.isChecked()
            or any(c.currentText() != 'multisite' for c in self._syn_type_widgets.values())
            or any(abs(v - 1.0) > 1e-9 for v in self.get_scale_config().values())
        )
