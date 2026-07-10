"""
examples_gui.py  --  Launch pad for cnmodel examples.

Presents a dropdown of all example scripts.  For scripts that accept
command-line arguments the parameters are exposed in a pyqtgraph
ParameterTree; clicking Run starts the example in a subprocess so it
can open its own plot window without blocking or polluting this session.
"""

import argparse
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pyqtgraph as pg
import pyqtgraph.parametertree as ptree
from pyqtgraph.Qt import QtCore, QtWidgets

# ── paths ──────────────────────────────────────────────────────────────────────

EXAMPLES_DIR = Path(__file__).parent.resolve()
CNMODEL_DIR  = EXAMPLES_DIR.parent / "cnmodel"

def _mech_dir() -> Path:
    machine = platform.machine()
    return CNMODEL_DIR / "mechanisms" / ("arm64" if machine == "arm64" else "x86_64")


def get_mechanism_list() -> List[str]:
    """Return the same mechanism list that test_mechanisms.py would present."""
    # Keep in sync with tdurs dict in test_mechanisms.py
    known = [
        "CaPCalyx", "hcno", "hcno_bo", "ihpkj", "ihpyr", "ihpyr_adj",
        "ihpyrlc", "ihsgc_apical", "ihsgc_basalmiddle", "ihvcn", "inav11",
        "jsrnaf", "ka", "kcnq", "kdpyr", "kht", "kif", "kir", "kis",
        "klt", "nacn", "nacncoop", "nacsh", "nap", "nappyr", "napyr",
        "pkjlk", "rsg",
        "GRC_NA", "GRC_KV", "GRC_KA", "GRC_KM", "GRC_KIR", "GRC_KCA", "GRC_CA",
    ]
    remap = {
        "nav11":             "inav11",
        "jsrnaf":            "jsrna",
        "hcno_bo":           "hcnobo",
        "ihsgc_apical":      "ihsgcApical",
        "ihsgc_basalmiddle": "ihsgcBasalMiddle",
    }
    md = _mech_dir()
    if not md.exists():
        return sorted(known)
    obj_stems = {f.stem for f in md.glob("*.o")}
    avail: List[str] = []
    for k in known:
        if k in obj_stems:
            avail.append(k)
        elif k in remap and remap[k] in obj_stems:
            avail.append(remap[k])
    return sorted(avail) if avail else sorted(known)


# ── parameter definition helpers ───────────────────────────────────────────────

def _p(name: str, argname: Optional[str], ptype: str, default: Any,
       tip: str = "", values: Optional[List] = None,
       positional: bool = False) -> Dict:
    """Compact factory for parameter-definition dicts."""
    d: Dict[str, Any] = dict(name=name, argname=argname, ptype=ptype,
                              default=default, tip=tip, positional=positional)
    if values is not None:
        d["values"] = values
    return d


# ── example catalogue ──────────────────────────────────────────────────────────
#
# Each dict describes one example script:
#   label       display name in the combo box
#   script      filename relative to EXAMPLES_DIR
#   description one-line description shown below the combo
#   params      list of parameter dicts (see _p() above)
#
# ptype may be: 'list' | 'float' | 'int' | 'bool' | 'str'
# For 'list' supply values=[...].
# positional=True  → no flag on the command line, value placed first.
# argname=='--clampmode' is handled specially (maps to --cc / --vc / --rmp).

_CELL_TYPES = [
    "bushy", "bushycoop", "tstellate", "tstellatenav11", "dstellate",
    "dstellateeager", "sgc", "cartwheel", "pyramidal", "pyramidalceballos",
    "octopus", "tuberculoventral", "mso", "granule",
]
_SYN_CELLS   = ["sgc", "bushy", "tstellate", "dstellate",
                 "tuberculoventral", "pyramidal"]
_POST_CELLS  = ["bushy", "tstellate", "dstellate", "octopus",
                "tuberculoventral", "pyramidal"]

EXAMPLES: List[Dict] = [
    # ── toy_model (no args) ────────────────────────────────────────────────────
    {
        "label":       "toy_model",
        "script":      "toy_model.py",
        "description": "IV curves for all implemented cell types (no parameters).",
        "params":      [],
    },
    # ── test_adex (no args) ────────────────────────────────────────────────────
    {
        "label":       "test_adex",
        "script":      "test_adex.py",
        "description": "Adaptive-exponential (AdEx) integrate-and-fire model (no parameters).",
        "params":      [],
    },
    # ── test_an_model ─────────────────────────────────────────────────────────
    {
        "label":            "test_an_model",
        "script":           "test_an_model.py",
        "simulator_aware":  True,
        "description": "Auditory-nerve model: fiber response to acoustic stimuli.",
        "params": [
            _p("Species",       "--species",   "list",  "cat",
               "AN model species", ["cat", "human", "human_glasberg"]),
            _p("Stimulus",      "--stimulus",  "list",  "tone",
               "Acoustic stimulus type", ["tone", "noise", "SAM", "clicks"]),
            _p("dB SPL",        "--dB",        "float", 30.,  "Sound level (dB SPL)"),
            _p("Fiber type",    "--fibertype", "list",  "hsr",
               "Spontaneous-rate group", ["hsr", "msr", "lsr"]),
            _p("CF (Hz)",       "--CF",        "float", 16000., "Characteristic frequency (Hz)"),
            _p("Mod depth (%)", "--dmod",      "float", 100.,   "SAM modulation depth (%)"),
            _p("Mod freq (Hz)", "--fmod",      "float", 200.,   "SAM modulation frequency (Hz)"),
            _p("Rate-intensity","--RI",        "bool",  False,  "Run rate-intensity series"),
        ],
    },
    # ── test_bushy_variation ──────────────────────────────────────────────────
    {
        "label":           "test_bushy_variation",
        "script":          "test_bushy_variation.py",
        "simulator_aware": True,
        "description": "Bushy cell IV (panel a) or sound-driven PSTH (panel d).",
        "params": [
            _p("Panel", None, "list", "a",
               "Figure panel: a = IV curves, d = sound-driven PSTHs",
               ["a", "d"], positional=True),
        ],
    },
    # ── test_ccstim (no args) ─────────────────────────────────────────────────
    {
        "label":       "test_ccstim",
        "script":      "test_ccstim.py",
        "description": "Current-clamp stimulus waveform shapes (no parameters).",
        "params":      [],
    },
    # ── test_cells ────────────────────────────────────────────────────────────
    {
        "label":       "test_cells",
        "script":      "test_cells.py",
        "description": "IV / VC curves for an individual cell type.",
        "params": [
            _p("Cell type",      "celltype", "list", "granule",
               "Cell type to test", _CELL_TYPES, positional=True),
            _p("Species",        "species",  "list", "mouse",
               "Animal species", ["guineapig", "cat", "rat", "mouse"],
               positional=True),
            _p("Model sub-type", "--type",   "str",  "",
               "Model sub-type variant (leave blank for default)"),
            _p("Temperature (°C)", "--temp", "float", 30.0,
               "Bath temperature (°C)"),
            _p("Morphology",  "-m",    "list", "stick",
               "Morphological configuration", ["point", "waxon", "stick"]),
            _p("Na channel",  "--nav", "list", "std",
               "Sodium channel variant", ["std", "jsrna", "nav11", "nacncoop"]),
            _p("TTX",         "--ttx", "bool", False, "Block Na channels (TTX)"),
            _p("Pulse type",  "-p",    "list", "step",
               "CC pulse shape", ["step", "pulse"]),
            _p("Durations",   "-d",    "str",  "[10, 100, 50]",
               'Pulse timings as "[delay_ms, dur_ms, post_ms]"'),
            _p("Clamp mode",  "--clampmode", "list", "cc",
               "Recording mode", ["cc", "vc", "rmp"]),
        ],
    },
    # ── test_circuit (no args) ────────────────────────────────────────────────
    {
        "label":       "test_circuit",
        "script":      "test_circuit.py",
        "description": "Simple circuit connectivity test (no parameters).",
        "params":      [],
    },
    # ── test_decorator (no args) ──────────────────────────────────────────────
    {
        "label":       "test_decorator",
        "script":      "test_decorator.py",
        "description": "Channel decorator test (no parameters).",
        "params":      [],
    },
    # ── test_mechanisms ───────────────────────────────────────────────────────
    {
        "label":       "test_mechanisms",
        "script":      "test_mechanisms.py",
        "description": "Plot kinetics of one compiled ion-channel mechanism.",
        "params": [
            _p("Mechanism", "mechanism", "list", "klt",
               "Mechanism to display", get_mechanism_list(), positional=True),
            _p("Export plot", "--export", "bool", False, "Save plot to file"),
        ],
    },
    # ── test_mso_inputs (no args) ─────────────────────────────────────────────
    {
        "label":       "test_mso_inputs",
        "script":      "test_mso_inputs.py",
        "description": "MSO cell synaptic input test (no parameters).",
        "params":      [],
    },
    # ── test_physiology ───────────────────────────────────────────────────────
    {
        "label":           "test_physiology",
        "script":          "test_physiology.py",
        "simulator_aware": True,
        "description": "Network physiology simulation with configurable target cell type.",
        "params": [
            _p("Cell type",     "--celltype",       "list",  "bushy",
               "Target cell type in the CN network",
               ["bushy", "tstellate", "dstellate", "pyramidal", "tuberculoventral"]),
            _p("Species",       "--species",        "list",  "mouse",
               "Species for cell parameters and connectivity",
               ["mouse", "guineapig", "rat"]),
            _p("CF (Hz)",       "--cf",             "float", 16000.,
               "Characteristic frequency of the target cell (Hz)"),
            _p("F min (Hz)",    "--fmin",            "float", 4000.,
               "Minimum test frequency (Hz)"),
            _p("F max (Hz)",    "--fmax",            "float", 32000.,
               "Maximum test frequency (Hz)"),
            _p("Oct. spacing",  "--octave-spacing",  "float", 0.25,
               "Octave spacing between test frequencies"),
            _p("dB min",        "--db-min",          "float", 20.,
               "Minimum sound level (dB SPL)"),
            _p("dB max",        "--db-max",          "float", 100.,
               "Maximum sound level (dB SPL)"),
            _p("dB step",       "--db-step",         "float", 10.,
               "Sound level step size (dB)"),
            _p("Execution",     "--execution",       "list",  "parallel",
               "Parallel or serial execution", ["parallel", "serial"]),
            _p("Stim type",     "--stim-types",      "list",  "tone_pip",
               "Stimulus type to simulate",
               ["tone_pip", "bandlimited_noise", "wideband_noise"]),
            _p("Noise oct. BW", "--noise-octave-width", "float", 1.0,
               "Octave bandwidth for band-limited noise (default 1.0)"),
            _p("Ignore cache",  "--ignore-cache",    "bool",  False,
               "Bypass cache and rerun all simulations"),
        ],
    },
    # ── test_populations ──────────────────────────────────────────────────────
    {
        "label":           "test_populations",
        "script":          "test_populations.py",
        "simulator_aware": True,
        "description": "Population-level synapse test between two cell types.",
        "params": [
            _p("Pre cell",  "pre_cell",  "list", "sgc",
               "Presynaptic population", _SYN_CELLS, positional=True),
            _p("Post cell", "post_cell", "list", "bushy",
               "Postsynaptic population", _SYN_CELLS, positional=True),
        ],
    },
    # ── test_sgc_input (no args) ──────────────────────────────────────────────
    {
        "label":       "test_sgc_input",
        "script":      "test_sgc_input.py",
        "description": "SGC → cochlear-nucleus synapse test (no parameters).",
        "params":      [],
    },
    # ── test_sgc_input_PSTH ───────────────────────────────────────────────────
    {
        "label":           "test_sgc_input_PSTH",
        "script":          "test_sgc_input_PSTH.py",
        "simulator_aware": True,
        "description": "AN-driven PSTH in a postsynaptic cochlear-nucleus cell.",
        "params": [
            _p("Cell type",    "cell",       "list",  "bushy",
               "Target cell type", _POST_CELLS, positional=True),
            _p("Stimulus",     "--stimulus", "list",  "tone",
               "Acoustic stimulus", ["tone", "SAM", "noise", "clicks"]),
            _p("Synapse type", "--type",     "list",  "simple",
               "Synapse model", ["simple", "multisite"]),
            _p("Repetitions",  "--nrep",     "int",   10,
               "Number of stimulus repetitions"),
            _p("dB SPL",       "--dB",       "float", 50.,
               "Sound level (dB SPL)"),
            _p("SR group",     "--sr",       "int",   2,
               "AN SR group (1 = high, 2 = medium, 3 = low)"),
            _p("Parallel",     "--parallel", "list",  "mp",
               "Execution mode", ["serial", "mp", "multiprocessing"]),
            _p("Temperature",  "--temp",     "float", 34.0,
               "Simulation temperature (°C)"),
        ],
    },
    # ── test_sgc_input_phaselocking ───────────────────────────────────────────
    {
        "label":           "test_sgc_input_phaselocking",
        "script":          "test_sgc_input_phaselocking.py",
        "simulator_aware": True,
        "description": "Phase-locking of a CN cell to auditory-nerve input.",
        "params": [
            _p("Cell type",    "--celltype",  "list",  "bushy",
               "Target cell type", ["bushy", "tstellate", "octopus", "dstellate"]),
            _p("Species",      "--species",   "list",  "mouse",
               "Animal species", ["guineapig", "mouse", "rat"]),
            _p("Stimulus",     "--stimulus",  "list",  "tone",
               "Acoustic stimulus", ["tone", "SAM", "clicks"]),
            _p("CF (Hz)",      "--CF",        "float", 16000.,
               "Carrier frequency (Hz)"),
            _p("dB SPL",       "--dB",        "float", 60.,
               "Sound level (dB SPL)"),
            _p("Mod depth (%)", "--dmod",     "float", 100.,
               "SAM modulation depth (%)"),
            _p("Mod freq (Hz)", "--fmod",     "float", 200.,
               "SAM modulation frequency (Hz)"),
            _p("Fiber type",   "--fibertype", "list",  "hsr",
               "AN fiber SR group", ["hsr", "msr", "lsr"]),
            _p("Rate-intensity","--RI",       "bool",  False,
               "Run rate-intensity series"),
        ],
    },
    # ── test_simple_synapses ──────────────────────────────────────────────────
    {
        "label":       "test_simple_synapses",
        "script":      "test_simple_synapses.py",
        "description": "Simple point-synapse between two cell types.",
        "params": [
            _p("Pre cell",  "pre_cell",  "list", "sgc",
               "Presynaptic cell type", _SYN_CELLS, positional=True),
            _p("Post cell", "post_cell", "list", "bushy",
               "Postsynaptic cell type", _SYN_CELLS, positional=True),
        ],
    },
    # ── test_sound_stim ───────────────────────────────────────────────────────
    {
        "label":           "test_sound_stim",
        "script":          "test_sound_stim.py",
        "simulator_aware": True,
        "description": "Sound stimulus generation tests.",
        "params":      [],
    },
    # ── test_sounds (no args) ─────────────────────────────────────────────────
    {
        "label":       "test_sounds",
        "script":      "test_sounds.py",
        "description": "Sound synthesis tests (no parameters).",
        "params":      [],
    },
    # ── test_synapses ─────────────────────────────────────────────────────────
    {
        "label":       "test_synapses",
        "script":      "test_synapses.py",
        "description": "Full multisite synapse test between two cell types.",
        "params": [
            _p("Pre cell",  "precell",  "list", "sgc",
               "Presynaptic cell type",
               ["sgc", "tstellate", "dstellate", "tuberculoventral"],
               positional=True),
            _p("Post cell", "postcell", "list", "bushy",
               "Postsynaptic cell type", _POST_CELLS, positional=True),
            _p("Synapse type",       "--type",       "list", "multisite",
               "Synapse model", ["simple", "multisite"]),
            _p("Convergence = 1",    "--convergence","bool", False,
               "Force convergence = 1 for all pairs"),
        ],
    },
    # ── test_synaptic_mechanisms ──────────────────────────────────────────────
    {
        "label":       "test_synaptic_mechanisms",
        "script":      "test_synaptic_mechanisms.py",
        "description": "AMPA / NMDA receptor open-probability kinetics.",
        "params": [
            _p("Receptor",    "--ampa",     "list",  "nmda",
               "Receptor type to test", ["ampa", "nmda"]),
            _p("Show plot",   "--plot",     "bool",  False, "Display plot"),
            _p("Retrieve",    "--retrieve", "bool",  False,
               "Load previously saved result without re-running"),
            _p("Filename",    "--filename", "str",   "test_syn_mech.p",
               "Output data filename"),
        ],
    },
]


# ── command builder ────────────────────────────────────────────────────────────

def build_command(script: str, params: List[Dict],
                  values: Dict[str, Any]) -> List[str]:
    """Translate ParameterTree values into a subprocess argument list."""
    positional: List[str] = []
    optional:   List[str] = []

    for p in params:
        val     = values.get(p["name"])
        argname = p["argname"]
        ptype   = p["ptype"]

        if p.get("positional"):
            positional.append(str(val))

        elif ptype == "bool":
            if val:
                optional.append(argname)

        elif ptype == "multistr":
            # Each whitespace-separated token becomes its own argument after the flag.
            tokens = str(val).split()
            if tokens:
                optional.append(argname)
                optional.extend(tokens)

        elif ptype in ("float", "int", "str", "list"):
            sval = str(val).strip()
            if not sval:
                continue  # leave unset optional args at their argparse default

            # --clampmode maps to one of three flags (mutually exclusive group)
            if argname == "--clampmode":
                optional.append("--" + sval)
            else:
                optional.append(argname)
                optional.append(sval)

    return [sys.executable,
            str(EXAMPLES_DIR / script)] + positional + optional


# ── GUI ────────────────────────────────────────────────────────────────────────

class ExamplesGUI(QtWidgets.QMainWindow):
    def __init__(self, initial_example: Optional[str] = None) -> None:
        super().__init__()
        self.setWindowTitle("cnmodel Examples")
        self.resize(500, 640)
        self._procs: List[subprocess.Popen] = []
        self._param_root: Optional[ptree.Parameter] = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        # ── global simulator selector ─────────────────────────────────────────
        sim_row = QtWidgets.QHBoxLayout()
        sim_row.addWidget(QtWidgets.QLabel("Simulator:"))
        self.simulator_combo = QtWidgets.QComboBox()
        self.simulator_combo.addItems(["py3", "matlab"])
        self.simulator_combo.setCurrentIndex(0)
        sim_row.addWidget(self.simulator_combo)
        sim_row.addStretch(1)
        layout.addLayout(sim_row)

        # ── example selector ──────────────────────────────────────────────────
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Example:"))
        self.combo = QtWidgets.QComboBox()
        for ex in EXAMPLES:
            self.combo.addItem(ex["label"])
        row.addWidget(self.combo, stretch=1)
        layout.addLayout(row)

        # ── description ───────────────────────────────────────────────────────
        self.desc = QtWidgets.QLabel()
        self.desc.setWordWrap(True)
        self.desc.setStyleSheet("color: #777; font-style: italic;")
        layout.addWidget(self.desc)

        # ── parameter tree ────────────────────────────────────────────────────
        self.tree = ptree.ParameterTree(showHeader=False)
        layout.addWidget(self.tree, stretch=1)

        # ── run / status / quit ───────────────────────────────────────────────
        bot = QtWidgets.QHBoxLayout()
        self.run_btn = QtWidgets.QPushButton("Run")
        self.run_btn.setMinimumHeight(34)
        bot.addWidget(self.run_btn)
        self.status_lbl = QtWidgets.QLabel("")
        self.status_lbl.setStyleSheet("color: #080;")
        bot.addWidget(self.status_lbl, stretch=1)
        self.quit_btn = QtWidgets.QPushButton("Quit")
        self.quit_btn.setMinimumHeight(34)
        bot.addWidget(self.quit_btn)
        layout.addLayout(bot)

        self.combo.currentIndexChanged.connect(self._on_example_changed)
        self.run_btn.clicked.connect(self._on_run)
        self.quit_btn.clicked.connect(QtWidgets.QApplication.quit)

        # Pre-select example if one was given on the command line
        start_idx = 0
        if initial_example is not None:
            labels = [ex["label"] for ex in EXAMPLES]
            if initial_example in labels:
                start_idx = labels.index(initial_example)
            else:
                print(f"Warning: unknown example '{initial_example}'; "
                      f"choices are: {labels}", file=sys.stderr)
        self.combo.setCurrentIndex(start_idx)
        self._on_example_changed(start_idx)

    # ── slots ──────────────────────────────────────────────────────────────────

    def _on_example_changed(self, idx: int) -> None:
        ex = EXAMPLES[idx]
        self.desc.setText(ex["description"])
        self._rebuild_tree(ex["params"])

    def _rebuild_tree(self, params: List[Dict]) -> None:
        children = []
        for p in params:
            d: Dict[str, Any] = {"name": p["name"], "tip": p.get("tip", "")}
            if p["ptype"] == "list":
                d["type"]   = "list"
                d["limits"] = p["values"]
                d["value"]  = p["default"]
            elif p["ptype"] == "bool":
                d["type"]  = "bool"
                d["value"] = bool(p["default"])
            elif p["ptype"] == "float":
                d["type"]     = "float"
                d["value"]    = float(p["default"])
                d["decimals"] = 6
            elif p["ptype"] == "int":
                d["type"]  = "int"
                d["value"] = int(p["default"])
            elif p["ptype"] == "multistr":
                # Space-separated list of tokens — rendered as a plain text field.
                d["type"]  = "str"
                d["value"] = str(p["default"])
            else:  # str
                d["type"]  = "str"
                d["value"] = str(p["default"])
            children.append(d)

        self._param_root = ptree.Parameter.create(
            name="params", type="group", children=children
        )
        self.tree.setParameters(self._param_root, showTop=False)

    def _current_values(self, params: List[Dict]) -> Dict[str, Any]:
        values: Dict[str, Any] = {}
        if self._param_root is None:
            return values
        for p in params:
            try:
                values[p["name"]] = self._param_root.child(p["name"]).value()
            except KeyError:
                values[p["name"]] = p["default"]
        return values

    def _on_run(self) -> None:
        idx = self.combo.currentIndex()
        ex  = EXAMPLES[idx]
        cmd = build_command(
            ex["script"], ex["params"], self._current_values(ex["params"])
        )
        if ex.get("simulator_aware"):
            sim = self.simulator_combo.currentText()
            cmd += ["--simulator", sim]
        self.status_lbl.setText(f"Starting {ex['label']}…")
        QtWidgets.QApplication.processEvents()
        try:
            proc = subprocess.Popen(cmd)
            self._procs.append(proc)
            self.status_lbl.setText(f"Running (PID {proc.pid})")
        except Exception as exc:
            self.status_lbl.setStyleSheet("color: #c00;")
            self.status_lbl.setText(f"Error: {exc}")

    def closeEvent(self, event) -> None:
        # Leave child windows alive when the launcher closes.
        event.accept()


# ── entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="cnmodel examples launcher — select and run examples interactively."
    )
    parser.add_argument(
        "example", nargs="?", default=None,
        help=("Example to pre-select on launch "
              "(e.g. toy_model, test_cells, test_mechanisms). "
              "Omit to start at toy_model."),
    )
    args = parser.parse_args()

    app = pg.mkQApp("cnmodel Examples")
    win = ExamplesGUI(initial_example=args.example)
    win.show()
    pg.exec()


if __name__ == "__main__":
    main()
