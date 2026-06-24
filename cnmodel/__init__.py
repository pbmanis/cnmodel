__author__ = "Paul B. Manis and Luke Campagnola"
__version__ = "0.58"
from pathlib import Path
try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass
import platform
import logging
logging.basicConfig(level=logging.INFO, format="[%(process)s] %(message)s")
import os
dirname = os.path.abspath(os.path.dirname(__file__))
libpath = Path(dirname)
if platform.system() == 'Windows':
    libpath2 = libpath / 'mechanisms'/ 'x86_64' / 'nrnmech.dll'
elif platform.system() == 'Darwin':
    libpath2 = libpath / 'mechanisms' / 'arm64' / 'libnrnmech.dylib'
else:
    libpath2 = libpath / 'mechanisms' / 'x86_64' / 'libnrnmech.so'
print("LIBPATH:", libpath2, libpath2.exists())
logging.info(f"cnmodel: Loading NEURON mechanisms from {libpath2}")
import neuron
try:
    neuron.h.MultiSiteSynapse  # already loaded — skip
    # print("cnmodel: NEURON mechanisms already loaded.")
except AttributeError:
    try:
        neuron.h.nrn_load_dll(str(libpath2))
        neuron.h.MultiSiteSynapse
        print("cnmodel: Successfully imported NEURON mechanisms!")
    except AttributeError:
        print("cnmodel: Failed to load NEURON mechanisms from", libpath)
        print("cnmodel: Checked path exists:", libpath2.exists())
        raise FileNotFoundError(f"Could not find NEURON mechanisms at {libpath}")
MODFILE_PATH = str(libpath2.parent)

# flag to allow unit tests to store / overwrite test results
AUDIT_TESTS = False



