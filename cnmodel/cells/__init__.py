"""
Cell definitions for models.

This class includes a number of different cell definitions and default
conductances for point models. 
"""

from .granule import *
from .bushy import *
from .tstellate import *
from .dstellate import *
from .cartwheel import *
from .pyramidal import *
from .sgc import *
from .octopus import *
from .tuberculoventral import *
from .msoprincipal import *
from .hh import *

__all__ = ["bushy", "tstellate", "dstellate", "cartwheel", "pyramidal", "sgc", "octopus",
            "tuberculoventral", "granule", "msoprincipal", "hh"]

from .cell import Cell

def cell_from_section(sec):
    return Cell.from_section(sec)

