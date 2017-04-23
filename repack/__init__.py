import sys, os

__all__ = ["parser", "repack", "utils", "constants", "cutils"]

from .repack import *
from . import utils
from . import constants
from . import cutils

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
