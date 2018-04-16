# Copyright (c) 2017-2018 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys, os

__all__ = ["parser", "repack", "utils", "constants"]

from . import utils
from . import constants
from . import VERSION as ver
from .repack import *

__version__ = "{:d}.{:d}.{:d}".format(ver.repack_VER,
                                      ver.repack_MIN, ver.repack_REV)

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
