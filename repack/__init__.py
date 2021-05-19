# Copyright (c) 2017-2021 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

__all__ = [
    'constants',
    'utils',
    ]

from . import constants
from . import utils
from .pack import *

from .VERSION import __version__

__all__ += pack.__all__


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
