#! /usr/bin/env python

# Copyright (c) 2017-2018 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys
import warnings
import time
from datetime import date

import repack as rep

if __name__ == "__main__":
  warnings.simplefilter("ignore", RuntimeWarning)
  if len(sys.argv) != 2:
    print("\n{:s}\n  Error: Wrong usage.\n{:s}\n".format(70*":", 70*":"))
    sys.exit(0)

  print("\n{:s}\n"
        "  repack: line-transition data compression.\n"
        "  Version {:s}.\n"
        "  Copyright (c) 2017-{:d} Patricio Cubillos.\n"
        "  repack is open-source software under the MIT license.\n"
        "{:s}\n\n".format(":"*70, rep.__version__, date.today().year, ":"*70))

  print("Starting: {:s}".format(time.ctime()))
  cfile =  sys.argv[1]
  args = rep.parser(cfile)
  rep.repack(*args)
  print("End: {:s}".format(time.ctime()))
