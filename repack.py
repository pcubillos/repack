#! /usr/bin/env python
import sys
import warnings

from repack import *

if __name__ == "__main__":
  warnings.simplefilter("ignore", RuntimeWarning)
  if len(sys.argv) != 2:
    print("\n{:s}\n  Error: Wrong usage.\n{:s}\n".format(70*":", 70*":"))
    sys.exit(0)

  cfile =  sys.argv[1]
  args = parser(cfile)
  repack(*args)
