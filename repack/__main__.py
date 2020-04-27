#! /usr/bin/env python

# Copyright (c) 2017-2020 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys
import warnings
import time
from datetime import date

import repack


banner = ':'*70
error_message = """\
  Error: Wrong usage.  For repacking:
      repack config_file.cfg
  For sorting (of MARVELized data):
      repack -sort config_file.cfg"""

def main():
    warnings.simplefilter("ignore", RuntimeWarning)
    if len(sys.argv) == 3:
        if sys.argv[1] != '-sort':
            print(f"\n{banner}\n{error_message}.\n{banner}\n")
            sys.exit(0)

    elif len(sys.argv) != 2:
        print(f"\n{banner}\n  Error: Wrong usage.\n{banner}\n")
        sys.exit(0)

    print(f"\n{banner}\n"
           "  repack: line-transition data compression.\n"
          f"  Version {repack.__version__}.\n"
          f"  Copyright (c) 2017-{date.today().year} Patricio Cubillos.\n"
           "  repack is open-source software under the MIT license.\n"
          f"{banner}\n\n")

    if len(sys.argv) == 2 and sys.argv[1] == '-v':
        return

    print(f"Start: {time.ctime()}")
    if len(sys.argv) == 3:
        repack.sort(sys.argv[2])
    else:
        repack.repack(sys.argv[1])
    print(f"End: {time.ctime()}")

if __name__ == "__main__":
    main()
