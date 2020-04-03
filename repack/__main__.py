#! /usr/bin/env python

# Copyright (c) 2017-2020 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys
import warnings
import time
from datetime import date

import repack as rep


banner = ':'*70

def main():
    warnings.simplefilter("ignore", RuntimeWarning)
    if len(sys.argv) != 2:
        print(f"\n{banner}\n  Error: Wrong usage.\n{banner}\n")
        sys.exit(0)

    print(f"\n{banner}\n"
           "  repack: line-transition data compression.\n"
          f"  Version {rep.__version__}.\n"
          f"  Copyright (c) 2017-{date.today().year} Patricio Cubillos.\n"
           "  repack is open-source software under the MIT license.\n"
          f"{banner}\n\n")

    print(f"Start: {time.ctime()}")
    cfile = sys.argv[1]
    rep.repack(cfile)
    print(f"End: {time.ctime()}")

if __name__ == "__main__":
    main()
