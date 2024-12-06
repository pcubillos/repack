import os
import bz2
import shutil
import subprocess

import numpy as np
import repack

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
os.chdir(ROOT+'tests')


def test_exomol_sort():
    states_file = 'data/1H2-16O__POKAZATEL.states.bz2'
    lbl_file = 'data/1H2-16O__POKAZATEL__34800-34900.trans.bz2'
    shutil.copy(lbl_file+'.original', lbl_file)
    os.remove('data/README_REPACK')

    subprocess.call('repack -sort exomol_sort.cfg'.split())

    e, g = repack.utils.read_states(states_file)
    with bz2.open(lbl_file, "rt") as file:
        lines = file.readlines()

    nlines = len(lines)
    up = np.zeros(nlines, int)
    lo = np.zeros(nlines, int)
    for i,line in enumerate(lines):
        up[i] = line[ 0:12]
        lo[i] = line[13:25]
    wn = e[up-1] - e[lo-1]

    assert np.all(np.ediff1d(wn)>=0)


def test_exomol_repack_sorted(capfd):
    subprocess.call('repack exomol_sort.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/1H2-16O__POKAZATEL__34800-34900.trans.bz2'.
Reading: 'data/1H2-16O__POKAZATEL__34800-34900.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       98.46%,      2,375/   154,580 lines.
  Flagging lines at 1000 K:
  Compression rate:       97.35%,      4,095/   154,580 lines.
  Total compression rate: 97.11%,      4,465/   154,580 lines.

With a threshold strength factor of 0.1,
kept a total of 4,465 line transitions out of 154,580 lines.

Successfully rewriten exomol line-transition info into:
  'H2O_exomol_0.26um_500-1000K_lbl.dat'.""" in capfd.out
