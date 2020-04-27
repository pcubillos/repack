import os
import subprocess
import pytest
import repack

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
os.chdir(ROOT+'tests')


def test_exomol_sort():
    if os.path.exists('data/1H2-16O__POKAZATEL__34800-34900.trans.sort.bz2'):
        os.remove('data/1H2-16O__POKAZATEL__34800-34900.trans.sort.bz2')
    subprocess.call('repack -sort exomol_sort.cfg'.split())
    assert os.path.exists('data/1H2-16O__POKAZATEL__34800-34900.trans.sort.bz2')


def test_exomol_repack_sorted(capfd):
    subprocess.call('repack exomol_sort.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/1H2-16O__POKAZATEL__34800-34900.trans.sort.bz2'.
Reading: 'data/1H2-16O__POKAZATEL__34800-34900.trans.sort.bz2'.
  Flagging lines at  500 K:
  Compression rate:       98.46%,      2,375/   154,580 lines.
  Flagging lines at 1000 K:
  Compression rate:       97.35%,      4,095/   154,580 lines.
  Total compression rate: 97.11%,      4,465/   154,580 lines.

With a threshold strength factor of 0.1,
kept a total of 4,465 line transitions out of 154,580 lines.

Successfully rewriten exomol line-transition info into:
  'H2O_exomol_0.26um_500-1000K_lbl.dat'.""" in capfd.out
