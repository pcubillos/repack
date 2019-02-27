import os
import sys
import pytest
import subprocess

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'

os.chdir(ROOT+'tests')


def test_exomol_single(capfd):
    subprocess.call(['python', '../repack.py', 'exomol_repack_single.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Reading: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       57.00%,         43/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       54.00%,         46/       100 lines.
  Total compression rate: 53.00%,         47/       100 lines.

Kept a total of 47 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files(capfd):
    subprocess.call(['python','../repack.py','exomol_repack_two_files.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Unzipping: '14N-1H3__MockBYTe__00200-00300.trans.bz2'.
Reading: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       57.00%,         43/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       54.00%,         46/       100 lines.
  Total compression rate: 53.00%,         47/       100 lines.

Reading: '14N-1H3__MockBYTe__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       80.00%,         20/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       77.00%,         23/       100 lines.
  Total compression rate: 75.00%,         25/       100 lines.

Kept a total of 72 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_isotopes(capfd):
    subprocess.call(['python','../repack.py','exomol_repack_two_isotopes.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Unzipping: '15N-1H3__MockBYTe-15__00100-00200.trans.bz2'.
Reading: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Reading: '15N-1H3__MockBYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       50.00%,        100/       200 lines.
  Flagging lines at  700 K:
  Compression rate:       50.00%,        100/       200 lines.
  Total compression rate: 48.00%,        104/       200 lines.

Kept a total of 104 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files_two_iso(capfd):
    subprocess.call(['python','../repack.py','exomol_repack_two_two.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Unzipping: '15N-1H3__MockBYTe-15__00100-00200.trans.bz2'.
Unzipping: '14N-1H3__MockBYTe__00200-00300.trans.bz2'.
Unzipping: '15N-1H3__MockBYTe-15__00200-00300.trans.bz2'.
Reading: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Reading: '15N-1H3__MockBYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       50.00%,        100/       200 lines.
  Flagging lines at  700 K:
  Compression rate:       50.00%,        100/       200 lines.
  Total compression rate: 48.00%,        104/       200 lines.

Reading: '14N-1H3__MockBYTe__00200-00300.trans.bz2'.
Reading: '15N-1H3__MockBYTe-15__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       72.00%,         56/       200 lines.
  Flagging lines at  700 K:
  Compression rate:       72.50%,         55/       200 lines.
  Total compression rate: 69.50%,         61/       200 lines.

Kept a total of 165 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call(['bzip2','-dk','14N-1H3__MockBYTe__00100-00200.trans.bz2'])
    subprocess.call(['bzip2','-dk','14N-1H3__MockBYTe.states.bz2'])
    subprocess.call(['python','../repack.py','exomol_repack_single_unzip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '14N-1H3__MockBYTe__00100-00200.trans'.
  Flagging lines at  500 K:
  Compression rate:       57.00%,         43/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       54.00%,         46/       100 lines.
  Total compression rate: 53.00%,         47/       100 lines.

Kept a total of 47 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out
    os.remove('14N-1H3__MockBYTe__00100-00200.trans')
    os.remove('14N-1H3__MockBYTe.states')


def test_exomol_single_chunks(capfd):
    subprocess.call(['python','../repack.py','exomol_repack_single_chunks.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
Reading: '14N-1H3__MockBYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K (chunk 1/2):
  Compression rate:       50.00%,         25/        50 lines.
  Flagging lines at  700 K:
  Compression rate:       50.00%,         25/        50 lines.
  Total compression rate: 48.00%,         26/        50 lines.

  Flagging lines at  500 K (chunk 2/2):
  Compression rate:       64.00%,         18/        50 lines.
  Flagging lines at  700 K:
  Compression rate:       58.00%,         21/        50 lines.
  Total compression rate: 58.00%,         21/        50 lines.

Kept a total of 47 line transitions.
Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out

