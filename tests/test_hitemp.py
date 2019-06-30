import os
import sys
import pytest
import subprocess

import numpy as np

#ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'

#os.chdir(ROOT+'tests')


def test_hitemp_single_zip(capfd):
    subprocess.call(['python', '../repack.py', 'hitemp_repack_single_zip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '02_3750-4000_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       76.13%,     51,032/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.99%,     68,431/   213,769 lines.
  Total compression rate: 66.28%,     72,084/   213,769 lines.

Kept a total of 72,084 line transitions.
Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out


@pytest.mark.skip
def test_hitemp_two_files(capfd):
    subprocess.call(['python','../repack.py','hitemp_repack_two.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '01_100-125_MockHITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       93.00%,          7/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       92.00%,          8/       100 lines.
  Total compression rate: 92.00%,          8/       100 lines.

Reading: '01_125-150_MockHITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       94.00%,          6/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       93.00%,          7/       100 lines.
  Total compression rate: 93.00%,          7/       100 lines.

Kept a total of 15 line transitions.
Successfully rewriten hitran line-transition info into:
  'H2O_hitran_050-100um_500-700K_lbl.dat' and
  'H2O_hitran_050-100um_500-700K_continuum.dat'.""" in capfd.out


@pytest.mark.skip
def test_hitemp_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call(['unzip', '01_100-125_MockHITEMP2010.zip'])
    subprocess.call(['python','../repack.py','hitemp_repack_single_unzip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '01_100-125_MockHITEMP2010.par'.
  Flagging lines at  500 K:
  Compression rate:       93.00%,          7/       100 lines.
  Flagging lines at  700 K:
  Compression rate:       92.00%,          8/       100 lines.
  Total compression rate: 92.00%,          8/       100 lines.

Kept a total of 8 line transitions.
Successfully rewriten hitran line-transition info into:
  'H2O_hitran_080-100um_500-700K_lbl.dat' and
  'H2O_hitran_080-100um_500-700K_continuum.dat'.""" in capfd.out
    # Teardown:
    os.remove('01_100-125_MockHITEMP2010.par')


@pytest.mark.skip
def test_hitemp_single_chunks(capfd):
    subprocess.call(['python','../repack.py','hitemp_repack_single_chunks.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '01_100-125_MockHITEMP2010.zip'.
  Flagging lines at  500 K (chunk 1/2):
  Compression rate:       86.00%,          7/        50 lines.
  Flagging lines at  700 K:
  Compression rate:       78.00%,         11/        50 lines.
  Total compression rate: 78.00%,         11/        50 lines.

  Flagging lines at  500 K (chunk 2/2):
  Compression rate:       94.00%,          3/        50 lines.
  Flagging lines at  700 K:
  Compression rate:       94.00%,          3/        50 lines.
  Total compression rate: 94.00%,          3/        50 lines.

Kept a total of 14 line transitions.
Successfully rewriten hitran line-transition info into:
  'H2O_hitran_080-100um_500-700K_lbl.dat' and
  'H2O_hitran_080-100um_500-700K_continuum.dat'.""" in capfd.out

