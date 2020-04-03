import os
import subprocess

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
os.chdir(ROOT+'tests')


def test_hitemp_single_zip(capfd):
    subprocess.call(['repack', 'hitemp_repack_single_zip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       76.13%,     51,032/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.99%,     68,431/   213,769 lines.
  Total compression rate: 66.28%,     72,084/   213,769 lines.

With a threshold strength factor of 0.01,
kept a total of 72,084 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out


def test_hitemp_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call(['unzip', '02_03750-04000_HITEMP2010.zip'])
    subprocess.call(['repack', 'hitemp_repack_single_unzip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '02_3750-4000_HITEMP2010.par'.
  Flagging lines at  500 K:
  Compression rate:       76.13%,     51,032/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.99%,     68,431/   213,769 lines.
  Total compression rate: 66.28%,     72,084/   213,769 lines.

With a threshold strength factor of 0.01,
kept a total of 72,084 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out
    # Teardown:
    os.remove('02_3750-4000_HITEMP2010.par')


def test_hitemp_two_files(capfd):
    subprocess.call(['repack', 'hitemp_repack_two.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       76.13%,     51,032/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.99%,     68,431/   213,769 lines.
  Total compression rate: 66.28%,     72,084/   213,769 lines.

Reading: '02_04000-04500_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       46.14%,     74,466/   138,258 lines.
  Flagging lines at  700 K:
  Compression rate:       18.99%,    111,996/   138,258 lines.
  Total compression rate: 18.73%,    112,367/   138,258 lines.

With a threshold strength factor of 0.01,
kept a total of 184,451 line transitions out of 352,027 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.2-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.2-2.6um_500-700K_continuum.dat'.""" in capfd.out


def test_hitemp_single_chunks(capfd):
    subprocess.call(['repack', 'hitemp_repack_single_chunks.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K (chunk 1/3):
  Compression rate:       83.40%,     11,830/    71,256 lines.
  Flagging lines at  700 K:
  Compression rate:       78.32%,     15,446/    71,256 lines.
  Total compression rate: 76.31%,     16,884/    71,256 lines.

  Flagging lines at  500 K (chunk 2/3):
  Compression rate:       77.32%,     16,163/    71,256 lines.
  Flagging lines at  700 K:
  Compression rate:       68.80%,     22,235/    71,256 lines.
  Total compression rate: 67.00%,     23,517/    71,256 lines.

  Flagging lines at  500 K (chunk 3/3):
  Compression rate:       67.65%,     23,050/    71,257 lines.
  Flagging lines at  700 K:
  Compression rate:       56.81%,     30,775/    71,257 lines.
  Total compression rate: 55.50%,     31,708/    71,257 lines.

With a threshold strength factor of 0.01,
kept a total of 72,109 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out

