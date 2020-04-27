import os
import subprocess

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
os.chdir(ROOT+'tests')


def test_hitemp_single_zip(capfd):
    subprocess.call('repack hitemp_repack_single_zip.cfg'.split())
    capfd = capfd.readouterr()
    assert """Reading: 'data/02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       76.11%,     51,071/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.97%,     68,468/   213,769 lines.
  Total compression rate: 66.26%,     72,118/   213,769 lines.

With a threshold strength factor of 0.01,
kept a total of 72,118 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out


def test_hitemp_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call('unzip data/02_03750-04000_HITEMP2010.zip -d data/'.split())
    subprocess.call('repack hitemp_repack_single_unzip.cfg'.split())
    capfd = capfd.readouterr()
    assert """Reading: 'data/02_3750-4000_HITEMP2010.par'.
  Flagging lines at  500 K:
  Compression rate:       76.11%,     51,071/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.97%,     68,468/   213,769 lines.
  Total compression rate: 66.26%,     72,118/   213,769 lines.

With a threshold strength factor of 0.01,
kept a total of 72,118 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out
    # Teardown:
    os.remove('data/02_3750-4000_HITEMP2010.par')


def test_hitemp_two_files(capfd):
    subprocess.call(['repack', 'hitemp_repack_two.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: 'data/02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       76.11%,     51,071/   213,769 lines.
  Flagging lines at  700 K:
  Compression rate:       67.97%,     68,468/   213,769 lines.
  Total compression rate: 66.26%,     72,118/   213,769 lines.

Reading: 'data/02_04000-04500_HITEMP2010.zip'.
  Flagging lines at  500 K:
  Compression rate:       46.12%,     74,492/   138,258 lines.
  Flagging lines at  700 K:
  Compression rate:       18.91%,    112,116/   138,258 lines.
  Total compression rate: 18.65%,    112,469/   138,258 lines.

With a threshold strength factor of 0.01,
kept a total of 184,587 line transitions out of 352,027 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.2-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.2-2.6um_500-700K_continuum.dat'.""" in capfd.out


def test_hitemp_single_chunks(capfd):
    subprocess.call('repack hitemp_repack_single_chunks.cfg'.split())
    capfd = capfd.readouterr()
    assert """Reading: 'data/02_03750-04000_HITEMP2010.zip'.
  Flagging lines at  500 K (chunk 1/3):
  Compression rate:       83.40%,     11,826/    71,256 lines.
  Flagging lines at  700 K:
  Compression rate:       78.32%,     15,451/    71,256 lines.
  Total compression rate: 76.31%,     16,883/    71,256 lines.

  Flagging lines at  500 K (chunk 2/3):
  Compression rate:       77.27%,     16,199/    71,256 lines.
  Flagging lines at  700 K:
  Compression rate:       68.75%,     22,271/    71,256 lines.
  Total compression rate: 66.94%,     23,554/    71,256 lines.

  Flagging lines at  500 K (chunk 3/3):
  Compression rate:       67.64%,     23,057/    71,257 lines.
  Flagging lines at  700 K:
  Compression rate:       56.82%,     30,771/    71,257 lines.
  Total compression rate: 55.50%,     31,706/    71,257 lines.

With a threshold strength factor of 0.01,
kept a total of 72,143 line transitions out of 213,769 lines.

Successfully rewriten hitran line-transition info into:
  'CO2_hitran_2.5-2.6um_500-700K_lbl.dat' and
  'CO2_hitran_2.5-2.6um_500-700K_continuum.dat'.""" in capfd.out

