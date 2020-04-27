import os
import subprocess
import pytest

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
os.chdir(ROOT+'tests')


@pytest.mark.parametrize('missfile',
    ['exomol_repack_missing1.cfg',
     'exomol_repack_missing2.cfg'])
def test_missing(capfd, missfile):
    subprocess.call(['repack', missfile])
    capfd = capfd.readouterr()
    assert """File(s) not Found Error: These files are missing:
  data/14N-1H3__MiSSinG__00200-00300.trans.bz2""" in capfd.out


def test_exomol_single(capfd):
    subprocess.call('repack exomol_repack_single.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       74.01%,    269,648/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.35%,    286,925/ 1,037,545 lines.
  Total compression rate: 70.75%,    303,456/ 1,037,545 lines.

With a threshold strength factor of 0.01,
kept a total of 303,456 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files(capfd):
    subprocess.call('repack exomol_repack_two_files.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: 'data/14N-1H3__BYTe__00200-00300.trans.bz2'.
Reading: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       74.01%,    269,648/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.35%,    286,925/ 1,037,545 lines.
  Total compression rate: 70.75%,    303,456/ 1,037,545 lines.

Reading: 'data/14N-1H3__BYTe__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       83.19%,    178,663/ 1,062,896 lines.
  Flagging lines at  700 K:
  Compression rate:       81.90%,    192,354/ 1,062,896 lines.
  Total compression rate: 80.69%,    205,243/ 1,062,896 lines.

With a threshold strength factor of 0.01,
kept a total of 508,699 line transitions out of 2,100,441 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_isotopes(capfd):
    subprocess.call('repack exomol_repack_two_isotopes.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: 'data/15N-1H3__BYTe-15__00100-00200.trans.bz2'.
Reading: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: 'data/15N-1H3__BYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       76.73%,    282,278/ 1,212,878 lines.
  Flagging lines at  700 K:
  Compression rate:       75.20%,    300,854/ 1,212,878 lines.
  Total compression rate: 73.40%,    322,684/ 1,212,878 lines.

With a threshold strength factor of 0.01,
kept a total of 322,684 line transitions out of 1,212,878 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files_two_iso(capfd):
    subprocess.call('repack exomol_repack_two_two.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: 'data/15N-1H3__BYTe-15__00100-00200.trans.bz2'.
Unzipping: 'data/14N-1H3__BYTe__00200-00300.trans.bz2'.
Unzipping: 'data/15N-1H3__BYTe-15__00200-00300.trans.bz2'.
Reading: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: 'data/15N-1H3__BYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       76.73%,    282,278/ 1,212,878 lines.
  Flagging lines at  700 K:
  Compression rate:       75.20%,    300,854/ 1,212,878 lines.
  Total compression rate: 73.40%,    322,684/ 1,212,878 lines.

Reading: 'data/14N-1H3__BYTe__00200-00300.trans.bz2'.
Reading: 'data/15N-1H3__BYTe-15__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       85.06%,    184,349/ 1,234,214 lines.
  Flagging lines at  700 K:
  Compression rate:       83.84%,    199,462/ 1,234,214 lines.
  Total compression rate: 82.57%,    215,138/ 1,234,214 lines.

With a threshold strength factor of 0.01,
kept a total of 537,822 line transitions out of 2,447,092 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call(['bzip2','-dk','data/14N-1H3__BYTe__00100-00200.trans.bz2'])
    subprocess.call('bzip2 -dk data/14N-1H3__BYTe.states.bz2'.split())
    subprocess.call('repack exomol_repack_single_unzip.cfg'.split())
    capfd = capfd.readouterr()
    assert """Reading: 'data/14N-1H3__BYTe__00100-00200.trans'.
  Flagging lines at  500 K:
  Compression rate:       74.01%,    269,648/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.35%,    286,925/ 1,037,545 lines.
  Total compression rate: 70.75%,    303,456/ 1,037,545 lines.

With a threshold strength factor of 0.01,
kept a total of 303,456 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out
    os.remove('data/14N-1H3__BYTe__00100-00200.trans')
    os.remove('data/14N-1H3__BYTe.states')


def test_exomol_single_chunks(capfd):
    subprocess.call('repack exomol_repack_single_chunks.cfg'.split())
    capfd = capfd.readouterr()
    assert """Unzipping: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: 'data/14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K (chunk 1/2):
  Compression rate:       70.54%,    152,852/   518,772 lines.
  Flagging lines at  700 K:
  Compression rate:       68.85%,    161,620/   518,772 lines.
  Total compression rate: 67.12%,    170,554/   518,772 lines.

  Flagging lines at  500 K (chunk 2/2):
  Compression rate:       77.49%,    116,797/   518,773 lines.
  Flagging lines at  700 K:
  Compression rate:       75.85%,    125,306/   518,773 lines.
  Total compression rate: 74.38%,    132,903/   518,773 lines.

With a threshold strength factor of 0.01,
kept a total of 303,457 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out

