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
  14N-1H3__MiSSinG__00200-00300.trans.bz2""" in capfd.out


def test_exomol_single(capfd):
    subprocess.call(['repack', 'exomol_repack_single.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: '14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       73.90%,    270,761/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.20%,    288,406/ 1,037,545 lines.
  Total compression rate: 70.61%,    304,908/ 1,037,545 lines.

With a threshold strength factor of 0.01,
kept a total of 304,908 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files(capfd):
    subprocess.call(['repack', 'exomol_repack_two_files.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: '14N-1H3__BYTe__00200-00300.trans.bz2'.
Reading: '14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       73.90%,    270,761/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.20%,    288,406/ 1,037,545 lines.
  Total compression rate: 70.61%,    304,908/ 1,037,545 lines.

Reading: '14N-1H3__BYTe__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       83.14%,    179,161/ 1,062,896 lines.
  Flagging lines at  700 K:
  Compression rate:       81.82%,    193,215/ 1,062,896 lines.
  Total compression rate: 80.62%,    206,002/ 1,062,896 lines.

With a threshold strength factor of 0.01,
kept a total of 510,910 line transitions out of 2,100,441 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_isotopes(capfd):
    subprocess.call(['repack', 'exomol_repack_two_isotopes.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: '15N-1H3__BYTe-15__00100-00200.trans.bz2'.
Reading: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: '15N-1H3__BYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       76.63%,    283,467/ 1,212,878 lines.
  Flagging lines at  700 K:
  Compression rate:       75.06%,    302,443/ 1,212,878 lines.
  Total compression rate: 73.26%,    324,263/ 1,212,878 lines.

With a threshold strength factor of 0.01,
kept a total of 324,263 line transitions out of 1,212,878 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_two_files_two_iso(capfd):
    subprocess.call(['repack', 'exomol_repack_two_two.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Unzipping: '15N-1H3__BYTe-15__00100-00200.trans.bz2'.
Unzipping: '14N-1H3__BYTe__00200-00300.trans.bz2'.
Unzipping: '15N-1H3__BYTe-15__00200-00300.trans.bz2'.
Reading: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: '15N-1H3__BYTe-15__00100-00200.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       76.63%,    283,467/ 1,212,878 lines.
  Flagging lines at  700 K:
  Compression rate:       75.06%,    302,443/ 1,212,878 lines.
  Total compression rate: 73.26%,    324,263/ 1,212,878 lines.

Reading: '14N-1H3__BYTe__00200-00300.trans.bz2'.
Reading: '15N-1H3__BYTe-15__00200-00300.trans.bz2'.
  Flagging lines at  500 K:
  Compression rate:       85.02%,    184,887/ 1,234,214 lines.
  Flagging lines at  700 K:
  Compression rate:       83.77%,    200,355/ 1,234,214 lines.
  Total compression rate: 82.50%,    216,004/ 1,234,214 lines.

With a threshold strength factor of 0.01,
kept a total of 540,267 line transitions out of 2,447,092 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_033-100um_500-700K_lbl.dat' and
  'NH3_exomol_033-100um_500-700K_continuum.dat'.""" in capfd.out


def test_exomol_single_unzip(capfd):
    # Unzip files before repacking:
    subprocess.call(['bzip2','-dk','14N-1H3__BYTe__00100-00200.trans.bz2'])
    subprocess.call(['bzip2','-dk','14N-1H3__BYTe.states.bz2'])
    subprocess.call(['repack', 'exomol_repack_single_unzip.cfg'])
    capfd = capfd.readouterr()
    assert """Reading: '14N-1H3__BYTe__00100-00200.trans'.
  Flagging lines at  500 K:
  Compression rate:       73.90%,    270,761/ 1,037,545 lines.
  Flagging lines at  700 K:
  Compression rate:       72.20%,    288,406/ 1,037,545 lines.
  Total compression rate: 70.61%,    304,908/ 1,037,545 lines.

With a threshold strength factor of 0.01,
kept a total of 304,908 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out
    os.remove('14N-1H3__BYTe__00100-00200.trans')
    os.remove('14N-1H3__BYTe.states')


def test_exomol_single_chunks(capfd):
    subprocess.call(['repack', 'exomol_repack_single_chunks.cfg'])
    capfd = capfd.readouterr()
    assert """Unzipping: '14N-1H3__BYTe__00100-00200.trans.bz2'.
Reading: '14N-1H3__BYTe__00100-00200.trans.bz2'.
  Flagging lines at  500 K (chunk 1/2):
  Compression rate:       70.42%,    153,438/   518,772 lines.
  Flagging lines at  700 K:
  Compression rate:       68.69%,    162,411/   518,772 lines.
  Total compression rate: 66.97%,    171,335/   518,772 lines.

  Flagging lines at  500 K (chunk 2/2):
  Compression rate:       77.38%,    117,324/   518,773 lines.
  Flagging lines at  700 K:
  Compression rate:       75.71%,    125,996/   518,773 lines.
  Total compression rate: 74.25%,    133,574/   518,773 lines.

With a threshold strength factor of 0.01,
kept a total of 304,909 line transitions out of 1,037,545 lines.

Successfully rewriten exomol line-transition info into:
  'NH3_exomol_050-100um_500-700K_lbl.dat' and
  'NH3_exomol_050-100um_500-700K_continuum.dat'.""" in capfd.out

