# Copyright (c) 2017-2018 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

__all__ = ["parse_file", "read_pf", "read_states", "lbl", "wnbalance",
           "count", "read_iso", "get_exomol_mol"]

import sys
import os
import re
import zipfile
import struct
import itertools

import numpy as np
import scipy.constants as sc

from .. import constants as c

topdir = os.path.realpath(
            os.path.dirname(os.path.realpath(__file__)) + "/../..") + "/"


def fopen(filename, mode="r"):
  """
  Find out file compression format (if any) and open the file.

  Parameters
  ----------
  filename: String
     File name.

  Returns
  -------
  file: FILE pointer
  """
  if   filename.endswith(".bz2"):
    return open(filename.replace(".bz2", ""), mode)
  elif filename.endswith(".zip"):
    zfile = zipfile.ZipFile(filename, mode)
    # Zip files have to be unzipped to seek them:
    fname = zfile.extract(zfile.namelist()[0])
    zfile.close()
    return open(fname, mode)
  else:
    return open(filename, mode)


def parse_file(lblfile, dbtype):
  """
  Extract info from an Exomol line-transition filename.

  Parameters
  ----------
  lblfile: String
     An Exomol trans file.
  dbtype: String
     Database type (hitran, exomol, or kurucz).

  Returns
  -------
  suffix: Strings
     Suffix of the trans file.
  molecule: String
     The molecule name.
  isotope: String
     The isotope name.
  pffile: String
     Partition-function file.
  sfile: String
     States file.
  """
  root, file = os.path.split(os.path.realpath(lblfile))
  if dbtype == "exomol":
    # Auxilliary files:
    sfile = file.replace("trans", "states")
    if sfile.count("__") == 2:
      suffix = sfile[sfile.rindex("__"):sfile.index(".")]
      sfile = sfile.replace(suffix, "")
    else:
      suffix = ""
    sfile  = root + "/" + sfile
    pffile = sfile.replace("states", "pf").strip(".bz2")
    # Get info from file name:
    molecule, isotope = get_exomol_mol(file)

  elif dbtype == "hitran":
    hitempID = {"01":"H2O", "02":"CO2", "05":"CO", "08":"NO"}
    molID = file[0:2]
    molecule = hitempID[molID]
    suffix = file[file.find("_")+1:file.rfind(".par")]
    if suffix.find("-") > 0:
      irange = suffix[0:suffix.find("-")]
      suffix = irange.zfill(5)
    # Get them later from the input PF file (required for HITEMP):
    isotope = None
    pffile  = None
    sfile   = None

  elif dbtype == "kurucz":  # TiO is the only available molecule so far
    molecule = "TiO"
    suffix  = ""
    isotope = None
    pffile  = None
    sfile   = None
  return suffix, molecule, isotope, pffile, sfile


def read_pf(pffile, dbtype="exomol"):
  """
  Read an Exomol partition-function file.

  Parameters
  ----------
  pffile: String
     Partition-function file name.
  dbtype: String
     Database type (hitran or exomol).

  Returns
  -------
  temp: 1D float ndarray
     Tabulated list of temperaures.
  pf: 1D float ndarray
     Partition-function values for each temperature.
  isotopes: 1D float ndarray
     The isotope names [returned only for pyrat dbtype].
  """
  # Read partition-function file:
  with fopen(pffile) as f:
    lines = f.readlines()

  if dbtype == "exomol":
    # Alocate outputs:
    ntemp = len(lines)
    temp = np.zeros(ntemp, np.double)
    pf   = np.zeros(ntemp, np.double)
    # Extract data:
    for i in np.arange(ntemp):
      temp[i], pf[i] = lines[i].split()
    return temp, pf

  if dbtype == "pyrat":
    # Number of header lines (to skip later when reading the tabulated data):
    nskip = 0
    while True:
      line = lines[nskip].strip()
      # Skip blank/empty lines:
      if line == "" or line.startswith('#'):
        pass
      # Read isotopes:
      elif line == "@ISOTOPES":
        isotopes = np.asarray(lines[nskip+1].strip().split())
      # Stop when the tabulated data begins:
      if line == "@DATA":
        nskip += 1
        break
      nskip += 1

    # Number of isotopes:
    niso = len(isotopes)
    # Number of temperature samples:
    ntemp = len(lines) - nskip

    # Allocate output arrays:
    temp = np.zeros(ntemp, np.double)
    pf   = np.zeros((niso, ntemp), np.double)
    # Extract the data:
    for i in np.arange(ntemp):
      info = lines[nskip+i].strip().split()
      temp[i] = info[0]
      pf[:,i] = info[1:]
    return temp, pf, isotopes


def read_states(states):
  """
  Read an Exomol states file.

  Parameters
  ----------
  states: String
     An Exomol states filename.

  Returns
  -------
  elow: 1D float ndarray
     State energy (cm-1).
  g: 1D integer ndarray
     State total statistical degeneracy.
  """
  # Read states file:
  with fopen(states) as f:
    lines = f.readlines()
  nstates = len(lines)
  # Alocate outputs:
  elow    = np.zeros(nstates, np.double)  # State energy
  g       = np.zeros(nstates, int)        # State degeneracy (incl. ns)
  # Extract data:
  for i in np.arange(nstates):
    elow[i], g[i] = lines[i].split()[1:3]
  return elow, g


class lbl():
  def __init__(self, lblfile, dbtype, elow, g, iso):
    """
    Parameters
    ----------
    lblfile: String
       A line-transition file.
    dbtype: String
       Database type (hitran or exomol).
    elow: 1D float ndarray
       The states energy (cm-1).
    g: 1D float ndarray
       The states total statistical degeneracy.
    iso: Integer
       The isotope index for this file (for Exomol data).
    """
    self.lblfile = lblfile                    # The file name
    self.dbtype  = dbtype                     # Database type

    if dbtype == "kurucz":
      self.file  = fopen(lblfile, "rb")
      self.llen  = 16
      self.ratiolog = np.log(1.0 + 1.0/2000000)
      self.tablog   = 10.0**(0.001*(np.arange(32769) - 16384))
    else:
      self.file  = fopen(lblfile, "r")  # The actual file
      dummy = self.file.readline()
      self.llen = self.file.tell()      # Length of lines in file
    self.file.seek(0,2)
    self.nlines  = int(self.file.tell()//self.llen)  # Number of lines
    self.elow    = elow
    self.g       = g
    self.iso     = iso  # Isotope index

  def bs(self, val, lo, hi):
    """
    Wavenumber binary search on database.

    Parameters
    ----------
    val: Float
       Target wavenumber value (in cm-1).
    lo: Integer
       Initial wavenumber index where to search.
    hi: Integer
       Final wavenumber index where to search.

    Returns
    -------
    index: Integer
       The index of the closest wavenumber entry to val.
    """
    # Out of bounds:
    if val <= self.getwn(0):
      return 0
    if val >= self.getwn(self.nlines-1):
      return self.nlines-1

    # End case:
    if hi-lo <= 1:
      if np.abs(val-self.getwn(hi)) < np.abs(val-self.getwn(lo)):
        return hi
      return lo

    # Half it:
    mid = int(0.5*(hi+lo))
    if self.getwn(mid) > val:
      return self.bs(val, lo, mid)
    return self.bs(val, mid, hi)


  def getwn(self, index):
    """
    Extract wavenumber from file at the requested index.

    Parameters
    ----------
    index: Integer
       The wavenumber index in database to extract.
    Returns
    -------
    wn: Float
       The wavenumber (cm-1) at position index.
    """
    if   self.dbtype == "exomol":
      self.file.seek(index*self.llen)
      iup = int(self.file.read(12)) - 1
      self.file.read(1)
      ilo = int(self.file.read(12)) - 1
      return self.elow[iup] - self.elow[ilo]

    elif self.dbtype == "hitran":
      self.file.seek(index*self.llen + 3)
      return float(self.file.read(12))

    elif self.dbtype == "kurucz":
      # Note I'm reversing the indexing because this DB has decreasing wn:
      self.file.seek((self.nlines-index-1)*self.llen)
      # 4 = struct.calcsize("i")
      iw = struct.unpack('i', self.file.read(4))[0]
      return 1.0/(np.exp(iw * self.ratiolog) * c.nano)


  def read(self, chunk):
    """
    Read a chunk of line transitions.

    Parameters
    ----------
    chunk: Two-element tuple
       Initial and final indices of the chunk to read.

    Returns
    -------
    gf: 1D float ndarray
       Transition weighted oscillator strength (unitless).
    Elow: 1D float ndarray
       Transition lower-state energy (cm-1).
    wn: 1D float ndarray
       Transition wavenumber (cm-1).
    isoID: 1D integer ndarray
       isotope ID (for hitran dbtype).
    """
    # Go to beginning of chunk:
    self.file.seek(chunk[0]*self.llen)
    nlines = chunk[1] - chunk[0]
    # Extract info:
    if self.dbtype == "exomol":
      iup = np.zeros(nlines, int)
      ilo = np.zeros(nlines, int)
      A21 = np.zeros(nlines, np.double)
      for i in range(nlines):
          line = self.file.readline()
          iup[i] = line[ 0:12]
          ilo[i] = line[13:25]
          A21[i] = line[26:36]
      iup -= 1
      ilo -= 1
      # Compute values:
      wn   = self.elow[iup] - self.elow[ilo]
      gf   = self.g[ilo] * A21 * c.C1 / (8.0*np.pi*100*sc.c) / wn**2.0
      Elow = self.elow[ilo]
      isoID = np.tile(self.iso, np.size(wn))

    elif self.dbtype == "hitran":
      isoID = np.zeros(nlines,       int)
      wn    = np.zeros(nlines, np.double)
      A21   = np.zeros(nlines, np.double)
      Elow  = np.zeros(nlines, np.double)
      g     = np.zeros(nlines, np.double)
      for i in np.arange(nlines):
        line = self.file.readline()
        isoID[i] = line[  2: 3]
        wn   [i] = line[  3:15]
        A21  [i] = line[ 25:35]
        Elow [i] = line[ 45:55]
        g    [i] = line[155:self.llen]
      gf = g * A21 * c.C1 / (8.0*np.pi*100*sc.c) / wn**2.0
      isoID = (isoID - 1) % 10

    elif self.dbtype == "kurucz":
      iw   = np.zeros(nlines, int)
      ieli = np.zeros(nlines, np.short)
      ielo = np.zeros(nlines, np.short)
      igf  = np.zeros(nlines, np.short)
      i = 0
      while i < nlines:
        idx = self.nlines - 1 - (chunk[0]+i)  # Reverse indexing again
        self.file.seek(idx*self.llen)
        # 10 = struct.calcsize("ihhh")
        iw[i], ieli[i], ielo[i], igf[i] = struct.unpack('ihhh',
                                                        self.file.read(10))
        i += 1
      wn    = 1.0/(np.exp(iw * self.ratiolog) * c.nano)
      gf    = self.tablog[igf]
      isoID = np.abs(ieli) - 8950
      Elow  = self.tablog[ielo]

    return gf, Elow, wn, isoID


  def close(self):
    if self.lblfile.endswith(".zip"):
      # Remove unzipped files:
      with zipfile.ZipFile(self.lblfile, "r") as zfile:
        os.remove(zfile.namelist()[0])
    self.file.close()


def wnbalance(lbls, wnmin, wnmax, targetsize, zero=0, tol=0.01):
  """
  Binary search of wavenumber value (wntarget) such that there are
  targetsize+zero line-transitions with wavenumber < wntarget.

  Parameters
  ----------
  lbls: List of lbl objects
     Line by line objects.
  wnmin: Float
     Minimum wavenumber to consider
  wnmax: Float
     Maximum wavenumber to consider
  targetsze: Integer
     The desired number of transitions with wavenumber < wntarget.
  zero: Integer
     Zero offset number of transitions to discount.
  tol: Float
     Fractional tolerance level (w.r.t. targetsize) of to accept wntarget.

  Returns
  -------
  wntarget: Float
     The wavenumber that has targetsize transitions to the left.
  """
  # Middle point between wn boundaries:
  wntarget = 0.5*(wnmax + wnmin)
  # Number of transitions with wn < wntarget:
  ntarget  = count(lbls, wntarget)

  if np.abs(ntarget-zero - targetsize) < tol*targetsize:
    return wntarget

  elif ntarget-zero < targetsize:
    return wnbalance(lbls, wntarget, wnmax,    targetsize, zero, tol)
  else:
    return wnbalance(lbls, wnmin,    wntarget, targetsize, zero, tol)


def count(lbls, wntarget):
  """
  Count the number of line transitions with wavenumber smaller than wntarget.

  Parameters
  ----------
  lbls: List of lbl objects
  wntarget: Float
     The target wavenumber.

  Returns
  -------
  nwave: Integer
     The number of transitions with wavenumber smaller than wntarget.
  """
  nwave = 0
  for k in np.arange(len(lbls)):
    nwave += lbls[k].bs(wntarget, 0, lbls[k].nlines-1)

  return nwave


def read_iso(mol, iso, dbtype="exomol", isofile=topdir+"inputs/isotopes.dat"):
  """
  Read an isotopes info file.

  Parameters
  ----------
  mol: String
     Molecule name.
  iso: List of strings
     Molecule's isotope name.
  dbtype: String
     Database format for isotope names (hitran or exomol).
  isofile: String
     File containing the isotopic information.

  Returns
  -------
  iratio: List of floats
    Isotopic abudance fraction (unitless).
  imass: List of floats
    Isotopic mass (amu).
  """
  # Alocate outputs:
  iratio = np.zeros(len(iso))
  imass  = np.zeros(len(iso))

  if dbtype in ["exomol", "kurucz"]:  # Kurucz's TiO uses Exomol iso notation
    iiso = 2
  elif dbtype == "hitran":
    iiso = 1
  # Read info file:
  with open(isofile, "r") as f:
    lines = f.readlines()
  # Get values for our molecule/isotopes:
  for i in np.arange(len(lines)):
    if lines[i].startswith("#") or lines[i].strip() == "":
      continue
    info = lines[i].split()
    if info[0] == mol:
      if info[iiso] in iso:
        iratio[iso.index(info[iiso])] = info[3]
        imass [iso.index(info[iiso])] = info[4]
  return iratio, imass


def get_exomol_mol(dbfile):
  """
  Parse an exomol file to extract the molecule and isotope name.

  Parameters
  ----------
  dbfile: String
      An exomol line-list file (must follow ExoMol naming convention).

  Returns
  -------
  molecule: String
      Name of the molecule.
  isotope: String
      Name of the isotope (See Tennyson et al. 2016, jmosp, 327).

  Examples
  --------
  >>> filenames = [
  >>>     '1H2-16O__POKAZATEL__00400-00500.trans.bz2',
  >>>     '1H-2H-16O__VTT__00250-00500.trans.bz2',
  >>>     '12C-16O2__HITEMP.pf',
  >>>     '12C-16O-18O__Zak.par',
  >>>     '12C-1H4__YT10to10__01100-01200.trans.bz2',
  >>>     '12C-1H3-2H__MockName__01100-01200.trans.bz2'
  >>>    ]
  >>> for db in filenames:
  >>>     print(get_exomol_mol(db))
  ('H2O', '116')
  ('H2O', '126')
  ('CO2', '266')
  ('CO2', '268')
  ('CH4', '21111')
  ('CH4', '21112')
  """
  atoms = os.path.split(dbfile)[1].split("_")[0].split("-")
  elements = []
  isotope  = ""
  for atom in atoms:
      match = re.match(r"([0-9]+)([a-z]+)([0-9]*)", atom, re.I)
      N = 1 if match.group(3) == "" else int(match.group(3))
      elements += N * [match.group(2)]
      isotope  += match.group(1)[-1:] * N

  composition = [list(g[1]) for g in itertools.groupby(elements)]
  molecule = "".join([c[0] + str(len(c))*(len(c)>1)
                      for c in composition])

  return molecule, isotope
