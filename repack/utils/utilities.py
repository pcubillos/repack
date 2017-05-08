# Copyright (c) 2017 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys, os
import re
import zipfile, bz2

import numpy as np
import scipy.constants as sc

from .. import constants as c

topdir = os.path.realpath(
            os.path.dirname(os.path.realpath(__file__)) + "/../..") + "/"

__all__ = ["parse_file", "read_pf", "read_states", "read_lbl", "read_iso"]


def fopen(filename):
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
    return bz2.BZ2File(filename, "r")
  elif filename.endswith(".zip"):
    zfile = zipfile.ZipFile(filename, "r")
    # Zip files have to be unzipped to seek them:
    fname = zfile.extract(zfile.namelist()[0])
    zfile.close()
    return open(fname, "r")
  else:
    return open(filename, "r")


def fclose(filename):
  """
  Remove temporary files.

  Parameters
  ----------
  filename: String
     File name.
  """
  if filename.endswith(".zip"):
    # Remove unzipped files:
    with zipfile.ZipFile(filename, "r") as zfile:
      os.remove(zfile.namelist()[0])


def parse_file(lblfile, dbtype):
  """
  Extract info from an Exomol line-transition filename.

  Parameters
  ----------
  lblfile: String
     An Exomol trans file.
  dbtype: String
     Database type (hitran or exomol).

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
    s = file.split("_")[0].split("-")
    molecule = ""
    isotope  = ""
    for i in np.arange(len(s)):
      match = re.match(r"([0-9]+)([a-z]+)([0-9]*)", s[i], re.I)
      N = 1 if match.group(3) == "" else int(match.group(3))
      molecule += match.group(2) + match.group(3)
      isotope  += match.group(1)[-1:] * N
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


def read_lbl(lblfile, dbtype, elow=None, g=None):
  """
  Read an Exomol line-transition file.

  Parameters
  ----------
  lblfile: String
     An Exomol trans filename.
  dbtype: String
     Database type (hitran or exomol).
  elow: 1D float ndarray
     The states energy (cm-1).
  g: 1D float ndarray
     The states total statistical degeneracy.

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
  f = fopen(lblfile)
  # Calculate file size:
  lenline = len(f.readline())
  f.seek(0,2)
  nlines = int(f.tell()/lenline)
  f.seek(0)

  # Extract info:
  if dbtype == "exomol":
    iup = np.zeros(nlines, int)
    ilo = np.zeros(nlines, int)
    A21 = np.zeros(nlines, np.double)
    for i in np.arange(nlines):
      line = f.readline()
      iup[i] = line[ 0:12]
      ilo[i] = line[13:25]
      A21[i] = line[26:36]
    # Compute values:
    wn   = elow[iup-1] - elow[ilo-1]
    gf   = g   [ilo-1] * A21 * c.C1 / (8.0*np.pi*100*sc.c) / wn**2.0
    Elow = elow[ilo-1]
    isoID = None

  elif dbtype == "hitran":
    isoID = np.zeros(nlines,       int)
    wn    = np.zeros(nlines, np.double)
    A21   = np.zeros(nlines, np.double)
    Elow  = np.zeros(nlines, np.double)
    g     = np.zeros(nlines, np.double)
    for i in np.arange(nlines):
      line = f.readline()
      isoID[i] = line[  2: 3]
      wn   [i] = line[  3:15]
      A21  [i] = line[ 25:35]
      Elow [i] = line[ 45:55]
      g    [i] = line[155:lenline]
    gf = g * A21 * c.C1 / (8.0*np.pi*100*sc.c) / wn**2.0
    isoID = (isoID - 1) % 10

  fclose(lblfile)
  f.close()
  return gf, Elow, wn, isoID


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

  if dbtype == "exomol":
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
