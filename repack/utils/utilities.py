import sys, os
import re
import numpy as np
import scipy.constants as sc

from .. import constants as c

topdir = os.path.realpath(
            os.path.dirname(os.path.realpath(__file__)) + "/../..") + "/"

__all__ = ["parse_file", "read_pf", "read_states", "read_lbl", "read_iso"]

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
  if dbtype == "exomol":
    root, file = os.path.split(os.path.realpath(lblfile))
    # Auxilliary files:
    sfile = file.replace("trans", "states")
    if sfile.count("__") == 2:
      suffix = sfile[sfile.rindex("__"):sfile.index(".")]
      sfile = sfile.replace(suffix, "")
    else:
      suffix = ""
    sfile  = root + "/" + sfile
    pffile = sfile.replace("states", "pf")

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
    pass

  return suffix, molecule, isotope, pffile, sfile


def read_pf(pffile):
  """
  Read an Exomol partition-function file.

  Parameters
  ----------
  pffile: String

  Returns
  -------
  temp: 1D float ndarray
     Tabulated list of temperaures.
  pf: 1D float ndarray
     Partition-function values for each temperature.
  """
  # Read partition-function file:
  with open(pffile, "r") as f:
    lines = f.readlines()
  # Alocate outputs:
  ntemp = len(lines)
  temp = np.zeros(ntemp, np.double)
  pf   = np.zeros(ntemp, np.double)
  # Extract data:
  for i in np.arange(ntemp):
    temp[i], pf[i] = lines[i].split()
  return temp, pf


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
  with open(states, "r") as f:
    lines = f.readlines()
  nstates = len(lines)
  # Alocate outputs:
  elow    = np.zeros(nstates, np.double)  # State energy
  g       = np.zeros(nstates, int)        # State degeneracy (incl. ns)
  # Extract data:
  for i in np.arange(nstates):
    elow[i], g[i] = lines[i].split()[1:3]
  return elow, g


def read_lbl(lblfile, elow, g):
  """
  Read an Exomol line-transition file.

  Parameters
  ----------
  lblfile: String
     An Exomol trans filename.
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
  """
  with open(lblfile, "r") as f:
    # Calculate file size:
    lenline = len(f.readline())
    f.seek(0,2)
    nlines = int(f.tell()/lenline)
    # Extract info:
    f.seek(0)
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
  return gf, Elow, wn


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
