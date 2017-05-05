# Copyright (c) 2017 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import sys, os
import struct
import warnings

import scipy.interpolate as sip
import scipy.constants as sc
import numpy as np

# Config Parser changed between Python2 and Python3:
if sys.version_info.major == 3:
  import configparser
else:
  import ConfigParser as configparser

from . import utils     as u
from . import constants as c


def parser(cfile):
  """
  Parse input arguments from configuration file.

  Parameters
  ----------
  cfile: String
     A configuration file.

  Returns
  -------
  files: List of strings
     List with line-transition file names.
  outfile: String
     Output file root name.
  tmin: Float
     Minimum temperature to sample.
  tmax: Float
     Maximum temperature to sample.
  dtemp: Float
     Temperature sampling rate.
  wnmin: Float
     Minimum wavenumber (cm-1) to sample.
  wnmax: Float
     Maximum wavenumber (cm-1) to sample.
  dwn: Float
     Wavenumber sampling rate (cm-1).
  sthresh: Float
     Threshold tolerance level for weak/strong lines.
  pffile: String
     Input partition-function file (for HITRAN input dbtype).
  """
  config = configparser.SafeConfigParser()
  config.read([cfile])
  section = "REPACK"
  # Input line-transition files:
  lblfiles = config.get(section, "lblfiles")
  lblfiles = lblfiles.split()
  # Database type:
  db = config.get(section, "dbtype")
  # Output file:
  outfile  = config.get(section, "outfile")
  # Partition-function file:
  pffile = None
  if config.has_option(section, "pffile"):
    pffile = config.get(section, "pffile")

  # Temperature sampling:
  tmin  = config.getfloat(section, "tmin")
  tmax  = config.getfloat(section, "tmax")
  dtemp = config.getfloat(section, "dtemp")

  # Wavenumber sampling:
  wnmin = config.getfloat(section, "wnmin")
  wnmax = config.getfloat(section, "wnmax")
  dwn   = config.getfloat(section, "dwn")

  # Line-strength threshold:
  sthresh = config.getfloat(section, "sthresh")

  return lblfiles, db, outfile, tmin, tmax, dtemp, wnmin, wnmax, dwn, \
         sthresh, pffile


def repack(files, dbtype, outfile, tmin, tmax, dtemp, wnmin, wnmax, dwn,
           sthresh, pffile=None):
  """
  Re-pack line-transition data into lbl data for strong lines and
  continuum data for weak lines.

  Parameters
  ----------
  files: List of strings
     List with line-transition file names.
  dbtype: String
    Database type (exomol or hitran).
  outfile: String
     Output file root name.
  tmin: Float
     Minimum temperature to sample.
  tmax: Float
     Maximum temperature to sample.
  dtemp: Float
     Temperature sampling rate.
  wnmin: Float
     Minimum wavenumber (cm-1) to sample.
  wnmax: Float
     Maximum wavenumber (cm-1) to sample.
  dwn: Float
     Wavenumber sampling rate (cm-1).
  sthresh: Float
     Threshold tolerance level for weak/strong lines.
  pffile: String
     Input partition-function file (for HITRAN input dbtype).
  """

  # Temperature sampling:
  ntemp = int((tmax-tmin)/dtemp + 1)
  temperature = np.linspace(tmin, tmax, ntemp)

  # Wavenumber sampling:
  nwave = int((wnmax-wnmin)/dwn + 1)
  wnspec = np.linspace(wnmin, wnmax, nwave)
  continuum = np.zeros((nwave, ntemp), np.double)

  if dbtype not in ["hitran", "exomol"]:
    print("\n{:s}\n  Error: Invalid database, dbtype must be either hitran "
          "or exomol.\n{:s}\n".format(70*":", 70*":"))
    sys.exit(0)

  # Parse input files:
  nfiles = len(files)
  suff, mol, isot, pf, states = [], [],  [], [], []
  for i in np.arange(nfiles):
    s, m, iso, p, st = u.parse_file(files[i], dbtype)
    suff.append(s)
    mol.append(m)
    isot.append(iso)
    pf.append(p)
    states.append(st)

  if len(np.unique(mol)) > 1:
    print("\n{:s}\n  Error: All input files must correspont to the same "
          "molecule.\n{:s}\n".format(70*":", 70*":"))
    sys.exit(0)
  mol = mol[0]

  z = []  # Interpolator function for partition function per isotope
  if pffile is not None:
    # Read input partition-function file (if given):
    pftemp, partf, isotopes = u.read_pf(pffile, dbtype="pyrat")
    isotopes = list(isotopes)
    for j in np.arange(len(isotopes)):
      z.append(sip.interp1d(pftemp, partf[j], kind='slinear'))
  else:
    isotopes = list(np.unique(isot))
  niso = len(isotopes)

  # Isotopic abundance ratio and mass:
  iratio, imass = u.read_iso(mol, isotopes, dbtype)

  s = suff[:]  # Make a copy
  # File indices for each wavenumber set:
  wnset = []
  while np.size(s) != 0:
    suffix = s[np.argmin(s)]
    idx = []
    for i in np.arange(nfiles):
      if suff[i] == suffix:
        idx.append(i)
        s.remove(suffix)
    wnset.append(idx)

  iso = np.zeros(nfiles, int)
  if dbtype == "exomol":
    for i in np.arange(nfiles):
      iso[i] = isotopes.index(isot[i])
    # Read partition-function and states files:
    lblargs = []
    for j in np.arange(niso):
      i = isot.index(isotopes[j])
      # Partition function:
      if pffile is None:
        temp, part = u.read_pf(pf[i], dbtype)
        z.append(sip.interp1d(temp, part, kind='slinear'))
      # States:
      elow, degen = u.read_states(states[i])
      lblargs.append([elow, degen])

  elif dbtype == "hitran":
    lblargs = [[None, None]]  # Trust me

  # Turn isotopes from string to integer data type:
  isotopes = np.asarray(isotopes, int)

  # Set output file names:
  lbl_out  = "{:s}_{:s}_{:s}_lbl.dat".      format(mol, dbtype, outfile)
  cont_out = "{:s}_{:s}_{:s}_continuum.dat".format(mol, dbtype, outfile)
  # Output line-by-line file:
  lbl = open(lbl_out, "wb")

  # Read line-by-line files:
  for i in np.arange(len(wnset)):
    gf   = np.array([])
    Elow = np.array([])
    wn   = np.array([])
    iiso = np.array([], int)
    # Gather data from all files in this wavenumber range:
    for k in np.arange(len(wnset[i])):
      idx = wnset[i][k]
      j   = int(iso[idx])
      # Read the LBL file:
      print("Reading: '{:s}'.".format(files[idx]))
      gfosc, el, wnumber, isoID = u.read_lbl(files[idx], dbtype, *(lblargs[j]))
      if dbtype == "exomol":
        isoID = np.tile(j, np.size(wnumber))
      wnrange = (wnumber >= wnmin) & (wnumber <= wnmax)
      gf   = np.hstack([gf,   gfosc  [wnrange]])
      Elow = np.hstack([Elow, el     [wnrange]])
      wn   = np.hstack([wn,   wnumber[wnrange]])
      iiso = np.hstack([iiso, isoID  [wnrange]])
    # Sort by wavelength:
    asort = np.argsort(wn)
    gf   = gf  [asort]
    Elow = Elow[asort]
    wn   = wn  [asort]
    iiso = iiso[asort]

    # Low temperature line flagging:
    Z = np.zeros(niso)
    for j in np.arange(niso):
      Z[j] = z[j](tmin)
    s = (gf*iratio[iiso]/Z[iiso] *
         np.exp(-c.C2*Elow/tmin) * (1-np.exp(-c.C2*wn/tmin)))
    # Line Doppler width:
    alphad = wn/(100*sc.c) * np.sqrt(2.0*c.kB*tmin / (imass[iiso]*c.amu))
    # Line-strength sorted in descending order:
    isort  = np.argsort(s/alphad)[::-1]
    flag = np.ones(len(iiso), bool)
    print("  Flagging lines at {:4.0f} K:".format(tmin))
    # Flag strong/weak lines:
    u.flag(flag, wn, s/alphad, isort, alphad, sthresh)
    nlines = len(wn)
    print("  Compression rate:       {:.2f}%,  {:9,d}/{:10,d} lines.".
          format(np.sum(1-flag)*100./len(flag), np.sum(flag), nlines))

    # High temperature line flagging:
    Z = np.zeros(niso)
    for j in np.arange(niso):
      Z[j] = z[j](tmax)
    s = (gf*iratio[iiso]/Z[iiso] *
         np.exp(-c.C2*Elow/tmax) * (1-np.exp(-c.C2*wn/tmax)))
    alphad = wn/(100*sc.c) * np.sqrt(2.0*c.kB*tmax / (imass[iiso]*c.amu))
    isort  = np.argsort(s/alphad)[::-1]
    print("  Flagging lines at {:4.0f} K:".format(tmax))
    flag2 = np.ones(len(iiso), bool)
    u.flag(flag2, wn, s/alphad/np.sqrt(np.pi), isort, alphad, sthresh)
    print("  Compression rate:       {:.2f}%,  {:9,d}/{:10,d} lines.".
          format(np.sum(1-flag2)*100./len(flag2), np.sum(flag2), nlines))
    # Update flag:
    flag |= flag2
    print("  Total compression rate: {:.2f}%,  {:9,d}/{:10,d} lines.\n".
          format(np.sum(1-flag)*100./len(flag), np.sum(flag), nlines))

    # Store weak lines to continuum file as function of temp:
    for t in np.arange(ntemp):
      T = temperature[t]
      for j in np.arange(niso):
        Z[j] = z[j](temperature[t])
      # Line strength in cm molec-1
      s = (c.C3 * gf*iratio[iiso]/Z[iiso] *
             np.exp(-c.C2*Elow/T) * (1-np.exp(-c.C2*wn/T)))
      # Continuum opacity in cm-1 amagat-1 (instead of cm2 molec-1):
      u.continuum(s[~flag], wn[~flag], continuum[:,t], wnspec)
    # Store strong lines to LBL data file:
    indices = np.arange(nlines)[flag]
    for i in indices:
      lbl.write(struct.pack("dddi", wn[i], Elow[i], gf[i], isotopes[iiso[i]]))

  # Close LBL file:
  lbl.close()

  # Convert from cm2 molec-1 to cm-1 amagat-1:
  continuum *= c.N0
  # Save Continuum data to file:
  with open(cont_out, "w") as f:
    # Write header:
    f.write("@SPECIES\n{:s}\n\n".format(mol))
    f.write("@TEMPERATURES\n        ")
    for j in np.arange(ntemp):
        f.write(" {:10.0f}".format(temperature[j]))
    f.write("\n\n")
    # Write the data:
    f.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-1:\n")
    f.write("@DATA\n")
    for i in np.arange(nwave):
      f.write(" {:12.6f} ".format(wnspec[i]))
      for j in np.arange(ntemp):
        f.write(" {:10.4e}".format(continuum[i,j]))
      f.write("\n")

  print("Successfully rewriten {:s} line-transition info into:\n  '{:s}' and"
        "\n  '{:s}'.".format(dbtype, lbl_out, cont_out))

