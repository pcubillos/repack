# Copyright (c) 2017-2020 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

__all__ = [
    'parser',
    'repack',
    'sort',
    ]

import sys
import os
import struct
import subprocess
import configparser
import multiprocessing as mp

import numpy as np
import scipy.interpolate as sip
import scipy.constants as sc

from . import utils as u
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
    chunksize: Integer
        Maximum size of chunks to read.
    ncpu: Integer
        Number of parallel CPUs to use.
    """
    config = configparser.ConfigParser()
    config.read([cfile])
    section = "REPACK"
    # Input line-transition files:
    lblfiles = config.get(section, "lblfiles")
    lblfiles = lblfiles.split()
    # Database type:
    db = config.get(section, "dbtype")
    # Output file:
    outfile = config.get(section, "outfile")
    # Partition-function file:
    pffile = None
    if config.has_option(section, "pffile"):
        pffile = config.get(section, "pffile")
    # Max chunk size:
    chunksize = 5000000
    if config.has_option(section, "chunksize"):
        chunksize = config.getint(section, "chunksize")

    ncpu = 1
    if config.has_option(section, "ncpu"):
        ncpu = config.getint(section, "ncpu")
    ncpu = np.clip(ncpu, 1, mp.cpu_count()-1)

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
        sthresh, pffile, chunksize, ncpu


def worker(input, output):
    """
    Multiprocessing worker that extracts the line-transition info
    and flags the strong/weak lines between the requested indices.
    """
    for args in iter(input.get, 'STOP'):
        wn, gf, Elow, iiso, tmin, tmax, zmin, zmax, imass, \
            iratio, sthresh, idx = args

        # Low temperature line flagging:
        s = (gf*iratio[iiso]/zmin[iiso] *
             np.exp(-c.C2*Elow/tmin) * (1-np.exp(-c.C2*wn/tmin)))
        alphad = wn/(100*sc.c) * np.sqrt(2.0*c.kB*tmin / (imass[iiso]*c.amu))
        # Line-strength sorted in descending line-trength order:
        isort = np.argsort(alphad/s)
        flag = np.ones(len(iiso), bool)
        u.flag(flag, wn, s/alphad, isort, alphad, sthresh)

        # High temperature line flagging:
        s = (gf*iratio[iiso]/zmax[iiso] *
             np.exp(-c.C2*Elow/tmax) * (1-np.exp(-c.C2*wn/tmax)))
        alphad = wn/(100*sc.c) * np.sqrt(2.0*c.kB*tmax / (imass[iiso]*c.amu))
        isort = np.argsort(alphad/s)
        flag2 = np.ones(len(iiso), bool)
        u.flag(flag2, wn, s/alphad, isort, alphad, sthresh)

        output.put((flag, flag2, wn, gf, Elow, iiso, idx))


def repack(cfile):
    """
    Re-pack line-transition data into lbl data for strong lines and
    continuum data for weak lines.

    Parameters
    ----------
    cfile: String
        A repack configuration file.
    """
    banner = 70 * ":"
    # Parse configuration file:
    args = parser(cfile)
    files, dbtype, outfile, tmin, tmax, dtemp, wnmin, wnmax, dwn, \
        sthresh, pffile, chunksize, ncpu = args

    # Auto-detect sorted files:
    files = [f.replace('.trans','.trans.sort')
             if os.path.exists(f.replace('.trans','.trans.sort')) else f
             for f in files]

    missing = [file for file in files if not os.path.exists(file)]
    if len(missing) > 0:
        miss_list = '\n  '.join(missing)
        print(f"\n{banner}\n"
               "  File(s) not Found Error: These files are missing:\n"
              f"  {miss_list}"
              f"\n{banner}\n")
        sys.exit(0)

    # Grid sampling for continuum:
    if dwn != 0 and dtemp != 0:
        ntemp = int((tmax-tmin)/dtemp + 1)
        nwave = int((wnmax-wnmin)/dwn + 1)
    else:
        ntemp, nwave = 0, 0

    temperature = np.linspace(tmin, tmax, ntemp)
    wnspec      = np.linspace(wnmin, wnmax, nwave)
    continuum = np.zeros((nwave, ntemp), np.double)

    if dbtype not in ["hitran", "exomol", "kurucz"]:
        print(f"\n{banner}\n  Error: Invalid database ({dbtype}), "
              f"dbtype must be either hitran, exomol, or kurucz.\n{banner}\n")
        sys.exit(0)

    # Parse input files:
    nfiles = len(files)
    suff, mol, isot, pf, states = [], [], [], [], []
    for dfile in files:
        s, m, iso, p, st = u.parse_file(dfile, dbtype)
        suff.append(s)
        mol.append(m)
        isot.append(iso)
        pf.append(p)
        if st is not None:
            states.append(st)

    # Uncompress states:
    allstates = np.unique(states)
    sdelete, sproc = [], []
    for state in allstates:
        if state.endswith(".bz2"):
            proc = subprocess.Popen(["bzip2", "-dk", state])
            sproc.append(proc)
            sdelete.append(os.path.realpath(state).replace(".bz2", ""))

    if len(np.unique(mol)) > 1:
        print(f"\n{banner}\n"
               "  Error: All input files must correspond to the same molecule."
              f"\n{banner}\n")
        sys.exit(0)
    mol = mol[0]

    z = []  # Interpolator function for partition function per isotope
    if pffile is not None:
        # Read input partition-function file (if given):
        pftemp, partf, isotopes = u.read_pf(pffile, dbtype="pyrat")
        isotopes = list(isotopes)
        for pfvalue in partf:
            z.append(sip.interp1d(pftemp, pfvalue, kind='slinear'))
    else:
        isotopes = list(np.unique(isot))
    niso = len(isotopes)

    # Isotopic abundance ratio and mass:
    iratio, imass = u.read_iso(mol, isotopes, dbtype)
    if np.any(iratio==0)  or np.any(imass==0):
        raise ValueError('One or more isotopes have missing isotopic ratio '
                         'or mass information in isotopes.dat file.')

    s = suff[:]  # Make a copy
    # File indices for each wavenumber set:
    wnset = []
    while np.size(s) != 0:
        suffix = s[np.argmin(s)]
        idx = []
        for i in range(nfiles):
            if suff[i] == suffix:
                idx.append(i)
                s.remove(suffix)
        wnset.append(idx)
    nsets = len(wnset)  # Number of wavenumber sets:

    # Number of sets ahead to unzip:
    zbuffer = np.amin([2,nsets])
    tdelete, tproc = [], []
    for b in range(zbuffer):
        tdelete.append([])
        tproc.append([])
        for idx in wnset[b]:
            if files[idx].endswith(".bz2"):
                print(f"Unzipping: '{files[idx]}'.")
                proc = subprocess.Popen(["bzip2", "-dk", files[idx]])
                tproc[b].append(proc)
                tdelete[b].append(files[idx].replace(".bz2", ""))

    for proc in sproc:
        proc.communicate()

    iso = np.zeros(nfiles, int)
    if dbtype == "exomol":
        for i in range(nfiles):
            iso[i] = isotopes.index(isot[i])
        # Read partition-function and states files:
        lblargs = []
        for j in range(niso):
            try:
                i = isot.index(isotopes[j])
            except:
                lblargs.append([None, None, None])  # placeholder
                continue
            # Partition function:
            if pffile is None:
                temp, part = u.read_pf(pf[i], dbtype)
                z.append(sip.interp1d(temp, part, kind='slinear'))
            # States:
            elow, degen = u.read_states(states[i])
            lblargs.append([elow, degen, j])

    else:  #dbtype in ["hitran", "kurucz"]:
        lblargs = [[None, None, None]]  # Trust me

    # Turn isotopes from string to integer data type:
    isotopes = np.asarray(isotopes, int)

    # Set output file names:
    lbl_out  = f"{outfile}_lbl.dat"
    cont_out = f"{outfile}_continuum.dat"
    # Output line-by-line file:
    lblf = open(lbl_out, "wb")

    # Create queues and start worker processes:
    task_queue = mp.Queue()
    done_queue = mp.Queue()
    for i in range(ncpu):
        mp.Process(target=worker, args=(task_queue, done_queue)).start()

    total_lines = 0
    # Read line-by-line files:
    for i in range(nsets):
        # Make sure current files are uncompressed:
        for p in tproc[i]:
            p.communicate()
        # Uncompress following set:
        if zbuffer < nsets:
            tdelete.append([])
            tproc.append([])
            for idx in wnset[zbuffer]:
                if files[idx].endswith(".bz2"):
                    print(f"Unzipping: '{files[idx]}'.")
                    proc = subprocess.Popen(["bzip2", "-dk", files[idx]])
                    tproc[zbuffer].append(proc)
                    tdelete[zbuffer].append(files[idx].replace(".bz2", ""))
            zbuffer += 1

        # Gather data from all files in this wavenumber range:
        lbl, istart, nlines = [], [], []
        for k in range(len(wnset[i])):
            # Initialize lbl object (not reading yet):
            idx = wnset[i][k]
            j = int(iso[wnset[i][k]])
            print(f"Reading: '{files[idx]}'.")
            lbl.append(u.lbl(files[idx], dbtype, *lblargs[j]))
            # Find initial value in range:
            i0 = lbl[k].bs(wnmin, 0,  lbl[k].nlines-1)
            while i0 > 0 and lbl[k].getwn(i0-1) >= wnmin:
                i0 -= 1
            # Find final value in range:
            iN = lbl[k].bs(wnmax, i0, lbl[k].nlines-1)
            while iN < lbl[k].nlines-1 and lbl[k].getwn(iN+1) <= wnmin:
                iN += 1
            istart.append(i0)
            nlines.append(iN-i0+1)

        # Count lines, set target chunk size:
        nchunks = int(np.sum(nlines)/chunksize) + 1
        target  = np.sum(nlines)/nchunks
        chunk   = np.zeros((len(wnset[i]), nchunks+1), int)
        # First and last are easy:
        chunk[:,      0] = istart
        chunk[:,nchunks] = chunk[:,0] + nlines
        # Easy-case split if only one file:
        if len(wnset[i]) == 1:
            chunk[0] = np.linspace(chunk[0,0], chunk[0,-1], nchunks+1)
        else:  # Complicated case:
            wnchunk = np.linspace(wnmin, wnmax, nchunks+1)
            # Intermediate boundaries:
            for n in range(1, nchunks):
                zero = np.sum(chunk[:,n-1])
                wnchunk[n] = u.wnbalance(lbl, wnchunk[n-1], wnmax, target, zero)
                for k in range(len(wnset[i])):
                    chunk[k,n] = lbl[k].bs(
                        wnchunk[n], chunk[k,n-1], chunk[k,nchunks])

        # Proccess chunks:
        for n in range(nchunks):
            gf   = np.array([])
            Elow = np.array([])
            wn   = np.array([])
            iiso = np.array([], int)
            for k in range(len(wnset[i])):
                # Read the LBL files by chunks:
                gfosc, el, wnumber, isoID = lbl[k].read(chunk[k,n:n+2])
                gf   = np.hstack([gf,   gfosc  ])
                Elow = np.hstack([Elow, el     ])
                wn   = np.hstack([wn,   wnumber])
                iiso = np.hstack([iiso, isoID  ])
            # Sort by wavelength:
            asort = np.argsort(wn)
            gf   = gf  [asort]
            Elow = Elow[asort]
            wn   = wn  [asort]
            iiso = iiso[asort]

            zmin = np.array([z[j](tmin) for j in range(niso)])
            zmax = np.array([z[j](tmax) for j in range(niso)])

            args = (wn, gf, Elow, iiso, tmin, tmax, zmin, zmax,
                    imass, iratio, sthresh, n)
            task_queue.put(args)

        collect_wn = []
        collect_gf = []
        collect_Elow = []
        collect_iiso = []
        chunk_idx = []
        for n in range(nchunks):
            flag, flag2, wn, gf, Elow, iiso, idx = done_queue.get()

            # Low temperature line flagging:
            comment = f" (chunk {n+1}/{nchunks})" if nchunks > 1 else ""
            ntotal = len(flag)
            nkept = np.sum(flag)
            print(f"  Flagging lines at {tmin:4.0f} K{comment}:")
            print(f"  Compression rate:       {100.0*(1-nkept/ntotal):.2f}%,"
                  f"  {nkept:9,d}/{ntotal:10,d} lines.")
            # High temperature line flagging:
            nkept = np.sum(flag2)
            print(f"  Flagging lines at {tmax:4.0f} K:")
            print(f"  Compression rate:       {100.0*(1-nkept/ntotal):.2f}%,"
                  f"  {nkept:9,d}/{ntotal:10,d} lines.")
            # Update flag:
            flag |= flag2
            nkept = np.sum(flag)
            print(f"  Total compression rate: {100.0*(1-nkept/ntotal):.2f}%,"
                  f"  {nkept:9,d}/{ntotal:10,d} lines.\n")

            # Store weak lines to continuum file as function of temp:
            for t,T in enumerate(temperature):
                Z = np.array([z[j](T) for j in range(niso)])
                # Line strength in cm molec-1
                s = (c.C3 * gf*iratio[iiso]/Z[iiso] *
                     np.exp(-c.C2*Elow/T) * (1-np.exp(-c.C2*wn/T)))
                # Continuum opacity in cm-1 amagat-1 (instead of cm2 molec-1):
                u.continuum(s[~flag], wn[~flag], continuum[:,t], wnspec)
            # Store strong lines to LBL data file:
            collect_wn.append(wn[flag])
            collect_gf.append(gf[flag])
            collect_Elow.append(Elow[flag])
            collect_iiso.append(iiso[flag])
            chunk_idx.append(idx)

        for n in np.argsort(chunk_idx):
            for wn, gf, Elow, iiso in zip(collect_wn[n], collect_gf[n],
                    collect_Elow[n], collect_iiso[n]):
                lblf.write(struct.pack("dddi", wn, Elow, gf, isotopes[iiso]))

        for k in range(len(wnset[i])):
            lbl[k].close()
        # Delete unzipped sets:
        for f in tdelete[i]:
            os.remove(f)
        total_lines += np.sum(nlines)

    # Close LBL file:
    print(f"With a threshold strength factor of {sthresh},\n"
          f"kept a total of {lblf.tell()/struct.calcsize('dddi'):,.0f} "
          f"line transitions out of {total_lines:,.0f} lines.\n")
    lblf.close()

    for i in range(ncpu):
        task_queue.put('STOP')

    if ntemp != 0:
        # Convert from cm2 molec-1 to cm-1 amagat-1:
        continuum *= c.N0
        # Save Continuum data to file:
        with open(cont_out, "w") as f:
            # Write header:
            f.write(f"@SPECIES\n{mol}\n\n")
            f.write("@TEMPERATURES\n        ")
            for temp in temperature:
                f.write(f" {temp:10.0f}")
            f.write("\n\n")
            # Write the data:
            f.write("# Wavenumber in cm-1, opacity in cm-1 amagat-1:\n")
            f.write("@DATA\n")
            for i in range(nwave):
                f.write(f" {wnspec[i]:12.6f} ")
                for j in range(ntemp):
                    f.write(f" {continuum[i,j]:10.4e}")
                f.write("\n")
    cont_msg = f" and\n  '{cont_out}'" if ntemp != 0 else ""

    print(f"Successfully rewriten {dbtype} line-transition info into:\n"
          f"  '{lbl_out}'{cont_msg}.")

    # Delete unzipped set:
    for f in sdelete:
        os.remove(f)


def sort_worker(input, output):
    """
    Multiprocessing worker that reads the line transition wavenumbers
    between the requested indices.
    """
    for args in iter(input.get, 'STOP'):
        lblfile, rec_size, elow, istart, iend, idx = args
        with u.fopen(lblfile, "r") as file:
            nlines = iend - istart
            file.seek(istart*rec_size)
            iup = np.zeros(nlines, int)
            ilo = np.zeros(nlines, int)
            for i in range(nlines):
                line = file.readline()
                iup[i] = line[ 0:12]
                ilo[i] = line[13:25]
        wn = elow[iup-1] - elow[ilo-1]
        output.put((wn, idx))


def sort(cfile):
    """
    Sort the ExoMol .trans files by wavenumber for MARVELized .states files

    Parameters
    ----------
    cfile: String
        A repack configuration file.
    """
    banner = 70 * ":"
    args = parser(cfile)
    files, dbtype, outfile, tmin, tmax, dtemp, wnmin, wnmax, dwn, \
        sthresh, pffile, chunksize, ncpu = args

    if dbtype != "exomol":
        sys.exit(0)

    missing = [file for file in files if not os.path.exists(file)]
    if len(missing) > 0:
        miss_list = '\n  '.join(missing)
        print(f"\n{banner}\n"
               "  File(s) not Found Error: These files are missing:\n"
              f"  {miss_list}"
              f"\n{banner}\n")
        sys.exit(0)

    # Parse input files:
    nfiles = len(files)
    suff, mol, isot, pf, states = [], [],  [], [], []
    for dfile in files:
        s, m, iso, p, st = u.parse_file(dfile, dbtype)
        suff.append(s)
        mol.append(m)
        isot.append(iso)
        pf.append(p)
        if st is not None:
            states.append(st)

    # Uncompress states:
    allstates = np.unique(states)
    sdelete, sproc = [], []
    for state in allstates:
        if state.endswith(".bz2"):
            proc = subprocess.Popen(["bzip2", "-dk", state])
            sproc.append(proc)
            sdelete.append(os.path.realpath(state).replace(".bz2", ""))

    isotopes = list(np.unique(isot))
    niso = len(isotopes)

    zbuffer = np.amin([2, nfiles])
    tdelete, tproc = [], []
    for idx in range(zbuffer):
        if files[idx].endswith(".bz2"):
            print(f"Unzipping: '{files[idx]}'.")
            proc = subprocess.Popen(["bzip2", "-dk", files[idx]])
            tproc.append(proc)
            tdelete.append(files[idx].replace(".bz2", ""))

    for proc in sproc:
        proc.communicate()

    iso = np.zeros(nfiles, int)
    for i in range(nfiles):
        iso[i] = isotopes.index(isot[i])

    lblargs = []
    for j in range(niso):
        i = isot.index(isotopes[j])
        elow, degen = u.read_states(states[i])
        lblargs.append([elow, degen, j])

    # Turn isotopes from string to integer data type:
    isotopes = np.asarray(isotopes, int)

    # Create queues and start worker processes:
    task_queue = mp.Queue()
    done_queue = mp.Queue()
    for i in range(ncpu):
        mp.Process(target=sort_worker, args=(task_queue, done_queue)).start()

    zproc = []
    for i in range(nfiles):
        # Make sure current files are uncompressed:
        tproc[i].communicate()
        # Uncompress following set:
        if zbuffer < nfiles and files[zbuffer].endswith(".bz2"):
            print(f"Unzipping: '{files[zbuffer]}'.")
            proc = subprocess.Popen(["bzip2", "-dk", files[zbuffer]])
            tproc.append(proc)
            tdelete.append(files[zbuffer].replace(".bz2", ""))
            zbuffer += 1

        # Initialize lbl object (not reading yet):
        j = iso[i]
        print(f"Reading: '{files[i]}'.")
        lbl = u.lbl(files[i], dbtype, *lblargs[j])

        nlines = lbl.nlines
        chunksize = int(nlines/ncpu) + 1
        chunks = np.linspace(0, nlines, ncpu+1, dtype=int)

        for k in range(ncpu):
            args = lbl.lblfile, lbl.llen, lbl.elow, chunks[k], chunks[k+1], k
            task_queue.put(args)

        wn = [None] * ncpu
        for k in range(ncpu):
            w, idx = done_queue.get()
            wn[idx] = w
        all_wn = np.concatenate(wn)
        wn_sort = np.argsort(np.argsort(all_wn))

        lines = np.zeros(nlines, f'U{lbl.llen}')
        lbl.file.seek(0)
        for k in range(nlines):
            lines[wn_sort[k]] = lbl.file.readline()
        sort_file = lbl.lblfile.replace('trans.bz2', 'trans.sort')
        with open(sort_file, 'w') as f:
            f.writelines(lines)

        proc = subprocess.Popen(["bzip2", "-z", sort_file])
        zproc.append(proc)

        lbl.close()
        os.remove(tdelete[i])

    for k in range(ncpu):
        task_queue.put('STOP')

    # Delete unzipped set:
    for state in sdelete:
        os.remove(state)

    for proc in zproc:
        proc.communicate()

