# Copyright (c) 2017-2020 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

__all__ = [
    "kB",
    "amu",
    "e",
    "me",
    "N0",
    "nano",
    "C1",
    "C2",
    "C3",
    "ROOT",
    ]

import os

import scipy.constants as sc


# Boltzmann constant in erg K-1:
kB  = sc.k * 1e7
# Unified atomic mass in g:
amu = sc.physical_constants["unified atomic mass unit"][0] * 1e3
# Elementary charge in statcoulomb (cm3/2 g1/2 s-1):
e  = 4.803205e-10
# Electron mass in g:
me = sc.m_e * 1e3
# Amagat in molecules cm-3:
N0 = sc.physical_constants[
        "Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1e-6

# One nanometer in centimeters:
nano = 1e-7

# Other constructed constants:
C1  = 4.0 * sc.epsilon_0 * sc.m_e * sc.c**2 / sc.e**2 * 0.01  # cm-1
C2  = sc.h * (sc.c * 100.0) / sc.k                            # cm K-1
C3  = sc.pi * e**2 / (me * (100*sc.c)**2)                     # cm

ROOT = os.path.realpath(os.path.dirname(__file__) + "/../..") + "/"
