import pytest

import numpy as np

import repack.utils as u


@pytest.mark.parametrize('db, molecule, isotope',
   [('1H2-16O__POKAZATEL__00400-00500.trans.bz2',   'H2O', '116'),
    ('1H-2H-16O__VTT__00250-00500.trans.bz2',       'H2O', '126'),
    ('12C-16O2__HITEMP.pf',                         'CO2', '266'),
    ('12C-16O-18O__Zak.par',                        'CO2', '268'),
    ('12C-1H4__YT10to10__01100-01200.trans.bz2',    'CH4', '21111'),
    ('12C-1H3-2H__MockName__01100-01200.trans.bz2', 'CH4', '21112')])
def test_get_exomol_mol(db, molecule, isotope):
    mol, iso = u.get_exomol_mol(db)
    assert mol == molecule
    assert iso == isotope


def test_read_lbl():
    wn, elow, gf, iiso = u.read_lbl('CO2_hitran_2.2-2.6um_500-700K_lbl.dat')
    nlines = 184451
    assert len(wn)   == nlines
    assert len(elow) == nlines
    assert len(gf)   == nlines
    assert len(iiso) == nlines
    assert wn[ 0] == 3750.00024
    assert wn[-1] == 4499.9982
    assert elow[ 0] == 6971.2258
    assert elow[-1] == 6302.3546
    np.testing.assert_approx_equal(gf[ 0], 6.27118677e-06)
    np.testing.assert_approx_equal(gf[-1], 5.64550190e-09)
    np.testing.assert_equal(np.unique(iiso),
        np.array([626, 627, 628, 636, 638, 828]))

