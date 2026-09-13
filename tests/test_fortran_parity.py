"""
Parity tests against executed FORTRAN.

Every other test in this suite checks a translation against a reading of the
source or against independent physics.  These check it against numbers the
compiled original actually produced.

The values come from running Example Problem 2 through the build that
`tools/fortran_parity.py` produces.  That deck carries a `DUMP A` card, so
its listing contains the 195-element WINGD COMMON block exactly as WTGEOM
filled it.

These tests skip when the build is absent, so the suite still runs on a
machine without gfortran.  To enable them:

    python tools/fortran_parity.py setup
    python tools/fortran_parity.py run 2

Agreement is asserted at 1e-5 relative.  The listing prints six significant
figures, so that is the floor set by the output format rather than by the
arithmetic.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'tools'))

from pydatcom.geometry.wing import calculate_straight_exposed_geometry

LISTING = ROOT / 'build' / 'fortran' / 'ex2.ours'

# Set by the listing's print format, not by the arithmetic.
PRINT_PRECISION = 1.0e-5

pytestmark = pytest.mark.skipif(
    not LISTING.exists(),
    reason="compiled FORTRAN listing absent; run tools/fortran_parity.py "
           "setup && run 2 to enable execution parity tests")


@pytest.fixture(scope='module')
def wingd():
    """The WINGD COMMON block as the compiled original filled it."""
    from extract_common import extract
    return extract(LISTING, 'A')


@pytest.fixture(scope='module')
def ex2_wing_state():
    """Example Problem 2 case 1 wing, exactly as the deck specifies it.

    $WGPLNF CHRDTP=0.64, SSPNE=1.59, SSPN=1.59, CHRDR=2.90,
            SAVSI=55.0, CHSTAT=0.0, TYPE=1.0$
    """
    return {
        'wing_type': 1.0,
        'wing_chrdr': 2.90, 'wing_chrdtp': 0.64,
        'wing_sspn': 1.59, 'wing_sspne': 1.59,
        'wing_savsi': 55.0, 'wing_chstat': 0.0,
    }


@pytest.fixture(scope='module')
def translated(ex2_wing_state):
    return calculate_straight_exposed_geometry(ex2_wing_state)


# WTGEOM index -> key in the translated geometry.  Indices established by
# reading wtgeom.f; the values are whatever the executable printed.
_WTGEOM_MAP = [
    (3, 'area', 'exposed planform area'),
    (7, 'aspect_ratio', 'exposed aspect ratio'),
    (10, 'root_chord', 'exposed root chord'),
    (16, 'mac', 'exposed mean aerodynamic chord'),
    (27, 'taper_ratio', 'exposed taper ratio'),
    (31, 'mac_span_location', 'exposed MAC span station'),
    (62, 'tan_le', 'tangent of leading-edge sweep'),
    (161, 'mac_c4_theoretical', 'theoretical MAC quarter chord'),
    (195, 'mac_le_theoretical', 'theoretical MAC leading edge'),
]


@pytest.mark.parametrize("index,key,label", _WTGEOM_MAP,
                         ids=[f"A{i}_{k}" for i, k, _ in _WTGEOM_MAP])
def test_wtgeom_matches_executed_fortran(wingd, translated, index, key, label):
    """Each WTGEOM quantity matches the value the original computed."""
    assert index in wingd, f"A({index}) absent from the listing"
    expected = wingd[index]
    assert translated[key] == pytest.approx(expected, rel=PRINT_PRECISION), (
        f"{label}: A({index}) = {expected:.6E} from executed FORTRAN, "
        f"translation gives {translated[key]:.6E}")


def test_quarter_chord_sweep_block_matches(wingd, translated):
    """A(43) and A(44) are cos and tan of the exposed quarter-chord sweep.

    This pins the ANGLES block layout that the whole index mapping rests on:
    a 6-element record of [deg, rad, sin, cos, tan, test].
    """
    import numpy as np
    tan_c4 = translated['tan_c4']
    assert wingd[44] == pytest.approx(tan_c4, rel=PRINT_PRECISION)
    assert wingd[43] == pytest.approx(1.0 / np.sqrt(1.0 + tan_c4**2),
                                      rel=PRINT_PRECISION)
    assert wingd[40] == pytest.approx(np.degrees(np.arctan(tan_c4)),
                                      rel=PRINT_PRECISION)
    # sin/cos/tan must be mutually consistent within the record.
    assert wingd[42] == pytest.approx(wingd[43] * wingd[44],
                                      rel=PRINT_PRECISION)


def test_listing_carries_a_full_wingd_block(wingd):
    """Guards the extractor: WINGD is declared A(195)."""
    assert max(wingd) == 195
    assert len(wingd) == 195


def test_exposed_equals_theoretical_when_sspne_equals_sspn(translated):
    """Example Problem 2 has SSPNE == SSPN, collapsing the two planforms.

    A(16) and A(122) therefore agree in the listing, which is why this deck
    is a clean first parity target.
    """
    assert translated['mac'] == pytest.approx(translated['mac_theoretical'])
    assert translated['area'] == pytest.approx(
        translated['area_theoretical'])
