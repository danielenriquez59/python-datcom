"""
Regression tests for ANGDET, DELY and ARCLSS.

Source of truth: datcom-legacy/datcom_2000/angdet.f, dely.f, arclss.f.

AREA1 was already translated in utils/math_utils.py and is covered by
tests/test_area_translation.py.
"""

import numpy as np
import pytest

from pydatcom.utils.table_lookup import angdet, fig68
from pydatcom.geometry.section_params import (
    calculate_dely, calculate_arclss,
    _FIG_41340_24B_TAPER, _FIG_41340_24B_VALUE, _DELY_STATIONS,
)


# --------------------------------------------------------------------------
# ANGDET
# --------------------------------------------------------------------------

@pytest.mark.parametrize("mach", [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0])
def test_angdet_agrees_with_fig68(mach):
    """Two independently translated routines must give the same limit.

    ANGDET solves NACA TR 1135 equations 168 and 138 directly; FIG68
    returns the maximum attached wedge angle from its own tabulated cubic
    when handed a deflection beyond the limit.  They come from different
    source files, so agreement is a genuine cross-check.
    """
    limit, error = fig68(mach, 60.0)
    assert error == 2, "a 60 degree wedge must be past detachment"
    assert np.degrees(angdet(mach)) == pytest.approx(limit, abs=1e-6)


def test_angdet_grows_with_mach():
    values = [angdet(m) for m in (1.2, 2.0, 3.0, 5.0)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_angdet_stays_below_the_asymptotic_limit():
    """The maximum turn angle approaches about 45.6 degrees as Mach grows."""
    for mach in (1.5, 3.0, 10.0, 50.0):
        assert 0.0 < np.degrees(angdet(mach)) < 45.6


def test_angdet_rejects_subsonic_mach():
    for mach in (0.5, 1.0):
        with pytest.raises(ValueError, match="supersonic"):
            angdet(mach)


# --------------------------------------------------------------------------
# DELY
# --------------------------------------------------------------------------

def _section(thickness=0.06, points=201):
    x = np.linspace(0.0, 1.0, points)
    return x, thickness * (1.0 - (2.0 * x - 1.0)**2)


def test_dely_samples_the_source_stations():
    assert _DELY_STATIONS == (0.0015, 0.0600)


def test_dely_is_the_scaled_ordinate_difference():
    x, t = _section()
    result = calculate_dely(x, t)
    assert result['deltay'] == pytest.approx(
        (result['y_at_06'] - result['y_at_0015']) * 100.0)


def test_dely_grows_with_section_thickness():
    values = [calculate_dely(*_section(t))['deltay']
              for t in (0.03, 0.06, 0.12)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_dely_is_positive_for_a_normal_section():
    """The 6 percent station sits above the 0.15 percent one."""
    assert calculate_dely(*_section())['deltay'] > 0.0


def test_dely_rejects_bad_input():
    x, t = _section()
    with pytest.raises(ValueError):
        calculate_dely(x, t[:10])
    with pytest.raises(ValueError):
        calculate_dely([0.5], [0.01])


# --------------------------------------------------------------------------
# ARCLSS
# --------------------------------------------------------------------------

@pytest.mark.parametrize("taper,expected",
                         list(zip(_FIG_41340_24B_TAPER,
                                  _FIG_41340_24B_VALUE)))
def test_figure_41340_24b_source_coordinates(taper, expected):
    """Every entry of the source DATA is reproduced."""
    assert calculate_arclss(float(taper), 1.0, 1.0)['factor'] == \
        pytest.approx(expected, abs=1e-12)


def test_figure_41340_24b_peaks_in_the_middle():
    """The curve rises to 0.496 at taper 0.3 and falls to zero at both ends."""
    assert _FIG_41340_24B_VALUE[0] == 0.0
    assert _FIG_41340_24B_VALUE[-1] == 0.0
    assert max(_FIG_41340_24B_VALUE) == 0.496


def test_arclss_forms_the_two_derived_quantities():
    cosine = np.cos(np.deg2rad(30.0))
    result = calculate_arclss(0.3, cosine, reference=4.0)
    assert result['classified'] == pytest.approx((0.496 + 1.0) * cosine)
    assert result['ratio'] == pytest.approx(4.0 / result['classified'])


def test_arclss_rejects_a_vanishing_denominator():
    with pytest.raises(ValueError, match="vanishing"):
        calculate_arclss(0.3, 0.0, 1.0)
