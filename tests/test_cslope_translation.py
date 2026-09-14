"""
Regression tests for the CSLOPE translation.

Source of truth: datcom-legacy/datcom_2000/cslope.f.

The slope itself comes from TBFUNX, which has its own translation tests, so
these cover the surface selection, the sign convention, the per-station
supplied-value skip, and the no-tail early return.
"""

import numpy as np
import pytest

from pydatcom.geometry.cslope import calculate_cslope, _STATIONS


def _section(thickness=0.06):
    """A symmetric parabolic section: yu = t*(1-(2x-1)^2), yl = -yu."""
    x = np.linspace(0.0, 1.0, 41)
    yu = thickness * (1.0 - (2.0 * x - 1.0)**2)
    return x, yu, x, -yu


def _run(**kwargs):
    x_upper, y_upper, x_lower, y_lower = _section()
    params = dict(x_upper=x_upper, y_upper=y_upper,
                  x_lower=x_lower, y_lower=y_lower,
                  z_wing=0.0, z_tail=2.0)
    params.update(kwargs)
    return calculate_cslope(**params)


def test_six_chord_stations():
    assert _STATIONS == (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    assert _run()['slope_angles_deg'].shape == (6,)


def test_tail_above_reads_the_upper_surface():
    assert _run(z_tail=2.0)['surface'] == 'upper'


def test_tail_below_reads_the_lower_surface():
    assert _run(z_tail=-2.0)['surface'] == 'lower'


def test_level_tail_reads_the_lower_surface():
    """The source test is strict: ZH > ZW, so a level tail falls through."""
    assert _run(z_wing=1.0, z_tail=1.0)['surface'] == 'lower'


def test_no_horizontal_tail_leaves_every_station_unset():
    result = _run(has_horizontal_tail=False)
    assert result['surface'] is None
    assert np.all(np.isnan(result['slope_angles_deg']))
    assert np.all(np.isnan(result['slopes']))


def test_sign_convention_makes_a_symmetric_section_agree():
    """The lower surface slope is negated, so both give the same angles.

    For yl = -yu the raw derivatives are opposite; the source's sign flip on
    the lower branch is what makes the two readings agree.
    """
    upper = _run(z_tail=2.0)['slope_angles_deg']
    lower = _run(z_tail=-2.0)['slope_angles_deg']
    assert np.allclose(upper, lower)


def test_slope_angles_are_antisymmetric_about_midchord():
    """A symmetric parabolic section rises then falls by equal amounts."""
    angles = _run()['slope_angles_deg']
    assert angles[0] == pytest.approx(-angles[5])
    assert angles[1] == pytest.approx(-angles[4])
    assert angles[2] == pytest.approx(-angles[3])


def test_angles_are_the_arctangent_of_the_slopes():
    from pydatcom.utils.constants import RAD
    result = _run()
    assert np.allclose(result['slope_angles_deg'],
                       np.arctan(result['slopes']) * RAD)


def test_supplied_stations_are_left_untouched():
    """A station the user has already given is skipped, as the UNUSED test does."""
    supplied = [None, 1.0, None, None, 2.0, None]
    result = _run(supplied=supplied)
    assert np.all(np.isnan(result['slope_angles_deg'][[1, 4]]))
    assert np.all(np.isfinite(result['slope_angles_deg'][[0, 2, 3, 5]]))


def test_thicker_section_has_steeper_slopes():
    thin = calculate_cslope(*_section(0.03), z_wing=0.0, z_tail=2.0)
    thick = calculate_cslope(*_section(0.12), z_wing=0.0, z_tail=2.0)
    assert abs(thick['slope_angles_deg'][0]) > abs(thin['slope_angles_deg'][0])


def test_rejects_bad_surface_arrays():
    x_upper, y_upper, x_lower, y_lower = _section()
    with pytest.raises(ValueError, match="upper"):
        calculate_cslope(x_upper[:3], y_upper, x_lower, y_lower, 0.0, 2.0)
    with pytest.raises(ValueError, match="lower"):
        calculate_cslope(x_upper, y_upper, x_lower[:3], y_lower, 0.0, -2.0)
    with pytest.raises(ValueError, match="six"):
        _run(supplied=[None, None])
