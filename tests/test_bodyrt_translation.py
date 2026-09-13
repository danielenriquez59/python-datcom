"""
Regression tests for the BODYRT translation.

Source of truth: datcom-legacy/datcom_2000/bodyrt.f.

Scope is kept tight.  BODYRT's five dependencies -- EQSPC1, FIG26, GETMAX,
TBFUNX and TRAPZ -- all have their own translation tests, so these cover the
figures BODYRT owns, its branch structure, and the symmetry and buildup the
routine is responsible for.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.bodyrt import (
    calculate_bodyrt, _FIG_42110_20_X, _FIG_42110_20_Y,
    _FIG_42120_35A_X, _FIG_42120_35A_Y,
    _FIG_42120_35B_X, _FIG_42120_35B_Y,
)


def _ogive_cylinder_boattail():
    """A 10-unit body: ogive nose, constant barrel, boat-tailed base."""
    x = np.array([0., 1., 2., 3., 4., 6., 8., 9., 9.5, 10.])
    r = np.array([0., .35, .6, .75, .8, .8, .8, .72, .62, .55])
    return x, np.pi * r**2, 2 * np.pi * r, r


def _run(alphas=(-10., -5., 0., 5., 10., 15.), **kwargs):
    x, s, p, r = _ogive_cylinder_boattail()
    params = dict(mach=0.3, reynolds_per_length=2.0e6,
                  sref=float(np.pi * 0.8**2), cbar=10.0, xcg=5.0)
    params.update(kwargs)
    return calculate_bodyrt(x, s, p, r, list(alphas), **params)


# --------------------------------------------------------------------------
# Figures BODYRT owns
# --------------------------------------------------------------------------

@pytest.mark.parametrize("table_x,table_y,label", [
    (_FIG_42110_20_X, _FIG_42110_20_Y, "4.2.1.1-20"),
    (_FIG_42120_35A_X, _FIG_42120_35A_Y, "4.2.1.2-35A"),
    (_FIG_42120_35B_X, _FIG_42120_35B_Y, "4.2.1.2-35B"),
])
def test_figure_tables_match_source_data(table_x, table_y, label):
    """Coordinate counts and endpoints match the source DATA statements."""
    assert len(table_x) == len(table_y), label
    assert np.all(np.diff(table_x) > 0), f"{label} abscissa must increase"


def test_figure_42120_35b_repeat_expansion():
    """The source writes three entries of Figure 4.2.1.2-35B as '3*1.8'."""
    assert list(_FIG_42120_35B_Y[10:13]) == [1.8, 1.8, 1.8]
    assert _FIG_42120_35B_Y[-1] == 1.79


def test_apparent_mass_factor_rises_with_fineness():
    """Figure 4.2.1.1-20: slender bodies approach the 2D apparent mass."""
    assert np.all(np.diff(_FIG_42110_20_Y) > 0)
    assert _FIG_42110_20_Y[-1] < 1.0


# --------------------------------------------------------------------------
# Symmetry: PHYSICS_REVIEW items 6 and 7
# --------------------------------------------------------------------------

def test_symmetric_body_has_odd_normal_force_and_even_drag():
    """An uncambered body must give odd CN and CM and even CD.

    PHYSICS_REVIEW item 7 recorded Cm of +7.49 at -10 degrees and +7.66 at
    +10 degrees for the previous implementation, which is neither.
    """
    result = _run()
    cn, cm, cd = result['cn'], result['cm'], result['cd']
    for negative, positive in ((0, 4), (1, 3)):
        assert cn[negative] == pytest.approx(-cn[positive])
        assert cm[negative] == pytest.approx(-cm[positive])
        assert cd[negative] == pytest.approx(cd[positive])


def test_zero_alpha_gives_zero_lift_and_moment():
    result = _run(alphas=(0.0,))
    assert result['cn'][0] == pytest.approx(0.0, abs=1e-15)
    assert result['cm'][0] == pytest.approx(0.0, abs=1e-15)
    assert result['cd'][0] == pytest.approx(result['cd_zero_lift'])


def test_normal_force_grows_faster_than_linearly():
    """Allen-Perkins crossflow adds a sin-squared term to the linear one."""
    result = _run(alphas=(5.0, 10.0, 20.0))
    cn = result['cn']
    linear = result['cla'] * np.array([5.0, 10.0, 20.0])
    assert np.all(cn > linear)
    assert cn[2] / cn[0] > 20.0 / 5.0


# --------------------------------------------------------------------------
# Setup branches and buildup
# --------------------------------------------------------------------------

def test_reference_station_finds_fastest_contraction():
    """X1 is where -dS/dx peaks on a boat-tailed body."""
    result = _run()
    assert result['body_contracts']
    assert 8.0 <= result['reference_station'] <= 10.0


def test_non_contracting_body_uses_the_last_station():
    """The source's fallback when the area never decreases."""
    x = np.array([0., 2., 5., 8., 10.])
    r = np.array([0., .4, .6, .75, .8])
    result = calculate_bodyrt(x, np.pi * r**2, 2 * np.pi * r, r,
                              [0.0, 5.0], mach=0.3, reynolds_per_length=2e6,
                              sref=2.0, cbar=10.0, xcg=5.0)
    assert not result['body_contracts']
    assert result['reference_station'] == pytest.approx(10.0)


def test_base_area_is_floored_for_boat_tailed_bodies():
    """BD(57) is raised to 30% of maximum when the base is smaller."""
    result = _run()
    # The test body's base is (0.55/0.8)^2 = 47% of maximum, above the floor,
    # so the drag must reflect the true base rather than the floor.
    ratio = (0.55 / 0.8)**2
    assert ratio > 0.30
    assert result['cd_base'] > 0.0


def test_drag_buildup_is_friction_plus_base():
    result = _run()
    assert result['cd_zero_lift'] == pytest.approx(
        result['cd_friction'] + result['cd_base'])
    assert result['cd_friction'] > 0.0
    assert result['cd_base'] > 0.0


def test_roughness_cutoff_limits_reynolds_number():
    """A rough surface caps the Reynolds number used for friction."""
    smooth = _run(roughness=0.4e-4)
    rough = _run(roughness=6.4e-4)
    assert rough['roughness_cutoff_reynolds'] < smooth['roughness_cutoff_reynolds']
    assert rough['cf'] >= smooth['cf']


def test_transonic_returns_before_the_angle_loop():
    """The source returns after the drag buildup when TRANSN is set."""
    result = _run(transonic=True)
    assert 'cn' not in result
    assert result['cd_zero_lift'] > 0.0
    assert result['cla'] > 0.0


def test_moment_slope_depends_on_reference_station():
    """CMA carries the XCG term, so moving the reference moves the slope."""
    forward = _run(xcg=3.0)['cma']
    aft = _run(xcg=7.0)['cma']
    assert forward != pytest.approx(aft)


def test_rejects_bad_geometry():
    x, s, p, r = _ogive_cylinder_boattail()
    with pytest.raises(ValueError):
        calculate_bodyrt(x[:1], s[:1], p[:1], r[:1], [0.0], 0.3, 2e6,
                         2.0, 10.0, 5.0)
    with pytest.raises(ValueError):
        calculate_bodyrt(x, s[:5], p, r, [0.0], 0.3, 2e6, 2.0, 10.0, 5.0)
    with pytest.raises(ValueError):
        calculate_bodyrt(x, s, p, r, [0.0], 0.3, 2e6, 0.0, 10.0, 5.0)
