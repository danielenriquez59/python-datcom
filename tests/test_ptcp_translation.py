"""PTCP and ARCCOS translation tests.

PTCP's quadrature weights admit an exact check that does not depend on the
source at all: each ``AREA`` element is an integration rule over one station
interval, so its weights must sum to that interval's width, and each
``AMON`` element is a *first moment* rule, so its weights must sum to the
interval's integral of ``r dr``.  Both hold to the precision the source
prints its constants at, which verifies every one of the 47 weight
expressions independently of how they were transcribed.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.ptcp import (
    calculate_ptcp, TIP, ROOT, _STATIONS,
    _area_two_to_four, _area_five_to_ten_staggered,
    _area_five_to_ten_regular, _area_eleven_to_nineteen,
    _moment_ten_point_a, _moment_ten_point_b, _moment_eleven_to_nineteen,
)
from pydatcom.utils.constants import PI
from pydatcom.utils.legacy_numeric import arccos


# --------------------------------------------------------------------------
# ARCCOS
# --------------------------------------------------------------------------

@pytest.mark.parametrize('value', [-0.999, -0.75, -0.5, -0.25, 0.0, 0.25,
                                   0.5, 0.75, 1.0])
def test_arccos_matches_the_real_inverse_cosine_inside_the_unit_range(value):
    # To the source's ten-digit PI, which zero and the negative half use.
    assert arccos(value) == pytest.approx(float(np.arccos(value)), abs=1e-9)


def test_arccos_uses_the_source_pi():
    assert arccos(0.0) == PI / 2.0
    assert arccos(-0.5) == pytest.approx(PI - np.pi / 3.0, abs=1e-15)
    assert arccos(0.5) == pytest.approx(np.pi / 3.0, abs=1e-15)


def test_arccos_returns_zero_at_minus_one_where_pi_is_correct():
    """Source defect, confined to exactly a = -1.

    ``sqrt(1-a^2)/a`` is negative zero there, and the source adds ``pi``
    only when ``X .LT. 0.0``, which negative zero does not satisfy.  Values
    just inside pick the correction up normally, so the error is a jump at
    one point rather than a region.
    """
    assert arccos(-1.0) == pytest.approx(0.0, abs=1e-15)
    assert np.arccos(-1.0) == pytest.approx(np.pi)
    # Immediately adjacent arguments are continuous with pi, not with zero.
    assert arccos(-0.9999999) == pytest.approx(np.pi, abs=1e-3)


@pytest.mark.parametrize('value', [1.5, 2.0, 5.0, 100.0])
def test_arccos_continues_as_arccosh_above_one(value):
    """The true arccos is imaginary there; the source returns its magnitude."""
    assert arccos(value) == pytest.approx(float(np.arccosh(value)))


def test_arccos_below_minus_one_follows_the_source_expression():
    assert arccos(-2.0) == pytest.approx(
        float(np.log(abs(-2.0 + np.sqrt(3.0)))))
    assert arccos(-2.0) < 0.0


def test_arccos_rejects_non_finite_arguments():
    for value in (np.nan, np.inf):
        with pytest.raises(ValueError, match="finite"):
            arccos(value)


# --------------------------------------------------------------------------
# The quadrature weights, checked against what they must integrate
# --------------------------------------------------------------------------

def _weights(builder, count, size):
    """Recover each element's weight sum by feeding in a constant one."""
    return np.array(builder(np.ones(size)))[:count]


def test_area_weights_integrate_each_interval_width():
    """A constant distribution must give each panel its own station width."""
    # AREA(2) to AREA(4) span 0.1 each.
    np.testing.assert_allclose(_weights(_area_two_to_four, 3, 19),
                               0.1, atol=2e-6)
    for builder in (_area_five_to_ten_staggered, _area_five_to_ten_regular):
        np.testing.assert_allclose(_weights(builder, 6, 19), 0.1, atol=2e-6)
    # AREA(11) to AREA(19) span 1.0 each.
    np.testing.assert_allclose(_weights(_area_eleven_to_nineteen, 9, 19),
                               1.0, atol=2e-6)


def _first_moments():
    edges = np.concatenate(([0.0], _STATIONS))
    return (edges[1:]**2 - edges[:-1]**2) / 2.0


def test_moment_weights_integrate_each_interval_first_moment():
    """A constant distribution must give each panel its integral of r dr.

    This holds exactly for every element of both ten-point sets except the
    last of set A, which is pinned separately below, and for all nine of
    the unit-step elements.
    """
    expected = _first_moments()
    for builder in (_moment_ten_point_a, _moment_ten_point_b):
        np.testing.assert_allclose(_weights(builder, 9, 19),
                                   expected[:9], atol=2e-6)
    np.testing.assert_allclose(_weights(_moment_ten_point_b, 10, 19),
                               expected[:10], atol=2e-6)
    np.testing.assert_allclose(_weights(_moment_eleven_to_nineteen, 9, 19),
                               expected[10:19], atol=2e-6)


def test_last_moment_weight_of_set_a_misses_its_first_moment():
    """AMON(10) on the ten-point tip path sums to 0.095589, not 0.095.

    Recorded rather than corrected, and deliberately not called a defect.
    The coefficient 0.017389 on the last ordinate is shared with AREA(10)
    of the same block, so the pair was evidently derived together for the
    endpoint rather than as plain Newton-Cotes.  What is notable is that
    AREA(10) does hit its own invariant exactly, as does AMON(10) of set B
    over the same interval.
    """
    got = _weights(_moment_ten_point_a, 10, 19)[-1]
    assert got == pytest.approx(-.030795 + .108995 + .017389)
    assert got != pytest.approx(_first_moments()[9], abs=1e-5)
    # The area rule over the same interval is exact.
    assert _weights(_area_five_to_ten_staggered, 6, 19)[-1] == \
        pytest.approx(0.1, abs=2e-6)
    # So is the other set's moment rule over the same interval.
    assert _weights(_moment_ten_point_b, 10, 19)[-1] == \
        pytest.approx(_first_moments()[9], abs=2e-6)


def test_first_area_panel_is_a_special_endpoint_rule():
    """AREA(1) does not integrate to 0.1; it handles the leading edge.

    The tip and root forms differ from each other and both from the plain
    panel width, which is why the invariant above starts at AREA(2).
    """
    probe = np.ones(19)
    from pydatcom.aerodynamics.ptcp import _area_first
    tip = _area_first(probe, root=False)
    root = _area_first(probe, root=True)
    assert tip == pytest.approx(.113301 - .032975)
    assert root == pytest.approx(.130786 - .045340)
    assert tip != pytest.approx(0.1, abs=1e-3)
    assert root != pytest.approx(0.1, abs=1e-3)
    assert root > tip


def test_the_two_ten_point_moment_sets_are_genuinely_different():
    """They integrate the same thing but are not the same expressions."""
    probe = np.linspace(0.3, 0.9, 19)
    first = np.array(_moment_ten_point_a(probe))
    second = np.array(_moment_ten_point_b(probe))
    assert not np.allclose(first, second)


def test_the_two_area_sets_for_five_to_ten_are_genuinely_different():
    probe = np.linspace(0.3, 0.9, 19)
    assert not np.allclose(np.array(_area_five_to_ten_staggered(probe)),
                           np.array(_area_five_to_ten_regular(probe)))


# --------------------------------------------------------------------------
# Path selection
# --------------------------------------------------------------------------

@pytest.mark.parametrize('region,location,path,points', [
    (1, TIP, 'nineteen_point_tip', 19),
    (2, ROOT, 'nineteen_point_root', 19),
    (3, TIP, 'ten_point_tip', 10),
    (4, ROOT, 'ten_point_root', 10),
    (5, TIP, 'ten_point_tip', 10),
    (6, ROOT, 'ten_point_root', 10),
    (7, TIP, 'nineteen_point_tip', 19),
    (8, ROOT, 'nineteen_point_root', 19),
])
def test_region_and_location_select_the_path(region, location, path, points):
    result = calculate_ptcp(0.5, region, location, 0.45, 1.2,
                            tan_te=0.10, tan_hinge_line=0.20)
    assert result['path'] == path
    assert result['points'] == points
    assert result['stations'].size == points
    assert result['area'].size == points
    assert result['moment'].size == points


def test_regions_three_and_four_take_the_trailing_edge_generator():
    """Regions 5 and 6 take the hinge line instead."""
    trailing = calculate_ptcp(0.5, 3, TIP, 0.45, 1.2, tan_te=0.10,
                              tan_hinge_line=0.20)
    hinge = calculate_ptcp(0.5, 5, TIP, 0.45, 1.2, tan_te=0.10,
                           tan_hinge_line=0.20)
    assert trailing['generator'] == pytest.approx(0.10 / 1.2)
    assert hinge['generator'] == pytest.approx(0.20 / 1.2)
    assert trailing['centre_of_pressure'] != pytest.approx(
        hinge['centre_of_pressure'])


def test_nineteen_point_regions_ignore_both_generators():
    """CASENO outside 3 to 6 never forms TGEN."""
    base = calculate_ptcp(0.5, 1, TIP, 0.45, 1.2, tan_te=0.10,
                          tan_hinge_line=0.20)
    moved = calculate_ptcp(0.5, 1, TIP, 0.45, 1.2, tan_te=0.90,
                           tan_hinge_line=0.80)
    assert base['generator'] == 0.0
    assert base['centre_of_pressure'] == pytest.approx(
        moved['centre_of_pressure'])


def test_regions_one_two_seven_and_eight_agree_with_each_other():
    """They differ only in the ONLY10 test, which excludes all four."""
    values = [calculate_ptcp(0.5, region, TIP, 0.45, 1.2)['pressure_ratio']
              for region in (1, 2, 7, 8)]
    assert values[0] == pytest.approx(values[1]) == pytest.approx(values[2])
    assert values[0] == pytest.approx(values[3])


# --------------------------------------------------------------------------
# Behaviour
# --------------------------------------------------------------------------

def test_cumulative_sums_are_the_running_integrals():
    result = calculate_ptcp(0.5, 1, TIP, 0.45, 1.2)
    np.testing.assert_allclose(result['cumulative_area'],
                               np.cumsum(result['area']))
    np.testing.assert_allclose(result['cumulative_moment'],
                               np.cumsum(result['moment']))
    np.testing.assert_allclose(
        result['pressure_distribution'],
        result['cumulative_area'] / result['stations'])


def test_outputs_are_read_off_the_distributions_at_the_station():
    """At a tabulated station the lookup returns that entry exactly."""
    result = calculate_ptcp(0.4, 1, TIP, 0.45, 1.2)
    index = int(np.argmin(np.abs(result['stations'] - 0.4)))
    assert result['pressure_ratio'] == pytest.approx(
        result['pressure_distribution'][index])
    assert result['centre_of_pressure'] == pytest.approx(
        result['centre_distribution'][index])


def test_pressure_and_centre_of_pressure_stay_physical():
    for region in range(1, 9):
        for location in (TIP, ROOT):
            result = calculate_ptcp(0.5, region, location, 0.45, 1.2,
                                    tan_te=0.10, tan_hinge_line=0.20)
            assert 0.0 < result['pressure_ratio'] < 1.0
            assert 0.0 < result['centre_of_pressure'] < 1.0
            assert np.all(np.isfinite(result['pp']))


def test_tip_and_root_differ_on_every_region():
    for region in range(1, 9):
        tip = calculate_ptcp(0.5, region, TIP, 0.45, 1.2, tan_te=0.10,
                             tan_hinge_line=0.20)
        root = calculate_ptcp(0.5, region, ROOT, 0.45, 1.2, tan_te=0.10,
                              tan_hinge_line=0.20)
        assert tip['pressure_ratio'] != pytest.approx(
            root['pressure_ratio'])


def test_ten_point_tip_and_root_use_opposite_generator_signs():
    """Tip forms SUMA + TGEN*SUMM, root forms SUMA - TGEN*SUMM."""
    for location, sign in ((TIP, +1.0), (ROOT, -1.0)):
        result = calculate_ptcp(0.5, 3, location, 0.45, 1.2, tan_te=0.10)
        area = result['cumulative_area']
        moment = result['cumulative_moment']
        expected = ((area - moment) /
                    (area + sign * result['generator'] * moment))
        np.testing.assert_allclose(result['centre_distribution'], expected)


def test_nineteen_point_centre_of_pressure_has_no_generator_term():
    result = calculate_ptcp(0.5, 1, TIP, 0.45, 1.2)
    area = result['cumulative_area']
    moment = result['cumulative_moment']
    np.testing.assert_allclose(result['centre_distribution'],
                               area / (area + moment))


# --------------------------------------------------------------------------
# Rejections
# --------------------------------------------------------------------------

def test_zero_beta_is_rejected():
    with pytest.raises(ValueError, match="beta"):
        calculate_ptcp(0.5, 1, TIP, 0.45, 0.0)


@pytest.mark.parametrize('region', [0, 9, -1])
def test_region_outside_one_to_eight_is_rejected(region):
    with pytest.raises(ValueError, match="outside"):
        calculate_ptcp(0.5, region, TIP, 0.45, 1.2)
