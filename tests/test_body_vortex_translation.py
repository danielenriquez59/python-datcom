"""
Regression tests for the DATCOM Section 4.3.1.3 body vortex translation.

Source of truth: datcom-legacy/datcom_2000/bodowg.f, ali.f, getmax.f, and
the caller in wbclb.f that defines the gating and the argument mapping.
"""

import numpy as np
import pytest

from pydatcom.interactions.body_vortex import (
    getmax, ali, fig4313_13a, fig4313_13b, fig4313_14, fig4313_15,
    calculate_bodowg, body_vortex_lift_increment,
)
from pydatcom.aerodynamics.wing_body_tail import calculate_clwbt


# --------------------------------------------------------------------------
# Figure tables: bodowg.f DATA statements
# --------------------------------------------------------------------------

_FIG_13A = list(zip(
    [6.8, 7.2, 7.5, 8., 8.4, 9., 9.6, 10.4, 11.2, 12.1, 13.4, 15., 16.,
     17.1, 18., 20.],
    [20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8.5, 8.,
     7.7, 7.]))

_FIG_13B = list(zip(
    [0., .5, 1., 2., 2.5, 3.3, 4., 4.7, 5.5, 6., 7.],
    [.86, 1.02, 1.2, 1.52, 1.65, 1.8, 1.92, 2., 2.09, 2.14, 2.23]))

_FIG_14 = list(zip(
    [0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.],
    [.5, .57, .62, .66, .69, .71, .72, .735, .75, .755, .76, .765, .765,
     .765, .765]))

_FIG_15 = list(zip(
    [0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.],
    [.4, .5, .6, .7, .81, .91, 1.02, 1.13, 1.24, 1.35, 1.47, 1.6, 1.74,
     1.88, 2.02]))


@pytest.mark.parametrize("x,y", _FIG_13A)
def test_fig4313_13a_source_coordinates(x, y):
    assert fig4313_13a(x) == pytest.approx(y, abs=1e-12)


@pytest.mark.parametrize("x,y", _FIG_13B)
def test_fig4313_13b_source_coordinates(x, y):
    assert fig4313_13b(x) == pytest.approx(y, abs=1e-12)


@pytest.mark.parametrize("x,y", _FIG_14)
def test_fig4313_14_source_coordinates(x, y):
    """Includes the FORTRAN '4*.765' repeat-count expansion."""
    assert fig4313_14(x) == pytest.approx(y, abs=1e-12)


@pytest.mark.parametrize("x,y", _FIG_15)
def test_fig4313_15_source_coordinates(x, y):
    assert fig4313_15(x) == pytest.approx(y, abs=1e-12)


def test_separation_station_moves_forward_with_alpha():
    """Figure 4.3.1.3-13A: the vortex separates earlier at higher alpha."""
    values = [fig4313_13a(a) for a in np.linspace(6.8, 20.0, 80)]
    assert np.all(np.diff(values) <= 1e-12)


def test_vortex_strength_grows_downstream():
    """Figure 4.3.1.3-15: circulation accumulates with distance."""
    values = [fig4313_15(x) for x in np.linspace(0.0, 7.0, 80)]
    assert np.all(np.diff(values) >= -1e-12)


# --------------------------------------------------------------------------
# GETMAX
# --------------------------------------------------------------------------

def test_getmax_returns_maximum_and_station():
    assert getmax([0., 1., 2., 3.], [1., 5., 3., 2.]) == (1.0, 5.0, 1)


def test_getmax_keeps_first_occurrence():
    """The source adopts a later station only when strictly greater."""
    assert getmax([0., 1., 2., 3.], [1., 5., 5., 2.])[2] == 1


def test_getmax_rejects_empty_or_mismatched():
    with pytest.raises(ValueError):
        getmax([], [])
    with pytest.raises(ValueError):
        getmax([0., 1.], [1.])


# --------------------------------------------------------------------------
# ALI
# --------------------------------------------------------------------------

def test_ali_is_finite_for_a_valid_vortex():
    value = ali(0.9, 0.8, 5.0, 0.8, 0.5)
    assert np.isfinite(value)


def test_ali_rejects_degenerate_geometry():
    """A zero exposed semispan and a taper ratio of -1 are both singular."""
    with pytest.raises(ValueError):
        ali(0.9, 0.8, 1.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        ali(0.9, 0.8, 5.0, 0.8, -1.0)


def test_ali_rejects_vortex_on_the_body_surface():
    """h=0 with the vortex at the body radius divides by zero in the source."""
    with pytest.raises(ValueError, match="body surface"):
        ali(0.0, 0.8, 5.0, 0.8, 0.5)


def test_ali_handles_zero_height_away_from_the_surface():
    """The H=0 branch of the source is reachable and finite."""
    assert np.isfinite(ali(0.0, 2.0, 5.0, 0.8, 0.5))


def test_ali_decays_as_the_vortex_moves_away():
    """A vortex farther above the surface plane interferes less."""
    close = abs(ali(1.0, 1.2, 5.0, 1.0, 0.5))
    far = abs(ali(6.0, 1.2, 5.0, 1.0, 0.5))
    assert far < close


# --------------------------------------------------------------------------
# BODOWG
# --------------------------------------------------------------------------

def test_bodowg_is_inactive_below_six_degrees():
    """The source skips the whole vortex system below |alpha| = 6."""
    for alpha in (0.0, 2.0, 5.9, -5.9):
        result = calculate_bodowg(alpha, 60.0, 3.0, 6.0, 0.5)
        assert not result['active']
        assert result['ivbw'] == 0.0
        assert result['go2pav'] == 0.0


def test_bodowg_activates_above_the_cutoff():
    result = calculate_bodowg(14.0, 60.0, 3.0, 6.0, 0.5)
    assert result['active']
    assert result['go2pav'] > 0.0


def test_bodowg_uses_absolute_alpha():
    """Vortex effects are symmetric in the sign of the angle of attack."""
    positive = calculate_bodowg(14.0, 60.0, 3.0, 6.0, 0.5)
    negative = calculate_bodowg(-14.0, 60.0, 3.0, 6.0, 0.5)
    assert positive['ivbw'] == pytest.approx(negative['ivbw'])
    assert positive['go2pav'] == pytest.approx(negative['go2pav'])


def test_bodowg_inactive_when_vortex_has_not_reached_the_surface():
    """XD <= 0 means separation occurs aft of the surface; the term is off."""
    result = calculate_bodowg(12.0, 10.0, 1.0, 5.0, 0.5)
    assert result['xd'] <= 0.0
    assert not result['active']


def test_bodowg_rejects_nonpositive_body_radius():
    with pytest.raises(ValueError):
        calculate_bodowg(14.0, 60.0, 0.0, 6.0, 0.5)


def test_bodowg_strength_grows_with_alpha():
    """Farther downstream of separation means a stronger vortex."""
    strengths = [calculate_bodowg(a, 60.0, 3.0, 6.0, 0.5)['go2pav']
                 for a in (8.0, 12.0, 16.0, 20.0)]
    assert all(b > a for a, b in zip(strengths, strengths[1:]))


# --------------------------------------------------------------------------
# The WBCLB / CLWBT vortex term
# --------------------------------------------------------------------------

def _fin_state(sspne=3.0):
    """A body-dominated fin, r/s = (SSPN-SSPNE)/SSPN."""
    return {
        'htail_type': 1.0, 'htail_chrdr': 3.0, 'htail_chrdtp': 1.5,
        'htail_sspn': 6.0, 'htail_sspne': sspne,
        'htail_savsi': 0.0, 'htail_chstat': 0.25,
        'synths_xh': 60.0,
    }


def test_vortex_term_gated_below_one_third_span_ratio():
    """WBCLB applies the term only when FACT(1) >= 1/3."""
    result = body_vortex_lift_increment(_fin_state(sspne=5.4), 14.0, 0.07)
    assert result['ratio'] == pytest.approx(0.1)
    assert result['gated']
    assert result['increment'] == 0.0


def test_vortex_term_active_above_one_third_span_ratio():
    result = body_vortex_lift_increment(_fin_state(sspne=3.0), 14.0, 0.07)
    assert result['ratio'] == pytest.approx(0.5)
    assert not result['gated']
    assert result['increment'] != 0.0


def test_vortex_term_is_zero_at_low_alpha_even_when_ungated():
    """A large body still contributes nothing below the six-degree cutoff."""
    result = body_vortex_lift_increment(_fin_state(sspne=3.0), 4.0, 0.07)
    assert not result['gated']
    assert result['increment'] == 0.0


def test_vortex_term_reduces_surface_lift_at_high_alpha():
    """Body vortices induce downwash on the aft surface, cutting its load."""
    increments = [body_vortex_lift_increment(
        _fin_state(sspne=3.0), a, 0.07)['increment']
        for a in (8.0, 12.0, 16.0, 20.0)]
    assert all(value < 0.0 for value in increments)
    assert all(b < a for a, b in zip(increments, increments[1:]))


def test_clwbt_accepts_the_vortex_term():
    """The increment enters CLWBT inside the QOQI-scaled bracket."""
    without = calculate_clwbt(0.4, 0.3, 0.06, 0.0, 1.0, 0.9)
    with_vortex = calculate_clwbt(0.4, 0.3, 0.06, 0.0, 1.0, 0.9,
                                  vortex_increment=-0.05)
    assert (with_vortex['cl_total'] - without['cl_total'] ==
            pytest.approx(-0.05 * 0.9))
    assert without['vortex_increment'] == 0.0
