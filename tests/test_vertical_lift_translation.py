"""VTLIFT/VFLIFT branch-logic tests.

The figure tables themselves are covered by
``test_vertical_lift_figures_translation.py``; this file exercises the
branch structure that chooses among them, and pins the three source defects
so that a later "tidy-up" cannot quietly change results.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.vertical_lift import (
    calculate_vtlift, STRAIGHT_TAPERED)
from pydatcom.aerodynamics.vertical_lift_figures import (
    fig4132_61, _FIG_61_GRID, _FIG_61)
from pydatcom.utils.legacy_interp import interx


def _panel(**overrides):
    base = dict(area=40.0, aspect_ratio=1.5, taper_ratio=0.45,
                taper_ratio_exposed=0.45, chrdtp=4.0, chrdbp=6.0,
                sspnop=2.0, delta_y=1.8)
    base.update(overrides)
    return base


def _sweep(**overrides):
    base = dict(sweple_deg=35.0, cosle=0.819, tanle=0.700, tanleo=0.900,
                tantei=0.10, tanteo=0.30, swtei_deg=6.0, swteo_deg=17.0)
    base.update(overrides)
    return base


_SREF = 250.0


# --------------------------------------------------------------------------
# Entry conditions
# --------------------------------------------------------------------------

@pytest.mark.parametrize('mach', [0.0, 0.8, 1.0, np.nan])
def test_subsonic_mach_is_rejected(mach):
    """BETA = sqrt(M^2-1) is real only above Mach 1."""
    with pytest.raises(ValueError, match="supersonic"):
        calculate_vtlift(mach, _panel(), _sweep(), _SREF)


def test_nonpositive_reference_area_is_rejected():
    with pytest.raises(ValueError, match="reference area"):
        calculate_vtlift(2.0, _panel(), _sweep(), 0.0)


def test_non_straight_path_requires_spans():
    with pytest.raises(ValueError, match="SPANS"):
        calculate_vtlift(2.0, _panel(), _sweep(), _SREF, planform_type=3.0)


def test_zero_leading_edge_tangent_is_replaced_not_divided_by():
    """The source guards both TANLEO and TANLEI with .00001 before dividing."""
    result = calculate_vtlift(2.0, _panel(), _sweep(tanle=0.0), _SREF)
    assert np.isfinite(result['cna'])
    # BETA/1e-5 is enormous, so the leading edge reads as supersonic.
    assert result['supersonic_leading_edge'] is True


# --------------------------------------------------------------------------
# Straight tapered path
# --------------------------------------------------------------------------

def test_leading_edge_regime_switches_figures_at_bovert_one():
    """BOVERT = BETA/TANLE crossing 1 swaps figure 4.1.3.2-60A for -60B."""
    sweep = _sweep(tanle=1.0)
    subsonic_le = calculate_vtlift(1.2, _panel(), sweep, _SREF)
    supersonic_le = calculate_vtlift(3.0, _panel(), sweep, _SREF)
    assert subsonic_le['bovert'] < 1.0 < supersonic_le['bovert']
    assert subsonic_le['supersonic_leading_edge'] is False
    assert supersonic_le['supersonic_leading_edge'] is True
    assert subsonic_le['correction_figure'] == '4.1.3.2-60A'
    assert supersonic_le['correction_figure'] == '4.1.3.2-60B'


def test_rectangular_panel_takes_the_rectangular_branch():
    """TRATIP == 1 and zero LE sweep select figures 56G or the closed form."""
    panel = _panel(taper_ratio=1.0)
    sweep = _sweep(sweple_deg=0.0)
    result = calculate_vtlift(2.0, panel, sweep, _SREF)
    assert result['rectangular'] is True
    assert result['slope_figure'] in ('4.1.3.2-56G', 'closed form')
    # Either condition alone is not enough.
    assert calculate_vtlift(2.0, _panel(taper_ratio=1.0), _sweep(),
                            _SREF)['rectangular'] is False
    assert calculate_vtlift(2.0, _panel(), _sweep(sweple_deg=0.0),
                            _SREF)['rectangular'] is False


def test_rectangular_branch_switches_on_aspect_ratio_times_beta():
    """Below AR*BETA = 1 figure 56G is read; above it the closed form runs.

    AR is doubled on entry, so the switch is at 2*AVT(7)*BETA = 1.
    """
    sweep = _sweep(sweple_deg=0.0)
    low = calculate_vtlift(1.02, _panel(taper_ratio=1.0, aspect_ratio=0.05),
                           sweep, _SREF)
    high = calculate_vtlift(3.0, _panel(taper_ratio=1.0, aspect_ratio=1.5),
                            sweep, _SREF)
    assert low['slope_figure'] == '4.1.3.2-56G'
    assert high['slope_figure'] == 'closed form'
    # The closed form is the slender-body limit 4/BETA less a correction.
    assert high['bcna'] == pytest.approx(
        4.0 - 2.0 / (2.0 * 1.5 * np.sqrt(3.0**2 - 1.0)))


def test_normal_force_slope_falls_with_mach():
    """CNA scales roughly as 1/BETA across the supersonic range."""
    values = [calculate_vtlift(m, _panel(), _sweep(), _SREF)['cna']
              for m in (1.5, 2.0, 3.0, 4.0, 5.0)]
    assert all(a > b for a, b in zip(values, values[1:]))
    assert all(v > 0.0 for v in values)


def test_slope_scales_with_panel_area_and_inversely_with_reference():
    base = calculate_vtlift(2.0, _panel(), _sweep(), _SREF)['cna']
    doubled = calculate_vtlift(2.0, _panel(area=80.0), _sweep(),
                               _SREF)['cna']
    halved_ref = calculate_vtlift(2.0, _panel(), _sweep(), _SREF / 2.0)['cna']
    assert doubled == pytest.approx(2.0 * base)
    assert halved_ref == pytest.approx(2.0 * base)


def test_straight_path_reports_no_defect_flags():
    result = calculate_vtlift(2.0, _panel(), _sweep(), _SREF)
    assert result['path'] == 'straight_tapered'
    assert result['spans_reads_wing_block'] is False
    assert result['figure_56a_second_call_uses_truncated_grid'] is False
    assert result['cnt2_divisor_inverted'] is False


# --------------------------------------------------------------------------
# Three-component path
# --------------------------------------------------------------------------

def _cranked(mach=2.0, spans=9.0, **kwargs):
    panel = _panel(**kwargs.pop('panel', {}))
    sweep = _sweep(**kwargs.pop('sweep', {}))
    return calculate_vtlift(mach, panel, sweep, _SREF, spans=spans,
                            planform_type=3.0, **kwargs)


def test_three_components_sum_into_the_total():
    """CNA = RKL*(CNABW + CNAGLV + CNAE)/RAD."""
    result = _cranked()
    assert result['path'] == 'three_component'
    total = (result['rkl'] *
             (result['cnabw'] + result['cnaglv'] + result['cnae']) /
             57.2957795)
    assert result['cna'] == pytest.approx(total)


def test_extension_component_needs_four_degrees_of_trailing_edge_sweep():
    """ABS(SWTEI-SWTEO) < 4 zeroes CNAE outright."""
    near = _cranked(sweep=dict(swtei_deg=10.0, swteo_deg=12.0))
    far = _cranked(sweep=dict(swtei_deg=10.0, swteo_deg=20.0))
    assert near['cnae'] == 0.0
    assert near['cn1'] is None and near['cn2'] is None
    assert far['cnae'] != 0.0
    assert far['cn1'] is not None and far['cn2'] is not None


def test_sharp_and_round_leading_edges_pick_different_curves():
    """KSHARP unset is the source's UNUSED sentinel and means round."""
    sharp = _cranked(ksharp=1.0)
    round_le = _cranked(ksharp=None)
    assert sharp['round_leading_edge'] is False
    assert round_le['round_leading_edge'] is True
    assert sharp['clebw'] != round_le['clebw']
    assert sharp['cna'] != round_le['cna']


def test_glove_geometry_is_a_triangle_on_the_inboard_span():
    """CRGLV = TANLEI*SPANIN makes the glove a delta of aspect ratio 4/tan."""
    result = _cranked()
    spanin = result['spanin']
    tanlei = _sweep()['tanle']
    assert result['sglv'] == pytest.approx(tanlei * spanin**2)
    assert result['arglv'] == pytest.approx(4.0 / tanlei)


def test_basic_wing_branch_switches_on_the_leading_edge_tangents():
    """TANLEO < TANLEI takes the simple form; otherwise the delta subtracts."""
    simple = _cranked(sweep=dict(tanle=1.2, tanleo=0.5))
    subtracting = _cranked(sweep=dict(tanle=0.5, tanleo=1.2))
    assert simple['cnt2'] is None
    assert simple['figure_56a_second_call_uses_truncated_grid'] is False
    assert subtracting['cnt2'] is not None
    assert subtracting['figure_56a_second_call_uses_truncated_grid'] is True


# --------------------------------------------------------------------------
# Source defects, pinned
# --------------------------------------------------------------------------

def test_spans_is_taken_from_the_wing_block():
    """Source defect: SPANS is EQUIVALENCEd to WINGIN(3), not VTIN(3).

    Every other planform quantity in the routine comes from VTIN, the
    vertical panel's own block, and nothing in DATCOM stages panel geometry
    into /WINGI/ before the call.  SPANS drives the basic-wing area, aspect
    ratio and taper, so the wing's semispan sizes the vertical panel.
    """
    result = _cranked(spans=9.0)
    assert result['spans_reads_wing_block'] is True
    # It genuinely drives the geometry, so a different wing changes the fin.
    other = _cranked(spans=14.0)
    assert other['sbw'] != result['sbw']
    assert other['arbw'] != result['arbw']
    assert other['cna'] != result['cna']


def test_second_figure_56a_call_sees_a_truncated_first_grid():
    """Source defect: LGB(1) is left at 12 by the figure 4.1.3.2-62 lookup.

    LGB(1) is set to 23 at label 1110, overwritten with 12 for figure
    4.1.3.2-62, and never restored before the second 56A call.  That caps
    the first-variable axis at 1.0 instead of 30, so every supersonic
    leading edge clamps to the same ordinate.
    """
    result = _cranked(mach=2.0, sweep=dict(tanle=0.5, tanleo=1.2))
    assert result['figure_56a_second_call_uses_truncated_grid'] is True
    # Beyond the truncated grid the lookup clamps, so CNT2 stops responding
    # to Mach even though the untruncated grid runs to 30.
    high = [_cranked(mach=m, sweep=dict(tanle=0.5, tanleo=1.2))['cnt2']
            for m in (2.0, 3.0, 4.0)]
    assert high[0] == pytest.approx(high[1]) == pytest.approx(high[2])


def test_cnt2_divisor_is_inverted_relative_to_every_other_use():
    """Source defect: CNT2 divides by BETA when the LE is subsonic.

    The prevailing convention in the routine is ``CNTHRY=BCNA/TA`` with
    ``IF(SUPLE)CNTHRY=BCNA/BETA`` -- a supersonic leading edge divides by
    BETA.  CNT2 does the opposite.  The consequence is a basic-wing
    component that goes negative at high Mach, where the consistent
    convention keeps it positive.
    """
    result = _cranked(mach=4.0, sweep=dict(tanle=0.5, tanleo=1.2))
    assert result['cnt2_divisor_inverted'] is True
    assert result['cnt2'] != pytest.approx(result['cnt2_consistent'])
    # Source value is the larger one here, and oversubtracts.
    assert result['cnt2'] > result['cnt2_consistent']
    assert result['cnabw'] < 0.0


def test_defect_flags_are_absent_when_their_branches_do_not_run():
    """A panel whose outboard LE is shallower never reaches either defect."""
    result = _cranked(sweep=dict(tanle=1.2, tanleo=0.5))
    assert result['figure_56a_second_call_uses_truncated_grid'] is False
    assert result['cnt2_divisor_inverted'] is False
    assert result['cnt2'] is None and result['cnt2_consistent'] is None


# --------------------------------------------------------------------------
# The figure end modes, which the call sites supply
# --------------------------------------------------------------------------

def test_figure_61_extrapolates_quadratically_off_its_grid():
    """Its call site passes LX1L=2 and LX1U=2, not the INTERX default.

    Inside the grid the modes are irrelevant, which is why they can be
    dropped without any in-range test noticing.  Outside it they are not.
    """
    inside = 1.0
    assert fig4132_61(inside) == pytest.approx(
        interx(1, _FIG_61_GRID, [inside], [12], _FIG_61))
    for outside in (-0.5, 5.0, 30.0):
        clamped = interx(1, _FIG_61_GRID, [outside], [12], _FIG_61)
        assert fig4132_61(outside) != pytest.approx(clamped)
