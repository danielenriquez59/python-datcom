"""
Regression tests for the WBTCDO translation.

Source of truth: datcom-legacy/datcom_2000/wbtcdo.f.

The figure lookups go through the already-tested INTERX, so these cover the
four-stage chain, the anchor recovery that makes the fairing correct, and
the branch and early-return behaviour.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.wbtcdo import (
    calculate_wbtcdo, calculate_drag_divergence_mach,
    _PARMA, _XUNIT, _PARMB, _YUNIT, _PARMC, _ZUNIT, _PARMD, _TMD,
    _SUBSONIC_ANCHOR, _SUPERSONIC_ANCHOR, _CD11_FROM_CD14,
)

_GEOMETRY = dict(sweep_c4_deg=30.0, thickness_ratio=0.06,
                 aspect_ratio=6.0, taper_ratio=0.5,
                 configuration_type=2.0)
_LEVELS = dict(cd_mach06=0.020, cd_mach07=0.021, cd_mach14=0.030)


def _run(mach, **kwargs):
    params = dict(_GEOMETRY, **_LEVELS)
    params.update(kwargs)
    return calculate_wbtcdo(mach, **params)


# --------------------------------------------------------------------------
# Table shapes
# --------------------------------------------------------------------------

def test_figure_part_shapes_match_the_source_calls():
    assert len(_PARMA) == 16 and len(_XUNIT) == 56    # LGH (7,8), LIND 8
    assert len(_PARMB) == 20 and len(_YUNIT) == 50    # LGH (5,10), LIND 10
    assert len(_PARMC) == 6 and len(_ZUNIT) == 9      # LGH (3,3), LIND 3
    assert len(_PARMD) == 26 and len(_TMD) == 39      # LGH (3,13), LIND 13


def test_grid_padding_between_packed_columns():
    """Each PARM array pads to its LIND stride before the next grid."""
    assert _PARMA[:7] == [0., 10., 20., 30., 40., 50., 60.]
    assert _PARMA[7] == 0.                     # padding to stride 8
    assert _PARMB[:5] == [2., 3., 4., 6., 8.]
    assert _PARMB[5:10] == [0.] * 5            # padding to stride 10
    assert _PARMD[:3] == [1., 2., 3.]
    assert _PARMD[3:13] == [0.] * 10           # padding to stride 13


# --------------------------------------------------------------------------
# The four-stage chain
# --------------------------------------------------------------------------

def test_each_stage_feeds_the_next():
    result = calculate_drag_divergence_mach(**_GEOMETRY)
    for key in ('stage_a', 'stage_b', 'stage_c', 'md'):
        assert np.isfinite(result[key])
    assert 0.5 < result['md'] < 1.2


def test_divergence_mach_rises_with_sweep():
    """Sweeping the wing delays drag divergence."""
    values = [calculate_drag_divergence_mach(
        **dict(_GEOMETRY, sweep_c4_deg=s))['md'] for s in (0., 20., 40., 60.)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_divergence_mach_falls_with_thickness():
    """A thicker wing diverges earlier."""
    values = [calculate_drag_divergence_mach(
        **dict(_GEOMETRY, thickness_ratio=t))['md']
        for t in (0.03, 0.06, 0.10)]
    assert all(b < a for a, b in zip(values, values[1:]))


# --------------------------------------------------------------------------
# The fairing
# --------------------------------------------------------------------------

def test_subsonic_anchor_is_recovered_exactly():
    """CDO equals CD7 at Mach 0.7, which pins the whole subsonic branch."""
    assert _run(_SUBSONIC_ANCHOR)['cdo'] == pytest.approx(0.021, abs=1e-12)


def test_supersonic_anchor_is_recovered_exactly():
    """CDO equals CD11 at Mach 1.1."""
    result = _run(_SUPERSONIC_ANCHOR)
    assert result['cdo'] == pytest.approx(result['cd_mach11'], abs=1e-12)


def test_cd11_is_derived_from_cd14_when_absent():
    result = _run(0.8)
    assert result['cd_mach11'] == pytest.approx(_CD11_FROM_CD14 * 0.030)


def test_supplied_cd11_is_used_directly():
    result = _run(0.8, cd_mach11=0.055)
    assert result['cd_mach11'] == pytest.approx(0.055)


def test_branch_switches_at_the_divergence_mach():
    md = calculate_drag_divergence_mach(**_GEOMETRY)['md']
    assert not _run(md - 0.01)['supersonic']
    assert _run(md + 0.01)['supersonic']
    assert _run(md - 0.01)['anchor'] == _SUBSONIC_ANCHOR
    assert _run(md + 0.01)['anchor'] == _SUPERSONIC_ANCHOR


def test_subsonic_drag_rises_gently_towards_divergence():
    md = calculate_drag_divergence_mach(**_GEOMETRY)['md']
    values = [_run(m)['cdo'] for m in (0.6, 0.7, 0.8, md - 0.02)]
    assert all(b > a for a, b in zip(values, values[1:]))
    assert values[-1] < 0.03


def test_subsonic_exponent_is_near_unity():
    """Which is why the subsonic branch stays well behaved."""
    assert 0.5 < _run(0.8)['exponent'] < 2.0


def test_supersonic_exponent_is_negative_for_a_peaked_drag_curve():
    """Drag peaks near Mach 1.1 and falls by 1.4, so DCD11 is negative.

    The resulting negative exponent makes the supersonic branch diverge as
    the query Mach approaches MD from above. That is the source's formula,
    which is built to be read meaningfully above the divergence Mach.
    """
    assert _run(1.2)['exponent'] < 0.0
    md = calculate_drag_divergence_mach(**_GEOMETRY)['md']
    near = _run(md + 0.005)['cdo']
    far = _run(1.2)['cdo']
    assert near > far


def test_early_return_when_no_supersonic_level_is_known():
    """The source leaves CDO untouched; this reports None rather than stale."""
    result = _run(1.2, cd_mach14=None, cd_mach11=None)
    assert result['cdo'] is None
    assert 'supersonic drag level' in result['reason']


def test_missing_cd14_below_divergence_is_rejected():
    with pytest.raises(ValueError, match="CD14"):
        _run(0.8, cd_mach14=None, cd_mach11=None)
