"""
Regression tests for the CLRDER translation.

Source of truth: datcom-legacy/datcom_2000/clrder.f.

The figure lookups go through the already-tested INTERX, so these cover the
tables CLRDER owns, the Section 7.1.3.2 compressibility correction, the
sweep-bracket interpolation, and the panel increments.  The whole routine
is also checked against a compiled probe (tools/probes/clrder.py).
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.clrder import (
    calculate_clr_wing, calculate_clr_panel_increment, calculate_clrder,
    _FIG_71320_10_AR, _FIG_71320_10_TAPER, _FIG_71320_10_DEP,
    _FIG_71320_11_AR, _FIG_71320_11_TAPER, _FIG_71320_11_DEP,
    _SWEEP_GRID, _UNITS, _UNITI,
)
from pydatcom.utils.legacy_interp import interx


def _wing(**kwargs):
    params = dict(cl=[0.0, 0.2, 0.4, 0.6], aspect_ratio=6.0,
                  taper_ratio=0.5, sweep_c4_deg=30.0,
                  mach=0.3, dihedral_deg=0.0, twist_deg=0.0)
    params.update(kwargs)
    return calculate_clr_wing(**params)


# --------------------------------------------------------------------------
# Figure tables
# --------------------------------------------------------------------------

def test_figure_tables_have_source_shapes():
    """LIND and LDEP from the source INTERX calls."""
    assert len(_FIG_71320_10_AR) == 10 and len(_FIG_71320_10_TAPER) == 4
    assert len(_FIG_71320_10_DEP) == 40
    assert len(_FIG_71320_11_AR) == 9 and len(_FIG_71320_11_TAPER) == 4
    assert len(_FIG_71320_11_DEP) == 36


def test_figure_71320_10_reproduces_every_source_entry():
    """All 40 entries, with the aspect ratio varying fastest."""
    for j, taper in enumerate(_FIG_71320_10_TAPER):
        for i, ar in enumerate(_FIG_71320_10_AR):
            value = interx(2, [_FIG_71320_10_AR, _FIG_71320_10_TAPER],
                           [ar, taper], [10, 4], _FIG_71320_10_DEP,
                           lx1l=1, lx2l=1, lx1u=1, lx2u=1)
            assert value == pytest.approx(
                _FIG_71320_10_DEP[j * 10 + i], abs=1e-12), (
                    f"mismatch at AR={ar}, taper={taper}")


def test_figure_71320_11_sentinel_column_repeats():
    """The 99. taper column repeats the 0.4 column, flattening the table."""
    assert _FIG_71320_11_DEP[18:27] == _FIG_71320_11_DEP[27:36]
    assert _FIG_71320_11_TAPER[-1] == 99.0


def test_figure_71320_10_increases_with_aspect_ratio():
    """Each taper column rises monotonically with aspect ratio."""
    for j in range(4):
        column = _FIG_71320_10_DEP[j * 10:(j + 1) * 10]
        assert all(b >= a for a, b in zip(column, column[1:]))


# --------------------------------------------------------------------------
# Section 7.1.3.2 correction and sweep bracket
# --------------------------------------------------------------------------

def test_compressibility_correction_is_unity_at_zero_mach():
    """At M=0 the numerator and denominator of CON coincide exactly."""
    result = _wing(mach=0.0)
    assert result['compressibility_correction'] == pytest.approx(1.0, abs=1e-12)


def test_compressibility_correction_grows_with_mach():
    corrections = [_wing(mach=m)['compressibility_correction']
                   for m in (0.0, 0.3, 0.6, 0.85)]
    assert all(b > a for a, b in zip(corrections, corrections[1:]))


def test_sweep_bracket_selects_the_source_interval():
    """The source keeps the last grid point at or below the wing sweep."""
    for sweep, expected in ((0.0, 0), (10.0, 0), (15.0, 1), (44.0, 2),
                            (45.0, 3), (59.0, 3)):
        assert _wing(sweep_c4_deg=sweep)['sweep_index'] == expected


def test_sweep_beyond_the_grid_uses_the_final_interval():
    """A sweep past 60 degrees must not index past the coefficient arrays."""
    result = _wing(sweep_c4_deg=70.0)
    assert result['sweep_index'] == len(_SWEEP_GRID) - 2
    assert np.all(np.isfinite(result['clr']))


def test_unswept_wing_uses_the_first_interval_exactly():
    """At zero sweep the offset is the full 15-degree step."""
    # CON is exactly 1 only at zero Mach, not merely at zero sweep: the
    # numerator keeps an AR*(1-BEE^2) term that survives an unswept wing.
    result = _wing(sweep_c4_deg=0.0, mach=0.0)
    assert result['compressibility_correction'] == pytest.approx(1.0)
    unit = result['figure_10_unit']
    expected_zero = -_UNITI[0] - _UNITS[0] * unit
    assert result['clr_per_cl'] == pytest.approx(-expected_zero, rel=1e-9)


# --------------------------------------------------------------------------
# The assembled derivative
# --------------------------------------------------------------------------

def test_clr_is_linear_in_lift_coefficient():
    result = _wing(cl=[0.0, 0.2, 0.4, 0.6])
    steps = np.diff(result['clr'])
    assert np.allclose(steps, steps[0])


def test_clr_per_cl_grows_with_aspect_ratio():
    values = [_wing(aspect_ratio=ar)['clr_per_cl'] for ar in (4.0, 6.0, 8.0)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_dihedral_and_twist_enter_additively():
    """Both terms are independent of CL and add to the lift-dependent part."""
    plain = _wing(cl=[0.0])['clr'][0]
    dihedral = _wing(cl=[0.0], dihedral_deg=5.0)['clr'][0]
    twist = _wing(cl=[0.0], twist_deg=-2.0)['clr'][0]
    both = _wing(cl=[0.0], dihedral_deg=5.0, twist_deg=-2.0)['clr'][0]
    assert both == pytest.approx(dihedral + twist - plain)


def test_unswept_wing_has_no_dihedral_term():
    """The dihedral factor carries sin(sweep), so it vanishes when unswept."""
    assert _wing(sweep_c4_deg=0.0)['dihedral_factor'] == pytest.approx(0.0)


def test_forward_swept_wing_is_rejected():
    """The source leaves CLR untouched rather than computing it."""
    with pytest.raises(ValueError, match="forward-swept"):
        _wing(sweep_c4_deg=-10.0)


def test_supersonic_input_is_rejected():
    with pytest.raises(ValueError, match="subsonic"):
        _wing(mach=1.2)


# --------------------------------------------------------------------------
# Panel increments
# --------------------------------------------------------------------------

def test_panel_increment_matches_source_expression():
    alpha, cyb, lp, zp, b = 8.0, -0.008, 20.0, 6.0, 30.0
    result = calculate_clr_panel_increment([alpha], cyb, lp, zp, b)
    sin_a, cos_a = np.sin(np.deg2rad(alpha)), np.cos(np.deg2rad(alpha))
    along = lp * cos_a + zp * sin_a
    across = zp * cos_a - lp * sin_a
    assert result['dclr'][0] == pytest.approx(
        -2.0 * cyb * along * across / b**2, rel=1e-12)
    assert result['dclp'][0] == pytest.approx(
        2.0 * cyb * across * (across - zp) / b**2, rel=1e-12)


def test_panel_roll_damping_vanishes_at_zero_alpha():
    """At alpha=0 the across-arm equals ZP, so dCLP is zero."""
    result = calculate_clr_panel_increment([0.0], -0.008, 20.0, 6.0, 30.0)
    assert result['dclp'][0] == pytest.approx(0.0, abs=1e-15)


def test_panel_increment_rejects_bad_reference_length():
    with pytest.raises(ValueError):
        calculate_clr_panel_increment([0.0], -0.008, 20.0, 6.0, 0.0)


# --- The whole routine, checked against a compiled probe -------------------
# tools/probes/clrder.py runs its cases in one program, so the replay
# carries the saved horizontal-tail carryover factors AKHB and AKBH.

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures'
                     / 'probes' / 'clrder.json').read_text())


def _replay():
    state, out = {}, []
    for p in _PROBE:
        r = calculate_clrder(dict(p['inputs'], state=state))
        state = r['state']
        out.append(r)
    return out


_WHOLE = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_whole_routine_matches_compiled_clrder(case):
    for name, values in _PROBE[case]['outputs'].items():
        lo = 201 if len(values) == 180 else 361
        ours = [_WHOLE[case][name.lower()][lo + i]
                for i in range(len(values))]
        np.testing.assert_allclose(ours, values, rtol=1e-9, atol=1e-14,
                                   err_msg=name)


def test_ventral_fin_damping_uses_the_tail_derivative():
    """VF(J+280) is formed with the vertical tail's DCYBV."""
    c = _PROBE[0]['inputs']
    r = _WHOLE[0]
    lpf, zpf = c['stbh']['11'], c['stbh']['12']
    ca, sa = np.cos(c['alpha'][0] * .01745329), \
        np.sin(c['alpha'][0] * .01745329)
    arm = zpf * ca - lpf * sa
    expected = 2. * c['vt']['141'] * arm * (arm - zpf) / c['blref'] ** 2
    assert r['vf'][281] == pytest.approx(expected)
