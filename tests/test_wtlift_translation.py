"""
Regression tests for the WTLIFT and CLMXBS translations.

Source of truth: datcom-legacy/datcom_2000/wtlift.f, clmxbs.f.

Every output is checked against a compiled probe of the legacy routines
(tools/probes/wtlift.py), and every embedded table against a re-parse of
the FORTRAN DATA statements.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import wtlift as module
from pydatcom.aerodynamics.wtlift import calculate_clmxbs, calculate_wtlift

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'wtlift.json').read_text())

_OUTPUTS = {'CLA': 'cla', 'B43': 'alpha_clmax', 'B44': 'clmax',
            'A144': 'a144', 'A145': 'a145', 'A146': 'a146', 'A159': 'a159',
            'A160': 'a160', 'A171': 'a171', 'A172': 'a172'}


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    """Every quantity WTLIFT writes, against execution.

    A quantity the source leaves untouched on a branch is zero in the probe,
    whose driver clears every array before each case, and absent here.
    """
    probe = _PROBE[case]
    result = calculate_wtlift(**probe['inputs'])
    for tag, key in _OUTPUTS.items():
        expected = probe['outputs'][tag]
        if key in result:
            assert result[key] == pytest.approx(expected, rel=1e-9,
                                                abs=1e-12), tag
        else:
            assert expected == 0.0, f"{tag} set by the source but not here"


def test_probe_reaches_every_branch():
    """Guard the fixture: it keeps reaching each path through the routine."""
    results = [calculate_wtlift(**p['inputs']) for p in _PROBE]
    computed = [r for r in results if r['computed']]
    low = [r for r in computed if r['low_aspect_ratio']]
    assert any(not r['low_aspect_ratio'] for r in computed)
    assert {r['clmxbs']['part'] for r in low} == {'A', 'B'}
    assert any(r['a160'] <= 4.5 for r in low)
    assert any(r['a160'] > 4.5 for r in low)
    assert any('cranked_ratio' in r for r in computed)
    assert any('a171' in r for r in computed)
    assert any(not r['computed'] for r in results)
    # The curved planform's diagnostic is what the source actually printed.
    curved = [p for p in _PROBE if p['inputs']['planform_type'] == 4.0]
    assert curved and 'NO CLALPHA COMPUTATION' in curved[0]['outputs']['_text'][0]


@pytest.mark.parametrize("stem,names", [
    ('clmxbs', ['C1ABC', 'DYAG', 'CBASE', 'C2A', 'AMN', 'DE']),
    ('wtlift', ['TR', 'C2', 'SALE', 'DELTAY', 'CLL', 'DY', 'DACLL', 'DYA',
                'AMACH', 'SALE4', 'DCAR', 'C1ABCS', 'ACLMX', 'DMN', 'C1TAB',
                'DACL', 'ACLE', 'C1TABO', 'DACLO', 'BA', 'CLOVCL']),
])
def test_every_table_value_matches_the_source(stem, names):
    source = parse(stem)
    for name in names:
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


def test_table_sizes_match_their_declarations():
    """The DATA counts fill each declared array exactly, unlike WINGCL's."""
    assert module._CBASE.size == 19 * 12
    assert module._DE.size == 15 * 3
    assert module._CLL.size == 13 * 7
    assert module._DACLL.size == 13 * 4
    assert module._DCAR.size == 5 * 6 * 4
    assert module._DACL.size == 20 * 3
    assert module._DACLO.size == 10 * 10


def _straight(aspect_ratio=6.0, mach=0.3, sweep=0.0):
    return dict(
        planform_type=1.0,
        geometry=dict(area=150.0, aspect_ratio=aspect_ratio, taper_ratio=0.5,
                      sweep_le_deg=sweep, cos_le=np.cos(np.radians(sweep)),
                      tan_le=np.tan(np.radians(sweep)),
                      tan_c2=np.tan(np.radians(sweep)) * 0.8,
                      arclss_classified=1.3, arclss_ratio=3.0),
        section=dict(deltay=2.0, xovc=0.3, cla=0.1, clmax=1.4),
        flight=dict(mach=mach, beta=np.sqrt(1 - mach**2), alpha_zero_lift=-2.0),
        sref=150.0)


def test_lift_slope_approaches_the_section_slope_for_a_long_wing():
    """Helmbold-Diederich: CLa -> cla as A grows, and 2*pi*A/2 as A -> 0."""
    long = calculate_wtlift(**_straight(aspect_ratio=500.0))['cla']
    assert long == pytest.approx(0.1, rel=0.01)
    stubby = calculate_wtlift(**_straight(aspect_ratio=0.05))['cla']
    assert stubby == pytest.approx(np.pi * 0.05 / 2 * np.pi / 180, rel=0.01)


def test_mach_enters_only_through_the_swept_term():
    """DATCOM's kappa is the Mach-dependent section slope over 2*pi/beta, so
    ``A*beta/kappa`` reduces to ``2*pi*A/cla``: beta cancels there and
    survives only in ``tan^2/beta^2``.  Compressibility otherwise arrives
    through the section slope A(131), the per-Mach CLALPA input.
    """
    unswept = calculate_wtlift(**_straight())['cla']
    assert calculate_wtlift(**_straight(mach=0.6))['cla'] == pytest.approx(
        unswept, rel=1e-14)
    swept = calculate_wtlift(**_straight(sweep=40.0))['cla']
    assert swept < unswept
    # At a fixed section slope the swept term alone lowers CLa with Mach...
    fixed = _straight(sweep=40.0, mach=0.6)
    assert calculate_wtlift(**fixed)['cla'] < swept
    # ...while a Prandtl-Glauert section slope, cla0/beta, raises it.
    scaled = _straight(sweep=40.0, mach=0.6)
    scaled['section']['cla'] = 0.1 * np.sqrt(1 - 0.3**2) / np.sqrt(1 - 0.6**2)
    assert calculate_wtlift(**scaled)['cla'] > swept


def test_lift_slope_scales_with_reference_area():
    """CLa*SREF is a physical force slope, independent of SREF."""
    a = _straight()
    b = _straight()
    b['sref'] = 300.0
    assert (calculate_wtlift(**a)['cla'] * 150.0 ==
            pytest.approx(calculate_wtlift(**b)['cla'] * 300.0))


def test_high_aspect_ratio_stall_angle_is_consistent_with_clmax():
    """B(43) = CLMAX/CLa + alpha_0 + delta-alpha, the linear-lift estimate."""
    r = calculate_wtlift(**_straight())
    assert r['alpha_clmax'] == pytest.approx(
        r['clmax'] / r['cla'] - 2.0 + r['a144'])
    assert 10.0 < r['alpha_clmax'] < 25.0
    assert 0.8 < r['clmax'] < 1.6


def test_clmxbs_part_switch_is_at_035():
    a = calculate_clmxbs(1.0, 0.3, 1.0, 0.0, 0.35, 1.0, 1.0)
    b = calculate_clmxbs(1.0, 0.3, 1.0, 0.0, 0.36, 1.0, 1.0)
    assert (a['part'], b['part']) == ('A', 'B')
    # Column DELTAY=0 at C1ABC=1.0 in each part of Figure 4.1.3.4-23.
    assert (a['clmax_base'], b['clmax_base']) == (1.58, 1.37)
