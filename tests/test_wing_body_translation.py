"""
Regression tests for the subsonic wing-body buildup: WBAERO's wing pass and
WBDRAG, WBLIFT, WBCM and WBCM0.

Source of truth: datcom-legacy/datcom_2000/wbaero.f, wbdrag.f, wblift.f,
wbcm.f, wbcm0.f.

Every output is checked against a compiled probe of the whole chain
(tools/probes/wbaero.py), which runs the legacy WBAERO with BODOWG, GETMAX,
ALI and TABLEC linked in, and every table against a re-parse of the FORTRAN.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import wing_body as module
from pydatcom.aerodynamics.wing_body import (
    NOT_AVAILABLE, calculate_wbaero, calculate_wbcm0, calculate_wbdrag,
)
from pydatcom.interactions import carryover
from pydatcom.utils.legacy_numeric import tbfunx

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'wbaero.json').read_text())

_CURVES = [('CD', 'cd'), ('CL', 'cl'), ('CM', 'cm'), ('CN', 'cn'),
           ('CA', 'ca'), ('CLA', 'cla'), ('CMA', 'cma')]

# Where each translated WB factor sits in /WHWB/, by result group.
_WB_SLOTS = {
    'lift': {'kwb': 2, 'kbw': 3, 'wb4': 4, 'wb5': 5, 'kwb_incidence': 7,
             'kbw_incidence': 8, 'wb9': 9, 'wb10': 10, 'wb11': 11,
             'wb20': 20, 'wb21': 21, 'wb22': 22, 'wb23': 23},
    'moment': {'wb12': 12, 'wb13': 13, 'wb14': 14, 'wb15': 15, 'cm0': 16},
    'drag': {'cd0': 17, 'interference': 18, 'reynolds': 19},
}


def _results():
    return [calculate_wbaero(**p['inputs']) for p in _PROBE]


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    """Every BW curve, the WB factors and the body-vortex FACT entries."""
    probe = _PROBE[case]
    out = probe['outputs']
    result = calculate_wbaero(**probe['inputs'])
    for tag, key in _CURVES:
        np.testing.assert_allclose(result[key], out[tag], rtol=1e-9,
                                   atol=1e-12, err_msg=tag)
    for group, slots in _WB_SLOTS.items():
        for key, slot in slots.items():
            if key in result[group]:
                assert result[group][key] == pytest.approx(
                    out['WB'][slot - 1], rel=1e-9, abs=1e-12), key
    assert out['WB'][5] == pytest.approx(result['drag']['cd0'])  # WB(6)
    n = len(result['cl'])
    vortex = result['vortex']
    assert out['FACT'][0] == pytest.approx(vortex['ratio'], rel=1e-12)
    np.testing.assert_allclose(out['FACT'][1:1 + n], vortex['ivbw'],
                               rtol=1e-9, atol=1e-12)


def test_probe_reaches_every_branch():
    """Guard the fixture: each branch of the chain stays covered."""
    results = _results()
    ratios = [r['lift']['ratio'] for r in results]
    assert min(ratios) <= 0.30 and any(0.3 < x <= 0.8 for x in ratios)
    assert max(ratios) > 0.8
    assert any('source_defect' in r['lift'] for r in results)
    assert any('kwb_incidence' in r['lift'] and 'wb9' in r['lift']
               for r in results)
    # The ellipse fit (beta*A < 4) and both sides of TEMP0 = 1.
    arg = [p['inputs']['surface']['beta'] * p['inputs']['surface']['a7']
           for p in _PROBE]
    assert min(arg) < 4.0 <= max(arg)
    assert {r['moment']['wb15'] == 0.5 for r in results} == {True, False}
    assert {r['moment']['cm0_regression'] for r in results} == {True, False}
    assert any(NOT_AVAILABLE in r['cma'] for r in results)
    assert any(p['inputs']['experimental'] for p in _PROBE)
    assert any(any(r['vortex']['ivbw']) for r in results)
    # A lift curve that turns over, so the inverse lookup is unordered.
    assert any(np.any(np.diff(r['cl']) < 0) for r in results)


@pytest.mark.parametrize("stem,names", [
    ('wbdrag', ['X137', 'X237', 'Y37']),
    ('wblift', ['X10A', 'Y10A', 'X10B', 'Y10B', 'X12A1', 'Y12A1', 'X12A2',
                'Y12A2', 'X12C', 'Y12C', 'XA12', 'XB12', 'Y412B', 'Y412C']),
    ('wbcm', ['X38B', 'Y38B', 'X21C', 'Y21C']),
])
def test_every_table_value_matches_the_source(stem, names):
    source = parse(stem)
    for name in names:
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


def test_wblift_carryover_tables_agree_with_clwbt():
    """WBLIFT and CLWBT digitise Figures 4.3.1.2-10 and -12A separately;
    the two copies are identical."""
    np.testing.assert_array_equal(module._Y10A, carryover._FIG_431210_KWB)
    np.testing.assert_array_equal(module._Y10B, carryover._FIG_431210_KBW)
    np.testing.assert_array_equal(module._Y12A1, carryover._FIG_431212A_KKWB)
    np.testing.assert_array_equal(module._Y12A2, carryover._FIG_431212A_KKBW)


# --------------------------------------------------------------------------
# Behaviour
# --------------------------------------------------------------------------

def test_source_defect_reads_the_supplied_stale_factors():
    """Above d/b = 0.8 the lift is scaled by whatever K factors were left."""
    probe = next(p for p in _PROBE if p['inputs']['stale'])
    inputs = dict(probe['inputs'])
    base = calculate_wbaero(**inputs)
    inputs['stale'] = dict(inputs['stale'], kwb=2.0 * inputs['stale']['kwb'])
    changed = calculate_wbaero(**inputs)
    assert base['lift']['source_defect'] == 'ratio_above_0.8_reads_stale_kwb_kbw'
    assert not np.allclose(base['cl'], changed['cl'])


def test_wing_body_lift_exceeds_the_wing_alone_through_carryover():
    """With a small body, K_W(B)+K_B(W) > 1 lifts the combination above the
    isolated wing by roughly the carryover, plus the small body lift."""
    probe = _PROBE[0]
    result = calculate_wbaero(**probe['inputs'])
    wing = np.array(probe['inputs']['surface_alone']['cl'])
    k = result['lift']['kwb'] + result['lift']['kbw']
    assert k > 1.0
    attached = slice(1, 5)
    assert np.all(result['cl'][attached] > wing[attached])


def test_zero_lift_angle_is_where_the_buildup_lift_vanishes():
    """ALOWB is the inverse lookup of the wing-body CL at zero."""
    probe = _PROBE[0]
    result = calculate_wbaero(**probe['inputs'])
    alpha = np.array(probe['inputs']['alpha_deg'])
    a0 = result['moment']['alpha_zero_lift']
    cl_at, _ = tbfunx(alpha, result['cl'], a0, 1, 1)
    assert cl_at == pytest.approx(0.0, abs=1e-12)


def test_forces_rotate_consistently():
    """CN and CA are the body-axis rotation of CL and CD."""
    result = calculate_wbaero(**_PROBE[0]['inputs'])
    alpha = np.radians(_PROBE[0]['inputs']['alpha_deg'])
    np.testing.assert_allclose(result['cn'] * np.cos(alpha) -
                               result['ca'] * np.sin(alpha), result['cl'],
                               atol=1e-14)


def test_wbdrag_experimental_substitution_needs_both_components():
    base = calculate_wbdrag(0.5, 2e6, 40.0, 0.007, 0.009, 0.0015,
                            [0.0, 0.0], [0.001, 0.002])
    supplied = calculate_wbdrag(
        0.5, 2e6, 40.0, 0.007, 0.009, 0.0015, [0.0, 0.0], [0.001, 0.002],
        {'kbody': True, 'body_cd': [0.02, module.UNUSED],
         'wing_cd': [0.01, 0.01]})
    assert supplied['cd'][0] == pytest.approx(0.03)
    assert supplied['cd'][1] == base['cd'][1]


def test_wbcm0_returns_none_outside_its_range():
    args = dict(aspect_ratio=4.0, tan_le=0.5, tovc=0.06, nose_length=4.0,
                afterbody_length=4.0, taper_ratio=0.5,
                leading_edge_radius=0.008, twist=-0.02, ycm=0.01, cld=0.2,
                reynolds=5e6, tr=0.4, wing_height=0.3, vt=0.0, hd=0.5,
                body_radius_ratio=0.15, mach=0.5)
    assert calculate_wbcm0(**args) is not None
    assert calculate_wbcm0(**dict(args, aspect_ratio=6.5)) is None
    # Positive twist is outside: the upper bound is UNUSED, 1e-30.
    assert calculate_wbcm0(**dict(args, twist=0.01)) is None
    # Reynolds number is clamped to 8e6.
    assert (calculate_wbcm0(**dict(args, reynolds=5e7)) ==
            calculate_wbcm0(**dict(args, reynolds=8e6)))


# --------------------------------------------------------------------------
# The horizontal-tail and body-vertical passes
# --------------------------------------------------------------------------

_TAIL_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'wbaero_tail.json').read_text())


def _tail_pass(inputs):
    t = module.tail_body_inputs(inputs['surface'], inputs['tail_alone'],
                                inputs['synthesis'], inputs['bd63'],
                                inputs['cbarr'])
    return calculate_wbaero(inputs['alpha_deg'], t['surface'],
                            t['surface_alone'], inputs['body'],
                            t['synthesis'], inputs['flight'], inputs['cbarr'],
                            moment_cutoff=False)


@pytest.mark.parametrize("case", range(len(_TAIL_PROBE)))
def test_tail_and_body_vertical_passes_match_compiled_routine(case):
    probe = _TAIL_PROBE[case]
    inputs, out = probe['inputs'], probe['outputs']
    tail = _tail_pass(inputs)
    body = module.calculate_body_vertical(inputs['alpha_deg'], inputs['body'],
                                          inputs['vertical_cd0'])
    for tag, key in _CURVES:
        np.testing.assert_allclose(tail[key], out['BH' + tag], rtol=1e-9,
                                   atol=1e-12, err_msg='BH' + tag)
        np.testing.assert_allclose(body[key], out['BV' + tag], rtol=1e-9,
                                   atol=1e-12, err_msg='BV' + tag)
    for group, slots in _WB_SLOTS.items():
        for key, slot in slots.items():
            if key in tail[group]:
                assert tail[group][key] == pytest.approx(
                    out['HB'][slot - 1], rel=1e-9, abs=1e-12), key
    assert out['FACT'][0] == pytest.approx(tail['vortex']['ratio'])


def test_tail_pass_slopes_through_unavailable_moments():
    """Unlike the wing pass, the tail pass has no CMa cutoff."""
    probe = next(p for p in _TAIL_PROBE
                 if NOT_AVAILABLE in p['inputs']['tail_alone']['cm'])
    tail = _tail_pass(probe['inputs'])
    assert NOT_AVAILABLE in tail['cm']
    assert NOT_AVAILABLE not in tail['cma']


def test_stale_tail_arm_reaches_only_the_zero_lift_moment():
    """BD(63), a Mach late, changes HB(16) and nothing on the curves."""
    inputs = json.loads(json.dumps(_TAIL_PROBE[0]['inputs']))
    a = _tail_pass(inputs)
    inputs['bd63'] = 0.0
    b = _tail_pass(inputs)
    assert a['moment']['cm0'] != b['moment']['cm0']
    for key in ('cl', 'cm', 'cd', 'cla', 'cma'):
        np.testing.assert_array_equal(a[key], b[key])
