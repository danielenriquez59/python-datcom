"""
Regression tests for the CMALPH translation.

Source of truth: datcom-legacy/datcom_2000/cmalph.f.

Checked against a compiled probe of CMALPH with FWDXAC linked in
(tools/probes/cmalph.py): the moment curve, CMa, CM0, the whole C work
array, and the COMMON values it overwrites.  Every table is checked against
a re-parse of the FORTRAN.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import cmalph as module
from pydatcom.aerodynamics.cmalph import calculate_cmalph
from pydatcom.aerodynamics.moment import calculate_cmalph_zero_lift_moment
from pydatcom.geometry.wing import calculate_straight_exposed_geometry

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'cmalph.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    probe = _PROBE[case]
    out = probe['outputs']
    result = calculate_cmalph(**probe['inputs'])
    np.testing.assert_allclose(result['cm'], out['CM'], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(
        [result['cma'], result['cm0'], result['a170'], result['a38'],
         result['a62']], out['OUT'], rtol=1e-9, atol=1e-14)
    for index, value in result['c'].items():
        assert value == pytest.approx(out['C'][index - 1], rel=1e-9,
                                      abs=1e-12), f'C({index})'


def test_probe_reaches_every_branch():
    results = [calculate_cmalph(**p['inputs']) for p in _PROBE]
    inputs = [p['inputs'] for p in _PROBE]
    assert {r['nonlinear'] for r in results} == {True, False}
    assert any('fwdxac' in r for r in results)
    assert any(r.get('source_defect') for r in results)
    assert any(4 in r['c'] for r in results)                 # twisted
    assert any(i['section']['xac'] != module.UNUSED for i in inputs)
    assert any(r['a38'] == 0.0001 for r in results)          # zero tangent
    straight_low = [r for r, i in zip(results, inputs)
                    if i['planform_type'] == 1.0 and 10 in r['c']
                    and i['geometry']['a7'] <= i['geometry']['a125']
                    and i['geometry']['a38'] > 0]
    assert {r['c'][10] >= 1.0 for r in straight_low} == {True, False}
    kinds = {i['planform_type'] for i in inputs}
    assert {1.0, 2.0, 3.0, 4.0} <= kinds
    # The nonlinear path reaches both angle regimes and TANRAT.
    nonlinear = [(r, i) for r, i in zip(results, inputs) if r['nonlinear']]
    assert any(max(i['alpha_deg']) > math.degrees(r['c'][36])
               for r, i in nonlinear)
    assert any(46 in r['c'] for r, _ in nonlinear)


def test_every_table_value_matches_the_source():
    source = parse('cmalph')
    names = ['XCMOM', 'YCMOM', 'X31412', 'X11412', 'X21412', 'X322A',
             'X122A', 'X222A', 'Y22A', 'X322B', 'X122B', 'X222B', 'Y22B',
             'X211A', 'Y11A', 'X211B', 'Y11B', 'X111C', 'X211C', 'Y11C',
             'X112A', 'X212A', 'Y12A', 'X112B1', 'X212B1', 'Y12B1', 'X112B2',
             'X212B2', 'Y12B2', 'X213A', 'X113A', 'Y13A', 'X213B1', 'X113B1',
             'Y13B1', 'X213B2', 'X113B2', 'Y13B2']
    for name in names:
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)
    # Y41412 is three EQUIVALENCEd DATA arrays end to end.
    np.testing.assert_array_equal(
        module._Y41412, source['Y415A'] + source['Y415B'] + source['Y415C'])


def test_zero_lift_moment_agrees_with_the_earlier_partial_translation():
    """moment.py's constant-section CMO and CMALPH's C(5) must agree."""
    state = {'wing_type': 1.0, 'wing_chrdr': 8.0, 'wing_chrdtp': 3.0,
             'wing_sspn': 14.0, 'wing_sspne': 12.5, 'wing_savsi': 25.0,
             'wing_chstat': 0.25, 'wing_cmo': -0.03, 'wing_cmot': -0.02,
             'options_sref': 150.0, 'options_cbarr': 6.0}
    wing = calculate_straight_exposed_geometry(state)
    cos_c4 = 1.0 / math.sqrt(1.0 + wing['tan_c4']**2)
    geometry = {'a3': wing['area'], 'a7': wing['aspect_ratio'],
                'a10': wing['root_chord'], 'a16': wing['mac'],
                'a27': wing['taper_ratio'], 'a30': 3.0, 'a34': 30.0,
                'a37': 0.87, 'a38': wing['tan_le'], 'a43': cos_c4,
                'a124': 1.5, 'a125': 3.0, 'a173': 5.0}
    for mach in (0.2, 0.5, 0.7):
        result = calculate_cmalph(
            1.0, [0.0, 4.0], geometry,
            {'cmo': -0.03, 'cmot': -0.02, 'twista': 0.0, 'deltay': 1.5},
            {'cla': 0.07, 'cn': [0.0, 0.28], 'alpha_clmax': 14.0},
            {'mach': mach, 'beta': math.sqrt(1 - mach**2)}, 150.0, 6.0)
        assert result['cm0'] == pytest.approx(
            calculate_cmalph_zero_lift_moment(state, mach), rel=1e-12)


def test_linear_moment_is_cn_times_the_static_margin():
    """Above 6/A(124), CM = CN*(x_cg - x_ac)*c_r/cbar + CM0 exactly."""
    probe = next(p for p in _PROBE
                 if not calculate_cmalph(**p['inputs'])['nonlinear'])
    r = calculate_cmalph(**probe['inputs'])
    g = probe['inputs']['geometry']
    lever = (g['a173'] / g['a10'] - r['xac']) * g['a10'] / probe['inputs']['cbarr']
    np.testing.assert_allclose(
        r['cm'], np.array(probe['inputs']['lift']['cn']) * lever + r['cm0'],
        rtol=1e-12)
    assert r['cma'] == pytest.approx(lever * probe['inputs']['lift']['cla'])


def test_nonlinear_moment_is_continuous_at_the_reference_angle():
    """Below the reference angle the centre of pressure comes from the
    Figure 4.1.4.3-23/-24 increments; above it a line to the centroid
    starts from C(49), the same point.  The centre of pressure is
    continuous there."""
    probe = next(p for p in _PROBE if p['inputs']['planform_type'] == 1.0
                 and calculate_cmalph(**p['inputs'])['nonlinear'])
    inputs = json.loads(json.dumps(probe['inputs']))
    r = calculate_cmalph(**inputs)
    ref = math.degrees(r['c'][36])
    inputs['alpha_deg'] = [inputs['alpha_deg'][0], ref - 1e-7, ref + 1e-7]
    inputs['lift']['cn'] = [0.1, 1.0, 1.0]
    r = calculate_cmalph(**inputs)
    assert r['cm'][1] == pytest.approx(r['cm'][2], abs=1e-5)


# --------------------------------------------------------------------------
# Overlays M31O37 (wing) and M33O41 (tail)
# --------------------------------------------------------------------------

from pydatcom.aerodynamics.cmalph import (  # noqa: E402
    NOT_AVAILABLE, calculate_moment_overlay,
)

_OVERLAY_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'm31o37.json').read_text())
_CMALPH_KEYS = ('planform_type', 'alpha_deg', 'geometry', 'section', 'lift',
                'flight', 'sref', 'cbarr')


def _overlay(inputs):
    return calculate_moment_overlay(
        inputs['surface'], inputs['free_alpha_deg'], inputs['cd'],
        inputs['cl'], {k: inputs[k] for k in _CMALPH_KEYS},
        inputs['experimental'])


@pytest.mark.parametrize("case", range(len(_OVERLAY_PROBE)))
def test_moment_overlay_matches_compiled_routine(case):
    probe = _OVERLAY_PROBE[case]
    result = _overlay(probe['inputs'])
    for tag, key in [('CM', 'cm'), ('CN', 'cn'), ('CA', 'ca'),
                     ('CLA', 'cla'), ('CMA', 'cma')]:
        np.testing.assert_allclose(result[key], probe['outputs'][tag],
                                   rtol=1e-9, atol=1e-12, err_msg=tag)


def test_moment_overlay_probe_reaches_every_branch():
    seen = set()
    for p in _OVERLAY_PROBE:
        r = _overlay(p['inputs'])
        seen.add((p['inputs']['surface'], p['inputs']['experimental'],
                  NOT_AVAILABLE in r['cm'], r['cmalph']['nonlinear']))
    for surface in ('wing', 'tail'):
        assert (surface, False, True, False) in seen      # cut off
        assert (surface, True, False, False) in seen      # experimental
        assert (surface, False, False, True) in seen      # test skipped


def test_tail_keeps_its_first_slope_when_the_test_is_skipped():
    """M33O41's restore sits inside the skipped branch; M31O37's does not."""
    probe = next(p for p in _OVERLAY_PROBE
                 if p['inputs']['surface'] == 'tail'
                 and not p['inputs']['experimental']
                 and _overlay(p['inputs'])['cmalph']['nonlinear'])
    # A lift curve with curvature, so its TBFUNX slope is not CLa.
    inputs = dict(probe['inputs'])
    local = np.array(inputs['alpha_deg'])
    inputs['cl'] = list(inputs['lift']['cla'] * local * (1 + 0.01 * local))
    tail = _overlay(inputs)
    wing = _overlay(dict(inputs, surface='wing'))
    free = np.array(inputs['free_alpha_deg'])
    first_slope = np.polyfit(free[:3], inputs['cl'][:3], 2)
    assert wing['cla'][0] == inputs['lift']['cla']
    assert tail['cla'][0] != pytest.approx(inputs['lift']['cla'])
    assert tail['cla'][0] == pytest.approx(
        np.polyval(np.polyder(first_slope), free[0]))


def test_tail_limit_is_tighter_than_the_wing_limit():
    """7.5 percent against 90: the tail cuts the moment off earlier."""
    wing = next(_overlay(p['inputs']) for p in _OVERLAY_PROBE
                if p['inputs']['surface'] == 'wing'
                and NOT_AVAILABLE in _overlay(p['inputs'])['cm'])
    tail = next(_overlay(p['inputs']) for p in _OVERLAY_PROBE
                if p['inputs']['surface'] == 'tail'
                and NOT_AVAILABLE in _overlay(p['inputs'])['cm'])
    assert (np.sum(tail['cm'] == NOT_AVAILABLE) >=
            np.sum(wing['cm'] == NOT_AVAILABLE))


# --- CMALPO (tools/probes/cmalpo.py)

from pydatcom.aerodynamics.cmalph import calculate_cmalpo  # noqa: E402

_CMALPO = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                      'cmalpo.json').read_text())


def _cmalpo(c):
    return calculate_cmalpo(c['planform_type'], c['geometry'], c['section'],
                            c['mach'], c['first_mach'], c['cbarr'])


@pytest.mark.parametrize("case", range(len(_CMALPO)))
def test_cmalpo_matches_compiled_routine(case):
    c, o = _CMALPO[case]['inputs'], _CMALPO[case]['outputs']
    r = _cmalpo(c)
    dcmdcl, a62 = o['OUT']
    assert r['dcmdcl'] == pytest.approx(dcmdcl, rel=1e-12, abs=1e-15)
    assert r['a62'] == a62


def test_cmalpo_panel_weights_follow_the_first_mach_number():
    """The Mach-zero routine weights its panels with lift slopes formed at
    FLC(3); B(1) only reaches FWDXAC."""
    c = next(p['inputs'] for p in _CMALPO
             if p['inputs']['planform_type'] == 3.0 and
             p['inputs']['geometry']['a62'] > 0)
    base = _cmalpo(c)
    moved = _cmalpo(dict(c, first_mach=c['first_mach'] + 0.3))
    assert moved['c'][171] != base['c'][171]
    assert moved['dcmdcl'] != base['dcmdcl']
    assert _cmalpo(dict(c, mach=c['mach'] + 0.3))['dcmdcl'] == \
        base['dcmdcl']


def test_cmalpo_high_aspect_ratio_centre_is_its_own():
    """Above A(125) CMALPO uses (A(161)-(SSPN-SSPNE)*A(62))/A(10)."""
    c = _CMALPO[0]['inputs']
    g, s = c['geometry'], c['section']
    assert g['a7'] > g['a125']
    assert _cmalpo(c)['xac'] == pytest.approx(
        (g['a161'] - (s['sspn'] - s['sspne']) * g['a62']) / g['a10'])
