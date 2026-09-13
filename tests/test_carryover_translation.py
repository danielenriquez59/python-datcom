"""
Regression tests for the DATCOM Section 4.3.1.2 lift carryover translation.

Source of truth: datcom-legacy/datcom_2000/hbtran.f (Figure 4.3.1.2-10) and
wbclb.f (Figure 4.3.1.2-12A).  Beyond reproducing the source DATA, the
angle-of-attack factors are checked against the slender-body result
``K_W(B) + K_B(W) = (1 + r/s)**2`` (Pitts, Nielsen and Kaattari), which is an
independent constraint the tables were built to satisfy.
"""

import numpy as np
import pytest

from pydatcom.interactions.carryover import (
    body_semispan_ratio, fig4312_10, fig4312_12a,
    calculate_carryover_factors,
)
from pydatcom.aerodynamics.moment import calculate_total_pitching_moment


# Figure 4.3.1.2-10: hbtran.f DATA TFIG10 / DKWB10 / DKBW10
_FIG10_SOURCE = [
    (0.0, 1.0, 0.0), (.1, 1.08, .13), (.2, 1.16, .29), (.3, 1.26, .45),
    (.4, 1.36, .62), (.5, 1.46, .80), (.6, 1.56, 1.0), (.7, 1.67, 1.22),
    (.8, 1.78, 1.45), (.9, 1.89, 1.70), (1.0, 2.0, 2.0),
]

# Figure 4.3.1.2-12A: wbclb.f DATA X12A / Y12A1 / Y12A2
_FIG12A_SOURCE = [
    (0., 1., 0.), (.1, .97, .11), (.2, .95, .21), (.3, .94, .31),
    (.4, .94, .41), (.5, .94, .51), (.6, .94, .60), (.7, .95, .70),
    (.8, .96, .80), (.9, .98, .90), (1., .99, 1.0),
]


@pytest.mark.parametrize("ratio,kwb,kbw", _FIG10_SOURCE)
def test_fig4312_10_source_coordinates(ratio, kwb, kbw):
    """Every Figure 4.3.1.2-10 table entry is reproduced exactly."""
    result = fig4312_10(ratio)
    assert result['kwb'] == pytest.approx(kwb, abs=1e-12)
    assert result['kbw'] == pytest.approx(kbw, abs=1e-12)


@pytest.mark.parametrize("ratio,kkwb,kkbw", _FIG12A_SOURCE)
def test_fig4312_12a_source_coordinates(ratio, kkwb, kkbw):
    """Every Figure 4.3.1.2-12A table entry is reproduced exactly.

    The source writes four of the KKWB entries with the FORTRAN repeat
    count ``4*.94``; this pins the expansion.
    """
    result = fig4312_12a(ratio)
    assert result['kkwb'] == pytest.approx(kkwb, abs=1e-12)
    assert result['kkbw'] == pytest.approx(kkbw, abs=1e-12)


def test_slender_body_sum_identity():
    """K_W(B) + K_B(W) must equal (1 + r/s)**2 to chart accuracy.

    The largest departure over the table is 0.02 at r/s = 0.3, 0.4 and 0.9,
    which is chart-digitization noise rather than a transcription error.
    """
    for ratio, kwb, kbw in _FIG10_SOURCE:
        assert kwb + kbw == pytest.approx((1.0 + ratio)**2, abs=0.021)


def test_slender_body_identity_exact_at_endpoints():
    """The theorem is satisfied exactly where the chart is anchored."""
    for ratio in (0.0, 0.1, 0.6, 0.7, 1.0):
        result = fig4312_10(ratio)
        assert result['kwb'] + result['kbw'] == pytest.approx(
            (1.0 + ratio)**2, abs=1e-12)


def test_zero_body_gives_isolated_surface():
    """With no body the surface carries its own lift and nothing else."""
    assert fig4312_10(0.0) == {'kwb': 1.0, 'kbw': 0.0}
    assert fig4312_12a(0.0) == {'kkwb': 1.0, 'kkbw': 0.0}


def test_factors_are_monotonic():
    """Carryover grows with body size; KWB and KBW never decrease."""
    ratios = np.linspace(0.0, 1.0, 200)
    kwb = [fig4312_10(r)['kwb'] for r in ratios]
    kbw = [fig4312_10(r)['kbw'] for r in ratios]
    assert np.all(np.diff(kwb) >= -1e-12)
    assert np.all(np.diff(kbw) >= -1e-12)


def test_factors_are_clamped_outside_the_chart():
    """INTERX mode 0 clamps; the ratio is physically bounded by [0, 1]."""
    assert fig4312_10(1.5)['kwb'] == pytest.approx(2.0)
    assert fig4312_10(-0.5)['kbw'] == pytest.approx(0.0)


def test_body_semispan_ratio():
    """r/s = (SSPN - SSPNE)/SSPN, the source FACT(1)."""
    state = {'wing_sspn': 15.0, 'wing_sspne': 13.5}
    assert body_semispan_ratio(state) == pytest.approx(0.1)
    state = {'htail_sspn': 6.0, 'htail_sspne': 5.4}
    assert body_semispan_ratio(state, 'htail') == pytest.approx(0.1)


def test_body_semispan_ratio_rejects_bad_spans():
    with pytest.raises(ValueError):
        body_semispan_ratio({'wing_sspn': 0.0})
    with pytest.raises(ValueError):
        body_semispan_ratio({'wing_sspn': 10.0, 'wing_sspne': 12.0})


def test_calculate_carryover_factors_bundles_all_four():
    state = {'htail_sspn': 6.0, 'htail_sspne': 5.4}
    factors = calculate_carryover_factors(state, component='htail')
    assert factors['ratio'] == pytest.approx(0.1)
    assert factors['kwb'] == pytest.approx(1.08)
    assert factors['kbw'] == pytest.approx(0.13)
    assert factors['kkwb'] == pytest.approx(0.97)
    assert factors['kkbw'] == pytest.approx(0.11)


# --------------------------------------------------------------------------
# Integration with the aircraft moment
# --------------------------------------------------------------------------

def _aircraft_state():
    return {
        'wing_type': 1.0, 'wing_chrdr': 6.0, 'wing_chrdtp': 3.0,
        'wing_sspn': 15.0, 'wing_sspne': 13.5,
        'wing_savsi': 0.0, 'wing_chstat': 0.25,
        'htail_type': 1.0, 'htail_chrdr': 3.0, 'htail_chrdtp': 1.5,
        'htail_sspn': 6.0, 'htail_sspne': 5.4,
        'htail_savsi': 0.0, 'htail_chstat': 0.25,
        'synths_xw': 10.0, 'synths_zw': 0.0, 'synths_aliw': 0.0,
        'synths_xh': 35.0, 'synths_zh': 2.0, 'synths_alih': 0.0,
        'synths_xcg': 12.0,
        'options_sref': 135.0, 'options_cbarr': 4.6667,
        'flight_mach': 0.3,
    }


def test_moment_uses_translated_carryover():
    """The aircraft moment now reports the CLWBT carryover path."""
    result = calculate_total_pitching_moment(_aircraft_state(), 0.5, 5.0, 0.3)
    assert result['tail_method'] == 'legacy_clwbt_with_carryover'
    assert result['carryover']['ratio'] == pytest.approx(0.1)
    assert result['carryover']['kwb'] == pytest.approx(1.08)


def test_carryover_increases_tail_effectiveness():
    """A body-mounted tail carries more load than an isolated one."""
    state = _aircraft_state()
    with_body = calculate_total_pitching_moment(state, 0.5, 5.0, 0.3)
    # SSPNE == SSPN is the no-body limit: KWB=1, KBW=0.
    isolated = dict(state, htail_sspne=6.0)
    without = calculate_total_pitching_moment(isolated, 0.5, 5.0, 0.3)
    assert with_body['carryover']['kbw'] > 0.0
    assert without['carryover']['kbw'] == pytest.approx(0.0)
    assert abs(with_body['cm_tail']) > abs(without['cm_tail'])


def test_tail_remains_stabilizing_with_carryover():
    """Carryover must not reverse the sign of the tail's stability effect."""
    state = _aircraft_state()
    no_tail = {k: v for k, v in state.items() if not k.startswith('htail_')}

    def slope(config):
        low = calculate_total_pitching_moment(config, 0.30, 2.0, 0.3)
        high = calculate_total_pitching_moment(config, 0.45, 8.0, 0.3)
        return (high['cm_total'] - low['cm_total']) / 6.0

    assert slope(state) < slope(no_tail)
