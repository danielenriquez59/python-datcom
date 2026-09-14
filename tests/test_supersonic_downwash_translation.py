"""
Regression tests for the Figure 4.4.1-76 supersonic downwash translation.

Source of truth: datcom-legacy/datcom_2000/sdwa.f, sdwb.f, sdwc.f,
sdwd.f, sdwe.f.

Each routine is a chain of INTERX lookups, which have their own translation
tests, so these cover the table packing each routine owns, the stage
chaining, and the corner values of the first stage where the source DATA can
be read off directly.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.supersonic_downwash import (
    calculate_sdwa, calculate_sdwb,
    _XA1, _YA1, _XA2, _YA2, _XA3, _YA3,
    _XB1, _YB1, _XB2, _YB2, _XB3, _YB3,
)


# --------------------------------------------------------------------------
# Table shapes: the LIND and LDEP of the source INTERX calls
# --------------------------------------------------------------------------

def test_figure_76a_table_shapes():
    assert len(_XA1) == 10 and len(_YA1) == 24     # LGH (6,4), LIND 6
    assert len(_XA2) == 6 and len(_YA2) == 6       # LGH (2,3), LIND 3
    assert len(_XA3) == 10 and len(_YA3) == 10     # LGH (2,5), LIND 5


def test_figure_76b_table_shapes():
    assert len(_XB1) == 12 and len(_YB1) == 32     # LGH (8,4), LIND 8
    assert len(_XB2) == 8 and len(_YB2) == 16      # LGH (4,4), LIND 4
    assert len(_XB3) == 6 and len(_YB3) == 6       # LGH (2,3), LIND 3


def test_lind_stride_leaves_padding_between_grids():
    """XA2 packs a 2-entry grid into a stride of 3, so index 2 is padding.

    This is the behaviour that makes the LIND argument necessary: the second
    grid starts at offset LIND, not at the end of the first grid.
    """
    assert _XA2[:2] == [0.0, 3.0]
    assert _XA2[2] == 0.0            # unused padding
    assert _XA2[3:6] == [4.0, 8.0, 12.0]


# --------------------------------------------------------------------------
# First-stage corners, readable straight from the source DATA
# --------------------------------------------------------------------------

@pytest.mark.parametrize("x,z,expected", [
    (0.90, 0.0, 0.90), (2.4, 0.0, 2.34),
    (0.90, 0.5, 0.00), (2.4, 0.5, 1.00),
    (1.4, 0.3, 0.87),
])
def test_figure_76a1_corners(x, z, expected):
    assert calculate_sdwa(x, 0.0, z, 4.0)['stage1'] == pytest.approx(
        expected, abs=1e-12)


@pytest.mark.parametrize("x,z,expected", [
    (1.0, 0.0, 0.32), (2.4, 0.0, 1.77),
    (1.0, 0.5, 0.00), (2.4, 0.5, 1.20),
    (1.5, 0.3, 1.00),
])
def test_figure_76b1_corners(x, z, expected):
    assert calculate_sdwb(x, 0.0, z, 4.0)['stage1'] == pytest.approx(
        expected, abs=1e-12)


# --------------------------------------------------------------------------
# Chaining and behaviour
# --------------------------------------------------------------------------

def test_each_stage_feeds_the_next():
    """The reported intermediates must be the ones the chain actually used."""
    from pydatcom.utils.legacy_interp import interx
    result = calculate_sdwa(1.4, 0.3, 0.1, 8.0)
    stage1 = interx(2, _XA1, [1.4, 0.1], [6, 4], _YA1, lind=6,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XA2, [stage1, 8.0], [2, 3], _YA2, lind=3,
                    lx1l=0, lx2l=1, lx1u=1, lx2u=1)
    expected = interx(2, _XA3, [stage2, 0.3], [2, 5], _YA3, lind=5,
                      lx1l=0, lx2l=0, lx1u=1, lx2u=2)
    assert result['stage1'] == pytest.approx(stage1)
    assert result['stage2'] == pytest.approx(stage2)
    assert result['sdw'] == pytest.approx(expected)


def test_downwash_grows_with_streamwise_distance():
    for routine in (calculate_sdwa, calculate_sdwb):
        values = [routine(x, 0.2, 0.1, 8.0)['sdw'] for x in (1.0, 1.4, 2.0, 2.4)]
        assert all(b > a for a, b in zip(values, values[1:]))


def test_downwash_falls_as_the_station_rises_out_of_the_sheet():
    """The Z axis of the first-stage figure runs 0.0 to 0.5 downward."""
    for routine in (calculate_sdwa, calculate_sdwb):
        values = [routine(1.6, 0.2, z, 8.0)['sdw'] for z in (0.0, 0.1, 0.3, 0.5)]
        assert all(b <= a for a, b in zip(values, values[1:]))


def test_routines_report_their_taper():
    assert calculate_sdwa(1.5, 0.2, 0.1, 8.0)['taper'] == 0.00
    assert calculate_sdwb(1.5, 0.2, 0.1, 8.0)['taper'] == 1.00


def test_results_are_finite_across_the_table_ranges():
    for routine in (calculate_sdwa, calculate_sdwb):
        for x in np.linspace(0.9, 2.4, 7):
            for z in (0.0, 0.25, 0.5):
                for ab in (4.0, 8.0, 12.0):
                    value = routine(float(x), 0.2, z, ab)['sdw']
                    assert np.isfinite(value)


# --------------------------------------------------------------------------
# Figures 4.4.1-76C, -76D and -76E
# --------------------------------------------------------------------------

from pydatcom.aerodynamics.supersonic_downwash import (
    calculate_sdwc, calculate_sdwd, calculate_sdwe,
    _XC1, _YC1, _XC2, _YC2, _XC3, _YC3A, _YC3B, _XC4, _YC4,
    _XD1, _YD1, _XD2, _YD2, _XD3, _YD3A, _YD3B, _YD3C, _XD4, _YD4,
    _XE1, _YE1, _XE2, _YE2, _XE3, _YE3, _XE4, _YE4,
)


def test_figure_76c_table_shapes():
    assert len(_XC1) == 11 and len(_YC1) == 28     # LGH (7,4), LIND 7
    assert len(_XC2) == 12 and len(_YC2) == 12     # LGH (2,6), LIND 6
    assert len(_XC3) == 5 and len(_YC3A) == 5 and len(_YC3B) == 5
    assert len(_XC4) == 4 and len(_YC4) == 4       # LGH (2,2), LIND 2


def test_figure_76d_table_shapes():
    assert len(_XD1) == 9 and len(_YD1) == 20      # LGH (5,4), LIND 5
    assert len(_XD2) == 9 and len(_YD2) == 20      # LGH (5,4), LIND 5
    assert len(_XD3) == 2
    assert len(_XD4) == 8 and len(_YD4) == 15      # LGH (5,3), LIND 5


def test_figure_76e_table_shapes():
    assert len(_XE1) == 9 and len(_YE1) == 20      # LGH (5,4), LIND 5
    assert len(_XE2) == 7 and len(_YE2) == 10      # LGH (5,2), LIND 5
    assert len(_XE3) == 4 and len(_YE3) == 4       # LGH (2,2), LIND 2
    assert len(_XE4) == 8 and len(_YE4) == 15      # LGH (5,3), LIND 5


def test_figure_76c2_padding_between_packed_grids():
    """The source writes the stride padding as '4*0.' between the grids."""
    assert _XC2[:2] == [0.00, 2.50]
    assert _XC2[2:6] == [0.0, 0.0, 0.0, 0.0]
    assert _XC2[6:] == [2.70, 3.20, 5.40, 6.40, 10.7, 12.8]


def test_figure_76e_repeats_76d_blocks_verbatim():
    """sdwe.f re-declares XD1/YD1 and XD4/YD4 rather than sharing them."""
    assert _XE1 == _XD1 and _YE1 == _YD1
    assert _XE4 == _XD4 and _YE4 == _YD4


@pytest.mark.parametrize("x,z,expected", [
    (0.900, 0.0, 0.80), (2.40, 0.0, 2.10),
    (0.900, 0.5, 0.00), (2.40, 0.5, 1.00),
])
def test_figure_76c1_corners(x, z, expected):
    assert calculate_sdwc(x, 0.0, z, 4.0)['stage1'] == pytest.approx(
        expected, abs=1e-12)


@pytest.mark.parametrize("x,z,expected", [
    (1.0, 0.0, 1.68), (2.4, 0.0, 3.62),
    (1.0, 0.5, 0.00), (2.4, 0.5, 1.51),
])
def test_figure_76d1_corners(x, z, expected):
    assert calculate_sdwd(x, 0.0, z, 4.0)['stage1'] == pytest.approx(
        expected, abs=1e-12)


def test_sdwc_returns_two_taper_values():
    result = calculate_sdwc(1.8, 0.15, 0.1, 6.0)
    assert result['sdw'].shape == (2,)
    assert result['taper'] == (0.25, 0.50)


def test_sdwd_returns_three_taper_values():
    result = calculate_sdwd(1.8, 0.15, 0.1, 6.0)
    assert result['sdw'].shape == (3,)
    assert result['taper'] == (0.00, 0.25, 0.50)


def test_sdwc_and_sdwd_fall_with_taper():
    """More taper puts the aft surface in weaker downwash."""
    c = calculate_sdwc(2.0, 0.2, 0.1, 8.0)['sdw']
    d = calculate_sdwd(2.0, 0.2, 0.1, 8.0)['sdw']
    assert c[0] > c[1]
    assert d[0] > d[1] > d[2]


def test_sdwc_and_sdwd_share_their_first_two_stages():
    """Both split only at the third stage, so the shared stages are scalars."""
    for routine in (calculate_sdwc, calculate_sdwd):
        result = routine(1.8, 0.15, 0.1, 6.0)
        assert isinstance(result['stage1'], float)
        assert isinstance(result['stage2'], float)


def test_sdwe_first_stage_matches_sdwd():
    """The duplicated DATA blocks must give identical first stages."""
    assert calculate_sdwe(1.8, 0.0, 0.1, 4.0, 0.3)['stage1'] == pytest.approx(
        calculate_sdwd(1.8, 0.0, 0.1, 4.0)['stage1'])


def test_sdwe_takes_taper_as_an_input():
    """Taper enters SDWE at stage 3 rather than selecting between curves."""
    values = [calculate_sdwe(2.0, 0.2, 0.1, 4.5, t)['sdw']
              for t in (0.25, 0.375, 0.50)]
    assert all(np.isfinite(v) for v in values)
    assert values[0] != values[2]


def test_sdwcde_downwash_grows_over_the_rising_range():
    """Downwash grows with streamwise distance up to the table turnover.

    Figure 4.4.1-76C1 is not monotonic to its last station: the z = 0.1
    curve runs ... 1.70, 1.81, 1.80, turning over between x = 2.20 and
    x = 2.40. That is chart data, not a translation artefact, so the growth
    assertion stops short of the turnover and the turnover itself is pinned
    separately below.
    """
    for routine, extra in ((calculate_sdwc, ()), (calculate_sdwd, ()),
                           (calculate_sdwe, (0.375,))):
        values = []
        for x in (1.0, 1.6, 2.2):
            out = routine(x, 0.2, 0.1, 6.0, *extra)['sdw']
            values.append(float(np.atleast_1d(out)[0]))
        assert all(b >= a for a, b in zip(values, values[1:]))


def test_figure_76c1_turns_over_at_its_last_station():
    """The z = 0.1 curve of 76C1 peaks at x = 2.20 rather than x = 2.40."""
    assert _YC1[12] == 1.81 and _YC1[13] == 1.80
    peak = calculate_sdwc(2.20, 0.0, 0.1, 4.0)['stage1']
    last = calculate_sdwc(2.40, 0.0, 0.1, 4.0)['stage1']
    assert last < peak
