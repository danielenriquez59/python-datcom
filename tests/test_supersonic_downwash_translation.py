"""
Regression tests for the Figure 4.4.1-76 supersonic downwash translation.

Source of truth: datcom-legacy/datcom_2000/sdwa.f, sdwb.f.

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
