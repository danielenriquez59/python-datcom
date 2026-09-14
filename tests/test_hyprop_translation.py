"""
Regression tests for the HYPROP translation.

Source of truth: datcom-legacy/datcom_2000/hyprop.f.

Two source defects are covered here rather than reproduced: a figure whose
call reads past the end of its array, and a branch that assigns the
temperature reading to the pressure variable.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.hyprop import (
    calculate_hyprop, _FIGURES, _LEVELS, _LEVEL_ALTITUDE, _BAND,
    _SHORT_FIGURE, _SHORT_DECLARED, _SHORT_REQUIRED,
)

_PROPERTIES = ('pressure', 'temperature', 'mach', 'density')


def _complete_43a():
    """A stand-in completed Figure 6.3.1-43A, 2 by 11."""
    return list(np.linspace(0.9, 2.0, 22))


# --------------------------------------------------------------------------
# Table inventory
# --------------------------------------------------------------------------

def test_four_properties_each_have_six_altitude_variants():
    assert sorted(_FIGURES) == sorted(_PROPERTIES)
    for prop in _PROPERTIES:
        assert sorted(_FIGURES[prop]) == list('ABCDEF')


def test_every_figure_is_rectangular_except_the_short_one():
    """len(y) == len(x1)*len(x2) for all but Figure 6.3.1-43A."""
    for prop in _PROPERTIES:
        for level in 'ABCDEF':
            x1, x2, y = _FIGURES[prop][level]
            if (prop, level) == _SHORT_FIGURE:
                continue
            assert len(y) == len(x1) * len(x2), f"{prop} {level}"


def test_the_short_figure_is_short_by_the_documented_amount():
    """6.3.1-43A holds 6 values where its call needs 22."""
    x1, x2, y = _FIGURES[_SHORT_FIGURE[0]][_SHORT_FIGURE[1]]
    assert len(y) == _SHORT_DECLARED == 6
    assert _SHORT_REQUIRED == 22
    # Every sibling in that position holds the full 22.
    for prop in _PROPERTIES:
        if prop == _SHORT_FIGURE[0]:
            continue
        assert len(_FIGURES[prop]['A'][2]) == _SHORT_REQUIRED


def test_level_altitudes_run_in_fifty_thousand_foot_steps():
    assert _LEVELS == 'ABCDEF'
    assert _LEVEL_ALTITUDE == [50000., 100000., 150000., 200000., 250000.,
                               300000.]
    assert all(b - a == _BAND
               for a, b in zip(_LEVEL_ALTITUDE, _LEVEL_ALTITUDE[1:]))


# --------------------------------------------------------------------------
# The short-figure guard
# --------------------------------------------------------------------------

@pytest.mark.parametrize("altitude", [0.0, 30000., 50000., 60000., 100000.])
def test_bands_touching_level_a_are_refused(altitude):
    """The source would read past the end of Figure 6.3.1-43A there."""
    with pytest.raises(ValueError, match="6.3.1-43A"):
        calculate_hyprop(altitude, 8.0, 5.0)


def test_a_completed_table_unblocks_those_bands():
    result = calculate_hyprop(30000., 8.0, 5.0,
                              mach_figure_43a=_complete_43a())
    assert all(np.isfinite(result[prop]) for prop in _PROPERTIES)


def test_higher_bands_need_no_replacement():
    assert np.isfinite(calculate_hyprop(120000., 8.0, 5.0)['pressure'])


# --------------------------------------------------------------------------
# Band selection and interpolation
# --------------------------------------------------------------------------

@pytest.mark.parametrize("altitude,band,lower,upper", [
    (120000., 3, 'B', 'C'),
    (175000., 4, 'C', 'D'),
    (240000., 5, 'D', 'E'),
    (290000., 6, 'E', 'F'),
])
def test_band_brackets_follow_the_source(altitude, band, lower, upper):
    result = calculate_hyprop(altitude, 8.0, 5.0)
    assert result['band'] == band
    assert result['level_lower'] == lower
    assert result['level_upper'] == upper


def test_interpolation_fraction_is_inside_the_band():
    for altitude in (110000., 120000., 149000., 175000., 240000.):
        assert 0.0 <= calculate_hyprop(altitude, 8.0, 5.0)['fraction'] <= 1.0


def test_the_top_band_collapses_to_the_last_level():
    """Above 300,000 feet both brackets are F, so the fraction cannot bite."""
    result = calculate_hyprop(400000., 8.0, 5.0)
    assert result['level_lower'] == result['level_upper'] == 'F'
    top = calculate_hyprop(300000., 8.0, 5.0)
    for prop in _PROPERTIES:
        assert result[prop] == pytest.approx(top[prop], rel=1e-9)


def test_properties_rise_with_altitude_through_the_bands():
    values = [calculate_hyprop(alt, 8.0, 5.0)['pressure']
              for alt in (120000., 175000., 240000.)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_all_four_properties_are_returned():
    result = calculate_hyprop(175000., 8.0, 5.0)
    for prop in _PROPERTIES:
        assert prop in result and np.isfinite(result[prop])


def test_each_property_reads_its_own_figure():
    """The source's top branch assigns the temperature figure to P; this
    translation does not, so the four outputs must differ."""
    result = calculate_hyprop(290000., 8.0, 5.0)
    assert result['pressure'] != pytest.approx(result['temperature'])
