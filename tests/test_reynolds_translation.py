"""Source-contract checks for FLTCON RNNUB dimensionalization."""

import numpy as np
import pytest

from pydatcom.aerodynamics.calculator import AerodynamicCalculator
from pydatcom.utils.atmosphere import Atmosphere


def wing_state(**updates):
    state = {
        'wing_aspect_ratio': 6., 'wing_area': 80., 'wing_mac': 4.,
        'options_sref': 80., 'options_cbarr': 5.,
        'flight_mach': [.6, .9, 1.4],
    }
    state.update(updates)
    return state


def test_rnnub_is_per_unit_length_and_tracks_mach_index():
    calc = AerodynamicCalculator(wing_state(flight_rnnub=[2e6, 3e6, 5e6]))
    assert calc._estimate_reynolds(.6) == 8e6
    assert calc._estimate_reynolds(.9) == 12e6
    assert calc._estimate_reynolds(1.4) == 20e6
    # Explicit loop index is accepted only when it identifies this Mach.
    assert calc._estimate_reynolds(.9, 1) == 12e6
    assert calc._estimate_reynolds(.9, 2) == 12e6


def test_scalar_rnnub_and_body_length():
    wing = AerodynamicCalculator(wing_state(flight_rnnub=2e6))
    assert wing._estimate_reynolds(.6) == 8e6
    body = AerodynamicCalculator({
        'body_nx': 2, 'body_x': [3., 13.], 'body_s': [0., 1.],
        'body_r': [0., .5], 'flight_rnnub': 2e6,
    })
    assert body._estimate_reynolds(.6) == 20e6


def test_source_atmospheric_rnnub_equation_and_loop_one_indexing():
    altitudes = [0., 10000., 30000.]
    calc = AerodynamicCalculator(wing_state(flight_alt=altitudes, flight_loop=1))
    for index, (mach, altitude) in enumerate(zip([.6, .9, 1.4], altitudes)):
        atmosphere = Atmosphere.calculate(altitude)
        expected_per_length = (1.2527e6 * atmosphere['pressure'] * mach *
                               (atmosphere['temperature'] + 198.6) /
                               atmosphere['temperature']**2)
        assert calc._estimate_reynolds(mach, index) == pytest.approx(4 * expected_per_length)


def test_source_fallback_and_explicit_dimensionless_reynolds():
    state = wing_state()
    calc = AerodynamicCalculator(state)
    assert calc._estimate_reynolds(.6) == 20e6
    result = calc.calculate_at_condition(3., .6, reynolds=7.5e6)
    assert result['reynolds'] == 7.5e6


def test_mach_sweep_reports_each_dimensionless_reynolds_number():
    calc = AerodynamicCalculator(wing_state(flight_rnnub=[2e6, 3e6, 5e6]))
    result = calc.calculate_mach_sweep(3., np.array([.6, .9, 1.4]))
    np.testing.assert_array_equal(result['reynolds'], [8e6, 12e6, 20e6])


def test_invalid_characteristic_length_or_rnnub_is_explicit():
    with pytest.raises(ValueError, match='characteristic length'):
        AerodynamicCalculator({'wing_area': 1., 'flight_rnnub': [2e6]})._estimate_reynolds(.6)
    with pytest.raises(ValueError, match='RNNUB'):
        AerodynamicCalculator(wing_state(flight_rnnub=[-2e6]))._estimate_reynolds(.6)
