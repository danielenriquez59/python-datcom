"""ATMOS checks against its source tables and independent hydrostatic integration.

These exercise the 1962 legacy model, not a modern replacement atmosphere.
"""

import numpy as np
import pytest

from pydatcom.utils.atmosphere import Atmosphere


# Base/upper heights, molecular-scale temperatures, and base pressure copied
# from datcom-legacy/datcom_2000/atmos.f (feet, Rankine, lb/ft**2).
LOW_H = [-16404., 0., 36089., 65617., 104987., 154199., 170604.,
         200131., 250186., 291160.]
HIGH_Z = [295276., 328084., 360892., 393701., 492126., 524934.,
          557743., 623360., 754593., 984252., 1312336., 1640420.,
          1968504., 2296588.]
TEMPERATURE = [577.17, 518.67, 389.97, 389.97, 411.57, 487.17,
               487.17, 454.77, 325.17, 325.17, 379.17, 469.17,
               649.17, 1729.17, 1999.17, 2179.17, 2431.17, 2791.17,
               3295.17, 3889.17, 4357.17, 4663.17, 4861.17]
PRESSURE = [3711.0839, 2116.2165, 472.67563, 114.34314, 18.128355,
            2.3162178, 1.2321972, .38030279, .021671352, .0034313478,
            .00062773411, .00015349091, .000052624212, .000010561806,
            .0000077083076, .0000058267151, .0000035159854, .0000014520255,
            .00000039290563, .000000084030242, .000000022835256,
            .0000000071875452]
MOLECULAR_WEIGHT = [28.9644, 28.88, 28.56, 28.07, 26.92, 26.66,
                    26.4, 25.85, 24.7, 22.66, 19.94, 17.94, 16.84, 16.17]
RADIUS = 20890855.
GMRS = .018743418
GRAVITY = 32.1740484


def geometric(h):
    return RADIUS * h / (RADIUS - h)


@pytest.mark.parametrize("high,index", [(False, i) for i in range(9)] +
                         [(True, i) for i in range(13)])
def test_each_atmos_layer_against_hydrostatic_integral(high, index):
    """Integrate d(log P)/dz directly, independent of ATMOS closed forms."""
    heights = HIGH_Z if high else LOW_H
    j = index + 9 if high else index
    coordinate0, coordinate1 = heights[index:index + 2]
    coordinate = (coordinate0 + coordinate1) / 2
    z0 = coordinate0 if high else geometric(coordinate0)
    z = coordinate if high else geometric(coordinate)
    slope = (TEMPERATURE[j + 1] - TEMPERATURE[j]) / (coordinate1 - coordinate0)

    nodes, weights = np.polynomial.legendre.leggauss(32)
    sample_z = z0 + (nodes + 1) * (z - z0) / 2
    sample_coordinate = sample_z if high else RADIUS * sample_z / (RADIUS + sample_z)
    sample_t = TEMPERATURE[j] + slope * (sample_coordinate - coordinate0)
    integral = (z - z0) / 2 * np.dot(
        weights, GMRS * (RADIUS / (RADIUS + sample_z)) ** 2 / sample_t)
    expected_pressure = PRESSURE[j] * np.exp(-integral)
    tms = TEMPERATURE[j] + slope * (coordinate - coordinate0)
    molecular_weight = ((MOLECULAR_WEIGHT[index] + MOLECULAR_WEIGHT[index + 1]) / 2
                        if high else 28.9644)
    result = Atmosphere.calculate(z)
    assert result['pressure'] == pytest.approx(expected_pressure, rel=1e-11)
    assert result['temperature'] == pytest.approx(tms * molecular_weight / 28.9644)
    assert result['density'] == pytest.approx(GMRS * expected_pressure / GRAVITY / tms)
    assert result['cs'] == pytest.approx(49.022164 * np.sqrt(tms))

    # The sound/density outputs are logarithmic derivatives, unlike dp_dz/dt_dz.
    lower, upper = Atmosphere.calculate(z - .1), Atmosphere.calculate(z + .1)
    for name, derivative in [('pressure', 'dp_dz'), ('temperature', 'dt_dz'),
                             ('density', 'drho_dz'), ('cs', 'dcs_dz')]:
        numerical = (upper[name] - lower[name]) / .2
        if name in ('density', 'cs'):
            numerical /= result[name]
        assert result[derivative] == pytest.approx(numerical, rel=1e-6, abs=1e-12)


def test_rounded_geopotential_endpoint_uses_final_low_layer():
    # H(295276) = 291160.6702 ft, slightly beyond the last tabulated HG.
    z = 295276.
    h = RADIUS * z / (RADIUS + z)
    result = Atmosphere.calculate(z)
    assert result['temperature'] == pytest.approx(325.17)
    expected = .021671352 * np.exp(GMRS * (250186. - h) / 325.17)
    assert result['pressure'] == pytest.approx(expected)
    assert all(np.isfinite(value) for value in result.values())


def test_above_table_retains_last_high_layer_as_in_fortran():
    # This verifies legacy extrapolation, not physical accuracy above 700 km.
    lower = Atmosphere.calculate(2296588.)
    upper = Atmosphere.calculate(2296589.)
    assert upper['pressure'] < lower['pressure']
    assert upper['pressure'] == pytest.approx(lower['pressure'], rel=1e-5)
    tms = 4861.17 + 198. / 328084.
    em = 16.17 - .67 / 328084.
    assert upper['temperature'] == pytest.approx(tms * em / 28.9644)


@pytest.mark.parametrize('z', [geometric(h) for h in LOW_H] + HIGH_Z)
def test_tabulated_boundaries_are_finite(z):
    assert all(np.isfinite(value) for value in Atmosphere.calculate(z).values())
