"""FIG68 source contracts and independent theta-beta-M verification."""

import numpy as np
import pytest

from pydatcom.utils.table_lookup import fig68


def _deflection(mach, shock_angle):
    """Independent perfect-gas oblique-shock equation, angles in radians."""
    return np.arctan(
        2.0 / np.tan(shock_angle)
        * (mach**2 * np.sin(shock_angle)**2 - 1.0)
        / (mach**2 * (1.4 + np.cos(2.0 * shock_angle)) + 2.0)
    )


@pytest.mark.parametrize("mach", [1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0])
def test_fig68_satisfies_oblique_shock_equation_on_weak_branch(mach):
    # Locate the maximum independently, then bracket the weak branch.
    lower, upper = np.arcsin(1.0/mach), np.pi/2.0
    for _ in range(100):
        left = lower + (upper-lower)/3.0
        right = upper - (upper-lower)/3.0
        if _deflection(mach, left) < _deflection(mach, right):
            lower = left
        else:
            upper = right
    peak = (lower+upper)/2.0
    maximum = _deflection(mach, peak)

    for fraction in [0.01, 0.25, 0.5, 0.9, 0.999]:
        delta = maximum * fraction
        angle, error = fig68(mach, np.rad2deg(delta))
        assert error == 0
        assert np.arcsin(1.0/mach) <= np.deg2rad(angle) <= peak
        assert _deflection(mach, np.deg2rad(angle)) == pytest.approx(delta, abs=1e-8)

    # Detached output is the maximum wedge angle, not the shock angle.
    angle, error = fig68(mach, np.rad2deg(maximum) + 0.01)
    assert error == 2
    assert angle == pytest.approx(np.rad2deg(maximum), abs=2e-8)
    # The exact threshold takes the legacy label 1040 branch.
    peak_angle, error = fig68(mach, angle)
    assert error == 0
    assert peak_angle == pytest.approx(np.rad2deg(peak), abs=2e-6)


@pytest.mark.parametrize(
    "mach,delta,angle,error",
    [(0.8, 10., 0., 3), (0.8, -10., 0., 3),
     (1., -1., 90., 1), (1., 0., 90., 0), (1., 1., 0., 2),
     (2., -1., 30., 1), (2., 0., 30., 0),
     (2., 10., 39.313931844818875, 0),
     (3., 20., 37.763634148375765, 0)],
)
def test_fig68_branch_contract_and_reference_angles(mach, delta, angle, error):
    actual_angle, actual_error = fig68(mach, delta)
    assert actual_error == error
    assert actual_angle == pytest.approx(angle, abs=1e-7)
