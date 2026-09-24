"""
IDEAL: section ideal aerodynamic parameters by Weber thin-airfoil analysis.

Builds the five Weber coefficient sets that ``SLOPE`` consumes, on a fixed
32-station cosine-spaced grid:

- ``ST1NU`` source distribution in parallel flow
- ``ST2NU`` slope of the thickness distribution
- ``ST3NU`` vortex distribution in normal flow
- ``ST4NU`` vortex distribution due to camber
- ``ST5NU`` slope of the mean line

It then derives the ideal angle of attack, the zero-lift angle, the ideal
lift coefficient and the quarter-chord pitching moment.

Reference: datcom-legacy/datcom_2000/ideal.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.interpolation import asmint
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The source fixes the Weber grid at 32 stations.
_STATIONS = 32


def _weber_grid(count: int = _STATIONS):
    """The cosine-spaced station grid, one-based like the source.

    ``THNU`` runs to pi, so the last station sits at the leading edge where
    ``x = 0`` and ``sin(THNU) = 0``.
    """
    nu = np.arange(1, count + 1, dtype=float)
    mu = np.arange(1, count, dtype=float)
    theta_nu = nu * PI / count
    theta_mu = mu * PI / count
    x = 0.5 * (np.cos(theta_nu) + 1.0)
    return theta_nu, theta_mu, x


def calculate_ideal(x_section: Sequence[float],
                    thickness: Sequence[float],
                    camber: Sequence[float],
                    leading_edge_radius: float,
                    supersonic_section: bool = False,
                    cli: Optional[float] = None,
                    alpha_zero_lift: Optional[float] = None
                    ) -> Dict[str, object]:
    """Translate IDEAL: Weber coefficients and section ideal parameters.

    Args:
        x_section: Section chord stations.
        thickness: Half-thickness ordinate at each station, ``THN``.
        camber: Mean-line ordinate at each station, ``CAM``.
        leading_edge_radius: ``RHO``; the source forms ``A0 = sqrt(2*RHO)``.
        supersonic_section: When set, thickness is resampled with TBFUNX
            rather than the ASMINT spline, matching the source's test for a
            supersonic airfoil designation.
        cli: A user-supplied ideal lift coefficient.  When given, the source
            skips deriving ``AI``, ``ALO`` and ``CLI`` entirely, so
            ``alpha_zero_lift`` must be supplied for the moment.
        alpha_zero_lift: ``ALO``, required only alongside ``cli``.

    Returns:
        Dictionary with the five Weber coefficient arrays, the resampled
        thickness and camber, and the derived section parameters.

    Raises:
        ValueError: If the section arrays are mismatched, or ``cli`` is
            supplied without ``alpha_zero_lift``.

    Notes:
        The last station sits at the leading edge, where ``sin(THNU)`` is
        zero.  ``ST2NU``, ``ST4NU`` and ``ST5NU`` all divide by it, so
        their final entries are singular.  The source is single precision
        and lands on a large finite number there; in double precision the
        value is far larger still.  None of the three is consumed at that
        station downstream: ``SLOPE`` reads them only over stations 1 to
        N-1 and takes ``ST1NU`` and ``ST3NU`` at the leading edge instead.
        All three final entries are returned as NaN so the singularity
        cannot be used by accident.
    """
    x_section = np.asarray(x_section, dtype=float)
    thickness = np.asarray(thickness, dtype=float)
    camber = np.asarray(camber, dtype=float)
    if x_section.ndim != 1 or x_section.shape != thickness.shape or \
            x_section.shape != camber.shape or x_section.size < 2:
        raise ValueError("IDEAL needs matching section coordinate arrays")
    if cli is not None and alpha_zero_lift is None:
        raise ValueError(
            "a supplied CLI makes the source skip deriving ALO, so the "
            "zero-lift angle must be supplied with it")

    count = _STATIONS
    last = count - 1                       # the source's L = N-1
    theta_nu, theta_mu, x = _weber_grid(count)
    a0 = np.sqrt(2.0 * leading_edge_radius)

    # Resample thickness and camber onto the Weber stations.
    if supersonic_section:
        zt = np.array([tbfunx(x_section, thickness, station,
                              lower=0, upper=0)[0] for station in x])
    else:
        zt = np.asarray(asmint(x_section, thickness, x), dtype=float)
    zc = np.asarray(asmint(x_section, camber, x), dtype=float)

    sin_nu = np.sin(theta_nu)
    cos_nu = np.cos(theta_nu)
    sin_mu = np.sin(theta_mu)
    cos_mu = np.cos(theta_mu)

    st1 = np.zeros(count)
    st2 = np.zeros(count)
    st3 = np.zeros(count)
    st4 = np.zeros(count)
    st5 = np.zeros(count)

    with np.errstate(divide='ignore', invalid='ignore'):
        for i in range(count):                     # zero-based; source I-1
            for mu_index in range(last):           # zero-based; source J-1
                sign = (-1.0)**((mu_index + 1) - (i + 1))
                same = (i == mu_index)
                delta = cos_mu[mu_index] - cos_nu[i]

                if same:
                    s1 = count / sin_nu[i]
                    s2 = cos_nu[i] / sin_nu[i]**2
                else:
                    s1 = ((sign - 1.0) / count * 2.0 * sin_mu[mu_index] /
                          delta**2)
                    s2 = (-2.0 * sign * sin_mu[mu_index] /
                          (sin_nu[i] * delta))
                st1[i] += s1 * zt[mu_index]
                st2[i] += s2 * zt[mu_index]

                if same:
                    s3 = count / sin_nu[i]
                else:
                    s3 = ((sign - 1.0) / count * 2.0 * sin_mu[mu_index] /
                          delta**2 + 2.0 / count * (1.0 - sign) /
                          (sin_mu[mu_index] * delta))
                st3[i] += s3 * zt[mu_index]

                sign_mu = (-1.0)**(mu_index + 1)
                if same:
                    s4 = (count / sin_nu[i] -
                          2.0 * (sign_mu - 1.0) /
                          (count * sin_nu[i] * (1.0 - cos_nu[i])))
                    s5 = -cos_nu[i] / sin_nu[i]**2
                else:
                    s4 = (2.0 * (sign - 1.0) / (count * sin_nu[i]) *
                          (1.0 - cos_mu[mu_index] * cos_nu[i]) /
                          (cos_nu[i] - cos_mu[mu_index])**2 -
                          2.0 * (sign_mu - 1.0) /
                          (count * sin_nu[i] * (1.0 - cos_mu[mu_index])))
                    s5 = -2.0 * sign / delta
                st4[i] += s4 * zc[mu_index]
                st5[i] += s5 * zc[mu_index]

        st1[count - 1] += count * a0

        # ST3NU picks up a separate leading-edge term.
        for i in range(count):
            if i == count - 1:
                st3[i] += count / 2.0 * a0
            else:
                s3 = ((-1.0)**(i + 1) - 1.0) / (count * (1.0 + cos_nu[i]))
                st3[i] += s3 * np.sqrt(leading_edge_radius / 2.0)

    # ST2NU, ST4NU and ST5NU all divide by sin(THNU) and are singular at
    # the last station; SLOPE never reads any of them there.
    st2[count - 1] = np.nan
    st4[count - 1] = np.nan
    st5[count - 1] = np.nan

    result = {
        'x': x,
        'theta_nu': theta_nu,
        'theta_mu': theta_mu,
        'thickness': zt,
        'camber': zc,
        'st1': st1, 'st2': st2, 'st3': st3, 'st4': st4, 'st5': st5,
        'a0': float(a0),
        'method': 'legacy_ideal',
    }

    if cli is None:
        derived = _section_parameters(zc, theta_nu, st5, count)
        result.update(derived)
        alpha_zero_lift = derived['alpha_zero_lift']
    else:
        result.update({'cli': float(cli),
                       'alpha_zero_lift': float(alpha_zero_lift),
                       'ideal_alpha': None})

    # Quarter-chord pitching moment, always computed.
    moment_sum = float(np.sum(zc[:last] * cos_nu[:last]))
    xmuo = moment_sum * (-PI / count)
    result['cm_c4'] = float(2.0 * xmuo + PI / 2.0 * alpha_zero_lift / RAD)
    return result


def _section_parameters(zc, theta_nu, st5, count: int) -> Dict[str, float]:
    """The source's ideal angle, zero-lift angle and ideal lift."""
    last = count - 1
    index_n_minus_4 = count - 4             # one-based N-4
    camber_at_95 = zc[3]                    # ZCNU(4)
    camber_at_05 = zc[index_n_minus_4 - 1]  # ZCNU(N-4)
    cos_nu = np.cos(theta_nu)
    sin_nu = np.sin(theta_nu)

    dai05 = .3739 * camber_at_05 + .04745 * st5[count - 2]    # ST5NU(N-1)
    dai95 = -.3739 * camber_at_95 + .04745 * st5[0]
    da95 = -.7834 * camber_at_95 + .09518 * st5[0]

    t1 = 0.5 * (-camber_at_95 * cos_nu[3] / (sin_nu[3] / 2.0)**2 -
                camber_at_05 * cos_nu[index_n_minus_4 - 1] /
                (sin_nu[index_n_minus_4 - 1] / 2.0)**2)
    t2 = 0.5 * camber_at_95 / (1.0 - cos_nu[3])

    index_n_minus_5 = count - 5
    sum_ai = 0.0
    sum_al = 0.0
    for i in range(5, last + 1):            # source I = 5 .. L, one-based
        position = i - 1                    # zero-based
        if i <= index_n_minus_5:
            sum_ai += -zc[position] * cos_nu[position] / \
                (sin_nu[position] / 2.0)**2
        sum_al += zc[position] / (1.0 - cos_nu[position])

    ideal = 0.5 / count * (t1 + sum_ai)
    ideal_alpha = RAD * (dai05 + ideal + dai95)
    zero = -2.0 / count * (t2 + sum_al)
    alpha_zero = RAD * (da95 + zero)
    return {
        'ideal_alpha': float(ideal_alpha),
        'alpha_zero_lift': float(alpha_zero),
        'cli': float(2.0 * PI / RAD * (ideal_alpha - alpha_zero)),
    }
