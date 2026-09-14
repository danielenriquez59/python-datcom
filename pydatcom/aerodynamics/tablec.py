"""
TABLEC: regression coefficients for the wing-body zero-lift moment.

Nineteen coefficients, each interpolated over the same fourteen-point Mach
grid.  The source stores them as one flat 266-element array split across
four DATA blocks, with coefficient ``i`` occupying positions
``14*(i-1)+1`` through ``14*i``.

The coefficient values here were extracted from the source DATA statements
by parsing rather than by hand, since 266 numbers transcribed by eye is
exactly where a silent error would hide.

Reference: datcom-legacy/datcom_2000/tablec.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The source's fourteen-point Mach grid, ZM.
_MACH = [0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 2.5]

_CCM = [
    # coefficient 1
    -6.1496E-01, 3.0521E-02, 1.0449E-02, 2.4846E-02, 2.3233E-02, -2.9503E-02, -7.2508E-02, -4.7882E-02, -7.7024E-03, -9.0333E-03, -6.4029E-02, 1.4988E-02, -1.0257E-01, -9.7877E-02,
    # coefficient 2
    3.4113E-02, -2.6482E-02, 1.9975E-03, -5.5509E-03, -5.4868E-04, 5.6041E-02, 1.1805E-01, 8.7861E-02, 6.7409E-02, 8.1884E-02, 1.1973E-01, 6.6427E-02, 1.7653E-01, 1.6776E-01,
    # coefficient 3
    -4.3111E-03, -3.1117E-03, -1.1029E-03, -2.4707E-03, -2.7806E-03, 5.2306E-03, 9.1112E-03, 1.2725E-02, 3.7585E-03, 2.9019E-03, 1.2987E-02, 7.8970E-04, 2.2724E-02, 2.8670E-02,
    # coefficient 4
    -6.7607E-03, 4.0854E-04, 1.9192E-04, -1.4904E-03, -1.4790E-03, -2.1117E-03, -5.5601E-03, 5.8939E-04, -2.0553E-03, -3.5078E-03, 1.4900E-03, -6.5439E-03, -6.5667E-03, -4.3908E-03,
    # coefficient 5
    5.8626E-01, 4.3533E-02, 2.8208E-02, 3.1764E-02, 2.4266E-02, 2.4656E-02, 1.6164E-02, 4.9004E-02, -1.3639E-02, -2.7873E-02, -2.9892E-02, -2.2854E-02, 2.6084E-03, 2.0467E-03,
    # coefficient 6
    -2.1631E-03, 1.0412E-04, 7.9415E-04, 3.6833E-04, 4.0862E-04, 8.2755E-04, -1.9504E-03, -1.8951E-04, -8.4471E-04, -1.6508E-03, -1.3285E-03, -9.8949E-04, -3.6847E-04, -1.3417E-03,
    # coefficient 7
    -2.9718E-03, -2.3971E-03, -1.2200E-03, 3.7272E-04, -9.2113E-04, -1.6469E-03, 4.3904E-03, -1.4263E-03, -2.7075E-03, -7.0990E-04, -1.9031E-03, -2.4785E-03, 1.5647E-03, -7.7258E-04,
    # coefficient 8
    -1.3161E-03, 1.0412E-02, 5.1328E-03, -3.5006E-03, 9.1995E-04, -8.5456E-03, -5.0916E-02, -7.3144E-03, -3.5333E-04, -1.0550E-02, 9.3873E-03, -1.1328E-02, -2.6630E-02, -2.3909E-02,
    # coefficient 9
    3.4451E-03, -4.9864E-03, -2.9594E-03, 7.3299E-04, -1.0519E-03, 9.2825E-03, 3.9982E-02, 2.5290E-03, -3.0929E-04, 3.2408E-03, -4.9618E-03, 2.6880E-03, 1.0089E-02, 9.4275E-03,
    # coefficient 10
    -9.7535E-03, -2.3881E-03, 7.4312E-04, 1.2306E-03, 2.1653E-03, 3.7542E-03, 1.1703E-04, -1.1616E-04, 1.7132E-03, 5.7507E-04, 1.1079E-03, -2.3730E-05, -1.8213E-03, 5.7215E-04,
    # coefficient 11
    -3.9389E+00, 5.6516E-01, 1.2045E-01, -1.0408E-01, 9.8669E-02, -8.2906E-01, 1.0256E+00, -6.2097E-01, -1.5639E+00, -4.3130E-01, -8.3720E-01, -1.6698E-01, 9.7612E-01, 1.3090E+00,
    # coefficient 12
    -5.2560E-01, -5.1829E-01, -5.3503E-01, -4.9728E-01, -6.2200E-01, -3.1518E-01, -3.9402E-01, -3.6576E-01, -2.1108E-01, -1.0254E-01, -2.9314E-02, 4.5021E-01, -4.7537E+00, -5.3830E+00,
    # coefficient 13
    -4.3982E+00, -4.0930E+00, -4.2059E+00, -3.9283E+00, -4.8034E+00, -3.1686E+00, -3.4999E+00, -3.1159E+00, -2.5628E+00, -2.0192E+00, -1.7911E+00, 1.3537E-01, -1.5553E+01, -1.6996E+01,
    # coefficient 14
    3.6116E-02, -6.9271E-03, -1.3549E-03, 4.8901E-03, 1.8263E-03, 7.5362E-03, 3.5393E-03, 1.0248E-02, 1.1013E-02, 6.7985E-03, 1.7144E-03, 2.2793E-03, -1.0677E-03, 1.4219E-02,
    # coefficient 15
    3.4068E-02, 1.1685E-02, 5.3796E-03, 1.4643E-03, 2.8848E-03, -3.8225E-03, 3.0079E-02, 4.7945E-03, 1.6835E-03, 8.3637E-03, 3.4691E-03, 4.0738E-03, -4.2213E-03, -9.1883E-03,
    # coefficient 16
    -7.3114E-02, -1.0598E-03, 2.0454E-03, 2.7298E-03, 2.5683E-04, -3.5619E-03, 1.6422E-03, -6.5869E-04, 1.6923E-03, 1.3512E-03, 1.1693E-03, 3.2272E-03, 4.5941E-03, 2.6961E-03,
    # coefficient 17
    1.2683E+00, -1.7713E-02, -3.9115E-03, -9.2493E-03, 6.5923E-03, 7.4832E-03, 2.5312E-02, -3.7542E-02, 1.7427E-02, 2.0635E-02, 1.8158E-02, 1.6570E-02, -6.3330E-04, -1.6885E-02,
    # coefficient 18
    -1.6503E-01, -5.2268E-02, -4.3905E-02, -5.1001E-02, -1.0221E-01, -1.9949E-02, -5.3897E-02, 6.3401E-03, -1.3287E-01, -1.3673E-01, -8.0305E-02, -1.8294E-01, -1.1299E-01, -6.8530E-02,
    # coefficient 19
    -7.7192E-04, 2.5439E-04, -1.2391E-03, -1.1483E-03, -7.8181E-04, -1.0290E-03, -1.3250E-03, -5.1180E-04, -1.3416E-03, -2.3095E-03, -2.2940E-03, -1.4919E-03, 9.4537E-04, -9.9505E-04,
]

_COEFFICIENTS = 19
_MACH_POINTS = 14
_CCM = np.asarray(_CCM, dtype=float).reshape(
    (_COEFFICIENTS, _MACH_POINTS))


def calculate_tablec(mach: float) -> Dict[str, object]:
    """Translate TABLEC: the nineteen WBCMO regression coefficients.

    Args:
        mach: Free-stream Mach number.

    Returns:
        Dictionary with ``c``, the nineteen coefficients at that Mach, and
        the grid they were read on.

    Notes:
        Each coefficient is read with TBFUNX end mode 0 at both ends, so a
        Mach outside the 0.40 to 2.50 grid is clamped to its nearest
        endpoint rather than extrapolated.
    """
    values = np.array([
        tbfunx(_MACH, _CCM[index], float(mach), lower=0, upper=0)[0]
        for index in range(_COEFFICIENTS)])
    return {
        'c': values,
        'mach_grid': np.asarray(_MACH, dtype=float),
        'method': 'legacy_tablec',
    }
