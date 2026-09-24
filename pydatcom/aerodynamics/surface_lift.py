"""
Subsonic lifting-surface lift passes: overlays M15O17 and M16O20, and CLMCH0.

``M15O17`` (wing) and ``M16O20`` (horizontal tail) run the same three
routines on a surface's own COMMON blocks: CALCA0 for the zero-lift angle,
WTLIFT for the slope and maximum lift, and LIFTCF for the curve.

``CLMCH0`` runs that pass once at Mach zero before the Mach loop, with the
section's low-speed slope and maximum lift, and keeps the results the rest
of the program reads: the Mach-zero lift curve in ``B(3)`` onward, its slope
in ``B(48)``, and the zero-lift and stall angles in ``A(126)`` and
``A(127)``, which DWASH scales its effective aspect ratio by.

Reference: datcom-legacy/datcom_2000/m15o17.f, m16o20.f, clmch0.f and
the head of calca0.f.
"""

from typing import Dict, Optional, Sequence
import logging

from pydatcom.aerodynamics.calca0 import calculate_calca0
from pydatcom.aerodynamics.liftcf import calculate_liftcf
from pydatcom.aerodynamics.wtlift import calculate_wtlift
from pydatcom.utils.constants import UNUSED

logger = logging.getLogger(__name__)


def section_alpha_zero(swafp: float, alphai: float, cli: float,
                       section_cla: float) -> float:
    """CALCA0's starting ``A(134)``: the section zero-lift angle.

    ``AIN(10)`` when it is supplied, else ``ALPHAI - CLI/A(131)``: the
    design angle less the design lift over the section slope.
    """
    if swafp != UNUSED:
        return float(swafp)
    return float(alphai) - float(cli) / float(section_cla)


def calculate_surface_lift(planform_type: float,
                           alpha_deg: Sequence[float],
                           geometry: Dict[str, float],
                           section: Dict[str, float],
                           flight: Dict[str, float],
                           sref: float,
                           angle_state: Optional[Sequence[float]] = None
                           ) -> Dict[str, object]:
    """Translate M15O17 / M16O20: CALCA0, WTLIFT and LIFTCF in sequence.

    Args:
        planform_type: ``AIN(15)``.
        alpha_deg: ``B(23)`` onward, the surface's local angles.
        geometry: Every ``A`` entry the three routines read: CALCA0's
            ``a7``, ``a27``, ``a40``, ``a43``; WTLIFT's (see
            :func:`calculate_wtlift`); LIFTCF's (see
            :func:`calculate_liftcf`), of which ``a159`` and ``a160`` come
            from WTLIFT here and need not be given.
        section: ``swafp`` ``AIN(10)``, ``alphai`` ``AIN(20)``, ``cli``
            ``AIN(19)``, ``twista`` ``AIN(11)``, ``tovc`` ``AIN(16)``,
            ``camber`` ``AIN(64)``, ``deltay``, ``xovc``, ``sspne``, and the
            section ``cla`` ``A(131)`` and ``clmax`` ``A(132)``.
        flight: ``mach`` ``B(1)`` and ``beta`` ``B(2)``.
        sref: ``SREF``.
        angle_state: LIFTCF's ``A(147)`` to ``A(158)``; see
            :func:`calculate_liftcf`.

    Returns:
        Dictionary with ``alpha_zero_lift`` (``B(49)``), and the CALCA0,
        WTLIFT and LIFTCF results under ``calca0``, ``wtlift``, ``liftcf``.
    """
    g = {k: float(v) for k, v in geometry.items()}
    a134 = section_alpha_zero(section['swafp'], section['alphai'],
                              section['cli'], section['cla'])
    calca0 = calculate_calca0(a134, g['a27'], g['a7'], g['a40'], g['a43'],
                              float(flight['mach']),
                              twist_deg=float(section['twista']),
                              thickness_ratio=float(section['tovc']),
                              camber=bool(section['camber']))
    alpha_zero = calca0['alpha_zero_lift']
    wing_flight = dict(flight, alpha_zero_lift=alpha_zero)

    wtlift = calculate_wtlift(
        planform_type,
        {'area': g['a3'], 'aspect_ratio': g['a7'], 'taper_ratio': g['a27'],
         'sweep_le_deg': g['a34'], 'cos_le': g['a37'], 'tan_le': g['a38'],
         'tan_c2': g['a50'], 'arclss_classified': g['a124'],
         'arclss_ratio': g['a125'],
         'aspect_ratio_inboard': g.get('a5'), 'tan_c2_inboard': g.get('a74'),
         'aspect_ratio_outboard': g.get('a168'),
         'tan_c2_outboard': g.get('a98')},
        section, wing_flight, sref)
    result = {'alpha_zero_lift': alpha_zero, 'section_alpha_zero': a134,
              'calca0': calca0, 'wtlift': wtlift,
              'method': 'legacy_surface_lift'}

    # LIFTCF runs even when WTLIFT computed nothing: its curved-planform
    # method needs none of WTLIFT's results.
    liftcf = calculate_liftcf(
        planform_type, alpha_deg,
        {'area': g['a3'], 'aspect_ratio': g['a7'], 'cos_le': g['a37'],
         'arclss_factor': g['a123'], 'arclss_ratio': g['a125'],
         'a159': wtlift.get('a159'), 'a160': wtlift.get('a160'),
         'inboard_span': g.get('a23'),
         'aspect_ratio_inboard': g.get('a5'), 'tan_le': g.get('a62'),
         'planform_length': g.get('a29')},
        section, wtlift, wing_flight, sref, angle_state)
    result['liftcf'] = liftcf
    return result


def calculate_clmch0(planform_type: float,
                     alpha_deg: Sequence[float],
                     geometry: Dict[str, float],
                     section: Dict[str, float],
                     sref: float, transonic: bool = False,
                     angle_state: Optional[Sequence[float]] = None
                     ) -> Dict[str, object]:
    """Translate CLMCH0: the Mach-zero lift pass.

    Args:
        planform_type, alpha_deg, geometry, sref: As for
            :func:`calculate_surface_lift`.
        section: As for :func:`calculate_surface_lift`, except that the
            slope and maximum lift come from ``cla_mach`` (``CLALPA(1)``),
            ``clamo`` (``CLAMO``, used instead when ``transonic``) and
            ``clmaxl`` (``CLMAXL``).
        transonic: The source's ``TRANSN`` flag.
        angle_state: LIFTCF's angle state; see :func:`calculate_liftcf`.

    Returns:
        Dictionary with ``cl`` (``B(3)`` onward), ``cla`` (``B(48)``),
        ``alpha_zero_lift`` (``A(126)``), ``alpha_clmax`` (``A(127)``), and
        the full pass under ``surface``.

    Notes:
        The source also swaps ``CLMAX(1)`` for ``CLMAXL`` around the pass
        and restores it; none of the three routines reads ``CLMAX(1)``, so
        that has no effect here.
    """
    cla = float(section['clamo'] if transonic else section['cla_mach'])
    pass_section = dict(section, cla=cla, clmax=float(section['clmaxl']))
    surface = calculate_surface_lift(planform_type, alpha_deg, geometry,
                                     pass_section, {'mach': 0.0, 'beta': 1.0},
                                     sref, angle_state)
    wtlift = surface['wtlift']
    result = {'surface': surface, 'a131': cla,
              'a132': float(section['clmaxl']),
              'alpha_zero_lift': surface['alpha_zero_lift'],
              'method': 'legacy_clmch0'}
    if wtlift['computed']:
        result.update({'cl': surface['liftcf']['cl'], 'cla': wtlift['cla'],
                       'alpha_clmax': wtlift['alpha_clmax']})
    return result
