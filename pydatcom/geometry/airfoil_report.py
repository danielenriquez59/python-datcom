"""
The airfoil-section report: THEORY.

THEORY runs the section analysis (IDEAL, then SLOPE at zero Mach and at
each case Mach number below the crest critical Mach number) and prints the
section characteristics.

Reference: datcom-legacy/datcom_2000/theory.f
"""

from typing import Callable, Dict, List, Sequence, Tuple

from pydatcom.io.fortran_format import fortran_write
from pydatcom.utils.constants import UNUSED

_F1020 = ('(1H1,29X,75HAUTOMATED STABILITY AND CONTROL METHODS PER APRIL '
          '1976 VERSION OF DATCOM    )')
_TITLES = {'W': '(1H ,55X,23HWING SECTION DEFINITION)',
           'H': '(1H ,50X,34HHORIZONTAL TAIL SECTION DEFINITION)',
           'V': '(1H ,51X,32HVERTICAL TAIL SECTION DEFINITION)',
           'F': '(1H ,52X,30HVENTRAL FIN SECTION DEFINITION)'}
_F1070 = ('(1H0,33X,23HIDEAL ANGLE OF ATTACK =,F10.5,5H DEG.//'
          '30X,27HZERO LIFT ANGLE OF ATTACK =,F10.5,5H DEG.//'
          '33X,24HIDEAL LIFT COEFFICIENT =,F10.5//'
          '18X,39HZERO LIFT PITCHING MOMENT COEFFICIENT =,F10.5//'
          '29X,28HMACH ZERO LIFT-CURVE-SLOPE =,F10.5,6H /DEG.//'
          '36X,21HLEADING EDGE RADIUS =,F10.5,15H FRACTION CHORD//'
          '30X,27HMAXIMUM AIRFOIL THICKNESS =,F10.5,15H FRACTION CHORD//'
          '48X, 9HDELTA-Y =,F10.5,14H PERCENT CHORD//)')
_F1080 = ('(1H0,25X,5HMACH=,F7.4,19H LIFT-CURVE-SLOPE =,F10.5,'
          '6H /DEG.,5X,6H XAC =,F10.5)')
_F1090 = ('(1H0,36X,43H*** CREST CRITICAL MACH NUMBER EXCEEDED ***//'
          '35X,22H CREST CRITICAL MACH =,F10.5//'
          '46X,11H LOCATION =,F10.5,15H FRACTION CHORD//'
          '38X,19H LIFT-CURVE-SLOPE =,F10.5,6H /DEG.)')
_F1110 = ('(   1X,69H     AIRFOIL MODULE UNABLE TO CALCULATE SECTION CLALPA'
          ' FOR MACH NO = ,F7.3,50H  THEREFORE A VALUE MUST BE SPECIFIED BY '
          'THE USER./1X,123H     SEE USERS MANUAL FOR SECTION '
          'CHARACTERISTICS INPUT.  (NOTE% CLALPA AND CLMAX SHOULD REFERENCE '
          'FOOTNOTE 1 ON PAGE 29. )/)')


def theory(section: Dict[str, object], machs: Sequence[float],
           reynolds: Sequence[float], cbar: float,
           run_ideal: Callable[[], object],
           slope: Callable[[float, float], Tuple[float, float]]
           ) -> Dict[str, object]:
    """Translate THEORY.

    Args:
        section: The section words after IDEAL: ``ai``, ``alo``, ``cli``,
            ``cmco4``, ``rho``, ``tmax``, ``deltay``, ``mcc``, ``xc``,
            ``clcc``, and ``cla``/``xac`` (per Mach, ``UNUSED`` where not
            given), and ``surface`` (``NACA(6)``: W, H, V or F).
        machs, reynolds: The case Mach numbers and unit Reynolds numbers.
        cbar: The reference chord.
        run_ideal: IDEAL.  slope: SLOPE, ``(mach, Re) -> (cla, xac)``;
            ``cla`` is ``UNUSED`` when it cannot be formed.

    Returns:
        ``lines``, ``cla0``, ``cla`` and ``xac`` as filled, and ``exit``:
        the run stops (EXIT) when a slope could not be formed.
    """
    s = section
    run_ideal()
    cla0, _ = slope(0.0, 0.0)
    lines = fortran_write(_F1020)
    if s.get('surface') in _TITLES:
        lines += fortran_write(_TITLES[s['surface']])
    lines += fortran_write(_F1070, [s['ai'], s['alo'], s['cli'], s['cmco4'],
                                    cla0, s['rho'], s['tmax'], s['deltay']])
    cla, xac = list(s['cla']), list(s['xac'])
    past_crest_critical = slope_failed = False
    for mach_index, mach in enumerate(machs):
        reynolds_at_cbar = reynolds[mach_index] * cbar
        if mach >= s['mcc']:
            past_crest_critical = True
        if cla[mach_index] != UNUSED:
            continue
        if mach >= s['mcc']:
            cla[mach_index] = s['clcc']
        if past_crest_critical:
            continue
        cla[mach_index], xac[mach_index] = slope(mach, reynolds_at_cbar)
        if cla[mach_index] != UNUSED:
            lines += fortran_write(_F1080, [mach, cla[mach_index],
                                            xac[mach_index]])
        else:
            slope_failed = True
            lines += fortran_write(_F1110, [mach])
    if past_crest_critical:
        lines += fortran_write(_F1090, [s['mcc'], s['xc'], s['clcc']])
    return {'lines': lines, 'cla0': cla0, 'cla': cla, 'xac': xac,
            'exit': slope_failed}
