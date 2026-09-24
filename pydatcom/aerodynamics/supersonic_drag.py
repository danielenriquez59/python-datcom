"""
Supersonic wing drag: SUPDRG and its overlay M18O22.

SUPDRG forms the wing's supersonic zero-lift drag, skin friction plus wave
drag, and for a straight tapered wing the Figure 4.1.5.2-58 drag-due-to-
lift factor.  It works from the COMMON words, so it covers the cranked and
double-delta planforms that :func:`supersonic.calculate_supdrg_straight_wing`
and :func:`supersonic.calculate_supdrg_skin_friction` (which build the
straight wing from a state dictionary) do not.

The tables are those of ``supersonic.py``, pinned to the source DATA by
test.

Reference: datcom-legacy/datcom_2000/supdrg.f, m18o22.f
"""

import math
from typing import Dict, Mapping

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.aerodynamics.supersonic import (_FIG_415127_CEPT,
                                              _FIG_415127_MACH,
                                              _FIG_415258_ROUND,
                                              _FIG_415258_SHARP,
                                              _FIG_415258_X)
from pydatcom.utils.constants import UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.table_lookup import fig26


def calculate_supdrg(mach: float, win: Mapping[int, float],
                     a: Mapping[int, float], sref: float, roughness: float,
                     sbw: float, stale: Mapping[str, float]
                     ) -> Dict[str, float]:
    """Translate SUPDRG: supersonic wing zero-lift and induced drag.

    Skin friction from FIG26 at the MAC Reynolds number, capped by the
    Figure 4.1.5.1-27 roughness cutoff: on the inboard MAC for a straight
    wing, and on the inboard then the outboard MAC, area weighted, for the
    other planforms.  Wave drag for a sharp (``KSHARP`` given) or round
    leading edge, on the inboard panel's sweep and area for a straight
    wing and the outboard panel's with ``SLG(119)`` for the others.  For a
    straight wing, the Figure 4.1.5.2-58 factor ``DRAGC`` and its planform
    parameter ``P``.

    Args:
        mach: ``FLC(N+2)``.
        win: ``WINGIN`` words 1 (``CT``), 3 (``SPANS``), 5 (``CB``), 6
            (``CR``), 15 (planform, 1 for straight tapered), 62, 63
            (``LERI``, ``LERO``), 70 (``TCEFF``), 71 (``KSHARP``, UNUSED
            for a round edge).
        a: ``A`` words 1, 2, 3 (panel and total areas), 7, 10, 15, 17
            (inboard and outboard MACs), 18 (``SIGMA``), 61 (``COSLE``), 62
            (``TANLEI``), 85, 86 (``COSLEO``, ``TANLEO``), 129 (unit
            Reynolds number).
        sref: ``SR``; roughness: ``RUFF``.
        sbw: ``SLG(119)``, the area the non-straight wave drag uses.
        stale: ``SLG`` words and a local the routine can read unset:
            ``rach`` (the Mach number FIG26 is read at, set only with
            roughness), ``rlcoff`` (89), and ``cfi`` (84, read when the
            inboard and outboard MACs are equal).

    Returns:
        The ``SLG`` words set, by name: ``beta`` (1), ``bovert`` (2),
        ``cdw`` (79), ``cdo`` (80), ``cdf`` (87), ``cf`` (88), ``rnn``
        (90), ``rlcoff`` (89), for a cranked wing ``rni``, ``cfi``,
        ``rno``, ``cfo`` (86, 84, 85, 83), and for a straight wing
        ``dragc`` (81) and ``p`` (82); with ``a62`` and ``a86``, the
        sweep tangents as the routine leaves them.

    Notes:
        Kept as executed: a zero sweep tangent is set to 1e-5 for the
        call and restored only on the straight-wing path, so a cranked
        wing leaves 1e-5 in ``A(62)`` or ``A(86)``; with no roughness the
        FIG26 Mach number is whatever the last call left; and a wing whose
        inboard and outboard MACs are equal takes its inboard friction
        from the previous call's ``SLG(84)``.
    """
    win_f = {int(k): float(v) for k, v in win.items()}
    geom = {int(k): float(v) for k, v in a.items()}
    zero_sweep = 1.0e-5
    tanleo = geom[86] if geom[86] != 0.0 else zero_sweep
    tanlei = geom[62] if geom[62] != 0.0 else zero_sweep
    straight = win_f[15] == STRAIGHT_TAPERED
    beta = math.sqrt(mach ** 2 - 1.0)
    result: Dict[str, float] = {'beta': beta}
    rach = float(stale['rach'])
    rlcoff = float(stale['rlcoff'])

    def friction(cbar):
        nonlocal rach, rlcoff
        reynolds_mac = cbar * geom[129]
        if roughness != 0.0:
            arg = 12. * cbar / roughness
            rach = min(mach, 3.0)
            cept = tbfunx(_FIG_415127_MACH, _FIG_415127_CEPT, rach, 0, 0)[0]
            rlcoff = arg**1.0482 * 10.0**cept
            if rlcoff < reynolds_mac:
                reynolds_mac = rlcoff
        return reynolds_mac, fig26(reynolds_mac, rach)

    reynolds_mac, cf = friction(geom[15])
    if straight:
        cdf = cf * geom[3] / sref * 2.0
    else:
        cfi = float(stale['cfi'])
        if geom[15] != geom[17]:
            result['rni'], result['cfi'] = reynolds_mac, cf
            cfi = cf
            reynolds_mac, cf = friction(geom[17])
        result['rno'], result['cfo'] = reynolds_mac, cf
        cdf = (cfi * geom[1] + cf * geom[2]) / sref * 2.0
    result.update({
        'rnn': reynolds_mac, 'cf': cf, 'rlcoff': rlcoff, 'cdf': cdf,
    })

    if straight:
        sweep_tangent, panel_area, cos_le = tanlei, geom[3], geom[61]
        lerbw = win_f[62] * ((win_f[6] + win_f[1]) / 2.0)
    else:
        sweep_tangent, panel_area, cos_le = tanleo, sbw, geom[85]
        lerbw = win_f[63] * ((win_f[6] + win_f[1] + 2.0 * win_f[5]) / 4.0)
    bovert = beta / sweep_tangent
    sonic_leading_edge = bovert >= 1.0
    wave_denominator = beta if sonic_leading_edge else sweep_tangent
    ksharp, tceff = win_f[71], win_f[70]
    if ksharp != UNUSED:
        volume_term = ksharp * tceff ** 2 * panel_area / sref
        cdw = volume_term / wave_denominator
    else:
        blunt_le = (1.28 * mach ** 3 * cos_le ** 6 /
                    (1.0 + mach ** 3 * cos_le ** 3))
        le_drag = 2.0 * lerbw * (2.0 * win_f[3]) / (sref * cos_le)
        thickness_term = 16.0 * tceff ** 2 * panel_area / (3.0 * sref)
        cdw = blunt_le * le_drag + thickness_term / wave_denominator
    result.update({
        'bovert': bovert,
        'cdw': cdw,
        'cdo': cdf + cdw,
        'a62': tanlei,
        'a86': tanleo,
    })
    if not straight:
        return result
    rlw = win_f[1] + geom[18] * geom[10]
    p = geom[3] / (rlw * 2.0 * win_f[3])
    table = _FIG_415258_SHARP if ksharp != UNUSED else _FIG_415258_ROUND
    result['dragc'] = interx(
        1, _FIG_415258_X, [beta * win_f[3] / rlw], [5], table,
        lind=5, lx1u=1,
    )
    result['p'] = p
    result['a62'] = 0.0 if tanlei == zero_sweep else tanlei
    result['a86'] = 0.0 if tanleo == zero_sweep else tanleo
    return result


def m18o22_options(sref: float, cbarr: float, roughness: float,
                   a4: float, a122: float) -> Dict[str, float]:
    """M18O22's ``/OPTION/`` defaults before SUPDRG.

    The overlay repeats M02O02's geometry prelude (SETUP1, WTGEOM on each
    surface and the vertical-panel halving, see
    :func:`pydatcom.geometry.wtgeom.vertical_panel_adjustments`), then
    fills an unset reference area and chord from the wing and a roughness
    below 1e-10 with 1.6e-4.  Unlike M02O02 it leaves ``BLREF`` alone and
    tests the roughness against 1e-10 rather than UNUSED.  SUPDRG runs
    when the wing is present and the case is not transonic.
    """
    return {'sref': a4 if sref == UNUSED else sref,
            'cbarr': a122 if cbarr == UNUSED else cbarr,
            'roughness': 1.6e-4 if roughness < 1.0e-10 else roughness}
