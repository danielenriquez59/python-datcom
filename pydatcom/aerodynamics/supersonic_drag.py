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
    w = {int(k): float(v) for k, v in win.items()}
    g = {int(k): float(v) for k, v in a.items()}
    tanleo = g[86] if g[86] != 0.0 else .00001
    tanlei = g[62] if g[62] != 0.0 else .00001
    straight = w[15] == STRAIGHT_TAPERED
    beta = math.sqrt(mach**2 - 1.)
    r: Dict[str, float] = {'beta': beta}
    rach = float(stale['rach'])
    rlcoff = float(stale['rlcoff'])

    def friction(cbar):
        nonlocal rach, rlcoff
        rnn = cbar * g[129]
        if roughness != 0.0:
            arg = 12. * cbar / roughness
            rach = min(mach, 3.0)
            cept = tbfunx(_FIG_415127_MACH, _FIG_415127_CEPT, rach, 0, 0)[0]
            rlcoff = arg**1.0482 * 10.0**cept
            if rlcoff < rnn:
                rnn = rlcoff
        return rnn, fig26(rnn, rach)

    rnn, cf = friction(g[15])
    if straight:
        cdf = cf * g[3] / sref * 2.
    else:
        cfi = float(stale['cfi'])
        if g[15] != g[17]:
            r['rni'], r['cfi'] = rnn, cf
            cfi = cf
            rnn, cf = friction(g[17])
        r['rno'], r['cfo'] = rnn, cf
        cdf = (cfi * g[1] + cf * g[2]) / sref * 2.
    r.update({'rnn': rnn, 'cf': cf, 'rlcoff': rlcoff, 'cdf': cdf})

    if straight:
        ta, s, ca = tanlei, g[3], g[61]
        lerbw = w[62] * ((w[6] + w[1]) / 2.)
    else:
        ta, s, ca = tanleo, sbw, g[85]
        lerbw = w[63] * ((w[6] + w[1] + 2. * w[5]) / 4.)
    bovert = beta / ta
    ksharp, tceff = w[71], w[70]
    if ksharp != UNUSED:
        arg = ksharp * tceff**2 * s / sref
        cdw = arg / ta if bovert < 1. else arg / beta
    else:
        arg1 = 1.28 * mach**3 * ca**6 / (1. + mach**3 * ca**3)
        arg2 = 2. * lerbw * (2. * w[3]) / (sref * ca)
        arg = 16. * tceff**2 * s / (3. * sref)
        cdw = arg1 * arg2 + (arg / beta if bovert >= 1. else arg / ta)
    r.update({'bovert': bovert, 'cdw': cdw, 'cdo': cdf + cdw,
              'a62': tanlei, 'a86': tanleo})
    if not straight:
        return r
    rlw = w[1] + g[18] * g[10]
    p = g[3] / (rlw * 2. * w[3])
    table = _FIG_415258_SHARP if ksharp != UNUSED else _FIG_415258_ROUND
    r['dragc'] = interx(1, _FIG_415258_X, [beta * w[3] / rlw], [5], table,
                        lind=5, lx1u=1)
    r['p'] = p
    r['a62'] = 0.0 if tanlei == 0.00001 else tanlei
    r['a86'] = 0.0 if tanleo == 0.00001 else tanleo
    return r


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
