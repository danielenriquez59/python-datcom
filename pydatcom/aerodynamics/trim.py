"""
Flap drag and trim: DRAGFP, TRIMRT, TRIMR2 and their overlay M38O46.

- ``DRAGFP``: the flapped surface's induced-drag increment at every angle
  and deflection (the Section 6.1.7 series over the span-load
  coefficients the flap routines leave), and the minimum-drag increment
  of Figures 6.1.7-22/-23 with the Figure 6.1.7-24 induced factor.
- ``TRIMRT``: trim by a control on the lifting surface, interpolating the
  flap routines' increments at the deflection that zeroes the untrimmed
  moment.
- ``TRIMR2``: trim by an all-moving horizontal tail, solving for the tail
  lift that balances the wing-body moment and the incidence that gives it.

M38O46 runs TRIMR2 for an all-moving tail without symmetric flaps (and a
tail type of at most 2); otherwise DRAGFP for flap types 1 to 6 and
TRIMRT when trimming.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/dragfp.f, trimrt.f, trimr2.f,
m38o46.f
"""

import math
from typing import Dict, Mapping, Optional, Sequence

import numpy as np

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex, tlinex_flat

# DRAGFP: the Section 6.1.7 induced-drag series coefficients.
_TEC = np.array([
    4.444, 0.4794, 0.1602, 0.06946, 0.02782,
])
_TOTC = np.array([
    38.65, 19.52, 13.24, 10.17, 8.398, 7.278, 6.538, 6.046,
    5.732, 5.557,
])
_SPU = np.array([
    0.0, 0.1423, 0.2817, 0.4153, 0.5407, 0.6549, 0.7557, 0.8412,
    0.9097, 0.9595, 0.9898, 1.0,
])
_SC1 = np.array([
    0.09772, 0.1463, 0.3188, 1.12, 13.95,
])
_SC2 = np.array([
    0.0493, 0.1232, 0.231, 0.7268, 7.616, 7.046,
])
_SC3 = np.array([
    0.117, 0.1928, 0.5431, 5.246, 5.167,
])
_SC4 = np.array([
    0.06423, 0.1739, 0.4384, 4.084, 4.029, 0.2948,
])
_SC5 = np.array([
    0.1753, 0.3869, 3.398, 3.372, 0.3127,
])
_SC6 = np.array([
    0.106, 0.3628, 2.951, 2.945, 0.2984, 0.06003,
])
_SC7 = np.array([
    0.3799, 2.673, 2.658, 0.285, 0.07869,
])
_SC8 = np.array([
    0.2635, 2.533, 2.459, 0.2785, 0.08798, 0.02287,
])
_SC9 = np.array([
    2.565, 2.401, 0.2858, 0.098, 0.03618,
])
_SC10 = np.array([
    2.244, 2.487, 0.2883, 0.116, 0.0491, 0.01405,
])
# Figure 6.1.7-22 (split and plain flaps) and -23 (other types): the
# profile-drag increment against deflection and flap-chord ratio.
_X11722 = np.array([
    0.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0,
    45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0,
])
_X21722 = np.array([
    0.0, 0.1, 0.16, 0.2, 0.249, 0.28, 0.32, 0.36,
    0.42,
])
_Y11722 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.002, 0.003, 0.004, 0.007, 0.008, 0.012,
    0.016, 0.025, 0.0, 0.005, 0.008, 0.012, 0.016, 0.019,
    0.025, 0.032, 0.045, 0.0, 0.009, 0.015, 0.02, 0.027,
    0.032, 0.039, 0.048, 0.07, 0.0, 0.014, 0.024, 0.03,
    0.041, 0.048, 0.06, 0.074, 0.105, 0.0, 0.019, 0.033,
    0.043, 0.057, 0.069, 0.088, 0.11, 0.151, 0.0, 0.026,
    0.043, 0.057, 0.077, 0.092, 0.114, 0.14, 0.191, 0.0,
    0.033, 0.054, 0.071, 0.097, 0.117, 0.15, 0.187, 0.247,
    0.0, 0.04, 0.067, 0.088, 0.119, 0.143, 0.184, 0.232,
    0.305, 0.0, 0.048, 0.08, 0.104, 0.14, 0.17, 0.218,
    0.264, 0.337, 0.0, 0.057, 0.093, 0.122, 0.165, 0.2,
    0.254, 0.307, 0.388, 0.0, 0.065, 0.108, 0.14, 0.192,
    0.238, 0.303, 0.366, 0.462, 0.0, 0.075, 0.122, 0.16,
    0.22, 0.273, 0.342, 0.414, 0.521, 0.0, 0.083, 0.138,
    0.179, 0.25, 0.305, 0.38, 0.451, 0.564, 0.0, 0.09,
    0.15, 0.199, 0.28, 0.345, 0.433, 0.52, 0.658,
])
_X11723 = np.array([
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0,
])
_X21723 = np.array([
    0.0, 0.1, 0.23, 0.36, 0.4,
])
_Y11723 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.002, 0.004,
    0.009, 0.012, 0.0, 0.004, 0.01, 0.018, 0.023, 0.0,
    0.011, 0.029, 0.05, 0.058, 0.0, 0.02, 0.057, 0.101,
    0.117, 0.0, 0.031, 0.088, 0.171, 0.205,
])
# Figure 6.1.7-24A-C: the induced-drag factor K'.
_X11724 = np.array([
    4.0, 6.0, 12.0,
])
_X21724 = np.array([
    0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    0.9, 1.0,
])
_Y61724 = np.array([
    2.34, 1.74, 1.41, 0.96, 0.62, 0.38, 0.22, 0.11,
    0.03, 0.0, 3.2, 2.54, 2.04, 1.24, 0.79, 0.5,
    0.3, 0.15, 0.07, 0.0, 4.9, 3.8, 2.8, 1.7,
    1.06, 0.65, 0.37, 0.18, 0.06, 0.0,
])


def _sqrt(x: float) -> float:
    """SQRT as the source runs it: NaN, not an error, below zero."""
    return math.sqrt(np.float64(x)) if x >= 0.0 else math.nan


def _series(gppa):
    """DRAGFP's two sums over the eleven span stations: (TCB, TCL)/A31422."""
    q = [0.0] * 13                   # Q(1..12): EN1, EN2, SM1..SM10
    k = l_ = 0
    for span_station in range(1, 12):
        g = gppa[span_station]
        if span_station % 2 == 0:
            q[1] += _TEC[k] * g
            for q_index, table in ((3, _SC1), (5, _SC3), (7, _SC5), (9, _SC7),
                                   (11, _SC9)):
                q[q_index] += table[k] * g
            k += 1
        else:
            for q_index, table in ((4, _SC2), (6, _SC4), (8, _SC6), (10, _SC8),
                                   (12, _SC10)):
                q[q_index] += table[l_] * g
            l_ += 1
    en1 = gppa[1] * (5.5 * gppa[1] - q[1])
    clsm = totsm = 0.0
    for span_station in range(2, 12):
        mirror = 13 - span_station
        temp = _SPU[span_station - 1] * gppa[mirror]
        totsm = (temp * (_TOTC[span_station - 2] * gppa[mirror] -
                         q[span_station + 1]) + totsm)
        clsm += temp
    return en1 + 2.0 * totsm, gppa[1] + 2.0 * clsm


def calculate_dragfp(alpha_deg: Sequence[float], flap: Dict[str, object],
                     wing: Dict[str, float], tail: Dict[str, float],
                     sref: float) -> Dict[str, object]:
    """Translate DRAGFP: flap induced and minimum drag increments.

    Args:
        alpha_deg: The wing's local angles ``B(23)`` onward.
        flap: ``delta`` (``F(1..)``), ``ftype`` (``F(17)``), and the flap
            routines' words: ``gd1``, ``gd2``, ``gd3`` (``TCD(29..)``,
            ``(1..)``, ``(15..)``, twelve each), ``adave`` (``FCM(63..)``),
            ``eta`` (``FLP(1..5)``), ``rkb`` (``FLP(20..23)``), ``cfoca``
            (``FLP(61)``) and ``dcl`` (``WING(201..)``).
        wing: ``a3``, ``a7`` and ``b49`` (the wing's area, aspect ratio and
            ``B(49)``).
        tail: ``a3``, ``a7`` (``AHT(3)``, ``AHT(7)``, ``UNUSED`` without a
            tail), ``alpha`` (``BHT(23..)``) and ``epsilon``
            (``DWASH(21..)``).
        sref: ``SREF``.

    Returns:
        ``dcdi`` (``BODY(201..)``, angle fastest then deflection), and for
        flap types up to 5 ``delcdf`` (``TCD(49..)``), ``delcdm``
        (``WING(231..)``) and ``kprm`` (``TCD(47)``).

    Notes:
        With a horizontal tail present (``AHT(3)``, ``AHT(7)`` not
        ``UNUSED``) the tail's area and aspect ratio replace the wing's
        throughout, and the angle is the tail's local angle less the
        wing's ``B(49)`` and the downwash.  The lift of the same series is
        formed and never stored.  Kept.
    """
    delta = [float(d) for d in flap['delta']]
    ndelta = len(delta)
    ftype = float(flap['ftype'])
    aspect, scale = float(wing['a7']), float(wing['a3']) / sref
    if float(tail['a3']) != UNUSED:
        scale = float(tail['a3']) / sref
    if float(tail['a7']) != UNUSED:
        aspect = float(tail['a7'])
    gd1, gd2, gd3 = (np.asarray(flap[k], dtype=float)
                     for k in ('gd1', 'gd2', 'gd3'))
    span_load_slope = gd1 / RAD
    span_load_delta = gd3 - gd2
    pi_ar_over_22 = PI * aspect / 22.
    b49 = float(wing['b49'])
    effective_alpha = []
    for angle_index, alpha_deg_value in enumerate(alpha_deg):
        value = float(alpha_deg_value) - b49
        if float(tail['alpha'][angle_index]) != UNUSED:
            value = (float(tail['alpha'][angle_index]) - b49 -
                     float(tail['epsilon'][angle_index]))
        effective_alpha.append(value)
    nalpha = len(effective_alpha)
    cd_base = np.zeros(nalpha)
    cd_by_delta = np.zeros((ndelta, nalpha))
    for deflection_slot in range(1, ndelta + 2):
        adaved = (float(flap['adave'][deflection_slot - 2]) *
                  delta[deflection_slot - 2] / RAD
                  if deflection_slot > 1 else 0.0)
        for angle_index in range(nalpha):
            gppa = [0.0] + [
                span_load_slope[span_load_index] * effective_alpha[angle_index] -
                (span_load_delta[span_load_index] * adaved
                 if deflection_slot > 1 else 0.0)
                for span_load_index in range(11)]
            induced = pi_ar_over_22 * _series(gppa)[0]
            if deflection_slot == 1:
                cd_base[angle_index] = induced * scale
            else:
                cd_by_delta[deflection_slot - 2, angle_index] = induced * scale
    result: Dict[str, object] = {
        'dcdi': (cd_by_delta - cd_base[None, :]).ravel(),
        'alpha': effective_alpha,
        'scale': scale,
        'aspect': aspect,
        'method': 'legacy_dragfp',
    }
    if ftype > 5.0:
        return result
    eta = flap['eta']
    kprm = tlinex_flat(_X11724, _X21724, _Y61724, aspect,
                       float(eta[4]) - float(eta[0]), 2, 0, 0, 0)
    kb = sum(float(r) for r in flap['rkb'][:4])
    arg1 = kprm / (aspect * PI)
    delcdf, delcdm = [], []
    for deflection_index, deflection_deg in enumerate(delta):
        if ftype in (1.0, 5.0):
            value = tlinex_flat(_X11722, _X21722, _Y11722, abs(deflection_deg),
                                float(flap['cfoca']), 2, 0, 0, 0)
        else:
            value = tlinex_flat(_X11723, _X21723, _Y11723, abs(deflection_deg),
                                float(flap['cfoca']), 2, 0, 0, 0)
        delcdf.append(value)
        delcdm.append(
            value * kb * scale +
            arg1 * float(flap['dcl'][deflection_index])**2 / scale)
    result.update({'kprm': kprm, 'kb': kb, 'delcdf': delcdf,
                   'delcdm': delcdm})
    return result


def calculate_trimrt(alpha_deg: Sequence[float], epsilon: Sequence[float],
                     untrimmed: Dict[str, Sequence[float]],
                     alpha_clmax: float, flap: Dict[str, object],
                     cdi: Optional[Sequence[float]] = None
                     ) -> Dict[str, object]:
    """Translate TRIMRT: trim by a control on the lifting surface.

    At each angle, the deflection whose moment increment cancels the
    untrimmed moment, and the flap increments there.  The pass stops at
    the first angle past ``alpha_clmax`` (``TSTOP = 2``) or where the
    moment lies outside the increments' range (``TSTOP = 1``).

    Args:
        alpha_deg: ``FLC(23)`` onward.
        epsilon: ``DWASH(21)`` onward; the trim angle is ``alpha - eps``.
        untrimmed: ``cl``, ``cm``, ``cd``, as TRIMRT chooses them: the wing
            alone, the wing-body with a body, or the wing-body-tail
            (with the vertical panels' drag when present) with a tail.
        alpha_clmax: ``B(43)``, or ``BHT(43)`` with wing and tail.
        flap: ``delta`` (``F(1..)``), ``ftype`` (``F(17)``) and the flap
            routines' ``dcm`` (``WING(211..)``), ``dcl`` (``201..``),
            ``dclmax`` (``221..``), ``dcdmin`` (``231..``), ``chd``
            (``261..``).
        cdi: DRAGFP's ``dcdi`` (``BODY(201..)``).

    Returns:
        ``ntrim``, ``tstop`` (``None`` when the pass completes: the source
        leaves ``TRM(22)`` as it was), ``alpha`` (``TRM(1..)``), the
        untrimmed curves, and per trimmed angle ``deltat``, ``dclt``,
        ``clmaxt``, ``cdmint``, ``chdt``, ``cdit`` (``VT(201)``, ``(221)``,
        ``(241)``, ``(281)``, ``(321)``, ``(261)`` onward) as the flap type
        allows.

    Notes:
        The source uses one deflection fewer than it was given
        (``NDELTA = F(16)+.5`` then ``NDELTA = NDELTA - 1``) in every
        lookup, in both bundled sources.  Kept.
    """
    ftype = float(flap['ftype'])
    ndelta = len(flap['delta']) - 1
    delta = np.asarray(flap['delta'][:ndelta], dtype=float)
    series = {k: np.asarray(flap[k][:ndelta], dtype=float)
              for k in ('dcm', 'dcl', 'dclmax', 'dcdmin', 'chd') if k in flap}
    dcm = series['dcm']
    if dcm[0] > dcm[-1]:
        delta_sorted, dcm_sorted = delta[::-1], dcm[::-1]
    else:
        delta_sorted, dcm_sorted = delta, dcm
    alpha = np.asarray(alpha_deg, dtype=float) - np.asarray(epsilon,
                                                            dtype=float)
    nalpha = len(alpha)
    untrimmed_cm = [float(v) for v in untrimmed['cm']]
    y_cdi = None
    if cdi is not None and ftype <= 6.0:
        y_cdi = np.asarray(cdi, dtype=float)[:nalpha * ndelta].reshape(
            ndelta, nalpha).T
    trimmed = {k: [] for k in ('deltat', 'dclt', 'clmaxt', 'cdmint', 'chdt',
                               'cdit')}
    ntrim, tstop = nalpha, None
    for angle_index in range(nalpha):
        if alpha[angle_index] > alpha_clmax:
            ntrim, tstop = angle_index, 2.0
            break
        moment_to_cancel = -untrimmed_cm[angle_index]
        if moment_to_cancel < dcm_sorted[0] or moment_to_cancel > dcm_sorted[-1]:
            ntrim, tstop = angle_index, 1.0
            break
        trim_delta = tbfunx(dcm_sorted, delta_sorted, moment_to_cancel,
                            0, 0, ordered=False)[0]
        trimmed['deltat'].append(trim_delta)
        trimmed['dclt'].append(
            tbfunx(delta, series['dcl'], trim_delta, 0, 0)[0])
        if ftype <= 5.0:
            trimmed['clmaxt'].append(
                tbfunx(delta, series['dclmax'], trim_delta, 0, 0)[0])
        if ftype <= 2.0:
            trimmed['cdmint'].append(
                tbfunx(delta, series['dcdmin'], trim_delta, 0, 0)[0])
        if ftype == 1.0:
            trimmed['chdt'].append(
                tbfunx(delta, series['chd'], trim_delta, 0, 0)[0])
        if y_cdi is not None:
            trimmed['cdit'].append(float(tlinex(
                delta, alpha, y_cdi, trim_delta, alpha[angle_index], 0, 0, 0, 0)))
    return {
        'ntrim': ntrim,
        'tstop': tstop,
        'alpha': alpha,
        'utcl': list(untrimmed['cl']),
        'utcm': untrimmed_cm,
        'utcd': list(untrimmed['cd']),
        **{k: [float(x) for x in v] for k, v in trimmed.items()},
        'method': 'legacy_trimrt',
    }


def calculate_trimr2(alpha_deg: Sequence[float], epsilon: Sequence[float],
                     q_ratio: Sequence[float], wing_body: Dict[str, object],
                     tail: Dict[str, object], wbt: Dict[str, object],
                     position: Dict[str, float], sref: float,
                     cbarr: float) -> Dict[str, object]:
    """Translate TRIMR2: trim by an all-moving horizontal tail.

    At each angle the tail lift that balances the wing-body moment and the
    tail's own zero-lift moment and drag about the CG (the quadratic in
    ``CL_h`` with the tail's induced drag), the incidence that produces
    it through the carryover factors, and the trimmed lift, drag, moment
    and hinge moment.  The pass stops past the tail's stall angle
    (``TSTOP = 2``) or at an unavailable wing-body moment (``TSTOP = 3``);
    an incidence outside ``+-alpha_clmax`` marks the angle with -1000.

    Args:
        alpha_deg, epsilon, q_ratio: ``FLC(23..)``, ``DWASH(21..)``,
            ``DWASH(1..)``.
        wing_body: ``cl``, ``cd``, ``cm`` (``BW(21..)``, ``(1..)``,
            ``(41..)``).
        tail: ``alpha`` (``BHT(23..)``), ``cl`` and ``cd`` (``HT(21..)``,
            ``(1..)``), ``cla`` (``HT(101)``), ``alpha_clmax``, ``cd0``,
            ``cm0`` (``BHT(43)``, ``(46)``, ``(47)``), ``ar`` (``AHT(7)``),
            ``area`` (``AHT(3)``), ``sspn``, ``sspne`` (``HTIN(4)``,
            ``(3)``), ``e`` (``DHT(30)``).
        wbt: ``kwb``, ``kbw``, ``kkwb``, ``kkbw`` (``WBT(1)``, ``(2)``,
            ``(150)``, ``(151)``), ``clbh`` (``WBT(130..)``), ``cdov``
            (``WBT(66)``), and the vortex factors ``f46`` and ``f68``
            (``WBT(46..)``, ``WBT(68..)``).
        position: ``xba``, ``zba`` (``BD(63)``, ``(64)``), ``xcg``,
            ``hinax`` (``SYNA(1)``, ``(11)``), ``alih`` (``SYNA(8)``).

    Returns:
        ``ntrim``, ``tstop``, and per angle reached ``aliht`` (``HT(221)``),
        ``clhtrm`` (``TRM2(1)``, the trimmed tail lift in the free-stream
        axes, also ``HT(261)``), ``cdhtrm`` (``HT(241)``), ``cmhtrm``
        (``HT(281)``), ``hmtrm`` (``HT(301)``), ``hmunt`` (``HT(201)``),
        ``clwbt`` (``VT(221)``), ``cdwbt`` (``VT(201)``); entries of angles
        marked -1000 hold ``None`` where the source leaves them.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    eps = np.asarray(epsilon, dtype=float)
    tail_inputs = {k: v for k, v in tail.items()}
    sref_over_tail_area = sref / float(tail_inputs['area'])
    tail_alpha_eff = alpha - eps
    tail_alpha_with_incidence = tail_alpha_eff + float(position['alih'])
    xbar = float(position['xba']) / cbarr
    zbar = float(position['zba']) / cbarr
    ar, e = float(tail_inputs['ar']), float(tail_inputs['e'])
    aclmax = float(tail_inputs['alpha_clmax'])
    kwb, kbw = float(wbt['kwb']), float(wbt['kbw'])
    kk = float(wbt['kkwb']) + float(wbt['kkbw'])
    keys = ('aliht', 'clhtrm', 'cdhtrm', 'cmhtrm', 'hmtrm', 'hmunt',
            'clwbt', 'cdwbt')
    trimmed = {k: [] for k in keys}
    nalpha = len(alpha)
    ntrim, tstop = nalpha, 1.0
    for angle_index in range(nalpha):
        if tail_alpha_eff[angle_index] > aclmax:
            ntrim, tstop = angle_index, 2.0
            break
        cmwb = float(wing_body['cm'][angle_index])
        if cmwb == 2.0 * UNUSED:
            ntrim, tstop = angle_index, 3.0
            break
        tail_alpha_rad = tail_alpha_eff[angle_index] / RAD
        argl = (-math.cos(tail_alpha_rad) * xbar +
                math.sin(tail_alpha_rad) * zbar)
        argd = (math.sin(tail_alpha_rad) * xbar +
                math.cos(tail_alpha_rad) * zbar)
        if not -aclmax <= tail_alpha_with_incidence[angle_index] <= aclmax:
            for k in keys:
                trimmed[k].append(None)
            trimmed['aliht'][-1] = trimmed['clhtrm'][-1] = -1000.0
            continue
        q_at_angle = float(q_ratio[angle_index])
        tail_cl_num = 2. * (
            (cmwb + q_at_angle * float(tail_inputs['cm0'])) / (argl * q_at_angle) +
            float(tail_inputs['cd0']) * argd / argl)
        tail_cl_den = 1. + _sqrt(
            1. - 4. * (argd / (argl * PI * ar * e) * sref_over_tail_area) *
            (tail_cl_num / 2.))
        tail_cl = tail_cl_num / tail_cl_den
        tail_cd = (float(tail_inputs['cd0']) +
                   tail_cl**2 / (PI * ar * e) * sref_over_tail_area)
        brackt = (float(wbt['f68'][angle_index]) *
                  float(wbt['f46'][angle_index]) *
                  (float(tail_inputs['sspn']) - float(tail_inputs['sspne'])) /
                  float(tail_inputs['sspn']))
        clhp, clahp = tbfunx(tail_inputs['alpha'], tail_inputs['cl'],
                             tail_alpha_eff[angle_index], 1, 1)
        aliht = (
            (tail_cl - (kwb + kbw) * (clhp + clahp * brackt *
                                      tail_alpha_eff[angle_index])) /
            (kk * (clahp + brackt * float(tail_inputs['cla']))))
        eps_rad = eps[angle_index] / RAD
        alpha_rad = alpha[angle_index] / RAD
        dcl = (tail_cl * math.cos(eps_rad) - tail_cd * math.sin(eps_rad)) * q_at_angle
        dcd = (tail_cd * math.cos(eps_rad) + tail_cl * math.sin(eps_rad)) * q_at_angle
        argl = (math.cos(tail_alpha_rad) * xbar +
                math.sin(tail_alpha_rad) * zbar)
        argd = (math.cos(tail_alpha_rad) * zbar -
                math.sin(tail_alpha_rad) * xbar)
        arm = (float(position['hinax']) + float(position['xba']) -
               float(position['xcg'])) / cbarr
        values = {
            'aliht': aliht,
            'clhtrm': dcl,
            'cdhtrm': dcd,
            'cmhtrm': dcl * argl + dcd * argd,
            'hmtrm': (dcl * math.cos(alpha_rad) + dcd * math.sin(alpha_rad)) * arm,
            'hmunt': (float(tail_inputs['cl'][angle_index]) * math.cos(alpha_rad) +
                      float(tail_inputs['cd'][angle_index]) * math.sin(alpha_rad)) * arm,
            'clwbt': (dcl + float(wbt['clbh'][angle_index]) +
                      float(wing_body['cl'][angle_index])),
            'cdwbt': (float(wing_body['cd'][angle_index]) + dcd +
                      float(wbt['cdov'])),
        }
        for k in keys:
            trimmed[k].append(float(values[k]))
    return {'ntrim': ntrim, 'tstop': tstop, **trimmed,
            'method': 'legacy_trimr2'}


def m38o46_routines(trim: bool, htpl: bool, symfp: bool, ftype: float,
                    tail_type: float) -> Sequence[str]:
    """M38O46's routing: which of DRAGFP, TRIMRT and TRIMR2 run.

    ``tail_type`` is ``WINGIN(101)``.  The overlay also copies
    ``WING(251)`` to ``VT(301)`` whenever it does not return early.
    """
    if trim and htpl and not symfp:
        return [] if tail_type > 2.0 else ['TRIMR2']
    routines = ['DRAGFP'] if ftype <= 6.0 else []
    return routines + (['TRIMRT'] if trim else [])
