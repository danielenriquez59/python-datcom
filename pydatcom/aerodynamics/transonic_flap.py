"""
Transonic flap and control increments: TRNYRL and its overlay M40O50.

TRNYRL scales the Mach 0.6 flap results (from LIFTFP, LATFLP and GDELTA,
run beforehand at M = 0.6) to the transonic Mach number by the ratio of the
wing's (or tail's) lift-curve slope to its Mach 0.6 slope: the symmetric
flap's lift increment and lift-slope change, and the asymmetric flap's
rolling and yawing moments.  An all-moving horizontal tail's rolling
moment is formed directly, from GDELTA's span loading below Mach 1 and
Figures 4.3.1.2-12A1/A2 above.

Reference: datcom-legacy/datcom_2000/trnyrl.f, m40o50.f
"""

from typing import Dict, Mapping, Sequence

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx

_X12A1 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
_Y12A1 = [1.0, 0.97, 0.95, 0.94, 0.94, 0.94, 0.94, 0.95, 0.96, 0.98, 0.99]
_X12A2 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
_Y12A2 = [0.0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.6, 0.7, 0.8, 0.9, 1.0]
_X12222 = [0.165, 0.2, 0.25, 0.3, 0.35, 0.4]
_Y21222 = [0.83, 0.818, 0.779, 0.738, 0.68, 0.615]


def calculate_trnyrl(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate TRNYRL: transonic flap increments from Mach 0.6 results.

    Args:
        data: By name: ``mach`` (``FLC(IM+2)``), ``nalpha``, ``flap``
            (``F(16)`` count, ``F(17)`` type, ``F(18)`` control type,
            ``deltal``/``deltar`` ``F(19..)``/``F(29..)``), ``asyfp``,
            ``htpl``, ``wing_cla`` (``WING(101)``), ``tail_cla``
            (``HT(101)``), ``tra70``, ``trah70`` (the Mach 0.6 slopes),
            ``delcl6`` (``WING(201..)``), ``cfact`` (``FLP(71..)``),
            ``clrlm6`` (``HT(211..)``), ``cnym6`` (``BODY(201..)``),
            ``bd`` (GDELTA's ``TCD(43..46)``), ``sspn``, ``sspne``
            (``HTIN(4)``, ``(3)``), ``astrw`` (``A(7)``), ``depsda``
            (``DWASH(41)``), ``blref``.

    Returns:
        The words set: ``delcl`` (``WING(201..)``), ``claldl``
        (``WING(241..)``), ``clafs`` (``HT(211..)``), ``cnafs``
        (``BODY(201..)``), ``clrolt`` (``HT(211..)``), ``transl``, and
        the ``TRN`` words ``encepe``, ``yh``, ``etaqrs``, ``cldelc``,
        ``cldalc``, ``kbh``, ``khb``.

    Notes:
        Kept as executed: Figure 6.2.1.2-22 is read with four of its six
        points (the call passes ``NP=4``), so above a span ratio of 0.3
        the tail's rolling-moment factor is extrapolated rather than read;
        the lift-curve slope with a translating flap is
        ``CLAW*(1+CFACT)+CLAW``, the basic slope counted twice; and the
        all-moving tail's subsonic rolling moment uses ``pi*A/RAD`` with
        the wing's aspect ratio.
    """
    flap_cfg = data['flap']
    ndelta = int(float(flap_cfg['ndelta']) + .5)
    ftype, stype = float(flap_cfg['type']), float(flap_cfg['stype'])
    transl = ftype not in (1.0, 5.0, 6.0)
    cla = float(data['wing_cla'])
    cla_ratio = cla / float(data['tra70'])
    result: Dict[str, object] = {'transl': transl}
    if not data['asyfp']:
        if data['htpl']:
            cla = float(data['tail_cla'])
            cla_ratio = cla / float(data['trah70'])
        result['delcl'] = [float(v) * cla_ratio for v in data['delcl6'][:ndelta]]
        if transl:
            result['claldl'] = [cla * (1. + float(c)) + cla
                                for c in data['cfact'][:ndelta]]
        else:
            result['claldl'] = [UNUSED] * ndelta
        return result
    if stype != 5.:
        result['clafs'] = [float(v) * cla_ratio for v in data['clrlm6'][:ndelta]]
        if stype == 4.:
            ncrol = ndelta * int(data['nalpha'])
            result['cnafs'] = [v if v == UNUSED else v * cla_ratio
                               for v in (float(x) for x in
                                         data['cnym6'][:ncrol])]
        return result
    sspn, sspne = float(data['sspn']), float(data['sspne'])
    exposed_span = sspn - sspne
    tcd = [float(v) for v in data['bd']]
    encepe = ((.352 * tcd[0] + .503 * tcd[1] + .344 * tcd[2] + .041 * tcd[3]) /
              (.383 * tcd[0] + .707 * tcd[1] + .924 * tcd[2] + .5 * tcd[3]))
    delta_left = [float(v) for v in flap_cfg['deltal'][:ndelta]]
    delta_right = [float(v) for v in flap_cfg['deltar'][:ndelta]]
    cnah = float(data['tail_cla'])
    blref = float(data['blref'])
    result['encepe'] = encepe
    if float(data['mach']) < 1.:
        yh = encepe * sspne + exposed_span
        # The call passes NP=4: the table's last two points are unread.
        etaqrs = tbfunx(_X12222[:4], _Y21222[:4], exposed_span / sspn, 1, 2)[0]
        cldelc = (.5 * (1. - (PI * float(data['astrw']) / RAD) *
                        float(data['depsda'])) * etaqrs * yh * cnah) / blref
        result.update({'yh': yh, 'etaqrs': etaqrs, 'cldelc': cldelc,
                       'clrolt': [cldelc * (left - right)
                                  for left, right in zip(delta_left, delta_right)]})
    else:
        yh = .4 * sspne + exposed_span
        span_ratio = (sspn - sspne) / (2. * sspn)
        khb = tbfunx(_X12A1, _Y12A1, span_ratio, 0, 0)[0]
        kbh = tbfunx(_X12A2, _Y12A2, span_ratio, 0, 0)[0]
        cldalc = .35 * (khb + kbh) * cnah * yh / blref
        result.update({'yh': yh, 'khb': khb, 'kbh': kbh, 'cldalc': cldalc,
                       'clrolt': [cldalc * (left - right)
                                  for left, right in zip(delta_left, delta_right)]})
    return result


def m40o50_words(asyfp: bool) -> Dict[int, float]:
    """M40O50 after TRNYRL: without an asymmetric flap, ``WING(221..230)``
    are set UNUSED and ``WING(252..260)`` -UNUSED."""
    if asyfp:
        return {}
    out = {220 + j: UNUSED for j in range(1, 11)}
    out.update({250 + j: -UNUSED for j in range(2, 11)})
    return out
