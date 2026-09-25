"""
Control-tab hinge moments: CTABS, and the flap overlay M36O44 that runs it.

CTABS forms, for a control surface with a trim or control tab and a
spring in the control circuit, the hinge moments with the tab free and
locked, the gearing increment, the stick (control) force and, for the two
tab types that deflect, the tab deflection, at each angle of attack and
control deflection.

The result arrays are ``(20, 9)`` FORTRAN arrays EQUIVALENCEd into the
``BW``, ``BH``, ``BV``, ``BWH`` and ``BWHV`` blocks from word 201, so entry
``(I, J)`` is word ``200 + I + 20*(J-1)``.

Reference: datcom-legacy/datcom_2000/ctabs.f, m36o44.f
"""

from typing import Callable, Dict, Mapping

from pydatcom.utils.constants import UNUSED


def _words(block):
    return {int(k): float(v) for k, v in block.items()}


def calculate_ctabs(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate CTABS: control-tab hinge moments and control force.

    Args:
        data: By name: ``nalpha``, ``mach`` and ``pinf`` (``FLC`` Mach and
            free-stream pressure at this Mach index), ``alpha`` (degrees),
            ``sref``, ``cbarr``, ``f`` (the ``/FLAPIN/`` words 1-137 as a
            1-based ``{word: value}`` map: the deflections 1-9, the count
            16, and ``TTYPE``, the tab chords and spans, ``B1``-``B4``,
            ``D1``-``D3``, ``GCMAX``, ``KS``, ``RL``, ``BGR``, ``DELR`` at
            117-137), and the result blocks ``bw``, ``bh``, ``bv``,
            ``bwh``, ``bwhv`` as ``{word: value}`` maps.

    Returns:
        The five result blocks.

    Notes:
        ``CFC`` and ``FC`` share ``BW`` 201 on, so with a known free-stream
        pressure the control force coefficient overwrites the control
        force, as the source does.
    """
    f = _words(data['f'])
    blk = {n: _words(data[n]) for n in ('bw', 'bh', 'bv', 'bwh', 'bwhv')}
    sref, cbarr = float(data['sref']), float(data['cbarr'])
    alpha = [0.0] + [float(v) for v in data['alpha']]
    ttype, cfitc, cfotc, bitc, botc = (
        f[word] for word in range(117, 122))
    b1, b2, b3, b4, d1, d2, d3 = (f[word] for word in range(126, 133))
    gcmax, ks, rl, bgr, delr = (f[word] for word in range(133, 138))
    stc = (cfitc + cfotc) * (botc - bitc)
    ctc = 2.0 * (cfitc + cfotc - cfitc * cfotc / (cfitc + cfotc)) / 3.0
    ac = stc * ctc / (sref * cbarr)
    q = kq = 0.0
    pinf, mach = float(data['pinf']), float(data['mach'])
    if pinf != UNUSED:
        q = 0.7 * pinf * mach ** 2
        kq = ks / q
    if rl < 0.0:
        r1, r2 = 0.0, 1.0
    if rl == 0.0:
        d = b2 / (ac * d2) + kq * bgr / d2
        r1 = delr / d
        r2 = -r1 * (kq / d2)
    if rl > 0.0:
        d = rl + b2 / (ac * d2) - kq * (rl - bgr) / d2
        r1 = (rl + delr) / d
        r2 = -r1 * (kq / d2)
    itype = int(ttype + 0.5)
    for angle_index in range(1, int(data['nalpha']) + 1):
        ndelta = int(f[16] + 0.5)
        for deflection_index in range(1, ndelta + 1):
            word = 200 + angle_index + 20 * (deflection_index - 1)
            delt = f[deflection_index]
            chctf = ((b1 + d1 * b2 / d2) * delt +
                     (b3 - d3 * b2 / d2) * alpha[angle_index])
            chctl = b1 * delt + b3 * alpha[angle_index]
            dchcg = ((bgr * b2 + bgr * ac * d1 + bgr * bgr * ac * d2) * delt
                     + bgr * ac * d3 * alpha[angle_index])
            cfc = gcmax * (r1 * chctf + r2 * chctl + r2 * dchcg)
            blk['bh'][word], blk['bv'][word], blk['bwh'][word] = (chctf,
                                                                   chctl,
                                                                   dchcg)
            blk['bw'][word] = cfc
            if itype == 2:
                blk['bwhv'][word] = -(b1 * delt + b3 * alpha[angle_index]) / b4
            if itype == 3:
                blk['bwhv'][word] = cfc / (b4 * (r1 + r2) * gcmax)
            if q > 0.0:
                blk['bw'][word] = cfc / (q * sref * cbarr)
    return blk


def m36o44(transn: bool, flap_type: float, ctab: bool,
           run_liftfp: Callable[[], object], run_hinge: Callable[[], object],
           run_ctabs: Callable[[], object]):
    """M36O44: LIFTFP; then, unless transonic, HINGE for flap types 1 and 5,
    ``WING(252..260)`` set to ``-UNUSED``, and CTABS when requested.

    Returns the results of the routines that ran, by name, and the
    ``WING`` words set (empty when transonic).
    """
    ran = {'LIFTFP': run_liftfp()}
    if transn:
        return ran, {}
    if flap_type == 1.0 or flap_type == 5.0:
        ran['HINGE'] = run_hinge()
    wing = {250 + word_index: -UNUSED for word_index in range(2, 11)}
    if ctab:
        ran['CTABS'] = run_ctabs()
    return ran, wing
