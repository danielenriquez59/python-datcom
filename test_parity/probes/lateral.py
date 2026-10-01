"""
Probe the subsonic lateral pass: M29O35, then M17O21 (SUBLAT twice).

Inputs are kept as COMMON slots keyed by one-based index, so the fixture
reads like the source.

Run from the repository root: ``python test_parity/probes/lateral.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m29o35', 'm17o21', 'sublat', 'trapz', 'getmax', 'tbfunx',
            'quad', 'tlinex', 'tlin1x', 'tlin3x', 'tlin4x', 'glook', 'switch']

UNUSED = 1.0e-30
TYPES = {1.0: 1, 2.0: 2, 3.0: 3, 4.0: 4}

COMMONS = """\
      COMMON /IBODY/   PBODY,  BODY(400)
      COMMON /IWING/   PWING,  WING(400)
      COMMON /IHT/     PHT,    HT(380)
      COMMON /IVT/     PVT,    VT(380)
      COMMON /IVF/     PVF,    VF(380)
      COMMON /IBW/     PBW,    BW(380)
      COMMON /IBH/     PBH,    BH(380)
      COMMON /IBV/     PBV,    BV(380)
      COMMON /IBWH/    PBWH,   BWH(380)
      COMMON /IBWV/    PBWV,   BWV(380)
      COMMON /IBWHV/   PBWHV,  BWHV(380)
      COMMON /WINGI/  WINGIN(101)
      COMMON /VTI/    VTIN(154), TVTIN(8), VFIN(154)
      COMMON /HTI/    HTIN(154)
      COMMON /WINGD/  A(195), B(49)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108), STBH(135)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /VTDATA/ AVT(195), AVF(195)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,K,NOVLY
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
      COMMON /SYNTSS/ SYNA(19)
      COMMON /BDATA/  BD(762)
      COMMON /BODYI/  BODYIN(126)
      COMMON /FLGTCD/ FLC(160)
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
"""

ALPHA = [-4., 0., 4., 8., 12., 16.]


def surface_blocks(rng, kind, ar, sweep_c2, sspn, sspne, chord, taper):
    """A surface's input and A/B blocks, loosely self-consistent."""
    c4 = sweep_c2 + 4.0
    win = {1: chord * taper, 2: 0.0, 3: sspne, 4: sspn, 5: chord * 0.7,
           6: chord, 11: -2.0, 12: UNUSED, 13: 3.0, 14: UNUSED,
           15: float(kind)}
    a = {1: sspne * chord * 0.8, 3: sspne * chord * 1.3,
         4: sspn * chord * 1.4, 10: chord * 0.9, 16: chord * 0.7,
         2: 5.0, 23: sspne * 0.4, 25: taper + 0.1, 30: chord * 0.4,
         42: math.sin(math.radians(c4)), 43: math.cos(math.radians(c4)),
         44: math.tan(math.radians(c4)),
         49: math.cos(math.radians(sweep_c2)),
         62: math.tan(math.radians(c4 + 5.0)),
         68: math.tan(math.radians(c4)), 70: sweep_c2,
         73: math.cos(math.radians(sweep_c2)), 86: 0.3,
         94: sweep_c2 * 0.6, 97: math.cos(math.radians(sweep_c2 * 0.6)),
         118: taper, 119: 12.0, 120: ar, 131: 0.105, 161: chord * 0.3,
         163: ar * 0.6, 167: sspne * chord * 0.5, 168: ar * 1.4,
         169: taper * 0.8, 171: 0.06, 172: 0.065,
         122: chord * 0.8, 136: 1.2, 195: 1.1}
    return win, a


def case(rng, kind=1.0, ar=6.0, nonuniform=False, zuzl=False, bo=True,
         htpl=True, vtpl=True, vfpl=True, tvtpan=False, phiv=UNUSED,
         phif=UNUSED, transn=False, vt_kind=1.0, sweep_c2=20.0, zh=2.0):
    n = len(ALPHA)
    alpha = np.array(ALPHA)
    mach = 0.6
    beta = math.sqrt(1 - mach**2)
    ww, wa = surface_blocks(rng, kind, ar, sweep_c2, 15.0, 13.0, 7.0, 0.4)
    hw, ha = surface_blocks(rng, 1.0, 4.0, 10.0, 6.0, 5.2, 3.0, 0.5)
    if nonuniform:
        ww.update({12: 4.0, 14: 7.0})
    vtin = {1: 2.0, 2: 1.5, 3: 0.2, 4: 5.5, 5: 3.0, 6: 4.5, 15: vt_kind,
            21: 0.1}
    vfin = {1: 1.0, 2: 0.0, 3: 0.1, 4: 1.2, 5: 1.4, 6: 1.6, 15: 1.0,
            21: 0.1}
    avt = {4: 22.0, 50: 0.4, 62: 0.6, 86: 0.5, 118: 0.45, 120: 1.6,
           122: 3.0, 136: 2.1, 195: 1.9}
    avf = {4: 3.0, 50: 0.2, 62: 0.4, 86: 0.3, 118: 0.6, 120: 1.1,
           122: 1.2, 136: 0.6, 195: 0.5}
    x = np.linspace(0.0, 45.0, 10)
    r = 2.2 * np.sin(np.pi * np.clip(x / 48.0, 0, 1))**0.6
    body = {'x': list(x), 'r': list(r)}
    if zuzl:
        body.update({'zu': list(r * 1.1 + 0.1), 'zl': list(-r * 0.9)})

    def curves(cla, a0):
        cl = cla * (alpha - a0)
        return {'cl': list(cl), 'cm': list(-0.02 - 0.1 * cl),
                'cn': list(cl / np.cos(np.radians(alpha)))}

    return {
        'alpha': ALPHA, 'mach': mach,
        'flags': {'wgpl': True, 'bo': bo, 'htpl': htpl, 'vtpl': vtpl,
                  'vfpl': vfpl, 'tvtpan': tvtpan, 'transn': transn,
                  'vertup': True},
        'wingin': ww, 'a': wa, 'b': {2: beta},
        'wing': curves(0.075, -2.0),
        'wing_cl0': list(0.07 * (alpha + 2.0)),
        'htin': hw, 'aht': ha, 'bht': {2: beta},
        'ht': curves(0.05, -0.5), 'ht_cl0': list(0.048 * (alpha + 0.5)),
        'bw_cl': list(0.08 * (alpha + 2.2)), 'bh_cl': list(0.004 * alpha),
        'body_lat': {'cla': 0.004, 'cla_h': 0.004, 'cyb': -0.003,
                     'cnb': 0.0012, 'clb': list(0.0001 * alpha)},
        'vtin': vtin, 'vfin': vfin, 'avt': avt, 'avf': avf,
        'tvtin': [1.2, 3.0, 0.8, 30.0, 9.0, 8.0, 18.0, 2.5],
        'syna': {1: 21.0, 2: 16.0, 3: -0.5, 4: 1.0, 5: 0.0, 6: 38.0,
                 7: zh, 8: 0.0, 9: 36.0, 12: 34.0, 18: phiv, 19: phif},
        'body': body, 'bd1': 45.0, 'bd66': 0.3, 'rn': 2.5e6,
        'sref': 170.0, 'cbarr': 6.0, 'blref': 30.0,
    }


def cases():
    rng = np.random.default_rng(5121)
    return [
        case(rng),
        case(rng, nonuniform=True, zuzl=True),
        case(rng, ar=0.8),
        case(rng, kind=3.0, sweep_c2=-15.0),
        case(rng, kind=2.0, ar=2.5, sweep_c2=40.0),
        case(rng, bo=False),
        case(rng, htpl=False, vfpl=False),
        case(rng, tvtpan=True, phiv=10.0),
        case(rng, phif=15.0),
        case(rng, transn=True),
        case(rng, vt_kind=3.0, zh=4.5),
        case(rng, vt_kind=3.0, zh=3.0),
        case(rng, zh=-1.0),
    ]


def block(lines, name, values):
    for index, value in values.items():
        if int(index) == 15 and name in ('WINGIN', 'HTIN', 'VTIN', 'VFIN'):
            lines.append(f'      {name}(15)=WT({TYPES[float(value)]})')
        else:
            lines.append(assign(f'{name}({int(index)})', value))


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      DIMENSION WT(4)', '      INTEGER*8 ITRUE',
             '      EQUIVALENCE (ITRUE,RTRUE)', '      ITRUE=1',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        f = c['flags']
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,400')
        for arr in ('BODY', 'WING'):
            lines.append(f'         {arr}(K)=0.')
        for arr in ('HT', 'VT', 'VF', 'BW', 'BH', 'BV', 'BWH', 'BWV', 'BWHV'):
            lines.append(f'         IF(K.LE.380) {arr}(K)=0.')
        for arr in ('A', 'AHT', 'AVT', 'AVF'):
            lines.append(f'         IF(K.LE.195) {arr}(K)=0.')
        for arr in ('VTIN', 'VFIN', 'HTIN'):
            lines.append(f'         IF(K.LE.154) {arr}(K)=0.')
        for arr in ('STB', 'STBH'):
            lines.append(f'         IF(K.LE.135) {arr}(K)=0.')
        lines.append('         IF(K.LE.126) BODYIN(K)=0.')
        lines.append('         IF(K.LE.101) WINGIN(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append('         IF(K.LE.49) BHT(K)=0.')
        lines.append(f' {label} CONTINUE')
        for flag in ('wgpl', 'bo', 'htpl', 'vtpl', 'vfpl', 'tvtpan',
                     'transn'):
            lines.append(f"      {flag.upper()}=.{'TRUE' if f[flag] else 'FALSE'}.")
        lines.append('      SYNA(10)=' + ('RTRUE' if f['vertup'] else '0.'))
        lines.append(f"      NALPHA={len(c['alpha'])}")
        block(lines, 'WINGIN', c['wingin'])
        block(lines, 'A', c['a'])
        block(lines, 'B', c['b'])
        block(lines, 'HTIN', c['htin'])
        block(lines, 'AHT', c['aht'])
        block(lines, 'BHT', c['bht'])
        block(lines, 'VTIN', c['vtin'])
        block(lines, 'VFIN', c['vfin'])
        block(lines, 'AVT', c['avt'])
        block(lines, 'AVF', c['avf'])
        block(lines, 'SYNA', c['syna'])
        for k, v in enumerate(c['tvtin']):
            lines.append(assign(f'TVTIN({k + 1})', v))
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
            for arr, key, base in [('WING', 'wing', 20), ('HT', 'ht', 20)]:
                for off, comp in [(0, 'cl'), (20, 'cm'), (40, 'cn')]:
                    lines.append(assign(f'{arr}({base + off + 1 + j})',
                                        c[key][comp][j]))
            lines.append(assign(f'B({3 + j})', c['wing_cl0'][j]))
            lines.append(assign(f'BHT({3 + j})', c['ht_cl0'][j]))
            lines.append(assign(f'BW({21 + j})', c['bw_cl'][j]))
            lines.append(assign(f'BH({21 + j})', c['bh_cl'][j]))
            lines.append(assign(f'BODY({181 + j})', c['body_lat']['clb'][j]))
        lines.append(assign('BODY(101)', c['body_lat']['cla']))
        lines.append(assign('BODY(141)', c['body_lat']['cyb']))
        lines.append(assign('BODY(161)', c['body_lat']['cnb']))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['rn']))
        lines.append(assign('BD(1)', c['bd1']))
        lines.append(assign('BD(66)', c['bd66']))
        body = c['body']
        lines.append(f"      BODYIN(1)={len(body['x'])}.")
        for k, xv in enumerate(body['x']):
            lines.append(assign(f'BODYIN({2 + k})', xv))
            lines.append(assign(f'BODYIN({62 + k})', body['r'][k]))
            if 'zu' in body:
                lines.append(assign(f'BODYIN({82 + k})', body['zu'][k]))
                lines.append(assign(f'BODYIN({102 + k})', body['zl'][k]))
        if 'zu' not in body:
            lines.append('      BODYIN(82)=UNUSED')
        for name in ('sref', 'cbarr', 'blref'):
            lines.append(assign(name.upper(), c[name]))
        lines.append('      CALL M29O35')
        lines.append('      CALL M17O21')
        for arr in ('WING', 'HT', 'VT', 'VF', 'BW', 'BH', 'BV', 'BWH', 'BWV',
                    'BWHV'):
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{arr}',")
            lines.append(f"     1  {arr}(141),{arr}(161),"
                         f"({arr}(180+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,135ES25.16)') 'STB',(STB(J),J=1,135)")
        lines.append("      WRITE(6,'(A,135ES25.16)') 'STBH',"
                     "(STBH(J),J=1,135)")
        lines.append("      WRITE(6,'(A,3ES25.16)') 'DIH',WINGIN(12),"
                     "WINGIN(13),WINGIN(14)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    # SUBLAT's /FLOLOG/ puts REAL XX(14) where the other routines have 14
    # LOGICALs; promoted to 8 bytes it would shift TRANSN.  Declared
    # LOGICAL, the layout is the original 4-byte one.
    patches = {'sublat': [
        ('      LOGICAL FLTC,OPTI,BO,WGPL,TVTPAN,TRANSN,HTPL,VTPL',
         '      LOGICAL FLTC,OPTI,BO,WGPL,TVTPAN,TRANSN,HTPL,VTPL,XX')]}
    records = parse_records(run('lateral', driver(all_cases), ROUTINES,
                                layout_patches=patches))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('lateral', payload))


if __name__ == '__main__':
    main()
