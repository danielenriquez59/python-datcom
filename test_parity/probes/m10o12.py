"""
Probe overlay M10O12 (WGEOTL, WBTAIL and its closing pass) and save the
fixture.

Run from the repository root: ``python test_parity/probes/m10o12.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m10o12', 'wgeotl', 'wbtail', 'bodowg', 'getmax', 'ali',
            'tbfunx', 'quad', 'tlin3x', 'tlinex', 'tlin1x', 'glook', 'switch']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

COMMONS = """\
      COMMON /IVT/    PVT,    VT(380)
      COMMON /IVF/    PVF,    VF(380)
      COMMON /IBW/    PBW,    BW(380)
      COMMON /IBWH/   PBWH,   BWH(380)
      COMMON /IBWV/   PBWV,   BWV(380)
      COMMON /IBWHV/  PBWHV,  BWHV(380)
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /FLGTCD/ FLC(160)
      COMMON /WINGI/  WINGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,K,NOVLY
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      COMMON /WINGD/  A(195),B(49)
      COMMON /WHWB/   FACT(182),WB(39),HB(39)
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /BDATA/  BD(762)
      COMMON /WBHCAL/ WBT(155)
      COMMON /WHAERO/ C(51),D(55),CHT(51),DHT(55),DVT(55),DVF(55)
      COMMON /IWING/  PWING,WING(400)
      COMMON /VTDATA/ AVT(195)
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /IHT/    PHT,HT(380)
      COMMON /BODYI/ XNX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1               BNOSE,BTAIL,BLN,BLA,DS
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB,VERTUP
"""

UNUSED = 1.0e-30


def case(rng, alpha=None, canard=False, htpl=True, missing_cm=None,
         missing_cd=None, alih=0.0, q_loss=0.1):
    alpha = np.array(alpha if alpha is not None else
                     [-4., 0., 4., 8., 12., 16., 20.])
    n = len(alpha)
    aliw = 1.0
    tail_local = alpha + alih
    ht_zero = -0.5
    ht_cl = 0.05 * (tail_local - ht_zero)
    eps = 0.35 * (alpha + 1.5)
    deda = 0.35 + 0.01 * np.cos(alpha / 7.0)
    q = 1.0 - q_loss * np.exp(-(alpha / 8.0)**2)
    bw_cl = 0.075 * (alpha + 2.0) - 0.0006 * np.maximum(alpha - 12, 0)**2
    bw_cm = -0.02 - 0.012 * (alpha + 2.0)
    if missing_cm is not None:
        bw_cm[missing_cm:] = 2 * UNUSED
    bw_cd = 0.02 + 0.06 * bw_cl**2
    if missing_cd is not None:
        bw_cd[missing_cd] = -UNUSED
    x = np.linspace(0.0, 45.0, 12)
    radius = 2.2 * np.sin(np.pi * np.clip(x / 50.0, 0, 1))**0.5
    wing_sspn, wing_sspne = 15.0, 13.0
    return {
        'alpha_deg': list(alpha),
        'htpl': htpl,
        'twash': 3.0 if canard else 0.0,
        'tail': {'sspn': 6.0, 'sspne': 5.0, 'a3': 55.0, 'a7': 4.1,
                 'a10': 3.2, 'a27': 0.5, 'a62': 0.3, 'a161': 1.4,
                 'c6': rng.uniform(0.2, 0.3), 'local_alpha': list(tail_local),
                 'cd0': 0.006, 'cm0': -0.004, 'alpha_zero_lift': ht_zero},
        'tail_alone': {'cd': list(0.006 + 0.04 * ht_cl**2),
                       'cl': list(ht_cl),
                       'cla': list(0.05 + 0.001 * np.cos(alpha / 5.0))},
        'wing_body': {'cd': list(bw_cd), 'cl': list(bw_cl), 'cm': list(bw_cm),
                      'cla': list(0.075 + 0.001 * np.sin(alpha / 9.0)),
                      'cma': list(-0.012 + 0.0005 * np.sin(alpha / 6.0)),
                      'cd0': 0.019, 'kwb': 1.12},
        'downwash': {'qoqi': list(q), 'angle': list(eps),
                     'gradient': list(deda),
                     'vortex_span': list(26.0 - 0.1 * alpha)},
        'wing': {'sspn': wing_sspn, 'sspne': wing_sspne, 'cla': 0.07,
                 'local_alpha': list(alpha + aliw), 'alpha_zero_lift': -1.8,
                 'a120': 6.2, 'a38': rng.uniform(0.1, 0.6), 'a118': 0.45,
                 'a12': 2.0, 'a24': 22.0, 'a80': -0.08, 'beta': 0.87},
        'synthesis': {'xcg': 21.0, 'xh': 38.0, 'zh': 2.0, 'zcg': 0.0,
                      'alih': alih, 'aliw': aliw},
        'body': {'x': list(x), 's': list(np.pi * radius**2), 'bd70': 1.3},
        'sref': 150.0, 'cbarr': 5.2,
        'vertical': {'dvt20': 0.003, 'dvf20': 0.0005, 'vt1': 0.0031,
                     'vf1': 0.0006},
    }


def cases():
    rng = np.random.default_rng(1012)
    return [
        case(rng),
        case(rng, alih=-2.0, q_loss=0.2),
        case(rng, alpha=[-8., -2., 3., 6., 9., 14.]),
        case(rng, missing_cm=4),
        case(rng, missing_cd=2),
        case(rng, canard=True),
        case(rng, canard=True, alih=1.5),
        case(rng, htpl=False),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        t, ta, bw, dw, w, sy, bd, v = (c['tail'], c['tail_alone'],
                                       c['wing_body'], c['downwash'],
                                       c['wing'], c['synthesis'], c['body'],
                                       c['vertical'])
        alpha = c['alpha_deg']
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.380) BWH(K)=0.')
        lines.append('         IF(K.LE.380) BWHV(K)=0.')
        lines.append('         IF(K.LE.380) BWV(K)=0.')
        lines.append('         IF(K.LE.182) FACT(K)=0.')
        lines.append('         IF(K.LE.155) WBT(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(f'      NALPHA={len(alpha)}')
        for j, a in enumerate(alpha):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'BHT({23 + j})', t['local_alpha'][j]))
            lines.append(assign(f'B({23 + j})', w['local_alpha'][j]))
            for base, arr in [(0, ta['cd']), (20, ta['cl']),
                              (100, ta['cla'])]:
                lines.append(assign(f'HT({base + 1 + j})', arr[j]))
            for base, key in [(0, 'cd'), (20, 'cl'), (40, 'cm'), (100, 'cla'),
                              (120, 'cma')]:
                lines.append(assign(f'BW({base + 1 + j})', bw[key][j]))
            for base, key in [(0, 'qoqi'), (20, 'angle'), (40, 'gradient')]:
                lines.append(assign(f'DWASH({base + 1 + j})', dw[key][j]))
            lines.append(assign(f'FACT({82 + j})', dw['vortex_span'][j]))
        lines.append(assign('WINGIN(3)', w['sspne']))
        lines.append(assign('WINGIN(4)', w['sspn']))
        lines.append(assign('WINGIN(101)', c['twash']))
        lines.append(assign('HTIN(3)', t['sspne']))
        lines.append(assign('HTIN(4)', t['sspn']))
        for key, index in [('a3', 3), ('a7', 7), ('a10', 10), ('a27', 27),
                           ('a62', 62), ('a161', 161)]:
            lines.append(assign(f'AHT({index})', t[key]))
        lines.append(assign('CHT(6)', t['c6']))
        lines.append(assign('BHT(46)', t['cd0']))
        lines.append(assign('BHT(47)', t['cm0']))
        lines.append(assign('BHT(49)', t['alpha_zero_lift']))
        lines.append(assign('WB(2)', bw['kwb']))
        lines.append(assign('WB(17)', bw['cd0']))
        lines.append(assign('WING(101)', w['cla']))
        lines.append(assign('B(49)', w['alpha_zero_lift']))
        lines.append(assign('B(2)', w['beta']))
        for key, index in [('a120', 120), ('a38', 38), ('a118', 118),
                           ('a12', 12), ('a24', 24), ('a80', 80)]:
            lines.append(assign(f'A({index})', w[key]))
        for key, name in [('xcg', 'XCG'), ('xh', 'XH'), ('zh', 'ZH'),
                          ('zcg', 'ZCG'), ('alih', 'ALIH'), ('aliw', 'ALIW')]:
            lines.append(assign(name, sy[key]))
        lines.append(f"      XNX={len(bd['x'])}.")
        for k, (xv, sv) in enumerate(zip(bd['x'], bd['s'])):
            lines.append(assign(f'X({k + 1})', xv))
            lines.append(assign(f'S({k + 1})', sv))
        lines.append(assign('BD(70)', bd['bd70']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append(assign('DVT(20)', v['dvt20']))
        lines.append(assign('DVF(20)', v['dvf20']))
        lines.append(assign('VT(1)', v['vt1']))
        lines.append(assign('VF(1)', v['vf1']))
        lines.append('      CALL M10O12')
        for block in ('BWH', 'BWHV', 'BWV'):
            for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                               ('CA', 80), ('CLA', 100), ('CMA', 120)]:
                lines.append(f"      WRITE(6,'(A,30ES25.16)') '{block}{tag}',"
                             f"({block}({start}+J),J=1,NALPHA)")
        for tag, start in [('FACT42', 41), ('FACT82', 81), ('FACT102', 101),
                           ('FACT122', 121)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"(FACT({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,30ES25.16)') 'BD95',"
                     "(BD(94+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,155ES25.16)') 'WBT',(WBT(J),J=1,155)")
        lines.append("      WRITE(6,'(A,9ES25.16)') 'GEOM',BD(8),BD(30),BD(31),")
        lines.append("     1  BD(58),BD(63),BD(64),BD(84),BD(761),BD(762)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('m10o12', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('m10o12', payload))


if __name__ == '__main__':
    main()
