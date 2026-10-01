"""
Probe overlay M06O06 (BODYRT, BODYJM, the slope pass) and save the fixture.

This is also the first execution check of the BODYRT translation.

Run from the repository root: ``python test_parity/probes/m06o06.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m06o06', 'bodyrt', 'bodyjm', 'eqspc1', 'interx', 'fig26',
            'getmax', 'tbfunx', 'quad', 'trapz', 'tlin1x', 'tlinex',
            'tlin3x', 'glook', 'switch']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

COMMONS = """\
      COMMON /IBODY/ PBODY, BODY(400)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IJKDUM(4),NOVLY
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /FLGTCD/ FLC(160)
      COMMON /WINGD/  A(195), B(49)
      COMMON /BODYI/  XNX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1                BNOSE,BTAIL,BLN,BLA,DS,ITYPE,METHOD,ELLIP
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /BDATA/  BD(762)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      COMMON /EXPER/  KKK(104),KBODY
      REAL ITYPE, METHOD
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB,KBODY,VERTUP
"""

UNUSED = 1.0e-30
ALPHA = [-6., -2., 0., 2., 4., 8., 12., 16., 20.]


def body_shape(kind):
    if kind == 'ogive':
        x = np.array([0., 1., 2.5, 4., 6., 9., 14., 20., 26., 30.])
        rmax = 1.5
        r = np.where(x < 6.0, rmax * np.sqrt(np.clip(x / 6.0, 0, 1)), rmax)
    else:
        x = np.array([0., 2., 4., 7., 11., 16., 22., 28., 33., 36.])
        rmax = 2.0
        r = rmax * np.sin(np.pi * np.clip(x / 40.0, 0, 1) + 0.15) ** 0.7
        r[0] = 0.0
    return list(x), list(r)


def case(kind='ogive', method=1.0, ellip=UNUSED, alpha_zero=0.0,
         kbody=False, mach=0.5):
    x, r = body_shape(kind)
    r = np.array(r)
    return {
        'x': x, 's': list(np.pi * r**2), 'p': list(2 * np.pi * r),
        'r': list(r), 'alpha': ALPHA, 'alpha_zero': alpha_zero,
        'mach': mach, 'rn': 2.0e6, 'sref': 40.0, 'cbar': 4.0,
        'blref': 20.0, 'bd33': 14.0, 'xcg': 13.5, 'roughness': 1.6e-4,
        'method': method, 'ellip': ellip, 'kbody': kbody,
    }


def cases():
    return [
        case(),
        case(kind='boattail', alpha_zero=-1.0, mach=0.7),
        case(method=2.0),
        case(method=2.0, ellip=1.0, kind='boattail'),
        case(method=2.0, ellip=0.6),
        case(method=2.0, ellip=1.8, kind='boattail', mach=0.3),
        case(kbody=True),
        case(method=2.0, kbody=True),
        case(kind='boattail', mach=0.7),
        case(alpha_zero=-1.0),
        case(kind='boattail', alpha_zero=-1.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      M=1',
             '      SUBSON=.TRUE.', '      TRANSN=.FALSE.']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.400) BODY(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      KBODY=.{'TRUE' if c['kbody'] else 'FALSE'}.")
        lines.append(f"      NALPHA={len(c['alpha'])}")
        lines.append(f"      XNX={len(c['x'])}.")
        for k in range(len(c['x'])):
            for name, key in [('X', 'x'), ('S', 's'), ('P', 'p'), ('R', 'r')]:
                lines.append(assign(f'{name}({k + 1})', c[key][k]))
        lines.append('      ZU(1)=UNUSED')
        lines.append('      ZL(1)=UNUSED')
        lines.append(assign('METHOD', c['method']))
        lines.append(assign('ELLIP', c['ellip']))
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'BD({255 + j})', a + c['alpha_zero']))
        lines.append(assign('BD(81)', c['alpha_zero']))
        lines.append(f"      FLC(2)={len(c['alpha'])}.")
        lines.append(assign('BD(33)', c['bd33']))
        lines.append(assign('XCG', c['xcg']))
        lines.append(assign('B(1)', c['mach']))
        lines.append(assign('FLC(43)', c['rn']))
        for name in ('sref', 'cbar', 'blref'):
            target = {'cbar': 'CBARR'}.get(name, name.upper())
            lines.append(assign(target, c[name]))
        lines.append(assign('ROUGFC', c['roughness']))
        lines.append('      CALL M06O06')
        for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                           ('CA', 80), ('CLA', 100), ('CMA', 120),
                           ('CYB', 140), ('CNB', 160), ('CLB', 180)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"(BODY({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,6ES25.16)') 'BD',BD(57),BD(61),")
        lines.append('     1  BD(62),BD(76),BD(34),BD(35)')
        lines.append("      WRITE(6,'(A,2ES25.16)') 'MORE',BD(275),ELLIP")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('m06o06', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('m06o06', payload))


if __name__ == '__main__':
    main()
