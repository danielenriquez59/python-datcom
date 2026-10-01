"""
Probe overlay M04O04 (BODOPT and the body slope pass) and save the fixture.

Run from the repository root: ``python test_parity/probes/m04o04.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m04o04', 'bodopt', 'eqspce', 'interx', 'fig26', 'getmax',
            'tbfunx', 'quad', 'trapz', 'tlin1x', 'tlinex', 'tlin3x', 'glook',
            'switch']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

COMMONS = """\
      COMMON /IBODY/  PBODY, BODY(400)
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IJKDUM(4),NOVLY
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /WINGD/  A(195), B(49)
      COMMON /BODYI/  XNX,X(20),S(20),P(20),R(20),
     1                ZU(20),ZL(20),BNOSE,BTAIL,RLN,RLA,DS
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,XVF,
     1                SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /BDATA/  BD(762)
      COMMON /EXPER/  KKK(104),KBODY
      LOGICAL KBODY,VERTUP
"""

ALPHA = [-4., 0., 4., 8., 12., 16., 20., 28., 40.]


def body(kind):
    x = np.array([0., 2., 5., 9., 14., 20., 26., 31., 35., 38.])
    if kind == 'closed':
        half = 2.0 * np.sin(np.pi * x / 38.0)**0.8
        half[-1] = 0.0                 # closes exactly: S(NX) = 0
    else:
        half = np.where(x < 9.0, 2.0 * np.sqrt(x / 9.0), 2.0)
        half = np.where(x > 31.0, 2.0 - 0.1 * (x - 31.0), half)
    width = half * 1.2
    camber = 0.4 * np.sin(np.pi * x / 38.0) + 0.02 * x
    if kind == 'droop':
        camber = np.where(x < 9.0, -0.5 * (9.0 - x) / 9.0, 0.0) + 0.01 * x
    zu = camber + half
    zl = camber - half
    s = np.pi * width * half
    p = np.pi * (width + half)
    return x, s, p, width, zu, zl


def case(kind='open', mach=0.5, kbody=False):
    x, s, p, r, zu, zl = body(kind)
    return {'x': list(x), 's': list(s), 'p': list(p), 'r': list(r),
            'zu': list(zu), 'zl': list(zl), 'alpha': ALPHA, 'mach': mach,
            'rn': 2.0e6, 'sref': 60.0, 'cbar': 5.0, 'blref': 25.0,
            'xcg': 17.0, 'aliw': 1.5, 'roughness': 1.6e-4, 'kbody': kbody}


def cases():
    return [case(), case(kind='closed'), case(kind='droop', mach=0.7),
            case(kbody=True)]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      M=1']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.400) BODY(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      KBODY=.{'TRUE' if c['kbody'] else 'FALSE'}.")
        lines.append(f"      NALPHA={len(c['alpha'])}")
        lines.append(f"      FLC(2)={len(c['alpha'])}.")
        lines.append(f"      XNX={len(c['x'])}.")
        for k in range(len(c['x'])):
            for name in ('x', 's', 'p', 'r', 'zu', 'zl'):
                lines.append(assign(f'{name.upper()}({k + 1})', c[name][k]))
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('B(1)', c['mach']))
        lines.append(assign('FLC(43)', c['rn']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbar']))
        lines.append(assign('BLREF', c['blref']))
        lines.append(assign('ROUGFC', c['roughness']))
        lines.append(assign('XCG', c['xcg']))
        lines.append(assign('ALIW', c['aliw']))
        lines.append('      CALL M04O04')
        for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                           ('CA', 80), ('CLA', 100), ('CMA', 120),
                           ('CYB', 140), ('CNB', 160), ('CLB', 180)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"(BODY({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,4ES25.16)') 'BD',BD(81),BD(62),")
        lines.append('     1  BD(61),BD(10)')
        lines.append("      WRITE(6,'(A,30ES25.16)') 'CDL',"
                     "(BD(214+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,30ES25.16)') 'BLOC',"
                     "(B(22+J),J=1,NALPHA)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('m04o04', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('m04o04', payload))


if __name__ == '__main__':
    main()
