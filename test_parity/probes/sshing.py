"""
Probe SSHING (with DFLCON and PTCP), the supersonic hinge moments, and
save the fixture.  SSSYM, which SSHING calls last, is stubbed: it has its
own tests.

The cases run in one program, in order: SSHING's locals are saved between
calls, so the replay carries them.

Run from the repository root: ``python test_parity/probes/sshing.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['sshing', 'dflcon', 'ptcp', 'arccos', 'arcsin', 'tbfunx',
            'quad']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE SSSYM
      RETURN
      END
"""

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /FLAPIN/ F(69)
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(77)
      COMMON /HTI/    HTIN(131)
      COMMON /IWING/  PWING, WING(400)
      COMMON /POWR/   SPR(59)
      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF
      COMMON /IHT/    PHT, HT(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1       HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,
     2       HYPERS,SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,
     2        HYPERS,SYMFP,ASYFP,TRIMC,TRIM
"""


def surface(tanle=0.3, tante=-0.1, bstro2=0.0, scale=1.0):
    return ({1: 3.0 * scale, 2: bstro2, 3: 9.0 * scale, 4: 10.0 * scale,
             5: 5.0 * scale, 6: 8.0 * scale, 18: 0.04, 66: 0.035},
            {10: 7.5 * scale, 25: 0.375, 28: 0.5, 59: math.atan(tanle),
             62: tanle, 77: math.atan(tante), 80: tante, 86: 0.5,
             104: -0.05})


def case(mach=2.0, ftype=1.0, htpl=False, aloci=3.0, aloco=8.0, cfi=2.0,
         cfo=1.5, tanle=0.3, tante=-0.1, bstro2=0.0, thick=0.02):
    win, a = surface(tanle, tante, bstro2)
    htin, aht = surface(0.45, -0.05, 0.0, 0.5)
    f = {11: thick, 12: cfi, 13: cfo, 14: aloci, 15: aloco, 16: 3.0,
         17: ftype}
    for i, (ci, co) in enumerate(((2.3, 1.8), (2.5, 1.9), (2.8, 2.1))):
        f[39 + i], f[49 + i] = ci, co
    return {'mach': mach, 'htpl': htpl, 'win': win, 'a': a, 'htin': htin,
            'aht': aht, 'f': f, 'claw': 0.052, 'clah': 0.047,
            'sref': 150.0,
            'spr': [0.013 * (k + 1) for k in range(59)]}


def cases():
    return [
        case(),                                     # inboard, tapered
        case(aloco=10.0),                           # tip control
        case(cfo=2.0),                              # untapered
        case(tanle=0.0),                            # unswept leading edge
        case(mach=1.05),                            # subsonic edge: RETURN
        case(ftype=2.0),                            # translating
        case(ftype=3.0, bstro2=4.0, aloci=7.0),     # outboard panel
        case(htpl=True, aloci=1.0, aloco=3.0, cfi=1.0, cfo=0.8),
        case(ftype=5.0),                            # stale geometry
        case(mach=3.0, aloci=0.5, aloco=8.8),
        case(mach=1.4, aloci=4.0, aloco=6.0, cfi=3.0, cfo=2.8),
        case(mach=2.5, aloci=1.0, aloco=9.9, cfi=1.0, cfo=0.5),
        case(mach=1.6, aloci=2.0, aloco=3.0, cfi=4.0, cfo=3.5),
        # Found by sampling: region cases 4 and 5, and a subsonic edge.
        case(mach=1.37, aloci=4.81, aloco=8.14, cfi=4.27, cfo=4.14),
        case(mach=3.04, aloci=2.16, aloco=6.29, cfi=3.23, cfo=1.6,
             tanle=0.1),
        case(mach=1.24, aloci=4.14, aloco=5.22, cfi=3.57, cfo=2.0,
             tante=0.2),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      IM=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('WING(101)', c['claw']))
        lines.append(assign('HT(101)', c['clah']))
        for name, key in (('WINGIN', 'win'), ('A', 'a'), ('HTIN', 'htin'),
                          ('AHT', 'aht'), ('F', 'f')):
            for k, v in c[key].items():
                lines.append(assign(f'{name}({k})', v))
        for k, v in enumerate(c['spr']):
            lines.append(assign(f'SPR({k + 1})', v))
        for k in list(range(241, 244)) + [251, 261]:
            lines.append(f'      WING({k})=-7.')
        lines.append('      CALL SSHING')
        lines.append("      WRITE(6,'(A,59ES25.16)') 'SPR',SPR")
        lines.append("      WRITE(6,'(A,5ES25.16)') 'W',")
        lines.append('     1(WING(K),K=241,243),WING(251),WING(261)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('sshing', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('win', 'a', 'htin', 'aht', 'f'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('sshing', payload))


if __name__ == '__main__':
    main()
