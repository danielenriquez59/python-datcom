"""
Probe CLRDER, the rolling moment due to yaw rate, and save the fixture.

Run from the repository root: ``python test_parity/probes/clrder.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['clrder', 'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  WINGIN(77)
      COMMON /WBHCAL/ WBT(156)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108), STBH(135)
      COMMON /SUPWBB/  SWB(61),SHB(61)
      COMMON /IWING/  PWING,WING(400)
      COMMON /IVT/    PVT,VT(380)
      COMMON /IVF/    PVF,VF(380)
      COMMON /IBODY/  PBODY,BODY(400)
      COMMON /IBW/    PBW,BW(380)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IBH/    PBH, BH(380)
      COMMON /IBV/    PBV, BV(380)
      COMMON /IBWHV/  PBWHV,BWHV(380)
      COMMON /IBWV /  PBWV,BWV(380)
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
"""

UNUSED = 1.0e-30
FLAGS = ('bo', 'subson', 'transn', 'supers', 'wgpl', 'htpl', 'vtpl', 'vfpl')
# Each result block: its fill step and the words the routine touches.
BLOCKS = {'wing': (0.001, list(range(21, 41)) + list(range(361, 381))),
          'body': (0.002, list(range(201, 381))),
          'bw': (0.003, list(range(361, 381))),
          'bwhv': (0.004, list(range(361, 381))),
          'bwv': (0.005, list(range(361, 381))),
          'vt': (0.0005, [141] + list(range(201, 381))),
          'vf': (0.0007, [141] + list(range(201, 381))),
          'ht': (0.0011, list(range(201, 381))),
          'bh': (0.0013, list(range(201, 381))),
          'bv': (0.0017, list(range(201, 381)))}


def case(sweep=30.0, ar=6.0, taper=0.4, mach=0.5, nalpha=6, unused=(),
         **flags):
    f = {k: False for k in FLAGS}
    f.update(subson=True, bo=True, wgpl=True, htpl=True, vtpl=True,
             vfpl=True)
    f.update(flags)
    blocks = {n: {k: step * k for k in words}
              for n, (step, words) in BLOCKS.items()}
    for n, k in unused:
        blocks[n][k] = UNUSED
    return dict(f, nalpha=nalpha, mach=mach, blref=30.0,
                alpha=[-2.0 + 3.0 * k for k in range(nalpha)],
                a={7: ar, 27: taper, 41: math.radians(sweep)},
                wingin={11: -2.0}, stb={11: 18.0, 12: 3.0, 122: 4.0},
                stbh={11: 14.0, 12: -2.0},
                wbt={1: 1.1, 2: 0.2, 65: 0.15, 66: 1.05},
                shb={11: 0.18, 35: 1.08}, unused=list(unused), **blocks)


def cases():
    return [
        case(),
        case(sweep=0.0, vfpl=False),
        case(sweep=62.0, ar=12.0, taper=1.5, vtpl=False),
        case(sweep=-10.0),                               # forward swept
        case(wgpl=False),                                # stale CLR
        case(subson=False, transn=True),
        case(subson=False, supers=True),
        case(subson=False),                              # stale AKHB
        case(unused=[('ht', 221), ('body', 301), ('vt', 241),
                     ('vf', 241), ('body', 341)]),
        case(unused=[('ht', 201)], mach=0.8, nalpha=20),
        case(bo=False, sweep=45.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,380',
                  *[f'      {name.upper()}(K)=K*{step}'
                    for name, (step, _) in BLOCKS.items()],
                  f'{100 + n:5d} CONTINUE']
        for name, k in c['unused']:
            lines.append(assign(f'{name.upper()}({k})', UNUSED))
        for flag in FLAGS:
            lines.append(f"      {flag.upper()}="
                         f".{'TRUE' if c[flag] else 'FALSE'}.")
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + k})', v))
        lines.append(assign('BLREF', c['blref']))
        for name in ('a', 'wingin', 'stb', 'stbh', 'wbt', 'shb'):
            for k, v in c[name].items():
                lines.append(assign(f'{name.upper()}({k})', v))
        lines.append('      CALL CLRDER')
        for name, (_, words) in BLOCKS.items():
            lo = 201 if len(words) > 40 else 361
            lines.append(f"      WRITE(6,'(A,180ES25.16)') '{name.upper()}',"
                         f"({name.upper()}(K),K={lo},380)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('clrder', driver(all_cases), ROUTINES))
    print(save('clrder', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
