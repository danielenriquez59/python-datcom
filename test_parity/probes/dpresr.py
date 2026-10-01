"""
Probe DPRESR, the supersonic dynamic pressure at the tail, and save the
fixture.

Run from the repository root: ``python test_parity/probes/dpresr.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['dpresr', 'fig68', 'arcsin', 'arccos', 'mach2', 'tbfunx',
            'quad']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /SYNTSS/ SYNA(19)
      COMMON /SUPDW/  DWA(237)
      COMMON /WINGI/  WINGIN(100)
      COMMON /WINGD/  A(195)
"""

CONVEX = [8.0, 5.0, 2.0, -1.0, -4.0, -7.0]
WAVY = [3.0, 4.0, 1.0, 2.0, -2.0, -5.0]
FLAT = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
THIN = [2.0, 1.2, 0.4, -0.4, -1.2, -2.0]


def case(zj, alpha, slope, mach=2.0, zwake=0.3, dwangl=0.02, cr=8.0,
         a12=1.0, rl2=12.0):
    return {'zj': zj, 'zwake': zwake, 'alpha': alpha, 'dwangl': dwangl,
            'mach': mach, 'cr': cr, 'a12': a12, 'rl2': rl2,
            'slope': slope, 'mj_in': -3.0, 'qq_in': -7.0}


def cases():
    return [
        case(2.0, 2.0, CONVEX),               # LE shock, then expansions
        case(2.0, 10.0, CONVEX),              # LE expansion above
        case(-2.0, 5.0, CONVEX),              # LE shock below
        case(-2.0, -12.0, CONVEX),            # LE expansion below
        case(2.0, 1.0, WAVY),                 # shocks mid-chord
        case(2.0, 1.0, WAVY, mach=1.3),
        case(-2.0, 12.0, CONVEX, mach=1.25),  # detached, subsonic behind
        case(0.05, 4.0, FLAT),                # inside the last wave
        case(-0.05, 4.0, FLAT, mach=3.0),
        case(3.0, 0.0, THIN, mach=1.6),
        case(-3.0, 0.0, THIN, mach=1.6),
        case(1.5, 6.0, THIN, mach=2.5, dwangl=0.05),
        case(-1.5, 3.0, WAVY, mach=4.0, dwangl=-0.03),
        case(0.8, 15.0, CONVEX, mach=1.8),
        case(12.0, 2.0, CONVEX),              # past the waves, Z falling
        case(-14.0, 5.0, CONVEX),             # past the waves, Z rising
        case(-2.0, 20.0, CONVEX, mach=1.5),   # TE shock detached
        case(2.0, 5.0, CONVEX, mach=1.05),    # subsonic behind LE shock
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=231,236')
        lines.append(f' {3000 + n} DWA(K)=0.')
        lines.append(assign('DWA(1)', c['mach']))
        lines.append(assign('WINGIN(6)', c['cr']))
        lines.append(assign('A(12)', c['a12']))
        lines.append(assign('A(24)', c['rl2']))
        for k, s in enumerate(c['slope']):
            lines.append(assign(f'WINGIN({95 + k})', s))
        for name, key in [('ZJ', 'zj'), ('ZWAKE', 'zwake'),
                          ('ALPHJ', 'alpha'), ('DWANGL', 'dwangl'),
                          ('QQ', 'qq_in'), ('XMJ', 'mj_in')]:
            lines.append(assign(name, c[key]))
        lines.append('      CALL DPRESR(ZJ,ZWAKE,ALPHJ,DWANGL,QQ,XMJ)')
        lines.append("      WRITE(6,'(A,8ES25.16)') 'R',QQ,XMJ,")
        lines.append('     1(DWA(K),K=231,236)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('dpresr', driver(all_cases), ROUTINES))
    print(save('dpresr', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
