"""
Probe TRANWG and CLMXB1 and save the fixture.

Run from the repository root: ``python test_parity/probes/tranwg.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['tranwg', 'clmxb1', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'quad']


def case(ar=3.0, sweep=40.0, taper=0.3, a25=None, deltay=1.5, bu4=2.0,
         a160=3.0, xovc=0.3):
    tan = math.tan(math.radians(sweep))
    return {'a': {3: 120.0, 7: ar, 25: taper if a25 is None else a25,
                  27: taper, 58: sweep, 61: math.cos(math.radians(sweep)),
                  62: tan, 86: tan * 0.8},
            'deltay': deltay, 'sref': 150.0,
            'clmxb1': {'bu4': bu4, 'a160': a160, 'xovc': xovc}}


def cases():
    return [
        case(),                                   # subsonic LE at 1.3-1.5
        case(sweep=20.0),                         # supersonic LE
        case(sweep=40.0, ar=6.0, deltay=3.0),     # mixed along the anchors
        case(sweep=0.0, taper=1.0, ar=0.6),       # rectangular, A*beta <= 1
        case(sweep=0.0, taper=1.0, ar=4.0),       # rectangular, A*beta > 1
        case(sweep=0.0, taper=0.999, ar=4.0),     # nearly rectangular: 56A
        case(sweep=60.0, ar=1.5, taper=0.0, deltay=0.5, xovc=0.4, a160=6.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      COMMON /WINGD/  A(195)',
             '      COMMON /WINGI/  WINGIN(77)',
             '      COMMON /OPTION/ SW',
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        for index, value in c['a'].items():
            lines.append(assign(f'A({index})', value))
        lines.append(assign('WINGIN(17)', c['deltay']))
        lines.append(assign('WINGIN(18)', c['clmxb1']['xovc']))
        lines.append(assign('SW', c['sref']))
        lines.append(assign('A(160)', c['clmxb1']['a160']))
        lines.append('      CALL TRANWG(CNA,DCNA)')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'OUT',CNA,DCNA,A(62),A(86)")
        lines.append(assign('BU4', c['clmxb1']['bu4']))
        lines.append('      CALL CLMXB1(BU4,CLS,A,WINGIN)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'CLS',CLS")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('tranwg', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('tranwg', payload))


if __name__ == '__main__':
    main()
