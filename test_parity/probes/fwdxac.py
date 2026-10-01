"""
Probe FWDXAC: run the legacy routine on chosen points and save the fixture.

Run from the repository root: ``python test_parity/probes/fwdxac.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['fwdxac', 'tlin3x', 'tlinex', 'tlin1x', 'glook', 'switch', 'quad']


def points():
    """(ATNSWP, TAPER, FACTOR, MACH) covering every branch and end mode."""
    rng = np.random.default_rng(4142)
    out = []
    for mach in (0.6, 1.6):
        for low, high in ((0.0, 1.0), (1.0, 6.0)):
            for _ in range(12):
                out.append((rng.uniform(1.0, 6.0), rng.uniform(0.0, 1.0),
                            rng.uniform(low, high), mach))
    out += [
        # Grid nodes, both sides of the abscissa and both regimes.
        (3.0, .25, .4, .6), (3.0, .25, 1 / .4, .6),
        (3.0, .25, .4, 1.6), (3.0, .25, 1 / .4, 1.6),
        # Quadratic extrapolation along the curve parameter, both ends.
        (0.4, .3, .5, .6), (8.5, .3, .5, .6), (9.0, .6, 3.0, 1.6),
        # Clamping in taper and at the abscissa ends.
        (2.5, 1.4, .7, .6), (2.5, -0.2, .7, 1.6), (2.5, .4, 0.0, .6),
        (2.5, .4, 1.0, .6), (2.5, .4, 1.0, 1.6),
        # Negative inputs, which the routine takes by magnitude.
        (-3.3, .45, -0.8, .6), (-3.3, .45, -2.2, 1.6),
        # Mach exactly one selects the supersonic tables.
        (4.0, .3, .5, 1.0), (4.0, .3, 2.0, 1.0),
        # The SUPT6 row the -.56 entry sits in.
        (5.0, 1.0, 1 / .6, 1.6),
    ]
    return [tuple(float(v) for v in p) for p in out]


def driver(all_points) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      REAL MACH', '      DIMENSION ROUTE(2)',
             "      DATA ROUTE /4HPROB,4HE   /",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, (atnswp, taper, factor, mach) in enumerate(all_points):
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(assign('ATNSWP', atnswp))
        lines.append(assign('TAPER', taper))
        lines.append(assign('FACTOR', factor))
        lines.append(assign('MACH', mach))
        lines.append('      CALL FWDXAC(ATNSWP,TAPER,FACTOR,MACH,ROUTE,XAC)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'XAC',XAC")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_points = points()
    records = parse_records(run('fwdxac', driver(all_points), ROUTINES))
    payload = [{'inputs': dict(zip(('atnswp', 'taper', 'factor', 'mach'), p)),
                'xac': r['XAC'][0]} for p, r in zip(all_points, records)]
    print(save('fwdxac', payload))


if __name__ == '__main__':
    main()
