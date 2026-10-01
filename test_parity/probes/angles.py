"""
Probe ANGLES and save the fixture.

Each case loads a six-word record, which serves as the previous state,
calls ANGLES with one entry, and prints the record it leaves.

Run from the repository root: ``python test_parity/probes/angles.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

EPS = 1.52587890625e-5


def record(degrees):
    """A consistent record for an angle, as ANGLES would leave it."""
    r = math.radians(degrees)
    return [degrees, r, math.sin(r), math.cos(r),
            math.tan(r) if abs(math.cos(r)) > EPS else 1.048576e6, r]


def cases():
    rng = np.random.default_rng(66)
    fresh = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    out = []
    # Entry 1 and 2 over the full circle and beyond, from a fresh record.
    for deg in [-400., -270., -180., -135., -90., -45., -1., 0., 1e-4, 1.,
                30., 45., 89.9999, 90., 135., 180., 181., 270., 725.]:
        out.append((1, [deg] + fresh[1:]))
        out.append((2, [0.0, math.radians(deg)] + fresh[2:]))
    # Entries 3 to 6 from a fresh record, including the near-zero and
    # near-unity branches.
    for value in [-1.0, -0.99999, -0.7, -1e-6, 0.0, 1e-6, 0.3, 0.99999, 1.0]:
        out.append((3, [0., 0., value, 1., 0., 0.]))
        out.append((4, [0., 0., 0., value, 0., 0.]))
    for value in [-3e6, -5.0, -1.0, -1e-6, 0.0, 0.2, 1.0, 7.0, 3e6]:
        out.append((5, [0., 0., 0., 1., value, 0.]))
    for s, c in [(0., 0.), (3., 4.), (-3., 4.), (3., -4.), (-3., -4.),
                 (1., 0.), (0., -2.), (1e-7, 1.)]:
        out.append((6, [0., 0., s, c, 0., 0.]))
    # The stateful early return: a new angle within EPS of the last one
    # leaves the old trigonometric values in place.
    for deg in (10.0, -35.0, 60.0):
        old = record(deg)
        stale = old[:2] + [9., 9., 9.] + old[5:]
        out.append((1, [deg + 0.5 * math.degrees(EPS)] + stale[1:]))
        out.append((1, [deg + 2.0 * math.degrees(EPS)] + stale[1:]))
        out.append((3, [0., 0., math.sin(math.radians(deg)), 9., 9., 0.]))
        out.append((4, [0., 0., 9., math.cos(math.radians(deg)), 9.,
                        math.radians(deg)]))
    # Random records of all six entries, with random previous state.
    for _ in range(30):
        entry = int(rng.integers(1, 7))
        arg = list(record(rng.uniform(-179, 179)))
        arg[entry - 1] = (rng.uniform(-200, 200) if entry == 1 else
                          rng.uniform(-3, 3) if entry == 2 else
                          rng.uniform(-1, 1) if entry in (3, 4, 6) else
                          rng.uniform(-20, 20))
        out.append((entry, arg))
    return [(e, [float(v) for v in a]) for e, a in out]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', '      DIMENSION ARG(6)']
    for n, (entry, arg) in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        for i, v in enumerate(arg):
            lines.append(assign(f'ARG({i + 1})', v))
        lines.append(f'      CALL ANGLES({entry},ARG)')
        lines.append("      WRITE(6,'(A,6ES25.16)') 'ARG',ARG")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('angles', driver(all_cases), ['angles']))
    payload = [{'entry': e, 'arg': a, 'result': r['ARG']}
               for (e, a), r in zip(all_cases, records)]
    print(save('angles', payload))


if __name__ == '__main__':
    main()
