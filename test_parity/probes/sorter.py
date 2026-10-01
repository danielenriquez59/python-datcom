"""
Probe SORTER, the INTEGER table sort, and save the fixture.

Run from the repository root: ``python test_parity/probes/sorter.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import parse_records, run, save  # noqa: E402

ROUTINES = ['sorter']
MARK = -7


def cases():
    t = [[5, 3, 9, 1], [2, 8, 2, 7], [5, 1, 4, 6]]
    big = [[2000000000, 3], [1500000000, 9], [-4, 1]]
    return [
        {'a': t, 'irow': 0, 'icol': 1},
        {'a': t, 'irow': 2, 'icol': 0},
        {'a': t, 'irow': 3, 'icol': 3},
        {'a': t, 'irow': 0, 'icol': 0},
        {'a': [[1, 2], [3, 4]], 'irow': 1, 'icol': 1},     # already sorted
        {'a': [[4, -2, 4], [4, 0, -9], [-1, 5, 4]], 'irow': 1, 'icol': 1},
        {'a': big, 'irow': 0, 'icol': 1},                  # sentinel wraps
        {'a': [[7]], 'irow': 1, 'icol': 1},
        {'a': [[3, 1, 2]], 'irow': 1, 'icol': 1},
        {'a': [[3], [1], [2], [1]], 'irow': 1, 'icol': 1},
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', '      INTEGER A(40),B(40),IR(10),IC(10)']
    for n, c in enumerate(all_cases):
        rows, cols = len(c['a']), len(c['a'][0])
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {100 + n} K=1,40')
        lines.append(f'      B(K)={MARK}')
        lines.append(f'{100 + n:5d} CONTINUE')
        for j in range(cols):
            for i in range(rows):
                lines.append(f'      A({i + 1 + rows * j})={c["a"][i][j]}')
        lines.append(f'      IFLAG={MARK}')
        lines.append(f"      CALL SORTER(A,{rows},{cols},{c['irow']},"
                     f"{c['icol']},B,IR,IC,IFLAG)")
        size = rows * cols
        lines.append(f"      WRITE(6,'(A,{size}I12)') 'A',(A(K),K=1,{size})")
        lines.append(f"      WRITE(6,'(A,{size}I12)') 'B',(B(K),K=1,{size})")
        lines.append(f"      WRITE(6,'(A,{rows}I4)') 'IR',(IR(K),K=1,{rows})")
        lines.append(f"      WRITE(6,'(A,{cols}I4)') 'IC',(IC(K),K=1,{cols})")
        lines.append("      WRITE(6,'(A,I4)') 'FLAG',IFLAG")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('sorter', driver(all_cases), ROUTINES))
    print(save('sorter', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
