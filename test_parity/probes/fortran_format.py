"""
Check pydatcom.io.fortran_format against gfortran: write a battery of
formats and items with the compiled compiler and save the raw records.

Run from the repository root: ``python test_parity/probes/fortran_format.py``.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import run  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'fortran_format.json')

REALS = [0.0, -0.0, 1.0, -1.0, 123.456, -0.000123456, 0.5, 0.05, 9.99995,
         99999.5, 1.0e-30, -1.0e-30, 2.5e120, 6.02e23, 0.125, 1234567.0,
         -98.7654321, 0.0999999, 1.0e10]
CASES = [
    ('(1P,5E12.5)', REALS),
    ('(5E12.5)', REALS),
    ('(5E14.6)', REALS),
    ('(4E10.3)', REALS),
    ('(2P,4E13.4)', REALS),
    ('(-1P,4E13.4)', REALS),
    ('(5F10.3)', REALS),
    ('(5F8.0)', REALS),
    ('(6F6.2)', REALS),
    ('(1P,5F10.3)', REALS),
    ('(4G12.4)', REALS),
    ('(1X,I5,I3,I1,I8)', [12, -34, 7, 123456789, 5, 1000, 0, -1]),
    ('(3I4.3)', [5, -5, 0]),
    ('(A,A6,A2,1X,A4)', ['HELLO', 'AB', 'WXYZ', 'CD  ']),
    ('(1H ,5X,A1,1H(,I3,2H)=1P  ,E12.5)', ['X', 3, 1.5]),
    ('(5(5X,A1,1H(,I3,2H)=,1P,E12.5))',
     [v for i in range(7) for v in ('Q', i + 1, 1.5 * i - 2.0)]),
    ('(1X,2HAB,2(I3,F6.2),/,1X,3HEND)', [1, 2.5, 3, -4.25]),
    ('(1X,I3,3X)', [7]),
    ('(T10,A3,TL6,A2,TR2,I2)', ['ABC', 'XY', 42]),
    ('(1X,2L3,L1)', [True, False, True]),
    ("(1X,'IT''S',I3,' DONE')", [5]),
    ('(1X,I3,(2F6.1))', [9, 1.0, 2.0, 3.0, 4.0, 5.0]),
    ('(2X,3(I2,1X),F5.1)', [1, 2, 3, 4.5, 6, 7]),
    ('(1X,I3,5X,4HTEXT)', [8]),
    ('(8X,E11.4,F9.4,E11.4)', [3.14159, -0.001, 1.0e-5]),
]


def _lit(v):
    if isinstance(v, bool):
        return '.TRUE.' if v else '.FALSE.'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        s = repr(v)
        return s.replace('e', 'D') if 'e' in s else s + 'D0'
    return "'" + v.replace("'", "''") + "'"


def _split(line):
    """Break a long statement into fixed-form continuation lines."""
    out, first = [], True
    while line:
        room = 66
        out.append(('      ' if first else '     1') + line[:room])
        line, first = line[room:], False
    return out


def driver() -> str:
    lines = ['      PROGRAM PROBE']
    for n, (fmt, items) in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') '@@CASE',{n}")
        fq = fmt.replace("'", "''")
        stmt = f"WRITE(6,'{fq}') " + ','.join(_lit(v) for v in items)
        lines += _split(stmt)
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = run('fortran_format', driver(), [])
    records, current = [], None
    for line in out.split('\n'):
        if line.startswith('@@CASE'):
            current = []
            records.append(current)
        elif current is not None:
            current.append(line)
    if records and records[-1] and records[-1][-1] == '':
        records[-1].pop()
    payload = [{'format': f, 'items': i, 'records': r}
               for (f, i), r in zip(CASES, records)]
    FIXTURE.write_text(json.dumps(payload, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
