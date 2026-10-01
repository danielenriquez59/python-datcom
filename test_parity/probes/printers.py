"""
Probe DMPARY, PRCSID and SWRITE, and save their raw printed records.

Run from the repository root: ``python test_parity/probes/printers.py``.
SWRITE is driven with AUXOUT's own formats ``FOR`` and ``FOR1``.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import run  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'printers.json')
U = 1.0e-30

FOR = ['(   ', '1H0 ', ',0P ', 'F6.3', ',1X ', ',1X ', ',1X,', 'F8.2',
       ',1X ', ',1X ', ',1X,', 'F8.2', ',1X ', ',2X ', ',1P ', 'E10.',
       '4   ', ',3X ', ',0P ', 'F9.3', ',1X ', ',4X ', ',1P ', 'E10.',
       '4   ', ',9X ', ',0P ', 'F9.3', ',1X ', ',1X ', ',1X,', 'F8.3',
       ',1X ', '    ', ',1X,', 'F8.3', ',1X ', '    ', ',1X,', 'F8.3',
       ',1X ', '    ', ',1X,', 'F8.3', ',1X ', ')   ']
ICONT = [7] + [10] * 10
FOR1 = ['(1H ', ',10X', ',1X,', 'E11.', '4,1X', ',1X ', ',1X,', 'F6. ',
        '2,1X', ',1X ', ',1X,', 'F5. ', '2,1X', ',2X ', ',1X,', 'F9. ',
        '4,1X', ',5X ', ',1X,', 'E11.', '4,1X', ',1X ', ',1X,', 'E11.',
        '4,1X', ',1X ', ',1X,', 'E11.', '4,1X', ',3X ', ',1X,', 'E11.',
        '4   ', ')   ']
ICON = [13, 7, 7, 11, 13, 13, 13, 13]

DMP = [
    {'array': [1.5, -2.25e-7, 3.0e12, 0.0, 7.0, 8.5, -9.0], 'name': 'ALPH',
     'nlet': 1},
    {'array': [0.001, 12345.678], 'name': 'CL  ', 'nlet': 2},
    {'array': [1.0, 2.0, 3.0, 4.0, 5.0], 'name': 'ABCD', 'nlet': 4},
    {'array': [6.5], 'name': 'XYZ ', 'nlet': 3},
]
TITLES = ['  SAMPLE WING-BODY CASE', '', 'X', 'A' * 74,
          '   PADDED   TITLE    ']
SWR = [
    {'form': FOR, 'icont': ICONT, 'last': 11, 'columns':
     [[0.6, 0.9, 1.2], [U, 5000.0, -U], [800.0, 2 * U, 900.0],
      [2116.2, U, 1500.0], [518.7, 400.0, -U], [2.5e6, U, 3.1e7],
      [300.0, 300.0, 2 * U], [8.0, -U, 8.0], [30.0, 30.0, U],
      [12.5, 12.5, 12.5], [U, 0.0, -1.25]]},
    {'form': FOR1, 'icont': ICON, 'last': 8, 'columns':
     [[1.25e-3, U], [2 * U, 4.5], [-U, 0.25], [0.035, U], [U, 2 * U],
      [1.1e-4, -U], [2.2e-5, 3.3e-6], [U, 7.7]]},
]


def _real(v):
    s = repr(float(v))
    return s.replace('e', 'D') if 'e' in s else s + 'D0'


def _split(stmt):
    out, first = [], True
    while stmt:
        out.append(('      ' if first else '     1') + stmt[:66])
        stmt, first = stmt[66:], False
    return out


def driver() -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /FLOLOG/ LDM(10),HEAD,LDM2(14)',
             '      COMMON /CASEID/ IDCSE(74)',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND',
             '      LOGICAL LDM,HEAD,LDM2,NDMF,NAF',
             '      CHARACTER*74 TITLE',
             '      CHARACTER*4 FOR(46),FOR1(34)',
             '      DIMENSION ARR(10),ICONT(11),ICON(8),C(3,14)',
             '      UNUSED=1.D-30']
    for k, w in enumerate(FOR):
        lines.append(f"      FOR({k + 1})='{w}'")
    for k, w in enumerate(FOR1):
        lines.append(f"      FOR1({k + 1})='{w}'")
    for k, v in enumerate(ICONT):
        lines.append(f'      ICONT({k + 1})={v}')
    for k, v in enumerate(ICON):
        lines.append(f'      ICON({k + 1})={v}')
    n = 0
    for c in DMP:
        lines.append(f"      WRITE(6,'(A,I4)') '@@CASE',{n}")
        for k, v in enumerate(c['array']):
            lines.append(f'      ARR({k + 1})={_real(v)}')
        lines.append(f"      CALL DMPARY(ARR,{len(c['array'])},"
                     f"4H{c['name']},{c['nlet']})")
        n += 1
    for head in (True, False):
        for t in TITLES:
            lines.append(f"      WRITE(6,'(A,I4)') '@@CASE',{n}")
            lines.append(f"      HEAD=.{'TRUE' if head else 'FALSE'}.")
            lines += _split(f"TITLE='{t}'")
            lines.append(f'      DO {600 + n} J=1,74')
            lines.append(f"{600 + n:5d} IDCSE(J)=TRANSFER(TITLE(J:J)//'   ',"
                         'IDCSE(J))')
            lines.append('      CALL PRCSID')
            n += 1
            if not head:
                break
    for c in SWR:
        lines.append(f"      WRITE(6,'(A,I4)') '@@CASE',{n}")
        for i, col in enumerate(c['columns']):
            for r, v in enumerate(col):
                lines.append(f'      C({r + 1},{i + 1})=' + (
                    'UNUSED' if v == U else '-UNUSED' if v == -U else
                    '2*UNUSED' if v == 2 * U else _real(v)))
        rows = len(c['columns'][0])
        name, icn = ('FOR', 'ICONT') if c['form'] is FOR else ('FOR1', 'ICON')
        args = ','.join(f'C(1,{i + 1})' for i in range(14))
        lines.append('      NDMF=.FALSE.')
        lines.append('      NAF=.FALSE.')
        lines += _split(f"CALL SWRITE({c['last']},{name},{len(c['form'])},"
                        f"{icn},{rows},{args},NDMF,NAF)")
        lines.append("      WRITE(6,'(A,2L2)') '@@FLAGS',NDMF,NAF")
        n += 1
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = run('printers', driver(), ['dmpary', 'prcsid', 'swrite'])
    records, current = [], None
    for line in out.split('\n'):
        if line.startswith('@@CASE'):
            current = []
            records.append(current)
        elif current is not None:
            current.append(line)
    if records and records[-1] and records[-1][-1] == '':
        records[-1].pop()
    FIXTURE.write_text(json.dumps({'dmpary': DMP, 'titles': TITLES,
                                   'swrite': SWR, 'records': records},
                                  indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
