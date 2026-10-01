"""
Probe the namelist reader internals (FINDCH, SKIPBL, EXTRST, TODEC, FINDVN,
FORINT, FORLOG, FORREA, READCD) and save their raw output.

Run from the repository root: ``python test_parity/probes/namelist_io.py``.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import run  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'namelist_io.json')

NUMBERS = ['1.5', '-2.25E3', '+.125', '  7  ', '3.E-2', '1.0E+10',
           '12345678901', '0.123456789012', '4.5E', '1-2', '2.5Q',
           '.5E-3X', '1E22', '6.02E23', '1.5E-300', '--1', '1.2.3',
           '9.87654E 2', '']
FINDS = [('ALPHA=1.0', '=', 1), ('   B', 'B', 2), ('NONE', '$', 1),
         ('   X', ' ', 4)]
BLANKS = [('    X', 1), ('X', 1), ('', 3), ('  $', 2)]
NAMES = ['MACH', 'ALT', 'NALPHA', 'ALSCHD']
LOOKUPS = ['NALPHA', 'ALT', 'ALSCHD ', 'ALTX', 'MAC', 'ALSCHD(2)', '']
INTS = [[1, -22, 333, 4444, 55555, 666666, 7777777, 88888888, 9],
        [0], [-2147483647, 5]]
LOGS = [[True, False, True], [False] * 10]
REALS = [[1.0, -2.5, 0.0, 1.0e-30, 123456.7, 1.0e7, 9.9999999e6, 0.1,
          0.09999999, 3.14159265], [-7.0e-8]]
CARDS = ['  $FLTCON NMACH=2.0, MACH(1)=0.6,0.8$', '', 'X' * 90, 'LAST']


def _card(text):
    return text.ljust(80)[:80]


def _split(stmt):
    out, first = [], True
    while stmt:
        out.append(('      ' if first else '     1') + stmt[:66])
        stmt, first = stmt[66:], False
    return out


def _load(lines, name, text, label):
    lines += _split(f"CARD='{text}'")
    lines.append(f'      DO {label} J=1,80')
    lines.append(f"{label:5d} {name}(J)=TRANSFER(CARD(J:J)//'   ',{name}(J))")


def driver() -> str:
    lines = ['      PROGRAM PROBE',
             '      INTEGER KOL(80),STR(80),VN(40),LENVN(4),CH',
             '      LOGICAL FOUND,IEOF,LG(10)',
             '      DIMENSION IV(10),RV(10),NM(8)',
             '      CHARACTER*80 CARD',
             '      CHARACTER*8 NAMSTR']
    lab = [1000]

    def nxt():
        lab[0] += 1
        return lab[0]

    def case():
        lines.append(f"      WRITE(6,'(A)') '@@CASE'")

    for text, ch, k in FINDS:
        case()
        _load(lines, 'KOL', _card(text), nxt())
        lines.append(f"      CH=TRANSFER('{ch}   ',CH)")
        lines.append(f'      K={k}')
        lines.append('      CALL FINDCH(KOL,CH,K)')
        lines.append("      WRITE(6,'(I4)') K")
    for text, k in BLANKS:
        case()
        _load(lines, 'KOL', _card(text), nxt())
        lines.append(f'      K={k}')
        lines.append('      CALL SKIPBL(KOL,K)')
        lines.append("      WRITE(6,'(I4)') K")
    case()
    _load(lines, 'KOL', _card('ABCDEFGHIJ'), nxt())
    lines.append('      CALL EXTRST(KOL,3,7,STR)')
    lines.append("      WRITE(6,'(5A1)') (STR(J),J=1,5)")
    for text in NUMBERS:
        case()
        _load(lines, 'KOL', _card(text), nxt())
        lines.append('      CALL TODEC(KOL,ANS,IERR)')
        lines.append("      WRITE(6,'(ES25.16,I3)') ANS,IERR")
    packed = ''.join(NAMES).ljust(40)
    lines += _split(f"CARD='{packed}'")
    lines.append(f'      DO {nxt()} J=1,40')
    lines.append(f"{lab[0]:5d} VN(J)=TRANSFER(CARD(J:J)//'   ',VN(J))")
    for k, n in enumerate(NAMES):
        lines.append(f'      LENVN({k + 1})={len(n)}')
    for text in LOOKUPS:
        case()
        _load(lines, 'KOL', _card(text), nxt())
        lines.append('      CALL FINDVN(4,LENVN,KOL,VN,40,NVN,FOUND)')
        lines.append("      WRITE(6,'(I4,L2)') NVN,FOUND")
    lines.append("      NAMSTR='VARNAME '")
    lines.append(f'      DO {nxt()} J=1,8')
    lines.append(f"{lab[0]:5d} NM(J)=TRANSFER(NAMSTR(J:J)//'   ',NM(J))")
    for vals in INTS:
        case()
        for k, v in enumerate(vals):
            lines.append(f'      IV({k + 1})={v}')
        lines.append(f'      CALL FORINT(6,NM,IV,{len(vals)})')
    for vals in LOGS:
        case()
        for k, v in enumerate(vals):
            lines.append(f"      LG({k + 1})=.{'TRUE' if v else 'FALSE'}.")
        lines.append(f'      CALL FORLOG(6,NM,LG,{len(vals)})')
    for vals in REALS:
        case()
        for k, v in enumerate(vals):
            s = repr(v)
            s = s.replace('e', 'D') if 'e' in s else s + 'D0'
            lines.append(f'      RV({k + 1})={s}')
        lines.append(f'      CALL FORREA(6,NM,RV,{len(vals)})')
    lines.append("      OPEN(20,FILE='cards.txt',STATUS='UNKNOWN')")
    for c in CARDS:
        lines += _split(f"WRITE(20,'(A)') '{c}'")
    lines.append('      REWIND 20')
    _load(lines, 'KOL', _card('PREVIOUS'), nxt())
    # A read at end of file is left out: gfortran stops READCD's retry
    # with "Read past ENDFILE record".
    for _ in range(len(CARDS)):
        case()
        lines.append('      CALL READCD(20,KOL,IEOF)')
        lines.append("      WRITE(6,'(80A1,L2)') KOL,IEOF")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = run('namelist_io', driver(), [
        'findch', 'skipbl', 'extrst', 'todec', 'findvn', 'forint',
        'forlog', 'forrea', 'readcd'])
    records, current = [], None
    for line in out.split('\n'):
        if line.startswith('@@CASE'):
            current = []
            records.append(current)
        elif current is not None:
            current.append(line)
    if records and records[-1] and records[-1][-1] == '':
        records[-1].pop()
    FIXTURE.write_text(json.dumps({
        'finds': FINDS, 'blanks': BLANKS, 'numbers': NUMBERS,
        'names': NAMES, 'lookups': LOOKUPS, 'ints': INTS, 'logs': LOGS,
        'reals': REALS, 'cards': CARDS, 'records': records}, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
