"""
Probe NAMEW (namelist printing) and TOINT (integer decoding), and save the
raw records.

Run from the repository root: ``python test_parity/probes/namelist_write.py``.

NAMEW hands a REAL block to FORLOG as LOGICAL, which lines up only when
REAL and LOGICAL words are the same size, as in the source's own build, so
NAMEW is compiled without the REAL*8 promotion, on exactly representable
values.  TOINT runs in the usual REAL*8 build.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'namelist_write.json')

NAMES = ['MACH', 'ALT', 'WGPL', 'TRIM', 'SREF', 'DUMMY']
LENVN = [4, 3, 4, 4, 4, 5]
# Length (negative LOGICAL, zero reuses the last type) and position.
LISTS = [
    {'vdime': [3, 2, -2, 0, 1, 4], 'loc': [1, 4, 6, 8, 9, 0]},
    {'vdime': [9, 0, -1, 1, 0, 0], 'loc': [1, 0, 10, 0, 11, 12]},
]
REALS = [0.5, 1.25, -3.0, 1500.0, 2.0e4, None, None, 1.0, 250.0, None,
         0.125, -1024.0, 7.75, 64.0, 0.0625, 3.5]
LOGICALS = {6: True, 7: False, 10: True}
INTS = ['42', '-7', '3.6', '-2.5', '1E3', '12X', '  +15  ', '', '2147483',
        '0.49']


def _card(text):
    return text.ljust(80)[:80]


def namew_driver() -> str:
    lines = ['      PROGRAM PROBE',
             '      INTEGER VN(24),LENVN(6),VDIME(6),LOC(6),NL(6),KD',
             '      DIMENSION CB(16)',
             '      LOGICAL LB(16)',
             '      EQUIVALENCE (CB,LB)',
             '      CHARACTER*24 S',
             "      S='" + ''.join(NAMES).ljust(24) + "'",
             '      DO 1 J=1,24',
             "    1 VN(J)=TRANSFER(S(J:J)//'   ',VN(J))",
             "      S='FLTCON'",
             '      DO 2 J=1,6',
             "    2 NL(J)=TRANSFER(S(J:J)//'   ',NL(J))",
             "      KD=TRANSFER('$   ',KD)"]
    for k, n in enumerate(LENVN):
        lines.append(f'      LENVN({k + 1})={n}')
    for k, v in enumerate(REALS):
        if v is not None:
            lines.append(f'      CB({k + 1})={v!r}')
    for k, v in LOGICALS.items():
        lines.append(f"      LB({k})=.{'TRUE' if v else 'FALSE'}.")
    for c in LISTS:
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        for k in range(6):
            lines.append(f"      VDIME({k + 1})={c['vdime'][k]}")
            lines.append(f"      LOC({k + 1})={c['loc'][k]}")
        lines.append('      CALL NAMEW(KD,6,NL,6,VN,24,LENVN,6,VDIME,CB,16,'
                     'LOC)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def toint_driver() -> str:
    lines = ['      PROGRAM PROBE', '      INTEGER KOL(80)',
             '      CHARACTER*80 CARD']
    for n, text in enumerate(INTS):
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        lines.append(f"      CARD='{_card(text)[:40]}'")
        lines.append(f"      CARD(41:80)='{_card(text)[40:]}'")
        lines.append(f'      DO {100 + n} J=1,80')
        lines.append(f"{100 + n:5d} KOL(J)=TRANSFER(CARD(J:J)//'   ',"
                     'KOL(J))')
        lines.append('      CALL TOINT(KOL,IANS,IERR)')
        lines.append("      WRITE(6,'(2I12)') IANS,IERR")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def _records(out):
    records, current = [], None
    for line in out.split('\n'):
        if line.startswith('@@CASE'):
            current = []
            records.append(current)
        elif current is not None:
            current.append(line)
    if records and records[-1] and records[-1][-1] == '':
        records[-1].pop()
    return records


def main():
    saved = probe.FLAGS
    probe.FLAGS = saved.replace(' -fdefault-real-8', '').replace(
        ' -fdefault-double-8', '')
    try:
        namew = _records(probe.run('namelist_write', namew_driver(),
                                   ['namew', 'forlog', 'forint', 'forrea']))
    finally:
        probe.FLAGS = saved
    toint = _records(probe.run('toint', toint_driver(),
                               ['toint', 'findch', 'todec']))
    FIXTURE.write_text(json.dumps({
        'names': NAMES, 'lenvn': LENVN, 'lists': LISTS, 'reals': REALS,
        'logicals': LOGICALS, 'ints': INTS, 'namew': namew,
        'toint': toint}, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
