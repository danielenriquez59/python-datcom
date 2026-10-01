"""
Probe TESTOR, the namelist card checker, over a sequence of cards in one
program (its END, NDML, NDM and IS are saved between calls), and save the
unit 6 and unit 11 records.

Run from the repository root: ``python test_parity/probes/testor.py``.
Cards are 80 columns padded to 84 with ``#``, the probe's stand-in for
what follows the card in memory.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'testor.json')
NAMES = [('NMACH', 1), ('MACH', 20), ('ALSCHD', 20), ('STMACH', 1),
         ('TRIM', -1), ('FLAGS', -3), ('RNNUB', 20)]
# (card, start column, NAM, K, IER)
CARDS = [
    (' $FLTCON NMACH=2.0, MACH(1)=0.6,0.8,', 9, 2, 1, 0),
    ('  ALSCHD=-2.,0.,2.,4.,', 1, 2, 2, 0),
    ('  6.,8.,10.$', 1, 2, 2, 0),
    (' $FLTCON TRIM=.TRUE., FLAGS=2*.FALSE.,.TRUE.$', 9, 2, 1, 0),
    (' $FLTCON XYZ=1.0, NMACH 3.0, STMACH=1.,2.$', 9, 2, 1, 0),
    (' $FLTCON NMACH(1)=2, MACH(3)=0.5,0.6,0.7$', 9, 2, 1, 0),
    (' $FLTCON MACH=0.1,0.2', 9, 2, 1, 0),
    (' $FLTCON NMACH=1.0$', 9, 2, 1, 0),
    ('  MACH=0.3$', 1, 2, 2, 0),
    (' $BADNAM X=1.$', 9, 2, 1, 1),
    (' SOME TEXT', 1, 1, 1, 0),
    (' $FLTCON MACH(0)=1., ALSCHD(X)=2.$', 9, 2, 1, 0),
    (' $FLTCON RNNUB=' + ','.join(['1.E6'] * 13) + '$', 9, 2, 1, 0),
    (' $FLTCON TRIM=.T., FLAGS=.F.,.TRUE,.FALSE.$', 9, 2, 1, 0),
]


def _card(text):
    return text.ljust(80)[:80] + '####'


def driver():
    packed = ''.join(n for n, _ in NAMES)
    lines = ['      PROGRAM PROBE',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND',
             f'      INTEGER KOL(84),NAMES({len(packed)}),LEN({len(NAMES)}),'
             f'LDM({len(NAMES)}),KD',
             '      CHARACTER*84 CARD',
             f'      CHARACTER*{len(packed)} NAMSTR',
             "      KD=TRANSFER('$   ',KD)", '      KAND=KD',
             "      OPEN(11,FILE='unit11.txt',STATUS='UNKNOWN')",
             f"      NAMSTR='{packed}'",
             f'      DO 10 J=1,{len(packed)}',
             "   10 NAMES(J)=TRANSFER(NAMSTR(J:J)//'   ',NAMES(J))"]
    for k, (n, d) in enumerate(NAMES):
        lines.append(f'      LEN({k + 1})={len(n)}')
        lines.append(f'      LDM({k + 1})={d}')
    for n, (text, col, nam, k, ier) in enumerate(CARDS):
        card = _card(text)
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        lines.append("      WRITE(11,'(A)') '@@CASE'")
        lines.append(f"      CARD='{card[:42]}'")
        lines.append(f"      CARD(43:84)='{card[42:]}'")
        lines.append(f'      DO {100 + n} J=1,84')
        lines.append(f"{100 + n:5d} KOL(J)=TRANSFER(CARD(J:J)//'   ',"
                     'KOL(J))')
        lines.append(f'      L={col}')
        lines.append(f'      CALL TESTOR(KOL,L,{nam},{k},{ier},LEN,LDM,'
                     f'{len(NAMES)},NAMES)')
        lines.append("      WRITE(6,'(A,I4)') '@@L',L")
    lines += ['      REWIND 11', "      WRITE(6,'(A)') '@@UNIT11'",
              '   20 READ(11,\'(A)\',END=30) CARD',
              "      WRITE(6,'(A)') CARD", '      GO TO 20',
              '   30 STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = probe.run('testor', driver(),
                    ['testor', 'vname', 'lvalue', 'rvalue'],
                    stubs=probe.STUBS + NMTEST)
    main_part, unit11 = out.split('@@UNIT11\n')

    def split(text):
        cases, current = [], None
        for line in text.split('\n'):
            if line.startswith('@@CASE'):
                current = []
                cases.append(current)
            elif current is not None:
                current.append(line)
        return cases

    six = split(main_part)
    eleven = [[ln.rstrip() for ln in c if ln.strip()]
              for c in split(unit11)]
    records = []
    for c6, c11 in zip(six, eleven):
        l_line = [ln for ln in c6 if ln.startswith('@@L')][0]
        records.append({'unit6': [ln for ln in c6 if ln and
                                  not ln.startswith('@@')],
                        'unit11': c11, 'l': int(l_line.split()[1])})
    FIXTURE.write_text(json.dumps({'names': NAMES, 'cards': CARDS,
                                   'records': records}, indent=1))
    print(FIXTURE)


NMTEST = """\
      LOGICAL FUNCTION NMTEST(KOL,KEY,NCHAR)
      INTEGER KOL(1) , KEY(1)
      DO 1000 I=1,NCHAR
         IF(KOL(I).NE.KEY(I)) GO TO 1010
 1000 CONTINUE
      NMTEST=.TRUE.
      RETURN
 1010 NMTEST=.FALSE.
      RETURN
      END
"""

if __name__ == '__main__':
    main()
