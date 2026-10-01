"""
Probe VNAME, LVALUE and RVALUE, the namelist card syntax checks, and save
the fixture.

Run from the repository root: ``python test_parity/probes/namelist_check.py``.
Each card is 80 columns padded to 84 with ``#``, the probe's stand-in for
whatever follows the card in memory.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import parse_records, run, save  # noqa: E402

ROUTINES = ['vname', 'lvalue', 'rvalue']

# NMTEST as in nmlist.f, whose other routines pull in the whole reader.
NMTEST = """\n      LOGICAL FUNCTION NMTEST(KOL,KEY,NCHAR)
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

HEAD = """\
      PROGRAM PROBE
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      INTEGER KOL(84), NUMBER(15), NAMES(40), LEN(8)
      INTEGER BLANK, EQUAL, COMMA, KD
      LOGICAL FOUND, ARRAY
      CHARACTER*84 CARD
      CHARACTER*40 NAMSTR
      DATA NUMBER / 4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,
     1       4H5   ,4H6   ,4H7   ,4H8   ,4H9   ,4H+   ,4H-   ,
     2       4H.   ,4H*   ,4HE   /
      DATA BLANK / 4H     /, EQUAL / 4H=    /, COMMA / 4H,    /
      DATA KD / 4H$    /
      KAND=KD
"""

NAMES = ['ALSCHD', 'NALPHA', 'MACH', 'ALT', 'WGPL', 'CD']


def card(text):
    return text.ljust(80)[:80] + '####'


def cases():
    out = []
    for text, l in [('  NALPHA=9.0,', 3), ('  ALSCHD(3)=1.,', 3),
                    ('  ALSCHD(0)=1.,', 3), ('  ALSCHD(A)=1.,', 3),
                    ('  MACHX=2.', 3), ('  CD =2.', 3), ('  ALT(12)=1.', 3),
                    (' ' * 74 + 'MACH(7', 75), ('  WGPL$', 3)]:
        out.append({'routine': 'VNAME', 'card': card(text), 'l': l,
                    'names': NAMES, 'nf': 1})
    for text, l in [('.TRUE.,.FALSE.$', 1), ('3*.TRUE., .FALSE.,', 1),
                    ('*.TRUE.,2.FALSE.$', 1), ('.TRUE. .FALSE.', 1),
                    ('2**.TRUE.,', 1), ('', 1), ('XYZ', 1),
                    (' ' * 78 + '0*', 79), (' ' * 74 + '.TRUE.', 75)]:
        out.append({'routine': 'LVALUE', 'card': card(text), 'l': l,
                    'ndml': 0, 'nf': 0})
    for text, l in [('1.0,2.5,-3.E2,4.5E-1$', 1), ('3*1.5,2*-.25$', 1),
                    ('1,2.,3.E,+-4.$', 1), ('1. 2.,', 1), ('..5,1.E2E3,', 1),
                    ('*2.,3**4.,', 1), ('$', 1), ('E5,', 1), ('', 1),
                    (' ' * 76 + '1.5,', 77), ('-1.2E+03,2.E-4$', 1),
                    (' ' * 77 + '7.1', 78)]:
        out.append({'routine': 'RVALUE', 'card': card(text), 'l': l,
                    'ndml': 0, 'nf': 0})
    return out


def driver(all_cases) -> str:
    lens = ','.join(str(len(n)) for n in NAMES)
    head = HEAD.replace('      KAND=KD',
                        f'      DATA LEN / {lens},0,0 /\n      KAND=KD')
    lines = [head.rstrip('\n'), f"      NAMSTR='{''.join(NAMES)}'",
             "      DO 4 J=1,40",
             "    4 NAMES(J)=TRANSFER(NAMSTR(J:J)//'   ',NAMES(J))"]
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        text = c['card']
        lines.append(f"      CARD='{text[:40]}'")
        lines.append(f"      CARD(41:84)='{text[40:]}'")
        # A formatted A1 read would take a comma as a field separator.
        lines.append(f"      DO {5000 + n} J=1,84")
        lines.append(f"{5000 + n:5d} KOL(J)=TRANSFER(CARD(J:J)//'   ',"
                     "KOL(J))")
        lines.append(f"      L={c['l']}")
        lines.append(f"      NF={c['nf']}")
        if c['routine'] == 'VNAME':
            lines.append('      ARRAY=.FALSE.')
            lines.append('      NDMS=0')
            lines.append(f'      CALL VNAME(KOL,L,LEN,{len(NAMES)},NAMES,I,'
                         'FOUND,ARRAY,')
            lines.append('     1  NDMS,NF,NUMBER,BLANK,EQUAL)')
            lines.append("      IFD=0")
            lines.append("      IA=0")
            lines.append("      IF(FOUND) IFD=1")
            lines.append("      IF(ARRAY) IA=1")
            lines.append("      WRITE(6,'(A,6I6)') 'OUT',L,I,IFD,IA,NDMS,NF")
        else:
            lines.append('      NDML=0')
            lines.append(f"      CALL {c['routine']}(KOL,L,NDML,NF,BLANK,"
                         "COMMA,NUMBER)")
            lines.append("      WRITE(6,'(A,3I6)') 'OUT',L,NDML,NF")
    lines += ['      STOP', '      END', NMTEST]
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('namelist_check', driver(all_cases),
                                ROUTINES))
    print(save('namelist_check', [{'inputs': c, 'outputs': r}
                                  for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
