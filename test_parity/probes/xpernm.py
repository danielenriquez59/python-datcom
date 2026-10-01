"""
Probe XPERNM, the experimental-namelist card counts, and save the
fixture.  Each case writes its deck to unit 8 and runs XPERNM on it.

Run from the repository root: ``python test_parity/probes/xpernm.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import parse_records, run, save  # noqa: E402

ROUTINES = ['xpernm', 'tbtrn']

COMMONS = """\
      COMMON /EXPER/ KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1               KBODY, KWING, KHT, KVT, KWB, KDWASH(3),
     2               ALPOW, ALPLW, ALPOH, ALPLH
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      LOGICAL KDWASH
"""

DECK = [
    ' $FLTCON NMACH=1.0,MACH(1)=0.6$',
    ' $EXPR01 CLAB(1)=0.1,',
    '   CLAB(2)=0.2$',
    ' $SYNTHS XCG=1.0$',
    ' $EXPR02 CMAB(1)=0.0$',
    ' $EXPRX1 CMAB(1)=0.0$',
    ' $EXPR CMAB(1)=0.0$',
    'CASEID ONE',
]


def cases():
    return [
        {'deck': DECK, 'nlist': [3, 5, 1500, 7] + [0] * 96},
        {'deck': DECK[:1] + [' $EXPR A=1$'] + DECK[1:4],
         'nlist': [0] * 100},
        {'deck': [' $EXPR A=1$'] * 103, 'nlist': [2] * 100},
        {'deck': [' $FLTCON A=1$', ' $SYNTHS B=2$'], 'nlist': [9] * 100},
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      CHARACTER*1 CD', "      CD='$'",
             "      READ(CD,'(A1)') KAND"]
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      OPEN(8,FILE='deck{n}.txt',STATUS='REPLACE')")
        for card in c['deck']:
            lines.append(f"      WRITE(8,'(A)') '{card}'")
        lines.append('      CLOSE(8)')
        lines.append(f"      OPEN(8,FILE='deck{n}.txt',STATUS='OLD')")
        for k, v in enumerate(c['nlist']):
            if v:
                lines.append(f'      NLIST({k + 1})={v}')
            else:
                lines.append(f'      NLIST({k + 1})=0')
        lines.append('      KLIST=1')
        lines.append('      CALL XPERNM')
        lines.append('      CLOSE(8)')
        lines.append("      WRITE(6,'(A,101I8)') 'R',KLIST,NLIST")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('xpernm', driver(all_cases), ROUTINES))
    print(save('xpernm', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
