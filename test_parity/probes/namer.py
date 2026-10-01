"""
Probe NAMER, DATCOM's namelist reader, on the FLTCON namelist of XNAM1,
and save the fixture.

Run from the repository root: ``python test_parity/probes/namer.py``.

The shipped datcom_2000/namer.f does not compile: in its ``DATA CARET``
the caret has become a carriage return (``4H<CR>   /``).  The build copy restores ``4H^`` from
the earlier datcom.f.  Each error case stops the program, so each runs as
its own build.  A LOGICAL stored into the REAL*8 block fills the first
four bytes of its word, read back through an EQUIVALENCEd LOGICAL array.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'namer.json')
ROUTINES = ['namer', 'readcd', 'skipbl', 'findch', 'extrst', 'toint',
            'todec', 'reptct', 'findvn', 'tolog', 'sublog', 'subint',
            'subrea']
# The caret became a carriage return (``4H<CR>   /``), which reads as a line
# break.
PATCH = {'namer': [('DATA CARET / 4H\n   /', 'DATA CARET / 4H^   /')],
         # READCD's retry after END stops a gfortran run (Read past
         # ENDFILE); one attempt gives the end of file the source expects.
         'readcd': [('DO 1000 J=1,2', 'DO 1000 J=1,1')]}

# XNAM1's FLTCON: names, dimensions, locations; a 161-word block.
NAMES = ['NMACH', 'MACH', 'NALPHA', 'ALSCHD', 'RNNUB', 'HYPERS', 'STMACH',
         'TSMACH', 'TR', 'ALT', 'PINF', 'TINF', 'VINF', 'WT', 'GAMMA',
         'LOOP', 'ALPHA', 'GRDHT']
LEN = [5, 4, 6, 6, 5, 6, 6, 6, 2, 3, 4, 4, 4, 2, 5, 4, 5, 5]
LDM = [1, 20, 1, 20, 20, -1, 1, 1, 1, 20, 20, 20, 20, 1, 1, 1, 1, 20]
LOC = [1, 3, 2, 23, 43, 161, 94, 95, 96, 97, 74, 117, 137, 157, 158, 159,
       160, 0]
MARK = 0.5

GOOD = [
    [' $FLTCON NMACH=2.0, MACH=0.6,0.8, ALSCHD(3)=1.0,2.0,$'],
    [' $FLTCON NMACH=1.0,', '   MACH(1)=0.5, RNNUB=2*1.E6,',
     '   3.E6, HYPERS=.TRUE.$'],
    [' $OPTINS SREF=1.0$', ' $FLTCON LOOP=2.0, WT=1.5E4$'],
    [' $OPTINS SREF=1.0$'],
    [' $FLTCONX NMACH=3.0, ALT=5*1000.$'],
    ['', '  $FLTCON   GAMMA = 1.4 ,  TR= 1.0,HYPERS=.F.$'],
]
BAD = [
    [' $ NMACH=1.$'],
    [' $FLTCON NMACH 1.$'],
    [' $FLTCON MACH(2=1.$'],
    [' $FLTCON XYZ=1.$'],
    [' $FLTCON MACH(21)=1.$'],
    [' $FLTCON MACH(19)=3*1.$'],
    [' $FLTCON NMACH=ABC$'],
    [' $FLTCON MACH(0)=1.$'],
    [' $FLTCON HYPERS=1.0$'],
    [' $FLTCON RNNUB(20)=1.0$', ''],
    [' $FLTCON TINF(20)=2*5.$'],
    # A blank card between value continuations: the values after it are
    # read as a variable name.
    [' $FLTCON ALSCHD=-2.,0.,', '', '   2.,4.,  GRDHT=1.,2.$'],
]


def _card_lines(cards):
    return [c.ljust(80)[:80] for c in cards]


def driver(files, stop_ok=False) -> str:
    packed = ''.join(NAMES)
    lines = ['      PROGRAM PROBE',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND',
             f'      INTEGER VN({len(packed)}),LEN(18),LDM(18),LOC(18),'
             'NL(6),KD',
             '      DIMENSION CB(161)', '      LOGICAL LB(322),EOF',
             '      EQUIVALENCE (CB,LB)',
             f'      CHARACTER*{len(packed)} S', '      CHARACTER*80 C',
             "      KD=TRANSFER('$   ',KD)", '      KAND=KD',
             f"      S(1:40)='{packed[:40]}'",
             f"      S(41:)='{packed[40:]}'", f'      DO 1 J=1,{len(packed)}',
             "    1 VN(J)=TRANSFER(S(J:J)//'   ',VN(J))",
             "      S='FLTCON'", '      DO 2 J=1,6',
             "    2 NL(J)=TRANSFER(S(J:J)//'   ',NL(J))"]
    for k in range(18):
        lines += [f'      LEN({k + 1})={LEN[k]}', f'      LDM({k + 1})={LDM[k]}',
                  f'      LOC({k + 1})={LOC[k]}']
    for n, cards in enumerate(files):
        unit = 20 + n
        lines.append(f"      OPEN({unit},FILE='c{n}.txt',STATUS='UNKNOWN')")
        for c in _card_lines(cards):
            lines.append(f"      C='{c[:50]}'")
            lines.append(f"      C(51:80)='{c[50:]}'")
            lines.append(f"      WRITE({unit},'(A)') C")
        lines.append(f'      REWIND {unit}')
    for n, _ in enumerate(files):
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        lines += [f'      DO {100 + n} K=1,161', f'{100 + n:5d} CB(K)={MARK}D0']
        lines.append(f'      CALL NAMER(KAND,{20 + n},NL,6,VN,{len(packed)},'
                     'LEN,18,LDM,CB,161,')
        lines.append('     1           LOC,EOF)')
        lines.append("      WRITE(6,'(A)') '@@DATA'")
        lines.append("      WRITE(6,'(L2)') EOF")
        lines.append("      WRITE(6,'(161ES25.16)') CB")
        lines.append("      WRITE(6,'(L2)') LB(321)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def _parse(out):
    cases = []
    for chunk in out.split('@@CASE\n')[1:]:
        if '@@DATA\n' in chunk:
            text, data = chunk.split('@@DATA\n')
            rows = data.strip('\n').split('\n')
            cases.append({'text': [ln for ln in text.split('\n') if ln],
                          'eof': rows[0].strip() == 'T',
                          'block': [float(v) for v in rows[1].split()],
                          'hypers': rows[2].strip() == 'T'})
        else:
            cases.append({'text': [ln for ln in chunk.split('\n') if ln]})
    return cases


def main():
    good = _parse(probe.run('namer', driver(GOOD), ROUTINES,
                            layout_patches=PATCH))
    bad = []
    for n, cards in enumerate(BAD):
        bad += _parse(probe.run(f'namer_e{n}', driver([cards]), ROUTINES,
                                layout_patches=PATCH))
    FIXTURE.write_text(json.dumps({
        'names': NAMES, 'len': LEN, 'ldm': LDM, 'loc': LOC, 'mark': MARK,
        'good': GOOD, 'bad': BAD, 'good_records': good,
        'bad_records': bad}, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
