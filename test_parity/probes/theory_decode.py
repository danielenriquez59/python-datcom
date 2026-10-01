"""
Probe TOLOG, REPTCT and THEORY and save their raw records.

Run from the repository root: ``python test_parity/probes/theory_decode.py``.
THEORY runs with IDEAL and SLOPE stubbed: SLOPE returns
``CLA = 0.1 + 0.02*MACH`` and ``XAC = 0.25 + 0.01*MACH``, or ``UNUSED``
above Mach 0.7, so the failure path is reached.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'theory_decode.json')
U = 1.0e-30

LOGS = ['.TRUE.', '  .T.', '.FALSE.,', '.F.', 'TRUE', '.TRU', '.TRUEX',
        '.FALSEY', '.FALS.', '', '   .T', '  .F.  ']
REPTS = ['3*1.5', '12*.TRUE.', '*7', '0*2.', '2.5*3', 'NOSTAR', '1*',
         ' 4 *9.']
THEORY = [
    {'surface': 'W', 'machs': [0.2, 0.5, 0.8], 'cla': [U, 0.11, U],
     'mcc': 0.9},
    {'surface': 'H', 'machs': [0.3, 0.75, 0.6], 'cla': [U, U, U],
     'mcc': 0.72},
    {'surface': 'X', 'machs': [0.4], 'cla': [U], 'mcc': 0.95},
]
SECTION = {'ai': 1.25, 'alo': -2.5, 'cli': 0.15, 'cmco4': -0.045,
           'rho': 0.0158, 'tmax': 0.12, 'deltay': 2.9, 'xc': 0.42,
           'clcc': 0.123}
SLOPE_STUB = """\
      SUBROUTINE IDEAL
      RETURN
      END
      SUBROUTINE SLOPE(MACH,CLA,RENN,XAC)
      REAL MACH
      CLA=0.1D0+0.02D0*MACH
      XAC=0.25D0+0.01D0*MACH
      IF(MACH.GT.0.7D0) CLA=1.D-30
      RETURN
      END
"""


def _real(v):
    s = repr(float(v))
    return s.replace('e', 'D') if 'e' in s else s + 'D0'


def _load(lines, text, label):
    card = text.ljust(80)[:80]
    lines.append(f"      CARD='{card[:40]}'")
    lines.append(f"      CARD(41:80)='{card[40:]}'")
    lines.append(f'      DO {label} J=1,80')
    lines.append(f"{label:5d} KOL(J)=TRANSFER(CARD(J:J)//'   ',KOL(J))")


def decode_driver():
    lines = ['      PROGRAM PROBE', '      INTEGER KOL(80)',
             '      LOGICAL LANS', '      CHARACTER*80 CARD']
    for n, t in enumerate(LOGS):
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        _load(lines, t, 100 + n)
        lines.append('      CALL TOLOG(KOL,LANS,IERR)')
        lines.append("      WRITE(6,'(L2,I3)') LANS,IERR")
    for n, t in enumerate(REPTS):
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        _load(lines, t, 200 + n)
        lines.append('      CALL REPTCT(KOL,IREPT,IERR)')
        lines.append("      WRITE(6,'(2I6,1X,80A1)') IREPT,IERR,KOL")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def theory_driver():
    s = SECTION
    lines = ['      PROGRAM PROBE',
             '      COMMON /IBODY/  PB, NACA(80), BF(232), CBAR',
             '      COMMON /FLGTCD/ FLC(93)',
             '      COMMON / IBH /  PBH, THN(60),CAM(60),A0,XC,MCC,CLCC,'
             'XAC(20)',
             '      COMMON /IBWH/   PBWH,AI,ALO,CLI,ASEP,CMCO4,CLA0,CLA(20)',
             '      COMMON /IBWHV/  PBWHV, RHO,TMAX,DELTAY',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      REAL MCC', '      INTEGER NACA',
             '      UNUSED=1.D-30', '      CBAR=4.0D0']
    for k in ('ai', 'alo', 'cli', 'cmco4', 'rho', 'tmax', 'deltay', 'xc',
              'clcc'):
        lines.append(f'      {k.upper()}={_real(s[k])}')
    for c in THEORY:
        lines.append("      WRITE(6,'(A)') '@@CASE'")
        lines.append(f"      NACA(6)=TRANSFER('{c['surface']}   ',NACA(6))")
        lines.append(f"      MCC={_real(c['mcc'])}")
        lines.append(f"      FLC(1)={_real(len(c['machs']))}")
        for i, m in enumerate(c['machs']):
            lines.append(f'      FLC({i + 3})={_real(m)}')
            lines.append(f'      FLC({i + 43})={_real(1.0e6 * (i + 1))}')
            lines.append(f"      CLA({i + 1})={_real(c['cla'][i])}")
            lines.append(f'      XAC({i + 1})=-1.0D0')
        lines.append('      CALL THEORY')
        n = len(c['machs'])
        lines.append("      WRITE(6,'(A)') '@@VALUES'")
        lines.append(f"      WRITE(6,'({2 * n + 1}ES25.16)') CLA0,"
                     f"(CLA(K),K=1,{n}),(XAC(K),K=1,{n})")
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
    decode = _records(probe.run('decode', decode_driver(),
                                ['tolog', 'reptct', 'skipbl', 'findch',
                                 'extrst', 'toint', 'todec']))
    # CALL EXIT binds to gfortran's intrinsic, which ends the run; the
    # build copy calls a stub that records it instead.
    stubs = probe.STUBS + SLOPE_STUB + (
        '      SUBROUTINE EXITX\n'
        "      WRITE(6,'(A)') '@@EXIT'\n      RETURN\n      END\n")
    theory = _records(probe.run(
        'theory', theory_driver(), ['theory'], stubs=stubs,
        layout_patches={'theory': [('CALL EXIT', 'CALL EXITX')]}))
    FIXTURE.write_text(json.dumps({
        'logs': LOGS, 'repts': REPTS, 'theory': THEORY, 'section': SECTION,
        'decode': decode, 'records': theory}, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
