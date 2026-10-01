"""
Probe DECODE and AIRFOL (with the section routines, CORDSP and XYCORD)
and save the fixture.

Run from the repository root: ``python test_parity/probes/airfol.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import parse_records, run, save  # noqa: E402

ROUTINES = ['airfol', 'decode', 'coord1', 'coord4', 'coord5', 'coord6',
            'cord4m', 'cord5m', 'cordsp', 'xycord', 'sleq', 'arccos',
            'tbfunx', 'quad']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /IBODY/  PB,NACA(80)
      COMMON /IWING/  PW, X(60)
      COMMON / IHT /  PHT, XU(60),XL(60),YUU(60),YLL(60)
      COMMON / IVT /  PVT, YUN(60),YLN(60)
      COMMON / IBW /  PBW,L,I,J,K,II,JJ,KK,III,JJJ,KKK,LLL
      COMMON / IBH /  PBH, THN(60),CAM(60)
      COMMON /IBWHV/  PBWHV, RHO,T,DELTAY,XOVC,TOVC,ZM,ZP
      COMMON /FLOLOG/ DUM(29),PART
      COMMON /WINGI/  WGIN(100)
      COMMON /HTI/    HTIN(154)
      COMMON /VTI/    VTIN(154),TVTIN(8),VFIN(154)
      LOGICAL DUM,PART
      CHARACTER*80 CARD
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

CARDS = [
    'NACA-W-4-2412',
    'NACA-H-4-0012',
    'NACA-W-4-2412-63',
    'NACA-W-5-23012',
    'NACA-W-5-23112',
    'NACA-W-5-23012-33',
    'NACA-W-1-16-212',
    'NACA-W-6-65-212',
    'NACA-W-6-64A010',
    'NACA-W-6-63,2-415',
    'NACA-W-6-65-218 A=0.5',
    'NACA-W-S-1-30-5.0',
    'NACA-W-S-2-00-4.0',
    'NACA-W-S-3-25-4.5-30',
    'NACA-W-4-4415',
]


def cases():
    return [{'card': c} for c in CARDS]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      PART=.FALSE.']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} M=1,60')
        for arr in ('XU', 'XL', 'YUN', 'YLN', 'THN', 'CAM', 'X', 'YUU',
                    'YLL'):
            lines.append(f'         {arr}(M)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      CARD='{c['card']}'")
        lines.append("      READ(CARD,'(80A1)') NACA")
        for k in (16, 18, 62, 63, 70, 71):
            lines.append(f'      WGIN({k})=UNUSED')
            lines.append(f'      HTIN({k})=UNUSED')
        lines.append('      WGIN(15)=WT(1)')
        lines.append('      HTIN(15)=WT(1)')
        lines.append('      NA=0')
        lines.append('      CALL DECODE(NACA,NA)')
        lines.append("      WRITE(6,'(A,12I6)') 'D',NA,L,I,J,K,II,JJ,KK,"
                     "III,JJJ,")
        lines.append('     1KKK,LLL')
        lines.append('      CALL AIRFOL')
        for arr in ('X', 'XU', 'XL', 'YUN', 'YLN', 'THN', 'CAM'):
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{arr}',")
            lines.append(f'     1({arr}(M),M=1,L)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('airfol', driver(all_cases), ROUTINES))
    print(save('airfol', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
