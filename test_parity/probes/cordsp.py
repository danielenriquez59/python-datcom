"""
Probe CORDSP, the supersonic section coordinates, and save the fixture.

Run from the repository root: ``python test_parity/probes/cordsp.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['cordsp']

COMMONS = """\
      COMMON /IWING/  PW,X(60)
      COMMON /IHT/    PHT,XU(60),XL(60),YUU(60),YLL(60)
      COMMON /IVT/    PVT,YU(60),YL(60)
      COMMON /IBW/    PBW,L,I,J,K,II,JJ,KK,III,JJJ,KKK,LLL
      COMMON /WINGI/  WGIN(100)
      COMMON /HTI/    HTIN(154)
      COMMON /VTI/    VTIN(154),TVTIN(8),VFIN(154)
      COMMON /IBWHV/  PBWHV,RHO,TOC
      COMMON /IBODY/  PB,NACA(80)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DIMENSION WT(4)
      INTEGER SURF(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
      DATA SURF /4HW   ,4HH   ,4HV   ,4HF   /
"""

X = [0.0, .025, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6,
     .65, .7, .75, .8, .85, .9, .95, 1.0]
UNUSED = 1.0e-30
BLOCK = {'W': 'WGIN', 'H': 'HTIN', 'V': 'VTIN', 'F': 'VFIN'}


def case(surface='W', planform=1, **digits):
    d = {k: digits.get(k, 0) for k in ('i', 'j', 'k', 'ii', 'jj', 'kk',
                                       'iii', 'jjj', 'kkk', 'lll')}
    words = {15: planform, 16: UNUSED, 18: UNUSED, 62: UNUSED, 63: UNUSED,
             70: UNUSED, 71: 3.0}
    return {'surface': surface, 'digits': d, 'x': X + [0.0] * 38,
            'surface_in': words}


def cases():
    return [
        case(i=1, j=3, jj=0, kk=5),                           # wedge 30%
        case('H', i=2, jj=0, kk=6),                           # biconvex
        case('V', i=3, j=2, k=5, jj=0, kk=4, jjj=3),          # hexagonal
        case('F', i=1, j=5, jj=0, kk=3, iii=5, planform=3),
        case(i=4, jj=0, kk=4),                                # falls through
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        blk = BLOCK[c['surface']]
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      NACA(6)=SURF({'WHVF'.index(c['surface']) + 1})")
        for m, v in enumerate(c['x']):
            lines.append(assign(f'X({m + 1})', v))
        for name, v in c['digits'].items():
            lines.append(f'      {name.upper()}={v}')
        for k, v in c['surface_in'].items():
            if k == 15:
                lines.append(f'      {blk}(15)=WT({v})')
            else:
                lines.append(assign(f'{blk}({k})', v))
        lines.append('      CALL CORDSP')
        for arr in ('XU', 'XL', 'YU', 'YL', 'YUU', 'YLL'):
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{arr}',{arr}")
        lines.append(f"      WRITE(6,'(A,8ES25.16)') 'S',RHO,TOC,{blk}(16),")
        lines.append(f'     1{blk}(18),{blk}(62),{blk}(63),{blk}(70),'
                     f'{blk}(71)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('cordsp', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c['surface_in'] = {str(k): v for k, v in c['surface_in'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('cordsp', payload))


if __name__ == '__main__':
    main()
