"""
Probe VTAREA (with PTINT1 and AREA1), the vertical-panel area in the Mach
shadows, and save the fixture.

Run from the repository root: ``python test_parity/probes/vtarea.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['vtarea', 'ptint1', 'area1']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(100)
      COMMON /HTI/    HTIN(154)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(73)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC
      COMMON /VTDATA/ AVT(195), AVF(195)
      COMMON /VTI/    VTIN(154), TVTIN(8), VFIN(154)
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,VERTUP
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

MACH_INDEX = 2
VT_COMMON = {59: math.atan(0.8), 77: 0.2, 83: math.atan(1.0), 101: 0.3}


def case(mach=2.0, vertup=True, htpl=True, cranked=False, ventral=False,
         xv=28.0, zv=1.0, aliw=2.0, alih=0.0, xh=30.0, zh=1.0):
    if cranked:
        vtin = {1: 1.2, 2: 2.0, 3: 5.2, 4: 6.0, 5: 3.0, 15: 3.0}
        avt = {1: 12.8, 2: 4.2, 3: 17.0, 10: 5.0, 21: 4.0, 62: 0.8,
               86: 1.0}
    else:
        vtin = {1: 2.0, 3: 5.2, 4: 6.0, 5: 2.0, 15: 1.0}
        avt = {1: 18.2, 2: 0.0, 3: 18.2, 10: 5.0, 21: 6.0, 62: 0.8,
               86: 0.0}
    return {
        'mach': mach, 'vertup': vertup, 'htpl': htpl, 'ventral': ventral,
        'xv': xv, 'zv': zv, 'vtin': vtin, 'avt': avt,
        'wing': {'span': 15.0, 'spans': 13.0, 'a62': 0.6, 'a10': 6.0},
        'tail': {'span': 5.0, 'spans': 4.0, 'a10': 3.0, 'a16': 2.5,
                 'a30': 0.5, 'a62': 0.7},
        'syna': {2: 10.0, 3: 0.0, 4: aliw, 6: xh, 7: zh, 8: alih},
        'vt_common': VT_COMMON if not ventral else
        {59: 0.5, 77: 0.1, 83: 0.6, 101: 0.2},
        'stale134': 0.77,
    }


def cases():
    return [
        case(),
        case(mach=1.3),
        case(mach=3.5),
        case(mach=1.15, xv=20.0),
        case(vertup=False, zv=-0.5),
        case(htpl=False),
        case(cranked=True),
        case(cranked=True, mach=1.4),
        case(cranked=True, vertup=False, zv=-0.5, mach=2.5),
        case(ventral=True, vertup=False, zv=-0.8),
        case(alih=-3.0),                        # tail turned to wing
        case(xv=14.0, mach=1.6),
        case(xv=8.0, mach=1.25, zh=3.0),
        case(xv=40.0, mach=5.0, aliw=-4.0),
        case(xv=24.0, mach=1.05, zv=4.0),
        # Found by sampling: every shadow shape, above and below the body.
        case(mach=1.2, vertup=False, xv=2.84, zv=-3.82, aliw=4.18,
             zh=-0.22, alih=-3.0),
        case(mach=1.5, vertup=False, cranked=True, xv=12.58, zv=-3.25,
             aliw=3.54, zh=-3.33, alih=-3.0),
        case(mach=1.05, xv=13.97, zv=4.39, aliw=3.9, zh=-0.67, alih=-3.0),
        case(mach=3.0, vertup=False, xv=13.65, zv=0.876, aliw=6.12,
             zh=2.77, alih=-3.0),
        case(mach=1.5, vertup=False, cranked=True, xv=6.93, zv=-3.03,
             aliw=5.23, zh=0.878, alih=-3.0),
        case(mach=1.05, xv=26.85, zv=4.88, aliw=2.81, zh=-1.56, alih=-3.0),
        case(mach=2.0, vertup=False, cranked=True, xv=2.82, zv=-3.31,
             aliw=-3.33, zh=0.33),
        case(mach=1.2, cranked=True, xv=9.76, zv=-2.64, aliw=-7.35,
             zh=0.125),
        case(mach=1.2, vertup=False, cranked=True, xv=5.16, zv=-3.7,
             aliw=1.39, zh=-3.02, alih=-3.0),
        case(mach=1.05, xv=13.3, zv=0.966, aliw=-0.8, zh=-1.49),
        # A ventral fin whose shadow turns on /VTDATA/'s sweeps.
        case(mach=1.05, cranked=True, xv=13.23, zv=2.69, aliw=5.96,
             zh=-3.65, ventral=True),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             f'      I={MACH_INDEX}']
    for n, c in enumerate(all_cases):
        vin, av = ('VFIN', 'AVF') if c['ventral'] else ('VTIN', 'AVT')
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,154')
        lines.append('         VTIN(K)=0.')
        lines.append('         VFIN(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(assign('FLC(4)', c['mach']))
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(f"      VERTUP=.{'TRUE' if c['vertup'] else 'FALSE'}.")
        for k, v in c['vtin'].items():
            if k == 15:
                lines.append(f'      {vin}(15)=WT({int(v)})')
            else:
                lines.append(assign(f'{vin}({k})', v))
        lines.append(assign(f'{vin}({134 + MACH_INDEX})', c['stale134']))
        for k, v in c['avt'].items():
            lines.append(assign(f'{av}({k})', v))
        for k, v in c['vt_common'].items():
            lines.append(assign(f'AVT({k})', v))
        w, t = c['wing'], c['tail']
        lines.append(assign('WINGIN(4)', w['span']))
        lines.append(assign('WINGIN(3)', w['spans']))
        lines.append(assign('A(62)', w['a62']))
        lines.append(assign('A(10)', w['a10']))
        lines.append(assign('HTIN(4)', t['span']))
        lines.append(assign('HTIN(3)', t['spans']))
        for k in (10, 16, 30, 62):
            lines.append(assign(f'AHT({k})', t[f'a{k}']))
        for k, v in c['syna'].items():
            lines.append(assign(f'SYNA({k})', v))
        lines.append(assign('XV', c['xv']))
        lines.append(assign('ZV', c['zv']))
        lines.append(f'      CALL VTAREA({vin},{av},VERTUP,XV,ZV)')
        m = MACH_INDEX
        lines.append("      WRITE(6,'(A,5ES25.16)') 'R',")
        lines.append(f'     1{vin}({94 + m}),{vin}({134 + m}),'
                     f'{vin}({114 + m}),SYNA(6),SYNA(7)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('vtarea', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('vtin', 'avt', 'syna', 'vt_common'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('vtarea', payload))


if __name__ == '__main__':
    main()
