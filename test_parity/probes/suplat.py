"""
Probe SUPLAT and SUPLAH, the supersonic sideslip derivatives of the wing,
tail and body, and save the fixture.

Run from the repository root: ``python test_parity/probes/suplat.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['suplat', 'suplah', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad', 'trapz']

COMMONS = """\
      COMMON /OPTION/ SR,CRBAR,ROUGFC,BLREF
      COMMON /FLGTCD/ FLC(160)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /BDATA/  BD(762)
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF
      COMMON /SYNTSS/ SYNA(19)
      COMMON /WINGD/  A(195)
      COMMON /SUPWH/  SLG(282)
      COMMON /SBETA/  SLA(62)
      COMMON /SUPBOD/ SBD(227)
      COMMON /HTDATA/ AHT(195)
      COMMON /BODYI/  BODYIN(126)
      COMMON /WINGI/  WINGIN(100)
      COMMON /HTI/    HTIN(154)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IBW/    PBW, BWI(380)
      COMMON /IBH/    PBH, BH(380)
      COMMON /IBWH/   PBWH, BWH(380)
      COMMON /IBWV/   PBWV, BWV(380)
      COMMON /IBWHV/  PBWHV, BWHV(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

UNUSED = 1.0e-30
I = 1
ALPHA = [-4.0, 0.0, 4.0, 8.0]
BX = [0.0, 3.0, 8.0, 15.0, 25.0, 35.0, 40.0]
BR = [0.0, 0.9, 1.5, 1.8, 1.8, 1.6, 1.4]


def case(tail=False, mach=2.0, taper=0.3, tanle=0.5, ar=3.0, planform=1,
         bo=True, htpl=True, zh=0.5, dihedral=(UNUSED, 3.0, UNUSED)):
    le = math.atan(tanle)
    return {
        'tail': tail, 'alpha': ALPHA, 'mach': mach, 'rl': 2.0e6, 'i': I,
        'sref': 250.0, 'cbarr': 8.0, 'blref': 30.0, 'bo': bo, 'htpl': htpl,
        'syna': {1: 20.0, 2: 12.0, 3: -0.4, 4: 2.0, 6: 34.0, 7: zh, 8: -1.0},
        'win': {3: 13.0, 4: 15.0, 6: 9.0, 12: dihedral[0], 13: dihedral[1],
                14: dihedral[2], 15: planform},
        'a': {10: 8.5, 59: le, 61: math.cos(le), 62: tanle,
              74: 0.6 * tanle, 118: taper, 120: ar},
        'clpcty': 1.1, 'cnaw': 0.05, 'cn': [0.05 * a for a in ALPHA],
        'clab': 0.004, 'bd': 0.3,
        'body': {'x': BX, 'r': BR},
        'htin': {3: 4.0, 4: 5.5, 94 + I: 18.0, 114 + I: 20.0,
                 134 + I: 25.0},
        'aht': {30: 0.4, 62: 0.7},
    }


def cases():
    root3 = math.sqrt(3.0)
    return [
        case(),                                          # tapered
        case(taper=1.0, tanle=0.0),                      # rectangular
        case(taper=0.0, tanle=2.0, ar=1.5),              # delta
        case(planform=3),                                # not straight
        case(bo=False),
        case(zh=3.0),                                    # tail off body
        case(zh=0.0, dihedral=(4.0, 3.0, 6.0)),          # tail on the axis
        case(tail=True, taper=0.0, tanle=root3, ar=1.5),  # beta/tan = 1
        case(tail=True),
        case(tail=True, tanle=-0.3),                     # swept forward
        case(taper=0.0, tanle=root3, ar=1.5),            # SUPLAT at 1.0
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', f'      I={I}']
    for n, c in enumerate(all_cases):
        t = c['tail']
        a, inn, blk, body, slg, sla0 = (
            ('AHT', 'HTIN', 'HT', 'BH', 141, 31) if t else
            ('A', 'WINGIN', 'WING', 'BWI', 0, 0))
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,380')
        for arr in ('WING', 'HT', 'BWI', 'BH', 'BWH', 'BWV', 'BWHV'):
            lines.append(f'         {arr}(K)=0.')
        lines.append('         IF(K.LE.62) SLA(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, al in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', al))
            lines.append(assign(f'{blk}({61 + j})', c['cn'][j]))
        lines.append(assign(f'FLC({I + 2})', c['mach']))
        lines.append(assign(f'FLC({I + 42})', c['rl']))
        for name, key in (('SR', 'sref'), ('CRBAR', 'cbarr'),
                          ('BLREF', 'blref'), ('SBD(18)', 'clab')):
            lines.append(assign(name, c[key]))
        lines.append(assign('BD(89)' if t else 'BD(66)', c['bd']))
        lines.append(f"      BO=.{'TRUE' if c['bo'] else 'FALSE'}.")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        for k, v in c['syna'].items():
            lines.append(assign(f'SYNA({k})', v))
        for k, v in c['win'].items():
            if k == 15:
                lines.append(f'      {inn}(15)=WT({v})')
            else:
                lines.append(assign(f'{inn}({k})', v))
        for k, v in c['a'].items():
            lines.append(assign(f'{a}({k})', v))
        lines.append(assign(f'SLG({slg + 3})', c['clpcty']))
        lines.append(assign(f'SLG({slg + 7})', c['cnaw']))
        lines.append(f"      BODYIN(1)={len(c['body']['x'])}.")
        for k, (x, r) in enumerate(zip(c['body']['x'], c['body']['r'])):
            lines.append(assign(f'BODYIN({2 + k})', x))
            lines.append(assign(f'BODYIN({62 + k})', r))
        if not t:
            for k, v in c['htin'].items():
                lines.append(assign(f'HTIN({k})', v))
            for k, v in c['aht'].items():
                lines.append(assign(f'AHT({k})', v))
        lines.append('      LF=-1')
        lines.append(f"      CALL {'SUPLAH' if t else 'SUPLAT'}")
        lines.append("      WRITE(6,'(A,31ES25.16)') 'SLA',")
        lines.append(f'     1(SLA(K),K={sla0 + 1},{sla0 + 31})')
        for tag, arr in (('SURF', blk), ('BODY', body), ('BWH', 'BWH'),
                         ('BWV', 'BWV'), ('BWHV', 'BWHV')):
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{tag}',")
            lines.append(f'     1({arr}(K),K=141,200)')
        lines.append(f"      WRITE(6,'(A,5ES25.16)') 'W',{inn}(12),"
                     f"{inn}(13),{inn}(14),BD(1)")
        lines.append("      WRITE(6,'(A,I4)') 'LF',LF")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('suplat', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('syna', 'win', 'a', 'htin', 'aht'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('suplat', payload))


if __name__ == '__main__':
    main()
