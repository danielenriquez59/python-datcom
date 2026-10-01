"""
Probe SUPLAV and SUPLAF (with MASRAT), the supersonic vertical-panel
sideslip increments, and save the fixture.

Run from the repository root: ``python test_parity/probes/suplav.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['suplav', 'suplaf', 'masrat', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /IVT/    PVT, VT(380)
      COMMON /IVF/    PVF, VF(380)
      COMMON /IBW/    PBW, BWI(380)
      COMMON /IBWV/   PBWV, BWV(380)
      COMMON /IBWHV/  PBWHV, BWHV(380)
      COMMON /FLGTCD/ FLC(160)
      COMMON /OPTION/ SREF, CBARR, ROUGFC, BLREF
      COMMON /SYNTSS/ XCG, XW, ZW, ALIW, ZCG, XH, ZHH, ALIH, XV,
     1                VERTUP, HINAX, XVF, SC, ZV, ZVF, YV, YF,
     2                PHIV, PHIF
      COMMON /WINGI/  WINGIN(101)
      COMMON /VTI/    VTIN(154), TVTIN(8), VFIN(154)
      COMMON /HTI/    HTIN(154)
      COMMON /SBETA/  SLA(62)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /VTDATA/ AVT(195), AVF(195)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,K
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
      LOGICAL VERTUP
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

ALPHA = [-4.0, 0.0, 4.0, 8.0]
MACH_INDEX = 2


def case(ventral=False, zw=0.0, zh=1.5, vertup=True, htpl=True,
         straight=True, alih=0.0):
    return {
        'alpha': ALPHA, 'ventral': ventral, 'htpl': htpl,
        'straight': straight, 'blref': 22.0, 'span': 10.0, 'cnav': 0.045,
        'vtin': {3: 4.5, 4: 6.0},
        'sv': {'svwb': 0.8, 'svb': 0.6, 'svhb': 0.5},
        'avt': {3: 30.0, 30: 1.4, 31: 2.2, 62: 0.9},
        'htin': {3: 3.0, 4: 4.0}, 'aht': {30: 0.7, 62: 0.6},
        'position': {'xcg': 20.0, 'zw': zw, 'zcg': 0.4, 'zh': zh,
                     'alih': alih, 'xv': 38.0, 'zv': 0.9,
                     'vertup': vertup},
        'bwv': {141: 0.01, 161: -0.002,
                **{181 + j: 0.001 * j for j in range(len(ALPHA))}},
        'bwhv': {141: 0.02, 161: -0.003,
                 **{181 + j: 0.002 * j for j in range(len(ALPHA))}},
    }


def cases():
    out = []
    for ventral in (False, True):
        out += [
            case(ventral),                                   # mid, same
            case(ventral, zw=1.5, zh=-0.6),                  # opposite, mid
            case(ventral, zw=-0.9, zh=3.0, vertup=False),    # clamp
            case(ventral, zw=-1.5, zh=-1.5),
            case(ventral, htpl=False, straight=False),
            case(ventral, straight=False, zw=0.3),
            case(ventral, alih=3.0, zw=0.45),
            case(ventral, htpl=False, zw=0.7),
        ]
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             f'      I={MACH_INDEX}']
    for n, c in enumerate(all_cases):
        v = c['ventral']
        blk, vin, av, sla0 = (('VF', 'VFIN', 'AVF', 31) if v else
                              ('VT', 'VTIN', 'AVT', 0))
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,380')
        for arr in ('VT', 'VF', 'BWV', 'BWHV'):
            lines.append(f'         {arr}(K)=0.')
        lines.append('         IF(K.LE.62) SLA(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        p = c['position']
        lines.append(f"      VERTUP=.{'TRUE' if p['vertup'] else 'FALSE'}.")
        for name, key in [('XCG', 'xcg'), ('ZW', 'zw'), ('ZCG', 'zcg'),
                          ('ZHH', 'zh'), ('ALIH', 'alih'), ('XV', 'xv'),
                          ('ZV', 'zv')]:
            lines.append(assign(name, p[key]))
        lines.append(assign('BLREF', c['blref']))
        lines.append(assign('WINGIN(4)', c['span']))
        lines.append(f"      WINGIN(15)=WT({1 if c['straight'] else 3})")
        lines.append(assign(f'SLA({sla0 + 31})', c['cnav']))
        for k, val in c['vtin'].items():
            lines.append(assign(f'{vin}({k})', val))
        for key, start in [('svwb', 95), ('svb', 115), ('svhb', 135)]:
            lines.append(assign(f'{vin}({start + MACH_INDEX - 1})',
                                c['sv'][key]))
        for k, val in c['avt'].items():
            lines.append(assign(f'{av}({k})', val))
        for k, val in c['htin'].items():
            lines.append(assign(f'HTIN({k})', val))
        for k, val in c['aht'].items():
            lines.append(assign(f'AHT({k})', val))
        for k, val in c['bwv'].items():
            lines.append(assign(f'BWV({k})', val))
        for k, val in c['bwhv'].items():
            lines.append(assign(f'BWHV({k})', val))
        lines.append(f"      CALL {'SUPLAF' if v else 'SUPLAV'}")
        lines.append("      WRITE(6,'(A,62ES25.16)') 'SLA',(SLA(K),K=1,62)")
        for arr in (blk, 'BWV', 'BWHV'):
            tag = 'P' if arr == blk else arr
            lines.append(f"      WRITE(6,'(A,50ES25.16)') '{tag}',")
            lines.append(f'     1({arr}(K),K=141,190)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('suplav', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('vtin', 'avt', 'htin', 'aht', 'bwv', 'bwhv'):
            c[key] = {str(k): val for k, val in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('suplav', payload))


if __name__ == '__main__':
    main()
