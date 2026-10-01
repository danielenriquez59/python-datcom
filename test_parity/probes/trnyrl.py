"""
Probe M40O50 (TRNYRL), the transonic flap increments, and save the
fixture.

Run from the repository root: ``python test_parity/probes/trnyrl.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m40o50', 'trnyrl', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA,IG,IJKDUM(3),NOVLY
      COMMON /FLGTCD/ FLC(93)
      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF
      COMMON /FLAPIN/ F(69)
      COMMON /HTI/    HTIN(131)
      COMMON /IBODY/  PBODY, BODY(400)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IHT/    PHT, HT(380)
      COMMON /WINGD/  A(195)
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /SBETA/  SB(351)
      COMMON /SUPDW/  DWA(35),TCD(58)
      COMMON /CONSNT/ PI,DR,UNUSED,RADS
      COMMON /POWR/   PW(300)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2        SUPERS,SUBSON,TRANSN,HYPERS,
     3        SYMFP,ASYFP,TRIMC,TRIM
"""

UNUSED = 1.0e-30
NALPHA = 3


def case(asyfp=False, htpl=False, ftype=2.0, stype=1.0, mach=0.9,
         sspne=4.0):
    nd = 4
    cnym6 = [0.001 * (k + 1) * (-1) ** k for k in range(nd * NALPHA)]
    cnym6[5] = UNUSED
    return {
        'mach': mach, 'nalpha': NALPHA, 'asyfp': asyfp, 'htpl': htpl,
        'flap': {'ndelta': float(nd), 'type': ftype, 'stype': stype,
                 'deltal': [5.0, 10.0, 15.0, 20.0],
                 'deltar': [-5.0, -8.0, -12.0, -15.0]},
        'wing_cla': 0.075, 'tail_cla': 0.05, 'tra70': 0.068,
        'trah70': 0.047, 'delcl6': [0.1, 0.2, 0.3, 0.38],
        'cfact': [0.05, 0.1, 0.12, 0.15],
        'clrlm6': [0.002, 0.004, 0.0055, 0.007], 'cnym6': cnym6,
        'bd': [0.5, 0.45, 0.3, 0.12], 'sspn': 5.0, 'sspne': sspne,
        'astrw': 3.5, 'depsda': 0.35, 'blref': 30.0,
    }


def cases():
    return [
        case(),
        case(htpl=True, ftype=1.0),
        case(asyfp=True, stype=1.0),
        case(asyfp=True, stype=4.0),
        case(asyfp=True, stype=5.0, mach=0.8),
        case(asyfp=True, stype=5.0, mach=1.2),
        case(asyfp=True, stype=5.0, mach=0.95, sspne=4.6),
        case(asyfp=True, stype=5.0, mach=0.95, sspne=2.5),
        case(ftype=6.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RADS=57.2957795', '      IM=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,300')
        for arr in ('WING', 'BODY', 'HT', 'PW'):
            lines.append(f'         {arr}(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(f"      ASYFP=.{'TRUE' if c['asyfp'] else 'FALSE'}.")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(assign('FLC(3)', c['mach']))
        fl = c['flap']
        lines.append(assign('F(16)', fl['ndelta']))
        lines.append(assign('F(17)', fl['type']))
        lines.append(assign('F(18)', fl['stype']))
        for k in range(4):
            lines.append(assign(f'F({19 + k})', fl['deltal'][k]))
            lines.append(assign(f'F({29 + k})', fl['deltar'][k]))
            lines.append(assign(f'WING({201 + k})', c['delcl6'][k]))
            lines.append(assign(f'PW({104 + 71 + k})', c['cfact'][k]))
            lines.append(assign(f'HT({211 + k})', c['clrlm6'][k]))
            lines.append(assign(f'TCD({43 + k})', c['bd'][k]))
        for k, v in enumerate(c['cnym6']):
            lines.append(assign(f'BODY({201 + k})', v))
        for name, key in [('WING(101)', 'wing_cla'), ('HT(101)', 'tail_cla'),
                          ('SB(205)', 'tra70'), ('SB(313)', 'trah70'),
                          ('HTIN(4)', 'sspn'), ('HTIN(3)', 'sspne'),
                          ('A(7)', 'astrw'), ('DWASH(41)', 'depsda'),
                          ('BLREF', 'blref')]:
            lines.append(assign(name, c[key]))
        lines.append('      CALL M40O50')
        for tag, text in [('WING', '(WING(K),K=201,260)'),
                          ('HT', '(HT(K),K=211,220)'),
                          ('BODY', '(BODY(K),K=201,220)'),
                          ('TRN', '(PW(K),K=294,300)')]:
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{tag}',{text}")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('trnyrl', driver(all_cases), ROUTINES))
    print(save('trnyrl', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
