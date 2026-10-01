"""
Probe SUPHYW, the horizontal tail's supersonic CL-p, CY-p and CN-p, and
save the fixture.

Run from the repository root: ``python test_parity/probes/suphyw.py``.
The cases run in one program (INTEP3 keeps its chart pair between calls).
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402

ROUTINES = ['suphyw', 'intep3', 'tlip2x', 'tlip1x', 'yup', 'switch',
            'glook', 'quad', 'interx', 'tlin1x', 'tlinex', 'tlin3x']
COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /HTDATA/ A(195)
      COMMON /HTI/    WINGIN(77)
      COMMON /FLGTCD/ FLC(93)
      COMMON /BDATA/  BD(300),DYN(213)
      COMMON /SUPWH/  GR(141),SLG(141)
      COMMON /IBH/    PBW,BW(380)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /IHT/    PWING, WING(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2        SUPERS,SUBSON,TRANSN,HYPERS,
     3        SYMFP,ASYFP,TRIMC,TRIM
"""


def case(mach, tanle, taper, htpl=True, bo=False, sspne=5.0, aw=3.0):
    return {'mach': mach, 'tanle': tanle, 'taper': taper, 'htpl': htpl,
            'bo': bo, 'sspne': sspne, 'aw': aw,
            'beta': math.sqrt(mach ** 2 - 1.),
            'lamle': math.degrees(math.atan(tanle))}


CASES = [case(1.5, 2.0, 0.0), case(1.5, 2.0, 0.1), case(1.8, 2.5, 0.4),
         case(2.0, 3.0, 0.6), case(1.6, 2.2, 0.9), case(2.2, 0.8, 0.0),
         case(2.5, 0.6, 0.2), case(1.8, 0.5, 0.5), case(3.0, 0.4, 0.8),
         case(2.0, 0.9, 1.0), case(1.8, 2.0, 0.5, htpl=True, bo=True,
                                   sspne=3.0),
         case(2.2, 1.5, 0.35, bo=True, sspne=5.5, aw=4.0)]
ALPHA = [-2.0, 0.0, 4.0, 10.0]


def driver():
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654D0', '      UNUSED=1.D-30',
             '      RAD=57.2957795D0', '      IM=1', '      SR=320.0D0',
             '      BLREF=36.0D0', f'      NALPHA={len(ALPHA)}']
    lines += [assign(f'FLC({23 + k})', v) for k, v in enumerate(ALPHA)]
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,380', '      WING(K)=K*0.001D0',
                  '      IF(K.LE.213) DYN(K)=K*0.002D0',
                  f'{100 + n:5d} CONTINUE']
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(f"      BO=.{'TRUE' if c['bo'] else 'FALSE'}.")
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in {4: 60.0, 58: c['lamle'], 62: c['tanle'], 74: 0.5,
                     118: c['taper'], 120: c['aw']}.items():
            lines.append(assign(f'A({k})', v))
        lines += [assign('WINGIN(3)', c['sspne']), assign('WINGIN(4)', 6.5),
                  assign('WINGIN(6)', 5.0), assign('SLG(1)', c['beta']),
                  assign('SLG(3)', 0.93), assign('SYNA(1)', 22.0)]
        lines.append('      CALL SUPHYW')
        lines.append("      WRITE(6,'(A,10ES25.16)') 'DYN',(DYN(K),"
                     'K=204,213)')
        lines.append("      WRITE(6,'(A,61ES25.16)') 'HT',(WING(K),"
                     'K=280,340)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


WING = COMMONS.replace('/HTDATA/ A(195)', '/WINGD/  A(195)').replace(
    '/HTI/    WINGIN(77)', '/WINGI/  WINGIN(77)').replace(
    '/BDATA/  BD(300),DYN(213)', '/POWR/   DYN(213)').replace(
    '/SUPWH/  GR(141),SLG(141)', '/SUPWH/  SLG(141)').replace(
    '/IBH/    PBW,BW(380)', '/IBW/    PBW,BW(380)').replace(
    '/IHT/    PWING, WING(380)', '/IWING/  PWING, WING(400)')


def main():
    records = parse_records(probe.run('suphyw', driver(), ROUTINES))
    # SUPRYW, the wing's twin: its blocks, gated on WGPL.
    wing_driver = driver().replace(COMMONS.rstrip('\n'), WING.rstrip('\n'))
    wing_driver = wing_driver.replace('HTPL=', 'WGPL=').replace(
        'CALL SUPHYW', 'CALL SUPRYW')
    wing = parse_records(probe.run('supryw', wing_driver,
                                   ['supryw'] + ROUTINES[1:]))
    print(save('suphyw', {'cases': CASES, 'alpha': ALPHA,
                          'records': records, 'wing_records': wing}))


if __name__ == '__main__':
    main()
