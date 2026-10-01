"""
Probe SUPCMQ (wing) and SUPHMQ (horizontal tail), the supersonic and
transonic Cm-q, and save the fixture.

Run from the repository root: ``python test_parity/probes/supcmq.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402

ROUTINES = ['interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'quad', 'tlip3x', 'tlip2x', 'tlip1x', 'yup', 'tbfunx']
FLOLOG = """\
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2        SUPERS,SUBSON,TRANSN,HYPERS,
     3        SYMFP,ASYFP,TRIMC,TRIM
"""
COMMON = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
"""
WING = COMMON + """\
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  WINGIN(77)
      COMMON /POWR/   DYN(213)
      COMMON /IWING/  PW, WING(400)
      COMMON /SUPWH/  SLG(141)
      COMMON /SBETA/  STB(135),TRA(108)
""" + FLOLOG
TAIL = COMMON + """\
      COMMON /HTDATA/ A(195)
      COMMON /HTI/    WINGIN(77)
      COMMON /BDATA/  XBD(300),DYN(213)
      COMMON /IHT/    PW, WING(380)
      COMMON /SUPWH/  XSLG(141),SLG(141)
      COMMON /SBETA/  XSTB(243),TRA(108),STB(135)
""" + FLOLOG


def case(mach, tanle, taper, transn=False, straight=True, ar=3.0):
    return {'mach': mach, 'tanle': tanle, 'taper': taper, 'transn': transn,
            'straight': straight, 'ar': ar,
            'beta': math.sqrt(mach ** 2 - 1.) if mach > 1 else 0.0}


CASES = [case(1.8, 1.0, 0.5), case(2.5, 0.3, 0.2), case(1.4, 2.0, 0.0),
         case(1.4, 2.0, 0.1), case(1.4, 2.0, 0.5), case(1.4, 2.0, 0.9),
         case(1.6, 2.0, 0.5, straight=False), case(1.5, 0.0, 0.5),
         case(0.9, 2.0, 0.5, transn=True), case(0.8, 2.0, 0.5, transn=True),
         case(0.9, 1.0, 0.1, transn=True), case(0.95, 0.0, 0.0, transn=True)]
TRA = {6: 0.85, 12: 0.06, 16: 0.8, 17: 0.9, 18: 1.0, 19: 1.1, 20: 1.2,
       21: 0.07, 22: 0.075, 23: 0.08, 24: 0.078, 25: 0.072, 105: 0.45}


def driver(commons, call):
    lines = ['      PROGRAM PROBE', commons.rstrip('\n'),
             '      DATA WS/4HSTRA/,WD/4HDOUB/',
             '      PI=3.141592654D0', '      DR=0.01745329D0',
             '      UNUSED=1.D-30', '      RAD=57.2957795D0', '      IM=1',
             '      SR=320.0D0', '      CBARR=8.5D0']
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,213', '      DYN(K)=K*0.01D0',
                  '      IF(K.LE.195) A(K)=K*0.03D0',
                  '      IF(K.LE.141) SLG(K)=K*0.02D0',
                  '      IF(K.LE.108) TRA(K)=K*0.005D0',
                  f'{100 + n:5d} CONTINUE', '      WING(101)=0.065D0',
                  '      WING(201)=1.7D0', '      WING(221)=-9.0D0']
        lines.append(f"      TRANSN=.{'TRUE' if c['transn'] else 'FALSE'}.")
        lines.append(f"      WINGIN(15)={'WS' if c['straight'] else 'WD'}")
        lines.append(assign('WINGIN(69)', 0.1))
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in {3: 280.0, 7: c['ar'], 10: 11.0, 16: 7.8,
                     27: c['taper'], 29: 12.0, 50: 0.3, 62: c['tanle'],
                     67: 0.9, 68: 0.5, 173: 1.25}.items():
            lines.append(assign(f'A({k})', v))
        lines += [assign('SLG(7)', 3.1), assign('SLG(134)', 0.42),
                  assign('SLG(135)', -0.12)]
        for k, v in TRA.items():
            lines.append(assign(f'TRA({k})', v))
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,213ES25.16)') 'DYN',DYN")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,ES25.16)') 'CMQ',WING(221)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = {}
    for name, commons in (('supcmq', WING), ('suphmq', TAIL)):
        out[name] = parse_records(probe.run(name, driver(commons,
                                                         name.upper()),
                                            ROUTINES + [name]))
    print(save('supcmq', {'cases': CASES, 'tra': TRA, 'records': out}))


if __name__ == '__main__':
    main()
