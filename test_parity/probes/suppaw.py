"""
Probe SUPPAW (wing) and SUPPAH (horizontal tail), the supersonic CL-q,
and save the fixture.

Run from the repository root: ``python test_parity/probes/suppaw.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402
from supcmq import CASES, COMMON, FLOLOG  # noqa: E402

ROUTINES = ['interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'quad', 'tlip3x', 'tlip2x', 'tlip1x', 'yup']
WING = COMMON + """\
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  WINGIN(77)
      COMMON /IWING/  PW, WING(400)
      COMMON /POWR/   DYN(213)
      COMMON /SUPWH/  SLG(141)
""" + FLOLOG
TAIL = COMMON + """\
      COMMON /HTDATA/ A(195)
      COMMON /HTI/    WINGIN(154)
      COMMON /IHT/    PW, WING(380)
      COMMON /BDATA/  XBD(300),DYN(213)
      COMMON /SUPWH/  XSLG(141),SLG(141)
""" + FLOLOG


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
                  f'{100 + n:5d} CONTINUE', '      WING(201)=1.7D0']
        lines.append(f"      TRANSN=.{'TRUE' if c['transn'] else 'FALSE'}.")
        lines.append(f"      WINGIN(15)={'WS' if c['straight'] else 'WD'}")
        for k, v in {3: 280.0, 7: c['ar'], 10: 11.0, 16: 7.8,
                     27: c['taper'], 29: 12.0, 62: c['tanle']}.items():
            lines.append(assign(f'A({k})', v))
        lines += [assign('SLG(1)', c['beta'] or 0.7),
                  assign('SLG(7)', 3.1), assign('SLG(134)', 0.42),
                  assign('SLG(135)', -0.12)]
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,213ES25.16)') 'DYN',DYN")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,141ES25.16)') 'SLG',SLG")
        lines.append("      WRITE(6,'(A,ES25.16)') 'CLQ',WING(201)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = {}
    for name, commons in (('suppaw', WING), ('suppah', TAIL)):
        out[name] = parse_records(probe.run(name, driver(commons,
                                                         name.upper()),
                                            ROUTINES + [name]))
    print(save('suppaw', {'cases': CASES, 'records': out}))


if __name__ == '__main__':
    main()
