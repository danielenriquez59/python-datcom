"""
Probe M19O23, the supersonic body overlay, with SYPBOD and EXSUBT stubbed,
and save the fixture.

Run from the repository root: ``python test_parity/probes/m19o23.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402

LEGACY = Path(__file__).resolve().parent.parent.parent / 'datcom-legacy' / \
    'datcom_2000'
U = 1.0e-30
CASES = [
    {'mach': 2.0, 'tsmach': 1.4, 'alpha': [-2.0, 0.0, 2.0, 4.0, 8.0],
     'option': [U, U, U, U]},
    {'mach': 1.6, 'tsmach': 1.2, 'alpha': [0.0, 5.0, 10.0, 20.0],
     'option': [250.0, 8.0, 2.0e-4, 30.0]},
    {'mach': 1.1, 'tsmach': 1.2, 'alpha': [0.0, 4.0, 8.0],
     'option': [U, 9.0, U, 40.0]},
]
STUBS = probe.STUBS + """\
      SUBROUTINE SYPBOD(I)
      RETURN
      END
      SUBROUTINE EXSUBT
      RETURN
      END
"""


def _declarations():
    text = (LEGACY / 'm19o23.f').read_text()
    body = text.split('\n', 1)[1].split('      NOVLY=19')[0]
    return '\n'.join(ln[:72] for ln in body.splitlines()
                     if ln[:1] not in 'cC' and 'DATA ROUTID' not in ln)


def driver():
    lines = ['      PROGRAM PROBE', _declarations(),
             '      UNUSED=1.D-30', '      RAD=57.2957795D0', '      I=1']
    for n, c in enumerate(CASES):
        lines.append("      WRITE(6,'(A,I4)') 'CASE',{}".format(n))
        lines += [f'      DO {100 + n} K=1,400',
                  f'{100 + n:5d} BODY(K)=(MOD(K,37)-11)*0.0625D0']
        lines.append(f"      NALPHA={len(c['alpha'])}")
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('TSMACH', c['tsmach']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + k})', v))
        for name, v in zip(('SREF', 'CBARR', 'ROUGFC', 'BLREF'),
                           c['option']):
            lines.append(assign(name, v))
        lines += [assign('A(4)', 310.0), assign('A(122)', 7.5),
                  assign('WINGIN(4)', 16.0)]
        lines.append('      CALL M19O23')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'OPT',SREF,CBARR,ROUGFC,"
                     'BLREF')
        lines.append("      WRITE(6,'(A,200ES25.16)') 'BODY',(BODY(K),"
                     'K=1,200)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    records = parse_records(probe.run(
        'm19o23', driver(), ['m19o23', 'tbfunx', 'quad'], stubs=STUBS))
    print(save('m19o23', [{'inputs': c, 'outputs': r}
                          for c, r in zip(CASES, records)]))


if __name__ == '__main__':
    main()
