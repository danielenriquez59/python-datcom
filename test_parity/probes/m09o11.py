"""
Probe overlay M09O11 (DWASH, then DYPRLS) and save the fixture.

Reuses the DWASH probe's cases, adding what DYPRLS reads: the wing drag
``B(46)``, the exposed MAC ``A(16)``, the tail angle ``A(11)``, and the
``KEPSLN`` switch that takes the wake deflection from DWASH.

Run from the repository root: ``python test_parity/probes/m09o11.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402
import dwash  # noqa: E402

ROUTINES = ['m09o11', 'dyprls'] + dwash.ROUTINES

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

EXPER = """\
      COMMON /EXPER/ KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1               KBODY, KWING, KHT, KVT, KWB, KDWASH(3),
     2               ALPOW, ALPLW, ALPOH, ALPLH
      LOGICAL KDWASH
"""


def cases():
    out = []
    for n, c in enumerate(dwash.cases()):
        for kepsln in (False, True):
            c = dict(c, dyprls={'cdow': 0.006 + 0.001 * n,
                                'mac': 7.5, 'gamma': 0.02 + 0.01 * n,
                                'kepsln': kepsln})
            out.append(c)
    return out


def driver(all_cases) -> str:
    base = dwash.driver(all_cases)
    lines = []
    case = -1
    for line in base.splitlines():
        if line.startswith("      WRITE(6,'(A,I3)') 'CASE'"):
            case += 1
        if line == '      CALL DWASH':
            d = all_cases[case]['dyprls']
            lines.append(assign('B(46)', d['cdow']))
            lines.append(assign('A(16)', d['mac']))
            lines.append(assign('A(11)', d['gamma']))
            lines.append(f"      KDWASH(2)=.{'TRUE' if d['kepsln'] else 'FALSE'}.")
            lines.append('      KDWASH(1)=.FALSE.')
            lines.append('      CALL M09O11')
            continue
        lines.append(line)
        if line.startswith('      LOGICAL VERTUP'):
            lines.extend(EXPER.rstrip('\n').splitlines())
    lines.insert(-2, '')
    text = '\n'.join(lines)
    return text.replace(
        "      WRITE(6,'(A,ES25.16)') 'A20',A(20)",
        "      WRITE(6,'(A,ES25.16)') 'A20',A(20)\n"
        "      WRITE(6,'(A,30ES25.16)') 'QOQI',(DWASHI(J),J=1,NALPHA)")


def main():
    all_cases = cases()
    records = parse_records(run('m09o11', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('m09o11', payload))


if __name__ == '__main__':
    main()
