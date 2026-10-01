"""
Probe M46O56, the dynamic-derivative overlay, with its five routines
stubbed to announce themselves, and save the fixture.

Run from the repository root: ``python test_parity/probes/m46o56.py``.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

LEGACY = Path(__file__).resolve().parent.parent.parent / 'datcom-legacy' / \
    'datcom_2000'
FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'm46o56.json')
CALLEES = ['DYNBOD', 'DNPAWB', 'DNPWBT', 'SUBWBT', 'CLRDER']
# Block: (COMMON name, size, fill step).
BLOCKS = {'body': ('BODY', 400, 0.001), 'wing': ('WING', 400, 0.002),
          'ht': ('HT', 380, 0.003), 'vt': ('VT', 380, 0.004),
          'vf': ('VF', 380, 0.005), 'bw': ('BW', 380, 0.006),
          'bh': ('BH', 380, 0.007), 'bv': ('BV', 380, 0.008),
          'bwh': ('BWH', 380, 0.009), 'bwv': ('BWV', 380, 0.011),
          'bwhv': ('BWHV', 380, 0.012), 'power': ('POWER', 200, 0.013),
          'dwsh': ('DWSH', 60, 0.014), 'dyn': ('DYN', 213, 0.015),
          'dynh': ('DYNH', 213, 0.016)}
CASES = [
    {'flags': {'bo': True, 'wgpl': True, 'htpl': True, 'vtpl': True,
               'subson': True}, 'nalpha': 6, 'mach': 0.6},
    {'flags': {'bo': True, 'wgpl': True, 'vfpl': True, 'transn': True,
               'dmpcse': True}, 'nalpha': 4, 'mach': 0.95},
    {'flags': {'bo': True, 'wgpl': True, 'htpl': True, 'hypers': True,
               'dpivf': True, 'dpiht': True}, 'nalpha': 3, 'mach': 6.0},
    {'flags': {'bo': True, 'wgpl': True, 'hypers': True, 'dpivf': True,
               'dpdynh': True}, 'nalpha': 2, 'mach': 1.2},
    {'flags': {'wgpl': True, 'htpl': True, 'vtpl': True, 'subson': True,
               'dpibdy': True, 'dpidwh': True, 'dpipwr': True},
     'nalpha': 20, 'mach': 0.3},
]


def _declarations():
    """M46O56's own declarations, with its COMMON names."""
    text = (LEGACY / 'm46o56.f').read_text()
    body = text.split('\n', 1)[1]
    decl = body.split('      NOVLY=46')[0]
    return '\n'.join(ln[:72] for ln in decl.splitlines()
                     if ln[:1] not in 'cC')


def driver() -> str:
    lines = ['      PROGRAM PROBE', _declarations(),
             '      UNUSED=1.D-30', '      I=1']
    flags = sorted({k for c in CASES for k in c['flags']} |
                   {'bo', 'wgpl', 'htpl', 'vtpl', 'vfpl', 'subson',
                    'transn', 'hypers'})
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A)') '@@CASE'")
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(f"      FLC(3)={c['mach']!r}D0")
        for flag in flags:
            v = c['flags'].get(flag, False)
            lines.append(f"      {flag.upper()}=.{'TRUE' if v else 'FALSE'}.")
        for k, (name, size, step) in enumerate(BLOCKS.values()):
            label = 1000 + 20 * n + k
            lines.append(f'      DO {label} K=1,{size}')
            lines.append(f'{label:5d} {name}(K)=K*{step!r}D0')
        lines.append('      CALL M46O56')
        lines.append("      WRITE(6,'(A)') '@@BLOCKS'")
        for name, size, _ in BLOCKS.values():
            lines.append(f"      WRITE(6,'(A,{size}ES25.16)') '{name}',"
                         f"({name}(K),K=1,{size})")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def stubs() -> str:
    out = [probe.STUBS]
    for name in CALLEES:
        out += [f'      SUBROUTINE {name}',
                f"      WRITE(6,'(A)') 'CALL {name}'",
                '      RETURN', '      END']
    return '\n'.join(out) + '\n'


def main():
    out = probe.run('m46o56', driver(), ['m46o56', 'dmpary'],
                    stubs=stubs())
    cases = []
    for chunk in out.split('@@CASE\n')[1:]:
        text, blocks = chunk.split('@@BLOCKS\n')
        values = {}
        for line in blocks.strip('\n').split('\n'):
            fields = line.split()
            values[fields[0]] = [float(v) for v in fields[1:]]
        cases.append({'text': text.rstrip('\n').split('\n') if text.strip()
                      else [], 'blocks': values})
    FIXTURE.write_text(json.dumps({'cases': CASES, 'blocks': BLOCKS,
                                   'records': cases}))
    print(FIXTURE)


if __name__ == '__main__':
    main()
