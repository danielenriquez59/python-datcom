"""
Probe HYPFLP (with HYPROP, FIG68 and SIMUL2), the hypersonic flap
increments, and save the fixture.

The cases run in one program, in order: HYPFLP's locals ``PHE`` and
``CPI2`` are saved between calls.

Run from the repository root: ``python test_parity/probes/hypflp.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['hypflp', 'hyprop', 'fig68', 'arcsin', 'arccos', 'simul2',
            'tbfunx', 'quad', 'tlinex', 'tlin1x', 'glook', 'switch',
            'interx', 'tlin3x']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA
      COMMON /OPTION/ SR,CBAR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLAPIN/ F(16)
      COMMON /FLGTCD/ FLC(93)
      COMMON /IBODY/  PBODY, BODY(200)
      COMMON /IWING/  PWING, WING(200)
      COMMON /IHT/    PHT, HT(200)
      COMMON /IVT/    PVT, VT(200)
      COMMON /BDATA/  HYP(80)
      LOGICAL LAMNR
      EQUIVALENCE (LAMNR,F(15))
"""

# Altitudes stay between 100,000 and 250,000 feet: below, HYPROP reads
# its defective Figure 6.3.1-43A; above, the translation corrects a
# pressure/temperature typo the compiled source keeps.
ALPHA = [0.0, 5.0, 10.0, 15.0, -5.0, 25.0]


def case(laminar=True, mach=8.0, alt=120000.0, rl=1.0e6, xhl=10.0,
         cf=2.0, deltas=(5.0, 15.0, 25.0, 35.0)):
    f = {1: alt, 2: xhl, 3: 0.3, 4: cf, 16: float(len(deltas))}
    for k, d in enumerate(deltas):
        f[5 + k] = d
    return {'alpha': ALPHA, 'mach': mach, 'rl': rl, 'f': f,
            'laminar': laminar, 'sref': 200.0, 'cbar': 8.0}


def cases():
    return [
        case(),
        case(laminar=False),
        case(mach=12.0, alt=180000.0),
        case(laminar=False, mach=6.0, rl=5.0e6),
        case(mach=5.0, rl=2.0e5, xhl=2.0, cf=4.0, alt=240000.0),
        case(laminar=False, mach=10.0, rl=1.0e7, deltas=(0.0, 10.0, 40.0)),
        case(deltas=(1.0, 3.0, 30.0)),                  # laminar, attached
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,200')
        for arr in ('BODY', 'WING', 'HT', 'VT'):
            lines.append(f'         {arr}(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['rl']))
        for k, v in c['f'].items():
            lines.append(assign(f'F({k})', v))
        lines.append(f"      LAMNR=.{'TRUE' if c['laminar'] else 'FALSE'}.")
        lines.append(assign('SR', c['sref']))
        lines.append(assign('CBAR', c['cbar']))
        lines.append('      CALL HYPFLP')
        for tag, arr in (('DCN', 'BODY'), ('DCA', 'WING'), ('DCMCN', 'HT'),
                         ('DCMCA', 'VT')):
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{tag}',")
            lines.append(f'     1({arr}(K),K=1,60)')
        lines.append("      WRITE(6,'(A,80ES25.16)') 'HYP',HYP")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('hypflp', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c['f'] = {str(k): v for k, v in c['f'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('hypflp', payload))


if __name__ == '__main__':
    main()
