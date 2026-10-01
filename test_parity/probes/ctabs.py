"""
Probe CTABS, the control-tab hinge moments and control force, and save the
fixture.

Run from the repository root: ``python test_parity/probes/ctabs.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['ctabs']

COMMONS = """\
      REAL KS,KQ,MACH
      COMMON /IBW/    PBW,    BW(380)
      COMMON /IBH/    PBH,    BH(380)
      COMMON /IBV/    PBV,    BV(380)
      COMMON /IBWH/   PBWH,   BWH(380)
      COMMON /IBWHV/  PBWHV,  BWHV(380)
      COMMON /FLGTCD/ NNNNN,MMMMMM,MACH(20),ALPHA(20),RNNUB(20),NGH
     1                ,GRDHT(10),PINF(20),STMACH,TSMACH,TR,ALT(20)
     2                ,TINF(20),VINF(20),WT,GAMMA,NALT,LOOP
      COMMON /OPTION/ SREF, CBARR, ROUGFC, BLREF
      COMMON /FLAPIN/ F(137)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,II,NALPHA,IG,NF,LF,K
"""

UNUSED = 1.0e-30
BLOCKS = {'bw': 0.001, 'bh': 0.002, 'bv': 0.003, 'bwh': 0.004,
          'bwhv': 0.005}


def case(ttype=1.0, rl=-1.0, pinf=2116.0, ndelta=3, nalpha=4):
    f = {k: 0.0 for k in range(1, 138)}
    for j in range(1, ndelta + 1):
        f[j] = -10.0 + 7.0 * j
    f[16] = float(ndelta)
    tab = [ttype, 0.25, 0.2, 3.0, 9.0, 0.3, 0.22, 4.0, 10.0,
           -0.006, -0.004, -0.003, 0.009, 0.5, -0.8, 0.3, 1.4,
           25.0, rl, 1.5, 0.2]
    for k, v in enumerate(tab):
        f[117 + k] = v
    return {'nalpha': nalpha, 'mach': 0.3, 'pinf': pinf, 'sref': 200.0,
            'cbarr': 6.0, 'alpha': [-4.0 + 3.0 * i for i in range(nalpha)],
            'f': f, **{n: {k: s * k for k in range(201, 381)}
                        for n, s in BLOCKS.items()}}


def cases():
    return [
        case(),
        case(ttype=2.0, rl=0.0),
        case(ttype=3.0, rl=2.0),
        case(ttype=2.0, rl=2.0, pinf=UNUSED),
        case(ttype=3.0, rl=0.0, ndelta=9, nalpha=20),
        case(ttype=1.0, rl=-0.5, pinf=UNUSED, ndelta=1, nalpha=1),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      UNUSED=1.E-30', '      II=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=201,380',
                  *[f'      {name.upper()}(K)=K*{step}'
                    for name, step in BLOCKS.items()],
                  f'{100 + n:5d} CONTINUE']
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign('MACH(1)', c['mach']))
        lines.append(assign('PINF(1)', c['pinf']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'ALPHA({k + 1})', v))
        for k, v in c['f'].items():
            lines.append(assign(f'F({k})', v) if v else f'      F({k})=0.')
        lines.append('      CALL CTABS')
        for name in BLOCKS:
            lines.append(f"      WRITE(6,'(A,180ES25.16)') '{name.upper()}',"
                         f"({name.upper()}(K),K=201,380)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('ctabs', driver(all_cases), ROUTINES))
    print(save('ctabs', [{'inputs': c, 'outputs': r}
                         for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
