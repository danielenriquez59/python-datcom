"""
Probe CONV, the input unit conversion and scaling, and save the fixture.

Run from the repository root: ``python test_parity/probes/conv.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['conv']

COMMONS = """\
      COMMON /FLGTCD/ FLC(42), RN(20), NGH, GRDH(10), PINF(20), FLC1(3),
     1                ALT(20), TINF(20), VINF(20), WT
      COMMON /OPTION/ SREF, CBARR, ROUGFC, BLREF
      COMMON /SYNTSS/ SYNA(19)
      COMMON /BODYI/  XNX, X(20), S(20), P(20), R(20), ZU(20),
     1                ZL(20), BTY(2), BL(3)
      COMMON /WINGI/  WGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /VTI/    VTIN(154), TVTIN(8), VFIN(154)
      COMMON /POWER/  PWIN(29), LBIN(21)
      COMMON /FLAPIN/ F(138)
      COMMON /BDATA/  BD(762)
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      REAL LBIN
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
"""

UNUSED = 1.0e-30
# Array: (length, fill step).
ARRAYS = {'rn': (20, 1.1e5), 'alt': (20, 500.0), 'pinf': (20, 25.0),
          'tinf': (20, 3.0), 'vinf': (20, 40.0), 'grdh': (10, 2.5),
          'syna': (19, 1.25), 'x': (20, 1.5), 's': (20, 0.75),
          'p': (20, 0.5), 'r': (20, 0.2), 'zu': (20, 0.1),
          'zl': (20, -0.1), 'bl': (3, 4.0), 'wgin': (101, 0.3),
          'htin': (154, 0.2), 'vtin': (154, 0.15), 'tvtin': (8, 0.9),
          'vfin': (154, 0.35), 'pwin': (29, 0.45), 'lbin': (21, 0.55),
          'f': (138, 0.05)}
SCALARS = {'wt': 12000.0, 'rougfc': 0.00016, 'sref': 300.0, 'cbarr': 8.0,
           'blref': 36.0}
FLAGS = ('symfp', 'asyfp', 'trajet', 'hypef')


def case(idim=1, scale=1.0, x0=2.0, xnx=12.0, unused=(), **flags):
    c = {n: [step * k for k in range(1, length + 1)]
         for n, (length, step) in ARRAYS.items()}
    c['x'] = [x0 + v for v in c['x']]
    c['x'][0] = x0
    for name, k in unused:
        c[name][k - 1] = UNUSED
    c.update(SCALARS, xnx=xnx, idim=idim, scale=scale, x0=x0,
             bd={33: 11.0, 65: 14.0, 74: 3.0, 82: 5.0},
             unused=[list(u) for u in unused],
             **{k: bool(flags.get(k, False)) for k in FLAGS})
    return c


def cases():
    spots = [('rn', 3), ('syna', 11), ('syna', 2), ('htin', 100),
             ('vtin', 120), ('f', 40), ('pwin', 19), ('x', 5)]
    return [
        case(),                                          # shift only
        case(x0=0.0),                                    # nothing at all
        case(idim=2, symfp=True, unused=spots),
        case(idim=3, scale=0.1, asyfp=True, trajet=True),
        case(idim=4, scale=2.0, hypef=True, symfp=True),
        case(idim=1, scale=0.5, x0=-3.0, xnx=20.0, trajet=True),
        case(idim=3, scale=1.0, x0=0.0, asyfp=True, hypef=True,
             unused=spots),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      UNUSED=1.E-30']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        for i, (name, (length, step)) in enumerate(ARRAYS.items()):
            label = 1000 + 30 * n + i
            lines.append(f'      DO {label} K=1,{length}')
            lines.append(f'{label:5d} {name.upper()}(K)=K*{step!r}')
        x0 = c['x0']
        lines.append(f'      DO {900 + n} K=1,20')
        lines.append(f'      X(K)=K*1.5+({x0!r})')
        lines.append(f'{900 + n:5d} CONTINUE')
        lines.append(assign('X(1)', x0) if x0 else '      X(1)=0.')
        for name, k in c['unused']:
            lines.append(assign(f'{name.upper()}({k})', UNUSED))
        for k, v in SCALARS.items():
            lines.append(assign(k.upper(), v))
        lines.append(assign('XNX', c['xnx']))
        for k, v in c['bd'].items():
            lines.append(assign(f'BD({k})', v))
        for flag in FLAGS:
            lines.append(f"      {flag.upper()}="
                         f".{'TRUE' if c[flag] else 'FALSE'}.")
        lines.append(f"      CALL CONV({c['idim']},{c['scale']!r})")
        for name, (length, _) in ARRAYS.items():
            lines.append(f"      WRITE(6,'(A,{length}ES25.16)') "
                         f"'{name.upper()}',{name.upper()}")
        lines.append("      WRITE(6,'(A,5ES25.16)') 'SC',WT,ROUGFC,SREF,"
                     "CBARR,BLREF")
        lines.append("      WRITE(6,'(A,5ES25.16)') 'BD',BD(11),BD(33),"
                     "BD(65),BD(74),BD(82)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('conv', driver(all_cases), ROUTINES))
    print(save('conv', [{'inputs': c, 'outputs': r}
                        for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
