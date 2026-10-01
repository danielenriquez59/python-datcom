"""
Probe CMALPH (with FWDXAC) and save the fixture.

Run from the repository root: ``python test_parity/probes/cmalph.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['cmalph', 'fwdxac', 'tbfunx', 'quad', 'tlin3x', 'tlinex',
            'tlin1x', 'glook', 'switch']

GEOMETRY = ['a1', 'a2', 'a3', 'a5', 'a7', 'a10', 'a16', 'a21', 'a23', 'a26',
            'a27', 'a30', 'a34', 'a37', 'a38', 'a40', 'a43', 'a62', 'a86',
            'a124', 'a125', 'a133', 'a164', 'a166', 'a167', 'a168', 'a169',
            'a171', 'a172', 'a173']

UNUSED = 1.0e-30

ALPHA = [-6., -2., 0., 4., 8., 12., 16., 20., 24., 28., 34., 40.]


def case(rng, kind=1, a7=2.0, sweep_le=35.0, sweep_out=None, mach=0.4,
         twist=0.0, cmot=-0.01, xac=UNUSED, a124=1.4, a125=3.0,
         alpha_clmax=22.0, alpha=ALPHA):
    tan_le = math.tan(math.radians(sweep_le))
    tan_out = math.tan(math.radians(sweep_out if sweep_out is not None
                                    else sweep_le * 0.6))
    beta = math.sqrt(abs(1.0 - mach**2))
    cla = 0.04 + 0.01 * a7
    alpha = np.asarray(alpha)
    cn = cla * alpha * (1.0 + 0.02 * np.abs(alpha))
    g = {
        'a1': 90.0, 'a2': 40.0, 'a3': 130.0, 'a5': a7 * 0.7, 'a7': a7,
        'a10': 12.0, 'a16': 8.0, 'a21': 3.0, 'a23': 5.0, 'a26': 0.35,
        'a27': rng.uniform(0.1, 0.6), 'a30': 6.5, 'a34': sweep_le,
        'a37': math.cos(math.radians(sweep_le)), 'a38': tan_le,
        'a40': sweep_le - 5.0, 'a43': math.cos(math.radians(sweep_le - 5.0)),
        'a62': tan_le, 'a86': tan_out, 'a124': a124, 'a125': a125,
        'a133': 7.0, 'a164': 2.5, 'a166': 4.0, 'a167': 45.0,
        'a168': a7 * 1.6, 'a169': 0.4, 'a171': 0.05, 'a172': 0.06,
        'a173': 9.0,
    }
    return {
        'planform_type': float(kind), 'alpha_deg': list(alpha),
        'geometry': g,
        'section': {'cmo': -0.02, 'cmot': cmot, 'twista': twist,
                    'deltay': rng.uniform(0.5, 2.2), 'xac': xac},
        'lift': {'cla': cla, 'cn': list(cn), 'alpha_clmax': alpha_clmax},
        'flight': {'mach': mach, 'beta': beta},
        'sref': 160.0, 'cbarr': 8.5,
    }


def cases():
    rng = np.random.default_rng(4142)
    return [
        # Straight, high aspect ratio: linear, MAC quarter chord.
        case(rng, a7=6.0, a124=1.5),
        case(rng, a7=6.0, a124=1.5, xac=0.27),
        case(rng, a7=6.0, a124=1.5, twist=-3.0, cmot=0.0),
        # Straight, low aspect ratio: Figure 26A, 26B, and the nonlinear
        # method across the stall and the reference angle.
        case(rng, a7=2.0, sweep_le=35.0, mach=0.3),
        case(rng, a7=1.6, sweep_le=60.0, mach=0.6),
        case(rng, a7=2.5, sweep_le=0.0, mach=0.3),
        case(rng, a7=1.2, sweep_le=45.0, mach=0.8, twist=2.0,
             alpha_clmax=15.0),
        # Straight, forward swept: FWDXAC, then linear since A(34) < 0.
        case(rng, a7=2.0, sweep_le=-30.0),
        # Low aspect ratio but above 6/A(124): linear.
        case(rng, a7=2.8, a124=2.5),
        # Double delta and cranked: both panels aft, inboard forward with
        # the outboard both ways, and the aft-inboard forward-outboard case.
        case(rng, kind=2, a7=1.8, sweep_le=65.0, sweep_out=45.0),
        case(rng, kind=3, a7=5.0, a124=1.5, sweep_le=40.0, sweep_out=20.0),
        case(rng, kind=3, a7=2.2, sweep_le=-20.0, sweep_out=-15.0),
        case(rng, kind=3, a7=2.2, sweep_le=-20.0, sweep_out=25.0),
        case(rng, kind=2, a7=2.2, sweep_le=50.0, sweep_out=-10.0),
        case(rng, kind=4, a7=2.0, sweep_le=55.0, sweep_out=58.0, mach=0.9),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA',
             '      COMMON /CONSNT/ PI,DR,UNUSED,RAD',
             '      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF',
             '      DIMENSION A(195),B(49),C(51),WINGIN(100),WING(400),WT(4)',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=2']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,400')
        lines.append('         WING(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.100) WINGIN(K)=0.')
        lines.append('         IF(K.LE.51) C(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha_deg'])}")
        lines.append(f"      WINGIN(15)=WT({int(c['planform_type'])})")
        for key in GEOMETRY:
            lines.append(assign(f'A({key[1:]})', c['geometry'][key]))
        s = c['section']
        lines.append(assign('WINGIN(61)', s['cmo']))
        lines.append(assign('WINGIN(67)', s['cmot']))
        lines.append(assign('WINGIN(11)', s['twista']))
        lines.append(assign('WINGIN(17)', s['deltay']))
        lines.append(assign('WINGIN(73)', s['xac']))
        for j, (a, v) in enumerate(zip(c['alpha_deg'], c['lift']['cn'])):
            lines.append(assign(f'B({23 + j})', a))
            lines.append(assign(f'WING({61 + j})', v))
        lines.append(assign('WING(101)', c['lift']['cla']))
        lines.append(assign('B(43)', c['lift']['alpha_clmax']))
        lines.append(assign('B(1)', c['flight']['mach']))
        lines.append(assign('B(2)', c['flight']['beta']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append('      CALL CMALPH(A,B,C,WINGIN,WING)')
        lines.append("      WRITE(6,'(A,30ES25.16)') 'CM',"
                     "(WING(40+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,51ES25.16)') 'C',(C(J),J=1,51)")
        lines.append("      WRITE(6,'(A,5ES25.16)') 'OUT',WING(121),B(47),")
        lines.append("     1  A(170),A(38),A(62)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('cmalph', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('cmalph', payload))


if __name__ == '__main__':
    main()
