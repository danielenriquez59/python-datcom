"""
Probe LIFTCF and save the fixture.

Run from the repository root: ``python test_parity/probes/liftcf.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['liftcf', 'clmxbs', 'angles', 'tbfunx', 'quad', 'tlinex',
            'tlin1x', 'tlin3x', 'glook', 'switch']

GEOMETRY = [('area', 3), ('aspect_ratio_inboard', 5), ('aspect_ratio', 7),
            ('inboard_span', 23), ('planform_length', 29), ('cos_le', 37),
            ('tan_le', 62), ('arclss_factor', 123), ('arclss_ratio', 125),
            ('a159', 159), ('a160', 160)]

ALPHA = [-8., -4., -2., 0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.,
         24., 28., 35.]


def case(rng, kind, aspect_ratio, mach=0.4, alpha=ALPHA, **overrides):
    sweep = rng.uniform(5.0, 60.0)
    c = {
        'planform_type': float(kind),
        'alpha_deg': list(alpha),
        'geometry': {
            'area': rng.uniform(80., 300.), 'aspect_ratio': aspect_ratio,
            'aspect_ratio_inboard': aspect_ratio * rng.uniform(0.3, 0.8),
            'inboard_span': rng.uniform(3.0, 8.0),
            'planform_length': rng.uniform(10.0, 30.0),
            'cos_le': np.cos(np.radians(sweep)),
            'tan_le': np.tan(np.radians(sweep)),
            'arclss_factor': rng.uniform(0.0, 0.5),
            'arclss_ratio': 3.0,
            'a159': rng.uniform(0.2, 1.08),
            'a160': rng.uniform(0.5, 8.0),
        },
        'section': {'deltay': rng.uniform(1.0, 4.0),
                    'xovc': float(rng.choice([0.3, 0.4])),
                    'sspne': rng.uniform(8.0, 20.0)},
        'lift': {'cla': rng.uniform(0.03, 0.08),
                 'clmax': rng.uniform(0.7, 1.4),
                 'alpha_clmax': rng.uniform(12.0, 26.0)},
        'flight': {'mach': mach, 'beta': np.sqrt(abs(1.0 - mach**2)),
                   'alpha_zero_lift': rng.uniform(-3.0, 0.0)},
        'sref': rng.uniform(100., 400.),
        'angle_state': [0.0] * 12,
    }
    for key, value in overrides.items():
        for block in ('geometry', 'section', 'lift', 'flight'):
            if key in c[block]:
                c[block][key] = value
    return c


def cases():
    rng = np.random.default_rng(4133)
    out = []
    for _ in range(6):
        out.append(case(rng, 1, rng.uniform(3.5, 10.0)))
        out.append(case(rng, 1, rng.uniform(1.1, 2.9)))
    # Aspect ratio at or below one takes the unscaled CNAA90.
    out.append(case(rng, 1, 0.8))
    out.append(case(rng, 1, 1.0))
    # An angle exactly at the zero-lift angle, and past 90 degrees less
    # the stall angle.
    c = case(rng, 1, 6.0, alpha_zero_lift=-2.0)
    c['alpha_deg'] = [-2.0, 0.0, 30.0, 60.0, 85.0]
    out.append(c)
    # Angle state left by a previous call: the stall record matching the
    # new stall angle to within EPS, with stale trig values in it.
    c = case(rng, 1, 6.0, alpha_zero_lift=-1.0, alpha_clmax=15.0)
    stall = np.radians(16.0)
    c['angle_state'] = [16.0, stall, 0.3, 0.9, 0.31, stall + 1e-6,
                        4.0, np.radians(4.0), 0.07, 0.99, 0.07,
                        np.radians(4.0)]
    out.append(c)
    # Double delta: the dashed region (A < 3, M >= 0.7, beta*tan <= 7) and
    # the solid curve, with angles on both sides of 12 degrees.
    for mach, ar in ((0.8, 2.0), (0.9, 1.4), (0.4, 2.0), (0.8, 4.0)):
        out.append(case(rng, 2, ar, mach=mach, tan_le=1.2))
    # Cranked: the same regions, and a sweep putting the break angle off
    # the Figure 4.1.3.3-57 grid.
    for mach, ar, tan_le in ((0.8, 2.0, 1.2), (0.4, 5.0, 0.8),
                             (0.4, 5.0, 0.5), (0.85, 2.5, 3.0)):
        out.append(case(rng, 3, ar, mach=mach, tan_le=tan_le))
    # Curved.
    out.append(case(rng, 4, 5.0))
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG',
             '      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF',
             '      COMMON /CONSNT/ PI, DEG, UNUSED, RAD',
             '      DIMENSION A(195),B(49),AOUT(380),AIN(100),WT(4)',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV /",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,380')
        lines.append('         AOUT(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append('         IF(K.LE.100) AIN(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      AIN(15)=WT({int(c['planform_type'])})")
        lines.append(f"      NALPHA={len(c['alpha_deg'])}")
        for j, a in enumerate(c['alpha_deg']):
            lines.append(assign(f'B({j + 23})', a))
        for key, index in GEOMETRY:
            lines.append(assign(f'A({index})', c['geometry'][key]))
        for k, v in enumerate(c['angle_state']):
            lines.append(assign(f'A({147 + k})', v))
        lines.append(assign('AIN(3)', c['section']['sspne']))
        lines.append(assign('AIN(17)', c['section']['deltay']))
        lines.append(assign('AIN(18)', c['section']['xovc']))
        lines.append(assign('AOUT(101)', c['lift']['cla']))
        lines.append(assign('B(44)', c['lift']['clmax']))
        lines.append(assign('B(43)', c['lift']['alpha_clmax']))
        lines.append(assign('B(1)', c['flight']['mach']))
        lines.append(assign('B(2)', c['flight']['beta']))
        lines.append(assign('B(49)', c['flight']['alpha_zero_lift']))
        lines.append(assign('SREF', c['sref']))
        lines.append('      CALL LIFTCF(A,B,AIN,AOUT)')
        for tag, array, start in [('CL', 'AOUT', 20), ('CN', 'AOUT', 60),
                                  ('ALPHA', 'B', 22)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"({array}({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'B45',B(45)")
        lines.append("      WRITE(6,'(A,12ES25.16)') 'STATE',"
                     "(A(146+J),J=1,12)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('liftcf', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('liftcf', payload))


if __name__ == '__main__':
    main()
