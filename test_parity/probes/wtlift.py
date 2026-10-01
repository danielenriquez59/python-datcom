"""
Probe WTLIFT (and CLMXBS through it) and save the fixture.

Run from the repository root: ``python test_parity/probes/wtlift.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['wtlift', 'clmxbs', 'tbfunx', 'quad', 'tlinex', 'tlin1x',
            'tlin3x', 'glook', 'switch']

# Index of each planform in the source's WTYPE, and the Python constant.
TYPES = {1: 'STRA', 2: 'DOUB', 3: 'CRAN', 4: 'CURV'}

# (key, A-block index) for the geometry the driver loads.
GEOMETRY = [('area', 3), ('aspect_ratio', 7), ('taper_ratio', 27),
            ('sweep_le_deg', 34), ('cos_le', 37), ('tan_le', 38),
            ('tan_c2', 50), ('arclss_classified', 124),
            ('arclss_ratio', 125), ('aspect_ratio_inboard', 5),
            ('tan_c2_inboard', 74), ('aspect_ratio_outboard', 168),
            ('tan_c2_outboard', 98)]

OUTPUTS = [('CLA', 'AOUT(101)'), ('B43', 'B(43)'), ('B44', 'B(44)'),
           ('A144', 'A(144)'), ('A145', 'A(145)'), ('A146', 'A(146)'),
           ('A159', 'A(159)'), ('A160', 'A(160)'), ('A171', 'A(171)'),
           ('A172', 'A(172)')]


def case(rng, kind, aspect_ratio, arclss_ratio=3.0, **overrides):
    sweep = rng.uniform(0.0, 55.0)
    mach = rng.uniform(0.1, 0.7)
    c = {
        'planform_type': float(kind),
        'geometry': {
            'area': rng.uniform(80., 300.), 'aspect_ratio': aspect_ratio,
            'taper_ratio': rng.uniform(0.0, 1.0), 'sweep_le_deg': sweep,
            'cos_le': np.cos(np.radians(sweep)),
            'tan_le': np.tan(np.radians(sweep)),
            'tan_c2': np.tan(np.radians(sweep)) * rng.uniform(0.5, 0.9),
            'arclss_classified': rng.uniform(0.6, 1.5),
            'arclss_ratio': arclss_ratio,
            'aspect_ratio_inboard': aspect_ratio * rng.uniform(0.3, 0.8),
            'tan_c2_inboard': rng.uniform(0.0, 1.5),
            'aspect_ratio_outboard': aspect_ratio * rng.uniform(1.0, 2.0),
            'tan_c2_outboard': rng.uniform(0.0, 0.5),
        },
        'section': {'deltay': rng.uniform(1.0, 4.8),
                    'xovc': rng.choice([0.25, 0.3, 0.35, 0.4, 0.45]),
                    'cla': rng.uniform(0.09, 0.115),
                    'clmax': rng.uniform(1.0, 1.7)},
        'flight': {'mach': mach, 'beta': np.sqrt(1.0 - mach**2),
                   'alpha_zero_lift': rng.uniform(-4.0, 0.0)},
        'sref': rng.uniform(100., 400.),
    }
    for key, value in overrides.items():
        for block in ('geometry', 'section', 'flight'):
            if key in c[block]:
                c[block][key] = value
    return c


def cases():
    rng = np.random.default_rng(3413)
    out = []
    for kind in (1, 2, 3):
        for _ in range(8):
            out.append(case(rng, kind, rng.uniform(3.5, 10.0)))
        for _ in range(8):
            out.append(case(rng, kind, rng.uniform(0.5, 2.9)))
    # Low aspect ratio with A(160) on both sides of 4.5.
    out.append(case(rng, 1, 2.0, sweep_le_deg=10.0, tan_le=0.176))
    out.append(case(rng, 1, 2.8, sweep_le_deg=60.0, tan_le=1.732))
    # Both Figure 4.1.3.4-23 parts at a fixed configuration.
    out.append(case(rng, 1, 1.5, xovc=0.35))
    out.append(case(rng, 1, 1.5, xovc=0.36))
    # Extrapolation: sweep past 60 degrees, DELTAY off both grid ends,
    # Mach below the Figure 4.1.3.4-22 grid.
    out.append(case(rng, 1, 6.0, sweep_le_deg=68.0, tan_le=2.475))
    out.append(case(rng, 1, 6.0, deltay=0.8))
    out.append(case(rng, 1, 6.0, deltay=5.2))
    out.append(case(rng, 1, 6.0, mach=0.05, beta=np.sqrt(1 - 0.05**2)))
    out.append(case(rng, 1, 1.2, mach=0.05, beta=np.sqrt(1 - 0.05**2)))
    # A cranked wing beyond the Figure 4.1.3.2-52 grid, both ends.
    out.append(case(rng, 3, 12.0))
    out.append(case(rng, 3, 1.1, arclss_ratio=1.0))
    # Curved: only the panel slopes are set.
    out.append(case(rng, 4, 5.0))
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF',
             '      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA,IG',
             '      COMMON /CONSNT/ PI, DEG, UNUSED, RAD',
             '      DIMENSION A(195),B(49),AOUT(101),AIN(101),WT(4)',
             "      DATA WT /4HSTRA ,4HDOUB ,4HCRAN ,4HCURV /",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        # Clear every array, so no output can be a previous case's value.
        label = 1000 + n
        lines.append(f'      DO {label} K=1,195')
        lines.append('         A(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append('         IF(K.LE.101) AOUT(K)=0.')
        lines.append('         IF(K.LE.101) AIN(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      AIN(15)=WT({int(c['planform_type'])})")
        for key, index in GEOMETRY:
            lines.append(assign(f'A({index})', c['geometry'][key]))
        lines.append(assign('AIN(17)', c['section']['deltay']))
        lines.append(assign('AIN(18)', c['section']['xovc']))
        lines.append(assign('A(131)', c['section']['cla']))
        lines.append(assign('A(132)', c['section']['clmax']))
        lines.append(assign('B(1)', c['flight']['mach']))
        lines.append(assign('B(2)', c['flight']['beta']))
        lines.append(assign('B(49)', c['flight']['alpha_zero_lift']))
        lines.append(assign('SREF', c['sref']))
        lines.append('      CALL WTLIFT(A,B,AIN,AOUT)')
        for tag, source in OUTPUTS:
            lines.append(f"      WRITE(6,'(A,ES25.16)') '{tag}',{source}")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('wtlift', driver(all_cases), ROUTINES))
    payload = [{'inputs': c,
                'outputs': {k: (v if k == '_text' else v[0])
                            for k, v in r.items()}}
               for c, r in zip(all_cases, records)]
    print(save('wtlift', payload))


if __name__ == '__main__':
    main()
