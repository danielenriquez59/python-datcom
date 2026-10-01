"""
Probe DWASH: run the legacy routine on chosen cases and save the fixture.

Run from the repository root: ``python test_parity/probes/dwash.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['dwash', 'tbfunx', 'quad', 'tlinex', 'tlin1x', 'glook',
            'switch', 'trapz']

# DWASH's own COMMON declarations, copied so the layout is identical.
COMMONS = """\
      COMMON /FLGTCD/ FLC(95)
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195), B(49)
      COMMON /WHWB/   FACT(182), WB(39)
      COMMON /WINGI/  CHRDTP,SSPNOP,SSPNE,SSPN,CHRDBP,CHRDR,SAVSI,SAVSO,
     1                CHSTAT,SWAFP,TWISTA,SSPNDD,DHDADI,DHDADO,
     2                TYPE,
     3                TOVC,DELTAY,XOVC,CLI,ALPHAI,CLALPA(20),
     4                CLMAX(20),CMO,LERI,LERO,CAMBER,TOVCO,XOVCO,CMOT,
     5                CLMAXL,CLAMO,TCEFF,KSHARP,XAC(20),ARCL,RSVD(2),
     6                SLOPE(6),TWASH
      COMMON /HTI/    HTIN(154)
      COMMON /IDWASH/ PDWASH, DWASHI(60)
      COMMON /IWING/  PWING, WING(400)
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      LOGICAL VERTUP
"""


def wing_block(chrdr, chrdtp, sspn, sspne, tan_le):
    """The straight-wing A-block values DWASH reads, as WTGEOM forms them."""
    exposed_root = chrdr * (chrdtp / chrdr + (1 - chrdtp / chrdr) * sspne / sspn)
    area = sspne * (exposed_root + chrdtp)
    taper = chrdtp / exposed_root
    tan_c4 = tan_le + 0.25 * (chrdtp - exposed_root) / sspne
    sweep = np.degrees(np.arctan(tan_c4))
    return {
        'area': area, 'aspect_ratio': 4 * sspne**2 / area,
        'taper_ratio': taper, 'taper_ratio_theoretical': chrdtp / chrdr,
        'sweep_c4_deg': sweep, 'cos_c4': np.cos(np.radians(sweep)),
        'tan_c4': tan_c4, 'tan_le': tan_le,
    }


def lift_curve(ang, alpha_zero, slope=0.075, bend=0.0012):
    """A smooth wing lift curve with some curvature, on the SREF basis."""
    d = np.asarray(ang) - alpha_zero
    return slope * d - bend * d * np.abs(d) * 0.1


def base_case():
    alpha = [-4., -2., 0., 2., 4., 8., 12., 16.]
    aliw = 1.0
    geometry = wing_block(6.0, 3.0, 15.0, 13.5, 0.12)
    geometry.update({'alpha_zero_lift': -2.0,
                     'alpha_zero_lift_reference': -2.0,
                     'alpha_clmax_reference': 14.0,
                     'mac_c4_theoretical': 2.6})
    ang = [a + aliw for a in alpha]
    return {
        'alpha_deg': alpha,
        'wing_alone': {'alpha': ang, 'cl': list(lift_curve(ang, -2.0))},
        'wing': {'sspn': 15.0, 'sspnop': 0.0, 'sspndd': 0.0, 'chrdtp': 3.0,
                 'chrdr': 6.0, 'chrdbp': 3.0, 'dhdadi': 0.0, 'dhdado': 0.0,
                 'deltay': 2.5, 'twash': 0.0, 'sspne': 13.5},
        'wing_geometry': geometry,
        'synthesis': {'aliw': aliw, 'xw': 10.0, 'xh': 35.0, 'alih': 0.0},
        'tail': {'sspn': 6.0},
        'tail_geometry': {'tail_arm': 22.5, 'tail_height': 2.0,
                          'a22': 3.1, 'mac_c4_theoretical': 1.1},
        'sref': 150.0, 'kwb': 1.07,
    }


def cases():
    out = []

    # 0: trailing-edge separation, vortex gradient.
    out.append(base_case())

    # 1: leading-edge separation on a swept wing with a sharp section.
    c = base_case()
    c['wing']['deltay'] = 1.2
    c['wing_geometry'].update(wing_block(6.0, 3.0, 15.0, 13.5, 0.75))
    out.append(c)

    # 2: the Section 4.4.1 gradient.
    c = base_case()
    c['wing']['twash'] = 2.0
    out.append(c)

    # 3: cranked planform with dihedral; a wide vortex sheet crosses the
    # break, so the two-panel dihedral drop applies.
    c = base_case()
    c['wing'].update({'sspnop': 6.0, 'chrdbp': 4.0, 'sspndd': 12.0,
                      'dhdadi': 3.0, 'dhdado': 8.0})
    out.append(c)

    # 4: the zero-lift reference below the first tabulated angle, so the
    # lower grid points take the source's LEX=-1 linear scaling.
    c = base_case()
    c['alpha_deg'] = [0., 3., 6., 9.]
    c['synthesis']['aliw'] = 0.0
    c['wing_alone'] = {'alpha': [0., 3., 6., 9.],
                       'cl': list(lift_curve([0., 3., 6., 9.], -3.0))}
    c['wing_geometry'].update({'alpha_zero_lift': -3.0,
                               'alpha_zero_lift_reference': -3.0})
    out.append(c)

    # 5: an angle exactly at zero lift, taking the CLWJ=0 vortex span.
    c = base_case()
    c['alpha_deg'] = [-2., 0., 5.]
    c['synthesis']['aliw'] = 0.0
    c['wing_alone'] = {'alpha': [-2., 0., 5.],
                       'cl': list(lift_curve([-2., 0., 5.], -2.0))}
    out.append(c)

    # 6: B(49) differing from the Mach-zero A(126), separating the
    # integration grid from the lift lookup grid.
    c = base_case()
    c['wing_geometry']['alpha_zero_lift'] = -1.5
    out.append(c)
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        w, g, s, t = (c['wing'], c['wing_geometry'], c['synthesis'],
                      c['tail_geometry'])
        nalpha = len(c['alpha_deg'])
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f"      NALPHA={nalpha}")
        for j, a in enumerate(c['alpha_deg']):
            lines.append(assign(f"FLC({j + 23})", a))
        for j, (a, cl) in enumerate(zip(c['wing_alone']['alpha'],
                                        c['wing_alone']['cl'])):
            lines.append(assign(f"B({j + 23})", a))
            lines.append(assign(f"WING({j + 21})", cl))
        for name, key in [('SSPN', 'sspn'), ('SSPNOP', 'sspnop'),
                          ('SSPNDD', 'sspndd'), ('CHRDTP', 'chrdtp'),
                          ('CHRDR', 'chrdr'), ('CHRDBP', 'chrdbp'),
                          ('DHDADI', 'dhdadi'), ('DHDADO', 'dhdado'),
                          ('DELTAY', 'deltay'), ('TWASH', 'twash'),
                          ('SSPNE', 'sspne')]:
            lines.append(assign(name, w[key]))
        for index, key in [(3, 'area'), (7, 'aspect_ratio'),
                           (27, 'taper_ratio'),
                           (118, 'taper_ratio_theoretical'),
                           (40, 'sweep_c4_deg'), (43, 'cos_c4'),
                           (44, 'tan_c4'), (62, 'tan_le'),
                           (126, 'alpha_zero_lift_reference'),
                           (127, 'alpha_clmax_reference'),
                           (161, 'mac_c4_theoretical')]:
            lines.append(assign(f"A({index})", g[key]))
        lines.append(assign('B(49)', g['alpha_zero_lift']))
        lines.append(assign('A(24)', t['tail_arm']))
        lines.append(assign('A(12)', t['tail_height']))
        lines.append(assign('A(22)', t['a22']))
        lines.append(assign('AHT(161)', t['mac_c4_theoretical']))
        lines.append(assign('HTIN(4)', c['tail']['sspn']))
        lines.append(assign('BD(77)', s['aliw']))
        lines.append(assign('ALIW', s['aliw']))
        lines.append(assign('XW', s['xw']))
        lines.append(assign('XH', s['xh']))
        lines.append(assign('ALIH', s['alih']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('WB(2)', c['kwb']))
        lines.append('      CALL DWASH')
        for tag, array, start in [('ANGLE', 'DWASHI', 20),
                                  ('GRADIENT', 'DWASHI', 40),
                                  ('HEIGHT', 'FACT', 61),
                                  ('SPAN', 'FACT', 81)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"({array}({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'A20',A(20)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('dwash', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('dwash', payload))


if __name__ == '__main__':
    main()
