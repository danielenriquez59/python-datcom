"""
Probe TABLES (directly, over a grid) and WBCD (both halves) and save the
fixture.  WBCD's regression is WBCDL, which reads TBSUB, TBTRN and TBSUP
through TABLES.

Run from the repository root: ``python test_parity/probes/wbcd.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['wbcd', 'wbcdl', 'tables', 'tbsub', 'tbtrn', 'tbsup']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(95), TR
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195)
      COMMON /WHWB/   FACT(182), WB(39), HB(39)
      COMMON /IHT/    PHT,HT(380)
      COMMON /HTI/    HTIN(131)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(100)
      COMMON /IBH/    PBH,BH(380)
      COMMON /SYNTSS/ XCG, XW, ZW, ALIW, ZCG, XH, ZH, ALIH, XV,
     1                VERTUP, HINAX, XVF, SCALE, ZV, ZVF, YV, YF,
     2                PHIV, PHIF
      COMMON /IBW/    PBW, BW(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,VERTUP
      DIMENSION BTAB(16),WT(4)
      LOGICAL NDM
"""

UNUSED = 1.0e-30


def table_points():
    rng = np.random.default_rng(1616)
    pts = []
    for mach in (0.0, 0.3, 0.6, 0.75, 0.9, 0.93, 0.95, 1.0, 1.05, 1.1, 1.25,
                 1.5, 1.8, 2.5, 2.6, -0.1):
        for alpha in (0.0, 0.4, 3.0, 7.5, -9.2, 11.0, 12.0, 14.9, 15.0,
                      17.3, 18.0, 18.5):
            pts.append((mach, alpha))
    for _ in range(40):
        pts.append((rng.uniform(0.0, 2.5), rng.uniform(-18.0, 18.0)))
    return [(float(m), float(a)) for m, a in pts]


def surface(kind=1.0, sspn=15.0, sspne=13.0, twist=-2.0, ar=4.0):
    return {'type': kind, 'sspn': sspn, 'sspne': sspne, 'chrdr': 7.0,
            'tovc': 0.06, 'ler': 0.008, 'twista': twist, 'ycm': 0.01,
            'cld': 0.2, 'a38': 0.5, 'a56': -0.1, 'a118': 0.4, 'a120': ar,
            'a122': 5.0}


def case(mach=0.5, wing=None, tail=None, xmax=20.0, alpha=None):
    alpha = alpha or [-4., 0., 4., 8., 12., 16.]
    return {
        'alpha': alpha, 'mach': mach, 'rn': 2.0e6, 'tr': 0.4,
        'wing': wing or surface(), 'tail': tail or surface(
            sspn=6.0, sspne=5.3, ar=3.5),
        'xw': 16.0, 'xh': 38.0,
        'body': {'length': 44.0, 'x_max_area': xmax, 'max_diameter': 4.2},
        'wing_body': {'cd0': 0.021, 'cl': list(0.08 * (np.array(alpha) + 2))},
        'tail_body': {'cd0': 0.009, 'cl': list(0.01 * np.array(alpha))},
    }


def cases():
    return [
        case(),
        case(mach=0.95),
        case(mach=1.2, alpha=[-2., 2., 6., 10., 14., 16.]),
        case(xmax=10.0),                            # BD(2) ahead of the LE
        case(xmax=30.0),                            # BD(2) behind the TE
        case(wing=surface(kind=3.0)),               # not straight: skipped
        case(wing=surface(ar=7.0)),                 # out of range
        case(wing=surface(twist=2.0)),              # positive twist: out
        case(mach=0.95, alpha=[12.0, 2.0, 4.0]),    # first angle has no data
        # A tail placed so its nose and afterbody lengths fall inside the
        # regression's range; an aft tail never does.
        dict(case(tail=surface(sspn=6.0, sspne=4.0, ar=3.5)), xh=15.0),
    ]


def driver(points, all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1',
             '      DO 10 K=1,16', '   10 BTAB(K)=-7.', ]
    for n, (mach, alpha) in enumerate(points):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(assign('XM', mach))
        lines.append(assign('XA', alpha))
        lines.append('      CALL TABLES(BTAB,XM,XA,NDM)')
        lines.append("      WRITE(6,'(A,L2)') 'NDM ',NDM")
        lines.append("      WRITE(6,'(A,16ES25.16)') 'B',BTAB")
    offset = len(points)
    for n, c in enumerate(all_cases):
        label = 2000 + n
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{offset + n}")
        lines.append(f'      DO {label} K=1,380')
        lines.append('         BW(K)=0.')
        lines.append('         BH(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append('      WGPL=.TRUE.')
        lines.append('      HTPL=.TRUE.')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'BW({21 + j})', c['wing_body']['cl'][j]))
            lines.append(assign(f'BH({21 + j})', c['tail_body']['cl'][j]))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['rn']))
        lines.append(assign('TR', c['tr']))
        for block, win, a, key in [('wing', 'WINGIN', 'A', 'wing'),
                                   ('tail', 'HTIN', 'AHT', 'tail')]:
            s = c[key]
            lines.append(f"      {win}(15)=WT({int(s['type'])})")
            for index, name in [(4, 'sspn'), (3, 'sspne'), (6, 'chrdr'),
                                (16, 'tovc'), (62, 'ler'), (11, 'twista'),
                                (93, 'ycm'), (94, 'cld')]:
                lines.append(assign(f'{win}({index})', s[name]))
            for index in (38, 56, 118, 120, 122):
                lines.append(assign(f'{a}({index})', s[f'a{index}']))
        lines.append(assign('WB(17)', c['wing_body']['cd0']))
        lines.append(assign('HB(17)', c['tail_body']['cd0']))
        lines.append(assign('XW', c['xw']))
        lines.append(assign('XH', c['xh']))
        lines.append(assign('BD(1)', c['body']['length']))
        lines.append(assign('BD(2)', c['body']['x_max_area']))
        lines.append(assign('BD(85)', c['body']['max_diameter']))
        lines.append('      CALL WBCD')
        for arr in ('BW', 'BH'):
            for tag, start in [('CD', 0), ('CN', 60), ('CA', 80)]:
                lines.append(f"      WRITE(6,'(A,30ES25.16)') '{arr}{tag}',"
                             f"({arr}({start}+J),J=1,NALPHA)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    points, all_cases = table_points(), cases()
    records = parse_records(run('wbcd', driver(points, all_cases), ROUTINES,
                                EXTRA_STUBS))
    tables = [{'mach': m, 'alpha': a, 'ndm': r['_text'][0].split()[-1] == 'T',
               'b': r['B']} for (m, a), r in zip(points, records)]
    wbcd = [{'inputs': c, 'outputs': r}
            for c, r in zip(all_cases, records[len(points):])]
    print(save('wbcd', {'tables': tables, 'wbcd': wbcd}))


if __name__ == '__main__':
    main()
