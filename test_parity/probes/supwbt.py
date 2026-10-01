"""
Probe SUPWBT, the supersonic wing-body-tail, and save the fixture.

Run from the repository root: ``python test_parity/probes/supwbt.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['supwbt', 'intkbw', 'quadin', 'bodowg', 'getmax', 'ali',
            'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'tbfunx', 'quad']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /OVERLY/ NLOG,NMACH,I,NALPH
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /BDATA/  BD(762)
      COMMON /SUPWBB/ SWB(61)
      COMMON /SUPDW/  DWA(237)
      COMMON /SUPWH/  SLG(141),STG(141)
      COMMON /WBHCAL/ STP(156)
      COMMON /SUPBOD/ SBD(227)
      COMMON /WHWB/   FACT(182)
      COMMON /FLGTCD/ FLC(95)
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /SYNTSS/ SYNA(19)
      COMMON /BODYI/  NXX,XCOOR(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1                BNOSE,BTAIL,BLN,BLA,DS
      COMMON /WINGI/  WINGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /IBW/    PBW, BW(380)
      COMMON /IBWH/   PBWH, BWH(380)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IBWHV/  PBWHV, BWHV(380)
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /IWING/  PWING, WING(380)
      REAL NXX
      EQUIVALENCE (JDETCH,DWA(237))
"""

UNUSED = 1.0e-30
RAD = 57.2957795
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0, 16.0]
BODY_X = [0.0, 5.0, 12.0, 25.0, 38.0, 48.0]
BODY_S = [0.0, 5.0, 11.0, 13.2, 12.0, 9.0]


def case(n, mach=1.8, jdetch=-1, user=False, canard=False, alih=-1.5,
         tail_taper=0.4, tail_ar=2.5, tanle=None, xh=40.0, trunc=False,
         alpha=None, eps=None):
    rng = np.random.default_rng(3100 + n)
    alpha = list(alpha or ALPHA)
    na = len(alpha)
    al = np.array(alpha)
    eps = list(eps if eps is not None else 0.25 * al + 0.3)
    tanle = math.tan(math.radians(40.0)) if tanle is None else tanle
    tcl = list(0.05 * (al - 1.0))
    wcl = list(0.09 * (al + 0.5))
    if trunc:
        tcl[-1] = UNUSED
        wcl[-1] = UNUSED
    return {
        'alpha': alpha, 'mach': mach, 'sref': 250.0, 'cbarr': 7.0,
        'user': user,
        'position': {'xcg': 21.0, 'aliw': 1.5, 'zcg': 0.3, 'xh': xh,
                     'zh': 1.2, 'alih': alih},
        'dwash': {'q': list(0.95 - 0.002 * al),
                  'eps': eps,
                  'dedalp': list(0.25 + 0.01 * rng.uniform(size=na)),
                  'hmach': list(mach - 0.05 - 0.01 * rng.uniform(size=na)),
                  'jdetch': jdetch},
        'tail': {'span': 5.5, 'spans': 3.7,
                 'a': {3: tail_ar * 5.5 * 1.8, 7: tail_ar, 10: 3.8,
                       27: tail_taper, 34: 40.0, 38: tanle, 62: tanle,
                       161: 1.6},
                 'cd': list(0.008 + 0.0004 * al**2),
                 'cl': tcl,
                 'cla': list(0.05 + 0.001 * rng.uniform(size=na)),
                 'alpha': list((al - 0.5) / RAD),
                 'cd0': 0.007, 'xac': 0.45, 'xac136': 0.6},
        'wing': {'span': 12.0, 'spans': 10.2, 'type': 3.0 if canard else 1.0,
                 'a': {7: 2.8, 10: 8.5, 12: 1.4, 24: 22.0, 27: 0.3,
                       38: 1.0, 80: 0.3, 120: 2.8},
                 'cla': 0.058, 'xac': 0.5, 'xac136': 0.4, 'kwb': 1.12,
                 'kkwb': 0.95, 'kkbw': 0.2, 'cd0': 0.02},
        'wing_body': {'cd': list(0.02 + 0.0008 * al**2), 'cl': wcl,
                      'cla': 0.09, 'cma': -0.02},
        'body': {'x': BODY_X, 's': BODY_S, 'dn': 3.6, 'd1': 4.0,
                 'bd68': -0.4, 'bd78': 0.15},
        'stale': {'cd0v': 0.003, 'cd0vf': 0.001, 'kkbw': 0.33,
                  'kkwb': 0.77},
    }


def cases():
    return [
        case(0),
        case(1, user=True, jdetch=0),                 # experimental data
        case(2, jdetch=0),                            # returns at once
        case(3, jdetch=4),                            # detachment
        case(4, canard=True),
        case(5, canard=True, alpha=[-4.0, 0.0, 2.0, 8.0, 12.0],
             eps=[0.5, 0.8, 0.5, 2.0, 3.0], alih=-1.5),   # ALPAHT = 0
        case(6, alih=0.0),                            # stale KK
        case(7, tail_taper=0.0, tail_ar=0.5),         # triangular, Fig. 10
        case(8, tail_taper=0.0, tail_ar=2.0),         # triangular, INTKBW
        case(9, mach=1.15, tail_ar=1.0),              # TRINO <= 4
        case(10, xh=47.0),                            # DX <= -CR
        case(11, tanle=0.0),                          # TANLE snapped
        case(12, trunc=True),                         # NA1
        case(13, alih=UNUSED, user=True),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        label = 3000 + n
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        for arr in ('BW', 'BWH', 'HT', 'BWHV', 'WING'):
            lines.append(f'         IF(K.LE.380) {arr}(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.195) AHT(K)=0.')
        lines.append('         IF(K.LE.237) DWA(K)=0.')
        lines.append('         IF(K.LE.182) FACT(K)=0.')
        lines.append('         IF(K.LE.156) STP(K)=0.')
        lines.append('         IF(K.LE.141) SLG(K)=0.')
        lines.append('         IF(K.LE.141) STG(K)=0.')
        lines.append('         IF(K.LE.60) DWASH(K)=0.')
        lines.append(f' {label} CONTINUE')
        na = len(c['alpha'])
        lines.append(f'      NALPH={na}')
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        p = c['position']
        for index, key in [(1, 'xcg'), (4, 'aliw'), (5, 'zcg'), (6, 'xh'),
                           (7, 'zh'), (8, 'alih')]:
            lines.append(assign(f'SYNA({index})', p[key]))
        d = c['dwash']
        for j in range(na):
            lines.append(assign(f'DWASH({1 + j})', d['q'][j]))
            lines.append(assign(f'DWASH({21 + j})', d['eps'][j]))
            lines.append(assign(f'DWASH({41 + j})', d['dedalp'][j]))
            lines.append(assign(f'DWA({189 + j})', d['hmach'][j]))
        lines.append(f"      JDETCH={d['jdetch']}")
        t = c['tail']
        lines.append(assign('HTIN(4)', t['span']))
        lines.append(assign('HTIN(3)', t['spans']))
        for k, v in t['a'].items():
            lines.append(assign(f'AHT({k})', v))
        for j in range(na):
            lines.append(assign(f'HT({1 + j})', t['cd'][j]))
            lines.append(assign(f'HT({21 + j})', t['cl'][j]))
            lines.append(assign(f'HT({101 + j})', t['cla'][j]))
            lines.append(assign(f'STG({33 + j})', t['alpha'][j]))
        lines.append(assign('STG(80)', t['cd0']))
        lines.append(assign('STG(134)', t['xac']))
        lines.append(assign('STG(136)', t['xac136']))
        w = c['wing']
        lines.append(assign('WINGIN(4)', w['span']))
        lines.append(assign('WINGIN(3)', w['spans']))
        lines.append(assign('WINGIN(101)', w['type']))
        for k, v in w['a'].items():
            lines.append(assign(f'A({k})', v))
        for name, key in [('WING(101)', 'cla'), ('SLG(134)', 'xac'),
                          ('SLG(136)', 'xac136'), ('SWB(35)', 'kwb'),
                          ('SWB(2)', 'kkwb'), ('SWB(37)', 'kkbw'),
                          ('SWB(4)', 'cd0')]:
            lines.append(assign(name, w[key]))
        wb = c['wing_body']
        for j in range(na):
            lines.append(assign(f'BW({1 + j})', wb['cd'][j]))
            lines.append(assign(f'BW({21 + j})', wb['cl'][j]))
        lines.append(assign('BW(101)', wb['cla']))
        lines.append(assign('BW(121)', wb['cma']))
        b = c['body']
        lines.append(f"      NXX={len(b['x'])}.")
        for k, (x, s) in enumerate(zip(b['x'], b['s'])):
            lines.append(assign(f'XCOOR({k + 1})', x))
            lines.append(assign(f'S({k + 1})', s))
        for name, key in [('SBD(4)', 'dn'), ('SBD(5)', 'd1'),
                          ('BD(68)', 'bd68'), ('BD(78)', 'bd78')]:
            lines.append(assign(name, b[key]))
        st = c['stale']
        for name, key in [('STP(1)', 'cd0v'), ('STP(156)', 'cd0vf'),
                          ('STP(131)', 'kkbw'), ('STP(132)', 'kkwb')]:
            lines.append(assign(name, st[key]))
        flag = 1 if c['user'] else 0
        lines.append(f'      CALL SUPWBT({flag},{flag},0)')
        for tag, text in [
                ('STP', '(STP(J),J=1,156)'), ('BWH', '(BWH(J),J=1,140)'),
                ('BWHV', '(BWHV(J),J=1,140)'),
                ('FACT', '(FACT(J),J=102,141)'),
                ('BD', 'BD(2),BD(3),BD(58),BD(63),BD(64),BD(84),BD(761),'
                       'BD(762)'),
                ('AHT', 'AHT(62)')]:
            lines.append(f"      WRITE(6,'(A,156ES25.16)') '{tag}',")
            text = text.replace(',BD(761)', ',\n     2BD(761)')
            lines.append(f"     1{text}")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('supwbt', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('tail', 'wing'):
            c[key]['a'] = {str(k): v for k, v in c[key]['a'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('supwbt', payload))


if __name__ == '__main__':
    main()
