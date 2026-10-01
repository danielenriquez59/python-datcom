"""
Probe M13O15 (with PRPWEF) and save the fixture.

PRPWEF's ``/WINGI/`` names ``LERI`` and ``LERO``, implicitly INTEGER, ahead
of ``CMOT``, which it reads; in a double-precision build they would take one
word where the reals around them take two, so a layout patch declares them
REAL, which is the source's single-precision storage.

Run from the repository root: ``python test_parity/probes/prpwef.py``.
"""

import copy
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m13o15', 'prpwef', 'angles', 'zerang', 'tbfunx', 'quad',
            'tlinex', 'tlin1x', 'tlin3x', 'glook', 'switch']

LAYOUT = {'prpwef': [('      REAL KN\n',
                      '      REAL KN\n      REAL LERI,LERO,KSHARP\n')]}

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /IWING/  PWING,  WING(400)
      COMMON /IPOWER/ PPOWER, POWER(200)
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /FLGTCD/ FLC(93)
      COMMON /OPTION/ SREF, CBARR, ROUGFC, BLREF
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,SYN(11)
      COMMON /WINGI/  WIN(101)
      COMMON /HTI/    HTIN(131)
      COMMON /POWER/  PWR(28), CROT
      COMMON /WBHCAL/ WBT(155)
      COMMON /WINGD/  A(195), B(49)
      COMMON /WHAERO/ C(51), D(55), CHT(51), DHT(55), DVT(55)
      COMMON /POWR/   PW(285)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC
     1                ,HEAD,PRPOWR,JETPOW,LOASRT
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
      COMMON /BDATA/  BD(762)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /VTDATA/ AVT(195)
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC
     1        ,HEAD,PRPOWR,JETPOW,LOASRT,CROT
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0, 16.0]
POWER_ORDER = ['aietlp', 'nengsp', 'thstcp', 'phaloc', 'phvloc', 'prprad',
               'kn', 'bwapr3', 'bwapr6', 'bwapr9', 'nopbpe', 'bapr75']


def case(n, htpl=True, twista=-2.0, **power):
    rng = np.random.default_rng(1300 + n)
    alpha = np.array(ALPHA)
    p = {'aietlp': 2.0, 'nengsp': 1.0, 'thstcp': 0.15, 'phaloc': 14.0,
         'phvloc': 0.2, 'prprad': 4.0, 'kn': 80.0, 'bwapr3': 0.1,
         'bwapr6': 0.08, 'bwapr9': 0.05, 'nopbpe': 3.0, 'bapr75': 30.0,
         'yp': 0.0, 'crot': False}
    p.update(power)
    return {
        'alpha': ALPHA, 'sref': 300.0, 'cbarr': 6.0, 'htpl': htpl,
        'power': p,
        'wing': {'ct': 2.5, 'bst0o2': 8.0, 'bsto2': 15.0, 'bo2': 17.0,
                 'cb': 5.0, 'cr': 8.0, 'twista': twista, 'cmo': -0.03,
                 'cmot': -0.02, 'cl': list(0.075 * (alpha + 2.0)),
                 'cla': 0.075, 'cd0': 0.009,
                 'cdl': list(0.0004 * (alpha + 2.0)**2 + 1e-5),
                 'd10': 0.0031, 'd11': 0.0012, 'd12': 0.0036},
        'a': {1: 90.0, 2: 40.0, 3: 130.0, 10: 7.5, 16: 5.5, 23: 7.0,
              26: 0.5, 32: 0.3, 33: 9.5, 62: 0.6, 67: 0.85, 69: 0.55,
              86: 0.4, 91: 0.9, 106: 30.0, 112: 20.0, 120: 4.5,
              121: 6.5, 134: -2.0, 161: 1.8},
        'position': {'xcg': 22.0, 'xw': 18.0, 'zw': 0.0, 'aliw': 1.0,
                     'zcg': 0.5, 'xh': 45.0, 'zh': 1.5, 'alih': -1.0},
        'tail': {'bo2': 6.0, 'cr': 3.5, 'ct': 1.8, 'area': 60.0,
                 'xbarr': 1.2, 'dht10': 0.0028, 'dht11': 0.001,
                 'cl': list(0.05 * (alpha + 1.0) + rng.uniform(0, .01, 6))},
        'vertical': {'area': 45.0, 'dvt10': 0.003, 'dvt11': 0.001},
        'body': {'cf': 0.0025, 'wetted_area': 900.0},
        'dwash': {'q_ratio': list(0.95 - 0.003 * alpha),
                  'epsilon': list(0.35 * alpha + 1.0)},
        'stale': {'sih': 7.5, 'ytemp': 0.12, 'dlh': 3.0},
    }


def cases():
    return [
        case(0),
        case(1, nengsp=2.0, yp=11.0),                # outboard, 2 engines
        case(2, nengsp=2.0, yp=5.5, htpl=False),     # inboard, no tail
        case(3, htpl=False),                          # stale SIH, YTEMP
        case(4, kn=UNUSED, crot=True, nopbpe=6.0, bapr75=35.0),
        case(5, phaloc=21.0),                         # XBARP/CRP < 0.25
        case(6, twista=1.5),                          # positive twist
        case(7, prprad=10.5, phvloc=0.0),             # past the break
        case(8, prprad=6.5, nengsp=3.0, yp=0.0),      # several engines
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,285')
        lines.append('         PW(K)=0.')
        lines.append('         IF(K.LE.200) POWER(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        na = len(c['alpha'])
        lines.append(assign('FLC(2)', na))
        lines.append(f'      NALPHA={na}')
        p, w, a, pos = c['power'], c['wing'], c['a'], c['position']
        for k, name in enumerate(POWER_ORDER):
            lines.append(assign(f'PWR({k + 1})', p[name]))
        lines.append(assign('PWR(28)', p['yp']))
        lines.append(f"      CROT=.{'TRUE' if p['crot'] else 'FALSE'}.")
        for k, name in enumerate(['ct', 'bst0o2', 'bsto2', 'bo2', 'cb',
                                  'cr']):
            lines.append(assign(f'WIN({k + 1})', w[name]))
        lines.append(assign('WIN(11)', w['twista']))
        lines.append(assign('WIN(61)', w['cmo']))
        lines.append(assign('WIN(67)', w['cmot']))
        lines.append(assign('WING(101)', w['cla']))
        lines.append(assign('B(46)', w['cd0']))
        for index, value in a.items():
            lines.append(assign(f'A({index})', value))
        for name, key in [('D(10)', 'd10'), ('D(11)', 'd11'),
                          ('D(12)', 'd12')]:
            lines.append(assign(name, w[key]))
        for name in ('XCG', 'XW', 'ZW', 'ALIW', 'ZCG', 'XH', 'ZH', 'ALIH'):
            lines.append(assign(name, pos[name.lower()]))
        t, v, b = c['tail'], c['vertical'], c['body']
        for name, value in [('HTIN(4)', t['bo2']), ('HTIN(6)', t['cr']),
                            ('HTIN(1)', t['ct']), ('AHT(3)', t['area']),
                            ('AHT(161)', t['xbarr']),
                            ('DHT(10)', t['dht10']), ('DHT(11)', t['dht11']),
                            ('AVT(3)', v['area']), ('DVT(10)', v['dvt10']),
                            ('DVT(11)', v['dvt11']), ('BD(92)', b['cf']),
                            ('BD(93)', b['wetted_area']),
                            ('SREF', c['sref']), ('CBARR', c['cbarr']),
                            ('PW(204)', c['stale']['sih']),
                            ('PW(257)', c['stale']['ytemp']),
                            ('PW(249)', c['stale']['dlh'])]:
            lines.append(assign(name, value))
        for j in range(na):
            lines.append(assign(f'FLC({23 + j})', c['alpha'][j]))
            lines.append(assign(f'WING({21 + j})', w['cl'][j]))
            lines.append(assign(f'D({36 + j})', w['cdl'][j]))
            lines.append(assign(f'WBT({110 + j})', t['cl'][j]))
            lines.append(assign(f'DWASH({1 + j})', c['dwash']['q_ratio'][j]))
            lines.append(assign(f'DWASH({21 + j})', c['dwash']['epsilon'][j]))
        lines.append('      CALL M13O15')
        lines.append("      WRITE(6,'(A,285ES25.16)') 'PW',(PW(J),J=1,285)")
        lines.append("      WRITE(6,'(A,140ES25.16)') 'POWER',"
                     "(POWER(J),J=1,140)")
        lines.append("      WRITE(6,'(A,3ES25.16,I4)') 'OUT',PWR(7),BD(79),"
                     "FLC(2),NALPHA")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('prpwef', driver(all_cases), ROUTINES,
                                layout_patches=LAYOUT))
    payload = []
    for c, r in zip(all_cases, records):
        c = copy.deepcopy(c)
        c['a'] = {str(k): v for k, v in c['a'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('prpwef', payload))


if __name__ == '__main__':
    main()
