"""
Probe WBAERO's horizontal-tail and body-vertical passes and save the
fixture.  The wing pass has its own probe, test_parity/probes/wbaero.py.

Run from the repository root: ``python test_parity/probes/wbaero_tail.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import assign, parse_records, run, save  # noqa: E402
from wbaero import COMMONS, EXTRA_STUBS, ROUTINES, UNUSED  # noqa: E402


def case(rng, alpha=None, sspne=5.0, alih=0.0, bd63=-17.5,
         missing_cm=None):
    alpha = np.array(alpha if alpha is not None else
                     [-4., 0., 4., 8., 12., 16., 20.])
    local = alpha + alih
    zero = -0.4
    cl = 0.055 * (local - zero) * np.where(local > 14, 0.9, 1.0)
    cm = -0.004 - 0.01 * cl
    if missing_cm is not None:
        cm[missing_cm:] = 2 * UNUSED
    x = np.linspace(0.0, 44.0, 12)
    radius = 2.1 * np.sin(np.pi * np.clip(x / 48.0, 0, 1))**0.5
    body_zero = -0.3
    return {
        'alpha_deg': list(alpha),
        'surface': {
            'sspn': 6.0, 'sspne': sspne, 'chrdr': 3.5, 'twista': 0.0,
            'tovc': 0.08, 'ler': 0.006, 'ycm': 0.0, 'cld': 0.0,
            'a7': 3.8, 'a10': 3.1, 'a27': 0.5, 'a38': 0.4, 'a44': 0.33,
            'a62': 0.4, 'a80': -0.05, 'a118': 0.5, 'a120': 4.0,
            'a122': 2.6, 'a129': 2.0e6, 'a160': rng.uniform(0.5, 6.0),
            'a161': 1.3, 'a173': rng.uniform(-18.0, -15.0),
            'beta': 0.87, 'cm0': -0.004, 'alpha_zero_lift': zero,
            'clmax': 1.0, 'alpha_clmax': 16.0, 'local_alpha': list(local),
            'cd0': 0.006, 'cdl': list(0.04 * cl**2), 'c6': 0.26,
        },
        'tail_alone': {'cd': list(0.006 + 0.04 * cl**2), 'cl': list(cl),
                       'cm': list(cm), 'cn': list(cl / np.cos(np.radians(local))),
                       'ca': list(0.005 - 0.01 * cl), 'cla': 0.055},
        'body': {
            'cd': list(0.01 + 0.0001 * alpha**2),
            'cl': list(0.004 * (alpha - body_zero)),
            'cm': list(0.006 * (alpha - body_zero)),
            'cla': list(0.004 + 0.0001 * np.abs(alpha)),
            'cma': list(0.006 + 0.0002 * np.abs(alpha)),
            'cm0': 0.002, 'alpha_zero_lift': body_zero,
            'cd_friction': 0.009, 'cd_base': 0.0015,
            'cdl': list(0.0005 * alpha**2 / 10), 'x': list(x),
            's': list(np.pi * radius**2),
        },
        'synthesis': {'xcg': 21.0, 'xh': 38.0, 'zh': 1.5, 'zcg': 0.0,
                      'alih': alih},
        'flight': {'mach': 0.5, 'reynolds_per_length': 2.0e6, 'tr': 0.4},
        'cbarr': 5.2, 'bd63': bd63, 'vertical_cd0': 0.0034,
    }


def cases():
    rng = np.random.default_rng(707)
    return [
        case(rng),
        case(rng, alih=-2.0),
        case(rng, sspne=5.6),
        case(rng, missing_cm=4),
        case(rng, bd63=0.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             '      WGPL=.FALSE.', '      HTPL=.TRUE.', '      BO=.TRUE.',
             '      VTPL=.TRUE.', '      I=1', '      NF=-1',
             '      KBODY=.FALSE.', '      KHT=.FALSE.']
    for n, c in enumerate(all_cases):
        sf, ta, bd, sy, fl = (c['surface'], c['tail_alone'], c['body'],
                              c['synthesis'], c['flight'])
        alpha = c['alpha_deg']
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.400) BODY(K)=0.')
        lines.append('         IF(K.LE.380) HT(K)=0.')
        lines.append('         IF(K.LE.380) BH(K)=0.')
        lines.append('         IF(K.LE.380) BV(K)=0.')
        lines.append('         IF(K.LE.195) AHT(K)=0.')
        lines.append('         IF(K.LE.182) FACT(K)=0.')
        lines.append('         IF(K.LE.160) FLC(K)=0.')
        lines.append('         IF(K.LE.154) HTIN(K)=0.')
        lines.append('         IF(K.LE.55) DHT(K)=0.')
        lines.append('         IF(K.LE.51) CHT(K)=0.')
        lines.append('         IF(K.LE.49) BHT(K)=0.')
        lines.append('         IF(K.LE.39) HB(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f'      NALPHA={len(alpha)}')
        for j, a in enumerate(alpha):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'BHT({23 + j})', sf['local_alpha'][j]))
            lines.append(assign(f'BD({255 + j})', a + bd['alpha_zero_lift']))
            lines.append(assign(f'DHT({36 + j})', sf['cdl'][j]))
            lines.append(assign(f'BD({215 + j})', bd['cdl'][j]))
            for base, key in [(0, 'cd'), (20, 'cl'), (40, 'cm'),
                              (100, 'cla'), (120, 'cma')]:
                lines.append(assign(f'BODY({base + 1 + j})', bd[key][j]))
            for base, key in [(0, 'cd'), (20, 'cl'), (40, 'cm'), (60, 'cn'),
                              (80, 'ca')]:
                lines.append(assign(f'HT({base + 1 + j})', ta[key][j]))
        lines.append(assign('HT(101)', ta['cla']))
        lines.append(assign('FLC(3)', fl['mach']))
        lines.append(assign('FLC(43)', fl['reynolds_per_length']))
        lines.append(assign('FLC(96)', fl['tr']))
        for key in ('a7', 'a10', 'a27', 'a38', 'a44', 'a62', 'a80', 'a118',
                    'a120', 'a122', 'a129', 'a160', 'a161', 'a173'):
            lines.append(assign(f'AHT({key[1:]})', sf[key]))
        lines.append(assign('BHT(1)', fl['mach']))
        for key, index in [('beta', 2), ('alpha_clmax', 43), ('clmax', 44),
                           ('cm0', 47), ('alpha_zero_lift', 49)]:
            lines.append(assign(f'BHT({index})', sf[key]))
        for key, index in [('sspne', 3), ('sspn', 4), ('chrdr', 6),
                           ('twista', 11), ('tovc', 16), ('ler', 62),
                           ('ycm', 93), ('cld', 94)]:
            lines.append(assign(f'HTIN({index})', sf[key]))
        lines.append(assign('DHT(20)', sf['cd0']))
        lines.append(assign('CHT(6)', sf['c6']))
        lines.append(assign('BD(59)', bd['cd_friction']))
        lines.append(assign('BD(60)', bd['cd_base']))
        lines.append(assign('BD(62)', bd['cm0']))
        lines.append(assign('BD(63)', c['bd63']))
        lines.append(assign('BD(81)', bd['alpha_zero_lift']))
        lines.append(f"      XNX={len(bd['x'])}.")
        for k, (xv, sv) in enumerate(zip(bd['x'], bd['s'])):
            lines.append(assign(f'X({k + 1})', xv))
            lines.append(assign(f'S({k + 1})', sv))
        lines.append(assign('BD(1)', bd['x'][-1]))
        for key, name in [('xcg', 'XCG'), ('xh', 'XH'), ('zh', 'ZH'),
                          ('zcg', 'ZCG'), ('alih', 'ALIH')]:
            lines.append(assign(name, sy[key]))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append(assign('DVT(20)', c['vertical_cd0'] * 0.8))
        lines.append(assign('DVF(20)', c['vertical_cd0'] * 0.2))
        lines.append('      CALL WBAERO')
        for block in ('BH', 'BV'):
            for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                               ('CA', 80), ('CLA', 100), ('CMA', 120)]:
                lines.append(f"      WRITE(6,'(A,30ES25.16)') '{block}{tag}',"
                             f"({block}({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,39ES25.16)') 'HB',(HB(J),J=1,39)")
        lines.append("      WRITE(6,'(A,21ES25.16)') 'FACT',"
                     "(FACT(141+J),J=1,21)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('wbaero_tail', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('wbaero_tail', payload))


if __name__ == '__main__':
    main()
