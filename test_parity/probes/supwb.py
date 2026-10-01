"""
Probe M20O24 (SUPWB, SUPHB, SUPCM0 and the slope pass) and save the
fixture.  The vertical panels are left out (VRTCDO and VFCDO are stubbed,
their flags off).

Run from the repository root: ``python test_parity/probes/supwb.py``.
"""

import copy
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m20o24', 'supwb', 'suphb', 'supcm0', 'wbcm0', 'tablec',
            'intkbw', 'quadin', 'bodowg', 'getmax', 'ali', 'interx',
            'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch', 'tbfunx',
            'quad']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE VRTCDO(I)
      RETURN
      END
      SUBROUTINE VFCDO(I)
      RETURN
      END
"""

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /FLGTCD/ FLC(160)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /WINGI/  WIN(101)
      COMMON /HTI/    HIN(154)
      COMMON /BODYI/  NXX,XCOOR(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1                BNOSE,BTAIL,RLN
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195), B(49)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,L,NOVLY
      COMMON /SUPBOD/ SBD(229)
      COMMON /SUPWBB/ SWB(61), SHB(61)
      COMMON /SUPWH/  SLG(141), STG(141)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108), STBH(135)
      COMMON /WBHCAL/ STP(156)
      COMMON /IBODY/  PBODY, BODY(400)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IVT/    PVT, VT(380)
      COMMON /IVF/    PVF, VF(380)
      COMMON /IBW/    PBW, BW(380)
      COMMON /IBH/    PBH, BH(380)
      COMMON /IBV/    PBV, BV(380)
      COMMON /IBWV/   PBWV, BWV(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
      REAL NXX
      DIMENSION WT(4)
"""

UNUSED = 1.0e-30
ALPHA = [-2.0, 0.0, 2.0, 4.0, 8.0, 12.0]
BODY_X = [0.0, 5.0, 12.0, 25.0, 38.0, 48.0]
BODY_S = [0.0, 5.0, 11.0, 13.2, 12.0, 9.0]


def surface(n, span=12.0, spans=10.2, cr=9.0, ar=2.8, taper=0.3,
            sweep=45.0, kind=1, a173=3.0, tovc=0.04):
    rng = np.random.default_rng(2000 + n)
    alpha = np.array(ALPHA)
    tanle = math.tan(math.radians(sweep))
    cla = 0.045 + 0.005 * ar
    return {
        'in': {3: spans, 4: span, 6: cr, 11: -1.0, 15: kind, 16: tovc,
               62: 0.005, 93: 0.005, 94: 0.1},
        'a': {3: ar * span, 7: ar, 10: cr * 0.95, 27: taper, 34: sweep,
              38: tanle, 62: tanle, 80: 0.3, 118: taper, 120: ar,
              122: 0.7 * cr, 161: 0.4 * cr, 173: a173},
        'aero': {'cla': cla, 'cl': list(cla * (alpha + 0.5)),
                 'cn': list(cla * (alpha + 0.5) * 1.01),
                 'ca': list(0.01 + 0.0001 * alpha**2),
                 'cd0': 0.006, 'cdl': list(0.0005 * alpha**2),
                 'xac': float(rng.uniform(0.4, 0.7))},
    }


def case(n, mach=1.6, wing=None, tail=None, xw=18.0, xh=40.0, aliw=1.5,
         alih=-1.0, bnose=0.0, rln=12.0, htpl=True, wgpl=True):
    return {
        'alpha': ALPHA, 'mach': mach, 'wgpl': wgpl, 'htpl': htpl,
        'wing': wing or surface(n),
        'tail': tail or surface(n + 50, span=5.5, spans=3.7, cr=4.0,
                                ar=2.5, taper=0.4, sweep=40.0, a173=-18.0),
        'position': {'xw': xw, 'xh': xh, 'aliw': aliw, 'alih': alih,
                     'zw': -0.8, 'zcg': 0.3, 'xcg': 21.0},
        'body': {'x': BODY_X, 's': BODY_S, 'rln': rln, 'bnose': bnose,
                 'dn': 3.6, 'd1': 4.0, 'alpha0': UNUSED, 'cla': 0.0045,
                 'cd0': 0.03, 'sbd120': 4.1,
                 'cl': list(0.0045 * np.array(ALPHA)),
                 'cd': list(0.03 + 0.0002 * np.array(ALPHA)**2),
                 'cm': list(0.02 * np.array(ALPHA) / 57.3)},
        'sref': 250.0, 'cbarr': 7.0, 'tr': 0.4,
        'stale': {'kkwb': 0.9, 'kkbw': 0.15, 'hkkwb': 0.8, 'hkkbw': 0.1,
                  'bw41': -0.02, 'bh41': -0.01},
    }


def cases():
    return [
        case(0),
        case(1, mach=2.4),
        case(2, mach=1.2, wing=surface(2, ar=1.8, taper=0.0)),    # untapered
        case(3, mach=1.3, wing=surface(3, ar=1.0, taper=0.5,
                                       sweep=60.0)),              # Fig. 10
        case(4, bnose=1.0),                                       # blunt nose
        case(5, rln=22.0),                                        # RLAP < 0
        case(6, aliw=0.0, alih=UNUSED),                           # stale KK
        case(7, wing=surface(7, kind=3)),                         # cranked
        case(8, xw=36.0, xh=44.0),                                # TE past body
        case(9, mach=3.0, wing=surface(9, sweep=20.0, ar=3.5)),   # no DX
        case(10, htpl=False),
        case(11, xw=22.0),                                        # CM0 in range
        case(12, xw=47.5, xh=47.0),                               # DX <= -CR
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1',
             '      BO=.TRUE.', '      VTPL=.FALSE.', '      VFPL=.FALSE.']
    for n, c in enumerate(all_cases):
        label = 3000 + n
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        for arr in ('BODY', 'WING'):
            lines.append(f'         IF(K.LE.400) {arr}(K)=0.')
        for arr in ('HT', 'VT', 'VF', 'BW', 'BH', 'BV', 'BWV'):
            lines.append(f'         IF(K.LE.380) {arr}(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.195) AHT(K)=0.')
        lines.append('         IF(K.LE.61) SWB(K)=0.')
        lines.append('         IF(K.LE.61) SHB(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      WGPL=.{'TRUE' if c['wgpl'] else 'FALSE'}.")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(96)', c['tr']))
        pos = c['position']
        for index, key in [(1, 'xcg'), (2, 'xw'), (3, 'zw'), (4, 'aliw'),
                           (5, 'zcg'), (6, 'xh'), (8, 'alih')]:
            lines.append(assign(f'SYNA({index})', pos[key]))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        b = c['body']
        lines.append(f"      NXX={len(b['x'])}.")
        for k, (x, s) in enumerate(zip(b['x'], b['s'])):
            lines.append(assign(f'XCOOR({k + 1})', x))
            lines.append(assign(f'S({k + 1})', s))
        for name, value in [('BNOSE', b['bnose']), ('RLN', b['rln']),
                            ('SBD(4)', b['dn']), ('SBD(5)', b['d1']),
                            ('BD(81)', b['alpha0']), ('SBD(18)', b['cla']),
                            ('SBD(124)', b['cd0']),
                            ('SBD(120)', b['sbd120']),
                            ('BD(1)', b['x'][-1])]:
            lines.append(assign(name, value))
        for j in range(na):
            lines.append(assign(f'BODY({21 + j})', b['cl'][j]))
            lines.append(assign(f'BODY({1 + j})', b['cd'][j]))
            lines.append(assign(f'BODY({41 + j})', b['cm'][j]))
        for key, win, a, blk, slg in [('wing', 'WIN', 'A', 'WING', 'SLG'),
                                      ('tail', 'HIN', 'AHT', 'HT', 'STG')]:
            s = c[key]
            for index, value in s['in'].items():
                if index == 15:
                    lines.append(f"      {win}(15)=WT({int(value)})")
                else:
                    lines.append(assign(f'{win}({index})', value))
            for index, value in s['a'].items():
                lines.append(assign(f'{a}({index})', value))
            ae = s['aero']
            lines.append(assign(f'{blk}(101)', ae['cla']))
            lines.append(assign(f'{slg}(80)', ae['cd0']))
            lines.append(assign(f'{slg}(134)', ae['xac']))
            for j in range(na):
                lines.append(assign(f'{blk}({21 + j})', ae['cl'][j]))
                lines.append(assign(f'{blk}({61 + j})', ae['cn'][j]))
                lines.append(assign(f'{blk}({81 + j})', ae['ca'][j]))
                lines.append(assign(f'{slg}({53 + j})', ae['cdl'][j]))
        st = c['stale']
        for name, key in [('SWB(2)', 'kkwb'), ('SWB(37)', 'kkbw'),
                          ('SHB(2)', 'hkkwb'), ('SHB(37)', 'hkkbw'),
                          ('BW(41)', 'bw41'), ('BH(41)', 'bh41')]:
            lines.append(assign(name, st[key]))
        lines.append('      CALL M20O24')
        for tag, text in [('BW', '(BW(J),J=1,140)'), ('BH', '(BH(J),J=1,140)'),
                          ('BWV', '(BWV(J),J=1,140)'),
                          ('SWB', '(SWB(J),J=1,61)'), ('SHB', '(SHB(J),J=1,61)'),
                          ('X', 'SLG(134),STG(134),BD(66),BD(89),BD(83)'),
                          ('AB', '(BD(254+J),J=1,NALPHA)'),
                          ('CM0', 'TRA(73),TRAH(73)')]:
            lines.append(f"      WRITE(6,'(A,140ES25.16)') '{tag}',")
            lines.append(f"     1{text}")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('supwb', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = []
    for c, r in zip(all_cases, records):
        c = copy.deepcopy(c)
        for key in ('wing', 'tail'):
            c[key]['in'] = {str(k): v for k, v in c[key]['in'].items()}
            c[key]['a'] = {str(k): v for k, v in c[key]['a'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('supwb', payload))


if __name__ == '__main__':
    main()
