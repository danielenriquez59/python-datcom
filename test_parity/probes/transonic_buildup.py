"""
Probe the transonic buildup and save the fixture: TRANCM and TRHTCM (with
WBTRAN, HBTRAN and WBCM1), TRACM0, TRANCD, and WBTRA (with TRAWBT), called
in the order overlays 24, 25 and 35 run them.

Run from the repository root: ``python test_parity/probes/transonic_buildup.py``.
"""

import copy
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['trancm', 'trhtcm', 'wbtran', 'hbtran', 'wbcm1', 'tracm0',
            'wbcm0', 'tablec', 'trancd', 'wbcdl', 'tables', 'tbsub',
            'tbtrn', 'tbsup', 'wbtra', 'trawbt', 'tranac', 'tlin4x',
            'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'tbfunx', 'quad']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,KK
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /FLGTCD/ FLC(160)
      COMMON /SYNTSS/ XCG, XW, ZW, ALIW, ZCG, XH, ZH, ALIH, XV,
     1                VERTUP, HINAX, XVF, SCALE, ZV, ZVF, YV, YF,
     2                PHIV, PHIF
      COMMON /WINGI/  WINGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195), B(49)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108), STBH(135)
      COMMON /SUPWBB/ SWB(61), SWBH(61)
      COMMON /WHWB/   FACT(182), WB(39), HB(39)
      COMMON /SUPBOD/ SBD(229)
      COMMON /BODYI/  NXX, XCOOR(20)
      COMMON /IBODY/  PB, BODY(400)
      COMMON /IWING/  PW, WING(400)
      COMMON /IBW/    PBW, BW(380)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IBH/    PBH, BH(380)
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /IVT/    PVT, VT(380)
      COMMON /IVF/    PVF, VF(380)
      COMMON /IBV/    PBV, BV(380)
      COMMON /IBWH/   PBWH, BWH(380)
      COMMON /IBWV/   PBWV, BWV(380)
      COMMON /IBWHV/  PBWHV, BWHV(380)
      COMMON /EXPER/  KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1                KBODY, KWING, KHT, KVT, KWB,
     2                KDWASH(3), ALPO, ALPL
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB,KDWASH
      REAL NXX
      DIMENSION WT(4)
"""

UNUSED = 1.0e-30
ALPHA = [-2.0, 0.0, 2.0, 4.0, 8.0, 12.0]


def surface(sspn=15.0, sspne=13.0, chrdr=7.0, ar=4.0, taper=0.4,
            tanle=0.7, tovc=0.06, kind=1.0, twist=-2.0):
    """One surface: its input words and its ``A`` block."""
    c2 = math.atan(tanle) - 0.12
    return {
        'in': {3: sspne, 4: sspn, 6: chrdr, 11: twist, 12: sspn * 0.4,
               13: 3.0, 14: -1.0, 15: kind, 16: tovc, 62: 0.008,
               93: 0.01, 94: 0.2},
        'a': {3: ar * sspn * 0.9, 7: ar, 10: chrdr * 1.05, 11: 0.03,
              12: 2.5, 16: chrdr * 0.7, 24: 22.0, 27: taper,
              38: tanle, 40: math.degrees(math.atan(tanle)) - 3.0,
              44: tanle * 0.9, 56: -0.1, 62: tanle, 73: math.cos(c2),
              80: 0.2, 118: taper, 120: ar, 122: chrdr * 0.7,
              126: -2.0, 127: 14.0, 161: 1.2, 173: 3.5},
    }


def case(mach=0.95, wing=None, tail=None, bo=True, wgpl=True, htpl=True,
         mfb=0.9, xw=16.0, xh=38.0, length=44.0, kdwash=(False,) * 3):
    return {
        'mach': mach, 'mfb': mfb, 'alpha': ALPHA, 'bo': bo, 'wgpl': wgpl,
        'htpl': htpl, 'kdwash': list(kdwash),
        'wing': wing or surface(),
        'tail': tail or surface(sspn=6.0, sspne=5.3, chrdr=3.5, ar=3.5,
                                taper=0.5, tanle=0.5),
        'sref': 160.0, 'cbarr': 5.5, 'blref': 30.0,
        'xcg': 21.0, 'xw': xw, 'zw': -1.2, 'aliw': 2.0, 'xh': xh,
        'reynolds': 2.0e6, 'tr': 0.4,
        'body': {'length': length, 'x_max_area': 20.0, 'alpha0': UNUSED,
                 'max_diameter': 4.2, 'bd87': 4.0, 'sbd120': 4.3,
                 'cla': 0.0042, 'cma': 0.021},
        # Words the earlier transonic passes leave: WING/HT(101) and (1),
        # BW/BH(1) from TRSONI/TRSONJ, B(48) from the Mach-zero pass, the
        # buildup BW(101) used without a body, stale BW/BH(41), WB/HB(13).
        'cla_wing': 0.075, 'cd0_wing': 0.009, 'cla_tail': 0.05,
        'cd0_tail': 0.006, 'cd0_wb': 0.021, 'cd0_hb': 0.011, 'b48': 0.07,
        'bw101': 0.083, 'stale': {'bw41': -0.01, 'bh41': -0.002,
                                  'wb13': 0.31, 'wb14': 0.29,
                                  'hb13': 0.21, 'hb14': 0.19},
        'dwash': {'q_ratio': 0.93, 'deda': 0.41},
    }


def with_surface(base, key, **kw):
    c = copy.deepcopy(base)
    c[key] = surface(**kw)
    return c


def cases():
    tail = dict(sspn=6.0, sspne=5.3, chrdr=3.5, ar=3.5, taper=0.5, tanle=0.5)
    base = case()
    out = [
        base,
        case(mach=0.85),
        with_surface(case(mach=1.05), 'wing', tovc=0.09),   # refaired
        case(mach=1.0),                                     # BETA = 1e-7
        with_surface(case(mach=1.2), 'wing', taper=0.0),    # untapered
        with_surface(case(mach=0.95), 'wing', ar=2.0),      # Fig. 10 K_B(W)
        case(mach=1.2, xw=38.0, length=44.0),               # TE past body
        case(bo=False),
        case(wgpl=False, bo=False),
        with_surface(case(), 'tail', kind=3.0, **{k: v for k, v in
                                                  tail.items()}),
        with_surface(case(), 'wing', sspn=8.0, sspne=6.0),  # b_w < 1.5 b_h
        case(kdwash=(True, False, True)),
        with_surface(case(mach=1.1), 'wing', ar=7.0),       # out of range
        with_surface(case(mach=1.3), 'wing', tanle=0.0),    # unswept
        with_surface(case(mach=0.9), 'wing', tovc=0.11, ar=2.5, taper=0.2,
                     tanle=1.2),
    ]
    out[10]['tail'] = surface(**tail)
    # A large dihedral break, so the wake's effective span passes it.
    out[1]['wing']['in'][12] = 13.0
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1',
             '      NF=1']
    for n, c in enumerate(all_cases):
        label = 3000 + n
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.400) BODY(K)=0.')
        lines.append('         IF(K.LE.400) WING(K)=0.')
        for arr in ('BW', 'HT', 'BH', 'BWH', 'BWHV', 'BWV', 'BV', 'VT',
                    'VF'):
            lines.append(f'         IF(K.LE.380) {arr}(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.195) AHT(K)=0.')
        lines.append('         IF(K.LE.154) HTIN(K)=0.')
        lines.append('         IF(K.LE.101) WINGIN(K)=0.')
        lines.append('         IF(K.LE.108) TRA(K)=0.')
        lines.append('         IF(K.LE.108) TRAH(K)=0.')
        lines.append('         IF(K.LE.61) SWB(K)=0.')
        lines.append('         IF(K.LE.61) SWBH(K)=0.')
        lines.append('         IF(K.LE.60) DWASH(K)=0.')
        lines.append(f' {label} CONTINUE')
        for flag in ('bo', 'wgpl', 'htpl'):
            lines.append(f"      {flag.upper()}="
                         f"{'.TRUE.' if c[flag] else '.FALSE.'}")
        for j, k in enumerate(c['kdwash']):
            lines.append(f"      KDWASH({j + 1})="
                         f"{'.TRUE.' if k else '.FALSE.'}")
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['reynolds']))
        lines.append(assign('FLC(96)', c['tr']))
        for name, key in [('SREF', 'sref'), ('CBARR', 'cbarr'),
                          ('BLREF', 'blref'), ('XCG', 'xcg'), ('XW', 'xw'),
                          ('ZW', 'zw'), ('ALIW', 'aliw'), ('XH', 'xh'),
                          ('B(48)', 'b48'), ('WING(101)', 'cla_wing'),
                          ('WING(1)', 'cd0_wing'), ('HT(101)', 'cla_tail'),
                          ('HT(1)', 'cd0_tail'), ('BW(1)', 'cd0_wb'),
                          ('BH(1)', 'cd0_hb'), ('BW(101)', 'bw101')]:
            lines.append(assign(name, c[key]))
        for win, a, tra, key in [('WINGIN', 'A', 'TRA', 'wing'),
                                 ('HTIN', 'AHT', 'TRAH', 'tail')]:
            s = c[key]
            for index, value in s['in'].items():
                if index == 15:
                    lines.append(f"      {win}(15)=WT({int(value)})")
                else:
                    lines.append(assign(f'{win}({index})', value))
            for index, value in s['a'].items():
                lines.append(assign(f'{a}({index})', value))
            lines.append(assign(f'{tra}(4)', c['mach']))
            lines.append(assign(f'{tra}(6)', c['mfb']))
        st = c['stale']
        for name, key in [('BW(41)', 'bw41'), ('BH(41)', 'bh41'),
                          ('WB(13)', 'wb13'), ('WB(14)', 'wb14'),
                          ('HB(13)', 'hb13'), ('HB(14)', 'hb14')]:
            lines.append(assign(name, st[key]))
        b = c['body']
        lines.append('      NXX=3.')
        lines.append(assign('XCOOR(3)', b['length']))
        lines.append(assign('BD(1)', b['length']))
        for name, key in [('BD(2)', 'x_max_area'), ('BD(81)', 'alpha0'),
                          ('BD(85)', 'max_diameter'), ('BD(87)', 'bd87'),
                          ('SBD(120)', 'sbd120'), ('BODY(101)', 'cla'),
                          ('BODY(121)', 'cma')]:
            lines.append(assign(name, b[key]))
        lines.append(assign('DWASH(1)', c['dwash']['q_ratio']))
        lines.append(assign('DWASH(41)', c['dwash']['deda']))
        # Overlay 25 (M25O31), then TRANCD (end of overlay 24), then 35.
        lines.append('      IF(WGPL) CALL TRANCM')
        lines.append('      IF(HTPL) CALL TRHTCM')
        lines.append('      IF(BO) CALL TRACM0')
        lines.append('      CALL TRANCD')
        lines.append('      CALL WBTRA')
        for tag, text in [
                ('WTRA', '(TRA(J),J=71,107)'), ('HTRA', '(TRAH(J),J=71,108)'),
                ('WSWB', '(SWB(J),J=1,61)'), ('HSWB', '(SWBH(J),J=1,61)'),
                ('W', 'WING(121),BW(101),BW(121),WB(13),WB(14),WB(15),'
                      'A(173),BW(41)'),
                ('H', 'HT(121),BH(101),BH(121),HB(13),HB(14),HB(15),'
                      'AHT(173),BH(41)'),
                ('AB', '(BD(254+J),J=1,NALPHA)'),
                ('BWCD', '(BW(J),J=1,NALPHA)'), ('BHCD', '(BH(J),J=1,NALPHA)'),
                ('WBT', 'DWASH(1),DWASH(41),BWHV(101),BWHV(121),BWH(101),'
                        'BWH(121),TRA(2),TRA(11),TRA(26),TRA(69)')]:
            lines.append(f"      WRITE(6,'(A,70ES25.16)') '{tag}',")
            parts, line = [], ''
            for item in text.split('),'):
                item = item if item.endswith(')') else item + ')'
                if len(line) + len(item) > 60:
                    parts.append(line)
                    line = ''
                line += (',' if line else '') + item
            parts.append(line)
            lines += [f"     1{part}{',' if k < len(parts) - 1 else ''}"
                      for k, part in enumerate(parts)]
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('transonic_buildup', driver(all_cases),
                                ROUTINES, EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r}
               for c, r in zip(all_cases, records)]
    print(save('transonic_buildup', payload))


if __name__ == '__main__':
    main()
