"""
Probe WBAERO's wing pass (WBDRAG, WBLIFT, WBCM, WBCM0, BODOWG) and save the
fixture.

Run from the repository root: ``python test_parity/probes/wbaero.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['wbaero', 'wbdrag', 'wblift', 'wbcm', 'wbcm0', 'tablec',
            'bodowg', 'getmax', 'ali', 'tbfunx', 'quad', 'tlinex', 'tlin1x',
            'glook', 'switch']

# EXSUBT reads experimental namelists from a file and EXIT closes files;
# neither has numerical content.
EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
      SUBROUTINE EXIT
      RETURN
      END
"""

# WBAERO's COMMON declarations, copied verbatim, and BODOWG's /BODYI/.
COMMONS = """\
      COMMON /IBODY/   PBODY,  BODY(400)
      COMMON /IWING/   PWING,  WING(400)
      COMMON /IHT/     PHT,    HT(380)
      COMMON /IBW/     PBW,    BW(380)
      COMMON /IBH/     PBH,    BH(380)
      COMMON /IBV/     PBV,    BV(380)
      COMMON /FLGTCD/ FLC(160)
      COMMON /OPTION/ SREF, CBARR, ROUGFC, BLREF
      COMMON /SYNTSS/ XCG, XW, ZW, ALIW, ZCG, XH, ZH, ALIH, XV,
     1                VERTUP, HINAX, XVF, SCALE, ZV, ZVF, YV, YF,
     2                PHIV, PHIF
      COMMON /WINGI/  WINGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /WINGD/  A(195), B(49)
      COMMON /BDATA/  BD(762)
      COMMON /WHWB/   FACT(182), WB(39), HB(39)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /WHAERO/ C(51), D(55), CHT(51), DHT(55), DVT(55), DVF(55)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,L
      COMMON /EXPER/  KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1                KBODY, KWING, KHT, KVT, KWB, KDWASH(3),
     2                ALPOW, ALPLW, ALPOH, ALPLH
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      COMMON /BODYI/ XNX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1               BNOSE,BTAIL,BLN,BLA,DS
      LOGICAL  VERTUP,MDATA,KBODY,KWING,KHT,KVT,KWB,KDWASH
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
"""

UNUSED = 1.0e-30


def case(rng, alpha=None, sspn=15.0, sspne=13.0, aliw=0.0, beta_a=5.0,
         temp0_small=False, regression=True, missing_cm=None,
         zero_cn=None, experimental=False, stale=None, stall=18.0):
    alpha = np.array(alpha if alpha is not None else
                     [-4., 0., 4., 8., 12., 16., 20., 24.])
    n = len(alpha)
    local = alpha + aliw
    alpha_zero = -1.5
    body_alpha_zero = -0.5
    mach = 0.5
    beta = np.sqrt(1 - mach**2)
    a7 = beta_a / beta
    # A wing lift curve with a stall, so CL is not monotonic past it.
    d = local - alpha_zero
    cl_wing = 0.07 * d * np.where(local > stall, np.exp(-(local - stall) / 6.0), 1.0)
    cn_wing = cl_wing / np.cos(np.radians(local))
    if zero_cn is not None:
        cn_wing[zero_cn] = 0.0
    cm_wing = -0.03 - 0.02 * cl_wing
    if missing_cm is not None:
        cm_wing[missing_cm:] = 2 * UNUSED
    tan_le = 0.02 if temp0_small else 0.5
    x = np.linspace(0.0, 40.0, 12)
    radius = 2.0 * np.sin(np.pi * np.clip(x / 44.0, 0, 1))**0.5
    s = np.pi * radius**2
    c = {
        'alpha_deg': list(alpha),
        'surface': {
            'sspn': sspn, 'sspne': sspne, 'chrdr': 6.0,
            'twista': rng.uniform(-3.0, 0.0), 'tovc': 0.06, 'ler': 0.008,
            'ycm': 0.01, 'cld': 0.2,
            'a7': a7, 'a10': 5.5, 'a27': 0.45, 'a38': tan_le,
            'a44': tan_le * 0.8, 'a62': tan_le, 'a80': -0.1,
            'a118': 0.5, 'a120': a7 * 1.05 if regression else 7.5,
            'a122': 4.9, 'a129': 2.0e6, 'a160': rng.uniform(0.5, 11.0),
            'a173': rng.uniform(1.0, 4.0),
            'beta': beta, 'cm0': -0.03, 'alpha_zero_lift': alpha_zero,
            'clmax': 1.2, 'alpha_clmax': stall,
            'local_alpha': list(local), 'x_quarter_chord': 19.0,
            'cd0': 0.007, 'cdl': list(0.05 * cl_wing**2), 'c6': 0.27,
        },
        'surface_alone': {
            'cd': list(0.007 + 0.05 * cl_wing**2), 'cl': list(cl_wing),
            'cm': list(cm_wing), 'cn': list(cn_wing),
            'ca': list(0.006 - 0.01 * cl_wing), 'cla': 0.07,
            'cma': -0.0014,
        },
        'body': {
            'cd': list(0.01 + 0.001 * alpha**2 / 10),
            'cl': list(0.004 * (alpha - body_alpha_zero)),
            'cm': list(0.006 * (alpha - body_alpha_zero)),
            'cla': list(0.004 + 0.0001 * np.abs(alpha)),
            'cma': list(0.006 + 0.0002 * np.abs(alpha)),
            'cm0': 0.002, 'alpha_zero_lift': body_alpha_zero,
            'cd_friction': 0.009, 'cd_base': 0.0015,
            'cdl': list(0.0005 * alpha**2 / 10),
            'x': list(x), 's': list(s),
        },
        'synthesis': {'xcg': 21.0, 'xw': 16.0, 'zw': -0.8, 'zcg': 0.0,
                      'aliw': aliw},
        'flight': {'mach': mach, 'reynolds_per_length': 2.0e6, 'tr': 0.4},
        'cbarr': 5.2,
        'experimental': ({'kbody': True, 'kwing': False} if experimental
                         else None),
        'stale': stale,
    }
    if experimental:
        # One angle without supplied body data keeps the computed drag.
        c['body']['cd'][2] = UNUSED
    return c


def cases():
    rng = np.random.default_rng(4312)
    return [
        case(rng),
        case(rng, aliw=2.0),
        case(rng, sspne=14.2),                      # d/b below 0.3
        case(rng, sspn=15.0, sspne=4.0, aliw=1.0,   # d/b between 0.3 and 0.8
             regression=False),
        case(rng, sspne=1.5, aliw=1.5,              # d/b above 0.8
             stale={'kwb': 1.3, 'kbw': 0.4, 'wb4': 0.09, 'wb5': 0.03,
                    'kwb_incidence': 0.9, 'kbw_incidence': 0.3}),
        case(rng, beta_a=2.5),                      # the ellipse fit
        case(rng, beta_a=3.2, temp0_small=True),    # TEMP0 below one
        case(rng, regression=False),
        case(rng, missing_cm=5),
        case(rng, zero_cn=3),
        case(rng, experimental=True),
        case(rng, alpha=[2., 4., 6., 10.]),         # all positive lift
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             '      WGPL=.TRUE.', '      HTPL=.FALSE.', '      BO=.FALSE.',
             '      VTPL=.FALSE.', '      I=1', '      NF=-1']
    for n, c in enumerate(all_cases):
        sf, sa, bd, sy, fl = (c['surface'], c['surface_alone'], c['body'],
                              c['synthesis'], c['flight'])
        alpha = c['alpha_deg']
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,762')
        lines.append('         BD(K)=0.')
        lines.append('         IF(K.LE.400) BODY(K)=0.')
        lines.append('         IF(K.LE.400) WING(K)=0.')
        lines.append('         IF(K.LE.380) BW(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.182) FACT(K)=0.')
        lines.append('         IF(K.LE.160) FLC(K)=0.')
        lines.append('         IF(K.LE.101) WINGIN(K)=0.')
        lines.append('         IF(K.LE.55) D(K)=0.')
        lines.append('         IF(K.LE.51) C(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append('         IF(K.LE.39) WB(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f'      NALPHA={len(alpha)}')
        for j, a in enumerate(alpha):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'B({23 + j})', sf['local_alpha'][j]))
            lines.append(assign(f'BD({255 + j})',
                                a + bd['alpha_zero_lift']))
            lines.append(assign(f'D({36 + j})', sf['cdl'][j]))
            lines.append(assign(f'BD({215 + j})', bd['cdl'][j]))
            for base, key in [(0, 'cd'), (20, 'cl'), (40, 'cm'),
                              (100, 'cla'), (120, 'cma')]:
                lines.append(assign(f'BODY({base + 1 + j})', bd[key][j]))
            for base, key in [(0, 'cd'), (20, 'cl'), (40, 'cm'), (60, 'cn'),
                              (80, 'ca')]:
                lines.append(assign(f'WING({base + 1 + j})', sa[key][j]))
        lines.append(assign('WING(101)', sa['cla']))
        lines.append(assign('WING(121)', sa['cma']))
        lines.append(assign('FLC(3)', fl['mach']))
        lines.append(assign('FLC(43)', fl['reynolds_per_length']))
        lines.append(assign('FLC(96)', fl['tr']))
        for key, index in [('a7', 7), ('a10', 10), ('a27', 27), ('a38', 38),
                           ('a44', 44), ('a62', 62), ('a80', 80),
                           ('a118', 118), ('a120', 120), ('a122', 122),
                           ('a129', 129), ('a160', 160), ('a173', 173)]:
            lines.append(assign(f'A({index})', sf[key]))
        lines.append(assign('B(1)', fl['mach']))
        for key, index in [('beta', 2), ('alpha_clmax', 43), ('clmax', 44),
                           ('cm0', 47), ('alpha_zero_lift', 49)]:
            lines.append(assign(f'B({index})', sf[key]))
        for key, index in [('sspne', 3), ('sspn', 4), ('chrdr', 6),
                           ('twista', 11), ('tovc', 16), ('ler', 62),
                           ('ycm', 93), ('cld', 94)]:
            lines.append(assign(f'WINGIN({index})', sf[key]))
        lines.append(assign('D(20)', sf['cd0']))
        lines.append(assign('C(6)', sf['c6']))
        lines.append(assign('BD(59)', bd['cd_friction']))
        lines.append(assign('BD(60)', bd['cd_base']))
        lines.append(assign('BD(62)', bd['cm0']))
        lines.append(assign('BD(77)', sy['aliw']))
        lines.append('      BD(79)=COS(BD(77)/RAD)')
        lines.append(assign('BD(81)', bd['alpha_zero_lift']))
        lines.append(assign('BD(83)', sf['x_quarter_chord']))
        lines.append(assign('BD(87)', 2.0 * (sf['sspn'] - sf['sspne'])))
        lines.append(f"      XNX={len(bd['x'])}.")
        for k, (xv, sv) in enumerate(zip(bd['x'], bd['s'])):
            lines.append(assign(f'X({k + 1})', xv))
            lines.append(assign(f'S({k + 1})', sv))
        lines.append(assign('BD(1)', bd['x'][-1]))
        for key, name in [('xcg', 'XCG'), ('xw', 'XW'), ('zw', 'ZW'),
                          ('zcg', 'ZCG'), ('aliw', 'ALIW')]:
            lines.append(assign(name, sy[key]))
        lines.append(assign('CBARR', c['cbarr']))
        experimental = c['experimental'] or {}
        lines.append(f"      KBODY=.{'TRUE' if experimental.get('kbody') else 'FALSE'}.")
        lines.append(f"      KWING=.{'TRUE' if experimental.get('kwing') else 'FALSE'}.")
        for key, index in [('kwb', 2), ('kbw', 3), ('wb4', 4), ('wb5', 5),
                           ('kwb_incidence', 7), ('kbw_incidence', 8)]:
            lines.append(assign(f'WB({index})', (c['stale'] or {}).get(key, 0.0)))
        lines.append('      CALL WBAERO')
        for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                           ('CA', 80), ('CLA', 100), ('CMA', 120)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"(BW({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,39ES25.16)') 'WB',(WB(J),J=1,39)")
        lines.append("      WRITE(6,'(A,41ES25.16)') 'FACT',(FACT(J),J=1,41)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('wbaero', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('wbaero', payload))


if __name__ == '__main__':
    main()
