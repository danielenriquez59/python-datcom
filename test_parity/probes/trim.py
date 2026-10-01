"""
Probe DRAGFP, TRIMRT and TRIMR2 and save the fixture.  Each case runs
DRAGFP then TRIMRT, as M38O46 does, and TRIMR2 on its own blocks.

Run from the repository root: ``python test_parity/probes/trim.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['dragfp', 'trimrt', 'trimr2', 'tbfunx', 'quad', 'tlinex',
            'tlin1x', 'glook', 'switch']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IZ,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF
      COMMON /IDWASH/ PDWASH, DWASH(60)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IBODY/  PBODY, BODY(400)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IVT/    PVT, VT(380)
      COMMON /IBW/    PBW, BW(380)
      COMMON /IBWH/   PBWH, BWH(380)
      COMMON /IBWHV/  PBWHV, BWHV(380)
      COMMON /SUPWH/  FCM(282)
      COMMON /SUPDW/  DWA(35),TCD(58)
      COMMON /FLAPIN/ F(69)
      COMMON /POWR/   PW(104),FLP(189),TRM(22)
      COMMON /WINGD/  A(195),B(49)
      COMMON /HTDATA/ AHT(195),BHT(49)
      COMMON /FLGTCD/ FLC(93)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /BDATA/  BD(762)
      COMMON /HTI/    HTIN(154)
      COMMON /WBHCAL/ WBT(155)
      COMMON /WHAERO/ XXX(157), DHT(55)
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
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0]


def case(n, ftype=1.0, tail=False, delta=(-10.0, 0.0, 10.0, 20.0, 30.0),
         wgpl=True, htpl=False, bo=False, vtpl=False, aclmax=18.0,
         dcm_sign=-1.0, cm_scale=1.0, alih=-1.0, stop_cm=None,
         tail_aclmax=14.0):
    rng = np.random.default_rng(700 + n)
    nd, na = len(delta), len(ALPHA)
    alpha = np.array(ALPHA)
    eps = 0.3 * alpha + 0.5
    blocks = {}
    for name, cl, cm, cd in [
            ('WING', 0.07, -0.004, 0.01), ('BW', 0.08, -0.006, 0.015),
            ('BWH', 0.075, -0.012, 0.02)]:
        blocks[name] = {'cl': list(cl * (alpha + 2.0)),
                        'cm': list(cm_scale * cm * (alpha + 1.0)),
                        'cd': list(cd + 0.0006 * alpha**2)}
    blocks['BWHV'] = {'cd': list(0.024 + 0.0006 * alpha**2)}
    if stop_cm is not None:
        blocks['BW']['cm'][stop_cm] = 2.0 * UNUSED
    return {
        'alpha': ALPHA, 'epsilon': list(eps), 'q_ratio': list(
            0.95 - 0.004 * alpha),
        'flags': {'wgpl': wgpl, 'htpl': htpl, 'bo': bo, 'vtpl': vtpl},
        'flap': {'delta': list(delta), 'ftype': ftype,
                 'gd1': list(rng.uniform(0.2, 1.2, 12)),
                 'gd2': list(rng.uniform(0.0, 0.3, 12)),
                 'gd3': list(rng.uniform(0.2, 0.9, 12)),
                 'adave': list(rng.uniform(0.3, 0.7, nd)),
                 'eta': [0.1, 0.2, 0.3, 0.4, 0.75],
                 'rkb': list(rng.uniform(0.1, 0.3, 4)),
                 'cfoca': 0.25,
                 'dcl': list(0.02 * np.array(delta)),
                 'dcm': list(dcm_sign * 0.004 * np.array(delta) +
                             rng.uniform(0, 1e-4, nd)),
                 'dclmax': list(0.01 * np.array(delta)),
                 'dcdmin': list(0.0004 * np.abs(delta)),
                 'chd': list(-0.003 * np.array(delta))},
        'wing': {'a3': 300.0, 'a7': 5.0, 'b49': -1.5, 'b43': aclmax,
                 'b23': list(alpha + 0.7)},
        'tail': {'a3': 60.0 if tail else UNUSED,
                 'a7': 3.8 if tail else UNUSED,
                 'alpha': list(alpha + 1.0) if tail else [UNUSED] * na,
                 'bht43': tail_aclmax},
        'blocks': blocks,
        # TRIMR2's tail and carryover words.
        't2': {'alpha': list(alpha + 1.0), 'cl': list(0.05 * (alpha + 1.0)),
               'cd': list(0.008 + 0.0005 * alpha**2), 'cla': 0.05,
               'cd0': 0.007, 'cm0': -0.01, 'ar': 3.8, 'area': 60.0,
               'sspn': 7.0, 'sspne': 6.1, 'e': 0.85},
        'wbt': {'kwb': 1.08, 'kbw': 0.14, 'kkwb': 0.95, 'kkbw': 0.12,
                'clbh': list(0.002 * alpha), 'cdov': 0.004,
                'f46': list(rng.uniform(0.0, 0.2, na)),
                'f68': list(rng.uniform(0.0, 0.5, na))},
        'position': {'xba': -22.0, 'zba': 1.5, 'xcg': 20.0, 'hinax': 41.0,
                     'alih': alih},
        'sref': 300.0, 'cbarr': 7.0,
    }


def cases():
    return [
        case(0),
        case(1, ftype=3.0, tail=True, htpl=True, wgpl=True),
        case(2, ftype=2.0, bo=True, wgpl=False, dcm_sign=1.0),
        case(3, ftype=5.0, aclmax=6.0),                  # stall stop
        case(4, ftype=6.0, cm_scale=40.0),               # out of range
        case(5, ftype=1.0, htpl=True, vtpl=True, tail=True, alih=12.0),
        case(6, ftype=4.0, bo=True, stop_cm=3, tail_aclmax=30.0),
        case(7, ftype=1.0, tail_aclmax=9.0, alih=-3.0),
        case(8, ftype=2.0, aclmax=30.0, cm_scale=0.3),   # trims throughout
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      IZ=1']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,400')
        for arr in ('WING', 'BODY'):
            lines.append(f'         {arr}(K)=0.')
        for arr in ('HT', 'VT', 'BW', 'BWH', 'BWHV'):
            lines.append(f'         IF(K.LE.380) {arr}(K)=0.')
        lines.append('         IF(K.LE.22) TRM(K)=-9.')
        lines.append(f' {label} CONTINUE')
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for flag, value in c['flags'].items():
            lines.append(f"      {flag.upper()}=.{'TRUE' if value else 'FALSE'}.")
        f, w, t = c['flap'], c['wing'], c['tail']
        lines.append(assign('F(16)', len(f['delta'])))
        lines.append(assign('F(17)', f['ftype']))
        for k, d in enumerate(f['delta']):
            lines.append(assign(f'F({k + 1})', d))
            lines.append(assign(f'FCM({63 + k})', f['adave'][k]))
            for name, start in [('dcl', 201), ('dcm', 211), ('dclmax', 221),
                                ('dcdmin', 231), ('chd', 261)]:
                lines.append(assign(f'WING({start + k})', f[name][k]))
        for k in range(12):
            lines.append(assign(f'TCD({29 + k})', f['gd1'][k]))
            lines.append(assign(f'TCD({1 + k})', f['gd2'][k]))
            lines.append(assign(f'TCD({15 + k})', f['gd3'][k]))
        for k in range(5):
            lines.append(assign(f'FLP({1 + k})', f['eta'][k]))
        for k in range(4):
            lines.append(assign(f'FLP({20 + k})', f['rkb'][k]))
        lines.append(assign('FLP(61)', f['cfoca']))
        lines.append(assign('A(3)', w['a3']))
        lines.append(assign('A(7)', w['a7']))
        lines.append(assign('B(49)', w['b49']))
        lines.append(assign('B(43)', w['b43']))
        lines.append(assign('AHT(3)', t['a3']))
        lines.append(assign('AHT(7)', t['a7']))
        lines.append(assign('BHT(43)', t['bht43']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        for j in range(na):
            lines.append(assign(f'FLC({23 + j})', c['alpha'][j]))
            lines.append(assign(f'B({23 + j})', w['b23'][j]))
            lines.append(assign(f'BHT({23 + j})', t['alpha'][j]))
            lines.append(assign(f'DWASH({21 + j})', c['epsilon'][j]))
            lines.append(assign(f'DWASH({1 + j})', c['q_ratio'][j]))
            for name, b in c['blocks'].items():
                for key, start in (('cd', 0), ('cl', 20), ('cm', 40)):
                    if key in b:
                        lines.append(assign(f'{name}({start + 1 + j})',
                                            b[key][j]))
        nd = len(f['delta'])
        lines.append('      CALL DRAGFP')
        lines.append("      WRITE(6,'(A,200ES25.16)') 'DCDI',"
                     f"(BODY(200+J),J=1,{na * nd})")
        lines.append("      WRITE(6,'(A,10ES25.16)') 'DELCDM',"
                     f"(WING(230+J),J=1,{nd})")
        lines.append("      WRITE(6,'(A,10ES25.16)') 'DELCDF',"
                     f"(TCD(48+J),J=1,{nd})")
        lines.append("      WRITE(6,'(A,ES25.16)') 'KPRM',TCD(47)")
        lines.append('      CALL TRIMRT')
        for tag, arr, start in [('DELTAT', 'VT', 200), ('DCLT', 'VT', 220),
                                ('CLMAXT', 'VT', 240), ('CDMINT', 'VT', 280),
                                ('CHDT', 'VT', 320), ('CDIT', 'VT', 260),
                                ('UTCL', 'HT', 200), ('UTCM', 'HT', 240),
                                ('UTCD', 'HT', 220)]:
            lines.append(f"      WRITE(6,'(A,20ES25.16)') '{tag}',"
                         f"({arr}({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,22ES25.16)') 'TRM',(TRM(J),J=1,22)")
        # TRIMR2 on its own words.
        t2, wbt, p = c['t2'], c['wbt'], c['position']
        lines.append(f'      DO {label + 500} K=1,380')
        lines.append('         HT(K)=0.')
        lines.append('         VT(K)=0.')
        lines.append('         IF(K.LE.22) TRM(K)=-9.')
        lines.append(f' {label + 500} CONTINUE')
        for j in range(na):
            lines.append(assign(f'BHT({23 + j})', t2['alpha'][j]))
            lines.append(assign(f'HT({21 + j})', t2['cl'][j]))
            lines.append(assign(f'HT({1 + j})', t2['cd'][j]))
            lines.append(assign(f'WBT({130 + j})', wbt['clbh'][j]))
            lines.append(assign(f'WBT({46 + j})', wbt['f46'][j]))
            lines.append(assign(f'WBT({68 + j})', wbt['f68'][j]))
        for name, value in [('HT(101)', t2['cla']), ('BHT(46)', t2['cd0']),
                            ('BHT(47)', t2['cm0']), ('AHT(7)', t2['ar']),
                            ('AHT(3)', t2['area']), ('HTIN(4)', t2['sspn']),
                            ('HTIN(3)', t2['sspne']), ('DHT(30)', t2['e']),
                            ('WBT(1)', wbt['kwb']), ('WBT(2)', wbt['kbw']),
                            ('WBT(150)', wbt['kkwb']),
                            ('WBT(151)', wbt['kkbw']),
                            ('WBT(66)', wbt['cdov']), ('BD(63)', p['xba']),
                            ('BD(64)', p['zba']), ('SYNA(1)', p['xcg']),
                            ('SYNA(11)', p['hinax']),
                            ('SYNA(8)', p['alih'])]:
            lines.append(assign(name, value))
        lines.append('      CALL TRIMR2')
        for tag, arr, start in [('ALIHT', 'HT', 220), ('CDHTRM', 'HT', 240),
                                ('CMHTRM', 'HT', 280), ('HMTRM', 'HT', 300),
                                ('HMUNT', 'HT', 200), ('HT261', 'HT', 260),
                                ('CLWBT', 'VT', 220), ('CDWBT', 'VT', 200)]:
            lines.append(f"      WRITE(6,'(A,20ES25.16)') '{tag}',"
                         f"({arr}({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,22ES25.16)') 'TRM2',(TRM(J),J=1,22)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('trim', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('trim', payload))


if __name__ == '__main__':
    main()
