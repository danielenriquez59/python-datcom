"""
Probe the transonic second-level routines and save the fixture: WBCLB
(called directly with its argument list) and SETUP2 (with CLBCLC) run
through its whole Mach schedule, the blocks it reads changed between
calls as the passes in between would change them.

Run from the repository root: ``python test_parity/probes/second_level.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['wbclb', 'bodowg', 'getmax', 'ali', 'setup2', 'clbclc', 'tbfunx',
            'quad']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /BDATA/  BD(762)
      COMMON /BODYI/  XNX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1                BNOSE,BTAIL,BLN,BLA,DS
      COMMON /IWING/  PWNG, WING(400)
      COMMON /IHT/    PHT,  HT(380)
      COMMON /IBW/    PBW,  BW(380)
      COMMON /IBH/    PBH,  BH(380)
      COMMON /FLGTCD/ FLC(96)
      COMMON /WINGI/  WINGIN(100)
      COMMON /HTI/    HTIN(154)
      COMMON /WINGD/  A(195), B(49)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108), STBH(135)
      COMMON /HTDATA/ AHT(195), BHT(49)
      COMMON /WBHCAL/ WBT(155)
      COMMON /LEVEL2/ SECOND(23)
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART
      REAL KWB,KBW,KKWB,KKBW,MACH,MFB
      DIMENSION CLW(20),CLB(20),CDW(20),CDB(20),ALP(20),ALPB(20),
     1          WIN(100),CLWB(20),CDWB(20),CLBB(20),FACT(41)
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0, 16.0]
BODY_X = [0.0, 6.0, 14.0, 26.0, 38.0, 44.0]
BODY_S = [0.0, 6.0, 12.5, 13.8, 11.0, 6.5]


def wbclb_case(sspn=6.0, missing_cl=(), missing_parts=(), unused_wb=(0, 1, 2,
               3, 4, 5), drag=False, clbb_set=(), alpb_shift=1.5):
    cl_w = [0.02 + 0.07 * a for a in ALPHA]
    cl_b = [0.004 * a for a in ALPHA]
    for j in missing_parts:
        cl_w[j] = UNUSED
    clwb = [0.05 + 0.09 * a for a in ALPHA]
    for j in unused_wb:
        clwb[j] = UNUSED if j not in missing_cl else -UNUSED
    cdwb = [0.02 + 0.001 * a * a for a in ALPHA]
    if not drag:
        cdwb[1] = UNUSED
    cd_w = [0.008 + 0.0009 * a * a for a in ALPHA]
    cd_b = [0.012 + 0.0002 * a * a for a in ALPHA]
    cd_b[3] = UNUSED
    clbb = [UNUSED] * len(ALPHA)
    for j in clbb_set:
        clbb[j] = 0.003 * j
    return {
        'alpha': ALPHA, 'alpha_body': [a + alpb_shift for a in ALPHA],
        'sspn': sspn, 'body_x': BODY_X, 'body_s': BODY_S,
        'xc': 18.0, 'taper': 0.4, 'aliw': 1.5,
        'cl_w': cl_w, 'cla_w': 0.072, 'cl_b': cl_b, 'cla_b': 0.004,
        'cd_w': cd_w, 'cd_b': cd_b, 'kwb': 1.12, 'kbw': 0.21,
        'cla_wb': 0.089, 'clbn14': 0.0012, 'cnam14': 0.05,
        'clblfb': 0.0019, 'clamfb': 0.095, 'mach': 1.05, 'mfb': 0.92,
        'clwb': clwb, 'cdwb': cdwb, 'clbb': clbb,
    }


def wbclb_cases():
    return [
        wbclb_case(),
        wbclb_case(sspn=15.0),                        # r/s < 1/3
        wbclb_case(missing_parts=(1, 4)),             # linear fallback
        wbclb_case(drag=True, clbb_set=(2, 5)),
        wbclb_case(unused_wb=(0, 3), missing_cl=(3,)),
        wbclb_case(alpb_shift=7.5),                   # vortex beyond 6 deg
    ]


def setup2_case(bo=True, wgpl=True, htpl=True, mfbw=0.92, mfbh=0.97):
    return {'bo': bo, 'wgpl': wgpl, 'htpl': htpl, 'mfbw': mfbw,
            'mfbh': mfbh, 'mach': 1.05, 'flc2': 6.0, 'nalpha': 6,
            'win68': 0.11, 'win69': 0.084, 'hin68': 0.09, 'hin69': 0.066}


def setup2_cases():
    return [setup2_case(), setup2_case(bo=False),
            setup2_case(htpl=False, mfbw=0.97),
            setup2_case(wgpl=False), setup2_case(bo=False, htpl=False)]


def block_values(step, name):
    """IOM block words the pass before step ``step`` leaves (CL at 21..,
    CLB at 181..), with a zero-lift first angle to exercise CLBCLC."""
    rng = np.random.default_rng(100 * step + len(name))
    cl = [0.0] + list(rng.uniform(-0.3, 0.8, 5))
    clb = [UNUSED] + list(rng.uniform(-0.01, 0.01, 5))
    if name == 'BH' and step == 5:
        cl = [0.0] * 6                               # no usable angle
    return cl, clb


def driver(wcases, scases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    n = 0
    lines.append(f'      XNX={len(BODY_X)}.')
    for k, (x, s) in enumerate(zip(BODY_X, BODY_S)):
        lines.append(assign(f'X({k + 1})', x))
        lines.append(assign(f'S({k + 1})', s))
    for c in wcases:
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        n += 1
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for k in range(na):
            for arr, key in [('ALP', 'alpha'), ('ALPB', 'alpha_body'),
                             ('CLW', 'cl_w'), ('CLB', 'cl_b'),
                             ('CDW', 'cd_w'), ('CDB', 'cd_b'),
                             ('CLWB', 'clwb'), ('CDWB', 'cdwb'),
                             ('CLBB', 'clbb')]:
                lines.append(assign(f'{arr}({k + 1})', c[key][k]))
        lines.append(assign('WIN(4)', c['sspn']))
        lines.append(assign('WIN(3)', c['sspn'] - 2.0))
        for name, key in [('XC', 'xc'), ('TAPR', 'taper'), ('ALIW', 'aliw'),
                          ('CLAW', 'cla_w'), ('CLAB', 'cla_b'),
                          ('KWB', 'kwb'), ('KBW', 'kbw'),
                          ('CLAWB', 'cla_wb'), ('CLBN14', 'clbn14'),
                          ('CNAM14', 'cnam14'), ('CLBLFB', 'clblfb'),
                          ('CLAMFB', 'clamfb'), ('MACH', 'mach'),
                          ('MFB', 'mfb')]:
            lines.append(assign(name, c[key]))
        lines += [
            '      CALL WBCLB(NALPHA,CLW,CLAW,CLB,CLAB,CDW,CDB,CLAWB,ALP,',
            '     1  ALPB,WIN,XC,TAPR,ALIW,KWB,KBW,CLBN14,CNAM14,CLBLFB,',
            '     2  CLAMFB,MACH,MFB,KKWB,KKBW,FACT,CLWB,CDWB,CLBB,CLBCL)',
            "      WRITE(6,'(A,20ES25.16)') 'CLWB',(CLWB(J),J=1,NALPHA)",
            "      WRITE(6,'(A,20ES25.16)') 'CDWB',(CDWB(J),J=1,NALPHA)",
            "      WRITE(6,'(A,20ES25.16)') 'CLBB',(CLBB(J),J=1,NALPHA)",
            "      WRITE(6,'(A,4ES25.16)') 'K',KKWB,KKBW,CLBCL,FACT(1)",
            "      WRITE(6,'(A,20ES25.16)') 'IV',(FACT(1+J),J=1,NALPHA)",
            "      WRITE(6,'(A,20ES25.16)') 'GO',(FACT(21+J),J=1,NALPHA)",
        ]
    for c in scases:
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        n += 1
        for flag in ('bo', 'wgpl', 'htpl'):
            lines.append(f"      {flag.upper()}="
                         f"{'.TRUE.' if c[flag] else '.FALSE.'}")
        lines += ['      SUBSON=.FALSE.', '      TRANSN=.TRUE.',
                  '      SUPERS=.FALSE.', f'      DO {8000 + n} K=1,23',
                  f' {8000 + n} SECOND(K)=-5.']
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign('FLC(2)', c['flc2']))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('TRA(6)', c['mfbw']))
        lines.append(assign('TRAH(6)', c['mfbh']))
        for name, key in [('WINGIN(68)', 'win68'), ('WINGIN(69)', 'win69'),
                          ('HTIN(68)', 'hin68'), ('HTIN(69)', 'hin69')]:
            lines.append(assign(name, c[key]))
        lines.append('      NF=-1')
        # Up to seven calls; each loads what the pass before would leave.
        for step in range(1, 8):
            lines.append(f"      IF(NF.LT.-7) GO TO {7000 + n}")
            lines.append(f"      NSTEP={step}")
            lines.append(assign('WBT(67)', 0.01 * step))
            lines.append(assign('WBT(155)', -0.02 * step))
            lines.append(assign('BW(101)', 0.08 + 0.001 * step))
            lines.append(assign('BH(101)', 0.05 + 0.001 * step))
            for name in ('WING', 'HT', 'BW', 'BH'):
                cl, clb = block_values(step, name)
                for k in range(6):
                    lines.append(assign(f'{name}({21 + k})', cl[k]))
                    lines.append(assign(f'{name}({181 + k})', clb[k]))
            lines.append('      CALL SETUP2')
            lines.append("      WRITE(6,'(A,2I4,10ES25.16)') 'S',NF,"
                         "NALPHA,FLC(3),")
            lines.append('     1B(1),B(2),BHT(1),BHT(2),WINGIN(21),'
                         'WINGIN(41),HTIN(21),')
            lines.append('     2HTIN(41)')
            lines.append("      WRITE(6,'(A,14ES25.16)') 'SEC',"
                         "(SECOND(J),J=1,14)")
            lines.append("      WRITE(6,'(A,7ES25.16)') 'SEX',"
                         "(SECOND(J),J=17,23)")
            lines.append("      WRITE(6,'(A,3L2)') 'LOG',SUBSON,TRANSN,"
                         "SUPERS")
        lines.append(f' {7000 + n} CONTINUE')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    wcases, scases = wbclb_cases(), setup2_cases()
    stdout = run('second_level', driver(wcases, scases), ROUTINES)
    records = parse_records(stdout)
    # SETUP2 records hold repeated tags, one set per call; split them.
    raw = stdout.split('CASE')[1 + len(wcases):]
    setup2 = []
    for c, block in zip(scases, raw):
        calls, current = [], None
        for line in block.splitlines()[1:]:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'S':
                current = {'nf': int(fields[1]), 'nalpha': int(fields[2]),
                           'S': [float(v) for v in fields[3:]]}
                calls.append(current)
            elif fields[0] == 'LOG':
                current['LOG'] = [v == 'T' for v in fields[1:]]
            elif fields[0] in ('SEC', 'SEX'):
                current[fields[0]] = [float(v) for v in fields[1:]]
        setup2.append({'inputs': c, 'calls': calls})
    payload = {
        'wbclb': [{'inputs': c, 'outputs': r}
                  for c, r in zip(wcases, records[:len(wcases)])],
        'setup2': setup2,
        'blocks': {f'{step}{name}': block_values(step, name)
                   for step in range(1, 8)
                   for name in ('WING', 'HT', 'BW', 'BH')},
    }
    print(save('second_level', payload))


if __name__ == '__main__':
    main()
