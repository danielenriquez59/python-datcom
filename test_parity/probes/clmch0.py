"""
Probe CLMCH0 and M15O17 (CALCA0, WTLIFT, LIFTCF) and save the fixture.

Each case runs CLMCH0, the Mach-zero pass, and then M15O17 at a flight Mach
number with that Mach's section data, as the main program's loop does; the
second call inherits LIFTCF's angle state from the first.

Run from the repository root: ``python test_parity/probes/clmch0.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['clmch0', 'm15o17', 'calca0', 'wtlift', 'liftcf', 'clmxbs',
            'angles', 'tbfunx', 'quad', 'tlinex', 'tlin1x', 'tlin3x', 'glook',
            'switch']

STUBS = """\
      SUBROUTINE MESSGE(ROUT,MESS,A,B,C,D,E)
      DIMENSION ROUT(2),MESS(20),A(1),B(1),C(1),D(1),E(1)
      RETURN
      END
      SUBROUTINE M16O20
      RETURN
      END
"""

GEOMETRY = ['a3', 'a5', 'a7', 'a23', 'a27', 'a29', 'a34', 'a37', 'a38',
            'a40', 'a43', 'a50', 'a62', 'a74', 'a98', 'a123', 'a124',
            'a125', 'a168']

UNUSED = 1.0e-30


def case(rng, kind=1, a7=6.0, sweep=20.0, twist=0.0, camber=False,
         swafp=UNUSED, mach=0.5):
    tan_le = math.tan(math.radians(sweep))
    c4 = sweep - 3.0
    alpha = [-4., -2., 0., 2., 4., 8., 12., 16., 20., 24.]
    return {
        'planform_type': float(kind), 'alpha_deg': alpha,
        'geometry': {
            'a3': 150.0, 'a5': a7 * 0.6, 'a7': a7, 'a23': 5.0,
            'a27': rng.uniform(0.2, 0.7), 'a29': 18.0, 'a34': sweep,
            'a37': math.cos(math.radians(sweep)), 'a38': tan_le,
            'a40': c4, 'a43': math.cos(math.radians(c4)),
            'a50': math.tan(math.radians(sweep - 6.0)), 'a62': tan_le,
            'a74': 0.9 * tan_le, 'a98': 0.4 * tan_le,
            'a123': rng.uniform(0.1, 0.45), 'a124': 1.3, 'a125': 3.0,
            'a168': a7 * 1.4,
        },
        'section': {'swafp': swafp, 'alphai': 1.1, 'cli': 0.15,
                    'twista': twist, 'tovc': 0.09, 'camber': camber,
                    'deltay': rng.uniform(1.2, 3.0), 'xovc': 0.3,
                    'sspne': 13.0, 'cla_mach': 0.105, 'clamo': 0.1,
                    'clmaxl': 1.45, 'cla': 0.11, 'clmax': 1.35},
        'flight': {'mach': mach, 'beta': math.sqrt(1.0 - mach**2)},
        'sref': 170.0,
    }


def cases():
    rng = np.random.default_rng(1517)
    return [
        case(rng),
        case(rng, twist=-3.0),
        case(rng, camber=True, mach=0.7),
        case(rng, twist=2.5, camber=True),
        case(rng, swafp=-2.2),
        case(rng, a7=2.0, sweep=45.0),
        case(rng, kind=3, a7=5.0, sweep=35.0, mach=0.75),
        case(rng, kind=2, a7=2.0, sweep=60.0, mach=0.8),
        case(rng, kind=4, a7=5.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IG,IJK(3),NOVLY',
             '      COMMON /FLGTCD/ FLC(93)',
             '      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,',
             '     1  VTPL,VTSC,HEAD,PRPOWR,JEQPOW,LOASRT,TVTPAN,SUPERS,',
             '     2  SUBSON,TRANSN,HYPERS',
             '      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,',
             '     1  HEAD,PRPOWR,JEQPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,',
             '     2  HYPERS',
             '      COMMON /WINGD/  A(195),B(49)',
             '      COMMON /IWING/  PWING,WING(400)',
             '      COMMON /WINGI/  WINGIN(101)',
             '      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      DIMENSION WT(4)',
             '      INTEGER*8 ITRUE',
             '      EQUIVALENCE (ITRUE,RTRUE)',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      ITRUE=1',
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             '      TRANSN=.FALSE.']
    for n, c in enumerate(all_cases):
        g, s, f = c['geometry'], c['section'], c['flight']
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,400')
        lines.append('         WING(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.101) WINGIN(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append(f' {label} CONTINUE')
        # The main program sets NALPHA before CLMCH0, which only resets it
        # after its own pass (datcom.f, ahead of the CLMCH0 calls).
        lines.append(f"      FLC(2)={len(c['alpha_deg'])}.")
        lines.append(f"      NALPHA={len(c['alpha_deg'])}")
        lines.append(f"      WINGIN(15)=WT({int(c['planform_type'])})")
        for key in GEOMETRY:
            lines.append(assign(f'A({key[1:]})', g[key]))
        for j, a in enumerate(c['alpha_deg']):
            lines.append(assign(f'B({23 + j})', a))
        for index, key in [(10, 'swafp'), (20, 'alphai'), (19, 'cli'),
                           (11, 'twista'), (16, 'tovc'), (17, 'deltay'),
                           (18, 'xovc'), (3, 'sspne'), (21, 'cla_mach'),
                           (22, 'cla'), (69, 'clamo'), (68, 'clmaxl'),
                           (41, 'clmax'), (42, 'clmax')]:
            lines.append(assign(f'WINGIN({index})', s[key]))
        lines.append('      WINGIN(64)=' + ('RTRUE' if s['camber'] else '0.'))
        lines.append(assign('SREF', c['sref']))
        lines.append('      CALL CLMCH0(A,B,WINGIN(21),WINGIN(41),WINGIN(68),')
        lines.append('     1            WINGIN(69),WING,0)')
        lines.append("      WRITE(6,'(A,30ES25.16)') 'B3',(B(2+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,6ES25.16)') 'M0',B(48),A(126),A(127),")
        lines.append('     1  A(134),A(137),A(135)')
        lines.append("      WRITE(6,'(A,12ES25.16)') 'S0',(A(146+J),J=1,12)")
        # The Mach loop's second pass, with that Mach's section data.
        lines.append('      A(131)=WINGIN(22)')
        lines.append('      A(132)=WINGIN(42)')
        lines.append(assign('B(1)', f['mach']))
        lines.append(assign('B(2)', f['beta']))
        lines.append('      CALL M15O17')
        for tag, start in [('CL', 20), ('CN', 60)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"(WING({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,5ES25.16)') 'M1',WING(101),B(43),")
        lines.append('     1  B(44),B(49),B(45)')
        lines.append("      WRITE(6,'(A,12ES25.16)') 'S1',(A(146+J),J=1,12)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('clmch0', driver(all_cases), ROUTINES,
                                STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('clmch0', payload))


if __name__ == '__main__':
    main()
