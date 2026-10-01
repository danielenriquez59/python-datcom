"""
Probe CMALPO (with FWDXAC) and save the fixture.  The surfaces are the
CMALPH probe's, plus the words only CMALPO reads.

Run from the repository root: ``python test_parity/probes/cmalpo.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import assign, parse_records, run, save  # noqa: E402
from cmalph import case  # noqa: E402

ROUTINES = ['cmalpo', 'fwdxac', 'tbfunx', 'quad', 'tlin3x', 'tlinex',
            'tlin1x', 'glook', 'switch']

GEOMETRY = ['a1', 'a3', 'a5', 'a7', 'a10', 'a16', 'a23', 'a26', 'a27',
            'a38', 'a40', 'a43', 'a62', 'a74', 'a86', 'a98', 'a125', 'a161',
            'a164', 'a166', 'a167', 'a168', 'a169', 'a173']


def cmalpo_case(rng, first_mach=0.5, **kw):
    c = case(rng, **kw)
    g = c['geometry']
    g.update({'a74': math.tan(math.radians(g['a34'] * 0.8)),
              'a98': math.tan(math.radians(g['a34'] * 0.5)),
              'a161': 7.5})
    return {'planform_type': c['planform_type'],
            'geometry': {k: g[k] for k in GEOMETRY},
            'section': {'sspn': 15.0, 'sspne': 13.2, 'cla': 0.105,
                        'cmo': -0.02, 'cmot': c['section']['cmot']},
            'mach': c['flight']['mach'], 'first_mach': first_mach,
            'cbarr': c['cbarr'], 'sref': c['sref']}


def cases():
    rng = np.random.default_rng(4153)
    return [
        cmalpo_case(rng, a7=6.0, a124=1.5),                  # A7 > A125
        cmalpo_case(rng, a7=2.0, sweep_le=35.0, mach=0.3),   # 26A
        cmalpo_case(rng, a7=1.6, sweep_le=60.0, mach=0.6),   # 26B
        cmalpo_case(rng, a7=2.5, sweep_le=0.0, mach=0.3),    # A38 = 0
        cmalpo_case(rng, a7=2.0, sweep_le=-30.0),            # FWDXAC
        cmalpo_case(rng, kind=2, a7=1.8, sweep_le=65.0, sweep_out=45.0),
        cmalpo_case(rng, kind=3, a7=5.0, sweep_le=40.0, sweep_out=20.0,
                    first_mach=0.9),
        cmalpo_case(rng, kind=3, a7=2.2, sweep_le=-20.0, sweep_out=-15.0),
        cmalpo_case(rng, kind=3, a7=2.2, sweep_le=-20.0, sweep_out=25.0),
        cmalpo_case(rng, kind=2, a7=2.2, sweep_le=50.0, sweep_out=-10.0),
        cmalpo_case(rng, kind=4, a7=2.0, sweep_le=55.0, sweep_out=58.0,
                    mach=0.9, first_mach=1.3),
        cmalpo_case(rng, kind=3, a7=2.0, sweep_le=0.0, sweep_out=30.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF',
             '      COMMON /FLGTCD/ FLC(95)',
             '      DIMENSION A(195),B(49),WINGIN(101),DYN(213),WT(4)',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,213')
        lines.append('         DYN(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.101) WINGIN(K)=0.')
        lines.append('         IF(K.LE.49) B(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      WINGIN(15)=WT({int(c['planform_type'])})")
        for key, value in c['geometry'].items():
            lines.append(assign(f'A({key[1:]})', value))
        s = c['section']
        for index, key in [(4, 'sspn'), (3, 'sspne'), (21, 'cla'),
                           (61, 'cmo'), (67, 'cmot')]:
            lines.append(assign(f'WINGIN({index})', s[key]))
        lines.append(assign('B(1)', c['mach']))
        lines.append(assign('FLC(3)', c['first_mach']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append('      CALL CMALPO(A,B,WINGIN,DYN)')
        lines.append("      WRITE(6,'(A,2ES25.16)') 'OUT',DYN(21),A(62)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('cmalpo', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('cmalpo', payload))


if __name__ == '__main__':
    main()
