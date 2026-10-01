"""
Probe M14O16 (with LOARWB) and save the fixture.

LOARWB EQUIVALENCEs its REAL results onto ``/SUPDW/ LB(237)``, which is
implicitly INTEGER.  In the source's single precision both are one word,
but promoted reals are two, so neighbouring scalars (``LB(22)`` and
``LB(23)``, say) would overlap in the probe build; a layout patch declares
``LB`` REAL, which restores the original storage.

Run from the repository root: ``python test_parity/probes/loarwb.py``.
"""

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m14o16', 'loarwb', 'fig26', 'tbfunx', 'quad', 'tlinex',
            'tlin1x', 'glook', 'switch']

LAYOUT = {'loarwb': [('     1     KNBA,KLBA,LBIN,L,LOK',
                      '     1     KNBA,KLBA,LBIN,L,LOK,LB')]}

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /FLGTCD/ FLC(160)
      COMMON /OPTION/ SREFF,CBARR,ROUGFC,BLREF
      COMMON /SYNTSS/ SYNA(19)
      COMMON /SUPDW/  LB(237)
      COMMON /IBW/    PBW,BW(380)
      COMMON /POWER/  PWIN(29),LBIN(21)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      REAL LB, LBIN
      LOGICAL BLF, ROUNDN
      EQUIVALENCE (LBIN(14),BLF), (LBIN(17),ROUNDN)
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0]


def case(**kw):
    lbin = {1: 1.2, 2: 400.0, 3: UNUSED, 4: 30.0, 5: 1.5, 6: 0.02,
            7: 10.0, 8: 60.0, 9: 1200.0, 10: 25.0, 11: 30.0, 12: 4.0,
            13: 10.0, 14: False, 15: 35.0, 16: 20.0, 17: True, 18: 120.0,
            19: 12.0, 20: 30.0, 21: 38.0}
    c = {'alpha': ALPHA, 'lbin': lbin, 'mach': 0.6, 'reynolds': 2.0e6,
         'roughness': 1.6e-4, 'stale_xocrb': 0.037}
    for key, value in kw.items():
        if isinstance(key, str) and key in c:
            c[key] = value
    for index, value in kw.get('set', {}).items():
        c['lbin'][index] = value
    return c


def cases():
    return [
        case(),
        case(set={1: 0.0}),                              # zero-lift angle 0
        case(set={3: 25.0, 6: UNUSED, 17: False}),       # DELTEP, sharp nose
        case(set={14: True, 16: 12.0}),                  # BLF, blunter nose
        case(mach=1.4, reynolds=9.0e7, roughness=4e-3),  # cutoff Reynolds
        case(set={3: 45.0, 6: UNUSED, 1: -2.0, 5: 1.1}),
        case(set={6: 0.09, 16: 32.0, 18: 200.0, 19: 30.0}),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,380')
        lines.append('         BW(K)=0.')
        lines.append('         IF(K.LE.237) LB(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['reynolds']))
        lines.append(assign('ROUGFC', c['roughness']))
        lines.append(assign('LB(119)', c['stale_xocrb']))
        for index, value in c['lbin'].items():
            if index == 14:
                lines.append(f"      BLF=.{'TRUE' if value else 'FALSE'}.")
            elif index == 17:
                lines.append(f"      ROUNDN=.{'TRUE' if value else 'FALSE'}.")
            else:
                lines.append(assign(f'LBIN({index})', value))
        lines.append('      CALL M14O16')
        for tag, start in [('CD', 0), ('CL', 20), ('CM', 40), ('CN', 60),
                           ('CA', 80), ('CLA', 100), ('CMA', 120),
                           ('KYB', 140), ('KNB', 160), ('KLB', 180),
                           ('XCP', 200)]:
            lines.append(f"      WRITE(6,'(A,20ES25.16)') '{tag}',"
                         f"(BW({start}+J),J=1,NALPHA)")
        lines.append("      WRITE(6,'(A,237ES25.16)') 'LB',(LB(J),J=1,237)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'XCG',SYNA(1)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('loarwb', driver(all_cases), ROUTINES,
                                layout_patches=LAYOUT))
    payload = []
    for c, r in zip(all_cases, records):
        c = copy.deepcopy(c)
        c['lbin'] = {str(k): v for k, v in c['lbin'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('loarwb', payload))


if __name__ == '__main__':
    main()
