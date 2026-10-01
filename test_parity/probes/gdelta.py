"""
Probe GDELTA (with AGENR and SIMUL4), the flap spanwise loading, and save
the fixture.

GDELTA declares ``COMMON /FLOLOG/ X(20),ASYFP``: twenty REALs standing in
for twenty LOGICALs, which real-8 promotion doubles, so the build copy
makes the placeholder INTEGER to keep ``ASYFP`` in place.

Run from the repository root: ``python test_parity/probes/gdelta.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['gdelta', 'agenr', 'simul4', 'det4', 'tbfunx', 'quad',
            'tlinex', 'tlin1x', 'glook', 'switch']

PATCHES = {'gdelta': [('COMMON /FLOLOG/ X(20),ASYFP',
                       'COMMON /FLOLOG/ IX(20),ASYFP')]}

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLOLOG/ LG(20),ASYFP
      COMMON /HTDATA/ AHT(195)
      COMMON /SUPDW/  DW(35),TCD(58)
      COMMON /HTI/    HTIN(131)
      LOGICAL LG,ASYFP
      DIMENSION GD1(14),GD2(14),GD3(14),BOCH(4)
"""


def case(boch=(6.0, 6.5, 7.5, 9.0), sb=30.0, efi=0.2, efo=0.7,
         asyfp=False):
    return {'boch': list(boch), 'sb': sb, 'efi': efi, 'efo': efo,
            'asyfp': asyfp,
            'tail': {'tante': -0.1, 'tanle': 0.7, 'bsto2': 5.0,
                     'crh': 4.0}}


def cases():
    return [
        case(),
        case(sb=0.0, efi=0.0, efo=1.0),
        case(boch=(3.0, 3.2, 3.6, 4.4), sb=45.0, efi=0.35, efo=0.9),
        case(boch=(10.0, 11.0, 13.0, 16.0), sb=-20.0, efi=0.1, efo=0.5),
        case(sb=60.0, efi=0.6, efo=0.95),
        case(asyfp=True, sb=35.0),
        case(asyfp=True, sb=0.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      ASYFP=.{'TRUE' if c['asyfp'] else 'FALSE'}.")
        for k, v in enumerate(c['boch']):
            lines.append(assign(f'BOCH({k + 1})', v))
        t = c['tail']
        for name, key in [('AHT(80)', 'tante'), ('AHT(62)', 'tanle'),
                          ('HTIN(3)', 'bsto2'), ('AHT(10)', 'crh')]:
            lines.append(assign(name, t[key]))
        lines.append(assign('SB', c['sb']))
        lines.append(assign('EFI', c['efi']))
        lines.append(assign('EFO', c['efo']))
        lines.append(f'      DO {5000 + n} K=1,4')
        lines.append(f' {5000 + n} TCD(42+K)=0.')
        lines.append('      CALL GDELTA(GD1,GD2,GD3,EFI,EFO,BOCH,SB)')
        lines.append("      WRITE(6,'(A,46ES25.16)') 'R',GD1,GD2,GD3,")
        lines.append('     1(TCD(K),K=43,46)')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'B',BOCH")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('gdelta', driver(all_cases), ROUTINES,
                                layout_patches=PATCHES))
    print(save('gdelta', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
