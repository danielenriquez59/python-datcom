"""
Probe WTGEOM (after SETUP1's sweep records) and save the fixture.

Run from the repository root: ``python test_parity/probes/wtgeom.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['wtgeom', 'angles', 'zerang']
UNUSED = 1.0e-30


def case(chrdtp=3.0, sspnop=0.0, sspne=13.0, sspn=15.0, chrdbp=None,
         chrdr=7.0, savsi=25.0, savso=None, chstat=0.25, xovc=0.30,
         xovco=0.35, repeat=False):
    single = sspnop < 1e-29
    return {
        'ain': {1: chrdtp, 2: sspnop, 3: sspne, 4: sspn,
                5: chrdtp if chrdbp is None else chrdbp, 6: chrdr,
                9: chstat, 66: xovco},
        'savsi': savsi, 'savso': savsi if savso is None else savso,
        'xovc': xovc, 'repeat': repeat,
    }


def cases():
    return [
        case(),
        case(chstat=0.0),
        case(chstat=0.5, savsi=-20.0),
        case(chstat=1.0, savsi=10.0),
        case(chstat=0.3, xovc=0.3),                 # CHSTAT at max t/c
        case(chstat=0.42, savsi=45.0, chrdtp=0.0),  # pointed tip
        case(sspnop=6.0, chrdbp=4.0, savsi=40.0, savso=20.0),
        case(sspnop=6.0, chrdbp=4.0, savsi=-25.0, savso=-15.0, chstat=0.0),
        case(sspnop=8.0, chrdbp=3.5, savsi=65.0, savso=45.0, chstat=0.0,
             chrdtp=1.0, xovco=0.3),
        case(sspnop=5.0, chrdbp=4.5, savsi=30.0, savso=10.0, chstat=0.25,
             xovc=0.25),
        case(repeat=True),
        case(sspnop=6.0, chrdbp=4.0, savsi=40.0, savso=20.0, repeat=True),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD',
             '      DIMENSION A(195),AIN(77)',
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        if not c['repeat']:
            label = 1000 + n
            lines.append(f'      DO {label} K=1,195')
            lines.append('         A(K)=0.')
            lines.append('         IF(K.LE.77) AIN(K)=0.')
            lines.append(f' {label} CONTINUE')
        # The block as it stood, for a repeated call.
        lines.append("      WRITE(6,'(A,195ES25.16)') 'BEFORE',"
                     "(A(J),J=1,195)")
        for index, value in c['ain'].items():
            lines.append(assign(f'AIN({index})', value))
        lines.append(assign('A(106)', c['savsi']))
        lines.append(assign('A(112)', c['savso']))
        lines.append(assign('A(174)', c['xovc']))
        lines.append('      CALL ANGLES(1,A(106))')
        lines.append('      CALL ANGLES(1,A(112))')
        lines.append("      WRITE(6,'(A,12ES25.16)') 'REC',(A(105+J),J=1,12)")
        lines.append('      CALL WTGEOM(A,AIN)')
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',(A(J),J=1,195)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'AIN5',AIN(5)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('wtgeom', driver(all_cases), ROUTINES))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('wtgeom', payload))


if __name__ == '__main__':
    main()
