"""
Probe SDWASH, the supersonic downwash and dynamic pressure at the tail,
and save the fixture.  INFTGM is stubbed: its ``A`` words are set by the
driver (INFTGM has its own tests).

Run from the repository root: ``python test_parity/probes/sdwash.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402

ROUTINES = ['sdwash', 'sddvc', 'sdwa', 'sdwb', 'sdwc', 'sdwd', 'sdwe',
            'dpresr', 'fig68', 'arcsin', 'arccos', 'mach2', 'tbfunx',
            'quad', 'trapz', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch']

EXTRA_STUBS = STUBS + """\
      SUBROUTINE INFTGM
      RETURN
      END
"""

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /OPTION/ SREF,CBARRW,ROUGFC,BLREF
      COMMON /OVERLY/ NLOG,NMACH,IMACH,NALPHA,IG,NF
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /FLGTCD/ FLC(160)
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /SUPDW/  DWA(236), JDETCH
      COMMON /WINGI/  WINGIN(101)
      COMMON /HTI/    HTIN(154)
      COMMON /SUPWH/  SLG(141)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IDWASH/ PDWASH, QQINFY(20), DWANGL(20), DEPDA(20)
"""

ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0]
THIN = [2.0, 1.2, 0.4, -0.4, -1.2, -2.0]


def case(mach=2.0, sweep_le=45.0, sweep_te=0.0, sweep_c2=20.0,
         user=False, nf=0, alpha=None, gamma=0.0, aliw=0.0, zh=1.5,
         taper=0.3):
    alpha = alpha or ALPHA
    na = len(alpha)
    return {
        'alpha': alpha, 'mach': mach, 'sref': 70.0, 'nf': nf,
        'user': user,
        'position': {'xw': 10.0, 'zw': 0.0, 'aliw': aliw, 'xh': 30.0,
                     'zh': zh},
        'a': {3: 60.0, 11: gamma, 12: zh, 16: 6.0, 24: 20.0,
              58: sweep_le, 70: sweep_c2, 76: sweep_te, 118: taper,
              120: 3.0},
        'win': {1: 1.8, 4: 10.0, 6: 6.0,
                **{95 + k: s for k, s in enumerate(THIN)}},
        'tail': {'tanle': 0.8, 'span': 4.0, 'spandi': 1.0, 'dihei': 2.0,
                 'diheo': 5.0},
        'wing': {'cl': [0.05 * (a + 1.0) for a in alpha],
                 'cla': [0.05 - 0.0005 * abs(a) for a in alpha],
                 'cd0': 0.012},
        'dwangl': [0.3 * a + 0.1 for a in alpha],
        'qqinfy': [-5.0] * na,
    }


def cases():
    return [
        case(),                                          # ICASE 1
        case(sweep_te=-10.0, sweep_c2=-5.0),             # ICASE 2
        case(sweep_te=-5.0, sweep_c2=10.0),              # ICASE 3
        case(sweep_te=10.0, sweep_c2=25.0),              # ICASE 4
        case(sweep_le=70.0, mach=1.5),                   # supersonic LE
        case(user=True),                                 # epsilon given
        case(user=True, sweep_le=70.0, mach=1.5),
        case(mach=1.05, sweep_le=10.0),                  # stops, subsonic
        case(mach=1.05, sweep_le=70.0),                  # stops, VISDW
        case(nf=-2, gamma=0.03, aliw=1.0),
        case(mach=3.0, zh=-1.0, alpha=[-6.0, -2.0, 2.0, 6.0]),
        case(mach=1.2, sweep_le=30.0, alpha=[0.0, 3.0, 6.0, 9.0]),
        case(mach=1.05, sweep_le=10.0, gamma=0.05),      # stops, returns
        case(mach=1.05, sweep_le=70.0, gamma=0.1,
             alpha=[0.0, 4.0, 8.0]),                     # stops at once
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,236')
        lines.append('         DWA(K)=0.')
        lines.append('         IF(K.LE.20) DEPDA(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        lines.append(f"      NF={c['nf']}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'WING({21 + j})', c['wing']['cl'][j]))
            lines.append(assign(f'WING({101 + j})', c['wing']['cla'][j]))
            lines.append(assign(f'DWANGL({1 + j})', c['dwangl'][j]))
            lines.append(assign(f'QQINFY({1 + j})', c['qqinfy'][j]))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('SLG(80)', c['wing']['cd0']))
        for name, value in c['position'].items():
            lines.append(assign(name.upper(), value))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        for k, v in c['win'].items():
            lines.append(assign(f'WINGIN({k})', v))
        t = c['tail']
        for name, key in [('AHT(62)', 'tanle'), ('HTIN(4)', 'span'),
                          ('HTIN(12)', 'spandi'), ('HTIN(13)', 'dihei'),
                          ('HTIN(14)', 'diheo')]:
            lines.append(assign(name, t[key]))
        lines.append(f"      CALL SDWASH(1,{1 if c['user'] else 0})")
        lines.append("      WRITE(6,'(A,236ES25.16)') 'DWA',")
        lines.append('     1(DWA(K),K=1,236)')
        lines.append("      WRITE(6,'(A,60ES25.16)') 'IDW',")
        lines.append('     1(QQINFY(K),K=1,20),(DWANGL(K),K=1,20),')
        lines.append('     2(DEPDA(K),K=1,20)')
        lines.append("      WRITE(6,'(A,3I6)') 'INT',NALPHA,NF,JDETCH")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('sdwash', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = []
    for c, r in zip(all_cases, records):
        c['a'] = {str(k): v for k, v in c['a'].items()}
        c['win'] = {str(k): v for k, v in c['win'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('sdwash', payload))


if __name__ == '__main__':
    main()
