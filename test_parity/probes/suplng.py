"""
Probe SUPLNG, the supersonic wing lift, pitching moment and drag, and save
the fixture.

Run from the repository root: ``python test_parity/probes/suplng.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['suplng', 'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook',
            'switch', 'tbfunx', 'quad', 'angdet', 'arcsin', 'fig60b', 'fwdxac', 'm27o33']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,MIDX,NALPHA
      COMMON /FLGTCD/ FLC(95)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /OPTION/ SW,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /SUPWH/  SLG(141)
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  HTIN(100)
      COMMON /IWING/ PWING,HT(400)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1    HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,HYPERS
      LOGICAL       FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1   HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,HYPERS
      LOGICAL DETACH
      EQUIVALENCE (DETACH,SLG(94))
      DATA WT/4HSTRA/,WD/4HDOUB/
"""

UNUSED = 1.0e-30
MIDX = 1
ALPHAS = [-2.0, 0.0, 0.5, 2.0, 5.0, 12.0]


def _edge(tan):
    """Sweep in degrees, sine and cosine of an edge of slope ``tan``."""
    s = math.atan(tan)
    return math.degrees(s), math.sin(s), math.cos(s)


def case(mach, tanle, straight=True, ar=3.0, taper=0.3, deltay=1.0,
         alphai=0.0, alphas=ALPHAS, tanleo=None, sweplo=None, tantei=0.1,
         tanteo=0.1, a_extra=None, win_extra=None, cl_unused=()):
    a = [0.25 * k for k in range(1, 196)]
    win = [0.125 * k for k in range(1, 101)]
    slg = [0.03125 * k for k in range(1, 142)]
    wing = [0.0078125 * k for k in range(1, 201)]
    for j in cl_unused:
        wing[20 + j] = UNUSED
    sweep, sinle, cosle = _edge(tanle)
    tanleo = tanle if tanleo is None else tanleo
    sweepo, _, cosleo = _edge(tanleo)
    words = {1: 150.0, 3: 300.0, 4: 280.0, 5: 2.0, 7: ar, 10: 12.0,
             16: 8.0, 23: 10.0, 25: taper, 26: 0.5, 27: taper,
             58: sweep, 60: sinle, 61: cosle, 62: tanle,
             76: _edge(tantei)[0], 80: tantei,
             82: sweep if sweplo is None else sweplo, 85: cosleo,
             86: tanleo, 100: _edge(tanteo)[0], 104: tanteo, 166: 6.0,
             167: 100.0, 168: 3.0, 169: 0.3}
    words.update(a_extra or {})
    for k, v in words.items():
        a[k - 1] = v
    wwords = {1: 2.0, 2: 6.0, 3: 14.0, 4: 15.0, 5: 10.0, 17: deltay,
              71: 0.5, 72: UNUSED}
    wwords.update(win_extra or {})
    for k, v in wwords.items():
        win[k - 1] = v
    slg[79], slg[80], slg[81] = 0.02, 1.0, 0.5
    syna = [0.0] * 19
    syna[0], syna[1], syna[3] = 20.0, 10.0, alphai
    return {'midx': MIDX, 'nalpha': len(alphas), 'mach': mach,
            'alpha': list(alphas), 'sw': 320.0, 'cbarr': 9.0,
            'straight': straight, 'syna': syna, 'a': a, 'wingin': win,
            'slg': slg, 'wing': wing,
            'a_words': sorted(words), 'win_words': sorted(wwords)}


def cases():
    wide = [-0.5, 0.0, 0.5, 1.2, 2.0, 4.0, 8.0, 15.0]
    return [
        case(1.5, 1.2),                                  # subsonic edge
        case(1.5, 2.2, ar=2.0, alphai=1.0),              # TLE/1.92 > 1
        case(2.0, 0.5, deltay=5.0),                      # detached shock
        case(2.0, 0.5, deltay=0.5, alphas=wide),         # attached
        case(2.0, 0.5, deltay=1.739, alphas=wide),       # detaches at 1 deg
        case(2.0, 0.0, taper=1.0, win_extra={72: 0.25}),  # rectangular
        case(1.5, 0.0, ar=0.5, taper=1.0),               # AR*BETA < 1
        case(1.5, -0.6),                                 # forward sweep
        case(1.6, 1.5, straight=False, tanleo=0.8, sweplo=20.0,
             tantei=0.3, tanteo=-0.1, cl_unused=(1,)),   # glove, extension
        case(1.6, 0.9, straight=False, tanleo=1.2, tantei=0.1,
             tanteo=0.1, win_extra={71: UNUSED}),        # round edge
        case(3.0, 0.9, straight=False, tanleo=0.6, sweplo=10.0),
        case(1.5, 1.2, alphai=1.0,
             alphas=[-4.0 + 1.5 * k for k in range(20)]),
    ]


def driver(all_cases, call='SUPLNG') -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             f'      MIDX={MIDX}']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,200',
                  '      IF(K.LE.195) A(K)=K*0.25',
                  '      IF(K.LE.100) HTIN(K)=K*0.125',
                  '      IF(K.LE.141) SLG(K)=K*0.03125',
                  '      IF(K.LE.200) HT(K)=K*0.0078125',
                  f'{100 + n:5d} CONTINUE']
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign(f'FLC({MIDX + 2})', c['mach']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + k})', v))
        lines.append(assign('SW', c['sw']))
        lines.append(assign('CBARR', c['cbarr']))
        for k, v in enumerate(c['syna']):
            lines.append(assign(f'SYNA({k + 1})', v))
        for k in c['a_words']:
            lines.append(assign(f'A({k})', c['a'][k - 1]))
        for k in c['win_words']:
            lines.append(assign(f'HTIN({k})', c['wingin'][k - 1]))
        lines.append(f"      HTIN(15)={'WT' if c['straight'] else 'WD'}")
        for k in (80, 81, 82):
            lines.append(assign(f'SLG({k})', c['slg'][k - 1]))
        for k, v in enumerate(c['wing']):
            if v == UNUSED:
                lines.append(assign(f'HT({k + 1})', v))
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,141ES25.16)') 'SLG',SLG")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,200ES25.16)') 'HT',(HT(K),K=1,200)")
        lines.append("      IDET=0")
        lines.append("      IF(DETACH) IDET=1")
        lines.append("      WRITE(6,'(A,ES25.16)') 'DET',FLOAT(IDET)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('suplng', driver(all_cases), ROUTINES))
    print(save('suplng', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))
    records = parse_records(run('m27o33', driver(all_cases, 'M27O33'),
                                ROUTINES))
    print(save('m27o33', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
