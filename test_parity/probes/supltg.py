"""
Probe SUPLTG, the supersonic horizontal-tail lift, pitching moment and
drag, and its overlay M22O26, and save the fixtures.

Run from the repository root: ``python test_parity/probes/supltg.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import assign, parse_records, run, save  # noqa: E402
from suplng import UNUSED, case as wing_case  # noqa: E402

ROUTINES = ['supltg', 'm22o26', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad', 'angdet', 'arcsin',
            'fig60b', 'fwdxac', 'fig26']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,MIDX,NALPHA,IG,IJKDUM(3),NOVLY
      COMMON /FLGTCD/ FLC(95)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /OPTION/ SW,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /SUPWH/  SLG(141),STG(141)
      COMMON /HTDATA/ A(195)
      COMMON /HTI/    HTIN(154)
      COMMON /IHT/    PHT,HT(380)
      COMMON /EXPER/  KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1                KBODY, KWING, KHT, KVT, KWB, KDWASH(3),
     2                ALPOW, ALPLW, ALPOH, ALPLH
      LOGICAL DETACH
      EQUIVALENCE (DETACH,STG(94))
      DATA WT/4HSTRA/,WD/4HDOUB/
"""

MIDX = 1


def case(*args, rl=3.0e6, ruff=0.0004, htin=None, a=None, **kw):
    c = wing_case(*args, **kw)
    # The tail sits at SYNA(6) with incidence SYNA(8).
    syna = c['syna']
    syna[5], syna[7] = 35.0, syna[3]
    syna[3] = 0.0
    win = c.pop('wingin') + [0.125 * k for k in range(101, 155)]
    hw = {6: 11.0, 62: 0.01, 63: 0.008, 70: 0.06}
    hw.update(htin or {})
    for k, v in hw.items():
        win[k - 1] = v
    aw = {2: 140.0, 15: 7.0, 17: 5.0, 18: 0.4, 122: 7.5}
    aw.update(a or {})
    for k, v in aw.items():
        c['a'][k - 1] = v
    c.update(htin=win, rl=rl, ruff=ruff,
             a_words=sorted(set(c['a_words']) | set(aw)),
             win_words=sorted(set(c['win_words']) | set(hw)))
    return c


def cases():
    wide = [-0.5, 0.0, 0.5, 1.2, 2.0, 4.0, 8.0, 15.0]
    return [
        case(1.5, 1.2),
        case(1.5, 2.2, ar=2.0, alphai=1.0, htin={71: UNUSED}),
        case(2.0, 0.5, deltay=5.0, ruff=0.0),            # stale RACH
        case(2.0, 0.5, deltay=1.739, alphas=wide),
        case(3.5, 0.5, htin={71: UNUSED}),               # RACH capped at 3
        case(2.0, 0.0, taper=1.0, win_extra={72: 0.25}),
        case(1.5, -0.6, ruff=0.0),
        case(1.6, 1.5, straight=False, tanleo=0.8, sweplo=20.0,
             tantei=0.3, tanteo=-0.1),
        case(1.6, 0.9, straight=False, tanleo=1.2, htin={71: UNUSED}),
        case(3.0, 0.9, straight=False, tanleo=0.6, sweplo=10.0,
             a={17: 7.0}),                               # equal chords
        case(1.5, 1.2, alphai=1.0, rl=1.0e8,
             alphas=[-4.0 + 1.5 * k for k in range(20)]),
    ]


def driver(all_cases, call) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             f'      MIDX={MIDX}', '      KDWASH(2)=0']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,200',
                  '      IF(K.LE.195) A(K)=K*0.25',
                  '      IF(K.LE.154) HTIN(K)=K*0.125',
                  '      IF(K.LE.141) STG(K)=K*0.03125',
                  '      HT(K)=K*0.0078125',
                  f'{100 + n:5d} CONTINUE']
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign(f'FLC({MIDX + 2})', c['mach']))
        lines.append(assign(f'FLC({MIDX + 42})', c['rl']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + k})', v))
        lines.append(assign('SW', c['sw']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append(assign('RUFF', c['ruff']))
        for k, v in enumerate(c['syna']):
            lines.append(assign(f'SYNA({k + 1})', v))
        for k in c['a_words']:
            lines.append(assign(f'A({k})', c['a'][k - 1]))
        for k in c['win_words']:
            lines.append(assign(f'HTIN({k})', c['htin'][k - 1]))
        lines.append(f"      HTIN(15)={'WT' if c['straight'] else 'WD'}")
        for k in (80, 81, 82):
            lines.append(assign(f'STG({k})', c['slg'][k - 1]))
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,141ES25.16)') 'SLG',STG")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,200ES25.16)') 'HT',(HT(K),K=1,200)")
        lines.append("      IDET=0")
        lines.append("      IF(DETACH) IDET=1")
        lines.append("      WRITE(6,'(A,ES25.16)') 'DET',FLOAT(IDET)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    for name, call in (('supltg', 'SUPLTG(0)'), ('m22o26', 'M22O26')):
        records = parse_records(run(name, driver(all_cases, call),
                                    ROUTINES))
        print(save(name, [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
