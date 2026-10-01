"""
Probe TRSONI (with TRANWG and CLMXB1) and save the fixture, then run the
same configurations through TRSONJ (with TRNHT) on the tail's blocks.

Each configuration runs a sequence of Mach numbers without clearing the
``TRA`` work array in between, as the main loop does, so the wave-drag
points TRSONI never stores carry over exactly as executed.

Run from the repository root: ``python test_parity/probes/trsoni.py``.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['trsoni', 'tranwg', 'clmxb1', 'interx', 'tlin1x', 'tlinex',
            'tlin3x', 'glook', 'switch', 'quad', 'tbfunx', 'tranf', 'fig26']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IG,NF
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /FLGTCD/ FLC(93)
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195),B(49)
      COMMON /WINGI/  WINGIN(77)
      COMMON /SUPBOD/ SBD(227)
      COMMON /SBETA/  STB(135),TRA(108)
      COMMON /IBODY/  PB, BODY(400)
      COMMON /IWING/  PW, WING(400)
      COMMON /IBW/    PBW, BW(380)
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
      DIMENSION WT(4)
"""

ALPHA = [-2., 0., 2., 4., 8.]
MACHS = [0.85, 0.95, 1.05, 1.2, 1.3]


def config(kind=1.0, ar=5.0, sweep=35.0, taper=0.4, toc=0.06, xovc=0.3,
           bo=True, nf=0, a128=3.0, deltay=1.6):
    le = math.radians(sweep)
    c2 = math.radians(sweep - 6.0)
    c4 = math.radians(sweep - 3.0)
    return {
        'wing': {'type': kind, 'tovc': toc, 'deltay': deltay, 'xovc': xovc,
                 'clamo': 0.105, 'cd0': 0.0071},
        'a': {3: 140.0, 7: ar, 16: 6.5, 25: taper, 27: taper, 58: sweep,
              61: math.cos(le), 62: math.tan(le), 67: math.cos(c4),
              71: c2, 74: math.tan(c2), 86: math.tan(le) * 0.9,
              128: a128, 129: 2.0e6},
        'sref': 160.0, 'roughness': 1.6e-4, 'nf': nf,
        'body': ({'cla_06': 0.0042, 'cma_06': 0.021, 'cf_06': 0.0028,
                  'cd0_06': 0.014, 'cd_base_06': 0.0018,
                  'wetted_area': 380.0, 'base_area': 3.0, 'length': 42.0,
                  'base_diameter': 1.0, 'cla_14': 0.0046, 'cma_14': 0.026,
                  'max_diameter': 4.5, 'cd0_14': 0.032} if bo else None),
        'alpha': ALPHA, 'machs': MACHS,
    }


def configs():
    return [
        config(),
        config(ar=1.8, sweep=55.0, taper=0.2),        # low AR, A(160) > 4.5
        config(ar=1.5, sweep=20.0, taper=0.5, xovc=0.4),  # low AR, <= 4.5
        config(ar=8.0, sweep=15.0, toc=0.1, deltay=2.5),
        config(kind=3.0),                             # wing part skipped
        config(bo=False),
        config(nf=-1),
        config(ar=2.5, sweep=65.0, taper=0.3),        # low AR, A(160) > 4.5
    ]


# TRSONJ reads the tail's blocks; its /SBETA/ TRA starts at word 244.
TAIL_COMMONS = (COMMONS
                .replace('/WINGD/  A(195),B(49)', '/HTDATA/ A(195),B(49)')
                .replace('/WINGI/  WINGIN(77)', '/HTI/    WINGIN(154)')
                .replace('STB(135),TRA(108)', 'STB(243),TRA(108)')
                .replace('/IWING/  PW, WING(400)', '/IHT/    PW, WING(380)')
                .replace('/IBW/    PBW, BW(380)', '/IBH/    PBW, BW(380)'))


def driver(all_configs, tail=False) -> str:
    commons = TAIL_COMMONS if tail else COMMONS
    lines = ['      PROGRAM PROBE', commons.rstrip('\n'),
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    n = 0
    for c in all_configs:
        label = 5000 + n
        lines.append(f'      DO {label} K=1,400')
        lines.append('         BODY(K)=0.')
        lines.append('         IF(K.LE.380) WING(K)=0.')
        lines.append('         IF(K.LE.380) BW(K)=0.')
        lines.append('         IF(K.LE.227) SBD(K)=0.')
        lines.append('         IF(K.LE.195) A(K)=0.')
        lines.append('         IF(K.LE.108) TRA(K)=0.')
        lines.append(f' {label} CONTINUE')
        lines.append(f"      NF={c['nf']}")
        lines.append(f"      BO=.{'TRUE' if c['body'] else 'FALSE'}.")
        lines.append(f"      NALPHA={len(c['alpha'])}")
        w = c['wing']
        lines.append(f"      WINGIN(15)=WT({int(w['type'])})")
        for index, key in [(16, 'tovc'), (17, 'deltay'), (18, 'xovc'),
                           (69, 'clamo')]:
            lines.append(assign(f'WINGIN({index})', w[key]))
        for j, alpha in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', alpha))
        lines.append(assign('SR', c['sref']))
        lines.append(assign('RUFF', c['roughness']))
        for mach in c['machs']:
            lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
            n += 1
            # Inputs the pass overwrites are reloaded each Mach, as the
            # program's upstream routines would.
            for index, value in c['a'].items():
                lines.append(assign(f'A({index})', value))
            lines.append(assign('WING(1)', w['cd0']))
            b = c['body']
            if b:
                for name, key in [('BODY(101)', 'cla_06'),
                                  ('BODY(121)', 'cma_06'),
                                  ('BD(92)', 'cf_06'), ('BD(61)', 'cd0_06'),
                                  ('BD(60)', 'cd_base_06'),
                                  ('BD(93)', 'wetted_area'),
                                  ('BD(57)', 'base_area'),
                                  ('SBD(2)', 'length'),
                                  ('SBD(6)', 'base_diameter'),
                                  ('SBD(18)', 'cla_14'),
                                  ('SBD(110)', 'cma_14'),
                                  ('SBD(120)', 'max_diameter'),
                                  ('SBD(124)', 'cd0_14')]:
                    lines.append(assign(name, b[key]))
            lines.append(assign('FLC(3)', mach))
            lines.append(f"      CALL {'TRSONJ' if tail else 'TRSONI'}(1)")
            lines.append("      WRITE(6,'(A,5ES25.16)') 'W',WING(101),"
                         "WING(1),A(160),A(62),A(86)")
            lines.append("      WRITE(6,'(A,82ES25.16)') 'TRA',"
                         "(TRA(J),J=1,82)")
            lines.append("      WRITE(6,'(A,4ES25.16)') 'B',BODY(101),"
                         "BODY(121),BW(1),SBD(6)")
            lines.append("      WRITE(6,'(A,20ES25.16)') 'CDJ',"
                         "(BODY(J),J=1,NALPHA)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_configs = configs()
    for name, routines, tail in [
            ('trsoni', ROUTINES, False),
            ('trsonj', [r.replace('trsoni', 'trsonj').replace('tranwg', 'trnht')
                        for r in ROUTINES], True)]:
        records = parse_records(run(name, driver(all_configs, tail),
                                    routines))
        payload, n = [], 0
        for c in all_configs:
            runs = []
            for mach in c['machs']:
                runs.append({'mach': mach, 'outputs': records[n]})
                n += 1
            payload.append({'inputs': c, 'runs': runs})
        print(save(name, payload))


if __name__ == '__main__':
    main()
