"""
Probe overlays M31O37 (wing) and M33O41 (horizontal tail) and save the
fixture.  Each runs CMALPH, CACALC and the slope and linearity pass.

Run from the repository root: ``python test_parity/probes/m31o37.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import STUBS, assign, parse_records, run, save  # noqa: E402
import cmalph  # noqa: E402

ROUTINES = ['m31o37', 'm33o41', 'cacalc'] + cmalph.ROUTINES

EXTRA_STUBS = STUBS + """\
      SUBROUTINE EXSUBT
      RETURN
      END
"""

INCIDENCE = 1.0


def cases():
    out = []
    for n, base in enumerate(cmalph.cases()):
        for surface in ('wing', 'tail'):
            for experimental in (False, True):
                if experimental and n % 3:
                    continue
                c = dict(base)
                local = np.array(base['alpha_deg'])
                cla = base['lift']['cla']
                # A lift curve that turns over past 16 degrees, so the
                # slope departs from linear.
                cl = cla * local * np.where(local > 16.0,
                                            np.exp(-(local - 16.0) / 10.0),
                                            1.0)
                c.update({'surface': surface, 'experimental': experimental,
                          'free_alpha_deg': list(local - INCIDENCE),
                          'cl': list(cl),
                          'cd': list(0.008 + 0.05 * cl**2)})
                out.append(c)
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE',
             '      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,IJK(3),NOVLY',
             '      COMMON /CONSNT/ PI,DR,UNUSED,RAD',
             '      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF',
             '      COMMON /WINGD/  A(195),B(49)',
             '      COMMON /IWING/  PWING,WING(400)',
             '      COMMON /WINGI/  WINGIN(101)',
             '      COMMON /HTDATA/ AHT(195),BHT(49)',
             '      COMMON /IHT/    PHT,HT(380)',
             '      COMMON /HTI/    HTIN(154)',
             '      COMMON /WHAERO/ C(51),D(55),CHT(51)',
             '      COMMON /FLGTCD/ FLC(95)',
             '      COMMON /EXPER/ KLIST, NLIST(100), NNAMES, IMACH, MDATA,',
             '     1               KBODY, KWING, KHT',
             '      LOGICAL MDATA,KBODY,KWING,KHT',
             '      DIMENSION WT(4)',
             "      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/",
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=2']
    for n, c in enumerate(all_cases):
        tail = c['surface'] == 'tail'
        aa, bb, win, arr = (('AHT', 'BHT', 'HTIN', 'HT') if tail else
                            ('A', 'B', 'WINGIN', 'WING'))
        label = 1000 + n
        lines.append(f"      WRITE(6,'(A,I3)') 'CASE',{n}")
        lines.append(f'      DO {label} K=1,380')
        lines.append(f'         {arr}(K)=0.')
        lines.append(f'         IF(K.LE.195) {aa}(K)=0.')
        lines.append(f'         IF(K.LE.101) {win}(K)=0.')
        lines.append('         IF(K.LE.51) C(K)=0.')
        lines.append('         IF(K.LE.51) CHT(K)=0.')
        lines.append(f'         IF(K.LE.49) {bb}(K)=0.')
        lines.append(f' {label} CONTINUE')
        flag = 'TRUE' if c['experimental'] else 'FALSE'
        lines.append(f"      KWING=.{flag}.")
        lines.append(f"      KHT=.{flag}.")
        lines.append(f"      NALPHA={len(c['alpha_deg'])}")
        lines.append(f"      {win}(15)=WT({int(c['planform_type'])})")
        for key in cmalph.GEOMETRY:
            lines.append(assign(f'{aa}({key[1:]})', c['geometry'][key]))
        s = c['section']
        for index, key in [(61, 'cmo'), (67, 'cmot'), (11, 'twista'),
                           (17, 'deltay'), (73, 'xac')]:
            lines.append(assign(f'{win}({index})', s[key]))
        for j, a in enumerate(c['alpha_deg']):
            lines.append(assign(f'{bb}({23 + j})', a))
            lines.append(assign(f'FLC({23 + j})', c['free_alpha_deg'][j]))
            lines.append(assign(f'{arr}({1 + j})', c['cd'][j]))
            lines.append(assign(f'{arr}({21 + j})', c['cl'][j]))
            lines.append(assign(f'{arr}({61 + j})', c['lift']['cn'][j]))
        lines.append(assign(f'{arr}(101)', c['lift']['cla']))
        lines.append(assign(f'{bb}(43)', c['lift']['alpha_clmax']))
        lines.append(assign(f'{bb}(1)', c['flight']['mach']))
        lines.append(assign(f'{bb}(2)', c['flight']['beta']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append('      CALL M33O41' if tail else '      CALL M31O37')
        for tag, start in [('CM', 40), ('CN', 60), ('CA', 80), ('CLA', 100),
                           ('CMA', 120)]:
            lines.append(f"      WRITE(6,'(A,30ES25.16)') '{tag}',"
                         f"({arr}({start}+J),J=1,NALPHA)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('m31o37', driver(all_cases), ROUTINES,
                                EXTRA_STUBS))
    payload = [{'inputs': c, 'outputs': r} for c, r in zip(all_cases, records)]
    print(save('m31o37', payload))


if __name__ == '__main__':
    main()
