"""
Probe SYNDIM, SECI and SECO, the geometry and airfoil-section set-up
routines, and save the fixture.

Run from the repository root: ``python test_parity/probes/section_setup.py``.
Every COMMON word is filled with a marker (``word*step``) first, so words
a routine leaves alone are checked too.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

LEGACY = Path(__file__).resolve().parent.parent.parent / 'datcom-legacy' / \
    'datcom_2000'
FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'section_setup.json')
U = 1.0e-30


def _declarations(stem, drop=()):
    text = (LEGACY / f'{stem}.f').read_text()
    lines = []
    for ln in text.splitlines()[1:]:
        if ln[:1] in 'cC':
            continue
        ln = ln[:72].rstrip()
        stmt = ln[6:].lstrip()
        if ln[5:6] not in (' ', '') and lines:
            lines.append(ln)
            continue
        if not stmt.split('(')[0].split()[0] in (
                'COMMON', 'DIMENSION', 'EQUIVALENCE', 'LOGICAL', 'REAL',
                'INTEGER', 'DATA'):
            break
        if any(d in stmt for d in drop):
            continue
        lines.append(ln)
    return '\n'.join(lines)


def _real(v):
    s = repr(float(v))
    return s.replace('e', 'D') if 'e' in s else s + 'D0'


def _records(out):
    cases, current = [], None
    for line in out.split('\n'):
        f = line.split()
        if not f:
            continue
        if f[0] == 'CASE':
            current = {}
            cases.append(current)
        else:
            current.setdefault(f[0], []).extend(
                float(v.replace('D', 'E')) for v in f[1:])
    return cases


SYNDIM = [
    {'bd77': 5.0, 'sspn': 15.0, 'sspne': 13.2, 'a62': 0.45, 'bd33': 20.0,
     'bd65': 6.0, 'bd74': 1.5, 'xcg': 22.0, 'xh': 40.0, 'alih': -2.0,
     'htin3': 5.0, 'htin4': 6.5, 'aht62': 0.6},
    {'bd77': -3.0, 'sspn': 10.0, 'sspne': 10.0, 'a62': 0.0, 'bd33': 12.0,
     'bd65': 4.0, 'bd74': -0.5, 'xcg': 11.0, 'xh': 30.0, 'alih': 0.0,
     'htin3': 4.0, 'htin4': 4.0, 'aht62': 0.2},
]


def syndim_driver():
    lines = ['      PROGRAM PROBE', _declarations('syndim'),
             '      UNUSED=1.D-30', '      RAD=57.2957795D0']
    for n, c in enumerate(SYNDIM):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,195', '      A(K)=K*0.01D0',
                  '      AHT(K)=K*0.02D0', f'{100 + n:5d} CONTINUE',
                  f'      DO {200 + n} K=1,100', '      BD(K)=K*0.03D0',
                  f'{200 + n:5d} CONTINUE']
        for name, v in (('BD(77)', c['bd77']), ('SSPN', c['sspn']),
                        ('SSPNE', c['sspne']), ('A(62)', c['a62']),
                        ('BD(33)', c['bd33']), ('BD(65)', c['bd65']),
                        ('BD(74)', c['bd74']), ('XCG', c['xcg']),
                        ('XH', c['xh']), ('ALIH', c['alih']),
                        ('HTIN(3)', c['htin3']), ('HTIN(4)', c['htin4']),
                        ('AHT(62)', c['aht62'])):
            lines.append(f'      {name}={_real(v)}')
        lines.append('      CALL SYNDIM')
        lines.append("      WRITE(6,'(A,100ES25.16)') 'BD',(BD(K),K=1,100)")
        lines.append("      WRITE(6,'(A,2ES25.16)') 'ARM',A(173),AHT(173)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def _section(atype, n_st):
    """A 380-word section block."""
    a = [0.001 * k for k in range(1, 381)]
    a[0] = atype
    a[81] = float(n_st)
    return a


SECI = [
    {'a': _section(1.0, 12), 'straight': True,
     'typein': {1: 3.0, 2: 6.0, 3: 12.0, 4: 15.0, 5: 5.0, 6: 9.0}},
    {'a': _section(2.0, 20), 'straight': False,
     'typein': {1: 3.0, 2: U, 3: 12.0, 4: 15.0, 5: 5.0, 6: 9.0}},
    {'a': _section(U, 7), 'straight': True,
     'typein': {1: 2.0, 2: 4.0, 3: 10.0, 4: 14.0, 5: 0.0, 6: 8.0}},
    {'a': _section(3.0, 30), 'straight': False,
     'typein': {1: 2.5, 2: 5.0, 3: 11.0, 4: 16.0, 5: 4.0, 6: 10.0}},
]


def seci_driver():
    lines = ['      PROGRAM PROBE', _declarations('seci'),
             '      DIMENSION AS(380),TIN(162)',
             '      LOGICAL CAMB', '      DATA WD/4HDOUB/',
             '      UNUSED=1.D-30']
    for n, c in enumerate(SECI):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,60', '      X(K)=-K*1.0D0',
                  '      XU(K)=-K*2.0D0', '      XL(K)=-K*3.0D0',
                  '      YU(K)=-K*4.0D0', '      YL(K)=-K*5.0D0',
                  '      THN(K)=-K*6.0D0', '      CAM(K)=-K*7.0D0',
                  f'{100 + n:5d} CONTINUE',
                  f'      DO {200 + n} K=1,380', '      AS(K)=K*0.001D0',
                  f'{200 + n:5d} CONTINUE',
                  f'      DO {300 + n} K=1,162', '      TIN(K)=K*0.5D0',
                  f'{300 + n:5d} CONTINUE']
        lines.append(f"      AS(1)={_real(c['a'][0])}")
        lines.append(f"      AS(82)={_real(c['a'][81])}")
        for k, v in c['typein'].items():
            lines.append(f'      TIN({k})={_real(v)}')
        lines.append(f"      TIN(15)={'STRA' if c['straight'] else 'WD'}")
        lines.append('      CALL SECI(AS,CAMB,AT,TIN,LOUT)')
        lines.append("      WRITE(6,'(A,ES25.16,2I5)') 'OUT',AT,LOUT,LL")
        for name in ('NACA', 'X', 'YU', 'YL', 'THN', 'CAM', 'CLA'):
            size = {'NACA': 80, 'CLA': 20}.get(name, 60)
            lines.append(f"      WRITE(6,'(A,{size}ES25.16)') '{name}',"
                         f"({name}(K),K=1,{size})")
        lines.append("      WRITE(6,'(A,ES25.16)') 'CBAR',CBAR")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


SECO = [
    {'atype': 1.0, 'camber': 1.0, 'nmach': 3, 'unset': [16, 17, 18, 21, 22,
                                                        41, 72, 61, 62, 68,
                                                        69, 10, 93, 63]},
    {'atype': -2.0, 'camber': 0.0, 'nmach': 20,
     'unset': [19, 20, 23, 40, 45, 60, 80, 91, 63, 65, 66, 67, 10]},
]
SECTION = {'tovc': 0.12, 'deltay': 2.4, 'xovc': 0.3, 'cli': 0.2, 'ai': 1.5,
           'cmco4': -0.05, 'rho': 0.016, 'clmax0': 1.4, 'cla0': 0.1,
           'alo': -2.0, 'covc': 0.03,
           'cla': [0.1 + 0.001 * m for m in range(20)],
           'clmax': [1.4 - 0.01 * m for m in range(20)],
           'xac': [0.25 + 0.002 * m for m in range(20)]}


def seco_driver():
    s = SECTION
    lines = ['      PROGRAM PROBE', _declarations('seco')]
    lines += [f'      TOVC={_real(s["tovc"])}',
              f'      DELTAY={_real(s["deltay"])}',
              f'      XOVC={_real(s["xovc"])}', f'      CLI={_real(s["cli"])}',
              f'      AI={_real(s["ai"])}', f'      CMCO4={_real(s["cmco4"])}',
              f'      RHO={_real(s["rho"])}',
              f'      CLMAX0={_real(s["clmax0"])}',
              f'      CLA0={_real(s["cla0"])}', f'      ALO={_real(s["alo"])}',
              f'      COVC={_real(s["covc"])}', '      UNUSED=1.D-30']
    for m in range(20):
        lines += [f'      CLA({m + 1})={_real(s["cla"][m])}',
                  f'      CLMAX({m + 1})={_real(s["clmax"][m])}',
                  f'      XAC({m + 1})={_real(s["xac"][m])}']
    for n, c in enumerate(SECO):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,162', '      A(K)=K*0.25D0',
                  f'{100 + n:5d} CONTINUE']
        for k in c['unset']:
            lines.append(f'      A({k})=UNUSED')
        lines.append(f"      FLC(1)={_real(c['nmach'])}")
        lines.append(f"      CALL SECO(A,{_real(c['camber'])},"
                     f"{_real(c['atype'])})")
        lines.append("      WRITE(6,'(A,162ES25.16)') 'A',A")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = {
        'syndim': {'cases': SYNDIM, 'records': _records(probe.run(
            'syndim', syndim_driver(), ['syndim']))},
        'seci': {'cases': SECI, 'records': _records(probe.run(
            'seci', seci_driver(), ['seci']))},
        'seco': {'cases': SECO, 'section': SECTION, 'records': _records(
            probe.run('seco', seco_driver(), ['seco']))},
    }
    FIXTURE.write_text(json.dumps(out))
    print(FIXTURE)


if __name__ == '__main__':
    main()
