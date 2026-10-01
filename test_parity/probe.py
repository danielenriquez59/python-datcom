"""
Compiled-routine probes: run one legacy routine under a generated driver.

`fortran_parity.py` checks the whole program against the shipped listings.
A probe goes the other way: it links a single legacy routine and its callees
against a driver that loads chosen values into the routine's COMMON blocks,
calls it, and prints what it wrote back.  The output is saved as a JSON
fixture under ``tests/fixtures/probes/``, so the Python tests compare against
execution without needing a compiler.

Probes build with ``-fdefault-real-8 -fdefault-double-8``.  The legacy
arithmetic is single precision, which limits agreement with the double
precision translation to about 1e-6 relative; promoting it lets the tests
hold 1e-9, tight enough to expose a transposed index or a wrong constant
that single-precision noise would hide.  Every real literal is promoted
with it, so no DATA value changes.

Each routine has a script in ``test_parity/probes/`` that builds its cases and
calls :func:`run`; rerun it to regenerate the fixture.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fortran_parity import ROOT, _on_windows, _wsl_path, shell  # noqa: E402

LEGACY = ROOT / 'datcom-legacy' / 'datcom_2000'
PROBES = ROOT / 'build' / 'probe'
FIXTURES = ROOT / 'tests' / 'fixtures' / 'probes'

FLAGS = '-w -std=legacy -fno-automatic -O0 -fdefault-real-8 -fdefault-double-8'

# MESSGE only prints extrapolation diagnostics; stubbing it keeps the
# probe from pulling in the output system.
STUBS = """\
      SUBROUTINE MESSGE(ROUT,MESS,A,B,C,D,E)
      DIMENSION ROUT(2),MESS(20),A(1),B(1),C(1),D(1),E(1)
      RETURN
      END
"""


def fortran_real(value: float) -> str:
    """A real literal that survives fixed-form parsing at full precision."""
    return f"{float(value):.17E}".replace('E', 'D')


def assign(name: str, value: float) -> str:
    """One fixed-form assignment statement, continued if it is long."""
    text = f"{name}={fortran_real(value)}"
    lines = []
    while len(text) > 66:
        lines.append(text[:66])
        text = text[66:]
    lines.append(text)
    return '\n'.join(
        ('      ' if i == 0 else '     1') + part for i, part in enumerate(lines))


def run(name: str, driver: str, routines, stubs: str = STUBS,
        layout_patches=None) -> str:
    """Compile ``driver`` with the named legacy routines and return stdout.

    Args:
        name: Probe name, used for the build directory.
        driver: Fixed-form FORTRAN main program.
        routines: Legacy routine file stems, such as ``['dwash', 'tbfunx']``.
        stubs: Extra FORTRAN source, by default a silent MESSGE.
        layout_patches: ``{stem: [(old, new), ...]}`` textual replacements
            applied to a build copy of a routine.  Only for restoring a
            COMMON layout that real-8 promotion breaks: a REAL placeholder
            standing in for 4-byte LOGICAL or INTEGER words, as in
            ``COMMON /FLOLOG/ FLTC,OPTI,BO,XX(14),TRANSN``, grows to 8 bytes
            and shifts every later name.  Each replacement must match.

    Raises:
        RuntimeError: If compilation or execution fails.
        ValueError: If a layout patch does not match its routine.
    """
    # Fixed form ignores everything past column 72, silently, so an over-
    # long statement loses its tail and still compiles.
    for number, line in enumerate(driver.splitlines(), 1):
        if len(line) > 72 and line[:1] not in 'cC*!':
            raise ValueError(f"probe {name} driver line {number} runs past "
                             f"column 72: {line}")
    directory = PROBES / name
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'driver.f').write_text(driver)
    (directory / 'stubs.f').write_text(stubs)
    legacy = (_wsl_path(LEGACY) if _on_windows() else str(LEGACY))
    layout_patches = layout_patches or {}
    for stem, replacements in layout_patches.items():
        text = (LEGACY / f'{stem}.f').read_text(errors='replace')
        for old, new in replacements:
            if old not in text:
                raise ValueError(f"layout patch for {stem} does not match: "
                                 f"{old!r}")
            text = text.replace(old, new)
        (directory / f'{stem}_layout.f').write_text(text)
    sources = ' '.join(
        f"'{stem}_layout.f'" if stem in layout_patches
        else f"'{legacy}/{stem}.f'" for stem in routines)
    build = shell(f"gfortran {FLAGS} driver.f stubs.f {sources} -o probe",
                  cwd=directory)
    if build.returncode != 0:
        raise RuntimeError(f"probe {name} failed to build:\n{build.stderr}")
    result = shell('./probe', cwd=directory)
    if result.returncode != 0:
        raise RuntimeError(f"probe {name} failed to run:\n{result.stderr}")
    return result.stdout


def parse_records(stdout: str):
    """Split probe output into ``{tag: [floats]}`` records per case.

    Drivers print ``CASE n`` to start a case and then lines of the form
    ``TAG v1 v2 ...``; a tag repeated within a case extends its list.  Any
    other line is the routine's own printing, kept under ``'_text'``.
    """
    cases, current = [], None
    for line in stdout.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'CASE':
            current = {}
            cases.append(current)
            continue
        try:
            values = [float(f.replace('D', 'E')) for f in fields[1:]]
        except ValueError:
            current.setdefault('_text', []).append(line.strip())
            continue
        current.setdefault(fields[0], []).extend(values)
    return cases


def save(name: str, payload) -> Path:
    """Write a probe fixture and return its path."""
    FIXTURES.mkdir(parents=True, exist_ok=True)
    path = FIXTURES / f'{name}.json'
    path.write_text(json.dumps(payload, indent=1) + '\n')
    return path
