"""
Compiled-FORTRAN parity harness for PyDATCOM.

Builds the original Digital DATCOM from `datcom-legacy/datcom.f`, runs the
11 sample cases from AFFDL-TR-79-3032, and compares the numbers it produces
against the reference output shipped alongside them.  This turns "matches my
reading of the source" into "matches execution", which is the check that
source-comparison alone cannot provide.

The build is the oracle for translation work: once it reproduces the shipped
reference, any routine can be probed by running a case through it and
comparing against the Python translation.

Usage, from the repository root:

    # Build the executable and extract the sample cases
    python test_parity/fortran_parity.py setup

    # Verify the build reproduces every shipped reference output
    python test_parity/fortran_parity.py verify

    # Run one case and keep its output
    python test_parity/fortran_parity.py run 2

On Windows this drives WSL; on Linux and macOS it runs natively.  Requires
gfortran (`sudo apt-get install gfortran`).

Status with gfortran 15.2.0 (`-w -std=legacy -fno-automatic -O1`): the build
compiles clean and 8 of the 11 cases reproduce the shipped reference to
printed precision.  Every residual difference has been identified, and none
indicates a bad build:

- **Lift-to-drag at exactly zero lift** (case 1).  The reference prints
  `******`, a format overflow from dividing by zero lift; this build prints
  a finite number.  That single column accounts for both of case 1's
  outliers and its one ragged line.
- **Underflow residue.**  The reference carries values such as -5.244E-12
  where this build accumulates to exact zero.  Handled by `ABSOLUTE_FLOOR`.
- **Signed zero.**  gfortran prints -0.00000E+00 for 0.00000E+00.
- **Last-place rounding** on values printed to three or four significant
  figures, plus a few iteration-sensitive quantities differing by about
  1e-3 absolute (cases 3 and 7).

The important consequence for translation work: **the oracle is the
executable, not the shipped listing.**  When a Python routine is checked
against a quantity produced by this binary, compiler rounding appears on
both sides of the comparison and cancels.  The shipped reference output is
a one-time sanity check that the build itself is sound, which it is.
"""

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEGACY = ROOT / 'datcom-legacy'
BUILD = ROOT / 'build' / 'fortran'
CASES = BUILD / 'cases'
EXE = 'datcom.exe'

# The shipped reference output is single precision printed to about four
# significant figures, so agreement is limited by print width rather than by
# arithmetic.  One unit in the last printed place is 1e-4 relative.
TOLERANCE = 2.0e-3

# Values this small are zero for every practical purpose.  The reference
# listings contain underflow residue such as -5.244E-12 where this build
# accumulates to exact zero; comparing those relatively yields a difference
# of 1.0 for two numbers that are both physically nothing.
ABSOLUTE_FLOOR = 1.0e-9

_NUMBER = re.compile(
    r'[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+\.?\d*[EeDd][-+]?\d+')


def _on_windows() -> bool:
    return platform.system() == 'Windows'


def _wsl_path(path: Path) -> str:
    """Translate a Windows path to its /mnt/<drive> form for WSL."""
    text = str(path.resolve()).replace('\\', '/')
    if len(text) > 1 and text[1] == ':':
        return f"/mnt/{text[0].lower()}{text[2:]}"
    return text


def shell(command: str, cwd: Path = None, stdin: str = None,
          timeout: int = 900) -> subprocess.CompletedProcess:
    """Run a POSIX shell command, through WSL when hosted on Windows."""
    directory = _wsl_path(cwd) if (cwd and _on_windows()) else (
        str(cwd) if cwd else None)
    if directory:
        command = f"cd '{directory}' && {command}"
    argv = (['wsl', '-d', 'Ubuntu', '--', 'bash', '-c', command]
            if _on_windows() else ['bash', '-c', command])
    return subprocess.run(argv, input=stdin, capture_output=True, text=True,
                          timeout=timeout)


# Progress and banner text that this build writes but the shipped listings
# do not, or do only inconsistently between cases.  None of it is output
# data, and leaving it in shifts the line alignment.
_NOISE = (
    'THIS SOFTWARE AND ANY ACCOMPANYING',
    'IS RELEASED "AS IS"',
    'WARRANTY OF ANY KIND',
    'THIS SOFTWARE AND ANY ACCOMPANYING DOCUMENTATION,',
    'INCLUDING, WITHOUT LIMITATION',
    'MERCHANTABILITY OR FITNESS',
    'IN NO EVENT WILL THE U.S. GOVERNMENT',
    'DAMAGES, INCLUDING LOST PROFITS',
    'INCIDENTAL OR CONSEQUENTIAL DAMAGES',
    'USE, OR INABILITY TO USE',
    'ACCOMPANYING DOCUMENTATION, EVEN IF',
    'OF THE POSSIBILITY OF SUCH DAMAGES',
    'Preparing to start the big loop',
    'At 1000',
    'Return to main program from',
)


def _is_noise(line: str) -> bool:
    return any(marker in line for marker in _NOISE)


def _line_numbers(path: Path):
    """Numbers in a listing, grouped by line, with program chatter removed.

    Grouping by line keeps the comparison aligned.  A flat sequence drifts
    irrecoverably as soon as one listing emits a token the other does not,
    turning a single formatting artifact into thousands of false mismatches.
    """
    lines = []
    with open(path, errors='replace') as handle:
        for line in handle:
            if _is_noise(line):
                continue
            lines.append([
                float(match.group().replace('D', 'E').replace('d', 'e'))
                for match in _NUMBER.finditer(line)])
    return lines


def extract_numbers(path: Path):
    """Every number in a DATCOM output listing, in order."""
    return [value for line in _line_numbers(path) for value in line]


def _difference(x: float, y: float) -> float:
    """Relative difference, with signed zero and underflow treated as zero.

    gfortran prints -0.00000E+00 where the original listing has
    0.00000E+00, and the two builds disagree on underflow residue near
    1e-12.  Both are numerically identical for any purpose this harness
    serves, so both collapse to no difference.
    """
    if abs(x) < ABSOLUTE_FLOOR and abs(y) < ABSOLUTE_FLOOR:
        return 0.0
    return abs(x - y) / max(abs(y), 1e-12) if y else abs(x - y)


def compare(ours: Path, reference: Path, tolerance: float = TOLERANCE):
    """Compare two listings numerically, line by line.

    Returns a dict with the count compared, how many exceed the tolerance,
    the worst relative difference and where it occurred.  Lines whose token
    counts differ are reported separately rather than silently shifting the
    rest of the comparison.
    """
    a, b = _line_numbers(ours), _line_numbers(reference)
    count = exceeded = ragged = 0
    worst, location = 0.0, None
    for index in range(min(len(a), len(b))):
        left, right = a[index], b[index]
        if len(left) != len(right):
            ragged += 1
            continue
        for position, (x, y) in enumerate(zip(left, right)):
            count += 1
            delta = _difference(x, y)
            if delta > tolerance:
                exceeded += 1
            if delta > worst:
                worst, location = delta, (index + 1, position, x, y)
    return {
        'count': count,
        'ours_lines': len(a),
        'reference_lines': len(b),
        'aligned': len(a) == len(b),
        'ragged_lines': ragged,
        'exceeded': exceeded,
        'worst': worst,
        'worst_line': location[0] if location else None,
        'worst_ours': location[2] if location else None,
        'worst_reference': location[3] if location else None,
    }


def setup() -> int:
    """Extract the sample cases and build the original program."""
    source = LEGACY / 'datcom.f'
    if not source.exists():
        print(f"error: {source} not found", file=sys.stderr)
        return 1

    BUILD.mkdir(parents=True, exist_ok=True)
    CASES.mkdir(parents=True, exist_ok=True)

    archive = LEGACY / 'exlinux.zip'
    if not archive.exists():
        archive = LEGACY / 'exwin.zip'
    if archive.exists():
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(CASES)
        print(f"extracted {len(list(CASES.glob('*.inp')))} cases from "
              f"{archive.name}")
    else:
        print("warning: no exlinux.zip or exwin.zip; cannot verify",
              file=sys.stderr)

    print("compiling datcom.f ...")
    result = shell(
        f"gfortran -w -std=legacy -fno-automatic -O1 "
        f"'{_wsl_path(source) if _on_windows() else source}' -o {EXE}",
        cwd=BUILD)
    if result.returncode != 0:
        print(result.stdout[-4000:], file=sys.stderr)
        print(result.stderr[-4000:], file=sys.stderr)
        return result.returncode
    print(f"built {BUILD / EXE}")
    return 0


def run_case(number: int, keep: bool = True) -> Path:
    """Run one sample case and return the path to its output listing."""
    name = f"ex{number}.inp"
    source = CASES / name
    if not source.exists():
        raise FileNotFoundError(f"{source} not found; run setup first")
    shutil.copy(source, BUILD / name)

    produced = BUILD / 'datcom.out'
    if produced.exists():
        produced.unlink()
    shell(f"./{EXE}", cwd=BUILD, stdin=f"{name}\n")
    if not produced.exists():
        raise RuntimeError(f"case {number} produced no datcom.out")
    if keep:
        target = BUILD / f"ex{number}.ours"
        shutil.copy(produced, target)
        return target
    return produced


def verify(tolerance: float = TOLERANCE) -> int:
    """Run every sample case and report agreement with the reference."""
    if not (BUILD / EXE).exists():
        print("error: executable missing; run setup first", file=sys.stderr)
        return 1

    print(f"{'case':>5} {'status':>7} {'numbers':>8} "
          f"{'over tol':>9} {'ragged':>7} {'worst rel':>11}")
    failures = 0
    for number in range(1, 12):
        reference = CASES / f"ex{number}.out"
        if not reference.exists():
            continue
        try:
            ours = run_case(number)
        except (RuntimeError, FileNotFoundError) as error:
            print(f"{number:>5} {'ERROR':>7}   {error}")
            failures += 1
            continue
        stats = compare(ours, reference, tolerance)
        ok = (stats['aligned'] and stats['worst'] < 1.0e-2 and
              stats['ragged_lines'] == 0)
        failures += 0 if ok else 1
        note = '' if stats['aligned'] else (
            f"  lines ours={stats['ours_lines']} "
            f"ref={stats['reference_lines']}")
        print(f"{number:>5} {'OK' if ok else 'CHECK':>7} "
              f"{stats['count']:>8} {stats['exceeded']:>9} "
              f"{stats['ragged_lines']:>7} {stats['worst']:>11.3e}{note}")
    print()
    if failures:
        print(f"{failures} case(s) need review")
    else:
        print("all cases reproduce the shipped reference output")
    return 1 if failures else 0


def report(number: int, tolerance: float = TOLERANCE, limit: int = 20) -> int:
    """List the individual values of one case that exceed the tolerance."""
    reference = CASES / f"ex{number}.out"
    ours = BUILD / f"ex{number}.ours"
    if not ours.exists():
        ours = run_case(number)
    if not reference.exists():
        print(f"no reference for case {number}", file=sys.stderr)
        return 1

    a, b = _line_numbers(ours), _line_numbers(reference)
    print(f"{'line':>6} {'pos':>4} {'ours':>14} {'reference':>14} {'rel':>11}")
    shown = 0
    for index in range(min(len(a), len(b))):
        left, right = a[index], b[index]
        if len(left) != len(right):
            print(f"{index + 1:>6}  ragged: {len(left)} vs {len(right)} values")
            shown += 1
            continue
        for position, (x, y) in enumerate(zip(left, right)):
            delta = _difference(x, y)
            if delta > tolerance:
                print(f"{index + 1:>6} {position:>4} {x:>14.6g} "
                      f"{y:>14.6g} {delta:>11.3e}")
                shown += 1
            if shown >= limit:
                print(f"... stopping at {limit}")
                return 0
    if not shown:
        print("no values exceed the tolerance")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('setup', help='extract cases and build the executable')
    verify_parser = sub.add_parser('verify', help='check all 11 cases')
    verify_parser.add_argument('--tolerance', type=float, default=TOLERANCE)
    run_parser = sub.add_parser('run', help='run a single case')
    run_parser.add_argument('case', type=int)
    report_parser = sub.add_parser('report', help='detail one case')
    report_parser.add_argument('case', type=int)
    report_parser.add_argument('--tolerance', type=float, default=TOLERANCE)

    args = parser.parse_args()
    if args.command == 'setup':
        return setup()
    if args.command == 'verify':
        return verify(args.tolerance)
    if args.command == 'report':
        return report(args.case, args.tolerance)
    if args.command == 'run':
        path = run_case(args.case)
        print(f"wrote {path}")
        reference = CASES / f"ex{args.case}.out"
        if reference.exists():
            stats = compare(path, reference)
            print(f"compared {stats['count']} numbers, "
                  f"worst relative difference {stats['worst']:.3e}")
        return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())
