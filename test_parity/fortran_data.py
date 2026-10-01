"""
Parse numeric DATA statements out of a fixed-form FORTRAN source file.

Used to generate embedded tables for translations, so no chart value is
ever transcribed by hand.  Handles continuation lines, comment lines,
several arrays in one DATA statement, and repeat counts such as ``3*.9``.
Hollerith constants (``4HSTRA``) are skipped, since they carry no numbers.

    python test_parity/fortran_data.py wtlift DCAR CLL

prints the named arrays of ``datcom-legacy/datcom_2000/wtlift.f``.
"""

import re
import sys
from pathlib import Path

LEGACY = Path(__file__).resolve().parent.parent / 'datcom-legacy' / 'datcom_2000'


def _statements(path: Path):
    """Yield each statement with its continuation lines joined."""
    current = None
    for raw in path.read_text(errors='replace').splitlines():
        if not raw or raw[0] in 'cC*!':
            continue
        line = raw[:72].ljust(6)
        if line[5] not in ' 0' and current is not None:
            current += line[6:]
            continue
        if current is not None:
            yield current
        current = line[6:]
    if current is not None:
        yield current


def _values(text: str):
    """Expand one slash-delimited value list."""
    out = []
    for token in text.split(','):
        token = token.strip().replace(' ', '')
        if not token:
            continue
        count, _, value = token.rpartition('*')
        if re.fullmatch(r'\d+H.*', value, re.IGNORECASE):
            return None
        number = float(value.replace('D', 'E').replace('d', 'e'))
        out.extend([number] * (int(count) if count else 1))
    return out


def parse(stem: str) -> dict:
    """Every numeric DATA array in ``<stem>.f``, as ``{NAME: [values]}``."""
    arrays = {}
    for statement in _statements(LEGACY / f'{stem}.f'):
        match = re.match(r'\s*DATA\s*(.*)', statement, re.IGNORECASE)
        if not match:
            continue
        # NAME /values/ pairs, optionally comma-separated.  A list of
        # scalar names shares one value list, one value each.
        for names, body in re.findall(
                r'([A-Z][A-Z0-9]*(?:\s*,\s*[A-Z][A-Z0-9]*)*)\s*/([^/]*)/',
                match.group(1), re.IGNORECASE):
            values = _values(body)
            if values is None:
                continue
            names = [n.strip().upper() for n in names.split(',')]
            if len(names) == 1:
                arrays[names[0]] = values
            elif len(names) == len(values):
                arrays.update({n: [v] for n, v in zip(names, values)})
            else:
                raise ValueError(f"DATA {','.join(names)} lists "
                                 f"{len(values)} values")
    return arrays


def literal(values, per_line: int = 8, indent: str = '    ') -> str:
    """Format values as the body of a Python list literal."""
    lines = []
    for i in range(0, len(values), per_line):
        lines.append(indent + ', '.join(repr(float(v)) for v in
                                        values[i:i + per_line]) + ',')
    return '\n'.join(lines)


if __name__ == '__main__':
    data = parse(sys.argv[1])
    for name in sys.argv[2:] or sorted(data):
        print(f'{name} ({len(data[name])}):')
        print(literal(data[name]))
