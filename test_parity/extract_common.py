"""
Extract COMMON block dumps from a Digital DATCOM output listing.

The `DUMP` card makes the original program print named COMMON blocks as
`NAME( n)= value` records.  Example Problem 2 carries `DUMP A`, which prints
the 195-element WINGD array that WTGEOM fills, so a run of the compiled
original yields the exact intermediate geometry a translation must match.

Used by tests/test_fortran_parity.py to pin Python translations against
values produced by executed FORTRAN rather than by reading the source.

Usage:
    python test_parity/extract_common.py build/fortran/ex2.ours A
    python test_parity/extract_common.py build/fortran/ex2.ours A --index 16
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List

# Matches records such as "A( 16)= 2.01047E+00" and "BODY(121)= 3.41518E-03".
_RECORD = re.compile(
    r'\b([A-Z][A-Z0-9]*)\(\s*(\d+)\)=\s*([-+]?[\d.]+E[-+]?\d+|[-+]?[\d.]+)')

# The legacy "unused" sentinel; these entries carry no value.
UNUSED = 1.0e-30


def extract(path: Path, name: str, occurrence: int = 0) -> Dict[int, float]:
    """Read one COMMON block dump from a listing.

    A listing repeats the dump once per case or condition.  ``occurrence``
    selects which one, defaulting to the first.

    Args:
        path: Output listing to read.
        name: COMMON array name, for example ``'A'`` or ``'BODY'``.
        occurrence: Which printing of the block to return.

    Returns:
        Mapping of one-based index to value.

    Raises:
        ValueError: If the named block never appears, or the requested
            occurrence does not exist.
    """
    blocks: List[Dict[int, float]] = []
    current: Dict[int, float] = {}
    previous_index = None

    with open(path, errors='replace') as handle:
        for line in handle:
            found = [(int(m.group(2)), float(m.group(3)))
                     for m in _RECORD.finditer(line)
                     if m.group(1) == name]
            if not found:
                continue
            # A new block starts whenever the index sequence restarts.
            if previous_index is not None and found[0][0] <= previous_index:
                blocks.append(current)
                current = {}
            for index, value in found:
                current[index] = value
            previous_index = found[-1][0]

    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"COMMON block {name!r} not found in {path}")
    if occurrence >= len(blocks):
        raise ValueError(
            f"{path} has {len(blocks)} printing(s) of {name!r}; "
            f"occurrence {occurrence} requested")
    return blocks[occurrence]


def is_set(value: float) -> bool:
    """Whether a dumped entry holds a real value rather than the sentinel."""
    return abs(value - UNUSED) > 1e-40 and value != 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('listing', type=Path)
    parser.add_argument('name', help="COMMON array name, e.g. A")
    parser.add_argument('--occurrence', type=int, default=0)
    parser.add_argument('--index', type=int, default=None)
    args = parser.parse_args()

    try:
        block = extract(args.listing, args.name, args.occurrence)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    if args.index is not None:
        value = block.get(args.index)
        if value is None:
            print(f"{args.name}({args.index}) not present", file=sys.stderr)
            return 1
        print(f"{args.name}({args.index}) = {value:.6E}")
        return 0

    print(f"{args.name}: {len(block)} entries "
          f"({sum(1 for v in block.values() if is_set(v))} set)")
    for index in sorted(block):
        value = block[index]
        note = '' if is_set(value) else '   (unused/zero)'
        print(f"  {args.name}({index:>3}) = {value:>14.6E}{note}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
