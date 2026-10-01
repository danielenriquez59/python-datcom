"""Make the Fortran parity helpers in ``test_parity/`` importable by the tests."""
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(_ROOT / 'test_parity'))

LEGACY = _ROOT / 'datcom-legacy' / 'datcom_2000'

# Marks a test that reads the DATCOM Fortran source, which is not distributed.
requires_fortran = pytest.mark.skipif(
    not LEGACY.is_dir(), reason="datcom-legacy/ Fortran source not present")
