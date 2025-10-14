# PyDATCOM

A Python implementation of the USAF Digital DATCOM (Data Compendium) aerodynamic analysis tool.

## Overview

PyDATCOM is a modernized Python conversion of the FORTRAN-based Digital DATCOM program. It provides aerodynamic analysis capabilities for aircraft and missiles across subsonic, transonic, supersonic, and hypersonic flight regimes.

**Original Reference**: AFFDL-TR-79-3032 (NTIS ADA-086557)

## Features

- **Input Parsing**: Compatible with original DATCOM namelist format (.inp files)
- **State Management**: Component-based state dictionary system
- **Modular Architecture**: Organized by functional areas (geometry, aerodynamics, interactions)
- **Modern Python**: Type hints, logging, exception handling
- **Multi-format Output**: DATCOM, YAML, JSON, CSV support

## Installation

### Requirements

- Python 3.8+
- NumPy >= 1.19.0
- PyYAML >= 5.3.0

### Install Dependencies

```bash
pip install -r requirements.txt
```

## Project Structure

```
pydatcom/
├── config/          # Configuration and defaults
├── io/              # Input/output and state management
├── geometry/        # Airfoil, body, wing, tail geometry
├── aerodynamics/    # Lift, drag, moment calculations
├── interactions/    # Wing-body, downwash, power effects
├── utils/           # Math, atmosphere, interpolation utilities
└── core/            # Main engine and case runner
```

## Usage

### Parse Input File

```python
from pydatcom.io import NamelistParser, StateManager

# Parse DATCOM input file
parser = NamelistParser()
cases = parser.parse_file('input.inp')

# Convert to state dictionary
state = StateManager()
state.update(parser.to_state_dict(cases[0]))

# Access parameters
mach = state.get('flight_mach')
wing_area = state.get('options_sref')
```

### State Dictionary Naming Convention

All state variables use component prefixes for organization:

- `constants_*` - Physical/math constants (PI, DEG, RAD)
- `flight_*` - Flight conditions (Mach, alpha, Reynolds)
- `options_*` - Reference dimensions (SREF, CBARR, BLREF)
- `synths_*` - Synthesis parameters (XCG, XW, ZW)
- `body_*` - Body geometry
- `wing_*` - Wing data
- `htail_*` - Horizontal tail
- `vtail_*` - Vertical tail
- `aero_*` - Aerodynamic outputs (CL, CD, CM)
- `flags_*` - Control flags

## Development Status

### Phase 1: Core Infrastructure ✅
- [x] Namelist parser
- [x] State manager
- [x] Constants module
- [x] Configuration system

### Phase 2: Geometry (In Progress)
- [ ] Airfoil coordinate generation
- [ ] Body geometry
- [ ] Wing geometry
- [ ] Tail geometry

### Phase 3: Utilities (In Progress)
- [x] Math utilities (ARCSIN, ARCCOS, AREA)
- [x] Atmosphere (US Standard 1962)
- [ ] Interpolation functions
- [ ] Table lookup

### Phase 4-7: Aerodynamics & Integration (Planned)
- [ ] Lift calculations
- [ ] Drag calculations
- [ ] Moment calculations
- [ ] Wing-body interactions
- [ ] Main computational engine
- [ ] Output formatting

## Testing

Run basic tests:

```bash
python tests/test_parser.py
```

## FORTRAN to Python Conversion Notes

### Key Differences

1. **Array Indexing**: FORTRAN 1-based → Python 0-based
2. **COMMON Blocks**: → State dictionary with prefixed keys
3. **GOTO Statements**: → Structured if/else/while
4. **DATA Statements**: → YAML/JSON tables
5. **Type System**: Implicit typing → Explicit type hints

### Example Conversion

FORTRAN:
```fortran
COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
SUBROUTINE ARCSIN(A)
    IF (A .GT. 1.0) A = 1.0
    ARCSIN = ASIN(A)
    RETURN
END
```

Python:
```python
def arcsin(a: float) -> float:
    """Arc sine with bounds checking (FORTRAN ARCSIN, line 434)."""
    if a > 1.0:
        a = 1.0
    return np.arcsin(a)
```

## References

- USAF Stability and Control DATCOM (1978)
- AFFDL-TR-79-3032: Digital DATCOM User's Manual
- Public Domain Aeronautical Software (PDAS)

## License

Converted from public domain USAF Digital DATCOM source code.

## Disclaimer

THIS SOFTWARE IS RELEASED "AS IS". NO WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
IS PROVIDED CONCERNING THIS SOFTWARE, INCLUDING WARRANTIES OF MERCHANTABILITY 
OR FITNESS FOR A PARTICULAR PURPOSE.

