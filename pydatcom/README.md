# PyDATCOM

A modern Python implementation of the USAF Digital DATCOM aerodynamic analysis tool.

## Purpose

PyDATCOM translates the legacy FORTRAN-based Digital DATCOM (Data Compendium) program into contemporary Python. The original DATCOM provides comprehensive aerodynamic analysis for aircraft and missiles across all flight regimes, but its FORTRAN codebase presents maintenance and usability challenges.

This Python conversion maintains full compatibility with existing DATCOM input formats while providing:
- **Modern Python architecture** with type hints, exception handling, and logging
- **Modular design** for easier maintenance and extension
- **NumPy/SciPy integration** for improved performance and numerical accuracy
- **Enhanced usability** through Python's ecosystem and tooling

## What is DATCOM?

The USAF Digital DATCOM (AFFDL-TR-79-3032) is the standard engineering tool for estimating aerodynamic characteristics of aircraft configurations. It provides empirical correlations and analytical methods for:

- **Lift, drag, and moment coefficients** across subsonic, transonic, supersonic, and hypersonic regimes
- **Stability derivatives** for flight control analysis
- **Component buildup methods** for wing-body-tail configurations
- **Power effects** and high-lift device analysis

Originally developed in the 1970s, DATCOM remains the industry standard for preliminary aerodynamic design due to its extensive validation against experimental data.

## Python Architecture

The PyDATCOM codebase is organized into logical modules that mirror DATCOM's functional areas:

### Core Components

- **`config/`** - Configuration management and default parameters
- **`io/`** - Input parsing and state management systems
- **`geometry/`** - Geometric definitions for airfoils, bodies, wings, and tails
- **`aerodynamics/`** - Aerodynamic coefficient calculations by flight regime
- **`interactions/`** - Component interference and coupling effects
- **`utils/`** - Mathematical utilities, atmospheric models, and interpolation functions

### Key Design Principles

1. **State Dictionary Pattern**: All data flows through a centralized state dictionary with prefixed keys for organization
2. **Functional Decomposition**: Complex FORTRAN subroutines are broken into focused, testable functions
3. **Type Safety**: Full type hints throughout for better IDE support and error catching
4. **Exception Handling**: Graceful error handling with informative messages
5. **Logging**: Comprehensive logging for debugging and monitoring

### Input/Output Compatibility

PyDATCOM maintains full compatibility with the original DATCOM namelist format (.inp files) while adding support for modern formats:

- **Input**: Traditional DATCOM namelist format
- **Output**: DATCOM format, YAML, JSON, CSV

## Installation & Requirements

### Requirements
- Python 3.8+
- NumPy >= 1.19.0
- PyYAML >= 5.3.0

### Install

From the repository root:

```bash
pip install -e ".[dev]"
```

## Usage Example

```python
from pydatcom.io import NamelistParser, StateManager

# Parse traditional DATCOM input file
parser = NamelistParser()
cases = parser.parse_file('aircraft.inp')

# Initialize state management
state = StateManager()
state.update(parser.to_state_dict(cases[0]))

# Access aerodynamic parameters
mach = state.get('flight_mach')
wing_area = state.get('options_sref')
lift_coeff = state.get('aero_cl')
```

## State Dictionary Convention

All state variables use descriptive prefixes for organization:

| Prefix | Description |
|--------|-------------|
| `constants_*` | Physical and mathematical constants |
| `flight_*` | Flight conditions (Mach, alpha, Reynolds) |
| `options_*` | Reference dimensions and configuration |
| `synths_*` | Synthesis parameters (CG locations, etc.) |
| `body_*` | Body geometry parameters |
| `wing_*` | Wing geometry and characteristics |
| `htail_*` | Horizontal tail parameters |
| `vtail_*` | Vertical tail parameters |
| `aero_*` | Aerodynamic coefficient outputs |
| `flags_*` | Control and option flags |

## FORTRAN to Python Translation

The conversion addresses key differences between FORTRAN and modern Python:

### Architectural Changes
- **COMMON blocks** → Centralized state dictionary
- **GOTO statements** → Structured control flow
- **Implicit typing** → Explicit type hints
- **Fixed arrays** → Dynamic NumPy arrays
- **Subroutine calls** → Object-oriented methods

### Performance & Accuracy
- **Vectorization**: Array operations replace element-wise loops
- **Numerical libraries**: NumPy/SciPy replace custom math routines
- **Precision**: Consistent floating-point handling
- **Memory management**: Automatic garbage collection

## References

- **USAF Stability and Control DATCOM** (1978) - Original technical reference
- **AFFDL-TR-79-3032** - Digital DATCOM User's Manual
- **Public Domain Aeronautical Software** (PDAS) - Additional validation data

## License

Converted from public domain USAF Digital DATCOM source code. This Python implementation is released under the same public domain terms as the original FORTRAN code.

## Disclaimer

THIS SOFTWARE IS RELEASED "AS IS". NO WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, IS PROVIDED CONCERNING THIS SOFTWARE, INCLUDING WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

