"""
PyDATCOM - Python Digital DATCOM
A Python implementation of the USAF Digital DATCOM aerodynamic analysis tool.

Converted from FORTRAN source code (AFFDL-TR-79-3032).
"""

__version__ = "0.1.0"
__author__ = "Converted from USAF Digital DATCOM"

from pydatcom.io.state_manager import StateManager
from pydatcom.io.namelist_parser import NamelistParser

__all__ = ['StateManager', 'NamelistParser']

