"""
Global state management for DATCOM calculations.

This module manages the global state dictionary that replaces FORTRAN COMMON blocks.
All state variables use component-prefixed naming for organization and future refactoring.

Reference: Multiple COMMON blocks throughout datcom.f
"""

from typing import Any, Dict, Optional, List
import yaml
from pathlib import Path


class StateManager:
    """
    Manages global state dictionary with component-prefixed keys.
    
    Replaces FORTRAN COMMON blocks with a centralized state dictionary.
    Component prefixes organize variables by functional area for future refactoring.
    """
    
    def __init__(self):
        """Initialize state manager with default values."""
        self._state: Dict[str, Any] = {}
        self._initialize_defaults()
    
    def _initialize_defaults(self) -> None:
        """Initialize default values from constants and COMMON blocks."""
        # Import here to avoid circular imports
        from pydatcom.utils.constants import get_constants_dict
        
        # Constants (COMMON /CONSNT/)
        self._state.update(get_constants_dict())
        
        # Flight conditions (COMMON /FLGTCD/)
        self._state['flight_nmach'] = 1
        self._state['flight_mach'] = []
        self._state['flight_nalpha'] = 0
        self._state['flight_alpha'] = []
        self._state['flight_rnnub'] = []
        self._state['flight_alt'] = []
        self._state['flight_vinf'] = []
        self._state['flight_pinf'] = []
        self._state['flight_tinf'] = []
        
        # Options (COMMON /OPTION/)
        self._state['options_sref'] = None
        self._state['options_cbarr'] = None
        self._state['options_rougfc'] = 1.6e-4
        self._state['options_blref'] = None
        self._state['options_irun'] = 0
        
        # Synthesis parameters (COMMON /SYNTSS/)
        self._state['synths_xcg'] = None
        self._state['synths_xw'] = None
        self._state['synths_zw'] = None
        self._state['synths_aliw'] = None
        self._state['synths_zcg'] = None
        self._state['synths_xh'] = None
        self._state['synths_zh'] = None
        self._state['synths_alih'] = None
        self._state['synths_xv'] = None
        self._state['synths_zv'] = None
        self._state['synths_xvf'] = None
        self._state['synths_zvf'] = None
        self._state['synths_yv'] = None
        self._state['synths_yf'] = None
        self._state['synths_phiv'] = None
        self._state['synths_phif'] = None
        self._state['synths_vertup'] = False
        self._state['synths_hinax'] = None
        self._state['synths_scale'] = None
        
        # Body geometry (COMMON /BODYI/, /IBODY/, /BDATA/)
        self._state['body_nx'] = 0
        self._state['body_x'] = []
        self._state['body_s'] = []
        self._state['body_p'] = []
        self._state['body_r'] = []
        self._state['body_zu'] = []
        self._state['body_zl'] = []
        self._state['body_bnose'] = None
        self._state['body_btail'] = None
        self._state['body_bln'] = None
        self._state['body_bla'] = None
        self._state['body_ds'] = None
        self._state['body_itype'] = 2
        self._state['body_method'] = 1
        
        # Wing data (COMMON /WINGI/, /IWING/, /WINGD/)
        self._state['wing_data'] = {}
        self._state['wing_a'] = [0.0] * 195
        self._state['wing_b'] = [0.0] * 49
        
        # Horizontal tail (COMMON /HTI/, /IHT/, /HTDATA/)
        self._state['htail_data'] = {}
        self._state['htail_a'] = [0.0] * 195
        self._state['htail_b'] = [0.0] * 49
        
        # Vertical tail (COMMON /VTI/, /IVT/, /VTDATA/)
        self._state['vtail_data'] = {}
        self._state['vtail_a'] = [0.0] * 195
        self._state['vtail_vf'] = [0.0] * 195
        
        # Aerodynamic outputs (COMMON /WHAERO/)
        self._state['aero_cl'] = None
        self._state['aero_cd'] = None
        self._state['aero_cm'] = None
        self._state['aero_cn'] = None
        self._state['aero_ca'] = None
        
        # Control flags (COMMON /FLOLOG/)
        self._state['flags_fltc'] = False
        self._state['flags_opti'] = False
        self._state['flags_bo'] = False
        self._state['flags_wgpl'] = False
        self._state['flags_wgsc'] = False
        self._state['flags_synt'] = False
        self._state['flags_htpl'] = False
        self._state['flags_htsc'] = False
        self._state['flags_vtpl'] = False
        self._state['flags_vtsc'] = False
        self._state['flags_supers'] = False
        self._state['flags_subson'] = False
        self._state['flags_transn'] = False
        self._state['flags_hypers'] = False
        
        # Case control
        self._state['case_id'] = ""
        self._state['case_save'] = False
        self._state['case_dump'] = False
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get value from state dictionary.
        
        Args:
            key: State variable key
            default: Default value if key not found
            
        Returns:
            Value from state or default
        """
        return self._state.get(key, default)
    
    def set(self, key: str, value: Any) -> None:
        """
        Set value in state dictionary.
        
        Args:
            key: State variable key
            value: Value to set
        """
        self._state[key] = value
    
    def update(self, data: Dict[str, Any]) -> None:
        """
        Update multiple state values.
        
        Args:
            data: Dictionary of key-value pairs to update
        """
        self._state.update(data)
    
    def reset(self, keep_constants: bool = True) -> None:
        """
        Reset state to defaults.
        
        Args:
            keep_constants: If True, preserve constant values
        """
        if keep_constants:
            constants = {k: v for k, v in self._state.items() if k.startswith('constants_')}
            self._initialize_defaults()
            self._state.update(constants)
        else:
            self._initialize_defaults()
    
    def export_to_yaml(self, filepath: Path) -> None:
        """
        Export current state to YAML file.
        
        Args:
            filepath: Path to output YAML file
        """
        with open(filepath, 'w') as f:
            yaml.dump(self._state, f, default_flow_style=False, sort_keys=True)
    
    def import_from_yaml(self, filepath: Path) -> None:
        """
        Import state from YAML file.
        
        Args:
            filepath: Path to input YAML file
        """
        with open(filepath, 'r') as f:
            data = yaml.safe_load(f)
            self._state.update(data)
    
    def get_all(self) -> Dict[str, Any]:
        """
        Get entire state dictionary.
        
        Returns:
            Complete state dictionary
        """
        return self._state.copy()
    
    def get_component(self, prefix: str) -> Dict[str, Any]:
        """
        Get all state variables for a component.
        
        Args:
            prefix: Component prefix (e.g., 'wing', 'body', 'flight')
            
        Returns:
            Dictionary of component variables with prefix removed from keys
        """
        prefix_key = f"{prefix}_"
        return {
            k.replace(prefix_key, ''): v 
            for k, v in self._state.items() 
            if k.startswith(prefix_key)
        }

