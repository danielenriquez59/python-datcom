"""
Parser for DATCOM namelist input files (.inp format).

Parses the FORTRAN-style namelist format used by Digital DATCOM:
  $NAMELIST param=value, param=value, ... $

Reference: datcom.f lines 6887-6895 and INPUT subroutine
"""

import re
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class NamelistParser:
    """
    Parse DATCOM .inp format namelist files.
    
    Handles:
    - $NAMELIST ... $ syntax
    - Array notation: PARAM=val1, val2 or PARAM(1)=val1, val2
    - Multi-case input with CASEID, SAVE, NEXT CASE
    - Comments and formatting
    """
    
    # Known namelists from DATCOM
    KNOWN_NAMELISTS = {
        'FLTCON', 'OPTINS', 'SYNTHS', 'BODY', 'WGPLNF', 'WGSCHR',
        'HTPLNF', 'HTSCHR', 'VTPLNF', 'VTSCHR', 'VFPLNF', 'VFSCHR',
        'SYMFLP', 'ASYFLP', 'DEFLCT', 'GROUND', 'TRIM', 'DAMP',
        'PART', 'DERIV', 'DUMP', 'BUILD', 'PWRINP', 'JET',
        'HYPER', 'PROPWR', 'JETPWR', 'NACON', 'CONTAB', 'EXPDATA'
    }
    
    # Commands
    COMMANDS = {'CASEID', 'SAVE', 'NEXT CASE', 'DUMP', 'BUILD', 'DIM', 'DERI', 'PART'}
    
    def __init__(self):
        """Initialize parser."""
        self.cases: List[Dict[str, Any]] = []
        self.current_case: Dict[str, Any] = {'namelists': {}, 'commands': []}
    
    def parse_file(self, filepath: Path) -> List[Dict[str, Any]]:
        """
        Parse DATCOM input file.
        
        Args:
            filepath: Path to .inp file
            
        Returns:
            List of case dictionaries
        """
        with open(filepath, 'r') as f:
            content = f.read()
        
        return self.parse(content)
    
    def parse(self, content: str) -> List[Dict[str, Any]]:
        """
        Parse DATCOM input content.
        
        Args:
            content: Input file content string
            
        Returns:
            List of case dictionaries, each containing namelists and commands
        """
        self.cases = []
        self.current_case = {'namelists': {}, 'commands': [], 'caseid': ''}
        
        lines = content.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('!') or line.startswith('C'):
                i += 1
                continue
            
            # Check for CASEID
            if line.upper().startswith('CASEID'):
                caseid = line[6:].strip()
                self.current_case['caseid'] = caseid
                i += 1
                continue
            
            # Check for SAVE command
            if line.upper() == 'SAVE':
                self.current_case['commands'].append('SAVE')
                i += 1
                continue
            
            # Check for NEXT CASE command
            if line.upper().startswith('NEXT'):
                if self.current_case['namelists'] or self.current_case['commands']:
                    self.cases.append(self.current_case)
                self.current_case = {'namelists': {}, 'commands': [], 'caseid': ''}
                i += 1
                continue
            
            # Check for DUMP command
            if line.upper().startswith('DUMP'):
                dump_what = line[4:].strip()
                self.current_case['commands'].append(f'DUMP {dump_what}')
                i += 1
                continue
            
            # Check for namelist start
            if line.startswith('$'):
                namelist_content, i = self._extract_namelist(lines, i)
                if namelist_content:
                    namelist_name, params = self._parse_namelist(namelist_content)
                    if namelist_name:
                        # Merge with existing namelist if same name appears multiple times
                        if namelist_name in self.current_case['namelists']:
                            existing = self.current_case['namelists'][namelist_name]
                            for key, value in params.items():
                                if key in existing and (isinstance(existing[key], list) or isinstance(value, list)):
                                    if not isinstance(existing[key], list):
                                        existing[key] = [existing[key]]
                                    if not isinstance(value, list):
                                        value = [value]
                                    # Merge arrays by extending
                                    max_len = max(len(existing[key]), len(value))
                                    # Ensure lists are same length
                                    while len(existing[key]) < max_len:
                                        existing[key].append(None)
                                    # Update with new values
                                    for idx, v in enumerate(value):
                                        if v is not None:
                                            existing[key][idx] = v
                                else:
                                    existing[key] = value
                        else:
                            self.current_case['namelists'][namelist_name] = params
                continue
            
            i += 1
        
        # Add final case if exists
        if self.current_case['namelists'] or self.current_case['commands']:
            self.cases.append(self.current_case)
        
        return self.cases
    
    def _extract_namelist(self, lines: List[str], start_idx: int) -> Tuple[str, int]:
        """
        Extract complete namelist from lines.
        
        Args:
            lines: Input lines
            start_idx: Starting line index
            
        Returns:
            Tuple of (namelist_content, next_line_index)
        """
        namelist_lines = []
        i = start_idx
        
        # Get first line
        line = lines[i].strip()
        namelist_lines.append(line)
        
        # Check if namelist is complete on one line
        if line.count('$') >= 2:
            return ' '.join(namelist_lines), i + 1
        
        # Continue reading until we find the closing $
        i += 1
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('!'):
                i += 1
                continue
            
            namelist_lines.append(line)
            
            if '$' in line:
                return ' '.join(namelist_lines), i + 1
            
            i += 1
        
        return ' '.join(namelist_lines), i
    
    def _parse_namelist(self, content: str) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        Parse namelist content.
        
        Args:
            content: Namelist string
            
        Returns:
            Tuple of (namelist_name, parameters_dict)
        """
        # Extract namelist name and content
        match = re.match(r'\$(\w+)\s+(.*?)\$', content, re.DOTALL)
        if not match:
            logger.warning(f"Could not parse namelist: {content[:50]}...")
            return None, {}
        
        name = match.group(1).upper()
        param_str = match.group(2)
        
        # Parse parameters
        params = self._parse_parameters(param_str)
        
        return name, params
    
    def _parse_parameters(self, param_str: str) -> Dict[str, Any]:
        """
        Parse parameter assignments from namelist.
        
        Args:
            param_str: Parameter string
            
        Returns:
            Dictionary of parameter name -> value
        """
        params = {}
        
        # Remove extra whitespace but preserve structure
        param_str = re.sub(r'\s+', ' ', param_str).strip()
        
        # Split into individual assignments
        assignments = self._split_assignments(param_str)
        
        # Unindexed lists start at element 1, just as indexed lists do.
        # Keep a lone unindexed value scalar until a continuation establishes
        # that it is an array. INPUT's MACH/ALSCHD/ALT are examples of arrays
        # whose legacy input decks commonly omit the explicit (1).
        last_array_name = None
        last_array_idx = None
        
        for assignment in assignments:
            assignment = assignment.strip()
            if not assignment:
                continue
            
            if '=' in assignment:
                # New assignment
                parts = assignment.split('=', 1)
                param_name = parts[0].strip().upper()
                value_str = parts[1].strip()
                
                # Parse value
                value = self._parse_value(param_name, value_str)
                
                # Handle array indexing
                if '(' in param_name:
                    array_name, index_str = param_name.split('(', 1)
                    index_str = index_str.rstrip(')')
                    start_idx = int(index_str) - 1  # Convert to 0-based
                    
                    # Initialize array if needed
                    if array_name not in params:
                        params[array_name] = []
                    elif not isinstance(params[array_name], list):
                        params[array_name] = [params[array_name]]
                    
                    # Ensure array is large enough
                    if isinstance(value, list):
                        # Multiple values starting at index
                        for i, v in enumerate(value):
                            idx = start_idx + i
                            while len(params[array_name]) <= idx:
                                params[array_name].append(None)
                            params[array_name][idx] = v
                        last_array_idx = start_idx + len(value)
                    else:
                        # Single value at index
                        while len(params[array_name]) <= start_idx:
                            params[array_name].append(None)
                        params[array_name][start_idx] = value
                        last_array_idx = start_idx + 1
                    
                    last_array_name = array_name
                else:
                    if isinstance(params.get(param_name), list):
                        values = value if isinstance(value, list) else [value]
                        for idx, v in enumerate(values):
                            while len(params[param_name]) <= idx:
                                params[param_name].append(None)
                            params[param_name][idx] = v
                    else:
                        params[param_name] = value
                    last_array_name = param_name
                    last_array_idx = len(value) if isinstance(value, list) else 1
            else:
                # Continuation values for last array
                if last_array_name and last_array_idx is not None:
                    if not isinstance(params[last_array_name], list):
                        params[last_array_name] = [params[last_array_name]]
                    cont_values = self._parse_value('', assignment)
                    if isinstance(cont_values, list):
                        for v in cont_values:
                            idx = last_array_idx
                            while len(params[last_array_name]) <= idx:
                                params[last_array_name].append(None)
                            params[last_array_name][idx] = v
                            last_array_idx += 1
                    else:
                        idx = last_array_idx
                        while len(params[last_array_name]) <= idx:
                            params[last_array_name].append(None)
                        params[last_array_name][idx] = cont_values
                        last_array_idx += 1
        
        return params
    
    def _split_assignments(self, param_str: str) -> List[str]:
        """
        Split parameter string into individual assignments.
        
        Args:
            param_str: Parameter string
            
        Returns:
            List of assignment strings
        """
        assignments = []
        current = []
        depth = 0
        
        for char in param_str:
            if char == '(':
                depth += 1
                current.append(char)
            elif char == ')':
                depth -= 1
                current.append(char)
            elif char == ',' and depth == 0:
                if current:
                    assignments.append(''.join(current).strip())
                current = []
            else:
                current.append(char)
        
        if current:
            assignments.append(''.join(current).strip())
        
        return assignments
    
    def _parse_value(self, param_name: str, value_str: str) -> Any:
        """
        Parse parameter value with type inference.
        
        Args:
            param_name: Parameter name
            value_str: Value string
            
        Returns:
            Parsed value (int, float, bool, str, or list)
        """
        value_str = value_str.strip()
        
        # Handle multiple values (comma-separated)
        if ',' in value_str:
            values = [self._parse_single_value(v.strip()) for v in value_str.split(',') if v.strip()]
            return values
        
        return self._parse_single_value(value_str)
    
    def _parse_single_value(self, value_str: str) -> Any:
        """
        Parse a single value with type inference.
        
        Args:
            value_str: Value string
            
        Returns:
            Parsed value
        """
        value_str = value_str.strip()
        
        # Boolean
        if value_str.upper() in ('.TRUE.', 'TRUE'):
            return True
        if value_str.upper() in ('.FALSE.', 'FALSE'):
            return False
        
        # Try integer
        try:
            return int(value_str)
        except ValueError:
            pass
        
        # Try float (including scientific notation)
        try:
            # Handle FORTRAN exponential notation (e.g., 4.28E6)
            value_str = value_str.replace('E', 'e').replace('D', 'e')
            return float(value_str)
        except ValueError:
            pass
        
        # String (remove quotes if present)
        if value_str.startswith(("'", '"')) and value_str.endswith(("'", '"')):
            return value_str[1:-1]
        
        return value_str
    
    def to_state_dict(self, case: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert parsed case to state dictionary format.
        
        Args:
            case: Parsed case dictionary
            
        Returns:
            State dictionary with component-prefixed keys
        """
        state = {}
        
        for namelist_name, params in case['namelists'].items():
            prefix = self._get_prefix_for_namelist(namelist_name)
            
            for param_name, value in params.items():
                key = f"{prefix}_{param_name.lower()}"
                state[key] = value
        
        # Add case metadata
        state['case_id'] = case.get('caseid', '')
        state['case_save'] = 'SAVE' in case.get('commands', [])
        
        return state
    
    def _get_prefix_for_namelist(self, namelist: str) -> str:
        """
        Get state dictionary prefix for namelist.
        
        Args:
            namelist: Namelist name
            
        Returns:
            Prefix string
        """
        mapping = {
            'FLTCON': 'flight',
            'OPTINS': 'options',
            'SYNTHS': 'synths',
            'BODY': 'body',
            'WGPLNF': 'wing',
            'WGSCHR': 'wing',
            'HTPLNF': 'htail',
            'HTSCHR': 'htail',
            'VTPLNF': 'vtail',
            'VTSCHR': 'vtail',
            'VFPLNF': 'vfin',
            'VFSCHR': 'vfin',
        }
        
        return mapping.get(namelist, namelist.lower())

