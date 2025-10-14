"""
US Standard Atmosphere 1962 calculations.

Reference: FORTRAN ATMOS subroutine, datcom.f line 621
"""

import numpy as np
from typing import Dict, Tuple
import logging

logger = logging.getLogger(__name__)


class Atmosphere:
    """
    US Standard Atmosphere 1962 model with inverse square gravitational field.
    
    Agrees with COESA document within 1% up to 700 km (2,296,588 ft).
    Reference: datcom.f line 621
    """
    
    # Constants
    G0 = 32.1740484  # Sea level gravitational acceleration, ft/sec²
    WM0 = 28.9644    # Sea level molecular weight
    R0 = 20890855.0  # Earth radius, ft
    GMRS = 0.018743418  # G0*M0/R*, deg R/ft
    
    # Geopotential altitude breakpoints (ft)
    HG = np.array([
        -16404., 0.0, 36089., 65617., 104987., 154199., 
        170604., 200131., 250186., 291160.
    ])
    
    # Geometric altitude breakpoints (ft)  
    ZM = np.array([
        295276., 328084., 360892., 393701., 492126., 524934.,
        557743., 623360., 754593., 984252., 1312336., 1640420.,
        1968504., 2296588.
    ])
    
    # Molecular weights
    WM = np.array([
        28.9644, 28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85,
        24.7, 22.66, 19.94, 17.94, 16.84, 16.17
    ])
    
    # Temperature at breakpoints (deg R)
    TM = np.array([
        577.17, 518.67, 389.97, 389.97, 411.57, 487.17, 487.17,
        454.77, 325.17, 325.17, 379.17, 469.17, 649.17, 1729.17,
        1999.17, 2179.17, 2431.17, 2791.17, 3295.17, 3889.17,
        4357.17, 4663.17, 4861.17
    ])
    
    # Pressure at breakpoints (lb/ft²)
    PM = np.array([
        3711.0839, 2116.2165, 472.67563, 114.34314, 18.128355,
        2.3162178, 1.2321972, 3.8030279E-01, 2.1671352E-02,
        3.4313478E-03, 6.2773411E-04, 1.5349091E-04, 5.2624212E-05,
        1.0561806E-05, 7.7083076E-06, 5.8267151E-06, 3.5159854E-06,
        1.4520255E-06, 3.9290563E-07, 8.4030242E-08, 2.2835256E-08,
        7.1875452E-09
    ])
    
    @classmethod
    def calculate(cls, altitude: float) -> Dict[str, float]:
        """
        Calculate atmospheric properties at given altitude.
        
        Args:
            altitude: Geometric altitude in feet
            
        Returns:
            Dictionary with atmospheric properties:
            - cs: Speed of sound, ft/sec
            - dcs_dz: Sound derivative, 1/ft  
            - altitude: Geometric altitude, ft
            - pressure: Pressure, lb/ft²
            - dp_dz: Pressure derivative, lb/ft³
            - density: Density, slugs/ft³
            - drho_dz: Density derivative, 1/ft
            - temperature: Temperature, deg Rankine
            - dt_dz: Temperature derivative, deg R/ft
        """
        z = altitude
        
        # Calculate gravity at altitude
        g = cls.G0 * (cls.R0 / (cls.R0 + z))**2
        
        if z <= 295276.0:
            # Temperature linear with geopotential
            h = cls.R0 * z / (cls.R0 + z)
            
            # Find temperature region
            j = 0
            for i in range(1, len(cls.HG)):
                if cls.HG[i] >= h:
                    j = i - 1
                    break
            
            # Calculate temperature slope and value
            if j < len(cls.HG) - 1:
                elh = (cls.TM[j + 1] - cls.TM[j]) / (cls.HG[j + 1] - cls.HG[j])
            else:
                elh = 0.0
            
            tms = cls.TM[j] + elh * (h - cls.HG[j])
            elz = elh * g / cls.G0
            dmdz = 0.0
            em = cls.WM0
            
            # Calculate pressure
            if elh != 0.0:
                # Non-zero slope
                pressure = cls.PM[j] * (cls.TM[j] / tms)**(cls.GMRS / elh)
            else:
                # Zero slope (isothermal)
                pressure = cls.PM[j] * np.exp(cls.GMRS * (cls.HG[j] - h) / tms)
        
        else:
            # Temperature linear with Z (high altitude)
            j = 8
            k = 0
            for i in range(1, len(cls.ZM)):
                if cls.ZM[i] >= z:
                    j = i + 8
                    k = i - 1
                    break
            
            # Calculate temperature slope and value
            if k < len(cls.ZM) - 1:
                elz = (cls.TM[j + 1] - cls.TM[j]) / (cls.ZM[k + 1] - cls.ZM[k])
            else:
                elz = 0.0
            
            tms = cls.TM[j] + elz * (z - cls.ZM[k])
            
            if k < len(cls.WM) - 1:
                dmdz = (cls.WM[k + 1] - cls.WM[k]) / (cls.ZM[k + 1] - cls.ZM[k])
            else:
                dmdz = 0.0
            
            em = cls.WM[k] + dmdz * (z - cls.ZM[k])
            zlz = z - tms / elz if elz != 0.0 else z
            
            # Pressure equation for high altitude
            if elz != 0.0:
                exp_term = cls.GMRS / elz * (cls.R0 / (cls.R0 + zlz))**2 * (
                    (z - cls.ZM[k]) * (cls.R0 + zlz) / (cls.R0 + z) / (cls.R0 + cls.ZM[k]) -
                    np.log(tms * (cls.R0 + cls.ZM[k]) / cls.TM[j] / (cls.R0 + z))
                )
                pressure = cls.PM[j] * np.exp(exp_term)
            else:
                pressure = cls.PM[j]
        
        # Calculate speed of sound and derivative
        cs = 49.022164 * np.sqrt(tms)
        dcs_dz = 0.5 * elz / tms
        
        # Calculate density and derivatives
        density = cls.GMRS * pressure / cls.G0 / tms
        drho_dz = -(density * g / pressure + elz / tms)
        dp_dz = -density * g
        
        # Calculate temperature and derivative
        temperature = em * tms / cls.WM0
        dt_dz = (em * elz + tms * dmdz) / cls.WM0
        
        return {
            'cs': cs,
            'dcs_dz': dcs_dz,
            'altitude': z,
            'pressure': pressure,
            'dp_dz': dp_dz,
            'density': density,
            'drho_dz': drho_dz,
            'temperature': temperature,
            'dt_dz': dt_dz,
        }
    
    @classmethod
    def get_properties(cls, altitude: float) -> Tuple[float, float, float]:
        """
        Get basic atmospheric properties (convenience method).
        
        Args:
            altitude: Geometric altitude in feet
            
        Returns:
            Tuple of (temperature_R, pressure_psf, density_slug_ft3)
        """
        atm = cls.calculate(altitude)
        return atm['temperature'], atm['pressure'], atm['density']

