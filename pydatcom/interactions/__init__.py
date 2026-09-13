"""Aerodynamic interaction calculations - wing-body, wing-tail, downwash, power effects."""

from pydatcom.interactions.carryover import (
    body_semispan_ratio, fig4312_10, fig4312_12a,
    calculate_carryover_factors,
)

__all__ = [
    'body_semispan_ratio',
    'fig4312_10',
    'fig4312_12a',
    'calculate_carryover_factors',
]
