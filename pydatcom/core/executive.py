"""
Executive sequencing: MAIN05, MAIN06 and MAIN07, which choose the flap and
control overlays for a case.

Each takes the configuration flags and a ``call(name)`` callable that runs
the named overlay, and returns the overlays called in order.

Reference: datcom-legacy/datcom_2000/main05.f, main06.f, main07.f
"""

from typing import Callable, List, MutableSequence


def main05(symfp: bool, asyfp: bool, trim: bool, htpl: bool,
           f: MutableSequence[float],
           call: Callable[[str], object]) -> List[str]:
    """MAIN05: the subsonic flap and control overlays.

    ``f`` is ``/FLAPIN/`` as a 1-based list, edited in place: for a
    symmetric flap without a jet flap one more deflection is added at
    zero for the overlays and the count restored afterwards, but the zero
    is left in the added slot.
    """
    ran: List[str] = []

    def run(name):
        ran.append(name)
        call(name)

    if symfp:
        if int(f[74] + 0.5) >= 1:
            run('M55O67')
        else:
            f[16] = f[16] + 1.0
            f[int(f[16] + 0.5)] = 0.0
            for name in ('M36O44', 'M37O45', 'M38O46'):
                run(name)
            f[16] = f[16] - 1.0
    if asyfp:
        if f[18] == 5.0:
            run('M37O45')
        run('M52O64')
    if not (symfp or asyfp) and (trim and htpl):
        run('M38O46')
    return ran


def main06(symfp: bool, asyfp: bool, control_type: float, mach: float,
           call: Callable[[str], object]) -> List[str]:
    """MAIN06: the transonic and supersonic flap and control overlays.
    ``control_type`` is ``F(18)``."""
    ran: List[str] = []

    def run(name):
        ran.append(name)
        call(name)

    if symfp:
        run('M36O44')
    if asyfp:
        run('M52O64')
        if control_type == 5.0 and mach < 1.0:
            run('M37O45')
    run('M40O50')
    return ran


def main07(symfp: bool, asyfp: bool,
           call: Callable[[str], object]) -> List[str]:
    """MAIN07: the supersonic hinge-moment and control-derivative
    overlays."""
    ran: List[str] = []
    if symfp:
        ran.append('M41O51')
        call('M41O51')
    if asyfp:
        ran.append('M53O65')
        call('M53O65')
    return ran
