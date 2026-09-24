"""
Hypersonic transverse jet sizing: TRANJT.

TRANJT sizes a transverse reaction-control jet on a hypersonic body
(Section 6.3.2): at each point of a force-time schedule it reads the
local flow ahead of the jet (Figures 6.3.2-30 to -33), the interaction
amplification K0 for a laminar (Figure 6.3.2-35A-E, through INTER3) or
turbulent (Figure 6.3.2-40) layer, the jet thrust and total pressure
(Figure 6.3.2-43), then the throat size, the propellant rate and weight,
and the centre of pressure of the interaction.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/tranjt.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.utils.constants import RAD
from pydatcom.utils.legacy_numeric import inter3, simul2, trapz
from pydatcom.utils.legacy_tables import tlinex

# Figures 6.3.2-30 to -33: the local flow ahead of the jet (pressure,
# dynamic pressure, Mach number, Reynolds number) against Mach number and
# angle of attack.  Figure 30's angle grid repeats 4.04.
_X13230 = np.array([
    30.0, 20.0, 15.0, 10.0, 7.0, 5.0,
])
_X23230 = np.array([
    0.0907, 2.06, 4.04, 4.04, 5.98, 8.0, 10.0, 12.1, 14.0, 16.0, 18.0,
])
_Y63230 = np.array([
    1.41, 4.45, 9.83, 9.83, 18.7, 32.1, 47.7, 67.1, 89.4, 116.0, 148.0,
    1.39, 2.72, 5.24, 5.24, 9.54, 15.4, 22.6, 31.9, 42.4, 54.2, 66.7,
    1.37, 1.85, 4.07, 4.07, 6.08, 9.18, 13.6, 18.6, 24.6, 31.4, 38.4,
    1.38, 1.63, 2.58, 2.58, 4.06, 5.26, 7.28, 9.89, 12.1, 15.0, 18.4,
    1.38, 1.63, 1.91, 1.91, 2.98, 3.95, 4.65, 5.98, 7.32, 8.65, 10.5,
    1.37, 1.41, 1.46, 1.46, 1.95, 3.3, 3.54, 4.02, 4.71, 5.4, 6.33,
])
_X13231 = np.array([
    20.0, 15.0, 10.0, 7.0, 5.0,
])
_X23231 = np.array([
    0.0184, 1.99, 4.02, 5.99, 7.94, 9.93, 11.9, 13.9, 16.0, 17.9, 19.9,
])
_Y63231 = np.array([
    1.01, 1.87, 2.85, 3.69, 4.2, 4.55, 4.75, 4.82, 4.81, 4.77, 4.63,
    0.996, 1.66, 2.36, 3.05, 3.59, 3.98, 4.27, 4.43, 4.48, 4.47, 4.43,
    1.01, 1.37, 1.85, 2.31, 2.71, 3.09, 3.42, 3.66, 3.8, 3.88, 3.89,
    1.0, 1.24, 1.53, 1.84, 2.15, 2.43, 2.7, 2.9, 3.07, 3.18, 3.23,
    1.01, 1.15, 1.35, 1.53, 1.75, 1.94, 2.12, 2.29, 2.4, 2.49, 2.56,
])
_X13232 = np.array([
    30.0, 20.0, 15.0, 10.0, 7.0, 5.0,
])
_X23232 = np.array([
    0.00264, 2.04, 4.0, 6.02, 8.03, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
])
_Y63232 = np.array([
    30.0, 24.1, 19.1, 14.9, 11.9, 9.76, 8.28, 7.15, 6.15, 5.5, 4.84,
    20.0, 17.1, 14.5, 12.4, 10.5, 8.86, 7.71, 6.62, 5.93, 5.28, 4.8,
    15.1, 13.5, 12.1, 10.8, 9.45, 8.37, 7.27, 6.5, 5.75, 5.06, 4.54,
    10.1, 9.34, 8.69, 8.0, 7.31, 6.67, 6.14, 5.67, 5.1, 4.66, 4.32,
    6.95, 6.69, 6.25, 5.96, 5.7, 5.22, 4.93, 4.48, 4.22, 3.75, 3.45,
    5.03, 4.81, 4.64, 4.47, 4.21, 3.91, 3.74, 3.61, 3.35, 3.09, 2.92,
])
_X13233 = np.array([
    30.0, 20.0, 15.0, 10.0, 7.0, 5.0,
])
_X23233 = np.array([
    -0.0475, 1.96, 2.91, 3.45, 4.02, 5.35, 6.06, 6.28, 8.06, 10.1, 12.0, 14.0, 16.0, 18.0, 20.0,
])
_Y63233 = np.array([
    1.01, 1.8, 1.96, 1.98, 1.97, 1.87, 1.77, 1.72, 1.47, 1.2, 1.01, 0.84, 0.71, 0.61, 0.53,
    1.01, 1.57, 1.77, 1.85, 1.91, 1.99, 1.99, 1.98, 1.84, 1.64, 1.45, 1.26, 1.09, 0.949, 0.842,
    1.01, 1.43, 1.62, 1.72, 1.78, 1.91, 1.96, 1.97, 1.96, 1.85, 1.72, 1.56, 1.41, 1.27, 1.13,
    1.01, 1.29, 1.41, 1.48, 1.53, 1.69, 1.76, 1.77, 1.89, 1.93, 1.91, 1.84, 1.76, 1.66, 1.54,
    1.01, 1.22, 1.3, 1.35, 1.4, 1.5, 1.56, 1.57, 1.69, 1.8, 1.85, 1.88, 1.85, 1.82, 1.75,
    1.01, 1.15, 1.22, 1.26, 1.29, 1.36, 1.41, 1.42, 1.52, 1.62, 1.68, 1.74, 1.76, 1.76, 1.75,
])
# Figure 6.3.2-35A-E: K0 for a laminar layer, in five Reynolds-number
# decades (INTER3); part B reads part A's second grid.  Part E's second
# row starts 1.60 where the row runs near 9.6, and part C's third row
# holds a 3.79 among values near 2.8; both kept.
_X1235A = np.array([
    20.0, 10.0, 4.0,
])
_X2235A = np.array([
    0.004, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 6.0,
])
_Y3235A = np.array([
    3.4, 3.3, 3.19, 2.94, 2.67, 2.45, 2.3, 1.98, 1.74, 1.56, 1.5, 1.36, 1.21,
    3.17, 3.1, 3.0, 2.84, 2.68, 2.5, 2.37, 2.02, 1.72, 1.56, 1.51, 1.36, 1.21,
    2.95, 2.94, 2.9, 2.8, 2.69, 2.6, 2.52, 2.19, 1.85, 1.58, 1.5, 1.36, 1.21,
])
_X1235B = np.array([
    20.0, 10.0, 6.0, 4.0, 2.0,
])
_Y3235B = np.array([
    4.1, 4.01, 3.68, 3.05, 2.69, 2.47, 2.32, 2.0, 1.75, 1.56, 1.5, 1.32, 1.2,
    4.0, 3.88, 3.63, 3.19, 2.81, 2.58, 2.4, 2.04, 1.78, 1.57, 1.5, 1.32, 1.2,
    3.84, 3.8, 3.61, 3.22, 2.92, 2.69, 2.5, 2.11, 1.8, 1.6, 1.52, 1.32, 1.2,
    3.67, 3.61, 3.5, 3.3, 3.09, 2.89, 2.7, 2.2, 1.85, 1.6, 1.52, 1.32, 1.2,
    3.23, 3.2, 3.14, 3.04, 2.96, 2.88, 2.8, 2.41, 2.03, 1.73, 1.68, 1.41, 1.2,
])
_X1235C = np.array([
    20.0, 15.0, 10.0, 6.0, 4.0, 2.0,
])
_X2235C = np.array([
    0.0004, 0.001, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.8,
])
_Y3235C = np.array([
    5.68, 5.61, 5.32, 5.12, 4.93, 4.72, 3.8, 3.0, 2.67, 2.32, 2.0, 1.77, 1.59,
    5.5, 5.49, 5.2, 5.04, 4.88, 4.72, 3.95, 3.1, 2.67, 2.37, 2.0, 1.77, 1.59,
    5.36, 5.31, 5.1, 4.96, 4.71, 4.68, 4.0, 3.19, 3.79, 2.4, 2.0, 1.77, 1.59,
    4.94, 4.93, 4.82, 4.78, 4.7, 4.6, 4.14, 3.41, 2.98, 2.55, 2.12, 1.82, 1.6,
    4.5, 4.46, 4.4, 4.38, 4.32, 4.28, 4.0, 3.53, 3.15, 2.65, 2.18, 1.82, 1.6,
    3.6, 3.6, 3.59, 3.58, 3.55, 3.53, 3.48, 3.34, 3.2, 2.88, 2.36, 2.0, 1.72,
])
_X1235D = np.array([
    20.0, 15.0, 10.0, 6.0, 4.0, 2.0,
])
_X2235D = np.array([
    0.0004, 0.001, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.1, 0.4, 0.8,
])
_Y3235D = np.array([
    7.55, 7.4, 6.6, 6.1, 5.6, 5.2, 3.89, 3.12, 2.75, 2.4, 1.77, 1.57, 7.33,
    7.2, 6.49, 6.05, 5.65, 5.3, 4.0, 3.15, 2.75, 2.4, 1.77, 1.57, 6.89, 6.84,
    6.38, 6.01, 5.7, 5.4, 4.18, 3.26, 2.83, 2.48, 1.77, 1.57, 6.02, 6.02, 5.8,
    5.6, 5.43, 5.29, 4.39, 3.4, 2.99, 2.6, 1.9, 1.7, 5.2, 5.19, 5.03, 4.99,
    4.9, 4.8, 4.38, 3.6, 3.09, 2.6, 1.9, 1.7, 3.81, 3.83, 3.8, 3.78, 3.73,
    3.7, 3.6, 3.48, 3.3, 2.9, 2.0, 1.7,
])
_X1235E = np.array([
    20.0, 15.0, 10.0, 6.0, 4.0,
])
_X2235E = np.array([
    0.0004, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.1, 0.2, 0.8,
])
_Y3235E = np.array([
    10.0, 9.52, 8.8, 7.45, 6.39, 5.6, 5.1, 3.88, 3.02, 2.37, 1.99, 1.5, 1.6,
    9.29, 8.7, 7.6, 6.6, 5.8, 5.2, 3.95, 3.04, 2.37, 1.99, 1.52, 8.73, 8.58,
    8.2, 7.5, 6.9, 6.26, 5.62, 4.13, 3.21, 2.41, 2.01, 1.54, 7.12, 7.0, 6.9,
    6.6, 6.32, 6.0, 5.7, 4.47, 3.4, 2.58, 2.1, 1.56, 5.72, 5.71, 5.67, 5.6,
    5.48, 5.34, 5.2, 4.52, 3.57, 2.62, 2.2, 1.59,
])
# Figure 6.3.2-43: the jet total-pressure ratio.
_X13243 = np.array([
    1.0, 1.5, 1.75, 2.0,
])
_X23243 = np.array([
    2.0, 4.0, 6.0, 8.0, 10.0, 12.0,
])
_Y63243 = np.array([
    7.5, 11.4, 15.3, 18.9, 22.6, 26.3,
    21.2, 31.2, 41.8, 52.0, 63.1, 74.0,
    50.0, 73.1, 97.0, 121.0, 148.0, 171.0,
    192.0, 282.0, 376.0, 477.0, 573.0, 650.0,
])
# Figure 6.3.2-40: K0 for a turbulent layer.
_X13240 = np.array([
    20.0, 15.0, 10.0, 6.0, 4.0, 2.0,
])
_X23240 = np.array([
    4e-05, 0.0001, 0.0002, 0.0004, 0.002, 0.006, 0.01, 0.02, 0.04, 0.1,
])
_Y63240 = np.array([
    3.08, 2.92, 2.78, 2.62, 2.21, 1.92, 1.8, 1.65, 1.49, 1.28,
    3.02, 2.96, 2.88, 2.76, 2.4, 2.11, 1.98, 1.8, 1.63, 1.4,
    2.95, 2.94, 2.9, 2.85, 2.63, 2.37, 2.23, 2.05, 1.87, 1.64,
    2.75, 2.75, 2.75, 2.75, 2.7, 2.58, 2.48, 2.33, 2.17, 1.8,
    2.5, 2.5, 2.52, 2.54, 2.55, 2.53, 2.49, 2.41, 2.29, 2.04,
    2.23, 2.23, 2.22, 2.21, 2.21, 2.21, 2.2, 2.2, 2.2, 2.2,
])


def _grid(flat, n1, n2):
    return np.asarray(flat, dtype=float).reshape(n1, n2).T


def _tl(x1, x2, flat, q1, q2, l1=0, l2=0, u1=0, u2=0):
    return float(tlinex(x1, x2, _grid(flat, len(x1), len(x2)), q1, q2,
                        l1, l2, u1, u2))


_K0_LAMINAR = [(_X1235A, _X2235A, _grid(_Y3235A, 3, 13)),
               (_X1235B, _X2235A, _grid(_Y3235B, 5, 13)),
               (_X1235C, _X2235C, _grid(_Y3235C, 6, 13)),
               (_X1235D, _X2235D, _grid(_Y3235D, 6, 12)),
               (_X1235E, _X2235E, _grid(_Y3235E, 5, 12))]


def calculate_tranjt(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate TRANJT: hypersonic transverse jet sizing.

    Args:
        data: By name: ``mach``, ``rl`` (``FLC(IM+42)``, the unit
            Reynolds number), ``pinf`` (``FLC(IM+73)``), and the
            ``/FLAPIN/`` jet words: ``time`` (``F(1..)``), ``nt``
            (``F(11)``), ``fc`` (``F(12..)``, the control forces),
            ``alpha`` (``F(22..)``), ``me`` (exit Mach, ``F(32)``),
            ``isp`` (``F(33)``), ``span`` (``F(34)``), ``phe`` (``F(35)``),
            ``gp`` (the jet's gamma, ``F(36)``), ``cc`` (``F(37)``), ``l``
            (``F(38)``), ``laminar`` (``F(39..)``).

    Returns:
        The ``JET`` curves ``m1``, ``rl``, ``p1``, ``q1``, ``cfc``,
        ``cfcr``, ``k0``, ``k``, ``fj0``, ``pj0p1m``, ``p0j``, ``p0jt``,
        ``rate``, ``weight`` and ``xcp`` (1-based, as ``JET(1..147)``
        after the overlapping writes), and the ``JETA`` words ``qinf``,
        ``cf0``, ``veoa``, ``fjmax``, ``pjmax``, ``dt`` and ``k``.

    Notes:
        Kept as executed: ``XCP`` starts at ``JET(138)``, inside
        ``WEIGHT`` (``JET(131..140)``), so with more than seven schedule
        points the centres of pressure overwrite the last weights; and the
        suspect entries of Figure 6.3.2-35C and E are read as they stand.
    """
    mach, rln, pinf = (float(data['mach']), float(data['rl']),
                       float(data['pinf']))
    nt = int(float(data['nt']) + .5)
    fc = [float(v) for v in data['fc']]
    alpha = [float(v) for v in data['alpha']]
    lam = [bool(v) for v in data['laminar']]
    me, isp, span, phe = (float(data['me']), float(data['isp']),
                          float(data['span']), float(data['phe']))
    gp, cc, ell = float(data['gp']), float(data['cc']), float(data['l'])
    time = [float(v) for v in data['time']]
    qinf = 0.7 * pinf * mach**2
    arg = gp + 1.
    cf0 = (2. / arg)**(gp / (gp - 1.)) * arg
    names = ('m1', 'rl', 'p1', 'q1', 'cfc', 'cfcr', 'k0', 'k', 'fj0',
             'pj0p1m', 'p0j', 'p0jt', 'rate')
    c = {k: [0.0] * nt for k in names}
    veoa = 0.0
    for time_index in range(nt):
        p1pi = _tl(_X13230, _X23230, _Y63230, mach, alpha[time_index], 0, 0, 0, 1)
        q1qi = _tl(_X13231, _X23231, _Y63231, mach, alpha[time_index])
        c['m1'][time_index] = _tl(_X13232, _X23232, _Y63232, mach, alpha[time_index])
        r1ri = _tl(_X13233, _X23233, _Y63233, mach, alpha[time_index])
        c['p1'][time_index] = p1pi * pinf
        c['q1'][time_index] = q1qi * qinf
        c['rl'][time_index] = r1ri * rln * ell
        c['cfc'][time_index] = fc[time_index] / (c['q1'][time_index] * span * ell * 144.)
        c['cfcr'][time_index] = 1.268 * c['cfc'][time_index] / cf0
        if lam[time_index]:
            c['k0'][time_index] = inter3(c['m1'][time_index], c['cfcr'][time_index], c['rl'][time_index],
                                _K0_LAMINAR)
        else:
            c['k0'][time_index] = _tl(_X13240, _X23240, _Y63240, c['m1'][time_index],
                             c['cfcr'][time_index])
        veoa = math.sqrt((arg * me**2) / (2. + (gp - 1.) * me**2))
        ph = phe / RAD
        arg1 = 1. + gp * veoa * math.sin(ph) / arg
        arg2 = (veoa + 1. / veoa) * math.cos(ph) / 2.
        c['k'][time_index] = (c['k0'][time_index] - 1.) * arg1 + arg2
        c['fj0'][time_index] = fc[time_index] / c['k'][time_index]
        c['pj0p1m'][time_index] = _tl(_X13243, _X23243, _Y63243, veoa, c['m1'][time_index])
        c['p0j'][time_index] = c['pj0p1m'][time_index] * c['p1'][time_index]
    fjmax, pjmax = c['fj0'][0], c['p0j'][0]
    for time_index in range(nt):
        if c['fj0'][time_index] > fjmax:
            fjmax = c['fj0'][time_index]
        if c['p0j'][time_index] > pjmax:
            pjmax = c['p0j'][time_index]
    arg = cc * cf0 * span * 12.0
    dt = fjmax / (arg * pjmax)
    for time_index in range(nt):
        c['p0jt'][time_index] = c['fj0'][time_index] / (arg * dt)
        c['rate'][time_index] = c['fj0'][time_index] / isp
    weight = list(trapz(c['rate'], time[:nt], 0))
    jet = [0.0] * 148
    for name, start in (('m1', 1), ('rl', 11), ('p1', 21), ('q1', 31),
                        ('cfc', 41), ('cfcr', 51), ('k0', 61), ('k', 71),
                        ('fj0', 81), ('pj0p1m', 91), ('p0j', 101),
                        ('p0jt', 111), ('rate', 121)):
        for slot, v in enumerate(c[name]):
            jet[start + slot] = v
    for slot, v in enumerate(weight):
        jet[131 + slot] = v
    xcp = []
    for time_index in range(nt):
        r2 = c['m1'][time_index]**2
        g = 1. - 1. / c['k'][time_index]
        if lam[time_index]:
            rs = [.2 * c['rl'][time_index]]
            for _ in range(9):
                rs.append(.75 * rs[-1])
            rs.reverse()
            a3 = (2. * cf0 * dt / (ell * 12.)) * c['p0jt'][time_index] / c['p1'][time_index]
            hl, hln = [], []
            for reynolds_step in range(10):
                a = rs[reynolds_step] * (r2 - 1.)
                cp2 = 1.60 / a**0.25
                cx = 4.75 * cp2
                hl.append(a3 / (gp * cx * r2 + 2.))
                zi = (cp2 * c['q1'][time_index] + c['p1'][time_index]) / c['p1'][time_index]
                a1 = 5. * (zi - 1.) / (7. * r2 - 5. * (zi - 1.))
                a = math.sqrt((7. * r2 - (6. * zi + 1.)) / (6. * zi + 1.))
                hln.append((1. - rs[reynolds_step] / c['rl'][time_index]) *
                           a1 * a)
            rls, _ = simul2(rs, hl, hln)
            if rls == -1000.:
                rls = rs[9]
            cp2 = 1.60 / (rls * (r2 - 1.))**0.25
        else:
            if c['m1'][time_index] > 5.0:
                cp2 = 0.2257 - 0.0232 * c['m1'][time_index] + 0.0014 * r2 - \
                    0.00003 * r2 * c['m1'][time_index]
            else:
                cp2 = 0.41 + 0.481 * c['m1'][time_index] - 0.0509 * r2 + \
                    0.0061 * r2 * c['m1'][time_index]
        a = 1. - g / 2. * c['cfc'][time_index] / cp2
        xcp.append((1. - g) + g * a)
    for slot, v in enumerate(xcp):
        jet[138 + slot] = v
    out = dict(c)
    out.update({'weight': weight, 'xcp': xcp, 'jet': jet[1:],
                'qinf': qinf, 'cf0': cf0, 'veoa': veoa, 'fjmax': fjmax,
                'pjmax': pjmax, 'dt': dt, 'method': 'legacy_tranjt'})
    return out
