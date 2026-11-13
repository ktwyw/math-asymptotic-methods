"""Asymptotic evaluation of integrals with a large parameter.

Laplace's method for I(lam) = integral of g(t) exp(lam f(t)) dt, dominated by an
interior maximum t0 of f, with its first correction.
"""

from __future__ import annotations

import math

__all__ = ["laplace_method"]


def laplace_method(lam, f0, f2, g0, f3=0.0, f4=0.0, g1=0.0, g2=0.0, correction=True):
    """Laplace's method for an interior maximum.

    The integral is I(lam) = int g(t) exp(lam f(t)) dt, with f'(t0) = 0 and
    f''(t0) = f2 < 0.  Arguments are the values at t0 of f, f'', f''', f'''' and
    of g, g', g''.  Returns

        exp(lam f0) sqrt(2 pi / (lam a)) [g0 + c1 / lam]          (a = -f2)

    with  c1 = g2/(2a) + g1 f3/(2a^2) + g0 f4/(8a^2) + 5 g0 f3^2/(24 a^3),
    or the leading term alone if correction is False.  The relative error is
    O(1/lam^2) with the correction and O(1/lam) without it.
    """
    if f2 >= 0:
        raise ValueError("f''(t0) must be negative at a maximum")
    a = -f2
    lead = math.exp(lam * f0) * math.sqrt(2 * math.pi / (lam * a))
    if not correction:
        return lead * g0
    c1 = g2 / (2 * a) + g1 * f3 / (2 * a**2) + g0 * f4 / (8 * a**2) + 5 * g0 * f3**2 / (24 * a**3)
    return lead * (g0 + c1 / lam)
