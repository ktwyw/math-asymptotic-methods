"""High-precision reference values.

Asymptotic expansions are checked against "exact" numbers.  Double precision is not
always enough: the error of a good expansion can be 1e-20 or smaller, and an
integral of a sharply peaked or oscillating function may lose many digits in
floating point.  These helpers evaluate references with mpmath at a chosen number
of decimal digits and return them as floats (or as mpmath numbers if requested).
"""

from __future__ import annotations

import mpmath as mp
import numpy as np

__all__ = ["mp_quad", "mp_map", "mp_sum"]


def mp_quad(f, a, b, dps=30, as_float=True, **kwargs):
    """Integrate f from a to b with mpmath at `dps` decimal digits.

    f receives and returns mpmath numbers.  Infinite limits use mp.inf.  Extra
    keyword arguments go to mpmath.quad (for example maxdegree).  Splitting points
    can be passed by giving a list instead of (a, b) as `a` and b=None.
    """
    with mp.workdps(dps):
        points = a if b is None else [a, b]
        val = mp.quad(f, points, **kwargs)
        return float(val) if as_float else +val


def mp_map(f, xs, dps=30):
    """Evaluate f(x) for each x at `dps` digits, returning a float array."""
    with mp.workdps(dps):
        return np.array([float(f(mp.mpf(x))) for x in np.atleast_1d(xs)])


def mp_sum(f, n0, dps=30):
    """Sum f(n) for n = n0, n0 + 1, ... to infinity with mpmath (nsum)."""
    with mp.workdps(dps):
        return float(mp.nsum(f, [n0, mp.inf]))
