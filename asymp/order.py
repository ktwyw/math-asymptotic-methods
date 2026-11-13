"""Measuring how fast an error goes to zero.

An asymptotic result predicts not only a value but the size of its error: if the
first neglected term is of order eps**p, a log-log plot of the error against eps
has slope p.  These functions measure that slope, so that every expansion in the
book can be checked against a numerical reference.
"""

from __future__ import annotations

import numpy as np

__all__ = ["loglog_slope", "local_slopes", "relative_error"]


def relative_error(approx, exact):
    """|approx - exact| / |exact|, element by element."""
    approx = np.asarray(approx, dtype=float)
    exact = np.asarray(exact, dtype=float)
    return np.abs(approx - exact) / np.abs(exact)


def loglog_slope(x, err, fit_range=None):
    """Least-squares slope of log(err) against log(x).

    Parameters
    ----------
    x, err : array_like
        Small parameter and the (positive) error at each value.
    fit_range : (xmin, xmax), optional
        Use only the points with xmin <= x <= xmax, typically the smallest x
        where the asymptotic regime holds but rounding error does not yet matter.

    Returns
    -------
    slope, prefactor
        err is approximately prefactor * x**slope over the fitted range.
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    keep = (err > 0) & np.isfinite(err)
    if fit_range is not None:
        keep &= (x >= fit_range[0]) & (x <= fit_range[1])
    if keep.sum() < 2:
        raise ValueError("need at least two positive errors in the fitted range")
    slope, intercept = np.polyfit(np.log(x[keep]), np.log(err[keep]), 1)
    return slope, float(np.exp(intercept))


def local_slopes(x, err):
    """Slopes d log(err) / d log(x) between neighbouring points.

    A constant local slope shows that the error follows a single power law; a
    drift shows a crossover (for example to the next term, or to rounding error).
    Returns the geometric midpoints of x and the slopes.
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    lx, le = np.log(x), np.log(err)
    return np.exp(0.5 * (lx[1:] + lx[:-1])), np.diff(le) / np.diff(lx)
