"""Working with series, convergent or divergent.

Partial sums, optimal truncation of a divergent asymptotic series at its smallest
term, Shanks and Richardson acceleration, and Padé approximants.
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp
import numpy as np

__all__ = [
    "partial_sums",
    "Truncation",
    "optimal_truncation",
    "shanks",
    "richardson",
    "pade_approximant",
]


def partial_sums(terms):
    """Cumulative sums S_0, S_1, ... of a sequence of terms."""
    return np.cumsum(np.asarray(terms, dtype=float))


@dataclass
class Truncation:
    """Result of truncating an asymptotic series at its smallest term."""

    n_terms: int        # number of terms kept (indices 0 .. n_terms - 1)
    value: float        # the truncated sum
    smallest: float     # magnitude of the smallest term, the first one omitted
    terms: np.ndarray   # the terms examined


def optimal_truncation(terms):
    """Truncate a (divergent) asymptotic series just before its smallest term.

    The rule of thumb of optimal (least-term) truncation: keep all terms before
    the smallest one; the error is then of the order of that smallest term.
    The terms must be given far enough for the minimum to be passed.
    """
    terms = np.asarray(terms, dtype=float)
    mags = np.abs(terms)
    k = int(np.argmin(mags))
    if k == len(terms) - 1:
        raise ValueError("the smallest term is the last one given; supply more terms")
    return Truncation(n_terms=k, value=float(terms[:k].sum()), smallest=float(mags[k]),
                      terms=terms)


def shanks(seq):
    """One Shanks transformation of a sequence of partial sums.

    S'_n = (S_{n+1} S_{n-1} - S_n**2) / (S_{n+1} + S_{n-1} - 2 S_n),
    exact for a sequence of the form A + B q**n.  Returns a sequence two shorter.
    """
    s = np.asarray(seq, dtype=float)
    num = s[2:] * s[:-2] - s[1:-1] ** 2
    den = s[2:] + s[:-2] - 2 * s[1:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(den != 0, num / den, s[2:])
    return out


def richardson(f_h, f_2h, p, ratio=2.0):
    """Richardson extrapolation for a quantity with error C h**p.

    Returns the extrapolated value and the error estimate of f_h.
    """
    delta = (f_h - f_2h) / (ratio**p - 1.0)
    return f_h + delta, abs(delta)


def pade_approximant(coeffs, m, n):
    """[m/n] Padé approximant from the Taylor coefficients c_0, c_1, ...

    Returns a callable x -> P(x) / Q(x), with numerator degree m and denominator
    degree n; needs at least m + n + 1 coefficients.
    """
    coeffs = np.asarray(coeffs, dtype=float)
    if len(coeffs) < m + n + 1:
        raise ValueError("need at least m + n + 1 coefficients")
    p, q = mp.pade([mp.mpf(float(c)) for c in coeffs[: m + n + 1]], m, n)
    p = np.array([float(c) for c in p])[::-1]          # highest power first, for np.polyval
    q = np.array([float(c) for c in q])[::-1]
    return lambda x: np.polyval(p, x) / np.polyval(q, x)
