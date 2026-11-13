"""asymp: small, readable tools for checking asymptotic approximations.

The module accompanies the notebooks of math-asymptotic-methods.  It favours
clarity over generality: every function is short enough to be read alongside the
text.  It provides four things:

* order:      measure the order of an error from a log-log plot;
* series:     partial sums, optimal truncation, Shanks and Richardson
              acceleration, Padé approximants;
* reference:  high-precision reference values with mpmath;
* integrals:  Laplace's method with its first correction;

plus the plotting style used by every figure.
"""

from .integrals import laplace_method
from .order import local_slopes, loglog_slope, relative_error
from .plotting import panel_label, use_style
from .reference import mp_map, mp_quad, mp_sum
from .series import (
    Truncation,
    optimal_truncation,
    pade_approximant,
    partial_sums,
    richardson,
    shanks,
)

__version__ = "0.1.0"

__all__ = [
    "loglog_slope", "local_slopes", "relative_error",
    "partial_sums", "Truncation", "optimal_truncation", "shanks", "richardson", "pade_approximant",
    "mp_quad", "mp_map", "mp_sum",
    "laplace_method",
    "use_style", "panel_label",
]
