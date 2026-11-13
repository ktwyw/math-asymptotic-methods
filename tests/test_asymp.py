"""Tests of asymp against closed-form results."""

import math

import mpmath as mp
import numpy as np
import pytest

import asymp as am


# ---------------------------------------------------------------- order
def test_loglog_slope_recovers_power_law():
    x = np.logspace(-4, -1, 20)
    slope, pref = am.loglog_slope(x, 3.0 * x**2.5)
    assert slope == pytest.approx(2.5, abs=1e-12)
    assert pref == pytest.approx(3.0, rel=1e-10)


def test_loglog_slope_fit_range_excludes_crossover():
    x = np.logspace(-6, 0, 61)
    err = x**2 + 1e-10                     # power law, then a floor (rounding)
    slope, _ = am.loglog_slope(x, err, fit_range=(1e-3, 1e-1))
    assert slope == pytest.approx(2.0, abs=0.01)


def test_loglog_slope_needs_two_points():
    with pytest.raises(ValueError):
        am.loglog_slope([1.0, 2.0], [0.0, 0.0])


def test_local_slopes_constant_for_power_law():
    x = np.logspace(-3, 0, 10)
    _, s = am.local_slopes(x, x**3)
    assert np.allclose(s, 3.0)


def test_relative_error():
    assert am.relative_error(1.01, 1.0) == pytest.approx(0.01)


# ---------------------------------------------------------------- series
def test_partial_sums_geometric():
    s = am.partial_sums(0.5 ** np.arange(5))
    assert s[-1] == pytest.approx(1.9375)


def test_optimal_truncation_euler_series():
    # Euler's series sum (-1)^n n! x^n for x = 0.1: smallest term at n = 1/x = 10
    x = 0.1
    terms = [(-1) ** n * math.factorial(n) * x**n for n in range(30)]
    tr = am.optimal_truncation(terms)
    assert tr.n_terms in (9, 10)
    # exact value: int_0^inf e^-t / (1 + x t) dt
    exact = float(mp.quad(lambda t: mp.exp(-t) / (1 + x * t), [0, mp.inf]))
    assert abs(tr.value - exact) < 2 * tr.smallest


def test_optimal_truncation_requires_minimum():
    with pytest.raises(ValueError):
        am.optimal_truncation([3.0, 2.0, 1.0])


def test_shanks_exact_for_geometric_sequence():
    s = am.partial_sums(0.9 ** np.arange(6))
    assert am.shanks(s)[0] == pytest.approx(10.0, rel=1e-12)


def test_shanks_accelerates_log2():
    s = am.partial_sums([(-1) ** n / (n + 1) for n in range(12)])
    once = am.shanks(s)
    assert abs(once[-1] - math.log(2)) < 1e-2 * abs(s[-1] - math.log(2))   # about 500 times better


def test_richardson_removes_leading_error():
    f = lambda h: 1.0 + 0.5 * h**2                          # error exactly C h^2
    val, err = am.richardson(f(0.1), f(0.2), p=2)
    assert val == pytest.approx(1.0, abs=1e-14)
    assert err == pytest.approx(0.005)


def test_pade_reproduces_exp():
    coeffs = [1 / math.factorial(k) for k in range(9)]
    r = am.pade_approximant(coeffs, 4, 4)
    assert r(1.0) == pytest.approx(math.e, rel=1e-7)          # [4/4] error is 4e-8
    assert r(np.array([0.0, 0.5]))[0] == pytest.approx(1.0)


def test_pade_needs_enough_coefficients():
    with pytest.raises(ValueError):
        am.pade_approximant([1, 1, 0.5], 2, 2)


# ---------------------------------------------------------------- reference
def test_mp_quad_gaussian():
    val = am.mp_quad(lambda t: mp.exp(-t**2), -mp.inf, mp.inf, dps=40)
    assert val == pytest.approx(math.sqrt(math.pi), rel=1e-15)


def test_mp_map_and_sum():
    assert am.mp_map(mp.erfc, [0.0])[0] == pytest.approx(1.0)
    assert am.mp_sum(lambda n: 1 / n**2, 1) == pytest.approx(math.pi**2 / 6, rel=1e-14)


# ---------------------------------------------------------------- integrals
def test_laplace_method_stirling():
    # Gamma(lam+1) = lam^(lam+1) int exp(lam (ln t - t)) dt; f = ln t - t at t0 = 1
    for lam in (10.0, 50.0):
        I = am.laplace_method(lam, f0=-1.0, f2=-1.0, g0=1.0, f3=2.0, f4=-6.0)
        approx = lam ** (lam + 1) * I
        exact = math.gamma(lam + 1)
        assert approx / exact - 1 == pytest.approx(-1 / (288 * lam**2), rel=0.1)  # minus the next Stirling term


def test_laplace_method_bessel_i0():
    lam = 40.0
    I = am.laplace_method(lam, f0=1.0, f2=-1.0, g0=1 / (2 * math.pi), f4=1.0)
    exact = float(mp.besseli(0, lam))
    assert I / exact - 1 == pytest.approx(-9 / (128 * lam**2), rel=0.05)      # minus the next term of I0


def test_laplace_method_leading_only_and_bad_curvature():
    lam = 100.0
    lead = am.laplace_method(lam, 0.0, -2.0, 1.0, correction=False)
    assert lead == pytest.approx(math.sqrt(math.pi / lam))
    with pytest.raises(ValueError):
        am.laplace_method(lam, 0.0, 1.0, 1.0)
