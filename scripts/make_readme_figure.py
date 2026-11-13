"""Overview figure of README.md: four checks from the book, recomputed from scratch.

(a) Optimal truncation of the divergent Euler series (Chapter 1).
(b) Laplace's method for n!: the error falls like 1/n, and like 1/n^2 with the first
    correction (Chapter 6).
(c) The Airy function and its saddle-point forms on both sides (Chapter 8).
(d) Kramers' escape time against the exact double integral and a simulation (Chapter 10).
Output: docs/images/overview.png.   Run from the repository root.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import airy

try:
    import asymp as am
except ImportError:          # running from a clone without `pip install -e .`
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import asymp as am

am.use_style()
fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))

# (a) Euler series
ax = axes[0, 0]
for x, c in zip((0.05, 0.1, 0.2), ("C0", "C1", "C3")):
    exact = float(mp.exp(1 / x) * mp.e1(1 / x) / x)
    s = np.cumsum([(-1) ** n * math.factorial(n) * x**n for n in range(40)])
    err = np.abs(s - exact)
    ax.semilogy(np.arange(1, 41), err, "o-", ms=2.5, color=c, label=f"$x = {x:g}$")
    ax.axhline(err.min(), color=c, ls=":", lw=0.8)
ax.set_ylim(1e-12, 1e4)
ax.set_xlabel("number of terms $N$"); ax.set_ylabel("error of the partial sum of the Euler series")
ax.legend(loc="lower right", fontsize=8.5); am.panel_label(ax, "(a)", y=0.08)

# (b) Laplace's method for n!
ax = axes[0, 1]
ns = np.logspace(0.5, 2.5, 20)


def ratio(n, corr):
    I = am.laplace_method(n, f0=0.0, f2=-1.0, f3=2.0, f4=-6.0, g0=1.0, correction=corr)
    return math.exp(math.log(I) + (n + 1) * math.log(n) - n - math.lgamma(n + 1))


for corr, c, lab in ((False, "C0", "leading term"), (True, "C3", "with first correction")):
    e = np.array([abs(ratio(n, corr) - 1) for n in ns])
    s, _ = am.loglog_slope(ns, e, fit_range=(10, 300))
    ax.loglog(ns, e, "o-", ms=3, color=c, label=f"{lab} (slope {s:.2f})")
ax.set_xlabel("$n$"); ax.set_ylabel(r"relative error of $n!$ by Laplace's method")
ax.legend(loc="lower left", fontsize=8.5); am.panel_label(ax, "(b)", x=0.9)

# (c) Airy function
ax = axes[1, 0]
x = np.linspace(-12, 3, 1500)
ai = airy(x)[0]
xp, xm = x[x > 0.3], x[x < -0.3]
zp, zm = 2 / 3 * xp**1.5, 2 / 3 * (-xm) ** 1.5
ax.plot(x, ai, "k", lw=1.8, label=r"$\mathrm{Ai}(x)$")
ax.plot(xp, np.exp(-zp) / (2 * np.sqrt(np.pi) * xp**0.25), "C0--", lw=1.4, label="one saddle")
ax.plot(xm, np.sin(zm + np.pi / 4) / (np.sqrt(np.pi) * (-xm) ** 0.25), "C3--", lw=1.4, label="two saddles")
ax.set_ylim(-0.6, 0.75)
ax.set_xlabel("$x$"); ax.set_ylabel("value")
ax.legend(loc="upper left", fontsize=8.5); am.panel_label(ax, "(c)", x=0.9)

# (d) Kramers
ax = axes[1, 1]
u = lambda y: 3 * y**2 - 2 * y**3


def tau(lam):
    inner = lambda y: quad(lambda z: math.exp(-lam * u(z)), -1.5, y, epsrel=1e-12, limit=200)[0]
    return quad(lambda y: math.exp(lam * (u(y) - 1)) * inner(y), 0, 2, points=[1], epsrel=1e-11,
                limit=200)[0] * math.exp(lam)


lams = np.linspace(2, 25, 24)
t = np.array([tau(l) for l in lams])
k = math.pi * np.exp(lams) / (3 * lams)
ax.semilogy(lams, t, "k", lw=1.8, label="exact double integral")
ax.semilogy(lams, k, "C0--", lw=1.4, label=r"Kramers: $\pi e^{\lambda}/3\lambda$")
rng = np.random.default_rng(3)
sim_l, sim_t, sim_e = [], [], []
for l in (2, 3, 4, 5):
    xw, alive, tex, tt, dt = np.zeros(1500), np.arange(1500), np.zeros(1500), 0.0, 2e-3
    while alive.size:
        xa = xw[alive]
        xa += -l * 6 * xa * (1 - xa) * dt + math.sqrt(2 * dt) * rng.standard_normal(xa.size)
        xw[alive] = xa
        tt += dt
        out = xa >= 2
        tex[alive[out]] = tt
        alive = alive[~out]
    sim_l.append(l); sim_t.append(tex.mean()); sim_e.append(tex.std() / math.sqrt(tex.size))
ax.errorbar(sim_l, sim_t, yerr=sim_e, fmt="o", color="C3", ms=4, label="Brownian simulation")
ax.set_xlabel(r"barrier height $\lambda = \Delta U/k_BT$"); ax.set_ylabel(r"mean escape time $\tau D/L^2$")
ax.legend(loc="lower right", fontsize=8.5); am.panel_label(ax, "(d)")

fig.savefig("docs/images/overview.png", dpi=110)
print("wrote docs/images/overview.png")
