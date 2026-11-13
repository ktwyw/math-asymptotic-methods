"""Animation at the top of README.md: adding terms to a divergent series.

The Euler series sum (-1)^n n! x^n is the asymptotic expansion of
E(x) = int_0^inf e^{-t}/(1 + x t) dt as x -> 0. It diverges for every x > 0, yet a few
terms give E(x) accurately when x is small. Left: partial sums S_N(x) against E(x).
Right: the error at x = 0.1 as terms are added; it is smallest near N = 1/x.
Output: docs/images/truncation.gif.

Run from the repository root:  python scripts/make_readme_animation.py   (needs Pillow)
"""
import io
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from PIL import Image

try:
    import asymp as am
except ImportError:          # running from a clone without `pip install -e .`
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import asymp as am

am.use_style()
xs = np.linspace(1e-4, 0.5, 300)
E = np.array([float(mp.exp(1 / x) * mp.e1(1 / x) / x) for x in xs])
x0 = 0.1
E0 = float(mp.exp(1 / x0) * mp.e1(1 / x0) / x0)
coef = np.array([(-1) ** n * math.factorial(n) for n in range(26)], dtype=float)
err0 = np.abs(np.cumsum(coef * x0 ** np.arange(26)) - E0)
best = int(np.argmin(err0))

frames = []
for N in list(range(1, 26)) + [25] * 6:
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.9))
    S = np.polyval(coef[:N][::-1], xs)
    a1.plot(xs, E, "k", lw=2, label=r"$E(x) = \int_0^\infty \frac{e^{-t}}{1 + xt}\,dt$")
    a1.plot(xs, S, "C3", lw=1.6, label=f"partial sum, $N = {N}$")
    a1.axvline(x0, color="0.6", lw=0.8, ls=":")
    a1.set_xlim(0, 0.5); a1.set_ylim(0.4, 1.3)
    a1.set_xlabel("$x$"); a1.set_ylabel("value")
    a1.legend(loc="lower left", fontsize=9)
    a1.set_title(r"$E(x) \sim 1 - x + 2!\,x^2 - 3!\,x^3 + \cdots$", fontsize=10)

    a2.semilogy(np.arange(1, N + 1), err0[:N], "o-", color="C0", ms=4)
    a2.axhline(err0[best], color="0.6", lw=0.8, ls=":")
    if N > best:
        a2.semilogy([best + 1], [err0[best]], "o", color="C3", ms=8, mfc="none", mew=1.6)
        a2.text(best + 2.2, err0[best] * 0.08, f"smallest error: {best + 1} terms", color="C3", fontsize=9)
    a2.set_xlim(0, 26); a2.set_ylim(1e-5, 10)
    a2.set_xlabel("number of terms $N$"); a2.set_ylabel(r"error at $x = 0.1$")
    a2.set_title("the error falls, then grows", fontsize=10)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90)
    plt.close(fig)
    frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=64))

frames[0].save("docs/images/truncation.gif", save_all=True, append_images=frames[1:],
               duration=[350] * (len(frames) - 6) + [700] * 6, loop=0, optimize=True)
print("wrote docs/images/truncation.gif,", len(frames), "frames")
