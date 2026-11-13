# Asymptotic methods

*Asymptotic expansions of integrals and sums, with every result checked against
high-precision numerics.*

Most problems in physics and engineering have no closed-form solution, but many contain
a parameter that is large or small: a Reynolds number, a barrier height in units of
$k_BT$, a Fourier number, the argument of a special function. Asymptotic methods turn
that parameter into approximations that are often more accurate, and always more
revealing, than a numerical solution alone. They are also easy to get wrong: a series
can diverge, an expansion can fail in part of the domain, and a correction can be of
the wrong order.

Each chapter of this book is an executed Jupyter notebook. Every expansion is
derived and then compared with the exact answer computed to high precision, and the
comparison is made quantitative: the error should fall like the first omitted term,
with the predicted slope on a log-log plot. The notebooks share a small Python
module, `asymp`, written to be read as much as to be used.

* **Core** chapters assume calculus, ordinary differential equations and a first
  course in complex variables, at the level of a senior undergraduate.
* **Advanced** chapters are for graduate students and researchers.

Books and articles for further study are listed in [Suggested references](references.md).
Differential equations, with boundary layers, matching, WKB and multiple scales,
are the subject of a companion book,
[*Perturbation methods for differential equations*](https://github.com/ktwyw/math-perturbation-methods). A
second companion, [*Error analysis and propagation*](https://github.com/ktwyw/error-analysis),
treats the accuracy of numbers from the experimental side.

## Part I. Foundations

| Chapter | Track | Topics |
|---|:-:|---|
| [1. What an asymptotic expansion is](../notebooks/part1_foundations/01_what_is_an_asymptotic_expansion.ipynb) | core | Poincaré's definition; convergent vs asymptotic; the Euler series and optimal truncation; erfc in double precision |
| [2. Scaling: finding the small parameter](../notebooks/part1_foundations/02_scaling_and_small_parameters.ipynb) | core | nondimensionalisation; a projectile; regular and singular roots; dominant balance |
| [3. Order symbols, gauge functions and uniqueness](../notebooks/part1_foundations/03_order_symbols_and_gauge_functions.ipynb) | core | O, o and ~; gauge functions and uniqueness; logarithmic gauges (Lambert W); terms beyond all orders; non-uniformity |
| [4. Divergent series in practice](../notebooks/part1_foundations/04_divergent_series_in_practice.ipynb) | advanced | optimal truncation and e<sup>−A/ε</sup>; Dyson's argument; Borel–Padé summation; Padé; Shanks and Richardson |

## Part II. Integrals and sums

| Chapter | Track | Topics |
|---|:-:|---|
| [5. Integration by parts and Watson's lemma](../notebooks/part2_integrals_and_sums/05_integration_by_parts_and_watsons_lemma.ipynb) | core | integration by parts (erfc, E<sub>1</sub>); Watson's lemma; fractional powers |
| [6. Laplace's method](../notebooks/part2_integrals_and_sums/06_laplaces_method.ipynb) | core | Laplace's method with its first correction; Stirling's formula; I<sub>0</sub>; Beta integrals |
| [7. Stationary phase](../notebooks/part2_integrals_and_sums/07_stationary_phase.ipynb) | core | stationary phase; J<sub>0</sub> and the Hankel function; endpoint contributions; diffraction at a straight edge |
| [8. Steepest descent and saddle points](../notebooks/part2_integrals_and_sums/08_steepest_descent.ipynb) | advanced | steepest-descent paths; the Airy function on both sides; coalescing saddles; a Fourier transform with complex saddles |
| [9. Sums: Euler–Maclaurin and the Stirling series](../notebooks/part2_integrals_and_sums/09_sums_euler_maclaurin_and_stirling.ipynb) | core | Euler–Maclaurin; harmonic numbers and ζ(2); the Stirling series and its optimal truncation; the rotational partition function |

## Part III. Worked examples

| Chapter | Track | Topics |
|---|:-:|---|
| [10. Chemical physics: Kramers' escape rate](../notebooks/part3_examples/10_kramers_escape_rate.ipynb) | core | the exact mean first-passage time; two Laplace integrals; Kramers' formula and its correction; Brownian simulation; Poisson escape |
| [11. Heat transfer: transient conduction at short and long times](../notebooks/part3_examples/11_transient_heat_conduction.ipynb) | core | images vs Fourier modes; equal convergence at t = 1/π; one-term and semi-infinite approximations; beyond all orders at the centre |

## Where to start

| If you are… | Read |
|---|---|
| a senior undergraduate meeting the subject | Chapters 1–3, 5 and 6, then Example 11 |
| a graduate student in physics or chemistry | Part I, Chapters 6, 8 and 9, then Example 10 |
| an engineer | Chapters 1–3, 6 and 9, then Example 11 |
| a researcher who wants to check an expansion | Chapters 1, 3 and 4: error slopes, gauge functions and divergent series |

## Licence and citation

Code: MIT. Text and figures: CC BY 4.0. Please cite the material using the
metadata in `CITATION.cff` in the repository.
