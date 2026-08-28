# Exact p=11 profile and quartic-trace theta reconstruction

Date: 2026-08-28. This is Proposition 15.667. It gives an exact finite
profile computation of the ordinary and quartic-trace theta coefficients at
`p=11`, uniquely reconstructs both modular forms through exponent 800, and
tests the full scalar trace-conservation cone with rational primal/dual LP
certificates. It does **not** prove R1.

## Normalization erratum for Proposition 15.665

The first version of Proposition 15.665 used the degree-four values from
Propositions 15.633--15.635 at `u/2` for the scaled-norm 20 and 24 shells,
but combined them with the radial correction for the unscaled shell vector
`u`. Since the polynomial is quartic,

\[
 H(u)=16H(u/2).
\]

The corrected p=11 traces are

| scaled norm `e` | shell count `N_e` | harmonic trace | raw trace `tau_e` |
|---:|---:|---:|---:|
| 11 | 244 | `-3538/63` | `0` |
| 20 | 16,104 | `-85888/21` | `89792/11` |
| 24 | 14,762 | `-63684/7` | `7076` |
| 27 | 442,860 | `-527406/7` | `538752` |

The general identity `A_e=R_e-rho_e I`, positivity of `R_e`, and conserved
raw trace in Proposition 15.665 are unchanged. Only the two numerical audit
rows and modular affine anchors derived from them required correction.

The profile computation below independently obtains the corrected values
from coordinate fourth moments, without using the two old channel spectra.

## Exact finite profile reduction

The ten-dimensional glue-dual code has `11^10=25,937,424,601` words. Plane
translations and nonzero scalar conjugacy reduce it to 21,437,340 exact
weighted representatives in four highest-degree strata. For every
representative, the six direction polynomials are replaced by their eleven
value multiplicities. There are 604 quartic value-distribution types and
2,558,543 weighted sorted six-tuples after exact aggregation.

For one direction profile `a=(a_s)`, integer dynamic programming records

\[
 \sum_s a_s=b,\qquad \sum_s a_s^2=b+2k,
 \qquad \sum_s y_sa_s=h\pmod {11}.
\]

The 604 distributions fall into 13 affine output types under
`y -> alpha*y+beta`. Tables for those 13 representatives reconstruct all
604 tables exactly; sampled tables are independently recomputed without
using that affine reconstruction. The complete tuple weights reconstruct
all `11^10` glue words.

For the six line profiles with common sum `t`, the corresponding dual vector
satisfies

\[
 x_\infty={t\over2p},\qquad
 x_u={\sum_{j=1}^6 a_{j,t_j(u)}-t/2\over p},\qquad
 e=2p\lVert x\rVert^2=2\sum_{j,s}a_{j,s}^2-t^2.
\]

The exact convolution accumulates the ordinary coefficient `N_e` and the
common-sum moments

\[
 M_{2,e}=\sum_t t^2N_{e,t},\qquad
 M_{4,e}=\sum_t t^4N_{e,t}.
\]

Coordinate transitivity gives

\[
 \sum_{x\in X_e}\sum_i x_i^4={p^2+1\over16p^4}M_{4,e},
\]

and Proposition 15.665's projection formula therefore yields

\[
 \boxed{\tau_e={p^2+1\over4p^2(p^2-1)}
              \left(N_e e^2-M_{4,e}\right).}                 \tag{1}
\]

The independent tight-frame check is

\[
 \boxed{(p^2+1)M_{2,e}=2peN_e.}                              \tag{2}
\]

Five 31-bit primes admitting primitive eleventh roots were used. Their
product is

`31999921744068749461247094447450713426945936557`,

which exceeds the unrestricted bounds for the count, second moment, and
fourth moment at every exponent through 120; the largest is

`4129286721131182422173790859253254080`.

Thus CRT reconstruction is unique over the integers. The V100 accelerates
modular arithmetic, but no floating-point result enters the certificate.
All 51 nonempty exponents through 120 satisfy (2), and all previously known
ordinary coefficients and raw traces match. The four calibration rows are:

| `e` | `N_e` | `M_2,e` | `M_4,e` | `tau_e` from (1) |
|---:|---:|---:|---:|---:|
| 11 | 244 | 484 | 29,524 | `0` |
| 20 | 16,104 | 58,080 | 2,555,520 | `89792/11` |
| 24 | 14,762 | 63,888 | 5,134,272 | `7076` |
| 27 | 442,860 | 2,156,220 | 66,363,660 | `538752` |

## Two exact modular reconstructions

The ordinary scalar theta series lies in the previously exported
41-dimensional exact modular space. Its profile coefficients through
exponent 88 have rank 41, so they uniquely determine the modular form. The
half-cusp target is then an independent prediction, not an extra anchor.
The resulting form has nonnegative integer coefficients through exponent
800, and all 32 held-out profile coefficients at exponents 89--120 match.

After recomputing the three harmonic channel reductions with the corrected
unscaled anchors, their homogeneous 32-column q-row matrices agree exactly.
The dimension-weighted base gives the trace series. Its exact profile prefix
first reaches rank 32 at exponent 92, with pivot exponents

```
31,32,35,36,39,40,43,44,47,48,51,52,55,56,
59,60,63,64,67,68,71,72,75,76,79,80,83,84,87,88,91,92.
```

All 28 held-out trace coefficients at exponents 93--120 match. Combining
the reconstructed harmonic trace with the scalar series and the exact
radial shift gives nonnegative raw traces through exponent 800. The exact
dimension-weighted half-cusp trace target is

\[
 -{4428472046531859136727844588716
    \over11389083011948997399715094836557}.
\]

These are coefficient-determination theorems for the scalar and trace
series. They do not determine how each shell's positive raw mass is split
among the individual PSL constituents.

## Exact trace-conservation cone through exponent 800

For each of the four p=11 component cases, every exponent through 800 was
constrained by nonnegative raw constituent mass and the exact conserved
total `tau_e`. After exact row conditioning, QSopt_ex solved both endpoint
problems. The checked solution parser independently verifies every rational
primal inequality, every dual sign, every stationarity equation, and equality
of primal and dual objectives.

| distinguished component | exact feasible target interval (decimal display) |
|---|---:|
| circle-kernel principal | `[-651.3544743742135, 614.4812457395238]` |
| circle-low Weil | `[-888.4183023833558, 874.9202223929279]` |
| circle-low principal | `[-651.2814664914692, 614.5489396383163]` |
| circle-high principal | `[-651.2668649149203, 614.5624784180748]` |

The exact rational endpoints are in the proposition JSON and archived LP
reports. Comparing independently certified truncations at exponents 120 and
800, seven of the eight endpoints are identical. The sole late contraction
is the circle-low-Weil maximum, from `880.0044032421115` to
`874.9202223929279`.

This is a genuine contraction relative to the uncoupled coefficient cone,
but it is decisive bad news for this particular closure route: even 800
exact shell constraints leave broad two-sided intervals. Aggregate scalar
trace conservation alone does not close R1.

## Literature and OEIS duplicate/context checks

Targeted searches after the result used the corrected traces `89792/11` and
`7076`, the exact trace half-cusp numerator, the phrases “profile theta
series”, “common coordinate fourth moment”, and Paley/conference quartic
theta combinations. They found the established adjacent theories of
weighted theta series and shell designs (for example
[arXiv:2308.14309](https://arxiv.org/abs/2308.14309)), and the standard
complete-weight-enumerator route from codes to Construction-A theta series.
No search result stated this Paley glue-profile reduction, the two finite
rank thresholds, or the shellwise conserved-trace LP.

Individual exact OEIS searches for

`526682934302285507648`, `7486058827792258651296`,
`33436768468462410480`, `46157661685482511780`, and
`2041917404166462280936`

returned no entries. Earlier searches of five large ordinary coefficients
also returned no matches. These are duplicate and context checks only; no
integer-sequence novelty or priority claim is made.

## Artifacts and next attack

The compact theorem record pins SHA-256 hashes for the finite-field orbit
report, tuple and profile-table archives, five-prime moment report, scalar
and trace reconstructions, corrected exact q-rows, and both QSopt endpoint
reports. Reproducibility artifacts are backed up under

`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/`.

The next structurally stronger attack is channel-resolved. Let `O` be the
square-circle tensor operator, whose three broad eigenvalues are already
known. In addition to `tr(R_e)`, exact profile moments for
`tr(R_e O)` and `tr(R_e O^2)` would solve separately for the kernel, low,
and high raw masses on every shell. That adds information which the scalar
conservation cone provably lacks. A still finer alternative is to compute
twisted traces `tr(R_e pi(g))` for selected PSL classes.
