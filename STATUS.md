# Current mathematical status

Updated 2026-09-05. The original MO limit is OPEN.

There is no reviewed proof of convergence or nonconvergence. The value
`1/2` is unproved. These statements are independent of the status of the
optional Paley research program.

## What is established

`CORE.md` gives the route-neutral definitions and proofs. In particular,
`m_n` is nondecreasing, consecutive `alpha_n` differences tend to zero,
the limit-point set is an interval, and

`1/pi <= liminf alpha_n <= limsup alpha_n <= 1/2`.

The upper bound uses conference constructions, not their optimality.
Ratio-dense transfer and two-multiplier Dini amplification are valid
conditional tools; neither supplies its missing hypothesis automatically.

The new [selected-half restriction theorem](evidence/NOTE_2026-09-05_SELECTED_HALF_RESTRICTION.md)
also proves `m_(2n+1)>=2m_n` and `m_(2n)>=2m_n-(n-1)`.
These unconditional bounds use deterministic phase coloring and selected
subsets of actual larger-order minimizers. They improve ordinary
restriction monotonicity but still leave a normalized factor `sqrt(2)`;
they do not establish sharp proportional restriction or convergence.

## What changed at the reset

Residual (ii) is parked as an optional route-local open lemma. Its former
status as an unavoidable obstacle was a bookkeeping error. The same applies
to treating an E1/bi-tight conjunction as an if-and-only-if test for the
original limit. The gap-two condition is stronger than the asymptotic
optimality needed to transfer the value `1/2`; neither has been proved
necessary for convergence.

The active global proof registry is `src/original_mo_status.py`.
`src/e1_main_chain_status.py` retains separate optional-route diagnostics.
Historical proof notes, scoped local theorems, and exact certificates remain
available; unproved bridges have not been promoted or declared false.
See `ARTIFACTS.md` for the preserved branches and terminology.

## Fresh mathematical target

The [full-strength boundary identity](evidence/NOTE_2026-09-05_FULL_STRENGTH_BOUNDARY_LIKELIHOOD.md)
retains the exact pair-dependent Gaussian boundary law and the entire
posterior covariance integral, including its singular endpoint limit.
An analytic order-three actual minimizer disproves a coordinatewise
sign premise, not the averaged comparison; the offending boundary
context is absent at the endpoint itself. The full-strength weighted
comparison remains open. Complete independent proof reads passed.

The [fixed-strength weak-rounding theorem](evidence/NOTE_2026-09-05_FIXED_WEAK_GAUSSIAN_CUSHION.md)
extends the vanishing-strength exclusion to a NONEMPTY interval of fixed
positive strengths. Conditional Gaussian replacement, removal of the
conditional mean by convexity, and a variance interpolation give the
actual pressure floor
`[c sqrt(2t)K0-2log2-c^2 t arcsin(rho)/(2pi)]n-o(n)`.
When its excess over `cn` is positive, the actual centered law has
exponentially small probability of any successful internal host, even
one selected after the cross block. Uniform Gaussian concentration covers
every generating host, and `exp(o(n))` such proposals fail without an
independence assumption. This does not use low sign-law information and
does not exclude the singular full-strength `rho=1`. The complete proof
passed independent review; no new computation or sampling was used.

The [Gaussian-sign information theorem](evidence/NOTE_2026-09-05_GAUSSIAN_SIGN_INFORMATION_SCALE.md)
separates two dependent-rounding regimes. For any sign law with second
moment matrix `C`, its relative entropy from iid signs is at least
`||C-I||_F^2/(4||C||_op)`. Gaussian arcsine and Schur-product identities
extend this to singular Gaussian-sign laws. For the actual centered
construction `Sigma_rho=I+rho H/mu`, every fixed `rho>0`, including the
singular canonical `rho=1`, has `Omega(n)` information on norm-capped
hosts. Therefore the low-information exclusion below does not apply to
full-strength rounding. Conversely, `rho=o(n^(-1/2))` has `o(n)`
information on actual half-product minimizers and fails the finite-step
comparison in mean and probability at the stated fixed parameters.
These are independently reviewed all-orders results, not a proof of
success or failure of the full-strength law or of the original limit.

The [iid quenched cross-block theorem](evidence/NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md)
now controls the ACTUAL average log pressure, not only annealed proxies.
An explicit admissible Parisi control gives `P_SK>=4/(3sqrt(pi))`;
Gaussian covariance comparison and a direct Bernoulli replacement prove
a host-free pressure floor. At the stated sufficiently large fixed `c,t`,
even choosing the internal host after seeing the cross block cannot make
`exp(o(n))` iid-marginal proposals succeed. Any successful proposal law
needs relative entropy `Omega(n)` from independent fair cross signs.
This does not exclude unrestricted dependent selection or the earlier
full-strength Gaussian-sign law. These are independently reviewed
all-orders proofs, with no new simulation or finite census.

The [positive-cone cutoff theorem](evidence/NOTE_2026-09-05_POSITIVE_CONE_TRUNCATION.md)
proves that retaining only degrees `k<=o(N)` loses a positive normalized
pressure at fixed `c>pi log2`, including actual optimizers. It does NOT
exclude a sparse selection of extensive degrees. The
[extensive coefficient comparison](evidence/NOTE_2026-09-05_EXTENSIVE_COEFFICIENT_MOMENTS.md)
identifies the full coefficients with high even moments up to
dimension-uniform factors. The
[degree-selector identities](evidence/NOTE_2026-09-05_POSITIVE_DEGREE_SELECTOR.md)
retain minimization/maximization order; a genuine mixed exchange costs
only `O(log N)` in pressure, while a pure exchange has the stated spin
entropy cost. Their scalar coefficient-rate convergence would suffice
for the original limit, but no cross-order transport is proved.

The [exact optimized order-six profile](evidence/NOTE_2026-09-05_EXACT_OPTIMIZED_ORDER_SIX_PROFILE.md)
is an independently reviewed analytic theorem for every temperature and
every point of the balanced `3+3` path. A coefficientwise polynomial
comparison gives its global minimum over ALL complete order-six signings;
for positive cross weight, equality holds exactly at `A^2=5I`.
Its endpoints cross once as temperature varies, but its interior excursion
above the starting endpoint grows linearly in `c` as `c->infinity`.
This excludes a temperature-uniform bounded excursion, not a fixed-`c`
small-oh comparison as the order grows. The candidate polynomial is
quadratic only at this order; no all-orders higher-coefficient comparison
has been established. Independent NUKA integer certificates and V100
histograms agree exactly, with finite evidence recorded separately in
`evidence/original_mo_optimized_profile_mesh.json`.

The [finite-step rounding and annealing theorem](evidence/NOTE_2026-09-05_FINITE_STEP_ROUNDING_ANNEALING.md)
is proved and independently reviewed. It extends Gram--Schmidt rounding
to an integral finite-step pressure bound, keeping the exponential inside
the actual spin expectation. For fixed `c>0` and `0<t<=1`, every such
quadratic moment-generating-function proxy is at least `c^2 t n/4`.
Separately, the actual canonical Gaussian-sign annealed pressure is at
least `c^2 t n/(2pi)-o(n)`, uniformly over complete hosts. Against the
optimized paired endpoint `2R_n<=cn+o(n)`, these certificates fail at
`c>4/t` and `c>2pi/t`, respectively. This does not exclude a good selected
cross block, the average log pressure, the actual Gram--Schmidt law,
vanishing step sizes, or another order comparison. The finite mesh check
is recorded separately from the all-orders analytic proofs.

The [integral cross-block rounding theorem](evidence/NOTE_2026-09-05_INTEGRAL_CROSS_BLOCK_COVARIANCE_ROUNDING.md)
is proved and independently reviewed. For every complete signing with
actual Gibbs covariances at opposite temperatures, it selects an integral
cross block with phase-averaged squared bilinear cost
`min_B qbar(B)<=n^2-8(a_A')^2/(pi ||A||_op^2)`.
The Gaussian construction uses the exact negative spectral edge of a
centered tensor matrix and an entrywise arcsine identity. Separate
covariance-rounding bounds retain fixed coordinate squares and spectral
tails; a stated low-effective-rank hypothesis suffices for the local
endpoint comparison. That hypothesis is not established for minimizers.
The general radial comparison and the integrated order comparison both
remain unproved. These are integral construction bounds, not convergence.

The [near-minimizer counterfamily](evidence/NOTE_2026-09-05_NEAR_MINIMIZER_OPPOSITE_PHASE_COUNTERFAMILY.md)
is proved and independently reviewed for actual complete quadratic hosts.
For every fixed `c>0` and all sufficiently large `N`, it gives
`tr(A U A V)=Omega_c(N^(9/4))` at `beta=c/sqrt(N)`, with either norm excess
`Phi(A)-m_N=O_c(N^(11/8))` or half-product pressure excess
`a_A-R_N=O_c(N^(7/8))`. These are two potentially different families.
Thus even leading-order near-minimality does not imply the proposed
quadratic trace bound. The constructed hosts fail edge-local half-product
optimality; exact global-minimizer and actual edge-local trace bounds are
not refuted. This is not a cross-order comparison or a convergence result.

The [full-row cavity counterexample](evidence/NOTE_2026-09-05_FULL_ROW_CAVITY_COUNTEREXAMPLE.md)
shows that even the complete family of row-replacement inequalities does
not imply bounded tilted second or fourth moments for arbitrary symmetric
cavity measures. Its measure has full support and a strict minimizing row.
It has not been realized as an actual quadratic Ising cavity; it is not
a counterexample to a signing, a pressure comparison, or convergence.
This rules out proving the missing moment bound from row optimality alone.

The [adaptive-perturbation theorem](evidence/NOTE_2026-09-05_ADAPTIVE_PERTURBATION_CORRELATIONS.md)
extends the signed correlation estimates to actual noise-adaptive
pressure minima. With the perturbation held fixed during sign flips,
even bounded Gaussian noise preserves signed Frobenius diffuseness in mean
and probability, uniformly along the balanced profile, by a Boolean-energy
norm refinement. This does not assert bounded-noise row-operator control.
The independently reviewed
[Gaussian switch-measure theorem](evidence/NOTE_2026-09-05_OPTIMIZED_GAUSSIAN_SWITCH_MEASURE.md)
retains the negative optimizer-switch measure in the weak Hessian and
bounds its entire Gaussian-weighted trace by `O_c((epsilon+epsilon^2)N)`.
These remove two genuine obstacles to optimized Gaussian smoothing.
They do not turn a change of deterministic block coefficients into a
heat derivative; dividing by the noise scale can destroy the small error.
The exact mean-profile derivative retains its deterministic term, and
the stochastic sign-flip reset is the previous defect up to `O_c(1)`.
No new computation was used for these two analytic results.

The [global-optimizer variational theorem](evidence/NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md)
is proved and independently reviewed. Sparse rounding controls the whole
shrinking near-sign class; global norm minimizers admit one common
near-maximizing signed-state ensemble with small aggregate correlations.
Actual edge-local symmetric-pressure minima also have small signed Gibbs
correlations, uniformly along the specified two-block weight profile.
The exact optimized-path identity isolates a signed defect imbalance.
Its needed integral bound is unproved; endpoint equality is not required
for the correctly directed dyadic subadditivity implication. These are
structural theorems, not an all-orders comparison or a convergence proof.

The [induced-restriction theorem](evidence/NOTE_2026-09-05_INDUCED_OPTIMIZER_RESTRICTIONS.md)
is proved and independently reviewed. For `n -> infinity` with
`n^2 = o(log N)`, every `n`-vertex signing occurs in every source signing
of norm `O(N^(3/2))`. Nevertheless, a uniform restriction typically has
normalized norm at least `(2/3)*sqrt(2/pi) = 0.531923...`; for globally
optimal sources, `exp(o(n))` uniform-marginal samples cannot preserve the
source's leading constant with probability bounded away from zero.

This proves that smaller-order optimal restrictions exist, not that their
constant matches the source. The latter remains precisely an unproved
cross-order comparison. The theorem excludes only the stated sampling
method and scale, not biased selection or other order comparisons.

Compare different orders using the actual global-minimizer property,
without assuming a Paley model, a limiting value, or a prescribed lift.
Any claimed convergence proof must control the normalized optimum, not just
an arbitrary low-norm signing or a typical induced restriction.

No all-orders comparison with sufficient error control is established.
Finite exclusions and moment identities do not change that status.
