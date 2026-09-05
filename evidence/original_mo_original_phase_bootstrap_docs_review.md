# Independent review: original-phase norm-only spectral bootstrap

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS, without requested mathematical changes.

## Exact reviewed inputs

The reviewer read every line of both proof files:

- Prerequisite: `/tmp/original_mo_original_phase_spectral_moment.md`,
  262 lines, SHA-256
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- Supplement: `/tmp/original_mo_original_phase_spectral_moment_bootstrap.md`,
  176 lines, SHA-256
  `3736db69d904b5a63ade46b32f6fddcc0505019f45ef483110c3ee67b24c8915`.

The reviewer additionally read the published direct-normalization proof,
including its complete Section 3 argument for the cited estimate
`||A||_op^2 <= 8 Phi(A)`. In particular the real/complex interpolation
factor and zero-diagonal cube-polarization factor were checked, not
merely accepted from the supplement's citation.

The reviewer did not contribute to the derivation of the prerequisite
or this bootstrap supplement. This is an independent mathematical
review, distinct from the root's disclosed collaborative-author check.
The proof inputs were left unchanged.

## Prerequisite checks

The positive and negative covariance matrices in (4) are genuinely PSD
and have diagonal one. Each phase uses its own admissible variance.
The linear terms in (7) have the correct respective signs and no
missing factor of two. The arcsine remainder in (8) bounds only the
off-diagonal correlations before using the spectral Frobenius norm;
the extra diagonal padding is not spuriously added to that estimate.
This proves both one-sided inequalities (2) and their sum (3).

The smallest admissible phase variances are strictly positive. Their
sum is exactly `2q + osc_i(A|A|)_ii`, even if the diagonal imbalance
is a nonzero constant. The finite additive error in (9) is correctly
qualified in the source note: the source alone does not declare it
subleading under an arbitrary operator norm.

The local scalar-moment argument in Section 3 was checked in full.
Both positive and negative masses and their common first moment are
nonzero. The probability estimates, the monotonic substitution
`t >= P_2/a`, and the final bound on `N_2` have the asserted directions.
The lower row-moment bound follows by applying the same argument to
`-X`. The spectral support identity implies `ab >= q`, but does not
require each extreme separately to exceed `sqrt(q)`. The spectral
variance choices and the factor between one and two follow correctly.

The nuclear-norm proof uses the same strictly positive diagonal for
the two actual PSD phases. The signed arcsine difference has the right
orientation also on negative edges. The bound between the two expected
quadratic energies costs exactly `2 Phi(A)`. The unordered-pair
Cauchy--Schwarz calculation gives
`Phi(A) >= kappa n^2 q / (2 tr|A|)` with the printed constant.

## Supplement checks

Equations (2)--(3) correctly use `S_4^+ <= a S_3^+` and its negative
analogue. They are valid regardless of the sign of the coefficient,
and division is explicitly restricted to positive coefficients.
From `ab >= q` and `L=max(a,b)`, both extremes lie in `[q/L,L]`.
Consequently both spectral coefficients are at least
`gamma_sp = kappa - 2 delta L/q`, proving (5) when positive.

The original-norm hypothesis gives `L <= sqrt(8C) n^(3/4)`, so
`gamma_sp = kappa - O_C(n^(-1/4))` uniformly. The leading universal
coefficient in (8) is therefore justified without an unproved
conference-scale operator bound. It is not described as a strict
improvement of the universal coefficient `4/kappa`.

For the actual variances in (9), trace zero gives each phase nuclear
mass `S_1/2`. Cauchy--Schwarz over at most n eigenvalues and the maximum
versus average diagonal estimate yield
`v_pm^* >= 2 S_2^pm/n >= S_1^2/(2n^2)`.
This implies the exact common coefficient (11), and hence (12),
with the stated positivity condition.

The prerequisite nuclear inequality implies
`S_1 >= kappa sqrt(n) q/(2C)` under the norm cap. Substituting this
and the operator bound into (11) gives exactly the error bound (14),
including its factor `16 delta C^2 sqrt(8C)/kappa^2` and power
`n^(7/4)/q^2`. Thus the coefficient error is uniformly
`O_C(n^(-1/4))`.

Once `gamma_* >= kappa/2`, equation (12) bounds
`S_3/V^* <= 2 Phi(A)/kappa`. Replacing `gamma_*` by `kappa` therefore
costs at most `O_C(n^(-1/4)) Phi(A) = O_C(n^(5/4))` in the ORIGINAL
quadratic norm. This justifies (15)--(16) with the actual diagonal
oscillation denominator, not merely the weaker spectral denominator.

## Scope and procedure

The supplement supplies a new bootstrap of the reviewed source, not a
correction or a retroactive reinterpretation of its fourth-moment
estimate. It preserves the same original signing and does not require
source regularization, cross-SDP scalar diagonals, small diagonal
oscillation, favorable singular-vector alignment, or spectral symmetry.
No evaluated paired Gaussian upper or all-orders convergence is claimed.

This review was analytic. File reading, line counting, and SHA-256
verification were used. No signing search, numerical scalar evaluation,
simulation, solver, or other tool-based mathematical computation was
performed for this review.

