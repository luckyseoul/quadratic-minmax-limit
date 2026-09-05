# Collaborative-author verification: tensor mixture and sign defects

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no mathematical changes requested.

## Exact input and review status

All 351 lines of `/tmp/original_mo_cross_tensor_mixture_sign_defect.md`
were read completely. The reviewed SHA-256 is
`66ed4d02f4cdd7d323ae5f0c717993bc3453d59dcb3a1e00b2e2c64720d424da`.
The input was left unchanged.

This is COLLABORATIVE-AUTHOR VERIFICATION, not a globally independent
full-note review. The docs worker contributed the tensor-mixture
construction, the interior Taylor estimate, and the observation that
the general interpolation lemma sharpens the cap-free remainder.
The exact worker contributed the two-sided sign-defect improvement
and wrote the integrated proof. Each component was cross-checked,
but an independent complete reviewer is still a separate provenance
requirement for the integrated artifact.

## Mathematical checks

The canonical Gram has unit diagonal and gives the exact third and
fourth singular-moment identities. Its SDP gap is nonnegative, with
no inference that actual conditional optimality makes that gap small.
The optimal diagonal is positive and has balanced block traces.

The weighted residual estimate is cap-free. The orthogonal polar
factor transfers the column residual with exactly the printed weight
orientation, including deficient rank. Applying the scalar negative
part inequality separately to rows and columns and adding gives
`N_- <= g/2`, not `g` or `g/4` in the general case.
The scalar-optimal-diagonal spectral calculation gives the distinct
improvement `N_- <= g/4` with its stated additional hypothesis.

The odd-tensor lifts have squared norms `1-t`; adding the canonical
families with weight t makes actual unit vectors. Their cross Gram
is exactly the printed sine-plus-canonical mixture. The countable
tensor notation does not require an infinite-dimensional Gaussian,
since only a finite Gram needs realization.

For every fixed `t<1`, `sin(asinh(1-t))+t<1`. The entire Taylor
segment lies in that strictly interior interval, and the absolute
second derivative gives exactly the factor `t^2 M_t/2`. Both the
generic gap coefficient `(1+sec a)/2` and the scalar coefficient
`(3+sec a)/4` follow with the correct sign. The t=0 bound recovers
the elementary tensor-rounding constant with zero remainder.

The real/complex interpolation factor in Section 4 is correct.
Combining it with the canonical primal gives the normalized fourth
moment estimate (17), and the Boolean cap alone gives (18).
The absolute error is `O_(C,t)(n^(5/4))`; no fixed conference-scale
operator bound is needed or imported.

The normalized inequality and gap rearrangement in (19) have the
right directions. The limiting envelope keeps t fixed before taking
the large-n limit. Reparameterization by `t=1-sinh(a)` gives
`F_U(0)=1-U` and `F_U'(0)=1-U`, proving the asserted strict gain
when `c_0<=U<1`. The scalar denominator has the same endpoint
derivative. For `U>=1` the envelope is zero. The small-t endpoint
constants in (21) and their required order of limits are correct.

The shell qualification is important and correctly stated: the
ratio in this theorem is `beta/tau`; an inactive shell cannot
substitute `|c|` for beta. The note does not identify arbitrary
moment laws with complete signings, assert scalar duals generally,
or conclude the unresolved original convergence theorem.

## Procedure

This was an analytic full read with file hashing. No signing search,
scalar numerical evaluation, solver, simulation, or other tool-based
mathematical computation was performed for this verification.

