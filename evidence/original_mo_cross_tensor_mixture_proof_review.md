# Independent full proof review: cross tensor mixture and sign defects

2026-09-05. Reviewer: optimized_profile_proof.

Reviewed source:
`/tmp/original_mo_cross_tensor_mixture_sign_defect.md`, 351 lines,
SHA-256 `66ed4d02f4cdd7d323ae5f0c717993bc3453d59dcb3a1e00b2e2c64720d424da`.
Every line was read. **PASS, with no corrections requested.**

The reviewer did not contribute to the tensor-mixture derivation. Before
receiving the note, the reviewer independently derived the two sign-defect
constants as a preliminary audit; the complete source was then checked
independently, including its distinct general weighted-residual proof.

## Actual SDP and residuals

The canonical families are unit and have the printed cross Gram and
cubic objective. The SDP dual is actually optimal, its diagonal entries
are positive, and its two block traces balance. No small canonical gap
or scalar diagonal is inferred from signing optimality.

For `Q=D-H_B`, the contraction argument gives
`Q D^{-1} Q <= 2Q` and the weighted residual bound `4g`. The column
orientation is correct: with an orthogonal polar factor, `R_c=U R_2^T`,
so `||R_c D_c^{-1/2}||_F=||D_c^{-1/2}R_2||_F`, including when B is
rank deficient. Applying the elementary negative-part inequality to
both residual blocks and adding gives precisely `N_-<=g/2`.

In the scalar case, actual dual optimality gives `tau=nd`; mere scalar
feasibility would not suffice. Since every singular value is at most d,
the displayed squared residual is at most `(d/n)g`. This yields the
stronger `N_-<=g/4`, without asserting that the canonical gap is zero.

## Tensor mixture and finite-order constants

For fixed `0<=t<1`, the tensor lifts have squared norm `sinh(a)=1-t`.
The orthogonal canonical component contributes squared norm t, giving
genuine unit vectors and cross Gram `sin(aV)+tW`. Finite Gram realization
legitimizes ordinary finite-dimensional Gaussian sign rounding.

The radius `r=sin(a)+t` is strictly less than one. Every Taylor segment
lies in `[-r,r]`; the absolute second derivative is bounded by the
printed `M_t`. The first derivative is `sec(aV)` with no sign ambiguity,
because `|aV|<pi/2`. The sign-defect estimate has the correct direction
on negative entries and produces exactly `(1+sec(a))/2`, or
`(3+sec(a))/4` for a scalar optimal diagonal.

The quadratic Taylor error retains `tr|B|^4/n^2` with the correct factor
`kappa t^2 M_t/2`. At t=0 it vanishes and recovers the exact elementary
tensor-rounding constant; t=1 is explicitly not substituted.

## Norm-only error and the limiting curve

The real/complex interpolation proof of `||B||^2<=2 beta(B)` has the
correct factor two. Combined with the canonical-primal inequality,
it gives `tr|B|^4/(n^2 tau)<=L/n<=sqrt(2 beta)/n` exactly. Thus the
normalized error is `O_(C,t)(n^(-1/4))` and the absolute error is
`O_(C,t)(n^(5/4))` under only the Boolean cap.

The normalized inequalities and the domain `c_0<=U<=1/kappa` are
correct. For each fixed t the error vanishes before taking the
supremum; no uniform Taylor control at t approaching one is required.
After reparameterization, the quotient has endpoint value `1-U` and
right derivative `1-U`. Hence its envelope is strictly larger than
the cubic-only bound for `c_0<=U<1`. For `U>=1` every numerator is
negative and the zero truncation is correct.

At `U=c_0`, expanding a at t=0 gives the exact two constants in (21).
The order of these limits is explicitly valid. The excluded joint
profile concerns the stated limiting SDP/Boolean and canonical-gap
variables, not a proof that every remaining abstract profile is an
actual signing.

## Scope and procedure

The active-shell restriction is retained: beta may be replaced by a
shell magnitude only when that magnitude actually attains beta.
Nothing identifies an arbitrary actual SDP diagonal with a scalar one,
evaluates the full conditional Gaussian maximum, or proves convergence.

This review is analytic. File reading, line counting, and SHA-256
verification were used; no mathematical numerical computation, signing
census, simulation, or solver run was performed for this review.
