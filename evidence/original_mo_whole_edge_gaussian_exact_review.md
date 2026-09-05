# Independent complete review: whole-edge Gaussian reduction and actual optimality

2026-09-05. Reviewer: optimized_profile_exact. **PASS.**

## Exact source and complete-read scope

Reviewed source:
`/tmp/original_mo_whole_edge_source_preserving_gaussian_reduction.md`

Final length: 502 lines.
Final SHA-256:
`6b22fb3ab1cc878b08fe79b5b57e0e661eaaa792dfc67f850d35db9f1b68bead`.

The complete original 318-line proof was read at SHA-256
`ae95099f76a54e13cc8b7edff5636a1a99b734ede04937079120bf6cd10ff0ce`.
The complete added Section 7 and final scope were read at the 502-line
SHA-256 `83e519418578aca2d4e0db04f304ffe7c4e2d82656d717308efe2cbc76edcec0`.
The final change is exactly the third-line status header. Restoring the
preceding header reproduces that preceding whole-file hash. No formula
or proof text changed after the completed mathematical review.

The mean-preserving shifted-sign universality prerequisite, including
its smooth-activation and Gaussian Holder ingredients, was independently
read in full in this review sequence. Its required hypotheses are
checked explicitly in the present source.

## Whole-edge reduction

1. The tensor covariance `I-K tensor K/L^2` is PSD and bounded by `2I`.
   Symmetric off-diagonal compression has diagonal `1-1/L^2`, so the
   denominator `D=L^2-1` is necessary. Matrix entries and orthonormal
   basis coefficients have the correct relative factor `sqrt(2)`.
2. The `n=2` compression is genuinely zero and cannot be normalized.
   Its independent one-edge fallback is valid. For `n>=3`, the triangle
   calculation and `tr(K^2)=n(n-1)` give `L^2>=3`, hence `R<=3I`.
3. Adjacent and disjoint edge correlations, their signs, and the entire
   shifted even-Hermite decomposition in equation (8) are exact. The
   retained rank-one covariance is PSD, and diagonal Hermite mass is
   not counted twice.
4. The four-cycle matrix calculation is complete: the ordered-pair
   square has the stated nonnegative-entry formula, the Schur-map
   argument bounds `V^T V`, and compression gives `2Q`, not `Q`.
   Therefore `||Q||<=L(n-1)/2` with no source operator cap.
5. The adjacent/disjoint row counts give equation (12). Combining the
   odd, four-cycle, line-graph and diagonal remainders gives equation
   (13), including the displayed `40/sqrt(n)` bound and small orders.
6. The scalar retained Gaussian has the claimed actual quadratic-norm
   cost. Even the trivial source norm bound suffices for the stated
   `O(n)` estimate; the argument is not a rank-only norm assertion.
7. Gaussian convex order with the augmented Boolean states gives the
   `O(n^(5/4))` covariance-to-norm error, uniformly in arbitrary fixed
   internal offsets. Shifted-sign universality gives `O(n^(16/11))`.
   There are no matching diagonal entries to restore in Sections 1--6.
8. Every threshold outcome is a complete signing of the same order.
   Consequently exact original norm optimality implies equation (17)
   in its stated LOWER direction. No opposite-order conclusion is used.

## Full symmetric lift and the integrated actual constraint

1. The full symmetric covariance operator in equation (18) has the
   correct factor two. Its actual off-diagonal entries extend exactly
   the previously defined edge Gaussian. Equation (19) follows by
   contracting against `xx^T/2` and `yy^T/2`.
2. The only contribution of the added matrix diagonal to a Boolean
   quadratic energy is the common scalar `tr(H)/2`. Its variance is
   exactly equation (20), at most `n/2`. Independence from off-diagonal
   entries is unnecessary. The combined lift therefore costs at most
   `sqrt(n/pi)` in expected absolute norm.
3. A common centered scalar cancels from expected ONE-PHASE pressure,
   but not from the augmented absolute phase. The source correctly
   retains the augmented posterior for all subsequent calculations.
4. The two replica expressions in equation (22) use the same current
   posterior and its factors `sigma sigma'`. The covariance contractions
   in equation (23) are exact, nonnegative, and account for the formerly
   present off-diagonal-only `K^2` terms through the explicit lift.
5. Frobenius Cauchy--Schwarz gives equation (24) even for a signed,
   indefinite `Gamma`. Substitution gives equation (25), including both
   its negative current energy-square term and its overlap subtraction.
   This inequality has not been given the wrong direction after a
   negative coefficient multiplies a variance.
6. Gaussian covariance differentiation gives equation (26) with the
   displayed signs and factors: `k'=-2hk`, `v'=2hk-2sa`, and the source
   matrix and latent law remain fixed. No current moment is replaced
   by an uncoupled or source Gibbs moment.
7. The scalar coefficient derivatives are integrable, the finite-state
   posterior observables and covariance contractions are bounded at
   fixed order, and the Gaussian noise vanishes as `h` tends to infinity.
   Hence the endpoint identity and integral in equation (27) are valid.
8. Equations (17), (21), the finite maximum-term lower bound, and
   `P_K<=beta Phi(K)` give equation (28) for an EXACT ORIGINAL NORM
   minimizer. No pressure-minimizer assumption is inserted. Choosing
   `beta=n^(-5/11)` gives equation (29) with the stated uniform absolute
   error. There is no integration of a pointwise approximation error
   over an infinite interval.

## Scope gate

The source proves a universal whole-edge shifted-sign Gaussian reduction,
a universal negative-energy-square current-posterior variance upper, and
an integrated same-order constraint from actual original norm optimality.
It does not sign the derivative pointwise or identify the integrated
constraint with the missing doubled-order upper bound. That transport
step and the original MO convergence statement remain open.
