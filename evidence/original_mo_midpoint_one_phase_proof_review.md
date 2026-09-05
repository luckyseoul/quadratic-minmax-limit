# Independent complete-read review: midpoint and one-phase Gaussian notes

2026-09-05. Reviewer: optimized_profile_proof. Verdict: **PASS** for both
frozen sources below, with no mathematical corrections requested.

## Sources completely read

1. `/tmp/original_mo_universal_spectral_midpoint_covariance.md`
   (424 lines), SHA256
   `1fc6f5bbb69038b6ac4ed845d26e0724a0ceb0b5a9d96d01b4554a8e37e6f968`.
2. `/tmp/original_mo_one_phase_gaussian_variance_upper.md`
   (236 lines), SHA256
   `1646f57b060db7fdaf15c2cc8a8766806d2f00297c6749e236d8e814e467bae0`.

The complete final midpoint source was reread after its final three-line
documentation/LaTex update. The complete final one-phase source was also
reread after clarifying that W need not have rank one and defining Z_A,a_A.
These are full mathematical reads, not checks of summaries or hashes alone.

## Independently checked: universal midpoint

- The trace argument gives ab>=n-1. The full admissible alpha interval is
  [-b,a], including exclusion of any mu<0 alternative by the incompatible
  same-corner conditions.
- Expansion of the P tensor P and Q tensor Q factors yields exactly the
  stated correlation matrix. The extreme spectral corners give its exact
  operator norm. Differentiating the two corner values proves the unique
  operator-minimizing midpoint, with covariance operator norm at most two.
- The alternative contraction-matrix formula, and all constants in the
  uniform O(1/n) arcsine linearization, check directly. In particular
  2|alpha_*|^3 L<=mu_*^2 is valid.
- The Gaussian covariance derivative retains the actual posterior Hessian.
  The free-alpha derivative has the correct signs and factors; neither of
  its two nonnegative traces controls the sign of their difference alone.
- The correlated-sign theorem is applicable with the absolute covariance
  bound two for every complete source, including exact original-norm
  minimizers. No source operator regularization is needed here.
- The auxiliary choices c=n^(1/22), epsilon=n^(-1/11) give the asserted
  absolute O(n^(16/11)) expected maximum error. The argument is uniform in
  the arbitrary fixed internal energy, theta, and rho. The prior can depend
  on the auxiliary temperature but not on the cross disorder.
- The resulting one-sided m_(2n) inequality follows because every rounded
  block is an actual complete signing. The unproved Gaussian upper target,
  and the additional all-orders issue beyond an unspecified dyadic little-o
  error, remain expressly separate.

## Independently checked: one phase and actual variance

- Swapping x,y in the negative phase gives f(-Z^T), with the same law under
  the stated transpose invariance. The correlated-pair absolute-difference
  bound uses at most twice the standard deviation, so the displayed phase
  reduction constants are correct.
- Gaussian Poincare follows from the displayed covariance identity; the
  gradient norm and singular-covariance treatment are valid. The direct
  maximum version follows by the same Lipschitz argument.
- The PSD variance split discards only s w^T Sigma w. It retains the
  posterior mean matrix W and its independent-component subtraction.
- The tensor contraction yields the stated current-posterior energy product
  and alpha terms with the correct factors of two and n.
- Differentiating eta_t and the cross covariance gives the exact finite-step
  identity, including the beta^2/4 coefficient. The fixed-source covariance
  is not silently varied with t. The endpoint is the fixed-source
  half-product, not an original symmetric-pressure or norm minimum.
- The integrated sign is not proved and is not claimed. The universal
  midpoint covariance supplies the advertised absolute phase-reduction
  constants, but does not evaluate the current coupled posterior.

No computation, solver, finite census, or external theorem claim was needed
for this review. The original MO convergence problem remains open.
