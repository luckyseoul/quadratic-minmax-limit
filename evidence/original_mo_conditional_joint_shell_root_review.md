# Root complete review: conditional joint-shell Gaussian package

2026-09-05. Verdict: PASS for the three scoped analytic results below.
The original all-orders convergence question remains OPEN.

## Complete-read scope and roles

The root independently read the complete conditional joint-shell proof
(451 lines) and intrinsic repair proof (338 lines), including all scope
limitations. The original 358-line shell core and every added line were
reviewed separately and the combined final file was reread. The root is
the author of the direct normalization proof (265 lines); the complete
root read of that proof is an author check, NOT an independent review.
An independent exact-worker complete review of that proof is recorded
separately, including its final complex-pairing clarification.

Approved artifacts, imported with exactly the reviewed bytes:

- `NOTE_2026-09-05_CONDITIONAL_CROSS_JOINT_SHELL_UPPER.md`:
  `64d68bb2feaa59a8049d6bcc42f3ab94c845249c3088fa618916522412d0a68a`.
- `NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`:
  `e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`.
- `NOTE_2026-09-05_INTRINSIC_CROSS_JOINT_SHELL_REPAIR.md`:
  `1dcd9b1e76b00887e406e505113c854b80f0661bb3bd69283f6486fb59fa2d53`.

## Mathematical audit

1. The symmetric cross-edge compression has exchange operator
   `S_B(X)=B X^T B`, unit diagonal, and norm bound three. The
   L-squared-normalized field covariance is genuinely PSD, by the
   two-dimensional source compression and tensor compression argument.

2. The Gaussian increment excess retains the mixed product
   `(x B y')(x' B y)`. Its exact decomposition into the covariance
   form on `(x-x') tensor (y-y')` and `(d-e)^2` is correct. The
   independent cushion supplies its separate nonnegative increment.

3. The three reference means follow from actual pair correlations.
   Spectral-midpoint polarization, expected Hamming fractions and the
   half-normal threshold inequality need neither coordinate independence
   nor independent selection. The two-block refinement and all factors
   in the concentration remainder are correct.

4. The masked refill proof uses the generic finite observable theorem
   with zero internal-edge coordinates and the fixed deterministic
   internal energy in its prior. The exact Hermite remainder restricts
   to the cross coordinates. The rank-one noise is bounded by its
   actual bilinear norm, not dismissed because its rank is one.

5. Direct intrinsic normalization uses the fixed-row compression
   `I-b_i b_i^T` to prove `mu>=n-1`, and the maximum with two handles
   n=2. Its even-Hermite identity is algebraic on unused full edges;
   no nonexistent PSD covariance extension is invoked. The four-cycle
   estimate and odd row counts give the stated operator error.
   The real/complex interpolation constant, original-norm polarization,
   conditional norm cap and uniform sign comparison give the floor in
   its correct direction for actual conditional optimizers.

6. For the intrinsic shell upper, the rank-four repair is PSD by the
   explicit completed-square identity on `u ybar^T+xbar v^T`.
   The unrepaired M need not be PSD and is never separately sampled.
   The increment calculation remains algebraic until the repair makes
   an actual Gaussian comparison possible.

7. The nonuniform marginal threshold inequality has the correct
   direction. The normalized correlation matrix is treated using W
   and a PSD matrix N; its weighted first-order error, trace estimate,
   nonlinear arcsine error and baseline error yield O(L_A sqrt(n)).
   The cross error is O(n). All ratios remain controlled when k or
   the total noise tends to zero; no tail lower bound is used.

8. Trace-zero spectral-width and block-compression bounds transfer
   those reference errors to O(n^(-1/2)) error in the mismatch r.
   Clipping and the one-half Holder bound then give O(n^(5/4)) in
   the upper, including the possible endpoint singularity. The
   marginal inflation itself is only O(n). The polynomial shell
   count and both absolute-value phases are retained.

## Scope

The result is a valid conditional optimizer floor combined with a
genuine, exchange-preserving Gaussian joint-shell upper, including an
intrinsic-normalization repair. Its sharp leading evaluation remains
unproved. No full-order optimality of the conditional source, optimality
of a noisy posterior, pointwise pressure-derivative sign, or original
convergence theorem is asserted. No numerical run or finite census is
part of these proof reviews.
