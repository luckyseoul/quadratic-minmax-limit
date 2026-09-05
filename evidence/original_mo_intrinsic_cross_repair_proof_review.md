# Independent complete proof review: intrinsic joint-shell extension

Reviewer: optimized_profile_proof. Date: 2026-09-05.

Reviewed the complete 338-line source, from title through final scope:

`/tmp/original_mo_mu_joint_shell_extension.md`

SHA256:
`1dcd9b1e76b00887e406e505113c854b80f0661bb3bd69283f6486fb59fa2d53`

Result: **PASS. No mathematical corrections requested.**

## Checks completed

1. Intrinsic normalization. The row principal block of H is
   `I-b_i b_i^T`, so mu>=n-1; the maximum with two handles n=2.
   `R=I+H/mu` has unit diagonal and lies between zero and 2I.
   The proof does not replace mu+1 by a whole-source spectral square.

2. Rank-four positive repair. Expanded the covariance quadratic form
   on `U=u ybar^T+xbar v^T` independently. All terms in equation (5)
   match, including both exchange contractions, both completed squares,
   and the separate `mu v n` term. Only the kR part requires the repair.
   No decomposition into a Gaussian with indefinite covariance M and
   independent repair noise is used.

3. Full mixed-replica comparison. Equation (6) has the correct factor
   k/(4mu) for increment excess divided by two. The exchange remains
   `(x B y')(x' B y)` until its exact difference-square identity is used.
   Adding P supplies nonnegative increments and a genuinely PSD field
   covariance, so Gaussian comparison is valid.

4. Heteroscedastic Hamming bound. The variance lower bound a_0 has the
   correct direction in the truncated half-normal expectation. The
   threshold argument permits flip indicators depending on every field
   coordinate. Concavity of square root gives exactly the marginal
   inflation term in (10). The normalization of the spectral mismatch
   is 2n R_lambda and its inversion is the stated delta(r).

5. Reference-energy stability. Verified the complete normalized
   correlation calculation using W and N. The first-order weighted
   error, the trace bound on AN, the Frobenius bound on the nonlinear
   arcsine remainder, and the separate baseline arcsine error all have
   the asserted orders. The off-block repair is zero, giving the exact
   attenuated cross arcsine formula and its O(n) error.

6. Uniform tail behavior. All estimates use ratios bounded by
   1/(n ell), rather than a positive lower bound on k or a_0. They
   therefore remain uniform when the threshold tends to either tail.

7. Leading mismatch transfer. Block compressions of H_lambda and its
   zero trace give the coefficient bounds used in (18). The reference
   errors consequently give O(n^(-1/2)) mismatch error. The stated
   one-half Holder bound, clipping step, and O(n^(5/4)) loss are correct.

8. Actual Gaussian maximum and scope. The number of shells is bounded
   by the stated polynomial; both absolute-value phases are included.
   The cross maximum has Lipschitz constant n and covariance norm at
   most 2k+v, giving (21). Conditional cross optimality is not promoted
   to full-order optimality. The intrinsic sign-to-Gaussian floor is a
   separately named prerequisite, not silently replaced by the
   whole-source covariance theorem.

The conclusion is a valid intrinsic-covariance upper with a rank-four
PSD repair and a leading-order-equivalent mismatch formula. The source
does not claim that its remaining sharp evaluation is proved or that
the original MO limit is closed.

No computation, signing census, or solver result was used in this
review. The reviewed source was not modified.
