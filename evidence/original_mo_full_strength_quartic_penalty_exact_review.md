# Independent complete review: full-strength quartic penalty

2026-09-05. Reviewer: exact-proof agent.

Reviewed all 239 lines of
`/tmp/original_mo_full_strength_quartic_penalty.md`, SHA256
`be57ac246c34ef75bb1e1a04a8b5aada87b99d00ff666fface5f10b3c0965c9a`.

**PASS. No corrections required.**

The following were checked independently and analytically:

1. The actual canonical sign-covariance formula yields
   `E BB^T=E B^TB=nI-nk_1A`. The mixed trace is exactly
   `E tr(B^TAB A)=k_0[n(n-1)]^2`: both factors from A exclude the
   potentially exceptional coordinate coincidences.
2. For the four-coordinate Gaussian density interpolation, the squared
   score has expectation `(1/2)tr(Y_s^2)`, and the absolute expected
   second log-density derivative is at most `(3/2)tr(Y_s^2)`.
   Thus `integral|p_s''|<=2tr(Y_s^2)<=8||X||_F^2`; the first two Taylor
   terms vanish by independent odd sign factors. Integrating twice gives
   the stated constant 8 on the unordered correlation-square sum.
3. A distinct rectangle has four correlations of magnitude
   `|alpha|/mu` and two of magnitude `1/mu`. The bound with constant 36
   also covers failure of the small-correlation hypothesis. Repeated
   rows or columns give exactly `2n^3-n^2` contributions. These arguments
   remain valid when the full canonical Gaussian covariance is singular.
4. Direct block squaring gives the stated quartic expression with the
   negative `-k_0 S_2^2` contribution and the cubic term of either sign.
   The inequality `|(2/pi)arcsin u|<=|u|` follows from the chord bound
   for convex arcsine on `[0,1]`. It gives `k^2<=8C^2`.
5. Maximizing the quadratic in `sqrt(v)` in (8) gives exactly
   `2t(2-t)^2 k^2 n/(4-t)`. Combining this with the positive mixed
   contribution and the rectangle bound gives `(41+88C^2)tn`.
   The resulting penalty cost is `O(lambda n)`, not `o(n)` at fixed
   lambda, as the source explicitly states.
6. The exact diagonal of `H^2` is
   `(n-1)[(n-1)+2alpha^2]`. The pair Schur complement is therefore the
   stated positive rank-at-most-two subtraction, with trace at most
   `4q` and diagonal entries at most `4 max_{e!=f}|T_ef|^2`.
   Conditioning is valid at full strength because the selected
   two-coordinate principal covariance stays positive definite for
   sufficiently large n; its complement may correctly remain singular.

The note does not identify the pair-conditioned sign law with the
unconditioned law or transfer source optimality to the cross-tilted
posterior. No unresolved pressure comparison is used in these estimates.
