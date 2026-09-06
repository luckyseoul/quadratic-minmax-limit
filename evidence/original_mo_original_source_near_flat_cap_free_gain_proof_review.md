# Independent full-source review: cap-free actual near-flat gain

2026-09-06. Role: independent reviewer of the NEW cap-free extension.
I authored the earlier 612-line near-flat and 553-line adaptive-gain
prerequisites. Those earlier arguments are imports, not independently
reauthored here. I did not supply the new tail, projector, trace-split,
or smoothed-output arguments in the present extension before this review.

## Frozen source and complete-read scope

    /tmp/original_mo_original_source_near_flat_cap_free_gain.md
    325 lines
    SHA256 0dfa5f62baaa57850a661bbc98d33d32440c783cccb11eaf5446feffbd81f7d4

I directly read all 325 lines in one untruncated response, then checked
every consequential new link below. I also read the complete 61-line
author receipt, SHA256
6c37407a6058031a66e07115e859a1eddc876f17b6304468b1ff076eee584f4a.
The frozen source was not edited. Verdict: PASS, with no correction needed.

## Independently checked extension links

1. The actual identity tr(M^2)=n-1 and the limit law's second moment
   exactly equal to one imply tr(M_o^2)=o(n) at the fixed cutoff C=5/3.
   The bounded-cutoff integral is legitimate because neither cutoff is
   an atom. Cauchy--Schwarz gives ||M_o||_1=o(n). No weak-law operator
   convergence or fourth-moment control is used.

2. Full threshold projectors remain full, including outlying eigenvectors.
   The squared first spectral error and absolute second spectral error
   are controlled by the second-moment tails. This proves the stated
   Frobenius and nuclear estimates, respectively. Diagonal control then
   gives P_ii->rho and T_ii->0 in empirical mean square. The good/bad
   diagonal normalization changes a bounded frame by o(sqrt(n)) in
   Frobenius norm and o(n) in nuclear norm, without deleting source
   coordinates. The repaired correlation has bounded operator norm.

3. Schur multiplication preserves the uniform power bound. Its Hermite
   series converges in operator norm. Literal complete entries give the
   uniform O(sqrt(n)) odd-power trace error. First-moment tail control and
   Frobenius pairing recover the actual positive baseline 5*kappa/8.

4. P+T=rho*R0 gives ||M^2/lambda0^2-rho*R+M/lambda0||_1=o(n).
   Pairing this nuclear error against uniformly bounded Schur powers is
   legitimate. Their even scalar-power sum is at least n. This yields
   the full higher-chaos mean lower 1-kappa, uniformly before summing
   the Hermite weights. Actual complete row lengths give every stated
   local variance upper, despite unbounded ||M||op.

5. First-chaos alignment uses a bounded spectral part paired with
   ||R-R0||_1=o(n), and a positive outlier part of trace o(n) paired with
   ||R-R0||op=O(1). No product of uncontrolled operator norms appears.
   The full positive-projector pairing is also controlled by the same
   second-moment tails. Thus the empirical squared alignment tends to zero.

6. The unchanged local Gaussianization applies to ACTUAL rows of M:
   max_j|M_ij|<=n^(-1/2) and sum_j M_ij^2<=1 exactly. It needs the
   correlation-frame cap, not a source cap. The new alignment and row
   variance bounds supply precisely the old empirical scalar hypotheses.
   I refreshed all of Sections 6, 7, and 9 of the frozen 612-line source;
   their use here does not introduce a growing-dimensional field CLT.

7. Y=clip((I+M_b/C)X) has coordinate displacement -2*X_i*p_i with
   p_i=min(r_(b,i)/(2C),1). The deterministic quadratic identity proves
   its gain even when M_b has a nonzero diagonal. The 2-Lipschitz gain
   function and E||M_o X||^2=o(n) transfer that scalar gain to actual
   full local fields; the input outlier energy is o(n) by nuclear pairing.

8. At each FIXED epsilon, the full-coordinate smoothed output is odd
   and 2*sqrt(B)/epsilon-Lipschitz. Applying Gaussian Poincare to every
   linear functional gives its stated covariance bound, including for
   singular R. This pays the entire original output outlier energy by
   ||M_o||_1=o(n). No covariance claim about unsmoothed Y or sign(MX) is
   smuggled into the argument.

9. The standard-normal small-ball bound controls the input smoothing
   error; the clipped bounded-bulk map is 2-Lipschitz. The resulting
   energy cost is exactly bounded by 2*C*kappa^(1/4)*n*sqrt(epsilon).
   Independent rounding is invoked ONLY for full zero-diagonal M.
   Taking n first and epsilon second closes the gain transfer without
   changing the actual objective or assuming a tail convergence rate.

10. The zero-field indicator in the fixed-p lower is correctly handled.
    With p=1/10 and auxiliary C=5/3 the old scalar bounds yield the same
    16/3125 gain. The rational comparison to 2/5+3/1100 is correct.

## Imports, classification, and remaining scope

The two previously fully read prerequisites were rehashed here:

    612-line source: 7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298
    553-line source: 0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51

I also refreshed the full adaptive-update Section 7 of the latter.
My earlier authorship is an overlap in imported prerequisites, not a
contribution to the new cap-free extension under review.

PASS: the full actual near-flat empirical law alone implies
liminf Phi(A)/n^(3/2)>=5*kappa/8+16/3125>2/5+3/1100, with no actual
source operator cap. C=5/3 is only the auxiliary update cutoff.
Other laws, the separate all-law cap removal, and global MO remain OPEN.
Only this separate /tmp review receipt was written. No mathematical
job, canonical edit, frozen-source edit, commit, or backup was performed.
