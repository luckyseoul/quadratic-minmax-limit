# Independent complete review: correlated signs to matched Gaussian

2026-09-05. Reviewer: exact-proof agent.

Reviewed all 395 lines of the final frozen source
`/tmp/original_mo_correlated_sign_gaussian_free_energy.md`, SHA256
`2e6537d0b1e2c4d8a72cc920e3fee50600d82be32417ba77c733aaedabc141c7`.

The complete final reread includes the explicit clarification that only
the independent additive Gaussian covariance is diagonally regularized;
the original unit-diagonal correlation matrix remains unchanged.

**PASS. No mathematical corrections required.**

This review is analytical; no numerical job was run. In particular:

1. Differentiating `C_t=Cov(f_t(G))` by the scalar OU generator and
   Gaussian integration by parts yields exactly
   `C_t'[i,j]=-2Sigma[i,j] E[f_t'(G_i)f_t'(G_j)]`.
   Including the independent Gaussian covariance derivative therefore
   gives the covariance identity (5) with its stated minus sign.
2. The Gaussian covariance identity is applied conditionally on the
   SAME additive Gaussian `Z_t`. Thus all three derivatives in (7)
   belong to the ONE posterior at `f_t(G_s)+Z_t`. Its third cumulant is
   exactly `<u_i u_j u_l>`, where `u=v-<v>`.
3. The two differentiated factors in (5) agree after exchanging i,j.
   Their sum supplies the factor 2 in (7). The resulting contraction is
   exactly `<sum_i a_i u_i (Sigma D u)_i(Sigma D_tilde u)_i>`.
   Using `|u_i|<=2`, `||u||^2<=4m`, and `||Sigma||<=K` proves the
   constant 8 in (10), without an entrywise absolute-correlation sum.
   Integrating the resulting `16 exp(-4t)` bound gives constant 4 in
   the smooth comparison (11).
4. The Gaussian Holder estimate has a complete heat-flow proof in
   (13). Its derivative is nonnegative by `Sigma<=K I`; it remains
   valid for singular Sigma, either directly in a Gaussian factor
   representation or by correlation-preserving regularization.
5. The residual has absolute mean `(2/pi)arctan(epsilon)` and variance
   at most that value. The scalar Bernstein series and Holder yield
   (16), and the specified exponential-maximum substitution gives
   exactly `sqrt(2aL)+bL`. Therefore (18) has the displayed constants
   and avoids a dimension-m coordinatewise absolute-error loss.
6. OU-smoothed sign has covariance
   `(2/pi)arcsin_entry(q^2 Sigma)`. Positive Schur-power coefficients
   prove both covariance decreases are positive semidefinite. The
   unital positive Schur map contracts the operator norm, giving
   `||C-C_epsilon||<=K[1-(2/pi)arcsin(1/(1+epsilon^2))]`.
   The half-angle identity and the stated linear bound on epsilon
   are correct.
7. The final Gaussian covariance interpolation uses the actual
   posterior Hessian, whose trace is at most `gamma^2 m`. It gives
   the stated positive direction and constant in (21).
8. The exact first and second derivative norms of the smoothing yield
   the third term of (2). With `m=n^2`, `gamma=O(n^-1/2)`,
   `L=O(n)`, and `epsilon=n^-1/9`, the first and third terms are
   `O(n^17/18)`, the second is `O(sqrt(n))`, and the last is
   `O(n^8/9)`. All are subextensive.
9. Singular endpoint Gaussian covariances cause no difficulty: Gaussian
   integration by parts can be performed on an independent standard
   Gaussian factor representation, and positive diagonal regularization
   of the additive Gaussian covariance can be removed using bounded
   derivatives and the Lipschitz log-partition bound. OU smoothing
   tends to zero in L2 because the activation is centered and bounded.
10. The bipartite application correctly absorbs the fixed internal
    interaction into a prior. If several spin configurations yield the
    same vector `sigma(y tensor x)`, their prior masses may simply be
    aggregated. The support bound and uniformity in this prior remain
    valid.

The proved conclusion is a genuine correlated-sign-to-matched-Gaussian
quenched free-energy equivalence under a bounded covariance operator.
It does not estimate the matched Gaussian pressure against the
smaller-order optimum or prove the remaining order comparison.
