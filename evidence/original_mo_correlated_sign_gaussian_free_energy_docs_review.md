# Independent complete-read review: correlated signs versus matched Gaussian

2026-09-05. Verdict: **PASS**, with no mathematical corrections requested.

Reviewed every line of the 392-line proof
`/tmp/original_mo_correlated_sign_gaussian_free_energy.md`, SHA256
`b19b54340b94085b17debf59b9f76ea97a940e73beae07a5f04791f75eec15ec`.
The reviewer independently derived the central OU contraction and the
Gaussian Holder smoothing mechanism before reading the completed proof.

## Checked implications

1. The OU covariance derivative has the correct factor and sign:
   `C_t'=-2 Sigma circ E[f_t'(G)f_t'(G)^T]`. The independent Gaussian
   covariance increase exactly cancels the product of expectations,
   leaving the negative covariance in equation (5).

2. Conditioning on the SAME independent Gaussian Z_t when applying
   Gaussian covariance interpolation is valid. Its first covariance
   factor depends only on G. The resulting posterior in (7) is evaluated
   at `f_t(G_s)+Z_t`, while the first two activation factors remain at G.

3. The third derivative is one central third moment under that actual
   posterior. With `u=v-<v>`, the signed contraction in (9) factors
   exactly. The inequalities `|u_i|<=2` and `||u||_2<=2sqrt(m)` give
   the factor 8 in (10). Both differentiated activation factors must be
   counted, producing 16 before integrating `e^(-4t)`; the final
   smooth-comparison constant is therefore 4, as stated.

4. For the activation used here, positivity of the OU covariance
   decrease follows from the positive arcsine Schur-power expansion.
   It also follows directly from conditional covariance for the Gaussian
   OU coupling. Singular covariance causes no distributional difficulty;
   Gaussian integration identities hold on its Gaussian linear image,
   or follow by positive covariance approximation.

5. The Gaussian Holder proof retains the full signed covariance in
   `a^T Sigma a`, and its bound by `K||a||_2^2` gives the claimed
   monotonicity. Bounded univariate approximation also covers singular
   joint support.

6. The sign residual is odd, centered, and bounded by one, with second
   moment at most `2 epsilon/pi`. Bernstein's series bound and the
   exponential finite-maximum optimization give (17), including its
   linear entropy term. The prior weights do not enter this bound.

7. The smoothed covariance is `(2/pi)arcsin^circ(q^2 Sigma)`, including
   its diagonal. The covariance difference is positive semidefinite;
   its operator bound and the posterior Hessian trace produce the
   stated quenched covariance error, without Gaussian annealing.

8. With `m=n^2`, `gamma=O(n^(-1/2))`, and `log|V|=O(n)`, the choice
   `epsilon=n^(-1/9)` makes both the smooth OU error and sign-smoothing
   error `O(n^(17/18))`; the covariance-smoothing error is smaller.
   Deterministic internal energies are absorbed into a fixed finite
   prior, so the application is uniform in such energies and in the
   source chosen before drawing the cross disorder.

The theorem establishes the stated quenched sign-to-matched-Gaussian
comparison. It does not evaluate the remaining Gaussian model or prove
the optimized cross-order comparison. Those scope limits are retained.

## Final frozen version

The final 395-line theorem has SHA256
`2e6537d0b1e2c4d8a72cc920e3fee50600d82be32417ba77c733aaedabc141c7`.
The reviewer read the changed opening section: the header now records
completed review, and the singular-covariance paragraph correctly keeps
the latent unit-diagonal Sigma unchanged while regularizing only the
independent Gaussian covariance by delta I. This strengthens the
exposition without changing the proof. All remaining proof text is the
previously completely reviewed version. Verdict remains PASS.
