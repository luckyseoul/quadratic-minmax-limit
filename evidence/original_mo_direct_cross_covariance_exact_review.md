# Independent exact review: directly normalized cross covariance

2026-09-05. Reviewer: optimized_profile_exact. Final verdict: PASS.

## Artifact and complete-read scope

The complete 265-line draft
`/tmp/original_mo_direct_cross_covariance_normalization.md` was read and
independently checked, covering every displayed formula (1)--(19).

- Initial fully reviewed SHA-256:
  `d0c1fc8ae9459f99751978bd99a88906d1fdb9cc16b23189891dd1f75abe87a8`.
- Final approved SHA-256:
  `e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`.

The initial draft's complex-pairing sentence had an ambiguous plus sign
if interpreted as ordinary bilinear pairing. The author clarified the
pairing to `u^* T v`, making the stated plus sign exact. I reread this
change, checked its final hash, and verified that replacing this single
line (182) by its previous text recovers the initial hash exactly.
No mathematical result or bound changed in this delta.

## Independent checks

1. The operator `S_B(X)=B X^T B` is self-adjoint under the real
   Frobenius inner product, with coordinate kernel `B_il B_kj`.
   Consequently `H=A tensor A-S_B+I` is self-adjoint and has zero
   diagonal. Its fixed-row compression is exactly `I-b_i b_i^T`,
   proving `||H|| >= n-1`. The specified normalization therefore gives
   diagonal-one covariance `0 <= R_mu <= 2I`, including `n=2`.

2. On two disjoint cross edges, write `a=A_ik A_jl` and
   `b=B_il B_kj`. The compressed whole-edge four-cycle entry is
   `Q=-ab`, whereas the correlation numerator is `a-b`. Thus the
   nonzero correlation magnitude `2/mu` occurs exactly at `Q=1`.
   Together with magnitude `1/mu` on adjacent edges and zero diagonal
   remainder, this proves the entire even-Hermite identity (8),
   including its orientation gauge and retained rank-one term.
   This identity is algebraic and needs no covariance extension on
   unused edges.

3. The full-edge row count in (9) is
   `2(N-2)*1 + [(N-2)(N-3)/2]*8 = 2(N-2)(2N-5)`.
   With `u<=1/mu^2`, `w<=2/mu^2`, and the stated independently
   reviewed four-cycle bound, (10)'s first estimate follows.
   Substituting `N=2n` and `mu>=n/2` bounds its odd contribution by
   `128/n`, its `L_K` contribution by `8 L_K/n`, and its remaining
   even contribution by `48/n`. The printed `180/n` is safely loose.

4. The rank-one cost is exactly bounded by
   `sqrt(w) E|xi| n^2 <= 2 n^2/(sqrt(pi) mu)=O(n)`.
   Operator covariance error `delta` is handled in both directions
   by PSD domination after adding independent `delta I` noise.
   The finite maximum has at most `2^(2n+1)` augmented states and
   coefficient norm at most `n`, irrespective of deterministic
   internal energies. This proves the error `C n sqrt(1+L_K)`.

5. The shifted-sign theorem applies to the actual cross coordinates
   with latent covariance constant two. Deterministic internal
   energies enter the prior and do not change its uniform derivative
   estimates. Its error and the preceding covariance comparison give
   (5), uniformly in the listed variables, including all real `h`.

6. The clarified complex pairing bounds the complex cube-to-one
   norm by twice its real counterpart. Halfway interpolation with
   the entry bound proves `||T||_op^2 <= 2 beta_R(T)`; the supplied
   direct three-lines argument gives the same estimate. Independent
   sign rounding extends the zero-diagonal quadratic bound to the
   real cube, and polarization gives `beta_R(T)<=4 Phi(T)`.
   Thus `L_K^2<=8 Phi(K)`. Under the stated norm cap,
   `n sqrt(1+L_K)=O(n^(11/8))=o(n^(16/11))`; no conference-scale
   spectral hypothesis has been introduced.

7. Flipping only `x` verifies (16) exactly. There are `2^(2n-1)`
   rank-one sign coefficient matrices, already including both signs.
   The independent-cross-sign maximum estimate therefore gives
   (17) with the printed constant. An actual conditional minimizer
   inherits the resulting original norm cap.

8. Every threshold output is an admissible cross signing with the
   same fixed diagonal blocks. This proves the one-sided conditional
   floor (19). The absolute order-n minimizer cap follows from the
   analogous elementary independent-sign maximum estimate. Taking
   an infimum over deterministic thresholds is legitimate because
   the error is uniform; attainment is neither assumed nor needed.

## Status limitation

This is an unconditional Gaussian reduction on actual conditional
optimizers, with a correctly directed optimizer floor. It does not
evaluate the Gaussian quantity against `2 sqrt(2) Phi(A)`, does not
make the old uncorrected joint fields PSD at the new normalization,
and does not establish original all-orders convergence.

No numerical experiment, finite-order search, or unverified source
optimality transport claim is used in this review.
