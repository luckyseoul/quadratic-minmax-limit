# Independent complete review: shifted signs and the simple mean Gaussian

2026-09-05. Reviewer: optimized_profile_exact. **PASS for both complete proofs.**
This is a proof-read receipt, not an order-comparison theorem.

## Exact final sources read

1. `/tmp/original_mo_shifted_sign_gaussian_universality.md`, 317 lines,
   SHA-256 `a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`.
2. `/tmp/original_mo_shifted_threshold_covariance_reduction.md`, 338 lines,
   SHA-256 `74457650912a515eaf6a209b184e5c1404a13fc48c68464068871ebd61236680`.

The first proof was read completely at its preceding hash
`a2ad221c8e967db5c018e4578afff30e609f3d56b8be0ec2203114416668dd0a`.
The final change is exactly its third-line status header: restoring that
single line reproduces the preceding whole-file hash. Its formulas and
proof text are unchanged. The second proof was read completely both at
its preceding 307-line hash and again at its 310-line hash
`6c663c8d7c1081333cce780fc729ffbf0430c418253f4bc901ad8ce85be5bf76`.
That change correctly replaces the conditional prerequisite language by
the now-proved first theorem and its exact final source hash. The entire
final added tail through equation (21) was then read at the 338-line hash
above. Its unchanged earlier formulas retain the full-read PASS.

The smooth-activation comparison and Gaussian Holder prerequisites were
also reread in full in Sections 2 and 3 of the canonical
`NOTE_2026-09-05_CORRELATED_SIGN_GAUSSIAN_FREE_ENERGY.md`.
The smooth proof does not require an odd activation.

## Shifted-sign universality checks

- Conjugating the latent covariance by the threshold-orientation diagonal
  retains the correlation diagonal and the same operator bound. The
  deterministic mean is absorbed exactly into the fixed prior; the same
  deterministic normalizer occurs on both sides of the comparison.
- The mean-preserving OU smoothing is exactly
  `2 Phi((x+h sqrt(1+epsilon^2))/epsilon)-1-s_h`. Its first two
  derivative suprema are uniform in the threshold. Gaussian invariance
  makes its centered mean exactly zero.
- Differentiating equal-threshold bivariate Gaussian disagreement gives
  the displayed derivative in equation (5). It increases up to threshold
  zero and decreases thereafter. The correlation-one endpoint follows
  by continuity. This proves the threshold-uniform sign-noise bound.
- All positive Hermite degrees, including the even degrees, occur in the
  covariance. Every integer Schur power is positive semidefinite and
  bounded above by `K I`. Therefore the OU covariance decrease is PSD,
  with operator norm at most `K v_h(epsilon)`.
- The exact relation `v_h(epsilon)=2 D_{-h}((1+epsilon^2)^(-1))`
  gives the stated `2 sqrt(2) epsilon/pi` bound.
- The residual is centered and bounded by two, not one. Its variance
  is bounded by the same `v_h`. The resulting Bernstein scale is
  `2 K |gamma|/3`; Gaussian Holder and the finite-state maximum bound
  give the coefficient `2^(5/4)` in the first term of equation (1).
- The smooth-activation contraction retains the one current posterior
  and the mixed-time diagonal factors. Its bound uses signed matrix
  contractions before taking absolute values. The integrated constant
  is `4 |gamma|^3 m K^2 ||f''|| ||f'||^2`.
- Gaussian covariance restoration uses the actual positive-semidefinite
  posterior Hessian, whose trace is at most `gamma^2 m`. Singular latent
  covariances and singular independent matching noise are handled
  without changing the latent unit diagonal.
- The fixed-critical-temperature exponent is `17/18`. The full error
  formula, with auxiliary `c=n^(1/22)` and `epsilon=n^(-1/11)`, gives
  normalized expected-norm error `O(n^(-1/22))`, equivalently raw error
  `O(n^(16/11))`. This does not exchange two limits without control.
- On the actual source threshold law, the omitted cross-diagonal entries
  cost at most `n` for the signs and at most `n` in expected Gaussian
  norm. Their correlations with off-diagonal entries are immaterial to
  this Lipschitz restoration bound.

## Even-Hermite covariance and Boolean-norm checks

- Column vectorization gives the claimed different-row/different-column
  tensor correlation and the two shared-index correlations. Threshold
  reversal changes exactly the even Hermite coefficients' signs.
- Equations (10) and (11) reproduce every off-diagonal covariance entry,
  and their diagonals are zero. The residual independent variance is
  already accounted for in the linear-plus-independent base covariance.
- The retained matrix in equation (12) is PSD. Expanding it gives exactly
  equation (13); the nonvanishing even correction is not discarded in
  operator norm.
- Midpoint bounds imply the odd remainder estimate and then the complete
  residual bound `(2n+3)/mu^2+2/mu <= 18/n`, including `n=2`.
- The independent Gaussian realization in equation (15) has exactly the
  retained covariance and zero diagonal. Polarization on the cube gives
  the used bound `beta(A)<=4 Phi(A)` for the original quadratic norm.
  Gaussian maximum control of the row and column multipliers yields
  its actual expected Boolean-norm cost
  `O(Phi(A) sqrt(log(2n)/n))`; rank alone is not used as a norm bound.
- The operator-error comparison follows from Gaussian convex order and
  addition of independent isotropic noise. With the augmented absolute
  value states, its explicit cost is at most `sqrt(90 log(2)) n`.
  Arbitrary fixed internal energies are valid additive offsets.
- Removing and restoring both Gaussian models' cross diagonals costs
  at most `2n`. Thus the Gaussian comparison is uniform in the real
  threshold, source, internal interaction, and `|theta|<=1`.
- The final use of the separately proved first theorem is valid and
  makes the displayed combined shifted-sign comparison unconditional.
- The added independent-Rademacher cap uses at most `2^(n-1)` distinct
  unsigned Boolean configurations and their two energy signs. With
  `n(n-1)/2` independent edge signs, the exponential-maximum estimate
  is `sqrt(n^2(n-1) log(2))`, and a minimum is at most this expectation.
  Every realized threshold block is a complete order-`2n` signing.
  Thus the uniform comparison legitimately gives equation (21) for
  any exact order-`n` minimizer after taking the infimum over fixed real
  thresholds. The error is absolute and uniform, and no minimizer in
  the threshold variable or disorder-adaptive threshold is assumed.

## Scope gate

Together the two proofs compare the actual shifted-sign paired norm to
the simple model `s_h A + 2 phi(h) G + sqrt(1-s_h^2-4 phi(h)^2) W`,
with the stated subleading expected error on original norm-capped sources.
They also bound `m_(2n)` above by the infimum of these expected Gaussian
paired norms plus the absolute `O(n^(16/11))` error. They do not bound
that infimum by `2 sqrt(2) m_n`, sign the model's mean/noise
derivative, or prove the original all-orders MO convergence statement.
Those missing upper and order-comparison steps remain explicitly open.
