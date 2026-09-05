# Independent review: actual Gibbs response and exact-minimizer spectral bound

2026-09-05. Complete-read independent mathematical audit. **PASS** for both
proofs below. No numerical experiment, solver, census, repository edit, or
documentation scanner was run for this review.

## Reviewed inputs

1. `/tmp/original_mo_norm_cap_field_response.md`, all 354 lines.
   SHA256: `46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`.
2. `/tmp/original_mo_exact_halfproduct_subcritical_spectral.md`, all 330 lines.
   SHA256: `36ae5852e60d2df2a92ca2b47702b0c431f915e3d61c8fef249a5d21a655d541`.

## Field-response proof

- The explicit tensor construction has unit norm because
  `sinh(log(1+sqrt(2)))=1`; finite Gram realization and Gaussian sign rounding
  recover the stated real bilinear constant without assuming an
  infinite-dimensional standard Gaussian.
- Zero-diagonal cube interpolation gives the claimed bilinear bound.
  Both finite SDP duals are feasible and attained; their diagonal entries
  are nonnegative. The resulting positive diagonal augmentation has the
  stated trace bound and satisfies both required Loewner inequalities.
- With `M=J+D`, conjugation by `D^(-1/2)` correctly gives
  `M D^(-1) M <= 2M <= 4D`. Gaussian augmentation produces the actual Ising
  marginal. Its field-dependent posterior, not the prior Gaussian, satisfies
  `E[g^T D^(-1) g] <= n+4 tr D` uniformly in the full external field.
- Conditional independence and total variance give the lower variance
  bound. The signed and bounded-nonuniform constants follow from the stated
  weighted second-moment estimates and Jensen/Cauchy--Schwarz.
- In the unrestricted-outside-field extension, `K=2B/delta` leaves at least
  `delta n/2` moderate coordinates, whose average second moment is at most
  `2H^2+8B(1+4B)/delta^2`. Multiplication by the squared lower field amplitude
  gives (13); zero-field spin reversal and integration give (14).
  Arbitrarily large coordinates outside that set remain in the actual
  posterior, and no clipping or monotonicity assumption is used.

## Exact-minimizer spectral proof

- At the fixed original temperature `beta=c/sqrt(N)`, randomizing all edges
  incident to a deleted set gives the exact-minimizer upper budget
  `a_A-a_(A_T) <= c^2 |S|/2`. Jensen gives the lower bound. Neither argument
  assumes that the induced host is optimal at any order or temperature.
- Whole-block spin reversal gives phasewise product lower bounds. Combining
  them with the completion budget and a single extreme spin configuration
  yields `Phi(A_S) <= (c+log(2)/c) sqrt(N) |S|` for every induced set.
- Independent ternary rounding converts that hereditary bound to (5).
  The eigenvector argument then gives exactly
  `||v||_1 >= |lambda|^2/(2CN)` and `|lambda| <= sqrt(2C) N^(3/4)`.
  The moderate-coordinate count and its constants are valid.
- The sparse pinning probabilities are at most one half. The squared
  sampling-error identity uses independence and the exact squared column
  norms `N-1`. Two factor-four Markov bounds jointly leave positive
  probability, so the deterministic sparse realization exists. It has
  `o(N)` pins and preserves a positive density of moderate field coordinates
  on the complement.
- The actual complement interaction has the required field-theorem norm
  cap at the original temperature. Pinning opposite spin vectors in the two
  phases makes the external fields equal and cancels the internal pinned
  energy. The resulting positive extensive response contradicts the
  sublinear completion budget.

## Verified conclusion and scope

For every fixed `c>0`, every sequence of exact global half-product pressure
minimizers at `c/sqrt(N)` has `||A_N||_op=o(N^(3/4))`.

This proves neither `O(sqrt(N))` spectral flatness nor an original-MO
cross-order comparison. It does not establish the result for absolute
Boolean norm minimizers or arbitrary leading near-minimizers. The
half-product zero-temperature slope is half the energy width, not in
general the absolute Boolean norm. The reviewed texts preserve all these
distinctions correctly.
