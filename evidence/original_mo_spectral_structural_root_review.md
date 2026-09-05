# Root complete-proof review: spectral and actual-Gibbs structure

2026-09-05. **PASS** after complete mathematical reads of the following
proofs. No mathematical computation, solver, or census was run. The
original convergence problem remains open.

## Positive structural results

- `NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md`, SHA256
  `8a52b7e4f171cc2089a00a6fd288e041d52605f820e49ace419ddd5fe850bec8`:
  root completed the full 290-line proof read. The unshifted diagonal
  majorizer, same-order recompletion, simultaneous net/Boolean estimates,
  uniform original-optimum sandwich and exact oscillation reduction all
  check. The probability sum is below one already at order two. The
  one-sided sandwich gives the stated single-error oscillation bound;
  growing the operator cutoff with N does not itself prove convergence.
  Both independent complete reads also passed. This theorem directly
  concerns the original norm objective, unlike the pressure-only results.
- `NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md`, SHA256
  `46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`:
  the elementary tensor rounding constant, finite SDP dual attainment,
  diagonal majorizer, actual field-dependent Gaussian posterior moment,
  and all signed/nonuniform/moderate-field variance and pressure bounds
  check. Arbitrarily large omitted field coordinates stay in the posterior;
  no clipping or upper susceptibility bound is used.
- `NOTE_2026-09-05_EXACT_HALFPRODUCT_SUBCRITICAL_SPECTRAL.md`, canonical
  SHA256 `10dfe02b63aa3c4aa987ce48d4a3e660e90509b43e6a50a1a002ba9ecc1cc522`:
  the full exact completion budget, hereditary cap, eigenvector
  delocalization, sparse pinning probabilities, common Markov realization,
  and opposite-phase internal-energy cancellation check. The canonical
  file differs from independently reviewed source SHA256
  `36ae5852e60d2df2a92ca2b47702b0c431f915e3d61c8fef249a5d21a655d541`
  only by replacing the temporary prerequisite path with its canonical
  repository path. Root inspected that exact one-line change.
- `NOTE_2026-09-05_HALFPRODUCT_NEARMINIMIZER_STRUCTURE.md`, SHA256
  `dccc256d3b7119c666102e54cffe3a2026d31edc1bcd0c4366a15ce92c762f0f`:
  every unnormalized near-minimizer error is retained. The l1 truncation
  avoids amplifying that error, and the reviewed restricted interpolation
  bound supplies the delocalization contradiction. Sparse pinning still
  contradicts the sublinear budget. The direct high-field estimate gives
  both a sublinear exceptional set and sublinear field l1 mass; the
  moderate-field contradiction and final two-threshold limit are uniform
  over all permitted sparse sets and pins. Both the subcritical operator
  bound and sparse Boolean-energy robustness hold for all leading
  half-product near-minima, not merely exact minima.
- `NOTE_2026-09-05_NUCLEAR_SPECTRAL_BUDGET.md`, SHA256
  `ee8ad5ff3dbf9aa9e251c4190e98ee1671c9a2140c759ba6f768f8c9c03ef13d`:
  PSD absolute-matrix shifts, positive diagonals, Gaussian angular
  rounding, pairwise harmonic-mean estimate, and both Gibbs variational
  inequalities check. The effective-rank conclusions are lower bounds
  on bulk size; they do not bound the largest eigenvalue or transport orders.

## Full-strength canonical-law consequences

- `NOTE_2026-09-05_FULL_STRENGTH_SPECTRAL_DEFICIT.md`, SHA256
  `c37123564ec9bba8c8f16048a3ce0d1a40348990cd82d87070f223ce7aa51ad6`:
  the missing-variance repair, angular coupling, conditional biased-sign
  replacement, Gaussian variance penalty, and necessary deficit constant
  check. Singular covariances are allowed. Probability transfer is only
  vanishing probability, not an exponential bound for the original law.
- `NOTE_2026-09-05_FULL_STRENGTH_CONSTRUCTION_CAP_HOSTS.md`, SHA256
  `add44c2a7a33d68cf18d06abf699f82c02c4850b35154681f6dd26b7b92de83b`:
  the reused complete twin modules and Walsh base meet the asserted norm
  construction cap. Exact invariant eigenvalues, bounded filler spectrum,
  actual Gibbs centering and mixed-sector multiplicity give precisely
  the stated sublinear deficit. No exact- or leading-near-minimum premise
  is supplied by this construction-cap example alone.
- `NOTE_2026-09-05_FULL_STRENGTH_HALFPRODUCT_NEARMINIMIZERS.md`, canonical
  SHA256 `ad83095163cf8e969e542a6626382dabaa5adb5e2ffce4bfffea274c813b53e4`:
  the new subcritical minimizing base permits one dominating but
  subcritical module. Its uniform energy cost, comparison at the final
  temperature and same-temperature optimum monotonicity prove actual
  leading half-product near-minimality. The actual tensor has exactly two
  fixed negative directions, so its deficit is `2(1-r)` eventually. The
  mean and probability consequences are correctly scoped. The canonical
  file differs from independently reviewed source SHA256
  `349e20c8d57df46256989f13725965c14641cff60102a2299fb9882a0a374dfc`
  only by replacing the temporary spectral-deficit prerequisite path with
  its canonical path. Root inspected that exact reference-only change.

## Scope audit

Root also completed the entire covariance-floor corollary, canonical
`NOTE_2026-09-05_ACTUAL_GIBBS_COVARIANCE_FLOOR.md`, SHA256
`98727eabf74902a45ef429a94c65729a783175a9474d9f21dad9e0523a673d51`.
The actual total-covariance diagonal component, projection-tail estimate,
two successive PSD trace comparisons, and zero-diagonal correction all
check. Its source hash was
`0b283e677e16b42dae158b085040e9b651f89d49e8a608a0f4cf1fbb24d67b7d`;
only the prerequisite's temporary path changed to its canonical path.
The quadratic lower floor does not refute the quadratic radial target.

Half-product pressure approaches half the energy width, not in general
the absolute Boolean norm. None of these arguments silently identifies
half-product, symmetric-pressure, and norm minimizers. A fixed-positive-
fraction order comparison is still absent. The source near-minimizer
construction is compatible with both new structural conclusions. The
original global status must remain OPEN.

Independent complete reviews are preserved alongside this record in
`original_mo_field_response_subcritical_spectral_review.md`,
`original_mo_exact_halfproduct_subcritical_spectral_exact_review.md`,
`original_mo_halfproduct_nearminimizer_structural_review.md`, and
`original_mo_full_strength_spectral_deficit_exact_review.md`.
