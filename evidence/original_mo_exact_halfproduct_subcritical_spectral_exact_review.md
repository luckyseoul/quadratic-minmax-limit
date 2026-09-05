# Independent review: exact half-product minimizers are spectrally subcritical

2026-09-05. Reviewer: optimized_profile_exact.

Completed a full read of all 330 lines of
`/tmp/original_mo_exact_halfproduct_subcritical_spectral.md`, SHA-256
`36ae5852e60d2df2a92ca2b47702b0c431f915e3d61c8fef249a5d21a655d541`.

**PASS.** No corrections required.

Checked the complete proof chain:

- Exact global half-product optimality and random completion give the
  all-subset deletion budget, with the actual complement held fixed.
- Phasewise block superadditivity yields the hereditary Boolean cap
  with the correct original-order temperature normalization.
- Independent ternary rounding extends that cap to real cube points;
  the resulting large-eigenvector l1 lower bound and all constants are
  correct. The moderate-coordinate density estimate is valid.
- The sparse pinning probabilities are at most one half. Its exact
  field-error expectation uses independence and the complete signing's
  squared column norms. Two factor-four Markov bounds give a common
  deterministic realization; no sampling claim is substituted for proof.
- The positive-density moderate field survives removal of the sparse
  support and its small mean-square error. No coordinatewise clipping
  or global bounded-field assumption is introduced.
- The actual complement satisfies the field-response prerequisite on
  its own q-spin dimension. The prerequisite's full final proof was
  independently checked by the root and another reviewer.
- Opposite pins in the two phases produce the same external field and
  cancel the signed internal energy exactly. The resulting positive
  linear pressure increment contradicts the sublinear deletion budget.

The theorem proves `||A_N||_op=o(N^(3/4))` for every exact half-product
minimizer sequence at fixed c. It does not prove an O(sqrt(N)) bound,
the analogous assertion for symmetric-pressure or Boolean-norm minima,
or convergence of the original normalized optimum. The energy-width
versus absolute-norm distinction is correctly retained.

No computation, signing census, or simulation was run for this review.
