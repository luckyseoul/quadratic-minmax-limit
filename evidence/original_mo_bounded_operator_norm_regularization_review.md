# Independent review: same-order regularization for the original norm minimum

2026-09-05. **PASS**, after a complete independent read of all 290 lines of
`/tmp/original_mo_bounded_operator_norm_regularization.md`.

Input SHA256:
`8a52b7e4f171cc2089a00a6fd288e041d52605f820e49ace419ddd5fe850bec8`.

The prerequisite tensor-rounding and diagonal-SDP majorizer proof was
previously independently read in full and checked. This review used no
numerical experiment, solver, census, or repository edit.

## Verified points

- The two attained diagonal duals have nonnegative entries and total
  trace at most `(4 pi/kappa) Phi(A)`. Their sum majorizes both A and -A;
  no identity shift is required for this positive-semidefinite application.
- Thresholding at `K sqrt(N)` gives the literal cardinality bound in (1).
  Principal restriction then gives `||A_T||_op<=K sqrt(N)`, including
  empty or singleton T. An empty exceptional set is handled separately.
- Each retained filler quadratic form has squared Rademacher coefficient
  sum at most two. The quarter-net estimate gives the stated operator
  failure bound. The Boolean union bound separately gives the stated
  `2 sqrt(N e)` threshold.
- The two failure bounds sum to less than one already at N=2 and decrease
  thereafter. Therefore the SAME completion satisfies both estimates;
  no independence between success events is assumed.
- The completed matrix is a complete signing at exactly the original
  order. The retained principal Boolean norm is at most Phi(A) by
  independent spin averaging. The triangle inequality with the NEW filler
  proves the objective bound. No estimate on the discarded OLD incident
  energy, or on `Phi(A'-A)`, is needed or asserted.
- The epsilon choice gives exactly the stated objective loss and
  exceptional-set size. The all-edge filler establishes nonemptiness of
  every constrained class with L>=8 and the uniform auxiliary C=2 bound.
- The same-order optimum sandwich holds for every N>=2 and K>0. Its
  constants are uniform, so an arbitrarily slowly diverging prescribed
  K_N is legitimate.
- The oscillation error in (14) is one sandwich width, not two: both the
  limsup and liminf increase by quantities in the same interval [0,error].
  This proves the equivalence in (15), while convergence at any fixed K
  remains a genuinely unproved sufficient next step.

## Conclusion

The original absolute Boolean minimum, not merely a pressure surrogate,
admits same-order leading near-minimizers with normalized operator norm
growing arbitrarily slowly. For every fixed operator threshold the
approximation error is uniform in order and decreases as `O(K^(-1/2))`.
The result is an unconditional reduction, not a proof of the original
limit or a uniform energy-function approximation.
