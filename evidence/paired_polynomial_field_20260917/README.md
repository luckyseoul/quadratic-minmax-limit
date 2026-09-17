# Paired polynomial Gaussian field bound: provenance and review

The [analytic proof](../NOTE_2026-09-17_PAIRED_POLYNOMIAL_FIELD_LOWER.md)
gives `liminf alpha_n >= B = 0.3258407554555742... > 13/40`.
Convergence remains OPEN. This is a new paired-field estimate, not a new
run of the completed September 15 spectral-sign calculation.

## Completed mesh verification

On 2026-09-17, NUKA ran `python3 check_certificate.py` once after the
proof and implementation were read in full and both staged hashes matched.
Python 3.14.4, SymPy 1.14.0; exit 0; all nine checks passed. The exact
stdout is preserved in [result_nuka.json](result_nuka.json).

Selected capability: a serial exact-rational and symbolic-polynomial
calculation. The 16-thread node was idle at preflight; GPU computation or
worker fan-out would not accelerate this tiny task. Fresh staging:
`/tmp/mo-paired-field-20260917.sGSbLlqH`.

The decisive polynomial certificate is `4*b*y^2*(a-b*q)`, nonnegative
for the explicitly stated real domain. The conservative rational bound
is `361809/1111000`, exceeding `13/40` by `367/555500`. No signings,
families, or spin states were enumerated. The decimal above was evaluated
separately on NUKA from the closed-form B, not used in the proof.

SHA-256:

    proof note       1f4ddcc4eb7f1098ccc480165834fbd784b763ab1ecc7734d0361d7698383a3c
    checker          73fda9a3a7c8dd87831e9eec40cb22e6ce1ca958e7e10c3fe3250cfefb4840e9
    NUKA result      73b0b787df8def4686675cd4804563b62c33843545a32b794fe74a2731212c32
    Gaussian lemma   7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298
    regularization   8a52b7e4f171cc2089a00a6fd288e041d52605f820e49ace419ddd5fe850bec8

The last two are the exact source notes linked from proof Section 1.

## Whole-proof author review

Root derived and read the complete proof and checked the following points:

- Both frames are actual normalized `(I+sM)^2` covariances, and both
  updated energies remain the original `sM` quadratic energies.
- Taylor expansion is only on off-diagonal entries with arguments at
  most 3/4 in magnitude. The diagonal correlation-one singularity is not
  differentiated. The Frobenius remainder is `O_L(1)`, not `O_L(sqrt(n))`.
- The entire even sign covariance is retained. Its remainder after
  subtracting `cI+bM^2` is PSD by the positive Schur-power series.
- `r_i^2<=q(w_i-q^2)` uses `M_ii=0`; dropping the diagonal correction
  would lose the paired-root estimate needed for this bound.
- Both polynomial radicands are nonnegative. The even-convex root-pair
  argument controls the actual opposite cubic moments together.
- Error transfer through M costs `L^2` at fixed L and vanishes after
  averaging and square roots. Absolute field moments are justified by
  the imported uniform scalar Gaussianization, not by variance alone.
- The same actual phase energy is used in the mean update. The update
  pays with the original Boolean norm, not a source operator penalty.
- The cap is fixed while n tends to infinity, then sent to infinity
  using the existing uniform same-order objective approximation. No
  growing-cap uniformity or two-sided source-energy approximation is used.

This is author-reviewed analytic evidence plus exact scalar corroboration.
It is not independent human review or a formalization of the matrix proof,
Gaussianization, or regularization.

## Continuation boundary

Reuse these results; do not rerun the checks as new research. No cross-order
upper comparison is supplied. The remaining convergence task is to close
such a comparison (for example the integrated optimized-pressure defect),
not to tune p or the last decimal of B.

A separate all-real-domain solver check was prepared for Soulkiller but
has not run; its default Python lacks Z3. It is not part of this milestone
or counted as evidence. No research worker was left running there.
