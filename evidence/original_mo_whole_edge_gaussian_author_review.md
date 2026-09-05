# Author receipt: whole-edge Gaussian reduction and actual optimizer constraint

Author: optimized_profile_docs_gate.  Date: 2026-09-05.
This is an author record, not an independent review of my own proof.

Final source: 502 lines in
`/tmp/original_mo_whole_edge_source_preserving_gaussian_reduction.md`,
SHA-256
`6b22fb3ab1cc878b08fe79b5b57e0e661eaaa792dfc67f850d35db9f1b68bead`.

## Exact complete-read provenance

The initial 318-line covariance and norm proof had SHA-256
`ae95099f76a54e13cc8b7edff5636a1a99b734ede04937079120bf6cd10ff0ce`.
Both root and optimized_profile_exact independently complete-read that
source and returned mathematical PASS with no corrections.

The appended full symmetric lift and integrated original-optimizer
constraint produced the 502-line source with SHA-256
`83e519418578aca2d4e0db04f304ffe7c4e2d82656d717308efe2cbc76edcec0`.
The original covariance proof was mathematically unchanged.  Both
reviewers complete-read the entire added Section 7 and final scope,
and returned mathematical PASS with no corrections.

The final source differs from that reviewed 502-line source only in
the status header, changed from an independent-review draft to an
analytic theorem with independent complete proof reads passed.
No formula, estimate, hypothesis, or scope statement changed.

## Author checks

Checked the symmetric compression and its normalization, including
the fact that basis coefficients and actual off-diagonal entries
differ by sqrt(2).  The proposed compression has zero variance at
n=2, so that case explicitly uses an independent-edge fallback.
For every n>=3, the normalized covariance has operator norm at most
three, without an operator bound on the source.

The Hermite expansion retains every even order.  The three possible
absolute correlation magnitudes give the exact rank-one plus
four-cycle/line-graph decomposition.  The nonnegative-entry F-squared
identity, the Schur bound for V-transpose-V, and the compression to
2Q yield the source-independent four-cycle estimate.  Checked the
row constants, residual 40/sqrt(n), and its Gaussian expected-max
cost O(n^(5/4)).  Even the trivial source norm bound controls the
rank-one Gaussian term at O(n), so no source norm cap is needed.

The separately reviewed shifted-sign theorem therefore gives the
absolute O(n^(16/11)) whole-order expected ORIGINAL norm comparison,
uniformly in the deterministic threshold and arbitrary fixed
internal energy.  All matrix diagonals are zero at this stage.

For the full symmetric extension, checked the covariance factors,
the common trace-half scalar, variance at most n/2, and expected
absolute-norm error at most sqrt(n/pi).  A one-phase expected
pressure is unchanged exactly; the augmented calculation retains
the signed phase factors.  The replica matrix Gamma need not be PSD.
The Frobenius contraction bound is valid for this signed Gamma and
gives the actual energy-square and overlap subtraction in (25).

Checked the exact h derivative, including both covariance derivative
signs, all factors of two, and the fully coupled posterior.  Its
integral converges at infinity.  The integrated constraint (28)
uses ORIGINAL norm optimality plus the finite-state maximum bound,
not pressure optimality.  The chosen beta=n^(-5/11) gives the
absolute n^(16/11) error in (29).

The result is a genuine SAME-order optimizer constraint.  It does
not reverse its inequality, assert a pointwise derivative sign,
prove a doubled-order upper mapping, or close the original MO limit.
