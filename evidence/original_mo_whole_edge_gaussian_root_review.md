# Root independent complete review: whole-edge optimizer constraint

2026-09-05. Mathematical PASS; the original convergence problem remains OPEN.

Completely read the initial 318-line proof at SHA-256
`ae95099f76a54e13cc8b7edff5636a1a99b734ede04937079120bf6cd10ff0ce`.
Then completely read the added Section 7 and final scope in the 502-line
proof at SHA-256
`83e519418578aca2d4e0db04f304ffe7c4e2d82656d717308efe2cbc76edcec0`.
The final status-header-only change was read and mechanically checked:
restoring the old third line recovers that preceding whole-file hash.
The final proof has SHA-256
`6b22fb3ab1cc878b08fe79b5b57e0e661eaaa792dfc67f850d35db9f1b68bead`.
Its canonical imported copy has that same hash.

## Independent mathematical checks

The symmetric compression has diagonal `1-1/L^2` and must be
renormalized. Its operator bound is at most three at every order n>=3;
the n=3 triangle and the degenerate n=2 fallback were checked explicitly.
The actual matrix entry is the orthonormal symmetric coordinate divided
by sqrt(2). No extra factor is introduced in the copied sign entries.

The three possible absolute edge correlations give the exact entire
even-Hermite identity, not a degree-two truncation. The full four-cycle
matrix square has nonnegative entries. Its row sum is bounded using
the positive Schur map for K^2, and compression has the required factor
two. Independently checked the odd row counts, the residual `40/sqrt(n)`
operator bound, the rank-one Gaussian's actual original-norm cost, and
the exponents 5/4 and 16/11. The threshold-sign comparison invokes the
already completely reviewed quenched theorem with covariance constant
three; matching first and second moments alone is not the argument.

The root independently derived the full symmetric extension used in
Section 7. The trace-half correction is a SINGLE scalar common to every
Boolean state, with variance at most n/2, even though it may correlate
with the off-diagonal matrix. Its expected absolute cost is at most
sqrt(n/pi). The full kernel and signed augmented replica contractions
in (19) and (22)-(23) have the correct factors. The Frobenius inequality
in (24) does not require a PSD replica matrix. It gives the actual
energy-square and signed-overlap subtractions in (25).

Checked both covariance derivative signs and every factor in (26),
the finite-state bounds proving absolute integrability, and the endpoint
limit in (27). The derivation of (28) uses original NORM optimality,
not pressure optimality. The finite-maximum entropy error has precisely
the same 16/11 exponent at beta=n^(-5/11). The integral bound is uniform
over deterministic thresholds; no noisy posterior is called a minimizer.

## Scope

This is a universal whole-edge distributional reduction and an actual
same-order optimizer constraint with subleading error. It also proves
a useful current-posterior variance upper. A negative coefficient on
that variance reverses the inequality; it must not be used to infer an
unproved derivative sign. No sharp order upper or all-orders limit
follows without a further argument. Source, tests and global predicates
are unchanged; a documentation scanner is not mathematical verification.
