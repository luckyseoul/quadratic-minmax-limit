# Independent full review: Hadamard sparse-flip finite template

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact source and independence

The reviewer read every line of
`/tmp/original_mo_hadamard_sparse_flip_template.md`, 351 lines,
SHA-256 `0d2355f94734b4c1e950c1e05c6df75df38b5ce181ba7fce550a4245e11328ed`.
The reviewer did not contribute to this construction or its proof.
The source was left unchanged.

The entire statement of the published rectangular matrix Bernstein
Theorem 1.6 was checked directly in
[Tropp's author-hosted paper](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).
The dimension factor, centeredness, almost-sure operator bound, two
variance operators, and tail denominator agree with the imported
theorem. This is verification of the published theorem statement and
its application, not a claim to have independently replayed the full
matrix-concentration paper.

## Construction and concentration checks

The base matrix is symmetric, orthogonal after normalization, and
regular. Its tensor powers have the stated sign entries and row sums.
The block-constant subspace is reducing for the normalized background.
Consequently the deterministic model is exactly the direct sum of C
and an orthogonal complementary restriction, with the stated complete
singular multiset and operator norm.

The entry-flip probabilities are valid once `n>=gamma^2`. Entrywise
expectation gives exactly (7): the factor p in M cancels n=pm to
produce the block-isometry correction. The two different uses of J
are explicitly distinguished. Centered entries have bound two and
variance at most `4 gamma/sqrt(n)`. Each variance operator is diagonal,
and its norm is bounded by the corresponding maximum row or column
variance sum `4 gamma sqrt(n)`.

The Bernstein substitution uses threshold `n^(3/8)` and produces
the exact exponent and dimension factor in (10). The subsequent
exponent weakening has the correct direction. Its failure probability
tends to zero for a fixed template. The expected flip count and
Markov bound are correct. The union bound gives intersection
probability at least one half for all sufficiently large allowed n.

Selecting realizations on this event gives actual complete sign
matrices, with at most the printed number of changed entries and
operator error `epsilon_n=n^(-1/8)+gamma n^(-1/2)` after normalization.
This is a probabilistic existence proof of a deterministic sequence;
no efficient search or sampled realization is asserted.

## Spectrum and SDP checks

Uniform singular-value perturbation gives the flat empirical bulk
and the fixed template's finitely many non-bulk limits, including
values below one. Under the stated scalar-template hypothesis the
top normalized norm tends to q.

The vector SDP operator upper and its Lipschitz estimate both use
the correct factor n. Repeating each optimal template vector in
its block gives the exact matching primal value, since the sum
within each block of the deterministic model is `m C_ij`.
This proves (14) with the printed error.

The scalar diagonal formed from the actual operator norm is feasible.
Equations (13)--(14) give its objective gap at most
`2 epsilon_n n^(3/2)`. Thus scalar optimality is asymptotic, not an
unproved exact statement about each realized matrix or its Gram.

## Boolean-norm direction and scope

The completion maximum exists. Its equivalent absolute bilinear
form and its lower and upper bounds in (17) are correct. Decomposing
actual signs into block means and perpendicular components gives
Gamma as an UPPER via orthogonal Cauchy--Schwarz, with the correct
normalization. The operator perturbation costs exactly epsilon_n.

The block-constant witness attains `beta(C)/p`. The separate baseline
one is genuinely witnessed: the printed balanced tensor sign vector
lies in the block-mean-zero complement and is a negative eigenvector
of the background. Taking its negative for the other sign vector
gives the exact normalized value one, unaffected by the template.
Both lower witnesses and the liminf/limsup directions in (18)--(19)
are therefore justified.

No equality of the actual Boolean limit with Gamma is claimed.
A lower bound on Gamma alone cannot exclude an actual Boolean cap.
The construction does not supply an original internal source,
conditional optimizer, or source/joint-shell compatibility, and it
does not prove the original convergence theorem.

## Procedure

This review was analytic, with primary-source reading and file
hash/line-count checks. No mathematical script, matrix sampling,
solver, census, numerical optimization, or other mathematical run
was performed by the reviewer.

