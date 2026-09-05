# Root review: canonical-gap compatibility and actual metric stability

2026-09-05. Verdict: PASS for both complete analytic proofs and their
stated conditional scope. No mathematical computation was needed.

## Complete source and review record

The full-SDP gap source was read completely, all 303 lines:
`/tmp/original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
Root also read the entire 72-line author receipt, 93-line exact-worker
review and 137-line independent docs-worker review. Their source hashes,
constants and limitations agree. Exact-worker development discussion is
disclosed; the docs reviewer did not contribute to this derivation.

The metric source received a complete 248-line initial read, followed by
both final changes: the explicit final-cell representative and the named
canonical-gap dependency with its exact hash. Final 252-line source:
`/tmp/original_mo_diagonal_majorizer_metric_stability.md`, SHA256
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.
Root read the complete author receipt and the full 99-line proof-worker
review. The proof worker authored the separate gap prerequisite but not
the metric-stability derivation; this involvement is disclosed.

Root contributed direction on source compatibility and the numerical
reference/fixed-metric interpretation. These root reads are collaborative
verification, not a claim of independent authorship or of an independent
review of root's own contributions. The complete independent reviews
are recorded separately with their actual roles.

## Canonical gap and uniform actual-state comparison

The canonical unit-row frames use the literal identity diag(K^2)=N-1.
Their SDP objective is tr|K|^3/(N-1), so the full SDP gap is nonnegative.
Whitening the attained dual slack and squaring proves the weighted
residual bound with constant 4(N-1)g. The polar factor can be completed
orthogonally on the kernel without inverting K or |K|.

Root independently checked the commutator identity with its transpose
and matrix order. The first residual bounds its right inverse square
root, the other its left inverse square root; no polar/diagonal
commutation is assumed. Completeness then identifies the exact squared
weighted commutator norm with 2(S tr(D^(-1))-N^2). The row-square
bound supplies d_i>=(N-1)^2/S, giving the stated dispersion inequality.

The square-root dispersion bound and contraction of the actual T give
Phi(K-(S/N)T)<=S sqrt(delta) on the entire actual cube. Averaging unused
coordinates and flipping one whole block give the different principal
and cross constants: 2N sqrt(delta) and N sqrt(delta). The argument
does not assume that the two error diagonal blocks are opposites.

On an original-zero-source cell, the representative is selected inside
the final refined cell. The actual weighted field and its pure-cross
counterpart are both proved PSD. Their covariance operator difference
is at most 2kN sqrt(delta), giving exactly the stated finite Gaussian
maximum cost and its safe absolute-value augmentation. The actual
weighted cross matrix remains in the pure-cross field.

## All-shell trace and metric stability

The two pointwise square-root inequalities bound both D^(1/2) and
its reciprocal rescalings without a maximum-diagonal assumption.
The constant field diagonal supplies the three Frobenius quantities
in the nuclear-norm congruence estimate. Root checked its matrix
factorization and the resulting E_delta bound without commutation.

Cyclic traces give the exact rescaled natural-D traces. Retaining
F_eta-(1-|eta|)F_eta^2, whose spectral function is nonnegative and at
most 1/[4(1-|eta|)], proves the stronger cancellation estimate.
The original and weighted radius difference is at most 2sqrt(delta).
For 0<=delta<=1, both product-error terms and the second square-root
error are included in 3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|eta|).

The reference is explicitly a numerical functional using the same actual
positive covariance and contraction. Its representative weighted energy
need not be constant throughout a bin. The proved pointwise comparison,
not an invented exact scalar-I shell, justifies using that functional
in the finite actual-cell upper. The eta window is fixed first; no
uniform endpoint limit or indefinite scalar source covariance appears.

## Publication scope and next implication

The gap theorem requires a trace-optimal full-K diagonal. The metric
theorem does not require optimality until it invokes that gap bound.
A bounded original norm controls the trace scale, not the small-gap
hypothesis. Even with small gap, the actual weighted trace supremum
has not been evaluated. The complementary gap range remains untreated.
Neither convergence, nonconvergence, an identified limit value, nor the
sharp conditional dyadic upper follows from these proofs alone.

The scoped canonical updates retain all these limitations. Frozen
proofs and receipts are imported byte-identically. No source signing,
src module, test or global proof predicate is modified. Documentation
gates, hashes and backups verify publication provenance, not mathematics.
