# Independent full review: invariant Boolean algebra and actual lower bound

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact source, prerequisites, and independence

The reviewer read every line of
`/tmp/original_mo_hadamard_boolean_algebra_lower.md`, 363 lines,
SHA-256 `68ce3f2f2a8fa2280208a9e145f508b6c2b2520d81e83185f579aaac89838a5d`.
The reviewer did not contribute to the modified construction, its
Boolean-algebra argument, or its clipping/dilation consequences.
No source was changed.

Both complete prerequisites were independently read in the same review:

- The 351-line sparse-flip construction at SHA-256
  `0d2355f94734b4c1e950c1e05c6df75df38b5ce181ba7fce550a4245e11328ed`.
- The final 324-line Gamma/clipping theorem at SHA-256
  `bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`.

The complete primary rectangular Bernstein theorem statement and the
existing exact-rational clipping certificate were also checked, as
recorded in the separate prerequisite review receipts. No computation
was rerun for this additive note.

## Invariant algebra and grid checks

The base permutation is a symmetric involution commuting with F.
The changed Hadamard background remains symmetric orthogonal after
normalization and fixes both the constant vector and the printed
balanced sign vector. Every tensor character is therefore fixed.
Those characters span all functions of the binary pattern, proving
the full invariant Boolean algebra, not merely invariance of a few
linear sign vectors.

The number of binary patterns, their equal multiplicities, the count
of blocks, and each block size all match. Choosing the first `2ell`
bits gives exactly `p=4^ell` blocks. Their indicators lie in the fixed
algebra, so the compressed background is I_p rather than the previous
O_p. This distinction is essential and correctly retained.

Within each block the remaining patterns give the exact stated mean
grid. Independent choices in different blocks are still global
functions of the binary pattern and hence fixed by the background.
Nested plus-sets give product mean exactly `1-|a_i-b_i|`, including
opposite endpoint means. This proves the actual paired-sign overlap
formula and grid density as the order grows.

## Modified embedding and actual Lambda lower

The deterministic model acts as C on block constants and as an
orthogonal operator on the complementary subspace. The changed
flip rule has the correct `C-I_p` mean correction. The previously
checked row/column variance estimates and concentration argument
apply unchanged, with the new fixed template constant.

All spectral and SDP conclusions remain asymptotic; finite-order
scalar optimality is not asserted. The actual sign witnesses give
the precise normalized energy (13). Continuity on the fixed compact
cube makes the grid maxima tend to Lambda_I. Uniform operator error
then proves the actual liminf lower in (15), while orthogonal
completion still gives Gamma as the upper.

For every template, choosing equal means gives (16). The nested
Bernoulli covariance is bounded by the product of its standard
deviations; averaging and Cauchy--Schwarz prove `Lambda_I<=Gamma`.
When C is symmetric PSD, its quadratic inequality and the residual
arithmetic-geometric mean bound reduce Gamma exactly to the common
mean quadratic maximum. The matching choice gives equality, so
both functionals and the actual Boolean limit in (17) coincide.
No such equality is asserted for general nonsymmetric templates.

## Positive top frame and dilation checks

A symmetric template with a unit-row frame satisfying `CW=qW`
has the required positive quadratic top frame. This does imply the
bipartite scalar-SDP value, but the reverse implication is not
silently imported. Gaussian conditional expectation cancels the
mixed residual terms, and the possibly indefinite quadratic residual
is bounded below by minus q times its squared norm. Thus (19) is
a genuine actual Boolean lower through Lambda_I, with exactly
the printed `D_f(q)` and no completion-fluctuation penalty.

The clipping moment identities and increasing D use the already
reviewed exact P/phi enclosures. They give the strict lower constant
in (22) under the additional positive-frame hypothesis. This is
reuse of the same prior rational certificate, not a new scalar run.

For a general scalar-optimal bipartite template, matched frames give
the two equations (23). Symmetric dilation has the stated norm and
positive top frame. Tensoring with I_2 fits the power-of-four size
and is permutation-equivalent to two copies. Its frame, SDP value,
Boolean value, and preserved finite-template Boolean/SDP ratio have
the correct factors of two and four.

Crucially, this dilation changes the template and the actual large
sign family. The actual lower bound for the modified symmetric class
is not transferred back to the original nonsymmetric embedding or
to every complete signing with a similar singular bulk.

## Scope and procedure

The note supplies actual lower bounds for its specified modified
families, not only failure of a Gamma upper certificate. It does
not construct an original conditional optimizer, internal source,
joint shell, or intrinsic source compatibility, and it does not
prove the original convergence theorem.

All review reasoning was analytic. Only source reading, metadata
inspection, line counting, and hashing were performed. No sampling,
numerical matrix construction, scalar certificate rerun, optimization,
solver, or other mathematical computation was run by this reviewer.

