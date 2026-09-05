# Independent full review: scalar-template completion obstruction

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact final source and independence

The reviewer directly read every line of the FINAL source
`/tmp/original_mo_scalar_template_gamma_bound.md`, 324 lines,
SHA-256 `bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`.
This was a complete read of the final version, not a claim to have
read the old 317-line version followed only by a delta. The reviewer
did not contribute to this theorem or its derivation. No source was
changed during this review.

## Matched frames and analytic inequality

Equality in both operator/Frobenius Cauchy--Schwarz bounds gives
`CV=qU` and `C^T U=qV`, with the correct positive proportionality.
It follows that the two frame operators are identical, including
at deficient rank. The trace and top-eigenvalue bounds in (2.2)
are correct. Linear rounding along a top eigenvector yields (2.3)
with the required factor `q-1`.

The Gaussian cross-term expansion uses a residual orthogonal to the
first Gaussian chaos. Conditional expectation cancels both mixed
terms even for endpoint correlations. The operator residual lower
bound then gives exactly `q(2c^2-v)` after normalization; no positivity
of a general template or preservation of its optimal Gram by clipping
is assumed.

The completion lower bound for two squared-norm variables is valid
on the whole unit square. The centered even function H has exactly
the stated variance after removal of its quadratic projection.
Its Hermite expansion starts at degree four. The covariance series
has nonnegative terms at every signed correlation, and the bound
by `R rho^2` is valid, including the endpoints.

Equality of the two frame operators gives pointwise equality of the
unclipped Gaussian squared norms. Thus the entire quadratic part of
their clipped norm difference cancels exactly. The remaining cross
covariance is nonnegative; the two separate variances are each at
most `R s`. This proves the fluctuation factor `sqrt(R s/2)` with
no independence assumption and no unproved equality of clipped norms.

Combining this bound with (2.3) has the correct inequality direction.
Solving the resulting quadratic in `sqrt(Gamma-1)` gives (1.3)
under exactly its stated `D(q)>=0` hypothesis.

## Clipping constants and the existing exact certificate

Gaussian integration by parts gives all five moment identities in
(5.1), including the coefficient `k=P-2phi` and the subtraction
`2k^2` in R. The enclosure signs used to bound D and R are correct,
with positivity of the lower factors before squaring.

The alternating-series parities, integration of the exponential
series, and the two signs in the Machin combination are correct.
The density-scale squared comparisons give the claimed interval
in the correct direction, and multiplication by positive intervals
gives the P and phi bounds. These analytic remainder arguments
were checked independently of the certificate output.

The reviewer also read completely the existing 105-line certificate
source, all 204 lines of its result, and the execution metadata:

- Certificate source:
  `/tmp/original_mo_scalar_template_gamma_rational_certificate.py`,
  SHA-256
  `d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`.
- Existing remote result:
  `/tmp/original-mo-template-rational.66FpRG/result.json`, SHA-256
  `fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
- Execution metadata:
  `/tmp/original_mo_template_gamma_rational_check.json`, SHA-256
  `fe31b7d0d5fac990a409093dfa1612ed0ef5ef0a6a496de62ccd2dbc0761a306`.

The certificate uses exact Fraction arithmetic, and its 28 checks
correspond to the finite comparisons printed in Sections 5--6.
All 28 stored results report PASS; the recorded source hash matches
the independently hashed certificate bytes. The metadata records
one prior remote run on soulkiller and no local mathematical run.
The reviewer did NOT rerun it locally or remotely.

The contradiction for q at least 5/2 uses the strict D lower bound,
R upper bound, and squared fluctuation comparison correctly. The
same strict arithmetic indeed excludes `Gamma<=283/200`, so the
added uniform margin needs no new numerical comparison or run.

## Scope and procedure

The result is a lower bound on the finite-template COMPLETION
functional Gamma. Where Gamma is only an upper certificate for
an actual large signing's Boolean norm, its failure does not imply
that actual norm exceeds the target. The source preserves this
distinction and supplies no hidden attainment or source-compatibility
claim. It does not settle original convergence.

This was an independent analytic proof review plus read-only audit
of an existing exact-rational certificate. No mathematical computation,
new certificate execution, matrix construction, search, solver, or
optimization was performed by the reviewer.

