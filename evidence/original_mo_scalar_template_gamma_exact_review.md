# Full proof review: scalar-optimal finite-template Gamma bound

2026-09-05. Reviewer: `optimized_profile_exact`.

Conclusion: PASS. No correction requested.

Final source: `/tmp/original_mo_scalar_template_gamma_bound.md`

Final length: 324 lines.

Final SHA256:
`bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`

The reviewer read the entire 317-line predecessor, checked the reported
delta, and then independently reread ALL 324 lines of the final source.
The reviewer did not edit this proof. Earlier collaborative discussion
included the elementary clipped cross-term estimate and checks of
intermediate constants; this receipt records the subsequent complete
verification of the author's finished proof, including the matched-frame
Hermite cancellation and explicit nonlinear Gamma bound.

Verified mathematical points:

- Equality in the scalar-dual SDP operator bound gives both `CV=qU`
  and `C^T U=qV`, hence the exact common frame matrix `M=U^T U=V^T V`.
- The top-eigenvector linear projection is cube-admissible and proves
  `Gamma>=1+(q-1)lambda_max(M)/p` with the printed normalization.
- Gaussian conditional expectation eliminates the linear/residual
  cross terms for arbitrary bounded measurable odd f, including
  endpoint Gaussian correlations. The remaining operator estimate
  gives `E[a^T C b]/p>=q(2c^2-v)` with no hidden sign assumption.
- The quadratic component of `f(X)^2-f(Y)^2` cancels pointwise.
  The even Hermite remainder has covariance between zero and `R rho^2`.
  Thus its two frame averages have nonnegative cross covariance and
  their difference variance is at most `2R tr(M^2)/p^2`.
- The completion floor, its factor one half, the square-root penalty,
  and substitution of the linear bound give exactly equation (1.2).
- Solving the scalar inequality in `sqrt(Gamma-1)` gives (1.3) when
  `D(q)>=0`, including its zero and zero-variance endpoint cases.
- For unit clipping, all integration-by-parts moments, the variance
  remainder, monotonicity of D, and the rational endpoint arithmetic
  are correctly normalized. Both `Gamma>sqrt(2)` and the stronger
  `Gamma>283/200` for `q>=5/2` follow with the printed strict margins.
- The finite-sum alternating bounds and Machin identity have the right
  directions; every multiplication and square uses a positive quantity.
- The scope is exactly a bound on Gamma. The proof does not treat an
  upper certificate as a lower bound on an actual constructed matrix.

The reviewer also completely read the exact rational checker source
`/tmp/original_mo_scalar_template_gamma_rational_certificate.py`, SHA256
`d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`.
It implements the displayed finite sums and all endpoint comparisons
using integer arithmetic and `fractions.Fraction`, with no float input,
numerical integration, search, or optimization.

The existing result was read in full, without rerunning it:
`/tmp/original-mo-template-rational.66FpRG/result.json`, SHA256
`fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
It contains 28 successful checks, overall PASS, and the matching checker
source hash. The reviewer verified the result entries against the proof
and checker. The parent performed the single remote execution; this
review does not claim a second execution.

No mathematical computation or test was run by this reviewer. Local
tools only read proof/checker/result text, computed byte hashes, and
wrote this receipt. The original MO convergence question remains open.
