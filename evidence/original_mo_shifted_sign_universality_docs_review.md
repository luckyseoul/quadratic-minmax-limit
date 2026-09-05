# Independent complete-read review: mean-preserving shifted-sign universality

Reviewer: optimized_profile_docs_gate.  Date: 2026-09-05.

Reviewed all 317 lines of
`/tmp/original_mo_shifted_sign_gaussian_universality.md`, SHA-256
`a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`.
Also reread the complete smooth-activation and Gaussian Holder proofs
in Sections 2--3 of the previously reviewed correlated-sign note.
The final source differs from the complete-read draft only in its
status header; that final header was reread and no mathematics changed.

Result: mathematical PASS; no corrections requested.

The orientation conjugation preserves unit diagonal and the covariance
operator bound.  Its deterministic mean is absorbed exactly into the
prior, with the same normalization on both sides.  The OU smoothing
formula contains the necessary shift `h sqrt(1+epsilon^2)` and is
exactly mean-preserving; its derivative norms are uniform in h.

Independently differentiated the Gaussian threshold disagreement:
`D_a'=2 phi(a)[1-2 Phi(a sqrt((1-q)/(1+q)))]`.  Its signs prove that
the maximum occurs at zero, including q=1 by continuity.  The
identity `v_h=2 D_-h(exp(-2t))`, the arccos/arctan conversion, and
the constant `2 sqrt(2)/pi` in the uniform variance bound are correct.

All Hermite orders are included.  Positive Schur powers and their
operator bound prove the full PSD covariance decrease, both for the
initial smoothing and for every further OU smoothing.  No oddness
assumption has been imported into the smooth-activation lemma.

The residual is centered and bounded by two, not one.  The Bernstein
scale is therefore `2K|gamma|/3`; the Gaussian Holder argument and
finite maximum calculation give exactly the first constant `2^(5/4)`.
The existing smooth proof contracts the signed mixed-time covariance
matrices against one actual posterior before taking absolute values;
it applies here and gives the stated cubic term.  The final Gaussian
covariance interpolation keeps the actual Hessian and has the stated
positive sign and constant.  Singular covariance causes no defect.

The full temperature-dependent estimate yields the claimed expected
maximum exponent with `c=n^(1/22)`, `epsilon=n^(-1/11)`, uniformly
in the threshold and arbitrary fixed internal energy.  Restricting to
the off-diagonal cross coordinates and restoring the n diagonal
entries costs at most 2n in expectation, irrespective of correlations.

Thus the separately stated shifted-sign prerequisite in the threshold
covariance-reduction note is now supplied by a fully reviewed proof.
The simple mean Gaussian model's sharp actual upper comparison remains
open; neither theorem asserts an original-norm order comparison or
closure of the original MO convergence problem.
