# Independent proof review: restricted scalar moment-feedback diagnostic

2026-09-05. Reviewer: optimized_profile_proof.
Verdict: **PASS**, with the source-provenance clarification incorporated.

The reviewer read every line of the initial 286-line source,
`/tmp/original_mo_scalar_moment_feedback_diagnostic.md`, SHA-256
`c807036a2b8008ebfeb3678235af0035faa8636bba55a7c5d4fbd13f5ab7ac13`.
The reviewer then checked the complete explicit repair/provenance delta
in the final 321-line source, SHA-256
`cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`.
The evaluated branch proofs and unattainability warnings are unchanged.

The root, docs, and exact workers supplied the branch derivations.
The reviewer independently checked them; matching preliminary algebra
was an audit rather than a contribution to those derivations. A separate
optional top-mass lemma proposed during discussion is NOT in this source
and is NOT certified by this receipt.

## Functional, repair, and mixed constraints

The normalized scalar variables, source-feedback coefficient, and two
resolvent traces have the correct factors. An actual scalar optimum
requires `d=||B||`, not merely scalar feasibility. Its optimal Gram's
top-eigenspace support is correctly treated as an additional condition.

The final source explicitly uses the actual intrinsic normalization.
The triangle inequality gives `L_A^2<=mu+q^2 n+1`, hence
`tr P/a_0<=4+2q^2` from the cited canonical repair. Dividing the two
raw trace increments by `2n a_0/d` gives exactly the constants and
fixed-metric-gap powers in (4b). No original-norm or conference-scale
cap on A is required for this step. The actual repaired field is used
before taking limits; an indefinite uncorrected matrix is not called
a Gaussian covariance.

The full mixed-rounding family (5) has the correct scalar sign-defect
coefficient. Retaining the exact squared residual gives (6), because
its bracket is exactly `q^2 E[y(1-sqrt(y))^2]`. No unattainable
optimal-Gram condition is inferred from these moments alone.

## Strong-feedback branch

At the stated endpoint, `chi=2/c_0>2`. A rank-at-most-four PSD repair
allows at most four negative eigenvalues of the uncorrected field.
Therefore the limiting singular law has the claimed support restriction.
The first trace has the printed limit as eta approaches one, and its
integrand is concave. Jensen gives (9) with no reversed inequality.

The elementary inequalities `3/5<kappa<2/3` and `5/6<c_0<1` yield
the exact strict bound `9/20<1/2`. The source first sends n to infinity
at a fixed metric gap using (4b), and only then sends eta to one.
The infimum/limsup inequality has the correct direction. No uniform
finite-n endpoint estimate is assumed.

## Weak-feedback branch

The Dirac law has the required first moment and the printed canonical
ratio. Concavity of `asinh(1-t)` proves feasibility of the ENTIRE mixed
family for every t, including the stronger fourth-moment residual form.

The lower estimate on the full two-trace functional reduces to the
quadratic in (13). Completing its square gives the stated Delta, and
all rational inequalities in the certificate (15) check exactly. In
particular the hyperbolic-sine bounds, square-root brackets, terminating
rational products, and strict comparison of the two squares are valid.
Because the denominator is positive and at most one, Delta gives the
uniform strict bound in (14), including the continuous eta=1 endpoint.

At zero feedback the negative-metric branch has identical spectral
traces and a larger shell factor, so it does not improve the bound.
For small positive feedback, `I-chi Q >= (1-chi)I`; positivity of the
two trace test matrices gives (16) before taking square roots. This
also applies to the negative branch and proves the stated persistence.

## Crucial scope

The Dirac law lacks the top singular value required by an actual scalar
optimum. The final source says this explicitly and does not supply a
top-mass repair, an actual optimal Gram, a complete signing, or joint
source/cross compatibility. It also correctly restricts the obstruction
to the retained scalar-affine two-trace/Jensen moment relaxation.

Thus the failed certificate is not a lower bound on an actual Gaussian
shell width, an actual optimizer counterexample, or a no-go theorem for
the desired upper comparison or original convergence. The strong and
weak branches establish precisely the limited diagnostic stated.

This independent review is analytic. File reading, line counting, and
SHA-256 checks were used; no mathematical numerical computation,
enumeration, simulation, or solver run was performed.
