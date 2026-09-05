# Independent full review: cubic alignment and the scalar-template threshold

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no mathematical corrections requested.
The single remote five-comparison result has also been completely audited.

## Exact source and independence

The reviewer read all 293 lines of the initial source, independently
checked the complete Section 6 revision, and then directly read every
line of the resulting FINAL 338-line source:
`/tmp/original_mo_scalar_template_cubic_alignment.md`, SHA-256
`60037f67234fbca8c17ee90bf52c7f4346b24e5f18eb5f2c922ebbd2d9382c2a`.
The initial source hash was
`82a48cb014a76d4f1ebeb352c42fe1fff6583b3637882684e2408f11f9839e56`.
The reviewer did not contribute to the cubic theorem or its derivation.
No source edits were made by the reviewer.

Both mathematical prerequisites were also completely read and checked:

- The 324-line scalar-template completion proof, SHA-256
  `bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`.
- The 321-line scalar-moment feedback diagnostic, SHA-256
  `cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`.

The reviewer coauthored parts of the older scalar diagnostic, but had
no authorship role in this new cubic theorem or its Section 6 extension.
The present independent judgment concerns the new source and the actual
validity of the extension, not independent authorship of that prerequisite.

## Actual optimal frames and the cubic inequality

The scalar-optimality premise is essential and is stated explicitly:
the actual vector SDP equals p times q, with q the operator norm.
Equality in the operator and Frobenius bounds gives both matched-frame
relations and equality of the two frame operators. This remains valid
at deficient frame rank. The trace identities, s<=mu<=1, and the
top-eigenvector completion bound (3) have the correct normalizations.

For the cubic argument, U-transpose T equals V-transpose with precisely
the required orientation. Tensor flattening in (7) preserves Frobenius
norm and turns its squared norm into exactly the three fourth-moment
arrays in (8), including the factor two on the mixed array.
The sum of squared inner products in EACH array is p-squared times s.
Cauchy--Schwarz therefore supplies p-squared times s-squared in each
case, with no rank assumption or replacement by a formal spectral law.

The operator norm of U squared is p times mu, so the lower bound in
(9) is 4p s-squared divided by mu. The contraction estimate (10) is an
upper bound on the same quantity, not a lower bound. Its combination
with (9) gives exactly j_3>=-1+2s-squared/mu. All divisions use mu>0.

## Hermite clipping and the strengthened completion estimate

The clipped Gaussian cubic coefficient has the correct sign and size:
b_3=-2phi/sqrt(6), hence b_3-squared=2phi-squared/3. The truncated
fourth moment and cubic absolute-tail identities in the source give
E[G-cubed f(G)]=3P-2phi, as required.

Gaussian Hermite orthogonality gives the stated tensor-feature
objectives at every odd degree. The first feature attains q, every
remaining normalized feature is at least minus one by contraction,
and the cubic feature admits the extra positive term already proved.
Finite Hermite sums pass to the limit by Gaussian L2 continuity for
the fixed finite matrix. No independence of the features is assumed.

The prior completion fluctuation estimate applies to these same actual
optimal frames. Its exact degree-two cancellation and even-remainder
covariance bound give sqrt(Rs/2), not a larger or smaller factor.
Adding this estimate to the corrected cross objective proves (13).
The simultaneous use of the top-frame inequality (3) is legitimate.

## Exact threshold and the existing new checker result

The companion proof's previously verified Gaussian enclosures imply
the strict D lower bound and positive slope in q. Under the proposed
Gamma cap, substitution of mu<=x/(q-1) has the correct direction in
the positive cubic coefficient. The resulting lower coefficient is
strictly greater than 0.63 for q>=12/5. Every denominator is positive.

Putting t=sqrt(s) reduces the fluctuation-minus-cubic maximization to
an affine term minus a positive quartic. Its critical-point equation
and substitution give exactly the cube 27 R_0-squared/(1024 b_0).
The displayed rational comparison makes the maximum strictly less
than 0.067, so (13) forces x>0.4154>0.415. Thus (19) follows with the
stated strict uniform margin and no effective-rank substitution.

The reviewer read completely the 74-line exact Fraction checker
`/tmp/original_mo_scalar_template_cubic_rational_certificate.py`,
SHA-256
`bc6b5fc08120a2ed645d16ce5a4762919776853bd94d264bf5ce26d6531979ee`.
Its five comparisons match the new D lower bound, cubic coefficient,
quartic maximum, final contradiction, and endpoint threshold square.
It correctly reuses the prior P, phi, and R enclosures rather than
claiming to reestablish them. The reviewer has not run this checker.

The reviewer subsequently read every line of the 48-line execution
metadata and 51-line stored result, and independently hashed both:

- `/tmp/original_mo_cubic_alignment_rational_check.json`, SHA-256
  `d865982cd37f28d2138eeccebb76e57ae63a947761d240c32f4094f0f5c67842`.
- `/tmp/original-mo-cubic-rational.PilWFP/result.json`, SHA-256
  `4bc6760b06927a05c104123ff858546ddac1729306c9f7c94a8ad490bc91ad27`.

All five exact Fraction comparisons report PASS with the correct
relations and quantities. The proof and certificate hashes match the
reviewed sources. The metadata records one completed remote run on
soulkiller at 2026-09-05T21:31:01Z, exit zero, one worker, a 30-second
limit, empty stderr, and no remaining remote PID at its recorded check.
It explicitly records zero local mathematical runs and no rerun of
the previous 28-comparison Gaussian-enclosure certificate. The new
result is therefore an arithmetic supplement to the analytic review,
not an execution of a general theorem or matrix-realizability checker.
The reviewer did not rerun either certificate locally or remotely.

## Variable-u diagnostic and order of limits

Section 6 now explicitly stipulates the diagnostic energy sqrt(2).
It does not infer actual Boolean saturation from an upper certificate.
The active normalization gives q>=sqrt(2), u=sqrt(2)/q in (0,1], and
m=u-squared/2 in (0,1/2]. These are separate from the older diagnostic's
fixed Krivine endpoint, as the new source states.

Direct substitution into the general trace algebra of the older note
gives exactly both rational functions in (20a) and the coefficients
in (20b). That trace algebra did not use the old numerical value of u.
Because m<=1/2, the limiting Dirac denominator is uniformly separated
from zero up to eta=1. Consequently (20) follows by taking the endpoint
of this limiting expression, and the second term indeed vanishes.

The repair estimates (20c) follow from the older note's intrinsic
source trace bound for every bounded q, independently of the value
of u. Their use here fixes eta<1, evaluates the actual repaired PSD
field first, and sends n to infinity before taking eta to one.
The fixed-eta limsup/inf inequality justifies this order. No assertion
of a uniform finite-n endpoint estimate or unrepaired PSD is needed.

The new Gamma threshold implies q<12/5 only within the stated
certificate condition. Then u>5sqrt(2)/12>2-sqrt(2), with the latter
comparison equivalent to 578>576. Solving the scalar inequality for
(1-u)/(1-u-squared/2) gives exactly this threshold in (0,1]. Thus the
identified formal weak-Dirac diagnostic passes in the remaining
Gamma-certifiable range, without a claim of matrix realizability.

## Scope and procedure

The theorem bounds the finite-template COMPLETION functional Gamma.
It does not convert an upper certificate into an actual Boolean norm
lower bound, infer actual frame rank, or prove scalar optimality or
a Dirac bulk for arbitrary sign matrices. The original source/joint
comparison and original MO limit remain open exactly as stated.

This is an independent analytic proof review and read-only checker
audit. No matrix construction, mathematical program, new certificate
execution, solver, search, or optimization was run by the reviewer.
