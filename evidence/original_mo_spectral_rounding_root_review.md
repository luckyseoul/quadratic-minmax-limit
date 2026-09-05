# Root review: actual spectral constraints and tensor-mixture rounding

2026-09-05. Full-read mathematical PASS in the stated finite-dimensional
and limiting scopes. Original convergence remains OPEN.

## Complete reads and authorship

The root read every line of the following six complete source proofs:

- Cross singular-moment rounding, 168 lines,
  `6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756`.
- Original phase spectral moment, 262 lines,
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- Original phase norm-only bootstrap, 176 lines,
  `3736db69d904b5a63ade46b32f6fddcc0505019f45ef483110c3ee67b24c8915`.
- Cross SDP complementarity, 293 lines,
  `21deee54db15dfd264106592e866bf57c2e954fadf86f84f87870ce9690ade70`.
- Cross tensor mixture and sign defect, 351 lines,
  `66ed4d02f4cdd7d323ae5f0c717993bc3453d59dcb3a1e00b2e2c64720d424da`.
- Evaluated scalar moment-feedback diagnostic, final 321 lines,
  `cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`.

The root also completely read every imported review receipt, and directly
checked the actual interpolation and polarization prerequisite in Section 3
of `NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`.

The root contributed phase padding and the spectral bootstrap to the
original-source pair. Those checks are collaborative-author verification,
not independent reviews. The docs worker independently read both complete
proofs and their interpolation prerequisite; the exact worker independently
reviewed the first source. The cross singular-moment and SDP-complementarity
proofs also received root complete independent checks and their respective
proof-worker and exact-worker complete independent reviews.

The tensor-mixture proof combines the exact worker's residual argument and
the docs worker's tensor construction. Its proof-worker complete review is
independent. The root supplied review feedback during development and then
checked the complete integrated proof; this receipt does not count every
collaborator as an independent author-level reviewer.

The scalar diagnostic combines docs-worker, exact-worker and root
contributions. The root completely read both the initial 286-line source
and the final 321-line version, contributing the rational certificate and
the explicit repair-trace argument. Its proof-worker full-read-plus-delta
review is independent. An optional top-mass construction proposed by that
reviewer is intentionally not part of the reviewed diagnostic.

## Cross singular-moment audit

The filtered canonical Gram is feasible for every spectrum function
`0<=f(sigma)<=sigma`, after adding separate orthogonal padding coordinates.
Its arcsine signal is `sum sigma^2 f(sigma)/n`; the squared remainder is
`sum sigma^2 f(sigma)^2/n^2`. The endpoint-valid Taylor estimate, clipped
error, exact partial-isometry SDP value, and normalized finite-order
coefficient all check. A shell magnitude cannot replace beta unless it
actually attains beta. No spectrum is replaced by an arbitrary measure
without stating the additional relaxation.

## Original-source phase and bootstrap audit

The padded positive and negative spectral phases are actual PSD Gaussian
correlations. Each uses its own variance. Zero diagonal of A removes the
padding from the quadratic trace pairing, and the unordered-pair remainder
has exactly the factor `2(1-kappa) S4_pm/v_pm^2`.

The local spectral measure has mean zero and second moment n-1. Its common
positive/negative first moment gives the sharp printed row-moment bound;
the argument does not require both extreme eigenvalues individually to
exceed sqrt(n-1). The actual diagonal-oscillation denominator and its
spectral upper bound are correctly distinguished. The common-diagonal
phases `|A|+A` and `|A|-A` give the nuclear-norm inequality with the
stated unordered-pair constant.

The bootstrap keeps the exact coefficients before any division and proves
their eventual positivity from the original norm cap. For actual minimal
phase variances, trace zero and Cauchy--Schwarz give
`v_pm >= S1^2/(2n^2)`. The nuclear lower bound then makes the common
coefficient `kappa-O_C(n^(-1/4))` under only `Phi(A)<=C n^(3/2)`.
Replacing it by kappa costs `O_C(n^(5/4))` in the ORIGINAL quadratic norm.
No fixed normalized operator bound or replacement source is needed.

## Actual SDP and complementarity audit

Strict primal and dual feasibility give the attained diagonal optimum.
Reciprocal rescaling proves block-trace balance for every optimum. The
Schur inequalities use the literal nonzero complete signing entries and
give the printed reciprocal and lower-diagonal bounds.

The canonical stacked Gram has unit rows and objective `tr|B|^3/n`.
Its gap g is an actual SDP gap, not an assumed optimizer consequence.
The weighted residual constant is 4g, while the unweighted diagonal
variance retains its maximum-diagonal/operator factor. The orthogonal
polar argument works at deficient rank. Zero canonical gap is equivalent
to equal nonzero singular values and forces the unique scalar optimal
diagonal; scalar optimality alone is not its converse.

The near-minimal SDP diffusion bounds have their stated hypotheses.
The norm-only cubic bootstrap in Section 6 retains its positive initial
coefficient `2 kappa-1` and legitimately obtains a vanishing quartic
relative error without a separate operator-cap hypothesis.

## Tensor mixture and evaluated necessary curve

The two weighted residual blocks have the correct transpose orientation.
Applying the negative-part inequality to both and adding gives `N_-<=g/2`.
For a scalar ACTUAL optimum, the exact spectral residual gives `N_-<=g/4`.
Neither bound assumes a small gap or transfers scalarity to arbitrary D.

The odd-tensor lift and canonical direct sum have total squared norm one.
Their finite Gram realization makes Gaussian rounding legitimate. For
each fixed t<1 the entire Taylor segment stays inside the stated strict
radius; its remainder and the general/scalar sign-defect coefficients
are exact. Real-entry interpolation and the canonical primal bound give
the absolute `O_(C,t)(n^(5/4))` error under only the Boolean norm cap.

The necessary limiting envelope takes the n limit first at each fixed t.
Its derivative at the reparameterized zero endpoint is 1-U, proving a
strict improvement over the cubic-only curve when c0<=U<1. The two
former tensor-endpoint constants follow with the required order of limits.
No divergent Taylor constant is silently bounded uniformly at t=1.

## Evaluated scalar diagnostic

The actual scalar-dual optimum requires both the top singular value and
an optimal Gram supported on that top eigenspace. These conditions are
not implied by the retained empirical moments. The source normalization
is the actual intrinsic mu, and the stated interval for its limiting
feedback is only a formal bound, not a realization theorem.

The explicit source triangle inequality gives
`L_A^2<=mu+q^2 n+1`. The published positive repair therefore satisfies
`tr P/a0<=4+2q^2`. The normalized trace increments in (4b) vanish at
each fixed metric gap; no original-norm cap on A is needed for this
restricted scalar-q argument. This fills the otherwise unspecified
repair-trace comparison before evaluating any uncorrected limiting law.

At strongest feedback, rank-four repaired positivity restricts the
limiting singular support. The endpoint integrand is concave, and Jensen
with the literal Frobenius moment gives a squared bound below 9/20.
The n limit precedes the metric-endpoint limit.

At weak feedback, the formal Dirac law satisfies the entire tensor-mixture
family, including its exact fourth-moment refinement. Completing the
square gives a uniform positive gap above the target for the full
two-trace/Jensen functional. The terminating-rational comparisons in (15)
are exact, not floating-point evidence. The negative metric branch and
small-positive-feedback persistence are correctly included.

This is an insufficiency result for the stated moment-only relaxation,
not a lower bound on an actual Gaussian width. The missing top singular
value of the literal Dirac law and the absence of an actual signing,
optimal Gram or compatible joint shell are explicitly recorded.

## Scope and procedure

The root also directly checked the primary statement of Theorem 1.1 in
Braverman--Makarychev--Makarychev--Naor's strict improvement of Krivine's
bound, linked in HANDOFF. It already excludes the literal old endpoint
for actual matrices. This is an imported published theorem statement,
not a claim that its full proof was independently replayed here. The
diagnostic remains a test of its expressly retained moment constraints;
its elementary constant must not be confused with the exact Grothendieck
constant. Actual-realization work must move to admissible nearby ratios
or uniform statements instead of that known excluded endpoint.

All six proofs concern actual source or cross matrices and their stated
relaxations. They do not establish a scalar majorizer for every optimizer,
small odd-spectral diagonal oscillation, favorable source/cross alignment,
or realizability of every permitted abstract moment profile. The all-shell
Gaussian evaluation and original convergence remain open.

These derivations and reviews were analytic. No signing census, solver,
simulation, numerical scalar evaluation, or other mathematical tool run
was used for this package. File hashing and the separate changed-document
scanner establish provenance and claim consistency, not the theorems.
