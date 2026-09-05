# Independent full review: evaluated actual-measure pure-cross upper

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no mathematical corrections requested.

## Exact source reads and staged-versus-final provenance

The reviewer directly read all 300 lines of the initially frozen source
`/tmp/original_mo_small_gap_pure_cross_upper.md`, SHA-256
`b80a8ab8cb765d7795958e53d44d982506439217971094d25929969e8e9b9579`,
and every line of the fixed 81-line certificate, SHA-256
`10d76c46fbdf75d8b856d06bae07a3d6304c78ce2d5b17de225567435f63fdf2`.
The reviewer had no involvement in either derivation and made no source
or checker edits.

After the single remote run, the reviewer read the COMPLETE three-hunk
diff between the retained staged proof and the final proof. The final
312-line source has SHA-256
`035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
Its changes only name the unchanged checker hash, record the completed
execution and exact result, and identify the full original-phase
dependency. There are no mathematical or checker changes.

The staged `/tmp/original-mo-pure-cross-rational.w3EsHK/proof.md`
still has the original b80a8ab8 hash above. The execution metadata
correctly records THAT staged hash, not the later documentary version.
This receipt is a full initial read plus a complete final-delta review.

All four required analytic prerequisites had also received the reviewer's
complete independent reads: gap303 at `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`,
metric252 at `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`,
original-phase274 at `1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`,
and weighted-field381 at `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.

## Positive actual field and exact f, u, m normalization

The pure-cross field is M_0=wn(I-a u L_D), with u=c_D/n and a=k/w
when w>0. It is positive semidefinite by the actual contraction and
|u|<=1. The w=0 case is separate. The prior gap comparison applies
because the representative belongs to the FINAL refined original-zero
cell. The metric-stability proof applies to M_0 because it uses only
positivity, constant diagonal wn, and the actual contraction and radius.
It does not need the deleted internal blocks or a scalar-D source.

Symmetrizing the two spectral signs of each actual singular value gives
exactly A_s,t and B_s,t in (2.1). The combined trace has the additional
factor t. After insertion in the Boolean ellipsoid upper, its coefficient
therefore becomes t, not sqrt(t). Both trace contributions yield the
normalization 2n^(3/2) in (2.2). The unsymmetrized factors establish
their nonnegativity, including zero singular values and endpoints.

The empirical measure is that of the n squared singular values of the
ACTUAL W, including zeros. The exact identity B_ij-squared=1 gives
m=(sum_left d_i-inverse)(sum_right d_i-inverse)/n. It is not an
independently chosen first moment or a formal Dirac-law assumption.

The inverse square-root dispersion estimate and Cauchy--Schwarz prove
(2.4). Each half's inverse-diagonal sum is consequently
(n/dbar)(1+o(1)). Since eta_source>=1 bounds n/dbar-squared,
the resulting additive error in m is o(1). The exact finite conversion
n/dbar-squared=n/[eta_source-squared(2n-1)] gives (2.5).

Here f_n is c/n^(3/2), not c/N^(3/2). Uniform cross compatibility gives
u=c/(n dbar)+o(1). Since f_n tends to a positive value, squaring and
dividing proves m=u-squared/f_n-squared+o(1), with all factors of two
correct. This is an actual matrix relation, uniform over the indicated
cell choices, not an assumed independent spectral budget.

## Fixed metric: actual-measure curvature bounds

At t=3/5, subtraction of A_s from one gives exactly the stated alpha
and beta. Differentiating twice yields the affine numerator in Section 3.
Both endpoint formulas are correct and strictly positive on 0<=s<=2/3,
so A_s is concave. Its endpoint values then give the lower bound 25/64.

The second derivative of B_s has the stated two endpoint numerators,
both positive on the same range. Thus B_s is convex. Concavity permits
Jensen for A_s, while convexity puts B_s BELOW its endpoint chord.
The chord coefficient is exactly (297-375s)/128. These two upper
estimates hold simultaneously for every actual measure, even though
one measure need not attain both. Equation (3.5) follows without a
spectral ansatz or any optimization over matrices or measures.

## Monotonic diagnostic expression

Along s=kappa u and m=u-squared/2, expansion of 1-C(u) gives the
polynomial G and differentiated numerator (4.2) exactly. The competing
kappa terms have the stated lower bound using u-cubed<=u, while the
remaining terms are nonnegative. The coarse interval 3/5<kappa<2/3
therefore proves C'(u)<0 on the required interval.

Together with C(u)>=25/64, differentiation of the first square root
gives its derivative at most -3/16 for u<1. The polynomial V is at
least one. Its derivative is bounded by 189/256 using u>=kappa and
kappa-squared>9/25, treating a negative numerator directly. The
second square-root derivative is at most 189/1280. Their sum is
strictly negative. Continuity gives the endpoint at u=1, so U(kappa)
is the maximum on the entire stated diagnostic interval.

## Exact rational enclosure and the existing one-run result

The checker reuses the already-proved pi interval, with baseline checker
SHA-256 `d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`
and result SHA-256 `fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
It does not rerun the prior Machin computation.

Every interval direction in A_up and B_up is correct: the subtracted
numerator product uses lower positive endpoints, the added product uses
upper endpoints, and the positive denominator uses its lower endpoint.
The coarse kappa interval also keeps 1-r m_+ positive. The chord
coefficient stays positive, so its endpoint substitutions preserve the
upper direction. Both square-root bounds use positive rational endpoints.
The checker maps all eleven prescribed comparisons to the proof:
coarse interval, two positivity checks, two squared-root enclosures,
their exact sum, its strict squared margin, and four analytic constants.

The reviewer completely read the 25-line stored result and 77-line
execution metadata, and independently checked their hashes:

- `/tmp/original-mo-pure-cross-rational.w3EsHK/result.json`:
  `0ea064435322e698b8e33a4d9bce8ab29156e3cfe013c9885f1f35e205156e41`.
- `/tmp/original_mo_small_gap_pure_cross_rational_check.json`:
  `bd2d6eda56412fb4a0788bfc68388bdfbffacec3860039096d8c4b77919864af`.

All eleven stored checks report PASS. Their total is 17677/25000,
equal to 70708/100000, with exact squared margin 23671/625000000 below
one-half. This proves the strict rational upper used in (5.3).
The metadata records one bounded one-worker remote run on soulkiller
at 2026-09-05T22:18:45Z, exit zero, empty stderr, and an absent recorded
PID at the post-run check. It records zero local mathematical runs and
no baseline rerun. The checker hash matches the completely read source.
The reviewer did not execute either checker locally or remotely.

## Actual small-gap application and global limits

The newly named original-phase theorem supplies u>=kappa-o(1) ONLY
under its actual small-gap premise and the separately active original
pure-cross condition p=q_A=0, c=Phi(K)>=0. A norm cap or diagonal
majorization alone does not supply that condition.

The additional normalization Phi(K)/n^(3/2) tending to sqrt(2) gives
m=u-squared/2+o(1). At fixed t=3/5 the denominator stays uniformly
positive on the compact spectral support. Uniform continuity controls
the small u and m discrepancies, including u just below kappa by o(1)
and the endpoint u=1. The previous positive-field, metric, bin, selection,
and Gaussian-padding errors remain lower order. This proves (6.1) for
the indicated actual pure-cross cell, with its strict rational margin.

For other positive f, the theorem retains the general expression with
m=u-squared/f-squared+o(1). It does not assert that expression is at most
f/2. In particular the needed conditional target is 2sqrt(2)Phi(A),
which can be smaller than sqrt(2)n^(3/2) if the original source constant
is below one-half. Neither this smaller target, the other original
internal-energy cells, nor the complementary large-gap branch is
settled. No actual law is replaced by a formal Dirac measure, and no
original MO convergence claim is made.

No mathematical computation, signing construction, numerical SDP,
optimization, search, simulation, census, solver, or new test was run
by this reviewer. Only existing outputs were read and audited.
