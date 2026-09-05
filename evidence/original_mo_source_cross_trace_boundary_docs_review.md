# Full analytic review: actual coupling and formal all-threshold boundary

2026-09-05. Reviewer: optimized_profile_docs_gate. Analytic PASS for
the complete final source. No correction requested. This is independent
review of the new relaxation and all-threshold arguments, with my
earlier contribution to the actual coupling explicitly disclosed below.

## Complete frozen-source and prerequisite reads

I directly read ALL 444 lines of the final source
`/tmp/original_mo_source_cross_nuclear_trace_boundary.md`, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.
This was a complete end-to-end final read, not only review of the
new threshold section or reliance on an earlier synopsis.

The required proofs were also read completely and their hashes checked:

- Original nuclear and phase source, 262 lines:
  `original_mo_original_phase_spectral_moment.md`, SHA256
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- Actual full weighted normalization, 280 lines:
  `original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
  `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
- Cross singular-moment rounding, 168 lines:
  `original_mo_cross_singular_moment_rounding.md`, SHA256
  `6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756`.
- Actual pure-cross reference functional, 312 lines:
  `original_mo_small_gap_pure_cross_upper.md`, SHA256
  `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
- Shifted covariance and coefficients, 384 lines:
  `original_mo_diagonal_majorizer_cross_covariance.md`, SHA256
  `0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
- Weighted field and original drift, 381 lines:
  `original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
  `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.

I read the complete reused baseline checker, 105 lines, SHA256
`d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`,
and its complete 204-line result, SHA256
`fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
Both are canonical evidence/original_mo_scalar_template_gamma_rational_*
artifacts. Their existing pi enclosure is reused, not re-executed.
The new rational comparisons were checked analytically from the printed
fractions. No new mathematical computation, checker run, scan, solver,
optimization, or test was performed by this reviewer.

## Contribution and review provenance

Before the note was written, I independently checked the author's
direct Frobenius/nuclear transfer and supplied the directional
arithmetic-geometric-mean refinement now appearing as (2.4).
I am therefore not independent of development of that finite-coupling
step. I subsequently checked every line of its integrated argument.

I did not derive, author, or edit the formal block witness, its retained
moment checks, the centered all-metric lower bound, or the all-threshold
extension. The exact worker authored the formal calculation; root
authored Section 9's extension. My checks of those arguments below
are independent full-source reviews, not reliance on their reported
PASS judgments. The final source accurately records these roles.

## Actual finite coupling

Jensen has the correct direction: sum sqrt(r_i)>=N when
r_i=dbar/d_i and the arithmetic mean of d_i/dbar is one.
The literal off-diagonal sign-square Frobenius identity becomes an
upper bound only when nonnegative diagonal terms are added.
Its expansion gives exactly N^2(2delta+delta^2)/dbar^2.
The nuclear-norm Lipschitz bound on either n-dimensional internal
block, divided by n, has error
sqrt[(2N/dbar^2)(2delta+delta^2)], with N=2n.

The original nuclear bound gives the coefficient
kappa sqrt(m_0)(1-1/n)/(2alpha), where m_0=n/dbar^2.
The exact inverse-half-sum identity m=m_0 ell h and
ell+h=2(1+delta) imply sqrt(m_0)>=sqrt(m)/(1+delta).
This is a one-sided transfer requiring no lower bound on m_0.

The actual contraction implies A_L^2<=I-WW^T and its right-hand
analogue. Operator-monotonicity of square root, not commutation,
therefore gives both nuclear upper bounds in (2.5). The displayed
integral representation has the correct Loewner direction and covers
singular matrices by regularization. These steps prove the full finite
inequality (1.1), including the factor 2 multiplying alpha.

Feasibility gives dbar^2>=N-1, while the original source nuclear
bound gives alpha>=kappa sqrt(1-1/n)/2. Thus the uniform error
O(sqrt(delta)+1/n) is cap-free for 0<=delta<=1. No trace optimality,
small canonical gap, active cell, covariance substitution, or trimming
is needed for this actual matrix inequality.

## Formal block data and exactly retained constraints

The formal choice m=9kappa^2/16, r=2m, f=4/3, u=kappa and
alpha=2/5 satisfies u=f sqrt(m), with 0<m<1/4.
For a=sqrt(m/(1-m)), the stated full law has
mu_1=m+sqrt(m(1-m)), mu_2=2m, and mu_3=m(1+a).
Both full normalization terms in (3.1) are at most u as asserted.
The associated scalar moment label gamma=(1-a)/2 is positive;
it is not called an actual signing's canonical gap.

The commuting P,H model realizes the listed real-matrix algebraic
data: HW=0, H^2=a^2(I-P), W^2=P, and T_form^2<=I.
Its full, cross and internal laws are mutually consistent. In dyadic
rational-rank approximations, an even projection rank and balanced
complementary Hadamard eigenvalues give zero internal diagonals and
constant diagonal second moments exactly. Irrational m is treated as
a weighted trace model or limit, not an exact finite empirical rank.

This construction does NOT identify its entries with complete-signing
magnitudes, its Boolean norms with alpha/f, or an actual active state
with u. The note does not claim that every possible signing realization
of these asymptotic laws has been ruled out; none is supplied.

The source/cross nuclear constraint has the printed slack
1-(3/2)kappa^2>1/3. The internal normalized nuclear condition is
alpha>=kappa/[2sqrt(1-m)] and is satisfied. In source sqrt(n)
units, the internal nonzero eigenvalue is 1/sqrt(1-m), with
balanced signs and zero odd diagonal. Its common-variance cubic
phase condition is the same inequality; the row-variance bound is
weaker. The cross endpoint law has integral y^(3/2)dnu=m, so
its retained cubic condition is met with u=kappa. At a genuine
active original cross state, c=Phi(K) and beta(B)<=Phi(K) imply
c=beta(B); no such Boolean attainment is inferred for this model.

The formal test retains only the stated block/moment and scalar
inequalities. It does not assert additional optimal-frame, tensor,
entrywise, or active-state constraints merely by listing these traces.

## Exact reference functional and signed metrics

The noise is normalized by 2n^(3/2). Hence the desired original
comparison has target sqrt(2)alpha=2sqrt(2)/5, not f/2.
Both expectations in (6.2) evaluate the SAME endpoint law; no separate
Jensen/chord extremizers are substituted. The atom at one gives
A_s,t(1)=(1+s)/(1+t)^2. Multiplying the second term by 1-t
before taking a limit yields its finite continuous t=1 expression.

The two-vector triangle inequality has the correct lower-bound
direction. Its zero-atom coordinate satisfies a(t)^2>1/4.
The inequality sqrt(1+t^2-2st)>=1-st is valid because the squared
difference is t^2(1-s^2)>=0. The resulting affine ratio divided by
1+t takes its minimum at an endpoint, giving exactly (7.1).

I checked both integer comparisons deriving 63/100<kappa<16/25
from the reused pi interval. The subsequent rational products give
A_*^2>125721/250000>49/100, sqrt(kappa)>79/100,
1-s>369/625, C(kappa^2)>72901/125000>29/50, and
m>35721/160000>11/50. Thus

    U_nu(t)^2>40501/125000>8/25,

with exact squared margin 501/125000. This is an analytic rational
comparison, not a newly executed certificate or sampled metric test.

For negative metrics, the zero-atom first coefficient increases and
the unit-atom numerator difference is 2kappa(1-kappa)>0.
The second trace difference is the nonnegative expression printed
in Section 8. Continuity after cancellation covers both metric
endpoint limits; no singular inverse is asserted at the endpoints.

## Every shifted threshold, full original drift, and limits

The prerequisite coefficients are exactly w=1-z^2 and
k=kappa exp(-h^2), where z=|2Phi_Gauss(h)-1|.
The planar Gaussian square lies in the disk of radius sqrt(2)|h|,
so z^2<=1-exp(-h^2), w>=exp(-h^2), and k/w<=kappa.
This holds for every finite real h, including h=0; reflection covers
negative h. Consequently the standardized spectral coefficient lies
in 0<=s<=kappa^2.

The whole function C(s), not its individual summands, is decreasing
on this interval by the printed derivative. Its centered endpoint
therefore supplies the same strict 29/50 bound uniformly. In
particular the proof does not incorrectly assert A_*>7/10 at s=0.
The universal sqrt(kappa) in the completion-square term is unchanged
by the threshold; the field covariance contributes the outer sqrt(w).

For negative metrics at general s, the first unit-atom difference is
2(kappa-s)>0 and the second difference remains nonnegative.
Therefore U_(nu,s)>L_0=sqrt(40501/125000) holds for all the
stated s and signed metrics, including their continuous endpoints.

The ORIGINAL drift on the formal pure-cross face is exactly
|s_h|c/(2n^(3/2))=z f/2. The inequality
sqrt(1-z^2)>=1-z proves

    z f/2+sqrt(1-z^2)L_0>=min{f/2,L_0}>sqrt(2)alpha.

This retains the full drift instead of using only a vanishing noise
bound at large threshold. At infinite threshold the zero-noise limit
has drift f/2=2/3, still above target. The argument is uniform over
the threshold and signed-metric parameters of this SAME functional.
It makes no near-scalar metric-error endpoint interchange for actual
signings and no statement about every possible rounding function.

## Scope and disposition

The actual finite coupling remains a valid theorem. The formal example
shows failure of the listed drift-plus-ellipsoid CERTIFICATE on that
relaxation. An upper certificate larger than the target is not a lower
bound on the actual Gaussian width, the original conditional norm,
or the optimum over complete signings. No actual active optimizer or
counterexample has been constructed, and original convergence remains
open. The listed trace-only route requires additional actual information
or a different argument; no universal impossibility theorem is proved.

I edited no mathematical source or canonical file and ran no computation.
This review and the separately requested /tmp documentation proposal
are local review outputs, not claims of integration or publication.
