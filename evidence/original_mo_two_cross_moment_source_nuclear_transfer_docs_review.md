# Independent frozen-source review: two-cross-moment nuclear transfer

2026-09-06. Reviewer: optimized_profile_docs_gate.

## Complete reads and contribution boundary

I directly read the ENTIRE frozen 255-line source
`/tmp/original_mo_near_scalar_two_cross_moment_source_nuclear_transfer.md`,
SHA256
`32c1e47608c1dc06037ababb6c3c34934fdea5546df17527e94f6509562e6525`.
I also directly read its complete 52-line author receipt, SHA256
`435a409ff60e3e523bd4e4c88dbefe8de0b623d1f17ad07b4527b31d564bdd93`,
at the matching `_author_receipt.md` /tmp path. Both hashes and line
counts were checked directly.

Both named prerequisites were directly read in full and independently
reviewed in this immediately preceding review sequence; their unchanged
hashes were refreshed for this task:

- `original_mo_near_scalar_cross_endpoint_source_nuclear_transfer.md`,
  230 lines, `6a486df0fd46aa76259e3f02e3734eb2529162500f98f89af58e90562e6a2187`;
- `original_mo_all_law_adaptive_nuclear_gain.md`,
  553 lines, `0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51`.

I previously checked root's proposed finite Schatten inequality and the
exact worker's reported threshold and rectangle before this source was
written. Those were checks of their already supplied arguments; I supplied
no new derivation, correction, parameter choice, or proof step. This is
my subsequent full independent check of the frozen new source, not a
substitute receipt based on that outline alone. My older scalar209 and
source/cross AM-GM contributions remain disclosed in the prerequisite
reviews; independence of every historical premise is not claimed.

## 1. Finite noncommuting bound

From 0<=H^2<=I-Y and 0<=Y<=I, trace positivity gives
0<=z=tr(H^2Y)/n<=min(M,d). In particular z<=d uses the trace
pairing of the PSD matrices (I-Y)-H^2 and Y, not commutation.

The exact decomposition H=H(I-Y)+HY factors the two terms through
sqrt(I-Y) and sqrt(Y). Schatten Cauchy--Schwarz and the trace-norm
triangle inequality give the four squared Frobenius factors
n(M-z), n(1-v), nz, and nv. Thus (2.1) has precisely the stated
two square-root terms. No Boolean norm is substituted here.

For fixed M,v the concave right side increases up to z=vM and
decreases afterward, attaining sqrt(M) there. Maximizing subject to
z<=d therefore gives z_*=min(d,vM). This point is at most M,
so its square roots are defined. The formula extends continuously
to zero M and the endpoint values of v. This proves all of Section 2.

## 2. Full moments and the same original source

The complete cross sign squares give m_D=m_0 e_L e_R. The exact
reciprocal half-means satisfy 1+delta<=e_L e_R<=(1+delta)^2,
with the stronger lower obtained from 1/e_L+1/e_R<=2. Thus only
m_D->m>0 and delta->0 are needed to recover m_0->m and the scale.
No endpoint law is used in this step.

The restated common-label trimming retains q/n->1, handles delta=0
explicitly, and yields the same complete A_J on both halves. Feasibility
and congruence give the 3epsilon eigenvalue comparison and retained
operator cap, without an operator assertion on the untrimmed A.

Positive/negative square interlacing gives the normalized 2b error;
the square-function Lipschitz bound and tr(A_J^2)=q(q-1) give (3.4).
Consequently the FULL internal second moment tends to m. The block
trace inequality then gives m<=1/2. The actual contraction spectrum
also gives 0<=Delta_D<=m_D(1-m_D), proving the stated limit range.

Continuity of the finite bound gives S(m,m,Delta). Dividing by
sqrt(m) produces exactly the piecewise L: below Delta=m^2 the
formula is the two displayed roots, and at or above it the cap is one.
The two branches agree at the boundary.

Unnormalized nuclear norm decreases under principal compression;
the congruence comparison contributes 3epsilon per normalized trace.
Together with dbar/sqrt(q)->1/sqrt(m), this proves the source nuclear
cap. Independent unbiased extension proves Phi(A_J)<=Phi(A).
Feasibility bounds Phi(A)/n^(3/2), justifying the normalized o(1)
comparison. All input moments remain those of the FULL actual H and Y.

## 3. Exact certified threshold and interior rectangle

On 0<=Delta<=m^2, the angles beta and theta satisfy
0<=theta<=beta and L=cos(beta-theta). The starting value sqrt(1-m)
is at most 4/5 exactly when m>=9/25. In that range inversion gives
sqrt(Delta/m)<=(4sqrt(m)-3sqrt(1-m))/5, with a nonnegative right
side. Its square is precisely Delta_crit(m).

That threshold is strictly below m^2, so the constant branch L=1
cannot add solutions. Together with C_m<=5/3 and m<=1/2, this proves
the exact region R for the TWO CERTIFIED caps. It does not assert
necessity for excluding actual sources by other arguments.
The endpoint values zero at m=9/25 and 1/100 at m=1/2 are correct.

In the rectangle m>=2/5, Delta<=1/1600, one has Delta<m^2.
L increases in Delta and decreases in m, the latter because the
first squared term has derivative -1+Delta/m^2<0. Its worst corner
is therefore sqrt(1917/3200)+1/40, strictly below 4/5 since
(31/40)^2=1922/3200. This verifies the full rectangle, not just a point.

## 4. New monotone functional composition

Differentiating the imported F_C using partial_w Psi<=1/2 gives
exactly (5.2). For C>=1 and 0<u<=1, its negative term is at most
-kappa/2 and its positive term is less than (1-kappa)/2. Since
kappa>1/2, F_C is strictly decreasing on the whole required interval.

The retained source's nuclear moment is at most sqrt((q-1)/q) and
at least ((q-1)/q)/(||A_J||/sqrt(q)). Hence its subsequential limits
are positive and no larger than L. Apply the uniform all-law theorem
first at each FIXED C=C_m+eta, then let eta tend to zero. Continuity
and decreasing dependence on u give F_(C_m)(u)>=F_(C_m)(L).
The original norm comparison proves (1.4), with no moving-cap error
assumption. Inside R, the independently proved two-cap consequence
gives the displayed constant 2/5+7/55000 directly.

The actual moment pairs lie in a compact bounded set. If every
accumulation point belongs to R, every subsequence relevant to the
original liminf has a further positive-moment limit in R, to which
the preceding theorem applies. This proves the no-moment-convergence
version. Under liminf m_D>=2/5 and limsup Delta_D<=1/1600,
the same extraction and near-scalar full-second-moment argument force
m<=1/2, so all such limits lie in the proved rectangle. Neither
formulation assumes a full empirical cross or internal spectral law.

## Verdict and remaining scope

PASS for the complete frozen 255-line source, with no required correction.
The theorem strengthens the actual source transfer from a full endpoint
law to two moments, while retaining positive diagonal feasibility and
dispersion tending to zero. It does not establish those hypotheses for
arbitrary original minimizers, and R is sufficient rather than necessary
for actual exclusion. No source, covariance, cross block, active field,
or optimizer is replaced. Complementary cells and original MO convergence
remain OPEN. The frozen 553-line and 230-line sources are unchanged.

No mathematical program, checker, solver, numerical evaluation,
construction, or search was run. Tools were used only for complete
reads, hashes and line counts, and this /tmp receipt. I made no
canonical edit and performed no publication or backup operation.
