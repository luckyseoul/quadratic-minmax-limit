# Exact review: near-scalar transfer of the actual complete-cross gain

2026-09-05. Complete analytic PASS for the new transfer. Earlier
prerequisite contributions are disclosed below; this is not a claim
of independence from the entire prerequisite chain.

## 1. Final source and completed reads

I read all 364 lines of the final source end to end, including a complete
reread after its one-line monotonicity clarification:
`/tmp/original_mo_near_scalar_cross_spectral_gain.md`, SHA256
`ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9`.
No correction is requested.

The substantive prerequisite was completely read and previously audited
in this same task sequence, all 411 lines:
`original_mo_complete_cross_flat_spectral_gain.md`, SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
Its hash was refreshed during this transfer review.

I also refreshed the following complete sources and verified their hashes:

- Full-spectrum good-coordinate comparison, 280 lines:
  `original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
  `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
- Full compatibility source, 303 lines:
  `original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
  `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
- All-shell metric and energy stability, 252 lines:
  `original_mo_diagonal_majorizer_metric_stability.md`, SHA256
  `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.

The transfer rederives its needed good-coordinate and uniform-energy
steps, rather than importing trace optimality from the gap theorem.

## 2. Normalization and balanced auxiliary block

The contraction row-square estimate and trace Cauchy--Schwarz give
dbar^2>=N-1; the separate cap gives dbar^2<=C^2 N. Hence the
lower bound m_0>=1/(2C^2) and the stated upper bound are correct.

The literal cross sign squares give m=m_0 ell h exactly. Separate
half-wise Cauchy--Schwarz and the identity t_L+t_R=2 give ell h>=1;
arithmetic-geometric mean gives ell h<=(1+delta)^2. Therefore the
direction m_0<=m is justified and available for the later one-sided
coefficient estimate. No scalar approximation of the full W is used.

The bad-coordinate estimate follows from the exact dispersion sum.
Balancing the two good halves can at most double the total discarded
fraction, yielding theta<=2b_0. The resulting q satisfies q>=n/2,
and m'=a m_0/(1+eta)^2>=1/(9C^2). These factors use N=2n
correctly; the two halves are not assumed to lose equal original counts.

B_J is a complete square sign matrix even when the retained left and
right coordinate sets differ. Its operator bound follows from the
actual weighted subblock and the upper bound on retained diagonal
entries. This places the operator hypothesis only on the auxiliary
matrix. It imposes no bound on the original B's largest singular value.

Independent unbiased extension of both removed halves makes every
omitted bilinear term have mean zero. This proves beta(B)>=beta(B_J)
without padding by zero Boolean coordinates or asserting auxiliary
optimality. The submatrix is used solely in this norm-lower-bound role.

## 3. Return both moments to the original full weighted spectrum

The Hermitian dilation of W_J is a principal compression of the full
dilation of W after 2(n-q) coordinate deletions. Interlacing applied
to (max(x,0))^(2k), k=1,2, gives loss between zero and 2(n-q)
in its positive-power sum. These sums are exactly the cross singular-
power sums. Division by n gives the factor 2theta in (4.1).

The retained diagonal congruence factors lie in [1-eta,1]. Both Y
and W_J are contractions, and their operator-norm difference is at
most 2eta. Singular-value perturbation followed by the 2k-Lipschitz
bound gives 4k a eta after normalization by the ORIGINAL n.
Thus the factors a in (4.3)--(4.4) are indispensable and correct.

The first moment of Y equals m' by literal sign squares. The fourth
moment inequality retains a nonnegative positive part and relates
a v_2' to the FULL original v_2. There is no replacement of nu by
the smaller empirical measure in the theorem's final conclusion.

## 4. Apply the prerequisite and audit every error factor

The flatness deficit of B_J satisfies 1-epsilon'=v_2'/m'. Applying
the prerequisite and returning to the original scale multiplies its
right-hand side by exactly a(1+eta). Both the dimension ratio q/n
and scale ratio d'/dbar are retained in (5.2).

Using the positive fourth-moment lower bound gives the finite formula
(5.3). Since m'<=m, its coefficient can be lowered to kappa/m
BEFORE expanding the positive part. This is valid even when v_2
is below the moment-comparison error; it does not incorrectly magnify
a negative subtraction by the larger auxiliary coefficient.

The gain term is exactly g_kappa a^2 m_0/(1+eta). Its difference
from g_kappa m is at most
g_kappa(2delta+delta^2+2theta+eta), using the directional inverse-
mean upper bound and m_0<=1. This checks the second normalization
factor independently of the fourth-moment calculation.

At eta=delta^(1/3)<=1/8, theta<=9eta/4. Therefore the fourth-
moment error is at most 25eta/2 and the gain error is at most
6g_kappa eta. Since 1/m<=2C^2, the total dispersion coefficient
is precisely 25kappa C^2+6g_kappa as stated.

The prerequisite error is uniform in all complete auxiliary matrices
at m'>=1/(9C^2). Because q>=n/2 and a(1+eta)<=3/2, the tail
envelope R_C(n) controls it and tends to zero. No monotonicity or
unproved finite-order rate of the original error function is assumed.

At delta=0, D is exactly scalar. The original B directly satisfies
the prerequisite at dbar, so the no-trimming endpoint has no
dispersion error and is covered by the same error envelope.

## 5. Actual active-state transfer and remaining scope

The pointwise dispersion inequality controls ||Qz-z||, while
||Qz||=||z||=sqrt(N). Contractivity of the ACTUAL T gives the
uniform original quadratic-energy comparison. The block sign flip
bounds the cross pairing by its quadratic Phi norm and yields exactly
|u_D-c/(n dbar)|<=2sqrt(delta), with no trace-optimality premise.

The additional actual conditions p=q_A=0 and positive c=Phi(K)
force c=beta(B): one inequality comes from the attaining state,
the other from a block sign flip in the full paired objective. Thus
the weighted active-state inequality (7.2) follows with the extra
2sqrt(delta) loss. A representative is eligible only when selected
inside its final refined cell with those same original conditions.

The asymptotic exclusion of an endpoint cross law with u_D=kappa+o(1)
is valid under the fixed cap and delta->0, because m is bounded below.
It is a new actual-entry constraint, not a claim that the prior formal
trace relaxation had supplied an actual signing. The proof neither
forces delta->0 for all optimizers nor evaluates all original cells.

## 6. Provenance and verdict

The root proposed the balanced-good-coordinate transfer. The proof
worker authored its finite factors, positive-part comparison, and
explicit error constants. I checked the proposed formulas before
freeze and independently checked the entire final transfer above.
An optional stronger inverse-mean inequality I mentioned is not used
in this artifact. I made no edit to its source.

I contributed to the prerequisite 411-line Gaussian proof and robust
gain, authored the earlier 280-line full-spectrum and 252-line metric
notes, and reviewed the earlier compatibility result. Thus this review
discloses those dependencies rather than claiming independence from
their development. It is a complete check of the new transfer itself.

No mathematical computation, checker, solver, signing search,
simulation, local mathematical job, or remote mathematical job was
run. This review used source reads and SHA256 checks only.

Verdict: complete analytic PASS for the final source and its explicit
scope. The final moments remain those of the FULL actual weighted W.
Original all-cell inequality and convergence remain open.
