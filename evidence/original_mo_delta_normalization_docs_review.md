# Independent full review: near-scalar actual full-spectrum normalization

2026-09-05. Reviewer: optimized_profile_docs_gate. Analytic PASS for the
entire frozen source. No mathematical or documentary correction requested.

## Frozen objects and independence

I directly read all 280 lines of
`/tmp/original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.

I also directly read the complete required prerequisites:

- Original Gaussian-phase and nuclear moments, all 262 lines of
  `/tmp/original_mo_original_phase_spectral_moment.md`, SHA256
  `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- Actual weighted compatibility, all 303 lines of
  `/tmp/original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
  `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
- Arbitrary-feasible-diagonal metric stability, all 252 lines of
  `/tmp/original_mo_diagonal_majorizer_metric_stability.md`, SHA256
  `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.

These hashes were directly verified. I had no theorem-development,
proof-writing, or pre-derivation involvement in the new source. Earlier
independent reviews of its prerequisites do not constitute authorship.
I ran no mathematical computation, signing search, solver, simulation,
test, or numerical spectral evaluation, locally or remotely. This is a
complete analytic proof review, not an execution or publication receipt.

## Feasible diagonal, trimming, and source preservation

The new theorem permits ANY positive diagonal with D+-K PSD, provided
the independent source-scale cap S=tr D<=C N^(3/2) is supplied.
Trace optimality and a small canonical-primal gap are not required.
The contraction row-square inequalities imply
S>=q tr(D^(-1)), and Cauchy--Schwarz then gives S>=N sqrt(q).
Consequently r=q/dbar^2 lies in (0,1], and the stated cap bounds r
below by a positive C-dependent constant uniformly for N>=4.

The arithmetic-harmonic identity is exact. Outside
|t_i-1|<=epsilon, its summand is at least
epsilon^2/(1+epsilon), proving the printed bound on b.
The hypotheses M>=3N/4 and N>=4 imply M>=3 and
(M-1)/(N-1)>=2/3. These verify every size and denominator
condition in the finite trimming argument.

The principal signing K_I remains a complete signing. Compressing the
two PSD majorizations and bounding D_I above gives
||K_I||op<=(1+epsilon)dbar. This does not bound ||K/dbar||op.
Extending a fixed state on I by independent unbiased signs makes
the full quadratic expectation exactly its K_I energy. Therefore
Phi(K)>=Phi(K_I), including for a negative extremizing energy.
Trimming is used ONLY in this auxiliary lower-bound mechanism.
It does not alter the actual full weighted spectrum, source covariance,
cross block, or final cells used later.

## Full actual spectral-moment transfer

For a principal compression, ordered eigenvalue interlacing bounds
the loss of the positive j-th powers between zero and N-M;
the analogous negative absolute powers have the same bound.
The full actual T has spectrum in [-1,1], so division by N proves
0<=mu_j-N^(-1)tr|T_I|^j<=2b for each j=1,2,3.
The moments on the retained side of the theorem still belong to the
FULL original T, including all exceptional coordinates.

On I, Q=diag(sqrt(t_i)) satisfies ||Q-I||op<=epsilon.
Expanding QT_IQ-T_I gives the printed operator error at most
epsilon(1+sqrt(1+epsilon))<=3epsilon. Both matrices have spectra
in [-(1+epsilon),1+epsilon]. The scalar Lipschitz constants for
|x|^j and Weyl comparison consequently give exactly
C_j=3j(3/2)^(j-1), namely 3, 9, and 81/4.
This establishes both directions of the moment estimate (3.3).

Literal sign squares give E_2=a r_I. Direct subtraction yields
r-a r_I=b(N+M-1)/dbar^2, bounded by (5/2)b for N>=4.
Combining the two estimates gives |mu_2-r|<=5b+9epsilon.
No claim that exceptional full-source eigenvalues disappear in
operator norm is used in any of these comparisons.

## Common-variance Gaussian phases and exact finite constant

The coordinate spectral measure of the auxiliary complete signing
has mean zero, second moment q_I, and support in [-L_*,L_*].
Writing its positive and negative second moments as P_2,N_2,
the probability bounds and P_2<=L_*t imply
q_I P_2<=L_*^2 N_2. Applying the same argument with reversed sign
proves the exact row-imbalance bound (4.1).

Thus the common variance
v_*=2q_I L_*^2/(L_*^2+q_I) is at least q_I+|h_i|.
Both padded correlations are genuinely PSD and unit diagonal.
Their linear ORIGINAL quadratic-energy difference is exactly
kappa tr|J|^3/v_*. The sum of their arcsine remainder magnitudes
is at most 2rho tr|J|^4/v_*^2. Both expected energies are in
[-Phi(J),Phi(J)], proving (4.2) with its factor 2Phi(J).
No rectangular-to-quadratic transfer or unequal-variance cancellation
is hidden here; padding affects neither off-diagonal remainder bound.

After division by S=N dbar, the cubic coefficient is exactly

    kappa[(1+epsilon)^2+r_I]/[2r_I(1+epsilon)^2].

The lower estimate E_3>=[mu_3-2b-(81/4)epsilon]_+ is valid because
E_3 is nonnegative. The bound tr|K_I|^4<=M L_*^4 gives exactly

    rho a[(1+epsilon)^2+r_I]^2/[2dbar r_I^2]

for the retained fourth-moment error. This verifies the entire finite
inequality (4.3), including all powers and factors. The fixed cap,
r_I>=2r/3, and dbar>=sqrt(N-1) make this error O_C(N^(-1/2)).

## Nuclear phases and denominator direction

The PSD matrices |J|+-J share the positive diagonal ell_i.
A zero ell_i would force a zero row of |J|, contradicting
(J^2)_ii=q_I>0, so normalization is legitimate even for singular J.
The signed edgewise correlation difference is
2/sqrt(ell_i ell_j), and the arcsine slope is at least one.
Their two original expectations therefore imply
2Phi(J)>=2kappa sum_(i<j)(ell_i ell_j)^(-1/2).

Arithmetic-geometric mean and Cauchy--Schwarz give
sum_(i<j)(ell_i ell_j)^(-1/2)>=M^2 q_I/(2tr|J|).
This checks the exact coefficient kappa M^2 q_I/tr|J| in (5.1).
Dividing by S converts it to kappa a^2 r_I/E_1. The UPPER
bound E_1<=mu_1+3epsilon gives a valid lower bound on this ratio,
proving (5.2) in the correct denominator direction. Its denominator
is positive. The phase and nuclear inequalities hold simultaneously,
so taking their maximum introduces no extremizer assumption.

## Uniform limit, paired use, and remaining scope

For 0<delta<=1/512, epsilon=delta^(1/3)<=1/8 gives
b<=(9/8)delta^(1/3)<1/4 and satisfies the finite hypotheses.
The second-moment bound proves mu_2=r+O(delta^(1/3)).
Because mu_1>=mu_2 and r has a C-dependent positive lower bound,
the nuclear denominator is uniformly separated from zero for
sufficiently small delta. All remaining coefficient changes are
O_C(delta^(1/3)); mu_j<=1 controls their products. This proves

    U>=max{kappa(1+r)mu_3/(2r), kappa r/mu_1}
                              -O_C(delta^(1/3)+N^(-1/2)).

For delta=0 the untrimmed scalar argument applies directly. The
alternative fixed-epsilon limit followed by epsilon tending to zero
is controlled by the explicit finite bounds, not an operator-norm
assertion or an unproved endpoint exchange.

Only for EXACT scalar D does the note substitute
mu_3=r(1-gamma). This identity follows by direct normalization of
tr|K|^3; no such identity is asserted for near-scalar D with outliers.
The resulting positive-gap scalar lower bound is not a proof that
every positive-gap source has near-scalar D or that every such
source satisfies the desired conditional field upper.

For the literal paired source, the quoted uniform energy comparison
uses feasibility and delta alone. This is explicit in the full
metric prerequisite, even though the earlier gap theorem additionally
assumes trace optimality to obtain delta from gamma. Thus the bound
|u_D-2c/S|<=2sqrt(delta) is available here without small gamma.
On the separate active original face p=q_A=0, c=Phi(K), it transfers
the U lower bound to the actual within-final-cell u_D. The extra
O(sqrt(delta)) error is absorbed by O(delta^(1/3)). General shells
must retain c/Phi(K); saturation is not imposed universally.

The actual weighted source and cross measures remain to be coupled
in the field upper. The separate metric theorem retains its fixed
interior-parameter limit order. No full-width evaluation, positive-gap
closure, or original MO convergence follows from this normalization
alone. The source expressly retains those unresolved implications.

No source or canonical file was changed by this reviewer. This
receipt is the sole new /tmp output; no mathematical run was needed.
