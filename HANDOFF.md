# Handoff: original convergence problem

Updated 2026-09-05. Start with `CORE.md` and `STATUS.md`.
The original MO limit is OPEN; `L=1/2` is also OPEN.

## Preservation and reset

The previous residual worktree, including all 48 uncommitted files, is
preserved at `archive/2026-09-05-paley-research`:
`ad8c6920412af0b3c23629afe2a9e95060c5471e`.
The separate 22-file dirty main snapshot is preserved at
`archive/2026-09-05-main-local-edits`:
`c2e13218cceb7e1fb36de8f2625bf4c4a7c0a606`.
The active checkout is `/home/nick/quadratic-minmax-limit`, branch `main`.
See `ARTIFACTS.md` for the exact scope and replay of historical documents.

The canonical entry documents no longer enforce Paley residual (ii),
gap-two optimality, a conjectured value, or a particular amplification
construction. The original-problem status is separate from route-local
predicates. Valid local results retain their stated scope and evidence.

## Next mathematical work

The newer [all-law source-gain milestone](evidence/original_mo_all_law_source_gain_milestone.json)
contains the [uniform adaptive gain](evidence/NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md).
For ACTUAL complete symmetric zero-diagonal A, the two caps
`limsup ||A||op/sqrt(n)<=5/3`, `limsup tr|A|/n^(3/2)<=4/5` imply
`liminf Phi(A)/n^(3/2)>=35/88+3/1250=2/5+7/55000`.
No limiting spectral law or diagonal homogeneity is assumed; the update
probabilities adapt to actual local fields on the same Boolean source.
The [two-moment transfer](evidence/NOTE_2026-09-06_TWO_CROSS_MOMENT_SOURCE_NUCLEAR_TRANSFER.md)
uses ACTUAL K=[[A,B],[B^T,-A]], positive diagonal D=diag(D_L,D_R)
with D+/-K>=0, and `delta=tr(D)tr(D^(-1))/(2n)^2-1->0`.
Put W=D_L^(-1/2) B D_R^(-1/2), Y=WW^T, and retain the FULL moments
`m_D=tr(Y)/n`, `Delta_D=tr[Y(I-Y)]/n`.
If EVERY accumulation point (m,Delta) belongs to
`R={9/25<=m<=1/2, 0<=Delta<=m[4sqrt(m)-3sqrt(1-m)]^2/25}`,
a common large original principal source supplies both caps above;
its norm comparison transfers the same strict lower back to A.
No full cross or internal law, extra trace cap, optimal diagonal, or
active-face premise is required. A simple sufficient condition is
`liminf m_D>=2/5`, `limsup Delta_D<=1/1600`, with delta->0 retained.
The [endpoint-law transfer](evidence/NOTE_2026-09-06_CROSS_ENDPOINT_SOURCE_NUCLEAR_TRANSFER.md)
is preserved as the Delta=0 special case. The older near-flat theorem
below keeps its larger gap under its narrower law hypothesis.
The ACTUAL source operator cap and paired diagonal-dispersion/moment
premises are not established for arbitrary optimizers; the paired field
is not replaced. No mathematical run was used. Other source regions,
the all-cell implication, and original convergence remain OPEN.

The latest reviewed actual-source milestone (2026-09-06) is
`evidence/original_mo_original_source_strict_gain_milestone.json`.
Its full proof is
`evidence/NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md`
(SHA-256 `7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`).
For ACTUAL complete symmetric zero-diagonal signings A, assume
`limsup ||A||op/sqrt(n)<=5/3` and the FULL empirical eigenvalue law
of A/sqrt(n) tends to
`(9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4))`. Then, kappa=2/pi,

`liminf Phi(A)/n^(3/2)>=5kappa/8+16/3125>2/5+3/1100`.

This is an ORIGINAL same-order quadratic-norm lower bound. It uses a
robustly normalized positive spectral projector, uniform higher-chaos
mean variance at least 1-kappa-o(1), and trace-of-square first-chaos
alignment. The Gaussianization extension proves a joint limit of ONE
local field and ONE distinguished Gaussian input coordinate, uniformly
over rows; it does not claim a growing-dimensional joint field limit.
An actual independent-coordinate Boolean update with FIXED probability
1/10 improves the positive phase. Its penalty retains the actual 5/3
operator cap, not the limiting nonzero atom 5/4. Weak empirical flatness
does not imply exact finite-order flatness, a large exact kernel, or
absence of spectral outliers. The distribution-free scalar support is
`evidence/NOTE_2026-09-06_ORIGINAL_SOURCE_LOCAL_UPDATE_SCALAR_GAIN.md`
(SHA-256 `7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`).

The separate cap-free transfer is
`evidence/NOTE_2026-09-06_NEAR_SCALAR_INTERNAL_FLAT_LAW_TRANSFER.md`
(SHA-256 `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`).
For the ACTUAL paired K=[[A,B],[B^T,-A]], take ANY positive diagonal
D=diag(D_L,D_R) with D+/-K>=0. If
`delta=tr(D)tr(D^(-1))/(2n)^2-1->0` and the FULL actual law of
H_L=D_L^(-1/2) A D_L^(-1/2) tends to
`chi_*=(9/25)delta_0+(8/25)(delta_(3/4)+delta_(-3/4))`,
one common original principal A_J, q=|J| with q/n->1, has the
law at 0 and plus or minus 5/4 for A_J/sqrt(q),
limsup ||A_J||op/sqrt(q)<=5/3, and satisfies
`Phi(A_J)/q^(3/2)<=Phi(A)/n^(3/2)+o(1)`.
Completeness and the full second moment force dbar/sqrt(q)->5/3,
dbar=tr(D)/(2n); no separate trace cap or trace optimality is needed.
No cross-law or active-state condition is required. The auxiliary
source only lower-bounds the original norm and never replaces the
paired covariance, cross block, or active field.

Composition excludes the specified ACTUAL near-scalar internal-law
regime with Phi(A)/n^(3/2)->2/5 underlying the strengthened formal
profile below. The older formal certificate boundary remains correct
for its explicitly listed relaxation, which did not include this new
original-source entry constraint; no actual signing was supplied there.
The following older checkpoints are retained as history. Neither small
delta nor chi_* is proved for arbitrary candidate extremizers. Other
actual profiles and the remaining all-cell implication still require
work; unchanged trace-only scans do not establish those implications.
No mathematical run was used. Original convergence remains OPEN.

The latest paired milestone is
`evidence/original_mo_weighted_cross_gain_boundary_milestone.json`.
It completes the weighted transfer previously marked unpublished in
the older scalar-gain checkpoint below, which is retained as history.
Its two proof sources are
`evidence/NOTE_2026-09-05_NEAR_SCALAR_CROSS_SPECTRAL_GAIN.md`
(SHA-256 `ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9`)
and `evidence/NOTE_2026-09-05_STRENGTHENED_TRACE_PROFILE_BOUNDARY.md`
(SHA-256 `903ac72c78c60706fbcfef09e50abeda0a18fe05976e3efab89d65becdbfccf1`).

For the ACTUAL paired signing K, ANY positive feasible D is allowed
with the SEPARATE cap S=tr D<=C N^(3/2), N=2n. Put dbar=S/N,
delta=S tr(D^(-1))/N^2-1, and let nu be the FULL actual squared-
singular-value law of W_D, including zeros. With m=integral y dnu,
v_2=integral y^2 dnu and g_kappa=sqrt(kappa)-kappa, kappa=2/pi,
uniformly for 0<=delta<=1/512 the actual original cross norm satisfies

`beta(B)/(n dbar)>=kappa v_2/m+g_kappa m`
`-[25kappa C^2+6g_kappa]delta^(1/3)-R_C(n)`, R_C(n)->0.

Here m>=1/(2C^2). No trace optimality, small canonical gap, maximum
diagonal bound or global unweighted operator cap on B is assumed.
A balanced complete submatrix is only an auxiliary ORIGINAL-norm
lower bound; interlacing and congruence return its second and fourth
singular powers to the full original nu. It never replaces W_D,
the covariance, source, or active cell. The uniform marginal CLT
prerequisite and its tail envelope are retained; no finite-n rate is claimed.
For the SEPARATE actual active conditions p=q_A=0 and c=Phi(K),
the same lower holds for u_D=c_D/n with an additional 2sqrt(delta)
loss. Thus near-flat full weighted laws exclude u_D=kappa+o(1)
in this fixed-cap, delta-to-zero branch, without the earlier operator cap.

The paired FORMAL diagnostic then tests a changed profile:
alpha=2/5, f=4/3, u=4/5, m=9/25 and
nu=(16/25)delta_0+(9/25)delta_1. Its specified full/source/cross
conditions, INCLUDING the new entry gain, all pass. Nevertheless the
same reference functional has U_s(t)>71/125>2sqrt(2)/5.
Adding the full ORIGINAL drift z f/2 keeps the certificate above
the target for EVERY signed metric and shifted Gaussian threshold,
including endpoints. One global supporting-line inequality supplies
the squared target margin 41/15625; the kappa enclosure was reused.

This is a lower bound on the FORMAL UPPER certificate, not on actual
Gaussian width or a Boolean norm. No complete signing or actual active
optimizer realizing the profile is supplied. The transfer theorem is
not retracted, and the earlier formal obstruction retains its scope.
A next attempt must add source-entry, Boolean-active-state, frame or
conditional-optimality information, or change the upper argument;
unchanged trace-only threshold/metric scans do not supply that step.
Small delta for optimizers and the all-cell original upper remain open.
No new mathematical run was used. Original convergence remains OPEN.

The latest actual-entry restriction is
`evidence/original_mo_complete_cross_flat_spectral_gain_milestone.json`
and `evidence/NOTE_2026-09-05_COMPLETE_CROSS_FLAT_SPECTRAL_GAIN.md`
(SHA-256 `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`).
For an ACTUAL n by n complete sign matrix B, choose the SEPARATE
unweighted operator bound d>=||B||op and put
`m=n/d^2>=m_0>0`, `epsilon=1-tr[(B^T B)^2]/(n^2 d^2)`.
Uniformly in these data as n tends to infinity, with kappa=2/pi,

`beta(B)/(nd)>=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-o_(m_0)(1)`.

Here 0<=epsilon<=1-m. Exact flatness of all nonzero singular values
at d is equivalent to epsilon=0; exact or asymptotic flatness therefore
forces a strictly positive leading-order gain over kappa for fixed m_0.
The proof retains the actual complete-entry higher-Hermite variance gain
and supplies its necessary Gaussianization bridge: mixed contractions,
a Gamma fluctuation bound and a characteristic-function equation prove
a uniform scalar absolute-moment lemma. Uniform second moments and an
L2 Hermite-tail estimate justify the limit passage. This is not an
absolute-moment inference from variance alone or a growing-dimensional
joint column CLT. No finite-n error rate or mathematical run is claimed.

For the SEPARATE actual active conditions p=q_A=0 and c=Phi(K),
one has c=beta(B); at scalar scale d the same lower applies to u=c/(nd).
Consequently the earlier FORMAL flat cross endpoint with u=kappa is
excluded in this actual scalar, bounded-operator setting. The formal
trace-relaxation certificate obstruction remains valid on its own stated
relaxation; it is not an actual-signing counterexample.

The unweighted d>=||B||op and m>=m_0 hypotheses are NOT obtained here
from small dispersion delta or trace control of a diagonal majorizer.
Transfer to the actual near-scalar weighted cross law is the next
unpublished implication, being treated separately. This note neither
replaces W_D nor evaluates every active cell, and original convergence
remains OPEN. Preserve the older scoped proofs and their exact premises.

The latest actual coupling and formal certificate boundary are recorded in
`evidence/original_mo_source_cross_trace_boundary_milestone.json` and
`evidence/NOTE_2026-09-05_SOURCE_CROSS_NUCLEAR_TRACE_BOUNDARY.md`
(SHA-256 `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`).
For ANY positive feasible D, let N=2n, dbar=tr(D)/N,
delta=tr(D)tr(D^(-1))/N^2-1, alpha=Phi(A)/n^(3/2), and let nu be
the ACTUAL squared-singular-value law of W_D, with mean m. The finite
source/cross inequality is

`integral sqrt(1-y)dnu >= [kappa sqrt(m)/(2alpha)](1-1/n)/(1+delta)`
`-sqrt[(2N/dbar^2)(2delta+delta^2)]`, kappa=2/pi.

This needs no trace cap, optimal diagonal, small canonical gap, or active
cell. Its O(sqrt(delta)+1/n) consequence keeps the original matrices;
the direct nuclear transfer requires no auxiliary trimming.

The SAME note also identifies a FORMAL limitation of the listed trace/block
relaxation. At alpha=2/5, f=4/3, u=kappa and
nu=(1-m)delta_0+m delta_1, m=9kappa^2/16, the retained full/source/cross
moment inequalities and block contraction all pass. Nevertheless EVERY
shifted Gaussian sign threshold and signed ellipsoid metric, including
endpoint limits, has formal drift-plus-certificate value above the target
sqrt(2)alpha. The full ORIGINAL drift z f/2 is retained, with
z=|2Phi_Gauss(h)-1|; the certificate's noise term is bounded below by
sqrt(1-z^2)L_0, L_0^2=40501/125000>8/25.

These formal parameters are not realized by complete signings or actual
active Boolean states in this note. A lower bound on this UPPER certificate
is not a lower bound on actual Gaussian width or the original norm.
The next step within this route is an additional actual entry/active-state
constraint, or a different upper argument, beyond the listed trace data.
Do not repeat threshold/metric scans on the unchanged relaxation. The pi
enclosure was reused analytically; no new mathematical run was needed.
Original convergence remains OPEN; this is not a signing counterexample
or an impossibility theorem for all methods.

The latest delta-only normalization is
`evidence/original_mo_delta_normalization_milestone.json` and its proof
`evidence/NOTE_2026-09-05_NEAR_SCALAR_DIAGONAL_SPECTRAL_NORMALIZATION.md`
(SHA-256 `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`).
It permits ANY positive feasible D with the SEPARATE cap
`S=tr D<=C N^(3/2)`: neither trace optimality nor small canonical gap
is required. With dbar=S/N, r=(N-1)/dbar^2 and mu_j the j-th absolute
moment of the ACTUAL full weighted contraction T, it proves

`2Phi(K)/S>=max{kappa(1+r)mu_3/(2r),kappa r/mu_1}`
`-O_C(delta^(1/3)+N^(-1/2))`, with `mu_2=r+O(delta^(1/3))`.

The proof's principal signing is only an auxiliary ORIGINAL-norm lower
bound. Interlacing and congruence transfer its moments back to the FULL
actual T, including exceptional coordinates. It never replaces the
source, covariance or W_D. The finite phase and nuclear inequalities
are (4.3) and (5.2); no mathematical execution was needed.
For EXACT scalar D, mu_3=r(1-gamma), so the phase term remains useful
at positive gamma where the previous masked-gap estimate was vacuous.
That identity is not silently imposed on near-scalar D with outliers.

The same lower transfers to u_D=c_D/n only on the separately active
original face p=q_A=0, c=Phi(K)>=0, at error 2sqrt(delta). Small delta
is still a hypothesis, not established by source or conditional
near-minimality. The next use must couple actual full and cross spectra
at the actual norm scale, rather than enlarge scalar D and recover the
old attenuation loss. Smaller-normalization and all-cell width estimates
remain open, as does the original convergence problem.

The latest evaluated small-gap package is
`evidence/original_mo_small_gap_evaluation_milestone.json`:

- `NOTE_2026-09-05_FULL_SDP_GAP_ORIGINAL_PHASE_BOUND.md`
  (SHA-256 `1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`)
  constructs two ACTUAL unit-diagonal PSD Gaussian phases with the same
  coordinate normalization. Subtracting their ORIGINAL quadratic energies
  gives the coefficient kappa/2 directly, not through beta(K)<=4Phi(K).
  Here D is trace-optimal for the FULL SDP, S=tr D=tau(K), q=N-1,
  and gamma=(S-tr|K|^3/q)/S. Under a fixed original norm cap,
  `Phi(K)>=kappa S/2-O(N^(3/2)sqrt(gamma)+N^(5/4))`.
  Its finite bound retains the actual weighted mask loss. It is vacuous
  for gamma>=1/4; do not reuse it as a complementary large-gap gain.
  The conclusion `u=c_D/n>=kappa-o(1)` additionally needs the ACTUAL
  active pure-cross conditions p=q_A=0 and c=Phi(K), with N=2n.
- `NOTE_2026-09-05_SMALL_GAP_PURE_CROSS_UPPER.md`
  (SHA-256 `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`)
  evaluates the actual
  squared-singular-value measure of W_D, including zeros. The exact
  sign-square identity and compatibility give
  `m=integral y dnu(y)=u^2/f_n^2+o(1)`, where c=f_n n^(3/2).
  For standard centered signs, at fixed t=3/5, concavity of A and B's convexity bound
  EVERY actual measure by the algebraic expression (3.5); no Dirac law
  or common extremizing measure is assumed. For
  f_n tending to sqrt(2), gamma tending to zero and the active conditions
  above, the resulting cell upper is
  `limsup E max X_z/(2n^(3/2))<=17677/25000<1/sqrt(2)`.
  The monotonicity argument covers all u>=kappa up to vanishing errors.

Both proofs and their complete review provenance are in the manifest.
Exactly eleven new fixed Fraction comparisons passed one soulkiller run;
result SHA-256
`0ea064435322e698b8e33a4d9bce8ab29156e3cfe013c9885f1f35e205156e41`.
The squared strict margin is `23671/625000000`. The earlier pi enclosure
was reused, not recomputed. The replay artifact is
`python3 evidence/original_mo_small_gap_pure_cross_rational_certificate.py`;
do not rerun unchanged arithmetic merely for another receipt.

This evaluates a previously unevaluated ACTUAL-law diagnostic face, not
the supremum over all coupled original/weighted cells. The general
formula retains f_n; the desired bound is still `F<=2sqrt(2)Phi(A)`.
If Phi(A)/n^(3/2) is below 1/2, the f=sqrt(2) result alone does not meet
that smaller target. The next implications are a bound at the actual
smaller normalization, control of nonzero original internal energies,
and a genuinely complementary positive-gap argument. A norm cap or
optimizer label does not supply small canonical gap. Original all-orders
convergence remains OPEN; none of these route-specific premises is a
necessary condition imposed on every possible convergence proof.

The latest quantified compatibility package is
`evidence/original_mo_gap_compatibility_milestone.json`:

- `NOTE_2026-09-05_FULL_SDP_GAP_WEIGHTED_COMPATIBILITY.md`
  (SHA-256 `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`)
  uses the canonical primal of the LITERAL complete K and ANY attained
  trace-optimal same-diagonal majorizer D. With S=tr D, q=N-1 and
  g=S-tr|K|^3/q, weighted residual squares are at most 4qg. The exact
  squared Frobenius norm of the inverse-weighted commutator is
  `2(S tr(D^(-1))-N^2)`, giving
  `delta<=4Sg/(qN^2)`. Uniform cube rescaling then proves
  `Phi(K-(S/N)T)<=S sqrt(delta)` for T=D^(-1/2)KD^(-1/2).
  The ORIGINAL source energy errors are at most 2N sqrt(delta), and
  the cross error at most N sqrt(delta), after scaling by S/N.
  No maximum-diagonal bound, nonsingular K, or unique optimum is needed.
  Its original-zero-source corollary compares two individually PSD
  Gaussian fields at O(N^(3/2)delta^(1/4)) cost. The pure-cross field
  retains actual W_D,c_D; its width is NOT evaluated by this corollary.
- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_METRIC_STABILITY.md`
  (SHA-256 `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`)
  applies to ANY actual majorizer D, not only an optimum. Constant field
  diagonal and a nuclear-norm congruence estimate control both exact
  resolvent traces and their cancellation. For 0<=delta<=1,
  `|B_D-B_flat|<=3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|eta|)`.
  Its (5.1) transfers the all-ACTUAL-cell upper on each fixed compact
  eta window, with the original drift and old bin/selection errors.
  B_flat is a NUMERICAL reference, retaining actual PSD M_theta and
  actual contraction L_D. The representative's weighted c_D need not
  be constant through a bin, so this is not an exact scalar-I shell
  constraint throughout that bin. Every representative is chosen within
  its FINAL refined original/weighted cell.

Both paths above are under `evidence/`. Complete author, root and
independent reviews are recorded in the manifest. These two analytic
proofs needed no mathematical execution, signing census or optimization;
the previous arithmetic certificates were reused without rerunning them.
No source signing, src module, test or global predicate was changed.

For bounded `S/N^(3/2)`, a relative canonical gap g/S tending to zero
gives delta tending to zero. This is a CONDITIONAL actual regime, not a
property established for every exact original or conditional minimizer.
Do not replace K by a purported contraction K/(S/N), remove its rare
diagonal outliers without accounting for them, or invoke an indefinite
scalar tensor covariance. The metric window must be fixed before its
asymptotic limit; the error is not uniform at |eta|=1.

The live next implications are an evaluated actual weighted trace upper
in the small-gap range and a correctly normalized full-K ORIGINAL-norm
argument for the complementary range. Rectangular beta bounds can lose
an essential factor through beta(K)<=4Phi(K); they cannot silently be
read as a quadratic-norm rounding gain. Neither these route-specific
targets nor the new compatibility bounds prove original convergence.

The latest same-source package is
`evidence/original_mo_weighted_covariance_milestone.json`.
Its four analytic results are:

- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_CROSS_COVARIANCE.md`
  (SHA-256 `0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`)
  chooses D with D+-K positive for the LITERAL complete block signing.
  The exact cross correlation is `R_D=I+Q(A tensor A-S_B+I)Q`,
  `q_ij=1/sqrt(d_i d_(n+j))`; its operator norm is below three.
  A norm cap supplies `tr D=O(N^(3/2))` and local correlations O(1/N)
  without trimming or changing A,-A,B. The complete weighted Hermite
  decomposition gives an O(n) retained-profile cost; the resulting
  conditional ORIGINAL-norm Gaussian floor has O(n^(16/11)) error.
  The separated even series is used only when epsilon<=1/2; the
  remaining bounded orders are treated separately. Dropping only the
  independent Gaussian variance padding costs O(n).
- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_WEIGHTED_SHELL_UPPER.md`
  (SHA-256 `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`)
  proves exact positive weighted linear fields. Actual representatives
  and width-1/n bins handle real weighted energies with O(sqrt(n))
  comparison error; selecting all original/weighted cells costs
  O(n sqrt(log n)). The same D majorizes the ORIGINAL cross form H_B,
  so D-eta H_B has exact shell radius `tr D-2eta c` using unweighted c.
  Equation (4.6) is the full all-actual-cell upper, with exact two-trace
  field width (4.5). Its weighted feedback c_D and internal p_D,q_D
  remain distinct from the original drift `(p-q)/2+s c`.
- `NOTE_2026-09-05_SCALAR_TEMPLATE_CUBIC_ALIGNMENT.md`
  (SHA-256 `60037f67234fbca8c17ee90bf52c7f4346b24e5f18eb5f2c922ebbd2d9382c2a`)
  proves `j_3>=-1+2s^2/mu` for ACTUAL matched optimal frames.
  For finite `tau(C)=p q`, this strengthens the Gamma UPPER-certificate
  exclusion to `q>=12/5 => Gamma(C)>283/200>sqrt(2)`.
  It crosses q=1+sqrt(2), the particular weak-Dirac diagnostic barrier.
  Section 6 separately stipulates leading energy sqrt(2), derives its
  variable-u expression, and takes n to infinity at fixed eta before
  eta tends to one. It does not infer actual Boolean saturation from
  a Gamma cap or prove arbitrary nonsymmetric attainability.
- `NOTE_2026-09-05_TENSOR_DEFLATION_FIXED_CAP_RATE.md`
  (SHA-256 `22febfa722afb3e18878f23f8e140895da90a3eb41fe0179356b08232d44f27a`)
  constructs actual complete sources under each FIXED norm cap C>1/2.
  Their tensor positive-part Gaussian repair has a uniform-rate lower
  of order n^(3/2)/sqrt(K), for both tensor signs and for the symmetric
  zero-diagonal quadratic norm. The amplitude is fixed before K, then
  n grows. No assertion is made for C=1/2, exact source minimizers,
  adaptively coupled slack and K, or the full A tensor A-S_B+I repair.

All four note paths are under `evidence/`. Each complete proof has an
independent review, with author/root collaboration disclosed in the
milestone manifest. The cubic note adds exactly five rational comparisons,
run once on soulkiller; the previous 28 clipping comparisons are reused.
The replay artifact is
`python3 evidence/original_mo_scalar_template_cubic_rational_certificate.py`.
Do not rerun unchanged arithmetic for another receipt. There were no
signing censuses, numerical optimizations, src-module or predicate changes.

The live unresolved implication is a sharp upper evaluation of the
weighted-shell note's (4.6) on ACTUAL coupled cells, using source and
conditional optimality. In particular its (5.3)-(5.6) retain the explicit
weighted/unweighted Delta_B and internal Delta_A discrepancies; setting
them to zero is not justified by the trace cap. The conditional Gaussian
floor does not supply this upper evaluation. The new covariance avoids
scalar attenuation and the generic deflation loss, but original all-orders
convergence and the proposed sharp dyadic inequality remain OPEN.

The new actual-sign realization package is
`evidence/original_mo_hadamard_template_milestone.json`. Its three complete
proofs passed independent reviews:

- `NOTE_2026-09-05_HADAMARD_SPARSE_FLIP_TEMPLATE.md`
  (SHA-256 `0d2355f94734b4c1e950c1e05c6df75df38b5ce181ba7fce550a4245e11328ed`)
  constructs actual cross sign matrices with flat singular bulk, finite
  outliers and asymptotically scalar-optimal SDP duals. The finite-template
  completion Gamma is an UPPER on their Boolean norm, not an attained value.
- `NOTE_2026-09-05_SCALAR_TEMPLATE_GAMMA_BOUND.md`
  (SHA-256 `bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`)
  uses matched optimal frames and exact quadratic Hermite cancellation
  to prove Gamma>283/200>sqrt(2) when q>=5/2 AND tau(C)=p q for
  the finite template. This excludes that upper certificate, not the
  actual matrices' Boolean cap.
- `NOTE_2026-09-05_HADAMARD_BOOLEAN_ALGEBRA_LOWER.md`
  (SHA-256 `68ce3f2f2a8fa2280208a9e145f508b6c2b2520d81e83185f579aaac89838a5d`)
  modifies the background to fix an entire Boolean algebra. Dense exact
  block-mean grids give a genuine actual lower Lambda_I. For PSD templates
  Lambda_I=Gamma and the actual normalized Boolean norm converges to it.
  With a symmetric POSITIVE top frame, actual liminf>1.524049912 for
  q>=5/2. Bipartite scalar SDP optimality alone is not that hypothesis;
  symmetric dilation changes the template and the actual sign family.

All three note paths above are under `evidence/`. The published rectangular
Bernstein theorem is stated and applied explicitly. One bounded soulkiller
run verified 28 exact rational comparisons for the clipping constants;
the later actual lower reuses that result without another run. No matrices
were sampled or enumerated. The live missing implication remains actual
low-norm nonsymmetric attainability and source/conditional-shell compatibility,
or a general upper using those constraints. These constructions are not
conditional optimizers and do not establish original convergence.

The preceding evaluated moment frontier is
`evidence/NOTE_2026-09-05_SCALAR_MOMENT_FEEDBACK_DIAGNOSTIC.md`
(SHA-256 `cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`).
At the old scalar endpoint, strongest source feedback plus rank-four
repaired positivity gives normalized squared upper `<9/20<1/2`.
The explicit fixed-metric repair-trace bounds require bounded scalar q,
not a conference-scale cap on A. Weak feedback still leaves a uniformly
positive gap for a formal Dirac moment law, even with the entire mixed
rounding family and its exact fourth-moment refinement. This is only a
counterexample to sufficiency of that MOMENT RELAXATION. The law does not
supply the top singular value and optimal Gram of an actual scalar dual,
an actual complete signing, or compatible source and joint-shell data.
Do not turn failure of an upper certificate into an actual width lower.

The literal old Krivine endpoint is already excluded for actual matrices
by [Braverman--Makarychev--Makarychev--Naor, Theorem 1.1](https://web.math.princeton.edu/~naor/homepage%20files/GroKri.pdf).
The explicitly defined `K_G=pi/(2 asinh(1))` in these diagnostic notes is
the elementary Krivine bound, not the exact real Grothendieck constant.
The retained-moment insufficiency theorem is still valid; it is not a
claim of insufficiency after imposing every published constraint. Actual
realizability work must move off that excluded exact endpoint and examine
admissible nearby ratios or a uniform range, not repeat its realization.

The new actual-matrix constraints supporting this evaluation are:

- `NOTE_2026-09-05_CROSS_SINGULAR_MOMENT_ROUNDING.md`: filtered cubic
  rounding with an exact quartic error and a cap-free clipped variant.
- `NOTE_2026-09-05_CROSS_SDP_COMPLEMENTARITY.md`: every optimal diagonal
  is block-balanced; its canonical-primal gap controls weighted residuals.
  Zero gap is equivalent to equal NONZERO singular values and forces a
  scalar optimum. A scalar optimum alone is not the converse.
- `NOTE_2026-09-05_ORIGINAL_PHASE_SPECTRAL_MOMENT.md` and
  `NOTE_2026-09-05_ORIGINAL_PHASE_MOMENT_BOOTSTRAP.md`: actual positive
  and negative source phases yield cubic and nuclear constraints. The
  bootstrap retains `2(n-1)+osc diag(A|A|)` with `O_C(n^(5/4))` ORIGINAL
  norm error under only `Phi(A)<=C n^(3/2)`, without replacing A.
- `NOTE_2026-09-05_CROSS_TENSOR_MIXTURE_SIGN_DEFECT.md`: actual canonical
  negative sign mass is at most half its SDP gap, or one quarter for a
  scalar optimal dual. Tensor mixing yields an explicit limiting curve
  strictly stronger than the cubic-only constraint below beta/tau=2/pi.
  Its `O_(C,t)(n^(5/4))` error requires t fixed before the n limit.

All paths in that list are under `evidence/`. Each complete proof has an
independent review; exact hashes, source aliases and reviewer roles are in
`evidence/original_mo_spectral_rounding_milestone.json`. No numerical
mathematics or signing census was used. The live next step is to exploit
actual optimal-Gram/coordinate and source-cross compatibility, or obtain
a sharper upper using it; optimizing the same weak-feedback moment-only
functional again does not address its now-proved insufficiency. The
original all-orders convergence question remains OPEN.

The underlying Gaussian upper is
`evidence/NOTE_2026-09-05_BOOLEAN_ELLIPSOID_SHELL_UPPER.md`
(SHA-256 `ede1b62a26a636179d918ba84a48d122ab013c38175bdb9cd164bcfd8bfeb9aa`).
For actual PSD C, positive P and nonnegative diagonal E<=P on a Boolean
shell z^T P z=q, it proves the exact completion-square remainder (3),
not merely a uniform sphere-to-cube multiplier. The diagonal-metric
limit is the exact cube width. Equations (12)-(17) completely evaluate
the weaker diagonal-affine specialization, including negative parameters
and singular-metric limiting infima. The stronger two-trace expression
(18) retains a trace with BOTH signs; do not substitute separate upper
bounds into it without checking the combined expression. The actual
cushioned scalar diagnostic improves the old bound but still exceeds
the desired leading constant. An indefinite reference is not a covariance.

The accompanying fixed-internal-block tool is
`evidence/NOTE_2026-09-05_CROSS_ONLY_OPERATOR_REGULARIZATION.md`
(SHA-256 `27d9ab77768e8b7afa2d48d041cf3fe6bf3b66e8b16e481ca12abcf906a28d4f`).
Its exact loss is the two selected A-cut norms plus `2n sqrt(s)`, where
s counts exceptional cross rows and columns. For `||A||<=K_A sqrt(n)`,
the normalized loss is at most `(2+sqrt(2)K_A)sqrt(Lambda C/K)` and
the new cross operator cap is `(K+8)sqrt(n)`. The source A,-A is unchanged.
Near-source selection and subsequent cross regularization have separate
slacks and cross cap `O(epsilon_A^(-4)epsilon_B^(-2))`. Do not transfer
exact optimizer properties to either regularized object or overlook
competition between this cap loss and the evaluated Gaussian gain.

Both complete proofs have independent reviews. Their provenance is in
`evidence/original_mo_boolean_ellipsoid_milestone.json`; the two elementary
scalar-check offloading exceptions are recorded there, not relabelled
as remote checks. The actual all-shell leading comparison remains open.

The current joint-shell package has three independently reviewed proofs:

- `evidence/NOTE_2026-09-05_CONDITIONAL_CROSS_JOINT_SHELL_UPPER.md`
  (SHA-256 `64d68bb2feaa59a8049d6bcc42f3ab94c845249c3088fa618916522412d0a68a`)
  proves the exchange-preserving cushioned field upper and the separate
  masked-cross conditional floor, with raw error `O(n^(16/11))`.
- `evidence/NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`
  (SHA-256 `e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`)
  uses the intrinsic cross operator `H=A tensor A-S_B+I`,
  `mu=max(2,||H||)`, and `R_mu=I+H/mu`. The entire threshold covariance
  correction is controlled; `||K||^2<=8 Phi(K)` and an elementary
  conditional norm cap give the same `O(n^(16/11))` Gaussian floor for
  actual cross optimizers over ANY original exact minimizer A.
- `evidence/NOTE_2026-09-05_INTRINSIC_CROSS_JOINT_SHELL_REPAIR.md`
  (SHA-256 `1dcd9b1e76b00887e406e505113c854b80f0661bb3bd69283f6486fb59fa2d53`)
  repairs the intrinsic linear fields by a rank-four PSD correction.
  It proves the genuine upper and retains the leading joint-mismatch
  formula with error `O(n^(5/4))`, uniformly even in vanishing-noise tails.

Use actual attainable `(x^TAx,y^TAy,x^TB_*y)` shells and conditional
optimality to sharpen/evaluate this upper. Do not drop the independent
cushion or mixed exchange term, declare the unrepaired intrinsic field
PSD, or replace conditional optimality by full-order optimality. The
needed leading comparison `F_A^*<=2sqrt(2) Phi(A)+o(n^(3/2))` is not proved.
Even proving a little-o dyadic inequality would not by itself settle
the original all-orders problem. These are optional analytic tools,
not newly mandatory proof architecture. Provenance and backup coverage:
`evidence/original_mo_conditional_joint_shell_milestone.json`.

The new whole-source variational tool is
`evidence/NOTE_2026-09-05_WHOLE_EDGE_SOURCE_PRESERVING_GAUSSIAN_REDUCTION.md`
(SHA-256 `6b22fb3ab1cc878b08fe79b5b57e0e661eaaa792dfc67f850d35db9f1b68bead`).
It uses all UNORDERED original edges and the normalized symmetric
compression `R=(L^2 I-T)/(L^2-1)`, with `T(X)=KXK` compressed to that
edge basis. For n>=3, `0<=R<=3I` for EVERY complete source; at n=2 the
compression is zero and the theorem uses an explicit independent fallback.
The entire even-Hermite correction is handled by a rank-one term and a
four-cycle operator bound. The resulting expected ORIGINAL whole-order
norm error is absolute `O(n^(16/11))`, uniform over deterministic h.

Its full symmetric Gaussian lift removes the diagonal at expected cost
at most `sqrt(n/pi)` and retains the exact augmented replica matrix
`Gamma=<sigma xx^T>`, not a positive-semidefinite substitute. Equations
(25)-(29) prove a negative-current-energy-square variance upper and an
integrated constraint for actual ORIGINAL norm minimizers, with error
`O(n^(16/11))` at `beta=n^(-5/11)`. Both independent complete reads passed.
These are SAME-order constraints. A valid mapping to the required order
upper remains unproved; do not reverse the lower inequality, assume
opposite diagonal blocks for a full optimizer, or treat a shifted
disorder-dependent posterior as another minimizer.

The threshold-optimized extension is now
`evidence/NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md`
(SHA-256 `a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`)
and `evidence/NOTE_2026-09-05_SHIFTED_THRESHOLD_COVARIANCE_REDUCTION.md`
(SHA-256 `74457650912a515eaf6a209b184e5c1404a13fc48c68464068871ebd61236680`).
The mean-preserving OU proof is uniform in all real thresholds and keeps
every Hermite order and the actual posterior. The even covariance term
is not discarded: its PSD low-rank part has actual Gaussian Boolean-norm
cost `O(Phi(A) sqrt(log(2n)/n))`; the remaining operator error is `O(1/n)`.
Thus ANY exact original minimizer A satisfies the proved one-sided bound
`m_(2n)<=inf_h E Phi([[A,Z_h],[Z_h^T,-A]])+D n^(16/11)`,
where `Z_h=s_h A+2 phi(h)G+sqrt(1-s_h^2-4 phi(h)^2)W`,
`s_h=2 Phi_Gauss(h)-1`, and G has the universal midpoint covariance.
The threshold is fixed before drawing disorder, not chosen adaptively.
Both complete proofs have independent full-read PASS receipts.
The remaining target is an actual evaluated Gaussian upper bound, not
another covariance identity or an unsupported derivative sign.

The underlying zero-threshold rounding reduction is
`evidence/NOTE_2026-09-05_UNIVERSAL_SPECTRAL_MIDPOINT_GAUSSIAN_REDUCTION.md`
(SHA-256 `1fc6f5bbb69038b6ac4ed845d26e0724a0ceb0b5a9d96d01b4554a8e37e6f968`).
For EVERY complete source A with extreme eigenvalues `a,-b`, freely choose
`alpha=(a-b)/2`, `mu=(a^2+b^2)/2`. The exact covariance has operator
norm `(a+b)^2/(a^2+b^2)<=2`, and its arcsine-linearization remainder
is bounded by `(1-2/pi)(2/(n-1)^2+1/(n-1))`. The generic quenched
theorem therefore gives an absolute `D n^(16/11)` expected ORIGINAL
paired-norm error for all sources and all fixed internal energies.
No source regularization is needed: for ANY exact original minimizer A,
`m_(2n)<=E Phi([[A,Z],[Z^T,-A]])+D n^(16/11)`.
The full alpha domain, operator-optimal midpoint and actual-pressure
derivative are proved; operator optimality is not pressure optimality.

The immediate Gaussian upper-bound tools are
`evidence/NOTE_2026-09-05_GAUSSIAN_ENERGY_SHELL_UPPER.md`
(SHA-256 `8bd3507b722d13077cdb47e8eaa47024b8e95144900226ae4e38272795c5c728`)
and `evidence/NOTE_2026-09-05_ONE_PHASE_GAUSSIAN_VARIANCE_UPPER.md`
(SHA-256 `1646f57b060db7fdaf15c2cc8a8766806d2f00297c6749e236d8e814e467bae0`).
The first retains exact source-energy shells and gives a quantitative
one-block width deficit; its central-shell two-field comparison is still
too weak. The second removes the augmented phase at subleading expected
cost and retains the coupled posterior in the actual variance derivative.
The live target is an evaluated Gaussian order upper bound on actual
original minima. Neither a favorable integral sign nor a sufficient
cross-order inequality has been proved. An unspecified little-o dyadic
inequality alone would still not establish convergence.

The following Gibbs-generated-law results remain valid separately. They
are not prerequisites for the freely chosen universal midpoint law.

The new comparison is
`evidence/NOTE_2026-09-05_CORRELATED_SIGN_GAUSSIAN_FREE_ENERGY.md`
(SHA-256 `2e6537d0b1e2c4d8a72cc920e3fee50600d82be32417ba77c733aaedabc141c7`).
Root and both independent complete proof reads passed. For any bounded
Gaussian covariance operator, `n^2` correlated Gaussian signs have the
same quenched critical pressure as their matched Gaussian, up to
`O(n^(17/18))`. The entire posterior and singular endpoint are retained.
The key third central moment is contracted BEFORE taking absolute
values; Gaussian Holder and sign smoothing control the nonsmooth limit.
Covariance matching alone was not the proof.

For canonical sources with `||A||op<=K sqrt(n)`,
`evidence/NOTE_2026-09-05_CANONICAL_COVARIANCE_GAUSSIAN_LINEARIZATION.md`
(SHA-256 `44188dde396587f1d148e01857365b44d1bddbe83d81dc5b085ccee0cdff9854`)
has an exact disjoint-support tensor identity. Its Gaussian covariance
remainder has operator norm `O_K(1/n)`, hence pressure cost `O_c,K(1)`.
The resulting simpler covariance is `(2/pi)Sigma+(1-2/pi)I`.
The remaining live implication is an UPPER comparison of this actual
Gaussian paired pressure against the appropriate optimized smaller-order
endpoint. No such bound has been proved. The new reduction also applies
to sufficiently slowly growing caps (fixed-c error
`O(K_n^4 n^(17/18))`), hence to leading ORIGINAL norm near-minimizers
provided by same-order regularization. It does not assert that every
unregularized or quartic-penalized minimizer has bounded operator norm.

The direct original-norm consequence is
`evidence/NOTE_2026-09-05_EXPECTED_PAIRED_NORM_GAUSSIAN_EQUIVALENCE.md`
(SHA-256 `bff778718c0f357598c035edba4598f2ed67b1c49359c668958afe1c39207df3`).
It compares EXPECTED maximum absolute paired energies with normalized
error `O_K(n^(-1/22))`, with constant `O(1+K^4)`. The auxiliary choices
`c=n^(1/22)` and `epsilon=n^(-1/11)` use the full explicit bounds;
the source covariance-generating temperature is held fixed separately.
Regularizing actual original-norm minima at threshold `n^(1/99)` gives
both objective loss and Gaussian-reduction error `O(n^(-1/198))`.
The next target is therefore a Gaussian doubled-norm upper comparison
on these genuine near-minimizers, not an identification of pressure
surrogates or pointwise closeness of individual cross outcomes.

Work directly with the global optima `m_n`. A genuine advance must compare
orders or otherwise control their normalized oscillation. A construction
checked only on a selected low-norm example need not extend to actual
minimizers. A theorem for every signing under a proved norm cap does apply
to minimizers; neither distinction may be hidden in a hypothesis.

The new same-order reduction is
`evidence/NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md`.
Its SHA-256 is `8a52b7e4f171cc2089a00a6fd288e041d52605f820e49ace419ddd5fe850bec8`;
root and both independent complete proof reads passed.
For every complete signing with `Phi(A)<=C N^(3/2)`, a diagonal SDP
majorizer, vertex trimming and one jointly good random recompletion
give a complete signing at the SAME order with
`||A'||_op<=(K+8)sqrt(N)` and normalized norm increase at most
`2sqrt(Gamma C/K)`, where `Gamma=4pi/log(1+sqrt(2))`.
This applies directly to ORIGINAL norm minima. Bounded-operator
constrained minima therefore approximate the actual normalized minimum
uniformly as the bound increases; an arbitrarily slowly diverging bound
admits leading norm near-minimizers. This is one-sided objective control,
not small `Phi(A'-A)`. The missing implication is a useful order
comparison in this controlled class. Neither bounded operator norm nor
typical restriction has been proved to supply one.

The next regularized comparison is now explicit in
`evidence/NOTE_2026-09-05_QUARTIC_PENALIZED_PROFILE_IDENTITY.md`
(SHA-256 `ad393709abb35ed760986b102e1b86ab4d23c80261efec04f35d03104c821013`).
For the actual minimum of `F_c(M)+lambda tr(M^4)` on the balanced
profile, all edge flips are admissible. Every row obeys
`E_i+8lambda(M^4)_ii+8lambda sum_j M_ij^4<=c^2 d+16lambda d^2`.
Thus the fourth diagonal moments are uniformly bounded. Tensor rounding
and the Boolean norm cap give `sum|Gamma_e|=O_c(N^(3/2))`, uniformly
for `0<lambda<=1`; the diagonal SDP majorizer also bounds `tr|M|^5`.
The exact identity is
`G_N(1)-G_N(0)=c^2/4+lambda(5-9/N)-integral J_N^lambda+O_c(sqrt(N))`.
The error is uniform over `0<lambda<=1`. Each actual penalized flip
gap is nonnegative, has bounded row sums and is `O_(c,lambda)(N^(-1/2))`,
but the mixed weighted gap integral still has no proved favorable sign.
Do not substitute a permutation average for the selected envelope
derivative. Its zero-cross endpoint is twice the penalized HALF-PRODUCT
minimum and is at most twice the penalized symmetric minimum; equality
with the latter is not needed or claimed.

The pressure approximation uses one and the same recompletion in
`evidence/NOTE_2026-09-05_SPECTRAL_REGULARIZATION_PRESSURE_PROFILES.md`
(SHA-256 `2f9f63f603fcae42a952fbae53a2301eaa6b95bbe7bac2e35bcab8997d28d7d7`).
It controls both actual phases for all c in a prescribed compact interval,
with normalized cost `Gamma C c^2/(2K)+O(log(N)/N)`, while retaining
the operator and norm bounds. At fixed c, quartic penalized minima
therefore approximate original symmetric minima within
`O_c(lambda^(1/3))`. Vanishing regularized oscillation would suffice;
it has not been established. The whole-row and multi-edge variational
constraints are now proved in `NOTE_2026-09-05_QUARTIC_PROFILE_ROW_RESET.md`.
The weighted signed force kernel is controlled in
`NOTE_2026-09-05_QUARTIC_FORCE_KERNEL_BOUNDS.md`; the weighted row-tilt
identity does not assert a fourth moment at its endpoint. Independent
coefficient refills in `NOTE_2026-09-05_QUENCHED_BIASED_COEFFICIENT_REFILL.md`
retain the full quenched posterior and exact quartic correction, with
`O(sqrt(N))` replacement error even over all edges. The separate actual
canonical cross law adds at most `(41+88C^2)lambda t n` to the paired
quartic penalty. None of these same-order finite variations supplies
the missing Gaussian endpoint inequality; no new signing census is needed.

The actual-Gibbs structural proofs are
`evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md`
(SHA-256 `46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`),
`evidence/NOTE_2026-09-05_EXACT_HALFPRODUCT_SUBCRITICAL_SPECTRAL.md`
(SHA-256 `10dfe02b63aa3c4aa987ce48d4a3e660e90509b43e6a50a1a002ba9ecc1cc522`),
and its strengthening
`evidence/NOTE_2026-09-05_HALFPRODUCT_NEARMINIMIZER_STRUCTURE.md`
(SHA-256 `dccc256d3b7119c666102e54cffe3a2026d31edc1bcd0c4366a15ce92c762f0f`).
A Boolean energy cap gives a positive extensive response to any field
with a positive density of moderate nonzero coordinates, even with
unbounded outside coordinates; the ACTUAL Gaussian posterior is retained.
For EVERY leading half-product near-minimizer at fixed `c/sqrt(N)`,
approximate optimality, eigenvector truncation and sparse pinning prove
`||A||_op=o(N^(3/4))`. Deleting ANY `o(N)` vertices also changes the
full Boolean energy uniformly by `o(N^(3/2))`. Every gap is retained.
These are not exact-minimizer-only properties. Half-product pressure
approaches half the energy WIDTH, not necessarily the absolute norm;
no original-norm transfer or fixed-fraction comparison is inferred.
Complete root and independent reads passed; see
`evidence/original_mo_spectral_structural_root_review.md`.

The singular full-strength criterion is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_SPECTRAL_DEFICIT.md`:
`V_r=tr[-rI-H/mu]_+=o(N)` implies the stated quenched mean failure
and vanishing success probability, not an exponential original-law tail.
The construction-cap example is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_CONSTRUCTION_CAP_HOSTS.md`.
The new subcritical base strengthens this to actual leading HALF-PRODUCT
near-minimizers in
`evidence/NOTE_2026-09-05_FULL_STRENGTH_HALFPRODUCT_NEARMINIMIZERS.md`
(SHA-256 `ad83095163cf8e969e542a6626382dabaa5adb5e2ffce4bfffea274c813b53e4`).
One reused sparse module costs `o(N)` pressure and gives exactly
`V_r=2(1-r)` eventually. This is not an exact-minimum, original-norm
near-minimum or unrestricted selected-outcome exclusion. The separate
`evidence/NOTE_2026-09-05_NUCLEAR_SPECTRAL_BUDGET.md` gives
`Phi(A)>=N^2(N-1)/(pi tr|A|)` and a linear nuclear effective rank under
the relevant objective caps, not spectral flatness. No new mathematical
computation was used in these results.

The actual covariance corollary
`evidence/NOTE_2026-09-05_ACTUAL_GIBBS_COVARIANCE_FLOOR.md` gives a
positive diagonal component of trace at least `chi_c N` in each actual
phase covariance. Thus rank-`o(N)` truncations cannot have `o(N)` tails,
and every integral full cross block has `qbar>=chi_c^2 N^2` under a
fixed Boolean cap. This retires the conditional sublinear-rank/tail
escape at fixed c, not the radial upper comparison, which can be
quadratic as well. The constant is not uniform at zero temperature.

The new unconditional selected-restriction theorem is
`evidence/NOTE_2026-09-05_SELECTED_HALF_RESTRICTION.md`
(SHA-256 `c8a9aa0b8c44fb14f444955fbe3eec8cba8e7f19c01fb8eeb2596418d3416c02`).
Its complete independent root read passed. An explicit odd cycle of
disjoint subsets gives a half-norm restriction at order `2n+1`; a
complementary-phase exchange argument gives boundary error `(n-1)/2`
at order `2n`. Applied to actual global norm minima, these prove
`m_(2n+1)>=2m_n` and `m_(2n)>=2m_n-(n-1)`.
This improves the old fixed-partition/monotonicity estimate, but the
normalized comparison still has factor `sqrt(2)`, not `1+o(1)`.
The missing leading factor is the issue, not the linear boundary error.
Do not confuse selected half-norm restrictions with typical tiny-n
restrictions at the sharper source-normalized threshold.

The exact full-strength identity is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_BOUNDARY_LIKELIHOOD.md`
(SHA-256 `8703433f6118f00dd589d711e9541f558489caa3d13059f8e71405333401fdb2`).
Root and independent complete reads passed. The derivative of the actual
planted log likelihood is a weighted sum of mixed finite differences
under PAIR-DEPENDENT Gaussian boundary laws, and its integral is valid
through singular `rho=1`. An exact actual order-three minimizer refutes
the coordinatewise sign premise; its counterexample context disappears
from the endpoint support. This does not refute the weighted average or
the full-strength finite-step comparison. Retain the boundary support,
full posterior and possible cancellations; the negative prior trace
alone does not control this integral. No computation was used.
The complete root review for both notes is
`evidence/original_mo_boundary_restriction_root_review.md`.

The fixed-strength strengthening is
`evidence/NOTE_2026-09-05_FIXED_WEAK_GAUSSIAN_CUSHION.md`
(SHA-256 `5df7258c4cf99deac09eaeb4a206e1270ffa7add1e49e176b70a4a232eb54d12`).
Root and independent complete reads passed. For ANY latent Gaussian
correlation matrix `S`, conditional independent-sign replacement and
convexity give the actual quenched floor
`[c sqrt(2t) K0-2log2-c^2 t arcsin(rho)/(2pi)]n-o(n)` for covariance
`(1-rho)I+rho S`. Thus a NONEMPTY interval of fixed positive strengths
is excluded at suitable fixed `c,t`, even though its information is
not `o(n)`. Precisely, the gap
`Delta_rho=c(sqrt(2t)K0-1)-2log2-c^2 t arcsin(rho)/(2pi)` must be positive.
The actual centered latent law has `||S||_op<=4n-3` for EVERY generating
host, so the heat-martingale and conditional bounded-difference bounds
make success exponentially rare. Even `exp(o(n))` proposals with these
marginals, including legitimate pre-draw host mixtures, fail; proposals
need not be independent and the internal host may be selected afterward.
This is not a mixture entropy lower bound or an exclusion of `rho=1`.
The complete independent review is
`evidence/original_mo_fixed_weak_gaussian_cushion_exact_review.md`.
No new computation is used. Do not repeat weak fixed-strength sampling
at these parameters or extend this result outside its explicit gap.

The preceding dependent-rounding information theorem is
`evidence/NOTE_2026-09-05_GAUSSIAN_SIGN_INFORMATION_SCALE.md`
(SHA-256 `5846e981204f03230bbfd415443824d1a320840d56b6163267e37ee1b8e5e566`).
Complete proof reads passed. Every sign law satisfies
`D(Q||iid)>=||C-I||_F^2/(4||C||_op)`, with the SECOND-MOMENT matrix
`C=E bb^T`, not a centered covariance absent a mean-zero hypothesis.
For any Gaussian correlation matrix, including singular ones, arcsine
and the Schur product give
`D(sign N(0,Sigma)||iid)>=||Sigma-I||_F^2/(pi^2||Sigma||_op)`.
For the actual centered tensor `H`, `mu=-lambda_min(H)` and
`Sigma_rho=I+rho H/mu`, this implies `Omega(n)` discrete information
at every fixed `rho>0` on norm-capped hosts, including singular `rho=1`.
The proof uses `||A||_op^2<=16Phi(A)` and retains `mu` in the operator
denominator before combining the ratio. No Gaussian determinant upper
bound is used at the singular endpoint. Thus full-strength canonical
rounding is OUTSIDE the low-information exclusion, not proved successful.
On actual half-product minimizers, `rho=o(n^(-1/2))` instead gives
`o(n)` information and is excluded in mean and with substantial success
probability by the following quenched theorem. Strengths outside the new
cushion criterion, the actual Gram--Schmidt law and unrestricted selected
cross blocks remain open. A conditional-law entropy lower
bound must not be extended to arbitrary mixtures over hosts.
No computation is used; the independent general entropy audit is
`evidence/original_mo_entropy_covariance_review.md`.

The preceding iid all-orders cross-block result is
`evidence/NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md`
(SHA-256 `97e1aeb3ac25c2570072d9f0ebdb0c4387f739ed3c005ec7b43d30409dd7ade4`).
Root and independent full reads passed. An explicit Gaussian martingale
control in the sourced zero-temperature Parisi formula proves
`P_SK>=K0=4/(3sqrt(pi))>1/sqrt(2)`. A host-free pure-cross pressure lower
bound, Gaussian covariance interpolation, and direct fixed-temperature
Bernoulli replacement give `E F>=(c sqrt(2t) K0-2log2)n-o(n)`.
Against `2R_n<=cn+o(n)`, the gap is positive when
`Delta=c(sqrt(2t)K0-1)-2log2>0`. Bounded differences prove an exponentially
small iid probability of ANY good internal host at such a cross block,
so `exp(o(n))` proposals with iid matrix marginals cannot succeed even
when dependent across proposals. A successful law must have relative
entropy at least `(Delta^2/(c^2 t)+o(1))n` from iid signs.
Thus iid QUENCHED and `o(n)`-information selection are excluded at those
fixed parameters; arbitrary dependent selection is not. The exact
planted-channel identity retains the reverse relative entropy and the
full actual Gibbs prior. No numerical SK constant or new computation is
used. Do not rerun iid samples or confuse this theorem with the earlier
annealed Gaussian-sign obstruction.

The independently reviewed all-orders coefficient results are
`evidence/NOTE_2026-09-05_POSITIVE_CONE_TRUNCATION.md`
(SHA-256 `632adeb92932db37ba1ac218621eb3f7d1b8bd24e8461273abf74a379d79d304`),
`evidence/NOTE_2026-09-05_EXTENSIVE_COEFFICIENT_MOMENTS.md`
(SHA-256 `b07772332265dea635c59a7d293562feedb5c57cb7b66d7850f77c1ffbd4107e`),
and `evidence/NOTE_2026-09-05_POSITIVE_DEGREE_SELECTOR.md`
(SHA-256 `20dae4c37ece2f5c5808595c54941de1b10a241d03b63c4431c76dc373849875`).
The exact central-factorial coefficients are nonnegative. A cutoff at
`k<=K_N=o(N)` loses extensive pressure for fixed `c>pi log2`, even on
actual norm or symmetric-pressure minimizers; the latter require their
separately justified norm cap. In any fixed positive `k/N` band, the
coefficients are within dimension-uniform multiplicative constants of
`E|Q_A|^(2k)/(2^k(2k)!)` for EVERY complete signing. Convergence of the
separately optimized coefficient rates at unbounded fixed `k/N` values
would imply convergence of `alpha_N`, but that transport is still open.
One selected extensive degree per signing already captures log pressure
to `O(log N)`; do not misstate the cutoff theorem as excluding sparse
degree selection. Mixed minimax is legitimate with its quantifiers;
moving a pure minimum through the coefficient sum is not. No new
census or simulation is used by these proofs. The complete independent
review and correction record is `evidence/original_mo_coefficient_quenched_review.md`.

The new fixed-order analytic theorem is
`evidence/NOTE_2026-09-05_EXACT_OPTIMIZED_ORDER_SIX_PROFILE.md`
(SHA-256 `a1469b34118da1bf971c7d53ad0fb8c50525f588a42bf7c4f2dda9b132966fd4`).
Root and independent complete proof reads passed. For ALL `u>=v>=0`,
the minimum of `E cosh(u I+v C)` over complete order-six signings is
`cosh(v)*(3X^2+3Y^2+2Y-4)/4`, where `X=cosh(2u),Y=cosh(2v)`.
For `v>0`, minimizers are exactly `A^2=5I`; there are twelve after
first-row-positive switching normalization. The proof compares every
coefficient on `X=1+p+q,Y=1+q`, not a finite grid. Its success uses the
candidate's exactly quadratic polynomial: low-moment comparisons do NOT
control higher positive candidate coefficients at larger orders.
Along `u=c sqrt((2-t)/6),v=c sqrt(t/6)`, the optimized endpoints cross
exactly once for positive `c`, and
`f6(c,2/17)-f6(c,0)=((sqrt(17)-4)/sqrt(3))*c-log(2)+o(1)` as `c->infinity`.
Thus no temperature-uniform bounded interior excursion holds even at an
actual global optimum. This does not refute a fixed-`c` small-oh order
comparison or convergence. The finite-temperature maximum is not asserted
to occur at `2/17`; the exact left derivative at `t=1` is positive.
Do not rerun the same catalog/grid or treat fourth moments alone as an
all-orders extension. A new argument must control the actual higher
coefficients, selected finite-step pressure, or another order comparison.

The latest finite-step theorem is
`evidence/NOTE_2026-09-05_FINITE_STEP_ROUNDING_ANNEALING.md`
(SHA-256 `058cdd3e17972be45a664b21e720fafd44c744194c2e8f3bcb37e818d474ee0a`).
The full proof passed root and independent reviews. Its Gram--Schmidt
bound retains `log E_nu exp((gamma^2/2) v^T G^-1 v)`; replacing this
log moment generating function by the average quadratic form is invalid.
For ALL `G>0,diag G<=1`, the resulting proxy has floor `c^2 t n/4`.
A separate Gaussian entropy-tilt proof gives the ACTUAL canonical
Gaussian-sign annealed floor `c^2 t n/(2pi)-o(n)`, uniformly over hosts
and admissible centering, with no covariance operator-norm hypothesis.
Here `c,t` are fixed as `n` grows. Since `2R_n<=cn+o(n)`, the respective
annealed certificates cannot give the needed small-oh finite-step
comparison when `c>4/t` or `c>2pi/t`. The Gaussian claim concerns
`log E_B exp F_B`, not `E_B F_B` or `min_B F_B`. The actual Gram--Schmidt
law is not excluded. The preceding local second-moment theorem survives.
Do not optimize the same quadratic proxy or substitute annealing for
selected-outcome control again; shrinking steps and different laws remain
outside the stated obstruction. No finite sample proves this theorem.

The preceding integral construction is
`evidence/NOTE_2026-09-05_INTEGRAL_CROSS_BLOCK_COVARIANCE_ROUNDING.md`
(SHA-256 `c02bcc4d73ca58ba701b80a1fd73fa1c54f928effd5a62fe77daa9925c7d5c01`).
Root and two independent agents checked the complete proof. With actual
opposite-temperature covariances `U,V`, define
`qbar(B)=(tr(B^T U B V)+tr(B^T V B U))/2`, allowing all `n^2` entries
of `B` to be independent choices of signs, including its diagonal.
The complete-sign Gaussian construction proves
`min_B qbar(B)<=n^2-8(a_A')^2/(pi ||A||_op^2)`.
Its sharper form uses the exact negative spectral edge of the centered
tensor matrix; scalar entrywise arcsine is justified by disjoint entry
types, not matrix functional calculus. The host is fixed during rounding.
General Gram--Schmidt covariance rounding additionally gives integral
spectral-tail bounds and retains fixed coordinate squares through diagonal
shifts. The conference-form scalar-shift optimum is not a limitation
theorem for all diagonal shifts or all rounding methods.
The general comparison `min_B qbar(B)<=2a_A'/beta+o(n^2)` remains open;
the sufficient low-effective-rank case is not established for minimizers.
Even that endpoint bound would not control the integrated balanced path.
No new computation is used by these analytic theorems.

The preceding analytic method check is
`evidence/NOTE_2026-09-05_NEAR_MINIMIZER_OPPOSITE_PHASE_COUNTERFAMILY.md`
(SHA-256 `8130ca8c0af67d9976f71f086a79607d0b7b640b1e5c35ba6eb08d87e81324f7`).
Its complete proof passed root and independent review. For fixed `c>0`
and every sufficiently large `N`, paired modules can be planted into an
arbitrary old signing while changing every spin energy by at most
`O_c(N^(11/8))`. Conditional entropy costs only the new vertices;
Rademacher averaging retains the full Gibbs feedback. A simultaneous
thermal/operator event and exact even-module covariance decomposition give
`tr(A U A V)=Omega_c(N^(9/4))`. Choosing an old norm minimizer gives norm
excess `O_c(N^(11/8))`; separately, choosing a half-product minimizer at
the SAME raw `beta=c/sqrt(N)` gives pressure excess `O_c(N^(7/8))`.
Both are leading-order near-minimizers, not merely correct-scale hosts.
They need not be the same family and are not edge-local half-product
minima. Do not extend the counterexample to exact minimizers, or treat
this moment condition as necessary for convergence. No computation is
used in the all-orders proof.

The preceding analytic method check is
`evidence/NOTE_2026-09-05_FULL_ROW_CAVITY_COUNTEREXAMPLE.md`
(SHA-256 `54de76afacf34c7443ece9f5a34c42ef32d741b6fa381a7f5b9412675a1b331f`).
Root and independent full reads passed; no computation was used. A positive,
exchangeable, even arbitrary cavity has a strict minimizing sign row up to
global reversal, yet its actual row-tilted second and fourth moments grow
at least as `sqrt(d)` and `d^(3/2)` at fixed critical row scale. Every subset
replacement and the complete row-noise hierarchy holds. Thus these
inequalities alone cannot prove the desired bounded tilted moments.
No actual quadratic-host realization was supplied. Do not relabel this
as an Ising counterexample or make bounded moments mandatory for convergence.

An independently derived local endpoint calculation retains actual phases.
Let `A` minimize `a_A(beta)=(log Z_+(A)+log Z_-(A))/2`, where
`beta=c/sqrt(n)`, and let `U,V` be its two phase covariance matrices.
The earlier coherent admissible choice in the paired family `A,-A`,
with cross block `B=A+D` and a fair independent signing of its diagonal, gives
`f_(2n)'(0+)<=-beta a_A'(beta)/2+(beta^2/4)(tr(A U A V)+n)`.
The exact derivative minimizes over ALL active block pairs and cross
signings; pairing with the negative is only an admissible upper bound.
The integral construction above allows other cross blocks; this coherent
trace condition is not required for that enlarged choice.
The trace comparison that would bound this derivative above by `o(n)` is unproved.
Even such an endpoint derivative bound would not by itself control the
whole interpolation. The formula is analytic; the fixed-order numerical
check below neither proves nor refutes an asymptotic small-oh comparison.

The preceding two analytic results are
`evidence/NOTE_2026-09-05_ADAPTIVE_PERTURBATION_CORRELATIONS.md`
(SHA-256 `054063ac00e2fda45b676fc9a257cb901f43627e83cf31ce0e061b7c8816bb5f`)
and `evidence/NOTE_2026-09-05_OPTIMIZED_GAUSSIAN_SWITCH_MEASURE.md`
(SHA-256 `9c5090ddf4e1b43222716182ce5de3c51216ad87cb6f216d4bcd3ea70571fa0a`).
Both complete proofs passed root and independent review; they have no
new finite-check coverage. For arbitrary fixed additive edge noise `E`,
edge-local sign optimality gives `sum |Gamma_e|<=4k tanh(beta)+||E||_2^2`,
with rowwise and balanced-profile versions. The physical edge flip is
`-2A_e`, not `-2(A_e+E_e)`. The sharper alternative is
`sum |Gamma_e|<=2k tanh(beta)+Phi(E)`, also groupwise. Thus Gaussian noise
`E=epsilon G` preserves signed Frobenius diffuseness whenever `beta->0`
and `epsilon=o(sqrt(N))`, including bounded noise and arbitrary edge-local
adaptive selections. At fixed critical `c` and bounded `epsilon`, the
balanced-path squared-correlation error is `O_(c,epsilon)(sqrt(N))`.
No corresponding bounded-noise row-operator improvement is claimed.

For `Psi(G)=min_A log E_(sigma,x) exp(sum u_e(A_e+epsilon G_e) sigma x_i x_j)`,
the weak Hessian is the selected smooth Gibbs covariance minus a PSD
switch measure `M`. Its standard-Gaussian-weighted trace is at most
`epsilon^2 ||u||_2^2+epsilon sqrt(2(N+1)log 2)||u||_2`, hence `o(N)`
for critical profiles and vanishing `epsilon`. The proof keeps optimizer
adaptation throughout; an exact order-two cusp shows why switches cannot
simply be dropped. Mixed variance paths need a bounded relative velocity
for this trace estimate. Changing deterministic weights still produces
`sum u'_e A_e Gamma_e`; a heat identity does not remove that term.
The physical-noise-coordinate switch bound is divided by `epsilon^2`.
Do not claim that bounded covariance transport supplies different
deterministic block endpoints or a cross-order comparison.

The coordinatewise switch identity also gives, for GLOBAL Gaussian
envelope minima, `E L_g<=2K_g tanh(u_g)+epsilon^2 u_g K_g`. This is
distinct from the Boolean-energy bound for arbitrary edge-local choices.
The exact changing-profile derivative is
`sum u'_e E[A_e Gamma_e]+epsilon^2 sum u_e u'_e E(1-Gamma_e^2)`
minus `sum (u'_e/u_e) m_e`, where `m_e` is the Gaussian-weighted
diagonal switch mass. The first term remains uncontrolled, even when
monotone variances make switching favorable. In the noiseless case the
signed normalized edge-flip-gap sum differs from the previous `D_N(t)`
by at most `c^4`. A martingale sign-flip generator therefore reproduces
that defect; it does not supply its sign or another independent obstacle.

The preceding analytic result is
`evidence/NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md`
(SHA-256 `96d2675bd0cf1ee48e962b2974a2b8649afc487454a3912044ee1e737c53a9a5`).
Its full integration passed two independent reads. It proves uniform
sparse near-flat rounding, a common diffuse near-maximizing ensemble for
global norm minimizers, and actual signed Gibbs diffuseness for pressure
minima. Along the specified balanced two-block path, the exact formula is
`f_N(1)-f_N(0)=c^2/4-integral D_N+E_N`, with
`|E_N|<=c^3 sqrt(N)/2+c^4/6`. The shared-sign endpoint is exactly twice
the minimum half log-product of the one-sided partition functions, hence
`f_N(0)<=2P_(N/2)(c)`. A lower bound `integral D_N>=-o(N)` would give
dyadic pressure subadditivity; it has not been proved. Do not require
endpoint equality for that direction, and do not treat a bare dyadic
small-oh estimate as all-orders convergence. No pressure minimizer is
silently identified with a norm minimizer.

The added local corollaries give a rowwise signed operator bound and a
bounded cavity exponential normalizer. Slow cooling of pressure minima
produces asymptotically norm-optimal sources with actual near-maximizing
Gibbs ensembles. The missing unsigned star-fluctuation estimate does not
follow from these signed bounds or from the cavity normalizer alone.
These corollaries are analytically reviewed, not new finite-check coverage.

The previous fresh analytic result is
`evidence/NOTE_2026-09-05_INDUCED_OPTIMIZER_RESTRICTIONS.md`
(SHA-256 `ab65d46bb48627170344219850131aa77ed9cbe9d152e143346a7fec71d42409`).
Root and independent review checked the full proof. For `n -> infinity`,
`n^2=o(log N)`, the complete induced signing law has explicit total-variation
control. Every smaller signing occurs, but uniform restrictions have
typical normalized norm at least `(2/3)*sqrt(2/pi)>1/2`.
The failure result extends to `exp(o(n))` samples with uniform marginals,
even when dependent. Do not omit the growth/scale or marginal hypotheses.
Do not confuse existence of an `m_n`-optimal restriction with existence of
one matching the source constant; the latter comparison remains open.

Do not automatically resume residual (ii), the old equation (33), a skew
ansatz, or a finite-prime census. They are optional archived avenues, not
the definition of progress. Before revisiting one, name the changed
premise and the implication for the original question.

## Verification

The order-six optimized-profile runs are recorded in
`evidence/original_mo_optimized_profile_mesh.json`. A single NUKA exact
integer/Fraction run covered all 1,024 switching-normalized signings and
all 64 spin states, produced 23 joint signatures and 22 nonnegative
polynomial difference certificates, and passed 3,397 checks in `0.184`
seconds. An independent V100 run produced all joint histograms and the
prescribed 455 floating-point profiles in `2.938` seconds. A separate
NUKA comparison of the already stored outputs, without re-enumeration or
pressure replay, matched all 20,480 histogram entries exactly and passed
3,285 checks. Floating-point near-minimizer tolerances do not classify
exact ties or derivatives. All three runs exited 0 once; worker absence
was verified. Full result JSONs, the GPU array archive, reviewed sources,
exact commands and raw/preserved hashes are retained. The analytic theorem
does not depend on any of these computations and concerns only order six.

The new finite-step mesh run and independent replay are recorded in
`evidence/original_mo_finite_cross_mesh.json`. Soulkiller's V100 evaluated
8,192 canonical Gaussian-sign cross blocks, 8,192 independent blocks,
and the two coherent references on ONE fixed order-six conference host,
using all 4,096 spin pairs at 20 prescribed `(c,t)` profiles. This is
neither a new host-minimizer census nor an exhaustive cross-block search.
NUKA independently enumerated all 16 order-two cross blocks and the
four-point Gaussian support, then replayed the selected order-six GPU
pressures with full direct spin sums. Sample minima are only upper
bounds; sample log means need not approximate rare-event annealed values.
The V100 run completed once in `4.249` seconds, exit 0. NUKA passed
1,105 order-two formula checks and 480 checks on 160 selected order-six
cases; maximum CPU/GPU pressure difference was `7.11e-15`. This replay
checks the pressure and endpoint, not the GPU's `qbar` values. Both
workers exited normally and absence was verified. At all 20 order-six
profiles the best Gaussian sample is index 1067; exact signed-permutation
algebra identifies it with the known `A-I` construction for every
temperature and step; the complete identity is in
`evidence/NOTE_2026-09-05_SAMPLED_CROSS_BLOCK_ORBIT.md`
(SHA-256 `98923ba2cf14f71b71511b7896734028a48d1c029866fcc88592c40d820da1aa`).
Thus this sample found no new noncoherent winner,
not an exhaustive proof of optimality. Do not enlarge the same sample
without a changed mathematical premise.
The finite results are not all-orders evidence, and no larger sample
or unchanged successful run is required for a cleaner receipt.

The new fixed-order opposite-phase probe ran on soulkiller's V100
(`2.434` seconds, exit 0) and independently on NUKA CPU (`0.603` seconds,
22,528 checks, exit 0). It examined 1,024 switching-normalized order-six
signings at exactly `c=0.5,1,2,4,8`; no larger census was run. CPU used
64 spin states and covariance traces, while CUDA used 32 antipodal representatives
and direct squared bilinear moments. Their candidate values agree within
`2.85e-14` in `T`; all five profiles have the same 12 numerical minimizing
signings. Floating-point comparisons are not a rigorous optimizer
classification. The positive finite virial gaps at the tested `c=1,2,4,8` do not refute
an asymptotic small-oh allowance, and this check is not evidence for the
all-orders planted theorem. Both workers exited normally and absence was
verified afterward. Exact commands, source hashes, results, tolerances,
and cleanup receipts are in `evidence/original_mo_opposite_phase_n6_mesh.json`.

One new soulkiller run of
`scripts/original_mo_weighted_pressure_n4_check.py` passed 7,110 formula
checks, exit 0, with one CPU worker. Its scope was exactly order four,
64 signings, 16 spin states, and six prescribed weight/temperature
profiles. This is a finite regression of the new pressure identities,
not a larger-order census or theorem certificate. The exact command,
reviewed proof input, input hashes, full log, and live preflight are in
`evidence/original_mo_weighted_pressure_regression.json`. There was no
rerun of unchanged mathematics. The all-orders claims rest on the proofs.

The reviewed technical reset replay on soulkiller passed 40 tests in
44.83 seconds: the new global registry, both independence directions,
legacy route aliases, and three existing wrapper regressions. It ran in
`/tmp/original-mo-reset-replay.W2zk3m` with explicit files and one worker.
The separate documentation replay checks the final entry documents and
retained proof scopes. Its result, both exact commands, input manifests,
and log hashes are recorded in
`evidence/original_mo_route_reset_regression.json`.

The earlier diagonal work's receipt is
`evidence/original_mo_diagonal_regression.json`: 65 technical tests and
17 documentation tests passed across two runs after three missing staging
inputs were supplied. It is not a verification of the present reset.

No convergence claim may be accepted merely by toggling a Boolean or by
closing optional-route checkboxes. A complete reviewed proof is required.
