# Current mathematical status

Updated 2026-09-05. The original MO limit is OPEN.

There is no reviewed proof of convergence or nonconvergence. The value
`1/2` is unproved. These statements are independent of the status of the
optional Paley research program.

## What is established

`CORE.md` gives the route-neutral definitions and proofs. In particular,
`m_n` is nondecreasing, consecutive `alpha_n` differences tend to zero,
the limit-point set is an interval, and

`1/pi <= liminf alpha_n <= limsup alpha_n <= 1/2`.

The upper bound uses conference constructions, not their optimality.
Ratio-dense transfer and two-multiplier Dini amplification are valid
conditional tools; neither supplies its missing hypothesis automatically.

The new [selected-half restriction theorem](evidence/NOTE_2026-09-05_SELECTED_HALF_RESTRICTION.md)
also proves `m_(2n+1)>=2m_n` and `m_(2n)>=2m_n-(n-1)`.
These unconditional bounds use deterministic phase coloring and selected
subsets of actual larger-order minimizers. They improve ordinary
restriction monotonicity but still leave a normalized factor `sqrt(2)`;
they do not establish sharp proportional restriction or convergence.

## What changed at the reset

Residual (ii) is parked as an optional route-local open lemma. Its former
status as an unavoidable obstacle was a bookkeeping error. The same applies
to treating an E1/bi-tight conjunction as an if-and-only-if test for the
original limit. The gap-two condition is stronger than the asymptotic
optimality needed to transfer the value `1/2`; neither has been proved
necessary for convergence.

The active global proof registry is `src/original_mo_status.py`.
`src/e1_main_chain_status.py` retains separate optional-route diagnostics.
Historical proof notes, scoped local theorems, and exact certificates remain
available; unproved bridges have not been promoted or declared false.
See `ARTIFACTS.md` for the preserved branches and terminology.

## Fresh mathematical target

The [cross-spike/bulk milestone](evidence/original_mo_cross_spike_bulk_upper_milestone.json)
and [whole-profile proof](evidence/NOTE_2026-09-06_WHOLE_PROFILE_CROSS_BULK_UPPER.md)
give a strict ORIGINAL-target upper for an entire actual pure-cross class.
Take K=[[A,B],[B^T,-A]], positive diagonal D with D+/-K>=0,
S=tr D=O(n^(3/2)), and `delta=S tr(D^(-1))/(2n)^2-1->0`.
The evaluated final cells have ORIGINAL p=q_A=0 and common c=x^TBy.
Set W=D_L^(-1/2) B D_R^(-1/2), u_n=2c/S, m_n=tr(W^TW)/n.
If liminf u_n>=7/8 and every fixed R>9/10 has only o(n) singular
values of W above R, then
`limsup [E max_cell X_z-2sqrt(2)Phi(A)]/n^(3/2)`
`<=2/sqrt(5)-2sqrt(2)/pi<0`, for the centered-sign base Gaussian process.
The proof uses the completeness bound limsup m_n<=1/2 and CORE's actual 1/pi
source lower, not alpha=2/5. No full spectral-law limit, trace optimality,
or convergence of the largest singular value is needed.
The [actual-radius metric](evidence/NOTE_2026-09-06_ACTUAL_CROSS_RADIUS_SHELL_UPPER.md)
and [low-rank spike allocation](evidence/NOTE_2026-09-06_LOW_RANK_CROSS_SPIKE_MASS_UPPER.md)
retain the original source, covariance, and cell radius. The latter pays
for Gaussian projection and conditioning while retaining Boolean spike mass;
its general weighted-field theorem is not restricted to zero-source cells.
Its zero-source weak-middle-law corollary gives
limsup E max_cell X_z/(2n^(3/2))<=14/25 at x^TWy/n->4/5 without radius convergence.
The wider bulk<=9/10, liminf u_n>=4/5, limsup m_n<=2/5 region is diagnostic at alpha=2/5.
Evaluated corollaries retain the fixed trace cap, delta->0, and ORIGINAL
p=q_A=0; existing cell/padding errors are controlled, not discarded.
The bulk/high-u/dispersion premises are not established for all conditional
optimizers. Other cells and original convergence remain OPEN.

The [actual central-cell companion](evidence/NOTE_2026-09-06_ACTUAL_CENTRAL_CELL_LINEAR_WIDTH_BOUNDARY.md)
proves that the maximum of actual linear-field cell widths, divided by
2n^(3/2), tends to sqrt(kappa*w). No trace cap, dispersion, or optimality
is needed. Thus the CENTERED linear-field upper cannot reach the target
at alpha<=1/2; a shifted-sign version necessarily needs w<=2alpha^2/kappa.
This is an architecture boundary, NOT a lower on the Gaussian cross
process or an obstruction to other proofs. The next upper work concerns
biased central/active-cell competition or a sharper cross-process comparison.

The [all-law source-gain milestone](evidence/original_mo_all_law_source_gain_milestone.json)
and [adaptive proof](evidence/NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md)
prove for ACTUAL complete symmetric zero-diagonal A that
`limsup ||A||op/sqrt(n)<=5/3` and `limsup tr|A|/n^(3/2)<=4/5` imply
`liminf Phi(A)/n^(3/2)>=2/5+7/55000`, without a limiting spectral law.
The [two-moment transfer](evidence/NOTE_2026-09-06_TWO_CROSS_MOMENT_SOURCE_NUCLEAR_TRANSFER.md)
needs ACTUAL paired K=[[A,B],[B^T,-A]], positive diagonal
D=diag(D_L,D_R) with D+/-K>=0, and
`delta=tr(D)tr(D^(-1))/(2n)^2-1->0`.
For W=D_L^(-1/2) B D_R^(-1/2), Y=WW^T, retain the FULL actual
moments `m_D=tr(Y)/n` and `Delta_D=tr[Y(I-Y)]/n`.
If every accumulation point satisfies `9/25<=m<=1/2` and
`0<=Delta<=m[4sqrt(m)-3sqrt(1-m)]^2/25`, the same original-source
gap follows via a common large principal source. In particular
`liminf m_D>=2/5`, `limsup Delta_D<=1/1600` suffice with delta->0.
No full cross/internal law, separate trace cap, diagonal optimality,
or active state is assumed. The [endpoint transfer](evidence/NOTE_2026-09-06_CROSS_ENDPOINT_SOURCE_NUCLEAR_TRANSFER.md)
remains the Delta=0 case; the older near-flat theorem keeps its larger
gap at its narrower law. The source operator cap and paired dispersion/
moment premises remain conditional, and the paired field is unchanged.
The region is sufficient for exclusion, not necessary; prior formal
certificates retain their listed scope. No mathematical run was used.
Other regions, the all-cell implication, and original convergence remain OPEN.

The [original-source strict-gain milestone](evidence/original_mo_original_source_strict_gain_milestone.json)
adds a reviewed ACTUAL source-entry restriction beyond the older formal
relaxation below. For complete symmetric zero-diagonal A with
`limsup ||A||op/sqrt(n)<=5/3` and FULL empirical law of A/sqrt(n)
tending to `(9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4))`, it proves
`liminf Phi(A)/n^(3/2)>=5kappa/8+16/3125>2/5+3/1100`, kappa=2/pi.
The proof improves an actual positive Gaussian phase by a FIXED 10%
independent-coordinate Boolean update on the original source. Its
uniform Gaussianization is joint only for one local field and one
distinguished input coordinate, not all fields simultaneously. The
stronger higher-chaos mean and trace-of-square alignment are retained.
The update penalty uses the actual 5/3 operator cap; the weak-law atom
5/4 is not an operator bound or an exact finite-order flatness claim.

The separately reviewed internal-law transfer needs only an ACTUAL
paired K=[[A,B],[B^T,-A]], ANY feasible positive diagonal
D=diag(D_L,D_R) with D+/-K>=0, dispersion
`delta=tr(D)tr(D^(-1))/(2n)^2-1->0`, and FULL actual H_L law tending to
`chi_*=(9/25)delta_0+(8/25)(delta_(3/4)+delta_(-3/4))`, where
H_L=D_L^(-1/2) A D_L^(-1/2). A common original principal source
supplies the 5/3 cap and the law at 5/4; completeness determines the
scale without a separate trace cap, trace optimality, cross law, or
active-state premise. Its norm comparison transfers the same strict
lower to the original Phi(A)/n^(3/2). No paired covariance, cross
block, or active field is replaced.

Thus the specified ACTUAL near-scalar internal-law regime at objective
2/5 is excluded. The older strengthened FORMAL certificate calculation
is still valid for its explicitly retained relaxation, which omitted
this new source-entry information; it was never an actual-signing
counterexample. Its checkpoint below is preserved as history. Small
dispersion and this internal law are not established for arbitrary
extremizers; other profiles and the all-cell implication remain open.
No mathematical run was used, and original convergence remains OPEN.

The [weighted cross gain and strengthened formal boundary](evidence/original_mo_weighted_cross_gain_boundary_milestone.json)
complete the weighted transfer marked unpublished in the historical
scalar-gain checkpoint below. For ANY positive feasible D of the ACTUAL
paired K, with N=2n and the SEPARATE cap S=tr D<=C N^(3/2), let
dbar=S/N, delta=S tr(D^(-1))/N^2-1 and m,v_2 be the first two
moments of the FULL actual squared-singular-value law of W_D.
With kappa=2/pi and g_kappa=sqrt(kappa)-kappa, for delta<=1/512,
`beta(B)/(n dbar)>=kappa v_2/m+g_kappa m`
`-[25kappa C^2+6g_kappa]delta^(1/3)-R_C(n)`, uniformly R_C(n)->0.
No trace optimality, small canonical gap or global unweighted operator
bound on B is needed. The balanced complete block is only an auxiliary
norm lower; the original full measure and covariance remain unchanged.
On the SEPARATE actual p=q_A=0, c=Phi(K) face, u_D=c_D/n has
the same bound with another 2sqrt(delta) loss. Since m>=1/(2C^2),
near-flat full weighted laws force a positive leading gain above kappa
when delta tends to zero; small delta is still a hypothesis.

The accompanying FORMAL profile alpha=2/5, f=4/3, u=4/5,
m=9/25 and nu=(16/25)delta_0+(9/25)delta_1 passes the specified
strengthened full/source/cross inequalities, including the new gain.
Its SAME drift-plus-ellipsoid reference certificate nevertheless exceeds
71/125>sqrt(2)alpha for every signed metric and shifted threshold,
including endpoints; the full original drift is kept. The supporting-
line proof has squared target margin 41/15625 and uses no new run.
These are formal trace data, not a complete signing, a Boolean norm,
an actual active optimizer or a lower bound on actual Gaussian width.
Thus the gain alone does not close this strengthened relaxation.
Additional actual source/active/frame information or another upper
argument is still needed, along with the remaining all-cell work.
The older proofs are preserved; original convergence remains OPEN.

The [actual complete-cross gain](evidence/original_mo_complete_cross_flat_spectral_gain_milestone.json)
adds an actual-entry restriction beyond the retained cubic cross bound.
For an ACTUAL complete n by n sign matrix B and a SEPARATELY given
d>=||B||op, let m=n/d^2>=m_0>0 and
epsilon=1-tr[(B^T B)^2]/(n^2 d^2). Uniformly as n tends to infinity,
`beta(B)/(nd)>=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-o_(m_0)(1)`.
Here kappa=2/pi, and epsilon=0 is exactly flatness of the nonzero
singular values at d. Exact or asymptotic flatness forces a strict
leading-order gain over kappa when m stays bounded below.
A self-contained uniform marginal Gaussianization argument, including
all mixed contractions and the absolute-moment/tail passage, converts
the actual higher-Hermite variance gain into this norm bound. It does
not use variance alone, a joint limit of all columns, or a mathematical run.

On an ACTUAL pure-cross active state p=q_A=0, c=Phi(K)=beta(B),
the scalar ratio u=c/(nd) has the same bound. Thus the earlier formal
flat endpoint u=kappa is excluded only in this actual scalar,
bounded-operator context. Neither small diagonal dispersion nor a trace
cap is assumed to supply the required unweighted operator bound.
The weighted transfer is a separate next implication, not a published
consequence here; the all-cell comparison and original convergence
remain OPEN. The earlier formal-relaxation obstruction is preserved.

The [actual source/cross coupling](evidence/original_mo_source_cross_trace_boundary_milestone.json)
proves a finite cap-free inequality for every positive feasible D and the
ACTUAL cross law nu. With alpha=Phi(A)/n^(3/2) and m=integral y dnu,
its consequence as diagonal dispersion delta tends to zero is
`integral sqrt(1-y)dnu >= kappa sqrt(m)/(2alpha)-O(sqrt(delta)+1/n)`.
No optimality, small canonical gap, active cell, or source replacement is
needed for that coupling; its exact finite error is retained in the proof.

The same note gives a FORMAL trace/block relaxation at alpha=2/5,
f=4/3, u=kappa and m=9kappa^2/16 that passes the specified full/source/
cross moment bounds. Its drift-plus-ellipsoid certificate exceeds the
correct target sqrt(2)alpha for every shifted Gaussian sign threshold
and signed metric, including endpoints. The full original drift remains.
The formal data do not supply an actual signing or an active Boolean
optimizer, and the result is not a lower bound on actual Gaussian width.

The current trace-only certificate therefore needs additional actual
entry/active-state information or a different upper argument. No new
mathematical run was used; the existing pi enclosure was reused.
This scoped obstruction does not resolve the original inequality or
original convergence, which remain OPEN.

The [delta-only normalization](evidence/original_mo_delta_normalization_milestone.json)
extends the original-norm spectral lower to positive canonical gap.
For ANY feasible positive D with the SEPARATE cap S=tr D<=C N^(3/2),
let r=(N-1)/(S/N)^2 and mu_j be absolute moments of the ACTUAL full
weighted contraction. As delta tends to zero it proves
`2Phi(K)/S>=max{kappa(1+r)mu_3/(2r),kappa r/mu_1}`
`-O_C(delta^(1/3)+N^(-1/2))`, and `mu_2=r+O(delta^(1/3))`.
Its auxiliary principal submatrix only lower-bounds the original norm;
the actual full spectral measure, covariance and cross block are retained.
Trace optimality and small canonical gap are not required, but small
diagonal spread is NOT inferred from near-minimality. The paired use
retains the separate active original-zero-cell conditions. This reviewed
analytic bound still needs actual full/cross spectral coupling and the
all-cell width estimate; it does not prove original convergence.

The [evaluated small-gap package](evidence/original_mo_small_gap_evaluation_milestone.json)
adds a correctly normalized ORIGINAL-phase bound and an actual-law upper.
For D trace-optimal for the FULL SDP, S=tr D=tau(K), two positive
phases of the same full signing give, under a fixed original norm cap,
`Phi(K)>=kappa tr(D)/2-O(N^(3/2)sqrt(gamma)+N^(5/4))`, kappa=2/pi.
The finite bound is vacuous for gamma>=1/4. Small gap and the SEPARATE
active original conditions p=q_A=0, c=Phi(K) imply `u=c_D/n>=kappa-o(1)`.
For the ACTUAL squared-singular-value measure of W_D,
`m=u^2/f_n^2+o(1)`, c=f_n n^(3/2). For standard centered signs, the
concavity and convexity bounds at the fixed
metric t=3/5 evaluate its two-trace upper for every actual measure, not
a formal Dirac choice. When f_n tends to sqrt(2), the active small-gap
cell satisfies `limsup E max X_z/(2n^(3/2))<=17677/25000<1/sqrt(2)`.
Eleven predetermined exact Fraction comparisons passed once on soulkiller;
the earlier pi certificate was reused without a rerun.

The smaller-f formula, other original internal-energy cells, and positive
canonical-gap sources are still unresolved. In particular the intended
conditional bound `F<=2sqrt(2)Phi(A)` is not replaced by an absolute
sqrt(2)n^(3/2) bound when the original normalized minimum is below 1/2.
The two fully reviewed proofs are conditional advances, not a proof of
original convergence or a claim that all optimizers have small gap.

The [canonical-gap compatibility package](evidence/original_mo_gap_compatibility_milestone.json)
now controls the actual weighted/unweighted discrepancy in a quantified
small-gap regime. For ANY trace-optimal full-signing diagonal D, with
`S=tr D=tau(K)` and `g=S-tr|K|^3/(N-1)`, it proves
`delta=S tr(D^(-1))/N^2-1 <=4Sg/((N-1)N^2)` and
`Phi(K-(S/N)D^(-1/2)KD^(-1/2))<=S sqrt(delta)`.
This gives uniform actual source/cross energy compatibility without
a maximum-diagonal bound. On an ORIGINAL zero-internal-energy cell,
the weighted field is within `O(N^(3/2)delta^(1/4))` of a separately
positive pure-cross field that still uses the actual weighted cross block.
For `0<=delta<=1`, the all-shell metric theorem also bounds the exact
two-trace expression's change by
`3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|eta|)`.
Its scalar-I reference is a compared numerical functional retaining
the ACTUAL covariance and contraction, not a scalarized source law.
The eta window is fixed before taking the small-gap/large-order limit.
Both proofs passed complete independent reviews; no mathematics was rerun.

A norm cap bounds the normalized diagonal trace but does NOT here prove
the canonical gap small for original or conditional optimizers. Even in
the small-gap range, the remaining actual weighted trace supremum has
not been evaluated. That evaluation and a correctly normalized argument
for the complementary gap range are live targets, not established
conditions necessary for every convergence proof. Original convergence
remains OPEN.

The [same-source weighted covariance package](evidence/original_mo_weighted_covariance_milestone.json)
replaces scalar tensor normalization by a diagonal majorizer of the
literal complete block signing. Its Gaussian covariance is positive,
has unit diagonal and operator norm below three, and leaves A,-A and
the cross drift B unchanged. Under an original norm cap, the full
weighted Hermite correction and shifted-sign comparison give
`O(n^(16/11))` expected ORIGINAL norm error. This applies to every
actual conditional cross optimizer without a separate source operator cap.
The weighted linear-field theorem then proves positive comparison fields,
a polynomial-cell reduction for real weighted energies, and an explicit
two-trace upper using the SAME diagonal and the ORIGINAL cross-energy
shell. The live missing step is to evaluate that upper on actual coupled
original/weighted cells. Its displayed Delta_B and internal Delta_A terms
are not known to vanish or have a favorable sign.

Two additional reviewed results constrain the routes into that comparison.
For finite scalar-optimal templates `tau(C)=p q`,
[cubic frame alignment](evidence/NOTE_2026-09-05_SCALAR_TEMPLATE_CUBIC_ALIGNMENT.md)
strengthens the completion-certificate exclusion to
`q>=12/5 => Gamma(C)>283/200>sqrt(2)`. This crosses the restricted
weak-Dirac diagnostic barrier but does not turn Gamma into an actual
Boolean lower or infer its separately stipulated energy.
The [fixed-cap deflation construction](evidence/NOTE_2026-09-05_TENSOR_DEFLATION_FIXED_CAP_RATE.md)
gives actual source signings under EVERY fixed cap C>1/2 for which
the tensor positive-part Gaussian repair costs at least a positive
multiple of `n^(3/2)/sqrt(K)`, including after symmetric restriction
and diagonal removal. It rules out faster uniform fixed-cap repair
rates, not such rates at exact minimizers or for the full coupled
cross covariance. Five new rational comparisons were offloaded once;
the previous 28-comparison clipping certificate was reused, not rerun.
These four analytic results do not prove original convergence.

The [actual Hadamard-template package](evidence/original_mo_hadamard_template_milestone.json)
now realizes a flat singular bulk and finite outliers with asymptotic
scalar SDP optimality in complete cross sign matrices. Its finite Gamma
completion is only a Boolean UPPER certificate. Matched-frame rounding
proves Gamma>283/200 for q>=5/2 and scalar-optimal finite templates
tau(C)=p q, without turning that failure into an actual Boolean lower.
A modified invariant Boolean algebra supplies a
separate actual lower Lambda_I, identifies the actual limit exactly for
PSD templates, and excludes the small Boolean cap for its symmetric
positive-top-frame family at q>=5/2. Symmetric dilation preserves the
finite-template ratio but changes the actual family. All three complete
proofs passed independent reviews; their clipping constants share one
28-comparison remote exact-rational certificate. General nonsymmetric
attainability, conditional optimality and source compatibility remain
unproved. No conclusion about original convergence follows.

The [spectral rounding package](evidence/original_mo_spectral_rounding_milestone.json)
now constrains the actual original source, cross singular moments, and
optimal SDP diagonal. The original-phase bootstrap needs only an original
norm cap, with `O(n^(5/4))` norm error. Tensor-mixture rounding gives a
strictly stronger necessary relation between the Boolean/SDP ratio and
the canonical spectral gap, also without a separate operator-cap premise.
The [evaluated scalar diagnostic](evidence/NOTE_2026-09-05_SCALAR_MOMENT_FEEDBACK_DIAGNOSTIC.md)
puts the old strongest-feedback case below target after imposing repaired
positivity. Weak feedback still admits a formal moment law for which the
entire metric-optimized two-trace/Jensen bound misses target. That law is
not a constructed signing, optimal Gram or shell, and this is not an actual
Gaussian-width counterexample. The next information needed by this route
is actual optimizer/coordinate and source-cross compatibility beyond these
moments, or a stronger upper that uses it. All six proofs passed independent
reviews; neither this diagnostic nor the moment constraints prove convergence.

The [Boolean ellipsoid upper](evidence/NOTE_2026-09-05_BOOLEAN_ELLIPSOID_SHELL_UPPER.md)
now retains an explicit expected Boolean-coordinate penalty in the
Gaussian completion-square bound. It gives a stronger shell bound, recovers
the exact cube width at a diagonal metric, and completely evaluates its
diagonal-affine scalar specialization, including endpoint infima and the
actual independent Hermite cushion. The evaluated scalar diagnostic
improves the prior comparison but still misses the sharp leading target;
the exact covariance traces and attainable joint shells remain relevant.
The [cross-only regularization theorem](evidence/NOTE_2026-09-05_CROSS_ONLY_OPERATOR_REGULARIZATION.md)
keeps A,-A literally fixed and retains the two internal boundary penalties
incurred by cross-strip deletion. It gives a vanishing cap error for a
bounded-operator source, with explicit near-optimizer slacks. It does not
give the same uniform assertion for every exact original minimizer, and
its cap loss cannot automatically be absorbed into the current Gaussian
gain. Both complete proofs passed independent reviews. These are scoped
analytic advances, not a sharp all-shell comparison or convergence proof.

The [conditional joint-shell upper](evidence/NOTE_2026-09-05_CONDITIONAL_CROSS_JOINT_SHELL_UPPER.md)
now compares the actual Gaussian cross proposal with linear fields while
retaining the full exchange covariance and independent Gaussian cushion.
Its joint mismatch deficit uses all three attainable energies, with no
independence assumption on the maximizing signs. A separate masked-cross
proof gives the correctly directed conditional-optimizer Gaussian floor.
The [direct cross normalization](evidence/NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md)
uses `mu=max(2,||A tensor A-S_B+I||)` and gives `0<=R_mu<=2I`.
For any actual conditional cross optimizer over fixed A,-A, its Gaussian
floor has error `O(n^(16/11))` whenever `Phi(A)=O(n^(3/2))`; in particular
this holds for every exact original order-n minimizer. No conference-scale
operator bound is assumed. The [intrinsic joint-shell repair](evidence/NOTE_2026-09-05_INTRINSIC_CROSS_JOINT_SHELL_REPAIR.md)
restores positive comparison fields by an explicit rank-four PSD term
and preserves their leading joint-mismatch upper at error `O(n^(5/4))`,
uniformly in the threshold. All three complete proofs passed independent
reviews. The live step is a sharp evaluation using attainable shells and
conditional optimality; `F_A^*<=2sqrt(2) Phi(A)+o(n^(3/2))` is NOT proved.
Neither that optional dyadic target nor these reductions alone establishes
original all-orders convergence. The global status remains OPEN.

The [whole-edge optimizer constraint](evidence/NOTE_2026-09-05_WHOLE_EDGE_SOURCE_PRESERVING_GAUSSIAN_REDUCTION.md)
now rounds every unordered edge of the source at its ORIGINAL order.
The correctly normalized symmetric compression of `I-K tensor K/L^2`
has covariance operator at most three for every source at order at least
three; order two has an explicit independent-edge fallback. A universal
four-cycle bound controls the entire even-Hermite correction, giving
absolute expected ORIGINAL norm error `O(n^(16/11))`, uniform in the
threshold, without any source norm or operator cap.
For exact ORIGINAL norm minima this supplies a same-order Gaussian lower
constraint and an integrated CURRENT-posterior inequality with the same
subleading error. The full symmetric lift accounts for the diagonal as
one common scalar; the variance upper retains a negative source-energy
square and signed overlap subtraction. No pressure-minimum substitution
is used. The missing step is a valid order-transport upper comparison;
the same-order inequality cannot be reversed to supply it.

The [shifted-sign comparison](evidence/NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md)
now preserves arbitrary deterministic threshold means and controls ALL
Hermite orders, uniformly in the threshold. Its expected ORIGINAL norm
error is still `O(n^(16/11))` under a bounded latent covariance operator.
The [shifted-threshold covariance reduction](evidence/NOTE_2026-09-05_SHIFTED_THRESHOLD_COVARIANCE_REDUCTION.md)
retains the even-Hermite correction as explicit PSD low-rank Gaussian
noise, bounds its Boolean-norm cost, and reduces the actual threshold
law to `Z_h=s_h A+2 phi(h)G+sqrt(1-s_h^2-4 phi(h)^2)W`.
For ANY exact original minimizer A this proves
`m_(2n)<=inf_h E Phi([[A,Z_h],[Z_h^T,-A]])+D n^(16/11)`
with an absolute constant. The infimum is over deterministic thresholds
chosen before the disorder. It has not been evaluated at the required
leading constant; no sign of the mean/noise derivative is asserted.
Both complete proofs passed independent reviews. The original limit
remains OPEN.

The new [universal spectral-midpoint reduction](evidence/NOTE_2026-09-05_UNIVERSAL_SPECTRAL_MIDPOINT_GAUSSIAN_REDUCTION.md)
removes the source operator-norm hypothesis for a freely chosen rounding
law. If the source eigenvalue extremes are `a,-b`, choosing
`alpha=(a-b)/2` and `mu=(a^2+b^2)/2` gives `0<=Sigma<=2I` for EVERY
complete signing. The exact Gaussian covariance remainder is `O(1/n)`
with an absolute constant. Thus the normalized expected ORIGINAL paired
norm error is `O(n^(-1/22))` uniformly over all sources, including ANY
exact minimizer of `m_n`, with no regularization. In particular,
`m_(2n)<=E Phi([[A,Z],[Z^T,-A]])+D n^(16/11)` for every source A.
This is a proved one-sided Gaussian reduction, not an evaluated upper
bound in terms of `m_n` or a convergence theorem. The freely selected
midpoint is not asserted to be generated by a source Gibbs temperature.

The [energy-shell upper bound](evidence/NOTE_2026-09-05_GAUSSIAN_ENERGY_SHELL_UPPER.md)
reduces this actual Gaussian maximum to constrained one-block widths,
with an absolute `O(n sqrt(log n))` remainder for the midpoint law.
Its explicit energy-mismatch width deficit does not close the central
shell comparison. The [one-phase reduction](evidence/NOTE_2026-09-05_ONE_PHASE_GAUSSIAN_VARIANCE_UPPER.md)
costs only `O(sqrt(n))` pressure or `O(n)` raw expected norm, and its
Gaussian differential upper bound retains the CURRENT posterior energy
product and mean-overlap subtraction. The required integrated upper
comparison remains unproved. All three proofs passed complete independent
reads; no computation or signing census was used.

The new [correlated-sign Gaussian comparison](evidence/NOTE_2026-09-05_CORRELATED_SIGN_GAUSSIAN_FREE_ENERGY.md)
proves a genuine quenched-pressure equivalence, with error
`O(n^(17/18))=o(n)` for `n^2` cross coefficients at critical scaling.
It retains the full Gibbs posterior, permits singular latent covariance,
and requires only a fixed covariance-operator bound. The proof includes
the sign-smoothing estimate, not merely covariance matching.
For bounded-operator source signings, the
[exact canonical covariance linearization](evidence/NOTE_2026-09-05_CANONICAL_COVARIANCE_GAUSSIAN_LINEARIZATION.md)
then replaces the matched Gaussian by covariance
`(2/pi)Sigma+(1-2/pi)I` at an additional `O(1)` pressure cost.
The estimates also cover sufficiently slowly growing operator caps,
so the same-order regularization supplies genuine ORIGINAL-norm leading
near-minimizers to which the Gaussian reduction applies. The remaining
Gaussian endpoint upper comparison is unproved. No cross-order theorem
or source symmetric/half-product identification is asserted.

The [expected ORIGINAL paired norm corollary](evidence/NOTE_2026-09-05_EXPECTED_PAIRED_NORM_GAUSSIAN_EQUIVALENCE.md)
uses the explicit temperature-dependent errors, not a fixed-temperature
limit exchange. Its normalized expected maximum-absolute-energy error
is `O_K(n^(-1/22))`, uniformly over fixed internal energies, with constant
`O(1+K^4)`. Same-order regularization with threshold `n^(1/99)` gives
genuine original-norm near-minimizers for which both the normalized
objective loss and the Gaussian-reduction error are `O(n^(-1/198))`.
This is an expectation comparison; the doubled Gaussian norm is not yet
bounded by the smaller-order optimum at the required leading constant.

The [same-order spectral regularization](evidence/NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md)
works directly for the ORIGINAL norm objective. Every norm-capped host
has a complete signing at the same order with operator norm at most
`(K+8)sqrt(N)` and normalized norm increase at most
`2sqrt(Gamma C/K)`, where `Gamma=4pi/log(1+sqrt(2))`.
Operator-constrained optima thus approximate the actual norm minimum
uniformly as K increases. This is one-sided objective approximation,
not uniform perturbation control. Order transport in that class remains
unproved; the theorem does not itself imply convergence.

The [same-filler pressure extension](evidence/NOTE_2026-09-05_SPECTRAL_REGULARIZATION_PRESSURE_PROFILES.md)
preserves both actual phase pressures over an entire compact temperature
interval. The [quartically penalized profile](evidence/NOTE_2026-09-05_QUARTIC_PENALIZED_PROFILE_IDENTITY.md)
then permits every sign flip and approximates the original symmetric
pressure minimum within `O_c(lambda^(1/3))` after normalization.
Actual row optimality bounds every diagonal fourth moment; the norm
majorizer supplies a bounded fifth moment and signed Gibbs control.
The balanced identity has the exact penalty remainder `lambda(5-9/N)`
and an `O_c(sqrt(N))` error uniform for `0<lambda<=1`.
Its mixed internal/cross edge-gap integral remains unproved. Bounded gap
row sums and vanishing individual gaps do not establish its sign.
These are independently reviewed analytic results, not order transport
or identification of symmetric and half-product pressure minima.

The [whole-row and multirow reset budgets](evidence/NOTE_2026-09-05_QUARTIC_PROFILE_ROW_RESET.md)
now retain both the actual deleted law and the full tilted law; their
tail and relative-entropy statements are distinct. The
[quartic-force bounds](evidence/NOTE_2026-09-05_QUARTIC_FORCE_KERNEL_BOUNDS.md)
give a uniform operator bound for the actual weighted SIGNED Gibbs
kernel, not either unsigned phase covariance. The
[weighted row-tilt identity](evidence/NOTE_2026-09-05_WEIGHTED_ROW_TILT_FOURTH_MOMENT.md)
does not assert an endpoint fourth-moment bound. Independent
[quenched coefficient refills](evidence/NOTE_2026-09-05_QUENCHED_BIASED_COEFFICIENT_REFILL.md)
have `O(sqrt(N))` Bernoulli/Gaussian replacement error even over all
edges, with the exact quartic correction retained. Separately, the
[actual canonical full-strength law](evidence/NOTE_2026-09-05_FULL_STRENGTH_QUARTIC_PENALTY.md)
adds only `O_C(lambda n)` to the paired quartic penalty. None of these
finite variations bounds the remaining Gaussian endpoint pressure.

The [actual-Gibbs response](evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md)
and [leading half-product near-minimizer structure](evidence/NOTE_2026-09-05_HALFPRODUCT_NEARMINIMIZER_STRUCTURE.md)
prove `||A||_op=o(N^(3/4))` and uniform `o(N^(3/2))` Boolean-energy
change on deletion of ANY `o(N)` vertices. The actual complement phases,
unrestricted outside fields and all near-minimizer errors are retained.
These are independently reviewed all-orders theorems. The objective is
half-product pressure, whose zero-temperature slope is half the energy
width, not necessarily the absolute norm. No fixed-fraction comparison
or original-norm identification is supplied.

The [full-strength near-minimizer example](evidence/NOTE_2026-09-05_FULL_STRENGTH_HALFPRODUCT_NEARMINIMIZERS.md)
has actual leading half-product near-minimality and exact spectral
deficit `2(1-r)`. Its canonical singular cross law fails in mean and has
vanishing success probability, not an asserted exponential tail. This
does not settle the law on exact minima, on original-norm near-minima,
or for unrestricted selected outcomes. The [nuclear budget](evidence/NOTE_2026-09-05_NUCLEAR_SPECTRAL_BUDGET.md)
separately forces a linear-size spectral bulk, not operator flatness.
No new numerical experiment was used.

The [full-strength boundary identity](evidence/NOTE_2026-09-05_FULL_STRENGTH_BOUNDARY_LIKELIHOOD.md)
retains the exact pair-dependent Gaussian boundary law and the entire
posterior covariance integral, including its singular endpoint limit.
An analytic order-three actual minimizer disproves a coordinatewise
sign premise, not the averaged comparison; the offending boundary
context is absent at the endpoint itself. The full-strength weighted
comparison remains open. Complete independent proof reads passed.

The [fixed-strength weak-rounding theorem](evidence/NOTE_2026-09-05_FIXED_WEAK_GAUSSIAN_CUSHION.md)
extends the vanishing-strength exclusion to a NONEMPTY interval of fixed
positive strengths. Conditional Gaussian replacement, removal of the
conditional mean by convexity, and a variance interpolation give the
actual pressure floor
`[c sqrt(2t)K0-2log2-c^2 t arcsin(rho)/(2pi)]n-o(n)`.
When its excess over `cn` is positive, the actual centered law has
exponentially small probability of any successful internal host, even
one selected after the cross block. Uniform Gaussian concentration covers
every generating host, and `exp(o(n))` such proposals fail without an
independence assumption. This does not use low sign-law information and
does not exclude the singular full-strength `rho=1`. The complete proof
passed independent review; no new computation or sampling was used.

The [Gaussian-sign information theorem](evidence/NOTE_2026-09-05_GAUSSIAN_SIGN_INFORMATION_SCALE.md)
separates two dependent-rounding regimes. For any sign law with second
moment matrix `C`, its relative entropy from iid signs is at least
`||C-I||_F^2/(4||C||_op)`. Gaussian arcsine and Schur-product identities
extend this to singular Gaussian-sign laws. For the actual centered
construction `Sigma_rho=I+rho H/mu`, every fixed `rho>0`, including the
singular canonical `rho=1`, has `Omega(n)` information on norm-capped
hosts. Therefore the low-information exclusion below does not apply to
full-strength rounding. Conversely, `rho=o(n^(-1/2))` has `o(n)`
information on actual half-product minimizers and fails the finite-step
comparison in mean and probability at the stated fixed parameters.
These are independently reviewed all-orders results, not a proof of
success or failure of the full-strength law or of the original limit.

The [iid quenched cross-block theorem](evidence/NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md)
now controls the ACTUAL average log pressure, not only annealed proxies.
An explicit admissible Parisi control gives `P_SK>=4/(3sqrt(pi))`;
Gaussian covariance comparison and a direct Bernoulli replacement prove
a host-free pressure floor. At the stated sufficiently large fixed `c,t`,
even choosing the internal host after seeing the cross block cannot make
`exp(o(n))` iid-marginal proposals succeed. Any successful proposal law
needs relative entropy `Omega(n)` from independent fair cross signs.
This does not exclude unrestricted dependent selection or the earlier
full-strength Gaussian-sign law. These are independently reviewed
all-orders proofs, with no new simulation or finite census.

The [positive-cone cutoff theorem](evidence/NOTE_2026-09-05_POSITIVE_CONE_TRUNCATION.md)
proves that retaining only degrees `k<=o(N)` loses a positive normalized
pressure at fixed `c>pi log2`, including actual optimizers. It does NOT
exclude a sparse selection of extensive degrees. The
[extensive coefficient comparison](evidence/NOTE_2026-09-05_EXTENSIVE_COEFFICIENT_MOMENTS.md)
identifies the full coefficients with high even moments up to
dimension-uniform factors. The
[degree-selector identities](evidence/NOTE_2026-09-05_POSITIVE_DEGREE_SELECTOR.md)
retain minimization/maximization order; a genuine mixed exchange costs
only `O(log N)` in pressure, while a pure exchange has the stated spin
entropy cost. Their scalar coefficient-rate convergence would suffice
for the original limit, but no cross-order transport is proved.

The [exact optimized order-six profile](evidence/NOTE_2026-09-05_EXACT_OPTIMIZED_ORDER_SIX_PROFILE.md)
is an independently reviewed analytic theorem for every temperature and
every point of the balanced `3+3` path. A coefficientwise polynomial
comparison gives its global minimum over ALL complete order-six signings;
for positive cross weight, equality holds exactly at `A^2=5I`.
Its endpoints cross once as temperature varies, but its interior excursion
above the starting endpoint grows linearly in `c` as `c->infinity`.
This excludes a temperature-uniform bounded excursion, not a fixed-`c`
small-oh comparison as the order grows. The candidate polynomial is
quadratic only at this order; no all-orders higher-coefficient comparison
has been established. Independent NUKA integer certificates and V100
histograms agree exactly, with finite evidence recorded separately in
`evidence/original_mo_optimized_profile_mesh.json`.

The [finite-step rounding and annealing theorem](evidence/NOTE_2026-09-05_FINITE_STEP_ROUNDING_ANNEALING.md)
is proved and independently reviewed. It extends Gram--Schmidt rounding
to an integral finite-step pressure bound, keeping the exponential inside
the actual spin expectation. For fixed `c>0` and `0<t<=1`, every such
quadratic moment-generating-function proxy is at least `c^2 t n/4`.
Separately, the actual canonical Gaussian-sign annealed pressure is at
least `c^2 t n/(2pi)-o(n)`, uniformly over complete hosts. Against the
optimized paired endpoint `2R_n<=cn+o(n)`, these certificates fail at
`c>4/t` and `c>2pi/t`, respectively. This does not exclude a good selected
cross block, the average log pressure, the actual Gram--Schmidt law,
vanishing step sizes, or another order comparison. The finite mesh check
is recorded separately from the all-orders analytic proofs.

The [integral cross-block rounding theorem](evidence/NOTE_2026-09-05_INTEGRAL_CROSS_BLOCK_COVARIANCE_ROUNDING.md)
is proved and independently reviewed. For every complete signing with
actual Gibbs covariances at opposite temperatures, it selects an integral
cross block with phase-averaged squared bilinear cost
`min_B qbar(B)<=n^2-8(a_A')^2/(pi ||A||_op^2)`.
The Gaussian construction uses the exact negative spectral edge of a
centered tensor matrix and an entrywise arcsine identity. Separate
covariance-rounding bounds retain fixed coordinate squares and spectral
tails; a stated low-effective-rank hypothesis suffices for the local
endpoint comparison. That hypothesis is not established for minimizers.
The general radial comparison and the integrated order comparison both
remain unproved. These are integral construction bounds, not convergence.

The [near-minimizer counterfamily](evidence/NOTE_2026-09-05_NEAR_MINIMIZER_OPPOSITE_PHASE_COUNTERFAMILY.md)
is proved and independently reviewed for actual complete quadratic hosts.
For every fixed `c>0` and all sufficiently large `N`, it gives
`tr(A U A V)=Omega_c(N^(9/4))` at `beta=c/sqrt(N)`, with either norm excess
`Phi(A)-m_N=O_c(N^(11/8))` or half-product pressure excess
`a_A-R_N=O_c(N^(7/8))`. These are two potentially different families.
Thus even leading-order near-minimality does not imply the proposed
quadratic trace bound. The constructed hosts fail edge-local half-product
optimality; exact global-minimizer and actual edge-local trace bounds are
not refuted. This is not a cross-order comparison or a convergence result.

The [full-row cavity counterexample](evidence/NOTE_2026-09-05_FULL_ROW_CAVITY_COUNTEREXAMPLE.md)
shows that even the complete family of row-replacement inequalities does
not imply bounded tilted second or fourth moments for arbitrary symmetric
cavity measures. Its measure has full support and a strict minimizing row.
It has not been realized as an actual quadratic Ising cavity; it is not
a counterexample to a signing, a pressure comparison, or convergence.
This rules out proving the missing moment bound from row optimality alone.

The [adaptive-perturbation theorem](evidence/NOTE_2026-09-05_ADAPTIVE_PERTURBATION_CORRELATIONS.md)
extends the signed correlation estimates to actual noise-adaptive
pressure minima. With the perturbation held fixed during sign flips,
even bounded Gaussian noise preserves signed Frobenius diffuseness in mean
and probability, uniformly along the balanced profile, by a Boolean-energy
norm refinement. This does not assert bounded-noise row-operator control.
The independently reviewed
[Gaussian switch-measure theorem](evidence/NOTE_2026-09-05_OPTIMIZED_GAUSSIAN_SWITCH_MEASURE.md)
retains the negative optimizer-switch measure in the weak Hessian and
bounds its entire Gaussian-weighted trace by `O_c((epsilon+epsilon^2)N)`.
These remove two genuine obstacles to optimized Gaussian smoothing.
They do not turn a change of deterministic block coefficients into a
heat derivative; dividing by the noise scale can destroy the small error.
The exact mean-profile derivative retains its deterministic term, and
the stochastic sign-flip reset is the previous defect up to `O_c(1)`.
No new computation was used for these two analytic results.

The [global-optimizer variational theorem](evidence/NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md)
is proved and independently reviewed. Sparse rounding controls the whole
shrinking near-sign class; global norm minimizers admit one common
near-maximizing signed-state ensemble with small aggregate correlations.
Actual edge-local symmetric-pressure minima also have small signed Gibbs
correlations, uniformly along the specified two-block weight profile.
The exact optimized-path identity isolates a signed defect imbalance.
Its needed integral bound is unproved; endpoint equality is not required
for the correctly directed dyadic subadditivity implication. These are
structural theorems, not an all-orders comparison or a convergence proof.

The [induced-restriction theorem](evidence/NOTE_2026-09-05_INDUCED_OPTIMIZER_RESTRICTIONS.md)
is proved and independently reviewed. For `n -> infinity` with
`n^2 = o(log N)`, every `n`-vertex signing occurs in every source signing
of norm `O(N^(3/2))`. Nevertheless, a uniform restriction typically has
normalized norm at least `(2/3)*sqrt(2/pi) = 0.531923...`; for globally
optimal sources, `exp(o(n))` uniform-marginal samples cannot preserve the
source's leading constant with probability bounded away from zero.

This proves that smaller-order optimal restrictions exist, not that their
constant matches the source. The latter remains precisely an unproved
cross-order comparison. The theorem excludes only the stated sampling
method and scale, not biased selection or other order comparisons.

Compare different orders using the actual global-minimizer property,
without assuming a Paley model, a limiting value, or a prescribed lift.
Any claimed convergence proof must control the normalized optimum, not just
an arbitrary low-norm signing or a typical induced restriction.

No all-orders comparison with sufficient error control is established.
Finite exclusions and moment identities do not change that status.
