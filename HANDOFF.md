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
