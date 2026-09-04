# Research continuity rules

This repository is a long-running proof project.  Avoid rediscovering old
branches or turning completed computations into new work.

## Original-problem gate (overrides the finite Paley front)

MathOverflow 413935 asks whether `alpha_n` converges; a proof need not identify
the value.  Proposition 6.3 reduces this direct route to two Dini-summable
amplification estimates for `H(n)=m_n^(2/3)`, at multipliers 2 and 3;
polynomial saving is unnecessary. Propositions 6.4--6.5 give the exact
equal-endpoint skew diamond for a Hadamard doubling lift; its hereditary
endpoint conditions are automatic. Proposition 6.6 closes that diamond
outside the exact Hamming-central/joint-energy residue (6.20). Exact Paley
conference maximizers can satisfy every inequality in (6.20) strictly, so
Hamming geometry, joint-energy deficit, and hereditary extremality alone
cannot empty that residue. Any further doubling attack must use the fact
that `A` globally minimizes over signings or a genuinely finer `A`-dependent
construction of `R`. Proposition 6.7 also gives an exact equal-endpoint
tetrahedral tripling frame. Its single-skew distance-product and
distinguished-endpoint shields do not close the unshielded tetrahedral
diamond. Proposition 6.8 supplies a genuinely different `1:2` composition
using independently optimal orders `n` and `2n`; its bi-balanced Hadamard
cross block closes every pair with `k_A k_B <= n^2/100` and leaves exactly
(6.42)--(6.43). Fixed-anchor signature refinement cannot close that residue:
states alternating on every retained pair force the *existing spectral
bound* into (6.42)--(6.43), even for any Dini-admissible mildly growing
anchor list. A successful `1:2` attack must correlate the actual Hadamard
bilinear value with the two internal energies or replace the tiled cross
block. Proposition 6.9 kills the uniform
signed-Eulerian free-energy target for every fixed temperature `c>0`, not
only `c=2`; do not list `c=3` as viable. A growing `c_n -> infinity` route
is soft-max equivalent to the original problem unless it supplies genuinely
uniform new information. Fixed-physical-temperature Fekete subadditivity,
annealed centering, edge-noise monotonicity, and cavity monotonicity do not
control the critical diagonal `beta=c/sqrt(n)` and are not a live substitute.
Proposition 6.10 gives the correctly optimized critical-pressure gate:
convergence of `s_n(c)` on an unbounded set of `c` would prove the MO limit.
But common-raw-temperature block interpolation has a nonzero `c^2/8`
per-spin equal-split defect and shifts `c` to `c/sqrt(2)`.  Ordinary signed
cut/graphon convergence also cannot determine this pressure: conference and
deterministic quasirandom signings both cut-converge to zero, while the
conference limsup is strictly below the quasirandom limit for every
`0<c<1`. Proposition 6.10a also disproves the literal conference-product
lower curve at every `c>0` already at `n=5`; its universal entropy fallback
has large-temperature slope only `1/pi`. Do not retry that determinant,
entropy, hypercontractive, or Lee--Yang target as a route to slope `1/2`,
and do not retry Fekete,
Guerra--Toninelli, or ordinary graphon compactness without a new
speed-`n^2` optimized lower tail or a stronger second-order limit object.
Do not reopen endpoint selection on either frame, an independent
skew budget, a finite pair census, the implication from three pairwise
diamonds to the tetrahedral diamond, or a statewise random-skew union bound.
Do not present the disk surrogate as an equivalent target: its asymptotic
form would prove a stronger
`1/sqrt(2pi)` lower bound, while only its zero-error form is disproved.
Do not launch another finite-prime, residue, orbit, or cell census as work on
the original question.  The long residual-(ii) section below is continuity
guidance only for a deliberately selected Paley route.

Proposition 6.5c records a second, nonduplicate multiplier-two route:
`[[A,C],[C^T,-A]]` has the exact opposite-diagonal diamond
`max(|Q_A(x)-Q_A(y)|+|x^T C y|)`.  For symmetric cross blocks it is an exact
hybrid-slice interpolation within additive `n`.  Do not retry ordinary
holomorphic complexification: it forces `C_0=A`, reduces to QPSK clique
flips, and the optimal order-four signing already violates the zero-loss
`sqrt(2)` bound. Proposition 6.5n strengthens this to an infinite family of
complete signings with `Phi(A)=Theta(n^(3/2))` and a leading-order violation.
Therefore even optimal-order scale is insufficient; a surviving coherent
argument must use the exact global-minimizer property `Phi(A)=m_n` (or a
quantitatively near-minimal leading constant). A live noncoherent use must
choose a genuinely
`A`-dependent cross block and prove statewise cancellation.  The exact
four-label expansion has one preferred live subproblem: for
`R=A circ S`, prove that a tournament can make every signed outgoing
half-cut have absolute energy at most
`Phi(A)/sqrt(2)+o_Dini(n^(3/2))`.  It is individually possible on every
fixed face; simultaneity is the entire issue.  See
`evidence/NOTE_2026-09-02_COMPLEXIFICATION_OPPOSITE_DIAGONAL_AUDIT.md` and
`evidence/NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md`.

The 2026-09-02 direct calculator adds one exact all-orders invariant.  With
`mu_k=E|sum_(1..k) eps_i|` and `K_n` the maximum total `L^1` influence of a
normalized complete equimodular homogeneous quadratic,
`m_n=n mu_(n-1)/K_n`; hence the MO limit exists iff `K_n` converges.  The
same audit writes the doubling diamond as the slope-one two-half envelope and
as a restricted fourth-phase norm.  Its complete `n=5,...,8` orientation
models are finite illustrations only: do not extend that census.  Exact
Clifford flattening and independent random-skew first moments are now recorded
dead ends.  Even the ideal central-energy assignment has first-moment
threshold `sqrt(log(2)/2)=0.588705...`, above the optimal upper scale
`1/2+o(1)`; do not retry it with a partition-function refinement.  So are
ordinary zero-padding, random padding,
constant-block graph blow-up, generic balanced tensor/Grothendieck relaxation,
and weighted compactness followed by pseudorandom sign filling: the first
three have the wrong scale, and every dense filler has leading-order switching
norm.  The unrestricted weighted influence constants converge, but they are
not the flat constants `K_n` and normalized flat forms escape their natural
compactification.  A live use must prove an all-orders comparison internal to
the flat class, an `A`-dependent integral cut construction, or a genuinely new
joint Bernoulli lower-tail theorem.  Start at
`evidence/NOTE_2026-09-02_ORIGINAL_MO_TWO_HALF_GEOMETRY.md` and
`evidence/NOTE_2026-09-02_SHARP_INFLUENCE_TENSORIZATION_AUDIT.md`.

Proposition 6.5a gives the sharper exact orientation form
`K(A,R)=max_U Phi(A^(F_S(U)))`, where `S=A circ R` is a tournament and
`F_S(U)` is the outward-oriented half of the cut of `U`.  Thus the live
doubling problem is to orient the edges so that every directed half-cut
neighbor of a globally optimal `A` has norm at most
`sqrt(2)Phi(A)+o_Dini(n^(3/2))`.  Minimality supplies only the reverse lower
bound.  A transitive tournament turns this into prefix total variation, not
range; do not claim prefix switching or one-edge minimality controls it.  See
`evidence/NOTE_2026-09-02_ORIENTATION_STRUCTURE_AUDIT.md` and
`evidence/NOTE_2026-09-02_ORDERED_SKEW_PREFIX_HALF_CUT.md`.  Single-arc
local descent can be blocked tautologically by tied active cuts, and the
fractional minimax has the trivial exact optimizer `R=0`; see
`evidence/NOTE_2026-09-02_HALF_CUT_LOCAL_MINIMAX_OBSTRUCTION.md`.
Proposition 6.5d gives the exact nonlinear replacement: bad state pairs form
a fractional cover of tournaments in decomposable bivector space.  Do not
retry its exact affine relaxation, displayed elliptope, or normalized
single-row even-moment certificate.  The LP value is exactly zero,
covariance is subcritical from `n=16`, and that moment blind range grows
linearly with `n`.  Proposition 6.5j also blocks the full degree-four
preordering in the exact squared-row encoding from `n=45`, including affine
SOS localizers, all pairwise cross-row products, and every identity among the
instantiated Pluecker rows.  Proposition 6.5l extends the obstruction to
every fixed raw polynomial degree at sufficiently large order, and through
half-degree `(1-o(1))log(n)/log(log(n))`.  Do not retry that encoding at
another fixed raw degree.  Quotient-degree or differently lifted encodings,
substantially growing raw degree, an `A`-dependent nonuniform functional, or
direct nonlinear rounding remain open; see
`evidence/NOTE_2026-09-02_BIVECTOR_ENERGY_LAYER_MINIMAX.md`,
`evidence/NOTE_2026-09-02_BIVECTOR_DEGREE4_PREORDERING_NO_GO.md`, and
`evidence/NOTE_2026-09-02_BIVECTOR_GROWING_DEGREE_PREORDERING_NO_GO.md`.
Proposition 6.5e proves the outgoing-half target is already sharp from the
other side.  Its signed-regular arcsine bound gives leading lower constant
`1/(sqrt(2)pi)` and a positive correction equal to the off-diagonal square
mass of `A^2-R^2` and `AR-RA`.  Do not seek a fixed improvement below
`1/sqrt(2)`: dyadic recurrence independently forbids it.  The correction
forces an approximate commuting/equal-square skew mate only along a
near-`1/pi` subsequence; it is not globally coercive under the current
`alpha<=1/2` bound.  See
`evidence/NOTE_2026-09-02_SIGNED_REGULAR_ARCSINE_RIGIDITY.md`.
Proposition 6.5f controls all constraints incident with a prescribed
finite-anchor family: a signature-cell tournament costs `Ln`, absorbed at
the critical threshold for `L<=0.2636965...sqrt(n)`.  Its global skew norm
may still be quadratic, so the theorem remains conditional on an open
vertex-cover condition.  A new attack must control the non-anchor
high-difference layer or prove that exact condition in
`evidence/NOTE_2026-09-02_FINITE_ANCHOR_SIGNATURE_TOURNAMENT.md`.
Proposition 6.5k replaces signature counting by an exact weighted
Banaszczyk rounding and shields an arbitrary positive linear fraction
`k<((3-2sqrt(2))/(25pi))n` of prescribed anchors.  This is a real
all-orders capacity gain, but the dangerous graph depends on the resulting
orientation and `Gamma(R)` is uncontrolled.  Do not infer or claim the
same-`R` vertex-cover condition; see
`evidence/NOTE_2026-09-02_BANASZCZYK_WEIGHTED_ANCHOR_ROUNDING.md`.
Proposition 6.5g retires the generic scalar-defect spectral conversion.  Random
orientation already produces the approximate commuting/equal-square mate
required near `1/pi`, but the optimal trace/Frobenius-to-spectral bridge
loses `pi/2` even at zero defect.  Do not retry that conversion without new
statewise or special spectral information;
see `evidence/NOTE_2026-09-02_RANDOM_SKEW_MATE_SECOND_MOMENT.md`.
Proposition 6.5m separately rules out the exact commuting-conference
shortcut.  For every symmetric signing `A` and skew signing `R` of even
order, each diagonal entry of `AR-RA` is `2 mod 4`, so
`||AR-RA||_F^2>=4n` and exact commutation is impossible.  In
orthogonal-design terminology the desired relation is anti-amicability;
amicable symmetric/skew pairs anticommute.  This `O(n)` obstruction does not
exclude the approximate `o(n^4)` mate required by Proposition 6.5e.  See
`evidence/NOTE_2026-09-02_CONFERENCE_COMMUTING_MATE_NO_GO.md`.
Proposition 6.5h gives the exact first-moment criterion for a uniform random
orientation and proves that criterion diverges exponentially throughout the
optimal-signing regime.  This kills the literal union bound, not the random
orientation distribution; do not turn it into a probabilistic impossibility
claim.  Proposition 6.5i then reverses the arcsine proof itself: near the
`1/pi` floor its positive and negative Gaussian outputs are near-extremal and
Hamming-central, and any sharp outgoing orientation forces the calculator's
balanced two-half saddle.  Thus absence of the central opposite-energy layer
cannot close the gate.  This is necessary structure, not an orientation;
see `evidence/NOTE_2026-09-02_GAUSSIAN_SATURATION_CENTRAL_SADDLE.md`.
The multiplier `sqrt(2)` cannot be weakened on an eventual dyadic tail:
Proposition 6.5b shows that any fixed smaller multiplier with normalized
vanishing error would contradict the uniform positive lower bound.  Random
vertex orders also have a recorded all-orders variance obstruction; do not
repeat order sampling without a theorem controlling hereditary switched
cross-degree energy and the joint ordering tail.  Finally, fixed-real-part
Hermitian interlacing controls only one spectral edge; its matching-polynomial
and spectral-radius steps fail after adding `A`.  See
`evidence/NOTE_2026-09-02_RANDOM_ORDER_HALF_CUT_VARIANCE_GATE.md` and
`evidence/NOTE_2026-09-02_FIXED_REAL_HERMITIAN_INTERLACING_AUDIT.md`.
Do not recursively reuse the same real and imaginary halves through the
rank-one phase gluing `K tensor (A+iR)`: its restricted fourth-phase norm is
exactly multiplied by `4`, while the next doubling step permits only
`2sqrt(2)`, and paired-edge signs change the norm by at most `2n`. Any new
bisection attack must instead control the coupled cross-rectangle flips; see
`evidence/NOTE_2026-09-03_TWO_HALF_SELF_GLUING_OBSTRUCTION.md`.
The subsequent cross-rectangle Fourier audit proves an exact stability/Gram
rigidity theorem and exhibits an infinite Gram-perfect family at the spectral
maximum.  It therefore retires norm-only and Gram-only closures, but it does
**not** prove multiplier two: the statewise diagonal-payment inequality remains
open.  See `evidence/NOTE_2026-09-03_CROSS_RECTANGLE_FOURIER_STABILITY.md`.

## Before starting an attack

1. Read the current gate in `STATUS.md` and `HANDOFF.md`.
2. Search `evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md`, `src/`, `tests/`,
   and `evidence/` for the proposed object, parameter range, invariant,
   script family, and expected output.
3. State the one unresolved implication the proposed work would close and
   the existing result it strictly advances.
4. Do not launch the work if its distinguishing output is already recorded.

## Duplication gate

- A rerun is allowed only after identifying a concrete changed premise: a
  code defect, corrected theorem hypothesis, new parameter range, or new
  invariant.  Record that delta before running it.
- Do not rerun finite-prime, CP-SAT, MILP, SAT, orbit, OEIS, or literature
  searches merely to reconfirm a result.  Reuse their stored artifacts.
- If the canonical documents disagree, repair the contradiction before
  spending compute.
- Keep one active mathematical gate.  Closing a subcase does not authorize
  appending another search list; update the gate and reassess the proof cold.
- Prefer an argument that closes an infinite family over a wider finite
  census unless the finite census is the explicitly recorded gate.

## Type-I closure at Proposition 15.750

`type_I_multilevel_bad_case_ND_closed()` is True. For `p>=11`, the proof
uses isolated signed-PSL transport, square-direction Johnson rigidity,
central-Krawtchouk parity halving, and the sharp 15.688 lift floor. The
`p=5,7` bases are tracked exact integer Farkas identities verified without
SciPy or eigenshell caches; signed-PSL 2-transitivity handles every
distinguished edge. Do not reopen the `3A+B`, Aut_e, finite-LP, or small-prime
census routes as Type-I work. They remain incomplete historical mechanisms,
not live global gates. Residual (ii), E1, and the limit remain open.

## Current residual-(ii) gate after audited Proposition 15.771 and the post-15.761 exact reductions

Do not reopen `p=13,k=58` or `p=17,k=74`.  Proposition 15.742 combines the
`M_2=0` congruence with six multiplicative interval cuts and closes the
former row by the exact energy contradiction `667<707+26C`.  Proposition
15.743 extends the same common difference-Radon mechanism to the latter.
At `p=17`, the local mean formula and the global signed total are both
needed: together they force `P_L=4+k_L` in every hard direction.  Do not
silently assume that identity in a local cell model, and do not use the
normalized exact-star row `q=(2)^8` to obtain its own normalization.  First
glue the unspecialized exact-row sums to get `hT=18P_L-69`; common `hT`
makes the exact-star `P_L` common, and then `6P_L<=75` together with the
isolated-chart congruence `P_L=5 (mod 8)` forces `P_L=5`.  Only afterward
may one set `hT=21` and `q=(2)^8`.  The complete catalog of 698
translation-averaged nine-set cuts then makes hard excess one
infeasible, gives sharp excess-two and excess-three energies 70 and 119,
and forces every opposite row to `(-3)^8`, of energy 72.  The only partition
not already killed by the excess-one row satisfies
`119+9*72=767<1211+34C`, so it is impossible as well.

At `p=13,t=4,k=60`, Proposition 15.744 replays all residues and closes the
exceptional `u=3` profile by a six-root quartic contradiction after rebuilding
the edge-count-sensitive mass-14 models at `|H|=61`.  Its `b=10` premise is
also exact: a rank-78 restriction promotes contact-layer equality to the
pointwise complement triple, and a separate 1,716-variable punctured-lift
model excludes the two-unit `b=10` cell.  Do not replace that model with
Proposition 15.688, because the difference can be negative on the omitted
intersection layer.  Proposition 15.745
closes `u=0`: the 74-cut row bounds and common Radon energy give `C<=1` in
the sole difficult partition, while its seven parallel edges in six classes
give `C>=1`; its unique doubled parallel displacement then makes
the transverse multiplicities Boolean and bounds the elevated row in
`[-7,6]`, giving `695<719`.  Do not import the old `|H|=59` height-four
infeasibility, omit the collision-one sign audit, or call the full
`p=13,k=60` row closed.  At that stage its exact remaining residues were
`u=4,6`; Proposition 15.749 later closes `u=4`.

Proposition 15.746 completes the sharp mass-ten Boolean equality
classification needed inside `u=4`.  The pointwise all-positive `b=2`
quadrature must precede Proposition 15.688; only then is the lift globally
nonnegative, Boolean, and of support 330.  The exact 1,716-variable,
1,710-constraint model proves that the 78 omitted-pair and 286 all-equal
triple supports are exhaustive.  Do not rerun or broaden that catalog.
Common `hT` then forbids mixing: omitted pairs force `P=3`, all-equal triples
force `P=5`, and at least two opposite cells have mean 12.  At `P=3,Q=5`
the literal is impossible modulo six, leaving only a `b=0` mass-12 lift of
height one/support 396 or height four.  The seven hard directions also force
`F6=2hM6+hM2^3-3M2M4` identically zero, so that opposite cell must satisfy
`F6=0`.  In opposite normalization `N'_(2r)=(-h)M_(2r)`, encode this as
`2N'_6+(N'_2)^3+3N'_2*N'_4=0`; do not copy the hard-sign formula unchanged.
The `P=5,Q=3` branch retains the literal-or-lift dichotomy.  Its 22,308
patterns have full weighted feature ranks through degree six, ruling out an
analogous universal polynomial identity in `N2,N4,N6` at those degrees, not
every conceivable invariant.  Proposition 15.746 is an
exhaustive finite equality classification and proved open reduction, not a
close of `u=4`.

Proposition 15.747 supersedes the mass-12 gate.  Its weighted-cut second
moment excludes every Boolean mass-12 lift modulo seven, and exact
one-worker necessary-relaxation models exclude height four at `Q=3,5`.
Thus the omitted-pair `P=3` branch is closed and every minimum cell in the
all-equal-triple `P=5` branch is a `Q=3,b=12` literal. Proposition 15.748
uses the resulting common roots of `M2,M4,M6`: root count excludes `z>=5`,
exact interpolation excludes `z=4,3`, and `z=2` leaves 336 moment-level
survivors per hard sign. The only remaining opposite excess partition is
`(1,1,1,1,1)`. These survivors are necessary moment data, not common graph
realizations.

Proposition 15.749 closes that last `u=4` branch without a common-graph
census.  Every opposite `Q=4` row obeys all 74 translated-cut inequalities.
Two exact dual combinations give `-5<=q_a<=1` for every distance bin, so
the admissible row list has 522 elements and 492 moment triples.  Its
intersection with the nonroot evaluation alphabet of each 15.748 survivor
has 12 triples, all with fourth moment zero.  Five `Q=4` directions would
therefore add five roots to the two literal roots and force the binary
quartic `M4` to vanish identically, contradicting the hard fourth-moment
alphabet.  Reuse the pinned lists and hashes in
`src/e1_gmin_m4_prop15749.py`; do not rerun a graph, orbit, or cell census.

Proposition 15.751 closes `k=4p+6` for every prime `p>=13`: above height
one, paired cubes force a half-mean restriction of maximum at most three;
at height one, corrected transposition influences and the fixed four-bit
catalog exclude the required density. Proposition 15.752 extends the same
mechanism to scaled mass `p+9`, closing `k=4p+8` for every `p>=23` and the
contiguous band `4<=t<=(p-9)/2` for `p=1 mod 4` or
`4<=t<=(p-7)/2` for `p=3 mod 4`. These are proved infinite-family inputs;
do not reopen them with a prime, graph, orbit, slice, or cell census.

Proposition 15.753 closes the two sharp fifth-shell endpoints
`p=17,k=76` and `p=19,k=84`. It derives the complete A/B and A/C branch
lists before optimization, glues every row to the common signed edge total,
and uses all 698/2,338 translated cuts with exact one-worker energy
certificates. In the p17-A last partition, preserve the opposite sign:
`S4=-S2^2`, not `S4=S2^2`. Do not rerun either endpoint or import the
wrong-sign temporary maximum.

Proposition 15.754 closes the last `p=13,k=60` residue `u=6` without a graph,
orbit, coefficient-cell, or common-realization census.  Common row-sum glue
forces `P_L=4+k_L`, `hT=5`, and opposite count `Q=4`.  The seven partitions
of the hard excess five are exhausted by exact common forms: joint `U=hM2`
and `G=hM4-M2^2` energy tables for the two low-root partitions; a joint
`U,G,J6` coefficient/row join for the two four-root partitions; and the
identically-zero quartic `G` plus sharp collision-aware row energies for the
three high-root partitions.  Solver-free cut-equality and quartic-character
checks independently exclude the sharp equality cases.  Reuse the pinned
scripts, artifacts, and fail-when-wrong tests named by
`src/e1_gmin_m4_prop15754.py`; do not replace them with a graph reconstruction
or another row census.

Proposition 15.755 next gives the exact full-cube dangerous-spike reduction.
At `p>=11`, a shared maximizer has defect `2p` or at least `6p-12`, and both
values are attained by explicit all-prime Boolean families.  Do not collapse
the first shell to one-bit flips, claim a larger second gap, or count the
`A` and `B` hereditary cut intervals twice. Proposition 15.756 separately
proves that an arbitrary-boundary character/Weil cap is weaker than the
trivial cap at every even size at least four; two parallel affine lines are
sharp. Do not resume a D-only Weil, fibre-Parseval, or pair-deficit route.
Propositions 15.757--15.761 now make the common edge--Radon target exact.
The binary image has no obstruction beyond boundary and total parity, and
the recorded `p=1 mod 4` compact aggregate family passes it. Sharp
coefficient cancellation kills an
atom-count `l1` argument, and scalar Parseval has no gap on either infinite
local ray. Proposition 15.759 supplies every characteristic-`p` moment row;
15.760 proves those rows are sufficient for an unrestricted integral lift,
with cokernel `(Z/pZ)^S(p)` and no hidden linear congruence; and 15.761 proves
the exact full-target Moore--Penrose inequality also has strict room on both
rays. Do not retry parity, another linear/Smith obstruction, coefficient
`l1`, scalar or full-target Parseval, or Euclidean least norm on this family.

Proposition 15.762 proves a separate universal conference cube gap. For a
symmetric conference matrix of order `p^2+1`, `p>=5`, every Boolean vector
is a `+p` eigenvector or lies at least eight below the spherical ceiling;
the same holds for `-C`. Therefore a conference class with no Boolean `+p`
or `-p` eigenvector is already a complete `Phi_sph-8` certificate. Do not
rerun optimizers for the gaps two, four, or six. No such conference class is
presently constructed, so residual (ii), E1, and `L=1/2` remain OPEN.
The sporadic Peisert `G(23^2,2)=P**(529)` is already audited: its scalar
orbit is a 12-direction linear-OA/PN connection set, and explicit Boolean
`+23` and `-23` eigenvectors regularize it.  See
`evidence/NOTE_2026-09-04_SPORADIC_PEISERT_529_DEDUP.md`; do not rerun a
solver or rediscover this class.

Proposition 15.763 proves a signed affine-alias refinement of Proposition
15.755 by retaining cut signs and the exact active-state row sum. For
admissible odd `r`, it gives
`|H|>=oddceil((pr^2+1)m(m-1)/(2r(m-r))-pr^2+2)` with
`m=(p+1)/2+r`; at `r=1` this is `(p^2+11)/4`. If no noncrossing alias reaches
deletion eigenshell level two, replace `pr^2+1` by `pr^2+3`. Do not apply this
to every defect-`2p` point: the three-coordinate integral-eigenvector branch
is not classified as affine, and all-deletions minimality does not provide a
common affine coordinate system. Residual (ii), E1, and `L=1/2` remain OPEN.

Proposition 15.764 audits the previously implicit deletion-to-unit bridge.
For an all-deletions minimal four-gap set H, an odd `|H|` has a deletion at
signed eigenshell level two exactly when H has a signed level-three row. Frame
averaging and a level-five degree/anticommutator exclusion prove this whenever
odd `|H|<=5p`. The H-floor also freezes the distinguished edge on every
critical row; Paley `-C` is switching/permutation equivalent to `C`, and the
level-three bi-tight exclusion gives deletion size `k>=3p+1`. Thus this is the
full official residual-(ii) entry in that range. Even H cannot enter
residual-(ii) level two by parity; the
unclosed failure ranges begin at even `|H|>=4p+2` and odd `|H|>=5p+2`. Do not
claim that every minimal four-gap path is already covered by the four E1 units.
The formal max-of-affine countermodel in the proposition is a method barrier,
not a Paley construction. Residual (ii), E1, and `L=1/2` remain OPEN.

Proposition 15.765 kills the proposed affine classification of the first
defect shell.  At `p=11`, an exact Kiss--Somlai four-special-direction set
gives a nonaffine `D subset F_11^2`, `|D|=77`, with
`3+2 sum_(v in D) chi(v-u)=11(2*1_D(u)-1)` for every `u`; the corresponding
integral `+11` eigenvector has a unique coordinate `3`, while its Boolean
shadow has defect `22=2p`.  Four nonconstant line profiles prove that `D`
is not a union of parallel lines.  Do not retry the claim that every
first-shell point is affine.  This constructs no common all-deletions `H`,
so residual (ii), E1, and `L=1/2` remain OPEN.

Propositions 15.768--15.770 close the first two generic post-band layers;
do not reopen them with another prime, graph, orbit, slice, cell, or solver
census.  For `p=1 mod 4`, `p>=29`, the closed rows are
`t=(p-7)/2,k=5p-7` and `t=(p-5)/2,k=5p-5`.  For `p=3 mod 4`, `p>=31`,
they are `t=(p-5)/2,k=5p-5` and `t=(p-3)/2,k=5p-3`.

The reusable mechanism is exact.  Proposition 15.768 makes the
complement-triple equality pointwise, obtains coefficient offset two from
`5-1-1-1`, glues every hard direction to one common row sum, and forces a
local cell of scaled mass `p+15`.  Its height-at-least-two argument includes
`p=29` through the sharp dimension-free cube theorem
`E[g]=3/4 => max(g)<=6`; its height-one argument uses the already pinned
four-bit catalog.  Proposition 15.769 classifies the sharp `p-3` lift as an
omitted pair or all-equal triple.  With the two parity baselines the four
offsets are `2,3,4,5`; common-row normalization forces `P=offset`, and the
opposite ledger forces local mass `p+13`.  Proposition 15.770 carries those
classified equalities forward as `m-1` low rows and one high row, preserving
`hT`, and reuses the `p+15`, `p+13`, and sharp `p-1` local exclusions.  Do
not treat the fixed four-bit catalog as a new finite-prime census.

Proposition 15.769 also closes `p=23,t=9,k=110`.  Do not retain the old
local-witness wording as an open endpoint: the height-three, half-mean cube
equality is exhaustively `F_4/F_5`; additive `2x2` minors globalize those
forms on `J(23,12)`; and offset compatibility leaves only `F_5`.  The twelve
hard roots make the quartic and octic moment forms identically zero, whereas
the exact five-set replay checks all 33,649 subsets and finds 1,518 quartic
zeros, 2,024 octic zeros, and no simultaneous zero.  Proposition 15.770
also closes `p=23,t=10,k=112`: the one-row carry leaves eleven low hard
triangle-minus-star roots, still more than the octic degree, and at least
eight opposite mass-36 rows; the same offset test forces `P=4,Q=5,F_5` and
reuses the same disjoint quartic/octic certificate.  The accelerator vector
`[33649,1518,2024,0]` is an independent implementation check; the exact CPU
five-set/orbit certificate is authoritative.

Proposition 15.771 closes `p=23,t=11,k=114` for every boundary size.
The `b=4,20` contact equalities globalize by the elementary quadratic
fixed-weight kernel and mixed-coefficient pair sums. For `b=6,8,...,18`,
every slice point lies in a cross-boundary swap cube of dimension at least
five; even-half degree-two injectivity excludes the equality. The remaining
mean46 families have integral offsets4,5,6,7,8. The common row identity
must be applied before selecting an offset: it forces common `P`, then
`12P<=115` and `P=offset mod11` give `P=offset`. Both this all-one branch
and the exact quotient-zero branch force at least five opposite mass32
rows. At phase zero, floors leave only `b=0,2,22`; the last two require
pointwise subtraction of XOR or the omitted bit, leaving mass8 below the
sharp floor20. The `b=0` cell invokes the local15.752 mass`p+9` theorem.
The `u=9` two-unit carry retains ten roots and the same fixed five-set
certificate; `u=10` invokes the already proved mass22 local exclusion.
Do not reopen any p23 layer `t<=11` with another graph, slice, coefficient,
or five-set census. The next unclaimed p23 layer is `t=12,k=116`.
See `NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md` and its small-boundary
and middle-boundary proof notes. The four-node independent arithmetic
replay is `evidence/p23_third_post_band_mesh_replay.json`; it supplements
the analytic proof. Generic later layers and residual(ii) globally remain
open.

The global Mobius-incidence audit also proves that distinct target
directions, all signed fixed-word intersections, and physical ternarity
alone cannot close the remaining coupled box: for every branch prime it
constructs, probabilistically but exactly, pairwise-disjoint physical halves
with collision surplus at least the full endpoint demand. Do not retry a
standalone Bezout/incidence bound; it must be coupled to required physical
cancellations or the prescribed even-moment cells.

The live Paley target is ordered. Degree five and all odd rows are blind on
the antipodal compact rays. The all-prime odd--Radon follow-up now covers an
initial branch-C band. For every prime `p=4r+3` with `r>=7`, an opposite row
with `b` arbitrary compact atoms and `r-1` all-equal atoms has a central
aggregate signed edge chain whenever all odd global forms vanish and
`3b<=r+2`. For the deterministic balanced allocation, putting
`delta=t-(2r^2-4r-2)`, every opposite row is therefore central throughout
`0<=delta<=(2r+2)floor((r+2)/3)`. This is a structural reduction under the
zero-odd-form and balanced-allocation hypotheses, not an even-moment
exclusion; it says nothing about nonzero odd forms or unbalanced allocations.

The full balanced `p=4r+3`, `r>=7`, support geometry is also settled through
the cubic alternative.  With `h=2r+1` and support budget at most `3h-6`, no
realizable odd-Radon word can be supported on one maximal line.  Couvreur
peeling reduces any word containing `h` collinear points to one or two
maximal lines, and the two-line coefficient, parity, and capacity arguments
exclude every such pair for all `0<=b<=r`.  At the boundary `b=r`, where a
degree-three/degree-`h-2` complete intersection could have size `3h-6`, every
boundary cubic is excluded: the reducible, singular, and smooth cases all
fail for `p=4r+3>=31`.  These are proved support exclusions, not centrality
beyond the earlier `3b<=r+2` interval; the high-intersection irreducible-conic
alternative remains.

That conic alternative is now classified exactly.  A word containing at
least `p-3` points on an irreducible conic peels completely onto it; every
such high-intersection conic is triangle-tangent, and its dual weights reduce
to a constant orbit difference.  The star case is parity-impossible.  Every
nonstar survivor forces `q^3=1`, `q!=1`, equivalently `k^2=-3`, and hence
`p=7 mod 12`.  At `p=31,b=7,k=11` an explicit six-all-equal/seven-compact
edge witness realizes the 29 constant-conic edges and kills all 105 odd
channels, while its even syndromes are `F6=(11,19,10)` and
`F8=(12,11,23,6)`.  It is therefore an exact counterexample to odd-only
centrality, not a zero-degree-six/eight witness.

The entire `p=31,b=7,k=11` constant-conic fiber with simultaneous
degree-six and degree-eight syndromes zero is nevertheless `UNSAT`.  The
independent exact meet-in-the-middle certificate partitions total alignment
deficit three into `(3)`, `(2,1)`, and `(1,1,1)`, replays 230,314,710 maximal
completions and 17,076 exact edge hits, and finds zero degree-six/eight hits.
This is an exhaustive finite-fiber theorem only.  Do not rerun it, extrapolate
it to another prime or conic parameter, or call residual (ii) closed.

At `p=31`, the original all-prime theorem covered `68<=t<=116`; the frozen
component-packing upgrade below supersedes this by the guaranteed balanced
opposite-row band `68<=t<=164`. The independently audited finite
classification excludes the
specific one-compact plus six-all-equal row when its odd, degree-six, and
degree-eight global forms are all zero: 449 noncentered scaling orbits are
`UNSAT`, the unique centered orbit is covered by the earlier exact theorem,
and there are no `SAT` fibers. Every balanced profile with `69<=t<=99` has at
least one such row, so no profile in that band can have all of those global
forms zero. The stored `t=69` is provenance for the row certificate, not a
restriction on reusing its atom profile. Do not rerun the fiber search or
promote this balanced zero-form band to unbalanced allocations, nonzero or
coordinated global forms, or a residual-(ii) close.

Separate nonzero seven-channel Jacobian minors show that the unrestricted
four-compact and four-all-equal degree-six/eight atom maps are dominant over
the algebraic closure in every characteristic at least 11. Thus no universal
polynomial identity among those seven channels, and no purely algebraic
projective root-count argument, can close branch C. The resulting common
forms and labels exist only over an algebraic closure or finite extension;
they are not admissible `F_p` labels or form coefficients and prove neither
odd/higher-moment compatibility nor a Boolean lift. Remaining branches
include unbalanced zero-form allocations and whether actual `F_p` labels can
realize simultaneous nonzero even global forms coordinated across directions.

The nonlinear lift gate now has two exact, non-closing descriptions.  For a
signed integral lift `Rz=y`, the coordinate defect
`beta(z)=(||z||_2^2-H_y)/2` is a nonnegative integer, and a signed Boolean
lift exists exactly when `beta_R(y)=0`.  A lift minimizes the defect exactly
when `|2 z.g|<=||g||_2^2` for every element `g` of the **complete** Graver
basis.  No complete basis or defect-zero proof for the compact target is
known.  The explicit Type-P and same-square Type-K ridge circuits instead
generate `K_ridge` with
`p ker_Z R subset K_ridge subset ker_Z R` and
`ker_Z R/K_ridge=(Z/pZ)^nu_p`, where, for `m=(p-1)/2`, `d=p+1`,
`nu_p=d p m^2+m(m-1)(4m+1)/6`.
This gives an exact one-step saturation and fiber parametrization, plus
necessary ridge-descent inequalities, but the proper quotient proves that
the ridge system is not the complete Graver basis and does not turn an
integral lift into a Boolean one.

The equianharmonic component theorem is now frozen.  In the only surviving
constant-conic branch (`p=7 mod 12`), exact pairing components force the
necessary compact count `b>=(2r+7)/3`, for both constant signs.  Combining
this with the line/cubic/conic classification proves opposite-row centrality
under zero odd global forms for `3b<=2r+4`; for `p=11 mod 12` it proves
centrality for every `b<=r`.  At `p=31`, this is the guaranteed
`68<=t<=164` band.  The exact `p=43,b=9` odd-zero witness attains the
threshold but has nonzero `F6,F8`.  Do not claim the threshold is sufficient,
or extend the p31 band past `t=164` without an additional theorem.

The balanced hard rows are also structurally settled under zero odd forms:
each is a fixed unit star plus `e<=2r-2` compact atoms, and the compact
residual is central.  The whole hard row is not central.  Central inversion
therefore leaves only the fixed hard-star antisymmetric target.  Its exact
integral cokernel consists of the odd moments, and the direction-localized
Mobius trades give disjoint ternary lifts for arbitrary centers in all hard
directions.  Treat the **antisymmetric Boolean half as CLOSED**, but do not
call this a common graph: the live gate is the coupled symmetric half,
`s_e=1` on used trade orbits, `s_e in {0,2}` on unused nonfixed orbits, and
binary choices on fixed antipodal edges.

Do not retry the superseded hard-star support floor, equality-pencil,
total/fixed/parallel-count, or Euclidean-norm routes.  The constructive
Mobius lift already proves those diagnostics cannot obstruct the
antisymmetric target.  Also do not infer an even-moment obstruction from the
threshold components: at `b=(2r+7)/3`, exact `U,V` four-compact/two-cycle
trades preserve every odd edge channel, and their mixed seven-channel
degree-six/eight map has nonzero Jacobian
`2^32*3^26*5^2*7*2161`.  This is a dominance barrier, not a finite-field
zero-syndrome construction.

Canonical frozen notes are
`NOTE_2026-09-03_EQUIANHARMONIC_COMPONENT_PACKING.md`,
`NOTE_2026-09-03_HARD_ROW_COMPACT_ODD_RADON_CENTRALITY.md`,
`NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md`,
`NOTE_2026-09-03_HARD_STAR_ANTISYMMETRIC_SUPPORT.md`, and
`NOTE_2026-09-03_EQUIANHARMONIC_THRESHOLD_EVEN_BARRIER.md`, all under
`evidence/`.

The symmetric lattice theorem is also frozen. The plus-map is mod-two
surjective and has exact cokernel
`(Z/pZ)^[(h-1)(2h^2+5h+6)/6]`, precisely the even moment rows. This proves
only unrestricted integral and mod-two central lifts. Do not retry a hidden
Smith, parity, `p^2`, or other-prime obstruction, and do not promote either
lift to the restricted central Boolean box.

For the localized Mobius half, retain the exact formulas
`P_L=P_M=1`, `P_(L-M)=0`, `P_(L+mM)=1+eta(1+m)` and
`tau_t=eta(Q(e1-t^2(e1+e2)))`. The forced symmetric pair chain contributes
twice this parallel vector. Do not use the disjoint construction's count
`m(p-1)=|H|_max+1` as a contradiction: two arbitrary nonzero hard-star
trades can share exactly one origin orbit with opposite signs, reducing the
nonzero support from `2(p-1)` to `2(p-1)-2`. The remaining restricted symmetric fiber is still
open.

The all-active support theorem does prove `c>=p-1`, but only when every hard
center is nonzero in balanced zero-odd branch C. It excludes `c=p-2` via the
equality pencil, a prime-order Redei--Megyesi direction bound, and the
opposite-row parallel quota. Do not apply it when any center is zero or call
it a Boolean lift. Frozen anchors are
`NOTE_2026-09-03_INVERSION_SYMMETRIC_LATTICE.md`,
`NOTE_2026-09-03_MOBIUS_HALF_SYMMETRIC.md`, and
`NOTE_2026-09-03_ALL_ACTIVE_PENCIL_SUPPORT.md` under `evidence/`.

The restricted symmetric box has now been reduced further, without being
solved. In fixed/nonfixed coordinates
`R+=[[A,2B],[0,C]]`, and `A mod 2` is an isomorphism with explicit inverse
`a_[v]=g_(L_v)(0)+sum_L g_L(L(v)^2)`. The fixed binary word is therefore
forced exactly. Subtract it, divide only the fixed block, and retain the
Boolean variables on unused nonfixed orbits, the Hamming equation
`2 sum b=|H|-|U|-|a|`, and every exact direction-weight slice. Do not rerun
the already solved first fixed parity or replace these integer slices by
their parity.

For a used orbit `(a,[delta])`, the forced-word change is `Phi=0` for
`a||delta`, otherwise the `p`-point antipodal affine block
`{[delta+c a]}`. These blocks satisfy `M M^T=M^T M=I` over `F_2` and have
`p`-column lifts in `ker C`; hence the full unpunctured halved map
`D=(C,Phi)` is onto with rank `d h(h+1)`. Do not retry unrestricted halved
parity or the actual structured puncture: the grouped-uncertainty square
theorem and row-code gap now prove the latter onto throughout the balanced
zero-odd branch-C regime. Frozen proofs:
`NOTE_2026-09-03_GROUPED_UNCERTAINTY_SQUARE.md` and
`NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE_GAP.md` under `evidence/`.

Universal robustness for `|U|<=|Delta|` is false. Deleting
`X_(L,beta)={L(a)=0,L(delta)^2=beta}`, of size `p h=|Delta|-h`, drops rank.
The all-prime group-support theorem proves `d_row(D)=p h`, classifies these
rectangles as all minimum words, and gives an empty weight interval
`p h<wt<|Delta|`. Since every Hamming-extendable actual branch-C `U` has
`|U|<|Delta|` and the Mobius midpoint theorem excludes rectangle containment,
the structured puncture is onto. The live target is now the **integral**
zero-one equation with its prescribed Hamming and direction weights, not
another mod-two rank or all-halves-cover calculation. The one-difference-slice kernel is only the
local `A_(h-1)` whole-slab exchange lattice; do not promote local
connectivity to normality or existence.

Do not infer the Boolean lift from onto plus the scalar quota bounds. The
exact cardinality comparison in
`NOTE_2026-09-03_SYMMETRIC_QUOTA_CARDINALITY_BARRIER.md` proves that
some compatible syndromes with identical feasible quotas have no preimage;
the next argument must use the actual transverse target.

For the balanced all-active branch-C Mobius ansatz, direction parity now
strengthens the cancellation floor without a prime scan. Put
`m=(p+1)/2`, `s=(t+1) mod (p+1)`, and
`kappa=t_max-t+1+j`. Then `j>=2` for `5<=s<=m`, `j>=1` for
`s in {4,m+1}`, and `j>=0` otherwise. Thus do not retry the bare `j=0`
endpoint on `4<=s<=m+1`, or `j=1` on `5<=s<=m`; the latter is wholly
excluded. This is a quota/fixed-edge obstruction, not the remaining
transverse Boolean close. See
`NOTE_2026-09-03_MOBIUS_PARALLEL_PARITY_ENDPOINT.md`.

The remaining top residue now has a genuine construction-family close.  Let
`p=4r+3`, `m=(p+1)/2`, `t=t_max`, `s=m+2`, and `j=0`, and restrict to the
localized-Mobius family with an opposite fixed direction and the forced
auxiliary SDR of `m-2` hard plus two opposite directions.  For every
projective direction `K`, pair the half discrepancy
`q_(L,M,c)=boundary(E_(L,M,c))+1_(L=c)` with zero and one representative of
each nonzero antipodal pair in `ker K`; write the resulting sign as
`sigma_K=(-1)^g_K`.  In the coordinate `D(z)=M+z(L-M)`, the exact local rule is
`sigma_L=-1`, `sigma_(L-M)=epsilon_L`, and
`sigma_(D(z))=epsilon_(D(z^2))` for finite `z!=1`.  Pairing `z` with `-z`
gives the all-prime identity

`product_K sigma_K = -epsilon_M`.

The `m`-half product is `+1`, because `m` is even and the auxiliary signs
have product `+1`.  Any valid ternary reduction on a nonorigin inversion
orbit changes two kernel bits and preserves that product.  The auxiliary SDR
makes all origin orbits distinct, so no origin-pair correction is available;
the unique fixed antipodal edge changes one kernel bit and forces product
`-1`.  Thus some adaptive `K` is a degree-boundary contradiction for every
design in this family, including arbitrary higher nonorigin overlaps.  This
closes the top opposite-fixed `j=0` localized-Mobius family, not residual (ii)
or arbitrary lifts.  At `p=31`, an independent exhaustive replay checked
14,880 half options, 446,400 center instances, and 14,284,800 selector
pairings with zero failures.  The evidence file
`evidence/e1_gmin_m4_p31_top_mobius_boundary_parity_global_replay.json` has
SHA-256
`9ce413b62c65b7a541388f49bb390c690af44389556e6e8676ec138ef2cbc533`
and
`evidence/e1_gmin_m4_p31_top_mobius_boundary_parity_component_replay.json`
has SHA-256
`1e6cb15792a681f6e9b02290e60e0e14f16225282c56258c562564d1b0556650`
(all 2,969 searchable bounded-component designs witnessed; no universal
single `K`).  Replay with
`PYTHONPATH=src:. python -m pytest -q -n 0 tests/test_p31_top_mobius_boundary_parity.py`.

Do not transfer the `j=0` fixed-edge conclusion to `j=1`.  At `j=1` retain
both Hamming alternatives: `f=3,d=0`, including the origin-collision branch,
and the conditional `f=1,d=1` branch.  The next target is simultaneous full
degree-boundary and compact-atom synchronization across those branches, not
another fixed-`F` selector restart and not a QR/BCH detour.

Do not promote the complementary-profile endpoint sketch to a construction,
or reopen it as a route around the adaptive-selector obstruction.
Once its auxiliary scales are fixed it requires `M_i(x)^2=4j_i^2` on one
fixed-edge line, but every proved quota and moment condition permits the
nonzero centres `j_i` to vary independently. This shows nonautomatic
coherence only for a preassigned family; an adaptive center-dependent choice
remains open only in the `j=1`, hard-fixed, and arbitrary-lift branches. The
complete two-half calculation also
has clean one-overlap points (for example `q=r=2`), so the local four-candidate
intersection test gives no contradiction. Finally, the centered compact atom
`K(v,-v;0)` realizes a singleton fixed word and exactly one silent group
while satisfying the full common-moment system. These remain useful warnings
for the open `j=1`, hard-fixed, and arbitrary-lift branches, but they no longer
make the closed opposite-fixed top `j=0` family a live target. See
`NOTE_2026-09-03_MOBIUS_ENDPOINT_BARRIER.md`.

The older adaptive-centre matching question is partly solved, but is
superseded as a target for the closed opposite-fixed top `j=0` family. For
every branch prime `p=4r+3>=31`, every chosen
opposite fixed direction, and every list of nonzero hard centres, a common
fixed-edge magnitude and singleton signs give a perfect matching of the hard
targets with nonsquare complementary parameters. Each matched target pair
nevertheless forces its two auxiliary directions. The unresolved condition
is that all `m` forced directions be distinct and equal the prescribed set
of `m-2` hard and two opposite directions. In affine coordinates this is the
exact bijection equation recorded in
`NOTE_2026-09-03_ADAPTIVE_MOBIUS_PAIRING.md`; its necessary paired invariant
is `g_U^2=g_V^2` and `(U-V)^2/g_U^2=c^2/4` with one common ratio. Do not
confuse the proved target matching with this paired SDR. Even multiplicities
of the `g^2` values, equivalently a square fibre polynomial, are not enough:
the abstract set `{0,1,2,4}` with constant `g^2` has no common chord ratio.

The actual target atoms now give a second exact gate. At `j=0`,
`a_Y+Phi(U)=e_x` is equivalent to `z=c_U+ell+s_x`, and a direction with
`n_D` target triangles is feasible at the central fixed-word/odd-moment layer
exactly when the required row has weight at most `n_D` with the same parity.
Consequently `Lambda>=kappa_0+m+q` and, after separating zero-`Phi`
cancellations, `sigma>=kappa_z+m+q`. Two distinct halves share at most eight
nonzero `Phi` block types; an explicit disjoint ternary pair shares three, so
the tempting bound one is false. The scalar bounds leave positive room.
Those scalar inequalities remain method warnings, not a way to evade the new
degree-boundary contradiction.  For `j=1`, the live task is full
degree-boundary/atom synchronization, including the even moments and
nonfixed target cells, across both Hamming alternatives; see
`NOTE_2026-09-03_MOBIUS_FIXED_WORD_ATOM_COUPLING.md`.

The construction budget also applies across the whole all-active ray. If
`|U|=m(p-1)-2 kappa`, then `kappa>=t_max-t+1` is necessary. The disjoint
lift is therefore never extendable there. Two halves cancel at most two
orbits, sharply, but the two-cancellation locus is rigid
`q=r=1/2,A=B=3/2`, so there is no free greedy-pair parameter. The exact
first objective is
`|U|+|a_Y+sum_(O in U) Phi(O)|<=|H|`; passing it still leaves the divided
Boolean fibre. Frozen anchors are
`NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md`,
`NOTE_2026-09-03_SYMMETRIC_HALVED_MOD2.md`,
`NOTE_2026-09-03_SYMMETRIC_HALVED_MOBIUS_COVER.md`,
`NOTE_2026-09-03_SYMMETRIC_UNUSED_SLICE_EXCHANGE.md`, and
`NOTE_2026-09-03_MOBIUS_HALF_INTERSECTIONS.md` under `evidence/`.

The top opposite-fixed `j=0` localized-Mobius gate is closed by the adaptive
kernel product.  The live zero-odd continuation starts with the two `j=1`
Hamming branches and their coupled full boundary/atom equations, not the full
unsplit signed box.  For hard-fixed, unbalanced, nonzero-form, or arbitrary
non-Mobius branches the corresponding complete gates remain open. An
unrestricted integral or least-norm real lift is not a nonnegative simple
graph. No general common `0/1` graph has been constructed, and residual (ii)
remains open.

Focused replay for these post-15.761 records:
`PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. /home/nick/.venvs/mo-exact/bin/python -m pytest -q tests/test_prop15768.py tests/test_prop15769.py tests/test_prop15770.py tests/test_p23_post_band_moment_close.py tests/test_p23_second_post_band_moment_close.py tests/test_compact_ray_moment_gate.py tests/test_conic_odd_radon.py tests/test_p31_equi_zero68_mitm.py tests/test_signed_boolean_defect.py tests/test_ridge_kernel.py tests/test_equianharmonic_component_packing.py tests/test_hard_compact_odd_radon.py tests/test_hard_star_antisymmetric_support.py tests/test_inversion_antisymmetric_radon.py tests/test_equianharmonic_threshold_even_barrier.py tests/test_inversion_symmetric_lattice.py tests/test_mobius_half_symmetric.py tests/test_all_active_pencil_support.py tests/test_symmetric_fixed_edge_elimination.py tests/test_symmetric_halved_mod2.py tests/test_symmetric_halved_mobius_cover.py tests/test_symmetric_slice_exchange.py tests/test_mobius_half_intersections.py tests/test_p31_top_mobius_boundary_parity.py tests/test_main_chain_docs.py`.
Residual (ii), E1, `L=1/2`, and the original MO limit remain OPEN.

Consequently the fifth shell `k=4p+8` is closed for every prime `p>=13`.
With `q=(p-1)/2` and `k=4p+2t`, the exact residual-(ii) frontier is now:
critical `p=5,7`; `p=11,t>=3` (`k>=50`); `p=13,17,19,t>=5`;
`p=23,t>=12`; `p=1 mod 4,p>=29,t>=q-1`; and
`p=3 mod 4,p>=31,t>=q`.  The separate positive `p=7,z=7` branch also
remains open.  Propositions 15.743--15.749 and
15.753--15.754 are finite certificates or branch theorems, not an all-prime
row theorem.  The preferred generic front is therefore a structural version
of the common-energy/cut mechanism that survives when the number of distance
bins grows.  Another independent coefficient-cell catalog, one-direction
floor, halving heuristic, or longer complete-domain timeout does not advance
that gate.  Do not reopen any `p=13,k=60` residue, couple the former 336
`u=4` survivors to common graphs, or rerun the `u=6` common-form tables.
Before launching a finite p11/p13-later computation, identify the invariant
that could extend beyond that one row or explain why the finite row is a
genuine base obstruction.

## Result discipline

Label every result as exactly one of: proved theorem, exhaustive finite
certificate, open reduction, counterexample, or retracted claim.  Never
promote computation or a heuristic pattern into a theorem.

After a genuine advance, update the proposition-dedup audit and the canonical
status/handoff documents in the same commit.  Preserve failed routes and
counterexamples: they are part of the project memory, not clutter to delete.
