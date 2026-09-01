# Proposition and route de-duplication audit

**Date:** 2026-09-01

**Scope:** Propositions 6.3--6.9, every assigned proposition through
Proposition 15.751, the live predicate wiring, and the attack scripts present
during the audit

**Purpose:** prevent a reformulation, solver/backend change, longer timeout,
or finite-prime rerun from being mistaken for a new mathematical attack

## Verdict

The duplication concern was correct.

### Proposition 15.751: generic `t=3` is closed; do not census it again

The branch-B cell `4p E[B]=p+7` is empty for every `p=1 mod 4,p>=29`.
A dimension-free half-mean cube theorem excludes height at least two after
paired-cube averaging. Correctly normalized transposition influences force
the height-one Boolean case to be a six-coordinate junta; direct slice
symmetrization and cube influence reduce it to four coordinates. The fixed
four-bit catalog has 222 tables and fourteen profiles, none of the target
density. Four accelerators and an exact scalar replay agree on hash
`63c9daf2b117b540a5199b1b007cb4e6997ba01704fbc6017efaaa9735859396`.
Together with prior `p=13,17` certificates and branch-A/C arithmetic,
`k=4p+6` is closed for every prime `p>=13`. Do not launch another generic
`t=3` cell, graph, orbit, slice, or prime census. The live large-prime layers
start at `t=4`; the global residual predicate remains False.

### Original-question correction and two-ray gate (Propositions 6.3--6.9)

The binding plan had silently strengthened MathOverflow 413935 from existence
of the limit to existence plus identification of its value.  That made the
Paley E(1) architecture look mandatory.  It is not.  Proposition 6.3 proves
that Dini-summable amplification at only multipliers 2 and 3 for
`H(n)=m_n^(2/3)` forces convergence: iteration controls the semigroup
`{2^a 3^b}`, its multiplicative gaps tend to one, and monotonicity fills the
gaps.  Errors `O(n/(log n)^(1+epsilon))` in `H` already suffice; an all-pairs
composition theorem and polynomial saving are unnecessary.  Another finite
Paley residue census does not advance this direct gate.

The same review retracts Section 10's claimed universal `0.282` additive
doubling floor: it improperly used a lower bound inside a triangle upper
bound. The exact two-block identity leaves coupled design live. Proposition
6.4 puts the all-Hadamard two-cloud subclass into an exact four-state normal
form. Proposition 6.5 then chooses an equal-endpoint skew frame: both
endpoints are the same optimal signing `A`, so every hereditary cut bound is
automatic. Endpoint selection is therefore retired. The live condition is
the mixed-state diamond (6.13). Proposition 6.6 proves it, with Dini error,
outside the explicit Hamming-central/joint-energy residue (6.20). Do not
replace that residue by an independent skew budget, the stronger disk
surrogate as though it were equivalent, or a finite pair census. The
asymptotic disk is not disproved; only its zero-error form is.

Proposition 6.7 supplies the corresponding endpoint correction for tripling.
Its tetrahedral three-cloud frame has all four projective endpoints equal to
`A` and differs from its exact minimax by at most the `3n` internal-cloud
term. The single-skew choice reduces the ray to the cyclic three-state
diamond and proves distance-product and distinguished-endpoint spectral
shields. **Status: PROVED exact reduction; OPEN on the unshielded
tetrahedral diamond.** It does not close multiplier three and must not be
recorded as a convergence theorem. Exact order-two data rule out deriving it
from three pairwise diamonds, while the actual low-energy bulk rules out the
literal statewise random-skew union bound at the required constant. Neither
no-go excludes a correlated or `A`-dependent three-skew construction.

### Post-audit direct corrections (Propositions 6.8--6.9)

Proposition 6.8 supplies the first independent `1:2` composition reduction.
It is not nested doubling: the order-`2n` block is independently optimal and
the rectangular cross block is paired simultaneously against its positive
and negative extrema and those of the order-`n` block. A Paley-Hadamard tile
then proves the exact cross estimate
`|x^T C y|<=4sqrt(q k_A k_B)+6n` and closes every pair with
`k_A k_B<=n^2/100`. The remaining conjunction (6.42)--(6.43) is the exact
two-state gate. A fixed-anchor signature refinement is authorized because it
preserves the spectral core with only `O(n)` border; an enumerative pair
census is still not authorized.

Proposition 6.9 retracts the claim, repeated in the 2026-08-29 cold notes,
that the signed-Eulerian target remains viable at `c=3`. Symmetric conference
signings obey
`E exp(+-c Q_C/sqrt(n))<=cosh(c sqrt(1-1/n))^(n/2)`, so the required fixed-
temperature lower bound fails by a linear margin for every fixed `c>0`.
The only unexcluded version has `c=c_n->infinity` and would require uniform
control of its error. Fixed-`c` SOS, character, and shell variants are now
duplicates of a disproved target.

### Post-audit correction (15.720)

This audit itself missed a semantic dependency failure: it grouped routes by
proposition and formula, but accepted Proposition 15.55's final kernel claim.
That claim is false. If `R=G-(n/2)P1`, then
`ker R=span{1}+ker G`, and Proposition 15.56 already exhibits `n-2`
star-difference vectors in `ker G`. Therefore the spectral floor cannot close
bi-tight through 15.167, and GLOBAL QVAR/R1 are not E(1) acceptance gates.

Proposition 15.720 supplies the valid replacement without a new small-prime
run: `ker(Gsum)=scheme+cross` forces a bi-tight degree congruence modulo
`(p^2-1)/2`, excluding the required levels 2 and 3 for every prime `p>=5`.
After Proposition 15.750 below, exactly one mathematical gate remains in the
four-unit E1 ledger: non-Walsh residual (ii).

A second semantic check caught a nearby downstream misuse before commit:
15.274/15.585 invoked the bi-tight result on one-sided `S≡±4` tight covers.
That implication is invalid. Bi-tight level 4 is indeed excluded by 15.720,
and Proposition 15.402 explicitly constructs one-sided Max-minus-tight
level-4 covers: unions of four parallel square-direction lines. Combining
15.402 with the already-proved k=1 cylinder classification in 15.272 G /
15.588 E gives a Max-plus score at most `0` for every member of this family
(`0` at `p=5`, `-2p` at `p=7`, and `-4p` from `p>=11`). Thus the family is
not residual-compatible, but generic one-sided-tight emptiness is false.
The only live level-4 target is a one-sided tight cover that also has
`s_+=2`. The former
15.274 E dichotomy and 15.585 A `min_+=2` conclusion are retracted. Parameter
or solver searches based on either conclusion are not authorized gates.

A third semantic check found that the finite Type-I LP implementation in
15.408 E and 15.410 C did not encode its displayed inequality.  The old row
`Fm + 3*f_e[:,None]` evaluates to `S+3k f_e` because `1^T x=k`; the intended
bad-case inequality is `S+3f_e<=0`, equivalently `Fm x<=-3f_e`.  The shared
row builder and tests now check that equality algebraically.  One corrected
run from the existing eigenshell caches remains infeasible at both `p=5` and
`p=7`, so the finite conclusions survive, but the old solver statuses were
not evidence for them.  These runs do not create a general route and must not
be extended to another prime.

A later exactification and uniform audit now gives Proposition 15.750 and
closes Type I without extending that finite search. The `p=5,7` numerical
duals were used only to select supports; tracked positive integer multipliers
now verify `A^T lambda=0` and `b^T lambda<0` from a regenerated Paley matrix.
The `p=7` certificate is a full nonnegative-cone certificate, while `p=5`
also uses the cardinality side condition. For every `p>=11`, isolated-chart
square-direction rigidity forces the signed mass profile, a nonsquare
direction has scaled lift mass four or six, and exact parity halving makes
the sharp 15.688 floor contradictory. Signed-PSL 2-transitivity normalizes
the distinguished edge. This is a proved all-prime theorem, not a census;
do not reopen the old `3A+B`, Aut_e, or floating-LP routes.

A fourth semantic check found the main source of the all-finite residual
duplication. Proposition 15.267's signed PSL action can move any selected
odd-boundary vertex to infinity while preserving the relative flip set up to
permutation, its size, and both shell-separation inequalities. Applying that
normalization before the existing infinity-present theorems gives Proposition
15.721: for every prime `p>=17`, 15.669 excludes every total boundary size
`6<=|D|<=p-3`, and 15.674 excludes `|D|=p-1`. Together with the old
`0/2/4` closures, every `|D|<=p-1` is impossible. The first unresolved
general shell is `|D|=p+1`, normalized to infinity plus `p` finite points.
Accordingly, the boundary-close role of the first/second all-finite campaigns
in 15.675--15.712 was redundant. Their internal lemmas remain available;
15.676 is still load-bearing on pair-deficit equality at `|D|=p+1`, and
15.690--15.691 are independent optional no-go results.

A fifth semantic check resolved the first new `p+1` subbranch without a
profile campaign. Proposition 15.722 tracks the signed phase cocycle across
all charts. Outside pair slack one is impossible; slack zero forces a
Miquelian circle with exact type alignment. Proposition 15.724 then reuses
three existing inputs—isolated-vertex counting, the xnor congruence from
15.672/15.673, and the sharp lift floor from 15.688—to exclude that full
circle. Do not launch another full-circle, `R=0`, conic, or circle-orbit
search at `|D|=p+1`. The same proposition now excludes every positive slack
through `max(3,floor(sqrt(p)-5/2))`; the active branch lies beyond that
cutoff at this stage of the audit. Proposition 15.726 below supersedes that
active cutoff without changing the historical 15.722 result.

A sixth semantic check found that `excess != 2` is not a universal profile
rule in the middle odd-fibre range. Proposition 15.723 proves the exclusion
by paired cubes except at the genuine cells
`(p,b,phase)=(17,5,1),(17,11,0)`, both realized by explicit integral
quadratics. Any script that deletes those cells by a blanket condition is a
relaxation bug until audited. Endpoint baseline uses may still be valid, but
must cite their pointwise baseline factorization rather than the blanket.

The tracked backward audit has the following verdicts:

| prior use | verdict after 15.723 |
|---|---|
| 15.674 and `tests/test_prop15674.py` | This odd-profile use reaches `p=17`. The code now retains `(b,phase)=(5,1),(11,0)` explicitly. At residue zero either cell needs quotient two, while all other directions need at least one and the total quotient sum is only `m`; at interior residues all `m` directions still need at least one but the sum is below `m`. The theorem and its four arithmetic rows are unchanged. |
| 15.676 | Its conic profiles use only `b=1,3,p`; neither exceptional cell occurs. |
| 15.675 and 15.679--15.683 | Their stated conclusions survive a parameter-aware replay. Proposition 15.723 handles reduced parity rank at least five for `p>=19`; reduced ranks three and four must be retained, but they alter only over-budget rows in these propositions. Their exact historical row ledgers are being regenerated. |
| 15.678 | **OPEN_RETRACTED_REDUCTION.** The corrected census has 108 compatible profiles spanning 47 arc profiles. The retained geometry excludes 14 arc profiles, leaving 94 compatible profiles uncovered. The old “exactly two profiles, both arcs” endpoint claim is false. Proposition 15.721 independently closes this all-finite boundary as a gate. |
| 15.684 | **OPEN_RETRACTED_REDUCTION.** Restoring the admissible phase-zero residue `u_0=9` gives scaled mass 18 and an explicit slack-zero profile, so the old positive-residue exclusion and total `1,247 -> 203` endpoint reduction are false. The exact residue-zero census and its reductions remain useful, and Proposition 15.721 independently closes this boundary as a gate. |
| 15.688 | The sharp lift theorem and its direct residue-zero census survive. The corrected generic ledger restores `u_0=7`, of scaled mass 14, but the sharp floor 16 excludes it before the residue-zero census; the final p=19 fourteen-profile block is unchanged. |
| 15.700--15.712 | **Corrected replay completed.** Propagated census IDs give `2503 -> 2219 -> 1744 -> 1481 -> 1368 -> 1228 -> 1215 -> 1213 -> 1020 -> 869 -> 321 -> 19 -> 14 -> 0`. In detail, 15.700 excludes 284 and sends slack zero `286 -> 2`; 15.701--15.704 exclude `475,263,113,140`; 15.705 is **PARTIAL/OPEN** and removes only 13 historical Orbiter targets, leaving 74 slack-16 rows; 15.706--15.712 exclude `2,193,151,548,302,5,14`, with 15.709 absorbing all 74 leftover slack-16 rows. The final nineteen- and fourteen-profile blocks are unchanged, and 15.712 still closes the endpoint. The whole all-finite size-16 ladder is superseded as an active gate by 15.721. |
| 15.724 | Its endpoint `b=2` use is the pointwise baseline factorization `A=2B` and the sharp 15.688 support floor, not the disputed middle shortcut. |

The generic `scripts/infinity_plus_p_quantized_dp.py` now routes every
two-unit excess through 15.723's parameter-aware
`floor_excess_admissible()` classifier.  In particular it retains the two
real `p=17` exception cells and the still-unproved reduced-size-three/four
cells; the former blanket filter is retired.  This repair does not turn the
script into a proof or reopen a live acceptance gate: 15.721 independently
supersedes all of the affected all-finite boundary campaigns.

A seventh semantic check retracted Proposition 15.725's attempted
parabola-plus-internal family close. Its 2,381-case finite phase-zero census
is exact, but the `p>=53` character-curve bounds are asserted rather than
proved, the admissible singular locus `4*a*nu+1=0` is untreated, and the
opposite product sign is unchecked. Do not cite 15.725 as an all-prime or
two-orientation exclusion.

An eighth adversarial check repaired four proof-certificate defects in the
new 15.722--15.724 chain without changing its valid conclusion.  The signed
Möbius cocycle now handles affine maps `c=0` separately instead of assigning
the impossible multiplier `chi(0)` at infinity; all finite-field APIs reject
odd composite moduli.  In 15.723 the far-contact active-coordinate minimum is
correctly `k-1` for odd `k` and `k-3` for even `k`; both remain at least five
in the stated range.  In 15.724 the imported two-coordinate baselines are
honestly XOR/XNOR rather than both XNOR: the needed congruence survives
because the sign parameter drops out of `(p-1)c=I+P_d-4`.  Independent
symbolic checks of the paired-cube operator, quadrature weights, gap
factorizations, and full-circle `(u,x,y)=(4,4,3)` arithmetic found no further
closure-affecting gap.  A subsequent exact finite-geometry check strengthens
the valid result: outside slack `R=0` is closed by 15.724, and 15.722 excludes
every positive `R<=max(3,floor(sqrt(p)-5/2))`.  The `R=2,3` cases use the
classified complete `(p-1)`/`(p-2)` arcs; the prime-dependent interval uses
an inclusion-minimal deletion to an arc and off-conic secant counting. Only
slack beyond that cutoff remained open at this stage.

A ninth semantic check gives Proposition 15.726 and strictly advances that
positive-slack gate without a finite-prime campaign.  For an outside
`p+1`-point set of slack `R`, let `T` be an inclusion-minimal deletion to an
arc `A` and put `t=|T|`.  The exact occupancy identity gives `1<=t<=R`,
minimality gives every `z in T` an `A`-secant, and the total number of these
deleted-point/secant incidences is at most `R`.  Hence each
`s_A(z)<=R-t+1`.  The arc has size `p+1-t` and tangent parameter `t+1`.
When `3R<=p-4`, Ball--Lavrauw's odd-order tangent envelope applies and has
dual degree `2(t+1)`: its size hypothesis follows from
`p+1-t-(2t+4)=p-3-3t>=p-3-3R>=1`.  Every deleted point lies on at least
`p-1+t-2R>2(t+1)` tangents.  Its dual line would therefore be a component of
the envelope, contradicting the nonzero tangent-polynomial value at an
`A`-secant through that point.  Thus for every prime `p>=17`,
`1<=R<=floor((p-4)/3)` is impossible.  Together with 15.724 at `R=0`, any
positive survivor must have
`R>=floor((p-1)/3)`.  This narrows but does not close the `p+1` shell:
residual (ii), multi-level Type I, and `L` remain open.

A tenth semantic check gives Proposition 15.727.  At the first integer left
by 15.726, `R=floor((p-1)/3)`, choose a minimum-cardinality deletion `T` to
an arc `A`.  The tangent-envelope bound excludes every `|T|<R`; equality in
the slack incidence count then forces `|T|=R`, every deleted point to have
arc-secant index one, and every rich line to be a pairwise `D`-disjoint
trisecant or 4-secant.  Hence `c_1(A)>=R`.  Exhaustive published arc
classifications and exact representative audits contradict this at
`p=17,19,23,29`, moving their first possible positive slacks to `6,7,8,10`.
At this stage the first endpoint not excluded was `p=31,R=10`, in the
disjoint rich-block normal form. Propositions 15.733--15.734 later supersede
that frontier.

An eleventh semantic check gives Proposition 15.728. At `p=31,R=10`, the
exact identity `sum_d b_d=72`, the two type budgets, common-residue
quantization, and the sharp integral-lift floor force one Paley type to have
scaled means `{30^15,62}`. At least fourteen of its directions have `b_d=2`.
If `y` is the number of 4-secants in 15.727, at least `4+y` of those
directions are nonrich and have fibre profile `(14,2,15,0,0)`. This is a
**proved necessary normal form, not an endpoint exclusion**. It changes no
acceptance predicate and does not authorize extrapolation from a `p=31`
profile to all primes.

A twelfth semantic check gives Proposition 15.729. In any remaining 15.727
endpoint block decomposition, retain three points on one rich block and two
on every other rich block. Exactly `R-1` points are deleted. The retained
affine set is a `(p+2-R,3)`-arc with exactly one trisecant; deleting any two
points of that trisecant leaves a `(p-R)`-arc whose deleted points are two
distinct extensions on one tangent. For `p=3R+1` the two sizes are
`2R+3,2R+1`; for `p=3R+2` they are `2R+4,2R+2`. This is a **proved
all-prime necessary reduction, not endpoint closure**. At the 15.729 stage,
the next implication was exclusion or classification of these near-extremal
affine unique-trisecant configurations while retaining their common
disjoint-block completion `D`; 15.730--15.732 supersede that formulation.

A thirteenth semantic check gives Proposition 15.730. Every maximum arc
contained in the same endpoint completion retains all singleton points and
exactly two points on each rich block. Hence there are exactly `3^x6^y`
maximum repairs. The complement of every repair is an `R`-arc all of whose
points have repair-secant index one. The full two-colour projective line
census, all point signatures, and the direction refinement
`b_d=p+1-2(sigma_d+tau_d+m_d)` follow exactly. A 4-secant supplies bases
with three displayed co-tangent extensions inside `D`. This is a **proved
simultaneous necessary normal form, not endpoint closure**.

A fourteenth semantic check gives Proposition 15.731. The squared
Segre-normalized tangent products on a repair's dual lines agree at every
node and glue constructively to a degree-`2(R+1)` plane form. In residue
`p=3R+2` the normalized envelope is unique; in residue `p=3R+1` its lifts
form an affine line directed by the product of all repair dual lines.
Adjacent repairs have a coherently normalized transition quotient of degree
two or three, respectively. This is a **proved algebraic refinement, not
endpoint closure**. The exact open implication is exclusion of the common
completion under the residual direction, phase, and lift constraints. A
nontrivial repair-cycle identity was the proposed next attack at this stage.

A fifteenth semantic check gives Proposition 15.732.  Put
`P_A=product_(u in A)L_u` and `Theta_A=P_A^2 Phi_A`.  Every coherently
normalized adjacent transition becomes the exact potential difference
`Theta_A'-Theta_A=P_(A intersect A')^3Q_(a,z)`.  Thus every additive
closed-walk circulation, and every linear functional of it, vanishes
identically; the cubic-pencil gauge adds another vertex coboundary.  The
edge quotient does have a nonzero gauge-invariant first jet at the rich-line
dual point, but its square character is repair-coloured and is not determined
by the recorded residual parity set.  A nonrich near-pairing direction has
at most `R+2` repair tangents against envelope degree `2R+2`, and trisecant
pair masks cannot recover the full block product modulo squares.  This is a
**proved method barrier, not endpoint narrowing**. Bare cycle evaluations
and tangent-component counts are duplicate/dead routes. Proposition 15.733
subsequently uses the simultaneous p31 coefficients, and 15.734 supersedes
the first-jet front at `k=4p`.

A sixteenth semantic check gives Proposition 15.733. At `p=31,R=10`, the
fifteen phase-one mean-30 directions have one common parallel count. The
exact `b=2` and `b=30` baselines have incompatible coefficient offsets, so
all fifteen have `b=2`. Finite-edge counting and the sharp lift floor force
the opposite type to `b=0`; the one remaining hard direction would then need
the impossible odd-fibre count `b=42>30`. This is a **proved symbolic
endpoint exclusion**, not a configuration census.

A seventeenth semantic check gives Proposition 15.734 and supersedes the
whole boundary-shell campaign at the `k=4p` endpoint. A graph with
`|H|=4p+1` has an isolated vertex for `p>=13`. Signed PSL transport sends it
to infinity while preserving residual compatibility, leaving `I=0` and an
all-finite even boundary of arbitrary size. Every directional `b_d` is
therefore even. The exact hard-type budget and phase-one floors leave only
three homogeneous baseline branches, with coefficient offsets `4,5,3`.
They force `(P,s)=(4,4),(5,5),(3,3)` and an opposite phase-zero `b=0` lift
of scaled mass `8,6,8`, all below `p-3`. Hence **every `k=4p`
residual-(ii) candidate, at every boundary size, is impossible for every
prime `p>=13`**. At `p=11` the same reduction lands at the sharp equality
`8=p-3`, which was the exact small-prime frontier at this stage and is
resolved in 15.736--15.737. Proposition 15.735 separately extends the
uniform argument two layers; the global residual-(ii) predicate remains open.

An eighteenth semantic check gives Proposition 15.735. For
`k=4p+2t`, `t=1,2`, the same isolated chart exists and each type has exact
budget `2m(m+t)`. The hard baseline offsets still force
`(P,s)=(4,4),(5,5),(3,3)`. The opposite parallel-count surplus remains
strictly below the `m` opposite directions, so one direction has scaled mean
`8`, `6`, or `8`, below every admissible phase-zero lift at `p>=13`.
Therefore **all boundary sizes at `k=4p,4p+2,4p+4` are impossible for every
prime `p>=13`**. The proof stops honestly at `t=3`: branch B has surplus
exactly `m`, all opposite directions may have `Q=3` and mean `p+7`, and a
nonnegative integral quadratic realizes that local mean at `p=17`. At
`p=13,t=3` an additional exact `b=10,a=20` hard branch survives. Extending
the one-direction floor or halving calculation to another layer is therefore
a duplicate dead route unless a simultaneous invariant is added.

A nineteenth semantic check gives Proposition 15.736 and replaces the former
external-classification premise at `p=11` by an **exhaustive finite
certificate**. On all 462 points of `J(11,6)`, the pair-monomial evaluation
space has rank 55 modulo 101. A deterministic family of eight-term third
differences supplies 407 independent annihilators, so its real nullspace is
exactly that quadratic space. An exact 462-variable CP-SAT model with support
84 and no-goods for 55 omitted-pair plus 165 all-equal-triple supports is
infeasible. Thus those 220 supports are exhaustive. Proposition 15.688 is
checked explicitly at the interface: scaled mass eight is below the
`H>=2` floor twelve, forcing a Boolean lift of support 84. This kills the
hard-`b=2` branch but, by itself, leaves the simultaneous all-equal-triple
branch open.

A twentieth semantic check gives Proposition 15.737. In each of the first
three `p=11` layers, at least three exact hard `4-z_j` baselines remain. Their
signed edge coefficients are stars, so the homogeneous binary quadratic
moment
`M_H(L)=sum_e chi(e)(L(u)-L(v))^2` vanishes at at least three projective
linear forms and hence vanishes identically. The only catalog survivor in an
opposite minimum direction is an all-equal triple, whose signed coefficients
form a triangle. Its moment is
`2(r^2-r+1)`, never zero because its discriminant `-3=8` is nonsquare modulo
11. Hence **`p=11,k=44,46,48` are impossible for every boundary size**.
Together with 15.735, the first three residual shells are closed for every
prime `p>=11`; critical `p=5,7`, `p=11,k>=50`, and
`p>=13,k>=4p+6` remain open.

A twenty-first semantic check gives Proposition 15.738, an **exhaustive
finite residual-cell certificate** on `J(13,7)`.  The changed premise is the
corrected exceptional target offset two, which forces opposite mass-14 cells
at `Q=0,6`; this is not a rerun of an old endpoint census.  Proposition
15.688 reduces their maximum to one or four, and exact necessary residual
models exclude height four at both parallel counts.  Rank 78 of the pair
evaluation space and 1,638 independent third differences reduce height one
to a support-462 Boolean catalog.  An anchored no-good CP-SAT certificate
proves that the 78 selected pairs, 156 oriented mixed pairs, and 858 signed
mixed triples are exhaustive.  Their coefficient offsets are `6,4,4`, so
only a selected pair survives at `Q=0,6`.  This classifies the forced local
cell but does not close a residual row by itself.

A twenty-second semantic check gives Proposition 15.739.  In the exceptional
`p=13,t=3,u=3` row, the corrected offset restricts the seven hard
complement-triple cells to common parallel count `P=2` or `8`; exact opposite
accounting forces a Proposition 15.738 selected-pair cell.  The sign-safe
binary quartic `G=2hM_4-M_2^2` vanishes on seven hard projective directions,
hence identically, while the opposite selected pair evaluates to
`-3(i-j)^4!=0` in `F_13`.  Thus this exceptional row is closed.  The generic
row is not: for `p=1 mod 4`, `p>=17`, exact stars force even moments through
degree `(p-9)/2` and conditional cuts force alphabet `{-1,0,1,2,3}`, but no
infeasibility theorem follows.  At `p=17`, a two-intersection stabilizer
average additionally gives cut range `[-26,-12]`, or equivalently a
`{0,...,7}`-valued quadratic of total mass 8,580; this is a sharper open
classification target, not closure.  At generic `p=13`, an explicit elevated local
cell has normalized `S_2=0,S_4=5` (and global `M_4=5h`).  Another
independent-direction floor, moment, backend,
or timeout is a duplicate dead route unless it adds common-graph or genuine
cross-direction constraints.

The bounded p17 mesh audit is terminal for the present exact model.  Five
complete-domain CP-SAT encodings with all 24,310 cuts and `M_2=M_4=0`
returned `UNKNOWN` after 700--1,200 seconds without an incumbent; a lazy-cut
encoding returned `UNKNOWN` before its first separation round; five Z3
variants returned no solver status.  No witness and no infeasibility
certificate exists.  Changing seed, backend, split encoding, conditioned-set
depth, or timeout is not a new attack.  The next p17 premise must classify or
otherwise exploit the proved `{0,...,7}` value range and total mass 8,580.

A twenty-third semantic check gives Proposition 15.740 and terminates the
finite local-catalog continuation at generic `p=13,t=3`.  The hard quotient
partitions are `1^6 4`, `1^5 2 3`, and `1^4 2^3`.  In the first two, five
exact stars force both `M_2` and `M_4` to vanish.  Six cyclic distance sums
then obey exact sum, `l1`, moment, and translation-averaged cut inequalities.
Nine cut vectors eliminate all 32,313 bounded aggregate rows; an independent
14-variable CP-SAT model is infeasible.  This is a necessary-relaxation
exclusion, so it soundly removes the five- and six-exact partitions.  It
does not close `p=13,k=58`: `1^4 2^3` remains.  Do not extend this into
another local catalog, seed, backend, or orbit census.  The changed premise
required for further p13 work is simultaneous realization by one 59-edge
graph of four exact stars, three elevated hard cells, and seven opposite
cells.  The affine binary-Radon reconstruction is already Proposition
15.692 and must be imported rather than renumbered.

A twenty-fourth semantic check gives Proposition 15.741 and supplies the
required common-graph premise without extending the local catalog.  Four
exact stars force the orientation-independent tensors `M_2=T_3=0` and
`U_4=lambda M_4`, with `M_4!=0`.  The common 59-edge graph then maps to 84
nonnegative integer displacement multiplicities whose 98-by-84
difference-Radon transform has Gram matrix `13I+2J-G`, an exact inverse, and
nonstar energy `707+26C`.  Six translated cuts imply only `C<=11`.
Explicit matched local cells and a strict fractional transform witness show
respectively that the independent cellwise scalar consequences and the bare
linear difference transform do not close the row.  They do not witness one
common global quartic.  The changed premise is layered: first add integrality,
the full quartic value code, and all cut/parity lifts to the 84-class system;
then, only if needed, add midpoint placement, simplicity, and full fibre-pair
compatibility in the binary edge lift.

A twenty-fifth semantic check gives Proposition 15.742 and closes that
specific p13 gate without either lift.  Jointly imposing the already-proved
`M_2=0` congruence and the six multiplicative interval-cut inequalities on
one integral six-bin row leaves exactly 30 elevated rows of maximum energy
31 and 24 opposite rows of maximum energy 82.  Direct bounded enumeration,
an independent slack-coordinate reconstruction, and independent one-worker
CP-SAT exclusions agree.  Thus the ten nonexact rows have energy at most
`3*31+7*82=667`, contradicting the common-graph identity
`707+26C>=707`.  This closes the generic `p=13,t=3` partition and, with
15.739, all of `p=13,k=58`.  It does not close later p13 layers, generic p17
and above, residual (ii), Type I, or the limit.

A twenty-sixth semantic check gives Proposition 15.743 and closes the first
generic resonance at `p=17,k=74`.  On an exact `k=1` hard star, equality of
the local and common row sums gives `hT=18P-69`.  Since `hT` is common, the
exact stars share `P`; there are at least six, so `6P<=75` gives `P<=12`, and
the literal congruence `P≡5 (mod 8)` forces `P=5`.  Thus `hT=21` and the exact
coefficient row is `q=(2)^8`.  For any hard quotient `k`, the common
difference-Radon row has off-bin sum `21-P`, whereas the directional cell
identity gives `17(P-3)-18k`; hence `P=4+k`.  Thus the three quotient
partitions `1^8 4`, `1^7 2 3`, and
`1^6 2^3` have fixed parallel profiles, total hard parallel count 48, and
nine opposite `Q=3` rows.

All 24,310 nine-sets give exactly 698 distinct translated cyclic-distance cut
vectors.  Deterministic one-worker CP-SAT models use broad coordinate domains
from the exact `l1` bounds and impose only the row sum, `l1`,
`M_2=M_4=0`, and all 698 cut inequalities.  With no prior energy upper bound,
they exclude excess one outright, excess-two energy at least 71,
excess-three energy at least 120, and opposite energy at least 73.  The
opposite coordinates have fixed sum `-24`, so equality in Cauchy makes
`(-3)^8` the unique energy-72 row.  The two partitions containing excess one
are already impossible.  For `1^8 4`, the nonstar upper bound is `119+9*72=767`,
contradicting its exact common-Radon energy `1211+34C>=1211`; the other two
nonstar baselines are 1251 and 1287.
This closes exactly `p=17,k=74`, not any `p>=17,t>=4` layer (beginning
with `p=17,k>=76`). The branch-B `t=3` range formerly left here is closed
later by 15.751. Residual (ii) and the limit remain open; Type I is closed by
15.750.

A twenty-seventh semantic check gives Proposition 15.744.  The complete
`p=13,t=4` residue sieve leaves `u in {0,3,4,6}`.  In `u=3`, the hard
quotient profile is six exact complement triples and one elevated row.  The
phase-one `b=10` equality is upgraded from its averaged quadrature by a
rank-78 restriction to the three contact layers.  Its two-unit alternative
is not treated as a globally nonnegative lift: a separate 1,716-variable
model with 1,638 third-difference identities, the exact punctured lower
bounds, and `sum B=66` is infeasible.  The
parallel congruence and common signed-total equations leave `P=2,8`; either
ledger forces an opposite phase-zero mass-14 cell at `Q=6` or `Q=0`.  The
old Proposition 15.738 height-four infeasibility is not reused: its
`|H|=59` l1 bound is changed here.  Both necessary models are rebuilt at
`|H|=61`, with `sum|W|<=61-Q`, and exact one-worker CP-SAT makes both
infeasible.  Only then is 15.738's edge-independent Boolean support catalog
imported, leaving a selected pair.  Six exact hard roots make
`G=2hM_4-M_2^2` identically zero, while the opposite selected pair gives
`-3(i-j)^4!=0`.  This closes exactly `p=13,t=4,u=3`, not the full `k=60`
row.

A twenty-eighth semantic check gives Proposition 15.745 and closes
`p=13,t=4,u=0` without a displacement census.  Exact-star/common-row gluing
forces `hT=17` and `P_L=4+k_L`; the opposite parallel profile is `3^6,4`.
All 74 translated cuts under the forced moment congruences make an opposite
`Q=3` row infeasible when at least five stars are exact, and give the broad
relaxation maxima 31, 96, 76, and 111 needed for the two remaining
partitions.  Three partitions die rowwise; `(1,1,1,1)` dies by
`691<721+26C`.  For `(2,1,1)`, the row bound gives `C<=1`, while its seven
parallel edges in six classes give `C>=1`.  The equality case puts the sole
duplicate in the elevated zero bin.  A direct audit of all 84 nonzero
functional bins in `F_13^2` proves the transverse sign split `6+/7-`, hence
`q_a in [-7,6]`; the elevated maximum becomes 66 and
`695<719`.  Together, 15.744--15.745 leave exactly `u in {4,6}` at
`p=13,k=60`.  Residual (ii), Type I, and the limit remain open.

A twenty-ninth semantic check gives Proposition 15.746.  The all-positive
`b=2` quadrature is used pointwise before Proposition 15.688, making each
hard `u=4` lift Boolean of support 330 on `J(13,7)`.  The exact classifier
uses 1,716 Boolean variables, all 1,638 third-difference identities, the
support and anchor equations, and all 70 anchored no-goods: 1,710 constraints
in total.  Exact infeasibility, independently replayed with one worker,
proves that the complete catalog is the 78 omitted-pair and 286 all-equal
triple supports.  Their full and anchored SHA-256 digests are
`4edf1fe1b9c73f05598b667dba121f064807c68421a4df2c8db7090a3e3ff35f`
and
`84ce6099dcca66f7cc2792dc60bcbb378672f2e9cac2b19e02812f2f20563c7a`.

The catalog offsets and common signed total forbid mixing: omitted pairs
force `P=3`, all-equal triples force `P=5`, and `hT=14P-61`.  The opposite
excesses sum to five, so at least two opposite cells have mean 12.  At
`P=3,Q=5`, the literal is excluded modulo six; the cell is a `b=0` mass-12
lift of height one/support 396 or height four.  The seven hard omitted-pair
directions force the sextic `F6=2hM6+hM2^3-3M2M4` identically zero, after
checking all `78^2=6,084` baseline/lift overlaps.  Thus each forced
`P=3,Q=5` cell satisfies `F6=0`.  The `P=5,Q=3` branch retains the
literal-or-lift dichotomy.  Across its 22,308 pair/triple patterns, the
weighted even-moment feature ranks are full through degree six, so there is
no analogous universal polynomial identity in `N2,N4,N6` at those degrees;
this does not rule out a different invariant.
Under the opposite-cell normalization `N'_(2r)=(-h)M_(2r)`, the constraint
is `2N'_6+(N'_2)^3+3N'_2*N'_4=0`; copying the hard-sign formula would give
the wrong mixed-term sign.
This is an exhaustive finite equality classification and proved open
reduction, not a close; `u in {4,6}` remains exact.

A thirtieth semantic check gives Proposition 15.747. For every integral
parallel count, a Boolean mass-12 lift would satisfy
`-7D2+84E2+182Q^2-1428Q+2598=0`; its residual is one modulo seven, so the
Boolean case is impossible. Two exact 169-variable, 3,526-constraint
necessary-relaxation models also exclude height four at `Q=3,5`. This closes
the omitted-pair `P=3` branch and forces every minimum all-equal-triple
`P=5,Q=3` cell to be the exact literal. The field sextic is not needed.

A thirty-first semantic check gives Proposition 15.748. The forced literals
are common roots of `M2,M4,M6`. Root count excludes at least five literals;
exact interpolation against the 69-element hard moment alphabet excludes
four and three for both hard signs. With two roots, 1,554 `(M2,M4)`
candidates and 2,688 allowed `N6` vectors per sign leave exactly 336
moment-level survivors per sign. Hence only excess partition `(1^5)`
remains. These records are necessary moment data, not common 61-edge graph
constructions; Proposition 15.749 closes them without graph reconstruction.

A thirty-second semantic check gives Proposition 15.749.  For each of the
five surviving opposite `Q=4` cells, all 74 translated cuts imply
`c.q<=-78`, alongside `sum q=-13` and `sum |q|<=57`.  Exact rational dual
identities for a coordinate and its negative give
`-52/9<=q_a<=26/15`; multiplicative symmetry and integrality reduce every
coordinate to `[-5,1]`.  Deterministic enumeration gives 522 rows and 492
distinct `(N2,N4,N6)` triples.  Their intersection with the 48-element
nonroot evaluation alphabet of either hard sign has 12 triples, all with
`N4=0`.  Thus the five `Q=4` directions and the two literal directions are
seven roots of the binary quartic `M4`, forcing `M4=0`; this contradicts the
nonzero hard fourth-moment alphabet.  Hence `p=13,t=4,u=4` is closed and the
exact `p=13,k=60` remainder is `u=6`.  This is an exhaustive finite
aggregate certificate and proved branch theorem, not a global residual-(ii)
close.

The Bartoli--Storme ``unique-trisecant ceiling'' previously recorded after
15.729 is **RETRACTED**. Under the corollary's other hypotheses, including
`d>3+2sqrt(q)` and existence of the configuration, its threshold is the upper
endpoint for which the associated hyperplane arrangement is second-smallest;
it is not a nonexistence theorem for larger unique-trisecant 3-arcs.
Likewise, Ball--Lavrauw's size hypothesis
is used only for their explicit interpolation formula. Proposition 15.731's
below-threshold envelope is derived from their scaled tangent lemma and the
new line-gluing argument, not by applying that theorem outside its scope.

The independent `p=31` complete-22-arc check is a finite historical side
route only.
Eleven of the twelve classified classes have sourced representatives and
exact audited `c_1<=2<10`; the twelfth representative remains unavailable.
Propositions 15.733--15.734 close the endpoint without that missing
representative, so completing the class list cannot advance a live gate. The
exact replay is
`scripts/p31_complete_22arc_public_audit.py`, with scope ledger
`evidence/NOTE_2026-08-31_p31_public_11_of_12_arc_audit.md` and deterministic
certificate `evidence/p31_complete_22arc_public_11_audit.json`.

The same replay exposed one stale exact-boundary diagnostic:
`p17_slack20_boundary_cryptominisat.py` still expected the 78 profiles and
69 signatures produced by the retracted blanket filter.  It now consumes the
corrected 193-profile block at census indices 1364--1556, deduplicates it to
184 signatures, and recomputes the full reflection ledger.  Proposition
15.707's algebraic exclusion already removes all 193 rows; the solver remains
an optional independent audit, not a proof dependency.

1. The final 300-second positive-`p=7,z=7` CP-SAT run repeated an existing
   exact full-torsion model for the same case.  Only the timeout changed.
2. Several long proposition blocks are different coordinates for the same
   unresolved scalar or relaxation.  In particular, most of Props.
   15.83--15.160 are the optional Path-C/Hypothesis-H residual in different
   forms, and much of 15.321--15.560 is the same unnamed `Q_tau`/class-function
   mixture under successive small-prime fits.
3. The old bounded acceptance AND is now exposed only as
   `e1_bounded_residual_split_closed()`. The corrected
   `e1_closed_general()` is the global gate and returns `False`, matching the
   still-open residual-(ii) and Type-I units. Historical prose that says the
   old AND is `True` does not close the current gate.
4. The first audit restored **GLOBAL QVAR** to what it then treated as a
   spectral-floor acceptance unit. The 15.720 correction above supersedes
   that conclusion: the whole spectral unit is no longer load-bearing.

No new computation should be launched from an attractive formula or script
name until this file is checked first.

## Coverage and numbering

- Propositions 15.1--15.82 are written directly in `solution.md` and related
  early modules.
- There are 665 source-backed proposition modules from 15.83 through 15.750.
- The labels 15.537, 15.583, and 15.584 have no proposition module.  They are
  unassigned labels, not unreviewed propositions; later source headers mention
  those numbers only as historical range/state markers.
- Therefore every assigned proposition through 15.750 was included in this
  audit.  The grouped ledger below is by shared mathematical route rather than
  a 736-row restatement of the assigned propositions.

## Authoritative acceptance chain

The public theorem is gated consistently by the corrected global
`e1_closed_general()` Boolean and `four_e1_units_closed()["closed"]`; both are
currently `False`. The historical bounded `True` is available only through
`e1_bounded_residual_split_closed()` and is not a global theorem predicate.

| unit | exact live content | status after audit of 15.751 |
|---|---|---|
| required bi-tight levels 2 and 3 | 15.720 degree congruence using 15.272/15.207 | **TRUE** |
| residual (ii) | non-Walsh multi-level Max-minus for every even `k>=4p` | **OPEN** — 15.734--15.737 close the first three shells; 15.751 closes `k=4p+6` for every `p>=13`; 15.744--15.749 close `u=0,3,4` at `p=13,k=60`. The exact remainder includes critical `p=5,7`, `p=11,k>=50`, `p=13,k=60,u=6` and later p13 layers, every `p>=17,t>=4` layer (beginning with `p=17,k>=76`), and positive `p=7,z=7`. The global predicate stays false. |
| Type I | the multi-level `3A+B>0` bad case | **TRUE** (15.750) — isolated-chart rigidity and parity halving close every prime `p>=11`; exact integer Farkas identities close `p=5,7` |
| Lemma D | every good-line triple and its Fejer two-plane amplitudes | **TRUE** (15.276) |

The spectral floor remains an interesting optional problem, but it has no
valid downstream role in the current E(1) proof. The shortest honest work map
is now residual (ii), the sole false E(1) unit, followed by the final
implication audit.

The positive `p=7,z=7` catalog is one finite residual subbranch, not a fourth
top-level front.

## Complete proposition-range account

| propositions | durable content | effect on the current gates / duplication rule |
|---|---|---|
| **15.1--15.19** | conference spectral calculus, switching, cube moments, the sandwich, and finite small-`n` results | The global fourth-moment/spectral-defect shell is asymptotically vacuous (15.16, 15.19).  Do not restart it. |
| **15.20--15.82** | Hamming/Max-Lipschitz reductions, matching and cover structure, finite `n=6,7,8,10` results, and early `m4`/resolvent forms | No live all-prime gate closes.  Perfect-matching exhaustion, continuous-Gamma/SDP transfer, generic Gaussian domination, structure-free projector bounds, and type6/cross-ratio pinning are already insufficient or false. |
| **15.83--15.160** | Path-C/Hypothesis-H/`16N`, `delta`/ED4/FFT/`R4`/cumulant/Gegenbauer dictionaries, and finite `p=5,7` certificates | This is mostly one optional residual written in many coordinates.  H is not an independent proof of its own residual (15.90).  Path C is not a current acceptance gate.  Class-key, raw PGL, one-line Aut, fixed type lists, Delsarte/moment LP, ULC, Jensen, pole, and Chebyshev routes are quarantined. |
| **15.161--15.240** | conditional spectral majorization; exact Gsum/kernel/`mu`/`delta` identities; old bounded residual-(ii) close | 15.179 and 15.236--15.237 close affine and even `k<=4p-2`, not the live `k>=4p` range.  The old Gsum scale, affine exhaustiveness, unsigned permanent, and the much looser 15.217 `delta` room must not be relabelled as R1. |
| **15.241--15.272** | exact residual-(i) reductions and the `k=1 union k=3` Veronese span | **15.272 genuinely closes only the two-level Type-I/residual-(i) slice.**  Aut-Schur and `k=3`-only span are false; Gsum and the cotangent pairing are unused. |
| **15.273--15.320** | Lemma D; `Aut dot F=Z`; character-pair, Gauss, Jacobi, torus, Kloosterman, and `Q_tau` floor reductions | **15.276 closes Lemma D.**  15.278 reduces the spectrum to `F`, but no proposition proves the floor.  Small-`p` orbit formulas, familywise floors, coarse Q-types, AP/QR0 generation, and interpolation are not `p`-laws. |
| **15.321--15.400** | increasingly refined `Q_tau`, class-function, occupancy, Jacobi, circle, LP, PSD, and floor models | Every floor statement remains open.  Two-point fits, low-degree names, occupancy LPs, Cauchy--Schwarz, PSD, and pointwise floor arguments are already recorded as insufficient or false. |
| **15.401--15.480** | further `Q_tau`/nonlinear-orbit names and finite `p=5,7` Type-I diagnostics | Almost all claims are finite-prime identities or killed extrapolations.  Aut-orbit-size guesses, Gauss/Jacobi/CM interpolation, type-count extrapolation, and character kernels do not name the full mixture. |
| **15.481--15.560** | more `A_full/Q_tau` reductions plus finite `p=5` residual slices | The finite `nF` exclusions do not close general residual (ii) or Type I.  Type-index Gram, one-dimensional Johnson, Max-minus Fourier support, and Aut-e inversion are insufficient or dead. |
| **15.561--15.589** | final class-function no-gos, exact profile classification, and the QVAR decomposition | 15.589 closes QVAR only for `k=1,...,6`. Per-stratum `k>=7` is false and is **not** the leftover. GLOBAL QVAR is now optional. Separately, 15.585 A relied on one-sided level-4 tight emptiness and is retracted; only its `{2,4,6}` mass calculation survives. |
| **15.590--15.628** | degree-four SoS countermechanism, exact R1/`delta` hierarchy, Walsh/W1/W2 investigation | Degree-four SoS cannot force Type I (15.590).  Character/PSD-only and fixed-channel R1 routes are insufficient.  **15.628 closes Walsh, W1, and W2 for all odd primes**, but explicitly leaves the non-Walsh 5+-level/even-`k>=4p` branch. |
| **15.629--15.668** | complete low R1 shells, modular-data no-go, nonlinear shell positivity, finite boundary closures, and exact `p=11` theta/channel work | Strong R1 is true at `p=11` by full census.  Scalar trace and broad square-circle conserved-mass cones through exponent 800 still admit sub-six targets and cannot prove general R1.  Props. 15.643--15.666 close the two-point, size-four, size-six, and finite `p=7` size-eight residual branches; rerunning their old solvers is duplication. |
| **15.669--15.712** | uniform residual ranges, infinity-plus-`(p-2)`, all-finite endpoint campaigns, and optional no-gos | After 15.721, the all-finite boundary-close role of 15.675--15.712 is superseded by signed transport into 15.669/15.674. Do not regenerate any first/second all-finite rows, including the former open `p=23` ledger. 15.676 remains load-bearing at total boundary `p+1`; 15.690--15.691 and reusable internal lemmas retain their independent content. |
| **15.713--15.719** | positive `p=7` infinity-plus-seven reductions | 15.713--15.717 close `z=0,1,2,3`.  15.718--15.719 identify and stabilize projected `z=7` semigroup supports but remove no source boundary.  All 56 actual `z=7` line boundaries remain open, and the semigroup/quotient route is terminated. |
| **15.720--15.721** | degree-congruence bi-tight close; signed boundary normalization | 15.720 closes the required bi-tight levels. 15.721 proves `|D|>=p+1` for every residual candidate at `p>=17` and identifies strict deficit in the normalized infinity-plus-`p` shell as the first general residual branch. Neither closes Type I or residual (ii). |
| **15.722--15.724** | exact phase cocycle; outside-pair slack; paired-cube floor-plus-two repair; full-circle exclusion | 15.722 identifies slack zero with an aligned Miquelian circle, excludes `R=1,2,3`, and more generally excludes `1<=R<=floor(sqrt(p)-5/2)` by minimal arc deletion plus the prime-field conic threshold. 15.724 excludes the circle; at that stage the `p+1` branch lay beyond `max(3,floor(sqrt(p)-5/2))`. 15.723 proves the middle floor-plus-two shortcut except for the explicit cells `(17,5,1)` and `(17,11,0)`, which every later profile audit must retain. |
| **15.725** | finite parabola-plus-internal census and attempted all-prime character bound | **RETRACTED as a family close.** The finite phase-zero census is retained; the all-prime character sums and opposite orientation are open. It changes no gate. |
| **15.726** | minimal arc deletion plus the Ball--Lavrauw dual tangent envelope | **PROVED historical narrowing.** For every prime `p>=17`, it excludes `1<=R<=floor((p-4)/3)` at `|D|=p+1`; 15.734 now supersedes the whole shell as a live endpoint gate. |
| **15.727** | endpoint tangent-envelope equality, disjoint rich blocks, and published arc classifications | **PROVED historical narrowing.** Equality forces disjoint trisecant/4-secant blocks and was excluded at `p=17,19,23,29`; 15.734 now supersedes every remaining endpoint case. |
| **15.728** | exact `p=31` odd-fibre sum, Paley-type residue budget, and nonrich `b=2` directions | **PROVED necessary normal form.** It fed the direct 15.733 close and is superseded as a live gate by 15.734. |
| **15.729** | retain-three/retain-two block deletion and co-tangent extension reduction | **PROVED all-prime necessary reduction, not endpoint closure.** Every remaining endpoint gives an affine `(p+2-R,3)`-arc with exactly one trisecant and a `(p-R)`-arc with two extensions on one tangent; 15.730--15.732 subsequently sharpen this to the full repair ensemble, tangent transitions, and exact-cycle barrier. |
| **15.730** | all `3^x6^y` maximum repairs, complementary arcs, and exact two-colour projective/directional census | **PROVED historical endpoint structure.** Every repair has an `R`-arc index-one complement; 4-secants give three displayed co-tangent extensions. Superseded as a live gate by 15.734. |
| **15.731** | direct tangent-envelope gluing and adjacent-repair transition quotient | **PROVED historical endpoint structure.** The normalized envelope is unique or a line-product pencil, with quadratic/cubic swap data. Superseded as a live gate by 15.734. |
| **15.732** | exact cleared repair cycles, rich-direction first jet, and natural phase-bridge audit | **PROVED method barrier.** Every cleared circulation is a coboundary and the natural phase bridges fail. Its algebra remains valid, but 15.734 makes a phase bridge unnecessary at `k=4p`. |
| **15.733** | simultaneous exact `p=31` baseline coefficients | **PROVED symbolic endpoint exclusion.** Closes `p=31,R=10` without a finite configuration census; subsumed by 15.734. |
| **15.734** | isolated-chart hard-type baselines at `|H|=4p+1` | **PROVED endpoint theorem.** Every boundary size is impossible for every prime `p>=13`; its p11 sharp-equality frontier is subsequently resolved by 15.736--15.737, and its larger-layer limitation is extended by 15.735. |
| **15.735** | isolated-chart type budgets and exact surplus arithmetic at `k=4p+2,4p+4` | **PROVED two-layer extension.** With 15.734, closes the first three residual shells for every `p>=13` and every boundary size. The `t=3` branch-B surplus equals `m`, so no later layer is claimed. |
| **15.736** | exact `J(11,6)` quadratic evaluation space and sharp support-84 Boolean catalog | **EXHAUSTIVE FINITE CERTIFICATE.** Rank 55 plus 407 independent third-difference identities make the 220 omitted-pair/all-equal-triple supports exhaustive. The hard-`b=2` p11 branch closes; simultaneous all-equal triples remain for 15.737. |
| **15.737** | signed star/triangle coefficient patterns and the binary quadratic moment over `F_11` | **PROVED p11 three-layer theorem.** Closes `k=44,46,48` for every boundary size. Together with 15.735, the first three shells are closed for every prime `p>=11`. |
| **15.738** | exact `J(13,7)` mass-14 height-four models and support-462 Boolean catalog | **EXHAUSTIVE FINITE CERTIFICATE.** Height four is infeasible at `Q=0,6`; the 1,092 Boolean supports are exhausted, and offsets leave only a selected pair. This is a local-cell theorem, not residual closure by itself. |
| **15.739** | corrected exceptional p13 ledger, sign-safe quartic, and generic even-moment/cut reduction | **PROVED BRANCH THEOREM AND HISTORICAL OPEN REDUCTION.** The exceptional `p=13,t=3,u=3` row is empty. Proposition 15.742 closes its generic p13 complement, 15.743 completes p17, and 15.751 closes the former generic `p>=29,t=3` remainder. |
| **15.740** | cyclic distance aggregation and translation-summed cuts in the generic p13 row | **PROVED BRANCH SPLIT WITH EXHAUSTIVE FINITE CERTIFICATE.** Nine inequalities eliminate all 32,313 aggregates when at least five hard stars are exact. Only `1^4 2^3` remains, subsequently closed by the common-graph energy certificate in 15.742. |
| **15.741** | common-graph cubic/quartic endpoint tensors and the 84-class difference-Radon transform | **PROVED OPEN REDUCTION AND METHOD BARRIER.** The four exact stars force `M_2=T_3=0`, `U_4=lambda M_4`, and `M_4!=0`; the transform gives an exact nonnegative-integer inverse and energy `707+26C` with `C<=11`. Its `M_2` and exact energy identities are the inputs completed by 15.742. |
| **15.742** | joint `M_2` congruence, six interval-dilate cuts, and common-graph energy | **EXHAUSTIVE FINITE CERTIFICATE.** The sharp elevated/opposite row energies are 31 and 82, so the ten nonexact rows have energy at most 667, contradicting `707+26C>=707`. The generic four-exact p13 branch and, with 15.739, all of `p=13,k=58` are closed. |
| **15.743** | cross-direction `P=4+k`, all 698 translated-cut vectors, and the p17 common-Radon energy | **EXHAUSTIVE FINITE CERTIFICATE.** Broad-domain one-worker CP-SAT excludes excess one, excess-two energy at least 71, excess-three energy at least 120, and opposite energy at least 73 without a prior energy cap; fixed sum `-24` then makes `(-3)^8` the unique opposite row of energy 72.  The only quotient partition not killed rowwise has `767<1211<=1211+34C`, closing `p=17,k=74`. |
| **15.744** | full p13 t4 residue sieve, b10 contact-layer and punctured-lift certificates, changed H61 mass-14 models, and six-root quartic | **PROVED BRANCH THEOREM WITH EXHAUSTIVE LOCAL CERTIFICATES.** Rank 78 makes the exact b10 cell pointwise and the separate punctured model excludes its two-unit alternative; the H59 height-four result is not imported, and both H61 models are rebuilt and infeasible.  The selected-pair survivor contradicts the identically zero quartic, closing exactly `p=13,t=4,u=3`. |
| **15.745** | full 74-cut p13 rows, exact Radon partition energies, and collision-one transverse sign audit | **EXHAUSTIVE FINITE AGGREGATE CERTIFICATE.** Three partitions die rowwise, one by `691<721+26C`, and the equality case of the last gives the sharpened `695<719` contradiction.  With 15.744, the exact `p=13,k=60` remainder is `u in {4,6}`. |
| **15.746** | exact support-330 Boolean catalog, branchwise `u=4` edge ledger, and omitted-pair sextic | **EXHAUSTIVE FINITE EQUALITY CLASSIFICATION AND PROVED OPEN REDUCTION.** Exact infeasibility makes the 78 omitted pairs and 286 all-equal triples exhaustive.  Their offsets force uniform `P=3` or `P=5` and at least two opposite mean-12 cells; at `P=3,Q=5` the cell is a `b=0` mass-12 lift satisfying `F6=2hM6+hM2^3-3M2M4=0`.  Neither `u=4` family is closed, so `u in {4,6}` remains exact. |
| **15.747** | general mass-12 Boolean cut second moment and exact height-four models at `Q=3,5` | **PROVED BRANCH EXCLUSION WITH EXHAUSTIVE FINITE CERTIFICATES.** The Boolean branch is impossible modulo seven and both height-four relaxations are infeasible. The `u=4` omitted-pair `P=3` branch is closed; minimum `P=5,Q=3` cells are forced literals. |
| **15.748** | common literal-root interpolation against the 69-element hard moment alphabet | **EXHAUSTIVE FINITE INTERPOLATION CERTIFICATE AND PROVED OPEN REDUCTION.** Cases with at least three literals are empty for both hard signs. Exactly 336 moment-level two-root survivors per sign remain, forcing opposite excess partition `(1^5)` without constructing a common graph. |
| **15.749** | all 74 translated-cut inequalities, exact coordinate duals, and the `Q=4` moment intersection | **EXHAUSTIVE FINITE AGGREGATE CERTIFICATE AND PROVED BRANCH THEOREM.** The exact box contains 522 rows and 492 moment triples. Its 12-point intersection with each survivor evaluation alphabet has zero fourth moment, so five `Q=4` roots plus two literal roots force the nonzero hard quartic to vanish. This closes exactly `p=13,t=4,u=4`. |
| **15.750** | isolated-chart rigidity, parity halving, and exact small-prime Farkas bases | **PROVED ALL-PRIME THEOREM.** Closes multi-level Type I for every prime `p>=5`. |
| **15.751** | half-mean cube height theorem, corrected transposition influences, slice-to-cube symmetrization, and fixed four-bit catalog | **PROVED INFINITE-FAMILY THEOREM WITH FIXED EXHAUSTIVE CERTIFICATE.** Closes generic branch B for `p>=29` and hence `k=4p+6` for every prime `p>=13`. The global residual predicate stays false. |

## Exact duplicated run

The files

- `/tmp/p7_z7_compact_cpsat_real_smoke.json`, and
- `/tmp/p7_z7_direct_semigroup_case0_300s.json`

use the same case `orbit0_leaf780_branchA`, the same target hash
`e7deeecdfad87ce61615d7e86ff4ef247c59c8be3e4ab27417c7710cd20ff3f1`,
and the same full `F_3^6 x F_7^21` model:

- 280 `L` variables and 27 quotient variables (307 total);
- 112 Johnson-kernel, 8 mass, and 27 modular equations
  (147 constraints total).

The first run returned `UNKNOWN` after 0.2 seconds.  The later run returned
`UNKNOWN` after 300.195 solver seconds, 95,634 conflicts, and 1,179,303
branches.  This was a timeout extension, not a new formulation or theorem.

## Route blacklist

Do not reopen the following without a genuinely new mathematical input that
is absent from the cited proposition chain.

### Spectral/R1

- global fourth-moment defect or spectral-shell separation (15.16, 15.19);
- Path-C/H/`16N` under another equivalent scalar name (15.83--15.160);
- character, representation, trace, or autocorrelation PSD alone (15.690);
- scalar or broad-channel `p=11` conserved-mass LPs (15.641, 15.667--15.668);
- per-stratum `k>=7` QVAR, pointwise QVAR, or separate profile bounds
  (15.589 and the explicit counterexamples);
- another small-prime Jacobi/Gauss/CM interpolation of `Q_tau`.

Any future optional R1 proof must use information absent from the abstract
spectra: the Boolean rank-one identity and the exact full Max-plus orbit
mixture. Any optional QVAR proof must couple the mixed-`k` ensemble. Neither
is a current E(1) gate.

### Type I and residual

- Gsum disjoint lower bound, Aut-Schur, `k=3`-only span, or the cotangent
  pairing;
- affine/two-level residual work already closed by 15.179, 15.236--15.237,
  and 15.272;
- generic one-sided level-4 tight-cover emptiness: false by the explicit
  four-line family in 15.402. The family itself is harmless by the 15.272 G /
  15.588 E cylinder witness; only residual-compatible one-sided tightness is
  a live target;
- Walsh/W1/W2 search after 15.628;
- size-four, size-six, finite `p=7` size-eight, `p=17`, or `p=19` solver
  reruns after their exact closures;
- any all-finite first/second-shell profile work from 15.675--15.712:
  15.721 transports those sizes into the already-closed infinity ranges;
- any boundary-size, outside-slack, endpoint-repair, tangent-envelope, or
  phase-bridge campaign in the first three residual shells for `p>=11`:
  Propositions 15.734--15.737 close every boundary at
  `k in {4p,4p+2,4p+4}`. The p31 class list, repair graph, tangent
  transitions, p11 sharp catalog, and binary moment are completed inputs,
  not new search invitations. At `t=3`, a changed premise must couple
  directions or add a global invariant; the one-direction mean/floor route
  is blocked by an explicit local lift;
- treating a projected/parity/semigroup survivor as a feasible graph;
- trying to reach the general residual with an `L2` bound on `delta`
  (15.595 proves the scale loses from `p>=11`).

### Positive `p=7,z=7`

- another projection dimension, seed, backend, encoding, or timeout for the
  Johnson-semigroup/quotient route;
- completing the cancelled odd `k=6` shard merely for coverage;
- interpreting target presence as capped lift or binary-edge feasibility.

The route is terminal unless a new separating invariant is proved.

## Working-tree quarantine

No file is deleted by this audit.  The following existing untracked groups
are retained as history/data but are removed from the active attack queue.

| files/globs | disposition |
|---|---|
| `scripts/w1_*`, `scripts/w2_*`, `scripts/walsh_*` and matching evidence | superseded by the general Walsh close, 15.628 |
| `scripts/p5_*size_four*`, `scripts/residual_size_four_*` | size-four branch closed by 15.652--15.656 |
| `scripts/p5_*size_six*`, `scripts/p7_size6_*`, `scripts/residual_boundary_size_six_*` | size-six branches closed by 15.657--15.661 |
| finite `p=7` size-eight/fixed-boundary scripts | finite branch closed by 15.662--15.666 |
| `scripts/p17_*`, old `p19`, `p23`, and `p31` endpoint diagnostics | superseded as first-three-shell boundary gates by 15.734--15.737; retain only as historical certificates or reusable finite-geometry work |
| `scripts/r1_p11_*`, broad theta/channel scripts | useful finite data, but the aggregate cones are certified insufficient and `p=11` itself is already R1-positive |
| `scripts/p7_infinity7_positive_z7_*` and matching `/tmp` artifacts | terminal semigroup/quotient campaign; not active |
| generic `residual_boundary_*` parity/projected models | necessary-condition diagnostics only; infeasibility can close a specified finite branch, feasibility cannot close or witness the graph problem |
| `scripts/infinity_plus_p_quantized_dp.py` | diagnostic only until every `excess != 2` use retains the two Proposition 15.723 `p=17` equality cells |

This quarantine was materialized on 2026-08-31.  The useful truth-state
corrections and twelve modules required by tracked certificate generators
were retained on `main`; the superseded snapshot was moved to the remote branch
`archive/dirty-worktree-2026-08-31`.  See
`NOTE_2026-08-31_DIRTY_WORKTREE_TRIAGE.md` for the exact disposition.  The
archive branch is historical, not an active attack queue.

## Stale or contradictory records

Central live documentation and the canonical JSON summaries are corrected
alongside this audit.  Former payloads are preserved under explicit
`historical_retracted` / `historical_pre_15723` names rather than left at the
canonical paths.  Treat the following as quarantined:

- The former 15.678 and 15.684 payloads are preserved as
  `e1_gmin_m4_prop15678.historical_retracted.json` and
  `e1_gmin_m4_prop15684.historical_retracted.json`.  They record the false
  endpoint closes and are not theorem evidence.
- The former bulky 15.700--15.712 payloads are preserved as
  `e1_gmin_m4_prop157NN.historical_pre_15723.json`.  Their counts descend from
  the retracted blanket floor-plus-two filter.  The canonical JSON files now
  carry the corrected replay summaries and point back to those historical
  payloads.

- `e1_closed_general()` is now the corrected global predicate and is `False`,
  in parity with `four_e1_units_closed()["closed"]`. The old bounded `True`
  survives only as `e1_bounded_residual_split_closed()`; never use that alias
  as an E(1) or global residual-(ii) close.
- Prop. 15.104's title claims a general `16N` proof, while its own final
  packaging retracts the proposed proof.  It has neither focused test nor
  evidence JSON.
- Prop. 15.45's displayed general average uses the `p=5` value `-1/15`;
  15.47 contains the corrected general threshold
  `-(p-2)/(p(2p-1))`.
- Prop. 15.69's heading overstates the all-prime scope of
  `lambda_max(T)=4p`; the body only certifies the relevant finite primes.
- Props. 15.72 and 15.73 record conflicting non-closing `p=7` type6 probe
  values.  Neither is an active input.
- Prop. 15.74's `p=7` table value conflicts with its formula; the corrected
  budget is `3/68`, and the finite conclusion remains below that budget.
- Historical 15.167--15.171, 15.272, `HINGE_GRAPH_15272.md`, and old evidence
  JSON may say residual (ii), E(1), or `L` is closed.  Those statements use
  obsolete scope or a false spectral-floor premise.

## Artifact coverage gaps

These gaps are inventory facts, not requests to regenerate data:

- no focused test: 15.83, 15.84, 15.85, 15.104, 15.202, 15.208;
- no canonical `e1_gmin_m4_prop*.json`: 15.104, 15.208, 15.278, 15.718,
  15.719.

Props. 15.718--15.719 instead have focused tests, prose certificates, and
hashed external artifacts.  Missing canonical JSON alone is not a reason to
rerun them.

## Mandatory preflight for future attacks

Before spending mesh/GPU time:

1. state which one of the two live fronts the attack can close;
2. name the exact proposition whose limitation it overcomes;
3. in any sparse layer, first audit the isolated-vertex count and transport an
   isolated vertex to infinity; Propositions 15.734--15.737 already close
   every boundary in the first three shells for `p>=11`;
4. do not reopen `k in {4p,4p+2,4p+4}` for `p>=11`, `k=4p+6` for
   `p>=13`, `p=13,k=58`,
   `p=17,k=74`, or the `p=13,k=60` residues `u=0,3,4`: Proposition 15.742 closes the p13 t3 row by the exact
   six-dilate/common-energy contradiction, and Proposition 15.743 closes the
   p17 row by the 698-vector/common-Radon certificate. Propositions
   15.744--15.745 close the two p13 t4 residues after, respectively,
   rebuilding the H61 mass-14 premise and auditing the collision-one
   transverse signs, and 15.749 closes `u=4` by its exact translated-cut
   moment intersection. The live residual ranges are critical `p=5,7`,
   `p=11,k>=50`, `p=13,k=60,u=6` and later p13 layers, every
   `p>=17,t>=4` layer (beginning with `p=17,k>=76`). Proposition 15.751
   closes the former generic branch-B `t=3` range; no further cell, graph,
   slice, or prime census there is authorized. Retain both 15.723
   floor-plus-two exceptions in any independent profile DP. Do not reopen the
   exceptional `p=13,t=3,u=3` row: 15.738--15.739 close it by an exact
   mass-14 catalog and quartic contradiction. Do not reopen the generic p13
   `t=3` aggregate: 15.740 leaves exactly `1^4 2^3`, 15.741 supplies its
   common-graph energy identity, and 15.742 proves the resulting integral
   row system empty. Propositions 15.746--15.749 complete the sharp mass-10
   classification and close both `u=4` branches. Do not rerun them or couple
   the 336 moment survivors to graphs. The next finite p13 residue is `u=6`,
   not a broad census or a rerun of `u=0,3,4`;
5. search tracked files, untracked scripts, `/tmp` artifact names and hashes,
   git history, GitHub, MathOverflow, literature notes, and OEIS when number
   patterns are involved;
6. write the acceptance condition before launching;
7. stop after the declared gate if the result is only `UNKNOWN`, a necessary
   survivor, a finite-prime fit, or another equivalent relaxation.

This preflight is the replacement for extending a finished line with one
more solver, projection, shell, channel, orbit, or timeout.
