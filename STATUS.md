# Status (2026-08-29)

**15.710 PROVED 208 of the 227 remaining `p=17,s=16` profiles are
impossible, not the endpoint**: every row has nine rigid phase-one `b=16`
directions. A genuine floor-rigid phase-zero `b=0` anchor excludes 176 rows
by forcing
`I=60` and gauge sum 14 while parallel nonnegativity requires at least 15.
Rigid `b=16` anchors in both phases exclude 32 more by forcing `I=68` and
gauge sum 16 while nonnegativity requires at least 17. Exactly nineteen
profiles remain, with residue split `(0,0):5,(7,0):9,(8,0):5`. The
endpoint, residual (ii), R1, Type I, and the limit remain OPEN.

**15.709 PROVED every remaining `p=17,s=16` profile with `u_1=8` is
impossible, not the endpoint**: all 280 such rows retain at least eight rigid
phase-one `b=2` directions. The 66 `(u_0,u_1)=(0,8)` rows retain rigid
phase-zero `b=0` and fail the global gauge comparison. The 214 `(8,8)` rows
retain rigid phase-zero `b=16` and fail the unique-even-fibre identity. The
exact p17 ledger drops from 507 to 227 profiles; every survivor has `u_1=0`
and pair slack at least 96. The endpoint, residual (ii), R1, Type I, and the
limit remain OPEN.

**15.708 PROVED all 54 pair-slack-twenty-four `p=17,s=16` profiles are
impossible, not the endpoint**: all 45 `(u_0,u_1)=(0,8)` rows retain rigid
phase-zero `b=0` and phase-one `b=2` directions; the global-sign identity
forces `I=68`, `g_0+g_1=16`, but parallel nonnegativity requires
`g_0>=9,g_1>=8`. The nine `(8,8)` rows instead force `I=4`. Every row retains
at least two rigid phase-zero `b=16` directions, and summing the exact cells
incident with the unique even fibre gives `N_j=delta_j-15z_j-I<=-3` for a
nonnegative edge count. Exactly 507 p17 profiles remain, all of pair slack at
least 28. The endpoint, residual (ii), R1, Type I, and the limit remain OPEN.

**15.707 PROVED all 78 pair-slack-twenty `p=17,s=16` profiles are
impossible, not the endpoint**: all 69 `(u_0,u_1)=(0,8)` rows retain at
least three rigid phase-zero directions with `b=0` or `2`, while phase one
retains at least eight rigid `b=2` directions. The 15.706 global-sign
comparison again forces `I=68`, contradicting boundary size 16. The nine
`(8,8)` rows have two or three undetermined directions; repair plus the
already-audited conic, complete-13, and complete-14 secant-index bounds
exclude repair depths one through five. Exactly 561 p17 profiles remain, all
of pair slack at least 24. The endpoint, residual (ii), R1, Type I, and the
limit remain OPEN.

**15.706 PROVED both remaining pair-slack-zero `p=17,s=16` profiles are
impossible, not the endpoint**: every mean allocation retains a rigid `b=2`
direction of each quadratic type. Summing its 136 coefficient cells and
comparing with the global finite-edge Paley-sign sum gives
`17I=4+72(g_++g_-)`. Thus `I=68`; the one remaining finite edge makes the
affine odd boundary have size 66, 68, or 70, never 16. This proof is
solver-free. Exactly 639 p17 profiles remain, all of pair slack at least 20.
The endpoint, residual (ii), R1, Type I, and the limit remain OPEN.

**FINITE-FIELD SIGN AUDIT:** three exploratory edge-lift scripts used integer
subtraction on encoded `F_{p^2}` elements. They now use componentwise field
subtraction and a full-edge canonical-conference regression test. The old
raw CP-SAT shards supporting 15.696 are invalid. The corrected archive closes
all twenty logical shards with 22 exact `INFEASIBLE` files; the hard
`022/I=28` case is partitioned over all three possible elevated phase-zero
roles. No old shard or timeout is treated as evidence.

**Cold theorem audit completed; L remains OPEN.** The actual asymptotic gate
is only `Phi(C_p)-m_{p^2+1}=o(p^3)` on a ratio-dense Paley tail
(Proposition 15.20e). The exact gap-2 four-unit gate is stronger. A general
random-plus-edge-descent construction kills local/product Paley stability:
correct-scale edge-local minima can remain `Theta(n^2)` from the signed
Paley orbit, and the product-frame second moments are signing-independent.
The surviving replacements are closest-global all-subsets witnesses, the
full-Max+ dilation-energy form of strong R1, and a signed even-Eulerian
free-energy bound above the fractional-moment barrier. The original `c=2`
free-energy target is false by a linear-margin fractional-moment
construction; `c=3` remains viable. Character orthogonality normalizes the
dilation energy exactly but cannot upper-bound it, and PSD/autocorrelation
surrogates admit explicit violations. See
`evidence/STRATEGY_2026-08-29_COLD_REVIEW.md`.

**15.694 PROVED the exact equality normal form for all four pair-slack-20
`p=19,s=16` profiles, not their exclusion**: every witness must split as an
11-arc `A` and a 5-arc `D`. Each deleted point lies on exactly one `A`-
secant, and line equality permits only `(a,d)=(0,0..2),(1,0..1),(2,0..2)`.
Thus the bad lines are either five triples, one quadruple plus three
triples, or two quadruples plus one triple. Adding any two undetermined
infinity points gives a 13-arc with `c1>=7` when `t=4` and `c1>=8` when
`t=5`; the exhaustive classified maximum is 9. This is a strict search-
space reduction but leaves all seven profiles and every top-level gate OPEN.

**15.693 PROVED all seven pair-slack-16 `p=19,s=16` profiles are
impossible, not residual (ii)**: repair deletes at most four points. With
at most three deletions, two undetermined infinity points produce an arc of
size at least 15 and hence a conic contradiction. With four deletions, the
repaired 12-arc plus two infinity points is a complete 14-arc. Slack equality
forces the four deleted points, and every unused undetermined infinity
point, to have secant index one. This gives at least five such outside
points, while Al-Zangana's exhaustive `PG(2,19)` classification gives at
most four for every 14-arc. Exactly seven profiles remain, with pair-slack
histogram `{20:4,24:1,28:1,32:1}`. The slack-20 rows now require exactly
five repair deletions. The endpoint and every top-level gate remain OPEN.

**15.692 PROVED the binary affine-Radon normal form for the fourteen
remaining `p=19` profiles, not the endpoint**: over `F_2`, the affine
line-point incidence matrix satisfies `A^T A=I+J`. On even point words it
is an isomorphism onto the direct sum of the even directional blocks, with
inverse `x=A^T r`. Hence there are no hidden linear compatibility equations:
the exact remaining condition is `wt(A^T r)=16`. Every profile passes the
mod-four inverse-weight test, and exact `{4,6,8}`-supported distributions
show that its first two stripe-count moments cannot force positive odd
density. The endpoint and every top-level gate remain OPEN.

**15.688--15.689 PROVED the sharp integral quadratic-lift floor, completed
the residue-zero census, and reduced the live `p=19,s=16` endpoint to 14
high-slack profiles, not residual (ii)**:
paired-cube quarter-integrality separates `H=1` from `H>=2`; combining the
latter branch with the exact 15.642 stabilizer weights gives
`4p E[B] >= p-3`, sharply attained by `(1-x_i)(1-x_j)`. At `p=19`, this
excludes phase-zero residues `2,3,4,6`. The minimum residue-zero profile
pair has impossible slack 34 but is not the complete row. Exact completion
leaves 143 phase-labelled profiles. The `PG(2,19)` complete-arc spectrum,
repair, and retained-conic-secant bounds exclude all 129 profiles of slack
at most twelve. Exactly 14 remain with histogram
`{16:7,20:4,24:1,28:1,32:1}`. The endpoint and all top-level gates remain
OPEN.

**15.687 PROVED all 68 pair-slack-20 `p=23,s=20` profiles are
impossible, not residual (ii)**: their undetermined-direction counts are
`2^2,3^36,4^30`. In the 66 rows with at least three directions,
overlapping infinity-point pairs either extend to the same conic, forcing
three collinear conic points, or give a complete 17-arc with five deleted
points of secant multiplicity one. The classified maximum is one. The two
two-direction rows close by the same complete-arc obstruction or the
five-point conic-core floor 24. Exactly 133 arithmetic profiles remain at
`p=23`, all of slack at least 24. The endpoint remains OPEN.

**15.686 PROVED the unique pair-slack-16 `p=23,s=20` profile is
impossible, not residual (ii)**: it has one undetermined direction. Fewer
than four repair deletions give an 18-arc after adjoining that infinity
point and hence the conic-core contradiction. In the four-deletion branch,
the repaired 16-arc plus the infinity point must be a complete 17-arc.
No deleted point can use the infinity point on a secant, and slack equality
forces all four to have secant multiplicity one. Proposition 15.685's
exhaustive five-class certificate shows that the maximum is one. Thus 201
exact arithmetic profiles remain at `p=23`, all of slack at least 20; the
endpoint remains OPEN.

**15.685 PROVED the unique pair-slack-12 `p=23,s=20` profile is
impossible, not residual (ii)**: the repair lemma produces an arc after at
most three deletions. A repaired arc of size at least 18, or an incomplete
17-arc, gives the conic-core contradiction from 15.684. Thus the hard branch
is a complete 17-arc `A` plus three points. Secant-line slack gives
`slack(S)>=4 sum mu_A(x)`; completeness and total slack 12 force all three
points to have secant multiplicity one. Five explicit complete-arc
representatives have pairwise-distinct projective-invariant multiplicity
histograms and exhaust Coolsaet--Sticker's five classified classes. Their
numbers of multiplicity-one points are `0,0,1,0,0`. The profile is
impossible, reducing `p=23` from 203 to exactly 202 arithmetic profiles;
the endpoint remains OPEN.

**15.684 PROVED every positive residue at the `p=23,s=20` next all-finite
endpoint is impossible and reduced residue zero from 1,247 exact profiles
to 203, not residual (ii)**: paired-cube value floors exclude scaled mass
12 and all height-at-most-three cases at mass 16. At height four, equality
in the stabilizer identity forces vanishing on the middle shell; its exact
23-dimensional restriction kernel factors as `(t-6)V_1`, and a
two-replacement affine identity contradicts integrality. For residue zero,
Segre's tangent envelope excludes all 363 arc profiles. The exhaustive
complete-arc classification of `PG(2,23)`, a slack-to-arc repair bound, and
an off-conic secant count exclude another 681 profiles. Exactly 203
arithmetic profiles remain, so `p=23` is reduced but OPEN. The other
endpoints `p=17,19`, later sizes, and the infinity-present remainder remain
OPEN.

**15.683 PROVED the `p=41,s=34` next all-finite endpoint is impossible,
not residual (ii)**: Proposition 15.681 removes all positive residues and
the exact residue-zero census leaves seven 34-arcs and two one-triple
near-arcs. In every arc row, eight exceptional directions have at least 28
tangents. Ball--Lavrauw's polynomial form of Segre's tangent envelope makes
their dual pencil lines double components of the degree-18 envelope,
leaving a conic; the three exceptional secant edges then force at least
three point-pencil lines into that conic. For a near-arc, deleting a triple
point while preserving the exceptional pair gives a 33-arc with high
tangent counts `33^7,31`; the eight double components leave a quartic, two
point-pencil factors leave a conic, and three further point pencils are
again forced. Both branches contradict degree. Only `p=17,19,23` remain at
this second boundary. Later sizes and the infinity-present remainder remain
OPEN.

**15.682 PROVED the `p=31,s=26` next all-finite endpoint is impossible,
not residual (ii)**: Proposition 15.681's integral paired-cube floor removes
all positive residues. The exact residue-zero ledger has fourteen
phase-labelled profiles: eleven are 26-arcs with at least three
undetermined directions, and three have one 3-secant and at least five.
Deleting one triple point in the latter case gives a 25-arc. Coolsaet's
complete classification has no complete arc in `PG(2,31)` of size 23
through 31, so adjoining two undetermined infinity points produces a
27-/28-arc contained in a conic. Three such infinity points force a conic
through three collinear points, a contradiction. Proposition 15.683
subsequently closes `p=41`; the three endpoints `p=17,19,23` remain open at
this boundary. Later sizes and the infinity-present remainder remain OPEN.

**15.681 PROVED the `p=29,s=24` next all-finite endpoint is impossible,
not residual (ii)**: paired-cube averaging applies directly to every
nonzero nonnegative integer-valued quadratic and gives scaled mass at least
`(p+1)/2` or `(p-1)/2` according to `p mod 4`. The floors
`14,16,18,20` remove every positive-residue row at `p=29,31,37,41`. At
`p=29`, pair-slack divisibility leaves 24-arcs with at least four
undetermined directions or a one-triple near-arc with six. Exact
`PGL(2,29)` complement-orbit counts match Coolsaet--Sticker's exhaustive
10 and 5 classes of 25- and 26-arcs, so all are conic-contained; three
undetermined infinity points then give a contradiction. Propositions 15.682
and 15.683 subsequently close `p=31,41`. The three endpoints `p=17,19,23`
remain open at this boundary. Later sizes and the infinity-present remainder
remain OPEN.

**15.680 PROVED the `p=37,s=30` next all-finite endpoint is impossible,
not residual (ii)**: the exact pair ledger leaves only phase-zero residues
`u=2,3,4,5`, each forcing a quotient-zero `b=0` quadratic lift of scaled
mass at most ten. Proposition 15.642 excludes the first three. At mass ten,
stabilizer averaging makes the lift `{0,1,2}`-valued; the exact degree-four
slice-distance floor excludes value two. A self-contained paired-cube
restriction gives density at least `17/74` for a nonzero Boolean quadratic
on `J(37,19)`, contradicting the required `5/74`. Proposition 15.681
subsequently closes `p=29`, 15.682 closes `p=31`, and 15.683 closes `p=41`;
the same boundary at `p=17,19,23`, later sizes, and the infinity-present
remainder remain OPEN.

**15.679 PROVED the next all-finite boundary is impossible for every prime
`p>=43`, not residual (ii)**: at the second even size above `3(p-1)/4`,
phase one is rigid and exact phase-zero quotient arithmetic leaves only
common residues `2<=u<=7`. Every such row forces a quotient-zero `b=0`
direction of scaled mean at most 14. The degree-two slice-distance floor is
strictly larger than 14 from `p=59`; exact pair/lift ledgers close the only
smaller in-scope primes `43,47,53`. Propositions 15.680--15.683 separately
close `p=37,29,31,41`; the three endpoints `p=17,19,23`, later sizes, and
the infinity-present remainder remain OPEN.

**15.678 PROVED the exceptional `p=17,s=14` first all-finite survivor is
impossible, completing this boundary size for every prime `p>=17`, not
residual (ii)**: exact residue arithmetic leaves two pair-equality 14-arc
profiles with common secant distribution `{7:6,6:8,1:1,0:3}`. Adjoining
two of their three undetermined infinity points gives a 16-arc. Sticker's
exhaustive `PG(2,17)` classification has a unique 16-arc class, represented
by conic-minus-two; the third infinity point is off that conic and must lie
on at least four surviving secants, a contradiction. The finite
classification is an explicit external dependency.

**15.675/15.677 PROVED the same first survivor for every prime `p>=19`**:
15.675 closes `p=3,5 mod 8` by positive pair gaps and 15.677 closes the outer
classes from `p=23` by the zero-quotient lift contradiction. Proposition
15.678 supplies the formerly exceptional residue-zero endpoint. Proposition
15.679 closes the next size from `p=43`, 15.680 closes its `p=37` endpoint,
15.681 closes `p=29`, 15.682 closes `p=31`, and 15.683 closes `p=41`; three
smaller endpoints, subsequent all-finite sizes, infinity-present remainder,
residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.674 PROVED the entire infinity-plus-`(p-2)` boundary shell is
impossible for both signs and every prime `p>=17`, not residual (ii)**:
all intermediate odd-fibre floors lie strictly above `p+1`. Exact type-sum
residues leave only an all-`p+1` form or `m-1` mean-`p-1` baselines plus one
mean-`2p` exception, so each type contains at most one intermediate
direction. Pair deficit excludes two `b=1` baseline types, while two
complementary baseline types determine at most two directions and are
collinear. The forced mixed pair has exactly 15.673's four arithmetic rows;
three close uniformly and the sole `p=17` endpoint again has `75>57`.
The infinity-plus-`p` shell, large all-finite boundaries, residual (ii),
Type I, R1, global QVAR, and L remain OPEN.

**15.673 PROVED every endpoint-only infinity-plus-`(p-2)` boundary is
impossible for both signs and every prime `p>=17`, not residual (ii)**:
Propositions 15.671--15.672 close the collinear case. For arbitrary endpoint
profiles `b_d in {1,p-2}`, same-type mean quantization and the four-unit
minimum lift cost leave four two-count arithmetic rows. Pair-deficit equality
would be a `(p-2)`-arc with three undetermined directions; Segre's odd-order
`p`-arc theorem forces those three collinear infinity points onto one conic.
Divisibility and `I<=p-2+2E` close three rows. The sole remaining endpoint is
`p=17,(x,y)=(0,7)`, where the exact inter-fibre `l1` minimum is 75 but only 57
transverse edges exist. Intermediate odd-fibre counts, all-finite large
boundaries, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.670 PROVED/certified every finite `p=11` size-eight boundary is
impossible, not residual (ii)**: affine similarity reduces every finite
eight-set to one of `C(119,6)=3,470,108,187` sets containing field points
zero and one. The verifier audits every scalar/translation direction action
and both phase signs. Complete V100/CUDA and RX 9070 XT/HIP replays have
identical cost-pair histograms, zero survivors, and exact minimum larger type
cost 76 against budget 72; an independent CPU prefix agrees entry by entry.
At `p=11`, infinity plus nine and finite size at least ten remain. Residual
(ii), Type I, R1, global QVAR, and L remain OPEN.

**15.669 PROVED uniform non-Walsh boundary ranges, not residual (ii)**: for
every odd `p>=17`, all-finite `6<=s<=3(p-1)/4` and infinity-present odd
finite `5<=s<=p-4` are impossible. Exact small-prime extensions close
`p=11` infinity plus seven and `p=13` finite eight / infinity plus seven or
nine. The larger floor-plus-pair profiles are only relaxed count profiles.
Proposition 15.670 subsequently closes finite `p=11` size eight.

**15.666 PROVED/certified every finite `p=7` size-eight boundary for both
product signs, not the separate infinity-plus-seven profile or residual
(ii)**: the 1,419,432 finite floor survivors per sign left by 15.664 have
23,892,792 exact mean-allocation leaves. Conditioned mod-seven and mod-three
omission scans intersect in 181,104 leaves. Exact 22-row local, all-triple,
and four-positive joins reduce these to 124,745, 78,126, and 62,892. A
single-catalog filter rejects 3,777, and the complete lossless base-seven
meet-in-the-middle join rejects the other 59,115, leaving zero. CPU/CUDA
prefixes agree at every stage, and three older full multi-characteristic
joins independently reject representative residual classes. The nonsquare
anti-isometry transfers the result to `c_H=+1`. Thus the complete finite
`C(49,8)` census is closed for both signs. Infinity plus seven, residual
(ii), Type I, R1, global QVAR, and L remain OPEN.

**15.664 PROVED/certified the four-allocation `p=7` size-eight stratum for
both product signs, not the full size-eight case or residual (ii)**: exact
allocation reconstruction partitions the 24,983,238 survivors left by
15.663. The dominant class contains 2,245 ordered profiles, 23,563,806
boundaries, and exactly four mean allocations per boundary, for 94,255,224
leaves. One quadratic type has floor sum 24, and each leaf raises one of its
four directions by eight. For every raised direction, 112 exact mod-seven
dependencies vanish on its complete 35-column score block. A 22-row
conditioned V100 sweep checks all 450,978,066 ranks and leaves 1,191
allocation leaves; all 135 mod-seven dependencies leave 1,176. NUKA
independently rebuilds the matrix, kernels, catalogs, and candidates. The
1,176 survivors are exactly the `4*7*42` affine-line-plus-off-line-point
family. Every member has two mod-seven and 756 mod-three catalog rows, with
empty intersection, so zero exact catalog choices remain. The nonsquare
anti-isometry transfers the exclusion to `c_H=+1`. Exactly 1,419,432 finite
size-eight floor survivors per sign remained at this stage; Proposition
15.666 subsequently closes all of them. The infinity-plus-seven profile,
residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.663 PROVED/certified the forced-floor `p=7` size-eight stratum for
both product signs, not the full size-eight case or residual (ii)**: exact
profile reconstruction selects 2,016 ordered profiles and 83,770,008
nonconic boundaries per sign whose two quadratic type floors are `(32,32)`.
The exact type mean sums force every directional mean to its floor, leaving
at most one 36-row catalog. An exhaustive V100 scan checks all 450,978,066
boundaries, leaves 526 candidates under an eight-row mod-seven projection,
and rejects all 526 with the complete 135 dependencies. NUKA independently
rebuilds the `282 x 1225` score matrix, rank-147 row space, dependencies,
catalogs, and all candidate failures. A nonsquare anti-isometry bijects the
stratum between signs. With 15.662, the open size-eight floor scope falls
from 108,754,569 to 24,983,238 boundaries per sign; Proposition 15.664
subsequently reduces that remainder to 1,419,432 and 15.666 closes it. The
infinity-plus-seven profile, residual (ii), Type I, R1, global QVAR, and L
remain OPEN.

**15.662 PROVED/certified the minimum-eight-odd-secant/conic `p=7`
size-eight subbranch for both product signs, not the full size-eight case or
residual (ii)**: complete CUDA censuses check all 450,978,066 finite
eight-point boundaries per sign. Exactly 6,174 attain eight odd secants and
are affine conics by Segre's theorem; 4,851 fail the exact floor. The 1,323
floor survivors form 32 stabilizer orbits. On `c_H=-1`, the 25 saturated
orbits have 600 exact mean allocations, split into 355 initial CP-SAT
exclusions, six long CP-SAT exclusions, and 239 multi-prime catalog joins.
The seven exceptional orbits have 1,260 allocations, split into 172 initial
exclusions, 662 ordinary projected V100 joins, and 426 high-direction-
eliminating joins. Independent audits leave zero allocations. A nonsquare
Paley anti-isometry maps all 1,323 survivors bijectively to `c_H=+1` and
transfers the exclusion. Its 108,753,246-case nonconic remainder is reduced
by 15.663--15.664 to 1,419,432 finite floor survivors per sign and closed by
15.666. The infinity-plus-seven profile, residual (ii), Type I, R1, global
QVAR, and L remain OPEN.

**15.661 PROVED/certified the `p=7` six-finite branch and therefore every
size-six boundary for odd `p>=5`, not residual (ii)**: exact floors leave
3,856,300 of `C(49,6)=13,983,816` boundaries and 80,704 square-semilinear
orbits. Simultaneous mod-three/mod-seven catalog joins reject all 160,745
elevation cases in 80,519 ordinary orbits. Compact exact high-mean models,
930 complete mean allocations, and 120 final low-catalog joins reject the
remaining 185 orbits. Independent NUKA and V100 floor sweeps have the same
survivor hash; NUKA reproduces the ordered orbit catalog and ordinary
exhaustion. A nonsquare anti-isometry transfers the opposite product sign.
With 15.657--15.660, every size-six boundary is closed. Proposition 15.662
subsequently closes the `p=7` minimum-eight-odd-secant/conic subbranch for
both signs, 15.663--15.664 reduce the finite size-eight floor remainder to
1,419,432 per sign, and 15.666 closes it. Residual (ii), Type I, R1, global
QVAR, and L remain OPEN.

**15.660 PROVED/certified every `p=5` size-six boundary, not residual
(ii)**: four exact catalogs are rebuilt from definitions for both product
signs and both infinity bits. Signed symmetry and complete coarse SCIP
batches reduce all survivors to six classes. Independent layered audits
close classes 0, 881, 2529, 3032, 4731, and 4939 while reconstructing every
finite quotient; all six have zero unresolved or feasible leaves. Together
with 15.657--15.659, size six then remained only for six finite points at
`p=7`; Proposition 15.661 subsequently closes that branch. Boundaries of
size at least eight, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.659 PROVED/certified the negative-product `p=7`
infinity-plus-five boundary exclusion, not all size-six boundaries**:
phase one and the exact type budget force exactly one scaled-mean-14
direction in each quadratic type. Independent V100 and Soulkiller floor
sweeps agree on 83,496 survivors among all `C(49,5)=1,906,884` boundaries,
and serial and GPU-seeded quotients agree on 1,750 square-semilinear orbits.
Left-kernel affine spans reject 2,205 of 2,230 elevation cases; exact testing
of all 32,400 catalog pairs in the remaining 25 cases leaves zero survivors.
NUKA and Soulkiller independently reproduce both modular stages. Together
with 15.658 this closes both `p=7` infinity-plus-five signs. Proposition
15.660 subsequently closes every `p=5` size-six case, and 15.661 closes six
finite points at `p=7`. Boundaries of size at least eight, residual (ii),
Type I, R1, global QVAR, and L remain OPEN.

**15.658 PROVED/certified the positive-product `p=7`
infinity-plus-five boundary exclusion, not all size-six boundaries**:
all eight directions have phase zero and scaled mean eight. The complete
`J(7,4)` classification is unique for odd-fibre sizes 1, 3, and 5, so each
of the `C(49,5)=1,906,884` finite boundaries fixes all 280 affine score
right sides. The 135 left-null dependencies of the common `282 x 1225`
system over `F_7` reject every boundary. An exact V100 integer sweep and
an independent NUKA/NumPy sweep agree on the full direction-mask histogram
and zero survivors. Proposition 15.659 subsequently closes the
negative-product infinity branch, 15.660 closes every `p=5` size-six case,
and 15.661 closes six finite points at `p=7`. Boundaries of size at least
eight, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.657 PROVED every six-point boundary for odd `p>=11`, not residual
(ii)**: exact positive quadrature extends the 15.652 floor table through
six odd fibres. If `s` finite boundary points have `b_d` odd fibres in
direction `d`, unique pair directions give
`sum_d(s-b_d) <= s(s-1)`. The resulting floor cost exceeds the total
affine slack budget by `p^2-9p+10` for infinity plus five points from
`p=11`, and by `p^2-12p+7` for six finite points from `p=13`. At `p=11`,
the separate quadratic-type budgets require deficits 20 and 18, exceeding
the geometric budget 30. Thus size six is closed for `p>=11`; Propositions
15.658--15.661 subsequently close `p=5,7`. Boundaries of size at least
eight, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.656 PROVED/certified every `p=5` size-four boundary, not residual
(ii)**: quotienting either complete eigenshell by antipodes gives 130 score
rows with normalized edge-column sum 26 and total slack 78. Edge count, the
distinguished edge, and 130 bad-count rows form a `132 x 325` system of rank
67 over `F_5`. Exact bounded lift syndromes exclude 712 orbit cases; the
only mod-five timeout is independently `INFEASIBLE` modulo seven. A signed
nonsquare Paley anti-isometry transfers the no-infinity negative-sign
exclusion to the positive sign. A fresh structural audit reconstructs all
four boundary classifications, shell ranks, parity/lift masses, and the
489-orbit sign-transfer bijection. Thus all 1,202 floor-surviving orbit/sign
cases, covering 26,450 boundary/sign cases, are excluded. With 15.632 this
closes `p=5`; with 15.652--15.655 every size-four boundary is closed for
every odd `p>=5`. Proposition 15.657 subsequently closes size six for
`p>=11`. Proposition 15.660 subsequently closes `p=5`, and 15.661 closes
the six-finite `p=7` branch. Larger boundaries remain OPEN.

**15.655 PROVED/certified the complete unsaturated `p=7` four-finite
exclusion, not residual (ii)**: the 280 exact affine score equations plus
edge count and the distinguished edge form a common `282 x 1225` integer
system. Its rank over `F_7` is 147, giving 135 exact left-null dependencies.
Complete syndrome joins reject all 1,716,742,440 catalog tuples in all 2,408
cases covering 23,520 boundaries and 518 orbits. An independent
implementation rebuilds the matrix, dependencies, catalog right sides, and
coverage and again finds zero survivors. The 15.654 nonsquare anti-isometry
transfers the result between product signs. Thus 15.653--15.655 close every
`p=7` size-four boundary; 15.656 subsequently closes `p=5`. Larger boundaries,
residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.654 PROVED/certified the doubly saturated `p=7` four-finite
exclusion; 15.655 subsequently closes the unsaturated remainder**: complete
Johnson-space
classification leaves one phase-zero and 36 phase-one saturated `b=4`
slacks. Exactly 58,800 boundaries per product sign have both type costs 32;
they form 1,225 square-semilinear orbits, whose exact coefficient models are
all infeasible with zero unknown and zero feasible case. A signed nonsquare
Paley anti-isometry fixes the distinguished edge, exchanges eigenshells,
preserves normalized scores, and flips the product sign, so one orbit sweep
covers both signs. Its 23,520-boundary unsaturated complement is now closed
by 15.655. The `p=5` branch is subsequently closed by 15.656. Larger boundaries, residual (ii),
Type I, R1, global QVAR, and L remain OPEN.

**15.653 PROVED/certified the `p=7` infinity-plus-three-finite exclusion,
not all size-four boundaries**: 15.652 already closes `c_H=-1`. For
`c_H=+1`, every direction saturates its type budget. An exact rank-21
Johnson-space calculation leaves one three-odd-fibre slack among 630 sparse
corrections. The resulting coefficient model reduces all 18,424 finite
triples to 416 square-semilinear boundary orbits, all infeasible with zero
unknown and zero feasible case. Four finite points at `p=7` are now closed
by 15.654--15.655. The `p=5` branch is subsequently closed by 15.656. Larger boundaries, residual
(ii), Type I, R1, global QVAR, and L remain OPEN.

**15.652 PROVED the complete four-point boundary exclusion for every odd
prime `p>=11`, not residual (ii)**: exact positive quadrature gives the
parity floors for zero through four odd affine fibres. Four finite boundary
points create only six pair-collision directions; infinity plus three finite
points creates only three. In either case the separate quadratic-type
budgets cannot be met for `p>=11`. With 15.632 and 15.650--15.651, the first
still-open boundary size there is at least six. Propositions 15.653--15.655
subsequently close the exceptional `p=7` size-four cases, and 15.656 closes
`p=5`. Larger boundaries, residual (ii), Type I, R1, global QVAR,
and L remain OPEN.

**15.651 PROVED/certified the complete positive-product
infinity-plus-point exclusion, not residual (ii)**: exact additive
coefficient equations close all seven `p=5` arithmetic cases. Fibrewise
`l1` profiles and a uniform-edge-type capacity argument close `p=11,13`.
At `p=7`, the rigid type split has 2250 surviving five-stars and 56
square-semilinear orbits per populated type; all 112 fixed-star models are
infeasible. The alternative all-eight-`kd=1` profile has three exhaustive
normalizations, also infeasible. With 15.643, `c_H=+1` is closed for every
odd `p>=5`; with 15.650, both signs of `D={infinity,v}` are closed. Other
boundary shapes, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.649--15.650 PROVED/certified the complete negative-product
infinity-plus-point exclusion, not residual (ii)**: 15.649 classifies all
1764 balanced `p=7` exceptional lifts and excludes all 6076 star-orbit
representatives. 15.650 reduces `p=5` to 24 arithmetic profiles and 33
placement orbits, all infeasible. Together with 15.647--15.648, `c_H=-1`
is closed for every odd `p>=5`. Other boundary shapes and all top-level
gates remain open.

**15.648 PROVED/certified the `p=11,13` negative two-point cases and four
`p=7` profiles, not residual (ii)**: at `p=13`, a zero-baseline direction
has prescribed inter-fibre `l1` norm at least 48 but only 44 transverse
edges, closing both orientations analytically. At `p=11`, exact CP-SAT
infeasibility covers both count orientations and all three proven
type-preserving exceptional-pair orbits. Direct all-pair certificates at
`p=7` exclude `(x,y)=(0,3),(0,6),(3,0),(6,0)`. The negative two-point
remainder is now only `p=5` and balanced `p=7 (3,3)`. Positive-product
small primes, other boundaries, residual (ii), Type I, R1, global QVAR, and
L remain OPEN.

**15.647 PROVED negative-product branch exclusion for every odd `p>=17`,
not residual (ii)**: the exact signed directional mean gives
`a_d-a_e=(p+1)(P_d-P_e)` within each quadratic type. Proposition 15.642
guarantees a baseline direction per type for all `p>=7`; since each type's
total lift excess is `p+1`, this forces exactly one exception of excess
`p+1` per type without the asymptotic slice theorem. Baseline coefficient
divisibility then makes both baseline parallel counts multiples of
`(p-1)/2`. For `p>=17` both must vanish, leaving a `4p-1` infinity star that
two finite edges cannot repair. Together with 15.643, both product signs of
`D={infinity,v}` are now closed for `p>=17`. The primes `5,7,11,13`, other
boundaries, residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.646 PROVED negative-product branch exclusion for sufficiently large
primes; superseded by 15.647**: in 15.644's normal form, summing the exact
inter-fibre identity gives signed transverse-edge sum zero in every baseline
direction. The exceptional splits `(U,V)=(3,1)` and `(1,3)` give global
finite-edge signed sums `+2` and `-2`; after removing two parallel edges,
a negative or positive baseline direction instead has transverse sum `+4`
or `-4`. Both are contradictions. The threshold remains qualitative because
15.644 uses an asymptotic slice-distance theorem. Small primes, the positive
branch at `p=5,7,11,13`, other boundary profiles, residual (ii), Type I, R1,
global QVAR, and L remain OPEN.

**15.645 PROVED fibre rigidity, no longer needed for branch exclusion**: inside 15.644's
large-prime negative normal form, every baseline direction meets the
`2p-1` infinity-neighbor set in the ideal profile `(1,2,...,2)` relative to
the boundary fibre, or in a profile obtained after one unit transfer. This
follows from an exact integral `l1` minimization; two transferred units
already exceed the `2p` transverse-edge budget for `p>=7`. Simultaneous
classification as two affine lines remains open but is bypassed by 15.646;
all top-level gates remain open.

**15.644 PROVED asymptotic normal form, not branch exclusion**: the
near-optimal degree-two polynomial-distance theorem forces exactly one
nonbaseline direction of each quadratic type in the `D={infinity,v}`,
`c_H=-1` branch for all sufficiently large odd primes. Exact mean,
divisibility, `l1`, boundary, and sign-product constraints then force
`2p-1` infinity edges, `2p+2` finite edges, exactly two parallel finite
edges in every baseline direction, and exceptional parallel counts `1,3`
(the negative-type count is odd). Proposition 15.646 now excludes this normal
form. Residual (ii), Type I, R1, global QVAR, and L remain OPEN.

**15.643 PROVED branch exclusion for `p>=17`, not residual (ii)**: in the
`D={infinity,v}`, `c_H=+1` branch, 15.642's pointwise baseline makes every
signed inter-fibre edge matrix additive. Divisibility forces each finite
parallel count to be a multiple of `(p-1)/2`; its exact `l1` edge budget
then excludes every populated direction for odd `p>=17`. The sole arithmetic
endpoint is an all-infinity star, whose boundary is not two. The cases
`p=5,7,11,13`, negative-product primes below the qualitative 15.644
threshold, other boundary profiles, residual (ii), Type I, R1, global QVAR,
and L remain OPEN.

**15.642 PROVED rigidity/sparsity, not residual (ii)**: a nonzero
nonnegative integer-valued quadratic on `J(p,(p+1)/2)` has an exact
stabilizer mass floor; combining it with the all-parameter degree-two slice
distance lemma gives support at least `(p^2-1)/(16p(p-2))`. For residual
boundary `D={infinity,v}`, `c_H=+1` forces every affine slack pointwise to
its parity baseline `x_s`. For `c_H=-1`, each quadratic direction type has
at most three nonbaseline directions, uniformly in `p>=5` (two at `p=7`).
This sharply reduces but does not exclude that boundary. Residual (ii),
Type I, R1, global QVAR, and L remain OPEN.

**15.641 PROVED route kill, not R1**: at `p=11`, the justified pre-second
infinity rows and geometric cusp gaps have rank 29 in the 66-dimensional
Kohnen space; the complete second-shell row raises this to 30, leaving a
36-dimensional kernel. An exact 21-coordinate rational witness annihilates
all of those rows and the second shell while its half-cusp R1 target is one.
Thus the current linear modular shell/cusp data do not determine or sign the
target. Additional exact data or nonlinear theta positivity could still
prove R1. Global QVAR, R1, residual (ii), Type I, and L remain OPEN.

**15.632 PROVED reduction + branch kill, not residual (ii)**: for every odd
affine separator `H`, the scaled integer slacks satisfy the type-split identity
`sum_{eps_d=tau} a_d=(p+1)(|H|-3p)/2`.  The odd-degree boundary fixes each
slack's Johnson-slice parity; exact quadratic majorization gives
`a_d >= 2 ceil(p M(p,b_d,eta_d))` and the corresponding separate budget in
each quadratic direction type.  At `|H|=4p+1`, an empty/Eulerian boundary
exceeds one type budget by `(p^2-1)/2`, so that branch is empty for every odd
prime.  A direct corrected p=5 affine solution survives with
`a=(12,4,0,6,10,4)` and boundary infinity plus one affine line, so nonempty
boundary profiles, the full non-Walsh shell, residual (ii), Type I, R1, and L
remain OPEN.

**15.631 PROVED structure, not R1**: for every Max+ odd vector `y0 in L`
and every `u=Pz in L*`,
`<u,y0> = 2p||u||^2 (mod 2)`.  Thus Poisson summation on `y0+2L` has
the scalar radial phase `(-1)^(2p||u||^2)`, independent of the Max+ vector
and of all other glue data.  The first dual shell is exactly `+-Pe_i` at
norm `1/2`; every other nonzero dual vector has norm at least `(p-1)/p`.
For the R1 degree-four harmonic, the phased first dual shell contributes
`||W||_F^2/[8(d+2)]`.  Higher dual-shell harmonic sums remain uncontrolled,
so R1 and global QVAR remain OPEN.

**15.630 PROVED structure, not R1**: the adjacent ETF lattice
`L*=P Z^n` has exact minimum `1/2`, with complete minimum shell
`{+-P e_i}` and kissing number `2(p^2+1)` for every odd prime.  The proof
uses the circle-frame energy identity, integer balancing for nonzero common
profile sum, and the Prop. 15.629 glue as a projective Reed--Solomon MDS code
plus Newton identities at common sum zero.  This recovers the finite
`(5,10)` and `(13,26)` ETF results of Böttcher et al. and settles the standard
Paley `(25,50)` minimum-shell question (the first `?` in their table).  It
does not separately assert a basis chosen from the minimal vectors.  It is the ordinary shell of
`L*`, not the odd-coset shell of `y0+2L`; R1 and global QVAR remain OPEN.

**15.629 PROVED structure, not R1**: for
`L=ker_Z(C-pI)` and the square-circle lattice `A`, the profile coefficient
kernels give `L/A = direct_sum_{d=1}^{m-2} ker(sum c_j t_j^d)`.  Hence
`[L:A]=p^C(m-1,2)`, `det(L)=2p^(m^2)`, `L*=P Z^n`,
`L*/L = Z/2 plus (Z/p)^(m^2)`, and `level(L)=4p` for every odd prime.
This makes the R1 shell a precise level-`4p` coset-theta problem; it does not
bound its degree-four harmonic coefficient.  R1 and global QVAR remain OPEN.

**15.628 PROVED close**: eligible nonsquare GQR circles span the fixed-edge
slice, and arbitrary affine halfspaces realise every eligible circle as a
U-difference.  Thus `dir(U)=H_0 cap ker(e_0+e_infinity)` and W1, W2, and
Walsh 15.406 E hold for every odd prime (p=3 direct rank 4).  The explicit
p=19 witness supersedes the generic-solver UNKNOWN.  This does not close
5+-level / even-`k>4p` or the other E1 leftovers; L remains OPEN.

**15.627 PROVED kill/split**: octic+(a mod 8) does not rescue
the linear stay box. Split-involution class is W2-nonempty at
p=31 (76 hits). Counting identity OPEN.

**15.626 PROVED kill/split**: bounded box |u,v,w|≤4 |k|≤8 of
ua+vb+wi+k is not W1 on p=a²+64c² (61 primes; not all of Z⁴).
W2 t=-2 at p=17, not a p-law (p=31).

**15.625 PROVED split**: W1 for (2/p)_4=−1 via d=−(p−1)/8.
Residual W1: p=a²+64c² (601…). CLASS exhaustive. W2 p-law OPEN.

**15.624 PROVED split**: inversion −1/x is Max− not in U. Named
W2 at p=11: π=x/(m(x+2)). Two Auts cover p=5,7,11,13 not 17+.

**15.623 PROVED split**: W1 for p≡73 or 97 (mod 120) via d=−3.

**15.622 PROVED split**: W1 for p≡17 (mod 24) via d=−2. Named
W2 at p=5: switched x/(x−1) clears Φ3.

**15.621 PROVED split**: W1 for every p≡5 (mod 8) via d=−1.

**15.620 PROVED kill/split**: s_N is not a W1 p-law (ε=0 at p=29).
Stay translates still hit ε=1. χ_p-pullback misses Φ3 at p=5.

**15.619 PROVED split**: odd_QNR(s_N)=0 for every p≡1 (biquadratic
fiber count (p−1)/4 ± b/2).

**15.618 PROVED split**: Φ=ε on W_0 (odd D-index on QR xor QNR);
scale Φ((D−I)γ)=2p−1≡1. s_N=f∘φ a p-law. 1_M coprime to g.

**15.617 PROVED split**: w∈(f)R iff f|content. 15.616 Walsh p=11
withdrawn. s_N nsq-class stay-sum has ε=1 at p=5,13,17 (W1 p≡1
construction; ε not a p-law). W2 generic at p=5,7. Walsh ∀p OPEN.

**15.616 PROVED** (not a leftover close): f(D)(z+Dz)≠0 for every
irred of g at p=5,7,11. Walsh at p=11 (W1 from 15.614 + W2).
Walsh ∀p / leftover 2 still OPEN.

**15.615 PROVED kill** (not a leftover close): ε(two-fiber)=1 iff
p≡1 is false (p=17). W2 named-pool misses g at p=11. leftover+splus
at p=5 k=20 already empty (15.528). residual_ii / Walsh OPEN.

**15.614 PROVED split** (not a leftover close): W1 for every p≡3
(mod 4) via W-lift of z+Dz; ε=p(p−1)/2 (mod 2). W1 for p≡1 and
W2 open (named vectors killed by g at p=11). Walsh OPEN.

**15.613 PROVED split** (not a leftover close): named halfspace-anti
z∈U; ε(y+Dy) constant on U; ε=odd Krylov sum vs 1_M. The ε-bit
(1 iff p≡3; T_α at p≡1) is certified not a p-law (Fable BLOCK).
Walsh OPEN.

**15.612 PROVED split** (not a leftover close): Walsh ⇔ W1 ∧ W2
on Aut-invariant ideals of W_0. CLASS of maximal Aut-invariant
ideals is a p-law (Fable PASS). W1 certified p=3,5,7, Fable BLOCK
as p-law (Frob one-point fails at p=5). Walsh OPEN.

**15.611 PROVED** (not a leftover close): W ≅ F2[M] ≅ F2[X]/(X^N+1)
as D-modules (Ind of regular F_p^×). W_0 is the unique D-invariant
hyperplane. dim ker((D−I)^2)∩W_0=2 is a p-law (a=v_2(N)≥2).
Walsh is F2[D]-ideal generation. Fable xhigh PASS. residual_ii False.

**15.610 PROVED** (not a leftover close; uniqueness KILL): W_0=extra^⊥∩W;
ker((D−I)^k)∩W_0 is I-invariant (char 2). ker((D−I)^2) has dim 2
(certified p=3,5,7,11), so Aut({0,∞}) uniqueness for Walsh is DEAD.
Walsh spanning OPEN. Fable xhigh BLOCK on irreducibility.

**15.609 PROVED** (not a leftover close): H0'=H0^⊥=rowspan(S);
opposite-type circles are never tangent, so I(H0)=H0 for every
odd p. Span of all square circles is rowspan(S). Walsh OPEN.

**15.608 PROVED** (not a leftover close): χ-type of an F_p-subline
is independent of basepoint; two PSL-orbits of circles (square vs
nsq); I∈PSL preserves each orbit. 1∈dir(U) by antipodes. I(H0)=H0
certified, not a p-law. Walsh spanning OPEN.

**15.607 PROVED** (not a leftover close): W is irreducible as a
G_aff^□-module for every odd p. F_p^× ⊂ M preserves each W^H and
transits the Φ_p-factors (fail: W^H simple as C_p-module at p=7).
Hence dir(affine_span(Max−))=H0. Walsh 15.406 E (xor-slice) OPEN.

**15.606 PROVED** (not a leftover close): square F_p-line averages
vanish on W; nsq averages are orthogonal projectors of rank p−1
summing to I, so W=⊕ W^H. Squares permute the nsq summands
transitively. If 2 is a primitive root mod p, each W^H is simple
and W is G_aff^□-irreducible. Not a p-law at p=7 (Φ_7 splits);
spin still fills W. Walsh OPEN.

**15.605 PROVED** (not a leftover close): Paley A of order p² satisfies
A²=A over F2 (q≡1 (mod 8); fail Paley-13). H0=⟨1⟩⊕W with W the
translate-span of the 15.604 extra vector, dim (q−1)/2. W is a
G_aff^□-module generated by extra, W^G=0. Irreducibility of W OPEN.
Walsh OPEN.

**15.604 PROVED** (not a leftover close): 1_QR ∈ H0 iff p≡1 (mod 4);
1_QNR ∈ H0 iff p≡3 (mod 4), via χ_q|_{F_p^×}≡1 and 15.598 line
sums. ker(D−I)∩H0 = ⟨1, extra⟩ dim 2; D^N=I. Restriction H0→F2^{QR}
is not onto (census). H0/⟨1⟩ irreducible OPEN. Walsh OPEN.

**15.603 PROVED** (not a leftover close): rank(S')=n/2; H0 ∩ H0'=⟨1⟩
(all affine line sums ⇒ constant on AG(2,p)); H0+H0'=even-weight.
Heart splits as two (q−1)/2 pieces. Irreducibility OPEN. Walsh OPEN.

**15.602 PROVED** (not a leftover close): G_aff^□ permutes rows of S;
unique 1-dim invariant subspace of H0 is ⟨1⟩. Inversion permutes the
0-pencil only. H0/⟨1⟩ irreducible OPEN. Walsh 15.406 E OPEN.
residual_ii False.

**F̂ Paley-field square/norm KILLED** (`src/e1_gmin_qvar_fhat_norm.py`,
not a leftover close): F̂=(λ−6)q² has odd valuation at a residue-degree-2
prime in Q(ζ_{q−1}) and Q(√p*) at p=5 (13) and p=7 (409). Bochner
F̂=|A|² in the Paley character field is false. F̂≥0 / leftover 1 still
OPEN.

History/references: `evidence/HISTORY_AND_REFERENCES.md` (MO/X/Paata/pre-internet). Not a close.

No leftover flag flipped. `phi_F_ge_6` imports
`leftover1_qvar_and_principal_proved()` from
`src/e1_gmin_leftover1_qvar_principal.py` (unnumbered hinge). That
predicate is False until **GLOBAL** QVAR (mixed-k; not per-stratum
k≥7) **and** principal ||δ||² / R1 are proved for all p≥5. Per-stratum
k≥7 is false at (41,7). Fail-eqs: drop 16 in QVAR; 32↦16 in λ_exc; drop
q+5 in V_sph; n-14↦n-6 in B_min. Uniform Paley E[S²]<20+12/p is false
(p=5 leftover witness 25.17>22.4). Leftover 2/3 False. Gsum / pairing
/ Aut-Schur False. L OPEN. Live `e1` is still the old AND.

**15.597 Theorem A* PROVED** (not a leftover close): Φ_part = λ̄ I on Z
for every prime p≥5. Leftover 1 ⇔ Φ_δ ⪰ −(2n+20)/(n−6) I. Corollaries:
tr(Φ_δ)=0; 0 ≤ λ_min(Φ) ≤ λ̄; remaining gap is [0,6). Distinct from
15.108 Theorem A*.

**15.598 PROVED** (not a leftover close): square-direction affine
F_p-lines cut Max−, \(\sum_{\{\infty\}\cup L} y=0\). Walsh ∀p reduces
to spanning of the xor-hyperplane of affine_span(Max−). residual_ii
stays False.

**15.599 PROVED** (not a leftover close): rank(SSᵀ)=n/2−1 and
rank(S)≤n/2, so rank(S)∈{n/2−1,n/2}. Equality n/2 certified p=3..37.
p=11 dim-60 was the y_∞=+1 half; antipodes restore n/2. Aut_e
irreducibility is false. Walsh spanning open. residual_ii stays False.

**15.600 PROVED** (not a leftover close): rank(S)=n/2 for every odd
prime, via \(1\in\ker S\cap(\ker S)^\perp\). dim H0=n/2 is a theorem.
Walsh spanning of the xor-hyperplane still open. residual_ii False.

**15.601 PROVED** (not Walsh, not a leftover close): pencil through 0
gives 1_QR ∈ rowspan(S) (p≡3) or rowspan(S)+span{ℓ} (p≡1). 15.406 E
stays OPEN; single-orbit Aut_e spanning is false at p=5.

`tr(Phi^2)=4||M||_F^2-3n^2+2n^2(n-1)/p^2` is a proved identity, not a
bound (TECHNICAL_NOTES §4). Leftover 3 next target: p=5 finite from C;
p≥7 `|μ|≤2/n` on |κ|=1 would close Type I (strictly stronger than L;
census 109/2863 vs 2/50 at p=7). Not imported. `|μ|≤|f4|` false at p=7.

Remaining Max+-free estimates isolated in
`evidence/NOTE_2026-08-21_remaining_general_p_estimates.md` and HANDOFF
2026-08-21: leftover 1 is QVAR k≥7 plus principal δ-room; leftover 2 is
leftover+splus at k=4p (Walsh interior 4-level only; Paley ES2 majorant
dead); leftover 3 is `|μ|≤L` or, for p≥7, `|μ|≤2/n`. p=13 orbits are
not a close.

HEAD 15.585. Sandwich and Paley \(\rho=1\) are proved. \(L=\tfrac12\) is not. Three leftovers remain (Lemma D is True). Ensemble \(Q_\tau\) and \(D=\lvert H_+\rvert/(2p)\) are unnamed in \(p\); \(\phi_F\) is not imported. 15.564: \(F=-2(3p^2+2)\). 15.573–15.575, 15.581–15.582: exclusive mix and named \(\mu_{1d}\), \(\mu_{k=3}\) give \(p=5\) \(Q_{++}/q^2=48/13<26/7\); the mix is not a general \(Q_\tau\). 15.578: \(2\chi\) fourth moment \(4p/(p-1)\); \(p=7\) occupancy mix is not a \(p\)-law. 15.585: leftover+\(s_+\) at \(k=4p\) forces \(\min_+=2\). 15.559 / 15.565 / 15.577 / 15.580: Aut\(_e\), Max± of \(C\), 1D Johnson, and Galois+\(F\) do not kill multi-level Type I. Pointwise \(Q^{++}\le4q^2\) fails on a positive fraction of Max+ at \(p=5,7\). Live `e1` is still the old AND. 15.558: \(J_{\mathrm{all}}=(1/8)\sum\chi_\Omega^{|\varepsilon|}S_\varepsilon\) and \(J_T\) via \(S_\Omega\); \(G_3=G_{02}+q^{-3}\sum_T Q_{3,T}J_T\); \(Q_{3,02}=-4N(2p^2+1)/p\) from \(F=-2(3p^2+2)\) certified \(p=5..23\). \(Q_{3,T}\) still live except generic. 15.559: Aut_e DEAD as a name of \(A_{\mathrm{full}}\) (inversion mixes full-Ω at p=7). 15.560: leftover+splus empty nF=0 at p=5 \(k=26,28,30\). 15.550: \(S(\lambda)=\mathrm{Kl}(1,\lambda^2/4)\) on every odd \(q\) (fail \(\pm G\), drop \(1/4\), \(\mathrm{Kl}(1,\lambda)\)). 15.553: Term0 of \(K_\lambda^{\mathrm{all}}\) is \(N(1+p^{-1}\sum_{\delta\neq0}\chi(\delta)e_{-\xi}(\delta)\mathrm{Kl}(1,r^2\xi^2\delta^2/16)^2)\); Term0/\(N=5,101,197,485\) at \(p=5,7,11,13\); fail drop \(\chi\); fail \(\mathrm{Kl}\mapsto G\); fail drop \(1/16\) at \(p\ge7\) (\(16=1\) in \(\mathbb F_5\)). Term0 is not the bulk (\(650\) vs \(249050\)). 15.551: \(\hat z\)-support is a Galois union of \(\Omega\)-lines; 1-line iff \(r\in\mathbb F_p\). 15.552: leftover+splus empty nF=9,11–13 at p=5 \(k=22\); nF=10 TLE. 15.549: \(K_\lambda=K_\lambda^{\mathrm{all}}-Nq(3q-5)\) with \(K_\lambda^{\mathrm{all}}\) the unrestricted Paley 4-linear (\(\chi(0)=0\)). Fail drop the Gauss collision; fail subtract only \(2Nq(q-2)\). \(K_\lambda^{\mathrm{all}}\) Ω-bulk not a p-law. 15.548: \(H_+\) 2-point \(G(a,b)=(N/p)\chi(b-a)\) on \(\mathbb F_q\); \(|\kappa|=3\) layer of the y-first 4-linear is \(K_3=(K_1+K_\lambda+K_{1-\lambda}+K_{\lambda(1-\lambda)})/4\). At p=5 \(K_3=83350\), \(K_\lambda=21550\), \(K_{1-\lambda}=K_{\lambda(1-\lambda)}=111550\). Fail \(G\equiv0\); fail drop the \(\lambda(1-\lambda)\) channel (\(55462.5\neq83350\)). The \(\chi(\lambda)\chi(1-\lambda)=+1\) bucket (100150) is not \(I_3\). \(K_\lambda\) unnamed in \(p\). \(16pA\) / \(Q_\tau\) still open. 15.547: leftover+splus at p=5 \(k=22\) (even \(>4p\)) is empty for nF=0,3,4,5,6,7,8,14 (HiGHS Infeasible, \(S\ge2\); nF=8 in 814s/119211 nodes). leftover-only official nF=3 exists (min\(_+=-8\)). leftover-only nF=1,2 and \(\ge15\) empty. nF=10 TLE and other even \(k>4p\) stay open. Fail: leftover-only empty. residual_ii_k_eq_4p_empty stays False. **Three leftovers remain** (Lemma D already True). 15.545: \(\mathrm{NUM\_SUM}=n_{1d}Q_{1d}^{++}+M_{\mathrm{NL}}\) and \(M_{\mathrm{NL}}=16pA-n_{1d}Q_{1d}^{++}\) hits live \(1280/3\), \(61936/3\). Fail drop 1D; fail drop \(3k\); fail \(Q_{1d}^{\mathrm{sub}}\). \(16pA\) still only live at p=5,7. 15.546: p=7 \(|\mu_{\mathrm{full}}|\le(4p+1)/(15p^2)=29/735\), \(|\mu_{\mathrm{part}}|\le(p-2)/(3p^2)=5/147\). Mix \(109/2863<|T|=5/91\), sharp. Fail: \(|\mu_{\mathrm{full}}|=1/p^2\); fail drop 15 in the full majorant (bound \(283/2863>|T|\)). Aut_e \(G>T\) on the p=7 ensemble. \(S(\rho)\) does not pin the three full magnitudes. Open for \(p\ge11\). 15.544: p=5 \(|\mu_{\mathrm{full}}|=1/p^2\) on all 1800 \(|\kappa|=1\) (sign not \(\kappa/p^2\)). Mix \(|\mu|\le3/65<|T|\). Fail: drop full-Ω (\(|\mu|=T\)). p=7 \(\mu_{\mathrm{full}}\in\{13,17,29\}/(15p^2)\), not a \(p\)-law. Aut_e \(p\ge11\) open. 15.543: Type+ 1D 3-point \(\mu_{1d}=\kappa/(p(p-2))\) on \(|\kappa|=1\). Fail \(\kappa/p^2\); fail \(T\) at p=7. \(|\mu_{1d}|\le|T|\) for \(p\ge5\), equality iff \(p=5\). Does not close Aut_e (free \(\mu\) unnamed). \(A_{\mathrm{full}}\) still not a \(p\)-law. 15.542: ns \(\mu\)-half-net count equals \(|H_+|\) at p=5 and p=7 (130, 5726; converse of 15.305 C, not p=5-only). Fail \(n_{\mathrm{hn}}=n_{1d}\); fail extra half-nets at p=7. Hence \(n_{\mathrm{free}}=(n_{\mathrm{hn}}-n_{1d})/q=c_{\mathrm{eq}}(p-1)\) live. Not a theorem for \(p\ge11\). 15.541: \(c_{\mathrm{eq}}=\lfloor(4Ap(p^2-5)-n_{1d}(p^2-9))/(p^2(p-1)(p^2-9))\rfloor\). Fail ceil; fail drop \(p^2-9\) at p=7 (\(20\neq19\)). Names the Hoffman endpoint, not live \(c\) for \(p\ge11\). 15.539: \(n_{\mathrm{free}}=c(p-1)\Leftrightarrow D=D_{\mathrm{lattice}}(c)\). Live pin \(c=c_{\mathrm{eq}}\) only at p=5,7. Fail \(u=c\); fail \(c_{\min}\) at p=7 (\(108\neq114\)). Not a theorem for \(p\ge11\). 15.507 is \(p\equiv1\) only. Floor / Type I multi-level / residual (ii) \(k\ge4p\) still **OPEN**. 15.535: twisted \(\sigma(x)=a\bar x+b\), \(a=\pm v^{1-p}\), \(N(a)=1\), \(b+a\bar b=0\), freely pairs \(T_{\mathrm{ns}}\)-orbits at p=7 (all 14) and at no p=5 (0/10; fix\(\in\{2,30\}\), inv\(\in\{2,6\}\)). Fail: \(a=+\) free at p=5. \(D=n_{\mathrm{orb}}/2\) is a free \(\langle T_{\mathrm{ns}},\sigma\rangle\) count only for \(p\equiv3\). 15.533: Aut\(_\infty\) involutions do not freely pair \(T_{\mathrm{ns}}\)-orbits. \(n_{\mathrm{inv}}(x\mapsto-x)=2\cdot3^{C((p-1)/2,2)}\) (15.509 Fix; 6,54). Frob \(n_{\mathrm{inv}}=2,218\); Frob\(\circ(-\mathrm{id})=22,218\). Fail \(n_{\mathrm{inv}}=0\); fail \(n_{\mathrm{inv}}=n_{1d}/p\) at p=7. \(D=n_{\mathrm{orb}}/2\) is not a named pairing. 15.532: \(T_v\) is free on \(H_+\) iff \(\chi_q(v)=-1\). \(\mathrm{Fix}(T_v)=\binom{p}{m}\) on square lines, 0 on nonsquare lines (Type− ⊄ \(H_+\), 15.451 C). Fail: \(T_1\) free. Hence \(|H_+|=p\,n_{\mathrm{orb}}(T_{\mathrm{ns}})\) and live \(n_{\mathrm{orb}}=2D\). The pairing that names \(D=n_{\mathrm{orb}}/2\) is still open. 15.531: 15.527 axis-only \(n_R\) is not translation-invariant (p=7: 16 on min-key vs 20 on any member). Lin-form affine occupancy is all of \(n_{\mathrm{free}}\) (4=4, 114=114) by 15.305 C; fail \(n_X^{\mathrm{lin}}=98\). Random \(k\)-subsets are lin-affine at rate \(0.183\)/\(0.012\), not 1. Occupancy splits cannot name \(n_{\mathrm{free}}\). 15.527: free \(H_+\) orbits split affine-R plus leftover, \(n_{\mathrm{free}}=n_R+n_X=(4,0)\) at p=5 and \((16,98)\) at p=7. Fail: \(n_R=p-1\) at p=7; fail \(n_X=0\) at p=7. \(n_R,n_X\) unnamed. 15.528: leftover+s₊=2 at p=5 \(k=20\) is empty for nF∈[7,20] (HiGHS Infeasible, \(S\ge2\), 1739s, 493277 nodes). leftover-only nF=8 exists. Combined with nF=0..6, leftover+splus is empty for all nF at p=5 \(k=20\). Fail: leftover-only empty. even \(k>4p\) with far stays open. 15.526: \(\mathbb F_q^+\) splits \(H_+\) into Type+ 1D (size \(n_{1d}\), size-\(p\) orbits) and free NL orbits (size \(q\), \(n_{\mathrm{free}}=4,114\)). Hence \(D=n_{1d}/(2p)+(p/2)n_{\mathrm{free}}\). Fail: \(n_{\mathrm{free}}=0\); fail \(n_{\mathrm{free}}=p-1\) at p=7; fail \(|H_+|=n_{1d}\). \(n_{\mathrm{free}}\) unnamed. 15.525: type-index even-character Gram \(G_{\tau\sigma}=((q-1)/2)(|\tau\cap\sigma|+|\tau\cap(-\sigma)|)-2|\tau||\sigma|\) is named (certified p=5..19). Live \(\delta\) is not an eigenvector (p=5 merged wedge \(17792\neq0\)). Catalog p-laws miss \(D=N(a+bi)\) at both p=5,7 (\((p-1)/2\) names \((2,3)\); nothing names \((3,20)\)). \(|J(\chi,\psi)|^2=q\neq D\). Fail: drop the minus intersection; claim the p=5 wedge vanishes; claim \(D=\dim V_+\). Ensemble \(Q_\tau\) unnamed. \(\phi_F\) not imported. 15.540: \(A_{1d}(r)=-4p^3/(p-2)\) on \(r\in\mathbb F_p\) (\(1+r\neq0\)), else 0. Fail: same value off \(\mathbb F_p\). \(A_{\mathrm{part}}(\mathbb F_p)=0\), \(A_{\mathrm{part}}(\mathrm{off})=-2p^3/3\) (p=7). \(A_{\mathrm{full,dbl}}=0\) at p=5, \(-4p^3/15\) at p=7 — not a \(p\)-law. \(A_{\mathrm{free,dbl}}=-4p^3/19\) at p=7 is D-free as a fraction but 19 is unnamed. Still den 13,409. 15.538: \(\sum_{\xi\in\Omega}\psi(\alpha\xi)=(-1+\chi_\Omega\chi(\alpha)G)/2\) (\(\alpha\neq0\)). Type+ 1D \(A_{\mathrm{dbl}}=-4p^3/(p-2)\) (Johnson \(E[\hat\varepsilon(1)^2\hat\varepsilon(-2)]=-2(p+1)/(p-2)\)). Fail: drop \(\chi(\alpha)\); \(-4p^3/(p+1)\). \(A_{\mathrm{free,dbl}}=0\) at p=5, \(-1372/19\) at p=7, so \(A\) is not 1D-only. Still den 13,409. 15.536: on \(\Omega\)-triples \(I=p\,S(r/(1+r))\). Doubles are CM: \(I_{\mathrm{dbl}}=p\,S(-1)=p(2p-a_p^2)\), equals \(2p^2\) iff \(p\equiv3\pmod4\). Fail: \(I=GS\) (p=5 dbl \(-30\neq30\)); \(I_{\mathrm{dbl}}=2p^2\) at p=5 (\(30\neq50\)). Generic \(S\) on distinct triples is not a \(p\)-law (\(-10\), \(-2\), split \(\{-10,6,22\}\) at p=11). \(\Delta_{\mathrm{conn}}\) still has den 13,409. 15.534: \(I(\eta,\theta)=G(\chi)\chi(\theta)S(-\eta/\theta)\) for \(\theta\neq0\); \(I(\eta,0)=-\chi(\eta)G\); \(S(0)=S(1)=-1\). Fail: drop \(\chi(\theta)\); \(G\equiv p\); \(I\) a single \(p\)-law (\((1,1)\) is \(-30\) at \(p=5\), \(98=2p^2\) at \(p=7\)). 3-\(\chi\) types do not pin \(S\) (p=5 type \((-1,-1,+1)\) has \(S\in\{-6,2\}\)). Live \(\Delta_{\mathrm{conn}}\) still has den 13,409. 15.529: \(\Delta_{\mathrm{conn}}=2G/((p^2-1)p^6)\sum_{\Omega}AI\) equals \(-328/65\), \(-1144/2863\). Fail: drop \(q^{-3}\). A still has den 13,409. 15.523: \(\sum_{x\neq y\neq0}\chi(x)\chi(y)\kappa\psi_\xi(x-y)=3p\,\chi_p(-1)\chi(\xi)\). Fail: drop one pairing. `type_I_multilevel_bad_case_ND_closed` stays False. 15.522: at p=5,7 every off-pm1 15.290 type has \((Q_\tau/q^2)\cdot D=a^2+2b^2\) in \(\mathbb Z[\sqrt{-2}]\). Fail: the same in \(\mathbb Z[i]\) (p=5 \(++\) gives 48, not a sum of two squares). \((a,b)\) unnamed (`10p-46` interpolates \(a_{++}\) at \{5,7\}). Ensemble \(Q_\tau\) still unnamed as a Gauss/Jacobi formula. \(\phi_F\) not imported. 15.528: leftover+s₊=2 at p=5 \(k=20\) is empty for nF∈[7,20] (HiGHS Infeasible, \(S\ge2\), 1739s, 493277 nodes). leftover-only nF=8 exists. Combined with nF=0..6, leftover+splus is empty for all nF at p=5 \(k=20\). Fail: leftover-only empty. The 0.4s \(S\equiv2\) harvest is not this certificate. even \(k>4p\) with far stays open. 15.524: nF=7 Infeasible. 15.521: nF=4,5,6 Infeasible. 15.520: p=11 leftover-\(2p\) Fejer equals \(904/45\) by Johnson partition + 15.519. Fail: \(832/45\) (integer-hat2 rounding); fail: leftover-\(2p\) at p=7; fail: same constant at p=19. Ensemble \(Q_\tau\) unnamed. \(\phi_F\) not imported. 15.519: QR0 \(Q=2(p+1)\) on \(\mathbb F_p^\times\) for \(p\equiv3\pmod4\). 15.518: p=11 affine leftover includes stab \(=2p\). 15.517: p=5 \(H_+\) \(W\)-hist is the two-orbit law \(\{20,24,36,40\}=15:50:50:15\), so \(\mathrm{Var}(W)=\mathrm{Var}_{2\mathrm{orb}}=660/13\). Fail: same 4-point support or \(\mathrm{Var}=\mathrm{Var}_{2\mathrm{orb}}\) at p=7. Does not prove \(\mathrm{Var}\le\mathrm{Var}_{2\mathrm{orb}}\) for \(p\equiv1\). 15.516: ensemble \(L_{\mathrm{ns\_mix}}=\mu_+\mu_-\). \(Q_\tau\) unnamed. \(\phi_F\) not imported. Live `e1` is still the old AND. \(L=\tfrac12\) is **not settled**.

| Claim | Status | Reference |
|-------|--------|-----------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) | CLOSED | `solution.md` |
| \(\rho=1\) on Paley \(n=p^2+1\) | CLOSED | `PROOF_rho_eq_1.md` |
| Bi-tight majorization algebra (15.167) | CANDIDATE | `prop15167.py` |
| Type I freeness-fail ND (two-level Max−) | **CLOSED** (that slice) | 15.272 \(k=1\cup k=3\) \(\Rightarrow G_+\succ0\Rightarrow\) 15.249 dual-eq empty |
| Residual (ii), affine + even \(k\le4p-2\) | **CLOSED** (that range) | 15.179 + 15.236 + 15.237 |
| Residual (ii), even \(k\ge4p\) | **OPEN** | multi-level Max− leftover (15.274); p=5 \(k=20\) leftover+splus empty all nF (15.521/524/528); p=5 \(k=22\) leftover+splus empty nF=0,3–9,11–14 (15.547+15.552); nF=10 TLE; 15.560 nF=0 empty at \(k=26,28,30\) |
| Type I, multi-level Max− | **OPEN** | 15.543 \(\mu_{1d}\); 15.544 p=5 mix; 15.546 p=7 mix \(109/2863<|T|\); 15.551 Galois line-union; 15.559 Aut_e DEAD as \(A_{\mathrm{full}}\) name; `type_I_multilevel_bad_case_ND_closed=False` |
| \(\lambda_{\min}(\Phi)\ge6\) | **OPEN** | 15.279; 15.550–15.582: \(S=\mathrm{Kl}\), \(\mu_{1d}\), \(\mu_{k=3}\) named; \(p=5\) \(Q_{++}=48/13\); \(A_{\mathrm{full}}\) / \(D\) / ensemble \(Q_\tau\) unnamed; `phi_F_ge_6=False` |
| E(1) / \(L=\tfrac12\) | **OPEN** | needs the four leftovers above |
| Path-C residual / 16N | OPEN optional | not required for denseness path |

### Fatal gaps (honest)

1. Residual **(i)** two-level dual-eq: **CLOSED** by 15.272 (not Aut-Schur, not Gsum LB). Multi-level Type I is **OPEN**. `gsum_disj_lb` remains False and is unused.
2. Residual **(ii)** for affine + even \(k\le4p-2\): **CLOSED** (15.179 + 15.236 + 15.237). Residual (ii) for even \(k\ge4p\): **OPEN**. Do not call residual (ii) globally closed.

### Remainder progress (15.172–194)

- **15.179 residual (ii) affine branch CLOSED:** dual two-level freeness-fail affine \(\Rightarrow S_H\equiv3\Rightarrow k=3p-1\); impossible for \(k\ge3p\); fail-eq empty under bi-tight.  
- **15.193 residual (ii) exhaustiveness still False (not required):** freeness-fail does **not** force \(S\in\{2,4\}\) and \(f_e=3-S\). Affine + even \(k\le4p-2\) is ND-closed (15.179/236/237). Even \(k\ge4p\) is **OPEN**.
- **15.194 row negative-mass (residual i):** proved sufficient lemma \(N_e<4-6/p\Rightarrow\) dual-eq empty (independent box min). **Census:** pure \(N_e\) target **fails** at \(p=3,5,7\) (\(N_e=16/3,384/65,\approx13>\) thr). Refined **row+mass LP** max \(\kappa_e\approx1.19,1.17<\) need at \(p=5,7\) (blocks; evidence only). General Max+-free row+mass / ker-box still OPEN.  
- **15.195 mass-corrected dual-eq:** dual-eq needs mass-min \(\sum a\kappa\le-2(2-\alpha)\). Criterion proved. **Census:** p=5 exact min\(-30/13>-254/65\) blocks; p=7 blocks; p=3 does not. PSD/G+ floors too weak for worst-case \(a\). General Max+-free mass-min bound still OPEN.  
- **15.196 spectral row energy:** \(a\ge-2\) (PSD); \(Q_e\le 2(n+\lambda_2(n-2))/(n-1)-4\). **If \(Q_e\le10\)** (or \(\lambda_2\le6+5/(n-2)\)), mass-min blocks dual-eq for all primes \(5\le p\le47\) (Fraction check). **Census:** \(Q\approx8.17,6.69\) at p=5,7; \(\lambda_2=88/13,2160/409\). Spectral ub alone too weak at p=5. **OPEN:** Max+-free \(Q_e\le10\) or \(\lambda_2\) thr.  
- **15.197 min-distance + Q(K₄):** **Proved Max+-free** \(d_H\ge p+1\) on Max± (\(y\neq z\)) ⇒ \(|y\cdot z|\le(p-1)^2-2\). **Proved** \(Q_e=K_4/(n(n-1))-2(3n-2)/(n-1)\) with \(K_4=E[(y\cdot z)^4]\) on Max+². **Proved** \(Q_e\le10\Leftrightarrow K_4\le16n^2-14n\). Census K₄ below thr at p=5,7; above at p=3. **OPEN:** Max+-free \(K_4\le16n^2-14n\).  
- **15.198 Wick_hi residual-(i) reductions:** **Proved Max+-free** \(\mathrm{Wick}_{hi}=12n^2+48n\le16n^2-14n\) for all primes \(p\ge5\); \(K_4\le\mathrm{Wick}_{hi}\Rightarrow Q\le6+58/(n-1)<10\). **Proved** partition \(K_4=12n^2-48n+C_{\mathrm{diag}}+24\sum\eta_S^2\) with \(C_{\mathrm{diag}}=4n(11n-14)/p^2\). Exact K₄ at p=7: \(5218435600/167281\). Census \(K_4\le\mathrm{Wick}_{hi}\) at p=3(=),5,7. **OPEN:** Max+-free \(K_4\le\mathrm{Wick}_{hi}\) (same scalar as 15.95 \(\sum M^2\le\mathrm{Wick}_{hi}\)).  
- **15.199 frame identity:** **Proved** \(K_4=4n^2+E[\langle X_c,Z_c\rangle_F^2]\) with \(\|X_c\|_F^2=n(n-2)\) constant; Wick ⇔ \(E[R]\le8(n+6)/(n-2)\). **Proved dead:** λ_max-only majorant (need 32/3 at p=5 but λ_max(𝒞)≈13.54). **OPEN:** same K₄/|μ|/ker-box hinge.  
- **15.200 C−2/n∈ker(Gsum):** **Proved Max+-free** \(\kappa_{ij}=C_{ij}-2/n\in\ker(\mathrm{Gsum})\) via 15.189 sum identity (no Max+ enum). Dual-eq = free-\(e\) ker-box max \(\kappa_e\ge2-\alpha\). Census free-\(e\) LP: p=3 feasible (~2.3≥1.8), p=5 blocked (~0.81<1.95). **OPEN:** max \(\kappa_e<2-\alpha\) general (or ≤(3/2)·scheme).  
- **15.201 α(n−2) sufficient free-e bound:** **Proved** α(n−2)<2−α for all primes p≥5 (⇔ p(p²−3p+1)>0). So free-e max κ_e≤α(n−2)=2·scheme_max closes dual-eq. **Proved** |μ|≤2/n ⇒ Farkas (Gsum≥−4/n beats μ_*). Census: free-e max ≤α(n−2) at p=5,7 (not p=3); p=7 max|μ|=109/2863<2/n. **OPEN:** free-e max≤α(n−2) general (or |μ|≤2/n / K₄≤Wick_hi).  
- **15.202 free-e dual form; scheme⊕cross ⊆ ker:** **Proved** free-e dual S* form (lower-box); scheme⊕cross ⊆ ker(Gsum) dim n−1+n²/4. **Certified** ker=sc at p=3,5,7; free-e max 14/5, 369/455, 11736/19775. **OPEN:** free-e max≤α(n−2) general; ker=sc general.  
- **15.203 dual construction D(C):** **Proved Max+-free** algorithm (scheme dual → Comm∩zero-diag∩reg-deg → scale → nonneg repair); α(n−1)<2−α for p≥5. **Certified** cost_D<2−α at p=5..23 (blocks dual-eq **if** ker=sc). **OPEN:** cost_D(p)<2−α for all p≥5 + ker=sc (or |μ|/K₄). Residual (i)/E1/L still OPEN — no predicate flip.  
- **15.204 ker PSD characterization:** **Proved Max+-free** κ∈ker(Gsum)⇔ yᵀ(κ⊙C)y≡0 on Max± (via κᵀGκ=E q²). Free-e on ker ≤ free-e on sc (= iff ker=sc). **Proved** 8/p<2−α for p≥5 (so sum_ne≤(4/3)n sc-dual blocks when ker=sc). **Certified** S*_sc=123/7,3912/113 (p=5,7 dual LP); cost_D·p≈6.5–7.7 and sum_ne≤(4/3)n for p≤19; ker Q₊=scheme_image at p=3,5. **OPEN:** ker=sc general / S*≤n−2 general / cost_D≤8/p general / |μ| / K₄. No predicate flip.  
- **15.205 M_cand / free-e ratio thresholds:** **Proved Max+-free** M_cand=(p−2)/(p(2p+3))≤1/(2p) for p≥3 (⇒ |μ|≤M_cand closes residual-(i) Farkas). **Proved** r(p)=(7p²−3p+4)/(5p²−3p+2)≤3/2 for p≥5 (⇒ free-e≤r·scheme closes dual-eq with 15.192). **Certified** free-e/scheme=r(p) at p=5,7; |μ|≤M_cand at p=5 (sharp), p=7. **OPEN:** |μ|≤M_cand general (15.74 companion) or free-e≤r·scheme general or K₄/ker paths. No predicate flip.  
- **15.206 local n₃ degree:** **Proved Max+-free** (any conference C²=(n−1)I): n₃(a,S)=(n−6)/4 on |κ|=1 centres, (n−6)/4−1 on |κ|=3; Paley n₃=(p²−5)/4, n₁=(3p²−7)/4. Proof: pair sums ⇒ N(e)+N(−e)=(n−6)/4; centre matches ±e iff |κ|=3. **Certified** full p=3,5 and sample p=7. Crude n₃·O(1) **insufficient** for S₃≤joint_budget (signed cancel needed). Residual (i) still OPEN — no predicate flip.  
- **15.207 ker=sc reduction:** **Proved Max+-free** scheme ++ formula Vpᵀ(DfC+CDf)Vp=2p VpᵀDf Vp; scheme⊆ker Q₊; M⊥scheme⇔B constant diag (Tr0⇒zero diag); Wick E_Wick[q²]=8‖B‖_F² on zero-diag ++; **ker=sc ⇔ G₊≻0 on 𝒲₊₊⁰** (and Max− twin). **Certified** G₊ min eig 8 (p=3), 40/13 (p=5) on 𝒲₊₊⁰; ker Q₊=scheme at p=3,5. **OPEN:** G₊≻0 on 𝒲₊₊⁰ for all p≥5, then free-e/sc bound. No predicate flip.  
- **15.208 reverse degrees unconditional:** **Proved** d₃=p²−5 on |κ|=1 (from 15.206); reverse d₁⁽³⁾=3(p²−1), d₃⁽³⁾=p²−9; L_abs≤1/(2p). **Certified** ‖T‖₂=4p at p=5,7; actual μ satisfies master. **OPEN:** |μ|≤M_cand / free-e / K₄ / G₊ PD general. No predicate flip.  
- **15.209 G₊ spectrum p=7 + free-e star saturation:** **Proved** dim 𝒲₊₊⁰=n(n−6)/8; implication G₊≻0⇒ker=sc⇒(free-e_sc≤r·scheme)⇒dual-eq empty. **Certified** G₊≻0 on 𝒲₊₊⁰ at p=3,5,7 with min eig 8, 40/13, **1536/409**; full positive moment spectra; |Max+|=12,260,11452; free-e optimum saturates all stars at −α (p=5,7) and stars-fixed LP recovers free-e max=r·scheme. **OPEN:** G₊≻0 general (Fourier on zero-diag ++), or |μ|/K₄. Residual (i)/E1/L still OPEN — no predicate flip.  
- **15.210 D(C) cost budget:** **Proved** 8/p<2−α (recall); implication (ker=sc)∧(cost_D≤8/p)⇒dual-eq empty. **Certified** cost_D·p∈[6.56,7.68]<8 for p∈{5..23}; cost_D≤8/p and <2−α throughout. **OPEN:** cost_D≤8/p general + ker=sc, or |μ|/K₄. No predicate flip.  
- **15.211 G₊ PSD free + λ_* candidate:** **Proved** G₊≽0 on 𝒲₊₊⁰ (Gram restriction); PD ⇔ min Rayleigh>0; λ_*=8(n−6)/n>0 for p≥5 would prove PD; α(n−1) dual target. **Certified** min Rayleigh ≥ λ_* at p=3,5,7 (**sharp** at p=5: 80/13); Comm proj of scheme dual has regular degree and sum_ne→n−1 (diag repair is the only inflation). **OPEN:** E[q²]≥λ_*‖B‖_F² general, or |μ|/K₄. No predicate flip.  
- **15.212 Veronese spanning ⇔ G₊ PD:** **Proved** G₊≻0 on 𝒲₊₊⁰ ⇔ rank{yyᵀ−S:y∈Max+}=n(n−6)/8; Gram spectrum identity G/N = Rayleigh; dim const-diag ++ =1+D. **Certified** full Veronese rank at p=3,5,7 (65=65, 275=275) ⇒ ker=sc there. **OPEN:** rank=D for all p≥5 (PSL/Fourier), then free-e_sc bound. No predicate flip.  
- **15.213 Veronese/λ_* attack:** **Proved** Max+ is spherical 2-design in V₊ (E[yyᵀ]=2P₊); λ_*>0 sufficiency for ker=sc; **4-design shortcut dead** at p≥5 (p=5 Rayleigh has 3 distinct eigs {80,144,176}/13). sc-dual sum_ne≤(4/3)n ⇒ free-e_sc≤8/p<2−α (recall). **Certified** Op=E[q²]−λ_*‖B‖²≽0 at p=3,5,7 (sharp ker dim n at p=5); D(C) sum_ne/n≤1.28<4/3 on p≤19. **OPEN:** E[q²]≥λ_* general or rank=D general (character-sum full Max+ blocked by multi-orbit). No predicate flip.  
- **15.214 Master/δ structure:** **Proved** T=Tᵀ; master solutions μ=μ_mn+δ with μ_mn=16(16p²I−T²)^+κ, δ∈E_{±4p}; **κ⊥E_{±4p}** (compatibility); resolvent form when ±4p∉σ(T) forces δ=0; matching Gsum PSD ⇒ |μ|≤1/2 (weak); Parseval ∑m_A²=2ⁿ/|Max+|. Particular majorant ≤1/(2p) recalled. **OPEN:** control δ on |κ|=1 so |μ|≤1/(2p) for all p≥5 (δ≠0 at p=5; pointwise |μ|≰|μ_part| at p=7). No predicate flip.  
- **15.215 K₄/Wick Tκ–η_part:** **Proved Max+-free (any conference):** star values on |κ|=1,3 (64-exhaust); ∑κ²=n(n−1)(n−2)(n−5)/8; n₃=n(n−1)(n−2)(n−6)/96; ‖Tκ‖₂²=6n(n−1)(n−2)(n−6); **Tφ=(n+2)star**; **T³κ=4(n+14)Tκ** (⇒ Tκ on λ²=4(p²+15) for Paley). Min-norm η_part has closed norm 5(p²−1)(p²+1)(p²+3)/(6p²(p²−5)) **strictly below** B_wick for p≥5. Full η=η_part+δ with δ∈E_{4p}.
- **15.216 Residual (i) K₄ thr path OPEN (fatal gap):** **Proved** ∑star κ=0; κ∈V_λ⊕ker(T); R_ke=128p/(3p²+17)≤2p; crude≤thr-η-budget; Q≤10⇒dual-eq empty (mass-min). **FATAL GAP (old):** R_ke≤2p⇏R≤2p via naive convex combo (AI-test BLOCK).  
- **15.217 Φ identity repairs Rayleigh criterion:** **Proved Max+-free** Φ(y)=∑κ∏y = n(n−1)(n−2)/8 constant on Max+ (Q²=2n(n−1)+8Φ; r=3 dies by C²/row-sum). Hence ⟨m₄,κ⟩=Φ and **R≤2p ⇔ ‖m₄‖₂² ≤ n(n−2)/4 ⇔ ‖δ‖₂² ≤ (p−1)(p+1)(p²+1)(3p²−47)/(24(p²−5))**. Conditional: that bound ⇒ η≤crude ⇒ K₄≤thr ⇒ dual-eq empty. **OPEN:** ‖m₄‖₂² ≤ n(n−2)/4 for all p≥5 (census OK at p=5: ≈143≤156). Path-C δ²≤room_hyp/24 is sufficient but open. Residual (i)/type_I/E1/L remain **OPEN**.  
- **15.218–224 structure (2026-08-11/12):** Max+⊥Max−; Gsum=2m₁μ₄; Tr(Gsum²); free-e row+ones dual **if** G_min≥−1/p; R-matrix; m4₂↔K₄; D-energy; ⟨D,TD⟩+CS ‖δ‖₂²≤k₂/|Max+|; gap_wick closed + Path C⇒Wick; T²φ=4(n+2)(φ−2κ) Max+-free; μ_part solves master always. Predicates stay **False**.  
- **15.225 constrained P± / μ_part=μ_mn (2026-08-12):** **Proved Max+-free** P±=(4pI±T)/(8p) on E_{±4p}; unique δ=δ₊+δ₋ split; **φ⊥E_{±4p}** (and star⊥E) via Tφ identity + κ⊥E; **μ_part=μ_mn** for all primes p≥5 (upgrades LSQR census); φ-averaged particular ≤1/(2p); room_δ=(p³−2p²−13p+18)/(2p²(p²−5))>0 with |δ|≤room_δ ⇒ |μ|≤1/(2p). **OPEN:** |δ_S|≤room_δ on |κ|=1 (or Aut-SOS / Path C / k₂ gap / K₄ / free-e+ker=sc). Residual (i)/E1/L still **OPEN** — no predicate flip. Evidence: `src/e1_gmin_m4_prop15225.py`, `evidence/e1_gmin_m4_prop15225.json`.  
- **15.226 Farkas-sharp thr / G-δ / cross / ‖μ_part‖₂² (2026-08-12):** **Proved** μ_f=(2p−3)/(p(3p−2))=−μ_*/2 with 1/(2p)<μ_f; maj≤μ_f so room_δ^f>room_δ>0; **∑ m₄⁺ m₄⁻ = n(n−2)/8** (y·z=0⇒e₄); ∑μ² identity; **δ∈E_{±4p}^G** (G-invariant, no dim≤1); **‖μ_part‖₂²=(p−1)(p+1)(p²+1)(3p²+1)/(24(p²−5))** closed. **OPEN:** |μ|≤μ_f or Aut-SOS on E^G for |δ|≤room_δ^f. No predicate flip. Evidence: `src/e1_gmin_m4_prop15226.py`.  
- **15.227 free-e dual sharp at μ_* (2026-08-12):** **Proved** at G_min=μ_* the row+ones dual has a=p(3p−2)/(6(p²−1)), sum_ne=p(p²+1)/3−1, **cost≡need** (zero slack / Farkas-calibrated). At G_min=−1/p: **need−cost=(p−4)/(2p+1)>0**. Monotone: cost(G)<need iff G>μ_*. One-sided form: residual-(i) ⇔ max|μ₄|≤1/(2p) (strict sublevel of max|μ|<μ_f). Aut-SOS on E^G restated; bound **OPEN**. No predicate flip. Evidence: `src/e1_gmin_m4_prop15227.py`.  
- **15.228 boolean diamond + Tμ_part (2026-08-12):** **Proved** |μ|+|ν|≤1 (ν=Tμ/(4p)); **Tμ_part=−8/(p²−5)·star**; on |κ|=1 star=0⇒ν=δ_{+}−δ_{-} and diamond |μ_part+δ|+|ν|≤1; L² orth ⟨Tμ_part,Tδ⟩=0; diamond L² ‖δ‖ bound (weak). Aut-SOS+diamond structure; bound **OPEN**. No predicate flip. Evidence: `src/e1_gmin_m4_prop15228.py`.  
- **15.229 size-3 Cy vanishes + functional eq (2026-08-12):** **Proved Max+-free** size-3 contribution to E[∏(Cy)_i] is 0 on |κ|=1 (C² collision calculus + 64-exhaust s=0); **(p⁴−1)μ+2φ=R̄₄** with R̄₄ the size-4 non-self-image remainder; conditional |R̄₄|≤(p⁴−1)/(2p)−4(p−2)⇒|μ|≤1/(2p). **OPEN:** bound R̄₄. No predicate flip. Evidence: `src/e1_gmin_m4_prop15229.py`.  
- **15.230 R_part≤budget + Cy≡δ unify (2026-08-12):** **Proved** R_part=(p⁴−1)μ_part+2φ has |R_part|≤R_part_max with closed form; **R_part_max≤B** for all primes p≥5 with slack (p−1)(p+1)(p²+1)(p³−2p²−13p+18)/(2p²(p²−5))>0 (same cubic as maj≤1/(2p)); **R̄₄=R_part+(p⁴−1)δ** and (B−R_part_max)/(p⁴−1)=room_δ. Particular Cy image is R-safe; residual-(i) = δ-correction only. **OPEN:** |R̄₄|≤B / |δ|≤room_δ. No predicate flip. Evidence: `src/e1_gmin_m4_prop15230.py`.  
- **15.231 permanent form of R̄₄; crude |per| dead (2026-08-12):** **Proved Max+-free** R̄₄(S)=∑_{T≠S} per(C[S,T]) μ_T=(Per μ)_S−μ_S on |κ|=1 (size-4 Cy = permanent operator). **Proved** crude ∑|per| majorant is dead for residual-(i) (p=5: 47500≫B=50.4; all p≥5: binom(n,4)·24≫B) — signed cancellation required (Jacobi/cycle before abs or Aut-SOS). Certified only: Per φ=(2n+1)φ on |κ|=1 for p=3,5,7 (not general). **OPEN:** |R̄₄|≤B. No predicate flip. Evidence: `src/e1_gmin_m4_prop15231.py`.
- **15.232 intersection split / k=3 per + C² pairing (2026-08-12):** **Proved Max+-free** R̄₄=∑_{k=0}^{3} R^{(k)}; k=3 permanent 11-term closed form (column-linear in the new vertex, |per|≤11); C² pairing ∑_{t∉S} C_{ut}C_{vt}=(n−1)1_{u=v}−∑_{s∈S} C_{su}C_{sv}. **Proved dead:** unsigned layers k=0,1,2 each exceed B for all primes p≥5. **Proved:** k=3 unsigned |R^{(3)}|≤44(p²−3)≤B for all primes p≥89 (f(p)=p⁴−88p³−8p²+280p−1; f(89)>0, f'>0); unsigned k=3 dead for 5≤p<89. Residual (i)/E1/L still **OPEN** — no predicate flip. Evidence: `src/e1_gmin_m4_prop15232.py`.
- **15.233 signed k=2 (2026-08-12):** **Proved Max+-free** 14 surviving perms, |per|≤14; per2 / κ grouped closed form; bilinear in the two new columns (C⊗C contraction). Unsigned |R^{(2)}|≤42(p²−3)(p²−4) still exceeds B for all p≥5. No predicate flip. Evidence: `src/e1_gmin_m4_prop15233.py`.
- **15.234 signed k=1 (2026-08-12):** **Proved Max+-free** 18 survivors, |per|≤18; Laplace along the shared column per=∑_{x≠a} C_{xa} per_3(C[S\{x},new]); trilinear in the three new vertices. Unsigned |R^{(1)}|≤12(p²−3)(p²−4)(p²−5) exceeds B for all p≥5. No predicate flip. Evidence: `src/e1_gmin_m4_prop15234.py`.
- **15.235 signed k=0 (2026-08-12):** **Proved Max+-free** full 24-perm permanent; S₄ cycle-type split 1+6+3+8+6; 4-cycle inverse pairing (P_4cyc even, |P_4cyc|≤6). Unsigned 24·C(n−4,4) still exceeds B for all p≥5 (n=p²+1). Layer types named; residual (i)/E1/L still **OPEN**. Evidence: `src/e1_gmin_m4_prop15235.py`.
- **15.236 residual (ii-b) ND CLOSED (2026-08-13):** Even \(k\le4p-2\) ⇒ max_Max− \(S\ge-2\); dichotomy weak-ND or dual-bad two-level; dual-bad empty by 15.50 slopes + \(v_p(|Max_\pm|)=1\). `residual_ii_b_ND_closed=True`. Evidence: `src/e1_gmin_m4_prop15236.py`.
- **15.237 residual (ii-a) ND CLOSED (2026-08-13):** Dual-bad \(U=\{S=-4,f_e=-1\}\) is a 0-1 pair-span; L²=L ⇒ star or triangle; constants / pair-slices / triangle 3-equals cannot be U. `residual_ii_a_ND_closed=True`. `residual_ii_full_closed` = affine ∧ (ii-a) ∧ (ii-b). Residual (i) / E1 / L still **OPEN**. Evidence: `src/e1_gmin_m4_prop15237.py`.  
- **15.238 Per-eigenrelations ⇒ Per μ_part Cy-FE (2026-08-13):** **Proved Max+-free** coefficient identity: if Per κ=p⁴κ−6φ and Per φ=(2n+1)φ, then Per μ_part=p⁴μ_part+2φ. Superseded as unconditional by 15.239. Evidence: `src/e1_gmin_m4_prop15238.py`.  
- **15.239 Per φ / Per κ Max+-free eigenforms (2026-08-13):** **Proved** for any symmetric conference: Per φ=(2n+1)φ and Per κ=(n−1)²κ−6φ on every 4-set (inj+IE / matching collision calculus + C²). Paley: (n−1)²=p⁴ ⇒ **Per μ_part=p⁴μ_part+2φ** unconditional; on |κ|=1 Cy-FE ⇒ Per δ=p⁴δ (structure). **OPEN:** |μ|≤1/(2p). No predicate flip. Evidence: `src/e1_gmin_m4_prop15239.py`.  
- **15.240 maj≤2/n + envelope criterion (2026-08-13):** **Proved** maj≤2/n for all primes p≥5 (cubic p³−3p²−5p−9). **Proved** criterion: envelope |μ|≤max(|μ_part|,|f4|) on |κ|=1 ⇒ |μ|≤2/n≤1/(2p) ⇒ residual-(i) Farkas. Free-e_sc budget recall (needs ker=sc). **OPEN:** envelope hypothesis / ker=sc (λ_*) / M_cand. No predicate flip. Evidence: `src/e1_gmin_m4_prop15240.py`.  
- **15.241 hull criterion + halfspace dead (2026-08-13):** **Proved** hull criterion: min(0,μ_part,f4)≤μ≤max(0,μ_part,f4) on |κ|=1 ⇒ residual-(i) (⇔ max-abs envelope when μ_part f4≤0; strictly stronger same-sign coherence when μ_part f4>0). **Proved** segment convexity on [μ_part,f4]. **Proved dead:** halfspace-only character sums (m4_H exceeds 1/(2p) at p=5,7; full Max+ needs non-halfspace cancellation). **OPEN:** hull/envelope hypothesis general. No predicate flip. Evidence: `src/e1_gmin_m4_prop15241.py`.  
- **15.242 Rayleigh spectrum p=3,5,7 (2026-08-13):** **Proved** λ_* sufficiency recall + dim 𝒲₊₊⁰ + Wick residual form. **Certified** full E[q²] spectrum on 𝒲₊₊⁰: p=3 single 16; p=5 eigenvalues {80,144,176}/13 with mult {n,n,d}, min=λ_* sharp; p=7 eigenvalues k/409 (409=|Max+|/(4p)) with mult {n,2n,n,n,d}, min=3072/409>λ_*. Min eigenspace mult=n at p=5,7. **OPEN:** E[q²]≥λ_* for all p≥5 (Aut-isotypic / orbit SOS). No predicate flip. Evidence: `src/e1_gmin_m4_prop15242.py`.  
- **15.243 ∑κ_C κ_B identity proved (2026-08-13):** **Proved Max+-free** (any conference, zero-diag B=P₊BP₊): ∑_S κ_C(S)κ_B(S)=(n+1)/4‖B‖_F² (parallel=1/4 via W=C⊙B row-sums; cross=n/4 via S_dist=n‖B‖² from Tr(CBCB)−corr). Upgrades 15.89.2 census→proved. **Proved** E[q²]=(8+4/p²)‖B‖²+8∑ρ κ_B and λ_* ⇔ ∑ρ κ_B≥−6/n−1/(2p²) (sharp p=5). **OPEN:** ρ residual lower bound. No predicate flip. Evidence: `src/e1_gmin_m4_prop15243.py`.  
- **15.244 ∑φ κ_B = −n/4 + μ_part rewrite (2026-08-13):** **Proved Max+-free** ∑_S φ(S)κ_B(S)=−(n/4)‖B‖_F² (per-star Z=xxᵀ⊙B, x=C_r·; σ₁=0 via CB=pB; Parallel_r=‖B‖²/4−(n/2)ρ_r; sum_r). **Proved** ∑ μ_part κ_B=(p²+3)/(4(p²−5))‖B‖² and γ(p)=∑(μ_part−κ_C/p²)κ_B/‖B‖²=(3p²+5)/(2p²(p²−5))>0. **Proved** λ_* ⇔ ∑(m4−μ_part)κ_B ≥ −6/n−1/(2p²)−γ(p). μ_part already overshoots Wick; Max+ residual m4−μ_part must not be too negative. **OPEN:** that residual (or hull/|μ|/K₄). No predicate flip. Evidence: `src/e1_gmin_m4_prop15244.py`.  
- **15.245 Z-frame / Op average (2026-08-13):** **Proved Max+-free** Z_y=yyᵀ−2P₊∈𝒲₊₊⁰ with ‖Z‖_F²=n(n−2), ⟨Z_y,Z_z⟩=(y·z)²−2n; yᵀBy=⟨B,Z_y⟩; Op=E[Z⊗Z]; scheme-image = ker Op on traceless ++ and traceless++=scheme⊕𝒲₊₊⁰ orthogonally; Tr(Op)=n(n−2), average Rayleigh=8(n−2)/(n−6)≥λ_* for p≥5. **OPEN:** Op≽λ_*I on 𝒲₊₊⁰ (CS on ρ dead; 4-design dead; s_max K₄ bound too crude). No predicate flip. Evidence: `src/e1_gmin_m4_prop15245.py`.  
- **15.246 edge Cov=Op/2 + R-path budget (2026-08-13):** **Proved Max+-free** edge features f_{ij}=y_i y_j: ‖f‖₂²=|E|, 1ᵀf=p, μ=C_edge/p, ‖f−μ‖₂²=n(n−2)/2; Z lifts f−μ; Rayleigh_Op=2·Rayleigh_Cov on 𝒲₊₊⁰ so Op≽λ_* ⇔ Cov≽4(n−6)/n I. **Proved** 15.217 ∑η²≤n(n−2)(n−5)/(8(n−1)) is **strictly weaker** than Wick_hi eta budget for all p≥5. **OPEN:** Cov gap / ∑η² bound / |μ|. No predicate flip. Evidence: `src/e1_gmin_m4_prop15246.py`.  
- **15.247 m4_part + room_δ^R (2026-08-13):** **Proved Max+-free** m4_part=aκ+bφ+z star solves (4pI−T)m=4κ/p in span{κ,φ,star} with a,b as μ_part, z=−2p/D; ‖m4_part‖₂²=‖μ_part‖₂²+z²‖star‖₂²; R≤2p ⇔ ‖m4−m4_part‖₂²≤room_δ^R=delta_room_for_R>0. **Proved** α·(3/2)(n−1)<2−α (dual budget). **OPEN:** ‖δ‖₂²≤room_δ^R. Evidence: `src/e1_gmin_m4_prop15247.py`.  
- **15.248 Comm scheme-dual We/sum_ne closed forms (2026-08-13):** **Proved Max+-free** after Comm-proj of scheme dual: We=½−1/(2p²(p²−2)); sum_ne^Comm=(n/2)/We−1=p²(p²+1)(p²−2)/(p⁴−2p²−1)−1 < (3/2)(n−1) for all primes p≥5 (restricted triple sum −C_e(n−2); Sp expansion). **OPEN:** full D(C) inflation (LS+nonneg) still ≤(3/2)(n−1), and ker=sc. Evidence: `src/e1_gmin_m4_prop15248.py`.




- Avg disj Gsum \(=2/(n-3)\); \(G_0\) PSD; **15.176** μ_* / −1/p sufficiency.  
- **15.177–178:** |μ₄| hinge form; star identity; dual-eq \(n_d\le1\) kill; \(n_d=2\) wedge kill for \(p\ge7\).  
- **15.180:** dual-eq Q_pairs=30−6p−24/p<0; open dual-eq core is \(n_d\ge2\) after PSD/score filters.  
- **15.181:** vertex-star \(\sum_{e\ni i}f_e=\pm p\) ⇒ \(\sum_{e\ni i}\mathrm{Gsum}_{ef}=2\); Max+-free κ-counts through edge \(n_3=(n-2)(n-6)/8\), \(n_1=3(n-2)^2/8\); matching \(n_d=2\) PSD floor −4 (harder than wedge); p=5 dual-eq k-sparse linear box empty (census).  
- **15.182:** dual-equality normal form \(x=\alpha\mathbf{1}-2e_*+\kappa\) with \(\alpha=6/(pn)\), \(\kappa\in\ker(G_+)\cap\ker(G_-)\), \(\kappa_e=2-\alpha\), box on \(\alpha+\kappa\); particular solution Max+-free; p=5 ker-box LP infeasible (census). General ker-box obstruction open.  
- **15.183:** Max+⊥Max− (any symmetric conference); \((G_+G_-)_{ee}=n/(2p^2)\); ker(Gsum)⊥1 automatic; matching \(n_d=2\) PSD \(a+b\ge-2\sqrt{2+c}\); binary dual-eq form. p=5 max \(\kappa_e\) under ker+box ≈0.811<2−α.  
- **15.184:** Max+-free \(T^2\kappa=-24\varphi+48\kappa\) on \(|\kappa|=1\) (C² reduction + 48-labeling); \(|T^2\kappa|=24|\varphi-2\kappa|\).  
- **15.185–187:** Paley \(|\varphi|\le2(p-2)\) all odd primes (Auer–Top supersingular\(\Rightarrow\)double residue + Hasse ladder). Global \(T^2\kappa=-24\varphi+48\kappa\) on **all** 4-sets (any conference, 64-labeling). \(T^2\varphi=4(p^2+3)(\varphi-2\kappa)\) (Paley census p=3,5,7). \(\mu_{\mathrm{part}}=[(p^2-1)\kappa-2\varphi]/(p^2(p^2-5))\) solves master when \(T^2\varphi\) holds; majorant \(\le1/(2p)\) for \(p\ge5\).  
- **15.188:** Target correction — \(|\mu_{\mathrm{actual}}|\not\le|\mu_{\mathrm{part}}|\) pointwise (p=7 census: ~29k violations). Viable sufficient bound \(|\mu|\le2/n=2/(p^2+1)\) (\(\le1/(2p)\) for \(p\ge5\)). p=5: actual\(=(4\kappa-\varphi)/(pn)\), max \(3/65\). p=7: max \(109/2863<1/14\) and \(<2/n\).  
- **15.189 (Max+-free):** \(1^\top y=(p+1)y_\infty\) on Max+; \(E_\pm[y_iy_j]=\pm C_{ij}/p\); tight frame; adjacent Gsum=0 from π; G+ 6×6 PSD ⇒ \(|\mu|\le1-2/p\) on \(|\kappa|=1\) (**too weak** for residual i).  
- **15.190 (Max+-free):** scheme-ker \(\kappa_{ij}=f_i+f_j\) (\(\sum f=0\)) lies in \(\ker(Gsum)\); scheme-ker max \(\kappa_e=\alpha(n-2)/2=3(n-2)/(pn)<2-\alpha\) for all \(p\ge3\). Full dual-eq ker-box empty at p=5,7 (census: max \(\kappa_e\approx0.811,0.593<2-\alpha\)); **feasible at p=3**. 3-point moments vanish (certified p=3,5).  
- **15.191 (Max+-free partial):** Derangement permanent of \(C[S,S]\) equals 1 on \(|\kappa|=1\) (64-exhaust); star-sum \(\sum_s\prod_{i\neq s}C_{is}=0\) on \(|\kappa|=1\); Cy-expansion size1+size2 \(=-2\varphi\) (any conference \(C^2\) + Paley \(\pi\)); envelope \(|4\kappa-\varphi|/(pn)\le2/n\le1/(2p)\) for \(p\ge5\). **Correction:** \(|\mu|\le|f_4|\) fails at p=7 (many classes; f4 not a pointwise majorant); viable target remains \(|\mu|\le2/n\) (census p=5,7) or \(\le1/(2p)\).  
- **15.192 (Max+-free):** Gsum diag\(=2\); row sum\(=n\); avg disj Gsum\(=2/(n-3)\). Aut_e averaging: dual-eq feasible iff Aut_e-invariant dual-eq feasible. \((3/2)\cdot\)scheme-max \(<2-\alpha\) for all \(p\ge5\). **Census Aut_e ker-box:** p=3 feasible (\(\max\kappa_e=14/5\)); p=5,7 empty (\(\max=369/455\), \(11736/19775\); ratios to scheme \(41/28\), \(163/113\), both \(<3/2\)).  
- **Still open (residual i):** Max+-free \(\max\kappa_e\le(3/2)\cdot\)scheme-max (or any bound \(<2-\alpha\)) for all \(p\ge5\), **or** \(|\mu|\le2/n\) (or \(\le1/(2p)\)) on \(|\kappa|=1\).  
- `gsum_disj_lb_proved_general()=false`; residual (ii) closed only for affine + \(k\le4p-2\); **E1/L OPEN**.

### Short package

`evidence/share/denseness_path_package.md`

### Required opens

1. \(\lambda_{\min}(\Phi)\ge6\) without treating \(G_{u,\mathrm{disj}}\) as a Gram (it is not PSD).
2. Residual (ii) for even \(k\ge4p\).
3. Type I when Max− is not two-level.
4. Lemma D: \(k=3\) existence and the 2-plane amplitude model.

**Non-required:** Path-C / 16N / Hypothesis H / 15.193 exhaustiveness.

**Current:** 15.272 span written; E(1)/\(L\) not settled. See the denseness-path package caveats.

**Residual-(i) attack (2026-08-13 late, not a close):** exact 4-point sums at \(p=3,5,7\) (\(\mu=\mu_{\mathrm{part}}\), \(f_4\), linear+CR-split); \(\nu\equiv0\); envelope holds. New dead: 4×4 Gram; affine-quadratic level sets; CR-class master+diamond LP (\(\max|\mu|=1\)); IP-valency \(K_4\) (not regular at \(p=7\)); global \(t(p)\) mix. Viable: general envelope, \(F(\lambda)\) matching those sums, \(G_+\succ0\) then Comm-repair dual, \(|\mu|\le2/n\), \(K_4\le n(15n-22)\), non-low-degree dual-eq Farkas.

### 15.249 (this resume) — algebraic free-e dual D_alg
- **Proved Max+-free (Paley/Weil):** cost_D < 2−α for free-e over scheme⊕cross for all primes p≥5 via Comm+Comm(diag) repair + Weil |Q|≤2p far bound.
- Closed forms: We_alg, sum_ne¹, stars>0; m_p≥1−2p; t_ub=2(2p−1)/den.
- **OPEN:** ker=sc (or ‖δ‖₂²≤room / |μ|). Predicates still **False**.

### 15.250 — R-path ⇔ E[s⁴]≤15n²−22n
- Odd Max+ moments vanish; Es4 expansion Max+-free.
- R≤2p ⇔ fourth-moment bound on Max+. Census p=5,7 holds.
- **OPEN:** prove Es4 bound general. Predicates **False**.

### 15.251 — Cy-identity; |μ₄|≤2/n path
- (p⁴−1)m₄+2φ=Ext on |κ|=1 Max+-free structure.
- |m₄⁺|,|m₄⁻|≤2/n ⇒ residual-i. Census p=5,7 holds (m₄⁺=m₄⁻).
- **OPEN:** |μ₄|≤2/n general. Predicates **False**.

### 15.252 — Master/T²/Ext residual-(i) criteria
- **Proved Max+-free:** T²μ=16(p²μ−κ) pointwise (master rewrite); |ρ|≤L_abs=(p−2)/(2p²)⇒|μ|≤1/(2p); |Ext|≤2p²−4p+6 (under |φ|≤2(p−2))⇒|μ|≤2/n.
- Triangle on |16κ+T²μ| alone is **vacuous** for bounding |μ| (equals 16p²|μ|).
- **Census p=5,7:** envelope OK; Ext≤uniform maj (p=7: 71.37≤76); |ρ|≤L_abs (p=5: 0.055≤0.060 sharpish); |μ|≰|f4| at p=7 (117600 viol) while max-abs envelope holds.
- **OPEN:** envelope / Ext maj / |ρ|≤L_abs general (or Es4 / ker=sc). Predicates **False**.

### 15.253 — Wick-reflection residual-(i) criteria (preferred)
- **Proved Max+-free:** |ρ_f4|≤L_abs for all primes p≥5 (maj (3p²−8p+1)/(p²n); cubic g=p³−8p²+17p−4>0).
- **Proved criterion:** |ρ|≤|ρ_f4| ⇔ μ∈[f4, f4^♯] with f4^♯=2κ/p²−f4 (Wick reflection) ⇒ residual-(i) Farkas.
- **Proved:** hull ⇒ |ρ|≤L_abs via triple majorant max(1/p²,|ρ_part|,|ρ_f4|).
- **Census p=5,7:** reflection viol=0; t=(μ−f4)/ρ_f4 ∈[−1.70,−0.10]⊂[−2,0]. Pure |μ|≤|f4| remains false at p=7.
- **OPEN:** prove reflection hyp |ρ|≤|ρ_f4| (t∈[−2,0]) general. Predicates **False**.


### 15.254 — T m₄± closed forms; Paley C∼−C
- **Proved Max+-free (any conference with π):** (Tm₄⁺)=4p m₄⁺−4κ/p; (Tm₄⁻)=−4p m₄⁻+4κ/p; Tμ=2p(m₄⁺−m₄⁻); ν=½(m₄⁺−m₄⁻).
- **Proved Paley:** D Pᵀ C P D=−C (P=mult-by-nonsquare, D_∞=−1); Max−=monomial image of Max+; μ=½[m₄⁺+χ_D m₄⁺∘π⁻¹].
- On |κ|=1: m₄⁺=m₄⁻ ⇔ m₄⁺(S)=χ_D(S)m₄⁺(π⁻¹S); then Tμ=ν=0.
- **(κ,φ) under π:** ∞∉S ⇒ (κ,φ) fixed; ∞∈S ⇒ (κ,φ)→(−κ,−φ); f4(S)=χ_D f4(π⁻¹S).
- **Census p=5,7:** T formulas exact; m₄⁺=m₄⁻ and Tμ=0 on all |κ|=1.
- **OPEN:** m₄⁺=m₄⁻ general on |κ|=1; reflection (15.253). Predicates **False**.

### 15.255 — m₄± ∞-expansions; m₄-eq ⇔ vanishing ∑μ
- **Proved Max+-free:** m₄⁺(S)=κ/p²+T₊/p, m₄⁻(S)=κ/p²−T₋/p on S∋∞; m₄⁺=m₄⁻ ⇔ ∑_{j∉S}μ({j,a,b,c})=0; Aut 3-transitive reduces to triples.
- Reflection on ∞-sets ⇔ |∑μ|≤p|ρ_f4| when equality holds.
- **Census p=5,7:** expansions exact; ∑μ=0. Predicates **False**.

### 15.256 — Global ∑μ and extension (κ,φ) sums
- **Proved Max+-free:** ∑_S μ(S)=−p²(p²−1)/12 (via E[Q²]=p²=binom(n,2)+6∑μ; equiv. E[(1ᵀy)⁴]=n+3n(n−1)+24∑μ with E=p⁴+6p²+1).
- **Proved (normalised Paley):** on S={∞,a,b,c}, ∑_{j∉S}κ(T_j)=3−κ0² (=2 on |κ0|=1); ∑_{j∉S}φ(T_j)=p²+P with P=(κ0²−3)/2 (=p²−1 on |κ0|=1).
- **Proved structure:** f4-extension sum (9−p²)/(pn)≠0 for p≥5 — local vanishing needs |κ|=3 cancellation (cannot hold if μ=f4 on all extensions).
- **Census p=5,7:** ∑μ and extension sums exact (1800 / 14112 triples).
- **OPEN:** vanishing ∑μ / reflection / residual-(i). Predicates **False**. Global ∑μ and (κ,φ) extension sums alone do not force local vanishing.

### 15.257 — Seidel 4×4 spectrum / resolvent; Schur setup
- **Proved Max+-free (all 4×4 Seidel):** |κ|=1 ⇒ σ={±1,±√5}, χ=λ⁴−6λ²+5; |κ|=3 ⇒ {3,−1³} or {−3,1³}.
- **Proved:** on |κ|=1, det(pI−C_S)=(p²−1)(p²−5); resolvent (pI−A)^{-1}=[p³I+p²A+p(A²−6I)+(A³−6A)]/det; same det for −pI−A.
- **Structure:** Schur y_S=(±pI−C_S)^{-1}e ⇒ ∏y_S multilin in external field; m₄⁺−m₄⁻ = difference of Max± averages (OPEN identity on |κ|=1).
- **OPEN:** Schur average identity / reflection / residual-(i). Predicates **False**.

### 15.258 — Orientation theorem C_S ∼_sp −C_S ⇔ |κ|=1
- **Proved Max+-free:** A ∼_sp −A (signed permutation) ⇔ |κ(A)|=1.
  - |κ|=3: σ(A)≠σ(−A) ⇒ not similar.
  - |κ|=1: exhaustive conjugacy, always with det D=+1 (4-product invariant under local conjugator).
- Matches census: ν=0 exactly on |κ|=1, nonzero on |κ|=3 (p=3,5,7).
- **OPEN:** lift local ∼_sp to m₄⁺=m₄⁻ (permutation need not extend to Aut(C); pure switching insufficient). Reflection / residual-(i) still OPEN. Predicates **False**.

### 15.259 — Transport lemma; π-covariance reduction
- **Proved Max+-free:** Local orientation extends to C'=D Pᵀ C P D with C'_S=−C_S; Φ(y)=D Pᵀ y bijection Max⁺(C)→Max⁺(C'); det D=+1 ⇒ m₄⁺(C,S)=m₄⁺(C',S).
- **Proved:** m₄⁺=m₄⁻ on |κ|=1 ⇔ m₄⁺(S)=χ_D(S) m₄⁺(π⁻¹S) (finite: m₄⁺(S)=m₄⁺(πS); ∞: m₄⁺(S)=−m₄⁺(πS)).
- **Proved:** π preserves (κ,φ) on finite 4-sets; flips sign on ∞-sets.
- **Census p=5,7:** transport exact; π-covariance on finite+∞ |κ|=1 samples; m₄⁺=m₄⁻.
- **OPEN:** π-covariance general (Aut-orbit of S vs πS not automatic: type-class sizes at p=5 do not divide |PGL(2,25)|). Residual-(i) still OPEN. Predicates **False**.

### 15.260 — ∑μκ and ∑μφ identities; layer counts
- **Proved Max+-free:** ∑_S μ κ = n p²(p²−1)/8 (via E[R²]=binom(n,2)+2∑κμ).
- **Proved Max+-free:** ∑_S μ φ = n ∑μ = −n p²(p²−1)/12 (via ∑_r(Cy)_r⁴=np⁴ expansion).
- **Proved:** n₁,n₃ from ∑κ² and binom(n,4).
- **Census p=5,7:** identities exact; max|μ| on |κ|=1 is 3/65, 436/11452 <1/(2p).
- **OPEN:** L^∞ hinge (π-covariance / reflection / |μ|≤1/(2p)). Predicates **False**.

### 15.261 — Exact Max± configuration multiplicities
- **Proved Max+-free:** mult_±(α)=(N/16)(1±qf(α)/(2p)+pr(α)m₄±(S)) via hypercube Fourier + odd vanishing + pairwise π (any 4-set).
- **Proved:** mult≥0 ⇒ |m₄±|≤1−q_max/(2p); on |κ|=1 q_max≤8 ⇒ |m₄|≤1−4/p (too weak for 1/(2p) when p≥5).
- **Proved:** on |κ|=1 with orientation Q, mult_+(α)=mult_-(Qᵀα) ∀α ⇔ ν=0.
- **Census:** formula exact p=5,7; Q-match and ν=0 both hold on |κ|=1.
- **OPEN:** ν=0 / π-covariance / reflection / residual-(i). Predicates **False**.

### 15.262 — Change-of-var identity; matched Max± particulars; μ_part-extension 0
- **Proved Max+-free:** m₄⁺(S)=χ_D(S) m₄⁻(πS) for every 4-set (change of variables y'_i=D_i y_{π^{-1}i}: Max+≅Max−).
- **Proved:** on finite 4-sets ν(S)=−ν(πS) (anti-invariance under π; π²∈Aut).  Does **not** alone force ν=0 on |κ|=1.
- **Proved Max+-free:** Max± particulars in span{κ,φ,star} share (a,b)=μ_part coeffs and have opposite star coeffs z_∓=±2/(p(p²−5)); on |κ|=1 both reduce to μ_part (star=0), so ν=0 ⇔ δ₊=δ₋ on |κ|=1.
- **Proved:** ∑_j μ_part(T_j)=0 on every |κ₀|=1 ∞-triple (uses 15.256 ∑κ=2, ∑φ=p²−1).  Hence m₄⁺=m₄⁻ ⇔ ∑(μ−μ_part)=0 on extensions.
- **Certified** identity/anti-ν/μ_part-ext at p=5,7; particular algebra for p≤13.
- **OPEN:** ν=0 / envelope |μ|≤max(|μ_part|,|f4|) / reflection. Predicates still **False**.
- Evidence: `src/e1_gmin_m4_prop15262.py`, `evidence/e1_gmin_m4_prop15262.json`.


### 15.263 — Same Ext equation for m₄±; Ext[ν]=(p⁴−1)ν on |κ|=1
- **Proved Max+-free:** on |κ|=1, (p⁴−1)m₄⁻+2φ=Ext[m₄⁻] with the **same** Ext as Max+ (Cy=−py ⇒ ∏(Cy)=p⁴∏y; Per=1; size1+size2=−2φ under π=−C/p; size3=0).
- **Proved:** both m₄± solve the same inhomogeneous Ext equation ⇒ **(p⁴−1)ν = Ext[ν]** on |κ|=1.
- Combined with 15.262: ν|_{|κ|=1} is π-odd in ker(Ext−(p⁴−1)I) (coupled to |κ|=3).
- **Certified** same Ext / eigenrelation at p=5 (8 |κ|=1 samples, machine precision); ν=0 there.
- **OPEN:** prove π-odd kernel of Ext−(p⁴−1)I is trivial on |κ|=1 (spectral gap), **or** envelope / reflection / |μ|. Predicates still **False**.
- Evidence: `src/e1_gmin_m4_prop15263.py`, `evidence/e1_gmin_m4_prop15263.json`.


### 15.264 — Scheme Ext eq for μ_part; Ext[δ]=(p⁴−1)δ; ExtΠ=ΠExt
- **Proved Max+-free (scheme algebra):** on |κ|=1, Ext[μ_part]=(p⁴−1)μ_part+2φ via Per κ=p⁴κ−6φ, Per φ=(2n+1)φ and a,b of μ_part (coeff identity −6a+2nb−b(p⁴−1)=2 for all primes p≥5).
- **Proved:** Ext[δ]=(p⁴−1)δ on |κ|=1 for δ=μ−μ_part (subtract true Ext eq for μ).
- **Proved Paley:** ExtΠ=ΠExt (Per commutes with Π via UᵀCU=−C / permanent signs); Ext preserves π-parity.
- **Structure:** ν|_{|κ|=1} ∈ ker(Ext−(p⁴−1)I) ∩ π-odd (coupled to |κ|=3); census {0} at p=5,7.
- **Certified** Fraction identity p=5..200; numeric Ext[μ_part] at p=5.
- **OPEN:** π-odd kernel triviality / envelope / reflection. Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15264.py`, `evidence/e1_gmin_m4_prop15264.json`.


### 15.265 — T anticommutes with Π; even/odd master split
- **Proved Paley:** C_{πi,πj}=−D_i D_j C_{ij} from D Pᵀ C P D=−C.
- **Proved:** TΠ+ΠT=0 on functions of 4-sets (reindex Johnson neighbours; sign product −1).
- **Proved:** even/odd split Tμ=4pν, Tν=4pμ−4κ/p; recovers (16p²I−T²)μ=16κ.
- **Proved:** on |κ|=1, ν=0 ⇔ Tμ=0 ⇔ Tδ=0 (Tμ_part=0 there via star=0).
- **Certified** TΠ+ΠT=0 at p=5,7; master split exact at p=5; ν=0 on |κ|=1 samples.
- **OPEN:** Tδ=0 / Ext π-odd gap / envelope / reflection. Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15265.py`, `evidence/e1_gmin_m4_prop15265.json`.


### 15.266 — Type-averages of π-odd vanish; type-T ±4p evecs agree on star=0 (p=5)
- **Proved Paley:** π-odd functions have zero average on every finite (κ,φ)-class (π preserves (κ,φ)).  In particular type-averages of ν on finite |κ|=1 vanish.
- **Proved:** Tδ=0 on finite |κ|=1 ⇔ within-type π-odd variation of ν vanishes.  Point-transitivity reduces the rest to ∞-sets.
- **Proved (p=5 type-quotient):** T on the 12 (κ,φ,star) classes has simple eigenvalues ±4p, and the two evecs **agree** on star=0; type-constant odd E is 0 on |κ|=1.
- **Certified** 12×12 eig at p=5 (star0 restrictions agree).  At p=7 the 14×14 type-T has extreme eigs ±16, **not** ±28=±4p (E_{±4p} is purely within-type).  Census ν̃=0 at p=5,7.
- **OPEN at 15.266:** within-type odd variation (closed later by 15.268). Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15266.py`, `evidence/e1_gmin_m4_prop15266.json`.

### 15.267 — CR dichotomy; PSL fusion S∼πS; signed Aut; ε hinge
- **Proved Max+-free:** χ(λ)=T₂T₃, χ(1−λ)=T₁T₂; |κ|=1 ⇔ χ(λ)χ(1−λ)=−1; V4 dets put Stab_PGL(S) in the nonsquare coset, so S∼_PSL πS on finite |κ|=1.
- **Proved:** signed PSL Aut D_x=χ(cx+d); m₄ transform with a sign ε (convention corrected in 15.268).
- **OPEN at 15.267:** ε=+1 was stated for the g⁻¹ convention (wrong slot). Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15267.py`, `evidence/e1_gmin_m4_prop15267.json`.

### 15.268 — Pairing-pole square; ε(m_σ∘τ,S)=+1; ν=0 on all |κ|=1
- **Proved Max+-free:** pairing involution swapping a↔b, c↔d has pole r=(ab−cd)/(a+b−c−d) when Δ≠0, and ∏(s−r) is a square ⇒ ε(τ,S)=+1. Δ=0 ⇒ affine pairing ⇒ χ(det)=+1.
- **Proved:** on finite |κ|=1 a pairing τ has χ(det)=−1; g=m_σ∘τ ∈ PSL, g(S)=πS, ε(g,S)=+1.
- **Proved:** correct Aut convention (Uy)_i=D_i y_{g(i)} is orthogonal, so m₄⁺(S)=ε(g,S) m₄⁺(g(S)).
- **Proved:** m₄⁺(S)=m₄⁺(πS)=m₄⁻(S) on finite |κ|=1 ⇒ ν=0; PSL 3-transitive transport ⇒ ν=0 on every |κ|=1 four-set.  Closes the 15.266 within-type π-odd obstruction.
- **Certified** polynomial square identity; p=5 all 9900 finite |κ|=1 have ε=+1 and (Max+ cache) m₄⁺(S)=m₄⁺(πS); χ(−1)=1.
- **OPEN:** |μ|≤1/(2p) (ν=0 does not kill even δ: p=5 μ=f4≠μ_part; p=7 within-type even span). Residual (i)/E1/L still **OPEN**. Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15268.py`, `evidence/e1_gmin_m4_prop15268.json`.

### 15.269 — Additive Fourier support of Max+; Wick 3-point; \(\kappa_3\) criterion
- **Proved (Paley):** \(y\in\mathrm{Max}+\), \(y_\infty=+1\) \(\Rightarrow\hat z\) supported on \(\{0\}\cup\Omega\), \(\hat z(0)=p\), \(\sum\lvert\hat z\rvert^2=p^4\). Two-point: \(E[\hat z(\xi)\hat z(\eta)]=2p^2\) if \(\eta=-\xi\in\Omega\), else 0.
- **Proved:** on \(\infty\)-sets, \(\mu=\mathrm{Wick}+\kappa_3\) with \(\mathrm{Wick}=\kappa/p^2-2/p^3\); \(\lvert\mathrm{Wick}\rvert\le(p+2)/p^3\le1/(2p)\) for all primes \(p\ge5\). Room for \(\kappa_3\): \((p^2-2p-4)/(2p^3)\).
- **OPEN:** \(\lvert\kappa_3\rvert\le\) that room (or another listed residual-(i) hinge). `residual_i_closed_via_269()=False`. Predicates **False**.
- Evidence: `src/e1_gmin_m4_prop15269.py`, `evidence/e1_gmin_m4_prop15269.json`.

### G+ / inversion-T (scratch 2026-08-14; not a shipped close)
- **Proved Max+-free (scratch, see 08-14 session handoff):** no cuspidals in \(\mathcal W_{++}^0\); G-span(\(F\))=all; affine disks span \(F_{\mathrm{aff}}\); \(\hat z_{\mathrm{inv}}=0\) off \(\Omega\); separable formula \(\hat z(\omega t^2)=p1_{x=0}+p1_{y=0}+G_p[\chi_p(-2)J(2x^2)+\chi_p(-2\,\mathrm{ib})J(2\,\mathrm{ib}\,y^2)]\).
- **Certified:** affine+inv-\(T\) spans \(F\) at \(7\le p\le23\); \(N\neq0\) on every even K-character at every prime \(7\le p\le79\).
- **OPEN (blocks residual i):** prove \(N(\varphi)\neq0\) for all \(p\ge7\). Do not flip predicates on the census.
