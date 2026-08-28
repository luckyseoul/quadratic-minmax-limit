# Handoff: min-max ±1 quadratic form

**Date:** 2026-08-28 (15.680 additionally closes the `p=37,s=30` next all-finite endpoint; no general flag flipped)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit
**HEAD:** on `main`. Working brain is ALWAYS main.

**No leftover flag flipped.** Leftover 1/2/3 False. L OPEN. Aut-Schur /
Gsum / pairing False. `e1_closed_general` True only by the old incomplete
wiring. The full infinity-plus-point boundary is closed. Every four-point
boundary is additionally closed for `p>=11`, and infinity plus three finite
points is closed at `p=7`. The doubly saturated four-finite profiles at
`p=7` are also closed for both signs, and 15.655 closes all 23,520
unsaturated boundaries (518 orbits) per sign. Thus every `p=7` size-four
case is closed. Proposition 15.656 closes all 1,202 floor-surviving `p=5`
orbit/sign cases by complete-shell modular syndromes and sign transfer.
Every size-four boundary is now closed for odd `p>=5`; 15.657 closes size
six for `p>=11`, 15.658--15.659 close both `p=7` infinity-plus-five
signs, 15.660 closes every `p=5` size-six boundary, and 15.661 closes the
six-finite `p=7` branch. Every size-six boundary is now closed for odd
`p>=5`. Proposition 15.662 additionally closes both signs of the finite
`p=7` size-eight minimum-odd-secant/conic subbranch: all 32 floor-surviving
conic orbits are excluded. Proposition 15.663 next excludes the disjoint
83,770,008-boundary forced-floor nonconic stratum for both signs. Proposition
15.664 excludes another 23,563,806 boundaries per sign having exactly four
mean allocations. Proposition 15.666 excludes every one of the last
1,419,432 finite boundaries per sign. Thus finite `p=7` size eight is closed
for both signs. This is not all of the projective size-eight case: the
separate infinity-plus-seven profile remains open. p=13 orbits /
mesh k=6 are not a close. cpu44
stays hard-closed.

Proposition 15.669 now excludes every all-finite even boundary with
`6<=s<=3(p-1)/4` and every infinity-present boundary with odd finite
`5<=s<=p-4` for every odd `p>=17`. It also closes `p=11` infinity plus seven
and, at `p=13`, eight finite points and infinity plus seven or nine. The
first profiles outside these ranges survive only the floor-and-pair
relaxation and are not actual graphs. Proposition 15.670 additionally closes
every finite `p=11` size-eight boundary by a losslessly normalized exact
census. Residual (ii) therefore remains open; the first `p=11` profiles are
now infinity plus nine and finite size at least ten.

Propositions 15.671--15.672 close both signs of the collinear
infinity-plus-`(p-2)` geometry. Proposition 15.673 closes every realization
whose directional odd-fibre counts all lie in `{1,p-2}`, for both signs and
all primes `p>=17`. Same-type mean residues leave four arithmetic rows;
Segre's `p`-arc theorem kills pair-deficit equality, and the unique `p=17`
endpoint has exact inter-fibre norm `75>57`. Profiles with intermediate
odd-fibre counts are then closed by Proposition 15.674: exact type sums allow
an intermediate count only as the unique high direction of its type, which
leaves the same four arithmetic rows. The next infinity-present shell has
`p` finite points.

## 15.680 — `p=37,s=30` next all-finite endpoint closed

`src/e1_gmin_m4_prop15680.py` and
`evidence/NOTE_2026-08-28_p37_next_all_finite_endpoint.md`. The exact pair
ledger leaves only `u_0=2,3,4,5`; each forces a quotient-zero `b=0`
direction and a nonzero nonnegative integral quadratic lift. Proposition
15.642 excludes scaled means four, six, and eight.

At the remaining scaled mean ten, stabilizer averaging forces lift values
into `{0,1,2}`. The degree-two support floor bounds the density of value two
by `2/1295`, while the degree-four floor for `B(B-1)` is
`1938/441595`, strictly larger. Thus the lift is Boolean. A paired-cube
restriction proves every nonzero Boolean quadratic on
`J(p,(p+1)/2)` has density at least `(p-3)/(4p)`; at `p=37`, this is
`17/74`, contradicting the required `5/74`.

Thus `p=37,s=30` is impossible. The same boundary at
`p=17,19,23,29,31,41`, later all-finite sizes, strict infinity-plus-`p`,
residual (ii), R1, QVAR, Type I, and the limit remain open.

## 15.679 — next all-finite boundary closed from `p=43`

`src/e1_gmin_m4_prop15679.py` and
`evidence/NOTE_2026-08-28_next_all_finite_boundary_p43.md`. At the second
even all-finite size above `3(p-1)/4`, phase one has its unique residue
`u_1=m-1`. Exact phase-zero quotient arithmetic excludes `u_0=0,1`, every
interior residue at least eight, and the final four residues, leaving only
`2<=u_0<=7`. Each forces a quotient-zero `b=0` direction with scaled mean at
most 14. Proposition 15.642's degree-two slice-distance floor exceeds 14
for `p>=59`; exact pair ledgers plus lift floors close `p=43,47,53`.

Thus this entire next boundary is impossible for every prime `p>=43`.
Proposition 15.680 separately closes `p=37`; the smaller endpoints
`17,19,23,29,31,41`, later all-finite sizes, strict infinity-plus-`p`,
residual (ii), R1, QVAR, Type I, and the limit remain open.

## 15.678 — exceptional `p=17` first all-finite survivor closed

`src/e1_gmin_m4_prop15678.py` and
`evidence/NOTE_2026-08-28_p17_first_all_finite_survivor_exclusion.md`.
At `p=17,s=14`, phase one has only residue eight. The six-unit lift floor
excludes phase-zero `u_0=2`, the coefficient `l1` ledger excludes `u_0=3`,
and all residues at least four exceed the pair budget. For `u_0=0`, pair
slack is divisible by four and leaves exactly two profiles. Both attain pair
equality and are 14-arcs with secant distribution
`{7:6,6:8,1:1,0:3}`.

Adjoining any two of the three undetermined infinity points gives a 16-arc.
Sticker's exhaustive classification has one PGL class of 16-arcs in
`PG(2,17)`; conic-minus-two is a representative, so every 16-arc lies on a
conic. The third infinity point is off that conic, but after deleting four
conic points every off-conic point still lies on at least four secants of the
14-set. Contradiction. Together with 15.675/15.677, the first all-finite
survivor is closed for every prime `p>=17`. Later sizes remain open.

## 15.677 — first all-finite survivor closed from `p=19`

`src/e1_gmin_m4_prop15677.py` and
`evidence/NOTE_2026-08-28_first_all_finite_survivor_complete_from_p19.md`.
For the `p=1,7 mod 8` classes from `p=23`, the exact type-mean quotient
ledger leaves `u_0=2`, plus `u_0=3` in the first class. Since
`sum k_d=m-u_0<m`, a phase-zero direction has quotient zero. Its mean is
four or six, below the phase-zero floor `p+1` for every nonzero even fibre
count, so it has `b=0` and factors pointwise as `A_d=2B_d` with `B_d`
nonzero. Proposition 15.642 gives `4p E[B_d]>=8`, a contradiction.

Together with 15.675, this excludes the first even all-finite size above
`3(p-1)/4` for every prime `p>=19`. This proposition does not claim the
smaller `p=17` endpoint because it has an additional `u_0=0` residue row;
Proposition 15.678 now closes that endpoint separately. Proposition 15.679
closes the next size from `p=43`, and 15.680 closes its `p=37` endpoint;
six smaller endpoints, subsequent sizes, strict deficit in
infinity-plus-`p`, residual (ii), R1, QVAR, Type I, and the limit remain
open.

## 15.675 — first all-finite survivor closed in two modulo-eight classes

`src/e1_gmin_m4_prop15675.py` and
`evidence/NOTE_2026-08-28_first_all_finite_survivor_half_close.md`. At the
first even `s>3(p-1)/4`, exact same-type mean residues sharpen 15.669's
floor-only deficit minimization. Phase one uniquely has `m-1` directions at
`b=2` and one at `b=s`; phase zero uniquely minimizes at residue four using
quotient weights zero, one, and two at `b=0,2,s`. The total deficit minus
the pair budget is

```text
p mod 8:       1          3          5          7
gap:        -(p-1)/4   (p+1)/2   (p-1)/2   -(p-7)/4.
```

Therefore the first survivor is excluded for all primes `p>=19` in the
middle two classes. The outer two classes are an exact negative-margin route
boundary, not a graph witness. Local tests and an independent NUKA run agree.

## 15.674 — the entire near-line shell is impossible

`src/e1_gmin_m4_prop15674.py` and
`evidence/NOTE_2026-08-28_full_near_line_shell_complete.md`. Every odd-fibre
floor is at least `p-1`, while every intermediate floor is strictly above
`p+1`. Writing same-type means as

```text
a_d=2u+(p+1)k_d,   sum k_d=m-u,
```

excludes all interior residues. Residue zero gives only all-`p+1`
baselines; residue `p-1` gives `m-1` low baselines and one arbitrary
mean-`2p` exception. Two `b=1` baseline types exceed the pair-deficit budget;
two complementary baseline types determine at most two directions and are
collinear. The forced mixed pair has 15.673's offsets `0,1,2,1`, coefficient
congruences, support contradictions, and exact `p=17` norm close. Therefore
all odd-fibre profiles on infinity plus `p-2` finite points are excluded for
both signs and every prime `p>=17`. This is one full shell, not residual (ii).

## 15.673 — every endpoint-only near-line boundary is impossible

`src/e1_gmin_m4_prop15673.py` and
`evidence/NOTE_2026-08-28_endpoint_near_line_complete.md`. With
`q=(p-1)/2,m=q+1`, all same-type directional means share an even residue
modulo `p+1`. The endpoint floors and Proposition 15.642's four-unit lift
cost leave baseline counts `x,y`, at most one exception per type, and
`x+y<=7`. The four phase/residue rows give

```text
q|y,x-1;   q|x,y-1;   q|y,x+1;   q|x,y+1.
```

The first two rows violate the boundary support inequality at their unique
candidates, the third has none, and the fourth has none above `q=8`. At
`p=17,q=8`, the sole candidate has `I=5,E=64,P_d=7`; its complementary
baseline matrix has exact `l1` minimum 75 and transverse capacity 57.
The geometric equality case uses the classical odd-order `p`-arc conic
theorem, as detailed in the note. This is a uniform endpoint subbranch
close, not general residual (ii).

## 15.670 — every finite p=11 size-eight boundary is impossible

`src/e1_gmin_m4_prop15670.py` and
`evidence/NOTE_2026-08-28_p11_size_eight_boundary_exclusion.md`. Every finite
eight-set has an affine-similarity image containing field points `0,1`.
The exact pointed-set identity reduces exclusion from `C(121,8)` sets to all
`C(119,6)=3,470,108,187` normalized sets. The verifier audits all nonzero
scalars, translations, and projective directions, including the nonsquare
type swap and corresponding `c_H` phase transfer.

Complete V100/CUDA and RX 9070 XT/HIP replays test both signs on every
normalized set. Both full cost-pair histograms agree exactly, both signs have
zero survivors, and the minimum larger type cost is 76 against the exact
budget 72. Independent CPU combinations code matches every histogram entry
on a 100,000-set prefix. This closes finite `p=11` size eight only. The
permanent archive is
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-p11-size-eight-boundary/`;
its manifest SHA-256 is
`197616fc71490427822491afff651dfb32f7752627ebbcdc7b26ddc619de11de`.

## 15.669 — uniform non-Walsh boundary-range exclusion

`src/e1_gmin_m4_prop15669.py` and
`evidence/NOTE_2026-08-28_uniform_boundary_range_exclusion.md`. For every odd
`p>=17` and `5<=b<=p-5`, an explicit positive degree-two quadrature on the
parity-one contact nodes matches the first two hypergeometric moments. Hence
both parity-majorant values are exactly one and every scaled directional
floor is `2p` throughout the full middle range.

With no infinity, the two quadratic direction types have opposite phases.
An exact saving/deficit knapsack forces total pair deficit at least

```text
(p+1)s/4 - 1 + (p-1)(s-2)/2.
```

Subtracting the geometric budget `s(s-1)` gives one quarter of
`h_p(s)=s(3p+3-4s)-4p`, positive on the full interval
`6<=s<=3(p-1)/4` by endpoint concavity. With infinity, phase zero forces all
directions to have `b=1`; phase one forces at least `(p-1)/2` such directions
per type. The latter proof explicitly handles the complemented `b=p-4`
six-unit saving when `p=1 mod 4`. This excludes odd finite
`5<=s<=p-4`.

Exact rational floor/count-profile programs also give gaps 30 and 18 for
`p=11` infinity plus seven, gap 4 for `p=13` eight finite, and phase gaps
`(42,30)` and `(40,24)` for `p=13` infinity plus seven and nine. The first
relaxed survivors of 15.669 are eight finite / infinity plus nine at `p=11`,
but 15.670 subsequently closes the finite-eight branch. At `p=13` they are
ten finite / infinity plus eleven, followed by the first even
`s>3(p-1)/4` without infinity for `p>=17`, and `s=p-2` with infinity. A
survivor is only a directional count profile; incidence realizability and
the full residual graph constraints are still open.

## 15.668 — exact p=11 broad-channel theta and finite strong R1

`src/e1_gmin_m4_prop15668.py` and
`evidence/NOTE_2026-08-28_p11_broad_channel_theta.md`. A fourth marked
Legendre-convolution statistic refines the complete p=11 glue-profile census
into the kernel, low, and high eigenspaces of the square-circle operator.
Five-modulus CRT recovers all three raw shell-mass prefixes exactly through
exponent 120. Their common 32-dimensional affine modular spaces reach full
rank at exponent 92; all 28 held-out coefficients per channel match, and the
three forms reconstruct through exponent 800 with nonnegative masses and
exact aggregate conservation.

Eight channelwise QSopt_ex endpoints have independently checked exact
rational primal and dual certificates. Every certified interval still
contains a target mapping to `Phi<6`, so broad-channel conservation is now a
proved-insufficient relaxation. It does not refute R1. Independently, the
complete finite census gives

```text
||delta||^2 = 1382747375360/583792784981 < n/12 = 61/6,
strong margin = 27314875631681/3502756709886 > 0.
```

Thus strong R1 is rigorously true at `p=11`, but general R1 remains open. The
next R1 route must be finer character-resolved/PSL transport or a uniform
Mellin--Parseval inequality, not another broad square-circle aggregate. The
33-file permanent archive is
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1-broad-channel/`;
its manifest SHA-256 is
`d1ef69b9af7007c0d2f09a3a5ea8a014cde62d9ed6109175cf4a6496d06b3f07`.

## 15.666 — every finite p=7 size-eight boundary closed for both signs

`src/e1_gmin_m4_prop15666.py` and
`evidence/NOTE_2026-08-27_p7_size_eight_complete_exclusion.md`. The exact
post-15.664 partition contains 1,419,432 finite boundaries and 23,892,792
mean-allocation leaves per sign, in allocation-count strata 11, 16, 24, and
44. Conditioned omission scans leave 458,822 leaves modulo seven and
2,671,872 modulo three; intersecting the same leaf triples leaves 181,104.
Exact 22-row local, all-triple, and four-positive catalog joins reduce these
to 124,745, 78,126, and 62,892. A single-catalog filter rejects 3,777 of the
last leaves, and a complete meet-in-the-middle join rejects the other 59,115,
leaving zero.

The 22 base-seven digits are packed losslessly because `7^22 < 2^64`; the
join uses exact key equality, not a probabilistic hash. CPU and CUDA prefixes
agree at every stage, including 512 leaves in the final join, and three
representative cases independently fail the older full multi-characteristic
engine. A separate counter confirms zero hash-partition capacity rejections.
The final audit SHA-256 is
`428b9604e21738d9b063f0edee8a42b31d471ecd56800e4366af8ed1d7a49eaa`.
The nonsquare anti-isometry from 15.662--15.664 transfers the zero result to
the other product sign. Raw records are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-complete/`.

Scope is finite boundaries only. Infinity plus seven, residual (ii), Type I,
R1, global QVAR, and L remain open.

## 15.665 — conserved positive quartic mass on every R1 dual shell

`src/e1_gmin_m4_prop15665.py` and
`evidence/NOTE_2026-08-27_r1_conserved_quartic_shell_mass.md`. For a complete
dual shell `X_s`, put `b_x=Pi_Z(xx^T)` and
`R_s=sum_x b_x tensor b_x`. Then `R_s` is positive semidefinite and its
degree-four harmonic operator is exactly

```text
A_s = R_s - rho_s I,   rho_s = 2 N_s r_s^2/[d(d+2)].
```

The diagonal-map Gram inverse is closed, so one scalar trace-harmonic theta
series gives `tau_s=tr(R_s)`. On every multiplicity-free PSL constituent,

```text
q_(s,c) >= 0,   sum_c dim(c) q_(s,c) = tau_s,
q_(s,c) <= tau_s/dim(c).
```

The trace polynomial is also the coordinate-transitive orbit of one zonal
quartic, reducing its exact PARI computation from 122 coordinate fourth
powers to one. Exact `p=11` checks reproduce the four proved shell operators.
This supplies the nonlinear coupling absent from 15.641, but is not an R1
close. Propositions 15.667--15.668 subsequently reconstruct the exact scalar
and broad-channel series through exponent 800 and prove both resulting cones
insufficient. An all-prime character/transport inequality remains necessary.

## 15.664 — p=7 size-eight four-allocation stratum closed for both signs

`src/e1_gmin_m4_prop15664.py` and
`evidence/NOTE_2026-08-27_p7_size_eight_four_allocation_exclusion.md`.
The 24,983,238-case remainder after 15.663 splits by exact mean-allocation
count as `4:23,563,806`, `11:154,056`, `16:1,194,816`, `24:1,176`, and
`44:69,384` boundaries per sign. Every four-allocation boundary has one
type-floor sum 24 and one 32; its four leaves each raise one deficient-type
direction by eight. For each raised direction, a 112-dimensional exact
mod-seven dependency subspace vanishes on its full 35-column score block.
NUKA materializes 22 conditioned rows per direction, and the V100 directly
checks all 450,978,066 boundary ranks and 94,255,224 selected leaves. It
leaves 1,191 projected leaves and 1,176 full mod-seven survivors. Those are
exactly all `4*7*42=1,176` boundaries formed by one affine line of the
deficient direction type plus one point off the line. NUKA independently
rebuilds all candidates and both complete left kernels: each geometric
survivor has two mod-seven and 756 mod-three catalog rows, with empty
intersection. Thus zero exact catalog choices survive, and the nonsquare
anti-isometry transfers the exclusion to the other sign. The V100 result
SHA-256 is
`96cfe751a6c0f6bbcd86a1ef799c25847653f8db907414c7b85da576e02efe47`;
the independent audit SHA-256 is
`8129b608ec2e09967e10a7da7b38a8e20584450772ac7aab6c1c8a984a370e67`.
Exactly 1,419,432 finite size-eight floor survivors per sign remained at
this stage; Proposition 15.666 subsequently closes all of them. Raw records
are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-four-allocation/`.

## 15.663 — p=7 size-eight forced-floor stratum closed for both signs

`src/e1_gmin_m4_prop15663.py` and
`evidence/NOTE_2026-08-27_p7_size_eight_forced_floor_exclusion.md`. The two
complete floor files contain 2,016 ordered profiles and 83,770,008 nonconic
boundaries per sign with exact type-floor sums `(32,32)`. Exact type means
force every direction to its floor, so at most one complete catalog has more
than one row (36 rows). A direct-rank V100 pass exhausts all 450,978,066
boundaries: 526 pass its eight-row dependency projection and zero pass all
135 mod-seven dependencies. NUKA independently rebuilds the `282 x 1225`
score matrix, its rank-147 row space, all catalogs, and all 526 failures. The
nonsquare anti-isometry explicitly permutes the directions and fibres and
transfers the exclusion to `c_H=+1`. The GPU result SHA-256 is
`6143d4eb269861b3d380c53262b534e0a54a9645c9bbe7c29d9327200ae30535`;
the independent audit SHA-256 is
`7adaa5e76bf4f5e128c82ec219650b390c8c087d3aed2a44857f9da7939a9c53`.
Together with 15.662, this left exactly 24,983,238 finite size-eight floor
survivors per sign; Proposition 15.664 reduces them to 1,419,432 and 15.666
subsequently closes the remainder. Raw records are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-forced-floor/`.

## 15.662 — p=7 size-eight conic subbranch closed for both signs

`src/e1_gmin_m4_prop15662.py` and
`evidence/NOTE_2026-08-27_p7_size_eight_conic_exclusion.md`. Complete CUDA
floor censuses check all 450,978,066 finite eight-point boundaries per sign.
The 6,174 minimum-eight-odd-secant boundaries are affine conics; 4,851 fail
the floor and 1,323 survivors reduce to 32 stabilizer orbits. For
`c_H=-1`, 25 saturated orbits contribute 600 exact allocations, excluded as
355 initial CP-SAT, six long CP-SAT, and 239 catalog-join certificates. The
seven exceptional orbits contribute 1,260 allocations, excluded as 172
initial, 662 ordinary projected, and 426 high-direction-omission
certificates. Independent aggregate and component audits leave zero. The
nonsquare anti-isometry maps the full conic survivor set bijectively to
`c_H=+1`, closing both signs. The aggregate audit SHA-256 is
`85f927f41b3ffc9afe1a101584e95ed852709ca6e861b439d8da1715008640a9`.
The full floor census has 108,754,569 survivors per sign, so the next target
was the 108,753,246 nonconic remainder—not another conic orbit. Proposition
15.663 subsequently removes 83,770,008 of those, 15.664 removes another
23,563,806, and 15.666 closes the last 1,419,432 per sign. Raw conic records
are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-conic/`.

## 15.661 — p=7 six-finite and all size-six boundaries closed

`src/e1_gmin_m4_prop15661.py` and
`evidence/NOTE_2026-08-27_p7_six_finite_exclusion.md`. Exact floors leave
3,856,300 boundaries and 80,704 square-semilinear orbits. Simultaneous
mod-three/mod-seven catalog signatures reject all 160,745 elevation cases
in 80,519 ordinary orbits. Of 185 deep orbits, compact exact high-mean
models close 92; the 93 timeouts split into 930 exact mean allocations, of
which 810 close directly and the last 120 close by complete low-catalog
hash joins. NUKA independently reproduces the V100 survivor hash, ordered
orbit catalog, profile histogram, and ordinary exhaustion. The global audit
is true, and a nonsquare anti-isometry transfers the other product sign.
Raw records are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-six-finite/`.
With 15.657--15.660, every size-six boundary is closed for odd `p>=5`.
Larger boundaries and the top-level predicates remain open.

## 15.660 — every p=5 size-six boundary closed

`src/e1_gmin_m4_prop15660.py` and
`evidence/NOTE_2026-08-26_p5_size_six_global_exclusion.md`. Four exact
catalogs are rebuilt from definitions for both signs and infinity bits.
Signed symmetry and complete coarse SCIP batches leave six classes:
`0,881,2529,3032,4731,4939`. Independent layered audits reconstruct every
finite quotient and close all six with zero unresolved or feasible leaves.
The from-definitions global audit has SHA-256
`d6650a9f71043dce2902e157b56f988305470b911eaccb130522dd2f55b3bbd8`.
Raw class records are archived under
`/mnt/storage/e1work/maxplus_p13/p5_size6_circle_attack_2026-08-26/`; all
97 class-881 artifact hashes were rechecked after transfer. With
15.657--15.659, only six finite points at `p=7` then remained at boundary
size six; Proposition 15.661 subsequently closes that branch. Larger
boundaries and the top-level predicates remain open.

## 15.659 — negative p=7 infinity-plus-five branch closed modulo seven

`src/e1_gmin_m4_prop15659.py` and
`evidence/NOTE_2026-08-26_p7_size_six_negative_infinity_mod7.md`. Phase one,
the exact type budget, and congruence modulo eight force exactly one
scaled-mean-14 direction per quadratic type. Independent V100 and
Soulkiller floor sweeps agree on 83,496 survivors among all 1,906,884
boundaries. A serial NUKA enumeration and a GPU-seeded quotient agree on
1,750 square-semilinear orbits. Affine-span filtering rejects 2,205 of
2,230 elevation cases; exact comparison of 32,400 catalog pairs in the 25
remaining cases leaves zero survivors. NUKA and Soulkiller independently
reproduce both modular stages. Raw records and implementation are archived
at `/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p7-size-six-negative-infinity/`.
Together with 15.658 this closes both infinity-plus-five signs at `p=7`.
Proposition 15.660 subsequently closes every `p=5` size-six branch, and
15.661 closes six finite points at `p=7`; larger boundaries remain open.

## 15.658 — positive p=7 infinity-plus-five branch closed modulo seven

`src/e1_gmin_m4_prop15658.py` and
`evidence/NOTE_2026-08-26_p7_size_six_positive_infinity_mod7.md`. Phase
zero and the exact type budget force scaled mean eight in every direction.
The unique `J(7,4)` slacks for `b=1,3,5` determine all 280 affine score
right sides from the finite five-set. The common `282 x 1225` edge system
has rank 147 over `F_7`, hence 135 dependency checks. A raw V100 integer
kernel rejected all `C(49,5)=1,906,884` boundaries in 2.83 seconds; an
independent NUKA/NumPy sweep reproduced the same mask histogram and zero
survivors in 4.47 seconds. Raw records and implementation are archived at
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p7-size-six-positive-infinity/`.
Proposition 15.659 subsequently closes the negative infinity sign. Six
finite points at `p=7` remain open; Proposition 15.660 subsequently closes
every `p=5` size-six branch.

## 15.656 — every p=5 four-point boundary closed by full-shell syndromes

`src/e1_gmin_m4_prop15656.py` and
`evidence/NOTE_2026-08-26_p5_four_point_full_shell.md`. Each complete `p=5`
eigenshell has 130 antipodal score rows, normalized edge-column sum 26, and
total slack 78. Edge count, the distinguished edge, and all bad-edge counts
give a `132 x 325` matrix of rank 67 over `F_5`, hence 65 exact left-null
dependencies. Bounded lift systems exclude 712 direct orbit cases modulo
five. The only timeout, no-infinity negative orbit 164 with representative
`[2,3,12,13]`, is independently infeasible modulo seven. A signed
nonsquare Paley anti-isometry transfers the 489 no-infinity negative orbits
to the positive sign. A fresh structural audit reconstructs all four orbit
sources, both shell ranks, every parity/lift mass, and the sign-transfer
bijection. All 1,202 floor-surviving orbit/sign cases and 26,450
boundary/sign cases are closed, with zero unknown. The archive is
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p5-four-point/`
and has SHA256
`d5db5e82389ebb0bfcb23e80da5e2322b1d65e74aa8f3804d25275793b7380da`.

## 15.655 — unsaturated p=7 four-finite profiles closed modulo seven

`src/e1_gmin_m4_prop15655.py` and
`evidence/NOTE_2026-08-26_p7_four_finite_unsaturated_mod7.md`.  For a fixed
boundary/elevation case, the 280 exact affine bad-edge counts, total edge
count, and distinguished edge give a `282 x 1225` integer system. Its rank
over `F_7` is 147, so every catalog right side must satisfy 135 left-null
dependencies. Exact complementary-syndrome joins cover all one- and
two-catalog patterns: 1,716,742,440 tuples across 2,408 cases, with zero
compatible tuple. Those cases cover all 518 unsaturated orbits and 23,520
boundaries for `c_H=-1`. An independent implementation rebuilds the matrix
from sign products, dependencies by incremental span reduction, catalogs
from interpolated target coefficients, and all coverage keys; it again finds
zero survivors. The 15.654 nonsquare anti-isometry transfers the result to
`c_H=+1`. The raw certificate, audit, and orbit source are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p7-four-point/`.

## 15.654 — saturated p=7 four-finite profiles closed

`src/e1_gmin_m4_prop15654.py` and
`evidence/NOTE_2026-08-26_p7_four_finite_saturated_exclusion.md`. Exact
Johnson evaluation gives one phase-zero and 36 phase-one saturated `b=4`
slacks. Of 82,320 boundary-only survivors per sign, exactly 58,800 have
both type costs equal to 32; they form 1,225 square-semilinear orbits. All
1,225 exact fixed-boundary coefficient models are infeasible. Multiplication
by a nonsquare, with switching at infinity, is an exact Paley anti-isometry
that fixes the distinguished edge, swaps eigenshells, preserves normalized
scores, and flips the product sign, so the sweep covers both signs. An
independent re-enumeration/audit passes every coverage and model-scope
check. The remaining 23,520 unsaturated boundaries (518 orbits) per sign
are subsequently closed by 15.655.

## 15.653 — p=7 infinity plus three finite points closed

`src/e1_gmin_m4_prop15653.py` and
`evidence/NOTE_2026-08-26_p7_infinity_three_exclusion.md`. The negative
product sign is already impossible by 15.652. In the positive sign, every
direction has scaled mean eight. The three-odd-fibre parity pattern has 630
possible mass-four corrections inside the 35-point Johnson slice; exact
rank-21 evaluation leaves uniquely `A=(|X cap B|-2)^2`. Sparse coefficient
equations then reduce all 18,424 finite triples to 416 square-semilinear
orbits. All 416 fixed-boundary models are infeasible; a fresh orbit
reclassification and independent audit report zero missing, malformed,
unknown, or feasible case. Four finite points at `p=7` are subsequently
closed by 15.654--15.655; all `p=5` size-four cases are subsequently closed
by 15.656.

## 15.657 — every six-point boundary closed for p>=11

`src/e1_gmin_m4_prop15657.py` and
`evidence/NOTE_2026-08-26_size_six_boundary_exclusion.md`. Exact positive
quadrature extends the parity floors through `b=6`. For `s` finite boundary
points, unique pair directions give the global deficit inequality
`sum_d(s-b_d)<=s(s-1)`. Combining its budgets 30 and 20 with the affine
slack floors excludes infinity plus five finite points for every `p>=11`
and six finite points for every `p>=13`; a separate exact type-profile
argument closes `p=11`. Proposition 15.660 subsequently closes `p=5`, and
15.661 closes the six-finite `p=7` case. Boundaries of size at least eight
remain open.

## 15.652 — every four-point boundary closed for p>=11

`src/e1_gmin_m4_prop15652.py` and
`evidence/NOTE_2026-08-26_four_point_boundary_exclusion.md`. Exact positive
quadrature proves the small-fibre floor table `(0,2p)`, `(p+1,p-1)`, and
`(2p-6,2p)`. Four finite points have only six pair-collision directions;
infinity plus three finite points has only three. Those counts contradict
the separate direction-type budgets for every odd prime `p>=11`. Combined
with the empty and two-point exclusions, the first open boundary size there
was six; 15.657 subsequently closes size six for `p>=11`. At `p=7`, 15.653 subsequently closes the
infinity-present positive-product case as well. The unsaturated four-finite
profiles are subsequently closed by 15.655. The exceptional `p=5`
size-four shapes are subsequently closed by 15.656.

## 15.651 — positive two-point remainder closed for p=5,7,11,13

`src/e1_gmin_m4_prop15651.py` and
`evidence/NOTE_2026-08-26_complete_positive_two_point.md`. Exact additive
coefficient equations and fibrewise `l1` profiles close all seven `p=5`
arithmetic cases and the nonzero-`k0` finite cases. At `k0=0`, one empty
direction forces a uniform finite-edge type, excluding `p=11,13` by type
capacity. At `p=7`, the remaining type split has 2250 five-stars and 56
square-semilinear orbits per populated type; all 112 exact fixed-star models
are infeasible. The alternative all-eight-`kd=1` profile has three exhaustive
normalizations, also all infeasible. With 15.643 this closes `c_H=+1` for
every odd `p>=5`; with 15.650, both product signs of `D={infinity,v}` are
closed. Other boundary shapes, residual (ii), R1, and L remain open.

## 15.649--15.650 — negative two-point remainder completely closed

Proposition 15.649 classifies all 1764 balanced `p=7` exceptional lifts and
certifies all 6076 star-orbit representatives infeasible. Proposition 15.650
then reduces `p=5` to 24 arithmetic profiles and 33 placement orbits, all
infeasible. Together with 15.647--15.648, `c_H=-1`,
`D={infinity,v}` is closed for every odd `p>=5`.

## 15.648 — negative two-point remainder reduced to p=5 and one p=7 profile

`src/e1_gmin_m4_prop15648.py` and
`evidence/NOTE_2026-08-25_finite_negative_two_point_exclusions.md`. At
`p=13`, the prescribed zero-baseline inter-fibre matrix has exact minimum
`l1=48>44` transverse edges, closing both orientations algebraically. At
`p=11`, both count orientations are CP-SAT infeasible on all three proven
type-preserving exceptional-pair orbits. At `p=7`, all-pair sweeps exclude
the four unbalanced profiles `(0,3),(0,6),(3,0),(6,0)`. Raw result hashes
and the storage archive are recorded in the proposition evidence. Only
`p=5` and balanced `p=7 (x,y)=(3,3)` remain in the negative two-point branch.

## 15.647 — both two-point product branches closed for p>=17

`src/e1_gmin_m4_prop15647.py` and
`evidence/NOTE_2026-08-25_all_prime_negative_product_exclusion.md`. In the
negative branch, the exact signed mean gives
`a_d-a_e=(p+1)(P_d-P_e)` for same-type directions. The all-prime sparsity
bound leaves a baseline in each type for every `p>=7`, so every nonzero lift
excess is a multiple of `p+1`. Its type total is exactly `p+1`; hence there
is exactly one exception per type, with no asymptotic theorem. Baseline
divisibility forces both baseline parallel counts to be multiples of
`(p-1)/2`. For `p>=17`, edge counts make both zero, leaving `4p-1` infinity
edges and only two finite edges, incompatible with boundary `{infinity,v}`.
Thus `c_H=-1` is closed for all odd `p>=17`; 15.643 already closes `c_H=+1`
there. The remaining two-point primes are `5,7,11,13`.

## 15.646 — negative-product infinity-point boundary excluded asymptotically

`src/e1_gmin_m4_prop15646.py` and
`evidence/NOTE_2026-08-25_negative_product_boundary_exclusion.md`. In every
baseline direction of the 15.644 normal form,
`K_st=-eps_d(a_s+a_t)` with `sum_s a_s=0`, so the signed sum of all
transverse finite edges is zero. The exceptional splits `(U,V)=(3,1)` and
`(1,3)` make the global finite-edge signed sum `+2` and `-2`, respectively.
Removing the two parallel edges in a negative baseline in the first case
leaves `+4`; removing the two positive parallel edges in the second leaves
`-4`. Both contradict zero. Hence `D={infinity,v}`, `c_H=-1` is empty for
all sufficiently large odd primes. The threshold is qualitative, and other
boundary profiles and all top-level gates remain open.

## 15.645 — baseline fibres are ideal or one-transfer

`src/e1_gmin_m4_prop15645.py` and
`evidence/NOTE_2026-08-25_baseline_fibre_profiles.md`. In every baseline
direction of 15.644, put `w_s=n_s+1_{s=j}` and `a_s=w_s-2`. Then
`sum a_s=0` and `sum_{s<t}|a_s+a_t|<=2p`. If `A` is the total positive
mass of the integral vector `a`, an exact sign decomposition gives a lower
bound exceeding `2p` whenever `A>=2` and `p>=7`. Hence `a=0` or `a` has one
`+1`, one `-1`, and all other entries zero. This is the precise two-line
fibre fingerprint. Its simultaneous geometric classification remains open
but is no longer needed for this branch because 15.646 bypasses it.

## 15.644 — asymptotic normal form of the negative two-point branch

`src/e1_gmin_m4_prop15644.py` and
`evidence/NOTE_2026-08-25_negative_product_asymptotic_normal_form.md`.
The near-optimal degree-two slice-distance theorem implies that, for all
sufficiently large odd primes, the `D={infinity,v}`, `c_H=-1` branch has
exactly one nonbaseline direction of each quadratic type. Combining the
forced exceptional mean with baseline divisibility and signed `l1` bounds
leaves the unique arithmetic profile

```
I=2p-1 infinity edges; E=2p+2 finite edges;
P_d=2 in every baseline direction;
exceptional parallel counts (positive,negative)=(1,3) or (3,1).
```

This proposition alone is not an exclusion. Proposition 15.646 now excludes
its complete normal form directly by signed transverse-sum bookkeeping.

## 15.643 — positive-product infinity-point boundary excluded for p>=17

`src/e1_gmin_m4_prop15643.py` and
`evidence/NOTE_2026-08-25_positive_product_boundary_exclusion.md`. In the
`D={infinity,v}`, `c_H=+1` branch, every directional slack is pointwise
baseline by 15.642. The signed inter-fibre edge matrix is therefore additive.
If `I` counts infinity edges and `P_d` finite edges parallel to direction
`d`, exact divisibility and the global edge count give

```
q=(p-1)/2,  I=5+q*k0,  P_d=q*kd,  sum_d kd=8-k0.
```

The directional `l1` budget gives
`abs(k0+p*kd-8)<=8-k0-kd`; a populated direction would satisfy
`2*k0+(p+1)*kd<=16`, impossible for `p>=17`. The remaining all-infinity
star has boundary size `4p+2`, not two. Thus this branch is empty for every
odd `p>=17`. At the time this left four small positive-product cases; they
are now closed by 15.651. Other boundary profiles remain open.

## 15.642 — nonzero quadratic mass; two-vertex boundary rigidity

`src/e1_gmin_m4_prop15642.py` and
`evidence/NOTE_2026-08-25_quadratic_lift_mass.md`. A stabilizer moment
certificate gives an exact elementary mass floor for any nonzero
nonnegative integer-valued quadratic on the middle slice. The exact
all-parameter degree-two distance lemma of Amireddy--Behera--Srinivasan--
Sudan additionally gives

```
Pr(B != 0) >= C(p-4,m-2)/C(p,m) = (p^2-1)/(16p(p-2)).
```

For residual boundary `D={infinity,v}`, the positive edge-product branch
has zero budget beyond parity, hence `A_d=x_{s_d(v)}` pointwise in every
direction. In the negative edge-product branch, each quadratic type has at
most three directions with a nonzero lift, uniformly for `p>=5` (two at
`p=7`). This is a major reduction from a growing exceptional set, but not a
boundary exclusion. The next live attack is the repeated additive-rank-two
constraint on the signed inter-fibre edge matrices in all baseline
directions. Exact p=5 CP-SAT probes remained `UNKNOWN` and are not evidence.

## 15.641 — current modular shell/cusp data do not determine R1

`src/e1_gmin_m4_prop15641.py` and
`evidence/NOTE_2026-08-25_p11_modular_independence.md`. At `p=11`, the
relevant Kohnen space has dimension 66. The justified infinity coefficients
`0..19` and geometric cusp gaps have rank 29; adding the complete second
dual shell at coefficient 20 gives rank 30 and leaves a 36-dimensional
kernel. The second-shell and target rows have joint rank two.

An exact 21-coordinate rational vector gives the stronger certificate

```
A_known w = 0,   c_second w = 0,   c_target w = 1.
```

Its binary payload is backed up under
`/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25/` with SHA-256
`5bdf184e653079c361f6ee1a2178dd3f4e9b051d9da6625cc3a4910a93b441e7`.
This closes coefficient determination from the current linear modular data,
not R1: additional shells/cusp constraints or nonlinear positivity specific
to lattice theta series remain possible.

## 15.640 — exact quartic saddle at scaled norm 3p-6

`src/e1_gmin_m4_prop15640.py` and
`evidence/NOTE_2026-08-25_scaled_norm_3p_minus_6_harmonic_saddle.md`.
The square-circle complement words through any point satisfy the new exact
frame identity

```
sum_{S contains i} w_S w_S^T = p^2 (P - 2(Pe_i)(Pe_i)^T).
```

After moving `i` to infinity, this is just `(p+1)/2` mutually orthogonal
parallel classes, each with Gram matrix `p(pI-J)` and frame eigenvalue
`p^2`; their dimensions sum to `d-1`.

For admissible `W`, `F=||W||_F^2`, the two complete families from 15.639
have quartic sums

```
negative triples:  2(p-3)(p+1) F
point--circles:     8(p-2) F + 2(p-5)/p^3 sum_S (w_S^T W w_S)^2.
```

Combining these with the radial terms in `H_W` and the 15.634 circle tensor
spectrum gives three complete shell eigenvalues:

```
-(p^4+2p^3-69p^2+136p+26)/(p^2+5),
 (p^4-14p^3+89p^2-196p+24)/(p^2+5),
 (p^4-10p^3+69p^2-176p-76)/(p^2+5).
```

They have signs `(-,+,+)` for every `p>=11`; positive-coefficient
expansions at `p=11` prove this uniformly. At `p=11` they are
`-582/7, 258/7, 426/7` with multiplicities `1220,305,244`. The odd norm
phase and evaluation at `x/2` multiply the spectrum by `-1/16`, giving
shadow signs `(+,-,-)`. Thus this shell reinforces the circle kernel but
cancels both circle-image channels. It is exact multi-scale R1 information,
not a bound on intervening or later shells, so R1 and L remain open.

## 15.639 — complete shell at scaled norm 3p-6 for every p>=11

`src/e1_gmin_m4_prop15639.py` and
`evidence/NOTE_2026-08-25_first_nonminimal_odd_scaled_shell.md`. Every
vector at the first possible nonminimal odd scaled norm `s=3p-6` has odd scaled
coordinates `r=2px`; the identity `sum r_i^2=2ps<9(p^2+1)` forces a unit
coordinate. Signed Paley transport makes its common profile sum `t=1`.
Equality in the 15.635 profile bound leaves only one active profile or
`R-2` active profiles.

The one-profile branch is the square-circle equality case; its only other
multiplicity pattern is exactly the exceptional profile killed in 15.636.
In the dense branch, the degree-two and degree-three moment differences
satisfy

```
Q2=-2AB,   Q3=-3AB(2mu+A+B).
```

Their two inactive roots force `Q2|Q3`. Adding the minimum vector encoded by
the resulting linear quotient turns the vector into a common-sum-two vector
on the complete second, third, or empty `2(p+3)` shell. Props. 15.633,
15.636, and 15.638 leave exactly two families:

```
Pz,  z supported on a signed triangle with all three conference edges -1;
+/-(Pe_i+w_S/p),  (w_S)_i=-1.
```

They are disjoint by their scaled-coordinate signatures (three versus one
coordinates of magnitude `p-2`). The complete signed count is

```
p^2 (p-1) (p+7) (p^2+1) / 6.
```

At `p=11`, exact NUKA `qfminim` through bound 28 returned `473,970` vectors
and maximum norm 27. Subtracting the complete count `31,110` through the
third shell leaves `442,860`, exactly the formula. This is the fourth shell
only at `p=11,13`; for `p>=17`, unclassified even candidates can intervene.
Its harmonic operator is supplied by 15.640. The intervening and later theta
tail remains open, so R1 and L are still open.

## 15.638 — the first post-third even candidate shell is empty

`src/e1_gmin_m4_prop15638.py` and
`evidence/NOTE_2026-08-25_empty_post_third_even_candidate_shell.md`. Proposition
15.637 removed `t=0` at scaled norm `2(p+3)`. Balancing leaves only
`|t|=2,p-1,p+1`, with total profile-energy excess respectively `4,4,2`.

For `t=2`, unsigned-pair cubic and quartic moment recurrences eliminate
all one-exception patterns by factored defects or Newton identities. With
two exceptions, the cubic forces both to be doubled points. Then
`D=2q2-q1^2` is square on the selected projective half and has two roots
there. The character sum of `ND` would have magnitude `p-3`, while the
smooth genus-one curve `Y^2=ND` gives at most `2sqrt(p)` by Hasse, impossible
for `p>=11`. For `t=p-1` and `p+1`, subtracting the all-one profile reduces
to common sum one; delta moment identities and Newton exclude every defect
pattern. Hence no dual vector has scaled norm `2(p+3)` for any `p>=11`.

This is a genuine uniform even-shell gap, but not the fourth norm or R1:
the odd-phase floor remains separate (and is 27 at `p=11`), while the next
nonempty shell and complete later harmonic theta tail remain unknown.
Exhaustive local-pattern tests pass, and an independent
finite-field audit finds no forbidden half-conic through `p=43`.

## 15.637 — first post-third zero-common-sum channel excluded

`src/e1_gmin_m4_prop15637.py` and
`evidence/NOTE_2026-08-25_one_profile_next_energy_gap.md`. At the first
possible even energy after the third shell, `E=p+3`, the profile MDS mass
bound initially leaves only active counts `h=1,R-1,R`. For `h=1`, integer
energy permits exactly one doubled entry with full field support or two
doubled entries with three omitted roots. The first case collapses to
`T=(X-alpha)^((p+1)/2)` and makes the positive/negative root supports
intersect. In the second case, the reverse formal square root
`S=sqrt(N/D)` has a long zero coefficient block. The differential equation

```
2 N D S' = (N' D - N D') S
```

descends through that block and forces `S=1`, hence `N=D`; repeated roots
would again be omitted roots. Therefore `h=1` is uniformly impossible for
`p>=11`.

The dense branches fall to low-degree moment recurrences. A signed pair's
first three moments obey

```
4 q1 q3 - 3 q2^2 - q1^4 = 0.
```

For `h=R`, its `k>=5` ordinary-pair zeros force this binary quartic to be
an identity, but the unique energy-four profile has nonzero factored
defect. For `h=R-1`, divisibility by the inactive direction reduces the
same relation to a cubic. The remaining six-unit energy-six case also uses
the degree-four pair recurrence; four matched power sums and Newton
identities would make two disjoint root multisets equal. Thus **no
zero-common-sum profile** inhabits energy `p+3`. The nonzero common sums
`2,p-1,p+1` were the boundary of 15.637 and are closed by 15.638 above.
Exact finite scouts independently agree through `p=127` for the
one-profile branch.

## 15.636 — complete third dual shell for every p>=11

`src/e1_gmin_m4_prop15636.py` and
`evidence/NOTE_2026-08-25_complete_third_dual_shell.md`. Proposition
15.635 left one possible zero-profile equality case. Its two degree-`m`
root polynomials differ by a constant and, after normalizing the two omitted
field points to `0,1`, would make

```
x^(p-1) + a(x+...+x^(p-2)) + b
```

a square with `a(a-1)!=0`. Reversing the square root creates a coefficient
gap in `(1+(a-1)y)^(m+1)(1-y)^m`. Comparing high Hasse derivatives at its
two roots forces `a-1=1/3` and `(a-1)^2=1/5`, hence `p|4`, impossible.
Therefore the complete third shell for every `p>=11` is exactly
`+/-P(e_i+C_ij e_j)`, with signed count `p^2(p^2+1)` and negative scalar
harmonic operator from 15.635. The first three dual shells are now complete;
R1 still needs a uniform bound on the fourth-and-later tail.

## 15.635 — third dual norm; p=11 audit

`src/e1_gmin_m4_prop15635.py` and
`evidence/NOTE_2026-08-25_third_dual_norm.md`.  With scaled dual norm
`s=2p||u||^2`, the next norm after `p` and `2(p-1)` is `2(p+1)` for every
`p>=11`.  The profile/MDS proof also gives the new odd-phase gap
`s>=3p-6` outside the minimum shell.  The signed point-pair family
`+/-P(e_i+C_ij e_j)` attains the third norm, has size `p^2(p^2+1)`, and
its complete harmonic contribution is the negative scalar

```
-(p^2+4p-3)/(4(p^2+5)) ||W||_F^2.
```

At `p=11`, exact saturated-dual `qfminim` through scaled norm 24 returned
31,110 vectors.  Subtracting the proved first/second counts `244+16,104`
leaves `14,762`, exactly the point-pair count, so the complete p=11 third
shell is classified.  For `p>11`, the norm and pair orbit are proved but
shell completeness was not claimed in 15.635; it is now supplied by 15.636.
The bare count is OEIS A071253/A069187; the lattice interpretation is the
content. R1 and the limit remain open.

## 15.634 — square-circle spectrum; complete second R1 shell is negative

`src/e1_gmin_m4_prop15634.py` and
`evidence/NOTE_2026-08-25_square_circle_operator.md`.  For the square
`F_p`-sublines, let `M` be point--circle incidence and let `A` join two
circles meeting in two points.  The exact all-prime operator identity is

```
A^2+pA = (p^2-1) M M^T/8 + (p-1)^2(p+1) J/8.
```

It gives `Spec(A)={k_2, (p-1)^2/4, -p, 0}` with multiplicities
`1,n-1,n(p-1)/4,n(p-3)/4`.  Projecting the signed-complement tensors
`w_S w_S^T` to the admissible harmonic space gives Gram spectrum
`0^n`, `[p^3(p-1)]^[n(p-1)/4]`, and
`[p^3(p+1)]^[n(p-3)/4]`.  Therefore the complete signed norm-`(p-1)/p`
harmonic shadow shell has the three explicit eigenvalues

```
-(p+2)(p^2-4p+1)/(4p(p^2+5)),
-(p^3-3p^2-19p+9)/(8p(p^2+5)),
-(p^3-5p^2-19p-1)/(8p(p^2+5)).
```

All are negative for every `p>=11`.  This is an exact cancellation channel,
so first-shell-only positivity is dead.  The first-two-shell truncation can
be positive only in an `O(1/log p)` Gaussian window; no tail estimate is yet
proved.  Exact construction audits pass at `p=3,5,7,11,13`.  R1, global
QVAR, and the final limit remain open.

## 15.633 — complete second dual shell

`src/e1_gmin_m4_prop15633.py` and
`evidence/NOTE_2026-08-25_dual_second_shell.md`.  For every `p>=5`, every
dual vector of norm `(p-1)/p` is, up to sign, exactly one of

```
P(e_i-C_ij e_j),                    i<j,
w_S/p,                              S a square F_p-subline.
```

The union is disjoint and the signed shell count is
`p(p+1)(p^2+1)`; at `p=3` the two descriptions overlap and the count is
`30`.  The proof uses integral circle profiles, MDS/Newton equality cases,
and a half-conic rigidity lemma.  Exact `qfminim` and CUDA numerator audits
match at `p=3,5,7,11`; OEIS and literature searches found no duplicate of
the formula-plus-lattice classification.  This classifies one complete
shell, not the later dual tail.

## 15.632 — type-split affine slack budget; Eulerian boundary empty

`src/e1_gmin_m4_prop15632.py` and
`evidence/NOTE_2026-08-25_affine_slack_parity_budget.md`.  For an odd
separator `H`, every affine direction gives a nonnegative integer quadratic
slack `A_d` on `J(p,(p+1)/2)`.  With `a_d=2p E A_d`, direct edge moments give

```
sum_{eps_d=tau} a_d = (p+1)(|H|-3p)/2,  tau=+-1,
```

and every `a_d` is even.  If `D` is the odd-degree boundary, its odd fibres
give the exact parity `A_d(x)=sum_{B_d}x+eta_d (mod 2)`.  Averaging under
`Sym(B_d)xSym(B_d^c)` reduces the sharp magnitude floor to an exact
three-variable LP `M(p,b,eta)`, hence
`a_d>=2 ceil(p M(p,b_d,eta_d))`, separately budgeted in each direction type.
At residual `|H|=4p+1`, `D=empty` forces one type to cost `p(p+1)` against
budget `(p+1)^2/2`, a gap `(p^2-1)/2`; the Eulerian branch is therefore
impossible for every odd prime.

Do not soft-close: the corrected p=5 affine model has an integral witness
with `a=(12,4,0,6,10,4)`, all pointwise slacks nonnegative, and boundary
infinity plus an affine line.  Nonempty profiles and the full shell remain
open.  `evidence/NOTE_2026-08-25_pbss_cross_audit.md` records a read-only
cross-audit: no PBSS theorem transfers directly; its useful R1 idea is an
adversarial-cancellation test followed by a multi-Gaussian theta-window LP.

## 15.628 — eligible GQR circles + affine completions close W2 and Walsh

`scripts/w2_affine_circle_close.py` and
`evidence/NOTE_2026-08-24_w2_gqr_circle_route.md`.  Nonsquare circles
meeting `{0,infinity}` evenly span `H_0 cap ker(e_0+e_infinity)` for every
odd prime: tangent pencils reduce the quotient to a connected bipartite
Cayley graph, and Katz's `t=-2` Soto--Andrade bound proves connectivity.
For every `T subset F_p` with `|T|=(p+1)/2`, the affine sign vector
`h_T(infinity)=1`, `h_T(a+b omega)=+1 iff b in T` satisfies `C h_T=p h_T`
(square-line character sum `-1`).  Nonsquare dilation turns every
`0 in T` into a completion of the standard nonsquare circle.  For an
outside pair with transverse coordinates `b,d`, prescribe
`s_T(b)s_T(d)=chi(u-v)`; this is always extendible for `p>=5`.  Hence every
eligible circle is an actual U-difference and
`dir(U)=H_0 cap ker(e_0+e_infinity)`.  The `p=3` slice is exact rank 4 and
W2 is vacuous.  Therefore W1, W2, and Walsh 15.406 E hold for every odd
prime.  Exact p=19 witness for the old MILP pair `[2,340]` uses
`T={0,1,2,3,4,5,6,7,8,18}`; the solver UNKNOWN is superseded.  This closes
the Walsh slice only; 5+-level / even-`k>4p` and other E1 leftovers remain
open.

## 15.627 — octic linear box empty; split-involution class W2 at p=31

`src/e1_gmin_m4_prop15627.py`. Fable rank 1 as a piecewise linear
stay on (2/p)_8 (and (octic, a_G mod 8)) is empty in the 15.626
box. Rank 2: every switched split involution is Max−; class size
p(p+1)/2. p=17: 17 W2 / 49 in U / 153. p=31: 76 W2 / 146 in U /
496, including x/(x−1), though t=-2 fails. W2 p-law = counting
identity, not another named t. k=6 mesh skipped. leftover flags
untouched.

## 15.626 — Fable a,b-windows empty; W2 t=-2 through p=17

`src/e1_gmin_m4_prop15626.py`. Residual W1: the bounded box
`ua+vb+wi+k` (|u,v,w|≤4, |k|≤8) is empty as an upper-half stay
on 61 primes p=a²+64c² ≤19441 (not a reduction of all of Z⁴).
Named ±a,±b,a/2,ib MIXED
(ε(-a)=1 at 601, 0 at 1201). Split involution t=-2 is W2 at
p=17 (and 7,23,41,47) but not a p-law: p=31 is in U with
gcd(c,g)≠1. t=i is W2 at p=17, dead at p=41. W1 residual /
W2 p-law / Walsh / leftover 2 OPEN.

## 15.625 — Fable Walsh consult: d=−(p−1)/8 iff (2/p)_4=−1

`src/e1_gmin_m4_prop15625.py`. Deep_review PASS-WITH-NOTE. Quarter
d=−(p−1)/4 is identically ε=0 on p≡1 (mod 8). Eighth-interval:
ε=1 iff (2/p)_4=−1 (Barrucand–Cohn). Hits 241,409; miss 601.
CLASS exhaustive ⇒ W1∧W2=Walsh, no generation gap. Residual W1:
p=a²+64c². W2: named Auts are conjugate split involutions.

## 15.624 — Fable W2 strategy: inversion misses U; named W2 at p=11

`src/e1_gmin_m4_prop15624.py`. Fable: close W2 via full switching Aut,
not x/(x−τ); no more W1 χ-tower. Inversion −1/x is Max−, swaps {0,∞},
never in U. nuka PGL(2,11): 12 W2 hits, first (1,0,5,10) =
π(x)=x/(m(x+2)), m=(p−1)/2. Covers p=11. Disjunction with
x/(x−1) covers p=5,7,11,13, not 17,19,23. W2 p-law OPEN. Next W1
per Fable: quarter-interval + h(−4p), not another modulus.

## 15.623 — W1 for p≡73 or 97 (mod 120) via d=−3

`src/e1_gmin_m4_prop15623.py`. Six-point count: χ(5)=−1 ⇒ ε=1.
OpenAI PASS. W1 remaining: p≡1 or 49 (mod 120). nuka: p=11
π=x/(x−τ) over F_q has 1 eigen and 0 in U — pole family too thin
for W2 at p=11. leftover 2 open.

## 15.622 — W1 for p≡17 (mod 24); named W2 through the p=5 Φ3-gate

`src/e1_gmin_m4_prop15622.py`. d=−2: four-point count, ε=1 iff
χ(3)=−1 among p≡1 (mod 8), i.e. p≡17 (mod 24). OpenAI PASS.
W1 remaining: p≡1 (mod 24). Named W2 at p=5: Paley Aut
π(x)=x/(x−1) with switching y_k=χ(k−1) z(π k); Cy=−py, y∈U,
gcd(c(z+y),Φ3)=1. Same at p=7,13; y∉U at p=11. leftover 2 open.

## 15.621 — W1 for p≡5 (mod 8); PGL(2,q)·z is Φ3-dead

`src/e1_gmin_m4_prop15621.py`. Named stay d=−1, a=−λ^{-1}:
S Δ (S−1)={m,−1}, |QR∩Δ| odd iff χ(2)=−1 iff p≡5 (mod 8).
OpenAI math_review PASS. W1 now: all p≡3 (15.614) and all p≡5
(mod 8). Remaining W1: p≡1 (mod 8). W2: full PGL(2,q) orbit of
named z, xor z, is Φ3-dead at p=5 (86400 in-U, 0 coprime); z xor
ensemble still 72/155 gcd1, so the other endpoint is not a PGL
image of z. leftover 2 open.

## 15.620 — s_N is not a W1 p-law; χ_p-pullback misses Φ3

`src/e1_gmin_m4_prop15620.py`. ε(s_N)=0 at p=29 (f·n_odd^{QR} even).
Stay F_p-translates of z still have some a with ε=1 at p=5,13,17,29,37;
which a is not a p-law (α=(p+1)/2 fails at p=17). n_odd^{QR} odd iff
χ_p=1 off 0 (certified). Fable W2 deep_review: z-span Φ3-dead at p=5;
χ_p-pullback was the live named candidate — empirically Φ3|c and ε=0.
W2 needs a Max- object outside the halfspace-anti D-module. leftover 2
open. OpenAI unused this unit; Fable one deep_review only.

## 15.619 — odd_QNR(s_N)=0 is a p-law

`src/e1_gmin_m4_prop15619.py`. Off-0 nsq fibers: n_ζ=(p−1)/4 ± b/2
for p=a²+b², both even (Gauss |J|=√p + Z[i] UFD: Re ∈ {±a,±b}
forced to ±b). 0-fiber already even. So odd_QNR of any φ-pullback
vanishes, including s_N. f·n_odd^{QR}≡1 still certified only.
W1 p≡1 / W2 / Walsh / leftover 2 open. OpenAI referee math_review
PASS-WITH-NOTE (Tr=+Tr(α); J associated to a±bi). Fable unused.

## 15.618 — Φ=ε p-law; s_N is a φ-pullback; 1_M coprime to g

`src/e1_gmin_m4_prop15618.py`. Fable xhigh math_review PASS-WITH-NOTE
on Claim A: ε(w)=odd_QR(w)+odd_QNR(w) on W_0 for every odd p. Scale
is Φ((D−I)γ)=|supp γ ∩ (QR∪QNR)|=2p−1≡1 (not a p≡3 transfer from
z+Dz). s_N=f∘φ is a p-law (halfspace symmetric differences on
φ-fibers). odd_QNR=0 and f·n_odd^QR=1 certified p=5,13,17 (every
QNR fiber has even odd-index count); not yet a p-law. 1_M has
content X+1, gcd(X+1,g)=1, so the membership test is sound;
(D−I)γ is not a U-difference. Named U-diffs miss the p=5 Φ3-gate.
W1 for p≡1 / W2 / Walsh / leftover 2 stay open. Fable directions:
cyclotomic numbers of order 2 for the s_N parities; per-factor
W2 witnesses glued by primitive idempotents of g; AP-supported
U-diff with geometric Mattson–Solomon.

## 15.617 — NSQ stay-sum; correct W2 test; Walsh p=11 withdrawn

`src/e1_gmin_m4_prop15617.py`. Fable: f(D)w≠0 is not w∉(f)R;
w∈(f)R iff f divides the γ-content. z+Dz has Φ3|c at p=5,7 and
misses Φ5/Φ15 at p=11, so 15.616 Walsh-at-p=11 is withdrawn.
W2 is generic (72/154 U-diffs at p=5). W1 for p≡1: s_N = F2-sum
of z+T_a z over nsq a∈F_p^× that stay; ε(s_N)=1 at p=5,13,17
(Zolotarev class). Construction is a p-law; ε-value not yet.
Fable strategy: DFT/Mattson–Solomon for W2; Zolotarev class for
W1 p≡1 (QR sum fails, nsq sum works). leftover 2 OPEN.

## 15.616 — W2 via z+Dz; Walsh at p=11 (withdrawn by 15.617)

`src/e1_gmin_m4_prop15616.py`. The Krylov-gcd test that put z+Dz in
ker g was wrong. Correct test: w∈(f)R iff f(D)w=0. For named z,
f(D)(z+Dz)≠0 for every irred f of g at p=5,7,11. Fail: Φ3(D)(z+Dz)=0
at p=5. At p=11 (≡3) W1 is 15.614 and W2 is this, so I_U=W_0 by
15.612 CLASS: Walsh at p=11. Not ∀p. W1 for p≡1 and leftover 2 open.

## 15.615 — two-fiber W1-1 killed; W2 named-pool miss; L2 not closed

`src/e1_gmin_m4_prop15615.py`. Three Fable xhigh queries. (1) Two-fiber
ε=1 iff p≡1 is **false** at p=17 (ε=0); Fable's (p+1)/2 guess dies
there. Named z still has stay diffs with ε=1 at p=5,13,17, no uniform
a. (2) W2: named z+Dz / two-fiber / Frob / stay of this z miss g at
p=11. Need new U-differences. (3) leftover 2: Fable's p=5 k=20
min_+=2 scan is already 15.528 (empty leftover+splus). residual_ii
False. Walsh OPEN.

## 15.614 — W1 for p≡3; named D-spans miss W2 (split)

`src/e1_gmin_m4_prop15614.py`. Fable: ε((1+D)v)=g(1). Lift
v=x+a e_0+b 1_QR of named z-bits into W (φ=L∘σ^{−1} nsq ker ⇒
wt=(p−1)/2 on every square line). Then ε(z+Dz)=v_0=p(p−1)/2
(mod 2)=1 iff p≡3. W1 is a p-law on that class. For p≡1 the
two-fiber {φ=(p−1)/2}∪{φ=p−1} is in W_0; ε=1 certified, not
proved. At p=11, g divides ann(z+Dz) and ann(two-fiber): these
vectors miss every g-orbit. Walsh / residual_ii False.

## 15.613 — W1 named Max- in U; one ε-bit per p mod 4 (split)

`src/e1_gmin_m4_prop15613.py`. Paley halfspace-anti z (15.254 of ρ=1
h) lies in U for every odd p. ε(y+Dy) is CONSTANT on U (affine:
(D-I) kills ⟨1⟩ and W_0 maps into ker ε). Under W≅F2[M],
γ=1_{H∪(1+H)} is 1_M, and ε(w)=∑_{k odd} g_k for w=g(D)γ.
Census: that U-constant is 1 iff p≡3 (mod 4) (p=3,5,7,11).
For p≡1, α=(p+1)/2 and σ with L(σ^{-1})=p-2 (exists, 15.604)
gives T_α z∈U (stay is a p-law: −L(σ^{-1}α)=1∈S) and
ε(z+T_α z)=1 at p=5,13. Fable xhigh: construction/stay/constancy
PASS; ε-value BLOCK. W1/Walsh/residual_ii stay open.

## 15.612 — Walsh ⇔ W1 ∧ W2; CLASS of Aut-invariant ideals (proved split)

`src/e1_gmin_m4_prop15612.py`. dir(U) is Aut({0,∞})-invariant, so
I_U=dir(U)∩W_0 is an ideal of R=F2[X]/h. Walsh ⇔ I_U=W_0.
Maximal proper Aut-invariant ideals are (X+1)R=(D−I)W_0 and
(f_O)R for each ⟨I,Frob⟩-orbit of irred factors of g
(Fable xhigh PASS on CLASS). Hence Walsh ⇔ W1 ∧ W2.
W1 (some U-difference has (X+1)-valuation 0) is certified p=3,5,7
and **not** a p-law (Fable BLOCK): Frob of one U-point has odd
ε at p=3,7 and even at p=5; translation-stay fills W_0 at p=5,7
but not p=3. W2 vacuous at p=3, implied by W1 at p=5,7, first
live at p=11 (orbits {Φ_3},{Φ_5},{two Φ_15 quartics}).
Walsh / residual_ii stay False.

## 15.611 — W ≅ F2[X]/(X^N+1); ker2 dim 2 is a p-law (proved)

`src/e1_gmin_m4_prop15611.py`. Even nsq-line invariants W^H ≅ F2[F_p^×]
(regular; restriction to F_p^×, f(0)=∑f). M transits nsq lines with
Stab=F_p^×, so W ≅ Ind F2[F_p^×] ≅ F2[M] ≅ F2[X]/(X^N+1) as
D-modules. Unique D-invariant hyperplane is W_0=im(D−I). p odd ⇒
4|N=2^a m with a≥2, so dim ker((D−I)^2)∩W_0=2 for every odd p
(upgrades 15.610 C from certified to a theorem). Fail: minpoly(D)
degree <N; fail: W^H simple as C_p at p=7; fail: ker² dim 1.
Fable xhigh PASS (0.93). Walsh is still F2[D]-ideal generation by
all U-differences. residual_ii False.

## 15.610 — Aut({0,∞}) uniqueness for Walsh is DEAD (proved kill)

`src/e1_gmin_m4_prop15610.py`. W_0={w∈W: w_0=0}=extra^⊥∩W.
I D I^{−1}=D^{−1}; in char 2, D^{−1}−I=D^{−1}(D−I), so each
ker((D−I)^k)∩W_0 is I-invariant. 4|N for odd p; ker((D−I)^2)∩W_0
has dim 2 (certified p=3,5,7,11), a proper Aut({0,∞})-submodule
strictly larger than ⟨extra⟩. Same role as Aut_e reducible at p=5:
one U-difference outside ⟨extra⟩ does not force Walsh. Fable xhigh
BLOCK on irreducibility of W_0/⟨extra⟩ (unipotent flag is a p-law
shape; dim 2 certified). Walsh is now F2[D]-ideal generation by
the full U-difference set. residual_ii False.

## 15.609 — I(H0)=H0 for every odd p (proved)

`src/e1_gmin_m4_prop15609.py`. Square and nsq ∞-circles have different
directions, so they meet in two points and H0'=H0^⊥=rowspan(S).
Distinct F_p-sublines meet in 0, 1, or 2 points; tangency (|∩|=1)
forces equal χ-type (PGL-normalize a flag to (P¹(F_p), ∞); the
other circle is then parallel to F_p, hence square). Mixed-type
pairs therefore meet in 0 or 2 points. I sends off-0 square
∞-circles into ker S'=rowspan(S), so I(H0)=H0. Fable xhigh
deep_review PASS. Walsh spanning still OPEN. Fail: a square–nsq
tangent pair; fail: I maps a square row into rowspan(S').

## 15.608 — square/nsq PSL-orbits of F_p-sublines; 1∈dir(U) (proved)

`src/e1_gmin_m4_prop15608.py`. Changing the basepoint of an F_p-subline
replaces the direction b by b^{−1}, so χ is well-defined. Setwise
Stab of P¹(F_p) in PGL(2,q) is PGL(2,p) (unique Möbius on three
points), inside PSL(2,q). All circles are one PGL-orbit; PSL
normal of index 2 with Stab⊂PSL splits that into two equal
PSL-orbits. PSL preserves χ-type (PSL_∞ is square dilations). I(z)=1/z lies in PSL
and preserves each orbit. I(H0)=H0 is 15.609 (this unit only certified it).
U is antipode-closed, so 1∈dir(affine_span(U)). Walsh spanning of
V/⟨1⟩ still OPEN. Fail: one PSL-orbit; fail: χ flips with basepoint.

## 15.607 — W irreducible as G_aff^□-module, all odd p (proved)

`src/e1_gmin_m4_prop15607.py`. F_p^× ⊂ M (15.598 B) preserves every
F_p-line through 0, hence each W^H. On F_q/H ≅ F_p it is Aut(C_p),
which transits the irreducible factors of Φ_p (the orbits F_p^×/⟨2⟩).
So W^H is irreducible for C_p ⋊ F_p^× even when 2 is not a primitive
root. A G_aff-submodule meets some W^H in 0 or all of it, and M
transits the summands, so W is irreducible. Antipodes put ⟨1⟩ in dir(affine_span(Max−)),
and G_aff^□ permutes Max−, so that direction is H0. Fail: W^H
simple as a C_p-module at p=7
(two cubics; F_p^×-span of either kernel is all of W^H). Walsh is
still the xor-slice; residual_ii stays False.

## 15.606 — nsq line averages split W; M transits (proved)

`src/e1_gmin_m4_prop15606.py`. π_H=∑_{h∈H} T_h on an F_p-line through
0. Square H: π_H=0 on W (sums on square affine lines). Nonsquare H:
W^H = all even H-invariants, dim p−1 (a square line meets every
H-coset once). The (p+1)/2 nsq projectors are orthogonal and sum
to I_W, so W=⊕ W^H. M=(F_q^×)² acts transitively on nsq directions
(Singer: F_q^×/F_p^× ≅ C_{p+1}, squares the index-2 subgroup).
If 2 is a primitive root mod p, Φ_p is irreducible over F2, each
W^H is a simple C_p-module, and every G_aff^□-submodule of W is
0 or W. Not a p-law: p=7 has ord_7(2)=3; G_aff-spin of either
cubic kernel still fills W (census). Walsh / residual_ii unchanged.
Fail: square π_H of rank p−1; fail: M mixes square with nsq.

## 15.605 — Paley A²=A over F2; H0=⟨1⟩⊕W (proved)

`src/e1_gmin_m4_prop15605.py`. Paley graph of order q=p² is srg with
(q−1)/4 even, so A²=A over F2; P=A+I+J likewise. Fail: Paley of
order 13 (q≡5 (mod 8)). The 15.604 extra vector (1_QR or 1_QNR)
has affine translates spanning a complement W of ⟨1⟩ in H0, dim
N=(q−1)/2: ker S_aff has dim N by 15.600, the 15.601 pencil writes
Paley neighborhoods as F2-sums of square lines (p≡1) or forces
Af=0 then Pf=f (p≡3). W is G_aff^□-invariant, generated by extra,
and has no trivial submodule. Irreducibility of W OPEN. Walsh /
residual_ii unchanged. Backend: serial F2; GPU unused.

## 15.604 — 1_QR ∈ H0 iff p≡1 (mod 4); ker(D−I)∩H0 dim 2 (proved)

`src/e1_gmin_m4_prop15604.py`. Square-line counts (15.598 A,B): a
0-line in square direction has L\\{0}⊂QR so |L∩QR|=p−1 even; an
off-0 line has Σ χ_q = −1 so |L∩QR|=(p−1)/2 even iff p≡1 (mod 4).
Thus 1_QR ∈ H0 iff p≡1, 1_QNR ∈ H0 iff p≡3; never both
(e_0+e_∞ ∉ H0). Square dilation D of order N=(q−1)/2 has ambient
ker(D−I)=⟨e_∞,e_0,1_QR,1_QNR⟩, and the two line types cut this
to a 2-space ⟨1, extra⟩. Restriction H0→F2^{QR} is **not**
surjective (ker 2,3,6,12 at p=3,5,7,11) — not a p-law; do not
name minpoly(D)=X^N+1 from a false onto map. Irreducibility of
H0/⟨1⟩ OPEN. Walsh / residual_ii unchanged. Fail: swap
congruences; fail: dim ker(D−I)∩H0 = 1 or 4. Backend: serial F2;
GPU unused.

## 15.603 — H0 ∩ H0' = ⟨1⟩, H0 + H0' = even-weight (proved)

`src/e1_gmin_m4_prop15603.py`. S' = nonsquare-direction ∞-lines;
rank(S')=n/2 by the same radical as 15.600. If x ∈ H0 ∩ H0' then every
affine F_p-line has the same sum, so f is constant on AG(2,p) (p odd:
the p+1 lines through a point). Intersection ⟨1⟩. Then
H0+H0'=even-weight. The PSL-heart E/⟨1⟩ splits as
(H0/⟨1⟩)⊕(H0'/⟨1⟩) of dim (q−1)/2 (Mortimer: heart not absolutely
simple for PSL(2,q), q odd). Irreducibility of each summand OPEN.
Walsh / residual_ii unchanged. Fail: H0=H0'; fail: non-constant
common kernel vector. Backend: serial F2; GPU unused.

## 15.602 — G_aff^□ on H0=ker S (proved; Walsh still OPEN)

`src/e1_gmin_m4_prop15602.py`. Translations, square dilations, and
Frobenius permute the rows of S (fail: nonsquare dilation). Unique
1-dimensional G_aff^□-invariant subspace of H0 is ⟨1⟩ (translations
force constancy on F_q; even weight). Inversion permutes the
square-pencil through 0, not every row; H0-invariance certified
p=3,5,7,11 (not a p-law). If H0/⟨1⟩ is irreducible then
dir(affine_span(Max−))=H0; that irreducibility is OPEN (cyclic-full
at p=5..31; p=3 hits ⟨1⟩ as expected). Walsh 15.406 E and residual_ii stay False.

## F̂ is not a Paley-field square or field-norm (leftover 1 still OPEN)

`src/e1_gmin_qvar_fhat_norm.py`. F̂=(λ−6)q² on even ψ∉{1,χ}.
At p=5 every named λ/13 gives F̂ with v_13=−1; 13 has residue
degree 2 in Q(ζ_{24}) and is inert in Q(√5); min 1250/13 is not a
Q-square. At p=7 every named λ/409 has v_409=−1; f=2 in Q(ζ_{48})
and inert in Q(√−7); 409 does not divide p(p²−1)(p²+1). Fail:
Q-square; fail: Bochner F̂=|A|² with A in the Paley character field
(valuations would be multiples of f). Positivity of F̂ stays OPEN.
Backend: serial Fractions (inherently); GPU unused. Residual flags
untouched.

## Walsh is NOT closed (single-orbit proof BLOCKED)

A proposed general-p Walsh proof (Aut_e-orbit of one U-difference spans
W because Aut_e-submodules would be PSL-submodules) was falsified by
both referees. Aut_e is the edge stabilizer, not PSL; it is reducible;
at p=5 a single orbit has dir ≤11 < 12=dim W. Do **not** flip 15.406 E
or residual_ii. Repair must use the full set of U-differences, not one
orbit.

## 15.601 — QR indicator in rowspan(S) or S+ℓ (proved)

`src/e1_gmin_m4_prop15601.py`. Pencil of the (p+1)/2 square-direction
lines through 0: Sᵀw = ((p+1)/2)(e_0+e_∞)+1_QR. Hence
p≡1 (mod 4) ⇒ 1_QR+ℓ ∈ rowspan(S); p≡3 ⇒ 1_QR ∈ rowspan(S).
On H0, QR·x equals ℓ(x) or 0. Aut_e-invariant extra duals of the
xor-slice are empty. Walsh spanning still open.

## 15.600 — rank(S)=n/2 for every odd prime (proved)

`src/e1_gmin_m4_prop15600.py`. The 15.599 Gram gap is closed.

- \(1\in\ker S\cap(\ker S)^\perp\): \(S1=0\) (\(|v_L|=p+1\) even) and \(1\)
  is a parallel-class sum of rows.
- Over \(\mathbb F_2\), \(\mathrm{rank}(SS^\top)=\mathrm{rank}(S)-\dim(K\cap K^\perp)\).
  Radical at least 1 plus 15.599 A,B forces \(\mathrm{rank}(S)=n/2\).
  Fail: \(\mathrm{rank}(S)=n/2-1\).
- \(\dim H_0=n/2\) is now a theorem, not a census. PΓL-cyclic modules of
  \(H_0\) are full at p=3,5,7 (8/8); Aut_e still reducible. Walsh is
  still spanning of the xor-hyperplane. residual_ii stays False.

Backend: identities (serial, inherently); rref cross-check p=5,7,11.

## 15.599 — rank pin, antipodes, Aut_e reducible (Walsh still open)

`src/e1_gmin_m4_prop15599.py`. Max-free F2 geometry of the square-line
incidence S.

- **Proved:** rank(SSᵀ)=n/2−1 (block-diagonal J−I, (p+1)/2 classes);
  class-sums = 1 ⇒ rank(S)≤n/2; hence rank(S)∈{n/2−1, n/2}.
- **Certified p=3..37:** rank(S)=n/2, so dim H0=n/2.
- **Antipodes:** y↦−y preserves Max− and the xor. The p=11 eps1 half
  (y_∞≡+1) has dim 60; with complements dim H=61=n/2 and dim U=60.
  15.596 compared a half-ensemble to n/2.
- **Killed:** Aut_e-irreducibility on H0 / xor-slice (proper cyclic
  modules at p=5,7,11). Line-flip of a square-line block is not Max−
  (exterior determines the signs on S). A single Aut_e-orbit of a U-point
  spans the slice at p=7, not at p=5.
- Walsh spanning still open. residual_ii stays False.

Referees (suggest_direction, both houses): do not more-census p=11
without antipodes; do not reopen Aut_e-irreducibility. ProcessPool
W=11 for ranks; GPU unused (F2).

## 15.598 — square-direction affine lines cut Max− (proved; Walsh still open)

`src/e1_gmin_m4_prop15598.py`. Max+-free character sums, then Cy=−py.

- Jacobi: \(\sum_{x\in\mathbb F_p}\chi_p(x(x+\delta))=-1\) for \(\delta\neq0\).
- \(\chi_{p^2}(z)=\chi_p(N(z))\). Off an affine line \(L=a+\mathbb F_p b\),
  \(\sum_{i\in L}C_{ij}=-\chi(b)\). Square direction \(\Rightarrow\)
  \(\sigma_{\mathrm{out}}=0\) and on-\(S\) row-sums \(=p\).
- Hence every Max− vector has \(\sum_{k\in S}y_k=0\) on
  \(S=\{\infty\}\cup L\), \(L\) square-direction. So
  \(\langle x,1_S\rangle=(p+1)/2\), even for \(p\equiv3\pmod4\), odd for
  \(p\equiv1\pmod4\). Fail: Max+ (many values); fail: nonsquare direction
  (not identically 0).
- Pair-slice \(U\) is the xor-hyperplane \(x_i+x_j=c\) of \(\mathbb F_2^n\).

Walsh (15.406 E) is now: affine_span(\(U\))=H\cap{ℓ=c} with
H=affine_span(Max−). Certified dim \(U\)=dim \(H\)−1 at p=3,5,7 (full)
and p=11 (full rank(\(B_U\))=60, 200k-sample dir(\(H\))=60, **not**
\(n/2=61\); 15.596's comparison to \(n/2\) assumed the p=5,7 value).
Square-line indicators with \(\chi(\mathrm{dir})=+1\) span the dual at
p=3,7 and their differences span it at p=5. Spanning for general p is
the remaining Walsh step. residual_ii stays False (leftover-only / 5+).

Scripts: `scripts/walsh_gf2_dual.py`, `walsh_subline_dual.py`,
`walsh_affine_line_dual.py`. Evidence jsons alongside. p=11 full
37,457,112×132 line popcounts: const_even=66, mixed=66, const_odd=0,
matching \(\frac12 p(p+1)\). Backend: field Jacobi serial (inherently
sequential, \(p<80\)); Max− inner products numpy; p=11 stream one
process over mmap (bandwidth, not 86-way GE). GPU unused (bit geometry).

Do **not** reopen Paley \(E_-[S^2]<20+12/p\) (false). Do not flip
residual_ii from Walsh alone.

## Since `5ce0258` (do not lose)

Literature (shipped): `evidence/HISTORY_AND_REFERENCES.md`. MO 413935 /
X prize / Paata SPbU+Volberg / Littlewood 4/3 / Paley 1933 /
Goethals–Seidel 1967 / Blei 34 / DMP / Talagrand. No classical
existence theorem. Do not reopen BH, typical-\(A\), Bowlin, or
“Paley \(\Rightarrow\lim=1/2\)”.

Wiring: `test_pairing_open_and_flag_imported` now ORs theorem I (shipped
flag already did). Dump JSON has equal-density \(=\bar\lambda\).
`src/e1_gmin_qvar_bool6.py` stays **uncommitted**.

GLOBAL QVAR still False. leftover-1 = global QVAR \(\land\) R1.

**Killed this stretch (scratch; not a new identity file):**

- **k=1 mass does not lift the floor.** \(E|Z|^2\ge n_{1d}S^2/|\mathrm{Max}+|\)
  beats QVAR iff \(|\mathrm{Max}+|\le n_{1d}(q-1)/3\). At p=5 that bound
  is 240 vs live 260: \(30\cdot900/260<225/2\). At p=7 it is worse
  (\(n_{1d}=140\), \(|\mathrm{Max}+|=11452\)). Named k=1+k=3 counts still
  need an upper bound on \(|\mathrm{Max}+|\). For \(p\ge13\) the bulk is
  the **top stratum** \(k=m=(p+1)/2\ge7\), not k=1. (p=41 k=7 is thin
  vs \(n_{1d}\); p=13 k=7 **is** the top stratum.)
- **Reconstruction inner product is the pairing.** With
  \(h_j=(\sigma_j-\varepsilon)/2\), \(F=\sum_j h_j\circ t_j\),
  \(F_w=\sum_j w_j h_j\circ t_j\): \(\sum_x F F_w=pZ_\psi\). Cross-direction
  terms vanish (affine plane: each pair of parallels from different
  square directions meets once; \(h_j\) mean-zero). \(\sum F^2=pS\) is
  the already-named \(\sum a_L=S\). Does not sign \(\hat F(\psi)\).
  Same as “B-weighted 15.588 tautology” in the kill list below.
- **\(\chi*1_D=p1_D-\frac{p-1}{2}1\) is linear** (15.317 A). It names
  2-point / Ω-support, not 4-point. L2 of the Fourier form
  \(f*f=-2p f\) on Ω is an identity, not a floor.
- **Independent \(\widehat N(\xi)\) would miss QVAR.** Named
  \(\mathrm{Var}(\widehat N)=q^2/4\) (15.305 A). Independence gives
  \(E[\Delta^2]=(q-1)q^2/8=2q^2(q-1)/16\), below the floor
  \(3q^2(q-1)/16\). Nyquist pairing needs **positive** same-sign
  frequency correlation (deficit on \(K\), not a Wick/independent lag).
- **Wrong quantifier:** \(\min_y Z^2=0\) (attained). QVAR is the
  **uniform** average on Max+, not a pointwise SOS of \(Z^2-\)threshold.
  Dirac on a \(Z=0\) vector is a different measure.
- **Orbitwise QVAR is false** (already in NOTE 2026-08-20): p=7 PSL
  orbit of size 1176 has \(Z_\psi=0\). Mixing \(k\) inside one orbit is
  allowed; mixing **orbits** is required. No named orbit sizes.
- Type A (2+2 on two square lines): Wick isotropic; off-diagonal
  Nyquist of Type A is **negative**; QVAR is carried by the
  same-direction block of \(\Pi_i=2a_i-p(p-1)/2\). Paley-type mix of
  Type A (sign of \(n_{++}^{\mathrm{same}}-n_{++}^{\mathrm{opp}}\))
  **reverses** p=7 vs p=11,13 — not a p-law. \(m_4\) splits inside
  Paley type (2–5 means). Occupancy 4th moments have Max+ denominators
  (409 at p=7).
- 3-point of \(1_D\) is **not** constant on Max+ (15.468 C, 16
  fingerprints at p=5). Not Wick. Gale–Ryser scratch used a reversed
  conjugate; those rates are **invalid**.
- Occupancy-energy probe on p=5/7 caches used the **wrong D-slice**
  (dropped infinity; \(k=9\) vs 10 / \(20\) vs 21). Discard.

Do **not** add another equivalent identity with `inequality_proved=False`.
Next constraint is still simultaneous Boolean ridge / Gauss 4-distinct
pairing on the size-weighted mix, or an odd-coset 4-harmonic that mixes
\(k\) **and** orbits. Do not split \(\lambda=0\). Do not require each
k-stratum. H/I Nyquist of occupancy covers \(p\equiv3\) only.

## Class-function route (d512824–HEAD) — leftover 1, not a close

Plan: `evidence/PLAN_2026-08-22_class_function_route.md`.  Scripts:
`gamma_class_function.py`, `gamma_class_p7_gate.py`,
`gamma_conjugacy_classes.py`, `gamma_class_parameter.py`,
`gamma_ae_fourier.py`.  Pointwise constituent energy is dead at p=7.
`Γ_δ` quantization is a p=5 artifact.  `Γ(−g)=Γ(g)` is a theorem.
Order+fix does not determine `Γ`.

Step 3–4: on PSL, elliptic `Γ=0`, involution `Γ=2(n−2)` (p=5,7),
unipotent `Γ` is determined by `λ_exc`.  Split `Γ` is a function of
`τ=tr²/det`; the p=5 O(1)-in-p split recipes **fail** at p=7.
`A_e` is the principal series `ρ(α_k)` with `4|k` (count identity
`(q−9)/8` for all `q=p²`; inner product 0/1 at p=5,7).  Fourier
inversion `λ_c=⟨Γ,χ_c⟩` reconstructs `Γ` with 0 mismatches on all of
PSL at both primes.  Principal floor binding value is p=5 `λ_4=80/13`.
p=11 stored `Φ` has min principal `8.054`, `λ_exc=8.664`.  No flag
flipped.  leftover 1 remains `λ_exc≥6` and `λ_k≥6` for those `k`.

Step 6: Aut-orbit values of L named at p=5,7
(`scripts/aut_orbit_L_and_lambda_fit.py`).  Binding leftover orbit is
the unique min (p=5 order-12 (−−); p=7 order-8 (−−)).  **Killed:**
L(i)=L(ω₃) as a p-law (true at p=5,7, false at p=11 sample);
p-independent cosine of λ(k) (p=5+7 does not predict p=11).  p=5
only: λ(k)=8+8/13−(64/13)cos(πk/6).  Uniform L≥L_min is not a floor
proof (binding character eats the slack).  Aut leftover dofs stay
n_orb−2.

Step 7: Boolean cubic on Ω **does not cut leftover Aut-dofs of Q**
(`scripts/boolean_cubic_orbit_relations.py`).  It is Fourier of
\(z_x^2=1\): `2ẑ(0)ẑ+∑_R ẑ(tξ)ẑ((1−t)ξ)=0` to 1e-12; `2pẑ+∑B` only
on y_∞=+1.  Squared form is the tautology 4p² Q(r)=E[|∑B|² u(r·)].
Off-diagonal bilinear Gram Γ is not linear in Aut-orbit Q (maxerr
757 / 1226 at p=5,7) and splits Aut-ratio buckets at p=7.  ρ(R)
misses leftover orbits.  Do not reopen as a floor argument (15.279 T
already: Boolean rewrite returns M).

Step 8: Boolean 4-point of V_+ is Aut-constant (42/42 p=5, 128/128
p=7) and **not** a function of {κ,φ,star} nor of (κ,CR,star)
(`scripts/m4_aut_orbit_vplus.py`).  Finite m₄ values are odd over
N/4 (ten at p=5, twenty at p=7).  ⟨m₄,κ_A⟩ lives on finite 4-sets
only; binding pairing 0.180 at p=5.  Character-sum of those orbit
values is 15.48 still open.  Do not add an identity file.  No flag
flipped.

Step 9: q-dependent split-Γ formula **fails the p=11 gate**
(`scripts/split_gamma_dilation_ansatz.py`).  Dilation t=−1 has
ensemble Γ=2(n−2) at p=5,7, not pointwise ±2p.  Cosine
a_j=A+B/q predicts p=11 λ_min≈5.6<6.  Re J(χ,α_k) is not linear
in λ−8 even at p=5,7.  Paley+order grouping of Γ(t) is a p=5,7
artifact.  No identity file.  No flag flipped.

Step 10: involution E[s²] splits as 2-point mass 2+2q and
4-distinct mass 2(n−2).  Mean m₄ on {a,−a,b,−b} is 2/(q−3), not
constant.  2-point cannot prove Γ(−1)=2(n−2).  Kloosterman/Bessel
are not a formula for Γ(t).  No flag flipped.

## 15.597 Theorem A* (proved — not a census)

`src/e1_gmin_m4_prop15597.py`. Distinct from 15.108's old "Theorem A*"
(16N). For every prime \(p\ge5\), on \(Z\):

\[
\Phi_{\mathrm{part}}=\bar\lambda\,I,\qquad
\bar\lambda=8(n-2)/(n-6).
\]

Closed-form contractions (lemmas L1–L3 on \(W\in Z\)): \(\sum\kappa t=(n+1)\|W\|^2/4\),
\(\sum\phi t=-n\|W\|^2/4\), \(\sum\star t=-p\|W\|^2\). Direct quadratic-form
checks at \(p=5,7,11,13\) match. The particular solution is spectrally
**invisible**: all deviation of \(\Phi\) is \(\Phi_\delta\).

**Sharpened leftover 1** (equivalence, inequality still open):

\[
\lambda_{\min}(\Phi)\ge6
\iff
\Phi_\delta\succeq -\frac{2n+20}{n-6}\,I
\quad\text{on }Z.
\]

The \(n/12\) R1 bound is the multiplicity-floor op-norm form of the same
statement. Global QVAR is the exceptional isotype of \(\Phi_\delta\);
principal room is the rest.

**Corollaries (proved, buy nothing toward 6):** \(\mathrm{tr}(\Phi_\delta)=0\),
so \(\lambda_{\min}(\Phi_\delta)\le0\). \(\Phi\) is Gram, so \(\Phi\succeq0\)
and \(0\le\lambda_{\min}(\Phi)\le\bar\lambda\). Target 6 sits strictly inside
that window. The remaining content of leftover 1 is pushing the proven
lower bound from 0 up to 6; nothing short of a genuine bound on \(\delta\)
crosses it. Do not re-derive Gram \(\ge0\) or tracelessness as a floor
argument.

Do **not** add another equivalent identity with `inequality_proved=False`.
A* removes \(m_4^{\mathrm{part}}\) from the spectral problem; it does not
bound \(\delta\in\ker(4pI-T)\).

**Killed (2026-08-22):** Aut-invariant 4-point master equation
\(Tm=4pm-4\kappa/p\) plus \(|m_4|\le1\) cannot prove QVAR.
`src/e1_gmin_qvar_box_master.py`: permutation Aut quotient of \(T\) is
exact (\(T\kappa=-6\star\)); at \(p=5\), \(\dim E_{4p}^{\mathrm{Aut}}=2\)
and \(\min\langle m,\kappa_{A_\psi}\rangle=-285/4<0\), while true Max+
pairing is \(+14.13\). Same kill at \(p=7\) (ker 7, LP min \(\approx-2708\)).
Need a constraint outside linear 4-point theory (15.589 I).
Checked and redundant/useless for recovering the sign: \((1^\top y)^4=(1+p)^4\)
is already in the master affine space; \(P_{E_{4p}}\kappa_{A_\psi}\) is
indefinite (both signs at \(p=5\)); \(\ell^1\) CS overshoots (\(-172\) vs need
\(-6.75\)). Next constraint is simultaneous Boolean ridge / 6-point coupling,
not another 4-point rewrite.
Linear 6-point contractions + box still miss: p=5 min \(-101/4\); p=7
15.590 joint deg-6 ker 4, box min \(-10633/8\). SOS-4 along the deg-4
kernel is feasible at pairing \(-45/4\). True pairings positive. Need
strictly stronger than SOS-4 / linear 6-point (overlapping Boolean
support or SOS-6). Local Boolean-6 (each 6-set a \(\{\pm1\}^6\) moment)
is **not a p-law**: p=5 min \(=+27/4\) equals the particular pairing;
p=7 HiGHS-IPM min \(\approx-172.75<0\) while true pairing is positive
(uncommitted `src/e1_gmin_qvar_bool6.py`, do not treat p=5 positivity
as general).

**Equal-density is the unproved ordering, not a weaker sufficient.**
Mean \(Q\) (or deficit) equal on fourth-powers-off-\(\pm1\) vs the
complementary squares gives \(B=4q^2(q-1)/(q-5)\) and
\(\lambda_{\mathrm{exc}}=8(q-1)/(q-5)=\bar\lambda\). That is theorem F
(exceptional above the \(\Phi\)-mean), strictly stronger than QVAR.
Singer/OA Nyquist pairing does not cover \(p=13\equiv1\pmod4\) (\(m\) odd).

Do not reopen as QVAR proofs (already killed / tautological):
- Unique \(G\)-invariant 4-harmonic: \(\dim\mathrm{Harm}_4(W_e)^G=2,3,6\)
  at \(p=5,7,11\) (`theorem_K_harm4_not_one_dimensional`). \(E|Z|^2-V_{\mathrm{sph}}\)
  is not a multiple of leftover-3.
- Pointwise SOS of \(Z^2-\mathrm{threshold}\): \(Z=0\) is attained.
- Spectrum of real \(A_\psi|_{V_+}\) is \(\{0\}\cup\{\pm\sigma\}\)
  (\(\sigma=p/4\) at \(p=7\); \(\sqrt5/2\) at \(p=5\)). Converts QVAR to
  an imbalance \(E[(t-s)^2]\ge3(p^2-1)\) for \(p\equiv3\), same floor.
- Frozen ridge 2-design: \(F=Ky\) linear, \(\|F\|^2=pS\), frozen \(F\in V_+\)
  gives \(E[\langle y,F\rangle^2]=(p+1)S<3pS\). The gap is \(y\)–\(F(y)\)
  correlation (RidgeAD), still unbound.

## What shipped (wiring, not a close)

Leftover-1 conjunct is now **GLOBAL QVAR**, not per-stratum k≥7:

| File | Role |
|---|---|
| `src/e1_gmin_global_qvar.py` | Mixed-k floor identities A–I, P. `global_qvar_proved_general()` is A.proved and (P or G or H or I `inequality_proved`). All inequalities False. No handwritten True. |
| `src/e1_gmin_leftover1_qvar_principal.py` | `leftover1 = global_qvar AND principal_delta_room`. Does **not** AND `qvar_k_ge_7`. |
| `src/e1_gmin_qvar_k_ge_7.py` | Per-stratum k≥7 is **false** at (41,7) E=0 Cy=py and (13,7) pointwise. Flag False. Not a leftover-1 import. |
| `src/e1_gmin_r1_principal_pge11.py` | R1 L² recorded, **not** proved. `r1_l2_bound_for_p_ge_11` False. Interpolant 4/(p−3)² killed. |
| `src/e1_gmin_m4_prop15278.py` | `phi_F_ge_6` imports leftover1 AND. |

Gating tests: `tests/test_global_qvar.py`, `tests/test_leftover1_qvar_principal.py`,
`tests/test_qvar_k_ge_7.py`, `tests/test_r1_principal_pge11.py`.

## Live obstruction

Prove \(\Phi_\delta\succeq -(2n+20)/(n-6)\,I\) on \(Z\) (15.597). Split as:

- exceptional isotype = GLOBAL QVAR: \(\hat F(\psi)\ge 0\) / Gauss 4-distinct
  pairing of \(m_4\) on the full Max+ mixture (all \(k\), \(\lambda=0\) unsplit);
- rest = R1 / \(\|\delta\|^2\le n/12\).

Import `inequality_proved` only if that sign is actually proved.

Do **not** add another equivalent identity with `inequality_proved=False`.
Do not set leftover 1 True until GLOBAL QVAR **and** R1. Do not set
`qvar_k_ge_7` True. Leftover 2/3 / L stay their own units.

Killed as proofs: per-stratum k≥7, CS on \(\langle\delta,\kappa_A\rangle\),
two-level occupancy as a p-law (fails p=11 k=4/5), B-weighted 15.588
tautology (\(\langle F,F_w\rangle=pZ\)), energy-only (15.589 I), Delsarte,
cyclotomy \(E[E_0 E_r]\), 2-point fit-as-proof, pointwise SOS (Z=0
attained), k=1 mass lift, independent \(\widehat N\), \(\chi*1_D\) as
4-point, Type A Paley-type factorisation, orbitwise QVAR, linear 4-point
+ box, SOS-4, linear 6-point, Boolean-6 as a p-law.

Scratch Aut/Torb probes on nuka (`/tmp/qvar-nuka/`,
`/tmp/grok-goal-f38dc225339a/implementer/probe_nuka_*.py`) are **not
shipped**. Torb vs \(T|_V\) still undiagnosed; do not import as a theorem.

## Mesh (out of QVAR scope)

PR #4 merged `mesh/k6-p13-enum` **into main** (`99cbf09`), so
`scripts/maxplus_profile_enum/` now lives on main. The live enum campaign
is still not a leftover close. cpu44 was hard-closed (2072 `orb*.npy`
kept; stop flag still set). Do **not** restart cpu44 unless the user
names cpu44 in the same turn. Worktree
`/home/nick/quadratic-minmax-limit-k6-mesh` may lag main.

## Compute

nuka (5700X3D) for serial / vcache. lucky is DNS only. Never 86 workers
on Orin. Soft-stop mesh unless the user hard-closes a named node.

## Suggested skills

`agent-cost-optimization`, `graph-engineered-completion`,
`use-available-compute`, `verification-before-completion`, `handoff`,
`scientific-critique`, `grill-me`, `self-refine-loop`, `research`,
`arxiv`, `litreview`. Referees only if the user names them.

---

**Date:** 2026-08-21 (R1 L² unit recorded, **not** proved; **no flag flipped**)

`tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2` is identity (I) in
`TECHNICAL_NOTES.md` §4 / `METHOD.md`: leftovers 1 and 3 are moments of
one four-point tensor. It is proved (index split + E[(y·z)^2]=2n). It
does **not** close leftover 1: CS on dim Z, tr Φ, tr(Φ²) alone cannot
get λ_min≥6, and the exceptional block still needs GLOBAL QVAR (mixed-k;
per-stratum k≥7 is false).

Leftover 3 next sufficient target, not imported: p=5 is already a finite
from-C theorem (`type_I_p5_through_e_3AB_positive`). For every prime
p≥7, `|μ|≤2/n` on |κ|=1 is strictly stronger than `|μ|≤L=(p-2)/(2p^2)`
and would close Type I (`2/n < L` iff p≥7). Census slack at p=7 is
tight: 109/2863 vs 2/50 (~5%). Do **not** use `|μ|≤|f4|` (false at p=7,
15.191 I). `|μ|≤maj` remains false at p=7 and p=11. Triangle
`|R̄₄|+2|φ|` is too weak for `|μ|≤L` even at p=5. L2 conversion of
3A+B is rejected. p=5 μ=f4 on each (κ,φ) class; p=7,11 split inside
(κ,φ).

Leftover 2: leftover+splus empty all nF at p=5 k=20 (15.528);
leftover-only is not residual (ii). Walsh cannot flip leftover 2
(interior 4-level only). Uniform Paley E_-[S²]<20+12/p is false.

R1 (`‖P_{E_{4p}} m₄⁺‖² ≤ n/12`) is the binding leftover-1/3 face
(15.595). `src/e1_gmin_r1_principal_pge11.py` records exact measured
‖δ‖² vs R1 (exceeds at p=5,7; census-holds at p=11). p=5 measured
equals κ_hyp_δ (1536/65). The interpolant κ_hyp_δ·4/(p−3)² is **killed
as a retained δ-bound** (equality law false at p=7,11; Aut-dim
ν_G-ratio dies at p=7; no operator identity). `r1_l2_bound_for_p_ge_11`
is False. `principal_delta_room_moment_proved` imports that unit only.

`src/e1_gmin_leftover1_qvar_principal.py` is the Max+-free leftover-1
import. `phi_F_ge_6_proved_general` is `leftover1_qvar_and_principal_proved()`,
which is True only if **GLOBAL** QVAR (mixed-k) **and** the principal
δ-moment both hold for every prime p≥5. Both estimates stay False.
Identities A–D (QVAR iff, V_sph>threshold, D+ room formula,
⟨m4,κ_B⟩≥0 iff floor 6) have fail-eqs and tests. p=13 orbits are not
imported. Aut-Schur / Gsum / pairing False. L OPEN. Leftovers 2 and 3
False. Live `e1` is still the old AND.

---

**Date:** 2026-08-20 (current branch; use `git log -1` for the exact checkpoint)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit
**Current branch:** `codex/leftover-moment-attack` (use `git log -1` for the
current hash).  The latest continuation adds Prop 15.589 Theorems J--U,
including the all-prime `k=4,5,6` QVAR closures described below.  All 124
focused Prop 15.588/15.589 and k=7 probe tests pass.
**Statement:** [MathOverflow 413935](https://mathoverflow.net/questions/413935).
\(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\lvert\sum_{i<j}a_{ij}x_ix_j\rvert\).

## Binding status: exact PSL and high-stratum reduction

Work on `codex/leftover-moment-attack` gives the multiplicity-free
decomposition

`Z = W_e direct-sum ((p^2-9)/8 distinct degree-(p^2+1) principal series)`,

where `dim W_e=(p^2+1)/2`. Thus every Phi eigenvalue has multiplicity at
least `n=p^2+1` except one exceptional scalar. That scalar is exactly

`lambda_exc = 32 E|Z_psi|^2/[q(q-1)]`, `psi^2=chi`,

so its floor is the single quartic variance inequality
`E|Z_psi|^2 >= 3q(q-1)/16`. The remaining principal floor route is the
exception-removed delta room
`||delta||^2 <= n(n+10)^2/[6(n-14)(n-6)]`.  This is sharper than the old
room by the factor `(n-6)/(n-14)` and is equivalent to
`E[(y dot z)^4] <= 4n(3n^2-37n+2)/(n-14)`.  It applies after QVAR proves
`lambda_exc>=6`. Neither open inequality is proved generally; no flag is
flipped. See
`evidence/NOTE_2026-08-20_psl_and_stratum_floor_reduction.md` and Prop 15.589.

Two routes are now explicitly dead: restricted Phi does not have floor 6 on
every profile stratum, and at p=7 a full PSL orbit of size 1,176 has
`Z_psi=0`, so the quartic bound is not pointwise/orbitwise.

The pointwise repair also fails inside the first live stratum itself.  An
exact coupled coefficient/profile/Boolean CP-SAT model at `p=13,k=7` found a
translation-gauged Max+ representative with `Z_psi=-28-42i`, hence
`|Z_psi|^2=2548 < 10647/2`.  The stored support is independently checked
against `Cy=13y`, all seven active square directions, and the direct quartic
kernel.  Thus QVAR on `k>=7` must use the ensemble mixture even at its first
unresolved prime.  See `evidence/k7_p13_cpsat_{probe.py,witness.json}`.

Orbit averaging makes real progress on that counterexample.  Its signed-PSL
lift orbit is free and has size `4,826,640`, split evenly between the two
global signs.  Jellyfin's A380 evaluated all `2,413,320` epsilon-plus vectors:

```
E_orbit |Z_psi|^2 = 806468/85 > 10647/2,
lambda_exc(orbit) = 19088/1785 > 6.
```

The orbit contains `14,196`, `28,392`, and `2,370,732` vectors at activities
`k=5,6,7`.  Crucially, its `k=7` slice alone clears QVAR with exact moment
`1606124/167`.  Independent profile reconstruction finds 12,852 depressed
nonzero-quintic representatives, exactly 1,071 in each of the twelve scalar
classes.  This does **not** close `p=13,k=7`: after forbidding all 1,071
scalar-7 representatives, independent eight-worker CP-SAT runs on Nuka and
Jellyfin found a second-orbit seed with `Z_psi=-132-198i` and
`|Z_psi|^2=56628`.  That seed also has a free 4,826,640-vector signed orbit.
Its full epsilon-plus orbit has the same exact mean `806468/85`; its
2,384,928-vector `k=7` slice has moment `198692/21`, again above QVAR.  The
next finite task is therefore iterative orbit decomposition, not pointwise
minimization.  Both packed orbits and the second census are in the platter
backup named below.

The exceptional target has since narrowed again.  For `p=3 mod 4`, `Z_psi`
is exactly a signed sum of nonnegative directional profile energies whose
pointwise total is `p(p^2-1)/4`.  Combining that identity with the affine
profile classification, and using the Euler-product lower bound for
`L(2,chi_p)` when `p=1 mod 4`, proves `(QVAR)` on every `k=1` and `k=3`
stratum for every prime.  Theorems L--O close `k=4`, and Theorem P makes
`k=5` empty for every `p>=41`.  Exact finite sieves close every remaining
`k=5` case as well.  Exact quartic-profile energies make `k=6` empty for every
`p>=47`; coefficient sieves close every remaining finite case, while its
`p=11` case already clears by complete census.  Thus the exceptional scalar
is closed through `k=6` for every prime and remains only on `k>=7` from
`p=13`; the principal delta-variance target is sharpened as above.  See the
updated reduction note and Prop 15.589.

The high-activity coefficient attack now has a general normal form.  Whenever
the top degree `k-2` is nonzero, translation uniquely kills the full
two-dimensional degree-`k-3` level.  Every lower degree `d` has
`k-d-1` free coefficients, with any `d+1` direction coefficients serving as
invertible pivots.  For `k=7` this gives depressed quintics
`a*s^5+c*s^3+d*s^2+e*s+f`.  Exact probes at `p=13,17,19` find minimum profile
energies `1,3,4`; energy alone does not eliminate the stratum, while all
`1,36,120` direction subsets have the predicted kernel ladder.  The next
exceptional computation should therefore be a recursive coupled coefficient
sieve, not a Cartesian product (`evidence/k7_quintic_profile_probe.py`).

The same degree theorem now gives exact arithmetic on every genuine profile:
`a_L in 2p Z`.  With `b_L=a_L/(2p)` and `T=(p^2-1)/8`, one has
`sum b_L=T` and `Z_psi/(2p)=T (mod 2)`.  Thus QVAR is the integer
anti-concentration target `E|sum psi(L)b_L|^2 >= 3T/8`.  Parity alone is far
too weak, but this normalization is binding and explains the p=11 histogram.

There is also an exact lattice-coset reformulation.  Max+ is the first shell
of the odd coset `y0+2 ker_Z(C-pI)`, not the ordinary lattice's first shell
(the latter has explicit norm-`p+1` Baer-line vectors).  The radius-sphere
benchmark for the exceptional quartic moment exceeds QVAR by
`q(q-1)(q-11)/(16(q+5))`.  Thus the live exceptional target is equivalently a
lower bound on one degree-4 odd-coset harmonic coefficient; proving that
coefficient nonnegative would suffice, but is still open.

The full p=11 directional covariance also kills a tempting shortcut: although
the quartic direction is top for the complete mixture, it is bottom among the
nonzero modes on k=4 and not top on k=6.  Stratum invariance alone cannot prove
QVAR by a top-eigenmode argument.

The latest attack kills a broader profile-only shortcut for every
`p=3 mod 4`, `p>=7`.  There are artificial full-support energy ensembles with
the exact conserved total, cyclic directional symmetry, equal means, integer
energies, separately admissible line profiles in every direction, and the
actual divisibility `a_L in 2p Z`, but quartic variance zero or `4p^2`,
below QVAR.  The fake profiles can also satisfy the individual polynomial
degree bound from Prop 15.588.
Therefore those facts cannot prove the bound even in combination.  Any
surviving profile argument must use the cross-direction coefficient kernels
and simultaneous Boolean ridge reconstruction, or an equivalent coupling
among directions (Prop 15.589 I).

The p=11 k=4 pure-parabola census further shows that the live inequality is
not fixed-active-subsetwise.  Each of nine genuine balanced four-direction
families has normalized moment `E B^2=5<45/8`; six unbalanced families have
`E B^2=63` and rescue the count-weighted aggregate to `39/2`.  Therefore the
proof must mix projective direction configurations even before mixing profile
strata (Prop 15.589 J).

The full-support coefficient descent kills another tempting induction.  At
`p=7` the top-degree-zero class is empty and every nonzero class clears QVAR.
At `p=11`, however, all `2,090,880` top-degree-zero vectors have actual profile
degree exactly three and moment `E B^2=137/36<45/8`; each of their twelve
projective leading-coefficient classes also fails.  The ten degree-four
nonzero classes each have the identical moment `111483/14039>45/8` and rescue
the exact mixture to `114771/14903`.  Thus QVAR cannot be proved separately by
actual profile degree or leading-coefficient class; adjacent degree families
must be mixed in their exact ensemble proportions (Prop 15.589 K).

There is also one positive high-prime closure.  A centered-Fourier bound for
nonconstant quadratic line profiles, plus six exact two-character-class
checks at `p=41,43,47,53,59,61`, proves that every active quadratic profile
uses more than one quarter of the conserved profile energy for every
`p>=41`.  Four active profiles are therefore impossible: the `k=4` stratum is
empty for all primes `p>=41` (Prop 15.589 L).

The same centered-Fourier argument plus Weil's additive-character estimate
works at arbitrary profile degree: a `k>=4` stratum is empty whenever
`p>4k^2`.  Hence every surviving stratum is `k=1`, `k=3`, or
`k>=sqrt(p)/2`; the first two already satisfy QVAR.  The unresolved
exceptional scalar is therefore asymptotically a high-activity problem
(Prop 15.589 M).

For `p=19,23,31`, the exact energy partitions and degree-2/degree-1
coefficient kernels leave zero constant-compatible candidates on all
`210,495,1820` direction subsets, respectively.  Together with the `p>=41`
barrier, this proves that when `p=3 mod 4`, `k=4` exists only at `p=7,11`;
its QVAR moment clears the target at both.  Thus the `k=4` contribution is
fully closed in this congruence class, and for `p>=19` QVAR starts at `k=5`
(Prop 15.589 N).

In the complementary congruence class, the same sieve proves `k=4` empty at
`p=29,37`; `p>=41` was already eliminated.  It regenerates the complete
`p=13,17` families and gives exact quartic moments `8788` and `314432/3`, both
above QVAR.  Since `p=5` has only three square directions, QVAR is now proved
on `k=4` for every prime.  The exceptional target starts at `k=5`
(Prop 15.589 O).

For `k=5`, a zero cubic-kernel scalar is impossible by the degree-at-most-two
energy bounds.  With nonzero scalar, all five profiles are cubic; translating
the input reduces exact energy enumeration to depressed cubics
`a s^3+c s+d`.  Their minima eliminate every prime `41<=p<101`, including the
sole numerical exception `p=43`, where all 28 relevant types have `b=45` and
sum to `225`, not the required `T=231`.  The general activity barrier handles
`p>=101`.  Hence `k=5` is empty for every `p>=41` (Prop 15.589 P;
`evidence/k5_cubic_energy_barrier.{py,json}`).

The finite cubic coefficient sieves close four more primes.  At `p=29`, all
736,828,092 low-energy type tuples give zero coefficient candidates; at
`p=37`, all 9,348 admissible leading patterns fail.  The `p=31` stratum is
nonempty but has only 8,000 translation representatives and exact moment
`E B^2=16704/5>45`.  The existing complete `p=11` census gives
`E B^2=163/9>45/8`.  Consequently `k=5` remains open only at
`p=13,17,19,23` (Prop 15.589 Q;
`evidence/k5_p{29,31,37}_coefficient_sieve.{py,json}`).

The same complete sieve closes those last four cases.  At `p=13,17`, direct
Gaussian-integer quartic evaluation gives moments `297468/31` and
`1650768/29`, above thresholds `10647/2` and `15606`.  At `p=19,23`, the
signed-energy moments are `29417/65>135/8` and `8908/19>99/4`.  Consequently
QVAR is proved on `k=5` for every prime (Prop 15.589 R;
`evidence/k5_p{13,17,19,23}_coefficient_sieve.{py,json}`).

For `k=6`, translating the input depresses every genuine quartic profile to
`a s^4+c s^2+d s+e`.  Exact minima at every prime `47<=p<=139` use more than
one sixth of the conserved energy; the general `p>4k^2` theorem handles
`p>=149`.  A vanishing quartic scalar is separately impossible by the cubic
and quadratic energy bounds.  Hence `k=6` is empty for every `p>=47`
(Prop 15.589 S; `evidence/k6_quartic_energy_probe.py` and
`evidence/k6_quartic_energy_probe_{low,high}.json`).

Exact coefficient-kernel sieves close three further `k=6` primes.  At
`p=37`, 8,189,942,400 raw type tuples leave no candidate on any of 27,132
direction subsets.  The same conclusion holds on all 54,264 subsets at
`p=41` and all 74,613 subsets at `p=43`.  Consequently the finite `k=6`
problem is reduced to `p=13,17,19,23,29,31` (Prop 15.589 T;
`evidence/k6_p{37,41,43}_coefficient_sieve.{py,json}`).

The accelerated exact elimination then closes all six residual cases at
once.  The `p=13,17,19` aggregate quartic moments are respectively
`8896212/955`, `149941632/2879`, and `10591740/103`, all above QVAR.  The
`p=23,29,31` coefficient sieves are empty after scanning 71.207, 20.937, and
2.971 trillion raw tuples.  Together with the complete `p=11` census and the
high-prime energy barrier, QVAR now holds on `k=6` for every prime
(Prop 15.589 U; `evidence/k6_p{13,17,19,23,29,31}_coefficient_sieve.json`).

**Settled.** Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) (`solution.md`). Paley \(\rho=1\) on \(n=p^2+1\) (`evidence/PROOF_rho_eq_1.md`).

**Not settled.** \(L=\lim\alpha_n=1/2\). That needs E(1): Paley \(C\) is a \(\Phi\)-minimizer for every prime \(p\ge5\). Soft-close from sandwich plus denseness, without E(1), is not a proof.

Live `e1_closed_general` is True only by the old wiring (affine residual (ii) plus two-level Type I). That is not E(1) for all \(p\ge5\). Aut-Schur, Gsum disj LB, and the cotangent pairing stay unused and False.

## 2026-08-21 — remaining estimates isolated; no flag flipped

The p=13 k=7 signed-PSL census is not a leftover close. Three Max+-free
estimates remain; none is proved. See
`evidence/NOTE_2026-08-21_remaining_general_p_estimates.md`.

| Leftover | Predicate | Remaining estimate |
|---|---|---|
| 1 | `phi_F_ge_6_proved_general` | QVAR on k≥7 (all p≥13) **and** principal `\|\|δ\|\|^2 ≤ n(n+10)^2/[6(n-14)(n-6)]`. Crude `E[s^4]≤2n^3` is too weak. |
| 2 | `residual_ii_k_eq_4p_empty` / `multilevel_ND_k_ge_4p_proved` | leftover+splus at k=4p. Walsh = U spans xor-hyperplane of affine_span(Max−) (15.598 cuts H; spanning open). Paley ES2 majorant is false. Leftover-only with min_+<2 exists and is not residual (ii). |
| 3 | `type_I_multilevel_bad_case_ND_closed` | `\|μ\|≤(p-2)/(2p^2)` on \|κ\|=1, equivalently `\|R̄₄\| ≤ \|L\|(p^4-1)+4(p-2)`. `\|μ\|≤\|T\|` does not close; `\|μ\|≤maj` is false at p=7. |

Live dump (twice): all three leftovers False; Gsum False; pairing False.
`e1_closed_general` is still True only by the old incomplete wiring.
Aut-Schur unused. L=1/2 OPEN. Do not import a leftover from p=13 orbits.

## 2026-08-20 session — read this first

**No flag flipped. `Max+` at `p=11` is now fully enumerated**, breaking the wall
`fable.md` names as the likely common blocker to leftovers 1 and 3 ("Max+ is
enumerable only for `p<=7`"). Full derivations, exact data, and eight corrections
to claims made mid-session: `evidence/TECHNICAL_NOTES_2026-08-20_maxplus_p11.md`.
Read `evidence/METHOD_why_500_props_never_moved_a_flag.md` before writing a new
numbered proposition — it names the loop that produced most of the ~500 that
never moved a flag, and gives four rules to avoid repeating it.

**Leftover 1, current form.**  The older variance-plus-multiplicity route below
has been refined by Prop 15.589.  Multiplicity is now proved exactly: every
principal block has degree `n`, and the only smaller block is one exceptional
degree-`n/2` scalar.  The live route therefore has two separate targets:

1. exceptional QVAR on `k>=7` from `p=13`,
   `E|Z_psi|^2 >= 3q(q-1)/16`; and
2. the principal sufficient room
   `||delta||^2 <= n(n+10)^2/[6(n-6)^2]`.

The identity `tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2` remains a live
way to attack target 2, but bounding `||M||_F^2` alone is no longer the single
floor target because the exceptional block must be handled by QVAR.  The
finite `p=5` floor is already checked.  See the PSL reduction note and the
older `NOTE_leftover1_variance_multiplicity_route.md` for the variance history.

**Leftover 3**: `mu/L` slack grows `0.769 -> 0.746 -> 0.307` (i.e. **growing**
headroom) from `p=5,7,11`. No structural change, just a third confirming point
with more room, not less.

**k=4 stratum of Max+ terminates at p=19** (confirmed both by GPU and by
independent uncapped CPU DFS on all three `p=19` subsets, with a `p=17` positive
control matching exactly — see technical notes §6). Not itself a leftover, but was
this session's original target before the enumeration wall turned out to be the
more useful thing to attack.

**Correction to the "Bottom multiplicity is exactly n" claim two sections below**:
false at `p=11`, where it is `244 = 2n`. The weaker `mult >= n` form survives and
is what the leftover-1 route above uses.

**Repo housekeeping**: GPU code that made the p=11 run possible
(`gpu_inner.py`'s auto-splitting flip resolution, memory-pool capping) had been
sitting **uncommitted on a tmpfs RAM disk** for a full session — a reboot would
have destroyed it. It is committed now. Branch `prop15586-maxplus-gram-reduction`
(`1fa0301`) is an ancestor of `main`; nothing from it was lost.

**Conflict worth knowing about**: the "Do not commit... 15.496, or 15.530" line
further down was written by a prior session and left unexplained beyond being
grouped with other dead ends. Both are now tracked on `main` per an explicit
live instruction in the 2026-08-20 session, with 15.530's two failing tests
quarantined as `xfail(strict=True)` rather than left red. If you are the one who
wrote that original instruction and had a reason beyond "these are dead ends",
that reason is not recorded anywhere — leaving this note so it is not silently
overridden twice.

## Open for E(1)

| Item | Flag | Status |
|---|---|---|
| \(\lambda_{\min}(\Phi)\ge6\) on \(Z\) | `phi_F_ge_6_proved_general=False` | Open. Exceptional QVAR is closed through `k=6` and remains on `k>=7` from `p=13`; principal blocks retain the delta-variance target. |
| Residual (ii), even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | Open. Walsh, the boundary ranges of 15.669, and finite p=11 size eight in 15.670 are closed; the first larger floor-plus-pair profiles and the full graph constraints remain. |
| Type I, Max− not two-level \(\{-1,-3\}\) | `type_I_multilevel_bad_case_ND_closed=False` | Open. Remainder is \(A_{\mathrm{full}}\). |
| Lemma D | True | Closed. Do not unflip. |

**Next attack.**  Complete the `p=13,k=7` orbit decomposition: accumulate the
depressed scalar-7 representatives from each new free orbit, ask CP-SAT for a
solution outside their union, and compute each orbit's exact `k=7` moment on
Jellyfin until infeasibility certifies exhaustion.  The first orbit contributes
1,071 scalar-7 representatives and its `k=7` slice clears QVAR; a second orbit
seed is already stored.  In parallel, prove QVAR generally on `k>=7`, or prove
the equivalent odd-coset degree-four harmonic excess is at
least `-q(q-1)(q-11)/(16(q+5))`.  Do not use a pointwise/orbitwise floor,
restricted-stratum PSD, ordinary minimum-shell design, or “quartic is top on
every stratum”: each is now disproved.  Positivity, conserved total, cyclic
symmetry, full support, and coarse divisibility are also insufficient; a
profile proof must now exploit the cross-direction coefficient kernels and
simultaneous Boolean realizability, with active-direction configurations mixed
before taking the second moment.  Quantization of exceptional projection norms
does not help by itself: the full p=11 census has 37 shells, no zero shell, and
minimum `4304/15 < 366=3n`, so even the nonzero-shell pointwise repair is
false.  The exact `p=13,k=7` witness `|Z_psi|^2=2548<10647/2` kills that
repair directly in the first unresolved stratum.  The coupled CP-SAT model is
now available for an orbit/weighted complete `k=7` second-moment computation;
more pointwise minimization cannot prove QVAR.  In parallel, an upper bound on
`||M||_F^2` may close the principal delta room.  Import `phi_F_ge_6` only when
both block types are controlled generally, never from finite-p data.

## Floor (leftover 1)

**Current block decomposition.** `Z=W_e direct-sum principal series`; QVAR is
the exact exceptional condition and the delta room is the current sufficient
principal condition.  The formulas below are the older equivalent Fourier
description and remain useful, but “name all of Q(r)” is no longer the binding
next step.

Wick: \(Q(\pm1)=8q^2\), off-diagonal \(4q^2\). \(\delta=4-Q/q^2\). Floor \(S_\square\ge6q^2\Leftrightarrow\langle\delta,\psi\rangle\le2\).

Live ensemble (not a general proof): \(Q_{++}/q^2=48/13\) at \(p=5\), \(1544/409\) at \(p=7\). \(Q=8A/D\) with \(D=\lvert H_+\rvert/(2p)=13,409\) (\(2^2+3^2\), \(3^2+20^2\); not a polynomial in \(p\)).

Named pieces:

- \(S(\lambda)=\mathrm{Kl}(1,\lambda^2/4)\) (15.550).
- \(F=-2(3p^2+2)\), \(Q_{3,02}=-4N(2p^2+1)/p\) (15.564).
- \(n_{1d}=m\binom{p}{m}\), \(n_{k=3}=\binom{m}{3}(p-1)q\), \(A_{1d}=-4p^3/(p-2)\).
- \(\mu_{k=3}=96p^4 P(r)/(p^2-1)\) (15.574). \(\mu_{1d}=2p^4(p^2-3p-2)/(p-2)\) (15.575).
- Exclusive 1D / \(k=3\) / full mix reconstructs live \(Q\) on every Paley×norm type at \(p=5,7\) (15.573). It is not a \(p\)-identity: the 1D+\(k=3\)-only mix is \(4.68>4\) at \(p=7\); \(\mu_{k=3}/q^2>4\) at \(p\ge11\).
- At \(p=5\), \(n_{\mathrm{full}}=0\), so \(Q_{++}/q^2=48/13<26/7\) (15.581; 15.507 \(J_{N^*}=2\)).
- 1D 4-point vanishes for \(r\notin\mathbb F_p\); \(p=5\) \(Q_{N^*}/q^2=32/13\) (15.582).
- Pointwise \(Q_y^{++}\le4q^2\) is false: about \(23\%\) of Max+ at \(p=5,7\) have \(Q_y/q^2>4\) (max \(5.33\), \(16\)). Any identity that uses only \(z_i^2=1\) and \(Cy=py\) cannot force the ensemble bound.
- Paley×norm types split into many \(\langle\mathrm{Frob},\mathrm{inv}\rangle\) orbits at \(p\ge11\) (++sub has \((p-3)/2\) orbits of size \(2\)). Two-type constancy is certified only for \(p\le7\).
- \(\mathrm{Gal}(\mathbb F_q/\mathbb F_p)\) acts on \(H_+\) with orbits of size \(1\) or \(2\). Orbit masses are \(1/\lvert H_+\rvert\) or \(2/\lvert H_+\rvert\), i.e. they name \(D\).
- \(n_{\mathrm{full}}=\lvert H_+\rvert-n_{1d}-n_{k=3}\) is \(0\) at \(p=5\) and \(90q\) at \(p=7\). \(\mu_{\mathrm{full}}\) is not a single formula in \(p\) (15.578).

## Residual (ii)

Official class is leftover Max− together with \(s_+\ge2\). leftover-only (\(s_+=0\)) at \(p=5\), \(k=20\) exists and is not this class.

- leftover+\(s_+\) empty for all \(n_F\) at \(p=5\), \(k=20\) (15.528).
- 15.585: leftover+\(s_+\) at \(k=4p\) forces \(\min_+=2\); \(\{2,4,6\}\) cannot have \(1_{S=2}\) a plus pair-slice.
- 15.598: square-direction \(\infty\cup L\) forces \(\sum_S y=0\) on Max−. Walsh ∀p is spanning of the xor-slice of H.
- 15.669: for `p>=17`, all-finite `6<=s<=3(p-1)/4` and infinity-present `5<=s<=p-4` are impossible; exact `p=11,13` extensions are listed at the top of this file. 15.670 additionally closes every finite `p=11` size-eight boundary. Larger count profiles survive only the current relaxation.
- No identity that leftover+\(s_+\) is empty at every even \(k\ge4p\).

## Type I

Two-level Max− is closed (15.272). Multi-level is open. Dead as a multi-level kill: Aut\(_e\) (15.559), Max± of \(C\) (15.565), Type+ 1D Johnson (15.577), Galois support plus \(F\) (15.580), \(\lvert\mu\rvert\le\lvert L\rvert\) on \(\lvert\kappa\rvert=1\) (unsigned \(\lvert\nu_{\mathrm{part}}\rvert\) exceeds \(\lvert L\rvert\)). Remainder is \(A_{\mathrm{full}}\).

## Do not reopen

Occupancy / Aut-involution pairing of \(T_{\mathrm{ns}}\) / \(\bar n_0\) interpolants / half-net census as a \(p\)-law / Aut\(_e\) as a name of \(A_{\mathrm{full}}\) / \((p-5)/15\) / \(10p-46\) / \(16(p-4)/D\) / Paley type as a \(Q\)-constant (false at \(p=7\)) / exclusive mix as a general \(Q_\tau\) / pointwise Wick or Boolean collision as a proof of \(Q_{++}\le4q^2\) / Gsum as a Gram / Aut-Schur.

15.495 catalogs, 15.496, and 15.530 **are now committed** (2026-08-20, see note
at top of file) — this line originally said not to. Left visible rather than
deleted so the reversal is traceable.

## Files

| File | Role |
|---|---|
| `STATUS.md` | Claim table |
| `GOAL.md` | Acceptance for E(1) / \(L=1/2\) |
| `solution.md` | Sandwich; Main Theorem (limit) stays OPEN |
| `evidence/share/denseness_path_package.md` | Stand-alone path; § Caveats |
| `src/e1_gmin_m4_prop15598.py` | Square-direction affine lines cut Max− (proved); Walsh spanning open |
| `src/e1_gmin_m4_prop15599.py` | Square-line F2-rank pin {n/2−1, n/2}; antipodes; Aut_e reducible |
| `src/e1_gmin_m4_prop15600.py` | rank(S)=n/2 all odd p (radical ⟨1⟩); Walsh spanning still open |
| `src/e1_gmin_m4_prop15601.py` | QR in rowspan(S) or S+ℓ (pencil); 15.406 E still OPEN |
| `src/e1_gmin_m4_prop15602.py` | G_aff^□ permutes rows of S; unique 1-dim invariant ⟨1⟩ |
| `src/e1_gmin_m4_prop15603.py` | H0 ∩ H0'=⟨1⟩; H0+H0'=even-weight; heart splits |
| `src/e1_gmin_m4_prop15604.py` | 1_QR ∈ H0 iff p≡1 (mod 4); ker(D−I)∩H0 dim 2 |
| `src/e1_gmin_m4_prop15605.py` | Paley A²=A over F2; H0=⟨1⟩⊕ translate-span of extra |
| `src/e1_gmin_m4_prop15606.py` | W=⊕ nsq W^H; M transits; irred if 2 primitive root mod p |
| `src/e1_gmin_m4_prop15607.py` | W G_aff-irred all odd p (F_p^× mixes Φ_p); dir(Max−)=H0 |
| `src/e1_gmin_m4_prop15627.py` | Octic linear-box empty; split-involution class W2-nonempty at p=31 |
| `src/e1_gmin_m4_prop15626.py` | Bounded a,b,i stay box empty on residual W1; W2 t=-2 at p=17 not p-law |
| `src/e1_gmin_m4_prop15608.py` | Two PSL-orbits of F_p-sublines; 1∈dir(U); I(H0) via 15.609 |
| `src/e1_gmin_m4_prop15609.py` | Opposite-type never tangent; I(H0)=H0 for every odd p |
| `src/e1_gmin_m4_prop15610.py` | Aut({0,∞}) uniqueness for Walsh DEAD; unipotent flag I-invariant |
| `evidence/SESSION_HANDOFF_2026-08-18_leftovers.md` | Named identities 15.550–15.585 |
| `evidence/TECHNICAL_NOTES_2026-08-20_maxplus_p11.md` | p=11 enumeration, derivations, corrections (2026-08-20) |
| `evidence/NOTE_2026-08-20_psl_and_stratum_floor_reduction.md` | Binding PSL decomposition, QVAR, low-stratum theorem, odd-coset route, and killed routes |
| `evidence/quartic_profile_attack.py` | Direct quartic/profile-energy diagnostic |
| `evidence/maxplus_p11/directional_energy_covariance_p11.{py,json}` | Full p=11 directional covariance and top-mode counterexample |
| `evidence/maxplus_p11/k4_active_subset_quartic_p11.{py,json}` | Exact p=11 k=4 active-subset split: balanced families fail QVAR, aggregate clears |
| `evidence/k4_p3mod4_coefficient_sieve.{py,json}` | Exact p=19,23,31 coefficient-level emptiness certificate |
| `evidence/k4_p1mod4_closure.{py,json}` | Exact p=13,17 quartic moments and p=29,37 emptiness; completes all-prime k=4 QVAR |
| `evidence/k6_coefficient_sieve_fast.py` | Exact Numba quadratic-three/linear-two elimination with orbit sharding |
| `evidence/merge_k6_coefficient_shards.py` | Validates and merges complete k=6 shard certificates |
| `evidence/k6_p{13,17,19,23,29,31}_coefficient_sieve.json` | Residual finite k=6 moments/emptiness; completes all-prime k=6 QVAR |
| `evidence/k7_quintic_profile_probe.py` | Exact depressed-quintic lift minima and universal seven-direction kernel audit |
| `evidence/k7_p13_cpsat_probe.py`, `k7_p13_cpsat_witness.json` | Exact coupled p=13,k=7 model and independently checked pointwise-QVAR counterexample |
| `evidence/k7_p13_signed_psl_orbit.py`, `k7_p13_signed_psl_orbit.json` | Packed-bit signed-PSL traversal: first witness has a free 4,826,640-vector signed orbit |
| `evidence/k7_p13_orbit_quartic_xpu.py`, `k7_p13_orbit_quartic_xpu.json` | Exact A380 quartic/activity census; first orbit and its k=7 slice both clear QVAR |
| `evidence/k7_p13_extract_orbit_representatives.py`, `k7_p13_orbit_completeness.py`, `k7_p13_second_orbit_seed.json` | Extracts 1,071 representatives per scalar and searches for the next orbit; a second seed was found and independently checked |
| `evidence/exceptional_projection_shell_probe.{py,json}` | Full p=5,7 exceptional shells and an explicit p=11 nonzero shell below `3n`; kills the nonzero-shell pointwise repair |
| `evidence/maxplus_p11/exceptional_projection_shells_p11_xpu.{py,json}` | Full 37.46M-row p=11 exceptional shell census: 37 shells, no zero, minimum `4304/15`, exact mean |
| `evidence/NOTE_leftover1_variance_multiplicity_route.md` | Historical principal variance route and `||M||_F^2` target |
| `evidence/METHOD_why_500_props_never_moved_a_flag.md` | Read before writing a new numbered proposition |
| `evidence/maxplus_p11/` | Scripts + logs for the p=11 spectrum/moment computations |
| `LONG_HORIZON_GOAL.md` | Terminal states |

Large `.npy` arrays (Max+ at p=11, 4.5 GB each) are **not in git** — they live at
`/mnt/storage/e1work/maxplus_p11/` on soulkiller, verified by md5 against the
original computation. Scripts there have hardcoded `/tmp/e1work` paths; repoint
before rerunning.

The p=13 orbit attack is centrally backed up at
`/mnt/storage/e1work/maxplus_p13/orbit_attack_2026-08-20/`.  It contains the
111 MB packed first and second orbits, the 1,071 packed scalar-7
representatives, both independent second-orbit CP-SAT certificates, both
orbit metadata records, both A380 censuses, and a SHA-256 manifest.  The
first-orbit packed hash is
`7223169420a18477dbdf95f6c3685186fbc6a7a1916ac875d761b22800c01eb2`.
The second is
`a3ce4e19e68770b41951b4ba28153fd5ed23884d0bcd912eeb43c421fa0e31c3`.

Jellyfin (`192.168.1.191`) now has a validated Intel Arc A380 environment at
`/home/nick/.venvs/mo-intel`: `torch 2.13.0+xpu`, `dpctl 0.22.1`, and
`pyopencl 2025.1`.  Torch XPU, Level Zero, and OpenCL all see the A380.  During
the shell census the soulkiller repo and p=11 data were mounted read-only under
`/home/nick/mnt/soulkiller-{repo,e1work}` via SSHFS.

**Checkpoint entering the latest continuation:** `54fd110` (exact `2p`
arithmetic, normalized integer QVAR, and the general coarse-profile
countermechanism).  The latest commit(s) after that checkpoint add the p=11
fixed-active-subset counterexample in Prop 15.589 J; use `git log -1` for HEAD.
