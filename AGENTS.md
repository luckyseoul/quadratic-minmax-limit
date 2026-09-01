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
Do not reopen endpoint selection on either frame, an independent
skew budget, a finite pair census, the implication from three pairwise
diamonds to the tetrahedral diamond, or a statewise random-skew union bound.
Do not present the disk surrogate as an equivalent target: its asymptotic
form would prove a stronger
`1/sqrt(2pi)` lower bound, while only its zero-error form is disproved.
Do not launch another finite-prime, residue, orbit, or cell census as work on
the original question.  The long residual-(ii) section below is continuity
guidance only for a deliberately selected Paley route.

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

## Current residual-(ii) gate after Proposition 15.753

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

Residual (ii) is still open at critical `p=5,7`, at `p=11,k>=50`, in the
sole `p=13,k=60` residue `u=6` and all later p13 layers, and beyond
Proposition 15.752's band. Propositions 15.743--15.749 and 15.753 are finite certificates
or branch theorems,
not an all-prime row theorem.  The preferred generic front is therefore a
structural version of its common-energy/cut mechanism that survives when
the number of distance bins grows.  Another independent coefficient-cell
catalog, one-direction floor, halving heuristic, or longer complete-domain
timeout does not advance that gate.  The next finite p13 implication is to
attack `u=6`; do not reopen the closed `u=4` branch, couple the 336 survivors
to common graphs, launch a broad support-396 or mass-12 census, or reopen
the already-closed `u=0,3` rows.  Before launching a finite p11/p13-later
computation, identify the invariant that could extend beyond that one row or
explain why the finite row is a genuine base obstruction.

## Result discipline

Label every result as exactly one of: proved theorem, exhaustive finite
certificate, open reduction, counterexample, or retracted claim.  Never
promote computation or a heuristic pattern into a theorem.

After a genuine advance, update the proposition-dedup audit and the canonical
status/handoff documents in the same commit.  Preserve failed routes and
counterexamples: they are part of the project memory, not clutter to delete.
