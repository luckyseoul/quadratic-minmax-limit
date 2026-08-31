# Proposition and route de-duplication audit

**Date:** 2026-08-30

**Scope:** every assigned proposition through Proposition 15.727, the live
predicate wiring, and the current untracked attack scripts

**Purpose:** prevent a reformulation, solver/backend change, longer timeout,
or finite-prime rerun from being mistaken for a new mathematical attack

## Verdict

The duplication concern was correct.

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
After this correction, exactly two mathematical gates remain: multi-level
Type I and non-Walsh residual (ii).

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
The first endpoint not excluded here is `p=31,R=10`; it remains in the
disjoint rich-block normal form.  Larger slack and the top-level gates remain
open.

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
- There are 642 source-backed proposition modules from 15.83 through 15.727.
- The labels 15.537, 15.583, and 15.584 have no proposition module.  They are
  unassigned labels, not unreviewed propositions; later source headers mention
  those numbers only as historical range/state markers.
- Therefore every assigned proposition through 15.727 was included in this
  audit.  The grouped ledger below is by shared mathematical route rather than
  a 719-row restatement of the writeup.

## Authoritative acceptance chain

The public theorem is gated consistently by the corrected global
`e1_closed_general()` Boolean and `four_e1_units_closed()["closed"]`; both are
currently `False`. The historical bounded `True` is available only through
`e1_bounded_residual_split_closed()` and is not a global theorem predicate.

| unit | exact live content | status after audit of 15.727 |
|---|---|---|
| required bi-tight levels 2 and 3 | 15.720 degree congruence using 15.272/15.207 | **TRUE** |
| residual (ii) | non-Walsh multi-level Max-minus for every even `k>=4p` | **OPEN** — for `p>=17`, 15.721 moves the general boundary floor to `|D|=p+1`; 15.676 and 15.722--15.724 close pair equality and slack zero, while 15.726 excludes every positive outside slack through `floor((p-4)/3)`; 15.727 excludes equality at `p=17,19,23,29` and rigidifies the other endpoint cases, beginning at `p=31,R=10` |
| Type I | the multi-level `3A+B>0` bad case | **OPEN** — `|κ|=1` needs `G>T` (for example `|μ|≤|L|`), while `|κ|=3` independently needs `χ_d((2p-1)μ+(p-2)ν)>-(p-2)/p`; the particular term is safe but the δ remainder is open |
| Lemma D | every good-line triple and its Fejer two-plane amplitudes | **TRUE** (15.276) |

The spectral floor remains an interesting optional problem, but it has no
valid downstream role in the current E(1) proof. The shortest honest work map
is now the two multi-level remainders followed by the final implication audit.

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
| **15.726** | minimal arc deletion plus the Ball--Lavrauw dual tangent envelope | **PROVED narrowing, not shell closure.** For every prime `p>=17`, it excludes `1<=R<=floor((p-4)/3)` at `|D|=p+1`; any positive survivor must have `R>=floor((p-1)/3)`. Residual (ii), Type I, and `L` remain open. |
| **15.727** | endpoint tangent-envelope equality, disjoint rich blocks, and published arc classifications | **PROVED narrowing, not shell closure.** Equality forces `R` index-one points outside an arc and disjoint trisecant/4-secant blocks. This excludes the endpoint at `p=17,19,23,29`; the first unexcluded endpoint is `p=31,R=10`. Residual (ii), Type I, and `L` remain open. |

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
- any `|D|=p+1` full-circle or outside-slack campaign with
  `0<=R<=floor((p-4)/3)`: 15.722--15.724 close `R=0`, and 15.726 closes the
  stated positive interval exactly; 15.727 also closes the endpoint at
  `p=17,19,23,29`, so do not rerun those finite endpoint classifications;
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
| `scripts/p17_*`, old `p19`, and `p23` all-finite endpoint diagnostics | superseded as boundary gates by 15.721 signed transport; retain only as historical certificates or reusable finite-geometry work |
| `scripts/r1_p11_*`, broad theta/channel scripts | useful finite data, but the aggregate cones are certified insufficient and `p=11` itself is already R1-positive |
| `scripts/p7_infinity7_positive_z7_*` and matching `/tmp` artifacts | terminal semigroup/quotient campaign; not active |
| generic `residual_boundary_*` parity/projected models | necessary-condition diagnostics only; infeasibility can close a specified finite branch, feasibility cannot close or witness the graph problem |
| `scripts/infinity_plus_p_quantized_dp.py` | diagnostic only until every `excess != 2` use retains the two Proposition 15.723 `p=17` equality cells |

The pre-existing modified files
`evidence/e1_gmin_m4_prop15626.json`,
`evidence/e1_gmin_m4_prop15627.json`,
`scripts/p19_second_boundary_profile_cryptominisat.py`, and
`scripts/residual_affine_johnson_milp.py` are user work and remain untouched.

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
3. for a nonempty residual boundary, first use signed PSL transport to put a
   boundary point at infinity; do not open an all-finite profile campaign;
4. at `|D|=p+1`, skip slack zero and every positive outside slack through
   `floor((p-4)/3)`; 15.722--15.724 close zero and 15.726 closes that positive
   interval.  At equality, 15.727 closes `p=17,19,23,29`; for the remaining
   primes start from its disjoint trisecant/4-secant normal form.  The first
   unexcluded endpoint is `p=31,R=10`.
   Retain
   both 15.723 floor-plus-two exceptions in any profile DP;
5. search tracked files, untracked scripts, `/tmp` artifact names and hashes,
   git history, GitHub, MathOverflow, literature notes, and OEIS when number
   patterns are involved;
6. write the acceptance condition before launching;
7. stop after the declared gate if the result is only `UNKNOWN`, a necessary
   survivor, a finite-prime fit, or another equivalent relaxation.

This preflight is the replacement for extending a finished line with one
more solver, projection, shell, channel, orbit, or timeout.
