# Proposition and route de-duplication audit

**Date:** 2026-08-30

**Scope:** every assigned proposition through Proposition 15.719, the live
predicate wiring, and the current untracked attack scripts

**Purpose:** prevent a reformulation, solver/backend change, longer timeout,
or finite-prime rerun from being mistaken for a new mathematical attack

## Verdict

The duplication concern was correct.

1. The final 300-second positive-`p=7,z=7` CP-SAT run repeated an existing
   exact full-torsion model for the same case.  Only the timeout changed.
2. Several long proposition blocks are different coordinates for the same
   unresolved scalar or relaxation.  In particular, most of Props.
   15.83--15.160 are the optional Path-C/Hypothesis-H residual in different
   forms, and much of 15.321--15.560 is the same unnamed `Q_tau`/class-function
   mixture under successive small-prime fits.
3. Old `e1_closed_general()` wiring and several historical prose blocks still
   say `True`/`CLOSED` for a smaller, obsolete acceptance problem.  They do not
   close the current four-unit gate.
4. The shortened attack plan accidentally omitted **GLOBAL QVAR**.  This is
   not a new conjecture or wider scope: it has always been one conjunct of the
   original spectral-floor acceptance unit.

No new computation should be launched from an attractive formula or script
name until this file is checked first.

## Coverage and numbering

- Propositions 15.1--15.82 are written directly in `solution.md` and related
  early modules.
- There are 634 source-backed proposition modules from 15.83 through 15.719.
- The labels 15.537, 15.583, and 15.584 have no proposition module.  They are
  unassigned labels, not unreviewed propositions; later source headers mention
  those numbers only as historical range/state markers.
- Therefore every assigned proposition through 15.719 was included in this
  audit.  The grouped ledger below is by shared mathematical route rather than
  a 716-row restatement of the writeup.

## Authoritative acceptance chain

The public theorem is gated by `four_e1_units_closed()`, not by the legacy
`e1_closed_general()` value.

| unit | exact live content | status after 15.719 |
|---|---|---|
| spectral floor | **GLOBAL mixed-`k` QVAR** and principal R1 | **OPEN** |
| residual (ii) | non-Walsh multi-level Max-minus for every even `k>=4p` | **OPEN** |
| Type I | the multi-level `3A+B>0` bad case | **OPEN** |
| Lemma D | every good-line triple and its Fejer two-plane amplitudes | **TRUE** (15.276) |

The spectral-floor unit is a conjunction.  Strong R1
`||delta||^2 <= n/12` would also imply the weaker Type-I estimate, but it
does **not** prove GLOBAL QVAR and cannot close residual (ii).  Consequently
the shortest honest work map is:

1. spectral front: GLOBAL mixed-`k` QVAR plus strong R1;
2. residual front: the general non-Walsh multi-level remainder;
3. final implication audit.

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
| **15.561--15.589** | final class-function no-gos, exact profile classification, and the QVAR decomposition | 15.589 closes QVAR only for `k=1,...,6`.  Per-stratum `k>=7` is false and is **not** the leftover.  GLOBAL QVAR must retain the full mixed-`k`, unsplit ensemble; profile-by-profile, stratum-by-stratum, or leading-class proofs cannot establish it. |
| **15.590--15.628** | degree-four SoS countermechanism, exact R1/`delta` hierarchy, Walsh/W1/W2 investigation | Degree-four SoS cannot force Type I (15.590).  Character/PSD-only and fixed-channel R1 routes are insufficient.  **15.628 closes Walsh, W1, and W2 for all odd primes**, but explicitly leaves the non-Walsh 5+-level/even-`k>=4p` branch. |
| **15.629--15.668** | complete low R1 shells, modular-data no-go, nonlinear shell positivity, finite boundary closures, and exact `p=11` theta/channel work | Strong R1 is true at `p=11` by full census.  Scalar trace and broad square-circle conserved-mass cones through exponent 800 still admit sub-six targets and cannot prove general R1.  Props. 15.643--15.666 close the two-point, size-four, size-six, and finite `p=7` size-eight residual branches; rerunning their old solvers is duplication. |
| **15.669--15.712** | uniform residual ranges and exact endpoint closures | Props. 15.693--15.699 close the `p=19` endpoint; 15.700--15.712 close the `p=17` endpoint.  Remaining scope includes strict infinity-plus-`p`, the `p=23` next-boundary endpoint, later all-finite sizes, and the `p=7` infinity-plus-seven remainder.  Do not regenerate the closed `p=17/19` rows. |
| **15.713--15.719** | positive `p=7` infinity-plus-seven reductions | 15.713--15.717 close `z=0,1,2,3`.  15.718--15.719 identify and stabilize projected `z=7` semigroup supports but remove no source boundary.  All 56 actual `z=7` line boundaries remain open, and the semigroup/quotient route is terminated. |

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

A live R1 proof must use information absent from the abstract spectra: the
Boolean rank-one identity and the exact full Max-plus orbit mixture.  A live
QVAR proof must couple the mixed-`k` ensemble.

### Type I and residual

- Gsum disjoint lower bound, Aut-Schur, `k=3`-only span, or the cotangent
  pairing;
- affine/two-level residual work already closed by 15.179, 15.236--15.237,
  and 15.272;
- Walsh/W1/W2 search after 15.628;
- size-four, size-six, finite `p=7` size-eight, `p=17`, or `p=19` solver
  reruns after their exact closures;
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
| `scripts/p17_*` and old `p19` endpoint diagnostics | endpoints closed by 15.699 and 15.712 |
| `scripts/r1_p11_*`, broad theta/channel scripts | useful finite data, but the aggregate cones are certified insufficient and `p=11` itself is already R1-positive |
| `scripts/p7_infinity7_positive_z7_*` and matching `/tmp` artifacts | terminal semigroup/quotient campaign; not active |
| generic `residual_boundary_*` parity/projected models | necessary-condition diagnostics only; infeasibility can close a specified finite branch, feasibility cannot close or witness the graph problem |

The pre-existing modified files
`evidence/e1_gmin_m4_prop15626.json`,
`evidence/e1_gmin_m4_prop15627.json`,
`scripts/p19_second_boundary_profile_cryptominisat.py`, and
`scripts/residual_affine_johnson_milp.py` are user work and remain untouched.

## Stale or contradictory records

Central live documentation is corrected alongside this audit, but historical
artifacts are not rewritten en masse.  Treat the following as quarantined:

- `e1_closed_general()` can be `True` while all three open acceptance
  predicates are `False`; only `four_e1_units_closed()` is authoritative.
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
3. search tracked files, untracked scripts, `/tmp` artifact names and hashes,
   git history, GitHub, MathOverflow, literature notes, and OEIS when number
   patterns are involved;
4. write the acceptance condition before launching;
5. stop after the declared gate if the result is only `UNKNOWN`, a necessary
   survivor, a finite-prime fit, or another equivalent relaxation.

This preflight is the replacement for extending a finished line with one
more solver, projection, shell, channel, orbit, or timeout.
