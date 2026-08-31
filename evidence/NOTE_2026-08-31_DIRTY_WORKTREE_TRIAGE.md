# Dirty-worktree triage

**Date:** 2026-08-31

**Base:** `81e432d` (through Proposition 15.727)

**Status:** repository hygiene only; no theorem gate changes

The post-audit working tree contained 88 modified tracked files and 261
untracked files (9,406,009 bytes).  The useful corrections are retained on
`main`; the remaining snapshot is preserved on the remote branch
`archive/dirty-worktree-2026-08-31`.  Nothing in that branch is an active
acceptance gate.

## Retained on `main`

- The legacy gate tests now require the current honest state: full residual
  (ii), `E(1)`, and `L` are open.  They no longer tolerate or assert obsolete
  bounded closes.
- `NOTE_2026-08-21_leftover3_contraction_closure.md` now says explicitly that
  the recorded `mu` bounds address only the `|kappa|=1` Type-I half; the
  independent signed `|kappa|=3` inequality remains open.
- The deterministic 15.241--15.248 JSON summaries are synchronized with the
  proved residual-I dual-equality predicate, and the 15.429 summary is
  synchronized with its source theorem text.
- `residual_affine_johnson_milp.py` now loads SciPy only on its HiGHS/sparse
  paths, so its CP-SAT and SCIP modes do not require that unrelated optional
  dependency at import time.
- Six `p=7,z=7` modules are retained solely as historical reproduction
  dependencies.  Tracked 15.718/15.719 certificate generators import them
  directly or transitively:
  `p7_infinity7_global_modular_cpsat.py`,
  `p7_infinity7_positive_z7_mod7_projection.py`,
  `p7_infinity7_positive_z7_pointed_mod7.py`,
  `p7_infinity7_positive_z7_pointed_full_cpsat.py`,
  `p7_infinity7_positive_z7_pointed_compact_cpsat.py`, and
  `p7_infinity7_positive_z7_survivor_compact_batch.py`.
  Their presence does not reactivate the terminated semigroup/quotient route.
- `w2_pgl2q.py` is retained because Proposition 15.621 cites its exhaustive
  `86400/0` PGL scan; dead overwritten determinant scratch was removed before
  tracking it.
- Five GP scripts are retained as the missing upstream reproduction chain for
  Proposition 15.667's 45-dimensional half-gap cache:
  `r1_p11_scalar_kohnen_cache.gp`,
  `r1_p11_scalar_support400_cache.gp`,
  `r1_p11_scalar_half_gap_shard.gp`,
  `r1_p11_scalar_half_gap_assemble.gp`, and
  `r1_p11_scalar_half_gap_reduce.gp`.  Their machine-specific `/home/nick`
  paths were replaced by required environment inputs.  R1 remains optional
  and non-load-bearing after 15.720.

## Retired from `main`

- ~~Partial JSON regenerations made without their external `/tmp` catalogs~~:
  these dropped valid payloads and are not replacement certificates.
- ~~Elapsed-time, floating-point, newline-only, and pretty-print churn~~:
  these contain no mathematical change.
- ~~The generalized `p19_second_boundary_profile_cryptominisat.py` worktree
  version~~: it mixes `p=17,19,23` with p=19-specific core sizes, slot bounds,
  and inconsistent string/integer profile keys.  It is an exploratory draft,
  not a sound public interface.
- ~~Finite p=5, p=7, p=17, and p=19 backend/profile reruns~~: their targeted
  branches are already closed or superseded by signed transport and 15.727.
- ~~The `p=7,z=7` projection dump and duplicate orbit directory~~: one is a
  necessary-only survivor set, not a close; the other duplicates a tracked,
  newer completed certificate.
- ~~W1/W2/Walsh, R1/QVAR, circle/fable, and residual fixed-boundary scratch~~:
  these are killed, optional, necessary-only, or superseded routes already
  classified in the proposition de-duplication audit.

The archive branch is for forensic recovery only.  Do not merge it into
`main`, cite an archived solver status as a theorem, or extend one of its
campaigns without satisfying the repository's changed-premise gate.

## Live frontier after cleanup

Proposition 15.727 remains the last proved strict advance.  At
`|D|=p+1`, the first unexcluded endpoint is `p=31,R=10`; for general primes
the endpoint has the disjoint trisecant/4-secant normal form.  Larger slack,
the non-Walsh residual-(ii) remainder, the `|kappa|=3` multi-level Type-I
remainder, and therefore the quadratic-minmax limit remain open.
