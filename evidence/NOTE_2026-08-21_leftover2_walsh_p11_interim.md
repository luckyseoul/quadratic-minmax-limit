# Leftover 2: Walsh mechanism (15.406 Theorem C) — first evidence at p=11

> **RESOLVED — full-ensemble result landed.**  `closed = True`, EXACT,
> over all 37,457,112 points (not a sample).  See Prop 15.596
> (`src/e1_gmin_m4_prop15596.py`).  This is still only a FOURTH census
> point (p=3,5,7,11), not a general-p proof; leftover 2 stays False.

Date: 2026-08-21.  Target: `residual_ii_k_eq_4p_empty` / `multilevel_ND_k_ge_4p_proved`
(leftover 2), currently False.  Route: 15.406 Theorem E — the Walsh
containment `W_{U^c} ⊆ W_U` is certified only at p=3,5,7; general p≥11 is
open.  **No flag flipped. This note documents an in-progress
extension, not a completed result.**

## Correction to an initial misreading

My first attempt tested the wrong object: score-slices of random size-4p
edge subgraphs. That is unrelated to 15.406's actual construction, which
splits the *entire* Max− ensemble by the sign of ONE fixed edge
`f_e = C_{ij} y_i y_j` (i,j fixed), independent of any leftover-witness
graph. The random-G test is discarded; it says nothing about Theorem C.

## The real mechanism, and what it reduces to

`walsh_U_implies_Uc(p)` (already in `src/e1_gmin_m4_prop15406.py`) computes,
in 0/1 coordinates `x = (1−y)/2`: the GF(2) nullspace `N₀` of `B_U` (rows =
`x` for `y∈U`) plus a particular solution `x₁` with `B_U x₁ = 1`, then
checks every element of `N₀ ∪ (x₁+N₀)` — i.e. every Walsh character
constant on U — is *also* constant on U^c. At p=5,7 this holds
(`closed=True`) and `rank(B_U)` there equals `n/2` exactly, the same
dimension as the full GF(2) affine span of Max− itself (verified
independently: `dim W(Max±) = n/2` at both primes, matching the ternary
p-eigenvector lattice bound).

## p=11: sample result (fast, suggestive, NOT yet rigorous)

Using the stored `maxplus_p11_eps1.npy` (37.4M rows) with the 15.254 swap
transport to Max−, a random 6000-row sample gives, in 8 seconds:

    |U|=3232 |Uc|=2768  rank(B_U)=60  (n/2=61)
    ker_mixed=0  aff_solvable=True  aff_mixed=0
    CLOSED = True

This is the first p=11 data point ever computed for this mechanism, and
it's a clean pass — but it is a **sample** of 6000 out of 37.4M points,
not the complete ensemble. p=5 and p=7's certified results use the full
cached arrays (260 and 11452 points respectively, small enough to hold
entirely); nothing analogous exists yet at p=11.

## Full-ensemble verification: DONE (Prop 15.596)

First implementation attempt had a performance bug (a ~2.3 billion-call
pure-Python `bin().count('1')` loop in the consistency-check pass) that
would have taken hours; killed before completion. Rewritten with
`np.bitwise_count` (numpy ≥2.0) for a fully vectorized pass 2, packing
each 122-bit row as a `(lo:uint64, hi:uint64)` pair
(`scripts/walsh_theorem_c_p11_full.py`). Three streaming passes over all
37,457,112 points (740s total, no subsampling):

    rank(B_U) = 60   (stable from the first 500K rows onward)
    ker_dim   = 62   solvable = True
    ker_mixed = 0    aff_mixed = 0
    CLOSED = True

**Note the mechanism differs from p=5,7**: there rank(B_U) = n/2 exactly
(B_U alone spans the whole Max− flat, making the containment trivial).
At p=11, rank(B_U) = 60 < n/2 = 61 — Theorem C holds through the genuine
algebraic containment N₀∪(x₁+N₀) ⊆ constant-on-U^c, NOT because B_U has
full rank. Any future general-p proof attempted via "B_U has full rank"
would be false at p=11; do not use that route.

## What this would mean if it lands True

Extending Theorem C's containment to p=11 would be the first evidence
beyond the certified p≤7 range for leftover 2's Walsh route — a genuine
census point, not yet a general-p proof (per fable.md's acceptance bar,
a single further prime does not close the leftover; the mechanism would
still need to be proved for all p, likely via the same signed-orbit /
character-sum toolkit used for leftover 3 today, applied to the fixed-edge
split of Max− rather than to four-point moments).

## Scope reminder

This is entirely independent of `NOTE_2026-08-21_nu_convolution_reduction.md`'s
δ-hierarchy (props 15.590–15.595): 15.595 proved leftover 2 does *not*
reduce to the ‖δ‖² bound that closes leftovers 1 and 3. Leftover 2 remains
its own separate root of the program (R2).
