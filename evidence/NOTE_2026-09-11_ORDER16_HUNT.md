# Order-16 local-search campaign: corrected completion result, floors, and exact two-vertex scans

**Status:** draft session note (2026-09-11, **corrected same day**). The first
version of this note landed two completion-form `Phi = 30` witnesses; both were
invalid signings (off-diagonal zeros) and are **retracted** — see
§Corrections. No all-orders claim; independent review pending.

## Method

`scripts/hunt.py`: `Q_K(s) = sum_{i<j} K_ij s_i s_j` is **linear in K**; with
per-pair basis vectors `G_a(s) = s_i s_j` over all states, a k-flip move
evaluates exactly as `Q' = Q - 2 sum_{a in move} K_a G_a` on the current Q
vector. Iterated local search: kick 3-9 free pairs, exact 1-flip descent,
exact 2-/3-flip repair scans (chunked), random restarts, optional periodic
deep repairs (`--deep-every`, `--repairs-every`). Modes: general (all 120
pairs free) and fixed-block (pin `A8` and `-A8`; free block = the 64-entry
cross block B).

**Validity gate (added after the correction below):** every recorded witness
must have all off-diagonal entries equal to `+-1`; otherwise `record()`
rejects it, and seeds are validated before descent.

## Runs (order 16; ~44 GPU-min V100 + ~44 GPU-min RX 9070 XT)

```
general    V100  seed 101: 242 rounds, floor 30
general    9070XT seed 202: 244 rounds, floor 30
general    V100  seed 111: 241 rounds, floor 30
general    9070XT seed 212: 258 rounds, floor 30
completion V100  seed 823 (clean): 3727 rounds, floor 32 (reached t=3.7s)
completion 9070XT seed 824 (clean): 3083 rounds, floor 32 (reached t=6.4s)
completion smoke  seed 813 (clean): 179 rounds, floor 32
```

~4,150 deep ILS rounds after the retune at orders 17/18/19 additionally floor
at the bank records 32 / 33 / 39, and further runs at orders 20 / 21
(2026-09-12) floor at 42 / 44 (see session logs; no improvements). Every
order in 16–21 floors exactly at its record under this machinery.

## Corrections (same-day retraction of the completion claim)

- **Original claim (first version of this note):** two distinct completion-form
  `Phi = 30` witnesses (kick-seeds 313, 414), hence "the completion family
  floors at 30".
- **Defect:** both matrices contain **zero entries at off-diagonal positions**
  — the cross-block diagonal `(i, i+8)`, `i = 0..7` (and no others). They are
  not complete signings, so they witness nothing for `m_16`. The zeros entered
  through the fixed-block seeds built as `B = +/-A8` (A8 is a graph matrix with
  zero diagonal); zeros are **absorbing under sign flips**, so they survived
  every descent and repair. The smoke witness had the same pattern.
- **Detection:** the archived two-vertex ternary engine's `validate()` (all
  off-diagonal `|value| = 1`) refused both matrices as inputs; an independent
  zero-count check confirmed the pattern. The contaminated witness file is
  rewritten in place as a RETRACTED record
  (`evidence/m16_phi30_completion_witnesses_20260911.json`), matrices kept for
  audit.
- **Corrected result:** clean-seed runs with the validity gate give
  **completion-family floor 32** (`evidence/m16_completion_ils_clean_20260911.json`).
  The general-space results are unaffected: all general-mode witnesses pass the
  validity check (off-diagonals all `+-1`).

## Exact two-vertex scans (archived ternary engine)

`two_vertex_ternary_cuda.py` (verifier SHA `384cc06d...`), run on C14 (the
K15 top-left block; input SHA in the receipt):

- **C14 -> two vertices -> order 16: `phis = [21, 29, 32]`;** the exact
  complete two-vertex endpoint optimum is **32**
  (`evidence/two_vertex_exact_c14_20260911.json`, engine 0.011 s).
- **Structural deduction:** the order-14 core of the tie intermediate has
  `Phi = 27` (not 21): the order-15 source `S15` is a `Phi = 27` minimizer
  whose 14-core *already* has `Phi = 27`. Hence the tie intermediate is not in
  the C14 two-vertex family — consistent with, and explaining, the exact
  family optimum 32.
- The same engine's `validate()` is what caught the retracted matrices above.

## Consolidated order-16 picture (all values valid-signing-verified)

- Known `Phi = 30` objects (two distinct): the tie intermediate `K16_phi30`
  (`nuka15.tie.first`) and the `nuka16` regenerated source. Both pass the
  full sign checks (CPU + GPU).
- Family floors: C14 two-vertex (exact) 32; K15 one-vertex (exact) 32;
  A8-completion (clean ILS, ~7K rounds) 32; general space (ILS, ~1K rounds
  plus seeds) no descent below 30.
- `m_16 in [27, 30]`; the exact value is open. The 30-record objects resist
  all current descent machinery; no evidence of anything below 30 was found,
  and the two structured families tested exactly both floor at 32.

## Repro

```
scripts/hunt.py --seed K16_phi30.json --seed seed16_nuka16src.json \
    --seed best_n16.json --minutes 11 --kick-seed 111
scripts/hunt.py --fixed-block m8_hunt.json --minutes 13 --kick-seed 823
```

Session drivers: `~/scratch/hunt_run16.sh`, `hunt_run16b.sh`, `hunt_run16c.sh`,
`tvwork_scan.sh`; seeds and intermediate outputs in `~/scratch/`.
