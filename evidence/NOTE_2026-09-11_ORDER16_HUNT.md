# Order-16 local-search campaign: floor at 30 in the full space and in the m8-completion family

**Status:** draft session note (2026-09-11). Finite search evidence plus two new
verified witnesses. No all-orders claim; independent review pending.

## Method

New machinery, promoted to `scripts/hunt.py`. The scoring function
`Q_K(s) = sum_{i<j} K_ij s_i s_j` is **linear in K**: with per-pair basis
vectors `G_a(s) = s_i s_j` precomputed over all states, any move touching a
few entries evaluates exactly as `Q' = Q - 2 * sum_{a in move} K_a G_a` on the
current `Q` vector — no per-move re-enumeration. Iterated local search:

- kick: flip 3-9 random free pairs;
- exact descent over all single free-pair moves;
- exact 2-flip and 3-flip repair scans (all combinations, chunked);
- random restarts every 40 rounds; best-so-far track with entrywise re-scoring.

Order-16 deep round cost ~2.5-2.7 s on both GPUs. Two modes: **general** (all
120 pairs free) and **fixed-block** (pin `A8` top-left, `-A8` bottom-right;
search the 64-entry cross block — the completion family `[[A8,B],[B^T,-A8]]`).

Bug note: the first completion implementation scored only the cross block (the
pinned block's contribution to Q was missing), producing an impossible
`Phi = 20 < m_15 = 27`. Caught by that monotonicity contradiction, fixed, all
completion runs re-executed; the results below are from the fixed version, and
every witness was re-verified by an independent exact-enumeration kernel.

## Runs (order 16; ~44 GPU-min total across V100 + RX 9070 XT)

Seeds: `K16` tie intermediate (`Phi = 30`), `nuka16` source (`Phi = 30`), bank
`n=16` winner (`Phi = 32`); completion seeds include `B = +/-A8`.

```
general    V100  seed 101: 242 rounds, floor 30
general    9070XT seed 202: 244 rounds, floor 30
general    V100  seed 111: 241 rounds, floor 30
general    9070XT seed 212: 258 rounds, floor 30
completion V100  seed 313: 1511 rounds, floor 30 (reached t=0.6s)
completion 9070XT seed 414: 1616 rounds, floor 30 (reached t=2.2s)
```

About 4,112 deep ILS rounds; **no witness with `Phi <= 29` found in any
stream.**

## New witnesses (verified)

`evidence/m16_phi30_completion_witnesses_20260911.json` holds two distinct
completion-form `Phi = 30` objects (V100 seed 313, 9070 XT seed 414). Each
was verified by exact enumeration on CPU via an independent kernel
(`einsum`, not the search path), and each has top-left block exactly equal to
the order-8 minimizer `A8` (`Phi = 10`) and bottom-right equal to `-A8`.
Neither object equals the tie-machine `K16` intermediate (also re-verified at
`Phi = 30` this session). There are now **three distinct known order-16
objects with `Phi = 30`**, from two independent construction routes.

## Reading

- `m_16 in [27, 30]` unchanged. The doubling-instance target
  `m_16 <= 2 sqrt(2) m_8 = 28.28` was not met: no `Phi <= 28` witness exists in
  any found object, and ~4.1K rounds of deep ILS from the best known seeds did
  not descend below 30.
- The completion family at the order-8 minimizer floors at exactly 30 as well.
  Consistency of the record across an exact structured family and two
  independent free searches (plus the campaign's 32) supports 30 as a durable
  upper record even though the exact value remains in [27, 30].

## Repro

```
scripts/hunt.py --seed K16_phi30.json --seed seed16_nuka16src.json \
    --seed best_n16.json --minutes 11 --kick-seed 111
scripts/hunt.py --fixed-block m8_hunt.json --minutes 11 --kick-seed 313
```

Session drivers: `~/scratch/hunt_run16.sh`, `~/scratch/hunt_run16b.sh`;
seeds in `~/scratch/`.
