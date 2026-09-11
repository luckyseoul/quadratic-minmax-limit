# Exact one-vertex extension scans from the winner bank: m_19 <= 39

**Status:** draft session note (2026-09-11). Exact computations; the improving
witness is verified by full enumeration on CPU and replayed on both GPUs
(V100 CuPy, RX 9070 XT ROCm CuPy). Independent review pending. No all-orders
claim.

## Method

For a signing `A` of order `n`, the optimal one-vertex extension satisfies the
exact identity

    min_b Phi([[A, b],[b^T, 0]]) = min_b max_x ( |Q_A(x)| + |b.x| ),

with `b` ranging over the `2^(n-1)` representatives (`b_0 = +1`; global flip is
a gauge). The scan enumerates all representatives against all `2^(n-1)` states
`x` with a chunked GPU matmul; the minimizing witness is re-verified by full
enumeration of the extended signing. (Identity: separation of the two new
spins, `PROOF.md` eq. (1) in
`evidence/original_mo_growth_checkpoint_20260906/`.)

## Results (exact optimal extension values)

```
source order / Phi(A)  ->  extended order : value   (bank value at that order)
15   K15       27         16 : 32                   (bank 30)
16   K16 tie   30         17 : 34                   (bank 32)
17   bank      32         18 : 33                   (bank 33)  matches
18   bank      33         19 : 39                   (bank 41)  IMPROVES
19   new(39)   39         20 : 46                   (bank 42)
19   bank      41         20 : 42                   (bank 42)  matches
20   bank      42         21 : 44                   (bank 44)  matches
```

## The improvement

`m_19 <= 39`. Witness: the **all-ones one-vertex extension** of the order-18
winner (`Phi = 33`): `K19 = [[A18, 1],[1^T, 0]]`; `Phi(K19) = 39` verified by
exact enumeration (`2^18` states) on CPU and replayed on the V100 and
RX 9070 XT. Matrix: `evidence/m19_phi39_witness_20260911.json`. (The all-ones
extension pattern also produced the exact `m_15` minimizer,
`K15 = C14 + all-ones vertex`; see `NOTE_2026-09-11_GROWTH_CHAIN_M16.md`.)

## Scope

Finite bounds only: `m_19 <= 39` (was `<= 41`). No convergence claim; lower
bounds are unchanged.

## Repro

`~/scratch/ext_scan.py --a ORDER|file:PATH --bank NEXT` (session scratch);
witness JSONs `~/scratch/ext_scan_{15..20}.json`.

## Two-vertex scans (checkpoint ternary engine, V100)

The archived reviewed engine (`two_vertex_ternary_cuda.py`,
`original_mo_path_v100_20260906_Q86jOz/`, verifier SHA `384cc06d…`) computes
the exact complete two-vertex endpoint optimum in `O(3^n)` on the GPU. Runs
(GPU time 0.03–0.30 s each):

```
source / [source, best intermediate, endpoint] : next order : comparison
K16 (30)   [30, 36, 37] -> 18 : 37   (bank 33)
17 bank    [32, 33, 39] -> 19 : 39   (the one-vertex 39, corroborated)
18 bank    [33, 43, 46] -> 20 : 46   (bank 42)
19 bank    [41, 44, 44] -> 21 : 44   (bank 44)  matches
K19 new    [39, 48, 48] -> 21 : 48
```

No two-vertex growth improves the bank bounds at 20/21. The `m_19 <= 39`
record is reproduced independently from the order-17 winner via the
two-vertex route. Every run carries the engine's built-in integer field
verifier receipts (`result.json` in `~/scratch/tvwork/out_*/`).

Follow-up scans: the two-vertex order-19 endpoint (`second.txt`, a second
`Phi=39` witness) extends one-vertex to exactly 46 at order 20 (all-ones);
the order-18 intermediate (`first.txt`, `Phi=33`) extends two-vertex with
`[33, 43, 46]` -> 46 at order 20. Every current best source saturates at 46
for order 20: the bank's 42 (n=20) and 44 (n=21) are not reachable by this
growth family, which now also rules out the natural extension attack for
20/21.
