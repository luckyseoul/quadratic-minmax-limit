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
