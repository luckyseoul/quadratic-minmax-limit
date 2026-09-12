# Paley conference evaluations: C26, C30, C38 — and the explicit bound m_38 <= 109

**Status:** draft session note (2026-09-12). Exact complete-cube evaluations,
independently replayed on both GPUs (V100 CuPy, RX 9070 XT ROCm CuPy). Finite
bounds only; no asymptotic claim. Independent review pending.

## Method

Builder `paley_conference(q)` (prime-power capable), two-block reorder via
`char_split` (nonsquare multiplier; identity `A = C[T,T] = -C[Tc,Tc]`
entrywise and cross-block transpose verified), exact norm by chunked GPU
enumeration `Phi = max |sum((S K) S)/2|` over the `2^(n-1)` projective states
(float32; all values are integers, so the path is exact). Builder validation:
C14 -> 21 and C18 -> 33 reproduce the known values.

## Results

```
 q   order   Phi    context
 13   14     21     = m_14 (exact); block 9 = m_7
 17   18     33     = order-18 record; block 12 = m_9
 25   26     65     order-26 record is 61 (conference does not improve); block 30
 29   30     75     matches the recorded order-30 value; spectral 80.8
 37   38     109    NEW explicit bound m_38 <= 109; block 47; spectral 115.6
```

- **`m_38 <= 109`** (was: 115.58 general spectral conference bound; 121
  coherent-lift value). Verified twice — V100 (1243 s) and RX 9070 XT
  (1178 s) — identical value; conference Gram checks `C C^T = 37 I` pass on
  both. `alpha_38 <= 109/38^(3/2) = 0.4653` for this construction.
- Spectral comparisons: C26 attains its spectral bound exactly
  (65 = 26 sqrt(25)/2); C14, C18, C30, C38 are strictly below (21 < 25.2,
  33 < 37.1, 75 < 80.8, 109 < 115.6).
- The two-block structure holds for q = 25 and 37 (prime-power extension of
  the reviewed theorem). Block norms 9, 12, 30, 47 at block orders
  7, 9, 13, 19: the equality `block = m_k` seen at (7, 9) does **not** persist
  (30 > m_13 = 20; 47 > m_19 = 39).

## Completion-construction sweeps (negative)

About 1,800 candidates `B = sign(alpha A + beta G)` (grid plus random) for
`[[A, B],[B^T, -A]]` at sources n = 10, 11, 12, 14. Best-found / target
(2 sqrt(2) m_n) ratios: **1.36 (n=10), 1.14 (n=11), 1.18 (n=12), 1.25 (n=14)**
— no record improves at orders 20/22/24/28. The program's successful
completion constructions (C14, C18) are algebraic (conference-derived);
noise-generated cross blocks do not approach their quality. This corroborates
the earlier midpoint-law ceiling readings (R >= ~1.2, no decay).

## Repro

`~/scratch/paley_eval.py` (builder + kernel), `~/scratch/eval_c38.py`,
`~/scratch/comp_sweep.py`; logs `c38.log`, `c38_verify_nuka.log`, `sweep*.log`
in `~/scratch/`; evidence `evidence/paley_conference_evaluations_20260912.json`.
