# Multiplier-2 route: midpoint-law readings and completion-family data

**Status:** draft session note (2026-09-11, Codewhale). Finite Monte Carlo and
exact enumerations; scratch certifiers (paths in §3); independent review
pending. No all-orders claim. Complements
`NOTE_2026-09-11_PALEY_TWO_BLOCK_DECOMPOSITION.md`.

## 1. Midpoint-law (T2) readings

Object (`NOTE_2026-09-05_UNIVERSAL_SPECTRAL_MIDPOINT_GAUSSIAN_REDUCTION`,
eq. (26)): `E Phi([[A_n, Z],[Z^T, -A_n]])` with
`Z = mat(sqrt(kappa) G + sqrt(1-kappa) W)`, `G ~ N(0, Sigma*)`, `A_n` an exact
order-`n` minimizer; target `2 sqrt(2) * m_n`. Errors are MC standard errors
(<= 0.02 in `R`). `Z` row: the smoothed law; `B` row: `mat(sign G)`.

n  | m_n | target | EPhi_Z | R_Z   | EPhi_B | R_B
---|-----|--------|--------|-------|--------|-------
2  | 1   | 2.83   | 3.83   | 1.354 | 4.00   | 1.414
3  | 3   | 8.49   | 7.70   | 0.908 | 7.20   | 0.849
4  | 4   | 11.31  | 13.22  | 1.168 | 13.51  | 1.194
5  | 4   | 11.31  | 18.32  | 1.619 | 18.76  | 1.658
6  | 5   | 14.14  | 23.81  | 1.684 | 24.67  | 1.744
7  | 9   | 25.46  | 32.38  | 1.272 | 32.93  | 1.294
8  | 10  | 28.28  | 39.92  | 1.411 | 40.46  | 1.430
10 | 13  | 36.77  | 56.86  | 1.546 | 57.53  | 1.565
12 | 18  | 50.91  | 73.39  | 1.442 | 74.24  | 1.458
14 | 21  | 59.40  | 90.22  | 1.519 | 91.65  | 1.543
15 | 27  | 76.37  | 105.89 | 1.387 | —      | —
16 | 32* | 90.51  | 116.26 | 1.284 | —      | —
17 | 32* | 90.51  | 122.91 | 1.358 | —      | —

Rows `n >= 15` computed on the V100 (CuPy 14.2, the campaign's
`Q = sum((S @ K) * S, axis=1)/2` scoring pattern); the GPU path reproduces
the `n = 10` CPU reading per-sample to `1e-6`. Minimizers used: exact
class-scan minimizers for `n <= 8`; for `n = 10, 12` the family composites
`[[A5,B5],[B5^T,-A5]]`, `[[A6,B6],[B6^T,-A6]]` with `Phi = 13 = m_10`,
`18 = m_12`; `n = 14` Paley `C14` (`Phi = 21 = m_14`); `n = 15` the exact
minimizer `K15 = [[C14, 1],[1^T, 0]]` — the all-ones one-vertex extension,
`Phi = 27 = m_15` (rebuilt and verified this session,
`evidence/K15_exact_minimizer_20260911.json`); `n = 16, 17` campaign winner
matrices (`Phi = 32` each; `m_16 in [27, 32]` and `m_17 <= 32` unknown —
marked `*`: targets/R for those rows are stated against `Phi(A)`).
`rho = 1` (midpoint); `rho = 0.5` is slightly worse.

**Reading.** `R_Z` stays in `[1.4, 1.7]` for `n >= 5` with no downward trend
over a factor-`7` range of `n`. If the ratio persists asymptotically, the
midpoint law alone will not meet the target (shortfall ~50% at measured
orders); the shifted-threshold sweep (§1.1) shows this ceiling is intrinsic to
the one-sided Gaussian family, and the constructive completion route (§2) is
the live path.

### 1.1 Shifted-threshold sweep (same session)

The same proved reduction admits a deterministic threshold `h`:
`Z_h = s_h A + 2 phi(h) G + sqrt(1 - s_h^2 - 4 phi(h)^2) W`,
`s_h = erf(h/sqrt(2))`; at `h = 0` this is exactly the midpoint law. Sweeping
`h in [-2, 4]` plus the `h -> inf` limit (`Z = A`), with the same minimizers:

n  | h=0 (R_Z) | best finite h | best R | h=inf, Z=A (R)
---|-----------|---------------|--------|----------------
7  | 1.268     | +4.00         | 1.180  | 1.179
8  | 1.411     | +1.00         | 1.369  | 1.414
10 | 1.532     | +4.00         | 1.362  | 1.360

Tuning `h` recovers part of the overshoot (`n = 10`: `1.53 -> 1.36`) but the
family plateaus at `>= ~1.2` across the tested orders; no threshold reaches the
target. The `h -> inf` limit is the deterministic object `H (x) A`
(`H = [[1,1],[1,-1]]`) — **not a signing**: its `n` matching entries `(i, n+i)`
are zero. Charging those zeros at `+n` gives only the weak bounds
`m_20 <= 60`, `m_16 <= 48`, `m_14 <= 37`, all dominated by §2; do not infer a
`2 sqrt(2) Phi(A) + n` repair from this limit
(`NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md` disproves any
`o(n^{3/2})` repair on the whole class `Phi(A) = O(n^{3/2})`). The limit value
is within ~1% of the best finite-`h` value — the ceiling is intrinsic, not a
tuning artifact.

## 2. Completion-family data

- `n = 2..6`: with exact minimizers and searched completions,
  `min_B Phi([[A,B],[B^T,-A]]) = m_2n` exactly (4, 5, 10, 13, 18; rigorous
  modulo the recorded exact values `m_10 = 13`, `m_12 = 18`).
- `n = 7`: the proved `C14` decomposition supplies a minimizer block with
  `Phi = 9 = m_7` and `Phi(C14) = 21 = m_14` — exists-A exactness, proved.
- `n = 8`: basin-hopping over 64 free signs finds `Phi = 32`, hence
  `m_16 <= 32` (target `28.28`, ratio `1.131`); search-based, witness in
  scratch. Extended runs (structured starts: Sylvester `H16 - I` = 40.0,
  locally improved to 36; random starts 46-50) do not improve 32.
  **Update (same session):** the archived growth-checkpoint tie intermediate
  gives `m_16 <= 30` (verified CPU+GPU; `NOTE_2026-09-11_GROWTH_CHAIN_M16.md`),
  so with `m_15 = 27`, `m_16` lies in `[27, 30]`.
- `n = 9`: the `C18` decomposition (block `12 = m_9`) gives
  `m_18 <= 33 <= 2 sqrt(2) m_9 = 33.94`.
- Weak explicit bounds from rounded midpoint samples: `m_20 <= 52`,
  `m_24 <= 68`, `m_28 <= 84` (min over 300 samples; each a genuine signing).

## 3. Repro

Session scratch certifiers (promotion pending):
`~/scratch/scout_v1.py`, `~/scratch/scout_v2.py`, `~/scratch/scout_v2b.py`,
`~/scratch/scout_v3.py` (threshold sweep), `~/scratch/scout_v3b.py` (order-16
search), `~/scratch/scout_v5.py` (GPU, run with
`~/.venvs/mo-exact/bin/python`) and `~/scratch/prep_k15.py`; merged outputs in
`evidence/mo_multiplier2_scout_results_20260911.json` (deterministic seeds
`12345+n` / `777+n` / `100+s` / `200-403` / `777+n` (v5); MC standard errors
<= 0.9%). `n = 14` chunking: `2^27` doubled states, 128 chunks; v5 chunking
`2^21` states per chunk.
