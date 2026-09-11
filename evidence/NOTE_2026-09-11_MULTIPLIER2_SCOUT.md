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

Minimizers used: exact class-scan minimizers for `n <= 8`; for `n = 10, 12`
the family composites `[[A5,B5],[B5^T,-A5]]`, `[[A6,B6],[B6^T,-A6]]` with
`Phi = 13 = m_10`, `18 = m_12`; for `n = 14` the Paley `C14`
(`Phi = 21 = m_14`). `rho = 1` (midpoint); `rho = 0.5` is slightly worse.

**Reading.** `R_Z` stays in `[1.4, 1.7]` for `n >= 5` with no downward trend
over a factor-`7` range of `n`. If the ratio persists asymptotically, the
midpoint law alone will not meet the target (shortfall ~50% at measured
orders); the constructive completion route (§2) is the live path, and sharper
law variants (shifted thresholds) remain untested here.

## 2. Completion-family data

- `n = 2..6`: with exact minimizers and searched completions,
  `min_B Phi([[A,B],[B^T,-A]]) = m_2n` exactly (4, 5, 10, 13, 18; rigorous
  modulo the recorded exact values `m_10 = 13`, `m_12 = 18`).
- `n = 7`: the proved `C14` decomposition supplies a minimizer block with
  `Phi = 9 = m_7` and `Phi(C14) = 21 = m_14` — exists-A exactness, proved.
- `n = 8`: basin-hopping over 64 free signs finds `Phi = 32`, hence
  `m_16 <= 32` (target `28.28`, ratio `1.131`); search-based, witness in
  scratch.
- `n = 9`: the `C18` decomposition (block `12 = m_9`) gives
  `m_18 <= 33 <= 2 sqrt(2) m_9 = 33.94`.
- Weak explicit bounds from rounded midpoint samples: `m_20 <= 52`,
  `m_24 <= 68`, `m_28 <= 84` (min over 300 samples; each a genuine signing).

## 3. Repro

Session scratch certifiers (promotion pending):
`~/scratch/scout_v1.py`, `~/scratch/scout_v2.py`, `~/scratch/scout_v2b.py`;
merged outputs in `evidence/mo_multiplier2_scout_results_20260911.json`
(deterministic seeds `12345+n` / `777+n` / `100+s`; MC standard errors
<= 0.4%). `n = 14` chunking: `2^27` doubled states, 128 chunks.
