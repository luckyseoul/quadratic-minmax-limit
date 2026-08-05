# Min-max quadratic form of ±1 coefficients

MathOverflow [413935](https://mathoverflow.net/questions/413935) /
[X challenge](https://x.com/PI010101/status/2081070728422752329):

\[
m_n
=
\min_{a_{ij}=\pm1}
\max_{x_j=\pm1}
\Bigl|\sum_{1\le i<j\le n}a_{ij}x_i x_j\Bigr|,
\qquad
\alpha_n=\frac{m_n}{n^{3/2}}.
\]

## Status

**Main claim (2026-08-05):** \(\displaystyle L=\lim_n\alpha_n=\tfrac12\).

Path: sandwich + \(\rho=1\) denseness on Paley \(n=p^2+1\) + E(1) via
bi-tight empty (Prop 15.167) and freeness-fail ND (Props 15.170–15.171).
Residual / \(16N\) spectral package remains **optional open** (not required for \(L\)).

**Proved (sandwich):**
\[
\frac1\pi
\;\le\;
\liminf_n\alpha_n
\;\le\;
\limsup_n\alpha_n
\;\le\;
\tfrac12.
\]

**Also proved:** \(\rho=1\) for Paley conference matrices of order \(n=p^2+1\) (halfspace boolean eigenvector).

See `HANDOFF.md`, `solution.md` (Props 15.167–15.171), `evidence/share/paper/`, `x-cards/`.

## Files

| Path | Role |
|------|------|
| `HANDOFF.md` | Research handoff / resume entry point |
| `solution.md` | Full mathematical writeup |
| `src/e1_gmin_m4_prop15167.py` … `prop15171.py` | Bi-tight + E(1) residual ND modules |
| `src/minmax_quadratic.py` | Exact `m_n`, Paley, \(\Phi\), bounds, \(\rho=1\) evec |
| `tests/test_prop15167.py` … `test_prop15171.py` | Load-bearing E(1)/L tests |
| `x-cards/` | X summary + key-lemmas JPEGs |
| `evidence/share/` | Paper PDF/TeX + share assets |
| `evidence/` | Verification JSON and session notes |

## Quick check

```bash
python3 -m pytest tests/test_minmax.py -v
python3 -c "from src.minmax_quadratic import exact_m; print([exact_m(n) for n in range(2,9)])"
```

## Exact small values

| n | m_n | α_n (approx) |
|---|-----|--------------|
| 2 | 1 | 0.354 |
| 3 | 3 | 0.577 |
| 4 | 4 | 0.500 |
| 5 | 4 | 0.358 |
| 6 | 5 | 0.340 |
| 7 | 9 | 0.486 |
| 8 | 10 | 0.442 |
| 9 | 12 | 0.444 |
| 10 | 13 | 0.411 |

At \(n=10\), Paley (order \(p^2+1\), \(p=3\)) has \(\Phi=15>m_{10}\): conference is not exactly optimal.
Exact optima first appear at Hamming distance 5 from Paley, and the only 5-edge undercutters are 144 perfect matchings — see `evidence/N10_STRUCTURE.md`. Those 144 form one PΓL(2,9)-orbit (maximizer-drop criterion) — see `evidence/N10_MATCHING_CLASSIFY.md`.