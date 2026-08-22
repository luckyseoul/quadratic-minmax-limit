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

**Goal:** settle the limit (see **`LONG_HORIZON_GOAL.md`**). Not done until \(L\) is proved or disproved.

**Main claim:** \(\displaystyle L=\lim_n\alpha_n\) is **OPEN** (2026-08-16).

Sandwich and Paley \(\rho=1\) are proved. E(1) on \(n=p^2+1\) is **not**.
Four leftovers (`GOAL.md`): \(\lambda_{\min}(\Phi)\ge6\); residual (ii) for
even \(k\ge4p\); Type I when Max− is multi-level; Lemma D (writeup exists,
hostile check still due). Residual (ii) is closed only for the affine branch
and even \(k\le4p-2\) (15.179/236/237), not for the statement E(1) needs.
Soft-close forbidden. Package: **`evidence/share/denseness_path_package.md`**.

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

**Also proved:** \(\rho=1\) for Paley conference matrices of order \(n=p^2+1\).

See **`STATUS.md`**, `HANDOFF.md`, denseness package, `solution.md`.

## Files

| Path | Role |
|------|------|
| `HANDOFF.md` | Research handoff / resume entry point |
| `evidence/HISTORY_AND_REFERENCES.md` | MO/X/Paata education and pre-internet sources (not a close) |
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