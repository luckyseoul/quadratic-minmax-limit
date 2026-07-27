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

**Open:** existence of \(\lim\alpha_n\) (whether \(\liminf=\limsup\)).

Soft multipartite/Hadamard comparison cannot force existence. See `solution.md` §9–§17 for dead approaches.

**Handoff for later sessions:** start at **`HANDOFF.md`** (proved list, open blockers, numerics, resume playbook).

## Files

| Path | Role |
|------|------|
| `HANDOFF.md` | Research handoff / resume entry point |
| `solution.md` | Full mathematical writeup + obstruction analysis |
| `src/minmax_quadratic.py` | Exact `m_n`, Paley (prime + prime-power), \(\Phi\), bounds, \(\rho=1\) evec |
| `tests/test_minmax.py` | Pytest suite |
| `evidence/` | Durable tables, \(\rho=1\) proof note, verification JSON |

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
Exact optima first appear at Hamming distance 5 from Paley, and the only 5-edge undercutters are 144 perfect matchings — see `evidence/N10_STRUCTURE.md`. Those 144 form one PΓL(2,9)-orbit (maximizer-drop criterion) — see `evidence/N10_MATCHING_CLASSIFY.md`. Limit existence remains open.