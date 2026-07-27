# N10-C6: Hamming-6 undercutters of Paley \(C_{10}\) are 6-cycles

**Date:** 2026-07-27  
**Status:** Certified exhaustive at \(n=10\). Supports path-cycle dichotomy for \(k_\star\); existence of \(\lim\alpha_n\) remains **OPEN**.

## Statement (Theorem N10-C6)

Let \(C\) be the Paley conference matrix of order \(n=10\) (\(\Phi(C)=15\), \(m_{10}=13\)).
For every 6-edge set \(F\subset\binom{[10]}{2}\),

\[
\Phi(C\oplus F)<15
\quad\Longleftrightarrow\quad
F\text{ is a }6\text{-cycle and }\Phi(C\oplus F)=13.
\]

There are exactly **360** such 6-cycles (full scan of \(\binom{45}{6}=8{,}145{,}060\) edge sets).

Combined with Theorem N10-S (exactly **144** five-edge undercutters, all perfect matchings):

| \(k\) | # undercutters | structure | \(\Phi\) | \(\Delta_{\max}\) |
|------:|---------------:|:----------|-------:|------------------:|
| 5 | 144 | perfect matching | 13 | 1 |
| 6 | 360 | single 6-cycle | 13 | 2 |

In particular every undercutter of minimal or next-to-minimal cardinality is a **path/cycle graph** with \(k\le n\) and \(\Delta\le 2\).

## Code / tests

- Shipped: `src/n10_cycle_undercutters.py` (`classify_k6_undercutters`, `k5_are_matchings_count`)
- Tests: `test_n10_k6_undercutters_are_cycles`, `test_n10_k5_undercutters_are_matchings_count`

## Relevance to E(1) / \(k_\star\)

Prop 15.20b: \(m_n\ge\Phi(C)-2k_\star\). If \(k_\star\le n\) (path-cycle dichotomy for a closest \(\Phi\)-minimiser), then
\[
m_n\ge\tfrac12 n\sqrt{n-1}-n,
\qquad
\alpha_n\ge\tfrac12\sqrt{1-1/n}-n^{-1/2}\to\tfrac12,
\]
and denseness Prop 6.1–6.2 along \(n=p^2+1\) forces \(\lim\alpha_n=\tfrac12\).

At \(n=10\), \(k_\star=5\le n\) via matchings (C6s give another minimising family at \(k=6\)).

**General path-cycle / \(\Delta\le 2\) dichotomy for all \(n=p^2+1\) is not proved.**
High-degree undercutters exist at larger \(k\) (e.g. \(k=7\), \(\Delta=3\)) but they are not cardinality-minimal; stars never undercut at \(n=10\).

## Not established

- \(k_\star=O(n)\) for general \(p\)
- \(\lim\alpha_n\) exists

**Do not mark Main Theorem / HANDOFF settled.**
