# Prop 15.38: two-sided size-5 Max-covers on Paley \(C_{10}\)

**Date:** 2026-07-28  
**Status:** **Proved** by exhaustive enumeration of \(\binom{45}{5}\) edge sets.  
**Existence of \(\lim\alpha_n\) remains OPEN** (this is n=10 structure, not E(1)).

## Theorem

Let \(C\) be the Paley conference matrix of order \(n=10\), \(\Phi(C)=15\).  
Among all edge sets \(F\subset\binom{[10]}{2}\) with \(|F|=5\) that are **two-sided Max-covers**
(\(\min_{\mathrm{Max}_{+}}S_F\ge1\) and \(\max_{\mathrm{Max}_{-}}S_F\le-1\)):

| \(\Delta(F)\) | count | \(\Phi(C\oplus F)\) | undercuts? |
|-------------:|------:|--------------------:|:-----------|
| 1 | **144** | **13** | **yes** |
| 2 | 8730 | \(\{15,17,19\}\) | no |
| 3 | 7920 | \(\{17,19\}\) | no |
| 4 | 360 | \(19\) | no |

**Conclusion:** undercutters among two-sided \(k=5\) Max-covers are **exactly** the 144
perfect matchings (N10-S). Every two-sided cover with \(\Delta\ge2\) has \(\Phi\ge\Phi(C)\).

Total two-sided \(k=5\) Max-covers: \(144+8730+7920+360=17154\).

## Method

- Max\(\pm\) = boolean \(\pm3\)-eigenvectors of \(C\) (12 each up to the shipped SVD free-variable scan; full \(\pm\) closure as needed).
- Enumerate all \(\binom{45}{5}=1{,}221{,}759\) five-edge sets.
- Filter two-sided Max-covers via score tables on Max\(\pm\).
- Exact \(\Phi\) by full cube (\(2^9\)) using shipped `form_Q`.
- Parallel re-verification: `src/n10_twosided_k5_classify.py` (ProcessPool).

## Link to E(1)

Supports the **path-cycle / low-\(\Delta\) undercutter** pattern: at the matching
cardinality, high-degree two-sided Max-covers do **not** undercut.  
Does **not** by itself prove \(k_\star=O(n^{3/2})\) for general \(p\).

Evidence JSON: `e1_n10_twosided_k5_classify.json`.  
**Do not mark Main Theorem settled.**
