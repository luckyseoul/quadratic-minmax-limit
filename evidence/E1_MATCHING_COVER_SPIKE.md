# Matching Max-covers at \(n=26\): two-sided covers that do **not** undercut

**Date:** 2026-07-27  
**Status:** Certified numerics + parity lemma. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Correction

Earlier `E1_MATCHING_MAXCOVER.md` reported no Max-covering perfect matching at \(n=26\) (SA never exceeded \(\min S=-1\)). A stronger 2-swap SA **does** find Max-covers (\(\min_{\mathrm{Max}_{+}}S=1\)). Those covers are typically **two-sided** and still **fail to undercut** \(\Phi=65\).

## Parity lemma (proved)

For odd prime \(p\), \(n=p^2+1\) is even and \(n/2\) is **odd**. Hence for every perfect matching \(M\) and every \(x\in\{\pm1\}^n\),
\[
S_M(x)=\sum_{\{i,j\}\in M}C_{ij}x_ix_j
\]
is an odd integer. In particular \(S_M(x)\neq0\), so either \(\min_{\mathrm{Max}_{+}}S_M\ge1\) (Max-cover) or \(\min_{\mathrm{Max}_{+}}S_M\le-1\).

**Corollary.** If \(M\) is **not** a Max-cover then some \(y\in\mathrm{Max}_{+}\) has \(S_M(y)\le-1\), hence
\[
|Q_{C\oplus M}(y)|=\Phi-2S_M(y)\ge\Phi+2,
\]
so \(\Phi(C\oplus M)\ge\Phi+2\). Non-covering perfect matchings **strictly raise** \(\Phi\).

## Census (Paley \(p=5\), \(n=26\))

Script: `src/e1_matching_mins_sa.py`, `src/e1_n26_matching_cover_census.py`  
JSON: `e1_matching_minS_sa.json`, `e1_n26_matching_covers_census.json`, `e1_n26_matching_cover_example.json`

| Quantity | Value |
|----------|------:|
| Seeds (min-\(S_{+}\) SA, 25k 2-swaps) | 48 |
| Max-covering PMs found | **3** (seeds 1012, 1013, 1032) |
| Of which two-sided (\(\max_{\mathrm{Max}_{-}}S\le-1\)) | **3** |
| \(\max_{\mathrm{Max}_{\pm}}\|Q\|\) on those covers | **63** \(=\Phi-2\) |
| Exact MITM \(\Phi(C\oplus M)\) | **65** \(=\Phi(C)\) for all 3 |
| Undercuts of \(\Phi=65\) | **0** |

So at \(n=26\), two-sided matching Max-covers exist and drop every boolean \(\pm p\)-eigenvector to \(|Q|\le63\), but **non-eigenvectors spike by \(+2\)** and restore \(\Phi=\Phi(C)\).

Contrast \(n=10\): 144 two-sided matching covers achieve \(\Phi=13=m_{10}\) (global Max-determination, no harmful spike).

## Relevance to E(1)

- Matching dichotomy is subtler than “no matching cover for \(p\ge5\)”: covers exist, but may be blocked by spikes.
- Prop 15.26 (local maximality after matching flips) holds for both \(p=3\) and \(p=5\); the difference is **global** Max-determination (n=10 only among tested).
- Still need a general spike theorem for two-sided covers when \(p\ge5\), or \(k_\star=O(n^{3/2})\).

**Do not mark Main Theorem settled.**
