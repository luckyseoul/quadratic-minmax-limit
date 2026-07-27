# E(1) rigidity attack: \(k_\star=o(n^{3/2})\)

**Date:** 2026-07-27  
**Status:** Partial — criterion sharpened; full rigidity **not proved**. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Criterion (Prop 15.20b)

On \(n=p^2+1\) with Paley \(C\) (\(\rho=1\)):
\[
m_n\ge\Phi(C)-2k_\star,
\qquad
k_\star
=\min\bigl\{\text{best switching Hamming distance to }C
\text{ among }\Phi\text{-minimisers}\bigr\}.
\]
**E(1)** holds if \(k_\star=o(n^{3/2})\).

## Certified facts

### \(n=10\) (exact)

| Object | best-\(k\) to Paley | \(\Phi\) |
|--------|-------------------:|---------:|
| Paley \(C\) | 0 | 15 |
| Matching undercutters (144) | **5** \(=k_\star\) | 13 |
| Other \(\Phi=13\) matrices (SA) | up to **≥15** | 13 |

So the **minimising** Hamming distance is \(k_\star=5=O(n)\), even though some minimisers sit far from Paley after switching. Edge lip uses the **closest** minimiser: \(m_{10}\ge15-10=5\) (true; actual 13).

Campaign: `src/e1_rigidity_k_threshold.py` / `{SCRATCH}/campaign_rigidity_k.log`  
— for every \(K\in\{5,6,7,8,10,12,15\}\), SA found \(\Phi=13\) matrices with best-\(k\ge K\).

### Boolean maximiser counts (rho=1)

| \(p\) | \(n\) | \(\lvert\mathrm{Max}_+\rvert\) (approx) | \(n^{3/2}\) |
|------:|------:|------------------------------------------:|------------:|
| 3 | 10 | 12 (exact) | 32 |
| 5 | 26 | ~\(10^2\) (2M sample) | 133 |

Aut upper bound \(\lvert\mathrm{PGL}(2,p^2)\rvert\cdot2\sim n^3\) is **larger** than \(n^{3/2}\), so \(\lvert\mathrm{Max}\rvert\le\lvert\mathrm{Aut}\rvert\) alone does **not** give \(k_\star=o(n^{3/2})\).

### Necessary covering

Any undercutter (\(\Phi(A)<\Phi(C)\)) must satisfy \(S_F(y)\ge1\) for every positive maximiser \(y\) (integer arithmetic). Hence undercutting flip sets are covering designs for \(\mathrm{Max}_+\). Inclusion-minimal covers have size \(\le\lvert\mathrm{Max}_+\rvert\), but a global \(\Phi\)-minimiser need not be an inclusion-minimal cover (far \(n=10\) optima).

## What would finish rigidity

1. Prove \(k_\star=O(n)\) (or any \(o(n^{3/2})\)) for all large \(n=p^2+1\): e.g. every closest undercutter is a matching / bounded-degree flip of Paley.  
2. Or prove the matching dichotomy: \(m_n=\min\bigl(\Phi(C),\min_M\Phi(C\oplus M)\bigr)\) over perfect matchings \(M\) — then degree lip gives gap \(\le n\).  
3. Or permanent relative gap / non-existence pair.

## Sufficient claim (still open)

**Claim (path-cycle / \(k_\star\le n\) dichotomy).** On \(n=p^2+1\), there exists a \(\Phi\)-minimiser whose disagreement graph with Paley \(C\) (after switching) is a disjoint union of paths and cycles — equivalently \(k_\star\le n\) (since such graphs have \(|E|\le n\)).

**Weaker form (matching only):** \(k_\star\le n/2\) via a perfect-matching minimiser.

**If either form holds**, then by Prop 15.20b
\[
m_n\ge\Phi(C)-2k_\star\ge\tfrac12 np-2n,
\]
so \(\alpha_n\ge\tfrac12\sqrt{1-1/n}-2n^{-1/2}\to\tfrac12\). Combined with \(\limsup\alpha_n\le\tfrac12\) and denseness Prop 6.1–6.2 along \(n_k=p_k^2+1\),
\[
\lim_{n\to\infty}\alpha_n=\tfrac12.
\]

**Evidence for the claim (n=10, certified):**
- \(k_\star=5\le n\) via 144 perfect-matching undercutters (N10-S / N10-C).
- All **360** Hamming-6 undercutters are single **6-cycles** with \(\Phi=13=m_{10}\) (Theorem N10-C6; `evidence/N10_CYCLE_UNDERCUTTERS.md`).
- So both minimal and next-to-minimal undercutters are path/cycle graphs (\(\Delta\le2\), \(k\le n\)).
- Stars never undercut (deg \(3..8\): 0 hits). High-\(\Delta\) undercutters appear only at larger \(k\) (e.g. \(k=7\), \(\Delta=3\)) and are not cardinality-minimal.
- Far \(m_{10}\)-optima with best-\(k\in[8,16]\) exist; many do **not** contain a \(k\le6\) undercutter as a subgraph, and single-edge deletion need not reduce them to one. This does **not** affect \(k_\star=\min\). Extraction of a sparse undercutter from *every* minimiser fails; only *existence* of a sparse minimiser is required for E(1).
- At \(n=26\), exact-Φ SA finds no undercut of \(\Phi=65\), consistent with \(k_\star=0\).

**Obstacle:** \(\Phi\) is only 2-Lipschitz, so near-optimal \(\Phi\) does **not** force small Hamming distance (far optima exist). Need a structural reason that a *closest* \(\Phi\)-minimiser has \(k\le n\) (path-cycle / matching). General \(p\) open.

## Not established

- \(k_\star=o(n^{3/2})\) in general  
- matching dichotomy for all \(p\)  
- \(\lim\alpha_n\) exists  

**Do not mark Main Theorem / HANDOFF settled.**
