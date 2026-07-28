# E(1) residual after Prop 15.33–15.34

**Date:** 2026-07-27  
**Status:** Existence of \(\lim\alpha_n\) remains **OPEN**.

## What is proved (load-bearing)

1. **Sandwich** \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) and **denseness** (Prop 6.1–6.2).
2. **\(\rho=1\)** on \(n=p^2+1\) (halfspace boolean evecs).
3. **Non-covers cannot undercut** (Prop 15.33): if \(\min_{\mathrm{Max}_{+}}S_F\le-1\) then \(\Phi(C\oplus F)\ge\Phi+2\).
4. Therefore **matching non-undercut reduces to Max-covering perfect matchings only**.
5. Matching flip algebra \(D^2=I\), \(A=C-2D\) (Prop 15.34).
6. Spike criterion (Prop 15.30) is **sufficient** for non-undercut but **not necessary**
   (`e1_criterion_fail_no_undercut.json`: \(\max R=54<60\), non-cover, \(\Phi=75\)).
7. Matching flip **block algebra** (Prop 15.36): \(B=CD+DC\) always commutes with \(C,D\);
   \(B|_{V_+}=2p D_{++}\); \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\). SA Max-covers at \(p=5\) all have
   \(\lambda_{\min}(B)=-6\), \(\mathrm{op}=\sqrt{41}\), with ≥2 distinct \(D_{++}\) types — all non-undercut.

## What would finish matching non-undercut for \(p\ge5\)

Prove: every Max-covering perfect matching \(M\) satisfies \(\Phi(C\oplus M)\ge\Phi(C)\).

Routes:
- **Clique-flip** (Prop 15.31): existence of \((y,F)\) for every Max-cover — certified on all **11** SA covers at \(p=5\) (Prop 15.35), not proved for all.
- **Spectral** (Prop 15.36): prove Max-covers have \(\lambda_{\min}(B)=-6\) and that this forces \(\Phi(A)\ge\Phi(C)\) for \(p\ge5\).
- **Structural**: all tested covers are two-sided, inclusion min+max; forced \(S=1\) (proved).

## Shortest remaining proof edges (2026-07-28)

1. **Upgrade Prop 15.40 → \(m_n\ge\Phi(C)-2\)** for all Seidel \(A\) on \(\rho=1\) orders  
   (edge-minimal undercutters already have gap \(\le2\); need no deeper far undercut).  
   Then gap \(O(1)=o(n^{3/2})\) \(\Rightarrow\) E(1) \(\Rightarrow L=\tfrac12\) by denseness. **Do not claim without this step (F13).**

2. **Or** \(N_{\mathrm{flip}}\ge1\) for every Max-cover PM at all \(p\ge5\) (Prop 15.39 pattern)  
   + path-cycle / \(k_\star=O(n^{3/2})\).

3. **Or** \(k_\star=O(n^{3/2})\) directly.

## What would finish \(L=\tfrac12\)

Matching non-undercut for \(p\ge5\) **plus** one of:
- matching / path-cycle dichotomy, or
- \(k_\star=O(n^{3/2})\) (Prop 15.20d / Max-Lipschitz),

together with \(p=3\) gap \(O(1)\) and denseness.

**New levers:** Prop 15.38 (only \(\Delta=1\) undercuts at n=10 k=5); Prop 15.40 (edge-minimal gap \(\le2\)).

## Census at \(p=5\) (not a proof)

| Object | Result |
|--------|--------|
| Random PM Max-covers in 20k samples | **0** (extremely rare) |
| SA Max-covers checked (unique) | **11**, all two-sided, min+max, \(\Phi=\Phi(C)\), clique-flip, \(\mathrm{op}=\sqrt{41}\) |
| Criterion-fail non-cover | exists; does **not** undercut |
| Full census JSON | `e1_maxcover_full_census.json` |

**Do not mark Main Theorem settled.**
