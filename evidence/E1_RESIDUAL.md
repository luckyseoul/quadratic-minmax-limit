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

## What would finish matching non-undercut for \(p\ge5\)

Prove: every Max-covering perfect matching \(M\) satisfies \(\Phi(C\oplus M)\ge\Phi(C)\).

Routes:
- **Clique-flip** (Prop 15.31): existence of \((y,F)\) for every Max-cover — certified on all 8 SA covers at \(p=5\), not proved.
- **Spectral**: Max-covers at \(p=5\) have \(\|A\|_{\mathrm{op}}=\sqrt{41}\); upgrade to \(\Phi(A)\ge\Phi(C)\).

## What would finish \(L=\tfrac12\)

Matching non-undercut for \(p\ge5\) **plus** one of:
- matching / path-cycle dichotomy, or
- \(k_\star=O(n^{3/2})\) (Prop 15.20d / Max-Lipschitz),

together with \(p=3\) gap \(O(1)\) and denseness.

## Census at \(p=5\) (not a proof)

| Object | Result |
|--------|--------|
| Random PM Max-covers in 20k samples | **0** (extremely rare) |
| SA Max-covers checked (unique) | **11**, all two-sided, min+max, \(\Phi=\Phi(C)\), clique-flip, \(\mathrm{op}=\sqrt{41}\) |
| Criterion-fail non-cover | exists; does **not** undercut |
| Full census JSON | `e1_maxcover_full_census.json` |

**Do not mark Main Theorem settled.**
