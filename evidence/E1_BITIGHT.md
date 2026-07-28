# Bi-tight S≡±s covers and star obstruction

**Date:** 2026-07-28  
**Status:** Structural lemmas for no-descent residual; p=5 integral infeasibility + star force certified.  
**Existence of lim α_n remains OPEN.**

## Definition

A flip set \(H\) is **bi-tight of level \(s\)** if \(|H|=sp\) and
\[
S_H\equiv s\text{ on }\mathrm{Max}_{+},\qquad S_H\equiv -s\text{ on }\mathrm{Max}_{-}.
\]
(Always \(\mathbb E_{+}[S]=|H|/p=s\) and \(\mathbb E_{-}[S]=-s\), so constancy is the extra content.)

## Why it matters

1. **Deep tight gap-2 undercutters** that are two-sided are bi-tight (Prop 15.29 + \(S\ge2\), \(\mathbb E=2\) ⇒ \(S\equiv2\); similarly on Max−).
2. **Type I freeness equality** (\(k=2p-1\), \(f_e=2-S\)) produces a tight \(S\equiv2\) cover \(H=G\cup\{e\}\). Descent to \(\Phi\le\Phi-4\) on that edge requires \(\max_{\mathrm{Max}_{-}}S_H\le-2\), hence bi-tight \(H\).

Thus: **if bi-tight size \(2p\) is integrally impossible, Type I freeness-failure cannot cause descent, and deep tight undercutters do not exist.**

## Proved (Prop 15.45)

| Claim | Status |
|--------|--------|
| Wedge \(G=\pm1/p\); wedge \(G^{+}+G^{-}=0\) | Proved |
| **Stars never bi-tight** (need \(\sum(G^{+}+G^{-})=2(2-p)\neq0\)) | Proved, all \(p>2\) |
| Level-1 tight ⇒ star when \(g_{\min}>-1/p\) | Proved |
| Level-2 matching blocked when \(g_{\min}>-1/15\) | Proved |
| \(g_{\min}>-1/p\) at \(p=5\) (\(-3/65>-1/5\)) | Certified |
| Non-star size-\(p\) tight infeasible at \(p=5\) | MILP certified |
| Bi-tight levels 2–4 infeasible at \(p=5\) | MILP certified |
| Deep two-sided \(k=10,12\) infeasible at \(p=5\) | MILP certified |
| \(\min\max S_{-}=2\) for size-\(2p\) with \(s_{+}\ge2\) at \(p=5\) | MILP certified |

## At p=3

Bi-tight exists: undercutting 6-cycles (\(C_6\)) of Paley \(C_{10}\), with \(\Phi=13=\Phi-2\).
Here \(g_{\min}=-1/p\), so star-force fails; avg degree of level-2 is \(1.2>1\).

## Algebraic identity

For any tight \(S\equiv2\) cover of size \(2p\),
\[
\frac1{\binom{2p}{2}}\sum_{e\neq e'\in H}\mathbb E_{+}[f_ef_{e'}]
=
-\frac1{15}.
\]
Same average on Max− for bi-tight. Stars cannot meet the bi-tight \((G^{+}+G^{-})\) sum.

## Residual

- Lift \(g_{\min}>-1/p\) and bi-tight integral infeasibility to all \(p\ge5\).
- Deep **non-tight** gap-2 control (or prove none exist) for full \(m_n\ge\Phi-2\).
- At \(p=5\): Type I + deep-tight closed; residual is deep non-tight with \(k>12\).

Scripts: `src/e1_bitight_infeas.py`, `src/e1_star_bitight_obstruction.py`, `src/e1_deep_cover_hunt.py`.

**Do not mark Main Theorem settled.**
