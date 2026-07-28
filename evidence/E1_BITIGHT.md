# Bi-tight S≡±2 covers of size 2p

**Date:** 2026-07-28  
**Status:** Structural lemma for no-descent residual; p=5 integral infeasibility certified.  
**Existence of lim α_n remains OPEN.**

## Definition

A flip set \(H\) is **bi-tight of level 2** if \(|H|=2p\) and
\[
S_H\equiv 2\text{ on }\mathrm{Max}_{+},\qquad S_H\equiv -2\text{ on }\mathrm{Max}_{-}.
\]
(Always \(\mathbb E_{+}[S]=|H|/p=2\) and \(\mathbb E_{-}[S]=-2\), so constancy is the extra content.)

## Why it matters

1. **Deep tight gap-2 undercutters** that are two-sided are bi-tight (Prop 15.29 + \(S\ge2\), \(\mathbb E=2\) ⇒ \(S\equiv2\); similarly on Max−).
2. **Type I freeness equality** (\(k=2p-1\), \(f_e=2-S\)) produces a tight \(S\equiv2\) cover \(H=G\cup\{e\}\). Descent to \(\Phi\le\Phi-4\) on that edge requires \(\max_{\mathrm{Max}_{-}}S_H\le-2\), hence bi-tight \(H\).

Thus: **if bi-tight size \(2p\) is integrally impossible, Type I freeness-failure cannot cause descent, and deep tight undercutters do not exist.**

## Certified at p=5

| Object | Result |
|--------|--------|
| Fractional bi-tight (uniform \(x_e\)) | **Feasible**, \(\sum x=10\) |
| Integral bi-tight \(\|H\|=10\) | **Infeasible** (HiGHS MILP, full Max± constraints) |
| Integral tight Max+ only \(\|H\|=10\) | **Feasible**; example has \(\Phi=85>\Phi(C)=65\) (spikes up), \(\max_{\mathrm{Max}_{-}}S=10\) (not bi-tight) |

## At p=3

Bi-tight exists: undercutting 6-cycles (\(C_6\)) of Paley \(C_{10}\), with \(\Phi=13=\Phi-2\).

## Algebraic identity

For any tight \(S\equiv2\) cover of size \(2p\),
\[
\frac1{\binom{2p}{2}}\sum_{e\neq e'\in H}\mathbb E_{+}[f_ef_{e'}]
=
\frac{4-2p}{2p(2p-1)}
=
-\frac1{15}
\]
(independent of \(p\)). Same average on Max− for bi-tight.

## Residual

- Prove integral bi-tight infeasible for all \(p\ge5\) (not just \(p=5\)).
- Then Type I no-descent is unconditional for \(p\ge5\); deep tight undercutters vanish for \(p\ge5\).
- Still need deep **non-tight** gap-2 control (or prove none exist) for full \(m_n\ge\Phi-2\).

**Do not mark Main Theorem settled.**
