# Clique-flip pair counts on Max-cover matchings (p=5)

**Date:** 2026-07-28  
**Status:** Certified invariant on all 11 stored Max-covers; forall proof open.  
**Existence of \(\lim\alpha_n\) remains OPEN.**

## Observation

For every stored Max-covering perfect matching \(M\) at Paley \(n=26\), the number
\[
N_{\mathrm{flip}}(M)
=
\#\bigl\{(y,F):y\in\mathrm{Max}_{+},\;S_M(y)=1,\;
F\text{ transversal }W\text{-clique},\;\Sigma_F=(1+p)/2=3\bigr\}
\]
takes only two values:

| type (by \(N_{\mathrm{flip}}\)) | \(N_{\mathrm{flip}}\) | \(\#y\) with a flip | count in census |
|--------------------------------|----------------------:|--------------------:|----------------:|
| A | **120** | 96 | 6 |
| B | **24** | 24 | 5 |

**In particular \(N_{\mathrm{flip}}(M)\ge24>0\) on every stored Max-cover**, so clique-flip
(Prop 15.31) applies and \(\Phi(C\oplus M)=\Phi(C)\).

## Proof target

Prove \(N_{\mathrm{flip}}(M)\ge1\) (or \(\ge24\)) for **every** Max-covering perfect matching
at \(p=5\), then matching non-undercut at \(p=5\) by Props 15.33+15.31+15.35.

## What this does not prove

- Forall Max-cover matchings (only 11 stored)
- General \(p\ge7\)
- \(k_\star\) / non-matching undercutters / E(1) / existence of \(L\)

**Do not mark Main Theorem settled.**
