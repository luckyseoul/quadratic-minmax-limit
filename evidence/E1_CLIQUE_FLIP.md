# Clique-flip construction for matching Max-cover spikes

**Date:** 2026-07-27  
**Status:** Construction proved sufficient; design constants certified; forall-Max-cover-\(M\) existence of flip still open for general \(p\ge5\).  
**Existence of \(\lim\alpha_n\) remains OPEN.**

## Setup

Paley \(C\) of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), perfect matching \(M\), \(A=C\oplus M\).

For \(y\in\mathrm{Max}_{+}\) and \(F\subset[n]\),
\[
Q_C(y^{\oplus F})=\Phi-2p|F|+4\sum_{i<j\in F}C_{ij}y_iy_j.
\]
If \(F\) is a clique in the \(y\)-switched graph \(W_{ij}=C_{ij}y_iy_j\) (\(|F|=r\)) then
\[
Q_C(y^{\oplus F})=\Phi-2r(p-r+1).
\]
Also \(S_M(y^{\oplus F})=S_M(y)-2\sum_{i\in F}\chi_i(y)\) with \(\chi_i=C_{i\pi(i)}y_iy_{\pi(i)}\) (one term per matched edge meeting \(F\); require \(F\) transversal to \(M\)).

## Arithmetic: only \(r\in\{1,p\}\) reach the threshold with a full clique

Need \(Q_C(y^{\oplus F})\ge\Phi-2p\), i.e. \(r(p-r+1)\le p\). The quadratic \(r(p-r+1)\) equals \(p\) at \(r=1\) and \(r=p\) and is **strictly larger** for \(1<r<p\). Thus only \(r=1\) or \(r=p\) work for full-clique flips.

- **\(r=1\)** forces \(S_M(y)\in\{2-p,-2-p\}\subset(-\infty,0)\) for the image to have \(S_M=-p\). Impossible for Max-covers (\(S_M\ge1\) on \(\mathrm{Max}_{+}\)).
- **\(r=p\)** requires \(\sum_{i\in F}\chi_i=(S_M(y)+p)/2\) and \(k_{+}=(p+\Sigma)/2\in\mathbb Z\cap[0,p]\), i.e. \(S_M(y)+3p\equiv0\pmod4\).

| \(p\bmod4\) | admissible \(S_M(y)\) for \(r=p\) |
|------------|----------------------------------|
| \(1\) (e.g. 5,13) | \(S_M(y)\equiv1\pmod4\) |
| \(3\) (e.g. 3,7) | \(S_M(y)\equiv3\pmod4\) |

At \(p=3\), undercutting matchings have \(S_M(\mathrm{Max}_{+})\subseteq\{1,5\}\) (no value \(\equiv3\bmod4\) in range), so **\(r=p\) clique-flip is impossible** — consistent with undercuts existing.  
At \(p=5\), \(\mathbb E[S_M]=2.6<3\), so every Max-cover has some \(y\) with \(S_M(y)=1\equiv1\pmod4\), unlocking \(r=p\).

## Design constants (Paley \(p=5\), certified)

| Object | Count |
|--------|------:|
| Seidel-consistent \(p\)-sets (\(C_{ab}C_{ac}C_{bc}=1\) all triples) | **390** |
| Max\(_{+}\) extensions per consistent \(F\) (pattern \(\pm\psi\)) | **60** (exactly 30+30) |
| Transversal consistent \(F\) per random matching (50 samples) | min **236**, mean **244.6** (all \(>0\)) |
| SA Max-covers tested with clique-flip | all succeeded; MITM \(\Phi=\Phi(C)\) |

Scripts: `src/e1_clique_flip_covers.py`, evidence JSON `e1_clique_flip_covers.json`.

## Sufficiency theorem (proved)

**Proposition 15.31.** Let \(M\) be a perfect matching Max-cover of Paley \(\rho=1\) order \(n=p^2+1\). Suppose there exist \(y\in\mathrm{Max}_{+}\) and a transversal set \(F\) of size \(p\) such that:
1. \(S_M(y)=s_0\) with \(s_0+3p\equiv0\pmod4\) and \(1\le s_0\le p\);
2. \(F\) is a clique in \(W_{ij}=C_{ij}y_iy_j\);
3. \(\sum_{i\in F}\chi_i(y)=(s_0+p)/2\).

Then \(Q_C(y^{\oplus F})=\Phi-2p\), \(S_M(y^{\oplus F})=-p\), and \(Q_A(y^{\oplus F})=\Phi\), hence \(\Phi(C\oplus M)\ge\Phi(C)\).

## What remains for matching non-undercut at all \(p\ge5\)

Prove that for every Max-covering perfect matching \(M\), a pair \((y,F)\) as in Prop 15.31 exists (or the equivalent level-set form \(\max_{S_M=-p}Q_C\ge\Phi-2p\)).

For \(p=5\): \(s_0=1\) is forced; need transversal \(W\)-clique with \(\Sigma=3\). Certified on all SA-found covers; not proved for every cover.

For \(p\equiv3\pmod4\), \(p\ge7\): need \(s_0\equiv3\pmod4\) (e.g. \(s_0=3\)) on some Max\(_{+}\) vector of every cover.

## Path to \(L=\tfrac12\) (still blocked)

Matching non-undercut for all \(p\ge5\) + matching dichotomy + \(p=3\) gap \(O(1)\) \(\Rightarrow\) E(1) \(\Rightarrow L=\tfrac12\).  
Matching dichotomy and non-matching undercutters remain open. **Do not mark settled.**
