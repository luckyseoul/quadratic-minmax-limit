# Max-cover matching flip spectrum (Prop 15.36)

**Date:** 2026-07-27  
**Status:** Algebra proved; spectral pattern certified on all SA Max-covers;  
**forall Max-cover ⇒ Φ≥Φ(C) still OPEN.**  
**Existence of lim α_n remains OPEN.**

## Setup

Paley conference \(C\) of order \(n=p^2+1\), perfect matching \(M\), signed matching matrix
\(D\) (\(D_{ij}=C_{ij}\) on edges of \(M\), else 0), \(A=C-2D=C\oplus M\).
Eigenspaces \(V_\pm=\ker(C\mp pI)\), dim \(n/2\) each. Write
\[
D=\begin{pmatrix}D_{++}&D_{+-}\\D_{-+}&D_{--}\end{pmatrix}
\quad\text{on }V_+\oplus V_-.
\]

## Algebra (proved for every matching)

1. \(D^2=I\) (Prop 15.34).
2. \(B:=CD+DC\) always commutes with \(C\) and with \(D\).
3. \(B|_{V_+}=2p\,D_{++}\), \(B|_{V_-}=-2p\,D_{--}\).
4. \(\mathrm{tr}(D_{++})=n/(2p)\) always; for \(y\in\mathrm{Max}_+\subset V_+\),
   \(S_M(y)=\tfrac12 y^\top D y=\tfrac12 a^\top D_{++}a\) with \(a=V_+^\top y\), \(\|a\|^2=n\).
5. \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\).

## Certified at \(p=5\) (Max-covers only)

Every SA Max-covering perfect matching tested satisfies:

| Property | Value |
|----------|------:|
| \(\lambda_{\min}(B)\) | \(-6\) exactly |
| \(\|A\|_{\mathrm{op}}\) | \(\sqrt{41}=\sqrt{p^2+16}\) |
| spike \(\max R\) on \(S_M=-p\) | \(\ge60=p(m-1)\) (tight) |
| exact MITM \(\Phi(A)\) | \(65=\Phi(C)\) |
| clique-flip (Prop 15.31) | yes |

**At least two spectral types** of \(D_{++}\) occur among Max-covers:

- **Simple type:** \(\mathrm{spec}(D_{++})=\{-3/5,\,0^{(8)},\,(4/5)^{(4)}\}\);  
  \(S\) on Max± takes values \(\{1,5,9\}\) (e.g. hist \(1^{168},5^{80},9^{12}\)).
- **Mixed type:** \(\mathrm{spec}(D_{++})\) involves additional values  
  \(\{\pm0.4,\,0.25359,\,0.94641,\ldots\}\) still with \(\lambda_{\min}(D_{++})=-3/5\);  
  \(S\in\{1,5\}\) only (hist \(1^{156},5^{104}\)).

Both types share \(\lambda_{\min}(B)=-6\) and non-undercut. Random (non-cover) matchings have
generic \(B\)-spectra with \(\lambda_{\min}(B)\in[-10,-8]\) and \(\|A\|_{\mathrm{op}}\in[6.8,7.0]\).

## Contrast \(p=3\) (undercutting covers)

All 144 Max-cover undercutters share one spectrum class:
\(\|A\|_{\mathrm{op}}\approx3.933\), \(B\)-eigs involving golden-ratio values
\(\{\pm(\sqrt5-1),\pm(\sqrt5+1),6\}\). Non-covers have strictly larger op-norm.
So **Max-covers minimise \(\|A\|_{\mathrm{op}}\) among perfect-matching flips** at both \(p=3\) and \(p=5\);
the discrete cube max undercuts only at \(p=3\).

## What this does **not** prove

- Every Max-cover has \(\lambda_{\min}(B)=-6\) (certified, not proved).
- \(\lambda_{\min}(B)=-6\Rightarrow\Phi(A)\ge\Phi(C)\) for matching flips when \(p\ge5\).
- Clique-flip existence for every Max-cover (Prop 15.31 open quantifier).
- Matching dichotomy / \(k_\star\) / E(1) / existence of \(L\).

## Evidence

- `evidence/e1_maxcover_spectrum.json` — two spectral types, both non-undercut
- `src/e1_maxcover_spectrum.py` — SA census of \(D_{++}\) / op / clique-flip / MITM
- Prior: `e1_maxcover_full_census.json` (11 covers, all \(\mathrm{op}=\sqrt{41}\))

**Do not mark Main Theorem settled.**
