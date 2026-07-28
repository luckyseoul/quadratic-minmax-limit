# Continuous Γ-bound on Max-cover matchings (p=5)

**Date:** 2026-07-27  
**Status:** Certified on all 11 stored Max-covers; not a forall proof.  
**Existence of lim α_n remains OPEN.**

## Statement (certified pattern)

Let \(C\) be Paley of order \(n=26\) (\(p=5\)), \(M\) a Max-covering perfect matching, and
\(\Gamma(z)\) the Γ-pairing matrix of Prop 15.32 on relative signs \(z\in\{\pm1\}^{13}\).
Write \(\lambda_{\max}(z)=\lambda_{\max}(\Gamma(z))\). Then on the full level \(S_M=\langle c,z\rangle=-p\):

\[
\min_{S_M=-p}\lambda_{\max}(z)\;\ge\;9.38758\;>\;\frac{2p(m-1)}{m}=\frac{120}{13}\approx9.2308.
\]

Hence the **continuous** bound \(\tfrac m2\lambda_{\max}(z)\ge60=p(m-1)\) holds for **every**
\(z\) on the spike level (not merely some \(z\)).

| Cover class (by \(\min\lambda_{\max}\)) | \(\min\lambda_{\max}\) | cont min \(\tfrac m2\lambda\) |
|----------------------------------------|----------------------:|------------------------------:|
| Type A (several covers) | 9.38758666 | 61.0193 |
| Type B (several covers) | 9.392527 | 61.0514 |

All 11 stored Max-covers fall into these two classes. Random non-covers often have
some \(z\) with cont \(<60\) (e.g. criterion-fail example: cont min \(\approx57.8\), discrete
\(\max R=54\)).

## Discrete gap

On the best \(z\), \(\lambda_{\max}\approx10.55\), continuous \(\approx68.57\), discrete
\(\max R=60\) (tight threshold). Sign of the top eigenvector achieves \(R=60\).
Goemans–Williamson theory is **not** quite enough for a free proof: the SDP value is
\(\approx63.8\) (evec non-flat), and \(\alpha_{\mathrm{GW}}\cdot\mathrm{SDP}\approx56<60\).

## Link to clique-flip / non-undercut

Discrete \(\max R\ge60\) \(\Leftrightarrow\) spike criterion \(\Leftrightarrow\) \(\Phi(C\oplus M)\ge\Phi(C)\).
Clique-flip constructs a point on the level with \(R=60\). All 11 covers admit clique-flip
and have MITM \(\Phi=\Phi(C)\).

**Additionally:** every tested residue-\(1\bmod4\) matching (not only covers) had at least one
clique-flip from an \(S_M=1\) vector. Max-covers are forced residue \(1\) (Prop 15.35).

## Clique-flip concentration note

On the stored Max-cover, among 168 vectors with \(S_M=1\), only **24** admit a
\(\Sigma=3\) transversal \(W\)-clique. For many other \(S=1\) vectors one has
\(\mathbb E[\Sigma^2]=1\) with \(\Sigma\in\{\pm1\}\) only (no clique-flip). So a
per-vector second-moment argument cannot use an arbitrary \(S=1\) vector — only a
positive-density subset works. Residue-\(1\) non-covers in sampling often had
**every** \(S=1\) vector clique-flip-capable (higher \(\mathbb E[\Sigma^2]\)).

## What this does not prove

- \(\min_z\lambda_{\max}\ge120/13\) for every Max-cover (only certified samples)
- Discrete \(\max R\ge60\) from the continuous bound alone
- Clique-flip for every Max-cover
- Matching dichotomy / \(k_\star\) / E(1) / existence of \(L\)

Scripts: inline census in session; covers from `e1_maxcover_full_census.json` sources.
**Do not mark Main Theorem settled.**
