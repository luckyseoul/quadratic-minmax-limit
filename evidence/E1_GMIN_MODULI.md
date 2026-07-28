# g_min residual: m4 moduli + Tr(G²) pin (updated 2026-07-30)

**Status:** Prop 15.53 shipped; structural path to closed form; **L still OPEN**.  
**Code:** `src/e1_gmin_moduli.py` → `evidence/e1_gmin_moduli.json`.

## Pairing identity (Prop 15.53.1, proved)

On every 4-set with \(|\kappa|=1\), the three pairing \(C\)-products are a perm of \((1,1,-1)\) or \((-1,-1,1)\), so
\[
g_{\min}=-\max\bigl\{|m_4(S)|:|S|=4,\,|\kappa(S)|=1\bigr\}.
\]
Certified \(p=5,7\).

## Setup (refined classes)

Stratify 4-sets by pure \(C\)-invariants \((\mathrm{type}_6,\mathrm{ext}\text{-sum hist})\):
- \(\mathrm{type}_6\): \(S_4\)-canonical 6-tuple of edge signs;
- ext-sum hist: histogram of \(\sum_{v\in S}C_{rv}\) for \(r\notin S\).

At \(p=5\): **37** classes, all with **constant** \(m_4\) (vs Max+).  
Bare \(C\)-types alone are **not** constant — need the external refinement.

## Evec system + nullity 1 (cert \(p=5\))

Averaging \(p\,m_4=\sum_r C_{ir}m(\cdots)\) over each class → combinatorial \(A\mathbf m=\mathbf b\)
(RHS uses only 2-design \(m_2=C/p\)). At \(p=5\): \(\mathrm{rank}=36\), **nullity 1**:
\[
\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n.
\]

## Pin via \(\mathrm{Tr}(G^2)\)

Dot form (Max+ pairwise): \(K_{ab}=((y_a\cdot y_b)^2-n)/2\), \(\mathrm{Tr}(G^2)=\|K\|_F^2/N^2\).
\[
\mathrm{Tr}(G^2)=\tfrac14\bigl(\mathbb E[\mathrm{dot}^4]-2n\,\mathbb E[\mathrm{dot}^2]+n^2\bigr).
\]
**Proved Max+-free:** \(\mathbb E[\mathrm{dot}^2]=n+n(n-1)/p^2\) from \(\mathbb E[yy^\top]=I+C/p\).

**Still Max+-dependent:** \(\mathbb E[\mathrm{dot}^4]\) (or closed \(G\)-spectrum).

Substituting \(\mathbf m(c)\) into edge-\(\mathrm{Tr}(G^2)\) → quadratic in \(c\).  
**Certified \(p=5\):** selected root (larger \(g_{\min}\)) recovers \(g_{\min}=-3/65\).  
**F16:** do not pin by max \(g_{\min}\) under PSD over the whole line (continuum of rank-66 PSD points; max PSD \(g_{\min}\approx-0.040\neq-3/65\)).

## Spectral check / tight obstruction (Prop 15.55)

Nonzero spectrum of \(G\):
- \(p=3\): \(8\) (×5), \(n/2=5\) (×1) — **λ_max > n/2**
- \(p=5\): \(n/2=13\) (×1), \(88/13\) (×d), \(72/13\) (×2d), \(40/13\) (×2d), \(d=13\)
- \(p=7\): \(n/2=25\) (×1), then \(2160/409,2016/409,\ldots\) (mult d,2d,…)

**Prop 15.55:** if λ_max(G)=n/2 is simple, no Max+-tight size-2p cover exists (bi-tight empty).  
Certified p=5,7; open to prove for all primes p≥5. Code: `e1_gmin_tight_obstruction.py`.

\[
\mathrm{Tr}(G^2)=23509/13\ (p=5).
\]
A general-\(p\) eigenvalue formula would also close the Tr(G²) pin.

## What remains for \(g_{\min}\ge L(p)\)

1. Prove refined classes + nullity 1 for all primes \(p\ge5\).
2. Max+-free \(\mathbb E[\mathrm{dot}^4]\) or \(G\)-spectrum as a function of \(p\).
3. Select correct quadratic root (larger \(g_{\min}\) among Tr(G²) roots — not PSD-max).
4. Prove \(g_{\min}(p)\ge L(p)\) or \(>T(p)\).

## Related

- Prop 15.51–15.54 in `solution.md`.
- F15: plain Fréchet too weak; F16: PSD-max pin banned; **F18:** PGL+Frob+sign orbit of halfspace = 60/260 (not full Max+).
- `e1_gmin_pgl_orbit.json`, `e1_gmin_cbound.json`, `e1_gmin_abound.json`.
