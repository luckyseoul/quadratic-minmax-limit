# g_min residual: m4 moduli + Tr(G²) pin (2026-07-29)

**Status:** Structural path to closed form; **L still OPEN**.

## Setup

Refine 4-sets by pure \(C\)-invariants: \((\mathrm{CR},\kappa,\mathrm{triangle\text{-}type})\).
At \(p=5\) this yields **14** classes; at \(p=7\), **24**. On each class \(m_4\) is constant (certified vs Max+).

## Evec identities

Averaging \(p\,m_4=\sum_r C_{ir}m(\cdots)\) over each class produces a linear system
\[
(pI-M)\mathbf m=\mathbf b
\]
with \(M,b\) **combinatorial** (from \(C\) only). The system is consistent with true Max+ moments.

## Nullity 1

At \(p=5,7\): \(\mathrm{rank}(pI-M)=n_{\mathrm{var}}-1\). Solution space is
\[
\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n,\qquad (pI-M)\mathbf n=0.
\]
The sum-of-disj-\(G\) constraint is already in the system (\(s_{\mathrm{null}}=0\)).

## Pin via \(\mathrm{Tr}(G^2)\)

\[
\mathrm{Tr}(G^2)=E+\frac{2n_{\mathrm{wedge}}}{p^2}+6\sum_{\mathrm{classes}}n_A m_A^2
\]
(\(G_{ee'}^2=m_4^2\) on each disj pairing). Substituting \(\mathbf m(c)\) gives a **quadratic** in \(c\).

**Certified \(p=5\):** one root recovers exact \(m_4\) and \(g_{\min}=-3/65\); the other is spurious (worse \(g_{\min}\)). Correct pin at \(p=5\): \(\mathrm{Tr}(G^2)\) quadratic (or match full \(G\)-spectrum). **Do not** use “max \(g_{\min}\) under PSD” alone — see below.

**\(p=7\):** needs full (not sampled) transitions for a stable quadratic; not yet closed-form. Distance-homogeneous Max+ only at \(p=5\); \(p=7\) has ≥2 Max+ types.

## Spectral check of \(\mathrm{Tr}(G^2)\) at \(p=5\)

Nonzero spectrum of \(G\): \(n/2\) (×1), \(88/13\) (×\(d\)), \(72/13\) (×\(2d\)), \(40/13\) (×\(2d\)) with \(d=13\).
\[
\mathrm{Tr}(G^2)=\Bigl(\tfrac{n}2\Bigr)^2+\sum_j m_j\lambda_j^2=1808.3846\ldots
\]
matches direct Gram computation. A general-\(p\) eigenvalue formula would close the pin.

## PSD / rank scan (anti-thrash; 2026-07-29 late)

On the nullity-1 line \(G(c)\) at \(p=5\):
- \(\mathrm{rank}(G)=\binom{d}{2}-d+1=66\) holds on a **continuum** of \(c\) (with \(\lambda_{\min}\approx0\)); rank alone does **not** pin \(c\).
- Maximizing matrix \(g_{\min}\) subject to PSD gives \(\approx-0.040\) near \(c\approx-0.42\), **strictly better** than the true Max+ value \(-3/65\approx-0.04615\) at true \(c\approx-0.291\).
- **Conclusion:** “max \(g_{\min}\) among PSD \(G(c)\)” is **not** the true selection rule. True selection needs \(\mathrm{Tr}(G^2)\) / full spectrum (or another Max+-free identity). Do not reopen PSD-max as a pin.

Also: \(\mathrm{Tr}(G^2)\) formula must count **\(2\times\)** off-diagonal contributions; missing the factor produced wrong pins until fixed.

## What remains for \(g_{\min}\ge L(p)\)

1. Prove nullity 1 for all primes \(p\ge5\) on this refined class set.
2. Closed form for \(\mathrm{Tr}(G^2)\) (or full \(G\)-spectrum) as a function of \(p\), **or** another \(c\)-pin independent of Max+.
3. Select the correct quadratic root via spectrum / \(\mathrm{Tr}(G^2)\) (**not** PSD-max \(g_{\min}\)).
4. Prove \(g_{\min}(p)\ge L(p)\) (or \(>T(p)\)).

## Related

- Prop 15.51: \(g_{\min}\ge T\Leftrightarrow\min a\ge1/(2p-1)\).
- F15: plain Fréchet on conditional cov is too weak.
- Incomplete Aut-orbits of halfspace ≠ full Max+ (60/260 at \(p=5\)).
