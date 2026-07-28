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

**Certified \(p=5\):** one root recovers exact \(m_4\) and \(g_{\min}=-3/65\); the other is spurious (worse \(g_{\min}\)). Selection: e.g. maximize \(g_{\min}\) among PSD \(G\), or match \(\lambda_{\max}(G)=n/2\).

**\(p=7\):** needs full (not sampled) transitions for a stable quadratic; not yet closed-form.

## What remains for \(g_{\min}\ge L(p)\)

1. Prove nullity 1 for all primes \(p\ge5\) on this refined class set.
2. Closed form for \(\mathrm{Tr}(G^2)\) (or full \(G\)-spectrum) as a function of \(p\), **or** another \(c\)-pin independent of Max+.
3. Solve and prove \(g_{\min}(p)\ge L(p)\) (or \(>T(p)\)).

## Related

- Prop 15.51: \(g_{\min}\ge T\Leftrightarrow\min a\ge1/(2p-1)\).
- F15: plain Fréchet on conditional cov is too weak.
- Incomplete Aut-orbits of halfspace ≠ full Max+ (60/260 at \(p=5\)).
