# Full-strength failure on actual construction-cap hosts

2026-09-05. **An actual-host application of the spectral-deficit theorem;
not an exact- or leading-near-minimizer assertion.**

Choose fixed (c>0\), (0<r<1\), and (t=1\) such that the gap

\[
 g_r(c,1)=c(\sqrt2 K_0-1)-2\log2
                        -\frac{c^2}{2\pi}\arcsin r>0
\]

in the full-strength spectral-deficit note is positive. Such choices
exist by the conditional-noise theorem. Along an infinite sequence of
orders (N\), there are actual complete signings (A_N\) satisfying

\[
 \Phi(A_N)\le(1/2+o(1))N^{3/2},
\]

whose canonical full-strength Gaussian-sign cross-block law, formed from
their actual opposite-temperature Gibbs covariances, satisfies

\[
 \mathbb E F_{A_N,B}(1)-2a_{A_N}(c/\sqrt N)
 \ge[g_r(c,1)-o(1)]N.                                      \tag{1}
\]

Its probability of attaining its own paired endpoint plus (o(N)\)
tends to zero. In particular, full-strength canonical rounding need
not succeed even for actual complete hosts meeting the best established
leading norm construction cap. There is **no claim** that these hosts
satisfy \(\Phi(A_N)=m_N+o(N^{3/2})\),
\(a_{A_N}=R_N+o(N)\), edge-local minimality, or exact minimality.

## 1. Complete signings and the inherited filler bound

Take (L=2^k\), let (H_L\) be the symmetric Walsh matrix
\((H_L)_{u,v}=(-1)^{u\cdot v}\), and put
\(B_L=H_L-\operatorname{diag}H_L\). Direct character orthogonality
gives (H_L^2=LI\). Thus (B_L\) is a complete signing and

\[
 \|B_L\|_{\rm op}\le\sqrt L+1,\qquad
 \Phi(B_L)\le\tfrac12L(\sqrt L+1).                           \tag{2}
\]

Fix (K\ge100/r\). Define

\[
 \ell=2\lceil K\sqrt L/2\rceil,\qquad
 s=\lfloor L^{1/4}\rfloor,\qquad
 m=2s\ell,\qquad N=L+m.
\]

In particular, (N/L\to1\), (m=O_K(N^{3/4})\), and
\(\ell=K\sqrt L+O(1)\).

Use exactly the complete twin-module construction of the reviewed
near-minimizer counterfamily, with this old signing (B_L\). Each of
the (s\) modules has two communities of \(\ell\) vertices, with
constant community pattern
\(\begin{pmatrix}1&1\\1&-1\end{pmatrix}\) and zero diagonal;
twins are paired within each community. Pair-to-old edges are independent
signs times \((1,-1)\); intermodule pair-to-pair blocks are independent
signs times \(\begin{pmatrix}1&-1\\-1&1\end{pmatrix}\).
These rules specify a sign on every off-diagonal entry of (A_N\).

The filler annihilates every pair-even module subspace. On the old plus
pair-odd subspace it is

\[
 \mathcal F=\begin{pmatrix}0&\sqrt2R^T\\\sqrt2R&2W\end{pmatrix}.
\]

For a fixed unit vector, the squared sum of its independent Rademacher
quadratic coefficients is at most eight. The exponential-moment bound
is consequently \(\Pr(|v^T\mathcal Fv|>u)\le2e^{-u^2/16}\).
A (1/4\)-net has at most (9^N\) points and controls the operator
norm within a factor two. Thus, exactly as in that reviewed proof,

\[
 \Pr(\|\mathcal F\|_{\rm op}>16\sqrt N)
 \le2e^{-(4-\log9)N}<1
\]

for large (N\). Fix any realization in the complementary event.
No thermal event or unverified Gibbs approximation is required here.

The same elementary block bound gives, uniformly over spin configurations,

\[
 |Q_{A_N}(z,x)-Q_{B_L}(z)|
 \le\frac{m\ell}{2}+16N\sqrt m+8m\sqrt N
 =O_K(N^{11/8})=o(N^{3/2}).                                 \tag{3}
\]

The first term bounds the internal module energies; the other two
bound the old--new and new--new filler forms. Equations (2)--(3)
prove the stated construction cap on \(\Phi(A_N)\).

## 2. Exact spectral separation and the actual centering

Every module's two-dimensional community-constant subspace is invariant.
In its normalized constant basis the internal matrix is

\[
 \begin{pmatrix}\ell-1&\ell\\\ell&-\ell+1\end{pmatrix},
\]

with exact eigenvalues \(\pm M\), where
\(M=\sqrt{2\ell^2-2\ell+1}\). Each has total multiplicity (s\).
The other pair-even eigenvalues are (+1\) or \(-1\). On the
old plus pair-odd remainder, the internal new operator has norm one,
and the filler bound and (2) give operator norm at most

\[
 (\sqrt L+1)+1+16\sqrt N\le18\sqrt N
\]

for large (N\). Consequently, for sufficiently large (N\), the
extreme eigenvalues of the full signing are exactly (M,-M\), and
all other eigenvalues have magnitude at most (18\sqrt N\). Also
\(M\ge\ell\ge(K/2)\sqrt N\).

Form the canonical law using the actual opposite Gibbs phases at any
positive temperature. Their means have opposite signs, so the actual
centering obeys

\[
 |\alpha|\le\Phi(A_N)/N\le\sqrt N
\]

for large (N\). Because the two spectral extremes have equal magnitude,
the exact canonical normalization is \(\mu=M^2\), independently
of that actual centering. Thus the eigenvalues of (T=H/\mu\) are

\[
 \frac{xy-\alpha(x+y)}{M^2},
\]

where (x,y\) are eigenvalues of (A_N\). If at least one of the
pair is not a spectral extreme, their absolute value is at most

\[
 \frac{36}{K}+\frac2K\left(1+\frac{36}{K}\right)
 \le\frac{39}{K}<r.
\]

The two same-sign extreme sectors have (T\)-eigenvalues
\(1-2\alpha/M>0\) and \(1+2\alpha/M>0\); the two mixed extreme
sectors have eigenvalue exactly \(-1\). Each mixed sector has
dimension (s^2\). Therefore the full-strength spectral deficit is
exactly

\[
 \boxed{\quad
 \operatorname{tr}[-rI-T]_+=2(1-r)s^2=O(\sqrt N)=o(N).
 \quad}                                                       \tag{4}
\]

This uses actual Gibbs data only through the proved centering bound;
neither a conditional conference covariance model nor optimizer
independence has been substituted for the actual law.

## 3. The actual quenched consequence and its limits

Apply the full-strength spectral-deficit theorem with (4). It gives

\[
 \mathbb E\min_{A'} F_{A',B}(1)
 \ge[c+g_r(c,1)-o(1)]N.
\]

The construction cap meanwhile gives
\(2a_{A_N}(c/\sqrt N)\le2(c/\sqrt N)\Phi(A_N)
\le cN+o(N)\). These facts prove (1), even allowing the internal
host to be reselected after the cross block is seen.

The spectral-deficit theorem's probability transfer also applies.
Its coupling error here is (O_{c,K}(N^{3/4})\), so Markov's
inequality bounds the additional probability loss at a fixed linear
pressure threshold by (O_{c,K}(N^{-1/4})\). Together with the
repaired-law exponential tail this gives a vanishing success probability,
not a claimed exponential tail for the original singular law.

The old arbitrary minimizer in the earlier near-minimizer construction
has no proved (O(\sqrt L)\) spectral bound. Replacing (B_L\) by
such a minimizer does not justify the spectral separation above. Thus
the exact-minimizer full-strength problem, as well as an analogous
claim for leading-order near-minimizers, remains open.

No numerical job, signing census, or Gaussian simulation was run.
