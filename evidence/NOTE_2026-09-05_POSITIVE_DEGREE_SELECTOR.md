# Positive-degree selection: exact minimax scope

2026-09-05. **All-orders variational identities; cross-order comparison OPEN.**

For \(A\in\mathcal S_N\), \(\beta>0\), set
\(\epsilon=\binom N2\bmod2\), \(q=\cosh(2\beta)-1\), and use
the exact nonnegative coefficients

\[
 Z_A(\beta)=(\cosh\beta)^\epsilon\sum_{k=0}^{D}a_{A,k}q^k,
 \qquad D=\left\lfloor\tfrac12\binom N2\right\rfloor,
 \qquad d=D+1.
\]

Coefficients beyond a signing's actual degree are zero. Define

\[
 b_{A,k}=a_{A,k}q^k,\qquad
 S=\min_A\sum_k b_{A,k},\qquad
 M=\min_A\max_k b_{A,k},\qquad
 R=\max_k\min_A b_{A,k}.
\]

The constant coefficient is one, so all these quantities are positive.

## 1. One degree per signing, with the quantifiers retained

For every signing, its largest weighted term lies between \(1/d\)
times its positive sum and the full sum. Minimizing gives

\[
 M\le S\le dM,
 \qquad 0\le\log S-\log M\le\log d=O(\log N).            \tag{1}
\]

Thus one appropriately selected degree captures the optimized pressure
to \(O(\log N)\), provided that the signing is minimized before its
maximizing degree is chosen. This does not replace \(M\) by \(R\).

## 2. The legitimate mixed-strategy exchange

Let \(\pi\) range over probability measures on complete signings and
let \(\lambda\) range over probability vectors on the \(d\) degrees.
The finite bilinear minimax theorem gives the exact identity

\[
 V:=\min_\pi\max_k\mathbb E_\pi b_{A,k}
   =\max_\lambda\min_A\sum_k\lambda_k b_{A,k}.             \tag{2}
\]

This is the permitted exchange: it introduces mixed strategies, not
a single universal degree. The exact comparisons are

\[
 \frac Sd\le V\le M\le S,
 \qquad R\le V.                                          \tag{3}
\]

Indeed, for any \(\pi\),
\(\max_k\mathbb E_\pi b_{A,k}\ge d^{-1}\mathbb E_\pi\sum_kb_{A,k}
\ge S/d\). Taking pure \(\pi\)'s proves \(V\le M\). Equivalently,
the uniform degree distribution already attains \(S/d\) on the right
side of (2). These formulas prove an \(O(\log N)\) approximation
after the valid mixed exchange. They do not assert that an optimal
\(\lambda\) is supported on a single degree.

## 3. A pure-degree exchange with a temperature-uniform entropy loss

There is a separate, weaker estimate for \(R\) which suffices to
preserve zero-temperature slopes. Put \(m_N=\min_A\Phi(A)\), and
recall the exact coefficient formula

\[
 d_k(m,\epsilon)=
 \frac{\prod_{j=0}^{k-1}[m^2-(2j+\epsilon)^2]}{2^k(2k)!},
 \qquad a_{A,k}=\mathbb E_x d_k(|Q_A(x)|,\epsilon).
\]

On the nonnegative integers of fixed parity \(\epsilon\), each
\(d_k(m,\epsilon)\) is nondecreasing. It is zero below the degree
threshold, and above that threshold every factor is nonnegative and
increasing in \(m\). The two antipodal states at \(\Phi(A)\), and
\(\Phi(A)\ge m_N\), therefore give, for every signing and every degree,

\[
 a_{A,k}\ge2^{1-N}d_k(\Phi(A),\epsilon)
          \ge2^{1-N}d_k(m_N,\epsilon).                    \tag{4}
\]

The integer \(m_N\) has parity \(\epsilon\). Consequently

\[
\begin{aligned}
 R&\ge2^{1-N}\max_k d_k(m_N,\epsilon)q^k\\
  &\ge\frac{2^{1-N}}d\,
       \frac{\cosh(\beta m_N)}{(\cosh\beta)^\epsilon},
 \qquad
 S\le\frac{\cosh(\beta m_N)}{(\cosh\beta)^\epsilon}.       \tag{5}
\end{aligned}
\]

The upper bound evaluates the pressure on a norm minimizer and uses
\(|Q_A|\le m_N\) pointwise. Hence

\[
 0\le\log S-\log R\le(N-1)\log2+\log d.                 \tag{6}
\]

This is not an \(o(N)\) fixed-temperature exchange. It is a bound
uniform in \(\beta\), with the ordinary spin entropy as its extensive
cost. In particular, at \(\beta=c/\sqrt N\), the pure-degree selector
\(\log[(\cosh\beta)^\epsilon R]/(cN)\) differs from
\(\alpha_N=m_N/N^{3/2}\) by at most
\(\log2/c+O(\log N/(cN))\). Thus it has the same zero-temperature
slope in the regime of large \(c\), without asserting a common degree
that minimizes every coefficient or a saddle point among pure strategies.

## 4. Endpoint of the present deduction

The exact positive expansion permits a one-degree maximum with the
original minimization order, a genuine mixed-strategy exchange up to
\(O(\log N)\) in pressure, and a pure-degree exchange with the explicit
entropy loss (6). These statements hold for all orders and all complete
signings and retain the full coefficient system.

They do not compare \(R,S\), or the scalar optimized coefficients,
between different orders. The remaining target is transport or
convergence of the extensive coefficient rates identified in the
companion moment-comparison note. Neither bilinear minimax nor
positivity alone supplies that transport. No convergence claim follows
from the present identities.
