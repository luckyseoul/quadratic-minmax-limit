# Sublinear positive-cone truncation loses extensive pressure

2026-09-05. **All-orders method-scoped theorem; convergence OPEN.**

This addresses the precise polynomial expansion used by the exact
order-six balanced-profile proof. It does not use fourth moments or a
spectral-defect shell. The conclusion applies to actual norm minimizers
and to actual symmetric-pressure minimizers, not just an abstract moment
measure. It rules out an approximation step; it is not a convergence
theorem or an obstruction to extensive-degree coefficient transport.

Let \(A\in\mathcal S_N\), \(N\ge2\), and set

\[
 Q_A(x)=\sum_{i<j}A_{ij}x_ix_j,\quad
 \Phi(A)=\max_x|Q_A(x)|,\quad
 \beta=c/\sqrt N,\quad
 Z_A(c)=\mathbb E_x\cosh(\beta Q_A(x)),
\]

with uniform spins. Put \(\epsilon=\binom N2\bmod2\) and
\(q=\cosh(2\beta)-1=2\sinh^2\beta\). Every energy has parity
\(\epsilon\), so there is an exact finite expansion

\[
 Z_A(c)=(\cosh\beta)^\epsilon\sum_{k\ge0}a_{A,k}q^k,
 \qquad a_{A,0}=1,\qquad a_{A,k}\ge0.                      \tag{1}
\]

Define its degree-\(K\) truncation by retaining \(0\le k\le K\)
in (1), and denote it by \(Z_{A,\le K}(c)\).

**Theorem.** Fix \(c>0\), a finite constant \(C\), and any sequence
of nonnegative integers \(K_N=o(N)\). Uniformly over all actual complete signings
with \(\Phi(A)\le C N^{3/2}\),

\[
 0\le\frac1N\log Z_{A,\le K_N}(c)=o(1).                  \tag{2}
\]

Moreover, uniformly over the same class,

\[
 \liminf_{N\to\infty}\frac1N
       \log\frac{Z_A(c)}{Z_{A,\le K_N}(c)}
 \ge\frac c\pi-\log2.                                   \tag{3}
\]

In particular, for every fixed \(c>\pi\log2\), a sublinear-degree
truncation loses a strictly positive normalized pressure, and captures
at most an exponentially vanishing fraction of the partition function.
The statements include any sequence of norm minimizers and any sequence
of symmetric-pressure minimizers at this same fixed \(c\).

## Exact coefficient bound

For a nonnegative integer \(m\) of parity \(\epsilon\), define
\(H_{m,\epsilon}(q)\) through

\[
 (\cosh\beta)^\epsilon H_{m,\epsilon}(\cosh2\beta-1)
 =\cosh(m\beta).
\]

For even \(m\) this is \(T_{m/2}(1+q)\); for odd \(m\) it is
the polynomial \(R_{(m-1)/2}(1+q)\) used in the order-six proof.
Its exact coefficients are

\[
 H_{m,\epsilon}(q)=\sum_{k\ge0}d_k(m,\epsilon)q^k,
 \qquad
 d_k(m,\epsilon)=
 \frac{\prod_{j=0}^{k-1}[m^2-(2j+\epsilon)^2]}
      {2^k(2k)!}.                                       \tag{4}
\]

The empty product is one. The coefficients vanish after degree
\((m-\epsilon)/2\), because the product then contains a zero factor.
Every coefficient before that zero is nonnegative. One direct verification
of (4), avoiding any asymptotic approximation, is to differentiate the
defining hyperbolic identity. The resulting equation is

\[
 q(q+2)H''+[(1+\epsilon)q+1]H'
             -\frac{m^2-\epsilon^2}{4}H=0,\qquad H(0)=1.
\]

Its coefficient recurrence is

\[
 d_{k+1}=
 \frac{m^2-(2k+\epsilon)^2}{4(k+1)(2k+1)}d_k,
\]

which gives (4) and terminates at the stated degree. Averaging (4) with
\(m=|Q_A(x)|\) gives (1). Crucially,

\[
 0\le a_{A,k}q^k
 \le\frac{(\Phi(A)\sinh\beta)^{2k}}{(2k)!}.              \tag{5}
\]

This bound retains the exact central-factorial coefficient object; it
does not replace the entire pressure by a low moment.

## Uniform truncation estimate

Under the norm cap, \(\lambda_N:=\Phi(A)\sinh(c/\sqrt N)\le DN\)
for all sufficiently large \(N\), where \(D>0\) is a constant depending
only on \(c,C\). If \(1\le K=K_N=o(N)\), eventually \(2K\le DN\).
The terms \((DN)^{2k}/(2k)!\) then increase for \(0\le k\le K\).
Using \((2K)!\ge(2K/e)^{2K}\), (5) gives

\[
\begin{aligned}
 \log Z_{A,\le K}(c)
 &\le\epsilon\log\cosh\beta+
       \log(K+1)+2K\log\frac{eDN}{2K}\\
 &=o(N).                                                  \tag{6}
\end{aligned}
\]

The last assertion follows from
\((K/N)\log(N/K)\to0\). If \(K=0\), the logarithm is merely
\(\epsilon\log\cosh\beta=O_c(N^{-1})\). Since the constant
coefficient is one, all these logarithms are nonnegative. This proves
(2), uniformly over the claimed class.

## Full pressure and actual optimizer scope

There are at least two spin states, \(x\) and \(-x\), attaining
\(|Q_A(x)|=\Phi(A)\). Their contribution and
\(\cosh z\ge e^{|z|}/2\) imply the exact entropy bound

\[
 \log Z_A(c)\ge\beta\Phi(A)-N\log2.
\]

The all-orders Gaussian lower bound proved in `CORE.md` is
\(\Phi(A)\ge N\sqrt{N-1}/\pi\). Thus every complete signing satisfies

\[
 \frac1N\log Z_A(c)
 \ge\frac c\pi\sqrt{1-\frac1N}-\log2.                    \tag{7}
\]

Combining (6) and (7) proves (3). This uses an already proved norm lower
bound, not a new claim about its optimality or its possible limiting value.

For actual norm minimizers, the proved conference-construction upper
bound and ratio-dense transfer give
\(\Phi(A)=m_N\le(1/2+o(1))N^{3/2}\), so the norm cap in (2) applies.

For an actual symmetric-pressure minimizer \(A_N^{\rm pr}\), comparison
with a norm minimizer instead gives

\[
 \log Z_{A_N^{\rm pr}}(c)
 \le\frac c{\sqrt N}m_N
 \le(\tfrac c2+o(1))N.
\]

Using the exact entropy bound once more yields

\[
 \Phi(A_N^{\rm pr})
 \le\left(\frac12+\frac{\log2}{c}+o(1)\right)N^{3/2}.      \tag{8}
\]

This is also a bounded norm cap for fixed \(c>0\), and proves the
pressure-minimizer assertion. It is intentionally not identified with
the sharper cap for a norm minimizer.

## What remains

The omitted mass in (3) comes from the exact positive polynomial of an
actual optimum, not from cancellations between signed coefficients.
Therefore retaining only degrees \(k\le K_N=o(N)\) cannot support
an \(o(N)\)-accurate pressure comparison at these fixed temperatures.
This does not exclude a sparse selection of extensive degrees. Indeed,
for each fixed signing and profile, its single largest weighted degree
term approximates the logarithm of the full positive sum within
\(O(\log N)\), since the polynomial has \(O(N^2)\) terms.
The needed next object is the extensive-degree range \(k=\Theta(N)\)
and its behavior under changes of order. No comparison of those
coefficient rates between orders is established here. The original
convergence problem remains open.
