# Dimension-uniform extensive coefficient comparison

2026-09-05. **All-orders coefficient theorem and an open cross-order target.**

Use the exact normalized positive polynomial

\[
 Z_A(c)=(\cosh(c/\sqrt N))^\epsilon
          \sum_k a_{A,k}(\cosh(2c/\sqrt N)-1)^k,
 \qquad \epsilon=\binom N2\pmod2,
\]

for a complete signing \(A\in\mathcal S_N\). The coefficients do not
depend on \(c\). Their exact formula, proved in the companion truncation
note, is

\[
 a_{A,k}=\frac1{2^k(2k)!}
   \mathbb E\prod_{j=0}^{k-1}
               \bigl(Q_A(x)^2-(2j+\epsilon)^2\bigr).       \tag{1}
\]

For each integer energy, the product is nonnegative: it stays positive
up to the appropriate degree and then contains a zero factor. In
particular no cancellations between positive and negative products are
being suppressed in (1).

**Theorem.** Fix \(0<\theta_0\le\theta_1<\infty\). There exist constants
\(b>0\) and \(N_0\), depending only on this interval, such that for
all \(N\ge N_0\), all integers \(\theta_0N\le k\le\theta_1N\),
and all complete signings \(A\),

\[
 \boxed{
 b\,\frac{\mathbb E|Q_A|^{2k}}{2^k(2k)!}
 \le a_{A,k}\le
 \frac{\mathbb E|Q_A|^{2k}}{2^k(2k)!}.}                    \tag{2}
\]

There is no norm-cap or optimizer hypothesis. Thus in the actual
extensive-degree range, the exact positive coefficients and the full
even moments have the same exponential rate, uniformly in the signing.

## Proof of the lower bound

The upper bound follows termwise from (1). Put
\(L=\|Q_A\|_{2k}\). The two antipodal states attaining \(\Phi(A)\)
give

\[
 L\ge2^{-(N-1)/(2k)}\Phi(A)
   \ge \frac{2^{-1/(2\theta_0)}}{\pi\sqrt2}N^{3/2},       \tag{3}
\]

using the proved \(\Phi(A)\ge N\sqrt{N-1}/\pi\) and \(N\ge2\).
Set \(\delta=2^{-1/(2\theta_0)}/(2\pi\sqrt2)>0\). On the set
\(E=\{|Q_A|\ge L/2\}\), one has \(|Q_A|\ge\delta N^{3/2}\).
Also

\[
 \mathbb E[|Q_A|^{2k}\mathbf1_E]
 \ge(1-2^{-2k})\mathbb E|Q_A|^{2k}
 \ge\tfrac34\mathbb E|Q_A|^{2k}.                          \tag{4}
\]

For sufficiently large \(N\), uniformly on \(E\),

\[
 0\le r_j:=\frac{(2j+\epsilon)^2}{Q_A^2}
 \le\frac{4\theta_1^2}{\delta^2N}\le\frac12,
 \qquad
 \sum_{j<k}r_j\le\frac{4\theta_1^3}{\delta^2}.
\]

Since \(\log(1-r)\ge-2r\) for \(0\le r\le1/2\),

\[
 \prod_{j<k}(Q_A^2-(2j+\epsilon)^2)
 \ge |Q_A|^{2k}\exp(-8\theta_1^3/\delta^2)
 \quad\text{on }E.                                      \tag{5}
\]

Outside \(E\), the product remains nonnegative. Equations (1), (4),
and (5) prove (2), for example with
\(b=(3/4)\exp(-8\theta_1^3/\delta^2)\).

## Optimized coefficient rates and the missing comparison

Define the separately optimized scalar quantities

\[
 a^{\min}_{N,k}=\min_{A\in\mathcal S_N}a_{A,k},\qquad
 V_{N,k}=N^{-3/2}\min_{A\in\mathcal S_N}\|Q_A\|_{2k}.
\]

The same multiplicative comparison (2) holds after minimizing, so for
\(k=\lfloor\theta N\rfloor\), fixed \(\theta>0\), Stirling's formula
gives

\[
 \frac1N\log\bigl(N^{-k}a^{\min}_{N,k}\bigr)
 =2\theta\log\frac{eV_{N,k}}{2\sqrt2\,\theta}+o(1).        \tag{6}
\]

The error is uniform on compact positive \(\theta\) intervals. The
norm lower bound, the norm upper constructions, and the following
inequality keep \(V_{N,k}\) uniformly above zero and below infinity,
so replacing \(k/N\) by \(\theta\) in (6) is harmless:

\[
 2^{-(N-1)/(2k)}\alpha_N\le V_{N,k}\le\alpha_N,
 \qquad \alpha_N=m_N/N^{3/2}.                              \tag{7}
\]

The lower bound in (7) holds for every signing before minimizing; the
upper bound follows by evaluating the moment on a norm minimizer.

Consequently, convergence of the scalar coefficient rate on the left
of (6) for an unbounded set of fixed \(\theta\)'s would imply convergence
of \(\alpha_N\), without imposing a value. Indeed, (6) would give
\(V_{N,\lfloor\theta N\rfloor}\to v_\theta\); (7) then gives

\[
 v_\theta\le\liminf_N\alpha_N\le\limsup_N\alpha_N
       \le2^{1/(2\theta)}v_\theta.
\]

The uniform upper bound on \(\alpha_N\) and unboundedness of these
\(\theta\)'s force the lower and upper limits to agree.

No such cross-order convergence or transport of coefficient rates is
proved here. In particular, separate minimization of each coefficient
does not identify a single signing minimizing their weighted sum, and
the minimum must not be moved through that sum. The theorem resolves
the comparison between the new exact coefficient object and actual
extensive moments; it does not resolve the remaining order comparison.
