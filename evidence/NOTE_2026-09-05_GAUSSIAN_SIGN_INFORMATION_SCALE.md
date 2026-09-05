# Information scale of the actual Gaussian-sign cross-block law

2026-09-05. **All-orders, method-scoped theorem; convergence OPEN.**

This applies the new quenched iid obstruction to the actual centered
Gaussian-sign construction, rather than its annealed certificate. It
distinguishes weakly regularized rounding from fixed-strength rounding.
For actual pressure-minimizing hosts, sufficiently weak rounding has
\(o(n)\) information relative to independent cross signs and is excluded
by the quenched theorem. In contrast, every fixed positive strength has
\(\Omega(n)\) information for every norm-capped host. The latter claim
includes the original singular Gaussian construction, but it is a
lower bound on its sign-law entropy, not an invalid Gaussian determinant
comparison at a singular covariance.

## 1. Exact covariance and normalization

Fix a complete symmetric zero-diagonal signing \(A\) of order \(n\ge2\).
The host is held fixed while generating a cross block. Let \(U,V\) be
its actual zero-field Gibbs covariances at any opposite positive and
negative temperatures. Write

\[
 p=\operatorname{tr}(AU),\quad q=\operatorname{tr}(AV),\quad
 \alpha=\frac{p+q}{2n},\quad L=\|A\|_{\rm op},
\]
\[
 H=A\otimes A-\alpha(A\otimes I_n+I_n\otimes A).
                                                               \tag{1}
\]

Let the extremal eigenvalues of \(A\) be \(a>0\) and \(-b<0\).
The reviewed covariance-rounding note supplies the exact negative spectral
edge and upper bounds below; the lower estimate for \(\mu\) is derived here:

\[
 -b/2\le\alpha\le a/2,\qquad
 \mu=-\lambda_{\min}(H)=ab+\alpha(a-b)>0,
\]
\[
 \frac{n-1}{2}\le\mu\le L^2,
 \qquad \|H\|_{\rm op}\le2L^2,
 \qquad \operatorname{diag}H=0.                            \tag{2}
\]

Here the exact eigenvalues of \(H\) are
\(xy-\alpha(x+y)\), with \(x,y\) eigenvalues of \(A\).
Their negative minimum occurs at the two mixed corners \((a,-b)\)
and \((-b,a)\). The same-sign corners are nonnegative. Also
\(ab\ge n-1\), obtained by summing
\((a-x)(b+x)\ge0\) over the spectrum and using
\(\operatorname{tr}A=0\), \(\operatorname{tr}A^2=n(n-1)\).
These facts give the lower bound on \(\mu\) in (2).

Let \(m=n^2\), use column vectorization of cross blocks, and define

\[
 \Sigma_\rho=I_m+\rho H/\mu,
 \qquad G\sim N(0,\Sigma_\rho),\qquad
 \operatorname{vec}B=\operatorname{sign}G,
 \qquad 0\le\rho\le1.                                    \tag{3}
\]

Denote this sign law by \(Q_{A,\rho}\) and the independent fair sign
law by \(P_0\). The covariance is positive definite for \(\rho<1\),
and positive semidefinite and singular at \(\rho=1\). Its
diagonal is one in both cases. The original canonical construction is
precisely \(\rho=1\); \(\rho<1\) is its stated regularization.

The exact squared Frobenius norm is

\[
 \boxed{\|H\|_F^2=n^2(n-1)[(n-1)+2\alpha^2].}             \tag{4}
\]

Indeed, the square of the first tensor term contributes
\((\operatorname{tr}A^2)^2\). Its cross trace against the other two
terms vanishes because \(\operatorname{tr}A=0\). The square of their
sum contributes \(2n\operatorname{tr}A^2\). This proves (4).

## 2. One-sided determinant upper bound

For \(0\le\rho<1\), data processing under componentwise sign and the
Gaussian relative entropy formula give

\[
 D(Q_{A,\rho}\Vert P_0)
 \le-\tfrac12\log\det\Sigma_\rho
 \le\frac{\rho^2\|H\|_F^2}{4(1-\rho)\mu^2}.              \tag{5}
\]

No upper bound on the positive eigenvalues of \(H/\mu\) is needed
for the second inequality. In fact, if \(X=\rho H/\mu\), then
\(\operatorname{tr}X=0\), all its eigenvalues satisfy \(x\ge-\rho\),
and

\[
 x-\log(1+x)\le\frac{x^2}{2(1-\rho)}\qquad(x\ge-\rho).
\]

This follows by integrating \(x/(1+x)\), separately on the negative
and positive intervals. Summing over the eigenvalues proves (5).
In particular, a sufficient condition for \(o(n)\) sign-law entropy is

\[
 \frac{\rho^2}{1-\rho}
 \frac{n^2(n-1)[(n-1)+2\alpha^2]}{\mu^2}=o(n).             \tag{6}
\]

The factor \(\|H\|_F^2/\mu^2\) is a negative-edge-normalized squared
Frobenius quantity; it is not the usual stable rank unless \(\mu\)
also equals the full operator norm. For fixed \(0<\rho<1\), (6) is
equivalent to \(\mu^2\gg n^3+n^2\alpha^2\). The exact determinant
condition \(-\log\det\Sigma_\rho=o(n)\) is another sufficient
condition, not asserted equivalent to the Frobenius test in general.
Neither Gaussian upper bound is asserted at \(\rho=1\).

## 3. A general lower bound for the actual sign-law entropy

For any probability law \(Q\) on \(\{-1,1\}^m\), put
\(C=\mathbb E_Q[bb^T]\). Then

\[
 \boxed{D(Q\Vert P_0)\ge
       \frac{\|C-I_m\|_F^2}{4\|C\|_{\rm op}}.}            \tag{7}
\]

To prove this, expose coordinates in any fixed order and define
\(z_i=\mathbb E_Q[b_i\mid b_1,\ldots,b_{i-1}]\). The relative
entropy chain rule and the elementary binary entropy inequality give

\[
 D(Q\Vert P_0)\ge\tfrac12\sum_i\mathbb E_Q z_i^2.
\]

Let \(c_i=(C_{ij})_{j<i}\). For every real vector \(v\) on those
coordinates, conditional expectation and Cauchy--Schwarz imply

\[
 |v^Tc_i|^2
 =\left|\mathbb E_Q z_i\sum_{j<i}v_jb_j\right|^2
 \le\mathbb E_Q z_i^2\,v^TC_{<i,<i}v
 \le\mathbb E_Q z_i^2\,\|C\|_{\rm op}\|v\|_2^2.
\]

Taking the dual Euclidean norm and summing gives (7), because
\(\sum_i\|c_i\|_2^2=\tfrac12\|C-I_m\|_F^2\).

Now let \(Q\) be the sign law of any Gaussian correlation matrix
\(\Sigma\), including a singular one. The Gaussian sign identity is
entrywise:

\[
 C=\tfrac2\pi\arcsin\circ\Sigma.
\]

The scalar inequality \(|\arcsin r|\ge|r|\) gives
\(\|C-I_m\|_F^2\ge(4/\pi^2)\|\Sigma-I_m\|_F^2\).
Also

\[
 \|C\|_{\rm op}\le\|\Sigma\|_{\rm op}.                  \tag{8}
\]

For (8), expand \((2/\pi)\arcsin r\) in positive odd powers of
\(r\); their coefficients sum to one. Each odd Schur power is
\(\Sigma\circ R\), where \(R\) is a positive semidefinite correlation
matrix. Since \(\Sigma\preceq\|\Sigma\|_{\rm op}I_m\), positivity
of the Schur product gives
\(\Sigma\circ R\preceq\|\Sigma\|_{\rm op}I_m\).
The positive series converges entrywise, and hence in finite-dimensional
matrix norm, including at correlations of absolute value one.
This proves (8) and

\[
 \boxed{D(Q_\Sigma\Vert P_0)\ge
       \frac{\|\Sigma-I_m\|_F^2}{\pi^2\|\Sigma\|_{\rm op}}.}           \tag{9}
\]

This lower bound concerns the discrete sign law directly. It does not
use absolute continuity of a possibly singular Gaussian relative to
an independent Gaussian.

Applied to (3), writing \(h_+=\lambda_{\max}(H)\ge0\), it yields

\[
 D(Q_{A,\rho}\Vert P_0)
 \ge\frac{\rho^2\|H\|_F^2}
             {\pi^2\mu(\mu+\rho h_+)}
 \ge\frac{\rho^2 n^2(n-1)^2}
             {\pi^2(1+2\rho)L^4},
 \qquad 0\le\rho\le1.                                   \tag{10}
\]

Here (2) and (4) justify the second inequality.

## 4. Norm caps force extensive information at fixed strength

A complete signing satisfies the useful interpolation bound

\[
 L^2\le16\Phi(A).                                        \tag{11}
\]

For completeness, let \(K_{\mathbb R}\) be the real bilinear norm
\(\sup_{x,y\in[-1,1]^n}|x^TAy|\). Its maximum occurs at vertices.
For sign vectors \(x,y\), put \(u=(x+y)/2\), \(v=(x-y)/2\).
Then \(x^TAy=u^TAu-v^TAv\). The zero diagonal makes \(Q_A\)
multiaffine, so its absolute maximum on the real cube is \(\Phi(A)\).
Consequently \(K_{\mathbb R}\le4\Phi(A)\).

Writing complex vectors in real and imaginary parts gives the complex
operator norm \(\|A\|_{\ell_\infty\to\ell_1}\le4K_{\mathbb R}\).
The other endpoint norm is \(\|A\|_{\ell_1\to\ell_\infty}=1\).
Riesz--Thorin interpolation halfway between these two endpoints gives
\(\|A\|_{\ell_2\to\ell_2}^2\le4K_{\mathbb R}\), proving (11).
No optimizer property enters this interpolation step.

Thus, for every fixed norm cap \(\Phi(A)\le Cn^{3/2}\), (10) implies

\[
 \boxed{D(Q_{A,\rho}\Vert P_0)
 \ge\frac{\rho^2}{256\pi^2(1+2\rho)C^2}
          n(1-1/n)^2.}                                   \tag{12}
\]

At every fixed \(\rho>0\), this is \(\Omega_{C,\rho}(n)\).
In particular, fixed-strength canonical Gaussian-sign rounding on a
norm-capped host never meets the \(o(n)\) entropy premise. For
\(\rho<1\), data processing implies that the Gaussian determinant
entropy also cannot be \(o(n)\). At \(\rho=1\), only the discrete
lower bound is being used.

These assertions apply to each fixed-host law. A lower bound for each
conditional law is not automatically a lower bound for an arbitrary
mixture over hosts, and no such mixture assertion is made.

## 5. Actual minimizers and the weak-rounding regime

The norm cap applies to actual minimizers, without any conference
optimality hypothesis. Norm minimizers satisfy
\(\Phi(A)=m_n\le(1/2+o(1))n^{3/2}\) by the proved construction bound.

For clarity, also take an actual global half-product minimizer at
\(\beta=c/\sqrt n\), so \(a_A(\beta)=R_n(\beta)\), with fixed
\(c>0\). If \(q_{\max},q_{\min}\) are its extreme spin energies and
\(w=(q_{\max}-q_{\min})/2\), the two antipodal extremal states in
each one-sided partition function give

\[
 a_A(\beta)\ge\beta w-(n-1)\log2.
\]

The uniform energy mean is zero, so \(\Phi(A)\le2w\). Comparison
with a norm minimizer yields

\[
 \Phi(A)\le2m_n+\frac{2(n-1)\log2}{\beta}
 \le\bigl(1+2\log2/c+o(1)\bigr)n^{3/2}.                   \tag{13}
\]

Hence every sequence of these actual pressure minimizers is covered
by (12). It is not necessary to identify them with norm minimizers.

For actual opposite-temperature phase covariances, the two mean energies
have opposite signs. Therefore

\[
 |\alpha|
 =\frac{|\mathbb E_+Q_A+\mathbb E_-Q_A|}{n}
 \le\frac{\Phi(A)}n\le C\sqrt n.                         \tag{14}
\]

This holds at any positive temperature used to form \(U,V\), not just
the temperature at which the host was selected. Equations (2), (4),
and (5) give the uniform upper bound

\[
 D(Q_{A,\rho}\Vert P_0)
 \le\frac{\rho^2}{1-\rho}
     n^2\left(1+\frac{2\alpha^2}{n-1}\right)
 \le\frac{\rho^2}{1-\rho}(1+4C^2)n^2,
 \qquad0\le\rho<1.                                      \tag{15}
\]

In particular, on every such actual-minimizer family,

\[
 \rho_n=o(n^{-1/2})\quad\Longrightarrow\quad
 D(Q_{A_n,\rho_n}\Vert P_0)=o(n).                         \tag{16}
\]

There is deliberately no claim settling every intermediate vanishing
\(\rho_n\) scale between this sufficient condition and fixed strength.

## 6. New finite-step consequence, not an annealing argument

Use the complete proof in `NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md`.
It defines \(K_0=4/(3\sqrt\pi)>1/\sqrt2\) and, for fixed \(c>0\),
\(0<t\le1\),

\[
 \Delta(c,t)=c(\sqrt{2t}K_0-1)-2\log2.
\]

The host-free quenched bound in that note proves, for any proposal law
\(Q\) on cross blocks,

\[
 \mathbb E_Q\min_{A'}F_{A',B}(t)
 \ge[c\sqrt{2t}K_0-2\log2]n
       -c\sqrt{tnD(Q\Vert P_0)}-o_{c,t}(n),               \tag{17}
\]

where
\(F_{A',B}(t)=\log\mathbb E\cosh[\eta(Q_{A'}(x)-Q_{A'}(y))
+\gamma x^TBy]\), \(\eta=(c/\sqrt n)\sqrt{1-t/2}\), and
\(\gamma=(c/\sqrt n)\sqrt{t/2}\). The same note shows that
\(D(Q\Vert P_0)=o(n)\) gives probability tending to zero for the event
\(\min_{A'}F_{A',B}(t)\le2R_n(c/\sqrt n)+o(n)\), whenever
\(\Delta(c,t)>0\).

Combining this theorem with (16), weakly regularized actual
Gaussian-sign rounding at \(\rho_n=o(n^{-1/2})\) fails the desired
finite-step comparison both in mean and with any probability bounded
away from zero at those parameters. Its mean excess is at least
\((\Delta(c,t)-o(1))n\). The internal host in (17) may even be chosen
after seeing the cross block, because the lower bound is host-free.

This is a new application to a dependent sign law and to the actual
average log pressure. It is not the old Gaussian annealed floor and
does not use the quadratic Gram--Schmidt proxy. On the other hand,
fixed-strength Gaussian-sign rounding, including the singular canonical
law, has extensive information by (12) and is outside this exclusion.
Extensive information is not a proof of success; its remaining quenched
or selected-outcome comparison is still open. Neither an unrestricted
\(\min_B\), arbitrary dependent rounding, nor the original convergence
problem is retired by these results.
