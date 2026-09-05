# Finite balanced steps: rounding bounds and annealing obstructions

## Scope and notation

This is a follow-up to
`NOTE_2026-09-05_INTEGRAL_CROSS_BLOCK_COVARIANCE_ROUNDING.md`.
It proves a genuine finite-step Gram--Schmidt bound, a uniform obstruction
to that quadratic proxy, and a separate uniform obstruction for the
**actual annealed Gaussian-sign construction**. These are method-scoped
statements, not counterexamples to a selected integral cross block,
an optimized finite-step comparison, or convergence.

Let \(\mathcal S_n\) denote complete symmetric zero-diagonal signings,
and put
\[
Q_A(x)=\tfrac12x^TAx,\qquad Z_A(s)=\mathbb E_xe^{sQ_A(x)},\qquad
a_A(s)=\tfrac12\bigl(\log Z_A(s)+\log Z_A(-s)\bigr),
\]
\[
R_n(s)=\min_{A\in\mathcal S_n}a_A(s),\qquad
\Phi(A)=\max_x|Q_A(x)|,\qquad
m_n=\min_{A\in\mathcal S_n}\Phi(A).
\]
All spin expectations without a measure subscript are uniform.
Fix \(c>0\) and \(0<t\le1\), independently of \(n\), and write
\[
\beta=\frac c{\sqrt n},\qquad
\eta=\beta\sqrt{1-t/2},\qquad
\gamma=\beta\sqrt{t/2},\qquad
k=\frac{\gamma^2}{2}=\frac{c^2t}{4n}.
\]
For a fixed host \(A\in\mathcal S_n\) and an integral cross block
\(B\in\{-1,1\}^{n\times n}\), define
\[
F_{A,B}(t)=\log\mathbb E_{x,y}\cosh\!\left[
 \eta\bigl(Q_A(x)-Q_A(y)\bigr)+\gamma x^TBy\right].             \tag{1}
\]
All \(n^2\) entries of \(B\), including its diagonal, are free signs.
This is a log-cosh branch, not an order-\(2n\) half-product pressure.

Let \(\mu_{A,s}(x)=2^{-n}e^{sQ_A(x)}/Z_A(s)\), and set
\[
\nu_\eta=\tfrac12\left(
 \mu_{A,\eta}\otimes\mu_{A,-\eta}
 +\mu_{A,-\eta}\otimes\mu_{A,\eta}\right),\qquad v=y\otimes x.
\]
The change \(x\mapsto-x\) leaves both internal energies unchanged and
reverses \(x^TBy\). Consequently, for every fixed \(B\),
\[
F_{A,B}(t)=2a_A(\eta)+
 \log\mathbb E_{\nu_\eta}e^{\gamma v^T\operatorname{vec}B}.    \tag{2}
\]
For any probability distribution \(\mathcal D\) on cross blocks,
define its annealed certificate by
\[
\mathcal A_{A,\mathcal D}(t)=
 \log\mathbb E_{B\sim\mathcal D}e^{F_{A,B}(t)}.
\]
At least one outcome satisfies
\(F_{A,B}(t)\le\mathcal A_{A,\mathcal D}(t)\). However, this
logarithm of an averaged partition function is not
\(\mathbb E_BF_{A,B}(t)\), nor is it \(\min_BF_{A,B}(t)\).

## 1. An exact finite-step Gram--Schmidt bound

The external input is Harshaw, Sävje, Spielman, and Zhang,
*Balancing Covariates in Randomized Experiments with the Gram--Schmidt
Walk Design*, arXiv:1911.03071v8, generalized
[Theorem 6.6*, Supplement S3.5, PDF p.50](https://arxiv.org/pdf/1911.03071v8#page=50).
For arbitrary input columns of norm at most one, the zero-start algorithm
returns a sign vector \(\sigma\) with
\[
\mathbb E e^{u^TW\sigma}\le e^{\|u\|_2^2/2}
\quad\text{for every }u.
\]
Zero start gives mean zero by the martingale property, Lemma 6.1;
the arbitrary-vector input is stated in Supplement S3.1, PDF pp.32--33.

For any \(G\succ0\) with \(G_{ii}\le1\), use \(W=G^{1/2}\).
Its columns have norm at most one. Substituting
\(u=\gamma G^{-1/2}v\) yields
\[
\mathbb E_\sigma e^{\gamma v^T\sigma}
 \le e^{k v^TG^{-1}v}.
\]
Apply this bound before integrating the actual paired Gibbs measure
in (2). Selecting one sign outcome proves
\[
\boxed{\quad
\min_BF_{A,B}(t)\le\mathcal C_{A,G}(t):=
2a_A(\eta)+\log\mathbb E_{\nu_\eta}e^{k v^TG^{-1}v}.
\quad}                                                        \tag{3}
\]
This is valid for a non-infinitesimal step and an actual integral
cross block. The exponential of the quadratic form remains inside
the spin expectation; replacing it by its mean is not justified.

## 2. Every quadratic Gram--Schmidt proxy has a leading floor

Returning (3) to uniform spins gives the exact expression
\[
\mathcal C_{A,G}(t)=\log\mathbb E_{x,y}
 \cosh\!\left[\eta(Q_A(x)-Q_A(y))\right]e^{k v^TG^{-1}v}.
\]
Since \(\cosh\ge1\), Jensen's inequality and
\(\mathbb E_{x,y}vv^T=I_{n^2}\) give
\[
\mathcal C_{A,G}(t)\ge k\operatorname{tr}G^{-1}.
\]
For each coordinate, Cauchy--Schwarz implies
\(1\le G_{ii}(G^{-1})_{ii}\). Therefore every admissible \(G\)
satisfies \(\operatorname{tr}G^{-1}\ge n^2\), and
\[
\boxed{\quad
\mathcal C_{A,G}(t)\ge\frac{c^2t}{4}\,n.
\quad}                                                        \tag{4}
\]
This is uniform over every host \(A\) and every allowed \(G\), even
when \(G\) depends on \(A,n,c,t\). Optimizing the quadratic proxy
does not remove its floor. Equation (4) is not a lower bound on the
actual Gram--Schmidt annealed distribution, since (3) used an upper
bound on that distribution's moment generating function.

## 3. A lower bound for an actual Gaussian-sign moment generating function

We first prove an elementary estimate for any positive semidefinite
correlation matrix \(\Sigma\), including singular matrices. Let
\(z\sim N(0,\Sigma)\), \(b=\operatorname{sign}z\), and let
\(v\in\{-1,1\}^d\). Put \(w=\Sigma v\) and \(V=v^T\Sigma v\).
Then, for every real \(\gamma\),
\[
\boxed{\quad
\log\mathbb E e^{\gamma v^Tb}
\ge\frac{\gamma^2V}{\pi}
 -\frac{2\gamma^4}{3\pi^2}\sum_{i=1}^d|w_i|^3.
\quad}                                                        \tag{5}
\]

To prove this, exponentially tilt the Gaussian law by
\[
\frac{dP_\theta}{dP}(z)=
 \exp\!\left(\theta v^Tz-\frac{\theta^2V}{2}\right).
\]
This remains a probability density when \(\Sigma\) is singular.
Its mean is \(\theta\Sigma v\), its covariance remains \(\Sigma\),
and its relative entropy is \(\theta^2V/2\). The entropy variational
inequality (equivalently Jensen after changing measure) gives
\[
\log\mathbb E e^{\gamma v^Tb}
\ge\gamma\sum_i v_i\bigl(2\Phi(\theta w_i)-1\bigr)
 -\frac{\theta^2V}{2},
\]
where \(\Phi\) is the standard Gaussian distribution function.
Writing \(a_0=\sqrt{2/\pi}\), the identity
\((2\Phi(s)-1)'=a_0e^{-s^2/2}\) implies, for every real \(s\),
\[
\left|2\Phi(s)-1-a_0s\right|\le\frac{a_0|s|^3}{6}.
\]
Choose \(\theta=a_0\gamma\). The quadratic term is
\(\gamma^2V/\pi\), and the error coefficient is
\(a_0^4\gamma^4/6=2\gamma^4/(3\pi^2)\), proving (5).
No variance-sharp subgaussian approximation or cumulant truncation was
assumed: (5) is a rigorous lower bound for the actual sign law.

## 4. Uniform application to the complete-sign Gaussian construction

Let the extremal eigenvalues of \(A\) be \(a>0\) and \(-b<0\).
For the reviewed Gaussian construction, write
\[
H=A\otimes A-\alpha(A\otimes I+I\otimes A),\qquad
\mu=ab+\alpha(a-b),\qquad \Sigma=I+H/\mu.                     \tag{6}
\]
Actual opposite-temperature covariances give
\(\alpha=(\operatorname{tr}(AU)+\operatorname{tr}(AV))/(2n)
\in[-b/2,a/2]\). The construction has
\(\mu=-\lambda_{\min}(H)>0\), so \(\Sigma\succeq0\) with
diagonal one. The following argument is uniform over the entire
stated \(\alpha\)-interval. Thus the covariances used to choose
\(\alpha\) may be taken at \(\beta\), at \(\eta\), or at any
other positive temperature, while the host is held fixed.

We need two uniform parameter bounds. For every eigenvalue \(\lambda\)
of \(A\), \((a-\lambda)(b+\lambda)\ge0\). Summing, and using
\(\operatorname{tr}A=0\) and \(\operatorname{tr}A^2=n(n-1)\), gives
\(ab\ge n-1\). At the two endpoints of the allowed interval for
\(\alpha\), \(\mu\) equals \(b(a+b)/2\) and \(a(a+b)/2\).
Also, direct maximization of \(|\alpha|/[ab+\alpha(a-b)]\) on
that interval gives
\[
\mu\ge\frac{ab}{2}\ge\frac{n-1}{2},\qquad
\frac{|\alpha|}{\mu}\le\frac1{a+b}
 \le\frac1{2\sqrt{n-1}}.                                      \tag{7}
\]

Let \(\ell_n=\log n\), and consider independent uniform \(x,y\)
satisfying
\[
\max_i|(Ax)_i|,\ \max_i|(Ay)_i|\le3\sqrt{n\ell_n},\qquad
|x^TAx|,\ |y^TAy|\le n\ell_n.                                \tag{8}
\]
Each row field is a sum of \(n-1\) independent signs. Hoeffding's
inequality and a union bound give failure probability at most
\(4n^{-7/2}\) for the two field conditions. Furthermore,
\(\mathbb E(x^TAx)^2=2n(n-1)\); Chebyshev bounds the combined
energy failure probability by \(4/\ell_n^2\). Therefore the event
\(\mathcal T_n\) in (8) satisfies
\[
\Pr(\mathcal T_n)\ge1-\delta_n,\qquad
\delta_n=4\ell_n^{-2}+4n^{-7/2}=o(1),                           \tag{9}
\]
uniformly over \(A\).

For \(v=y\otimes x\), the coordinate formula for \(\Sigma v\)
and (7)--(8) imply, for \(n\ge3\),
\[
\|\Sigma v\|_\infty\le
 1+\frac{18n\ell_n}{n-1}
 +3\sqrt{\frac{n\ell_n}{n-1}}
 \le M_n:=42\ell_n.                                           \tag{10}
\]
Writing \(P_x=x^TAx\), \(P_y=y^TAy\), the quadratic form is exactly
\[
V=v^T\Sigma v
=n^2+\frac{P_xP_y-\alpha n(P_x+P_y)}{\mu}.
\]
In particular, throughout \(\mathcal T_n\),
\[
|V-n^2|\le D_n:=n^2\left(
 \frac{2\ell_n^2}{n-1}+\frac{\ell_n}{\sqrt{n-1}}\right)=o(n^2).
                                                                    \tag{11}
\]
Applying (5), with \(d=n^2\), gives the uniform bound
\[
\log\mathbb E_B e^{\gamma x^TBy}
\ge \frac{\gamma^2(n^2-D_n)}{\pi}
 -\frac{2\gamma^4n^2M_n^3}{3\pi^2}
\quad ((x,y)\in\mathcal T_n).                                 \tag{12}
\]
Here \(B=\operatorname{sign}z\) reshaped as a matrix, with precisely
the Gaussian law (6); it is not a substitute covariance proxy.

## 5. The actual Gaussian annealed floor and its consequence

The Gaussian-sign law is invariant under \(B\mapsto-B\). Thus,
for each fixed pair \((x,y)\),
\[
\mathbb E_B\cosh\!\left[\eta(Q_A(x)-Q_A(y))+\gamma x^TBy\right]
=\cosh\!\left[\eta(Q_A(x)-Q_A(y))\right]
 \mathbb E_Be^{\gamma x^TBy}
\ge\mathbb E_Be^{\gamma x^TBy}.
\]
Integrate (12) over \(\mathcal T_n\). For all sufficiently large
\(n\), so that \(\delta_n<1\), the actual annealed certificate obeys
\[
\mathcal A_{A,\mathrm{Gauss}}(t)
\ge\frac{\gamma^2(n^2-D_n)}{\pi}
 -\frac{2\gamma^4n^2M_n^3}{3\pi^2}+\log(1-\delta_n).
                                                                    \tag{13}
\]
The typical set has probability tending to one; no exponentially small
single-spin-state estimate is used. At the stated critical scale,
\(\gamma^2D_n=O_c(\sqrt n\log n+\log^2n)\) and
\(\gamma^4n^2M_n^3=O_c(\log^3n)\). Hence, uniformly over all
complete hosts and all allowed \(\alpha\),
\[
\boxed{\quad
\mathcal A_{A,\mathrm{Gauss}}(t)
\ge\frac{c^2t}{2\pi}\,n
 -O_c(\sqrt n\log n+\log^3n).
\quad}                                                        \tag{14}
\]
This requires no bound on \(\|\Sigma\|_{\rm op}\).

For comparison with the active order-\(n\) endpoint, the reviewed
all-orders upper bound in `CORE.md` is
\(m_n\le(1/2+o(1))n^{3/2}\). Since \(a_A(\beta)\le\beta\Phi(A)\),
\[
2R_n(\beta)\le2\beta m_n\le c n+o_c(n).                       \tag{15}
\]
The paired optimized path starts at \(2R_n(\beta)\). Combining
(4), (14), and (15) yields two different obstructions:
\[
\inf_{A,G}\mathcal C_{A,G}(t)-2R_n(\beta)
\ge\left(\frac{c^2t}{4}-c-o(1)\right)n,                       \tag{16}
\]
\[
\inf_A\mathcal A_{A,\mathrm{Gauss}}(t)-2R_n(\beta)
\ge\left(\frac{c^2t}{2\pi}-c-o(1)\right)n.                   \tag{17}
\]
The second infimum may also range over any admissible choice of
\(\alpha\) for each host. Thus the quadratic Gram--Schmidt proxy
cannot supply the desired finite-step small-oh comparison when
\(c>4/t\); direct annealing of the actual Gaussian-sign construction
cannot supply it when \(c>2\pi/t\).

These statements do not invalidate the reviewed endpoint second-moment
saving. They do not rule out selecting a good Gaussian-sign outcome,
controlling the average of \(F_{A,B}\) instead of the logarithm of the
average partition function, using another signing law, or giving another
order-comparison proof. The actual Gram--Schmidt law itself is not
excluded by (16). All asymptotic bounds here keep \(c,t\) fixed; no
claim is made for growing temperatures or vanishing step sizes.
