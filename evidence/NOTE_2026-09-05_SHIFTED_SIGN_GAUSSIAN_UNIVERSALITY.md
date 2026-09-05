# Mean-preserving shifted signs have the matched-Gaussian quenched pressure

2026-09-05. Analytic theorem; independent complete proof reads passed.
This extends the reviewed zero-threshold comparison to arbitrary fixed
threshold size, uniformly in that size. It retains the actual posterior.
It does not provide an upper bound comparing the optima at two orders.

## 1. Statement, including arbitrary deterministic threshold orientations

Let G be a centered Gaussian vector in R^m with correlation matrix Sigma,
where diag(Sigma)=1 and 0<=Sigma<=K I, K>=1. Fix signs d_i in {-1,1}
and any real h. Define

\[
 B_i=\operatorname{sign}(G_i+h d_i),\qquad
 s_h=2\Phi(h)-1,\qquad b_i=s_h d_i.
\]

Here Phi is the standard normal distribution function. Let Y be Gaussian
with mean b and covariance EXACTLY Cov(B). Singular covariance is allowed.
For a finite nonempty V subset of [-1,1]^m, a deterministic probability
measure pi on V, and real gamma, put

\[
 F(z)=\log\sum_{v\in V}\pi(v)e^{\gamma v^Tz},
 \qquad L=\log(2|V|),\qquad \kappa=2/\pi.
\]

The prior, orientations and threshold may depend on m and gamma, but
are fixed independently of the Gaussian disorder. Zero prior weights
can be discarded. For every 0<epsilon<=1,

\[
\begin{split}
 |\mathbb EF(B)-\mathbb EF(Y)|\le{}&
 2^{5/4}|\gamma|\sqrt{K m\epsilon L/\pi}
       +{2K|\gamma|\over3}L\\
 &+4\kappa^{3/2}e^{-1/2}K^2|\gamma|^3m\epsilon^{-4}
       +{\sqrt2 K\over\pi}\gamma^2m\epsilon.
                                                               \tag{1}
\end{split}
\]

Every constant is uniform in h, d and pi. In particular, m<=n^2,
|gamma|<=C_0/sqrt(n), L<=C_1 n and epsilon=n^(-1/9) give an
O_{K,C_0,C_1}(n^(17/18))=o(n) quenched pressure error.

The proof uses the smooth-activation estimate and Gaussian Holder
inequality proved in Sections 2--3 of
`NOTE_2026-09-05_CORRELATED_SIGN_GAUSSIAN_FREE_ENERGY.md`.
Their hypotheses are verified below; oddness is not one of the
hypotheses of the smooth-activation estimate.

## 2. Reduce to one centered scalar activation without a mean error

Let D=diag(d_i), G'=DG and R=D Sigma D. Then R is a correlation
matrix with 0<=R<=K I, and

\[
 B=D\{s_h\mathbf1+f_h(G')\},\qquad
 f_h(x)=\operatorname{sign}(x+h)-s_h.
                                                               \tag{2}
\]

The identity is almost sure because each Gaussian coordinate has
variance one. Its deterministic mean can be absorbed EXACTLY into the
prior: replace each observable v by Dv and its prior weight by a weight
proportional to pi(v) exp(gamma s_h d^T v). The same deterministic
normalizing constant is added to both pressures being compared.
Therefore it suffices to prove (1) for the centered vector f_h(G')
and its centered matched Gaussian, with an arbitrary fixed prior.
All transformed observables still lie in [-1,1]^m.

Let P_t be the one-dimensional Gaussian Ornstein--Uhlenbeck semigroup,
and choose

\[
 t_\epsilon={1\over2}\log(1+\epsilon^2),\qquad
 f_{h,\epsilon}=P_{t_\epsilon}f_h.
\]

The smoothing formula is

\[
 f_{h,\epsilon}(x)
 =2\Phi\left({x+h\sqrt{1+\epsilon^2}\over\epsilon}\right)-1-s_h.
                                                               \tag{3}
\]

Gaussian invariance of P_t gives E f_{h,epsilon}(Z)=0 exactly.
There is no threshold-dependent O(epsilon^2) mean discarded here.
The derivative norms, uniformly in h, are

\[
 \|f_{h,\epsilon}'\|_\infty=\sqrt\kappa\,\epsilon^{-1},
 \qquad
 \|f_{h,\epsilon}''\|_\infty
       =\sqrt\kappa\,e^{-1/2}\epsilon^{-2}.
                                                               \tag{4}
\]

## 3. Uniform threshold disagreement and covariance decrease

For a centered standard Gaussian pair (X,Y) of correlation q in [0,1],
let D_a(q)=P(sign(X-a) != sign(Y-a)). For q<1, differentiation of
the two-variable Gaussian distribution gives

\[
 {\partial\over\partial a}D_a(q)
 =2\phi(a)\left[1-2\Phi\left(a\sqrt{1-q\over1+q}\right)\right].
                                                               \tag{5}
\]

Indeed D_a(q)=2[Phi(a)-Phi_2(a,a;q)], and each of the two boundary
derivatives of Phi_2 is phi(a) Phi(a sqrt((1-q)/(1+q))). Equation (5)
is nonnegative when a<=0 and nonpositive when a>=0. Thus disagreement
is maximal at threshold zero. The planar Gaussian sign identity gives

\[
 D_a(q)\le D_0(q)={\arccos q\over\pi}.                    \tag{6}
\]

The endpoint q=1 follows by continuity, or directly. This argument
does not assume that the shifted sign has mean zero.

Expand the centered activation in the normalized Hermite basis:

\[
 f_h(x)=\sum_{j\ge1}a_j(h){\operatorname{He}_j(x)\over\sqrt{j!}},
 \qquad \sum_{j\ge1}a_j(h)^2=1-s_h^2\le1.
\]

The joint Gaussian generating function proves the Hermite covariance
identity and hence

\[
 C_h:=\operatorname{Cov}(f_h(G'))
       =\sum_{j\ge1}a_j(h)^2 R^{\circ j},\qquad
 C_{h,\epsilon}
       =\sum_{j\ge1}a_j(h)^2 e^{-2jt_\epsilon}R^{\circ j}.
                                                               \tag{7}
\]

All integer Schur powers appear; the even ones are not omitted.
Each is positive semidefinite and has operator norm at most K.
For the latter assertion, R<=K I and the positive Schur map for the
correlation matrix R^(circ(j-1)) give R^(circ j)<=K I for j>=2;
the case j=1 is the hypothesis. The series converge in operator norm
by the displayed square-summability. Consequently

\[
 0\preceq C_h-C_{h,\epsilon}\preceq K v_h(\epsilon)I,
 \qquad
 v_h(\epsilon)=\sum_{j\ge1}a_j(h)^2(1-e^{-2jt_\epsilon}).
                                                               \tag{8}
\]

For a Gaussian pair of correlation e^(-2t_epsilon), subtracting the
constant mean from both signs cancels in the variance decrease, so

\[
 \begin{split}
 v_h(\epsilon)
 &=2D_{-h}(e^{-2t_\epsilon})\\
 &\le {2\over\pi}\arccos{1\over1+\epsilon^2}
   ={4\over\pi}\arctan{\epsilon\over\sqrt{2+\epsilon^2}}
   \le {2\sqrt2\over\pi}\epsilon.
                                                               \tag{9}
 \end{split}
\]

The same Hermite expansion proves that C_{h,epsilon} minus the
covariance of any further OU smoothing is positive semidefinite.
This verifies the covariance-decrease hypothesis required by the
reviewed smooth-activation comparison.

## 4. Smoothing error at the finite-state quenched scale

Set r=f_h-P_(t_epsilon) f_h. Gaussian invariance gives E r(Z)=0,
and |r|<=2. Parseval and (9) give

\[
 \mathbb E r(Z)^2
 =\sum_{j\ge1}a_j(h)^2(1-e^{-jt_\epsilon})^2
 \le v_h(\epsilon)\le {2\sqrt2\over\pi}\epsilon.
                                                               \tag{10}
\]

The absolute range two is retained; the zero-threshold range-one
argument is not applicable at a shifted threshold.

For a centered variable U with |U|<=2 and variance at most v,
the exponential series and k!>=2*3^(k-2), k>=2, yield

\[
 \log\mathbb Ee^{zU}
 \le {v z^2\over2(1-2|z|/3)},\qquad |z|<3/2.
                                                               \tag{11}
\]

Apply the reviewed Gaussian Holder inequality for R<=K I to the
coordinate exponentials of r(G'_i). For any observable v in V,
the resulting Bernstein bound for gamma sum_i v_i r(G'_i) has
variance proxy K gamma^2 m v_h(epsilon) and scale 2K|gamma|/3.
Applying the exponential maximum bound to both signs and all |V|
observables, and optimizing its parameter, gives

\[
 \mathbb E\max_{v\in V}
       \left|\gamma\sum_i v_i r(G'_i)\right|
 \le |\gamma|\sqrt{2K m v_h(\epsilon)L}
                         +{2K|\gamma|\over3}L.
                                                               \tag{12}
\]

For example, this follows by inserting
z=sqrt(2L/A)/(1+B sqrt(2L/A)) into
L/z+A z/[2(1-Bz)], where A=K gamma^2 m v_h(epsilon) and
B=2K|gamma|/3. Zero variance or gamma=0 follows by continuity.
Pointwise log-sum-exp comparison bounds the pressure difference by
the maximum in (12). Equations (9) and (12) give the first two
terms of (1). No sum of coordinatewise absolute errors is used.

## 5. Smooth comparison and restoration of the matched covariance

For a bounded smooth centered activation f with the covariance
decrease verified in Section 3, the previously reviewed smooth
comparison is

\[
 |\mathbb EF(f(G'))-\mathbb EF(Y_f)|
 \le4|\gamma|^3mK^2\|f''\|_\infty\|f'\|_\infty^2,
                                                               \tag{13}
\]

where Y_f has covariance Cov(f(G')). That proof differentiates the
actual log partition along further OU smoothing plus independent
matched Gaussian covariance. The third derivative is the current
posterior central third moment. Its signed contraction is bounded
before taking absolute values, by
8m K^2 ||f_t''||_infinity ||f_t'||_infinity^2. Integrating the
resulting 16 e^(-4t) factor proves (13). Neither oddness nor a
zero-threshold Hermite parity assumption occurs in that argument.
Singular Gaussian covariance is handled by a linear Gaussian factor;
any auxiliary diagonal regularization is applied only to the
independent matching noise, not to the unit-diagonal latent law.

Apply (13) to f_(h,epsilon) and use (4). This gives the third term
of (1). Finally, (8)--(9) and Gaussian covariance interpolation give

\[
 0\le\mathbb EF(Y_h)-\mathbb EF(Y_{h,\epsilon})
 \le {\gamma^2m\over2}\|C_h-C_{h,\epsilon}\|
 \le {\sqrt2 K\over\pi}\gamma^2m\epsilon.
                                                               \tag{14}
\]

Here the Hessian is gamma^2 Cov_current(v), positive semidefinite
with trace at most gamma^2 m. The arbitrary fixed prior remains
inside this derivative. Combining (12)--(14) proves (1).

## 6. Quantitative expected ORIGINAL norm consequence

Suppose m<=n^2 and the coefficient states, after augmenting by the
sign of an absolute value, number at most 2^(2n+1). Allow arbitrary
fixed deterministic internal energies I_n and a cross coefficient
theta with |theta|<=1. For beta=c/sqrt(n), the associated cosh
pressure has a fixed prior proportional to exp(beta sigma I_n).
Apply (1) with |gamma|<=c/sqrt(n). For a fixed K this yields

\[
 |\mathbb EP_B-\mathbb EP_Y|
 \le D_K[cn\sqrt\epsilon+c\sqrt n
       +c^3\sqrt n\epsilon^{-4}+c^2n\epsilon],             \tag{15}
\]

uniformly over thresholds h, orientations and internal energies.
For each outcome the pressure is between beta times the maximum
absolute energy and that quantity minus (2n+1)log 2. Consequently

\[
 { |\mathbb E M(B)-\mathbb E M(Y)|\over n^{3/2}}
 \le D_K[\sqrt\epsilon+n^{-1/2}
          +c^2n^{-1/2}\epsilon^{-4}+c\epsilon+1/c].        \tag{16}
\]

Take c=n^(1/22) and epsilon=n^(-1/11). Every principal term in
(16) is n^(-1/22), so the expected raw maximum error is

\[
 |\mathbb E M(B)-\mathbb E M(Y)|\le D_K n^{16/11}.
                                                               \tag{17}
\]

The estimate uses the complete temperature-dependent error, not an
unjustified exchange of a fixed-temperature limit with zero temperature.

For a complete source A and its universal midpoint latent covariance,
apply the result on the n(n-1) off-diagonal cross coordinates, with
d_(ij)=A_ij and K=2. Their marginal covariance has operator norm at
most two. The omitted n cross-diagonal entries contribute at most n
for the signs and at most n in expected norm for the matched Gaussian,
regardless of their correlations with the other coefficients.
Absorb this 2n cost into the absolute constant in (17).

Thus (17) applies to the actual law B_(h,ij)=sign(G_ij+h A_ij),
including its zero-threshold diagonal, and the Gaussian Y_h with the
same mean and covariance. For I_n(x,y)=Q_A(x)-Q_A(y) and theta=1,
M is exactly the ORIGINAL paired maximum-absolute quadratic norm.
No source operator bound, regularization or half-product identification
is required for this sign-to-MATCHED-Gaussian step.

The separate even-Hermite covariance reduction is needed before
replacing this matched Gaussian by a simple mean-plus-noise model.
An upper comparison of that model against a smaller-order original
minimum remains an additional unresolved task. The original MO
convergence question is not settled by this theorem.
