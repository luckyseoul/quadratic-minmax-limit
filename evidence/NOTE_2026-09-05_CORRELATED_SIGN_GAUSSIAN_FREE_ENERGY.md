# Correlated Gaussian signs have the matched-Gaussian quenched free energy

2026-09-05. **Analytic theorem; independently complete-read and verified.**
The comparison retains the entire spin posterior. It is valid for singular
Gaussian covariances and needs only a uniform covariance operator bound,
not small individual correlations or independent coefficient refills.
Its application to an order comparison still requires the separate
matched-Gaussian pressure estimate.

## 1. Statement

Let G be a centered Gaussian vector in R^m with correlation matrix Sigma,
where

\[
 \operatorname{diag}\Sigma=1,\qquad
 0\preceq\Sigma\preceq K I,
 \qquad K\ge1.
\]

Let V be any finite nonempty subset of [-1,1]^m and let pi be any fixed
probability measure on V. Zero weights can be discarded. Define

\[
 \mathcal F(b)=\log\sum_{v\in V}\pi(v)e^{\gamma v^Tb},
 \qquad L=\log(2|V|),\qquad
 C=\mathbb E[\operatorname{sign}(G)\operatorname{sign}(G)^T]
   =\frac2\pi\arcsin^{\circ}\Sigma.
                                                               \tag{1}
\]

The superscript denotes entrywise arcsine, including the diagonal.
Let Y be centered Gaussian with covariance C. For every 0<epsilon<=1,

\[
 \begin{aligned}
 |\mathbb E\mathcal F(\operatorname{sign}G)
                         -\mathbb E\mathcal F(Y)|
 \le{}&2|\gamma|\sqrt{K m\epsilon L/\pi}
       +\frac{K|\gamma|}{3}L\\
 &+4(2/\pi)^{3/2}e^{-1/2}K^2|\gamma|^3m\epsilon^{-4}
       +\frac{\sqrt2}{\pi}K\gamma^2m\epsilon.
 \end{aligned}                                                \tag{2}
\]

In particular, if m=n^2, |gamma|<=C_0/sqrt(n), and L<=C_1 n, choosing
`epsilon=n^(-1/9)` gives

\[
 \boxed{\quad
 |\mathbb E\mathcal F(\operatorname{sign}G)
                         -\mathbb E\mathcal F(Y)|
 =O_{K,C_0,C_1}(n^{17/18})=o(n).
 \quad}                                                       \tag{3}
\]

All constants are uniform in the fixed prior pi. Thus pi may incorporate
an arbitrary deterministic internal Hamiltonian depending on n.
It must not depend on the Gaussian disorder G or Y.

## 2. Smooth activation comparison with the full posterior

First let f be a bounded smooth centered function with bounded first
and second derivatives, and assume that the covariance decrease in the
following construction is positive semidefinite. This last property
holds for every Gaussian OU smoothing, as verified below for the only
activation needed here.

Write P_t for the one-dimensional Ornstein--Uhlenbeck semigroup,

\[
 f_t(x)=P_tf(x)=\mathbb E f(e^{-t}x+\sqrt{1-e^{-2t}}Z),
 \qquad Z\sim N(0,1),
\]

and put `C_t=Cov(f_t(G))`, `C_f=C_0` within this section. Let Z_t be a
Gaussian vector independent of G with covariance C_f-C_t, and define

\[
 \Psi(t)=\mathbb E\mathcal F(f_t(G)+Z_t).
\]

For each t, only this joint distribution is needed; no pathwise
differentiation of a covariance square root is used. Gaussian integration
by parts for singular G holds directly by representing G as a linear
image of independent standard Gaussians. For covariance differentiation
of Z_t one may first use covariance C_f-C_t+delta I and then send
delta down to zero. This regularizes only the independent Z_t, leaving
Sigma and its unit diagonals unchanged. The bounded derivatives of
mathcal F justify the limit.

The generator identity `partial_t f_t=f_t''-x f_t'` and Gaussian
integration by parts give

\[
 (C_t')_{ij}=-2\Sigma_{ij}\,
                 \mathbb E[f_t'(G_i)f_t'(G_j)].             \tag{4}
\]

Applying the Gaussian OU generator to
`mathcal F(f_t(G)+z)` at fixed z, then accounting for the covariance
derivative of the independent Z_t, gives the EXACT identity

\[
 \Psi'(t)=-\sum_{i,j}\Sigma_{ij}\,
   \operatorname{Cov}\!\left(
       f_t'(G_i)f_t'(G_j),\,
       \mathcal F_{ij}(f_t(G)+Z_t)\right).                 \tag{5}
\]

The covariance in (5) includes G and Z_t, but its first variable depends
only on G. In particular it may be evaluated by first conditioning on
the same fixed value of Z_t in its two factors. No Gibbs derivative is
replaced by its prior value.

For completeness the Gaussian covariance identity used next is

\[
 \operatorname{Cov}(U(G),W(G))
 =\int_0^\infty e^{-s}
   \mathbb E[\nabla U(G)^T\Sigma\nabla W(G_s)]\,ds,
 \quad G_s=e^{-s}G+\sqrt{1-e^{-2s}}G',                       \tag{6}
\]

where G' is an independent copy. It follows by differentiating
`E[U(G)W(G_s)]` in s, or by Gaussian integration by parts and the OU
semigroup. The same identity holds conditionally on an independent
Z_t, which is then held FIXED while forming G_s.

Apply (6) to (5). Differentiating the two factors
`f_t'(G_i)f_t'(G_j)` produces two terms. They are equal after exchanging
i,j, using symmetry of Sigma and of the third derivative. Hence

\[
 \begin{aligned}
 \Psi'(t)=-2\int_0^\infty e^{-s}\mathbb E
 \sum_{i,j,l}\Sigma_{ij}\Sigma_{il}
 f_t''(G_i)f_t'(G_j)f_t'((G_s)_l)\,
 \mathcal F_{ijl}(f_t(G_s)+Z_t)\,ds.
 \end{aligned}                                                \tag{7}
\]

This formula explicitly keeps the mixed-time diagonal factors:
`f_t''(G_i), f_t'(G_j)` occur at G, whereas `f_t'((G_s)_l)` and the
ENTIRE posterior occur at G_s. The independent Gaussian Z_t is the
same in that posterior as in the conditional application of (6).

Fix G,G_s,Z_t in the integrand. Let angle brackets denote the actual
posterior on V at `b=f_t(G_s)+Z_t`, let `bar v=<v>`, and set
`u=v-bar v`. The exact third log-partition derivative is

\[
 \mathcal F_{ijl}(b)=\gamma^3\langle u_i u_j u_l\rangle.
                                                               \tag{8}
\]

Thus it is one third CENTRAL moment under one posterior. There is no
need to decouple replicas or move that posterior outside the Gaussian
expectation. Set

\[
 D=\operatorname{diag}(f_t'(G_j)),\qquad
 \widetilde D=\operatorname{diag}(f_t'((G_s)_j)),\qquad
 a_i=f_t''(G_i).
\]

The full contraction in (7), divided by gamma^3, equals exactly

\[
 \left\langle
  \sum_i a_i u_i(\Sigma D u)_i(\Sigma\widetilde D u)_i
 \right\rangle.                                             \tag{9}
\]

Since |u_i|<=2 and ||u||_2<=2sqrt(m), Cauchy--Schwarz and the covariance
operator bound yield, for EACH posterior configuration,

\[
 \begin{aligned}
 \left|\sum_i a_i u_i(\Sigma D u)_i
                          (\Sigma\widetilde D u)_i\right|
 &\le2\|a\|_\infty\|\Sigma Du\|_2
                          \|\Sigma\widetilde Du\|_2\\
 &\le8mK^2\|f_t''\|_\infty\|f_t'\|_\infty^2.
 \end{aligned}                                                \tag{10}
\]

The signed matrix contractions are retained until (10); an entrywise
absolute bound on the two factors Sigma would lose this estimate.

The semigroup derivative bounds are
`||f_t'||_infty<=e^(-t)||f'||_infty` and
`||f_t''||_infty<=e^(-2t)||f''||_infty`. Consequently (7)--(10) give

\[
 |\Psi'(t)|\le16|\gamma|^3mK^2e^{-4t}
                     \|f''\|_\infty\|f'\|_\infty^2.
\]

At t=0 the vector is f(G); as t tends to infinity it converges in law
to a centered Gaussian Y_f with covariance C_f. Boundedness of f and
the Lipschitz bound for mathcal F justify convergence of expectations.
Integration therefore proves

\[
 \boxed{\quad
 |\mathbb E\mathcal F(f(G))-\mathbb E\mathcal F(Y_f)|
 \le4|\gamma|^3mK^2\|f''\|_\infty\|f'\|_\infty^2.
 \quad}                                                       \tag{11}
\]

## 3. A self-contained Gaussian Holder inequality

Under `Sigma<=K I`, nonnegative one-coordinate functions satisfy

\[
 \boxed{\quad
 \mathbb E\prod_{i=1}^m q_i(G_i)
 \le\prod_{i=1}^m(\mathbb E q_i(Z)^K)^{1/K},
 \qquad Z\sim N(0,1).
 \quad}                                                       \tag{12}
\]

Here is a direct proof for bounded positive smooth functions; bounded
measurable functions follow by approximation. Let
`g_i(x,t)=P_t(q_i^K)(x)`, `a_i=partial_x log g_i`, and
`U_t(G)=prod_i g_i(G_i,t)^(1/K)`. The scalar semigroup equation and
Gaussian integration by parts give

\[
 \frac{d}{dt}\mathbb E U_t
 =\mathbb E\left[
  U_t\left(\frac1K\sum_i a_i^2
              -\frac1{K^2}a^T\Sigma a\right)\right]\ge0.
                                                               \tag{13}
\]

The initial value is the left side of (12), and the limiting value is
its right side. This proves (12), also for singular Sigma. Approximation
can be made in each univariate standard-Gaussian marginal, so singular
joint support introduces no additional discontinuity issue.

## 4. Controlling the sign-smoothing error at the correct scale

Choose

\[
 f_\epsilon(x)=2\Phi(x/\epsilon)-1,
 \qquad r_\epsilon(x)=\operatorname{sign}(x)-f_\epsilon(x).
\]

Both functions are odd, and `|r_epsilon|<=1`. Under a standard normal Z,
write `v_epsilon=E r_epsilon(Z)^2`. The sign-disagreement formula for
`Z` and `Z+epsilon Z'` gives

\[
 v_\epsilon\le\mathbb E|r_\epsilon(Z)|
 =\frac2\pi\arctan\epsilon\le\frac{2\epsilon}{\pi}.       \tag{14}
\]

For a centered random variable R bounded in absolute value by one,
with variance at most v, the elementary exponential series gives

\[
 \log\mathbb E e^{aR}
 \le\frac{v a^2}{2(1-|a|/3)},\qquad |a|<3.                \tag{15}
\]

Indeed, `E|R|^k<=v` for k>=2 and `k!>=2*3^(k-2)` suffice.
Apply (12) to exponential functions of r_epsilon, then (15). For any
v in V and every `s>=0` with `K s|gamma|<3`,

\[
 \log\mathbb E\exp\left(s\gamma\sum_i v_i r_\epsilon(G_i)\right)
 \le\frac{K s^2\gamma^2m v_\epsilon}
                   {2(1-Ks|\gamma|/3)}.                    \tag{16}
\]

The same bound holds for its negative. The exponential maximum bound
over these 2|V| variables, optimized in s, gives

\[
 \mathbb E\max_{v\in V}\left|\gamma\sum_i v_i r_\epsilon(G_i)\right|
 \le |\gamma|\sqrt{2K m v_\epsilon L}
                          +\frac{K|\gamma|}{3}L.           \tag{17}
\]

For example substitute
`s=sqrt(2L/a)/(1+b sqrt(2L/a))` with
`a=K gamma^2 m v_epsilon`, `b=K|gamma|/3` in
`L/s+a s/[2(1-bs)]`. Zero-variance cases follow by a limit.

Pointwise log-sum-exp comparison bounds
`|mathcal F(sign G)-mathcal F(f_epsilon(G))|` by the maximum in (17).
Combining (14)--(17) proves

\[
 |\mathbb E\mathcal F(\operatorname{sign}G)
                  -\mathbb E\mathcal F(f_\epsilon(G))|
 \le2|\gamma|\sqrt{K m\epsilon L/\pi}
                             +\frac{K|\gamma|}{3}L.       \tag{18}
\]

This is the step that permits a slowly vanishing epsilon. A sum of
coordinatewise absolute errors would instead give an inadequate
`O(|gamma| m epsilon)` bound.

## 5. Smoothing covariance and matched-Gaussian comparison

The function f_epsilon is OU-smoothed sign with contraction
`q=(1+epsilon^2)^(-1/2)`. Therefore

\[
 C_\epsilon=\operatorname{Cov}(f_\epsilon(G))
 =\frac2\pi\arcsin^\circ(q^2\Sigma).
                                                               \tag{19}
\]

Its further OU smoothing has covariance
`(2/pi) arcsin^circ(q^2 e^(-2t) Sigma)`. The arcsine power series has
positive coefficients and odd powers. Every Schur power of Sigma is
positive semidefinite, so `C_epsilon-C_t>=0` as required in Section 2.
Also every such Schur power has operator norm at most K: the Schur
map for a correlation matrix is positive and sends I to I, hence
contracts the operator norm on real symmetric matrices.

It follows that D_epsilon=C-C_epsilon is positive semidefinite and

\[
 \begin{aligned}
 \|D_\epsilon\|_{\rm op}
 &\le K\left[1-\frac2\pi\arcsin\frac1{1+\epsilon^2}\right]\\
 &=\frac{4K}{\pi}\arctan
                     \frac\epsilon{\sqrt{2+\epsilon^2}}
 \le\frac{2\sqrt2K}{\pi}\epsilon.
 \end{aligned}                                                \tag{20}
\]

For independent centered Gaussians Y_epsilon with covariance C_epsilon
and W with covariance D_epsilon, their sum has covariance C. Gaussian
covariance interpolation and the actual posterior Hessian give

\[
 0\le\mathbb E\mathcal F(Y)-\mathbb E\mathcal F(Y_\epsilon)
 \le\frac{\gamma^2m}{2}\|D_\epsilon\|_{\rm op}
 \le\frac{\sqrt2K}{\pi}\gamma^2m\epsilon.                 \tag{21}
\]

Here `mathcal F''=gamma^2 Cov_posterior(v)` is positive semidefinite
and has trace at most gamma^2 m. Thus (21) is a quenched comparison,
not an annealed moment-generating estimate.

Finally,

\[
 \|f_\epsilon'\|_\infty=\sqrt{2/\pi}\,\epsilon^{-1},
 \qquad
 \|f_\epsilon''\|_\infty
       =\sqrt{2/\pi}\,e^{-1/2}\epsilon^{-2}.
\]

Use these values in (11), then combine (18), (11), and (21).
This is exactly (2), and the stated choice of epsilon proves (3).

## 6. The actual bipartite finite-step pressure

For any fixed host A and any fixed deterministic internal term
`I(x,y)=eta[Q_A(x)-Q_A(y)]`, set

\[
 \mathscr F(b)=\log\mathbb E_{x,y}
                \cosh[I(x,y)+\gamma x^T\operatorname{mat}(b)y].
\]

Augment by sigma=+/-1. After subtracting the fixed internal partition
constant, this is exactly (1) with

\[
 v=\sigma(y\otimes x)\in\{-1,1\}^{n^2},\qquad
 \pi(\sigma,x,y)\propto e^{\sigma I(x,y)},
 \qquad |V|\le2^{2n+1}.
\]

The prior is fixed independently of the cross disorder. Consequently
(3) applies uniformly to this ACTUAL pressure whenever
`||Sigma||op<=K` and `gamma=O(n^(-1/2))`, with no bound needed on I.
The covariance may be the full-strength singular canonical covariance;
only its stated operator bound is used here.

The conclusion is a sign-to-MATCHED-GAUSSIAN quenched equivalence.
Any subsequent replacement of C by a simpler Gaussian covariance, or
any comparison of that Gaussian pressure with the optimized smaller-
order endpoint, is a separate step. No independence of the original
correlated signs, frozen-posterior substitution, or favorable endpoint
pressure inequality has been assumed.
