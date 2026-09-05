# Universal spectral-midpoint covariance for the original paired norm

Date: 2026-09-05.  Status: exact finite-dimensional covariance theorem;
the expected paired Gaussian-norm upper comparison is not proved here.

This note concerns a freely chosen Gaussian rounding law.  Its parameter
need not equal a parameter obtained from the original Gibbs profile.

## 1. Source, full parameter interval, and exact factorization

Let `n >= 2` and let `A` be a real symmetric `n x n` matrix with zero
diagonal and off-diagonal entries in `{ -1, 1 }`.  Write

\[
 a=\lambda_{\max}(A)>0,\qquad -b=\lambda_{\min}(A)<0,
 \qquad L=\|A\|_{\rm op}=\max(a,b).
\]

The trace identity and the positive semidefinite matrix
`(aI-A)(A+bI)` give

\[
 0\le\operatorname{tr}((aI-A)(A+bI))
   =n\{ab-(n-1)\},\qquad ab\ge n-1.                 \tag{1}
\]

For a real parameter `alpha`, set

\[
 H_\alpha=A\otimes A-\alpha(A\otimes I+I\otimes A),
 \qquad \mu_\alpha=ab+\alpha(a-b).                  \tag{2}
\]

For every `alpha in [-b,a]`, `mu_alpha > 0` and

\[
 \Sigma_\alpha=I+H_\alpha/\mu_\alpha                \tag{3}
\]

is a correlation matrix.  Indeed, put

\[
 P={A+bI\over a+b},\qquad Q={aI-A\over a+b}.
\]

Then `P,Q` are positive semidefinite, `P+Q=I`, and

\[
 \boxed{\displaystyle
 \Sigma_\alpha=u_\alpha P\otimes P+v_\alpha Q\otimes Q,
 \quad u_\alpha={(a+b)(a-\alpha)\over\mu_\alpha},
 \quad v_\alpha={(a+b)(b+\alpha)\over\mu_\alpha}.}   \tag{4}
\]

Both coefficients are nonnegative, and the diagonal of (3) is one.
Also `mu_alpha` lies between its positive endpoint values `b^2,a^2`.
Identity (4) follows by expansion: the coefficients of `A tensor A`,
`A tensor I + I tensor A`, and `I` are respectively
`1/mu_alpha`, `-alpha/mu_alpha`, and `1`.

The interval `[-b,a]` is the full admissible interval for formula (3).
For `mu_alpha > 0`, evaluating on the extreme eigenvector pairs
`(a,a)` and `(-b,-b)` requires respectively `alpha <= a` and
`alpha >= -b`.  If `mu_alpha < 0`, those same two requirements reverse
and become incompatible.  Formula (3) is undefined when `mu_alpha=0`.

Equivalently, the scalar eigenvalue of `H_alpha` at eigenvalues `s,t`
of `A` is `st-alpha(s+t)`.  It is bilinear on the rectangle
`[-b,a]^2`; its two mixed corners equal `-mu_alpha`, and the other
two corners plus `mu_alpha` are

\[
 (a+b)(a-\alpha),\qquad (a+b)(b+\alpha).
\]

Thus `lambda_min(H_alpha)=-mu_alpha` throughout the interval.
Moreover,

\[
 \|\Sigma_\alpha\|_{\rm op}=\max(u_\alpha,v_\alpha).\tag{5}
\]

The corners attaining `u_alpha,v_alpha` occur in the actual spectrum;
bilinearity shows that no interior eigenvalue is larger.

At the endpoints, (4) reduces to the single-sector correlations

\[
 \Sigma_{-b}=(I+A/b)\otimes(I+A/b),\qquad
 \Sigma_a=(I-A/a)\otimes(I-A/a).                    \tag{6}
\]

## 2. The operator-optimal midpoint has a universal bound

Choose

\[
 \boxed{\displaystyle
 \alpha_*={a-b\over2},\qquad
 \mu_*={a^2+b^2\over2},\qquad
 \Sigma_*=I+H_{\alpha_*}/\mu_*.}                    \tag{7}
\]

In fact `alpha_*` belongs to the smaller interval `[-b/2,a/2]`.
The two same-corner eigenvalues of `H_alpha*` are both `ab`;
the mixed corners are `-mu_*`.  Hence

\[
 \boxed{\displaystyle
 0\preceq\Sigma_*\preceq
 { (a+b)^2\over a^2+b^2}I\preceq2I,
 \qquad \operatorname{diag}\Sigma_*=1.}             \tag{8}
\]

No bound on `L/sqrt(n)` is assumed in (8).

This parameter uniquely minimizes (5) over the full interval:

\[
 u_\alpha'=-{a^2(a+b)\over\mu_\alpha^2}<0,
 \qquad v_\alpha'={b^2(a+b)\over\mu_\alpha^2}>0,     \tag{9}
\]

and the two functions meet exactly at `alpha_*`.
This minimizes the covariance operator norm, not necessarily the
expected paired norm or its pressure.

One alternative form makes positivity transparent.  If

\[
 d={a-b\over a+b},\qquad
 T={2A-(a-b)I\over a+b},
\]

then `-I <= T <= I`, `diag(T)=-d`, and

\[
 \Sigma_*={I+T\otimes T\over1+d^2}.                 \tag{10}
\]

## 3. Uniform exact arcsine linearization

Fix any `rho in [0,1]`.  Let

\[
 \Sigma_\rho=(1-\rho)I+\rho\Sigma_*,\qquad
 G\sim N(0,\Sigma_\rho),\qquad B=\operatorname{sign}(G),
 \qquad \kappa=2/\pi.
\]

The coordinates are indexed by ordered pairs `(i,j)`, so `B` is an
`n x n` cross block.  All marginal variances are one, even when the
joint covariance is singular.  Thus the convention for `sign(0)` is
irrelevant.  Let `C=Cov(B)`.  The Gaussian arcsine identity gives

\[
 C=\kappa\Sigma_\rho+(1-\kappa)I+E_\rho,             \tag{11}
\]

where, writing `r(t)=arcsin(t)-t`,

\[
 E_\rho=\kappa\left[
 r(\rho/\mu_*)A\otimes A
 -r(\rho\alpha_*/\mu_*)
       (A\otimes I+I\otimes A)\right].             \tag{12}
\]

To check (12), distinct cross coordinates with both indices different
have covariance `rho A_ik A_jl/mu_*`; those sharing exactly one index
have covariance `-rho alpha_* A_ik/mu_*` (or the analogous other
index).  The diagonal is exact in (11).

For `|t| <= 1`, the positive Taylor coefficients of the odd arcsine
series imply

\[
 \kappa|r(t)|\le(1-\kappa)|t|^3.                  \tag{13}
\]

Indeed, factor out `|t|^3` from the terms of degree at least three,
and sum their coefficients at one.  The relevant arguments are in
this interval by (1) and the following elementary bounds:

\[
 \mu_*\ge ab\ge n-1,\qquad
 L^2\le2\mu_*,\qquad \alpha_*^2\le\mu_*/2.
\]

Using `||A tensor I + I tensor A|| <= 2L` gives

\[
 \begin{split}
 \|E_\rho\|_{\rm op}
 &\le(1-\kappa)\rho^3
       {L^2+2|\alpha_*|^3L\over\mu_*^3}\\
 &\le(1-\kappa)\rho^3
       \left({2\over\mu_*^2}+{1\over\mu_*}\right)\\
 &\le(1-\kappa)
       \left({2\over(n-1)^2}+{1\over n-1}\right).
                                                               \tag{14}
 \end{split}
\]

Here `2|alpha_*|^3 L <= mu_*^2` proves the second term.
This is an absolute `O(1/n)` estimate for every complete source.

Consequently, for any fixed finite prior on vectors
`v in [-1,1]^(n^2)`, and

\[
 F(z)=\log\sum_v\pi(v)\exp(\gamma\langle v,z\rangle),
\]

if `Y~N(0,C)` and `Z~N(0,kappa Sigma_rho+(1-kappa)I)`, Gaussian
covariance interpolation yields

\[
 |\mathbb EF(Y)-\mathbb EF(Z)|
 \le{\gamma^2n^2\over2}\|E_\rho\|_{\rm op}.        \tag{15}
\]

For completeness, `Hess(F)=gamma^2 Cov_current(v)` is positive
semidefinite with trace at most `gamma^2 n^2`; interpolate between
the two positive semidefinite endpoint covariances.  No sign of
`E_rho` is required.  The prior may depend on `A` and `gamma`, but
is fixed independently of the Gaussian disorder being compared.

By (8), `||Sigma_rho|| <= 2`.  The separately proved correlated-sign
free-energy comparison therefore applies with its absolute covariance
constant `K_G=2`, without any regularization of the source `A`.
Together with (14)-(15), this removes the source operator-norm
hypothesis from the freely chosen midpoint Gaussian reduction.
It applies, in particular, to every exact original-norm minimizer.

## 4. Actual-pressure variational condition for free alpha

This last identity is independent of midpoint optimality.  Fix a
finite deterministic prior as above.  For `alpha in [-b,a]`, let
`Z_alpha` have covariance

\[
 C_\alpha=(1-\kappa)I+\kappa\Sigma_\alpha,
 \qquad f(\alpha)=\mathbb E F(Z_\alpha).
\]

The independent component makes these covariances positive definite.
Differentiation using (4) and (9) gives the exact current-posterior
identity

\[
 \boxed{\displaystyle
 f'(\alpha)={\kappa\gamma^2(a+b)\over2\mu_\alpha^2}
 \left[b^2\mathbb E\operatorname{tr}
       ((Q\otimes Q)\operatorname{Cov}_{\rm current}(v))
       -a^2\mathbb E\operatorname{tr}
       ((P\otimes P)\operatorname{Cov}_{\rm current}(v))\right].}
                                                               \tag{16}
\]

Both weighted traces are nonnegative, but their difference has not
been proved to have either sign.  An interior pressure minimizer
balances them; endpoint minimizers obey the corresponding one-sided
conditions.  Formula (16) does not make `alpha_*` a pressure minimizer.

For `v=x tensor y` with `x,y in {-1,1}^n`, write
`e=Q_A(x)=x^T A x/2` and `f=Q_A(y)`.  The cross-field variance is

\[
 \operatorname{Var}(\langle v,Z_\alpha\rangle)
 =n^2+{\kappa\over\mu_\alpha}
           \{4ef-2\alpha n(e+f)\}.                \tag{17}
\]

Optimizing this single variance is not the same as optimizing the
expected maximum over all pairs or the actual-posterior expression
(16).  In particular, the Gaussian expected paired-norm upper
comparison needed for order transport remains open.

## 5. Uniform expected ORIGINAL norm corollary, including exact minimizers

Keep any complete source `A`, the freely chosen midpoint covariance,
and any fixed `rho in [0,1]` from Section 3.  Let `W` be an independent
standard Gaussian vector in dimension `n^2`, and set

\[
 B=\operatorname{mat}(\operatorname{sign}G),\qquad
 Z=\operatorname{mat}(\sqrt\kappa\,G+\sqrt{1-\kappa}\,W).
                                                               \tag{18}
\]

Thus the real cross matrix `Z` has precisely the linearized covariance
in (11).  For ANY deterministic real internal energy `I_n(x,y)` and
ANY deterministic real `theta_n` with `|theta_n| <= 1`, define

\[
 \mathcal M_I(D)=\max_{x,y\in\{-1,1\}^n}
                  |I_n(x,y)+\theta_n x^TDy|.
\]

The internal energy and `theta_n` may depend on `A,n`, but are fixed
before the cross disorder is drawn.  There is an absolute constant
`D_0` such that, simultaneously for all these choices,

\[
 \boxed{\displaystyle
 |\mathbb E\mathcal M_I(B)-\mathbb E\mathcal M_I(Z)|
       \le D_0 n^{16/11}.}                          \tag{19}
\]

In particular, no hypothesis on `||A||/sqrt(n)` or on the internal
energy's size is present in this statement.

Here is the quantitative deduction from the separately proved
correlated-sign quenched theorem.  For `c>0`, put `beta=c/sqrt(n)`
and introduce the auxiliary pressure

\[
 P_D=\log\mathbb E_{x,y}\cosh\!\left(
          \beta[I_n(x,y)+\theta_n x^TDy]\right).
\]

Augment by a sign `sigma`, using the fixed prior proportional to
`exp(beta sigma I_n(x,y))`.  Its cross observables are
`sigma(y tensor x)` and there are at most `2^(2n+1)` underlying
states.  The prior may depend on `c` but not on the cross draw.
The effective cross coefficient obeys
`|gamma|=|beta theta_n| <= c/sqrt(n)`.

Let `Y` be Gaussian with the exact covariance of `sign(G)`.  For
`0<epsilon<=1`, the correlated-sign theorem with `K_G=2` gives an
absolute-constant estimate

\[
 |\mathbb E P_B-\mathbb E P_{\operatorname{mat}Y}|
 \le D_1\left[
 cn\sqrt\epsilon+c\sqrt n
   +c^3\sqrt n\,\epsilon^{-4}+c^2n\epsilon\right].   \tag{20}
\]

The universal bound (14), followed by (15), gives

\[
 |\mathbb E P_{\operatorname{mat}Y}-\mathbb E P_Z|
       \le D_2 c^2                                 \tag{21}
\]

with another absolute constant.  No Gaussian-sign replacement is
inferred from covariance matching alone: (20) is precisely the
separate quenched theorem that justifies that step.

For every real matrix `D`, the exact maximum-term bounds are

\[
 \beta\mathcal M_I(D)-(2n+1)\log2
          \le P_D\le\beta\mathcal M_I(D).           \tag{22}
\]

The two pressure-to-maximum defects lie in the same interval of
length `(2n+1)log 2`.  Taking expectations, using (20)-(22), and
dividing by `beta n^(3/2)=cn` proves

\[
 { |\mathbb E\mathcal M_I(B)-\mathbb E\mathcal M_I(Z)|\over n^{3/2}}
 \le D_3\left[
  \sqrt\epsilon+n^{-1/2}
  +c^2n^{-1/2}\epsilon^{-4}+c\epsilon+c/n+1/c\right].\tag{23}
\]

Take `c=n^(1/22)` and `epsilon=n^(-1/11)`.  The four principal
terms in brackets are all `n^(-1/22)`; the remaining terms are
`n^(-1/2)` and `n^(-21/22)`.  Since `3/2-1/22=16/11`, this proves
(19).  All constants are absolute because both input covariance
estimates (8) and (14) are universal.

For the ORIGINAL norm, use `I_n(x,y)=Q_A(x)-Q_A(y)` and `theta_n=1`.
For any real cross block `D`, set

\[
 \mathcal A_D=\begin{pmatrix}A&D\\D^T&-A\end{pmatrix},
 \qquad \Phi(\mathcal A_D)
       =\max_{z\in\{-1,1\}^{2n}}|z^T\mathcal A_Dz/2|.
\]

Then `mathcal M_I(D)=Phi(mathcal A_D)` exactly, so

\[
 \boxed{\displaystyle
 |\mathbb E\Phi(\mathcal A_B)-\mathbb E\Phi(\mathcal A_Z)|
    \le D_0 n^{16/11},\qquad
 { |\mathbb E\Phi(\mathcal A_B)-\mathbb E\Phi(\mathcal A_Z)|
       \over(2n)^{3/2}}=O(n^{-1/22}).}               \tag{24}
\]

Define `m_N=min_C Phi(C)` over complete symmetric zero-diagonal
signings of order `N`.  For every outcome, `mathcal A_B` is such a
complete signing of order `2n`; hence

\[
 \boxed{\displaystyle
 m_{2n}\le\mathbb E\Phi(\mathcal A_Z)+D_0n^{16/11}
 \quad\hbox{for every complete source }A.}          \tag{25}
\]

This holds, in particular, when `A=A_n` is ANY exact minimizer of
the original order-`n` norm.  No regularization and no replacement
of that exact minimizer by a near-minimizer is necessary.

Inequality (25) is a proved all-orders, one-sided Gaussian reduction.
The dyadic upper-comparison target

\[
 \mathbb E\Phi\begin{pmatrix}A_n&Z\\Z^T&-A_n\end{pmatrix}
       \le 2\sqrt2\,m_n+o(n^{3/2})                 \tag{26}
\]

remains unproved.  A dyadic comparison with an unspecified little-o
error would itself need further error control or another argument
to settle convergence across all orders.  Equations (19) and (24)
compare expectations, not
pointwise coupled outcomes or high-probability selected signings.
The freely chosen midpoint law is not asserted to equal the
Gibbs-generated law, nor does it identify a pressure minimum with
an original-norm minimum.  The original MO limit remains OPEN.
