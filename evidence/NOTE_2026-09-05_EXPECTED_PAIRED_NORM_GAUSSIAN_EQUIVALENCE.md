# Expected paired Boolean norm: full-strength signs and a linearized Gaussian

2026-09-05. **Quantitative expected-norm comparison; original limit OPEN.**
This is a consequence of the independently reviewed correlated-sign
quenched theorem and canonical covariance linearization. It concerns
expectations of maxima, not pointwise agreement of the two random matrices.
No computation is used.

## 1. Statement for the actual canonical law

Fix K<infinity. For each n>=2 let A be a complete symmetric zero-diagonal
signing with `||A||op<=K sqrt(n)`. Form the canonical covariance from
this fixed host and its ACTUAL opposite-phase Gibbs covariances at any
chosen opposite temperatures:

\[
 \alpha=\frac{\operatorname{tr}(AU)+\operatorname{tr}(AV)}{2n},
 \quad H=A\otimes A-\alpha(A\otimes I+I\otimes A),\quad
 \mu=-\lambda_{\min}(H),\quad\Sigma=I+H/\mu.
\]

Draw G with centered Gaussian covariance Sigma, and let
`B=mat(sign G)` with column vectorization. This is the actual
full-strength canonical law, including its possibly singular Gaussian
covariance. Let W be an independent standard Gaussian vector and put

\[
 Y^{\rm lin}=\sqrt{2/\pi}\,G+\sqrt{1-2/\pi}\,W,
 \qquad Z=\operatorname{mat}(Y^{\rm lin}).
                                                               \tag{1}
\]

Its covariance is exactly
`C_lin=(2/pi)Sigma+(1-2/pi)I`.

Let I_n(x,y) be ANY deterministic real internal energy, and let
`|theta_n|<=1`. These may depend on A and n, but are fixed before the
cross disorder is drawn. Define

\[
 \mathcal M_I(D)=\max_{x,y\in\{-1,1\}^n}
               |I_n(x,y)+\theta_n x^TDy|.
\]

There is a constant D_K, independent of the internal energy, source,
covariance-generating temperatures, and theta_n, such that

\[
 \boxed{\quad
 \left|\mathbb E\mathcal M_I(B)-\mathbb E\mathcal M_I(Z)\right|
 \le D_K n^{16/11}.
 \quad}                                                       \tag{2}
\]

One may take `D_K<=D(1+K^4)` for an absolute constant D. In particular,
the displayed estimate also applies when the source cap K varies with n.

In particular, for `I_n(x,y)=Q_A(x)-Q_A(y)` and theta_n=1, write

\[
 \mathcal A_D=\begin{pmatrix}A&D\\D^T&-A\end{pmatrix},
 \qquad
 \Phi(\mathcal A_D)=\max_{z\in\{-1,1\}^{2n}}
                         |\tfrac12 z^T\mathcal A_Dz|.
\]

Then `mathcal M_I(D)=Phi(mathcal A_D)` exactly, so (2) is the
original maximum-absolute Boolean quadratic norm comparison

\[
 \boxed{\quad
 \left|\mathbb E\frac{\Phi(\mathcal A_B)}{(2n)^{3/2}}-
       \mathbb E\frac{\Phi(\mathcal A_Z)}{(2n)^{3/2}}\right|
 =O_K(n^{-1/22}).
 \quad}                                                       \tag{3}
\]

The first block matrix is a complete signing. The second has real
Gaussian cross coefficients; Phi has the same Boolean-max definition
on both. No half-product pressure minimum is identified with a norm
minimum, and no exact optimizer assumption is needed for this result.

## 2. Uniform covariance hypotheses

The canonical identities give

\[
 |\alpha|\le\|A\|_{\rm op}/2,\qquad
 \mu\ge(n-1)/2,\qquad
 \|H\|_{\rm op}\le2\|A\|_{\rm op}^2.
\]

Consequently

\[
 0\preceq\Sigma\preceq K_G I,
 \qquad K_G=1+8K^2,\qquad\operatorname{diag}\Sigma=1.
                                                               \tag{4}
\]

This bound is uniform in the temperatures used to generate alpha.
Those temperatures are held fixed during the following auxiliary
zero-temperature comparison; no equality with its auxiliary temperature
is required.

Let C be the exact covariance of sign G. The reviewed tensor-support
linearization proves

\[
 C=C_{\rm lin}+E_n,\qquad \|E_n\|_{\rm op}=O_K(n^{-1}).    \tag{5}
\]

Its stronger n^(-1) rate uses the bounded operator norm of A and the
two exact tensor supports. A Boolean norm cap alone only gives the
weaker rate stated in the linearization note.

## 3. Pressure comparison with free smoothing and temperature parameters

For arbitrary c>0 set beta=c/sqrt(n), and use the auxiliary pressure

\[
 P_D=\log\mathbb E_{x,y}\cosh\left(
               \beta[I_n(x,y)+\theta_n x^TDy]\right).
\]

After augmenting by sigma and subtracting the fixed internal partition
constant, this is a log partition with coefficient observables
`sigma(y tensor x)` and a deterministic prior proportional to
`exp(beta sigma I_n(x,y))`. There are at most `2^(2n+1)` underlying
states and `|gamma|=|beta theta_n|<=c/sqrt(n)`. The prior can depend
on c, but remains independent of the cross draw. The correlated-sign
theorem is uniform over all such priors.

Let Y be Gaussian with the exact matched covariance C. The full
epsilon-dependent estimate in that theorem gives, for 0<epsilon<=1,

\[
 |\mathbb E P_B-\mathbb E P_{\operatorname{mat}Y}|
 \le D_K\left[
  cn\sqrt\epsilon+c\sqrt n+
  c^3\sqrt n\,\epsilon^{-4}+c^2n\epsilon\right].           \tag{6}
\]

Gaussian covariance interpolation, retaining the actual posterior
Hessian, and (5) then give

\[
 |\mathbb E P_{\operatorname{mat}Y}-\mathbb E P_Z|
 \le\tfrac12\gamma^2 n^2\|E_n\|_{\rm op}
 \le D_K c^2.                                             \tag{7}
\]

In particular the normalized contribution from (7), after division
by c n, is O_K(c/n). No Gaussian-sign distribution is replaced merely
because its covariance matches: that replacement is supplied by the
separate quenched theorem in (6).

## 4. Passage to the expected maximum and optimized exponents

The finite maximum-term bounds hold for every real cross matrix D:

\[
 \beta\mathcal M_I(D)-(2n+1)\log2
       \le P_D\le\beta\mathcal M_I(D).                    \tag{8}
\]

The two pressure-to-maximum defects therefore lie in the same interval
of length `(2n+1)log2`. Taking expectations in (8), using (6)--(7),
and dividing by beta n^(3/2)=c n proves

\[
 \begin{aligned}
 \frac{|\mathbb E\mathcal M_I(B)-\mathbb E\mathcal M_I(Z)|}
      {n^{3/2}}
 \le D_K\bigl[
  &\sqrt\epsilon+n^{-1/2}
    +c^2n^{-1/2}\epsilon^{-4}\\
  &+c\epsilon+c/n+1/c\bigr].
 \end{aligned}                                               \tag{9}
\]

Now choose

\[
 c=n^{1/22},\qquad \epsilon=n^{-1/11}=c^{-2}.
\]

The four principal terms in (9) are exactly

\[
 \sqrt\epsilon=c\epsilon=1/c=n^{-1/22},\qquad
 c^2n^{-1/2}\epsilon^{-4}
     =n^{2/22-11/22+8/22}=n^{-1/22}.
\]

The remaining terms are n^(-1/2) and c/n=n^(-21/22), hence smaller.
This proves (2), since `3/2-1/22=16/11`, and then (3).
The pressure errors at this growing c may be O(n); their vanishing
effect on the normalized maximum follows from the explicit division
by c in (9), not a claim of uniform fixed-c pressure asymptotics.

## 5. A genuine original-norm near-minimizer family

The proved same-order original-norm regularization applies to exact
minimizers of `m_n=min_A Phi(A)`, not merely to a pressure surrogate.
With `Gamma_0=4pi/log(1+sqrt(2))`, for every r_n>0 it supplies actual
complete signings A'_n satisfying

\[
 0\le\frac{\Phi(A'_n)-m_n}{n^{3/2}}
       \le2\sqrt{2\Gamma_0/r_n},\qquad
 \|A'_n\|_{\rm op}\le(r_n+8)\sqrt n.                       \tag{10}
\]

For clarity, the uniform polynomial dependence in (2) follows directly
from the preceding proof: `K_G=1+8K^2`, so the sign-smoothing term has
at most order 1+K dependence, its entropy term and covariance smoothing
have at most order 1+K^2 dependence, and the smooth OU term has at most
order 1+K^4 dependence. The exact tensor linearization gives

\[
 \|E_n\|_{\rm op}
 \le64(1-2/\pi)K^2/n^2+16(1-2/\pi)K^4/n,
\]

so its pressure contribution has the same polynomial bound. Absorbing
all smaller powers proves `D_K<=D(1+K^4)`.

Choose any r_n tending to infinity with
`(r_n+8)^4 n^(-1/22) -> 0`, for example a log-log sequence. Form the
actual canonical source covariance from A'_n and draw its full-strength
sign cross block B_n and linearized Gaussian cross block Z_n as above.
Then (2), (10), and the paired norm identity prove simultaneously

\[
 \Phi(A'_n)=m_n+o(n^{3/2}),
\]
\[
 \mathbb E\Phi\begin{pmatrix}A'_n&B_n\\B_n^T&-A'_n\end{pmatrix}
 =\mathbb E\Phi\begin{pmatrix}A'_n&Z_n\\Z_n^T&-A'_n\end{pmatrix}
                                                    +o(n^{3/2}).
                                                               \tag{11}
\]

As one explicit quantitative choice, take `r_n=n^(1/99)`. Both the
normalized objective loss in (10) and the normalized expected norm
discrepancy in (11) are then `O(n^(-1/198))`, since

\[
 r_n^{-1/2}=n^{-1/198},\qquad
 (r_n+8)^4n^{-1/22}=O(n^{4/99-1/22})=O(n^{-1/198}).
\]

The signings A'_n are leading ORIGINAL norm near-minimizers constructed
before drawing cross disorder. No bounded-operator assertion about
every exact optimizer is required or concluded.

## Scope

The conclusion compares EXPECTATIONS of maximum absolute paired Boolean
energies. It does not assert pointwise closeness under a coupling, a
high-probability bound for individual rounded outcomes, or optimality of
a selected cross signing. Evaluating the linearized Gaussian model and
comparing it with the optimized smaller-order endpoint remain separate
tasks. The original MO limit is still OPEN.
