# One-phase Gaussian reduction and an actual-posterior variance upper bound

2026-09-05. **Quenched reductions and finite-step differential inequality.**
The derivative retains the current coupled posterior throughout. Its
integrated comparison with the smaller-order optimum is not proved here.

## 1. Removing the augmented phase at subextensive cost

Let Z be a centered Gaussian n by n matrix. Its vectorized covariance C
may be singular, but is assumed invariant under matrix transposition.
Let I(x,y) be any fixed deterministic internal interaction satisfying
`I(y,x)=-I(x,y)`, and let gamma be real. Define

\[
 f(Z)=\log\mathbb E_{x,y}
               e^{I(x,y)+\gamma x^TZy},
 \qquad
 P(Z)=\log\mathbb E_{x,y}
               \cosh(I(x,y)+\gamma x^TZy).
\]

Here x,y are independent uniform Boolean vectors. Interchanging them
in the negative phase gives exactly

\[
 P(Z)=\log\frac{e^{f(Z)}+e^{f(-Z^T)}}2.                    \tag{1}
\]

The two random variables on the right have the same law, because the
Gaussian covariance is transpose-invariant and the Gaussian is centered.
The arithmetic-geometric mean inequality and the maximum bound give

\[
 0\le\mathbb E P-\mathbb E f
 \le\frac12\mathbb E|f(Z)-f(-Z^T)|
 \le\sqrt{\operatorname{Var}(f(Z))}.                     \tag{2}
\]

The gradient of f with respect to the cross vector is
`gamma <y tensor x>`, where the expectation is under the actual
Gaussian-dependent posterior. Thus its squared Euclidean norm is at
most `gamma^2 n^2`. Gaussian Poincare, including singular covariance,
therefore yields

\[
 \boxed{\quad
 0\le\mathbb E P-\mathbb E f
 \le|\gamma|n\sqrt{\|C\|_{\rm op}}.
 \quad}                                                     \tag{3}
\]

For completeness, Poincare here follows directly from the Gaussian
covariance identity

\[
 \operatorname{Var}f(Z)=\int_0^\infty e^{-u}
    \mathbb E[\nabla f(Z)^TC\nabla f(Z_u)]\,du,
\]

where `Z_u=e^(-u)Z+sqrt(1-e^(-2u))Z'` and Z' is independent.
The gradient bound proves (3). A Gaussian factor representation makes
the same identity valid for singular C.

In particular, `gamma=O(n^(-1/2))` and `||C||op=O(1)` give an
`O(sqrt(n))` pressure error, uniformly in the internal interaction.

There is also a direct maximum version. For any fixed antisymmetric
internal energy J and real theta, put

\[
 h(Z)=\max_{x,y}[J(x,y)+\theta x^TZy],\qquad
 M(Z)=\max_{x,y}|J(x,y)+\theta x^TZy|.
\]

Then `M=max(h(Z),h(-Z^T))`. The maximum is Lipschitz with constant
`|theta|n` in the cross vector. Applying (2) and Poincare, or a smooth
maximum approximation, proves

\[
 \boxed{\quad
 0\le\mathbb E M-\mathbb E h
 \le|\theta|n\sqrt{\|C\|_{\rm op}}.
 \quad}                                                     \tag{4}
\]

This is an expectation bound, not pointwise agreement of the phases.

## 2. A PSD split with the full current posterior

Let Sigma be any correlation matrix on the n^2 cross coordinates, and
write

\[
 C=(1-s)I+s\Sigma,\qquad 0\le s\le1.
\]

At any fixed Gaussian disorder value, use the one-phase posterior on
(x,y), and set

\[
 T=y\otimes x,\quad w=\langle T\rangle,
 \quad W=\langle xy^T\rangle,
 \quad V=\operatorname{Cov}(T).
\]

One has `||T||_2^2=n^2` and `||w||_2^2=||W||_F^2`. Therefore

\[
 \begin{aligned}
 \operatorname{tr}(CV)
 &=(1-s)(n^2-\|W\|_F^2)
       +s\bigl[\langle T^T\Sigma T\rangle-w^T\Sigma w\bigr]\\
 &\le(1-s)(n^2-\|W\|_F^2)
                         +s\langle T^T\Sigma T\rangle.
 \end{aligned}                                                \tag{5}
\]

Only the nonnegative term `s w^T Sigma w` has been discarded. In
particular, the mean of the rank-one observables is not set to zero;
the mean matrix W itself need not have rank one.

Specialize to a fixed complete symmetric zero-diagonal source A and

\[
 H=A\otimes A-\alpha(A\otimes I+I\otimes A),\qquad
 \Sigma=I+H/\mu\succeq0,
 \qquad \mu>0.
\]

No Gibbs-energy formula for alpha is needed for the following identity.
Direct tensor multiplication in (5) gives the pointwise posterior bound

\[
 \boxed{\quad
 \begin{aligned}
 \operatorname{tr}(CV)
 \le{}&n^2+\frac{s}{\mu}
  \left\langle
   (x^TAx)(y^TAy)-\alpha n(x^TAx+y^TAy)
  \right\rangle\\
 &-(1-s)\|W\|_F^2.
 \end{aligned}
 \quad}                                                     \tag{6}
\]

This is an upper bound on the actual quenched Gaussian variance term.
It keeps the favorable opposite-energy term and an explicit posterior
overlap subtraction. Both are evaluated after the cross coupling is
present, not in the original independent opposite-phase Gibbs prior.

## 3. The actual finite-step derivative

Fix beta>0 and fix A, alpha, mu, s before the interpolation. Suppose
the preceding Sigma is positive semidefinite. For 0<=t<=1 set

\[
 \eta_t=\beta\sqrt{1-t/2},\quad
 \gamma_t=\beta\sqrt{t/2},\quad
 F(t)=\mathbb E_Z\log\mathbb E_{x,y}
       e^{\eta_t[Q_A(x)-Q_A(y)]+\gamma_t x^TZy},
 \qquad \operatorname{Cov}(\operatorname{vec}Z)=C.
\]

Let angle brackets and W_t denote its CURRENT posterior. Gaussian
covariance differentiation and `eta_t'=-beta^2/(4eta_t)` give exactly

\[
 F'(t)=\frac{\beta^2}{4}\mathbb E
 \left[\operatorname{tr}(C V_t)
       -\frac1{\eta_t}\langle Q_A(x)-Q_A(y)\rangle\right].
                                                               \tag{7}
\]

The formula holds on (0,1), extends to the right derivative at zero,
and is integrable up to both endpoints. There are finitely many spin
states, eta_t stays positive, and the Gaussian derivatives have uniform
finite bounds at each fixed n. Singular C can again be handled by a
Gaussian factor representation or an independent positive diagonal
regularization.

Combining (6) and (7) proves the actual-posterior differential upper

\[
 \boxed{\quad
 \begin{aligned}
 F'(t)\le\frac{\beta^2}{4}\mathbb E\Bigg[
 &n^2+\frac{s}{\mu}
    \left\langle4Q_A(x)Q_A(y)
            -2\alpha n[Q_A(x)+Q_A(y)]\right\rangle\\
 &-(1-s)\|W_t\|_F^2
       -\frac1{\eta_t}\langle Q_A(x)-Q_A(y)\rangle
 \Bigg].
 \end{aligned}
 \quad}                                                     \tag{8}
\]

Define `Z_A(u)=E_x exp(u Q_A(x))` and
`a_A(u)=(log Z_A(u)+log Z_A(-u))/2`.
At t=0, the one-phase partition factors and

\[
 F(0)=\log Z_A(\beta)+\log Z_A(-\beta)=2a_A(\beta).
                                                               \tag{9}
\]

Integrating (7) gives the exact change from this endpoint; integrating
(8) gives its upper bound. Neither step by itself bounds that change
by o(n). In particular, no sign has been proved for the integral on
the right of (8).

## 4. Application to the Gaussianized sign proposal

The matched linearized Gaussian proposal has `s=2/pi` and
`C=(1-2/pi)I+(2/pi)Sigma`. Its tensor form is transpose-invariant,
so all preceding statements apply.

For the spectrally centered choice
`alpha=(a-b)/2`, `mu=(a^2+b^2)/2`, where a and -b are the two extreme
eigenvalues of A, the separate midpoint covariance factorization gives
`0<=Sigma<=2I`. Thus

\[
 \|C\|_{\rm op}\le1+2/\pi.
\]

Consequently the errors in (3)--(4) have absolute constants for EVERY
complete source A; no bounded-source-operator hypothesis is needed
for this centered proposal. The true Gaussian-sign-to-Gaussian
equivalence is supplied by the separately reviewed theorem, not by
the phase reduction in this note.

The remaining upper comparison requires control of the CURRENT energy
product and overlap terms in (8), possibly with an optimized covariance
choice. The source Gibbs covariances and source optimality cannot be
substituted for these coupled-posterior quantities. No original
symmetric minimum is identified with a half-product minimum here.
