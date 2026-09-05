# Canonical sign covariance and matched-Gaussian linearization

2026-09-05. **Exact covariance identity and a uniform Gaussian-to-Gaussian
pressure comparison.** The resulting error is `O(sqrt(n))` under a norm
cap, and `O(1)` under a fixed bounded-operator cap. This note does not
replace the Gaussian signs themselves by Gaussian coefficients.

## 1. Canonical law and its covariance

Let `n>=2`, let `A` be a symmetric zero-diagonal complete sign matrix,
and write its extreme eigenvalues as `a>0` and `-b<0`. Set
`L=max(a,b)=||A||op`. Let `U,V` be the actual zero-field Gibbs covariance
matrices of this fixed host at opposite inverse temperatures. Write

\[
 p=\operatorname{tr}(AU),\qquad q=\operatorname{tr}(AV),\qquad
 \alpha=\frac{p+q}{2n},\qquad
 H=A\otimes A-\alpha(A\otimes I+I\otimes A).
\]

Convexity of the two log partition functions and the spectral bounds give

\[
 0\le p\le na,\quad -nb\le q\le0,\quad
 -b/2\le\alpha\le a/2,\quad |\alpha|\le L/2.
                                                               \tag{1}
\]

Define

\[
 \mu=ab+\alpha(a-b),\qquad
 \Sigma_\rho=I+\rho H/\mu,\qquad 0\le\rho\le1.
\]

For completeness, the eigenvalues of `H` are
`xy-alpha(x+y)`, where `x,y` are eigenvalues of `A`. This bilinear
function is minimized on `[-b,a]^2` at a corner. Its same-sign corner
values are nonnegative by (1), and its two mixed values are `-mu`.
Thus `lambda_min(H)=-mu`.

Also `(x+b)(a-x)>=0` on the same interval, so
`x^2<=(a-b)x+ab`. Summing over the spectrum and using
`tr A=0`, `tr A^2=n(n-1)`, proves `ab>=n-1`. If `a>=b`, use
`alpha>=-b/2` in the definition of `mu`; if `a<=b`, use
`alpha<=a/2`. In either case,

\[
 \mu\ge ab/2\ge(n-1)/2>0.                                  \tag{2}
\]

Consequently `Sigma_rho` is positive semidefinite, with diagonal one,
including at full strength `rho=1`.

Let `Z~N(0,Sigma_rho)`, and put `B=sign(Z)` using column vectorization
of the `n by n` cross block. Each coordinate has variance one, so its
sign is defined almost surely even if the joint Gaussian is singular.
The two-coordinate Gaussian sign identity gives

\[
 C_\rho:=\mathbb E[BB^T]
 =I+\frac2\pi\left[
       \arcsin(\tau)A\otimes A
       -\arcsin(\tau\alpha)(A\otimes I+I\otimes A)
                     \right],\qquad \tau=\rho/\mu.          \tag{3}
\]

Here `B` inside `BB^T` denotes the vectorized block. All arcsines in
(3) are scalar. Indeed, for distinct tensor coordinates, either both
indices differ, giving covariance `tau` times a product of two signs,
or exactly one differs, giving `-tau alpha` times one sign. These
supports are disjoint, and oddness of arcsine gives (3). The diagonal
is one. This also proves that all displayed scalar arguments lie in
`[-1,1]`.

Set `kappa=2/pi` and `b_0=1-2/pi`. Equation (3) yields the exact
linearization

\[
 \boxed{\quad C_\rho=\kappa\Sigma_\rho+b_0 I+E_\rho,\quad}
                                                               \tag{4}
\]
\[
 E_\rho=\kappa\left[
 (\arcsin\tau-\tau)A\otimes A
 -(\arcsin(\tau\alpha)-\tau\alpha)
                            (A\otimes I+I\otimes A)\right].
                                                               \tag{5}
\]

## 2. Operator-norm error, including full strength

For `|u|<=1`, the positive-coefficient power series of arcsine gives

\[
 |\arcsin u-u|\le(\pi/2-1)|u|^3.                            \tag{6}
\]

Explicitly, `arcsin u=u+sum_{k>=1} a_k u^(2k+1)` with `a_k>0`;
each absolute power in the remainder is at most `|u|^3`, and the sum
of its coefficients is `pi/2-1`. Endpoint values follow by monotone
convergence. Thus (5), tensor-product operator norms, and the triangle
inequality prove

\[
 \boxed{\quad
 \|E_\rho\|_{\rm op}
 \le b_0\left(\frac\rho\mu\right)^3
                     (L^2+2|\alpha|^3L).
 \quad}                                                     \tag{7}
\]

There is no entry-count factor `n^2` in (7); the two exact tensor
supports must be retained.

If only `alpha^2<=Cn`, the trivial `L<=n-1` and (2) imply

\[
 \|E_\rho\|_{\rm op}
 \le8b_0\rho^3\left[
        \frac1{n-1}
       +\frac{2C^{3/2}n^{3/2}}{(n-1)^2}\right]
 =O_C(\rho^3 n^{-1/2}).                                    \tag{8}
\]

An original Boolean norm cap `Phi(A)<=K_0 n^(3/2)` supplies this
hypothesis with `C=K_0^2`: since `p>=0>=q` and both absolute energies
are at most `2Phi(A)`, one has `|alpha|<=Phi(A)/n`.

Under the stronger fixed cap `L<=K sqrt(n)`, (1)--(2) give instead

\[
 \boxed{\quad
 \|E_\rho\|_{\rm op}
 \le\frac{8b_0\rho^3}{(n-1)^3}
                  \left(K^2n+\frac{K^4n^2}{4}\right)
 =O_K(\rho^3 n^{-1}).
 \quad}                                                     \tag{9}
\]

All bounds are uniform over `rho in [0,1]`. No small-strength limit
or nonsingularity of the original canonical covariance is assumed.

## 3. Replacing the covariance of a matched Gaussian

Let `I(x,y)` be any fixed real internal interaction, with
`x,y in {-1,1}^n`. For a centered Gaussian cross vector `G_R` of
covariance `R`, define its quenched symmetric pressure

\[
 \mathcal F_I(R)=\mathbb E_{G_R}\log\mathbb E_{x,y}
          \cosh\left(I(x,y)+\gamma\langle x\otimes y,G_R\rangle\right).
                                                               \tag{10}
\]

The internal interaction and its coefficients may be arbitrary, but
are fixed independently of the Gaussian draw. Put

\[
 C_\rho^{\rm lin}=\kappa\Sigma_\rho+b_0 I.
\]

This covariance is positive definite, with smallest eigenvalue at
least `b_0>0`. Interpolate by
`R_s=C_rho^lin+s E_rho`, `0<=s<=1`. Both endpoints are covariance
matrices, so every `R_s` is positive semidefinite. For `s<1` it is
positive definite.

Introduce the actual augmented-spin observable
`T=sigma(x tensor y)`, where `sigma` is also a uniform sign. Its
squared Euclidean norm is exactly `n^2`. Gaussian differentiation of
the FULL log partition function gives

\[
 \frac{d}{ds}\mathcal F_I(R_s)
 =\frac{\gamma^2}{2}\mathbb E
           \operatorname{tr}\left(E_\rho\operatorname{Cov}_s(T)\right).
                                                               \tag{11}
\]

The covariance in (11) belongs to the current Gaussian-dependent
posterior, not the original internal Gibbs law. It is positive
semidefinite and has trace
`n^2-||<T>_s||_2^2<=n^2`. Hence

\[
 \boxed{\quad
 |\mathcal F_I(C_\rho)-\mathcal F_I(C_\rho^{\rm lin})|
 \le\frac{\gamma^2n^2}{2}\|E_\rho\|_{\rm op}.
 \quad}                                                     \tag{12}
\]

Integrate (11) on `s<1` and pass to the endpoint by continuity to
cover a possibly singular `C_rho`. Finite-spin log partitions are
Lipschitz in their Gaussian vector, so this passage is valid.

In the paired scaling `gamma=c sqrt(t/(2n))`, (12) becomes

\[
 |\mathcal F_I(C_\rho)-\mathcal F_I(C_\rho^{\rm lin})|
 \le\frac{c^2tn}{4}\|E_\rho\|_{\rm op}
 =\begin{cases}
 O_{c,C}(\sqrt n),&\alpha^2\le Cn,\\
 O_{c,K}(1),&\|A\|_{\rm op}\le K\sqrt n.
 \end{cases}                                                \tag{13}
\]

These estimates are uniform over all fixed internal interactions,
`t in [0,1]`, and `rho in [0,1]`.

## Scope

The vector of canonical Gaussian signs has covariance `C_rho`.
Equations (12)--(13) compare a Gaussian vector with THAT covariance
against a Gaussian vector with the linearized covariance. They do
not compare the signs to either Gaussian vector. Covariance agreement
alone supplies no such pressure equivalence. That separate genuine
sign-to-Gaussian comparison, and any subsequent order transport,
remain unproved by this note.
