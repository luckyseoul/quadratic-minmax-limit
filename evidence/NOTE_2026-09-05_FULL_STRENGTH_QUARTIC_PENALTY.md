# Quartic cost and boundary covariance of full-strength canonical rounding

2026-09-05. **Actual-law estimates; the boundary posterior integral is OPEN.**
The cross signs below have the canonical correlated Gaussian-sign law,
not an independent coefficient-refill law. No computation is used.

## 1. Canonical source and sign covariance

Let A be a complete symmetric zero-diagonal signing of order n>=2 with
`Phi(A)<=C n^(3/2)`. Use its ACTUAL opposite-phase Gibbs covariances U,V
at the chosen covariance-generating temperature and define

\[
 \alpha=\frac{\operatorname{tr}(AU)+\operatorname{tr}(AV)}{2n},
 \quad H=A\otimes A-\alpha(A\otimes I+I\otimes A),
 \quad\mu=-\lambda_{\min}(H),\quad T=H/\mu.
\]

The reviewed canonical construction gives

\[
 |\alpha|\le C\sqrt n,\qquad \mu\ge(n-1)/2,\qquad
 \Sigma=I+T\succeq0,\qquad\operatorname{diag}\Sigma=1.
                                                               \tag{1}
\]

Draw G with centered Gaussian covariance Sigma and set
`vec(B)=sign(G)`, with column vectorization. The host is fixed before
this draw. Singular Sigma is allowed. Set

\[
 k_0=\frac2\pi\arcsin(1/\mu),\qquad
 k_1=\frac2\pi\arcsin(\alpha/\mu).
\]

The exact entrywise Gaussian-sign covariance is

\[
 \mathbb E[\operatorname{vec}B\operatorname{vec}B^T]
 =I+k_0 A\otimes A-k_1(A\otimes I+I\otimes A).            \tag{2}
\]

In particular, writing S_j=tr(A^j),

\[
 \begin{aligned}
 \mathbb E BB^T=\mathbb E B^TB&=nI-nk_1A,\\
 \mathbb E\operatorname{tr}(B^TA^2B)
 =\mathbb E\operatorname{tr}(B^TBA^2)&=nS_2-nk_1S_3,\\
 \mathbb E\operatorname{tr}(B^TAB A)&=k_0 S_2^2.
 \end{aligned}                                               \tag{3}
\]

No posterior replacement is involved in these distributional identities.

## 2. A four-sign rectangle bound

For four centered unit-variance Gaussian variables with correlation
matrix I+X, let f be the product of their four signs. If
`max_(i!=j)|X_ij|<=1/6`, then `||X||_op<=1/2`. Interpolate the Gaussian
density p_s with covariance I+sX, 0<=s<=1. At s=0,
`E f=0` and `d(E f)/ds=0`: in the first density derivative every term
contains at most two Gaussian coordinates, leaving two independent
mean-zero sign factors.

For completeness, with `Y_s=(I+sX)^(-1/2)X(I+sX)^(-1/2)` and a standard
Gaussian z, the density score and its derivative have the laws

\[
 \frac{p_s'}{p_s}=\tfrac12(z^TY_sz-\operatorname{tr}Y_s),
 \qquad
 (\log p_s)''=\tfrac12\operatorname{tr}Y_s^2-z^TY_s^2z.
\]

Their Gaussian second-moment and first-absolute-moment bounds imply

\[
 \int|p_s''|\le2\operatorname{tr}Y_s^2\le8\|X\|_F^2.
\]

Twice integrating the density derivative therefore gives

\[
 |\mathbb E f|\le4\|X\|_F^2
                  =8\sum_{i<j}X_{ij}^2.                 \tag{4}
\]

Apply this to the four distinct coordinates at the corners of a
rectangle in B. Four pair correlations have magnitude |alpha|/mu,
and the other two have magnitude 1/mu. Thus (4) bounds the rectangle
expectation by `32(1+alpha^2)/mu^2` whenever its small-correlation
hypothesis holds. If it fails, the trivial absolute bound one is at
most `36(1+alpha^2)/mu^2`. Consequently, without a smallness assumption,

\[
 |\mathbb E B_{ij}B_{ik}B_{lk}B_{lj}|
 \le\frac{36(1+\alpha^2)}{\mu^2}
 \quad(i\ne l,\ j\ne k).                                \tag{5}
\]

The cases with a repeated row or column contribute exactly
`2n^3-n^2` to the fourth singular moment. Hence, with
`R=tr((BB^T)^2)`,

\[
 \boxed{\quad
 \mathbb E R\le2n^3-n^2+
             \frac{36n^2(n-1)^2(1+\alpha^2)}{\mu^2}.
 \quad}                                                       \tag{6}
\]

This includes the singular full-strength Gaussian law.

## 3. Exact expected paired quartic penalty

For 0<=t<=1 put `a=sqrt(2-t)`, `b=sqrt(t)`, and

\[
 M_t=\frac1{\sqrt{2n}}
        \begin{pmatrix}aA&bB\\bB^T&-aA\end{pmatrix},
 \qquad v=\operatorname{tr}(A/\sqrt n)^4=S_4/n^2.
\]

Block multiplication gives exactly

\[
 \begin{aligned}
 4n^2\operatorname{tr}M_t^4
 ={}&2a^4S_4+4a^2b^2(X+Y-Z)+2b^4R,\\
 X={}&\operatorname{tr}(B^TA^2B),\quad
 Y=\operatorname{tr}(B^TBA^2),\quad
 Z=\operatorname{tr}(B^TAB A).
 \end{aligned}
\]

For example, the off-diagonal block of the squared unnormalized
matrix is `ab(AB-BA)`, and its squared Frobenius norm is
`a^2 b^2(X+Y-2Z)`. Using (3) proves

\[
 \boxed{\quad
 \mathbb E\operatorname{tr}M_t^4
 =\frac{a^4}{2}v+
 \frac{a^2b^2}{n^2}(2nS_2-2nk_1S_3-k_0S_2^2)
                 +\frac{b^4}{2n^2}\mathbb ER.
 \quad}                                                       \tag{7}
\]

At t=0 this is exactly 2v. To bound its possible increase, discard
the nonpositive k_0 term. Cauchy--Schwarz gives
`|S_3|<=sqrt(S_2 S_4)=n sqrt(n(n-1)) sqrt(v)`.
Set `k=|k_1|sqrt(n-1)`. The inequality
`|(2/pi)arcsin u|<=|u|` gives, by (1), `k^2<=8C^2`.
Thus the internal quartic change plus the absolute cubic term is at most

\[
 -\frac{t(4-t)}2v+2t(2-t)k\sqrt{nv}
 \le\frac{2t(2-t)^2}{4-t}k^2n
 \le2tk^2n.                                               \tag{8}
\]

The first inequality maximizes a quadratic in sqrt(v); it also holds
at t=0 by continuity. The last uses `(2-t)^2<=4-t`.
Equations (1) and (6) also give

\[
 \frac{\mathbb ER}{2n^2}
 \le n+72(1+\alpha^2)\le(37+72C^2)n.
\]

Since `2a^2b^2 S_2/n<=4tn`, equations (7)--(8) establish

\[
 \boxed{\quad
 \mathbb E\operatorname{tr}M_t^4
 \le2\operatorname{tr}(A/\sqrt n)^4+(41+88C^2)tn.
 \quad}                                                       \tag{9}
\]

Therefore this actual full-strength cross law adds at most
`(41+88C^2)lambda t n` to the paired endpoint penalty. No row-fourth
bound or exact optimality was needed beyond the source norm cap.
For penalized minimizing sources whose cap C is uniform over
0<lambda<=1, this is a uniform O(lambda n) penalty cost. It is not an
o(n) bound at fixed lambda.

## 4. The actual pair-conditioned boundary covariance

The same source also gives a precise boundary estimate for
`Sigma_rho=I+rho T`, 0<=rho<=1. Every column of T has the SAME squared norm

\[
 q=\|Te_e\|_2^2
   =\frac{(n-1)[(n-1)+2\alpha^2]}{\mu^2}
   \le4+16C^2,                                           \tag{10}
\]

by direct expansion of the diagonal of H^2. Moreover,
`max_(e!=f)|T_ef|<=max(1,|alpha|)/mu=O_C(n^(-1/2))`.
For all sufficiently large n each principal pair covariance is
invertible uniformly in rho. For P={e,f} and its complement R, the
EXACT covariance conditioned on G_e=G_f=0 is

\[
 \begin{aligned}
 \Sigma^{ef}_{\rho,R}
 &=\Sigma_{\rho,RR}-K^{ef}_\rho,\\
 K^{ef}_\rho
 &=\rho^2[T_{Re},T_{Rf}]
   \begin{pmatrix}1&\rho T_{ef}\\\rho T_{ef}&1\end{pmatrix}^{-1}
   \begin{bmatrix}T_{eR}\\T_{fR}\end{bmatrix}.
 \end{aligned}                                               \tag{11}
\]

If `max_(e!=f)|T_ef|<=1/2`, the subtracted matrix is positive
semidefinite of rank at most two, with

\[
 \operatorname{tr}K^{ef}_\rho\le4q\le16+64C^2,
 \qquad\max_j(K^{ef}_\rho)_{jj}=O_C(n^{-1}).               \tag{12}
\]

These assertions remain valid at rho=1. The conditional covariance
can then be singular; its range is the actual linear boundary support
before thresholding. It is still pair-dependent. A bounded-trace,
rank-two correction does not identify its sign law with the
unconditioned law, nor identify the subsequently cross-tilted Gibbs
posterior with the original paired prior.

## Scope

The canonical full-strength proposal has controlled quartic cost,
and (11) explicitly retains the correct boundary covariance and
support. The new actual-source force bounds are not automatically
bounds for the cross-tilted posterior: that larger paired host need
not be a penalized global minimizer. No favorable bound for the
weighted posterior-covariance integral in the full-strength likelihood
identity has been proved here. Exact quartic-minimizer success or
failure, and the original convergence problem, remain OPEN.
