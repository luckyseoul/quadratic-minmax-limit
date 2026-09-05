# Shifted thresholds: the even-Hermite correction and a simple mean Gaussian

2026-09-05.  Status: proved Gaussian covariance-to-norm reduction and
shifted-sign expected-norm corollary.  The threshold-sign replacement
uses the separately proved and independently complete-read
mean-preserving shifted-sign universality theorem, cited in Section 5.

## 1. Setup and the unconditional Gaussian statement

Let `A` be a complete symmetric zero-diagonal signing of order `n>=2`,
and write `a=lambda_max(A)>0`, `-b=lambda_min(A)<0`.  Use the freely
chosen spectral midpoint

\[
 \alpha=(a-b)/2,\quad \mu=(a^2+b^2)/2,\quad
 \Sigma=I+\{A\otimes A-\alpha(A\otimes I+I\otimes A)\}/\mu.
                                                               \tag{1}
\]

The independently proved midpoint identities give

\[
 \operatorname{diag}\Sigma=1,\quad0\preceq\Sigma\preceq2I,
 \quad\mu\ge n-1,\quad L^2\le2\mu,\quad\alpha^2\le\mu/2,
 \qquad L=\|A\|_{\rm op}.                                  \tag{2}
\]

Let `G~N(0,Sigma)`, use column vectorization, and for ANY real `h` set

\[
 (B_h)_{ij}=\operatorname{sign}(G_{ij}+hA_{ij}),\quad
 s=2\Phi_{\rm Gauss}(h)-1,\quad
 k=4\phi(h)^2,\quad v=1-s^2-k.                             \tag{3}
\]

Here `phi` is the standard normal density.  The diagonal thresholds
are zero, because `A_ii=0`.  Thus `E B_h=sA`, including its diagonal.
The constant `v` is nonnegative, by the first Hermite coefficient
bound proved below.

Let `Y_h` be the Gaussian matrix with mean `sA` and centered covariance
EXACTLY `Cov(vec B_h)`.  Independently draw a standard Gaussian matrix
`W`, and define the simpler model

\[
 Z_h=sA+\sqrt{k}\,G+\sqrt v\,W.                            \tag{4}
\]

For arbitrary deterministic real `I_n(x,y)` and `|theta|<=1`, put

\[
 \mathcal M_I(D)=\max_{x,y\in\{-1,1\}^n}
               |I_n(x,y)+\theta x^TDy|.
\]

There is an absolute constant `C` such that, uniformly in all these
choices, including `A,h,I_n,theta`,

\[
 \boxed{\displaystyle
 |\mathbb E\mathcal M_I(Y_h)-\mathbb E\mathcal M_I(Z_h)|
 \le C\left[n+\Phi(A)\sqrt{\log(2n)\over n}\right].}       \tag{5}
\]

Here `Phi(A)=max_x |x^TAx/2|` is the ORIGINAL source norm.  In
particular, `Phi(A)<=C_0 n^(3/2)` makes the right side
`O_{C_0}(n sqrt(log n))=o(n^(3/2))`.  No source operator cap is needed.
Equation (5) itself concerns two Gaussian models.  Section 5 supplies
the threshold-sign replacement using the separate universality theorem.

## 2. Centered threshold Hermite coefficients and parity

For a scalar standard Gaussian `g`, define

\[
 f_h(g)=\operatorname{sign}(g+h)-s.
\]

Expand in the normalized probabilists' Hermite basis:

\[
 f_h(g)=\sum_{j\ge1}c_j(h){\operatorname{He}_j(g)\over\sqrt{j!}}.
\]

The mean is zero and Parseval gives

\[
 \sum_{j\ge1}c_j(h)^2=1-s^2\le1.
\]

Direct one-dimensional integration gives `c_1(h)=2 phi(h)`, hence
`k=c_1(h)^2` and `v=sum_{j>=2} c_j(h)^2>=0`.  Also

\[
 f_{-h}(g)=-f_h(-g),\qquad
 c_j(-h)=(-1)^{j+1}c_j(h).                                 \tag{6}
\]

For two standard Gaussian variables with correlation `t`, the Hermite
identity is `E[He_j(g)He_l(g')]=1_{j=l} j! t^j`.  It follows, for
example, by comparing coefficients in the joint Gaussian exponential
generating function.  The resulting covariance series converges
absolutely also at `t=+/-1`, by Cauchy-Schwarz and Parseval.

Define

\[
 o(t)=\sum_{\substack{j\ge3\\j\ \mathrm{odd}}}c_j(h)^2t^j,
 \qquad e(t)=\sum_{\substack{j\ge2\\j\ \mathrm{even}}}c_j(h)^2t^j.
\]

For `|t|<=1`,

\[
 |o(t)|\le|t|^3,\qquad0\le e(t)\le t^2.                   \tag{7}
\]

On the off-diagonal cross-coordinate set
`O={(i,j):i!=j}`, the threshold is either `+h` or `-h`.
By (6), its odd Hermite coefficients are unchanged and its even
coefficients acquire the sign `A_ij`.  Thus the covariance of two
distinct coordinates `z,z' in O`, whose Gaussian correlation is `t`,
is exactly

\[
 kt+o(t)+A_z A_{z'}e(t).                                  \tag{8}
\]

The variance at each coordinate in `O` is `1-s^2`.

## 3. Exact tensor supports and the low-rank even term

Let `Pi` project onto `O`, let `D_A=diag(vec A)`, and let `J` be the
`n x n` all-one matrix.  Thus `D_A^2=Pi` and `||D_A||<=1`.
All the covariance statements in this section are on `O`; equivalently
they may be embedded as matrices zero on its complementary coordinates.
Put

\[
 q=1/\mu,\quad r=\alpha/\mu,\quad u=e(q),\quad w=e(r).
\]

The linear-plus-independent base covariance on `O` is the restriction
of `C_0=k Sigma+vI`.  Distinct cross coordinates with both indices
different have correlation `q A_ik A_jl`; those sharing one index
have correlation `-r A_ik` or its counterpart.  Hence (8) gives

\[
 C_h|_O=C_0|_O+R_{\rm odd}+R_{\rm even},                   \tag{9}
\]

where

\[
 R_{\rm odd}=\Pi\left[o(q)A\otimes A
             -o(r)(A\otimes I+I\otimes A)\right]\Pi,       \tag{10}
\]
\[
 R_{\rm even}=D_A\left[u(J-I)\otimes(J-I)
              +w\{(J-I)\otimes I+I\otimes(J-I)\}\right]D_A.
                                                               \tag{11}
\]

Both remainders have zero diagonal; the residual variance `vI`
already accounts for all their diagonal Hermite mass.

Define the positive semidefinite matrix

\[
 \mathcal L=uD_A(J\otimes J)D_A
       +(w-u)_+D_A(J\otimes I+I\otimes J)D_A.              \tag{12}
\]

It has rank at most `2n+1`.  Expanding (11) shows

\[
 R_{\rm even}-\mathcal L
  =(u-2w)\Pi-(u-w)_+D_A(J\otimes I+I\otimes J)D_A.        \tag{13}
\]

By (2) and (7),

\[
 \|R_{\rm odd}\|\le{L^2+2|\alpha|^3L\over\mu^3}
                      \le2/\mu^2+1/\mu,
\]

while `u<=mu^(-2)`, `w<=alpha^2/mu^2`, and
`||J tensor I+I tensor J||=2n` give

\[
 \boxed{\displaystyle
 C_h|_O=(C_0|_O+\mathcal L)+R,\qquad
 \|R\|\le{2n+3\over\mu^2}+{2\over\mu}\le{18\over n}.} \tag{14}
\]

Every norm in this section is the operator norm.  The two covariance
matrices compared in (14) are positive semidefinite, even though
`R` need not have either sign.  The even correction is not discarded
in operator norm: its generally nonvanishing low-rank part is retained.

## 4. Its Gaussian norm cost is subleading on original norm-capped sources

The centered Gaussian matrix with covariance (12) can be realized as

\[
 U=\sqrt u\,\xi A+\sqrt{(w-u)_+}
                \{\operatorname{diag}(g)A+A\operatorname{diag}(g')\},
                                                               \tag{15}
\]

where `xi`, the `n` coordinates of `g`, and those of `g'` are
independent standard Gaussians.  It is zero on the diagonal.

Let `beta(A)=max_{x,y in {-1,1}^n}|x^TAy|`.  The original norm gives

\[
 \beta(A)\le4\Phi(A).                                    \tag{16}
\]

Indeed, a zero-diagonal quadratic form has absolute value at most
`Phi(A)` throughout `[-1,1]^n`, by independent sign rounding with
prescribed coordinate means.  Apply polarization to `(x+y)/2`
and `(x-y)/2`.  The bilinear maximum also bounds every pair of
vectors in that cube, by multilinearity.

It follows from (15), the standard Gaussian finite-maximum bound,
and `sqrt(u)<=1/mu`, `sqrt((w-u)_+)<=|alpha|/mu`, that

\[
 \begin{split}
 \mathbb E\beta(U)
 &\le4\Phi(A)\left[
 {\sqrt{2/\pi}\over\mu}
       +{2|\alpha|\over\mu}\sqrt{2\log(2n)}\right]\\
 &\le C\Phi(A)\sqrt{\log(2n)\over n}.                    \tag{17}
 \end{split}
\]

This is an actual Boolean-norm estimate, not a rank-only argument.

For clarity, the remaining operator-norm error in (14) costs only
`O(n)` in the expected maximum.  In general, if two centered Gaussian
covariances on `m` coordinates differ by at most `delta` in operator
norm, then for any fixed finite family of coefficient vectors of
norm at most `sqrt(m)` and arbitrary fixed additive offsets,

\[
 |\mathbb E\max_\ell(a_\ell+v_\ell\cdot X_1)
       -\mathbb E\max_\ell(a_\ell+v_\ell\cdot X_2)|
       \le\sqrt{2\delta m\log N},                         \tag{18}
\]

where `N` is the number of states.  To prove (18), use
`C_1 <= C_2+delta I` and Gaussian convex order, realized by adding
independent Gaussian noise with the PSD covariance difference.
Then couple `C_2+delta I` as `X_2+sqrt(delta) W` and bound the
increase of the maximum by `sqrt(delta) max_l v_l dot W`.
The exponential Gaussian maximum bound proves one direction;
interchange the two covariances for the other.  Singular covariances
cause no exception.  Gaussian means can be absorbed into the offsets.

For `mathcal M_I`, augment by the sign of the absolute value, so
`N<=2^(2n+1)` and `m<=n^2`.  With `delta=18/n`, (18) is at most
`sqrt(90 log 2) n`.  Adding the independent Gaussian (15) changes
the expected maximum by a nonnegative quantity at most `E beta(U)`:
the lower bound is conditional Jensen, and the upper bound is the
pointwise bilinear-norm Lipschitz estimate.  This proves (5) on `O`.

Finally restore the `n` cross-diagonal entries.  Their means are zero.
The matched Gaussian has variance one there and the simple Gaussian
has variance `1-s^2<=1`; their expected sums of absolute diagonal
entries are each at most `n`.  Replacing these entries by zero thus
costs at most `2n` altogether, irrespective of their correlations
with the other entries.  This proves the full statement (5).

## 5. Combined shifted-sign expected ORIGINAL norm reduction

The separately proved theorem
`NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md`
(final source SHA-256
`a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`)
supplies the absolute bound

\[
 |\mathbb E\mathcal M_I(B_h)-\mathbb E\mathcal M_I(Y_h)|
                      \le D n^{16/11},                    \tag{19}
\]

uniformly for `Sigma<=2I`.  Combining it with (5) gives

\[
 |\mathbb E\mathcal M_I(B_h)-\mathbb E\mathcal M_I(Z_h)|
 \le D n^{16/11}
       +C\left[n+\Phi(A)\sqrt{\log(2n)\over n}\right].    \tag{20}
\]

Equation (19) uses the separate mean-preserving quenched proof, not
an inference from matching first and second moments.  Both (19) and
the combined estimate (20) are unconditional proved comparisons.

For `I_n(x,y)=Q_A(x)-Q_A(y)` and `theta=1`, these are precisely the
expected ORIGINAL norms of the paired blocks `[A,D;D^T,-A]`.
Consequently, on an original norm-capped source, the even-Hermite
correction causes only a subleading `O(n sqrt(log n))` cost.
This includes every exact original minimizer family.  More explicitly,
let `m_N=min_C Phi(C)` over complete signings of order `N`, and choose
ANY exact minimizer `A_n` at order `n`.  The independent-Rademacher
exponential-maximum bound gives the
elementary cap `Phi(A_n)<=sqrt(log 2) n^(3/2)`.  Indeed, there are
`2^(n-1)` distinct Boolean energies up to global spin reversal,
each a sum of `n(n-1)/2` independent signs; the exponential bound
for their absolute maximum is at most
`sqrt(n^2(n-1) log 2)<=sqrt(log 2) n^(3/2)`.

For every deterministic `h`, the block `[A_n,B_h;B_h^T,-A_n]` is
a complete signing of order `2n`.  Therefore (20), with the original
internal energy, proves the actual one-sided optimized reduction

\[
 \boxed{\displaystyle
 m_{2n}\le\inf_{h\in\mathbb R}\mathbb E\Phi
     \begin{pmatrix}A_n&Z_h\\Z_h^T&-A_n\end{pmatrix}
       +D' n^{16/11}.}                                    \tag{21}
\]

The constant `D'` is absolute: the original norm cap makes the
additional error `O(n sqrt(log n))`, which is absorbed by
`n^(16/11)`.  The bound is uniform in h, so taking the infimum
needs no minimizing threshold to exist.  The threshold is chosen
deterministically before the disorder is drawn, not adaptively
from a realized Gaussian matrix.

This optimized Gaussian reduction does not evaluate its infimum or
provide the still-missing
upper comparison for the simple mean Gaussian model (4), or any
sign for its actual mean/noise interpolation derivative.  The
original MO limit remains OPEN.
