# A norm cap forces an extensive response to bounded dense fields

2026-09-05. **All-orders actual-Gibbs lower bounds.** These are not
upper bounds on susceptibility, not spectral flatness, and not a
cross-order comparison. The elementary tensor-rounding argument below
proves the semidefinite estimate used in the proof; no unreviewed
Grothendieck theorem is imported. No computation was run.

## 1. Statements

Let J be a real symmetric zero-diagonal n by n matrix, and write

\[
 \Phi(J)=\max_{x\in\{-1,1\}^n}\left|\tfrac12x^TJx\right|.
\]

Suppose Phi(J)<=b n, where b is independent of n. Put

\[
 \kappa=\log(1+\sqrt2),\qquad B=1+\frac{4\pi}{\kappa}b.
                                                               \tag{1}
\]

For a sign vector v and a real h, let the actual Ising law be

\[
 \mu_{J,hv}(x)\ \propto\ \exp\left(\tfrac12x^TJx+h v^Tx\right).
\]

Then

\[
 \boxed{\quad
 \operatorname{Var}_{\mu_{J,hv}}(v^TX)
 \ge\frac n2\exp\left[-2\sqrt{4h^2+8B(1+4B)}\right].
 \quad}                                                       \tag{2}
\]

This lower bound is uniform in J, v, and n under the stated norm cap.
The variance may be much larger than n. In particular, defining

\[
 \Psi_{J,v}(h)=\log\mathbb E_{\mu_{J,0}}e^{h v^TX},
\]

global spin reversal gives Psi'(0)=0, and integration of (2) gives

\[
 \boxed{\quad
 \Psi_{J,v}(h)
 \ge\frac{n h^2}{4}
       \exp\left[-2\sqrt{4h^2+8B(1+4B)}\right].
 \quad}                                                       \tag{3}
\]

There is also a nonuniform-field version. Let H>0, 0<s<=H^2, and let
w be any real vector satisfying

\[
 \|w\|_\infty\le H,\qquad \|w\|_2^2\ge s n.
\]

For every 0<=t<=1,

\[
 \boxed{\quad
 \operatorname{Var}_{\mu_{J,tw}}(w^TX)
 \ge\frac{s n}{2}
   \exp\left[-2\sqrt{2H^2+
                     \frac{8BH^4(1+4B)}{s^2}}\right].
 \quad}                                                       \tag{4}
\]

Consequently

\[
 \boxed{\quad
 \log\mathbb E_{\mu_{J,0}}e^{w^TX}
 \ge\frac{s n}{4}
   \exp\left[-2\sqrt{2H^2+
                     \frac{8BH^4(1+4B)}{s^2}}\right].
 \quad}                                                       \tag{5}
\]

All constants are explicit and positive at fixed b,H,s. The statements
apply equally to J and -J, with the same b. Their useful setting here
is J=beta A_T, an ACTUAL complement Hamiltonian, not a replacement by
the uniform spin law.

More generally, no global bound on the field is needed if it has a
positive-density set of moderate nonzero coordinates. Suppose a>0,
H>=a, 0<delta<=1, and

\[
 \#\{i:a\le|w_i|\le H\}\ge\delta n.
\]

There is no restriction on the other coordinates of w. Then, for
every 0<=t<=1,

\[
 \boxed{\quad
 \operatorname{Var}_{\mu_{J,tw}}(w^TX)
 \ge\frac{a^2\delta n}{2}
 \exp\left[-2\sqrt{2H^2+\frac{8B(1+4B)}{\delta^2}}\right],
 \quad}                                                       \tag{13}
\]

and

\[
 \boxed{\quad
 \log\mathbb E_{\mu_{J,0}}e^{w^TX}
 \ge\frac{a^2\delta n}{4}
 \exp\left[-2\sqrt{2H^2+\frac{8B(1+4B)}{\delta^2}}\right].
 \quad}                                                       \tag{14}
\]

The field is not clipped in these inequalities. Its arbitrarily large
coordinates remain in the actual Gibbs law throughout the proof.

## 2. Elementary tensor rounding and a diagonal majorizer

We first prove the finite matrix estimate needed for the augmentation.
For any real array K and finitely many unit vectors u_i,v_j,

\[
 \left|\sum_{ij}K_{ij}\langle u_i,v_j\rangle\right|
 \le\frac\pi{2\kappa}
       \max_{x_i,y_j\in\{-1,1\}}\left|\sum_{ij}K_{ij}x_i y_j\right|.
                                                               \tag{6}
\]

Indeed, define the Hilbert direct-sum vectors

\[
 U_i=\bigoplus_{k\ge0}
       \sqrt{\frac{\kappa^{2k+1}}{(2k+1)!}}\,
                          u_i^{\otimes(2k+1)},
\]
\[
 V_j=\bigoplus_{k\ge0}(-1)^k
       \sqrt{\frac{\kappa^{2k+1}}{(2k+1)!}}\,
                          v_j^{\otimes(2k+1)}.
\]

Their norms are one because sinh(kappa)=1, and
`<U_i,V_j>=sin(kappa <u_i,v_j>)`. The Gram matrix of these finitely
many vectors has a finite-dimensional realization. Apply a single
standard Gaussian to that realization and take the signs of its
projections, obtaining x_i,y_j. The planar Gaussian sign identity and
`0<kappa<1<pi/2` give exactly

\[
 \mathbb E x_i y_j
 =\frac2\pi\arcsin\left(\sin(\kappa\langle u_i,v_j\rangle)\right)
 =\frac{2\kappa}\pi\langle u_i,v_j\rangle.
\]

The absolute value of the expected bilinear form is at most its maximum
over sign pairs. This proves (6); no infinite-dimensional Gaussian
random variable is required.

Next, for a symmetric zero-diagonal J,

\[
 \max_{x,y\in\{-1,1\}^n}|x^TJy|\le4\Phi(J).                \tag{7}
\]

For completeness, zero diagonal and independent spin averaging imply
`|Q_J(z)|<=Phi(J)` for every z in [-1,1]^n. Apply this to
`(x+y)/2` and `(x-y)/2` and use
`x^TJy=2[Q_J((x+y)/2)-Q_J((x-y)/2)]`.

Consider separately the two finite semidefinite programs

\[
 \max\{\operatorname{tr}(\pm J X):X\succeq0,
                                      \operatorname{diag}X=\mathbf1\}.
\]

Every feasible X is a Gram matrix of unit vectors. Equations (6)--(7)
bound each optimum by `(2 pi/kappa) Phi(J)`. Their diagonal duals are

\[
 \min\left\{\sum_i d_i^{\pm}:
        \operatorname{diag}(d^{\pm})\mp J\succeq0\right\}.
\]

Strong duality applies: X=I is strictly positive definite and feasible,
and sufficiently large scalar diagonal matrices are strictly dual
feasible. The diagonal entries d_i^+ and d_i^- are nonnegative because
J has zero diagonal. The dual minima are attained; bounded objective
sublevels have bounded nonnegative diagonal coordinates.

Set

\[
 D=I+\operatorname{diag}(d^++d^-),\qquad
 \overline d=\frac1n\operatorname{tr}D.
\]

Then D is positive diagonal, its entries are at least one, and

\[
 \boxed{\quad D-J\succeq I,\quad D+J\succeq I,
 \qquad \overline d\le1+\frac{4\pi\Phi(J)}{\kappa n}\le B.
 \quad}                                                       \tag{8}
\]

This is a trace-controlled diagonal augmentation. It does not claim
that the largest diagonal entry, or ||J||_op, is bounded.

## 3. Actual Gaussian augmentation and weighted moments

Put M=J+D. Then M is positive definite, `diag M=diag D`, and

\[
 0\preceq M\preceq2D,\qquad
 M D^{-1}M\preceq2M\preceq4D.                              \tag{9}
\]

For the second inequality, conjugate the first by D^(-1/2); the resulting
PSD matrix has eigenvalues in [0,2], so its square is at most twice
itself, and conjugate back.

The standard Gaussian augmentation for the actual law mu_(J,f) can be
written as a joint law of X and g. Start with the centered Gaussian
of covariance M, weight it together with a uniform spin X by
`exp((f+g)^T X)`, and normalize. Integrating g gives the desired spin
law because
`X^T M X=X^T J X+tr D`. Conditional on g, the spins are independent
with means `tanh(f_i+g_i)`. Conditional on X,

\[
 g\mid X\sim N(MX,M).
\]

Thus, independently of the external field f,

\[
 \begin{aligned}
 \mathbb E[g^T D^{-1}g]
 &=\operatorname{tr}(D^{-1}M)
       +\mathbb E[X^T M D^{-1}M X]\\
 &\le n+4\operatorname{tr}D
 \le n(1+4B).                                               \tag{10}
 \end{aligned}
\]

Here `tr(D^-1 M)=n` uses the zero diagonal of J, and
`X^T D X=tr D` uses that the coordinates of X are signs. The expectation
in (10) is under the full, field-dependent augmented law, not the prior
Gaussian. This distinction is essential.

For any real vector w, conditional independence and total variance give

\[
 \operatorname{Var}_{\mu_{J,f}}(w^TX)
 \ge\sum_i w_i^2\,\mathbb E\operatorname{sech}^2(f_i+g_i).
                                                               \tag{11}
\]

## 4. Signed fields

Let `I={i:d_i<=2 overline d}`. Then |I|>=n/2. For f=h v with v a
sign vector, (10) implies

\[
 \sum_{i\in I}\mathbb E(hv_i+g_i)^2
 \le 2h^2 n+4\overline d\,n(1+4\overline d).
\]

The average second moment over these coordinates is therefore at most
`4h^2+8 overline d(1+4 overline d)`. Since
`sech^2 u>=exp(-2|u|)`, Jensen's inequality followed by Cauchy--Schwarz,
over both the coordinates and their actual augmented laws, gives

\[
 \sum_{i\in I}\mathbb E\operatorname{sech}^2(hv_i+g_i)
 \ge\frac n2\exp\left[-2\sqrt{4h^2+8B(1+4B)}\right].
\]

Combine this with (11), where v_i^2=1, to prove (2). For h>=0,
`Psi(h)=int_0^h(h-t)Psi''(t)dt`; the bound in (2) decreases with |t|.
This proves (3) for h>=0 and then for all h by evenness.

## 5. Bounded dense nonuniform fields

Suppose `||w||_infinity<=H` and `||w||_2^2>=s n`. Set

\[
 K=\frac{2BH^2}{s},\qquad I=\{i:d_i\le K\},\qquad
 q=\sum_{i\in I}w_i^2.
\]

At most `(tr D)/K` coordinates are outside I. Hence

\[
 q\ge sn-H^2\frac{\operatorname{tr}D}{K}\ge\frac{sn}{2}.
\]

For f=t w, 0<=t<=1, the weighted average second moment satisfies

\[
 \begin{aligned}
 \frac1q\sum_{i\in I}w_i^2\mathbb E(tw_i+g_i)^2
 &\le2H^2+\frac{2H^2K}{q}\mathbb E[g^TD^{-1}g]\\
 &\le2H^2+\frac{8BH^4(1+4B)}{s^2}.                         \tag{12}
 \end{aligned}
\]

Apply Jensen and Cauchy--Schwarz using coordinate weights w_i^2/q,
then multiply by q and use (11). This proves (4). The zero-field law
has global spin-reversal symmetry, so its mean of w^T X is zero.
Integrating the second derivative from t=0 to t=1 proves (5).

## 6. Moderate coordinates inside an otherwise unrestricted field

Let `I_0={i:a<=|w_i|<=H}`, so `|I_0|>=delta n`, and set

\[
 K=\frac{2B}{\delta},\qquad I=I_0\cap\{i:d_i\le K\}.
\]

At most `tr D/K<=delta n/2` coordinates fail the diagonal bound,
so `|I|>=delta n/2`. Under the ACTUAL law with full field f=t w,
0<=t<=1, equation (10) is still valid independently of all coordinates
of f. On I alone the fields are bounded by H, giving

\[
 \begin{aligned}
 \frac1{|I|}\sum_{i\in I}\mathbb E(tw_i+g_i)^2
 &\le2H^2+\frac{2K}{|I|}\mathbb E[g^TD^{-1}g]\\
 &\le2H^2+\frac{8B(1+4B)}{\delta^2}.
 \end{aligned}
\]

Apply the same unweighted Jensen/Cauchy--Schwarz estimate on I and
use `w_i^2>=a^2` there in (11). This proves (13). The initial law is
zero-field and spin-reversal symmetric, so integrating its log moment
generating function's second derivative gives (14). Large fields outside
I affect the posterior law but do not invalidate (10); no monotonicity
under coordinatewise clipping has been assumed.

## Scope

The result uses a Boolean energy cap and actual Gibbs structure. It does
not replace that structure by a generic symmetric cavity measure, whose
behavior need not satisfy this conclusion. It does not bound either
covariance operator from above, and it permits macroscopic magnetization
and susceptibility. The constants are positive but deliberately not
optimized. A separate argument is needed to turn this field response
into an exclusion for actual global minimizers or into order transport.
