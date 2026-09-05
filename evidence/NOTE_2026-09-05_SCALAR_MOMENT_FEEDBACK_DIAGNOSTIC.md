# Evaluated scalar-diagonal diagnostic: strong and weak feedback

2026-09-05. This note evaluates two limiting cases of the scalar-diagonal
pure-cross two-trace bound. It is deliberately a RESTRICTED diagnostic.
Scalar optimal diagonals are not assumed for actual conditional
optimizers in general. The weak-feedback example is a counterexample
to the sufficiency of the retained MOMENT RELAXATION, not an actual
signing, shell, covariance, or Gaussian-width counterexample.

Contributions: the exact worker supplied the strong-feedback evaluation;
the docs worker supplied the weak-feedback evaluation; the root
supplied an independent coarse rational check
and helped delimit the moment-only scope. All derivations are analytic.
No tool-based mathematical computation is used.

## 1. Precisely restricted parameters and functional

Write
\[
 \kappa=2/\pi,\qquad c_0=\operatorname{asinh}1,
 \qquad u=\kappa c_0=1/K_G,\qquad f=\sqrt2,
 \qquad q=f/u,\qquad m=q^{-2}=u^2/2.                 \tag{1}
\]
The letter q here is an operator-scale parameter, NOT the quantity
`n-1` used in the original-phase spectral note.

For the scalar-optimal-diagonal diagnostic, write
\[
 D=dI_{2n},\quad d=q\sqrt n,\quad \tau=nd,
 \quad c=\beta=f n^{3/2},\quad H_B=\begin{pmatrix}0&B\\B^T&0\end{pmatrix}.
\]
An actual scalar optimum must have `d=||B||_op`. Its optimal Gram
must additionally be supported on the top eigenspace of `H_B`.
These are real attainability conditions, not consequences of the
empirical moment constraints below.

Let `y=\sigma^2/d^2` denote the normalized squared singular values and
let E denote their empirical average or a limiting probability law.
The retained Frobenius and canonical-SDP moments are
\[
 0\le y\le1,\qquad Ey=m,\qquad
 R:=S/\tau=q^2E y^{3/2},\qquad S=\operatorname{tr}|B|^3/n. \tag{2}
\]
In the zero-threshold, zero-internal-energy joint-shell diagnostic,
the uncorrected covariance has the form
\[
 M=a_0(I-\chi Q),\qquad Q=H_B/d,\qquad
 \chi=\rho_{\rm src}u q^2,
 \qquad \rho_{\rm src}={\kappa n\over\mu+\kappa}.       \tag{3}
\]
Here `a_0=n(1+kappa/mu)` and the actual intrinsic normalization is
\[
 H_{\rm src}=A\otimes A-\mathcal S_B+I,\qquad
 \mathcal S_B(X)=BX^TB,\qquad
 \mu=\max(2,\|H_{\rm src}\|_{\rm op})\ge n-1.
\]
The actual covariance is a PSD repair `M+P`, with `P>=0` and
`rank(P)<=4`; M itself must not be assumed PSD at finite n.
The available normalization bound permits the FORMAL limiting range
`0<=rho_src<=kappa`, without proving each point attainable.
The weak endpoint means a limit, not a zero finite-n covariance
coefficient.

For `0<=eta<1`, the normalized two-trace/Jensen functional is
\[
 \begin{split}
 U_\eta(\chi,\nu)
 &=\eta\sqrt{(1-u)E A_\eta(y)}
       +(1-\eta)\sqrt\kappa\sqrt{E B_\eta(y)},\\
 A_\eta(y)
 &={1+(\eta^2-2\chi\eta-2\eta+\chi)y+\chi\eta^2y^2
       \over(1-\eta^2y)^2},\\
 B_\eta(y)
 &={1+(\eta^2-2\chi\eta)y\over(1-\eta^2y)^2}.
 \end{split}                                             \tag{4}
\]
The width normalization is `2n sqrt(a_0)`; the leading target at (1)
is `f/2=1/sqrt(2)`.

For clarity, (4) comes from the two spectral traces of
`P_eta=D-eta H_B` and `E_eta=(1-eta)D`. Their normalized resolvents are
\[
 T_\eta=E{1-\chi\eta y\over1-\eta^2y},\qquad
 R_\eta=E B_\eta(y),\qquad
 T_\eta-(1-\eta)R_\eta=\eta E A_\eta(y).
\]
These identities explain all factors in (4). They do not assert that
the uncorrected M is a covariance when it is indefinite.

Here is an explicit repair-trace bound, not an appeal to an unspecified
previous comparison. The canonical source is
`evidence/NOTE_2026-09-05_INTRINSIC_CROSS_JOINT_SHELL_REPAIR.md`,
equations (4) and (8), which give `P>=0` and
\[
 {\operatorname{tr}P\over a_0}
       \le2\left(1+{L_A^2\over\mu+1}\right).
\]
Since `||S_B||_op=L_B^2=d^2=q^2n`, the definition of intrinsic mu
and the triangle inequality imply
\[
 L_A^2\le\mu+q^2n+1,\qquad
 {\operatorname{tr}P\over a_0}\le4+2q^2,                \tag{4a}
\]
using `mu+1>=n`. No original-norm cap or conference-scale bound on A
is needed. The exact repaired traces, normalized by `2na_0/d`, obey
\[
 \begin{split}
 0\le\widetilde T_\eta-T_\eta
   &\le {2+q^2\over n(1-\eta)},\\
 0\le\widetilde R_\eta-R_\eta
   &\le {2+q^2\over n(1-\eta)^2}.                       \tag{4b}
 \end{split}
\]
Indeed `||P_eta^{-1}||<=1/[d(1-eta)]`, so the respective trace
increments are bounded by this inverse norm times `tr P`, and by
`d ||P_eta^{-1}||^2 tr P`, before normalization. Thus both increments
vanish for fixed eta and bounded q. The two-trace square-root expression
is applied to the actual PSD repaired field first. Passing its traces
to their limits is then justified by (4b) and continuity, even if M
was indefinite at finite n.

## 2. The complete retained mixed-rounding constraint

The scalar-diagonal tensor-mixture theorem gives, at leading order,
for every fixed `0<t<1`, with `a(t)=asinh(1-t)`,
\[
 {\beta\over\kappa\tau}
 \ge a(t)+tR
    -{t(\sec a(t)-1)\over4}(1-R).                        \tag{5}
\]
Equivalently it bounds `1-R` below by the scalar mixed envelope.
At (1), the left side is `c_0`.

One may retain the stronger fourth-moment residual bound instead:
\[
 {\beta\over\kappa\tau}
 \ge a(t)+tR
    -{t(\sec a(t)-1)\over4}
                   [1-2R+q^2Ey^2].                     \tag{6}
\]
The bracket is nonnegative since it equals
`q^2 E[y(1-sqrt(y))^2]`. It follows by retaining
`N_-/tau <= ||dB-B|B|||_F^2/(4nd tau)` in the scalar sign-defect
proof. Indeed the latter ratio is exactly
`[1-2R+q^2Ey^2]/4`.

Equations (5)--(6) are necessary constraints of the existing mixed
rounding argument, not a sufficient characterization of attainable
sign matrices or optimal Gram matrices.

## 3. Strong feedback repairs the old extremal diagnostic

Take `rho_src=kappa` in the limit. Equations (1)--(3) give
\[
                    \chi=2/c_0>2.                       \tag{7}
\]
Actual repaired positivity supplies an essential additional condition.
The number of negative eigenvalues of M cannot exceed four, because
adding a rank-at-most-four matrix makes it PSD. For positive chi,
those negative eigenvalues correspond precisely to singular values
with `chi sqrt(y)>1`. Hence every limiting empirical law obeys
\[
                 \operatorname{supp}\nu
                          \subset[0,\chi^{-2}]
                          =[0,c_0^2/4]\subset[0,1/4).     \tag{8}
\]

The limiting law is separated from one. Thus in (4), as `eta` increases
to one, the second term tends to zero and
\[
 U_1^2=(1-u)E{1-\chi y\over1-y}.
\]
The integrand is concave because its second derivative is
`2(1-chi)/(1-y)^3<0`. Jensen and (2) therefore give
\[
 \inf_{0\le\eta<1}U_\eta^2
 \le { (1-u)(1-\chi m)\over1-m}
 = { (1-u)(1-\kappa u)\over1-u^2/2}
                                      <{9\over20}<{1\over2}. \tag{9}
\]
Here is an entirely rational proof of the strict estimate. One has
`3/5<kappa<2/3`, and
\[
 c_0=\int_0^1{dx\over\sqrt{1+x^2}}
       >\int_0^1(1-x^2/2)\,dx=5/6,\qquad c_0<1.
\]
Thus `1/2<u<2/3` and `kappa u>3/10`. The numerator in (9) is below
`7/20` and its denominator is above `7/9`, proving the printed bound.

The order of limits is important. At finite n, evaluate the PSD repaired
trace functional at each FIXED eta first and use the explicit
rank-four trace bound (4b). Then pass to a limiting empirical law. Only
then take eta to one, using (8). The inequality
`limsup_n inf_eta F_n(eta) <= limsup_n F_n(eta)` for fixed eta
justifies the resulting upper bound. No uniform finite-n estimate
as eta tends to one is asserted.

Thus strong feedback and repaired positivity already remove the old
extremal scalar point from this diagnostic obstruction, even without
using the extra mixed moment restriction.

## 4. Weak feedback leaves an obstruction to the moment relaxation

Now take the weak-feedback limit `rho_src=0`, hence `chi=0`, and the
formal bulk law
\[
                              \nu=\delta_m.              \tag{10}
\]
It has the correct mean and
\[
                   R=q^2m^{3/2}=\sqrt m=u/\sqrt2.
\]
It satisfies the ENTIRE family (5), not merely its small-t endpoint.
Indeed `a(t)` is concave and `a'(0)=-1/sqrt(2)`, so
\[
 a(t)\le c_0-t/\sqrt2,\qquad
 a(t)+tR\le c_0-{t(1-u)\over\sqrt2}<c_0.                \tag{11}
\]
The terms subtracted in (5) and (6) are nonnegative. Thus (11) also
proves feasibility for the stronger fourth-moment version (6), for
every t simultaneously.

For this law the complete two-trace functional is
\[
 U_\eta={\eta\sqrt{(1-u)(1-2m\eta+m\eta^2)}
       +(1-\eta)\sqrt\kappa\sqrt{1+m\eta^2}
                       \over1-m\eta^2}.                 \tag{12}
\]
Put
\[
 b=1/\sqrt2,\qquad S_0=\sqrt\kappa,\qquad
 A_0=\sqrt{(1-u)(1-m)}.
\]
Since `1-2m eta+m eta^2 >= 1-m` and `sqrt(1+m eta^2)>=1`,
\[
 U_\eta\ge {\eta A_0+(1-\eta)S_0\over1-m\eta^2}.
\]
Subtracting b from this lower bound leaves the numerator
\[
 (S_0-b)-(S_0-A_0)\eta+bm\eta^2
 \ge\Delta:=(S_0-b)-{(S_0-A_0)^2\over4bm}.              \tag{13}
\]
The constant Delta is strictly positive, as verified below. Since
`0<1-m eta^2<=1`, (13) proves the uniform strict obstruction
\[
             \boxed{\quad U_\eta\ge b+\Delta>b
                  \quad\hbox{for all }0\le\eta\le1.\quad} \tag{14}
\]
At chi zero the negative-eta resolvent traces are the same as the
positive-eta traces, while the first shell factor replaces `1-u` by
`1+u`. Thus negative eta cannot improve (14).

### Exact coarse rational certificate for Delta

The classical bounds `3.14<pi<22/7` imply
`7/11<kappa<100/157`, hence
\[
 .7977<S_0<.7981,\qquad .7071<b<.7072.
\]
Also `.88<c_0<.89`. For example the cubic term of `sinh(.88)` is
below `.114`, and all subsequent term ratios are below `.04`, so
`sinh(.88)<.88+.114/.96=.99875<1`. Meanwhile the cubic term of
`sinh(.89)` exceeds `.11`, so `sinh(.89)>.89+.11=1`.
Consequently
\[
 .56<u<.567,\qquad .1568<m<.1608.
\]
All decimals here denote exact terminating rationals. They give
\[
 A_0^2>(.433)(.8392)=.3633736>(.6028)^2=.36336784,
\]
so `0<S_0-A_0<.1953`; positivity also follows from
`A_0^2<1-u<.44<kappa`. Finally
\[
 4bm(S_0-b)
 >4(.156)(.707)(.09)=.03970512
 >(.1953)^2=.03814209>(S_0-A_0)^2.                        \tag{15}
\]
This proves Delta positive by exact rational inequalities. No numerical
scalar calculation was run to establish it.

## 5. Positive weak feedback and the missing attainability condition

The pure law (10), with `q>1`, is NOT itself the full spectrum of an
actual scalar-optimal-D candidate: it has no singular value d. This
must not be hidden by calling it an actual optimizer counterexample.

The present note does not construct a top-singular-value repair or an
actual optimal Gram. Possible vanishing top outliers are not ruled out
by the bulk law, but no such realization is asserted here.

For the formal law (10), the obstruction persists for sufficiently
small positive chi.
For `0<=chi<1` and `0<=y<=1`, the covariance spectral factors
`1+-chi sqrt(y)` are bounded below by `1-chi`. Both positive trace
terms defining (4) are linear in these covariance factors before
taking square roots. Consequently
\[
 U_\eta(\chi,\nu)\ge\sqrt{1-\chi}\,U_\eta(0,\nu).        \tag{16}
\]
This holds on both eta branches. Equations (14) and (16) give a uniform
positive obstruction for the formal bulk law when source feedback is
sufficiently weak. Formal PSD support is satisfied in this case as well.

## 6. What the evaluation does and does not settle

The extremal strong-feedback case passes the target, by (9), for all
permitted limiting laws with the repaired-positivity restriction.
The weak-feedback case does not: the retained moment relaxation allows
laws for which even the eta-optimized full two-trace/Jensen functional
is strictly above the target, despite all currently retained mixed
rounding constraints.

No complete sign matrix realizing this formal diagnostic has been
constructed. In particular the actual optimal-Gram support condition, coordinate
structure, source/cross compatibility, and exact endpoint attainability
are not supplied by this moment calculation. The result is not a lower
bound on an actual Gaussian shell width and is not a counterexample to
the desired original-order upper comparison. It identifies a precise
limitation of this restricted moment-only proof route when feedback
is allowed to vanish.
