# Random tournaments: the sharp universal partition threshold

**Status:** proved all-orders first-moment bound and sharp asymptotic audit of
what degree two, `||Q_A||_2^2=binom(n,2)`, and `||Q_A||_infty=M` can supply.
The bound succeeds only above an explicit constant `2.00996...`, far outside
the optimal-signing regime `M/n^(3/2)<=1/2+o(1)`.  No finite signing census is
used.

Put

\[
 N={n\choose2},\qquad Q=Q_A,\qquad M=\|Q\|_\infty,
 \qquad c_0=\sqrt2-1.
\]

Choose a tournament `S` by orienting every edge independently and uniformly,
and put `R=A circ S`.  By the directed half-cut identity, the desired
zero-error doubling diamond is equivalent to

\[
 \Phi(A^{E_S^+(U)})\le\sqrt2 M
 \quad\hbox{for every }U\subseteq[n].               \tag{1}
\]

## 1. Exact annealed reduction

Fix `U` and a Boolean state `y`.  Let `h=|U|(n-|U|)` and let `I=I(U,y)` be
the contribution of edges not crossing `U`.  The random orientations split
the crossing edges independently, so

\[
 Q_{A^{E_S^+(U)}}(y)\ \buildrel d\over=\ I+S_h,     \tag{2}
\]

where `S_h` is a sum of `h` independent Rademacher signs.  If
`x_i=-y_i` on `U` and `x_i=y_i` off `U`, then

\[
 I={Q(x)+Q(y)\over2}.                                \tag{3}
\]

Let

\[
 L_A(t)=\mathbb E_y e^{tQ(y)}.
\]

Using the free eight-element symmetry of nontrivial `(U,y)` events, the two
one-sided Chernoff bounds, `h<=n^2/4`, and (3), the expected number of bad
event orbits is at most

\[
 {4^n\over8}e^{-\sqrt2\lambda M}
 (\cosh\lambda)^{n^2/4}
 \left(L_A(\lambda/2)^2+L_A(-\lambda/2)^2\right).   \tag{4}
\]

Thus any upper bound making (4) less than one produces a tournament
satisfying (1).

## 2. Best exponential-scale partition information from the stated data

The random variable `Q(y)` has mean zero, second moment `N`, and absolute
value at most `M`.  For every real `t`, the power series and
`|Q|^r<=M^(r-2)Q^2` for `r>=2` give

\[
 \boxed{
 L_A(t)\le
 1+{N\over M^2}\left(e^{|t|M}-1-|t|M\right).
 }                                                     \tag{5}
\]

Suppose `M=alpha n^(3/2)` with fixed `alpha>0`.  Equation (5) implies, for
fixed `s`,

\[
 \limsup_{n\to\infty}{1\over n}
 \log L_A(s/\sqrt n)\le |s|\alpha.                  \tag{6}
\]

Degree-two hypercontractivity does not improve (6) at speed `n`.  Indeed,

\[
 \Pr(|Q|\ge uM)
 \le\inf_{p\ge2}
 \left({(p-1)\sqrt N\over uM}\right)^p
 =\exp(-\Theta(\sqrt n))                             \tag{7}
\]

for fixed `u>0`; this is only a subexponential energy penalty.

More strongly, the right side of (6) is sharp among homogeneous degree-two
polynomials with exactly the two stated norms.  Let

\[
 k=\lfloor2\alpha\sqrt n\rfloor,qquad
 P_n(x)=a_n\sum_{1\le i<j\le k}x_ix_j,qquad
 a_n=\sqrt{{N\over {k\choose2}}}.
\]

Then `||P_n||_2^2=N` and

\[
 \|P_n\|_\infty=\sqrt{N{k\choose2}}
   =(\alpha+o(1))n^{3/2}.                            \tag{8}
\]

The maximum is attained whenever the first `k` signs agree, an event of
probability `2^(1-k)=exp(-O(sqrt(n)))`.  Hence, for every fixed `s>0`,

\[
 \lim_{n\to\infty}{1\over n}\log
 \mathbb E e^{sP_n/\sqrt n}=s\alpha.                \tag{9}
\]

Thus no universal partition theorem using only degree two, `L2`, and
`L-infinity` can replace `|s|alpha` in (6) by a smaller speed-`n` rate.
The example is deliberately not equimodular; an improvement for the MO
problem would have to use the complete unit coefficients, not merely the
three inputs audited here.

## 3. Optimized entropy-energy constant

Set `lambda=c/sqrt(n)` in (4).  Equations (5)--(6) and
`log cosh(lambda)=lambda^2/2+O(lambda^4)` give the optimized certificate
rate

\[
 \Psi_\alpha(c)
 =\log4-c(\sqrt2-1)\alpha+{c^2\over8}.              \tag{10}
\]

No other Chernoff scale improves this rate.  If `lambda sqrt(n)->0`, its
negative energy term is `o(n)` and cannot pay the `n log(4)` entropy.  If
`lambda sqrt(n)->infinity` while `lambda->0`, the positive
`lambda^2 n^2/8` term dominates the linear `lambda M` gain; if `lambda`
does not tend to zero, `log(cosh(lambda)) n^2/4` already has quadratic
order.  Hence `lambda=Theta(n^(-1/2))` is the only competitive scale.

Its minimizer and minimum are

\[
 c_*=4(\sqrt2-1)\alpha,
 \qquad
 \min_c\Psi_\alpha(c)
 =\log4-2(\sqrt2-1)^2\alpha^2.                     \tag{11}
\]

Consequently the optimized universal partition bound has negative
exponential rate precisely in the strict constant regime

\[
 \boxed{
 \alpha>{\sqrt{\log2}\over\sqrt2-1}
 =2.009964633673172\ldots .
 }                                                     \tag{12}
\]

The endpoint is also sufficient if `alpha_n=M/n^(3/2)` is eventually at
least the constant in (12).  To see the polynomial margin, use (5) directly
at `lambda=4c_0 alpha_n/sqrt(n)`.  Since
`rho=N/M^2=Theta(1/n)`, the two partition functions contribute the factor
`rho^2`; when `c_0^2 alpha_n^2=log2`, the exponential terms in (4) cancel
and the remaining bound is `O(n^(-2))`.  More generally the same choice
works whenever

\[
 n\bigl(\log2-c_0^2\alpha_n^2\bigr)-\log n
 \longrightarrow-\infty.                            \tag{13}
\]

This success range concerns arbitrary signings whose cube norm has that
large fixed `n^(3/2)` constant.  It is vacuous for the original minimizers,
for which `alpha_n<=1/2+o(1)`.

There are two different thresholds worth keeping separate.  Even under the
fictitious ideal assignment `I=0` for every central event, the exact
Rademacher moderate-deviation rate is

\[
 \log4-4\alpha^2,                                    \tag{14}
\]

so the independent first moment cannot possibly work below

\[
 \alpha_{\rm ideal}=\sqrt{{\log2\over2}}
 =0.5887050112577373\ldots .                         \tag{15}
\]

The much larger threshold (12) is what follows uniformly from the stated
energy information: those facts allow a subexponential population of states
with `|I|` near `M`, whose random allowance is only `(sqrt(2)-1)M`.

In large-deviation language, if event orbits with cut density `theta` and
`|I|/M` near `s` have energy cost `J_theta(s)` at speed `n`, their exact
first-moment exponent is

\[
 \boxed{
 \Gamma_\alpha=
 \sup_{0<\theta<1,\ 0\le s\le1}
 \left[
  \log2+H(\theta)-J_\theta(s)
  -{\alpha^2(\sqrt2-s)^2\over2\theta(1-\theta)}
 \right].
 }                                                     \tag{16}
\]

Here `H` is binary entropy with natural logarithms.  The criterion closes
when `Gamma_alpha<0` (the zero-rate boundary needs its polynomial
prefactor).  The three universal inputs above give no positive speed-`n`
cost for fixed `s`, so `J_theta(s)=0` is the only universally valid rate.
The worst point in (16) is then `s=1, theta=1/2`, yielding (11)--(12).
If all energy mass were instead at `s=0`, the same formula yields
(14)--(15).  At either threshold the central cut is uniquely extremal:
`H(theta)<=log2` and `theta(1-theta)<=1/4`.

## 4. Verdict

Using the full universal degree-two energy information does not rescue the
random-tournament route.  It exposes exactly why: at the exponential scale
needed to defeat `4^n` events, degree two plus `L2` and `L-infinity` gives
zero near-extremal entropy cost.  The route would need a new equimodular,
`A`-dependent theorem giving a positive `J_theta(s)` in (16).  Even such a
theorem cannot overcome the intrinsic ideal barrier (15), already above the
known optimal upper constant `1/2`.
