# Original MO problem as a sharp equimodular L1-influence constant

**Status:** proved exact reformulation, proved sharp two-state calculus
inequality, and proved six-variable counterexample to the literal conjectural
degree bound for general bounded functions.  None of these statements proves
that the original limit exists.

Write

\[
 Q_A(x)=\sum_{i<j}A_{ij}x_ix_j,\qquad
 M(A)=\|Q_A\|_{L^\infty(\{\pm1\}^n)},\qquad
 m_n=\min_A M(A),
\]

where `A` is symmetric, has zero diagonal, and has off-diagonal entries in
`{+1,-1}`.  For a function on the cube use the discrete derivative

\[
 D_i f(x)={f(x)-f(x^{\oplus i})\over2}
\]

and total `L1` influence

\[
 \operatorname{Inf}_1(f)=\sum_i\mathbb E|D_i f|.
\]

## 1. Exact equivalence, including the leading constant

Let

\[
 \mu_k=\mathbb E|\varepsilon_1+\cdots+\varepsilon_k|
 = {k\binom{k-1}{\lfloor(k-1)/2\rfloor}\over2^{k-1}}.
\]

For `f_A=Q_A/M(A)`, direct differentiation gives

\[
 D_i f_A(x)={x_i\over M(A)}\sum_{j\ne i}A_{ij}x_j.
\]

The absolute value of the last sum has exactly the law of the absolute value
of a sum of `n-1` independent Rademacher signs, independently of `A` and of
`i`.  Therefore

\[
 \boxed{\operatorname{Inf}_1(f_A)={n\mu_{n-1}\over M(A)}.}       \tag{1}
\]

Define the sharp constant in the **real, complete, equimodular degree-two
class** by

\[
 K_n=\max_A\operatorname{Inf}_1(Q_A/M(A)).
\]

Since the numerator in (1) is independent of the signing,

\[
 \boxed{K_n={n\mu_{n-1}\over m_n},\qquad
 {m_n\over n^{3/2}}={\mu_{n-1}/\sqrt n\over K_n}.}               \tag{2}
\]

In particular, since `mu_(n-1)/sqrt(n) -> sqrt(2/pi)` and the sharp constants
are bounded away from zero and infinity, the MathOverflow limit exists if and
only if `K_n` converges.  If `K_n -> K`, its value is

\[
 \lim_n {m_n\over n^{3/2}}={\sqrt{2/\pi}\over K}.                \tag{3}
\]

This is an exact equivalence, not a lower-bound application of influence
theory.  The difficulty has moved entirely into whether the sharp constants
of the non-nested equimodular classes converge.

The known interval

\[
 {1\over\pi}\le\liminf {m_n\over n^{3/2}}
 \le\limsup {m_n\over n^{3/2}}\le {1\over2}
\]

corresponds asymptotically to the reverse interval

\[
 \sqrt{8/\pi}\ \lesssim K_n\lesssim\sqrt{2\pi}
\]

in the appropriate liminf/limsup sense.

## 2. Exact best-response interpolation

The homogeneous-polynomial argument can be made exact in degree two.  Let
`f=Q_A/M(A)`.  For a Boolean `x`, choose

\[
 T(x)_i=\operatorname{sgn}((Ax)_i),
\]

with either sign allowed when the entry vanishes.  Put

\[
 a=f(x),\qquad b=\sum_i|D_i f(x)|,\qquad c=f(T(x)).
\]

The continuous multilinear extension of `f` is bounded by one on the full
cube.  Homogeneity therefore gives, for every real `t`,

\[
 h_x(t):=f(x+tT(x))=a+bt+ct^2,
 \qquad |h_x(t)|\le(1+|t|)^2.                                  \tag{4}
\]

The middle coefficient is exactly `b`, because

\[
 {x^TAT(x)\over M(A)}={\|Ax\|_1\over M(A)}
 =\sum_i|D_i f(x)|.
\]

Apply the upper bound in (4) at `t>0` and the lower bound at `-t`.  Optimizing
the two resulting bounds over `t` gives

\[
 \boxed{
 b\le 2+2\min\left\{
 \sqrt{(1-a)(1-c)},\sqrt{(1+a)(1+c)}
 \right\}.}                                                     \tag{5}
\]

Equivalently,

\[
 b\le2+2\sqrt{1+ac-|a+c|}
 \le2+2\sqrt{1-a^2}\le4.                                      \tag{6}
\]

Thus the exact quantity missing from one-variable Harris calculus is the
joint energy of a state and its synchronous best response.  Averaging (5)
gives the rigorous `A`-dependent route

\[
 K(A)\le 2+2\,\mathbb E_x
 \sqrt{1+f(x)f(Tx)-|f(x)+f(Tx)|}.                              \tag{7}
\]

There is also a coarser mean-energy consequence.  With
`eta_A=E_x f(Tx)`, averaging (4) first and then optimizing in `t` gives

\[
 \boxed{K(A)\le2+2\sqrt{1-|\eta_A|}.}                          \tag{8}
\]

Equations (7)--(8) are true for every signing.  They do not close the MO
problem: conference-like candidates can have `eta_A=0`, and (6) is least
informative where the normalized energy `f(x)` is near zero.  In fact

\[
 \mathbb E f(x)^2={\binom n2\over M(A)^2}=\Theta(1/n)
\]

at the relevant `M(A)=Theta(n^(3/2))` scale, so uniform states concentrate
precisely in that scalar blind spot.  A successful use of (7) would need new
control of the **push-forward energy** `f(Tx)`, not another one-point moment
bound for `f(x)`.

## 3. A six-variable theorem: the conjectural `Inf <= degree` ceiling is false

Consider

\[
A=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&-1&-1&1\\
1&1&0&1&-1&-1\\
1&-1&1&0&1&-1\\
1&-1&-1&1&0&1\\
1&1&-1&-1&1&0
\end{pmatrix}.                                                  \tag{9}
\]

Direct multiplication gives `A^2=5I`, while `Q_A(1)=5`.  For every Boolean
`x`,

\[
 |Q_A(x)|={1\over2}|x^TAx|\le {1\over2}\sqrt5\,\|x\|_2^2
 =3\sqrt5<7.
\]

The value `Q_A(x)` is an odd integer because it is a sum of fifteen signs.
Consequently `|Q_A(x)|<=5`, and the value at `1` shows `M(A)=5` exactly.
For `f=Q_A/5`, equation (1) and `mu_5=15/8` now give

\[
 \boxed{\|f\|_\infty=1,\quad \deg f=2,\quad
 \operatorname{Inf}_1(f)={6(15/8)\over5}={9\over4}>2.}          \tag{10}
\]

Filmus--Hatami--Keller--Lifshitz explicitly conjectured
`Inf_1(f)<=degree(f)` for every bounded real function after proving the
general `d^2` bound and a homogeneous `O(d log d)` bound.  Equation (10) is a
hand-checkable counterexample to that literal conjecture, already within the
dense equimodular quadratic class.  This note makes no priority claim about
the counterexample.

There is no normalization mismatch here.  In Section 2 of their paper they
write

\[
 f_i(x)={f(x)-f(x^{\oplus i})\over2}
        =x_i{\partial f\over\partial x_i}(x),
 \qquad \operatorname{Inf}_i[f]=\|f_i\|_1,
\]

and sum these quantities over `i`.  This is exactly the derivative and the
uniform probability normalization used above.  For (9), each of the six
summands is

\[
 \mathbb E|D_i(Q_A/5)|={1\over5}\mathbb E
 |\varepsilon_1+\cdots+\varepsilon_5|={1\over5}{15\over8}={3\over8},
\]

so their total influence is `6(3/8)=9/4`, whereas their conjectured ceiling
for this homogeneous degree-two function is `2`.

The distinction between pointwise and average gradient is important.  Their
published quadratic example makes the maximum pointwise gradient approach
four, but does not refute the conjecture about **total average** `L1`
influence.  Equation (10) does.

## 4. Why the established theories do not imply convergence of `K_n`

### General bounded-polynomial influence

The theorem `Inf_1(f)<=d^2` gives only `K_n<=4`, hence only

\[
 m_n\ge {n\mu_{n-1}\over4}
 \sim {1\over\sqrt{8\pi}}n^{3/2},
\]

which is weaker than the repository's universal `1/pi` lower constant.
The unrestricted sharp constants over **all** degree-two polynomials do have
a limit for an elementary monotonicity reason: adding a dummy variable nests
those classes, and the `d^2` theorem bounds their suprema.  The equimodular classes
in (2) are not nested: dummy variables, restrictions, and direct sums create
zero or unequal coefficients.  Therefore that monotonicity says nothing
about convergence of `K_n`.

### Boolean Sidon and Bohnenblust--Hille theory

Volberg's Section 6 introduces the analogous extremal `T_(n,d)` for
homogeneous Boolean polynomials whose degree-`d` Fourier coefficients have
constant modulus.  Its real square-free degree-two subproblem is exactly
`m_n`.  The Bohnenblust--Hille estimate gives the correct power law

\[
 T_{n,2}\gtrsim \binom n2^{3/4},
\]

but its dimension-free constant is not an asymptotic leading constant.
Likewise, the homogeneous Boolean-radius theorem of
Defant--Mastylo--Perez determines the scale for the larger arbitrary-
coefficient class, not the minimum over real complete equimodular signings.
Neither result supplies an amplification law, monotonicity, or a limiting
variational formula for (2).

### Grothendieck and decoupling

The influence identity also reads

\[
 n\mu_{n-1}=\mathbb E_x\|Ax\|_1
 \le\|A\|_{\infty\to1}.
\]

For symmetric `A`, polarization with
`u=(x+y)/2`, `v=(x-y)/2` gives

\[
 |x^TAy|=2|Q_A(u)-Q_A(v)|\le4M(A).                              \tag{11}
\]

Thus the elementary polarization estimate already reproduces `K_n<=4`.
Grothendieck and Rademacher-chaos decoupling replace (11) by comparisons with
other norms up to universal constants.  Those results establish the
`Theta(n^(3/2))` scale for complete unit weights, but their fixed constant
loss does not determine the leading extremal constant and gives no
cross-order inequality.

### Gaussian-chaos moment comparison

For an even quadratic cube function,

\[
 2^{-(n-1)/p}\|Q_A\|_\infty\le\|Q_A\|_p\le\|Q_A\|_\infty.     \tag{12}
\]

The lower estimate uses only the unavoidable antipodal pair of maximizers;
without quantitative information about the other states, cardinality alone
cannot improve its factor.  Hence a fixed-`p` invariance principle or
Rademacher/Gaussian-chaos comparison cannot uniformly recover `M(A)`; even
`p=Theta(n)` retains a fixed multiplicative loss in (12).  Recovering the
leading `L-infinity` constant by this route requires `p/n -> infinity`, outside
the fixed-moment regime of the standard chaos comparison theorems.  The
leading bivariate-Gaussian surrogate for the best-response fields similarly
produces an arcsine expression involving `(A^2)_(ij)/(n-1)` only at the level
of pair correlations; conference matrices have all these off-diagonal
correlations equal to zero.

## 5. Audited conclusion

The influence formulation is exact and useful nomenclature:

\[
 \text{MO convergence}\quad\Longleftrightarrow\quad
 \text{convergence of the real complete equimodular constants }K_n.
\]

It also yields the new exact best-response target (7) and the finite theorem
(10).  But none of the currently available degree-two bounded-polynomial,
Sidon/Bohnenblust--Hille, Grothendieck, or Gaussian-chaos results proves that
`K_n` converges.  A theorem that would genuinely advance the original
problem must add one of the missing pieces:

1. an asymptotically lossless operation relating the equimodular classes at
   different orders (equivalently, the multiplier-two/three gate already in
   the repository); or
2. a compact limiting object that retains `n^(3/2)` fluctuation information,
   together with continuity of the extremal norm; or
3. new `A`-dependent control of the best-response joint law in (7) strong
   enough to imply one of those cross-order statements.

General low-degree calculus alone supplies none of the three.

## Primary references checked

- Y. Filmus, H. Hatami, N. Keller, N. Lifshitz,
  [On the sum of the L1 influences of bounded functions](https://arxiv.org/abs/1404.3396),
  especially the conjecture in the Introduction, the exact convention in
  Section 2, and the bounds in Sections 3--4.
- A. Volberg,
  [An estimate of Sidon constant for complex polynomials with unimodular coefficients](https://arxiv.org/abs/2205.04936),
  especially Sections 4 and 6.
- A. Defant, M. Mastylo, A. Perez,
  [Bohr's phenomenon for functions on the Boolean cube](https://arxiv.org/abs/1707.09186),
  especially Theorem 3.1.
- A. Eskenazis, P. Ivanisvili,
  [Polynomial inequalities on the Hamming cube](https://arxiv.org/abs/1902.02406).
- S. Astashkin, K. Lykov,
  [Random unconditional convergence of Rademacher chaos in L-infinity and
  sharp estimates for discrepancy of weighted graphs and hypergraphs](https://arxiv.org/abs/2412.20107).
- P. Ivanisvili, T. Tkocz,
  [Comparison of moments of Rademacher chaoses](https://arxiv.org/abs/1807.04358).
