# Best-response calculus versus edge-flip minimality

Date: 2026-09-02

Status: **proved all-orders identities and a conference-family obstruction;
no comparison of `K_n` across dimensions follows.**  No finite census is used.

Let `A` be a symmetric zero-diagonal signing of `K_n`, put

\[
 Q_A(x)=\sum_{i<j}a_{ij}x_ix_j={1\over2}x^TAx,
 \qquad M=\max_x|Q_A(x)|,
\]

For the first-moment statements below, fix once and for all a tie value
`sgn_0(0) in {+-1}` (a separately fixed value at each coordinate would also
work), and define the synchronous best response

\[
 T_A(x)_i=\operatorname {sgn}_0((Ax)_i).                     \tag{1}
\]

The tie choice must not depend on the other input coordinates; an arbitrary
input-dependent tie rule would invalidate the exchangeability argument in
Proposition 1.  The exact interpolation
inequality involving `Q_A(x)` and `Q_A(T_Ax)` is recorded in
`NOTE_2026-09-02_ORIGINAL_MO_L1_INFLUENCE_CALCULUS.md`.  The issue here is
whether global minimality of `M` under coefficient-edge flips adds the
missing control.

## 1. The complete first-order law of the best-response map is universal

Put `r=n-1` and `mu_r=E|epsilon_1+...+epsilon_r|`.

### Proposition 1

For every signing `A` and every two distinct indices `i,j`,

\[
 \boxed{\quad
 \mathbb E_x[T_A(x)_i x_j]
   =a_{ij}{\mu_{n-1}\over n-1}.\quad}                        \tag{2}
\]

Consequently

\[
 \boxed{\quad
 \mathbb E_x\,T_A(x)^TAx=n\mu_{n-1}.\quad}                  \tag{3}
\]

**Proof.**  Fix `i` and absorb the row signs into independent Rademacher
variables `epsilon_k=a_ik x_k`, `k ne i`.  By exchangeability, the numbers

\[
 \theta=\mathbb E\left[\epsilon_j
       \operatorname {sgn}\left(\sum_{k\ne i}\epsilon_k\right)\right]
\]

are the same for all `j ne i`.  Summing them gives

\[
 (n-1)\theta
 =\mathbb E\left[\left(\sum_{k\ne i}\epsilon_k\right)
       \operatorname {sgn}\left(\sum_{k\ne i}\epsilon_k\right)\right]
 =\mu_{n-1}.
\]

This is unaffected by which fixed tie value was chosen, because the product
with the zero sum vanishes.  Since
`x_j=a_ij epsilon_j`, (2) follows.  Summing `a_ij` times (2) over ordered
pairs gives (3).  QED.

The left side of (3) is the unnormalized pointwise `L^1` gradient averaged
over the cube.  It is therefore *exactly the same for every vertex of the
coefficient hypercube*.  In particular, if `A^e` is obtained by flipping any
edge,

\[
 \mathbb E\,T_{A^e}(x)^TA^ex
 =\mathbb E\,T_A(x)^TAx.                                    \tag{4}
\]

Thus the natural averaged best-response action has identically zero discrete
variation under edge flips.  Global minimality cannot enter through its
first moment; after division by `M`, (3) is just the already-known identity
`K(A)=n mu_(n-1)/M`.

## 2. The only possible gain is genuinely joint push-forward energy

For a fixed `x`, abbreviate `T=T_A(x)` and `G=T^TAx=||Ax||_1`.  With
`u=(x+T)/2` and `v=(x-T)/2`, homogeneity gives

\[
 4Q_A(u)=Q_A(x)+Q_A(T)+G,
 \qquad
 4Q_A(v)=Q_A(x)+Q_A(T)-G.
\]

Both `u` and `v` lie in `{0,+-1}^n`; their quadratic values are conditional
expectations of Boolean values and hence have absolute value at most `M`.
Therefore

\[
 \boxed{\quad
 G+|Q_A(x)+Q_A(T_Ax)|\le4M.\quad}                            \tag{5}
\]

After averaging and normalizing `f=Q_A/M`,

\[
 K(A)+\mathbb E|f(x)+f(T_Ax)|\le4.                           \tag{6}
\]

Equation (2) determines the entire degree-one Fourier part of every output
coordinate of `T_A`, but supplies no lower bound for the second term in (6).
That term necessarily involves correlations between distinct output
coordinates and hence higher Fourier levels of the threshold functions.

There is also an exact obstruction at the extremal states used by edge-flip
arguments.  Here we make a pointwise best-response selection, distinct from
the fixed map used in Proposition 1.  If `Q_A(x)=M`, coordinatewise
maximality gives `x_i(Ax)_i>=0` for every `i`; selecting zero-field ties
toward this particular `x` gives `T=x`.  If `Q_A(x)=-M`, selecting its
zero-field ties oppositely gives `T=-x`.  In either case

\[
 G=2M,\qquad |Q_A(x)+Q_A(T)|=2M,                             \tag{7}
\]

so (5) is equality.  Exact maximizers and minimizers are fixed points or
antipodal two-cycles of the best-response dynamics, precisely where the
calculus has no slack.

## 3. What one-edge global minimality actually says

Let `e={i,j}`, put `s_e(x)=a_ij x_i x_j`, and let `A^e` denote the signing
with that edge flipped.  Then

\[
 Q_{A^e}(x)=Q_A(x)-2s_e(x).                                  \tag{8}
\]

If `A` is a global minimizer with value `M=m_n`, then
`M<=Phi(A^e)<=M+2`; parity makes `Phi(A^e)` equal to `M` or `M+2`.
More precisely, the condition `Phi(A^e)>=M` is equivalent to the existence
of a state `x` satisfying

\[
 \boxed{
 \begin{split}
 &s_e(x)=-1\quad\hbox{and}\quad Q_A(x)\ge M-2,\\
 &\hspace{27mm}\text{or}\\
 &s_e(x)=+1\quad\hbox{and}\quad Q_A(x)\le-M+2.
 \end{split}}                                                \tag{9}
\]

This follows immediately by solving
`|Q_A(x)-2s_e(x)|>=M` under `|Q_A(x)|<=M`.

Thus one-edge minimality is an outer-two-energy-layer covering condition.
It supplies witnesses, not probability mass.  When a witness is on the
outermost layer, (7) says that the best-response inequality is exactly
saturated there.  Nothing in (9) lower-bounds the uniform measure of these
layers, so it cannot be inserted into the average in (6) without a new
anti-concentration or multiplicity theorem.

## 4. Symmetric conference matrices kill the mean push-forward route

The obstruction persists one level beyond (3).

### Proposition 2

Let `C` be any symmetric conference signing of even order `n`, so
`C^2=(n-1)I`.  Its local fields are odd and never vanish.  For uniform `X`,

\[
 \boxed{\quad
 \mathbb E[T_C(X)_iT_C(X)_j]=0\quad(i\ne j),
 \qquad \mathbb E Q_C(T_CX)=0.\quad}                         \tag{10}
\]

In fact the output coordinates are pairwise unbiased and independent.

**Proof.**  Fix `i ne j`, set `a=C_ij`, and condition on every coordinate
except `X_i,X_j`.  Write

\[
 U=\sum_{k\ne i,j}C_{ik}X_k,
 \qquad V=\sum_{k\ne i,j}C_{jk}X_k.
\]

Both are even.  Averaging separately over `X_j` and `X_i` gives

\[
 \mathbb E_{X_i,X_j}[T_iT_j\mid(X_k)_{k\ne i,j}]
 =\sigma(U)\sigma(V),                                       \tag{11}
\]

where `sigma(t)=sgn(t)` for `t ne 0` and `sigma(0)=0`.
Orthogonality of rows `i,j` says

\[
 \sum_{k\ne i,j}C_{ik}C_{jk}=0.
\]

Hence among the `n-2` remaining positions the two signed rows agree and
disagree equally often.  After absorbing coefficients into the variables,

\[
 (U,V)=(S+R,S-R),
\]

where `S,R` are independent, identically distributed sums of `(n-2)/2`
Rademacher variables.  The right side of (11) is
`sigma(S^2-R^2)`.  Swapping `S` and `R` reverses its sign and preserves its
law, while ties contribute zero.  Its expectation is therefore zero.
Each `T_i` is itself unbiased, proving pairwise independence.  Summing the
first assertion against `C_ij` proves the second.  QED.

Consequently the coarser best-response bound

\[
 K(C)\le2+2\sqrt{1-|\mathbb E(Q_C(T_CX)/M)|}
\]

reduces exactly to the uninformative ceiling four on every symmetric
conference signing.  If the conference matrix is regular,
`C 1=sqrt(n-1) 1`, then `1` is also an exact maximizer and fixed point of
`T_C`; (5) is attained with equality at that state.  These statements hold
at every order at which such a matrix exists and use no small-order data.

Conference matrices are not asserted here to be global minimizers at every
order.  Rather, (10) proves that a universal best-response estimate cannot
distinguish the canonical flat-spectrum candidates.  To use global
minimality, one would need information strictly stronger than its one-edge
outer-layer witnesses (9).

## 5. Conclusion

The best-response map provides a clean diagnostic but not an all-orders
comparison:

1. its complete first-order joint law with the input is forced by
   equimodularity and is independent of the signing;
2. edge-flip global minimality only covers edges by states in the outer two
   energy layers, where exact outer-layer states saturate the calculus;
3. even the mean energy after one synchronous response is exactly zero for
   every symmetric conference matrix.

Therefore no inequality using the averaged action, degree-one response
correlations, or mean push-forward energy can compare `K_n` with `K_(n+1)`.
The surviving target is genuinely higher-order: an `A`-dependent lower bound
on the absolute joint term in (6), together with an outer-layer multiplicity
theorem derived from **global** rather than one-edge minimality.  Neither is
proved here, so the original limit remains open.
