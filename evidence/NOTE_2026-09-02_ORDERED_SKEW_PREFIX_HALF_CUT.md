# Ordered skew: prefix variation and the one-sided cut-flip obstruction

**Status:** proved two exact all-orders identities.  The ordered construction
does not by itself close the multiplier-two diamond.  It converts that
diamond into a uniform upper bound on structured neighbors of an optimal
signing, whereas global minimality supplies only the reverse lower bound.
No finite-order search is used here.

Let `A` be a symmetric zero-diagonal signing of order `n`, fix a linear order

\[
 v_1\prec v_2\prec\cdots\prec v_n,
\]

and orient the same edge signs by

\[
 R^\prec_{v_i v_j}=A_{v_i v_j}\quad(i<j),\qquad
 R^\prec_{v_j v_i}=-A_{v_i v_j}.                 \tag{1}
\]

Thus `R^prec` is a skew signing.  Write `A^E` for the symmetric signing
obtained by reversing the signs of `A` on an unordered edge set `E`.

## 1. Exact prefix-flip and summation-by-parts identities

Fix a Boolean state `y`.  Let `P_k={v_1,...,v_k}`, let `y^(k)` be obtained
from `y` by reversing the signs on `P_k`, and put

\[
 q_k=Q_A(y^{(k)})\qquad(0\le k\le n).
\]

Since `y^(n)=-y` and the form is quadratic, `q_n=q_0`.  At the `k`th step,
direct differentiation gives

\[
 {q_{k-1}-q_k\over2}
 =y_{v_k}(R^\prec y)_{v_k}.                       \tag{2}
\]

Consequently, for another Boolean state `x` and
`z_i=x_i y_i`,

\[
 \boxed{
 x^T R^\prec y
 ={1\over2}\sum_{k=1}^n z_{v_k}(q_{k-1}-q_k).}   \tag{3}
\]

Define the signed `y`-energy of the prefix cut by

\[
 C_k(y)=\sum_{\substack{i\in P_k\\j\notin P_k}}
          A_{ij}y_i y_j={q_0-q_k\over2},
 \qquad C_0=C_n=0.                                \tag{4}
\]

Abel summation in (3) gives the equivalent identity

\[
 \boxed{
 x^T R^\prec y
 =\sum_{k=1}^{n-1}(z_{v_k}-z_{v_{k+1}})C_k(y).}  \tag{5}
\]

In particular, maximizing over `x` (equivalently, choosing every `z_i`
independently) yields the exact total-variation formula

\[
 \boxed{
 \max_x|x^T R^\prec y|
 ={1\over2}\sum_{k=1}^n|q_k-q_{k-1}|.}           \tag{6}
\]

Thus the ordered skew is not controlled by the range of the prefix energies;
it is controlled by their total variation.  Even if `|q_k|<=M` for every
`k`, range information alone permits the right side of (6) to be as large as
`nM`.  A useful proof must control oscillation, not merely reuse the endpoint
bound.

There is also a cut form of (5).  If `S={i:x_i=-y_i}` and
`s_k=1_(v_k in S)`, split the boundary of `S` into the two directional parts

\[
\begin{aligned}
 E^+_\prec(S)&=\{\{v_i,v_j\}:i<j, v_i\in S, v_j\notin S\},\\
 E^-_\prec(S)&=\{\{v_i,v_j\}:i<j, v_i\notin S, v_j\in S\}.
\end{aligned}                                                     \tag{7}
\]

Writing `F_y` and `B_y` for the sums of `A_ij y_i y_j` on `E^+` and
`E^-`, respectively, (5) says

\[
 F_y-B_y=\sum_{k=1}^{n-1}(s_k-s_{k+1})C_k(y),
 \qquad x^T R^\prec y=2(B_y-F_y).                \tag{8}
\]

For a prefix or suffix `S`, one directional part is empty.  The difficulty
starts exactly when membership in `S` alternates along the order.

## 2. Exact half-cut norm identity

Let `I_y` be the contribution to `Q_A(y)` of edges not crossing `S`.  Then

\[
 Q_A(y)=I_y+F_y+B_y,\qquad
 Q_A(x)=I_y-F_y-B_y,\qquad
 x^T R^\prec y=2(B_y-F_y).                       \tag{9}
\]

The elementary identity
`|u|+|v|=max(|u+v|,|u-v|)` therefore gives

\[
\begin{aligned}
 {1\over2}\bigl(&|Q_A(x)+Q_A(y)|+|x^T R^\prec y|\bigr)\\
 &=|I_y|+|B_y-F_y|\\
 &=\max\bigl\{
 |Q_{A^{E^+_\prec(S)}}(y)|,
 |Q_{A^{E^-_\prec(S)}}(y)|
 \bigr\}.                                        \tag{10}
\end{aligned}
\]

The two one-sided cut flips in (10) are switching-equivalent: switching all
vertices of `S` reverses the full cut and takes
`A^(E^+_prec(S))` to `A^(E^-_prec(S))`.  They consequently have the same
`Phi` norm.  Since every pair `(x,y)` is represented by a unique choice of
`(S,y)`, taking the maximum proves the exact all-orders identity

\[
 \boxed{
 {1\over2}\max_{x,y}
 \bigl(|Q_A(x)+Q_A(y)|+|x^T R^\prec y|\bigr)
 =\max_{S\subseteq[n]}\Phi\!\left(A^{E^+_\prec(S)}\right).}     \tag{11}
\]

This is precisely the four-state norm `K(A,R)` of Proposition 6.5 for the
ordered choice (1).  Hence this construction proves the multiplier-two
diamond, with its required Dini error, if and only if one can choose an order
such that

\[
 \max_S\Phi\!\left(A^{E^+_\prec(S)}\right)
 \le\sqrt2\,M+{1\over2}n^{3/2}\Omega(n),
 \qquad M=\Phi(A).                               \tag{12}
\]

For `S=P_k`, `E^+_prec(S)` is the full cut, so the corresponding signing is
a vertex switching of `A` and its norm is exactly `M`.  Equation (12) asks
whether the same near-`sqrt(2)` stability holds uniformly for every
**one-sided** cut flip.

## 3. What global minimality does, and does not, provide

Suppose now that `A` is globally optimal, so `M=Phi(A)=m_n`.  Every
one-sided cut flip is another signing of order `n`; global minimality gives

\[
 \Phi\!\left(A^{E^+_\prec(S)}\right)\ge M.        \tag{13}
\]

The inequality needed in (12) is an upper bound.  Thus the raw variational
consequence of global minimality has exactly the wrong direction.  A large
value of a half-cut neighbor is entirely compatible with `A` being a global
minimizer.

The prefix formulation makes the same obstruction transparent.  If `y` is
a positive maximizer, switching by `y` puts `Q_A(1)=M`, and every cut energy
satisfies

\[
 0\le C(S)={M-Q_A(1-2\mathbf1_S)\over2}\le M.     \tag{14}
\]

This cut positivity follows from maximality of the state `y`; it does not use
global minimality among coefficient signings.  Moreover, the matrices
obtained by successively switching the prefixes `P_k` all have norm exactly
`M` for every `A`, optimal or not.  Applying global minimality along that
prefix chain is therefore tautological and gives no control of the total
variation in (6).  Moving instead to the half-cut matrices in (11) invokes
global minimality only through the reverse inequality (13).

Switching the input representative does not add another parameter to the
construction.  For every diagonal sign matrix `D`,

\[
 R^\prec(DAD)=D R^\prec(A)D,                      \tag{15}
\]

and the simultaneous change `(A,R,x,y) -> (DAD,DRD,Dx,Dy)` leaves the
objective invariant.  Permutations amount exactly to choosing the linear
order `prec`; simultaneous switchings are gauge.

## 4. Verdict

The ordered recipe has a rigorous and potentially useful reformulation, but
not a closure from minimality alone:

\[
 \text{ordered-skew diamond}
 \quad\Longleftrightarrow\quad
 \text{uniform one-sided cut-flip stability (12)}.
\]

This is finer than the already known undirected hereditary cut inequalities.
Those inequalities control `F_y+B_y`; the ordered skew depends on
`F_y-B_y`.  Equivalently, they bound the range of the prefix path, while the
ordered construction measures its total variation.  Any successful use of
global optimality must therefore prove a new structural theorem about global
minimizers -- namely an order for which all the structured neighbors in
(12) have norm at most `sqrt(2)M+o_Dini(n^(3/2))`.  The definition of `m_n`
and its immediate switching consequences do not supply that theorem.

The original MO limit and the multiplier-two ray remain open.
