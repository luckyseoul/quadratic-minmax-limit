# The square-circle operator and exact second-shell cancellation

Date: 2026-08-25. Proposition 15.634. This note proves the spectrum suggested
by the exact `p=5,7,11,13` shell computations. It diagonalizes the complete
second dual shell, but it does not control later shells and does not prove R1.

Let `X=P^1(F_{p^2})`, `n=p^2+1`, and let `B` be the PSL orbit of square
`F_p`-sublines. Its point--circle incidence matrix `M` has

\[
 |B|={p(p^2+1)\over2},\qquad
 M^TM={p^2-1\over2}I+{p+1\over2}J.                 \tag{1}
\]

Thus `M` has rank `n`. Two distinct circles in `B` meet in zero, one, or two
points. For a fixed circle, the respective valencies are

\[
 k_0={p(p-1)(p-3)\over4},\quad
 k_1=p^2-1,\quad
 k_2={p(p^2-1)\over4}.                             \tag{2}
\]

## The intersection-graph identity

Let `A` be adjacency by two-point intersection. Fix two circles `S,T` and
put `j=|S\cap T|`. Reducing by `PSL(2,p^2)` to the standard subline and the
three standard relative positions, the remaining circle is parametrized by
one norm equation. The nonzero discriminants occur in opposite pairs; the
elementary quadratic-character sum `sum_x eta(x(x+a))=-1` gives

\[
 \#\{U:|U\cap S|=|U\cap T|=2\}=
 \begin{cases}
 (p-1)^2(p+1)/8,&j=0,\\
 p(p^2-1)/8,&j=1,\\
 (p^3+p^2-9p-1)/8,&j=2.
 \end{cases}                                      \tag{3}
\]

The diagonal count is `k_2`. Equations (2)--(3) are exactly the matrix
identity

\[
 \boxed{A^2+pA={p^2-1\over8}MM^T
              +{(p-1)^2(p+1)\over8}J.}             \tag{4}
\]

On the point-incidence module, direct counting gives

\[
 AM={p^2-1\over4}J+{(p-1)^2\over4}M,
\]

with the first coefficient interpreted after separating the incident and
nonincident counts. Hence `A` has eigenvalues `k_2` on constants and
`(p-1)^2/4` on the nonconstant part of `col(M)`. On `ker(M^T)`, (4) becomes
`A(A+pI)=0`. Taking the trace fixes the remaining multiplicities:

\[
 \operatorname{Spec}(A)=
 \left\{
 k_2^1,\ \left({(p-1)^2\over4}\right)^{n-1},\
 (-p)^{n(p-1)/4},\ 0^{n(p-3)/4}
 \right\}.                                       \tag{5}
\]

## Projected signed complements

For each square circle choose the signed complement `w_S` satisfying
`Cw_S=pw_S`. Explicitly:

- for a square affine line `c+F_p g`, use
  `w_infinity=0` and `w_x=eta(det(g,x-c))`;
- for a finite circle `N(x-c)=r`, where
  `eta(r)=-eta(-1)`, use
  `w_infinity=1` and `w_x=eta(N(x-c)-r)`.

The same quadratic-character sum proves `Cw_S=pw_S`. For two square circles,

\[
 |w_S^Tw_T|=
 \begin{cases}
 2p,&|S\cap T|=0,\\
 p,&|S\cap T|=1,\\
 0,&|S\cap T|=2,
 \end{cases}\qquad \|w_S\|^2=p(p-1).              \tag{6}
\]

Let `b_S` be the projection of `w_Sw_S^T` to
`Z={W:PWP=W, diag(W)=0}`. The normal tensors `(Pe_i)(Pe_i)^T` have Gram
matrix

\[
 {p^2-1\over4p^2}I+{1\over4p^2}J,
\]

whose inverse is `4p^2 I/(p^2-1)-2J/(p^2-1)`. Since the coordinate-square
vector of `w_S` is the complement indicator of `S`, (6) gives every entry
of `G=(<b_S,b_T>)`. Substitution of (1), (4), and (5) yields

\[
 \boxed{\operatorname{Spec}(G)=
 0^n,\quad [p^3(p-1)]^{n(p-1)/4},\quad
 [p^3(p+1)]^{n(p-3)/4}.}                          \tag{7}
\]

This also proves that the only linear relations among the projected circle
tensors are the `n` point-incidence relations.

## Complete second harmonic shell

Proposition 15.633 writes the square-circle contribution as

\[
 {1\over8p^4}\sum_S(w_S^TWw_S)^2
 -{(p-1)^2\over4p(d+2)}\|W\|_F^2,
\]

and gives the point-pair scalar. Combining it with (7), the complete signed
norm-`(p-1)/p` shell has three eigenvalues on `Z`:

\[
\begin{aligned}
 \lambda_0&=-{(p+2)(p^2-4p+1)\over4p(p^2+5)},\\
 \lambda_-&=-{p^3-3p^2-19p+9\over8p(p^2+5)},\\
 \lambda_+&=-{p^3-5p^2-19p-1\over8p(p^2+5)}.
\end{aligned}                                     \tag{8}
\]

Their multiplicities are respectively

\[
 {n(p-1)(p-3)\over8},\qquad {n(p-1)\over4},\qquad
 {n(p-3)\over4}.                                  \tag{9}
\]

The three numerators in (8) are positive from `p=11` onward (the first is
already positive at `p=5`, the second at `p=7`), and the relevant cubics are
strictly increasing thereafter. Therefore:

\[
 \boxed{\text{The complete second harmonic shadow shell is negative
 definite for every odd prime }p\ge11.}            \tag{10}
\]

Its scaled norm `2(p-1)` is even, so the radial phase of Proposition 15.631
does not change this sign. The first dual shell is positive after its odd
phase, while the second shell cancels it in every channel. This rules out a
first-shell-only positivity argument. Higher shells or a genuinely
multi-scale modular/theta inequality are still required for R1.

More sharply, the first transformed shell has scalar coefficient
`1/(4(p^2+5))` and Gaussian weight `exp(-pi/(8t))`. The most negative
second-shell channel is `lambda_0`, with weight
`exp(-pi(p-1)/(4pt))`. Consequently the truncation to exactly these two
shells is positive semidefinite only if

\[
 \exp\!\left(-{\pi(p-2)\over8pt}\right)
 \le {p\over(p+2)(p^2-4p+1)},
\]

or equivalently

\[
 t\le {\pi(p-2)\over
  8p\log((p+2)(p^2-4p+1)/p)}.                    \tag{11}
\]

Thus even the two-shell truncation forces the single-Gaussian parameter
into an `O(1/log p)` window. This is a quantitative obstruction, not a tail
estimate: omitted shells can still change the sign.

Exact construction audits at `p=3,5,7,11,13` verify (1)--(7) entrywise. The
proof above is uniform; those computations are checks, not interpolation.
