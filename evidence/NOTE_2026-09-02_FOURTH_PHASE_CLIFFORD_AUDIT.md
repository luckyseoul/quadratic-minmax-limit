# Fourth-phase form of the doubling diamond and the Clifford parity barrier

**Status:** exact reformulation and exact obstruction.  The phase identity is
a new formulation of the live multiplier-two target.  Exact Clifford
flattening is impossible from order three onward; approximate or
phase-restricted flattening remains open.

Let `A` be a symmetric zero-diagonal sign matrix and `R` a skew
zero-diagonal sign matrix.  Put

\[
 Q_A(x)=\sum_{i<j}A_{ij}x_ix_j
 \quad\hbox{and}\quad
 \mathcal H=A+iR.
\]

The matrix `mathcal H` is Hermitian and every off-diagonal entry is one of
`+1+i,+1-i,-1+i,-1-i`.

## 1. The exact fourth-phase identity

For Boolean `x,y`, put

\[
 u={x+y\over2},\qquad v={x-y\over2},\qquad z=u+iv.
\]

Thus `z` ranges bijectively over `{+1,-1,+i,-i}^n` as `(u,v)` ranges over
ordered signed complementary supports.  Define

\[
 I=Q_A(u)+Q_A(v)={Q_A(x)+Q_A(y)\over2},
 \qquad C=-u^TRv={x^TRy\over2}.
\]

Symmetry of `A` and skewness of `R` give

\[
 z^*\mathcal H z=2(I+C),\qquad
 \bar z^*\mathcal H\bar z=2(I-C).                  \tag{1}
\]

Since
`max(|I+C|,|I-C|)=|I|+|C|` and the fourth-phase cube is closed under
conjugation, (1) proves the global identity

\[
 \boxed{
 \max_{x,y\in\{\pm1\}^n}
 \left(|Q_A(x)+Q_A(y)|+|x^TRy|\right)
 =
 \max_{z\in\{\pm1,\pm i\}^n}|z^*(A+iR)z|.
 }                                                       \tag{2}
\]

Consequently Proposition 6.5's multiplier-two problem is exactly a
Hermitian extension problem: prescribe the real part `A` and choose all
imaginary off-diagonal signs so that the fourth-phase quadratic norm in
(2) is at most `2 sqrt(2) Phi(A)` plus the admissible error.

The full operator norm gives the sufficient estimate

\[
 \max_z|z^*\mathcal H z|\le n\|\mathcal H\|_{op}.     \tag{3}
\]

It is generally stronger than needed because (2) samples only the
fourth-phase torus.

## 2. What exact Clifford flattening would do

If, hypothetically,

\[
 A^2=(n-1)I,\qquad R^2=-(n-1)I,qquad AR+RA=0,        \tag{4}
\]

then

\[
 \mathcal H^2=2(n-1)I,qquad
 \|\mathcal H\|_{op}=\sqrt{2(n-1)}.                 \tag{5}
\]

If in addition `A` has a Boolean spectral extremizer, then
`Phi(A)=n sqrt(n-1)/2`; equations (2)--(5) give the exact desired diamond

\[
 \max_{x,y}
 \left(|Q_A(x)+Q_A(y)|+|x^TRy|\right)
 \le 2\sqrt2\Phi(A).                                \tag{6}
\]

Thus the two halves really do assemble into one complex spectral ball.
The next subsection shows why equality cannot occur for sign matrices.

## 3. Exact anticommutation is parity-impossible

For `i != j`, define

\[
 s_{ij}=A_{ij}R_{ij}.
\]

Then `s_ji=-s_ij`.  Put

\[
 p_i=\prod_{k\ne i}s_{ik}\in\{\pm1\}.
\]

Consider `(AR+RA)_ij/2` modulo two.  For a fixed
`k notin {i,j}`, the two signs

\[
 A_{ik}R_{kj},\qquad R_{ik}A_{kj}
\]

are equal precisely when `s_ik != s_jk`.  Hence, if
`(AR+RA)_ij=0`, the Hamming distance

\[
 d_{ij}=|\{k\ne i,j:s_{ik}\ne s_{jk}\}|
\]

must be even.  On the other hand

\[
 (-1)^{d_{ij}}
 =\prod_{k\ne i,j}s_{ik}s_{jk}
 =-p_ip_j.                                           \tag{7}
\]

Thus exact anticommutation would force `p_i p_j=-1` for every pair
`i != j`.  Three indices make this impossible.  Therefore

\[
 \boxed{AR+RA\ne0\quad\hbox{for every such sign pair when }n\ge3.} \tag{8}
\]

This also rules out a Hermitian conference matrix whose upper-triangular
entries are all `(A_ij+iR_ij)/sqrt(2)` at every order `n >= 3`.

There is a quantitative parity floor.  Whenever `p_i=p_j`, (7) makes
`d_ij` odd, and therefore `(AR+RA)_ij` is a nonzero multiple of two.  The
number of same-`p` unordered pairs is minimized by splitting the `p_i` as
evenly as possible.  Consequently

\[
 \|AR+RA\|_F^2\ge
 \begin{cases}
 2n(n-2),&n\text{ even},\\
 2(n-1)^2,&n\text{ odd}.
 \end{cases}                                          \tag{9}
\]

The factor two in (9) counts both orientations of each matrix entry.
This is only an order-`sqrt(n)` operator-norm floor, so it does **not**
exclude approximate anticommutation at the accuracy relevant to a Dini
error.

## 4. Full spectral flattening is value-specific

For every sign pair, independently of anticommutation,

\[
 \|A+iR\|_F^2=2n(n-1),
 \qquad
 \|A+iR\|_{op}\ge\sqrt{2(n-1)}.                    \tag{10}
\]

Therefore a proof of the doubling target using only (3) would require

\[
 n\sqrt{2(n-1)}
 \le 2\sqrt2\,m_n+o(n^{3/2}),
\]

or

\[
 m_n\ge {n\sqrt{n-1}\over2}-o(n^{3/2}).             \tag{11}
\]

Together with the known upper envelope, (11) already forces limiting value
`1/2`.  Hence the unrestricted operator norm cannot supply a genuinely
value-free subadditivity proof.  The live information in (2) is the
**restricted fourth-phase norm**, or an `A`-dependent approximate Clifford
construction strong enough on that discrete torus without controlling the
whole unit sphere.
