# Circle geometry against the two-outlier boundary

**Date:** 2026-09-02
**Status:** proved boundary-intersection reduction; residual (ii) remains open

This note sharpens the positive-mismatch boundary conditions (28)--(29) in
`NOTE_2026-09-02_TWO_HALF_NEAR_PENCIL_REDUCTION.md`.  It is an all-prime
geometric argument, not a prime census.  It also records two uniform
exceptional mechanisms which show that boundary incidence alone cannot close
the remaining mismatch fibres.

Translate the pencil point to \(a=0\), identify the affine plane with
\(\mathbf F_{p^2}\), and write

\[
 D_0=\{z:\chi_{p^2}(z)=-1\}.
\]

Let \(U_1,U_2\) be the two square one-dimensional
\(\mathbf F_p\)-subspaces whose pencil lines were replaced, and let \(b\) be
the common point of the two replacement lines.  Thus

\[
 D=D_0\mathbin\triangle U_1\mathbin\triangle(b+U_1)
      \mathbin\triangle U_2\mathbin\triangle(b+U_2).       \tag{1}
\]

The point at infinity is not in \(D\).  Put \(s=7\) in the ordinary case
and \(s=9\) in the triple case.

For a positive-mismatch spike circle \(\Gamma\), let \(E\) be the original
three-set and \(N_\Gamma\) the mismatch set.  When \(\mu=1\),

\[
 P_\Gamma=\Gamma\setminus(E\cup N_\Gamma),\qquad
 |P_\Gamma|=p-3,
\]

and the strict boundary-cut lemma requires

\[
 |P_\Gamma\setminus D|\le 2s.                       \tag{2}
\]

The four omitted points \(E\cup N_\Gamma\) give the elementary lower bound

\[
 |P_\Gamma\setminus D|\ge |\Gamma\setminus D|-4.   \tag{3}
\]

## Square-circle normal forms

Because \(p\equiv1\pmod4\), a square Miquelian circle is one of the
following.

* An affine line \((c+u\mathbf F_p)\cup\{\infty\}\) with
  \(\chi_{p^2}(u)=1\).
* A finite norm circle

  \[
  \Gamma(c,R)=\{z:N(z-c)=R\},\qquad \eta(R)=-1,     \tag{4}
  \]

  where \(N=N_{p^2/p}\) and \(\eta\) is the quadratic character of
  \(\mathbf F_p\).

For (4), use a trace-zero \(\theta\in\mathbf F_{p^2}^{\times}\).  The
standard parametrization

\[
 t\longmapsto \lambda {t-\theta\over t+\theta},
 \qquad N(\lambda)=R,
\]

has determinant \(2\lambda\theta\).  Since
\(\eta(N(2\theta))=-1\), its determinant is square precisely when
\(\eta(R)=-1\).  This also agrees with the square-circle normal form in
Proposition 15.634.

## Affine circles through infinity

Let \(\Gamma=(c+u\mathbf F_p)\cup\{\infty\}\) be square.  If the affine
line does not pass through \(0\), then

\[
 \sum_{t\in\mathbf F_p}\eta(N(c+tu))=-1.            \tag{5}
\]

Indeed, \(N(c+tu)\) is a quadratic with square leading coefficient and
discriminant

\[
 (cu^p-c^pu)^2.
\]

The discriminant is a nonsquare in \(\mathbf F_p\) when
\(c/u\notin\mathbf F_p\), so the standard quadratic-character identity
gives (5).  There are consequently \((p+1)/2\) points of \(D_0\) on the
affine line and \((p+1)/2\) points of its projective completion outside
\(D_0\).

If \(\Gamma\) is not one of the four lines in (1), the four toggled affine
lines change at most four of its points.  Equations (2)--(3) then give

\[
 |P_\Gamma\setminus D|\ge {p+1\over2}-4-4
 ={p-15\over2}>18                                      \tag{6}
\]

for every \(p\ge53\).  If \(\Gamma=b+U_i\) is a replacement line, toggling
the whole line first leaves at least \((p+3)/2\) projective points outside
the boundary, and the two lines in the other direction can remove at most
two of them.  Hence

\[
 |P_\Gamma\setminus D|\ge {p-1\over2}-4
 ={p-9\over2}>18.                                      \tag{7}
\]

A square line through \(0\), other than \(U_1,U_2\), has no point in
\(D_0\).  The two original-line intersections at \(0\) cancel in (1), and
the two replacement-line intersections contribute at most two boundary
points.  It is therefore even farther from (2).

There are exactly two affine-line exceptions:

\[
 \boxed{\Gamma=U_1\cup\{\infty\}\quad\hbox{or}\quad
        \Gamma=U_2\cup\{\infty\}.}                 \tag{8}
\]

For example, on \(U_1\), the \(U_1\) toggle inserts all \(p\) affine
points, \(b+U_1\) is disjoint, and \(U_2\) and \(b+U_2\) remove their two
distinct intersection points.  Thus

\[
 |D\cap\Gamma|=p-2,
 \qquad |\Gamma\setminus D|=3.                     \tag{9}
\]

Consequently every other square affine circle is excluded by (2) in both
two-outlier geometries for every live prime \(p\ge53\), while (8) survives
the boundary-cardinality test uniformly.

## Finite circles centered at the pencil point

Let \(\Gamma=\Gamma(0,R)\).  By (4), \(\eta(R)=-1\), and hence

\[
 \boxed{\Gamma\subset D_0.}                        \tag{10}
\]

Moreover, \(\Gamma\) misses \(U_1,U_2\): on a square pencil line
\(N(tu)=t^2N(u)\) has square character and cannot equal \(R\).  Each of the
two replacement lines meets \(\Gamma\) in at most two points.  Therefore

\[
 \boxed{|\Gamma\setminus D|\le4.}                  \tag{11}
\]

This is a second uniform exceptional family.  At the level of (2) alone,
the four points \(E\cup N_\Gamma\) can cover all of
\(\Gamma\setminus D\), leaving \(P_\Gamma\subset D\).  This is only a
countermechanism to a boundary-cardinality proof: it does not construct a
Boolean three-spike completion or a common graph \(H\).

## Noncentral finite circles

Set

\[
 S(\Gamma)=\sum_{z\in\Gamma}\eta(N(z)).             \tag{12}
\]

Suppose first that \(0\notin\Gamma(c,R)\) and \(c\ne0\).  On the conic
\(N(z-c)=R\), the double cover

\[
 Y^2=N(z)                                           \tag{13}
\]

has four distinct geometric branch points: the two finite intersections
with the isotropic lines \(N(z)=0\), and the two points at infinity.  Thus
its smooth projective model has genus one.  Every
\(\mathbf F_p\)-point of the base conic is affine, and fibre counting in
(13) gives

\[
 \#(13)(\mathbf F_p)=p+1+S(\Gamma).
\]

Hasse's bound therefore yields

\[
 |S(\Gamma)|\le2\sqrt p.                            \tag{14}
\]

If \(0\in\Gamma(c,R)\) and \(c\ne0\), then \(R=N(c)\).  Every projective
direction except the tangent at \(0\) supplies exactly one further point
of \(\Gamma\), and the character of that point is the character of the
direction.  The tangent direction has character
\(-\eta(N(c))=1\), because the circle is square and
\(\eta(N(c))=-1\).  Removing that one square direction gives the exact
identity

\[
 S(\Gamma)=-1,qquad
 |\Gamma\setminus D_0|={p+1\over2}.                \tag{15}
\]

For a noncentral finite circle not through \(0\), (14) gives

\[
 |\Gamma\setminus D_0|\ge {p+1-2\sqrt p\over2}.
\]

Each of the four lines in (1) meets \(\Gamma\) in at most two points.
After allowing all eight boundary toggles and the four points omitted from
\(P_\Gamma\),

\[
 \boxed{|P_\Gamma\setminus D|
 \ge {p-23\over2}-\sqrt p.}                        \tag{16}
\]

The right side is larger than \(14\) once
\(p>(1+\sqrt{52})^2\), and larger than \(18\) once
\(p>(1+\sqrt{60})^2\).  In the live congruence class this proves:

\[
\begin{array}{c|c|c}
\text{geometry}&\text{all noncentral, non-through-0 finite circles excluded}&
\text{primes not decided by (16)}\\ \hline
\text{ordinary }(s=7)&p\ge73&p=53,61\\
\text{triple }(s=9)&p\ge89&p=53,61,73
\end{array}                                         \tag{17}
\]

For a noncentral circle through \(0\), the two original lines meet it at
\(0\) and at most one other point each; the two copies of \(0\) cancel in
(1).  The replacement lines contribute at most four more points.  Using
(15), at most six effective toggles, and then (3), gives

\[
 |P_\Gamma\setminus D|\ge {p+1\over2}-6-4
 ={p-19\over2}.                                    \tag{18}
\]

Thus these circles are excluded in the ordinary geometry for every
\(p\ge53\), and in the triple geometry for every \(p\ge61\).  Only the
triple \(p=53\) row is not decided by (18).

## What remains for positive mismatch

For \(\mu=1\), boundary geometry has reduced the infinite family to:

* the two removed pencil lines (8), for every live prime;
* the centered norm circles (10), for every live prime;
* noncentral, non-through-pencil finite circles only at \(p=53,61\) in the
  ordinary case and at \(p=53,61,73\) in the triple case;
* noncentral circles through the pencil point only for the triple
  \(p=53\) case.

The first two bullets are genuine uniform boundary-cardinality exceptions,
so no refinement which uses only \(|D\cap P_\Gamma|\) can close the branch.
One must use the circle signing, the boundary phase, another Boolean alias,
or more of the complete Max cut box.

For \(\mu=2\), the necessary condition is

\[
 |D\cap P_\Gamma|\ge {p-3\over2}-2s,
 \qquad P_\Gamma=\Gamma\setminus(E\cup N_\Gamma).  \tag{19}
\]

The same uniform exceptions show why (19) is not by itself an exclusion.
For a centered circle, (11) and the removal of five points give

\[
 |D\cap P_\Gamma|\ge p-8,
\]

and for a removed pencil line, (9) gives
\(|D\cap P_\Gamma|\ge p-7\).  Both are stronger than (19) throughout the
live range.  As above, this proves survival of the boundary-incidence
test, not existence of the required signed completion.
