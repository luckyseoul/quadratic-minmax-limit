# Hard-row compact-residual odd-Radon centrality

This note proves a structural fact about the **hard** rows of the balanced
branch-C compact survivor.  It is separate from the opposite-row theorem.
The atom profile is different, and none of the opposite-row line bounds is
silently reused.

Let

\[
 p=4r+3\ge31,\qquad h=(p-1)/2=2r+1,\qquad m=h-2=2r-1.
\]

A balanced branch-C hard row is

\[
       \text{one fixed unit star}+e\text{ compact atoms}.             \tag{1}
\]

The conclusion below concerns the compact residual after subtracting the
star.  It does **not** say that the whole hard-row edge chain is central.

## 1. The full hard ray has \(e\le2r-2\)

The branch-C parameter interval is

\[
 2r^2-4r-2\le t\le4r^2-2r-5,
\]

and the \(2r+2\) hard compact counts are the balanced allocation of
\(t+1\).  At the upper endpoint,

\[
 4r^2-2r-4=(2r+2)(2r-3)+2.                             \tag{2}
\]

Thus exactly two hard rows have \(e=2r-2\), the other \(2r\) have
\(e=2r-3\), and throughout the full balanced ray

\[
                         0\le e\le2r-2.                 \tag{3}
\]

For \(d\le p-2\), every summand in the degree-\(d\) moment of a translated
unit star is a polynomial in its running label of degree at most \(d\).
The field sums of all monomials of degrees \(0,\ldots,p-2\) vanish.
Consequently the star in (1) contributes zero to every odd form
\(d=3,5,\ldots,p-2\).  If all those global forms vanish, the compact
residual by itself is odd-zero.

Pair an edge with its negative and write its integer orbit difference as
\(n_E\).  The usual invariant coordinates give a word on
\(\Omega=H\times H\), where \(H\) is the set of nonzero squares,

\[
 W(U,D)=n_E D\sigma,
\]

orthogonal to all bivariate polynomials of total degree at most \(m\).
The compact-only occurrence budget is

\[
 N:=\sum_E |n_E|\le3e\le6r-6=3m-3.                    \tag{4}
\]

## 2. Couvreur leaves only line and conic branches

Choose a minimally \(m\)-linked subset of the support of a hypothetical
nonzero \(W\).  Couvreur's first three low-weight alternatives are:

1. \(m+2=h\) collinear points;
2. \(2m+2=p-3\) points on a conic, in the absence of the first case;
3. a \(3m\)-point cubic/degree-\(m\) complete intersection.

Here (4) is strictly below \(3m\), so the cubic alternative cannot occur.
It remains to exclude supports containing a line and supports containing an
irreducible conic.  A reducible conic needs a separate two-line check:
its minimally linked subset can have only \(h-1\) points on each component,
so neither component need contain the \(h\)-point line circuit.

The linked-set input is A. Couvreur, *The dual minimum distance of
arbitrary-dimensional algebraic-geometric codes*, J. Algebra 350 (2012),
Theorem 3.8 and Lemma 2.13.  The latter permits the standard algebraic-field
extension used after homogenizing the affine Cartesian points.

## 3. One and two maximal lines are impossible

The line-peeling lemma from the compact-ray audit applies because
\(N\le3m\): a support containing one full maximal line is either that line
or lies on the union of it with one further maximal line.  We now use the
hard compact-only budget (4), not the larger opposite-row budget.

On one horizontal or diagonal maximal line, the \(h\) nonzero orbit
differences represent all \(h\) projective residue classes.  Hence

\[
 \sum |n_E|\ge1+\cdots+h=\frac{h(h+1)}2>6r-6.          \tag{5}
\]

On a vertical line, all residues are one nonzero constant \(\kappa\).
Since \(N<3h\), its least absolute representative has size at most two.
Using the other integer lift even once costs at least

\[
 (h-1)+(p-1)=6r+2>6r-6,                                \tag{6}
\]

so every actual coefficient is the same integer \(\pm1\) or \(\pm2\).
For \(\pm1\), reduction modulo two and projection of label vertices by
\(x\sim-x\) is decisive: every compact atom becomes a triangle boundary
and is Eulerian, whereas a full fixed-sum matching has the two odd quotient
vertices \([0]\) and \([\sigma/2]\).  For \(\pm2\), the target needs
\(2h=4r+2\) aligned occurrences.  A compact atom has at most two aligned
occurrences on the two physical fixed-sum matchings, so all compact atoms
together supply at most

\[
 2e\le4r-4<4r+2.                                       \tag{7}
\]

For two lines from different maximal families, at least one nonvertical
component has \(h-1\) exclusive points in distinct projective classes.
This already costs

\[
 1+\cdots+(h-1)=\frac{h(h-1)}2>6r-6.                  \tag{8}
\]

For two lines of the same family, the affine-Cartesian dual factorization
leaves an affine-linear residual factor.  On a full vertical component a
nonconstant factor is injective in the other coordinate, with least
integer mass

\[
 (r+1)^2>6r-6.                                         \tag{9}
\]

If the factor is constant, both vertical components are full.  Their two
nonzero constant residues have total least absolute value at most two, so
both are units; an alternative lift costs at least \(4h-1\).  Modulo two,
the two fixed-sum graphs leave the distinct odd vertices
\([\sigma_1/2]\) and \([\sigma_2/2]\), again contradicting Euler parity.

On a full horizontal or diagonal component, one projective coefficient
class occurs at most three times: after clearing denominators the equality
is cubic.  If \(h=3q+s\), \(0\le s<3\), its exact least mass is

\[
 L_3(h)=\frac{3q(q+1)}2+s(q+1)>6r-6\qquad(h\ge15).     \tag{10}
\]

This excludes every line-containing support.

It also excludes the reducible-conic configuration that did not itself
contain a full line circuit.  A nonmaximal line meets \(\Omega\) in at most
\(r+1\) points.  Since the conic circuit has \(2h-2=4r\) points, two
nonmaximal components, or one maximal and one nonmaximal component, are too
small; both components must be distinct maximal lines.  For different
families (8) already uses only the \(h-1\) exclusive points.  For two
vertical lines, a nonconstant residual factor gives \(h-1=2r\) distinct
nonzero residues on a component, of least mass

\[
                       r(r+1)>6r-6.                    \tag{10a}
\]

For two horizontal or two diagonal lines, the cubic projective-fibre bound
on only \(h-1\) nonzero coordinates gives

\[
 L_3(h-1)>6r-6\qquad(r\ge7).                            \tag{10b}
\]

The constant two-vertical case remains full on both components and was
already excluded by quotient parity.  Thus every reducible-conic support is
excluded as well.

## 4. Every irreducible-conic branch is impossible

Suppose the support contains at least \(2m+2=p-3\) points on an
irreducible conic \(Q\).  Multiplying the test polynomials by an equation
of \(Q\) shows that the support outside \(Q\), if nonempty, is a nonzero
degree-\((m-2)\) dual word and therefore has at least \(m\) points.  But
(4) gives

\[
 |\operatorname{supp}(W)\setminus Q|
 \le(3m-3)-(2m+2)=m-5.                                 \tag{11}
\]

Thus the whole word is conic-supported.  The character-sum and tangent
classification proved in
`NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md` now applies verbatim:

\[
 U=u z^2,\qquad D=d(z-1)^2,
\]

and on \(z\in\mathbf F_p\setminus\{0,1\}\),

\[
 n(z)=\alpha+\frac{\beta}{z-1}.                        \tag{12}
\]

If \(\beta\ne0\), the exact least integer mass is at least
\(h(h-1)>N\).  If \(\beta=0\), a nonunit residue costs at least
\(2(p-2)>N\), as does one nonminimal integer lift.  The only survivor could
therefore be the constant integer word \(+1\) or \(-1\) on all \(p-2\)
target edge orbits.

The star conics \(k=\pm1\) are excluded by the same quotient-Euler parity
as above: their projected target graph has odd vertices \([0]\) and
\([1/2]\), while a sum of compact triangle boundaries is Eulerian.

For a nonstar conic, normalize the constant word to \(+1\) and use the
alignment score from the conic note.  In the nonequianharmonic case
\(q^3\ne1\), all compact atoms have score at most two except for at most two
explicit score-three types.  The repeated-candidate argument is important.
If \(T=e\), \(n_3\) is the number of score-three occurrences, and
\(A,R,O\) count aligned, reverse, and other occurrences, then

\[
 A-R=p-2,\qquad 2R+O=3e-(p-2).                         \tag{13}
\]

Non-score-three atoms contribute at least \(1+R_i\) to \(2R_i+O_i\).
If \(n_3\le2\), this gives \(2R+O\ge e-2\).  If \(n_3\ge3\), their
\(3n_3\) aligned occurrences lie on at most six target coordinates, so
\(R\ge3n_3-6\) and in fact \(2R+O\ge e\).  In both cases

\[
 p-2\le2e+2\le4r-2<4r+1=p-2,                           \tag{14}
\]

a contradiction.

In the equianharmonic case \(q^3=1\ne q\), the characteristic-uniform
eight-case table in
`NOTE_2026-09-03_EQUIANHARMONIC_COMPONENT_PACKING.md` proves that a
distinct-label compact atom has score at most two; all exceptional
characteristics are below 31.  Therefore

\[
 p-2\le2e\le4r-4<4r+1=p-2,                             \tag{15}
\]

again impossible.  Notice that the more delicate component-packing
threshold is unnecessary here: the hard residual has no all-equal atoms,
so the elementary per-compact score bound already has a strict margin.

Equations (11)--(15) exclude every irreducible-conic branch, including the
star, nonequianharmonic, and equianharmonic cases.

## 5. The exact conclusion and its boundary

Couvreur leaves no possible support for a nonzero compact-residual word.
Hence every orbit difference is zero modulo \(p\).  A single compact atom
can change one antipodal orbit difference by at most two, so

\[
 |n_E|\le2e\le4r-4<p=4r+3.                              \tag{16}
\]

Thus congruence zero is integer zero: the compact residual is centrally
symmetric as a signed integer edge chain.

The proved statement is exactly:

> For every prime \(p=4r+3\ge31\), every hard row in the balanced branch-C
> compact survivor, and every admissible labeling of its \(e\) compact
> atoms, vanishing of all odd global forms through degree \(p-2\) forces
> the compact residual to be centrally symmetric over the integers.

The fixed unit-star term remains.  No claim is made for nonzero odd global
forms, unbalanced hard allocations, the joint degree-six/eight system, a
common \(\mathbf F_p\) edge lift, the Boolean affine-box intersection, or
residual-(ii) itself.  All of those gates remain open.

The executable arithmetic/dependency certificate is
`src/e1_gmin_m4_hard_compact_odd_radon.py`; its focused replay is
`tests/test_hard_compact_odd_radon.py`.
