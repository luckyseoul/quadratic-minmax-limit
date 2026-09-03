# Hard-star antisymmetric support and the equality pencil

Date: 2026-09-03

Status: proved conditional support inequality and equality classification.
The later direction-localized Möbius construction solves the entire
antisymmetric ternary box; the coupled symmetric half and residual (ii)
remain open.

## 1. Input from hard-row compact centrality

Let \(p=4r+3\ge31\) be prime. In a balanced branch-C hard row, the compact
residual is centrally symmetric over the integers when all odd global forms
through degree \(p-2\) vanish. The only antisymmetric target left in that
row is therefore the fixed unit-star difference

\[
  \pm(S_j-S_{-j}).
\]

It is zero when \(j=0\). When \(j\ne0\), it has exactly \(2(p-2)\) nonzero
cells, every coefficient is \(\pm1\), and its squared norm is \(2(p-2)\).

Write \(J(e)=-e\) for antipodal negation of source edges. For a Boolean edge
set \(x\), let:

- \(c\) be the number of nonfixed \(J\)-orbits occupied on exactly one side;
- \(d\) be the number occupied on both sides; and
- \(f\) be the number of selected fixed antipodal edges.

Then

\[
 |H|=f+c+2d,
 \qquad
 |\{e:x_e\ne x_{-e}\}|=2c,
\]

and, for \(x^\pm=(x\pm Jx)/2\),

\[
 \|x^-\|^2={c\over2},
 \qquad
 \|x^+\|^2=|H|-{c\over2}.
\]

One active hard row partitions the source edges. Its \(2(p-2)\) nonzero
target cells consequently require \(2c\ge2(p-2)\). Thus

\[
             c\ge p-2.                                      \tag{1}
\]

This is sharper than the real least-norm estimate. If \(A\) hard rows have
nonzero centres, then the target squared norm is \(2A(p-2)\), while the
exact real least squared norm is only \(2A(p-2)/p^2<1\), even for
\(A=(p+1)/2\).

## 2. Exceptional-direction refinement

A nonfixed source orbit has midpoint \(a\ne0\) and difference class
\([\delta]\). Its image is a nonself transverse cell in an active direction
\(L\) unless \(L(a)=0\) or \(L(\delta)=0\). Let \(E\) be the sum, over the
\(c\) single source orbits, of the number of distinct active directions
among these two exceptions. Counting active nonself images gives

\[
              A(p-2)\le Ac-E.                               \tag{2}
\]

In particular,

\[
       c\ge p-2+\left\lceil {E\over A}\right\rceil.          \tag{3}
\]

Equality \(c=p-2\) forces \(E=0\), and every active row is a bijection from
the selected single source orbits to its \(p-2\) star-cell orbits.

## 3. Equality and a degree-eight resultant

Assume \(c=p-2\) and at least nine hard rows are active. Let
\(e=\{u,v\}\) and \(e'=\{u',v'\}\) be two selected source edges. In an
active row \(L\), both endpoint-square sets contain the exceptional value
\(j_L^2\). Therefore

\[
 R_{e,e'}(L)
 =\prod_{a\in\{u,v\}}\prod_{b\in\{u',v'\}}
      \bigl(L(a)^2-L(b)^2\bigr)=0.                           \tag{4}
\]

Since

\[
 L(a)^2-L(b)^2=L(a-b)L(a+b),
\]

equation (4) is a product of eight projective linear factors. A nonzero
homogeneous polynomial of degree eight on the projective line has at most
eight roots. Nine distinct active directions force one factor to vanish
identically. Hence \(e\) and \(e'\) share an endpoint modulo sign.

The projective endpoint pairs of all selected source orbits are now a
pairwise-intersecting family of two-subsets. Such a family either has one
common vertex or is contained in a triangle. A fixed projective endpoint
pair has only two physical edge-orbit lifts modulo simultaneous negation,
so the triangle alternative contains at most six source orbits. Since
\(p-2\ge29\), all selected orbits share one projective source vertex
\([P]\).

Finally, if \(L(P)^2\ne j_L^2\), the other endpoint of every selected edge
must lie on one of the two exceptional affine lines. With the projective
centre fixed this yields at most two target cell-orbits, not the required
\(p-2\) distinct cell-orbits. Consequently

\[
                    j_L^2=L(P)^2                            \tag{5}
\]

for every active row. Thus equality with at least nine active rows forces a
full projective pencil and coherent hard phases up to sign.

Rows with \(j_L=0\) are inactive and contribute no root to (4). The argument
therefore says nothing when fewer than nine hard rows have nonzero centres.

## 4. Why two opposite rows do not kill the pencil

A tempting next claim is that a full pencil cannot have zero antisymmetric
projection in two opposite directions. That claim is false without the
simultaneous active-row bijections.

Choose the two directions as the \(x\)- and \(y\)-coordinates. In translated
neighbour coordinates \(z=w+P\), take the signed five-cell path

\[
 (0,1)^+,(1,1)^-,(1,2)^+,(2,2)^-,(2,0)^+                  \tag{6}
\]

and, for \(i=1,\ldots,r-1\), the disjoint checkerboard rectangle on rows and
columns \(2i+1,2i+2\):

\[
 (a,a)^+,(a,b)^-,(b,a)^-,(b,b)^+ .                         \tag{7}
\]

There are \(5+4(r-1)=p-2\) distinct nonzero points. Every nonzero row and
column has signed sum zero; the zero row and column have sum one. Take
\(P=(2r+1,2r+1)\). Then \(P_x,P_y\ne0\), none of these points is \(0\) or
\(2P\), and the edges \(\{P,z-P\}\) give \(p-2\) distinct nonfixed
antipodal edge orbits.

For a direction \(M\), the zero fibre \(M(z)=0\) is invisible to the
antisymmetric edge projection: the edge and its antipode land in the same
self-antipodal cell and cancel. Thus (6)--(7) have exactly zero
antisymmetric Radon projection in both coordinate directions.

This construction is not a residual-(ii) counterexample: it does not
satisfy the active hard-row bijectivity constraints. It does prove that two
zero rows alone cannot be used to improve (1) to \(c\ge p-1\).

## 5. Scalar count methods at this gate

At the floor \(c=p-2\), both endpoints of the balanced branch-C ray admit
nonnegative integer ledgers

\[
 |H|=f+c+2d,
 \qquad
 P_D=f_D+c_D+2d_D,
\]

within the available fixed and nonfixed orbit capacities. The old
full-target norm bound also lies strictly below the available symmetric
half norm throughout the ray. These checks do not construct a common simple
graph; they show only that scalar total, fixed, parallel, and norm
bookkeeping cannot supply the missing contradiction.

Before the direction-localized construction was found, the equality case
reduced to a simultaneous finite-affine problem:

1. one set of \(p-2\) oriented pencil edges had to be bijective onto every
   active hard star-cell family;
2. its signed fibre sums had to vanish in every zero hard/opposite row; and
3. the unused double and fixed coordinates had to realize the prescribed
   symmetric target.

## 6. Subsequent resolution of the antisymmetric gate

The later theorem in
NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md gives, for every hard
direction and every nonzero centre, a direction-localized Möbius trade on
exactly \(p-1\) inversion orbits. Its Radon image is the required hard-star
difference in that direction and zero in every other direction. A sharp
greedy avoidance count makes the trades for all hard directions pairwise
disjoint. After the Paley column-sign transport their sum is a ternary
antisymmetric source and hence the inversion difference of a simple graph.

Thus the cross-directional antisymmetric incidence problem isolated above
is no longer live. The lower bound, exceptional-incidence refinement,
equality pencil theorem, and two-row counterexample remain valid structural
facts, but they are not needed to construct the antisymmetric lift.

The live problem is the coupled symmetric half. Every orbit used by a
Möbius trade has forced pair total \(s_e=1\); unused nonfixed orbits have
\(s_e\in\{0,2\}\), and fixed antipodal edges are separate binary variables.
No result here or in the Möbius theorem realizes those symmetric totals.
Residual (ii) therefore remains open.

## Reproduction

    PYTHONPATH=src python src/e1_gmin_m4_hard_star_antisymmetric_support.py
    PYTHONPATH=src pytest -q tests/test_hard_star_antisymmetric_support.py
