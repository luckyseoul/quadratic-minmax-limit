# All-centered-radii zonotope audit

**Status:** exact aggregate formula and rigorous method barrier; residual
(ii) remains open.

This note tests the proposed zonotope covector built from all centered
square norm circles at the lower branch-B endpoint

\[
 p=4r+1,\qquad r\ge13.
\]

It uses no graph census.  All calculations are over `F_p` or
`F_{p^2}`.

## The summed circle-pair tensor

Let `eta` be the quadratic character of `F_p`, put `N=N_{p^2/p}`, and let

\[
 \mathcal R_-:=\{R\in\mathbb F_p^\times:\eta(R)=-1\}.
\]

The centered square circles

\[
 \Gamma_R=\{z:N(z)=R\},\qquad R\in\mathcal R_-,                 \tag{1}
\]

are disjoint and partition

\[
 D_0=\{z:\eta(N(z))=-1\}.                                     \tag{2}
\]

There are `2r` circles, each of size `p+1`.

Undo the common diagonal gauge, so the finite Paley edge sign is

\[
 \kappa(u,v)=\eta(N(u-v)).                                    \tag{3}
\]

For `Gamma_R`, write its complement and sparse words as

\[
 w_R(z)=\eta(N(z)-R),\qquad c_R|_{\Gamma_R}\in\{\pm1\},
\]

with disjoint supports, and put `y_{R,+}=w_R+c_R`,
`y_{R,-}=w_R-c_R`.  The gauge-invariant half of their two edge features is

\[
 g_R(u,v):={\kappa(u,v)\over2}
 \bigl(y_{R,+}(u)y_{R,+}(v)+y_{R,-}(u)y_{R,-}(v)\bigr)
 =\kappa(u,v)(w_R(u)w_R(v)+c_R(u)c_R(v)).             \tag{4}
\]

The oriented sparse circle is a switched positive clique, so

\[
 \kappa(u,v)c_R(u)c_R(v)=1\qquad(u,v\in\Gamma_R).              \tag{5}
\]

For a finite edge put

\[
 A=N(u),\qquad B=N(v),\qquad
 T(A,B)=\sum_{t\in\mathbb F_p}\eta\bigl(t(A-t)(B-t)\bigr),     \tag{6}
\]

and let `Q(A,B)=p-1` if `A=B`, and `Q(A,B)=-1` otherwise.  Expanding

\[
 \mathbf1_{\mathcal R_-}(t)
 ={1-\eta(t)-\mathbf1_{t=0}\over2}
\]

gives the exact aggregate

\[
 \boxed{
 G(u,v):=\sum_{R\in\mathcal R_-}g_R(u,v)
 ={\kappa(u,v)\over2}
   \bigl(Q(A,B)-T(A,B)-\eta(AB)\bigr)
 +\mathbf1_{A=B,\ \eta(A)=-1}.}                              \tag{7}
\]

Thus the all-radius average is not a scalar function of the edge direction.
When `A` and `B` are distinct and nonzero, (6) is the trace of a Legendre
elliptic curve and varies with `B/A`; the direction sign (3) does not remove
it.  A Weil estimate leaves an error of order `sqrt(p)` per edge, hence an
order `p^(5/2)` error over `|H|` of order `p^2`, while the summed Max floor
below is only order `p`.  Consequently a Weil bound on (6) cannot separate
the live fibre.

Before (7) could literally be used as a zonotope covector, its edge tensor
would also have to lie in `im(R^t)`.  No such row-span identity follows from
(7); the cubic trace is the surviving endpoint-norm channel not reduced to
the directional quotas.  The next two obstructions show that even granting
such an identity would not make the aggregate useful with the presently
known quotas.

## Exact quota-capacity obstruction

Two specializations of (7) are elementary.

First suppose `A=B`.  From (6),

\[
 T(A,A)=\begin{cases}-\eta(A),&A\ne0,\\0,&A=0.\end{cases}
\]

If the edge direction is hard, `kappa(u,v)=+1`, equation (7) gives the
single value

\[
 \boxed{G(u,v)=2r\qquad(A=B,\ \kappa(u,v)=+1).}                \tag{8}
\]

For every hard direction there are exactly `p` such edges: with fixed
half-difference `delta`, equality of endpoint norms is the one linear
condition

\[
 N(a+\delta)=N(a-\delta)
 \quad\Longleftrightarrow\quad
 \operatorname{Tr}(a\delta^p)=0                               \tag{9}
\]

on the midpoint `a`.

Second take `u=0`.  If `B=N(v)` is nonsquare, then

\[
 T(0,B)=-\eta(B),\qquad \kappa(0,v)=\eta(B)=-1,
\]

and (7) gives

\[
 \boxed{G(0,v)=1\qquad(\eta(N(v))=-1).}                       \tag{10}
\]

Every opposite direction contains `p-1` distinct center edges of this
form.

At the lower endpoint there are `2r+1` directions of each type, with

\[
 P_U\in\{r+2,r+3\}\quad(U\text{ hard}),\qquad
 P_U=r\quad(U\text{ opposite}),                               \tag{11}
\]

and hence

\[
 H_+:=\sum_{U\ \mathrm{hard}}P_U=2r^2+5r+5,qquad
 H_-:=\sum_{U\ \mathrm{opposite}}P_U=2r^2+r.                    \tag{12}
\]

The capacities `p` in (9) and `p-1` in (10) exceed every quota in (11).
Choosing the required number of (8) edges in each hard direction and (10)
edges in each opposite direction produces one simple edge set: edges from
different directions cannot coincide.  It has the exact quotas and summed
pair score

\[
 \sum_{R\in\mathcal R_-}
   \bigl(S_H(y_{R,+})+S_H(y_{R,-})\bigr)
 =2\sum_{e\in H}G(e)
 =8r^3+24r^2+22r.                                             \tag{13}
\]

Even if both members of every pair had positive boundary phase, their
strongest combined floor would be `10`; over all `2r` radii this asks only
for `20r`.  The gap

\[
 (13)-20r=2r(4r^2+12r+1)>0.                                  \tag{14}
\]

Thus the projection to the all-radius score and the exact individual
direction quotas already contains simple points satisfying the strongest
summed Max floor.  No zonotope covector using only this aggregate and the
parallel-count rows can separate the branch.

## The partition of `D_0` is itself a countermechanism

The small surplus `s=7` or `9` does not make the all-radius sum negative.
There is an exact obstruction already at surplus zero and boundary `D_0`.

For one circle, the sparse word has equally many signs.  Indeed its infinity
coordinate is zero, so the infinity row of its `+p` eigenvector equation
gives

\[
 \sum_{z\in\Gamma_R}c_R(z)=0.
\]

Therefore each sign class has size `(p+1)/2=2r+1`.  By (5), hard chords join
equal signs and opposite chords join unequal signs.  Choose one opposite
edge between the sign classes, then match the remaining `2r` vertices
inside each class.  This is a perfect matching of `Gamma_R` with `2r` hard
edges and one opposite edge.

For `A=B` nonsquare, (7) gives

\[
 G(e)=2r\quad(e\text{ hard}),\qquad
 G(e)=-2r+2\quad(e\text{ opposite}).                           \tag{15}
\]

Hence this matching has aggregate half-score

\[
 2r(2r)+(-2r+2)=4r^2-2r+2.                                   \tag{16}
\]

Taking the disjoint union over all `2r` norm circles produces a perfect
matching of `D_0`, with boundary exactly `D_0` and surplus zero, whose full
summed pair score is

\[
 2(2r)(4r^2-2r+2)=16r^3-8r^2+8r.                             \tag{17}
\]

This is already much larger than `20r`.  Thus partitioning the base boundary
into centered circles supplies high positive capacity rather than an
aggregate contradiction.  The four near-pencil line toggles and the
`s=7/9` core may still matter to a covector retaining their exact locations,
but their information is lost in the all-radius sum.

## Verdict

The all-centered-radii tensor has the exact formula (7), but it retains a
nonconstant elliptic trace.  More decisively, (8)--(14) give exact simple
quota capacity, and (15)--(17) give a surplus-zero matching on the circle
partition of `D_0` with a very large admissible aggregate score.  Therefore
neither the raw all-radius sum nor that sum augmented only by hard/opposite
parallel quotas yields a zonotope separator.  A viable covector must retain
the two replacement-line locations or individual midpoint/fibre-pair rows;
averaging over the radii destroys precisely that information.
