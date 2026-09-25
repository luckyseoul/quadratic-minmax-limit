# Opposite-phase stable pairs have vertexwise bipartite dominance

2026-09-25. Direct continuation of the frozen one-vertex extension route.
This note strengthens the opposite-phase branch in
`NOTE_2026-09-25_BANASZCZYK_STABLE_RECURSION.md`.

Repository-wide searches for the vertexwise statement below, including
"opposite phase cross degree internal degree", "cross-degree dominates
internal", "signed cross degree nonnegative stable opposite", and
"opposite stable states vertexwise cut", found no earlier occurrence.
The aggregate cut identity from the Banaszczyk recursion is prior work.
The new content here is the per-vertex domination and its induced
`l1`-field consequences.

Let `A` be a complete symmetric zero-diagonal signing of order `n`,
put `F=Phi(A)`, and let `x` and `y` be stable in opposite
orientations. After multiplying `A` by the orientation of `x`, assume

    x is +stable,
    y is -stable,

with

    Q_A(x)=F-D_x,
    Q_A(y)=-F+D_y.                                         (1)

Switch by `x`:

    B_ij = A_ij x_i x_j.                                   (2)

Then the all-ones state is +stable for `B`. Put

    z_i=x_i y_i,

and let

    T={i:z_i=-1},    S={i:z_i=+1}.                          (3)

The state `z` is -stable for `B`.

For a vertex `i`, write `a_i` for its signed degree inside its own
side of the cut and `c_i` for its signed degree across the cut:

for `i in T`,

    a_i=sum_(j in T, j!=i) B_ij,
    c_i=sum_(j in S) B_ij,                                 (4)

and for `i in S` use the analogous definitions with `S,T`
interchanged.

## Theorem: vertexwise bipartite dominance

For every vertex,

    boxed: c_i >= |a_i| >= 0.                              (5)

In particular every signed cross-degree is nonnegative.

### Proof

First take `i in T`. Since the all-ones state is +stable,

    a_i+c_i >= 0.                                          (6)

For the state `z` one has

    (Bz)_i = -a_i+c_i,

and `z_i=-1`. Since `z` is -stable,

    z_i(Bz)_i <= 0,

hence

    a_i-c_i <= 0,
    c_i >= a_i.                                            (7)

Together (6)--(7) give

    c_i >= max(a_i,-a_i)=|a_i|.

For `i in S`,

    (Bz)_i=a_i-c_i,
    z_i=+1.

The same two stability inequalities give

    a_i+c_i>=0,
    a_i-c_i<=0,

and again `c_i>=|a_i|`. QED.

## Exact aggregate identities

Let

    C=sum_(i in T,j in S) B_ij                              (8)

be the signed cross energy. The two parent energies are

    Q_B(1)=q_T+q_S+C=F-D_x,
    Q_B(z)=q_T+q_S-C=-F+D_y.                               (9)

Therefore

    boxed: C=F-(D_x+D_y)/2,                                (10)

and

    q_T+q_S=(D_y-D_x)/2.                                   (11)

Summing (5) on either side gives the new field bounds

    boxed:
    ||B[T] 1_T||_1 <= C,
    ||B[S] 1_S||_1 <= C.                                   (12)

Consequently

    |q_T|<=C/2,
    |q_S|<=C/2.                                             (13)

Equation (12) is vertexwise information compressed into a principal-block
statement; it is strictly stronger than the previously recorded support-size
bound

    |T||S| >= C.                                           (14)

Indeed, (14) follows immediately because each `c_i` is at most the size
of the opposite side, while (5) additionally controls every internal signed
degree separately.

## Signed-cut interpretation

Moving one vertex across the cut changes the signed cut weight by twice
`a_i-c_i`. Equation (5) therefore implies that `T|S` is a one-vertex
local maximum of the signed cut weight. More strongly, the same cut also
dominates the absolute internal signed degree at every vertex, a consequence
of simultaneously having the all-ones state +stable and the cut state
-stable.

Thus the opposite-phase child in the first Banaszczyk projection is not an
arbitrary principal support, even though it does not inherit one-sided
principal stability: after switching, its support is one side of a
vertexwise-dominant signed cut.

## Direct connection to the one-vertex recursion

The first coordinate Banaszczyk transform creates an opposite-phase child

    |x_T.u_T| < (D_x+D_y)/2 + r - 1.                        (15)

After switching by `x`, the child normal is simply `1_T`. The latest
recursion note previously retained only the aggregate fact that both sides
of `T|S` are `Omega(sqrt(n))` near the edge. Equations (5) and (12)
supply additional inherited structure:

1. every cross row and column sum is nonnegative;
2. each such cross sum dominates the magnitude of the corresponding
   principal signed degree;
3. the total absolute local-field mass of the all-ones child state on
   either principal side is at most the parent cross energy
   `F-(D_x+D_y)/2`.

This closes no convergence claim by itself. It does, however, remove the
description of opposite-phase children as structurally uncontrolled. Any
next recursion step may now use the vertexwise inequalities (5), not merely
the support lower bound (14).
