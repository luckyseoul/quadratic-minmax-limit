# One-coordinate Banaszczyk projection and stable-pair recursion

2026-09-25. Direct continuation of the frozen one-vertex extension route.
This note does not prove convergence. It identifies exactly what happens to
the stable-state slabs after one coordinate Banaszczyk transform and connects
that transformation to the existing stable-pair factorization theorem.

## 1. Setup

Let A be a complete symmetric zero-diagonal signing of order n, F=Phi(A),
and let the active oriented stable constraints be indexed by pairs (sigma,x).
Write

    D_sigma(x)=F-sigma Q_A(x),

and give the constraint the half-width

    b_x=D_sigma(x)+r,                                       (1)

where r>1. The corresponding open slab is

    |x.z|<b_x.                                               (2)

Fix a coordinate i. Since (2) is unchanged when x is replaced by -x,
orient every active normal so that x_i=+1. Write z=(s,u), where s=z_i and
u contains the remaining coordinates, and write

    x=(1,x_-).

Then the x-slab has fiber

    |s+x_-.u|<b_x.                                          (3)

## 2. Exact projection after one Banaszczyk transform

For an open convex body K, Guo--Fang--Lu use the open Banaszczyk transform

    T_(e_i) K
      = ((K-e_i) intersect (K+e_i)) + (-2,2)e_i.            (4)

The last Minkowski summand does not change the projection onto e_i^perp.
Hence the projection of T_(e_i)K equals the projection of

    L=(K-e_i) intersect (K+e_i).                            (5)

For a fixed u, membership (s,u) in L means that both (s+1,u) and
(s-1,u) lie in every parent slab. By (3), for each parent x this is exactly

    |s+x_-.u|<b_x-1.                                        (6)

Thus a projected point u is feasible iff the open intervals in (6) have a
common point. Intervals on the line satisfy Helly with Helly number two, so
this is equivalent to every pair intersecting. Therefore the projection is
described EXACTLY by

    |(x_- - y_-).u| < b_x+b_y-2                            (7)

for every pair of active oriented parent normals x,y.

Let

    T=T(x,y)={j != i : x_j != y_j}.

Because x_i=y_i=1,

    x_- - y_- = 2 x_T                                      (8)

after padding x_T by zero off T. Dividing (7) by two gives the exact child
constraint

    boxed:
    |x_T.u_T|
      < (D_sigma(x)+D_tau(y))/2 + r - 1.                    (9)

Here sigma,tau are the stable orientations of the two parent states.
Equation (9) is the first-step compatibility law for the actual one-vertex
extension body. It is not a heuristic or an entropy estimate.

## 3. Same-phase parents recurse to a stable principal problem

Assume sigma=tau and x,y are both sigma-stable. The existing
stable-pair factorization theorem gives, on their disagreement block T,

    x_T is sigma-stable for A[T],                           (10)

and

    delta_T
      :=P_sigma(A[T])-sigma Q_(A[T])(x_T)
      <= (D_sigma(x)+D_sigma(y))/2.                         (11)

Consequently the exact child constraint (9) is implied by

    boxed:
    |x_T.u_T| < delta_T+r-1.                                (12)

Thus every SAME-PHASE pair created by the first Banaszczyk coordinate step
is dominated by a stable-state slab of a smaller principal signing, with
exactly one unit of extension buffer consumed.

This is the missing direct connection between the September-24 stable-pair
factorization and the one-vertex transform route. The factorization theorem
by itself was duplicate geometry; after (9), its inherited-stability clause
has an immediate extension-row use.

For exact same-phase maximizers D_x=D_y=0, (11) sharpens to delta_T=0:
their child normal is an exact one-sided stable maximizer of A[T], and its
new half-width is exactly r-1.

## 4. Opposite-phase parents must have a large disagreement cut

Now let x be +stable and y be -stable, with

    Q_A(x)=F-D_x,
    Q_A(y)=-F+D_y.                                         (13)

Again let T be their disagreement set and S=T^c. Internal pair products on
T and S are identical for x and y, while their T--S cross products have
opposite signs. Writing C for the cross energy of x,

    Q_A(x)=q_T+q_S+C,
    Q_A(y)=q_T+q_S-C.                                      (14)

Subtracting gives

    boxed:
    C=F-(D_x+D_y)/2.                                       (15)

Since C is a signed sum of exactly |T||S| unit terms,

    boxed:
    |T|(n-|T|)
      >= F-(D_x+D_y)/2.                                    (16)

Therefore dangerous opposite-phase pairs cannot have a tiny disagreement
side. If d=(D_x+D_y)/2<F and s=min(|T|,n-|T|), then

    s(n-s)>=F-d,                                            (17)

and in particular

    boxed:
    s >= (F-d)/n.                                          (18)

At the edge F=Theta(n^(3/2)) and d=O(n), (18) forces

    s=Omega(sqrt(n)).                                      (19)

Thus the first transform creates two qualitatively different child types:

1. same phase: hereditary stable principal slabs, at cost one unit of buffer;
2. opposite phase: no inherited stability is asserted, but the child support
   and its complement are both Omega(sqrt(n)) in the near-edge regime.

## 5. Why this matters for a Dini extension proof

A constant excess is already sufficient for convergence: if one can prove

    Phi(A_n extended by a)
      <= m_n+r_n+C                                          (20)

for one absolute C and one exact minimizer A_n at every sufficiently large
n, then rho_n=C and

    sum C/(n+1)^(3/2)<infinity.                             (21)

Equation (12) shows why a constant budget is plausible for same-phase
recursion: each genuine coordinate-elimination compatibility step spends
one raw unit, while its support may shrink.

Equation (16) supplies the complementary shrinkage mechanism for
opposite-phase pairs: at small deficit neither side can have fewer than
order sqrt(n) vertices.

What is NOT yet proved is a full iteration invariant for the entire
transformed body. The transform remains n-dimensional; (9) describes its
projection onto e_i^perp, not the whole transformed body. In addition,
same-phase pairs may have disagreement support very close to n after the
forced orientation x_i=y_i=1, and opposite-phase children do not inherit a
stable orientation on A[T]. Those are the two concrete obstructions that a
complete recursive proof must still resolve.

## 6. Exact general-normal formula

For later iterations it is useful to record the non-Boolean version.
Suppose two parent slabs have normals c,d with c_i,d_i>0,

    |c.z|<b_c,   |d.z|<b_d.                                (22)

Eliminating coordinate i from the shrunken fibers gives horizontal normal

    h=d_i c-c_i d,     h_i=0,                              (23)

with compatibility threshold

    boxed:
    |h.z|
      < d_i b_c+c_i b_d-2c_i d_i.                          (24)

For Boolean parent normals with c_i=d_i=1, (24) reduces to (7).

Equation (24) is the correct object for more than one transform step: after
the first step the normals are no longer generally Boolean. Any successful
iteration must control the growth of these coefficients and the associated
thresholds, rather than silently reusing the Boolean formula.
