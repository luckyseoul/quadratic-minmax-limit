# Field-width potential controls all non-Boolean Banaszczyk children

2026-09-25. Direct continuation of the frozen one-vertex extension route.
This note resolves the coefficient-growth concern left explicit in
`NOTE_2026-09-25_BANASZCZYK_STABLE_RECURSION.md`. It does NOT prove
convergence: the remaining issue is a lower/positivity control on the child
half-widths, not an upper control on transformed normals.

Duplication audit: repository searches for "field width invariant",
"Banaszczyk field norm", "transformed slab potential", "non-Boolean normals
invariant", and "compatibility child field" found no prior occurrence of the
potential or recurrence below. The general compatibility formula itself is
prior work from the Banaszczyk recursion note.

Let A be a complete symmetric zero-diagonal signing of order n and put

    F=Phi(A).

Pad every transformed normal by zeros on already eliminated coordinates, so
all normals continue to live in R^n and A acts on them.

For a symmetric slab

    |c.z| < b,                                               (1)

define its field-width potential

    P_A(c,b) := b + (1/2)||A c||_1.                          (2)

The potential is positively homogeneous.

## 1. Original stable slabs exactly saturate one common potential

Let x be stable in orientation sigma and let

    D=F-sigma Q_A(x).

Stability gives the exact identity

    ||A x||_1 = 2 sigma Q_A(x) = 2(F-D).                    (3)

The original one-vertex slab at desired increment r is

    |x.z| < D+r.                                             (4)

Therefore every original stable slab satisfies exactly

    boxed:
    P_A(x,D+r)=F+r.                                          (5)

Thus deficit disappears from this combined quantity.

## 2. General compatibility recurrence

Suppose two current parent slabs are

    |c.z|<b_c,    |d.z|<b_d,                                (6)

and the next eliminated coordinate is i. Negate either normal if needed so

    c_i>0,    d_i>0.                                        (7)

The exact general-normal compatibility formula from the previous recursion
note gives the child normal and threshold

    h = d_i c-c_i d,
    beta = d_i b_c+c_i b_d-2c_i d_i.                        (8)

Scale the child slab by 1/2:

    g=h/2,    b_g=beta/2.                                   (9)

Then triangle inequality gives

    ||A h||_1
      <= d_i ||A c||_1+c_i ||A d||_1.                      (10)

Consequently

    boxed:
    P_A(g,b_g)
      <= [d_i P_A(c,b_c)+c_i P_A(d,b_d)]/2
         -c_i d_i.                                         (11)

This is the basic all-depth recurrence.

## 3. Coefficients never blow up

Assume the parents have been normalized so that

    ||c||_infty<=1,    ||d||_infty<=1.                      (12)

Then 0<c_i,d_i<=1 and, coordinatewise,

    |h_j|
      <= d_i |c_j|+c_i |d_j|
      <= c_i+d_i
      <=2.

Hence the normalized child in (9) obeys

    boxed: ||g||_infty<=1.                                  (13)

Parents with zero coefficient at the eliminated coordinate do not enter a
compatibility pair; their slab simply persists unchanged. Therefore every
generated slab can be kept, at every depth, with infinity-normal norm at
most one.

This removes the coefficient-growth obstruction stated in Section 6 of the
previous recursion note.

## 4. Every genuine generation drops the potential by one

Suppose both parents satisfy

    P_A(c,b_c)<=L,
    P_A(d,b_d)<=L,                                          (14)

with L>=2 and the normalization (12). Equation (11) gives

    P_A(g,b_g)
      <= L(c_i+d_i)/2-c_i d_i.                              (15)

For 0<c_i,d_i<=1,

    L-L(c_i+d_i)/2
      = L[(1-c_i)+(1-d_i)]/2
      >= (1-c_i)+(1-d_i)
      >= 1-c_i d_i,                                        (16)

where L>=2 is used in the first inequality. Rearranging,

    boxed:
    L(c_i+d_i)/2-c_i d_i <= L-1.                            (17)

Therefore

    boxed:
    P_A(g,b_g)<=L-1.                                        (18)

The drop is independent of how small the nonzero pivot coefficients are.

## 5. Generation-depth invariant

Give each original stable slab generation 0. A persisted zero-pivot slab keeps
its generation. A compatibility child has generation

    gen(child)=1+min(gen(parent 1),gen(parent 2))             (19)

when using a common upper envelope by generation.

Starting from (5), induction with (18) gives, as long as the displayed right
side is at least two,

    boxed:
    gen(c,b)>=k
      => P_A(c,b) <= F+r-k.                                 (20)

In particular:

- every first-generation child has P_A<=F+r-1;
- every second-generation non-Boolean child has P_A<=F+r-2;
- arbitrary later non-Boolean children remain infinity-normalized and acquire
  one additional unit of potential drop for each genuine compatibility
  generation.

Because F=Theta(n^(3/2)) while there are only n coordinate eliminations,
the harmless technical condition L>=2 holds throughout the asymptotic regime.



## 6. Complementary lower bound: normalized width loses at most one

There is also a lower recurrence, after using slab homogeneity to normalize
every nonzero parent normal by

    ||c||_infty=||d||_infty=1.                              (21)

Write the corresponding normalized half-widths simply as

    w_c=b_c,    w_d=b_d,

and put

    p=c_i,    q=d_i,     0<p,q<=1.                          (22)

For the compatibility child,

    h=q c-p d,
    beta=q b_c+p b_d-2pq.                                  (23)

As before set g=h/2 and b_g=beta/2. If g is nonzero, let

    m=||g||_infty.

Because the parent infinity norms are one,

    m <= (p+q)/2.                                           (24)

If the weighted parent width

    W=(q w_c+p w_d)/(p+q)                                  (25)

is greater than one, then beta>0 because

    2pq/(p+q) <=1.                                         (26)

Renormalize the child slab by dividing by m. Its normalized half-width is

    w_child=b_g/m.

Using (24) and positivity,

    w_child
      >= [q w_c+p w_d-2pq]/(p+q)
      = W-2pq/(p+q)
      >= W-1.                                               (27)

Therefore

    boxed:
    w_child
      >= (q w_c+p w_d)/(p+q)-1.                            (28)

This is the complementary lower estimate to the upper field-width potential.

Original stable slabs have normalized width

    w_x=D_x+r.                                              (29)

Thus along any compatibility tree, normalized width is bounded below by the
pivot-weighted ancestral average of the original quantities D_x+r, minus one
for each genuine compatibility step in the corresponding weighted recursion.

In the equal-depth case this takes the clean form

    boxed:
    generation k child:
    w >= r + D_bar - k,                                     (30)

where D_bar is a convex combination of the stable deficits at its leaves.

More generally define recursively

    R_leaf=D_leaf,
    K_leaf=0,

and at a child with pivot weights q,p,

    R_child=(q R_1+p R_2)/(p+q),
    K_child=1+(q K_1+p K_2)/(p+q).                          (31)

Then, whenever the right side remains positive through the recursion,

    boxed:
    w_child >= r+R_child-K_child.                           (32)

Hence deep dangerous children can only be assembled from ancestors whose
pivot-weighted average stable deficit is small enough that it fails to pay
for the accumulated compatibility cost.

This converts the remaining positivity problem into a near-edge ancestry
problem rather than an uncontrolled coefficient problem.

## 7. Exact second-step specialization

First-generation normals are in {0,+-1}^n. If two such slabs participate in
the next genuine elimination, orient them so their pivot entries are both
one. Their normalized child is

    g=(c-d)/2,                                               (33)

so

    g_j in {0,+-1/2,+-1}.                                   (34)

If the parent potentials are at most F+r-1, then

    boxed:
    b_g+(1/2)||A g||_1 <= F+r-2.                            (35)

If the two parent normalized half-widths are w_c,w_d, then also

    boxed:
    w_child >= (w_c+w_d)/2-1.                               (36)

Thus the second transform has simultaneous upper field control and lower
width control, even though its normal is genuinely non-Boolean.

## 8. Relation to the opposite-phase dominance result

First-generation normals are in {0,+-1}^n. If two such slabs participate in
the next genuine elimination, orient them so their pivot entries are both
one. Their normalized child is

    g=(c-d)/2,                                               (21)

so

    g_j in {0,+-1/2,+-1}.                                   (22)

If the parent potentials are at most F+r-1, then

    boxed:
    b_g+(1/2)||A g||_1 <= F+r-2.                            (23)

Thus the second transform's non-Boolean coefficients are quantitatively
controlled without requiring the two parent constraints to share a switching
gauge.



For an opposite-phase first-generation child v=x 1_T with
d=(D_x+D_y)/2, the previous vertexwise-dominance note gives the sharper
direct field estimate

    ||A v||_1 <= 2(F-d).                                    (24)

Its half-width is

    b_v=d+r-1,                                              (25)

so (24)--(25) recover

    P_A(v,b_v)<=F+r-1.                                      (26)

The new recurrence shows that this is not an isolated opposite-phase
phenomenon: the same potential is the correct common quantity for every
compatibility child and it continues through non-Boolean generations.

## 9. What remains open

Equations (20) and (32) now give both an upper field-width control and a
lower normalized-width recursion. They guarantee positivity through every
branch for which the pivot-weighted ancestral deficit reserve pays the
accumulated compatibility cost.

The remaining obstruction is narrower: one must rule out a deep
compatibility branch with

    r+R_child-K_child <= 0,

that is, a branch built from an unusually low pivot-weighted average of
near-edge stable deficits. The coefficient-growth problem and the generic
non-Boolean-normal problem are removed; what remains is a near-edge ancestry
problem in the stable skeleton.
