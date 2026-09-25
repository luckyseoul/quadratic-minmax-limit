# Coordinate-variation one-vertex extension criterion

2026-09-25. Directly in scope under
`SCOPE_FREEZE_2026-09-25_ONE_VERTEX_ONLY.md`.

This note extracts a finite-direction strengthening of the Guo--Fang--Lu
directional-total-variation machinery and applies it directly to the exact
one-vertex extension body. It does **not** prove convergence, but it replaces
the reciprocal-square Komlos mass condition by a correlation-sensitive
convex-body criterion.

Reference: S. Guo, E. X. Fang, J. Lu, *Vector Balancing via Directional
Total Variation*, arXiv:2609.11189v1 (2026), especially Lemmas 3.2--3.4 and
4.1--4.3.

## 1. Exact extension body

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), and let T_A be the antipodal stable skeleton from
NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md. For [x] in T_A write

    D([x]) = F-E_A([x]).

Fix r>0 and M>1. Define the bounded open symmetric convex body

    K(A,r,M)
      = { z in R^n :
            |z_i|<M for every i,
            |z.x|<D([x])+r
              for every [x] with D([x])+r<=n }.             (1)

If a sign vector a lies in K(A,r,M), then every active stable class obeys

    |a.x|-D([x]) < r.

Every inactive class obeys the same inequality automatically because
|a.x|<=n<D([x])+r. Hence the exact stable-state reduction gives

    Phi(A extended by a) < F+r.                             (2)

Thus it is enough to prove that K(A,r,M) contains a vertex of the Boolean
cube.

## 2. Finite-coordinate version of Guo--Fang--Lu stability

For an integrable density rho and direction u, write V_u(rho) for its
directional total variation.

The Guo--Fang--Lu proof of finite-step stability actually yields the
following finite-direction specialization.

### Lemma

Let K be a nonempty bounded open convex set in R^n, and suppose K supports a
probability density rho in BV such that

    V_{e_i}(rho) <= 1/3,       i=1,...,n.                   (3)

Then for every coordinate direction e_j, the open Banaszczyk transform
T_{e_j}K is nonempty and supports another probability density rho' satisfying

    V_{e_i}(rho') <= 1/3,      i=1,...,n.                   (4)

### Proof

Use the geometric lift and fiber rearrangement of Guo--Fang--Lu,
Lemmas 4.1--4.3, with v=e_j.

Their exact rearranged-height identity needs only the single step bound

    V_{e_j}(rho) <= 1/3

to guarantee mean absolute height at least one.

Fiber rearrangement separately preserves every horizontal translation
variation, so after rearrangement

    V_{(e_i,0)}(R*) <= V_{e_i}(rho) <= 1/3                 (5)

for every i.

It remains only to extract one density on the height-one section with all
the coordinate bounds simultaneously. The proof of their prescribed-section
Lemma 3.4 already works for this finite list. For arbitrary simplex weights
alpha_i define

    H_alpha(xi)=sum_i alpha_i |e_i.xi|.

Equation (3.9) and the same Jensen/Minkowski-convexity argument give

    h_{H_alpha}(D_1) <= sum_i alpha_i V_{(e_i,0)}(R*)
                     <= 1/3.                               (6)

The finite-direction separation argument in the first part of their
Lemma 3.3 then supplies one probability density on D_1 for which every
V_{e_i} is at most 1/3. By Lemma 4.1, D_1=T_{e_j}K. QED.

The important point is that no variation estimate in an arbitrary unit
direction is needed. Only the n actual signing directions e_i are retained.

## 3. Coordinate-variation signing theorem

Iterate the lemma through e_1,...,e_n. Every transformed body remains open,
bounded, convex, nonempty and symmetric. At the final stage symmetry and
convexity put zero in the final body. Backward sign recovery through the
translate-containment lemma therefore gives signs a_i in {+-1} with

    a=(a_1,...,a_n) in K.                                  (7)

Consequently:

    THEOREM.
    If K(A,r,M) supports a probability density rho with

        max_i V_{e_i}(rho) <= 1/3,                          (8)

    then there is an actual incident sign row a such that

        Phi(A extended by a) < Phi(A)+r.                    (9)

For an exact minimizer and r=r_n this gives alpha_(n+1)<alpha_n.
For r=r_n+rho_n with sum rho_n/(n+1)^(3/2)<infinity it gives convergence.

This is a direct one-vertex certificate. Unlike the existing Xi_A(r)
corollary, it does not replace the stable constraints by independent
coordinate directions in constraint space; correlations among the stable
normals are retained in the geometry of K(A,r,M).

## 4. Equivalent finite Cheeger certificate

For alpha in the probability simplex let

    H_alpha(xi)=sum_i alpha_i |xi_i|.

The finite-direction separation argument also gives the equivalent
sufficient condition

    h_{H_alpha}(K(A,r,M)) <= 1/3
    for every alpha in Delta_(n-1).                         (10)

Indeed, (10) supplies one common density satisfying (8).

A more concrete, stronger-than-necessary test comes from the uniform density
on K. For a convex body K,

    V_{e_i}(1_K/vol K)
      = 2 vol_(n-1)(proj_(e_i^perp) K) / vol_n(K).          (11)

Therefore it is sufficient that

    vol_n(K)
      >= 6 vol_(n-1)(proj_(e_i^perp) K)
      for every i.                                         (12)

Equivalently, the mean length of a coordinate fiber of K, averaged over its
coordinate projection, is at least six in every coordinate direction.

Equation (12) is not asserted here for minimizers. It is a concrete new
geometric target whose verification would immediately produce the required
extension row.

## 5. Why the old reciprocal-square target is too lossy at low order

For the repository's exact order-15 minimizer
`evidence/K15_exact_minimizer_20260911.json` with F=27, exhaustive
projective enumeration gives 340 stable antipodal classes with deficit
histogram

    D=0:  66
    D=2:  66
    D=4:  39
    D=6: 117
    D=12: 13
    D=14: 39.                                               (13)

The neutral increment is

    r_15 = 27[(16/15)^(3/2)-1]
         = 2.7445120988729617... .                          (14)

The active reciprocal-square mass is

    Xi_A(r_15) = 14.141428337053116... ,                    (15)

whereas the Guo--Fang--Lu cube/Komlos corollary requires

    Xi_A(r_15) <= 1/(18 pi)
               = 0.01768388256576615... .                  (16)

Thus the scalar Xi certificate misses this exact low-order object by nearly
three orders of magnitude. This does not disprove an asymptotic Xi bound,
but it demonstrates that a successful extension argument should preserve
correlation among the dangerous states rather than pay for them as
independent reciprocal-square constraints.

The coordinate-variation body (1) does exactly that.

## 6. Next live implication

The project is now reduced to a concrete in-scope geometric question:

    prove (8), or the Cheeger form (10), for
    K(A_n,r_n+rho_n,M)

with a Dini-summable rho_n, uniformly for exact minimizers A_n and some
arbitrary fixed M>1 (or with M sent large after the estimate).

The strongest elementary subtarget is the coordinate-fiber estimate (12).
Failure of (12) does not kill the route because a nonuniform density can
still satisfy (8).
