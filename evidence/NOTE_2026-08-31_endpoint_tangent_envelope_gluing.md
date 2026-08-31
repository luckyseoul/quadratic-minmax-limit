# Endpoint tangent envelopes glue below the quoted interpolation threshold

**Date:** 2026-08-31
**Status:** proved structural reduction, not endpoint exclusion.

## Setup

Fix any maximum repair `A` from Proposition 15.730.  Write

`n=|A|=k=p+1-R`, `t=p+2-n=R+1`, `d=2t=2R+2`.             (1)

For `a in A`, let `f_a(X)` be the product of the `t` tangent-line forms at
`a`.  Scale these products as in Segre's lemma of tangents, so that

`f_a(b)=(-1)^(t+1) f_b(a)` for distinct `a,b in A`.       (2)

Ball--Lavrauw [44, Lemma 10] proves (2).  Their Theorem 11 gives one explicit
degree-`2t` interpolation when `|A|>=2t+2`; it does not say that no tangent
envelope exists below that sufficient hypothesis.

## Compatible sections on the dual lines

For a point vector `a`, write

`a*={Z:a dot Z=0}` and `L_a(Z)=a dot Z`.

The linear map `X -> X cross a` has kernel spanned by `a` and maps onto
`a*`.  Every tangent-line form at `a` vanishes on `a`, so

`h_a(X cross a)=f_a(X)^2`                                  (3)

is a well-defined homogeneous degree-`d` section on `a*`.

The dual lines of an arc have no triple point.  At the node `a* intersect
b*`, the two values in (3) are `f_a(b)^2` and `f_b(a)^2`, up to the harmless
sign from reversing a cross product.  Equation (2) and the even degree make
them equal.  Thus the sections `{h_a:a in A}` are compatible at every node.

## Elementary line-gluing lemma

Let `L_1,...,L_n` be projective lines with no triple intersection, and give
each `L_j` a homogeneous degree-`d` section, compatible at all pairwise
nodes.  These sections extend to a homogeneous degree-`d` polynomial on the
whole plane.

To prove this, suppose a polynomial already has the right restriction on
`L_1,...,L_(j-1)`.  Its error on `L_j` vanishes at the `j-1` distinct old
nodes.  If `j-1<=d`, divide that binary form by the restrictions of the old
line equations, extend the quotient off `L_j`, and add it times the product
of the old line equations.  If `j-1>d`, the degree-`d` error has more roots
than its degree and is zero.  This induction glues every line.

The difference of two extensions vanishes on every `L_j`, hence is divisible
by their product.  Therefore the kernel has dimension

`binom(d-n+2,2)` if `n<=d`, and `0` if `n>d`.               (4)

Applying the lemma to (3) gives a homogeneous tangent-envelope polynomial
`Phi_A(Z)` of degree `d`, whose restriction to every `a*` is `h_a`.

## Exact endpoint uniqueness

The two endpoint residues land on opposite sides of the kernel boundary:

| residue | `n=|A|` | `d` | fixed-normalization envelope space |
|---|---:|---:|---|
| `p=3R+1` | `2R+2=d` | `2R+2` | `Phi_0 + lambda product_(a in A)L_a` |
| `p=3R+2` | `2R+3=d+1` | `2R+2` | unique |

Changing the initial common normalization of all tangent products rescales
the prescribed sections and the envelope.  Thus, without fixing that
normalization, the first row is the projective pencil
`P(span{Phi_0,product L_a})` with the pure product point omitted, and the
second row is a unique projective envelope.

This applies simultaneously to every one of the `3^x 6^y` repairs in the
common endpoint completion.  There are no missing existence conditions to
recover from the one- or two-point deficit in the sufficient hypothesis of
Ball--Lavrauw Theorem 11.

## One-block swaps have quadratic or cubic transition data

Let two adjacent repairs differ on one rich block:

`A=C union {a}`, `A'=C union {z}`.                           (5)

For `u in C`, put `ell_(u,a)(X)=det(u,a,X)`.  The tangent products

`ell_(u,a) f_u^A` and `ell_(u,z) f_u^(A')`                  (6)

are tangent products for the smaller arc `C`.  The determinant sign together
with (2) gives the Segre normalization for `C`.  Hence the two families in
(6) differ by one common scalar, not a point-dependent collection of
scalars.  Rescale the primed family by that scalar.

On every common dual line `u*`, equations (3) and (6) now give

`L_z^2 Phi_(A') = L_a^2 Phi_A`.                            (7)

Consequently

`L_z^2 Phi_(A') - L_a^2 Phi_A = P_C Q_(a,z)`,

`P_C=product_(u in C)L_u`.                                 (8)

The left side has degree `d+2`.  If `p=3R+2`, then `|C|=d`, so `Q_(a,z)` is
quadratic.  The two envelopes are unique after normalization, and this
quadratic is correspondingly fixed up to their common scalar.

If `p=3R+1`, then `|C|=d-1`, so `Q_(a,z)` is cubic.  Replacing the two
envelopes by other members of their pencils changes this cubic by an element
of

`span{L_a^3,L_z^3}`.                                       (9)

Thus the invariant transition datum in this residue is the cubic class
modulo those two pure cubes.

The relative edge scalings can be chosen coherently.  Every repair contains
the same singleton set `S`, with `|S|=c+1+2y>=2`.  Fix `e in S` and a form
`rho_e` for the unique line through `e` avoiding `D`.  For every repair set

`f_e^A=rho_e product_(v in D\A) det(X,e,v)`.               (10)

The `R` complement factors and `rho_e` are exactly the `R+1` tangents at
`e`.  If a swap replaces `a` by `z`, (10) changes by the exact factor
`det(X,e,a)/det(X,e,z)`.  Normalizing every other tangent product relative
to this fixed base makes the scalar in (7) equal to one on every repair-graph
edge.  The factors telescope around a closed walk.  This removes a potential
scaling cocycle, but it does not remove the line-product lift kernel in the
`p=3R+1` residue.

## Correct next attack

The repair graph is a Cartesian product of one `K_3` for every trisecant
block and one `J(4,2)` for every 4-secant block.  Proposition 15.731 puts a
unique envelope and quadratic transition on every vertex/edge in the
`p=3R+2` case, and an envelope pencil with a cubic transition class on every
vertex/edge in the `p=3R+1` case.

The exact open implication is exclusion of the 15.727/15.730 common
completion under the direction, phase, and lift constraints. Deriving a
nontrivial cycle-compatibility identity from these degree-two or degree-three
edge coordinates is the preferred next attack, but it is not itself a proved
necessary condition yet. This is a more precise route than asking for a
nonexistent unique-trisecant size ceiling or for missing tangent-envelope
existence conditions.

No such cycle contradiction is proved here.  Endpoint equality, larger
slack, residual (ii), Type I, and the quadratic-minmax limit remain open.

## Artifacts

- `src/e1_gmin_m4_prop15731.py`
- `tests/test_prop15731.py`
- `evidence/e1_gmin_m4_prop15731.json`
