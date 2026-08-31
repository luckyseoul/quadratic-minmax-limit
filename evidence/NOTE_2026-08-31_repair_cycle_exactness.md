# Repair-cycle exactness and phase-bridge audit

**Date:** 2026-08-31
**Result status:** proved method barrier; no endpoint exclusion
**Scope:** the full Proposition 15.730 repair graph in both endpoint residues

**Later status:** Proposition 15.734 closes every `k=4p` residual boundary
for `p>=13` by an isolated-chart coefficient argument. The cycle no-go and
local-jet formulas below remain correct, but they are no longer a live
endpoint front.

## Question audited

Proposition 15.731 proposed composing its quadratic or cubic repair-swap
quotients around the `K_3` and `J(4,2)` factors of the repair graph.  The
target was a direction-sensitive cycle invariant which could be compared to
the residual phase and lift constraints.

The additive cycle route is now completely determined.  It has an explicit
vertex potential, so it cannot produce nonzero holonomy without importing
new data.

## Exact potential

For every repair `A`, define

`P_A=product_(u in A)L_u`, `Theta_A=P_A^2 Phi_A`.             (1)

On an oriented adjacent swap `A=C union {a}` to
`A'=C union {z}`, Proposition 15.731 gives

`L_z^2 Phi_A'-L_a^2 Phi_A=P_C Q_(a,z)`.                      (2)

Since `P_A=P_C L_a` and `P_A'=P_C L_z`, multiplying (2) by
`P_C^2` proves

`Theta_A'-Theta_A=P_C^3 Q_(a,z)`.                            (3)

Therefore, on every closed repair walk,

`sum P_C^3 Q_(a,z)=0`.                                      (4)

The cancellation occurs as a polynomial identity before applying any
linear functional.  Evaluation at a direction, extraction of a coefficient,
differentiation, and polarization all preserve zero.

The degree check is exact.  Write `p=3R+c`, `c in {1,2}`.  Then

- `k=|A|=2R+c+1`;
- `deg Phi_A=2R+2`;
- `deg Q=4-c`;
- `deg Theta_A=(2R+2)+2k=3(k-1)+(4-c)=deg(P_C^3Q)`.

In residue `c=1`, the envelope change

`Phi_A -> Phi_A+mu_A P_A`

changes `Theta_A` by `mu_A P_A^3` and changes an edge quotient by

`Q_(a,z) -> Q_(a,z)+mu_A' L_z^3-mu_A L_a^3`.                 (5)

This is again a vertex coboundary.  The cycle identity is gauge invariant
when the representatives are correlated by vertex; independently quotienting
every edge by its two cube directions only loses that correlation.

## Small cycles

Fix three points `i,j,k` on one rich block and the same repair choice `B`
off that block.  The triangle on repairs `A_ij,A_ik,A_jk` gives

`L_i^3 Q_(j,k)^i + L_k^3 Q_(i,j)^k + L_j^3 Q_(k,i)^j=0`.     (6)

Equation (6) applies to a trisecant `K_3` and to every triangular face of a
4-secant `J(4,2)`.  Chordless squares have the analogous four-term identity.
All are instances of (4).  They are useful consistency relations, but they
are not independent curvature.

An exact `p=5` probe was used only as an adversarial check, not as evidence
for the all-prime theorem.  It found actual nonzero quadratic quotient
triples whose cleared circulation is nevertheless identically zero, and it
showed that arbitrary fixed-direction evaluations can vary between zero and
nonzero.  The proposition and its deterministic certificate do not rely on
that finite computation.

## The surviving edge-local jet

Let `q` be the dual point of the rich line containing the exchanged points
`a,z` and the retained point `b`.  This rich line is a secant of both repairs,
so both envelope values at `q` are nonzero.  Taking the lowest homogeneous
part of (2) proves

`Q_(a,z)(q)=0`,                                             (7)

`P_(C-{b})(q) L_b j_q^1 Q_(a,z)`
`  =Phi_A'(q)L_z^2-Phi_A(q)L_a^2`.                          (8)

The first jet in (8) is nonzero.  If it vanished, two nonproportional local
line forms would have nonzero scalar squares equal.

With `l_r` the local class of `L_r`,
`Delta_rs=det(l_r,l_s)`, and

`K=Phi_ij(q) Delta_ij^2`,

the normalized formula is

`j_q^1 Q_i(j,k)=(K/P_B(q))`
`  (l_k^2/Delta_ik^2-l_j^2/Delta_ij^2)/l_i`.                 (9)

Here `K` is independent of the selected pair and is a nonzero square.  The
first-jet square character is consequently `chi_p(P_B(q))`.  Cubic gauge
changes vanish to order three, so the first jet is gauge invariant in
residue `c=1` as well.

This is genuine local structure, but it is not yet residual phase data.
`P_B(q)` is a repair-coloured product.  The residual parity set records the
symmetric difference of the `A`- and `T`-tangent fibres, so it fixes only the
product of the two colour characters, not the `A` character in (9).

## Two failed bridges

For a nonrich `b=2` direction, the `D`-fibre profile is

`((p-3)/2 empty, 2 singleton, (p-1)/2 double)`.              (10)

Deleting `R` points destroys at most `R` of the double fibres.  Hence a
repair has at least `(p-1)/2-R` direction secants and at most

`|A|-2((p-1)/2-R)=R+2` tangents.                             (11)

These are at most `R+2` distinct known roots on the dual direction-pencil
line. The envelope has degree `2R+2`, so (11) cannot force that line as a
component by distinct-root counting. At `p=31`
this is at most 12 roots against degree 22, for every one of the at least
`4+y` nonrich Paley-hard directions supplied by Proposition 15.728.

Products over repairs also lose the required colour.  Modulo squares, the
three pair selections on a trisecant are the masks `110,101,011`.  Their
span is the even-weight subspace of `F_2^3`, which does not contain the full
mask `111`. The full `1111` mask on a 4-secant is recoverable. Thus a
repair-product character cannot reconstruct a factor containing an
unselected trisecant block. For a jet on one rich block this applies whenever
another trisecant is present; a sole-trisecant target is the exact exception
to the mask argument and still lacks a proved colour-separating phase rule.

## Exact status and next front

Proposition 15.732 proves:

- every naturally cleared additive repair transition is an exact
  coboundary;
- the cubic-pencil ambiguity does not change that conclusion;
- the rich-direction quotient has a nonzero gauge-invariant first jet;
- the near-pairing directions are too tangent-sparse for a component
  root-count argument;
- the simplest product-over-repairs phase recovery loses trisecant parity.

It does **not** itself exclude an endpoint completion or change any top-level
gate. At the 15.732 stage, the next useful lemma had to add information absent
from the bare transition algebra. The two grounded possibilities were:

1. express the repair-coloured first-jet character in terms of coefficients
   of the signed residual lift; or
2. work directly with the simultaneous exact `p=31` identities
   `A_d=(1-x_s-x_t)^2` in at least fourteen same-Paley directions and
   exclude the common completion from that coefficient system.

Proposition 15.733 executes the second route at `p=31`; Proposition 15.734
then makes both routes unnecessary throughout the `k=4p,p>=13` endpoint.

## Artifacts

- `src/e1_gmin_m4_prop15732.py`
- `tests/test_prop15732.py`
- `evidence/e1_gmin_m4_prop15732.json`
