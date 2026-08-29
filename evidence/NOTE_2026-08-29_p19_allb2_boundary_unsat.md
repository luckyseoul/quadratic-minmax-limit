# p=19 all-b2 slack-twenty boundary exclusion

**Date:** 2026-08-29
**Proposition:** 15.698
**Status:** proved computationally; one profile excluded, endpoint open

The last slack-twenty profile is

```text
phase zero: {0:5,16:5}
phase one:  {2:10}.
```

Proposition 15.694 forces every witness to be an 11-arc core plus five
deleted points. Both sets are arcs, and each deleted point lies on exactly
one core secant. The exact native-XOR model imposes this repair structure,
the complete affine line parity transform `r=A*x`, its inverse `x=A^T*r`,
weight 16, and the full directional profile.

The normalization is lossless. Choose a phase-zero `b=0` direction. Every
line meets the boundary evenly, so a retained core point has a boundary
partner on its line. Translation and division by their square difference
send the pair to field elements zero and one while preserving the Paley
phase. Thus the model's fixed pair covers every orbit.

CryptoMiniSat 5.11.21 returned `UNSATISFIABLE` twice:

```text
nuka, 8 threads:              174.08 solver seconds
soulkiller registered ECC, 16 threads: 160.89 solver seconds
```

Each run used 1,184,892 clauses, 741 native XOR constraints, and 776 exact
cardinality constraints. Raw JSON and hashes are under
`evidence/p19_allb2_boundary_unsat/`.

This is infeasibility of the boundary itself, not failure of an attack
method. Therefore no affine edge lift exists. Proposition 15.697 leaves only
even infinity degrees `0,20,38`; nonsquare dilation flips both the direction
type and `c_H`, so both signs are excluded.

All p=19 slack-twenty profiles are now closed. Three profiles remain, with
slack histogram `{24:1,28:1,32:1}`. The p=19 endpoint, residual (ii), R1,
Type I, and the limit remain open.
