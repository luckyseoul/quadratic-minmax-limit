# Proposition 15.716: positive p7 infinity-plus-seven z=2 is closed

Every seven-point boundary with at least two undetermined directions is a
permutation graph for each pair of those directions.  Enumerating the 28
direction pairs and all `7!` permutation graphs gives 141,120 incidences:

```text
123,480 + 3(5,488) + 21(56) = 28(7!) = 141,120.
```

The three terms are the exact actual-boundary counts at `z=2,3,7`, with
multiplicity `C(z,2)`.  The full 2,352-element square affine-semilinear group
reduces them to 92, 10, and two orbits.  The `z=2` orbit sizes are
`18*588 + 52*1176 + 22*2352 = 123,480`.

Because these are unpointed translation orbits, the distinguished-edge
normalization is not imposed.  Deleting precisely that row leaves the
translation-equivariant 281-by-1,225 edge system, rank 146 modulo seven,
with 135 audited left dependencies and direction blocks beginning at row
offset one.

The exact `z=2` mean ledger has 1,232 leaves, not merely the 1,184
residue-zero leaves:

```text
192  one-q=2 high-direction leaves
992  two-q=1 complete-catalog leaves
 48  residue-four four-catalog leaves
----
1232  exact leaves.
```

The last 48 occur when both undetermined directions have the same quadratic
type.  Their `b=7` directions have mean four, their two ordinary same-type
directions mean twelve, and the opposite type mean eight.  The variable
catalogs have only 56 rows, or 62 for `b=3`.  This family was caught by an
independent coverage audit before any proposition was claimed.

`scripts/p7_infinity7_positive_z2_mod7_join.py` rejects the 192 high leaves
using the full 112-dimensional dependency space annihilating the entire high
direction block.  This is a necessary block-image relaxation larger than the
unknown high catalog, so rejection remains rigorous.  It joins the 992
two-catalog leaves exactly on all 135
coordinates and handles the 48 residue-four leaves by exact 2+2
meet-in-the-middle joins.  No affine-span relaxation is called exact.

Nuka processed all 92 orbits and all 1,232 leaves in 10.5 seconds.  Every
leaf fails: zero survive modulo seven.  Thus all 123,480 actual positive
`z=2` boundaries are excluded.  The positive remainder is 5,544 actual
boundaries in twelve orbits: 5,488 at `z=3` and 56 at `z=7`.  The projected
outer envelope falls from 492 to 212 profiles.  The positive endpoint,
negative branch, residual (ii), and every top-level gate remain open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found no
prior occurrence of this orbit/catalog exclusion.  OEIS contains the number
123,480 in an unrelated injective-function triangle; that numerical match
has no evidentiary role.  This is attribution context, not a formal priority
claim.
