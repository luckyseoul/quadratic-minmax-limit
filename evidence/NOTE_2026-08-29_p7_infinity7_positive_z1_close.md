# Proposition 15.715: positive p7 infinity-plus-seven z=1 is closed

Let `z` count undetermined affine directions of the seven-point finite part
of the boundary.  In the positive-product branch, a `z=1` boundary has one
direction with odd-fibre count `b=7`.  The two type budgets leave exactly
four mean allocations: either that direction rises from mean zero to eight,
or it remains at zero and one of the three ordinary directions of the same
quadratic type rises from mean eight to sixteen.

The variable Johnson-slice catalogs have 1,764 rows, except for `b=3` at
mean sixteen, which has 2,233.  The common normalized 282-by-1,225 edge
system has rank 147 modulo seven and 135 audited left-null dependencies.
The normalization loses no boundary: the infinity vertex has positive odd
degree, so one selected infinity edge exists, and an affine translation can
move its finite endpoint to zero while fixing infinity and preserving all
Paley signs and direction types.

For each elevated direction,
`scripts/p7_infinity7_positive_z1_mod7_gpu.py` selects 23 dependencies of
full rank on that direction block and packs their base-seven values
losslessly.  This is only a prefilter.  Every projected boundary is
reconstructed on the host, all four allocations are considered, and each
complete catalog is compared on all 135 dependencies.

Two complete V100 runs with 65,535 and 32,768 CUDA blocks independently scan
all `C(49,7)=85,900,584` seven-point finite boundaries.  Both find

```text
actual z=1 boundaries:             6,324,528
boundaries passing the projection:     1,326
survivors of all 135 dependencies:          0.
```

The sorted 1,326-rank set agrees by SHA-256 between launch geometries.  Thus
every positive `z=1` boundary is excluded before any edge solve.  Together
with Proposition 15.714, the actual positive scope falls from 6,453,552 to
129,024 boundaries, and the projected outer envelope falls from 792 to 492
profiles with `z=2,3,7`.  Those branches, the entire negative branch,
residual (ii), and every top-level gate remain open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found only
the surrounding finite-direction, Paley-49, and modular-incidence theory,
not this exact census or exclusion.  This is attribution context, not
mathematical evidence or a formal priority claim.
