# Complete p=7 size-eight conic-subbranch exclusion

Proposition 15.662 closes both product signs of the finite `p=7`
size-eight branch having the minimum possible eight odd secants. By Segre's
odd-order oval theorem, this is exactly the affine-conic subbranch. It does
not close all size-eight boundaries, residual (ii), Type I, R1, global QVAR,
or the limit.

## Exact scope

The two complete CUDA censuses each check
`C(49,8)=450,978,066` finite boundaries. For each product sign they find
108,754,569 exact directional-floor survivors in 5,152 ordered profiles.
Exactly 6,174 boundaries attain the minimum of eight odd secants. The
projective incidence count independently gives

`16,758 * 21 / 57 = 6,174`

nondegenerate conics disjoint from the fixed line at infinity. Segre's
theorem identifies every eight-point arc over the odd field `F_7` as a
conic. Of those 6,174 conics, 4,851 fail the floor and 1,323 survive. The
48-element square-semilinear stabilizer partitions the survivors into 32
orbits.

The full floor census is much larger than this conic branch. Subtracting the
1,323 conic survivors leaves exactly 108,753,246 nonconic floor survivors
per sign. That number is the honest open size-eight remainder.

## Saturated 25-orbit component

For `c_H=-1`, 25 of the 32 orbits are saturated. They cover 1,176
boundaries and have 24 exact type-compatible directional-mean allocations
each, hence 600 allocations total. Independent reconstruction gives the
following disjoint exclusion ledger:

- 355 allocations are infeasible in the first exact CP-SAT pass.
- Six more are infeasible in longer exact CP-SAT runs.
- The remaining 239 are excluded by simultaneous mod-three/mod-seven
  complete-catalog joins.

The largest reused catalog contains all 575,407 zero-parity mean-16 slacks.
Its `575407 x 35` array has SHA-256
`78874943f492ea6a041d516aaf90a04c0ccc30701416d877e6458d8ffdeaa00e`.
Every row has the required mean, parity, and 14 primitive integer
degree-two kernel equations. The independent replay leaves no deferred or
modularly consistent allocation.

## Exceptional seven-orbit component

The other seven orbits cover 147 boundaries and have 180 exact allocations
per representative, hence 1,260 total. The first exact pass excludes 172.
An ordinary three-block projected V100 join independently excludes 662 of
the 1,088 initially unknown leaves, with no supported survivor.

The final 426 leaves each have exactly one zero-parity direction with scaled
mean 20, 24, or 32. Directly materializing those giant catalogs is
unnecessary. The common mod-seven edge system has a `135 x 282` dependency
matrix. For each possible high direction `h in {0,2,5,7}`, the left kernel
of its 35-column block has dimension 112. A deterministic 66-row dependency
family is selected as three disjoint 22-row projections. Its rank is zero
on direction `h` and 14 on every other direction, so the entire high catalog
contributes the zero signature and is eliminated algebraically. The seven
remaining catalogs have at most 1,764 rows. Exact V100 joins then exclude
all 426 leaves, split `59,66,59,61,60,60,61` across orbit indices
`0,5,8,25,26,30,31`.

`scripts/p7_exceptional_omit_high_audit.py` independently reconstructs the
dependency combinations, all 84 projected catalogs, and every high-leaf
certificate. Together with the ordinary audit, it verifies the exact
partition `172+662+426=1260` and zero remaining exceptional allocations.
An independent planted CPU/GPU self-test also returns one match on each
backend.

## Opposite-sign transfer and aggregate audit

Multiplication by a nonsquare in `F_49`, together with switching only the
infinity coordinate, fixes infinity, finite zero, and the distinguished
edge and satisfies the exact signed conference anti-isometry
`S C_perm S = -C`. A finite boundary gives infinity even degree; because
the graph has 29 edges, the transformation flips its Paley product sign
while preserving the normalized score and eigenshell conditions. The
aggregate audit explicitly maps all 6,174 minimum boundaries and all 1,323
floor-surviving conics for `c_H=-1` bijectively onto the corresponding
`c_H=+1` sets. Thus one complete exclusion transfers to the other sign.

`scripts/p7_size8_conic_global_audit.py` independently checks both complete
floor files and their hashes, reconstructs all 32 stabilizer orbits,
reenumerates all 600 saturated allocations, validates the disjoint 355/6/239
certificate partition, checks the ordinary and high exceptional audits,
and verifies the exact sign-transfer bijection. Its output is
`evidence/p7_size8_cminus1/global_conic_audit.json`, SHA-256
`85f927f41b3ffc9afe1a101584e95ed852709ca6e861b439d8da1715008640a9`.
It records both `all_32_conic_orbits_both_signs_excluded=true` and
`closes_all_p7_size8=false`.

Raw records are archived at
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-conic/`.

## Literature and OEIS context check

[Müller's modern proof of Segre's theorem](https://arxiv.org/abs/1311.3082)
confirms the classification input that every oval over a finite field of odd
order is a conic. [Ball--Csajbók](https://arxiv.org/abs/1711.10876) studies
few odd secants for `q+2` points, and
[Di Benedetto--Solymosi--White](https://arxiv.org/abs/2001.06994) studies
directions in affine Galois planes and a Paley-graph application; neither
contains this residual allocation classification. The previously flagged
[arXiv:2305.03523](https://arxiv.org/abs/2305.03523) concerns Bellman
functions on planar nonconvex domains and does not overlap this finite
conic audit.

Exact OEIS API searches on 2026-08-27 return no sequence entry containing
108,754,569, 108,753,246, or 575,407. The input count 450,978,066 appears,
as expected, in binomial-coefficient sequences including A017765 because it
is `C(49,8)`. These are duplication/context checks only; no new-sequence
claim is made.
