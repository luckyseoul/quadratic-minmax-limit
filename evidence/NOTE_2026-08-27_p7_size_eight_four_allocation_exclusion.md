# Complete p=7 size-eight four-allocation exclusion

Proposition 15.664 closes, for both product signs, the finite `p=7`
size-eight stratum having exactly four admissible directional-mean
allocations. It eliminates 23,563,806 floor-surviving boundaries per sign.
It does not close all size-eight boundaries, residual (ii), Type I, R1,
global QVAR, or the limit.

## Exact scope

After the conic stratum of Proposition 15.662 and the forced-floor stratum
of Proposition 15.663, the complete floor census leaves 24,983,238
boundaries per sign. Reconstructing every exact directional-mean allocation
from the 5,152 ordered floor profiles gives this complete partition:

| allocations | ordered profiles | boundaries per sign |
|---:|---:|---:|
| 4 | 2,245 | 23,563,806 |
| 11 | 248 | 154,056 |
| 16 | 516 | 1,194,816 |
| 24 | 8 | 1,176 |
| 44 | 110 | 69,384 |

The boundary column sums to 24,983,238. The four-allocation stratum has
94,255,224 allocation leaves. Its odd-secant distribution is

| odd secants | boundaries per sign |
|---:|---:|
| 16 | 691,488 |
| 20 | 5,603,640 |
| 24 | 9,190,146 |
| 28 | 5,990,544 |
| 32 | 1,846,908 |
| 36 | 232,848 |
| 40 | 5,880 |
| 44 | 2,352 |

For `c_H=-1`, 17,298,078 boundaries have type-floor sums `(24,32)` and
6,265,728 have `(32,24)`. The two counts swap for `c_H=+1`.

## Four-leaf rigidity and raised-direction omission

Each quadratic type has exact scaled-mean sum 32. In this stratum one type
has floor sum 24 and the other has floor sum 32. The saturated type stays at
its floor. The deficient type has exactly four allocations, each obtained
by raising one of its four directions by eight and leaving the other seven
directions at their floors.

The common affine score equations form a `282 x 1225` integer matrix. Over
`F_7` it has rank 147 and a 135-dimensional left kernel. Restricting those
dependencies to any one direction's 35 score columns has rank 23. Therefore
the subspace which is identically zero on that complete direction block has
dimension `135-23=112`. The raised catalog has zero contribution to every
conditioned dependency, regardless of whether it has 1,764 or 2,233 rows.

`scripts/p7_size8_one_elevation_tables.py` deterministically selects 22
independent conditioned dependencies for each of the eight possible raised
directions and materializes every remaining floor-catalog contribution.
Every non-raised floor catalog is a singleton except the phase-one,
four-odd-fibre, scaled-mean-14 catalog, which has 36 rows; at most one such
catalog occurs in a leaf. NUKA built and audited the cache in 6.07 seconds.
The compressed table has SHA-256
`8427e4db27fd165dd8e21535434c81fe2f20349ee5e7c1878c63a187c01d040b`;
its summary has SHA-256
`7c990b5200ec1dbc1a82c5fe263cee85ef60e74fe0929ff4eb8f8fc71fab8692`.

## Complete V100 exhaustion

`scripts/p7_size8_one_elevation_gpu.py` directly un-ranks all
`C(49,8)=450,978,066` finite boundaries on Soulkiller's Tesla V100. It
reproduces the exact 23,563,806-boundary scope, both floor-pair counts, and
the complete odd-secant histogram. An independent `itertools` CPU prefix of
100,000 ranks agrees with CUDA on every count, histogram, and surviving
rank/direction pair.

The CUDA pass tests all 94,255,224 leaves in 1.55 seconds and leaves 1,191
projected leaves on 1,177 boundaries. Host reconstruction with all 135
mod-seven dependencies leaves 1,176 leaves on 1,176 boundaries. The full
run took 132.21 seconds, mostly in source reconstruction, the independent
CPU prefix, and exact candidate catalog rebuilding. The V100 record is
`evidence/p7_size8_four_allocation/four_allocation_cminus1_v100.json`,
SHA-256
`96cfe751a6c0f6bbcd86a1ef799c25847653f8db907414c7b85da576e02efe47`.

## Exact survivor geometry

All 1,176 mod-seven survivors have the same intrinsic form:

`boundary = one affine line L together with one point P not on L`.

The raised direction is normal to `L`. Its odd-fibre count is two; every
other direction has odd-fibre count six, for 44 odd secants in total. For
`c_H=-1`, the line normal belongs to one of the four negative quadratic
direction types. There are four such direction classes, seven parallel
lines in each class, and 42 points off each line, hence exactly

`4 * 7 * 42 = 1,176`.

The independently generated line-plus-point family agrees as a set of
lexicographic rank/direction pairs with all 1,176 full mod-seven survivors.
Every member has exactly two matching rows in its 1,764-row raised catalog.

## Independent two-modulus closure on NUKA

`scripts/p7_size8_one_elevation_audit.py` does not import the CUDA scanner.
It independently rebuilds the direction partitions, integer score matrix,
complete left kernels, catalogs, conditioned projections, lexicographic
ranks, and all 1,191 recorded candidates. The matrix ranks and dependency
dimensions are

| field | rank | left-kernel dimension |
|---:|---:|---:|
| `F_3` | 162 | 120 |
| `F_7` | 147 | 135 |

Fifteen projected leaves fail both complete systems. Each of the 1,176
line-plus-point leaves has 756 matching catalog rows modulo three and two
matching rows modulo seven. The row-index sets are disjoint for every leaf.
Catalog ordering describes exact integer Johnson-slice rows and is the same
in every characteristic, so any genuine solution would have to select one
row belonging to both sets. None does. Thus all 1,191 candidates, and hence
all 94,255,224 leaves, are infeasible.

The audit also regenerates all `4*7*42` geometric candidates and verifies
set equality with the mod-seven survivors. It completed in 22.83 seconds.
The independent record is
`evidence/p7_size8_four_allocation/independent_nuka_audit.json`, SHA-256
`8129b608ec2e09967e10a7da7b38a8e20584450772ac7aab6c1c8a984a370e67`.

## Opposite sign and honest remainder

The same audit reconstructs multiplication by the nonsquare field element
8. It fixes the distinguished edge, is a signed conference anti-isometry,
bijects the eight direction and fibre partitions, and reverses every
quadratic direction type. Simultaneously changing `c_H=-1` to `c_H=+1`
preserves the phase condition and swaps the two type-floor sums. It maps the
negative-type line-plus-point family bijectively to the positive-type
family. Therefore the complete exclusion transfers to `c_H=+1`.

The exact open size-eight floor scope per sign is now

`24,983,238 - 23,563,806 = 1,419,432`.

These remaining boundaries have 11, 16, 24, or 44 exact mean allocations.
No claim about them is made here.

The scripts, proposition record, proof note, conditioned tables, GPU record,
and independent audit are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-four-allocation/`.

## Literature and OEIS context check

[Ball--Csajbók, *On sets of points with few odd secants*](https://arxiv.org/abs/1711.10876)
is the adjacent finite-geometry literature already relevant to the conic
minimum. It does not provide this directional mean-allocation
classification or the two-characteristic catalog exclusion. Searches on
2026-08-27 for 94,255,224, 23,563,806, and 1,419,432 found no matching
mathematical result or relevant OEIS entry; visible exact-number hits were
unrelated financial, administrative, or game data. The intermediate 1,176
is not a new sequence phenomenon: it is also `C(49,2)`, a standard binomial
coefficient visible in [OEIS A011001](https://oeis.org/A011001), while its
role here is explained directly by `4*7*42`. This is a duplication/context
check, not an OEIS submission or a novelty claim.
