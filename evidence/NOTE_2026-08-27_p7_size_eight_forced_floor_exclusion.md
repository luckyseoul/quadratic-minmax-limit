# Complete p=7 size-eight forced-floor exclusion

Proposition 15.663 closes, for both product signs, the finite `p=7`
size-eight stratum whose two quadratic direction types both have exact floor
sum 32. It eliminates 83,770,008 nonconic floor survivors per sign. It does
not close all size-eight boundaries, residual (ii), Type I, R1, global QVAR,
or the limit.

## Exact scope and rigidity

The complete Proposition 15.662 floor files each check
`C(49,8)=450,978,066` boundaries and contain 108,754,569 survivors in 5,152
ordered odd-fibre profiles. Reading those profiles with the exact parity
floors selects 2,016 profiles and 83,770,008 boundaries per sign with type
floor sums `(32,32)`. Their odd-secant distribution is

| odd secants | boundaries per sign |
|---:|---:|
| 16 | 254,016 |
| 20 | 8,396,640 |
| 24 | 32,673,984 |
| 28 | 30,465,456 |
| 32 | 10,459,344 |
| 36 | 1,467,648 |
| 40 | 52,920 |

The histogram sums to 83,770,008. In particular, this stratum is disjoint
from the eight-odd-secant conics closed by 15.662.

The exact directional mean sums are also 32 in each quadratic type. Every
directional mean must therefore equal its parity floor. At these floor
means, every complete Johnson-slice slack catalog is a singleton except the
phase-one, four-odd-fibre, scaled-mean-14 catalog, which has 36 rows. Two
such catalogs would already cost 28 in one four-direction type, while each
of the other two directions costs at least six, exceeding the type budget
32. Thus a boundary has at most one nonsingleton catalog.

## Complete CUDA exclusion

The common affine score equations consist of edge count, the distinguished
edge, and 280 exact bad-edge equations, giving a `282 x 1225` integer
matrix. Its rank over `F_7` is 147, so its complete left dependency space
has dimension 135.

`scripts/p7_size8_forced_floor_gpu.py` directly un-ranks and checks all
450,978,066 boundaries on Soulkiller's Tesla V100. Eight deterministic
linear combinations of the 135 dependencies are used only as a rejection
prefilter. A true solution of all 135 equations must pass every such linear
combination, so this reduction loses no possible solution. The projected
sweep leaves 526 ranks. Host reconstruction then tests every allowed row of
the sole possible 36-row catalog against all 135 dependencies and leaves
zero survivors.

The CUDA pass itself took 0.197 seconds; the complete run took 121.33
seconds because it rebuilt and audited the exact dependency basis before
launch. The complete result is
`evidence/p7_size8_forced_floor/forced_floor_cminus1_v100.json`, SHA-256
`6143d4eb269861b3d380c53262b534e0a54a9645c9bbe7c29d9327200ae30535`.
The record explicitly keeps `closes_all_p7_size_eight=false` and
`closes_residual_ii=false`.

## Independent NUKA replay and opposite sign

`scripts/p7_size8_forced_floor_audit.py` does not import the CUDA scanner.
On NUKA it independently rebuilt the direction tables, the `282 x 1225`
matrix, its rank-147 row space, the complete 135-row left kernel, and every
Johnson catalog used by all 526 projected candidates. The dependency basis
has SHA-256
`33f13cc757a528c64f9beb12004838cd4ea63f36a1af4e42daabc8590d26cdb1`.
All 526 candidates pass the recorded eight-row projection and all 526 fail
the complete dependency system. The replay took 10.54 seconds.

The same audit reconstructs multiplication by the nonsquare field element
8. It gives the direction permutation `(4,7,1,0,2,3,5,6)`, reverses every
quadratic direction type, and permutes the seven fibres within each
direction. Odd-fibre counts are preserved. Simultaneously changing
`c_H=-1` to `c_H=+1` preserves the phase condition, while the two type sums
swap; hence the `(32,32)` stratum is mapped bijectively between signs. The
signed conference anti-isometry fixes the distinguished edge and preserves
the normalized feasibility problem, so the exclusion transfers.

The independent audit is
`evidence/p7_size8_forced_floor/independent_nuka_audit.json`, SHA-256
`7adaa5e76bf4f5e128c82ec219650b390c8c087d3aed2a44857f9da7939a9c53`.

## Honest remaining size-eight scope

Proposition 15.662 excluded 1,323 conic floor survivors per sign. Proposition
15.663 excludes a disjoint 83,770,008 nonconic survivors. Therefore the
remaining floor scope is exactly

`108,754,569 - 1,323 - 83,770,008 = 24,983,238`

boundaries per sign. No claim is made about those 24,983,238 cases.

The code, compact evidence, proposition record, and proof note are archived
under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-forced-floor/`.

## Literature and OEIS context check

[Ball--Csajbók](https://arxiv.org/abs/1711.10876) and the odd-secant
references already recorded for 15.662 concern geometric lower bounds and
the minimum/conic regime. They do not supply this `(32,32)` directional
floor classification or the exact modular-catalog exclusion. Exact web and
OEIS-domain searches on 2026-08-27 for 83,770,008 and 24,983,238 returned no
matching mathematical sequence or prior result; the visible exact-number
hits were unrelated database and listing identifiers. This is only a
duplication/context check, not an OEIS submission or novelty claim.
