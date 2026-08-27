# Complete finite `p=7` size-eight exclusion

Date: 2026-08-27. This is Proposition 15.666.

Proposition 15.666 closes every **finite** size-eight boundary at `p=7`, for
both Paley product signs. It does not close the distinct
infinity-plus-seven profile, residual (ii), Type I, R1, global QVAR, or the
final limit.

## 1. Exact scope after Proposition 15.664

The complete `C(49,8)=450,978,066` finite-boundary floor census and
Propositions 15.662--15.664 leave exactly 1,419,432 boundaries per sign.
Their exact mean-allocation counts are:

| allocations per boundary | boundaries | allocation leaves |
|---:|---:|---:|
| 11 | 154,056 | 1,694,616 |
| 16 | 1,194,816 | 19,117,056 |
| 24 | 1,176 | 28,224 |
| 44 | 69,384 | 3,052,896 |
| **total** | **1,419,432** | **23,892,792** |

`scripts/p7_size8_remaining_allocation_structure.py` reconstructs this
partition from the exact floor profiles. Each leaf raises a known support of
one to five direction catalogs. After those raised catalogs are omitted, at
most one 36-row floor catalog remains.

## 2. Exact conditioned omission scans

The common integer score system has shape `282 x 1225`. Its ranks and left
dependency dimensions are

| characteristic | rank | left dependencies |
|---:|---:|---:|
| 3 | 162 | 120 |
| 7 | 147 | 135 |

For every raised support, `scripts/p7_size8_multi_elevation_tables.py`
constructs the exact dependency subspace that vanishes on every raised
35-column direction block. Forty independent rows are selected and audited.
The minimum conditioned dimensions are 44 in characteristic three and 42 in
characteristic seven.

`scripts/p7_size8_remaining_gpu.py` then checks all `C(49,8)` boundaries and
all 23,892,792 leaves. Independent `itertools` CPU prefixes of 100,000 ranks
agree exactly with direct-rank CUDA in each characteristic.

| necessary test | survivor leaves by the four strata | total |
|---|---:|---:|
| mod 7 omission | 150,528; 14,294; 5,880; 288,120 | 458,822 |
| mod 3 omission | 363,384; 1,906,296; 0; 402,192 | 2,671,872 |
| same leaf passes both | 77,616; 0; 0; 103,488 | 181,104 |

The last line intersects boundary-rank/stratum/leaf triples, not modular
catalog witnesses. It is only a necessary reduction, so differing witnesses
between characteristics cannot cause a false exclusion.

## 3. Restoring local catalog compatibility

For variable support `V` and tested subset `T`, an exact mod-seven dependency
space is conditioned to vanish on `V\T`. A zero join for the catalogs in `T`
therefore rules out the complete leaf. Twenty-two independent dependency
digits are packed losslessly into one base-seven unsigned 64-bit integer:

`7^22 = 3,909,821,048,582,088,049 < 2^64`.

No probabilistic hash is used. The packed integer is the complete selected
signature, and every search compares that integer exactly.

`scripts/p7_size8_subset_catalog_gpu.py` applies three bounded stages:

| stage | input | survivors | CPU/GPU prefix |
|---|---:|---:|---:|
| local pairs / positive triples | 181,104 | 124,745 | 256 |
| all relevant catalog triples | 124,745 | 78,126 | 32 |
| four-positive catalog join | 78,126 | 62,892 | 64 |

Every prefix agrees exactly. These stages are reductions only; no closure is
claimed until the complete join below.

## 4. Single-filtered complete catalog join

`scripts/p7_size8_full_catalog_filtered_gpu.py` constructs, for every
variable direction `d`, dependencies that vanish on `V\{d}`. They filter the
rows of catalog `d` without making any choice in the other variable catalogs.
Each surviving exact row retains its signature under one common independent
22-row full-dependency projection.

A CUDA block then joins **all** variable catalogs of that leaf (three or five)
with a shared-memory meet-in-the-middle table. The table uses exact open
addressing on the full packed key. The observed maximum build side is 1,764,
below both the imposed 1,800-row cap and the 4,096-slot hash capacity; the
maximum probe side is 2,744. Thus no table overflow or dropped insertion is
possible. A separately instrumented rejection counter reports zero leaves
without an admissible hash partition, so none of the exclusions below is a
capacity fallback.

The complete accounting is

| outcome | leaves |
|---|---:|
| at least one isolated catalog has no allowed row | 3,777 |
| isolated rows exist, but the complete catalog join is empty | 59,115 |
| complete-join survivors | **0** |
| total | **62,892** |

The final run materializes 1,439,451 distinct
`(isolate signature, full signature)` pairs. An independent CPU
implementation agrees with CUDA on the first 512 candidates. The older
multi-characteristic full-join engine also independently returns zero on one
representative from each of the stratum-11, three-catalog stratum-44, and
five-catalog stratum-44 classes, with raw catalog sizes respectively
`(62,56,36,56,56)`, `(2233,2233,1764)`, and
`(56,1764,62,62,56)`.

Consequently no one of the 23,892,792 exact mean-allocation leaves can satisfy
the score equations for `c_H=-1`.

## 5. Opposite sign

Propositions 15.662--15.664 already audit multiplication by a nonsquare in
`F_49` as a signed conference anti-isometry. It bijects the complete finite
size-eight floor sets, reverses the quadratic direction types, and exchanges
`c_H=-1` with `c_H=+1` while preserving the phase and score constraints.
Therefore the zero-survivor result transfers to the opposite sign.

Every finite `p=7` size-eight boundary is excluded for both signs.

## 6. Evidence and archive

The compact committed chain is under `evidence/p7_size8_complete/`. The final
512-prefix capacity audit has SHA-256
`428b9604e21738d9b063f0edee8a42b31d471ecd56800e4366af8ed1d7a49eaa`.

The 121 MB mod-three survivor record and the other raw stage records are not
placed in Git. They are archived with pinned hashes at

`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-size-eight-complete/`.

## 7. Literature and OEIS context check

[Ball--Csajbók, *On sets of points with few odd secants*](https://arxiv.org/abs/1711.10876)
is the adjacent finite-geometry literature used in the preceding conic
branch. Searches of the arXiv finite-geometry and Paley-graph literature found
no directional mean-allocation or exact catalog exhaustion matching this
finite `F_49` computation. This is a context/duplication check, not a broad
novelty claim.

The input `450,978,066=C(49,8)` is already the eighth nontrivial entry of
[OEIS A017765](https://oeis.org/A017765). The intermediate numbers 181,104
and 62,892 occur in unrelated OEIS entries (for example A004143/A037972 and
A247345), with no connection to this catalog chain. Searches for 1,419,432
and 23,892,792 found no relevant sequence interpretation. No OEIS submission
is proposed.
