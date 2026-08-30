# Positive p=7 infinity-plus-seven z=7: global join and semigroup audit

## Scope and conclusion

This note audits the positive `p=7`, `z=7`, phase-zero, `c_H=1` branch.
It records a rigorous reduction and 348 rigorous pointed-case exclusions, but it
does **not** close `z=7`, the positive endpoint, or the main theorem.

The committed orbit source has 56 line boundaries in two orbits of size 28
(source orbit indices 102 and 103).  Their 2,160 exact-mean leaves, split into
the two disjoint pointed-star branches `A` and `B`, give 4,320 pointed leaf
cases.  The parent exact affine-hull sieve over moduli 3 and 7 gives

```text
4,320 pointed cases
-3,024 rigorous affine-hull projection rejections
------
1,296 parent survivors
```

Passing that sieve is only a necessary condition: high catalogs are replaced
by their full exact-zero-mean degree-two affine hulls, and a passing modular
right side is not yet a binary edge lift.

## Exact four-case symmetry and transfer

The four compact cases are the two line orbits crossed with pointed branches
`A/B`.  The exact audit proves that `A` and `B` have the same augmented compact
row space within each orbit and that the orbit-0 to orbit-1 maps are

```text
f_A(u) =  8u
f_B(u) = 32u + 24
sigma  = (4, 7, 1, 0, 2, 3, 5, 6)
```

Both affine maps send the first line to the second, induce the displayed
direction permutation, swap the two quadratic direction types and the common
residues 0 and 4, and transport the fixed-edge normalization, parity floors,
14 primitive Johnson equations per direction, exact means, and all compact
constraint families.  Constants and branch-fixed right sides are included in
the augmented rows.  The common augmented ranks and canonical RREF hashes are

| modulus | rank | canonical RREF SHA-256 |
|---:|---:|---|
| 3 | 120 | `4896adb1116b115941a784a3ae3bf5310f2dc73125a70b4a65f5e9e422197ba4` |
| 5 | 114 | `f4ec7f92643c521afe188f77d247b9ae2970d7d398e6e07ae1192a46b3c58810` |
| 7 | 135 | `10c0a4c13efbd4573bd75c2f24225a907eddeccbdbed37aead5ebbd3c174d609` |
| 11 | 114 | `5c9ed22d7a49521cb7323e0b707ce70dc319d4fceb49b745564b8f5b82d853f8` |

The leaf transport is a bijection: all 4,320 pointed leaves form 1,080
four-element classes.  Its permutation hash is
`cc5af4ca13a2d16713a2b84e98f693c9825e20d7dba5844c87c05207e0642bae`.
The 1,296 parent survivors are exactly 324 complete four-case classes, with no
partial class; the parent survivor-key hash is
`f756e0128e12c78ad7a17f85dd621e4e9ff00f0be80c06ac60a71f74045fc784`
and the ordered transfer-class hash is
`57b89949fb776704bb7e0d9e4a62bd9175e402741bd4917f07b7c537b16f7083`.

The transfer used by the global join is stronger than row-space equality.
Across the branch maps, the audit checks all 2,758 `S/M` occurrences as exact
catalog-row-set bijections, using one row permutation simultaneously modulo 3
and 7.  It separately checks 58 high-catalog transport domains.  Every high
direction transports the full rank-20 exact-zero-mean hull, and all 67
nonliteral anchor transports differ from the target anchor by an explicitly
verified hull vector.  That shift either cancels between the base and retained
catalog contribution or is annihilated by high-hull conditioning.  Therefore
the exact same-index mod-3/mod-7 global-join decision transfers throughout
each four-case class.

This certificate does **not** claim that a raw branch-specific binary edge
witness transfers between `A` and `B`; only the audited compact feasibility or
infeasibility statements transfer.

## Exact global mod-3/mod-7 join

One orbit-0/branch-`A` representative from each of the 324 classes was tested.
Every enumerable `U/S/M` direction uses one complete exact catalog row, with
the same row index in both characteristics.  Contributions are concatenated
as raw mod-3/mod-7 bytes, deduplicated exactly, and joined by exact byte
comparison; no digest-collision assumption enters a decision.  High
directions remain relaxed only to their audited full affine hulls.

With a six-million-state maximum projected side, all 246 budget-eligible
representatives were processed (236 with enumerable catalogs and 10 with no
enumerable direction).  The other 78 were recorded as explicit budget skips:

| representative decision | representatives | four-case transfer |
|---|---:|---:|
| rigorous empty global join | 87 | 348 |
| necessary-only join survivor | 159 | 636 |
| explicit side-budget skip | 78 | 312 |
| **total** | **324** | **1,296** |

Thus the rigorous new content is the exclusion of 87 representatives, or 348
of the 1,296 parent pointed survivors.  The remaining 636 survivors and 312
skips are not feasibility claims.  Internal certificates are

```text
case results          c34cc913c27910e3876e1b78aed0e9c8c2f42cb2f4368f95054bcd6ead1db7a7
rejections            03ef66e6b4529d05c0351762dd8a190eb9b400e3b050852fafa8f2ea16d0da78
survivor witnesses    b52fbafed16a04514ce7d403d562e39d12fde0ac6ba31845caf3e20b8c8a8dd0
transferred decisions 7b045433541f00e6b866e443ccd23187906550bfe77c1985ee512ad04a41f210
```

### Hostile audit

An independent read-only hostile audit passed the rejection path.  It traced
coverage, budget classification, exact raw-byte deduplication, meet-in-the-
middle signs, witness reconstruction, hashes, and the symmetry transfer; it
also compared 250 randomized small joins with direct brute force and sampled
real cases.  It found no false-negative route.  The built-in manufactured
tests additionally include positive and impossible joins and the trap in
which mod 3 and mod 7 each match separately but no shared catalog row matches
both.  The review certifies the 87 rejections, not any survivor.

The hostile review emitted no separate artifact.  Its initial provenance
hardening concern was addressed by the explicit rerun archived below, which is
bound to global-join source SHA-256
`86ca9a8055ba20f129f284b6e9001478880a56f20185b701c708000577f34bd8`
and the explicit symmetry certificate.

## Why moduli 5 and 11 add zero information here

The four-prime audit is marked as a bounded smoke run because it executes only
five final joins.  Its no-go preflight, however, exhausts all 159 mod-3/mod-7
join survivors and all 312 distinct complete-catalog signature records, so the
following conclusion does not depend on those five sample joins.

All four pointed systems were reconstructed directly at moduli 3, 5, 7, and
11.  For branches `A/B`, their ranks are respectively

```text
             mod 3   mod 5   mod 7   mod 11
branch A       162     168     147      168
branch B       169     175     154      175
dependencies   120     114     135      114
```

Rank-nullity, complete left-null identities, and manufactured consistent and
inconsistent right sides were checked in every case.  The exact integral
degree-two evaluation lattice is saturated of rank 21, and its zero-mean
sublattice is saturated of rank 20.  Modulo 5, the primitive equations alone
have a misleading nullity 21; Smith saturation proves that the extra
one-dimensional direction is non-liftable.  The liftable zero-mean space is
still exactly the rank-20 reduction of the integral lattice.

After conditioning on every high-direction set occurring among the 159
candidates, moduli 5 and 11 both retain all 114 rational dependencies and see
every high-hull image as zero.  Every candidate base has zero projected
mod-5/mod-11 component, and every complete `S/M` catalog contribution has zero
projected component there (201,936 catalog rows audited across 312 signature
records).  Consequently the four-prime unique-state counts equal the original
mod-3/mod-7 counts for every candidate, and the four-prime quotient is not
strictly stronger on this universe.  Further repetition of this same
mod-5/mod-11 compact test cannot reject a remaining case.

## Exact Johnson semigroup structure

Let `K` be the primitive 14-by-35 integer kernel cutting out the Johnson
degree-two functions on a direction fibre.  Its int64 SHA-256 is
`f90bafe2de158e7a4f08cd32e603631bfa0dab7fdf4fcae2bb013642294d7ccb`.
For parity floor `P_b`, every exact directional slack can be written

```text
A = P_b + 2L,       L in Z_{≥0}^35,       K L = 0,
P_7 = 0,            P_1(X) = 1_{u in X}.
```

Both floors lie in the Johnson degree-two space.  The score box
`0 <= A_X <= 13` is exactly `0 <= L_X <= 6`.  With grade
`g = sum_X L_X / 5`, define

```text
C_g = {L >= 0 : K L = 0, sum_X L_X = 5g, L_X <= 6}.
```

Then the `b=7` catalog at scaled mean `a` is `2 C_(a/4)`, while the
`b=1` catalog is `P_1 + 2 C_((a-8)/4)`, up to the already audited fibre
coordinate permutation.  Hence `U/S/M` are exactly grades `0/1/2`; the high
catalogs use grades `3,4,5,6,8`, and required grade 7 never occurs.

Normaliz computes the Hilbert basis of the uncapped semigroup

```text
S = {L in Z_{≥0}^35 : K L = 0}.
```

It has 896 distinct binary generators and no primitive generator above grade
3:

| grade | Hilbert-basis generators | exact identification |
|---:|---:|---|
| 1 | 56 | exactly the lifted `S56` catalog |
| 2 | 168 | the indecomposable part of `M1764` |
| 3 | 672 | a genuinely new primitive family |

More precisely, `M1764` is the disjoint union of 1,596 distinct sums of two
grade-1 generators and the 168 primitive grade-2 generators.  At grade 3,
30,436 distinct triple-`S` sums and 8,428 distinct grade-1 plus primitive-
grade-2 sums overlap in 1,680 rows, giving 37,184 decomposable rows; adjoining
the 672 new atoms gives the complete 37,856-row layer.  Therefore `S56` and
`M1764` alone do not generate the high catalogs, but the full 896-element
Hilbert basis does.

The Normaliz Hilbert series has denominator
`(1-t)^15 (1-t^2)^5 (1-t^6)` and gives these exact **uncapped semigroup**
layer counts:

| grade | rows |
|---:|---:|
| 0 | 1 |
| 1 | 56 |
| 2 | 1,764 |
| 3 | 37,856 |
| 4 | 575,407 |
| 5 | 6,496,938 |
| 6 | 57,232,105 |
| 7 | 410,200,367 |
| 8 | 2,474,264,653 |

The coordinate cap matters.  Because every generator is binary, every layer
through grade 6 automatically satisfies `L_X <= 6`, so the required layers
through grade 6 equal their capped catalogs.  At grade 8, arbitrary semigroup
sums can violate the cap—for example, eight copies of one grade-1 generator
have coordinate 8.  The grade-8 Hilbert-series count is therefore only an
uncapped outer count, and any support enumeration must track or filter the
cap.  This does not weaken generation of an actual capped target: in a
nonnegative Hilbert decomposition every summand and partial sum is
coordinatewise bounded by the final target.  It only forbids treating the raw
grade-8 layer as the exact catalog.

The trust boundary is explicit: Normaliz certifies Hilbert-basis completeness.
The audit independently reparses all 896 rows, verifies distinctness,
nonnegativity, binarity, `K L=0`, all grades, the complete `S56/M1764`
identifications, and the displayed Hilbert-series expansion.

## Torsion-support calibration and finite projected stabilization

For the orbit-0/branch-`A` 282-by-1,225 pointed matrix, 112 Johnson relations
plus two quadratic-type sum relations give 114 independent exact rational
left dependencies.  A nonzero mod-11 minor pins the rational rank at 168.
The complete characteristic-3 and characteristic-7 dependency dimensions are
120 and 135, so quotienting by the rational dependencies gives the effective
torsion quotient

```text
F_3^6 x F_7^21.
```

The four `H0_S0_M7` representative skips (leaf indices 405, 411, 413, and
414) are a clean calibration because all seven variable directions use the
complete `M1764` catalog and no high hull is present.  Exact projected-support
convolution retained all six mod-3 coordinates and selected mod-7 coordinates:

| projection | group size | tested projections | result |
|---|---:|---:|---|
| `k=2`: `F_3^6 x F_7^2` | 35,721 | 16 (4 per case) | every support saturated the full group; every target present |
| `k=3`: `F_3^6 x F_7^3` | 250,047 | 32 (8 per case) | every support saturated the full group; every target present |
| `k=4`: `F_3^6 x F_7^4` | 1,750,329 | 64 (16 per case) | every support saturated the full group; every target present |
| `k=5`: `F_3^6 x F_7^5` | 12,252,303 | 20 (5 per case) | every support saturated the full group; every target present |
| `k=6`: `F_3^6 x F_7^6` | 85,766,121 | 4 (1 per case) | every support saturated the full group; every target present |

The `k=2/3` runs use exact CPU support sets.  The `k=4/5/6` runs use the
self-audited exact CUDA open-addressing engine: hash collisions are resolved
by key comparison and probing, and no digest or probabilistic membership test
enters a decision.  No completed projection hit its state cap and no case was
rejected.  These are exact necessary-condition calibrations, not evidence of
binary feasibility.  Their useful conclusion is negative: coordinate-only
torsion projections through `k=6` saturate too quickly to separate these four
hard cases, so this `H0_S0_M7` projection family is closed as nonseparating.

### Finite stabilization lemma

There is also a positive structural consequence for the high catalogs.  In a
fixed finite projected group `G`, let `A_i` be the projected Hilbert generators
of grade `i` and define the uncapped exact-grade recurrence

```text
T_0 = {0},
T_g = union_(i=1,2,3) (T_(g-i) + A_i).
```

The maximum Hilbert-generator grade is three.  If four consecutive supports
are equal,

```text
T_r = T_(r+1) = T_(r+2) = T_(r+3) = S,
```

then `T_g=S` for every `g>=r`, and `S` is the subgroup of `G` generated by
all `A_i`.  Indeed, the recurrence for `T_(r+3)` shows `S+A_i` is contained
in `S` for every generator.  Translation is injective on finite `G`, so each
containment is equality.  Thus `S` is invariant under the generated subgroup;
conversely every element of `S` is a generator sum.  Induction gives permanent
stabilization.

The complete `k=3` and `k=4` support artifacts reproduce the stronger equality
of both raw and anchor-relative support hashes at grades 3, 4, 5, and 6 in all
eight directions.  The stabilized subgroup orders are

| projection | directions 0--4 | direction 5 | directions 6--7 |
|---|---:|---:|---:|
| `F_3^6 x F_7^3` | 147 | 3 | 3 |
| `F_3^6 x F_7^4` | 1,029 | 21 | 3 |

Because every Hilbert generator is binary, uncapped and capped catalogs agree
through grade 6.  Hence this subgroup statement is exact for every required
high catalog of grades 3--6, not merely an outer relaxation.  At grade 8 it
remains only an outer-support statement: projection has discarded the
coordinate data needed to enforce `L_X<=6`.  Missing projected membership can
still reject a grade-8 case, but projected presence cannot certify a capped
lift.

### Terminal case-join audit

The bounded case-join campaign did not turn the projected structure into a
rejection.  All 51 grade-three-only representatives survived without skips
under the eight selected `k=3` projections, the completed one-projection
`k=4` prefix, and the completed one-projection `k=5` prefix.  The `k=5` run
was a disjoint `26+25` V100/Orin shard cover.  The even-index `k=6` shard also
retained all 26 selected representatives without skips; its odd-index mate
was cancelled, so there is no complete `k=6` claim.

One final fixed-case gate removed projection selection entirely.  For
`orbit0_leaf780_branchA`, CP-SAT imposed the 112 Johnson kernel equations,
eight mass equations, coordinate bounds, and all 27 coordinates of the full
`F_3^6 x F_7^21` quotient in a 307-variable, 147-constraint model.  After
300.195 solver seconds, 95,634 conflicts, and 1,179,303 branches, it returned
`UNKNOWN`.  By the audited decision semantics, only `INFEASIBLE` would be a
rigorous quotient rejection and only an independently audited witness would
be a necessary survivor; `UNKNOWN` has no mathematical force.

These computations add no proposition and remove no boundary.  No further
seed, projection, encoding, or solver variant is scheduled: this
Johnson-semigroup/quotient attack is terminated as non-closing.

## Unresolved mathematical scope

After the global join, 237 representatives remain unresolved: 159
necessary-only survivors and 78 explicit budget skips, representing 948
pointed cases under the proved fourfold transfer.  Organized by the largest
high-catalog grade needed, they are

| catalog tier | representatives |
|---|---:|
| `H0_S0_M7` torsion calibration cases | 4 |
| grade-3-only high catalogs | 51 |
| maximum high grade 4 | 137 |
| maximum high grade 5 | 24 |
| maximum high grade 6 | 13 |
| cap-sensitive grade 8 | 8 |
| **total** | **237** |

Accordingly, the 56 positive `z=7` source boundaries remain open.  The
negative `p=7` branch, later-prime endpoint work, `R1`, the non-Walsh
multi-level remainder, and the main theorem are outside this certificate and
remain open.

## Artifact and source manifest

All hashes below are SHA-256 values read from the files named here.

```text
/home/nick/quadratic-minmax-limit/evidence/p7_infinity7_positive_zge2_orbits.json
  4cf79ffb19b0e6961f6f187efe513d0fc5503ff2f0c63317603f4b0205d7cd63

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/p7_inf7_z7_affine_hull_multimod_full.json
  fa1f1e52d24389a9863274cfa6d2b251a4b06e5e9f2b05e624a0b2587ec65f79

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/global-join/p7_z7_symmetry_global.json
  0df450d3be1b6897d2bf5e7ae55473ad8c53facffdab86256745d37d5e49cdf3

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/global-join/p7_z7_global_join_full_explicit.json
  5d3bb1a3385fe848932f3405032eb45501741afc8a3aef36d76b126f5859b93a

/tmp/p7_z7_global_catalog_join_fourmod_smoke.json
  e59f5609a74ee84ed2e020003ee912ad72b930c36c8d63d98621eeb29e965bcd

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_johnson_semigroup_rebuilt.in
  14b670086a068da3f07c3bc1fb6c8ece339ae9dd78a6b1e412219a286cf45a5f

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_johnson_semigroup.gen
  3b582d6a0e7c83cb8ed41a421e4950be2645ce2d3aa18dab17432628b787b789

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_johnson_semigroup_hs.inv
  867279bb23320ab7ea9f03dc41765b8ebddcfef7b041c4ee2a1373085d6a42dd

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_johnson_semigroup_hs.out
  92e8a4f648b9cc13395e392164dd0d7e7632e0027cd2954cf7915622a3e31bdf

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_johnson_semigroup_hilbert_audit.json
  09c9ffaeace3eab5fa3bcca99a78d2d63496c6a23bf7f1b5f162769395798a24

/tmp/p7_z7_torsion_smoke.json
  a54fad1c270eb20656e2455bc2fa65b1e0f65d41c748f1dd905e6334bb570c6f

/tmp/p7_z7_torsion_k3.json
  05e4694b7f9b09b7a4fce591b1ad084099a8bf8398ea99171711e421a16a1e0a

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_high_semigroup_k3.json
  e9cd5f24a9ec6961404b95e1240abb44a21d0ef4e5d78c9fc0324f643edce079

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_high_semigroup_k4.json
  8a3f58c2728421458821dc6207928836c085b49ea35018abd2c8587b7adacf53

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_high_semigroup_k5.json
  33dadc47dfb8c4016ef79e1ae04c6c6711c84c31ce972a2a191c6f1084fcf9fd

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_high_semigroup_k6.json
  32b660c5b47523c96fa5b81e2357c9752a5bb62749d5c84c9e3e94463562ffd1

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_torsion_gpu_k4_orin_sixteen.json
  581d47b0ac0ecdebc349d3486ea86ada593fd85b71bd70efb70438efda7a4b3d

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_torsion_gpu_k5_v100_five.json
  50cf489fceb527c55b676ef350b4cf969bbeb7320caa4d91aaa1cf06c413d8e4

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_torsion_gpu_k6_v100_fullcap.json
  8dec55818d3913b0c2e50f34341c8424b80b18fa3e5778d8803ef623d4fcf007

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_smoke.json
  c41fc908c58f27536223d9c824b63394c6f47ed3a19c8a23a970183ab14379ed

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_all51_k3_auto8.json
  bb9359af2e09e6a7700d74853a9e18be74e30077edde0c854127ab2b9638b7d3

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_all51_k4_prefix.json
  b5bdec818852fb2d5ec8039941dbcd25ae71bf792f81a3b55e19ad9520028e94

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_gpu_k5_shard0_v100.json
  738b0dc44c1cfd8b75a723e4841535ca01dabe93ec1053bf6987f1705e608667

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_gpu_k5_shard1_orin.json
  fa966371e4387bf5dc8b4295eebed7cccaeb84dcad10bf5785f92d4474745be6

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_semigroup_case_join_gpu_k6_shard0_v100.json
  f7cb5123237ecab4f1adf9b3bffd533f00d4e7852e461e49a9fb7edd5983ff11

/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7/semigroup/p7_z7_direct_semigroup_case0_300s.json
  5051cefee62a0ea6f4e7c38235df6bfa54a8f2a7ade80d06d046470ba95ab6f0
```

The archived Hilbert audit includes the explicit grade-3 decomposition counts
and grade-8 cap warning recorded above; the underlying Normaliz input,
generator file, and Hilbert series are unchanged.

Current local source hashes are

```text
scripts/p7_infinity7_positive_z7_pointed_affine_hull_multimod.py
  0ae46bf3a0ad64975b7e0ac55aea562c9efcdee43828cb205e4e7cabd196d8b8
scripts/p7_infinity7_positive_z7_compact_symmetry_audit.py
  42ce41934612e0cde03338744975f470ec1c507816b1ee29dea000d19d73950c
scripts/p7_infinity7_positive_z7_global_catalog_join.py
  86ca9a8055ba20f129f284b6e9001478880a56f20185b701c708000577f34bd8
scripts/p7_infinity7_positive_z7_global_catalog_join_fourmod.py
  a3bbb8b7dcc47d5acf89da888e214a468053b3f77ae61e68f66b88ff9cbabd3d
scripts/p7_infinity7_positive_z7_johnson_semigroup_normaliz.py
  465ecf21bfd586379959a292493e717453a81c672e8f40d61e381ab8fb4ed0c5
scripts/p7_infinity7_positive_z7_torsion_support_projection.py
  520f8ca73dd60659e63745c4aebe91c5dcb6bd9ae3b48e9ad98650a7cacaa8dd
scripts/p7_infinity7_positive_z7_torsion_support_gpu.py
  e3a12c05c992ccd7acb1b470e9dcd63e9e3836a3c2e0ae338ca40fbd3c4c9820
scripts/p7_infinity7_positive_z7_high_semigroup_support.py
  718bc28c3f62395b0c581d9773ba42bcef1f2b05d8f6e10fdd11b254703ce256
scripts/p7_infinity7_positive_z7_semigroup_case_join.py
  7fa61f570be56692e1bc6e53e3fe8cc2824893f5a350fd57210b054d6049d55b
scripts/p7_infinity7_positive_z7_semigroup_case_join_gpu.py
  c05311a062a5cf5426e6f49b1ee011947d484ed45ccbc5274e841f459135c89c
scripts/p7_infinity7_positive_z7_semigroup_case_cpsat.py
  c7c9098395d97b3dabc25ed7bd03b446c81c9afb278883926f225818c105f827
```
