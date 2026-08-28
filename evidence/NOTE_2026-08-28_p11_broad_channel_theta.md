# Exact p=11 broad-channel theta reconstruction

Proposition 15.668 refines the scalar shell trace of Proposition 15.667 into
the three eigenspaces of the square-circle operator.  It is an exact finite
`p=11` theorem and an exact negative result about one proof route.  It does
not prove the all-prime R1 inequality.

## 1. The marked profile statistic

For a quartic residue profile `a=(a_s)_{s in F_11}`, the ordinary theta
counter only retains its value histogram.  The square-circle trace also
requires

```text
U4(a) = sum_c (sum_s eta(s-c) a_s)^4.
```

Input-affine permutations preserve this integer.  Exhaustive canonicalization
reduces the `11^4=14,641` zero-constant quartics to 1,007 input-affine types.
Allowing an output-affine transformation reduces those to 20 canonical
dynamic programs; the recorded output scale and translation transport the
additive-character phase exactly.

The complete ten-dimensional glue code is then reduced as follows:

```text
11^10 codewords
  -> 21,437,340 translation/nonzero-scalar representatives
  -> 2,584,901 weighted sorted rich-profile tuples.
```

Their weights reconstruct exactly `11^10=25,937,424,601`.  The V100 output
was independently rebuilt on 268 sampled entries across all orbit strata.

For each profile the dynamic program records the ordinary count, the squared
profile excess, and `U4`.  Five independent 31-bit primes with primitive
eleventh roots give product

```text
31,999,921,744,068,749,461,247,094,447,450,713,426,945,936,557.
```

That exceeds the unrestricted integer bounds for all three marked quantities
through exponent 120.  CRT recovery is therefore deterministic integer
recovery, not a probabilistic residue check.  The ordinary coefficients agree
with Proposition 15.667, and the four already classified shell operators give
an independent calibration of all three channel masses.

## 2. Recovering the three raw masses

For a shell raw operator `R_e`, square-circle evaluation gives a vector `z`.
After projecting away the point-incidence columns, the circle Gram operator is

```text
G = p^3(p+1) I + 2p^2 A_2,
```

where `A_2` joins square circles meeting in two points.  Consequently the two
marked profile contractions `z^T z` and `z^T A_2 z`, together with the scalar
trace, recover exactly

```text
tr(R_e Pi_kernel),  tr(R_e Pi_low),  tr(R_e Pi_high).
```

At `p=11` the dimensions are respectively

```text
1220, 305, 244;     1220+305+244 = dim Z = 1769.
```

Every recovered mass through exponent 120 is nonnegative, as required by
`R_e >= 0`.

## 3. Modular reconstruction

Subtracting the universal radial term and dividing by the channel dimension
puts each broad average in the same 32-dimensional affine modular space used
by the channel q-row exports.  All three channels have the common pivot list

```text
31,32,35,36,...,87,88,91,92.
```

Thus the prefix through exponent 92 has full rank 32.  The remaining 28 exact
profile coefficients, exponents 93 through 120, were not used to solve for
the coordinates and all match.  The unique forms are then evaluated through
exponent 800.  All three reconstructed raw mass series remain nonnegative,
sum coefficientwise to the aggregate raw trace, and have transformed targets
whose dimension-weighted sum is exactly the aggregate target.

The three broad target averages are

```text
kernel = -10463154194187058501821423212 / 56945415059744986998575474182785
low    = -19210249628300203741452825212 / 56945415059744986998575474182785
high   =  -2883758999278296860880307324 / 11389083011948997399715094836557.
```

## 4. Exact endpoint certificates and the route limit

The LP now conserves raw mass on every shell separately in each broad channel
and also conserves each transformed broad target.  Constituents of the same
case are interchangeable, so a distinguished constituent plus the average of
its peers is an exact convex symmetry quotient for one endpoint objective.

QSopt_ex solved all eight rational endpoint problems.  Every reported primal
constraint and every dual stationarity equation was independently recomputed
over `Fraction`.  The resulting harmonic-target intervals are

| case | exact certified decimal interval |
|---|---:|
| circle-kernel principal | `[-522.933314, 508.493608]` |
| circle-low Weil | `[-382.405131, 392.954765]` |
| circle-low principal | `[-257.123541, 268.631365]` |
| circle-high principal | `[-219.228926, 219.228419]` |

These are much narrower than the aggregate intervals in Proposition 15.667,
but the Poisson conversion is exact and monotone, and every interval still
contains target values mapping below `Phi=6`.  Therefore the broad-channel
cone itself cannot certify the spectral floor or R1.  This is a certified
route limit, not evidence that R1 is false.

## 5. The independent finite p=11 theorem

The full Max+ census had already computed

```text
||delta||^2 = 1382747375360 / 583792784981
             = 2.368558520991...
```

whereas

```text
n/12 = 61/6,
exact R1 threshold = 22143/1682.
```

The exact strong margin is

```text
27314875631681 / 3502756709886 > 0.
```

Hence strong R1 is rigorously true at `p=11`.  This conclusion comes from the
complete finite census, not from the broad LP, and supplies no uniform
certificate for `p>=13`.  General R1, mixed-character QVAR, residual (ii),
Type I, and the limit all remain open.

The exact source record is
`evidence/e1_gmin_r1_principal_pge11.json`, SHA-256
`bc78840eadf843db041cc601f949b0f305a1a8027b3920bc95d37029e616a6f4`.

## 6. Literature and OEIS check

Tasaka's survey on [harmonic strength and weighted theta
series](https://arxiv.org/abs/2308.14309) confirms the standard bridge between
lattice-shell harmonic sums and modular-form coefficients.  Ozeki's work on
[Siegel theta series and association
schemes](https://doi.org/10.2206/kyushujm.68.053) is adjacent to the use of
orbital decompositions, but neither source supplies the channel transport
inequality needed here.  The previously requested Bellman-function paper
[arXiv:2305.03523](https://arxiv.org/abs/2305.03523) remains an analogy for a
possible concavity argument, not a theorem about these Paley shells.

Targeted OEIS searches for the large exact census values and the dimension
triple returned no matching sequence.  Generic “circle theta” results lead to
classical square-lattice sequences and are unrelated.  No novelty claim is
made from a negative search.

## 7. Reproduction and archive

The permanent archive is

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1-broad-channel/
```

It contains the rich type table, weighted tuple table, marked dynamic-program
table, exact CRT moment report, reconstructed forms, all eight LP/solution/log
triples, and the final report.  `SHA256SUMS` has 33 entries and verifies in
full; its own SHA-256 is

```text
d1ef69b9af7007c0d2f09a3a5ea8a014cde62d9ed6109175cf4a6496d06b3f07.
```

The principal commands, with the paths used in this run, are:

```bash
python scripts/r1_p11_channel_profile_types.py \
  --output /tmp/r1_p11_channel_profile_types.npz \
  --report /tmp/r1_p11_channel_profile_types.json

python scripts/r1_p11_channel_dual_tuple_gpu.py \
  --output /tmp/r1_p11_channel_dual_tuples.npz \
  --report /tmp/r1_p11_channel_dual_tuples.json

python scripts/r1_p11_channel_profile_tables.py \
  --types /tmp/r1_p11_channel_profile_types.npz --max-k 30 \
  --output /tmp/r1_p11_channel_profile_tables_k30.npz \
  --report /tmp/r1_p11_channel_profile_tables_k30.json

python scripts/r1_p11_channel_moments_gpu.py \
  --tuples /tmp/r1_p11_channel_dual_tuples.npz \
  --profiles /tmp/r1_p11_channel_profile_tables_k30.npz \
  --ordinary-moments /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_profile_theta_moments_e120_exact_20260828.json \
  --max-e 120 --output /tmp/r1_p11_channel_moments_e120.json

python scripts/r1_p11_broad_channel_reconstruct.py \
  --channel-moments /tmp/r1_p11_channel_moments_e120.json \
  --scalar-reconstruction /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_scalar_theta_reconstruct_moments_e88_to800_20260828.json \
  --trace-reconstruction /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_trace_reconstruct_e92_to800_20260828.json \
  --exact-row-directory /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_channel_unscaled_v3_20260828 \
  --output /tmp/r1_p11_broad_channel_reconstruct_e800.json

python scripts/r1_p11_broad_endpoint_qsopt.py \
  --broad-reconstruction /tmp/r1_p11_broad_channel_reconstruct_e800.json \
  --scalar-reconstruction /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_scalar_theta_reconstruct_moments_e88_to800_20260828.json \
  --exact-row-directory /mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/r1_p11_channel_unscaled_v3_20260828 \
  --coefficient-through 800 \
  --output-directory /tmp/r1_p11_broad_endpoint_all_e800
```
