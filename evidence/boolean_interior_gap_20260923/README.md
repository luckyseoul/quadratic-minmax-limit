# Finite interior-gap algebra receipt

The analytic statement and complete author proof are in
[the proof note](../NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md).
This receipt corroborates local algebra only; it is not a finite census,
an all-orders computational certificate, or independent proof review.

## Run

- Date: 2026-09-23 UTC.
- Baseline: main at `d9705fced9dbaad096b4b2a09599ae54105bf9af`.
- Host: NUKA, selected for a small serial exact-algebra check. No GPU or
  parallel worker pool was needed; the controller did not execute the test.
- Live preflight at 07:49:37 UTC: 16 threads, load 0.00/0.00/0.07,
  about 10 GiB available memory, SymPy 1.14.0. No competing math job
  appeared in the top-CPU process listing.
- Fresh staging: `/tmp/mo-interior-gap-20260923.FBxiMxek`.
- Both source hashes matched on the remote host before execution.
- Command, from that directory:

  ```sh
  PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 check_identities.py
  ```

- Result: 15 checks passed; process exited 0. There is no running job
  from this check. The emitted JSON is retained in `result.json`.

## Inputs and review limits

| Input | SHA-256 |
| --- | --- |
| Proof note | `0048e75faea73b4882e82b1ce485cecbd8391b48e349b802342aa1144b07ad1d` |
| Checker | `498af0d728f5edb4f47183d21825622414156cd9bc198b7eac75e36f884de92c` |
| Result JSON | `89257922a17dd1ee870d365813e81434cd5a5cdcf5f85584ce1f11185c0fb483` |

The checker covers the fourth-moment coefficients, split-size arithmetic,
quadratic optimizer value, symmetric perturbation identity, layer-cake
constant, mean-update scaling and rearrangement, and one retained C5
zero-field boundary point. C5's norm is reused, not re-enumerated.

The root author read the complete proof and checker before the run and
checked the selection argument, one-sided extrema, restrictions, k=0/1
exceptions, and normalization. This is author review, not a claim of
independent human or proof-assistant verification. Holder and the
all-orders probabilistic argument are proved in prose, not by this script.

The September 6 archive already recorded an asymptotic cubic one-sided
estimate. The new proof explicitly acknowledges it; the finite bound
and mean-update correction do not turn its old cubic posterior control
into the missing linear comparison.

This is a supporting milestone, not a major research or structural
checkpoint, so the standing backup rule does not require another full
large-drive snapshot. No old artifacts were deleted or overwritten.
Convergence and the post-update law estimates remain open.
