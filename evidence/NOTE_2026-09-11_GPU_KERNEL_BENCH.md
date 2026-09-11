# GPU kernel bench: V100 vs RX 9070 XT on the scout Phi kernel

**Status:** draft session note (2026-09-11). Timings from identical scripts on
identical objects; cross-engine float32 agreement recorded. Independent review
pending. No asymptotic claim; finite evidence only.

## Object

The scout scoring kernel (`scripts/bench_phi.py`):

    Phi(K) = max over the 2^(2n-1) projective states S of |sum((S @ K) * S)/2|

float32, chunked at 2^21 states per chunk, `Q = sum((S @ K) * S, axis=1) * 0.5`
(the campaign scoring pattern). Bench objects are the doubling
`K = [[A, Z],[Z^T, -A]]` with the fixed deterministic `Z` (seed 999), built from:

- `K15_minimizer.json` (order 15, `Phi = 27`) -> order-30 doubled, 2^29 states
- `K16_phi30.json` (order 16, `Phi = 30`) -> order-32 doubled, 2^31 states
- `K17_phi34.json` (order 17, `Phi = 34`) -> order-34 doubled, 2^33 states

## Timings (best of reps, seconds)

| node       | device            | runtime             | 2^29  | 2^31   | 2^33   |
|------------|-------------------|---------------------|-------|--------|--------|
| soulkiller | V100-SXM2-16GB    | cupy 14.2.0 (CUDA)  | 3.675 | 15.769 | 69.278 |
| nuka       | RX 9070 XT        | cupy 13.5.1 (ROCm)  | 3.477 | 14.808 | 66.544 |
| soulkiller | 88-thread CPU     | numpy 2.5.2         | 230.3 | -      | -      |

Raw rep lists: V100 `[4.069, 3.677, 3.675]`, `[16.159, 15.769, 15.782]`,
`[69.637, 69.278]`; 9070 XT `[5.463, 3.494, 3.477]`, `[16.383, 14.83, 14.808]`,
`[68.042, 66.544]`; CPU `[230.307]`. First-rep overhead is cupy kernel/context
warmup (amortized in campaign runs).

Cross-engine agreement: the `phi` values are identical to float32 print
precision on all three objects on both GPUs
(105.97763061523438, 118.06329345703125, 129.4598388671875). This doubles as an
independent hardware/software-stack re-verification of the K15, K16 and K17
readings ("replayed on the RX 9070 XT" is now literal with numbers).

## Reading

- On this kernel the RX 9070 XT with ROCm 10.0 is at parity with the V100
  (about 4-6% faster) at roughly 125-150 M states/s; the 88-thread CPU is
  ~60x slower.
- The 9070 XT's theoretical fp32 advantage does not materialize here: the
  kernel is launch/memory-bound (skinny 2^21 x 34 @ 34 x 34 matmuls), not
  FLOP-bound. One integer- or shape-heavy kernel class may still favor RDNA4;
  this bench closes the question only for the enumeration scoring pattern.
- Dispatch policy: use either GPU for enumeration-bound campaign work; split
  by node convenience and seed disjointness, not by speed. Orin and the A380
  were not benched for this pattern (edge GPUs, not on the campaign path).

## Repro

`python bench_phi.py --matrix K16_phi30.json --doubled --reps 3` on each
engine (`~/.venvs/mo-exact` on soulkiller, `~/.venvs/rocm72` on nuka);
objects in `~/scratch/`.
