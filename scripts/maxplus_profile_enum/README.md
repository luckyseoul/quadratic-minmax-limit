# Max+ profile enumeration toolkit (Prop 15.588)

Enumerates Max+ at any prime p through the polynomial line-sum profile
classification (15.588 Part I): a k-active element of Max+ is a tuple of
degree <= max(0, k-2) polynomial profiles mod p over the k active square
directions, with level-d coefficient vectors in the kernel of
sum_j c_j t_j^d = 0, plus row-sum, point-carry and +-p flip constraints.

The hard-coded paths point at the session scratchpad (`/tmp/e1work`); edit
`sys.path.insert` lines before reuse. Files:

| file | role |
|---|---|
| `kgen.py` | field context, square directions/coordinates, mod-p kernels |
| `kgen3.py` | `prep_subset`: per-subset kernels, u-lattice, class constant |
| `kgen4.py` | staged CPU enumerator (row filters, u-join, batched point tests) |
| `kgen5.py` | numba join `_prep_tables`/`_join_outer` + `process_outer` |
| `kgen6.py` | translation-gauged phases (top level lam != 0 => next level 0) |
| `flipsolve.py` | exact +-p flip assignment: propagation + pruned DFS (python) |
| `flipnb.py` | the same solver in numba (`flip_batch`), ~26us -> ~1us/cand |
| `gpu_inner.py` | V100 pipeline: numba generates packed candidate codes with probe prefilters; cupy tests batches (gathers, elementwise, matmul-with-ones reductions only — no NVRTC paths); flips resolved by `flipnb` |
| `dilation.py` | 120-element gauge group (square dilations x Frobenius): action on outer coefficient tuples, orbit/transversal machinery |
| `run_all.py` / `run_gauged.py` | multiprocess drivers |
| `k4_scan.py` | k=4 counts across primes (p=13: 168q with 28/35 subsets empty) |
| `k5_gauge_validate.py` | end-to-end check: dilation+translation gauge reproduces the ungauged k=5 total exactly |
| `spectrum.py` | Phi spectrum via svec(zz^T) in eigenspace coordinates (reproduces the exact p=5,7 spectra; D*lambda integral) |
| `moments.py` | 4-point moments / |kappa|=1 mu checks from an enumerated Max+ |

## Validation ledger (all exact)

- p=5: N=260; p=7: N=11452 — byte-identical to the cached enumerations.
- p=7 k=4 = 4410 through: CPU staged, numba join, GPU pipeline (3 paths).
- p=11 k=4 = 58080 (CPU and GPU paths agree; every vector passes Cy=py
  with residual exactly 0).
- p=11 k=5 = 1306800, twice: ungauged (run_all) and dilation+translation
  gauged (k5_gauge_validate) — identical totals.
- p=11 k=1 = 2772 = m*C(p,m); k=3 = 24200 = C(m,3)(p-1)q.
- p=13 k=4 = 168q, only 7 of 35 direction-quadruples nonempty.

## k=6 at p=11 (the remaining stratum)

Cost structure: ~94M probe-surviving candidates in the first top-stratum
outer alone; 13310 top-stratum outers, 161051 lower. The working plan is
dilation-gauged outers (~111 + ~1342 orbit reps) x the GPU inner
(~2-8 s/outer after `flipnb`), i.e. an hour-scale run, not a month-scale
one. `run_k6_gauged.py` (to be written from `k5_gauge_validate.py` +
`gpu_inner.py`) is the missing driver.
