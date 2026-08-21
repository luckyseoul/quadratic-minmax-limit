# p=13 / p=17 Max± census plan — for the |μ| ≤ 2/n attack (leftover 3)

Written 2026-08-21 as a handoff. Goal: extend the |κ|=1 four-point moment
census (`evidence/mu_kappa_phi_census.json`, currently p=5,7,11) to p=13 and
if feasible p=17, to (a) re-verify `|μ| ≤ 2/n` at new primes and (b) feed the
formula hunt for a general-p proof (15.191 K → leftover 3 → residual-(i)).

**What μ needs:** μ(S) = ½(E₊+E₋)[∏_{i∈S} y_i] over the FULL Max+ and Max−
sets — every stratum k=1..m, m=(p+1)/2. NOT just k=4. At p=13 that is
k ∈ {1,3,4,5,6,7} (k=2 empty at p=5,7,11; verify before assuming at p=13).

## Known counts / anchors

| p | m | full |Max+| (eps=+1 half) | status |
|---|---|---|---|
| 5 | 3 | 130 | done, ground truth |
| 7 | 4 | 5,726 (k: 140/1176/4410) | done, ground truth |
| 11 | 6 | 37,457,112 (k6 = 36,065,260) | done, verified Cy=py |
| 13 | 7 | unknown; k=4 = 168q = 28,392 (28/35 direction-quads EMPTY) | k=7 PSL-orbit data exists |
| 17 | 9 | unknown; k=4 = 216q (27/126 quads) | nothing else |

Growth is ~20-25x per k-stratum (p=11: 58k → 1.3M → 36M for k=4,5,6).
Extrapolation: p=13 full set plausibly 10⁸–10⁹ vectors (top stratum k=7
dominates), i.e. 20–200 GB int8 and days-to-weeks of V100 outers if done by
raw enumeration. p=17 (m=9) raw enumeration is likely OUT OF REACH.

## Route A (recommended): orbit-representative moment sums — no full enumeration

Max+ is invariant under the signed automorphism group G ⊇ signed-PSL(2,q)
(15.267/15.268 for the exact signed action; |PSL(2,169)| ≈ 2.4M). If the
stratum decomposes into orbits Orb(y_r) with known representatives y_r and
orbit sizes, then for any four-set S:

    Σ_{y ∈ Orb(y_r)} y_S = (|Orb(y_r)|/|G|) · Σ_{g∈G} ε(g,·) (g·y_r)_S

so the full moment sum over Max+ is a sum over (orbit reps × group elements),
never materializing the stratum. Cost per orbit: |G| ≈ 2.4M signed-permutation
applications restricted to 4 coordinates — GPU-trivial (gather + product).
Equivalently: for each rep, accumulate the rank-4 tensor contribution by
pushing S through g⁻¹; batch over sampled four-set class representatives
(one per (κ,φ,cross-ratio) class, ~O(q) classes) instead of all C(n,4) sets.

Existing assets for p=13 k=7: `/mnt/storage/e1work/maxplus_p13/`
(`k7_p13_signed_psl_orbit.npy` 111 MB, per-orbit seeds/jsons in
`orbit_attack_2026-08-20/`, SHA256SUMS). Confirm whether the k=7 orbit list
is COMPLETE (the union-completeness logs on NUKA:
`/tmp/e1work_nuka/k7_p13_union*_completeness_seed*.log` were the check).
Lower strata (k=1..6) at p=13 are small enough for Route B below; k=7 via
orbits. p=17: Route A only, and only if someone first produces the top-strata
orbit reps (CP-SAT seed + orbit closure, as in `k7_p13_cpsat_probe.py`).

## Route B: raw stratum enumeration (toolkit facts, corrected API)

Toolkit: `scripts/maxplus_profile_enum/` (committed at 5f0ac34). Key APIs —
note "subset" means a 4/5/6-tuple of DIRECTION indices, 0..m-1, not field
points (my first draft got this wrong; C((p+1)/2, k) subsets):

```python
import sys, itertools
sys.path.insert(0, "scripts/maxplus_profile_enum")   # also fix the hardcoded
        # sys.path.insert(0,'/tmp/e1work') inside the modules, or symlink.
from kgen import square_coords
from kgen3 import prep_subset
from kgen4 import enum_chunk
from multiprocessing import Pool

p = 13
dirs, forms, coords = square_coords(p)          # m = len(dirs) = (p+1)/2
m = len(dirs)
tasks = []
for sub in itertools.combinations(range(m), k):  # k = stratum
    ctx = prep_subset(p, list(sub), forms, coords)
    tasks.append((ctx, 0, ctx["outer_total"], 1))  # eps=+1
with Pool(workers) as pool:
    for sub, lo, hi, sols in pool.imap_unordered(enum_chunk, tasks):
        ...  # sols: list of int8 y-vectors over the q finite points
```

For k ≥ 5 use the dilation/translation gauge driver (`run_kgauged.py`,
currently p=11-hardcoded — parameterize p) and remember the k=6 lesson:
**the translation action on phase-T reps is NOT free**; expand then dedup
(numpy void-view), never assume ×q multiplicity. GPU inner loop:
`gpu_inner.py` — V100 is sm_70 and CUDA 13 dropped it, so NO NVRTC paths
(no cp.unique / axis reductions / bool .sum(); elementwise + gather +
matmul only). Chunk fp32 GEMMs ≤ 50k rows so int entries stay exact (< 2^24).

Time estimates (scale from p=11 logs in /mnt/storage/e1work/logs/):
k=4: minutes (35 subsets). k=5: hours (21 subsets × bigger outers).
k=6: ~1-2 days on V100. k=7: weeks raw — use Route A instead.

## Validation protocol (non-negotiable, from the p=11 campaign)

1. Reproduce ground truth first: the same pipeline at p=5 must give 130 and
   at p=7 must give 5726 (per-k {1:140, 3:1176, 4:4410}). Set equality, not
   totals, wherever an independent path exists.
2. Exact eigencheck on everything kept: Cy = py, residual 0, chunked int64.
3. Strata pairwise disjoint; k=2 emptiness verified, not assumed.
4. Closure spot-tests under dilation/Frobenius and translation (catches
   partially-enumerated orbits; CANNOT catch a missing whole stratum — only
   ground-truth reproduction covers that).
5. μ integrality anchor: max|μ| must clear denominators to an integer over
   p·D = N/4 (fable.md); at p=11 this was 213,924 exactly.
6. Estimate total runtime out loud before launching multi-hour jobs; batch
   GPU candidate tests from the start (per-candidate Python paths became the
   bottleneck three separate times at p=11).

## Machines

- soulkiller: 88 cores, 60 GB RAM, V100-16GB (sm_70; see CUDA caveat above).
- NUKA (192.168.1.192): 16 cores, 15 GB RAM ONLY — fine for DFS verification
  (`kgen3` CPU path), cannot hold large arrays; numpy lives in
  `/home/nick/.venvs/rocm72/bin/python3`.
- Store results under `/mnt/storage/e1work/maxplus_p13/` (canonical, 9.1T);
  /tmp is tmpfs and dies on reboot. Big .npy stay OUT of git.

## What this feeds

`evidence/NOTE_2026-08-21_remaining_general_p_estimates.md` leftover 3:
sufficient |μ| ≤ L, target |μ| ≤ 2/n (15.191 K, strictly stronger for p≥7).
Census slack so far: 0.769 (p=5) → 0.746 (p=7, TIGHT for 2/n: 0.952) →
0.307 (p=11). A p=13 point decides whether the 2/n ratio keeps falling —
and any new prime where it fails kills the 2/n route immediately, which is
also progress. No census flips the flag by itself (fable.md acceptance bar).
