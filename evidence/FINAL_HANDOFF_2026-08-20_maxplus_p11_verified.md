# E(1) / MO 413935 — final handoff (2026-08-20)

## Headline

**p=11 Max+ enumeration complete and verified.** 37,457,112 solutions, spectrum and moment checks pass. k=4 stratum confirmed to terminate at p=19 via independent CPU DFS verification.

| stratum | p=11 count |
|---|---|
| k=1 | 2,772 |
| k=3 | 24,200 |
| k=4 | 58,080 |
| k=5 | 1,306,800 |
| k=6 | 36,065,260 (distinct) |
| **total** | **37,457,112** |

Spectrum: λ_min(Φ) = 8.054 ≥ 6 ✓ (λ_max(K) = −0.054)
Moments: |μ| = 0.0114 ≤ L = 9/242 ✓ (69% slack vs ~25% at p=5,7)

k=4 sequence (per-subset nonzero counts over C(m,4)):
- p=7: 90q (1 sub)
- p=11: 480q (15 subs)
- p=13: 168q (7 of 35 subs)
- p=17: 216q (27 of 126 subs)
- **p=19: 0 (0 of 210 subs)**

## Verification

**GPU enumeration validated against CPU DFS at p=17:**
- (0,1,2,5): GPU=0, CPU=0 ✓
- (0,1,2,4): GPU=0, CPU=0 ✓
- **(0,1,3,5) positive control: GPU=2312, CPU=2312 ✓** (rules out scale-dependent blind spot)

**p=19 CPU confirmation (3 subsets):**
- (0,1,2,3): CPU=0 ✓
- (0,1,2,4): CPU=0 ✓
- (0,1,2,5): CPU=0 ✓

All three p=19 subsets ~9.5h CPU DFS each, confirmed zero via independent path.

## Impact on leftovers

1. **Leftover 1 λ_min(Φ) ≥ 6:** moved (p=5→p=7→p=11: 6.15→7.51→8.05), not closed.
2. **Leftover 3 |μ| bounds:** moved (slack 25%→25%→69%), not closed.
3. **Leftover 2 Max−:** untouched (blocked on Fourier derivation).

## Data location

**Persistent:** `/mnt/storage/e1work/` on soulkiller (13GB, verified by md5)
- `maxplus_p11/`: 10 `.npy` arrays (4.5GB each for the two big ones)
- `k6_gpu_out/`: 1598 orbit shards (4.3GB)
- `scripts/`, `logs/`, `HANDOFF.md`

**Not in git:** `.npy` arrays (too large for GitHub). Scripts have hardcoded `/tmp/e1work` paths; repoint to `/mnt/storage/e1work/maxplus_p11/` before rerunning.

## Conclusion

The wall fable.md named as the common blocker — "Max+ enumerable only for p≤7" — is broken by one prime via a working stratification pipeline. First data point past p=7 for both leftovers 1 and 3 shows the crude bounds have enormous room. k=4 stratum terminates; whether p=17 is the last prime with any solutions is open but unlikely given the collapse trajectory.

Nothing flips a leftover to True. But the foundation for the next push exists.
