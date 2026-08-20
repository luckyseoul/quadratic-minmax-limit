# E(1) / MO 413935 — session handoff (2026-08-19)

## 1. Headline results

### p=11 Max+ eigenspace fully enumerated and verified

| stratum | count |
|---|---|
| k=1 | 2,772 |
| k=3 | 24,200 |
| k=4 | 58,080 |
| k=5 | 1,306,800 |
| k=6 | 36,065,260 (distinct) |
| **total (eps=+1 half)** | **37,457,112** |

- `N = 2 x 37,457,112 = 74,914,224`
- `D = N/(4p) = 1,702,596` — **exact integer** (was non-integer before dedup; see §3)
- **Eigen-equation residual = 0** over all 37.4M vectors (`Cy = py`, chunked int64 check)
- Strata are pairwise **disjoint** (verified: zero cross-stratum overlap)

Artifact: `/tmp/e1work/maxplus_p11_eps1.npy` (37,457,112 x 122 int8, col 0 = y_inf)

### Completeness: validated against ground truth

The k-stratum pipeline (k=1 family + `prep_subset`/`enum_chunk` for k>=2)
**reproduces the known Nh exactly** at both primes where Max+ is independently known:

| p | per-k | total | expected Nh | match |
|---|---|---|---|---|
| 5 | {1:30, 2:0, 3:100} | 130 | 130 | YES |
| 7 | {1:140, 2:0, 3:1176, 4:4410} | 5726 | 5726 | YES |

`k=2` is **empty** at p=5, p=7, and (verified directly, 15 subsets) at **p=11** —
so omitting it from `assemble.py` is correct, and k=1,3,4,5,6 covers all k <= m=6.

This is the load-bearing completeness argument. **The symmetry tests below are
weaker than I first claimed**: closure under the automorphism group cannot detect
a *missing whole orbit* (an absent k-stratum would pass closure trivially). They
are a necessary check, not a sufficient one.

### Completeness corroborated by two independent symmetry tests

Max+ must be invariant under Aut(C). Both tested by sorted-key membership on
20,000 random solutions per group element:

- **Dilation/Frobenius group** (120 elements, sampled 12): 240,000 image tests, **0 failures**
- **Translation** (q=121 shifts, sampled 13): 260,000 image tests, **0 failures**

These rule out a *partially enumerated orbit* (the k=6 gauge bug's failure mode),
but NOT a missing whole stratum. Treat as corroboration of the ground-truth
validation above, not as an independent completeness proof.

### Spectrum check PASSES

```
dimZ from constraints: 1769   = n(n-6)/8 = 1769   [matches exactly]
lambda_min = 8.054448680      >= 6 ?  TRUE
lambda_max = 8.664378396
```

Full cluster decomposition (multiplicities sum to 244*4 + 122*6 + 61 = 1769 = dimZ):

| lambda | mult | D*lambda | err |
|---|---|---|---|
| 8.054448680 | 244 | 13713472.10 | 1.05e-01 |
| 8.108397112 | 244 | 13805324.49 | 4.89e-01 |
| 8.169097011 | 244 | 13908671.90 | 1.05e-01 |
| 8.219347110 | 244 | 13994227.51 | 4.89e-01 |
| 8.301818329 | 122 | 14134642.68 | 3.20e-01 |
| 8.335108219 | 122 | 14191321.91 | 8.69e-02 |
| 8.410919161 | 122 | 14320397.32 | 3.20e-01 |
| 8.428691246 | 122 | 14350656.00 | **5.59e-09** |
| 8.451605717 | 122 | 14389670.09 | 8.69e-02 |
| 8.637088305 | 122 | 14705472.00 | **7.45e-09** |
| 8.664378396 |  61 | 14751936.00 | **1.86e-09** |

**The `D*lambda` question — analysed, partly resolved.**

Exact identity (from fable.md's `spec(Phi) = nonzero spec(Ghat/N)` with `Ghat`
an **integer** matrix, and `D = N/(4p)`):

```
    D * lambda  =  eig(Ghat) / (4p)        so   D*lambda integral  <=>  4p | eig(Ghat)
```

Verified against fable.md's exact values: p=5 eigs {1600,2880,3520} all divisible
by 20; p=7 {86016,...,120960} all divisible by 28. At p=11 (4p=44) only 3 of 11
clusters are divisible. **So the p<=7 all-integer pattern is an arithmetic accident
of divisibility that stops at p=11 — not a defect in the enumeration.**

Three clusters are determined (D*lambda exact to ~1e-9), and their `eig(Ghat)` are
**smooth**, matching the p=5,7 character:

| eig(Ghat) | mult | factorization |
|---|---|---|
| 631,428,864 | 122 | 2^8 * 3 * 11 * 41 * 1823 |
| 647,040,768 | 122 | 2^8 * 3 * 11 * 191 * 401 |
| 649,085,184 |  61 | 2^8 * 3^3 * 11 * 8537 |

**The other eight are UNDETERMINED, not non-integer.** lambda carries ~6e-9..1e-7
error (float64 accumulation over 37.4M rows), which is +/-0.5 to +/-7 in
`eig(Ghat) = lambda*N`. Rounding them yields implausible values (one is a 9-digit
prime with multiplicity 244; two are odd) — the signature of noise, since a genuine
integer-Gram eigenvalue at this p should be smooth.

*Recorded so it is not repeated:* an attempt to argue the fractional parts cluster
on multiples of 1/44 was a **null result** — max deviation 0.0107 vs half-spacing
0.01136, so the deviations fill 94% of the available range and any value "fits".

**To settle:** `Z = {B sym : CB = pB, diag B = 0}` has a rational basis (C is
integer), and the exact 4-point moment matrix is already available as integers/Nh
from `moments_gpu.py`. Rebuilding Phi over that rational basis gives exact
eigenvalues, hence exact `eig(Ghat)` at p=11 — which is what fable.md asks for
("seek a p-formula for `eig(Ghat)`, not for `lambda`"). Est. 1-2 hrs.

### k=4 stratum appears to DIE at p=19

| p | q | total | nonzero subsets | per-subset |
|---|---|---|---|---|
| 7 | 49 | 90q | 1/1 | 90q |
| 11 | 121 | 480q | 15/15 | 9 @ 40q, 6 @ 20q |
| 13 | 169 | 168q | 7/35 | 24q |
| 17 | 289 | 216q | 27/126 | 8q |
| **19** | **361** | **0** | **0/210** | — |
| 23+ | | untested | | flip node cap blowup |

Evidence that p=19's zero is **real mathematics, not a bug**:

1. **Pipeline instrumentation** — p=19 generates 3.2M candidates and 385,136 raw
   solutions per subset; nothing dies silently. Each prime has a constant
   subset-independent "degenerate background" that the activity filter strips
   (p=13: 19,032 = 14,976 background + 4,056 genuine).
2. **Scale probe** — same code finds genuine solutions at p=17 (q=289): 8 of first
   40 subsets, each exactly 2,312 = 8.00q. Finds none at p=19 (q=361), 0/40. So it
   is not a scale breakdown.
3. **Direct activity measurement** — at p=13, 80/300 sampled solutions are active in
   all 4 directions; at p=19, **0/300**, with per-direction activity itself dropping
   from ~64% to ~25%.
4. **Independent CPU path** (`kgen3` DFS, no node cap, different implementation)
   reproduces the GPU exactly at p=13 on both a nonzero case (4,056 = 24.000q) and a
   zero case (0), and at p=17 sub (0,1,2,3) → 0.

**Status: p=19 CPU confirmation still running** on NUKA (`p19_verify_nuka_mp.py`,
5-way parallel over p=17 subs (0,1,2,4),(0,1,2,5) and p=19 subs (0,1,2,3),(0,1,2,4),(0,1,2,5)).
This is the last thing that could overturn the conclusion.

**Implication:** the original objective "extract the p-law for k=4 counts" is likely
void. The sequence 90 → 40/20 → 24 → 8 → 0 terminates. The live question becomes
*why* they vanish and whether p=17 is the last prime with any.

### Moments / mu-bound check PASSES

```
pairing consistency max diff: 0.00e+00        [all 3 pairings agree exactly]
max |mu| over |kappa|=1 four-sets: 0.011422344   (mu*Nh = 427848 exactly)
L = (p-2)/2p^2      = 0.037190083   |mu| <= L ?  TRUE
|T| = (p-2)/(p(2p-1)) = 0.038961039  |mu| < |T| ? TRUE
max |m4| over ALL four-sets: 0.033641675
through-e |kappa|=1 sets: 5400, min G = -0.011422344 > -|T| ? TRUE
```

`mu * (p*D) = 0.011422344 * 18,728,556 = 213,924` **exactly** — an integer
numerator over the denominator `p*D = N/4` that fable.md predicts. This is an
arithmetic check independent of whether the enumeration is complete.

## 1b. Impact on the three fable.md leftovers

**Nothing flips a flag to True. No route is killed.** But:

### Leftover 1 `phi_F_ge_6_proved_general` (lambda_min(Phi) >= 6) — MOVED

| p | lambda_min(Phi) | lambda_max(K) = 8 - lambda_min |
|---|---|---|
| 5 | 80/13 = 6.1538 | 48/26 = 1.846 (binds) |
| 7 | 3072/409 = 7.5110 | 200/409 = 0.489 |
| **11** | **8.0544** | **-0.054** |

fable.md: "leftover 1 is exactly `lambda_max(K) <= 2`". At p=11 it is **negative**.
First data point beyond p=7; fits the increasing trend.

**The wall named in fable.md is broken by one prime.** fable.md asserts "Max+ is
enumerable only for `p <= 7`" and "at p=11 the nullity is 61, so an exhaustive
sweep is `2^61 ~ 2.3e18`". Max+ at p=11 is now fully enumerated via the
polynomial-profile stratification (not brute force), validated against ground
truth at p=5,7. Since fable.md also says "getting Max+ moments at general `p` is
plausibly the single underlying problem", this is the load-bearing piece.

Does NOT close it: acceptance forbids "a census at p <= 7 standing in for general
`p`", and p=11 is still a census. fable.md's own route is "a crude bound closes
leftover 1 for `p >= 7`" — this shows that bound has enormous room.

### Leftover 3 `type_I_aut_e_3AB_positive_general` (3A+B > 0) — MOVED, route easier

| p | mu | L = (p-2)/2p^2 | mu/L |
|---|---|---|---|
| 5 | 3/65 | 3/50 | 0.769 |
| 7 | 109/2863 | 5/98 | 0.746 |
| **11** | **213924/18728556** | **9/242** | **0.307** |

fable.md notes ~25% slack at p=5,7 and "a crude bound suffices here too". At p=11
the slack is **69%** — the crude bound gets easier, not harder. Also
`max|m4|` over all four-sets: 21/65 (p=5), 327/2863 (p=7), **0.033641675** (p=11).

### Leftover 2 `multilevel_ND_k_ge_4p_proved` — UNTOUCHED

Max−; blocked on the Fourier derivation (see §2). The identified route is
transplanting the Part I polynomial-profile machinery to Max−; the p=11 Max+
apparatus would carry over once that derivation lands.

## 2. Still open

- p=19 independent CPU confirmation (running, NUKA)
- ~~Moments / mu-bound check~~ **DONE, passes** (see §1)
- Max− polynomial-profile analog — **blocked**. Empirical class-swap verified at
  p=5,7 (see `NOTE_maxminus_flat_marginal.md`) but the naive Fourier sign-flip
  derivation does NOT reproduce it. Needs the direction<->frequency duality tracked
  exactly. Do not re-derive from the shortcut; it gives the wrong answer.
- p=23+ k=4: `flip_batch` blows its 5M-node per-candidate budget. Probably moot if
  k=4 is genuinely empty past p=17, but unconfirmed.
- p-adic structure on deep runs — **shelved by user** (needs primes to 500M x bases
  to 500B; days of compute, may be unrelated).

## 3. Bugs found and fixed this session

| file | bug | fix |
|---|---|---|
| `assemble.py` | OOM: `set(map(tuple,...))` over 39.3M rows | numpy void-view dedup |
| `assemble.py` | **saved the un-deduplicated array** despite computing the distinct count | rebuild with dedup applied |
| `spectrum.py` | 34 GB float64 cast of the full array | chunked projection (`assemble2.py` pattern) |
| `moments.py` | materializes `Q` = 37.4M x 7381 = **1.1 TB** | `moments_gpu.py`: chunked `Q^T Q` on V100 |
| `run_kgauged.py` | k=6 gauge over-expansion: phase-T reps expanded by all q=121 translations assuming the action is **free**; it isn't for solutions with nontrivial translation stabilizer | dedup downstream (37,925,570 raw → 36,065,260 distinct, ~4.9%) |

**Precision note on `moments_gpu.py`:** entries of `Q` are ±1, so a per-chunk
`Q^T Q` has integer entries bounded by the chunk size (50k << 2^24) — each fp32
GEMM is *exact*, and V100 has no TF32 path to silently reduce precision. Chunks
accumulate into float64.

## 4. Corrections to earlier claims in this session

Recorded so the next person doesn't chase them:

- **`gpu_inner.py` "decode() inside recursion" was NOT a bug.** Decoding a slice
  decodes exactly those candidates. The "fix" was a no-op — proven by rerun:
  2687s vs 2684s, byte-identical result. `_resolve_flips_recurse` is harmless but
  its rewrite bought nothing.
- **p=19 was repeatedly called "corruption" without verification.** It is almost
  certainly a real zero (§1).
- **p=19's zero and the crash are two separate events**: p=19 *completed* with 0,
  then **p=23** crashed on the flip node cap.
- A `Cy == p*y` check that returned 0/300 "valid" was **my broken matrix
  construction** (q x q instead of the order-(q+1) conference matrix with a point
  at infinity), caught only because p=13 was included as a control.

## 5. Where things live

**CANONICAL (persistent): `/mnt/storage/e1work/` on soulkiller** — 9.1T drive, 13 GB used.
Data was originally in `/tmp/e1work/`, which is **tmpfs (RAM-backed)** and does not
survive reboot; it was copied out and verified (md5 match on both multi-GB arrays,
1598/1598 shards, 10/10 `.npy`).

```
/mnt/storage/e1work/
  maxplus_p11/     10 .npy incl. maxplus_p11_eps1.npy (4.5 GB), k6_p11_full.npy (4.6 GB)
  k6_gpu_out/      1598 orb*.npy shards (4.3 GB)
  scripts/         all .py
  logs/            all .log
  HANDOFF.md
```

The `.npy` arrays are **NOT in git** (4.5 GB each, far over GitHub's 100 MB limit).
Scripts in the repo still contain hardcoded `/tmp/e1work` paths — repoint them to
`/mnt/storage/e1work/maxplus_p11/` before rerunning.
- NUKA (192.168.1.192): `/tmp/e1work_nuka/`, venv `/home/nick/.venvs/rocm72/bin/python3`
- Repo: `~/quadratic-minmax-limit` (per user: commit **results/data only**, not scratch code)

### Machine notes
- **soulkiller**: 88 cores, 60 GB, Tesla V100-SXM2-16GB. This is where the session runs.
- **NUKA**: 16 cores, **15 GB only** — cannot hold the 34 GB float64 spectrum load; use
  soulkiller for anything loading the full Max+ array. numpy is NOT in system python
  there; use the rocm72 venv. `minmax_quadratic.py` had to be copied over manually.
