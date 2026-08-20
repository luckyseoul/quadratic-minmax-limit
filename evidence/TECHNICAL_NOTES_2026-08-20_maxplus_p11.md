# E(1) / MO 413935 — detailed technical notes, session 2026-08-19/20

Working record of a session that enumerated `Max+` at p=11, computed its spectral
and four-point moment data exactly, and reduced leftover 1 to a single blocking
estimate. Nothing here flips a leftover flag to `True`.

Throughout: `p` prime, `q = p^2`, `n = q + 1`, `C` the Paley conference matrix of
order `n`, `P = (I + C/p)/2`, `Max+ = {y in {+-1}^n : Cy = py}`, `N = |Max+|`,
`Z = {B symmetric : CB = pB, diag B = 0}`, `d = dim Z = n(n-6)/8`,
`Phi(B) = E_y[(y^T B y)^2]`, `K = 8I - Phi`, `Ghat_ab = <y_a,y_b>^2 - 2n`.

`Max+` is closed under `y -> -y`. All four-point moments are even, so every
statistic below is computed on the `eps=+1` half `H_+` (where `y_inf = +1`),
of size `Nh = N/2`. Write `D = Nh/(2p) = N/(4p)`.

---

## 1. Enumeration of Max+ at p=11

### 1.1 Why brute force is out

`fable.md`: the nullity at p=11 is 61, so a sweep is `2^61 ~ 2.3e18`. It also
asserts "Max+ is enumerable only for `p <= 7`". That is the wall this section
removes.

### 1.2 Method — polynomial-profile stratification

Solutions are stratified by `k`, the number of *active* directions among the
`m = (p+1)/2` square-class directions. For each `k`-subset of directions the
solver enumerates candidates consistent with a degree-`<= k-2` polynomial profile,
then resolves the remaining sign freedoms. `k` ranges `1..m`; at p=11, `m = 6`.

Two independent implementations exist and are used against each other throughout:

- **GPU** (`k4_scan_gpu.py` + `gpu_inner.py`): batched candidate generation on the
  V100, CPU flip-assignment. Has a 5M-node per-candidate cap (`flip_batch`).
- **CPU** (`kgen3.enum_chunk`): plain backtracking DFS, **no node cap**.

### 1.3 Results

| stratum | count at p=11 |
|---|---|
| k=1 | 2,772 |
| k=2 | **0** |
| k=3 | 24,200 |
| k=4 | 58,080 |
| k=5 | 1,306,800 |
| k=6 | 36,065,260 (distinct) |
| **Nh** | **37,457,112** |

- `N = 2 * 37,457,112 = 74,914,224`
- `D = N/(4p) = 1,702,596` — exact integer
- `d = dim Z = 1769 = n(n-6)/8` — matches the formula exactly
- Eigen-equation residual `Cy - py` is identically zero over all 37.4M vectors
  (chunked int64 check; float64 would need a 34 GB cast, see §8.3)
- Strata verified pairwise disjoint

Artifact: `maxplus_p11_eps1.npy`, shape `(37457112, 122)`, `int8`, column 0 is
`y_inf`.

### 1.4 Validation — the load-bearing argument

The stratification pipeline **reproduces the independently known `Nh` exactly** at
both primes where `Max+` is enumerable by other means:

| p | per-k counts | total | expected `Nh` | match |
|---|---|---|---|---|
| 5 | {1:30, 2:0, 3:100} | 130 | 130 | yes |
| 7 | {1:140, 2:0, 3:1176, 4:4410} | 5726 | 5726 | yes |

`k=2` is empty at p=5, p=7, and — checked directly over all `C(6,2)=15` subsets —
at **p=11**. So `k in {1,3,4,5,6}` covers every `k <= m` and `assemble.py`'s
omission of `k=2` is correct rather than an oversight. (This was very nearly
shipped as an assumption; it is now a check.)

### 1.5 Symmetry closure — weaker than it looks

`Max+` is invariant under `Aut(C)`. Tested by sorting row-bytes once and using
`searchsorted` for membership, on 20,000 random solutions per group element:

- Dilation/Frobenius group (order 120, 12 sampled): 240,000 image tests, 0 failures
- Translation (`q=121` shifts, 13 sampled): 260,000 image tests, 0 failures

**These do not prove completeness.** Closure under a group cannot detect a missing
*whole orbit* — an entirely absent `k`-stratum passes closure trivially, because
the group permutes directions and therefore preserves `k`. They rule out a
*partially* enumerated orbit (the failure mode of the `k=6` gauge bug in §8.5) and
nothing more. The ground-truth reproduction in §1.4 is the real evidence.

This correction matters: closure was initially reported as "strong evidence of
completeness", which overstated it.

---

## 2. Spectrum of Phi at p=11

```
dim Z from constraints: 1769   = n(n-6)/8 = 1769
lambda_min = 8.054448680       >= 6 ?  TRUE
lambda_max = 8.664378396
lambda_max(K) = 8 - lambda_min = -0.054448680   <= 2 ?  TRUE
```

Full cluster decomposition (multiplicities sum to `244*4 + 122*6 + 61 = 1769 = d`):

| lambda | mult | D*lambda | dev. from integer |
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
| 8.664378396 | 61 | 14751936.00 | **1.86e-09** |

Trend across primes:

| p | lambda_min(Phi) | lambda_max(K) |
|---|---|---|
| 5 | 80/13 = 6.15385 | 48/26 = 1.84615 (binds) |
| 7 | 3072/409 = 7.51100 | 200/409 = 0.48900 |
| 11 | 8.05445 | **-0.05445** |

### 2.1 The `D*lambda` integrality question — resolved

Since `spec(Phi) = nonzero spec(Ghat/N)` with `Ghat` an **integer** matrix, and
`D = N/(4p)`:

```
    D * lambda  =  eig(Ghat) * D / N  =  eig(Ghat) / (4p)
```

So `D*lambda` is an integer **iff `4p` divides `eig(Ghat)`**. Checked against
fable.md's exact values: at p=5 the eigenvalues `{1600, 2880, 3520}` are all
divisible by `4p = 20`; at p=7 `{86016, ..., 120960}` are all divisible by `28`.
At p=11, `4p = 44` divides only 3 of the 11. **The p<=7 "all integral" pattern is
an arithmetic accident of divisibility that stops at p=11, not a defect.**

The three determined clusters give new exact `eig(Ghat)` values, and they are
*smooth*, matching the p=5,7 character:

| eig(Ghat) | mult | factorization |
|---|---|---|
| 631,428,864 | 122 | 2^8 * 3 * 11 * 41 * 1823 |
| 647,040,768 | 122 | 2^8 * 3 * 11 * 191 * 401 |
| 649,085,184 | 61 | 2^8 * 3^3 * 11 * 8537 |

**The other eight are undetermined, not non-integer.** `lambda` carries roughly
`6e-9..1e-7` accumulated float64 error, which is `+-0.5` to `+-7` in
`eig(Ghat) = lambda * N`. Naively rounding produces implausible values — one is a
nine-digit prime with multiplicity 244, two are odd — whereas a genuine
integer-Gram eigenvalue here should be smooth. That mismatch is the tell.

**Recorded so it is not repeated:** an attempt to argue the fractional parts
cluster on multiples of `1/44` was a **null result**. Max deviation was 0.0107
against a half-spacing of 0.01136, so the deviations fill 94% of the available
range and *any* value "fits". Eleven-for-eleven agreement there means nothing.

To settle the remaining eight: `Z` has a rational basis (`C` is integer), and the
exact four-point moment matrix is available as integers over `Nh` (§4). Rebuilding
`Phi` over that rational basis gives exact eigenvalues, hence exact `eig(Ghat)` —
which is what fable.md asks for ("seek a p-formula for `eig(Ghat)`, not for
`lambda`"). Estimated 1-2 hours.

---

## 3. Four-point moments and leftover 3 at p=11

```
pairing consistency max diff: 0.00e+00
max |mu| over |kappa|=1 four-sets: 0.011422344     (mu * Nh = 427848 exactly)
L   = (p-2)/2p^2      = 0.037190083    |mu| <= L  ?  TRUE
|T| = (p-2)/(p(2p-1)) = 0.038961039    |mu| < |T| ?  TRUE
max |m4| over ALL four-sets: 0.033641675
through-e |kappa|=1 sets: 5400, min G = -0.011422344 > -|T| ? TRUE
```

"Pairing consistency" checks that the three pairings of each 4-set give identical
moments. It is exactly zero over all 8,783,390 four-sets — a strong structural
check that the moment tensor is genuinely a four-point tensor of a common measure.

Arithmetic confirmation independent of completeness: fable.md states
"denominators throughout are `p*D = N/4`". Here `p*D = 11 * 1,702,596 = 18,728,556
= N/4`, and

```
    mu * (p*D) = 0.011422344 * 18,728,556 = 213,924    exactly an integer
```

so `mu = 213924/18728556` in lowest-terms-over-`p*D` form.

Margin trend — the slack **grows**, so leftover 3's crude bound gets easier:

| p | mu | L = (p-2)/2p^2 | mu/L |
|---|---|---|---|
| 5 | 3/65 | 3/50 | 0.769 |
| 7 | 109/2863 | 5/98 | 0.746 |
| 11 | 213924/18728556 | 9/242 | **0.307** |

`max|m4|` over all four-sets: `21/65` (p=5), `327/2863` (p=7), `0.033641675` (p=11).

---

## 4. tr(Phi^2) from the four-point tensor — derivation

This is the concrete form of fable.md's remark that "leftovers 1 and 3 are moments
of one tensor".

### 4.1 Setup

`spec(Phi) = nonzero spec(Ghat/N)` gives `tr(Phi^2) = tr(Ghat^2)/N^2`, and

```
    tr(Ghat^2)/N^2 = E_{a,b}[ (<y_a,y_b>^2 - 2n)^2 ]
                   = E[<y_a,y_b>^4] - 4n E[<y_a,y_b>^2] + 4n^2
```

The centring constant is forced:

```
    E_{a,b}[<y_a,y_b>^2] = sum_{i,j} E_a[y_i y_j] E_b[y_i y_j]
                         = sum_{i,j} (2P)_{ij}^2 = ||2P||_F^2 = 4 tr(P^2) = 4 tr(P) = 2n
```

using that `P` is a projection and `tr(P) = (n + tr(C)/p)/2 = n/2` since
`diag(C) = 0`. **This is exactly why `Ghat` is centred at `-2n`.** Hence

```
    tr(Phi^2) = E[<y_a,y_b>^4] - 4n^2,     E[<y_a,y_b>^4] = sum_{ijkl} m4(ijkl)^2
```

where `m4(ijkl) = E_y[y_i y_j y_k y_l]`.

### 4.2 Collapsing onto the pair-moment matrix

Split `sum_{ijkl} m4^2` by index coincidences, using `y_i^2 = 1`:

- `i=j, k=l`: `m4 = 1`. Count `n^2`. Contributes `n^2`.
- `i=j, k!=l`: `m4 = E[y_k y_l] = (2P)_{kl} = C_{kl}/p`. Since `C` has zero
  diagonal and `+-1` off-diagonal, `sum_{k!=l} C_{kl}^2 = n(n-1)`. Contributes
  `n * n(n-1)/p^2`.
- `i!=j, k=l`: same by symmetry.
- `i!=j, k!=l`: equals `4 * sum_{i<j, k<l} m4^2 = 4 ||M||_F^2`, where
  `M[(ij),(kl)] = E[y_i y_j y_k y_l]` over `i<j`, `k<l`.

Therefore

```
    tr(Phi^2)  =  4 ||M||_F^2  -  3 n^2  +  2 n^2 (n-1) / p^2                (I)
```

### 4.3 Verification

`(I)` reproduces the known exact spectra **exactly**:

| p | tr(Phi^2) from spectrum | from identity `(I)` | match |
|---|---|---|---|
| 5 | 85248/13 = 6557.538462 | 6557.538462 | yes |
| 7 | 3545625600/167281 = 21195.626521 | 21195.626521 | yes |

(`tr(Phi) = n(n-2)` also confirmed: 624 at p=5, 2400 at p=7.)

### 4.4 Exact values at p=11

Computed from the **integer** pair-moment Gram `G = Q^T Q` (`Q` the `Nh x 7381`
matrix of `y_i y_j`, entries `+-1`), so no eigensolver and no float rounding:

```
    ||G||_F^2  = 47,738,086,747,745,638,464
    tr(Phi^2)  = 2440162570133760 / 20130785689  = 121215.466094
    tr(Ghat^2) = 680,278,281,952,170,147,840          (integer, as required)
               = 2^16 * 3^3 * 5 * 11^2 * 61 * 139 * 181 * 414061
```

The float eigendecomposition of §2 gave `121215.46` — agreement to **8 significant
figures by a wholly independent path**. The two computations validate each other.

For p-formula hunting: `tr(Ghat^2)` = `2^12*3^2*5^2*13*37` (p=5),
`2^14*3^2*5^2*7^2*11*1399` (p=7), and the p=11 value above. No formula apparent.

### 4.5 Numerical technique

Entries of `Q` are `+-1`, so a per-chunk `Q_c^T Q_c` has integer entries bounded by
the chunk size. At 50,000 rows that is far below `2^24`, so each fp32 GEMM is
**exact**, and the V100 has no TF32 path to silently reduce precision. Chunks
accumulate into int64/float64.

`||G||_F^2` overflows int64 if summed naively: entries reach `Nh = 3.75e7`, squares
reach `1.4e15`, and a full row of 7381 such reaches `1.04e19` against an int64 max
of `9.22e18`. Rows are therefore summed in halves and combined as Python ints.

---

## 5. Leftover 1 — reduction to a variance bound plus a multiplicity fact

### 5.1 The sharper target

`fable.md` notes `lambda_max(K) = 48/n` *exactly* at p=5. Since
`8 - 48/n = 8(n-6)/n`, the natural general statement is

```
    lambda_min(Phi)  >=  8(n-6)/n
```

and `8(n-6)/n >= 6  <=>  n >= 24  <=>  p >= 5`. So this single inequality would
close leftover 1 for **every** p at once. Status:

| p | lambda_min(Phi) | 8(n-6)/n | holds |
|---|---|---|---|
| 5 | 6.15385 | 6.15385 | equality |
| 7 | 7.51100 | 7.04000 | yes |
| 11 | 8.05445 | 7.60656 | yes |

(`spectrum.py` already prints `8(n-6)/n` as "candidate"; p=11 confirms it.)

### 5.2 Pure second-moment bounds provably cannot work

Suppose only `d`, `tr(K)`, `tr(K^2)` are known. If `lambda_max(K) = m`, the other
`d-1` eigenvalues sum to `tr(K) - m` and their squares to `tr(K^2) - m^2`, so
Cauchy-Schwarz forces `(tr(K)-m)^2 <= (d-1)(tr(K^2)-m^2)`. The largest `m`
satisfying this is the best possible bound from that data:

| p | needed (`48/n`) | best 2-moment bound | gap |
|---|---|---|---|
| 5 | 1.8462 | **22.0308** | 11.9x |
| 7 | 0.9600 | **15.0608** | 15.7x |

Leftover 1 needs `lambda_max(K) <= 2`. The bound is an order of magnitude too weak
**and the gap widens with p**. So no argument using only the first two spectral
moments and the dimension can close leftover 1. Something structural is required.

### 5.3 The multiplicity input

Observed multiplicity of `lambda_min`: `n` (p=5), `n` (p=7), **`2n`** (p=11). All
`>= n`. More generally every multiplicity at all three primes is a multiple of
`n/2`:

- p=5: 26, 26, 13 (sum 65 = d)
- p=7: 50, 100, 50, 50, 25 (sum 275 = d)
- p=11: 244x4, 122x6, 61 (sum 1769 = d)

This is what defeats the worst case in §5.2: an outlier eigenvalue cannot be
*alone*. It drags at least `n` copies with it, and they all pay into the centred
second moment.

### 5.4 The bound

If `lambda_min` has multiplicity at least `m0`, then it contributes
`m0 (mean - lambda_min)^2` to the centred second moment, so

```
    m0 (mean - lambda_min)^2  <=  sum_i (lambda_i - mean)^2  =  d * var
    =>   mean - lambda_min  <=  std * sqrt(d / m0)                          (II)
```

With `mean = tr(Phi)/d = 8(n-2)/(n-6)`, `d = n(n-6)/8`, `m0 = n`, requiring
`mean - lambda_min <= mean - 6` rearranges to the clean sufficient condition

```
    var  <=  32 (n+10)^2 / (n-6)^3                                          (III)
```

**`(III)` + `mult(lambda_min) >= n`  =>  `lambda_min >= 6`  =>  leftover 1.**

### 5.5 Status of (III)

| p | n | var (exact) | threshold (III) | slack | closes? |
|---|---|---|---|---|---|
| 5 | 26 | 8.725207 | 5.184000 | **-68.3%** | no |
| 7 | 50 | 0.909716 | 1.352367 | +32.7% | **yes** |
| 11 | 122 | 0.032134 | 0.357210 | **+91.0%** | **yes** |

`var * n` = 226.86, 45.49, 3.92 — empirically `var ~ n^-3.6`, while the threshold
decays only like `32/n`. The margin therefore widens rapidly with `p`, and the
crossover sits between p=5 and p=7.

This is exactly the shape fable.md predicts: *"the floor binds only at the smallest
prime — p=5 is a finite check and by p=7 there is 4x slack. A crude bound closes
leftover 1 for p >= 7."* At p=5 the finite check is immediate: `lambda_min = 80/13
= 6.1538 >= 6`.

### 5.6 What remains for leftover 1

1. **Prove `mult(lambda_min) >= n` for `p >= 7`.** Observed `n, n, 2n`. No proof.
   Likely accessible from the `Aut(C)`-representation structure, given that every
   multiplicity is a multiple of `n/2`.
2. **Prove `(III)` for `p >= 7`** — equivalently an upper bound on `tr(Phi^2)`,
   equivalently on `||M||_F^2` via identity `(I)`.

Item 2 is the single blocking quantity, and `(I)` converts it from a spectral
problem into a **four-point moment estimate on the same tensor leftover 3 needs**.
Attacking that tensor now serves both leftovers.

---

## 6. The k=4 stratum terminates at p=19

| p | q | total | nonzero subsets | per-subset |
|---|---|---|---|---|
| 7 | 49 | 90q | 1/1 | 90q |
| 11 | 121 | 480q | 15/15 | 9 @ 40q, 6 @ 20q |
| 13 | 169 | 168q | 7/35 | 24q |
| 17 | 289 | 216q | 27/126 | 8q |
| **19** | **361** | **0** | **0/210** | — |
| 23+ | | untested | | `flip_batch` node-cap blowup |

### 6.1 Evidence

1. **Pipeline instrumentation.** p=19 generates 3.2M candidates and 385,136 raw
   solutions per subset; nothing dies silently. Each prime has a constant,
   subset-independent "degenerate background" that the activity filter strips
   (p=13: 19,032 = 14,976 background + 4,056 genuine).
2. **Scale probe.** The same code finds genuine solutions at p=17 (8 of the first
   40 subsets, each exactly `2312 = 8.00q`) and none at p=19 (0/40). Not a scale
   breakdown.
3. **Activity measurement.** At p=13, 80/300 sampled solutions are active in all
   four directions; at p=19, **0/300**, with per-direction activity itself falling
   from ~64% to ~25%.
4. **Independent CPU DFS**, no node cap, different algorithm — see below.

### 6.2 CPU cross-check (NUKA, `kgen3.enum_chunk`)

| p | subset | GPU | CPU | time | role |
|---|---|---|---|---|---|
| 13 | (0,1,2,3) | 4056 | 4056 | 267s | positive control |
| 13 | (0,1,2,4) | 0 | 0 | 263s | negative |
| 17 | (0,1,2,3) | 0 | 0 | 6310s | negative |
| 17 | (0,1,2,5) | 0 | 0 | 6546.7s | negative |
| 17 | (0,1,2,4) | 0 | 0 | 6600.6s | negative |
| **17** | **(0,1,3,5)** | **2312** | **2312** | 6423.4s | **positive control** |
| 19 | (0,1,2,3) | 0 | 0 | 34329.4s | the question |
| 19 | (0,1,2,4) | 0 | 0 | 34237.0s | the question |
| 19 | (0,1,2,5) | 0 | 0 | 34551.8s | the question |

The p=17 positive control was added mid-session after noticing that **every**
p=17 subset originally queued was one the GPU already reported as zero. Three
zeros agreeing with three zeros cannot distinguish "genuinely zero" from "both
paths blind at this scale" — and that capability is precisely what the p=19
conclusion rests on. `(0,1,3,5)` returning exactly 2312 closes that gap.

**Caveat on independence.** Both paths share `square_coords`/`prep_subset` for
context setup; only the search differs. So this catches algorithm-specific bugs
(e.g. the node cap) but not a shared setup bug.

### 6.3 Cost scaling — why the CPU runs took so long

The GPU cost grew `2.7x` from p=13 to p=17; the uncapped CPU DFS grew **23.8x**
(265s -> 6310s) while `q` grew only `1.71x`. Branching factor and depth both grow
with `p` and compound through the recursion. p=19 then took ~34,300s per subset,
another `5.2x`.

### 6.4 Consequence

The original objective "extract the p-law for k=4 counts" is void: the sequence
`90q -> 480q -> 168q -> 216q -> 0` terminates. The live question is *why* they
vanish, and whether p=17 is the last prime with any.

---

## 7. Corrections to fable.md

1. **"Bottom multiplicity is exactly `n`, top exactly `n/2`"** (stated from p=5,7).
   At p=11 the bottom multiplicity is **244 = 2n**, so the exact-`n` claim **fails
   at p=11**. Top is 61 = n/2, consistent. The weaker `>= n` form survives and is
   what §5 uses — anyone relying on exact `n` should know.

2. **"Max+ is enumerable only for `p <= 7`"** — no longer true. p=11 is enumerated
   (§1), not by brute force but by polynomial-profile stratification.

3. **The `D*lambda`-integral pattern at p=5,7 does not persist.** It is divisibility
   of `eig(Ghat)` by `4p` (§2.1), and it fails for 8 of 11 clusters at p=11.

---

## 8. Corrections to claims made during this session

Recorded so they are not chased again.

### 8.1 The `gpu_inner.py` "decode() inside recursion" bug was not a bug
Decoding a slice decodes exactly those candidates. The "fix" was a no-op, proven
by rerun: 2687s vs 2684s, byte-identical output. `_resolve_flips_recurse` is
harmless but its rewrite bought nothing.

### 8.2 p=19 was repeatedly called "corruption" without verification
It is a genuine zero (§6). Separately: p=19 *completed* with 0 and then **p=23**
crashed on the node cap — two distinct events that were conflated.

### 8.3 A `Cy == p*y` check that returned 0/300 valid was my own broken matrix
Built `q x q` instead of the order-`(q+1)` conference matrix with a point at
infinity. Caught only because p=13 was included as a control.

### 8.4 Symmetry closure was oversold as proof of completeness
See §1.5. It cannot detect a missing whole orbit.

### 8.5 `run_kgauged.py` gauge over-expansion (real bug)
Phase-`T` representatives were expanded by all `q=121` translations assuming the
action is **free**. It is not, for solutions with nontrivial translation
stabiliser. Raw 37,925,570 -> 36,065,260 distinct (~4.9% duplicates), fixed by
downstream dedup.

### 8.6 `assemble.py` saved the un-deduplicated array
It computed the distinct count but wrote the raw array. This is what made `D`
non-integral before the fix.

### 8.7 The "fractional parts cluster on multiples of 1/44" argument was a null result
See §2.1. Deviations filled 94% of the half-spacing; any value fits.

### 8.8 An ETA of "~53 min per p=17 CPU subset" was wrong by 2x
Actual ~1.8h. The underlying error was extrapolating CPU cost using the GPU's
scaling exponent (§6.3).

---

## 9. Reproduction

### 9.1 Data locations

**Canonical (persistent): `/mnt/storage/e1work/` on soulkiller** (9.1T drive).
`/tmp/e1work/` is **tmpfs (RAM-backed)** and does not survive reboot; contents
were copied out and verified (md5 match on both multi-GB arrays, 1598/1598 shards,
10/10 `.npy`).

```
/mnt/storage/e1work/
  maxplus_p11/    maxplus_p11_eps1.npy (4.5 GB)  <- the Max+ set
                  k6_p11_full.npy (4.6 GB)
                  G_pairmoment_p11.npy           <- integer pair-moment Gram
                  phiZ_p11.npy, m4diag_p11.npy, k4/k5 arrays
  k6_gpu_out/     1598 orb*.npy shards (4.3 GB)
  scripts/        all .py
  logs/           all .log
```

`.npy` arrays are **not in git** (4.5 GB each vs GitHub's 100 MB limit). Scripts
retain hardcoded `/tmp/e1work` paths — repoint to
`/mnt/storage/e1work/maxplus_p11/` before rerunning.

### 9.2 Key scripts

| script | purpose | cost |
|---|---|---|
| `assemble2.py` | build `maxplus_p11_eps1.npy` from strata, dedup | minutes |
| `closure_test.py` | symmetry closure checks (§1.5) | ~1 min |
| `spectrum.py` | `Phi` spectrum, `lambda_min` (§2) | ~31 min |
| `moments_gpu.py` | four-point moments, `mu` bounds (§3) | ~9 min |
| `trphi2_p11.py` | exact `tr(Phi^2)`, `tr(Ghat^2)` (§4.4) | ~10 min |

### 9.3 Hardware notes

- **soulkiller**: 88 cores, 60 GB RAM, Tesla V100-SXM2-16GB. All GPU work here.
- **NUKA** (192.168.1.192): 16 cores, **15 GB only** — cannot hold the 34 GB
  float64 spectrum load; use soulkiller for anything touching the full array.
  `numpy` is absent from system python; use `/home/nick/.venvs/rocm72/bin/python3`.
  `minmax_quadratic.py` had to be copied over manually.
- **CUDA 13 dropped sm_70**, so GPU JIT is broken for this V100: `numba.cuda`
  fails, and any CuPy path needing NVRTC (axis reductions, `cp.unique`, `.sum()`
  on bool, `RawKernel`) fails compiling `cuda_fp4.hpp`. Precompiled cuBLAS,
  cuSOLVER, elementwise ops and boolean masking work; substitute a matvec against
  a ones-vector for axis reductions.

### 9.4 Memory traps hit

- `moments.py` as written materialises `Q` at `37.4M x 7381` = **1.1 TB**. Must be
  chunked (`moments_gpu.py`).
- `spectrum.py` originally cast the full array to float64 = **34 GB**. Chunked.
- `assemble.py` used `set(map(tuple, ...))` over 39.3M rows and OOM'd. Replaced
  with numpy void-view dedup.

---

## 10. Open problems, in priority order

1. **Bound `||M||_F^2` (hence `tr(Phi^2)`, hence `var`) at general `p`.** Single
   blocking quantity for leftover 1 via `(III)`; same tensor leftover 3 needs.
   Identity `(I)` is the bridge.
2. **Prove `mult(lambda_min) >= n` for `p >= 7`.** Observed `n, n, 2n`. Every
   multiplicity is a multiple of `n/2` at all three primes — suggests an
   `Aut(C)`-representation argument.
3. **Exact `eig(Ghat)` at p=11 for the remaining 8 clusters** via a rational basis
   for `Z` plus the exact integer moment matrix. ~1-2 hrs. Feeds the p-formula
   hunt fable.md asks for.
4. **Leftover 2 (`multilevel_ND_k_ge_4p`)** — untouched. Max-; blocked on the
   Fourier derivation. The empirical class-swap is verified at p=5,7 (see
   `NOTE_maxminus_flat_marginal.md`) but the naive sign-flip derivation does **not**
   reproduce it; the direction<->frequency duality must be tracked exactly. Do not
   re-derive from the shortcut. The p=11 Max+ apparatus should transplant once the
   derivation lands.
5. **Why does k=4 vanish at p=19?** Structural explanation absent. Also unknown
   whether p=17 is the last prime with any k=4 solutions.
6. **p=23+ k=4** — `flip_batch` exceeds its 5M-node budget. Probably moot given (5).

---

## 11. Honest summary

Nothing flips a leftover flag to `True`. Two of three leftovers moved:

- **Leftover 1** now has an explicit route with a single named blocking estimate,
  and the crude bound fable.md anticipated is identified and verified at p=7,11.
- **Leftover 3** gained its first data point past p=7, with the slack *growing*
  from ~25% to 69%.
- **Leftover 2** is untouched.

The enumeration wall fable.md names as the common blocker is broken by one prime,
with a pipeline validated against ground truth rather than a one-off script. That
is the load-bearing change; everything in §2-§5 follows from having `Max+` at p=11.

By fable.md's own acceptance bar — a predicate returning `True` via a real import
from a unit that goes `False` under perturbation, never a census standing in for
general `p` — **this session did not close a leftover.** p=11 is still a census.
