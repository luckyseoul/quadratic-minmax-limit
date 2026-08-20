# Leftover 1 reduces to a variance bound plus a multiplicity fact

**Summary.** `lambda_min(Phi) >= 6` follows from `var(spec Phi) <= 32(n+10)^2/(n-6)^3`
together with `mult(lambda_min) >= n`. Both hold at p=7 and p=11; p=5 fails the
variance test and must stay a finite check (where `lambda_min = 80/13 >= 6` directly).
This is the "crude bound closes leftover 1 for p >= 7" that `fable.md` anticipates,
made explicit.

## The argument

Elementary. If `lambda_min` has multiplicity at least `m0`, then it contributes
`m0 * (mean - lambda_min)^2` to the centred second moment, so

```
    m0 * (mean - lambda_min)^2  <=  sum_i (lambda_i - mean)^2  =  d * var
    =>  mean - lambda_min  <=  std * sqrt(d / m0)
```

With `mean = 8(n-2)/(n-6)`, `d = n(n-6)/8`, and `m0 = n`, demanding
`mean - lambda_min <= mean - 6` rearranges to

```
    var  <=  32 (n+10)^2 / (n-6)^3                     (*)
```

`(*)` is sufficient for `lambda_min >= 6`, i.e. for leftover 1.

## Status at the three enumerable primes

| p | n | var (exact) | threshold `(*)` | slack | closes? |
|---|---|---|---|---|---|
| 5 | 26 | 8.725207 | 5.184000 | −68.3% | no |
| 7 | 50 | 0.909716 | 1.352367 | +32.7% | **yes** |
| 11 | 122 | 0.032134 | 0.357210 | **+91.0%** | **yes** |

`var * n` = 226.86, 45.49, 3.92 — empirically `var ~ n^-3.6`, while the threshold
decays only like `32/n`. The margin therefore widens rapidly; the crossover sits
between p=5 and p=7, exactly where `fable.md` places it.

## Why this is not just the standard second-moment bound

Without the multiplicity input, the best possible bound from `d`, `tr(Phi)`,
`tr(Phi^2)` alone (Cauchy-Schwarz on the remaining `d-1` eigenvalues) gives
`lambda_max(K) <= 22.03` at p=5 and `15.06` at p=7, against a requirement of `<= 2`
— and that gap *widens* with p (11.9x -> 15.7x). **So no bound using only the first
two spectral moments and the dimension can ever close leftover 1.** The multiplicity
fact is doing the essential work: it forbids a lone far outlier, because any outlier
drags at least `n` eigenvalues with it.

## Second moment is computable from leftover 3's tensor

Derived and verified exactly at p=5 and p=7 against the known spectra:

```
    tr(Phi^2)  =  4 ||M||_F^2  -  3 n^2  +  2 n^2 (n-1) / p^2
```

where `M[(ij),(kl)] = E_y[y_i y_j y_k y_l]` on `i<j`, `k<l` — leftover 3's own
four-point object. Route: `spec(Phi) = nonzero spec(Ghat/N)` gives
`tr(Phi^2) = E_{a,b}[(<y_a,y_b>^2 - 2n)^2]`, then split `sum_{ijkl} m4^2` by index
coincidences. (`E[<y_a,y_b>^2] = ||2P||_F^2 = 2n`, which is why `Ghat` is centred
at `-2n`.) This is the concrete form of fable.md's "leftovers 1 and 3 are moments
of one tensor", and it computes `tr(Phi^2)` with no eigensolver.

New exact values at p=11 (from the integer pair-moment Gram, `||G||_F^2 =
47738086747745638464`):

```
    tr(Phi^2)  = 2440162570133760 / 20130785689 = 121215.466094
    tr(Ghat^2) = 680278281952170147840          (integer, as required)
               = 2^16 * 3^3 * 5 * 11^2 * 61 * 139 * 181 * 414061
```

The float eigendecomposition gave 121215.46 — agreement to 8 significant figures
by a wholly separate path, so the two computations validate each other.

For comparison: `tr(Ghat^2)` = `2^12*3^2*5^2*13*37` (p=5),
`2^14*3^2*5^2*7^2*11*1399` (p=7). No p-formula is apparent yet.

## What remains

1. **Prove `mult(lambda_min) >= n` for p >= 7.** Observed n, n, 2n at p=5, 7, 11.
2. **Prove `var <= 32(n+10)^2/(n-6)^3` for p >= 7**, equivalently an upper bound on
   `tr(Phi^2)`, equivalently on `||M||_F^2` via the identity above.

Item 2 is the single blocking quantity and it is now a four-point moment estimate,
not a spectral one — the same tensor leftover 3 needs.

## Correction to fable.md

fable.md states "Bottom multiplicity is exactly `n`, top exactly `n/2`" (checked at
p=5,7 only). At p=11 the **bottom multiplicity is 244 = 2n, not n** — so the exact
claim fails at p=11. The weaker `>= n` form survives and is what the argument above
uses. Top multiplicity is 61 = n/2, consistent.

Full p=11 multiplicities: 244x4, 122x6, 61 (sum 1769 = dim Z = n(n-6)/8). All are
multiples of n/2 = 61 at all three primes.
