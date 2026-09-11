# Structural Paley lift family `B=C+I`: closed form refuted, family pinned near 1/2

This extends `evidence/paley_structural_lift_20260911/`, which recorded three
exact orders of the family

\[
K=\begin{pmatrix}C&B\\B^T&-C\end{pmatrix},\qquad B=C+I,
\]

for `C` a Paley conference signing of order `n=q+1`, and described the
mechanism as "scalable". Two questions were left open there: whether the
family has a formula, and whether it is asymptotically useful. Both are
answered here, negatively.

## 1. A fourth exact order, and the closed form it suggests

The previous scorer (`scripts/paley_two_block_score.cpp`) hard-codes the
quadratic character mod a **prime** `q`, so it cannot reach `q=9`. The new
`scripts/paley_structural_lift_family.py` builds `C` over any prime power
`q=1 (mod 4)` via explicit `GF(q)` arithmetic, with `C=C^T`, zero diagonal
and `CC^T=qI` asserted at construction. Full enumeration gives

| q | lifted order N | exact Phi(K) | Phi/N^(3/2) |
| ---: | ---: | ---: | ---: |
| 5 | 12 | 18 | 0.4330 |
| 9 | 20 | 40 | 0.4472 |
| 13 | 28 | 70 | 0.4725 |
| 17 | 36 | 108 | 0.5000 |

`q=5,13,17` independently reproduce the earlier artifact's values through a
different (field-generic) construction path. `q=9` is new.

All four points satisfy `Phi = N(N+12)/16` exactly. This is a trap: the fit
was read off the same four points it explains, and it implies quadratic
growth `Phi ~ N^2/16`, hence a divergent normalized value.

## 2. The closed form is false

`N=36` is the largest order this family admits to exact enumeration
(`2^35` states); the next prime power is `q=25`, at order `N=52`. The
spectral bound of section 3 gives `Phi(K)<=203.07` there, while the fit
predicts `208`. The fit is therefore false at `q=25` -- the very first order
beyond the four it was read off -- and at every larger order, where its
quadratic growth runs away from the true `O(N^(3/2))` ceiling:

| q | N | fit `N(N+12)/16` | spectral bound | fit admissible |
| ---: | ---: | ---: | ---: | --- |
| 17 | 36 | 108 | 118.37 | yes (exact value 108) |
| 25 | 52 | 208 | 203.07 | **no** |
| 29 | 60 | 270 | 250.59 | **no** |
| 49 | 100 | 700 | 531.51 | **no** |

Independently, a witness search (`scripts/paley_structural_lift_search.py`,
multi-restart tabu on the exact integer form, every reported value replayed
as an explicit `+-1` witness, validated by recovering all four exact optima
above) reaches only `182` at `q=25`. Independent deep shards agree on a
stable value at each of the three orders retested, in every case far under
the fit:

| q | N | fit | deep shards (seeds) | stable value |
| ---: | ---: | ---: | ---: | ---: |
| 25 | 52 | 208 | 4 x (3000 restarts, 30000 sweeps) | 182 |
| 29 | 60 | 270 | 6 x (1200 restarts, 20000 sweeps) | 222 |
| 53 | 108 | 810 | 6 x (1200 restarts, 20000 sweeps) | 534 |

These are still lower bounds, not exact norms; their agreement across seeds
is stability evidence only.

The fit agreeing on exactly the four enumerable orders and failing at the
first order past them is the characteristic shape of a pre-asymptotic
coincidence, not evidence for a family formula.

## 3. Exact spectral bound for the family

The lift decomposes as `K=M+J` with `M=[[C,C],[C,-C]]` and `J=[[0,I],[I,0]]`,
and direct multiplication gives the identity

\[
K^2=(2q+1)I+2\,(C\oplus C),
\]

verified symbolically and numerically. Since `C` has spectrum `+-sqrt(q)`,

\[
\|K\|=\sqrt{2q+1+2\sqrt q},\qquad
\Phi(K)\le \tfrac12 N\|K\| = (q+1)\sqrt{2q+1+2\sqrt q}.
\]

So `Phi(K)/N^(3/2) <= sqrt(2q+1+2 sqrt q)/(2 sqrt(2q+2))`, which exceeds
`1/2` for every `q` and decreases to `1/2`. Two consequences:

- The family's growth is `O(N^(3/2))`, not `N^2`: the section-1 closed form
  is false for all large `q`, and the exact agreement at `N<=36` is a
  small-order coincidence.
- The spectral bound by itself never certifies a constant below `1/2`, so it
  cannot improve the existing `limsup alpha_n <= 1/2`.

## 4. Where the family actually sits

Witness lower bounds were computed on the four-node local mesh (soulkiller
88 threads, jellyfin 16, NUKA 16, orin 6) for all prime powers `q=1 (mod 4)`
up to `q=197`, i.e. orders up to `N=396`. Normalized values are flat:

| q | N | witness | witness/N^(3/2) | spectral bound / N^(3/2) |
| ---: | ---: | ---: | ---: | ---: |
| 61 | 124 | 654 | 0.4736 | 0.5287 |
| 97 | 196 | 1294 | 0.4716 | 0.5233 |
| 137 | 276 | 2156 | 0.4702 | 0.5199 |
| 169 | 340 | 2938 | 0.4686 | 0.5181 |
| 197 | 396 | 3708 | 0.4705 | 0.5168 |

The complete per-order table is in `raw/`. Every tested order satisfies
`0.468 <= witness/N^(3/2)`, with no downward trend over a factor of eleven in
order, while the exact small-order values *rise* from `0.433` to `0.500`.

Because `m_N <= Phi(K)` is the only direction this construction supplies, and
`Phi(K)` is at least the witness value, the family cannot yield an upper
bound constant below about `0.47` at any tested order. It is therefore not a
route to a value below `1/2` by more than a few percent, and not a route to
anything near `1/pi`.

## Interpretation boundary

Section 3 is an exact all-orders bound for this one family. Section 4's
values are rigorous lower bounds on `Phi(K)` obtained by heuristic search:
they are not exact norms, and a larger true `Phi` only strengthens the
negative conclusion. Nothing here bounds `m_N` from below, decides the
original convergence question, or applies to cross blocks other than
`B=C+I`. The failure of one algebraic family is not evidence about the
value of the limit.
