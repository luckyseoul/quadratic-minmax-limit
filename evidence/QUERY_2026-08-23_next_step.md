# Integration and next-step query (2026-08-23)

## Where the two roots actually stand

**R1 (leftovers 1 + 3, unified at 15.594).**  `Φ_part = λ̄I` is proved for
all p (15.597), so every open bit is `Φ_δ`, and leftover 1 ⟺
`Φ_δ ⪰ −(2n+20)/(n−6)·I` ⟺ `‖δ‖² ≤ n/12`.  Recast as a class function:
`Γ(g) = tr(Φ·π(g))`, `λ_c = ⟨Γ, χ_c⟩`, `Γ = λ̄ψ_Z + Γ_δ` with ψ_Z known.
Step 5 re-identified leftover 1 as 15.279 `Φ|_F : R̂_rest ≥ −2q²` — QVAR on
every even character, **binding minimum at k=8**, not QVAR's k=30.

**R2 (leftover 2) = W1 ∧ W2**, no generation gap once CLASS is exhaustive
(15.612).  W1 is being chipped away by covering congruences; W2 is an
existence question.

| W1 status | mechanism |
|---|---|
| p ≡ 5 (mod 8) | d=−1, ε=1 iff χ(2)=−1 (15.621) |
| p ≡ 17 (mod 24) | d=−2 four-point count (15.622) |
| p ≡ 73, 97 (mod 120) | d=−3 six-point count (15.623) |
| (2/p)₄ = −1 | d=−(p−1)/8, ε ≡ h(−8p)/4, Barrucand–Cohn (15.625) |
| **(2/p)₄ = +1 — OPEN** | linear box in (a,b,i,k) killed (15.627 A) |

## New this session: the interval family is exhausted too

Scan over the residual class `(2/p)₄=+1` (27 primes below 4000; ε as
defined in 15.625 `_eps_d`):

```
  p    p%16    (a,b)   (2/p)_8   eps(-(p-1)/8)   eps(-(p-1)/16)
  73     9    (3, 8)      --           0              n/a
  89     9    (5, 8)      --           0              n/a
 257     1   (1, 16)      +1           0               0
 337     1   (9, 16)      +1           0               1
 881     1  (25, 16)      +1           0               1
1217     1  (31, 16)      -1           0               0
1249     1  (15, 32)      -1           0               1
1553     1  (23, 32)      -1           0               0
```

* `ε(−(p−1)/4) = 0` on every residual prime — the quarter-interval is dead
  everywhere (confirms 15.625's setup).
* `ε(−(p−1)/8) = 0` on **every** residual prime — the eighth-interval is
  dead exactly off its (2/p)₄=−1 domain.  Clean confirmation of 15.625.
* **The halving pattern does NOT continue.**  `−(p−1)/16` exists only for
  p ≡ 1 (mod 16) — i.e. half the residual class — and is MIXED there
  (0,1,1,0,1,0), and mixed within each octic half too.  So it is neither a
  p-law nor octic-determined, consistent with 15.627 A's "MIXED on each
  octic half".
* No single named stay among {−1,−2,−3,−4,−a,−b,−(a+b),−(p±1)/4,8,16} is
  ε=1 across the residual class.

**Structural obstruction now named:** the residual class splits by
`b (mod 16)` — p ≡ 9 (mod 16) has b ≡ 8 (mod 16), p ≡ 1 (mod 16) has
b ≡ 0 (mod 16) — and the sixteenth-interval stay only *exists* on the
second half.  Any W1 completion must therefore use **two different stays**
on the two halves, or an invariant that is not an interval count.

## Query 1 (recommended, most tractable): close W2 by counting

W2 is an *existence* statement, and the data says the witnesses are a
growing positive fraction, not a near-miss:

| p | split-involution class | in U | W2 hits | hit rate |
|---|---|---|---|---|
| 17 | 153 | 49 | 17 | 11% |
| 31 | 496 | 146 | 76 | 15% |

The class is the conic `α²+βγ=1` (with ±), so `p(p+1)/2` points — a curve,
not a search space.  **Query: is the W2 condition cut out on that conic by
polynomial non-vanishing / coprimality conditions of bounded degree?**  If
yes, Weil gives

    #W2 = c·p² + O(p^{3/2}),   c > 0,

hence W2 ≠ ∅ for all p beyond an explicit bound, with p = 5, 11, 17, 31
already censused.  That closes W2 outright.  15.627 B already names the
missing piece as a *counting identity* ("some conjugate always coprime") —
this is exactly the shape a Weil bound proves, and the 11% → 15% growth is
evidence the main term is genuinely positive rather than marginal.

## Query 2 (W1 residual): two stays, or leave the interval family

Given the exhaustion above, the concrete question is whether the octic
layer of the tower supplies the next class-number identity, matching
15.625's `ε = N1+N4 ≡ N1−N4 = h(−8p)/4`:

* on p ≡ 1 (mod 16): is `ε(−(p−1)/16) ≡ h(−16p)/8 (mod 2)`, or tied to a
  Western/Reuschle octic-residuacity criterion in terms of `c` where
  p = a² + 64c²?
* on p ≡ 9 (mod 16): the sixteenth interval does not exist — what is the
  substitute stay?

**Warning for whoever takes this:** the residual class has only 27 primes
below 4000, and a criterion fitted on ~8 of them is exactly the
pre-asymptotic trap that killed the C=110 hypothesis on 2026-08-21.  Fit on
p < 2000, then PREDICT p ∈ (2000, 10⁴) before claiming anything.

## Not recommended

Re-deriving anything in the kills table of `README.md`, and in particular
another *linear* stay form on (a,b,i,k) — 15.627 A closed that box over 61
residual primes.
