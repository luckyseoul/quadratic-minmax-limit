# W2: the coprimality model — fingerprint-confirmed in both directions

Date: 2026-08-23.  Target: W2 (15.612), the odd-factor half of Walsh, which
with W1 closes leftover 2's 15.406 E.  **No flag flipped; W2 stays open.**

## Reframing (structural, from the PIR)

R = F2[X]/h is a principal ideal ring, so I_U = (γ) for a unique monic
γ | h, Aut-invariant: γ = (X+1)^i·∏_O f_O^{e_O}.  Walsh ⟺ γ=1;
W1 ⟺ i=0; W2 ⟺ all e_O=0; deg γ = codim I_U.  Hence **one U-difference
with content coprime to g settles every orbit at once** — the per-orbit
structure of W2 is a red herring for the search; the object of study is
the DENSITY of unit-content differences.  (`_switched`'s W2 test already
checks exactly this: gcd(content, f)=1 for all f | g.)

## The model and the census

Random-element heuristic: Pr[content coprime to g] = ∏_{f|g}(1−2^{−deg f}).
Census of the switched split-involution class (rates conditioned on in-U —
the class-level 11%/15% figures conflate three filters):

| p | m | model | measured in-U rate | direction predicted/observed |
|---|---|---|---|---|
| 17 | 9 | 0.738 | 17/49 = 0.347 | — |
| 31 | 15 | 0.618 | 76/146 = 0.521 | ↑ / ↑ |
| 41 | 105 | 0.458 | 63/237 = **0.266** | **↓ / ↓** |
| 47 | 69 | 0.749 | 205/318 = **0.645** | **↑ / ↑** |

p=41 and p=47 were run AFTER the model was fixed (this session,
`scripts/w2_class_scan_general.py`, 60-way ProcessPool over `_switched`).
The measured rate zigzags 0.35 → 0.52 → 0.27 → 0.65 and the model called
BOTH reversals in advance: down at p=41 because 7 | m brings two cubic
factors, up at p=47 because m=69 has only three orbit-factors.  A naive
trend-reading would have been wrong twice.  Measured/model ratio stays in
[0.47, 0.86]: a persistent deficit (differences are not perfectly random)
but bounded away from collapse, with the variation structure entirely the
factorization of m.  Also uniform across all four primes: EVERY class
element is Max− under switching (n_eigen = n_class), and the in-U fraction
drifts slowly (0.320, 0.294, 0.275, 0.282).

## The remaining gap (proof-shaped)

w_t mod f is NOT low-degree in the conic parameter t — the difference bits
carry χ(γj+δ) switches, so the residue is a mixed character sum and the
"few roots" closure does not apply.  What remains is a second-moment
estimate: for each irreducible f | g,

    |B_f| := #{t : f | w_t}  ≲  2^{−deg f}·|T|·(1+o(1)),

via Σ_t |w_t(ζ_f)|² (Weil-class input).  Then Σ_f |B_f| < |T| leaves a
surviving unit-content witness, closing W2 at that p (small p by census).
Risk: Φ_3 (present for every p≠3) contributes the largest bad set and the
factor count grows with ω(m), so the union bound needs genuinely tight
constants, not O(·).

## Next steps

1. γ(p) exactly at p ≤ 41 (one gcd chain over a spanning set of
   U-differences).  If γ = (X+1)^c uniformly, W2 is subsumed by the named
   conjecture γ-form.
2. Per-factor breakdown of the p=41 misses (which f kills each) against
   the per-factor model 2^{−deg f}; mismatch localizes the structure a
   proof must use.
3. The second-moment estimate itself.
