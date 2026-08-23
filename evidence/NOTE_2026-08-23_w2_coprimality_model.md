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

---

## Per-factor breakdown (same day): the union-bound proof shape is DEAD

`scripts/w2_perfactor_breakdown.py` records, for every in-U class element,
which orbit-product f_O shares a factor with the content.  Two corrections
to the section above, both against my own earlier reading.

**Correction 1 — model bug.**  The W2 test is `gcd(content, f_O) = 1`, i.e.
the bad event is "divisible by SOME irreducible in O", with probability
`1 − ∏_{f∈O}(1 − 2^{−deg f})` — not `2^{−deg f_O}` ("divisible by all"),
which is what I first compared against.  Only matters where an orbit-product
is reducible; the corrected table below uses the true irreducible content
(p=41's deg-72 orbit is six deg-12 irreducibles; p=47's deg-66 is
11,11,22,22).

**Correction 2 — "prime-power localization" is dead.**  I had guessed the
anomaly sat on Φ_{ℓ²}-type factors.  At p=41, m=105 = 3·5·7 is squarefree,
yet the deg-6 and deg-12 factors are inflated 17–29×.  The inflation is not
a prime-power effect.

| p | orbit (irred degs) | model | measured | ratio |
|---|---|---|---|---|
| 17 | [2] | 0.2500 | 0.6531 | 2.61 |
| 17 | [6] | 0.0156 | 0.4694 | **30.0** |
| 31 | [2] | 0.2500 | 0.3562 | 1.42 |
| 31 | [4] ×3 | 0.0625 | 0.055–0.089 | 0.88–1.42 |
| 41 | [2] | 0.2500 | 0.6414 | 2.57 |
| 41 | [3] ×2 | 0.1250 | 0.3840 | 3.07 |
| 41 | [4] ×3 | 0.0625 | 0.110–0.131 | 1.76–2.09 |
| 41 | [6] ×2 | 0.0156 | 0.2700 | **17.3** |
| 41 | [12]×6 | 0.0015 | 0.0422 | **28.8** |

p=31 fits the independence model almost exactly (0.88–1.42); p=17 and p=41
do not.  So "differences behave like random elements of R" is FALSE in
general — it is a good aggregate predictor and a bad per-factor one.

### The kill

Summing the **exact measured** bad-rates per orbit:

| p | Σ_O P[bad_O] | true union (miss rate) | union bound |
|---|---|---|---|
| 17 | 1.123 | 0.653 | **FAILS** |
| 31 | 0.582 | 0.479 | ok |
| 41 | 2.363 | 0.734 | **FAILS** |

The bad events are **strongly positively correlated** — at p=41 the sum
overshoots the true union by 3.2×.  Therefore the proof shape proposed in
the consult ("bound each |B_f| by a second-moment/Weil estimate, sum, get
< |T|") **cannot close W2 — not even given perfect per-factor bounds.**
Any union-bound-over-orbits route is dead at p=17 already.

### What survives, and the revised gap

* Unit-content elements remain a healthy fraction at every censused prime:
  0.347, 0.521, 0.266, 0.645 (p = 17, 31, 41, 47).  No drift toward 0.
* The aggregate Euler product still tracks the *direction* of change
  correctly (it called both reversals out of sample) because it is driven
  by the factor-degree profile of m; it systematically **over**estimates
  the unit-content rate (ratios 0.47, 0.84, 0.58, 0.86) precisely because
  the bad events correlate positively.
* **Revised gap:** a bound on the union / on the unit-content indicator
  *directly* — e.g. a second moment on the count of unit-content witnesses
  over the class — rather than factor-by-factor.  Inclusion–exclusion over
  ~9 orbits with these correlation sizes is not obviously tractable either.

W2 stays OPEN.  No flag flipped.
