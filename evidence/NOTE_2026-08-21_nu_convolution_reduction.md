# Leftover 3 → ONE scalar: the w-line convolution and the fiber budget

> **FALSIFIED SAME DAY — read §7 before §5.**  The scalar hypothesis stated
> in §5 (`M ≤ 110/p⁴`) is **false**: a rigorous data-free lower bound gives
> `M·p⁴ ≥ 138.39` at p=17.  The scaling itself was wrong — `M ≳ 8/p³`, not
> `C/p⁴`.  §§1–4 (the identities, the fibers, the budget inequality) stand;
> §5's reduction does not.

Date: 2026-08-21.  Supersedes the "signed ν required" scoping in
`NOTE_2026-08-21_signed_nu_dataset.md` (see §0).  Scripts:
`/mnt/storage/e1work/leftover3_mu/{line_convolution,fiber_budget}.py`.
Data: `evidence/fiber_budget_dataset.json`, `evidence/nu_of_w_dataset.json`.
**No flag flipped.**  Leftover 3 is now reduced to a single scalar estimate.

## 0. Correction to the previous note

The earlier claim "no magnitude-only estimate on ν can close leftover 3"
was proved only inside the **affine-parametrization** hypothesis set
(keep true |ν| per free δ-orbit, flip signs).  That set discards the fiber
and character-sum structure.  The convolution below retains it and is a
magnitude-only estimate that DOES close for p≥11.  The sign-flip
counterexamples remain valid against the affine route; they do not apply
here.  Scope the earlier statement accordingly.

## 1. Why a "closed form for ν(w)" is the wrong target

At p=11 the fiber value is ν̂ = −2866/425649 in lowest terms, and 425649
carries a large prime factor of N.  No Weil-type character-sum expression
has such a denominator: a closed form would have to be for the integer
sums, i.e. the Max+ counting problem again.  Abandon "find ν(w) exactly".

## 2. The convolution identities (new, exact, Max-free)

Frame S={∞,0,1,w}; κ(w)=1+χ(w)+χ(1−w); φ(w)=−a_q(Legendre E_w).  For all
l ∈ F_q∖{0,1} (verified with **zero violations**, all 189 instances at
p=5,7,11):

    (A)   p·ν(l) = Σ_{w≠l} χ(l−w)·μ(w)
    (B)   p·μ(l) = Σ_{w≠l} χ(l−w)·ν(w) + κ(l)·N/p
    (C)   Σ_w μ(w) = Σ_w χ(w)μ(w) = Σ_w χ(1−w)μ(w) = N

So (μ,ν) is a conference-eigenpair on the punctured w-line, and **ν is
supported only on the both-squares locus** {χ(w)=χ(1−w)=1} (=|κ|=3).
(B) is the lever: it expresses μ on |κ|=1 as κ/p² plus a convolution of
the *small-support* function ν.

## 3. Fibers are data-free

ν is constant on the orbits of ⟨anharmonic S₃, Frobenius⟩ acting on the
locus.  Predicted-vs-measured fiber sizes:

| p | predicted | measured (from ν data) |
|---|---|---|
| 5 | 2,3 | 2,3 |
| 7 | 2,3,6 | 2,3,6 |
| 11 | 2,3,6,6,12 | 2,3,6,6,12 |

Exact match — so fibers, and hence the character sums S_f(l)=Σ_{w∈f}χ(l−w),
are computable at ANY p without touching Max±.

## 4. The budget, and the complete case analysis

From (B) with |κ(l)|=1, writing M := max_f |ν̂_f| and A(p) := max_{|κ(l)|=1}
Σ_f |S_f(l)|:

    |μ(l)| ≤ 1/p² + (1/p)·M·A(p),   and |μ|≤2/n holds as soon as
    M·A(p) ≤ (p²−1)/(p(p²+1)).

| range | mechanism | status |
|---|---|---|
| p = 5, 7 | census; p=5 is a from-C theorem (15.275 L), p=7 max\|μ\|=109/2863 < 1/25 | **done** |
| p = 11,13,17,19,23 | computed A(p) = 11,15,33,41,51; slack 1.08, 1.32, 1.34, 1.51, 2.16 | **closes, given M** |
| p ≥ 29 | trivial \|S_f\|≤\|f\|, Σ\|f\|=\|locus\|≤(p²+3)/4 ⇒ need C ≤ 4p³(p²−1)/((p²+1)(p²+3)), RHS ≥ 115 at p=29 and monotone increasing | **closes, given M** |

Slack grows monotonically (0.35, 0.43 at p=5,7 — hence the census there —
then 1.08 → 4.35 by p=47), and asymptotically slack ~ 4p/C.

## 5. The single remaining estimate  — ***FALSIFIED, see §7***

> ~~**Leftover 3 follows from:  max over fibers f of |ν̂_f| ≤ C/p⁴ with
> C = 110**~~ — retracted.  The premise is false for p ≥ 17.

Retained only as a record of what was tried.  The original reasoning: with
`M := max_f|ν̂_f|` bounded by `C/p⁴`, §4's case analysis (census at p=5,7;
computed `A(p)` at p=11..23; trivial tail at p≥29) closes leftover 3.  The
empirical basis was `M·p⁴ = 96.2, 109.0, 98.6` at p=5,7,11 — three points,
flagged at the time as carrying no evidence per the METHOD note, and indeed
carrying none.  See §7.

## 6. Also settled

The joint degree-4+6 kernel at p=7 projects onto the (μ₄,δ₄) block with
rank 2 = the degree-4 kernel dimension: degree 6 adds nothing to the
4-point sector at p=7.  Degree escalation stays dead (see the K₄ table in
`NOTE_2026-08-21_leftover3_contraction_closure.md`).

## 7. Falsification (same day, data-free)

### The pinned functional gives a rigorous lower bound on M

The degree-4 system pins exactly one linear functional of the ν-vector:
`Σ_f c_f ν_f = V`, with the annihilator `c` and the value `V` both computed
from the equivariant system alone — **no Max± data**.  Hence

    M = max_f |ν̂_f|  ≥  |V| / Σ_f |c_f|  =:  LB(p)          (rigorous)

`LB` is normalization-independent (`V` and `c` scale together).  Computed by
`scripts/frame_line_system.py`; validated against the four-set implementation
at p=5,7,11,13 (kernel dims 1,2,4,6 and LB·p⁴ = 50.00, 62.36, 91.79, 107.17
reproduced exactly).

| p | LB·p⁴ | allowed·p⁴ | LB/allowed | **LB·p⁴/p** |
|---|---|---|---|---|
| 5 | 50.00 | 38.46 | 1.300 | 10.00 |
| 7 | 62.36 | 47.04 | 1.326 | 8.91 |
| 11 | 91.79 | 119.02 | 0.771 | 8.34 |
| 13 | 107.17 | 144.74 | 0.740 | 8.24 |
| 17 | **138.39** | 147.85 | 0.936 | 8.14 |
| 19 | **154.13** | 166.37 | 0.926 | 8.11 |
| 23 | **185.76** | 237.67 | 0.782 | 8.08 |
| 29 | **233.39** | 286.25 | 0.815 | 8.048 |

(Self-consistency: LB exceeds the budget exactly at p=5,7 — the two primes
§4 hands to census — and falls below from p=11 on.  Not built in.)

### Two conclusions

1. **`M ≤ 110/p⁴` is false.**  LB·p⁴ = 138.39 > 110 at p=17, and grows.
2. **The scaling was wrong, not just the constant.**  `LB·p⁴/p` →
   8.34, 8.24, 8.14, 8.11, 8.08, 8.048: flat.  So `M ≳ 8/p³`.  No constant
   `C` bounds `M·p⁴`, so the §4 tail argument (which needs a uniform `C`)
   has no premise.

### Does the per-prime budget survive?

Not formally dead: `LB < allowed` at every p ≥ 11, because `allowed·p⁴/p`
also scales like p (range 8.70–11.13).  But the margin is gone.  At the
binding primes p=17,19 the gap is ~7%, and the measured LB-to-true gap at
p=11 is **7.4%** (true/LB = 98.58/91.79 = 1.074).  Propagating it:

    estimated true/allowed  =  1.005 (p=17),  0.995 (p=19)

i.e. the budget most likely **fails at p=17**.  Deciding it needs true `M`
at p=17 — Max± at p=17, out of reach.  A route needing an estimate accurate
to 7% is not one a crude bound closes.

### Methodological note

The three-point fit at p=5,7,11 sat in the pre-asymptotic regime (`M·p⁴`
appears flat there only because `p` has not yet grown).  The note hedged
that it was "recorded only as scaling, NOT fitted" and cited the METHOD
warning — the hedge was insufficient.  The correct action was to compute a
fourth point *before* writing §5.  The frame-line reduction that makes this
a 30-second computation per prime was one refactor away and was built two
hours later.  **Compute the next point before naming the estimate.**

### What survives

- §2 convolution identities (A),(B),(C): exact, 189/189 verified.
- §3 fibers = anharmonic S₃ × Frobenius orbits; ν-death is the V₄ pairing
  mechanism (15.268 generalized: ν dies on a fiber iff some pairing has
  χ(det)=−1, i.e. iff |κ|=1 — verified set-wise).
- §4's budget *inequality* (only the hypothesis feeding it is dead).
- `scripts/frame_line_system.py`: ~120× faster than the four-set code,
  reaches p=31 in ~2 min where the four-set path could not reach p=17.
