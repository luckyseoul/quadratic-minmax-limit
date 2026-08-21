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

## 8. Reassessment against the ACTUAL leftover-3 target (L, not 2/n)

§§4–7 measured the budget against `2/n` — the 15.191 convenience target.
Leftover 3 itself needs only the weaker `|μ| ≤ L = (p−2)/2p²` on |κ|=1
(15.275/15.268: on |κ|=1 ν=0, 3A+B>0 ⟺ G>T, and |μ|≤L < |T| suffices;
`2/n < L` for p≥7).  Same identity (B), same A(p), weaker target:

    |μ| ≤ 1/p² + (1/p)·M·A(p) ≤ L   ⟸   M ≤ (p−4)/(2p·A(p)) ≈ 1/(6p).

| p | allowed_L·p⁴ | rigorous LB·p⁴ | true·p⁴ | truth/allowed |
|---|---|---|---|---|
| 5 | 20.8 | 50.00 | 96.2 | 4.62 — census (15.275 L) |
| 7 | 73.5 | 62.36 | 109.0 | 1.48 — census (109/2863 < 5/98) |
| 11 | 423.5 | 91.79 | 98.6 | **0.233** |
| 13 | 659.1 | 107.17 | — | ~0.17 (est from LB+7.4%) |
| 17 | 967.7 | 138.39 | — | ~0.15 |
| 19 | 1254.7 | 154.13 | — | ~0.13 |
| 23 | 2266.4 | 185.76 | — | ~0.09 |
| 29 | 3586.6 | 233.39 | — | ~0.07 |
| 31 | — | 249.30 (·/p = 8.042) | — | — |

The §7 falsification does NOT touch this: the requirement `M ≲ 1/(6p)` is
three orders above the true scale `M ~ 8.6/p³`, the rigorous LB sits far
below it (a lower bound cannot threaten an upper-bound requirement it does
not reach), the margin grows like p²/50, and the budget fails at exactly
the two primes (5, 7) that census already covers.

**Status per the pre-asymptotic-fit discipline:** the inequality chain is
proved arithmetic (identity (B) exact; A(p) data-free); `M ≤ 1/(6p)` is an
OPEN estimate whose truth is verified only at p=5,7,11 (holding at 11, the
only budget-relevant prime with known truth).  Candidate proof route, so
far untried: split ν = (1/p)χ⋆μ_part + (1/p)χ⋆(μ−μ_part); the μ_part term
is Max-free (15.247 A) and Weil-bounded at O(1/p²); the residual needs
only a crude O(1) line-L¹ bound — the first leftover-3 route with room
for crude bounds.

Open verification items before any flip claim:
1. Trace the |μ|≤L ⇒ type_I sufficiency WIRING in code (15.275 imports).
2. Extend the rigorous LB to p=37..47 and check LB < allowed_L throughout.
3. The Weil bound on χ⋆μ_part (the real work).

## 9. Sufficiency wiring traced; the complete reduction (both class types)

Traced in `src/e1_gmin_m4_prop15275.py` (T_bitight, L_abs_gmin,
pairing_G_values, three_AB_kappa3_from_mu_nu, and the 3A+B identity
lemma).  Signs: T(p) = −(p−2)/(p(2p−1)) and L(p) = −(p−2)/2p² are both
NEGATIVE; "L>T" and "|L|<|T|" are the same statement.

`type_I_aut_e_3AB_positive_general` needs, on every Aut_e far class:

* **|κ|=1 classes.**  The three pairing G's are {+μ,+μ,−μ} (up to the
  class sign), and 3A+B>0 ⟺ G>T.  Worst pairing needs −|μ| > T, i.e.
  |μ| < |T|; |μ| ≤ |L| closes it strictly.
* **|κ|=3 classes.**  3A+B = 2n_O[(p−2) + p·χd·w]/(p²−1),
  w = (2p−1)μ + (p−2)ν.  Frame fact: on a κ=+3 four-set EVERY edge has
  χd = +1 (all pairs involve ∞, {0,1} with χ(−1)=1, or χ(w), χ(1−w) —
  all +1 on the locus).  So the condition is just w > −(p−2)/p.
* **κ=−3 four-sets do not exist** (theorem, via sharp 3-transitivity:
  frame κ = 1+χ(w)+χ(1−w) ∈ {3,1,−1}).  No third case.

Both conditions route through identity (B).  For (ii), the convolution at
locus points gives μ(l) ≥ 3/p² − M·A₃/p, and with ν ≥ −M:
w > −(p−2)/p ⟸ M ≤ c₂(p) := [(p−2)/p + 3(2p−1)/p²] / [(2p−1)A₃/p + (p−2)],
where A₃(p) = max over locus l of Σ_f|S_f(l)|.  Combined, data-free:

| p | A₁ | A₃ | c₁ (κ=1) | c₂ (κ=3) | binding·p | truth M·p |
|---|---|---|---|---|---|---|
| 5 | 3 | 4 | 3.33e-2 | 1.65e-1 | 0.167 | 0.880 — census |
| 7 | 7 | 10 | 3.06e-2 | 6.41e-2 | 0.214 | 0.318 — census |
| 11 | 11 | 16 | 2.89e-2 | 3.39e-2 | 0.318 | **0.074** |
| 13 | 15 | 28 | 2.31e-2 | 1.99e-2 | 0.259 | ~0.051 |
| 17 | 33 | 46 | 1.16e-2 | 1.17e-2 | 0.197 | ~0.030 |
| 23 | 51 | 74 | 8.10e-3 | 7.05e-3 | 0.162 | ~0.016 |
| 29 | 85 | 120 | 5.07e-3 | 4.32e-3 | 0.125 | ~0.010 |
| 37 | 125 | 192 | 3.57e-3 | 2.67e-3 | 0.099 | ~0.006 |
| 47 | 217 | 310 | 2.11e-3 | 1.65e-3 | 0.077 | ~0.004 |

Budget fails at exactly the census primes 5, 7 and closes from p=11 with
margin growing like ~p/2.4 (4.3× at p=11 → 20× at p=47).  A(p) grows like
~p²/10 because the per-fiber triangle costs the FIBER COUNT (~p²/24, each
|S_f| = O(1) by Weil); hence binding ~ 3.6/p² vs truth ~ 8.6/p³.

### The complete reduction of leftover 3 (current sharpest form)

> **leftover 3 ⟸ census at p=5,7 (done) + `M(p) = max_f|ν̂_f| ≤
> binding(p) ≈ 3.6/p²` for all p ≥ 11**, binding(p) data-free
> (table above to p=47; general p needs the A(p) growth made rigorous,
> a per-fiber Weil bound |S_f| ≤ 2√|f|+O(1) plus the fiber count).

Status: OPEN estimate; truth known only at p=5,7,11 (holds at 11 with
4.3× margin); rigorous LB ~ 8.05/p³ sits a factor ~p/2.2 below the
requirement at every computed prime — no computable quantity currently
threatens it, unlike the falsified §5 form.  Proof target: M = O(1/p²)
with constant < 3.6 via ν = (1/p)χ⋆μ_part + (1/p)χ⋆(μ−μ_part)
(μ_part Max-free from 15.247 A; Weil on the first term; the second needs
line-L¹ = O(p) — crude-tolerant).

## 10. First step of the proof route: the μ_part split, evaluated on data

Split ν = ν_part + ν_res with ν_part := (1/p)·K(μ_part), where μ_part is
the Max-free 15.247 A particular solution transported to the line:
μ_part(w) = [(p²−1)κ(w) − 2φ(w) − 8p·1_{κ=3}]/(p²(p²−5))  (star = +4 on
the locus in the frame, 0 off it).  ν_part is a FINITE EXPLICIT character
sum — no Max± input.  Exact evaluation against the ν data:

| p | max\|ν_part\|·p³ | max\|ν_res\|·p³ | true M·p³ | residual/binding | sign of ν_part on locus |
|---|---|---|---|---|---|
| 7 | 11.45 | 9.21 | 15.57 | 0.878 | all negative |
| 11 | 11.38 | 3.97 | 8.96 | **0.103** | all negative |

- **The uniform negativity of ν is carried by the explicit term**:
  K(μ_part) < 0 on the whole locus at both primes.
- ν_part·p³ is stable (~11.4) — the explicit term captures the true scale.
- The open residual uses 10% of its budget at p=11 and is falling.
- Caveat: the naive leading term −2q/(p²(p²−5)) is NOT the full closed
  form (K(μ_part)/(−2q/D) = 2.86, 5.45 at p=7,11 — growing); the complete
  evaluation of K(μ_part) is classical character-sum algebra (Jacobsthal +
  one elliptic-type double sum), still to be done.  Two data points on the
  residual decay — per §7's lesson, NOT extrapolated.

Refined open estimate (current sharpest form of leftover 3):

> |ν_res(l)| = |ν(l) − (1/p)K(μ_part)(l)| ≤ binding(p) − |ν_part|(p)
> on the locus, for p ≥ 11 — with ν_part explicit and binding(p)
> data-free.  Measured usage at p=11: 10%.

Rigorous LB tail (frame_line): p=37 gives LB·p⁴ = 297.09, LB·p⁴/p = 8.03
(seventh point on the flat line); LB/binding ≈ 0.06 — no tension.
