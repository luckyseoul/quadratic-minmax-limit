# Leftover 3 → ONE scalar: the w-line convolution and the fiber budget

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

## 5. The single remaining estimate

> **Leftover 3 follows from:  max over fibers f of |ν̂_f| ≤ C/p⁴ with
> C = 110**, where ν̂ = ν/N on the both-squares locus.  Everything else in
> §4 is proved or data-free-computable.

**Honest status of the hypothesis:** M·p⁴ = 96.2, 109.0, 98.6 at
p=5,7,11 — three points.  Per `METHOD_why_500_props_never_moved_a_flag.md`
a three-point fit carries no evidence, and this note does NOT claim it.
It is named here because it is now the *only* missing input, it is a
magnitude bound on a single scalar (not a closed form, not a census), and
C is generous: any C ≤ 115 works from p=29 up, and the p=11..23 rows
tolerate C up to ~8.1e-3·11⁴ = 119 (p=11 is the binding case).

Testing it out-of-sample at p=13 needs only max|ν| on the 41-point locus
— Route A of `p13p17enum.md` (orbit-representative moments), NOT a full
Max± census.  That is the next experiment, and it is the one that decides
whether this route lives.

## 6. Also settled

The joint degree-4+6 kernel at p=7 projects onto the (μ₄,δ₄) block with
rank 2 = the degree-4 kernel dimension: degree 6 adds nothing to the
4-point sector at p=7.  Degree escalation stays dead (see the K₄ table in
`NOTE_2026-08-21_leftover3_contraction_closure.md`).
