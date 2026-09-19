# Assessment 2026-09-18: cross-transfer / C_cross sweep — duplication check (option C)

**Purpose.** The handoff `HANDOFF_2026-09-18_cross_transfer.md` (section 5 & 7)
requires re-evaluating duplication against `solution.md §10 (bilinear floor)` and
the shield note *before continuing*, and then picking one of {A, B, C}. This
document records that check and the decision.

**Status:** COMPLETE. This is a routing decision, not a mathematical claim.

## 1. Data-integrity finding (critical)

The measurement scripts listed in the handoff (§6) **do not exist on soulkiller**:

  find / -name "mo_*.py"   ->  (no hits)

and there are no result data files (`.npy`/`.npz`/sweep logs) dated after
2026-09-17 other than the handoff doc itself. The only 2026-09-18/19 file in
the project is `HANDOFF_2026-09-18_cross_transfer.md`.

Consequence: the empirical findings in §3 of the handoff (A-term p95s; C_cross
values; the 2n transfer max 0.898 at n=13) exist **only as prose in the handoff
doc**. They cannot be re-run or independently reproduced. This materially weakens
their status as evidence and rules out a cheap "just re-run" for option A.

## 2. What the REAL frontier is (from the two-ray note + solution.md)

From `evidence/NOTE_2026-09-01_ORIGINAL_LIMIT_TWO_RAY.md` (two-ray theorem) and
`solution.md` Propositions 6.5-6.10 (through 15.774):

* **Two-ray convergence criterion (proved).** α_n converges iff we have
  Dini-summable *doubling* (multiplier two) and *tripling* (multiplier three /
  1:2 split) estimates:
      m_{2n}^{2/3} <= 2 m_n^{2/3} + 2 n η(n),   m_{3n}^{2/3} <= 3 m_n^{2/3} + 3 n η(n)
  (equivalently a doubling bound plus the 1:2 split H(3n) <= H(n)+H(2n)).
* **Multiplier two** is narrowed (Prop 6.5, 6.6) to the **explicit residue
  (6.20)** — pairs (x,y) with h_x,h_y>ρn, d_H(x,y)(n-d_H(x,y))>ρn²,
  h_x h_y>(ρ/4)n², |Q_A(x)+Q_A(y)|>2√2·M-n^(3/2), ρ=a², a=(√2-1)/π.
  Prop 6.6 builds a skew R shielding every pair *except* those in (6.20).
* **Multiplier three / 1:2 split** is OPEN (Prop 6.7 tripling diamond;
  Prop 6.8 1:2 reduction, residual 6.42-6.43).
* **The live gate:** close residue (6.20), then prove the 1:2 split.

## 3. §10 bilinear floor (the referenced shield note) — what it does/does not say

`solution.md §10` gives the exact two-block identity

    Φ(S) = max_{x,y} ( |Q_{A1}(x)+Q_{A2}(y)| + |x^T B y| )          (10.1)

so internal and cross energies cannot *cancel* in the absolute maximum. §13
(13.3 table) records the separate-norm cross contributions: i.i.d. random gives
√log2, Hadamard (k=2) pays 1/(2√2)≈0.354, etc. **Key scope sentence:**
"(10.1) does **not** rule out the multiplier-two target in Proposition 6.3. ...
A coupled state/profile construction remains live."

So §10 does **not** foreclose the 2n-transfer route; it only kills *uncoupled*
separate-norm bounds. The shield note (RG2_EQUAL_ENDPOINT_PALEY_SHIELD.md)
gives the exact remaining multiplier-two diamond:

    |Q_A(x)+Q_A(y)| + |x^T R y|  <=  2√2·M + n^(3/2) Ω(n)       (exact, up to Dini tail)

and confirms the residue (6.20) is the open part.

## 4. Is the sweep a duplicate?

The "2n object" measured in the sweep is

    max_{x,y} | x^T Z y + Qx - Qy |   (real Gaussian structured Z, conference A)

which is the **actual 2n norm of a specific conference-based 2n construction**,
normalized by n^(3/2). That is *precisely* the quantity the two-ray / multiplier-
two criterion bounds. So the sweep is **not off the frontier** (contrary to the
handoff's own §1 claim); it re-measures the *same frontier object* (multiplier
two / 2n transfer) that residue (6.20) targets, but by a different method:

* residue (6.20)      = **analytic** route (construct R, prove the diamond).
* the sweep            = **empirical** route (measure the 2n norm of one construction).

These are complementary views of the *same* target, not a duplicate sub-task. The
handoff's worry ("re-measure the same object") is correct in that it is the same
frontier object, but it is *not* a no-op duplicate: it gives a numeric value for
the transfer constant.

## 5. Why it still cannot advance the frontier as-is

1. **Not a proof.** 0.898 < 0.9216 at a single n=13 is a measurement, not a
   bound on m_{2n} for all n, and not a closure of residue (6.20) or the
   doubling estimate.
2. **Tiny n, thin margin.** 0.898 vs 0.9216 is a 2.5% gap; the sweep itself
   notes it "could close at larger n."
3. **Code lost.** No `mo_*.py` and no result data remain (see §1). Re-doing the
   n=13 measurement requires rebuilding the whole apparatus from prose.
4. **Option A (verify at n>13) is not feasible "as described."** The handoff
   itself says n>13 is intractable by full cube and needs a *different* bound;
   with the code gone, A becomes "build a new n>13 2n-transfer bound from
   scratch" — a new project, not a continuation.

## 6. The one useful artifact the sweep leaves

The measured 2n transfer ≈ 2√2·B (0.898..0.9216, B=0.32584) at n=13 is a
**concrete target value** for the analytic closure of residue (6.20): the skew
R construction in Prop 6.6 must ultimately match a 2n transfer near 2√2·B ≈
0.9216. That is useful context for the main-line work, even though it is not
evidence on its own.

## 7. Decision

Per the handoff's own section 5 ("check §10 + shield note before continuing")
and section 7 ("pick ONE"):

  C (re-assess duplication)   -> DONE.  (this document)
  A (re-verify n>13 + 3n)    -> NOT recommended: code lost, n=13 tiny, thin
                               2.5% margin, would be a new project.
  B (return to the frontier)  -> **RECOMMENDED.** Close residue (6.20) + the
                               1:2 split. Record the sweep as a signal (this
                               file) and do not re-run it.

**Recommended next step: B.** Resume the main-line multiplier-two work at the
explicit residue (6.20) from Prop 6.6 (the "exact two-block identity + residue
is open" state in §13 / the shield note), with the n=13 value 2n-transfer≈
2√2·B as a sanity-check target, not a claim.

---
Assessment authored by: Bonsai agent (nuka), 2026-09-18, working from
soulkiller `/home/nick/quadratic-minmax-limit/`.
Sources read: NOTE_2026-09-01_ORIGINAL_LIMIT_TWO_RAY.md,
NOTE_2026-09-01_RG2_EQUAL_ENDPOINT_PALEY_SHIELD.md, solution.md §10/§13 &
Prop 6.5-6.6, STATUS.md.
