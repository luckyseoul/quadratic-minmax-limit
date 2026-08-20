# Method note: why ~500 propositions never moved a flag

Written for whoever picks this up next. This is not about any single proposition;
it is about the loop that keeps producing them.

## The failure loop

The repo contains hundreds of numbered propositions. Many carry a docstring that
says, in the author's own words, some version of:

> "Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing / ... flags."

and often also:

> "Not proved for p >= 11."

Two concrete examples, both of which sat **uncommitted** on disk until 2026-08-20:

- **15.496** — Fourier of the 3-linear `A_r` is affine in `Re J(chi,psi)` and the
  15.302 weight `sigma`. Four tests pass. Self-declares that it flips nothing.
- **15.530** — `n_R = (p-3)^2` at p=3,5,7; self-declares "Not proved for p>=11".
  Its own test **fails**: `prove_A()` returns `proved=False`.

That is the loop:

1. Pick a substructure that is computable at p <= 7.
2. Verify it at p = 3, 5, 7.
3. Notice it does not generalise; write "not proved for p >= 11" in the docstring.
4. Ship it as a numbered proposition anyway.
5. Repeat.

Each step is locally reasonable. The aggregate is ~500 propositions and three
flags that never moved.

## Why small-p census cannot close anything

`fable.md`'s acceptance bar is explicit: *"never a census at p <= 7 standing in for
general p."* This is not pedantry. `Max+` is a different object at small p —
`k=4` solution counts run `90q, 480q, 168q, 216q` for p = 7, 11, 13, 17 and then
**terminate at zero** for p = 19 (confirmed 2026-08-20 by uncapped CPU DFS on all
three tested subsets, with a positive control at p=17 returning exactly 2312).

A structure fitted to p = 3,5,7 is fitted to the regime *before* that collapse.
15.530's `n_R = (p-3)^2` is a three-point fit. Three points determine a quadratic
exactly, so the fit carries no evidence at all.

## What counts as progress

Per `fable.md`: flip a flag, or kill a route by a general counter-mechanism.
Nothing else. Two additions that this session found useful:

3. **Reduce a leftover to a single named estimate.** Not another side fact, but a
   statement of the form "leftover X follows from bound Y", with Y isolated.
4. **Remove a shared blocker.** If one obstruction blocks several leftovers,
   attacking it beats attacking any leftover directly.

## Worked example: the 2026-08-19/20 session

For contrast, not self-congratulation. What was done, and why each part counted:

- **Removed the shared blocker.** `fable.md` states "Max+ is enumerable only for
  p <= 7" and calls general-p Max+ moments "plausibly the single underlying
  problem". `Max+` at p=11 is now fully enumerated (37,457,112 in the eps=+1 half)
  by polynomial-profile stratification, validated by reproducing the known `Nh`
  exactly at p=5 (130) and p=7 (5726). That is (4).
- **Killed a route.** No bound using only `dim Z`, `tr(Phi)` and `tr(Phi^2)` can
  close leftover 1: the optimal such bound gives `lambda_max(K) <= 22.03` at p=5
  and `15.06` at p=7 against a requirement of `<= 2`, and **the gap widens with
  p**. That is (2) — a general counter-mechanism, not a failed attempt.
- **Reduced a leftover to a named estimate.** Leftover 1 follows from
  `var <= 32(n+10)^2/(n-6)^3` together with `mult(lambda_min) >= n`. Both hold at
  p=7 (+32.7% slack) and p=11 (+91.0%). That is (3).
- **Connected two leftovers.** `tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2`,
  verified exactly against the known spectra at p=5 and p=7. Leftover 1's variance
  is a function of leftover 3's four-point tensor, so one estimate serves both.

None of this flipped a flag. The session did **not** close a leftover, and p=11 is
still a census. But each item either removed an obstruction, eliminated a class of
attempt, or named the remaining gap — rather than adding proposition 501.

## Rules to break the loop

1. **Before writing a proposition, state which flag it flips or which route it
   kills.** If the answer is neither, it is a lab note. Put it in `evidence/`,
   not in `src/` with a number.
2. **A verified-at-p<=7 result is data, not a proposition.** Record the numbers.
   Do not name it, number it, or give it a `proved` field.
3. **Do not ship a `proved=False` proposition as if it were a result.** 15.530
   does this. Its tests are now `xfail(strict=True)` so the work stays tracked
   without redding the suite; `strict` means a genuine fix surfaces as XPASS.
4. **Attack the shared blocker.** Leftovers 1 and 3 are moments of one tensor.
   Bounding `||M||_F^2` at general p serves both. Another Fourier identity for
   `A_r` serves neither.
5. **Use positive controls.** Agreement between two methods that both return zero
   proves nothing — they may share a blind spot. This session nearly shipped the
   p=19 result on three zero-vs-zero agreements before adding a p=17 subset where
   the answer was known to be nonzero (2312, matched exactly). Only then did the
   p=19 zeros mean anything.
6. **Commit working code.** The GPU fixes that made p=11 tractable sat uncommitted
   in a clone on `tmpfs` for a full day. A reboot would have destroyed the only
   copy of the thing that made the result possible.

## The single next target

Bound `||M||_F^2` — equivalently `tr(Phi^2)`, equivalently the spectrum variance —
at general `p`. It is the one blocking quantity for leftover 1 via rule (3)'s
reduction, and it is the same four-point tensor leftover 3 needs.

Known exact values to fit or bound against:

| p | tr(Phi^2) | tr(Ghat^2) |
|---|---|---|
| 5 | 85248/13 | 443289600 = 2^12*3^2*5^2*13*37 |
| 7 | 3545625600/167281 | 2779770470400 = 2^14*3^2*5^2*7^2*11*1399 |
| 11 | 2440162570133760/20130785689 | 680278281952170147840 = 2^16*3^3*5*11^2*61*139*181*414061 |

No p-formula is apparent. `fable.md` advises seeking one for `eig(Ghat)` rather
than for `lambda`, since the former are plain integers. Note also that only 3 of
11 eigenvalue clusters at p=11 have `4p | eig(Ghat)`, so the p<=7 pattern where
all of them do is an arithmetic accident that stops at p=11 — do not build on it.

An upper bound suffices. An exact formula is not required.
