# Route plan: prove the stronger value-specific claim \(L=1/2\)

> **Scope correction (2026-09-01).**  This is an optional Paley route plan,
> not the binding acceptance test for MathOverflow 413935.  The original
> problem is settled by proving existence without identifying the value.
> Proposition 6.3 reduces that direct route to Dini-summable amplification at
> multipliers 2 and 3; polynomial saving is more than is needed.  Proposition
> 6.5 gives the exact equal-endpoint diamond, and Proposition 6.6 proves it
> outside the explicit residue (6.20).  See
> `LONG_HORIZON_GOAL.md`.

Use this file as the `/goal` body only when deliberately selecting the
stronger Paley route. Example:

```
/goal complete GOAL.md in /home/nick/quadratic-minmax-limit
```

Do not attach a `--budget` flag. Do not `/goal clear` to start fresh.

## Goal kind
code-change

## Acceptance criteria

Current partial advance (15.775): the complete layer `t=q+3,k=5p+5`
is closed for every prime `p>=259201`, with all boundary sizes, and the
same argument excludes both signed shell floors `r=3,4,5` on
`46656*h^3<=p^3*(p-1)` when `h=r mod 2`. See
`evidence/NOTE_2026-09-04_EVENTUAL_FIRST_LAYER_CLOSE.md`.
These are proved infinite-family exclusions, not all-size localization:
the global residual-(ii) and minimal-four-gap bridge criteria below remain
open. No acceptance Boolean is changed by the bounded support theorem.

1. For every prime \(p\ge5\), E(1) is proved Max+-free on Paley \(n=p^2+1\). The corrected gate is:
   - **Required bi-tight levels 2 and 3:** **satisfied by 15.720.** From 15.272/15.207, a centered bi-tight indicator lies in `scheme+cross`; 15.720 forces all degrees into one residue class modulo \((p^2-1)/2\), and the handshake identity excludes both levels. It also excludes *bi-tight* level 4. Generic one-sided level-4 covers exist (15.402), while the known four-line family is outside residual by a 15.272/15.588 Max+ witness; only joint residual compatibility remains open. Proposition 15.167's spectral implication is retracted. Global QVAR, principal R1, and the spectral floor are no longer acceptance units.
   - **Residual (ii)** ND (or another Max+-free kill) for **even \(k\ge4p\)**, not only \(k\in[3p+1,4p-2]\). Affine two-level at \(k\ge3p\) is already dead (15.179).
   - **Restricted Type I when Max− is not two-level \(\{-1,-3\}\): satisfied by 15.750 for `k=3p-2` and `S_G=3-2*f_e on Max+`.** Isolated-chart rigidity and parity halving close \(p\ge11\); exact integer Farkas identities close \(p=5,7\). General odd-k level-one no-descent is not implied.
   - **Lemma D** fully written and checkable: every triple of good lines is a \(k=3\) Max+ (pointwise \(z\), three-line Fourier support, \(\hat z(0)=p\)), and each locked triple spans the 2-plane \(\{x+y+z=0\}\) on its edges (amplitudes must reconcile Fejer characters with convolution). **Satisfied by 15.276 and `A3_PROOF.md`;** census rank-2 at \(p=5,7,11\) is a check, not the general proof.
   - **Minimal-four-gap implication bridge:** **OPEN.** Proposition 15.764 proves only odd `|H|<=5p`. The no-contact ranges begin at even `|H|>=4p+2` and odd `|H|>=5p+2`; even contact at two also needs a general odd-k no-descent theorem or the actual hypotheses of a proved unit. Contact alone does not enter 15.750. The former 15.774 unconditional even bound is retracted; its r4 bound holds without a signed level-two row. See `evidence/NOTE_2026-09-04_GLOBAL_BRIDGE_TYPE_I_SCOPE_AUDIT.md`. Closing residual (ii) alone is not sufficient to assert E(1).
2. The legacy-named `four_e1_units_closed()` is True **only** by importing the valid bi-tight, residual-(ii), Type-I, Lemma-D, **and minimal-four-gap bridge** units (no handwritten `return True`). `e1_closed_general()` must agree with this expanded acceptance gate. Aut-Schur **False**, Gsum **False**, pairing **False**, and the open spectral floor stay as they are but are not acceptance gates.
3. `solution.md` Main Theorem (limit) states \(L=\lim\alpha_n=1/2\) **only after** 1–2, and `evidence/share/denseness_path_package.md` § Caveats is empty or records them as closed. Soft-close (sandwich + denseness without E(1)) fails.

## Verification plan
1. **gating:** Run the corrected main-unit dump **twice**. Both runs: required bi-tight levels 2/3 True, residual (ii) True, multi-level Type I True, Lemma D True, minimal-four-gap bridge True, `four_e1_units_closed=True`, and \(L\) CLOSED. The spectral floor, `gsum_disj_lb`, Aut-Schur, and pairing may remain False.
2. **gating:** Pytest of 15.272 / 15.270 / 15.249 / 15.167 / 15.236 / 15.237 / 15.170 plus new caveat-unit tests (`pytest -n W`, \(W\) from compute-budget). Capture `{SCRATCH}/pytest_e1.txt`. Exit 0. Acceptance units must go False if the \(k\ge4p\) or multi-level argument is wrong; the optional \(\lambda_{\min}\) route does not gate them.
3. **gating:** Independent cold read of `evidence/share/denseness_path_package.md`. Q1 (Type I and the required bi-tight levels) must be proved; Q2 (E(1) / \(L=1/2\)) remains blocked until residual (ii) and the minimal-four-gap bridge both close.
4. **evidence:** Diff of the expanded acceptance units + predicate wiring + package + `solution.md`. Census \(p\le7\) must not be the only reason a “general” flag is True (except the already-finite \(p=5\) Veronese rank 65).

## Non-goals
- Aut-Schur / Jacquet / PSL-span of \(k=3\) \(F\).
- Gsum disj LB; cotangent pairing \(1^\top K^{-1}v\).
- Envelope / reflection / \(K_4\le\mathrm{Wick}_{hi}\) / \(|\mu|\le2/n\).
- Path-C / \(16N\) / Hypothesis H / 15.193 exhaustiveness.
- Re-proving sandwich, denseness Prop 6.1, \(\rho=1\), Johnson Lemma E, cost_D polynomial, \(V_+\) Fourier iso.
- Unflipping \(L\) to “start honest” and then flipping it back. If a caveat dies, record it and switch listed attack — do not cycle the Main Theorem.

## Assumed scope
Repo `/home/nick/quadratic-minmax-limit/`. Hinge already shipped: 15.272 Lemmas B–C, E, F, G Singer on \(F\), 15.207, 15.249, and Lemma D in 15.276 / `A3_PROOF.md`. Route-local leftovers: package § Caveats, `HANDOFF.md`, `STATUS.md`. Writeup: `solution.md`. Statement: \(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}|\sum_{i<j}a_{ij}x_ix_j|\); within this stronger route the target value is \(L=1/2\).

## Implementation approach within this route
If this stronger Paley route is deliberately selected, attack residual
(ii) and the minimal-four-gap implication bridge exposed by Proposition 15.764;
Proposition 15.750 closed the Type-I multi-level remainder. QVAR/R1 work
is optional mathematics, not a prerequisite for this proof. Wire flags by
import after a unit is live and failing-when-wrong. Set `solution.md` \(L=1/2\)
only after criteria 1–2. Use compute only behind a predeclared mathematical
gate.  Otherwise the direct original-problem front is Proposition 6.3, not a
finite Paley residue.

## Task checklist
- [x] Replace the invalid 15.167 bi-tight arrow with a valid all-prime obstruction for required levels 2 and 3 (15.720; bi-tight level 4 is a non-load-bearing corollary).
- [ ] Close residual (ii) for even \(k\ge4p\)
- [ ] Close the all-size minimal-four-gap bridge, including even contact and general odd-k Type-I entry
- [x] Close the restricted `k=3p-2`, full-affine Type I bad case when Max− is multi-level (15.750)
- [x] Write Lemma D: every triple + 2-plane (15.276; `A3_PROOF.md`)
- [ ] Import those units so `four_e1_units_closed()` is True for the right reasons
- [ ] Package caveats cleared; `solution.md` states \(L=1/2\) only then
- [ ] Two 15.170 dumps + honesty pytest + cold read of the package

## Risks
- The residual (ii) multi-level leftover is real. A structure increment that leaves it open fails the goal.
- The historical four-unit ledger is not exhaustive without the Proposition 15.764 bridge; residual-(ii) closure alone must not flip E1.
- Historical `e1_closed_general` wiring is not acceptance; the corrected main dictionary is authoritative.
- Do not `/goal clear` to “start fresh.”
