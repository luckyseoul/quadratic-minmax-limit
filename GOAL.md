# Plan: Settle \(L=\lim\alpha_n=1/2\) by closing the two genuine remainders

Use this file as the `/goal` body. Example:

```
/goal complete GOAL.md in /home/nick/quadratic-minmax-limit
```

Do not attach a `--budget` flag. Do not `/goal clear` to start fresh.

## Goal kind
code-change

## Acceptance criteria
1. For every prime \(p\ge5\), E(1) is proved Max+-free on Paley \(n=p^2+1\). The corrected gate is:
   - **Required bi-tight levels 2 and 3:** **satisfied by 15.720.** From 15.272/15.207, a centered bi-tight indicator lies in `scheme+cross`; 15.720 forces all degrees into one residue class modulo \((p^2-1)/2\), and the handshake identity excludes both levels. It also excludes *bi-tight* level 4. Generic one-sided level-4 covers exist (15.402), while the known four-line family is outside residual by a 15.272/15.588 Max+ witness; only joint residual compatibility remains open. Proposition 15.167's spectral implication is retracted. Global QVAR, principal R1, and the spectral floor are no longer acceptance units.
   - **Residual (ii)** ND (or another Max+-free kill) for **even \(k\ge4p\)**, not only \(k\in[3p+1,4p-2]\). Affine two-level at \(k\ge3p\) is already dead (15.179).
   - **Type I** when Max− is **not** two-level \(\{-1,-3\}\): the 15.169 bad case \(f_e\equiv-1\) on \(\{S=-1\}\) must still force \(\Phi(H)\ge\Phi-2\) (or be empty). Dual-eq empty on the two-level law is not enough.
   - **Lemma D** fully written and checkable: every triple of good lines is a \(k=3\) Max+ (pointwise \(z\), three-line Fourier support, \(\hat z(0)=p\)), and each locked triple spans the 2-plane \(\{x+y+z=0\}\) on its edges (amplitudes must reconcile Fejer characters with convolution). **Satisfied by 15.276 and `A3_PROOF.md`;** census rank-2 at \(p=5,7,11\) is a check, not the general proof.
2. `four_e1_units_closed()` is True **only** by importing the valid bi-tight, residual-(ii), Type-I, and Lemma-D units (no handwritten `return True`). Legacy `e1_closed_general()` is not an acceptance predicate. Aut-Schur **False**, Gsum **False**, pairing **False**, and the open spectral floor stay as they are but are not acceptance gates.
3. `solution.md` Main Theorem (limit) states \(L=\lim\alpha_n=1/2\) **only after** 1–2, and `evidence/share/denseness_path_package.md` § Caveats is empty or records them as closed. Soft-close (sandwich + denseness without E(1)) fails.

## Verification plan
1. **gating:** Run the corrected main-unit dump **twice**. Both runs: required bi-tight levels 2/3 True, residual (ii) True, multi-level Type I True, Lemma D True, `four_e1_units_closed=True`, and \(L\) CLOSED. The spectral floor, `gsum_disj_lb`, Aut-Schur, and pairing may remain False.
2. **gating:** Pytest of 15.272 / 15.270 / 15.249 / 15.167 / 15.236 / 15.237 / 15.170 plus new caveat-unit tests (`pytest -n W`, \(W\) from compute-budget). Capture `{SCRATCH}/pytest_e1.txt`. Exit 0. New units must go False if the identity or the \(k\ge4p\) / multi-level / \(\lambda_{\min}\) argument is wrong.
3. **gating:** Independent cold read of `evidence/share/denseness_path_package.md`. Q1 (the two-level residual-i hinge and required bi-tight levels) must be proved; Q2 (E(1) / \(L=1/2\)) remains blocked until the two multi-level remainders close.
4. **evidence:** Diff of the four units + predicate wiring + package + `solution.md`. Census \(p\le7\) must not be the only reason a “general” flag is True (except the already-finite \(p=5\) Veronese rank 65).

## Non-goals
- Aut-Schur / Jacquet / PSL-span of \(k=3\) \(F\).
- Gsum disj LB; cotangent pairing \(1^\top K^{-1}v\).
- Envelope / reflection / \(K_4\le\mathrm{Wick}_{hi}\) / \(|\mu|\le2/n\).
- Path-C / \(16N\) / Hypothesis H / 15.193 exhaustiveness.
- Re-proving sandwich, denseness Prop 6.1, \(\rho=1\), Johnson Lemma E, cost_D polynomial, \(V_+\) Fourier iso.
- Unflipping \(L\) to “start honest” and then flipping it back. If a caveat dies, record it and switch listed attack — do not cycle the Main Theorem.

## Assumed scope
Repo `/home/nick/quadratic-minmax-limit/`. Hinge already shipped: 15.272 Lemmas B–C, E, F, G Singer on \(F\), 15.207, 15.249, and Lemma D in 15.276 / `A3_PROOF.md`. Binding leftovers: package § Caveats, `HANDOFF.md`, `STATUS.md`. Writeup: `solution.md`. Statement: \(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}|\sum_{i<j}a_{ij}x_ix_j|\) and the target limit is \(L=1/2\).

## Implementation approach
Attack only the two open units: residual (ii) and the Type-I multi-level remainder. QVAR/R1 work is optional mathematics, not a prerequisite for this proof. Wire flags by import after a unit is live and failing-when-wrong. Set `solution.md` \(L=1/2\) only after criteria 1–2. Use compute only behind a predeclared mathematical gate.

## Task checklist
- [x] Replace the invalid 15.167 bi-tight arrow with a valid all-prime obstruction for required levels 2 and 3 (15.720; bi-tight level 4 is a non-load-bearing corollary).
- [ ] Close residual (ii) for even \(k\ge4p\)
- [ ] Close Type I bad case when Max− is multi-level
- [x] Write Lemma D: every triple + 2-plane (15.276; `A3_PROOF.md`)
- [ ] Import those units so `four_e1_units_closed()` is True for the right reasons
- [ ] Package caveats cleared; `solution.md` states \(L=1/2\) only then
- [ ] Two 15.170 dumps + honesty pytest + cold read of the package

## Risks
- Residual (i)/(ii) multi-level leftovers are real. A structure increment that leaves either open fails the goal.
- Historical `e1_closed_general` wiring is not acceptance; the corrected main dictionary is authoritative.
- Do not `/goal clear` to “start fresh.”
