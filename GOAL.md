# Plan: Settle \(L=\lim\alpha_n=1/2\) by closing the remaining three of four units

Use this file as the `/goal` body. Example:

```
/goal complete GOAL.md in /home/nick/quadratic-minmax-limit
```

Do not attach a `--budget` flag. Do not `/goal clear` to start fresh.

## Goal kind
code-change

## Acceptance criteria
1. For every prime \(p\ge5\), E(1) is proved Max+-free on Paley \(n=p^2+1\). That means all four of the following, not a subset:
   - **Spectral floor \(\lambda_{\min}(\Phi)\ge6\)** on \(Z\): prove both current conjuncts, global mixed-\(k\) QVAR and principal R1. The argument must not treat \(G_{u,\mathrm{disj}}\) as a Gram (it has negative eigenvalues), split QVAR into false per-stratum claims, or substitute the finite `p=11` R1 census. Then 15.167 majorization empties bi-tight for all such \(p\).
   - **Residual (ii)** ND (or another Max+-free kill) for **even \(k\ge4p\)**, not only \(k\in[3p+1,4p-2]\). Affine two-level at \(k\ge3p\) is already dead (15.179).
   - **Type I** when Max− is **not** two-level \(\{-1,-3\}\): the 15.169 bad case \(f_e\equiv-1\) on \(\{S=-1\}\) must still force \(\Phi(H)\ge\Phi-2\) (or be empty). Dual-eq empty on the two-level law is not enough.
   - **Lemma D** fully written and checkable: every triple of good lines is a \(k=3\) Max+ (pointwise \(z\), three-line Fourier support, \(\hat z(0)=p\)), and each locked triple spans the 2-plane \(\{x+y+z=0\}\) on its edges (amplitudes must reconcile Fejer characters with convolution). **Satisfied by 15.276 and `A3_PROOF.md`;** census rank-2 at \(p=5,7,11\) is a check, not the general proof.
2. `four_e1_units_closed()` is True **only** by importing those four units (no handwritten `return True`). Legacy `e1_closed_general()` is not an acceptance predicate. 15.272 Johnson same-line, \(V_+\) Fourier iso, cost_D polynomial, Aut-Schur **False**, Gsum **False**, pairing **False** stay as they are.
3. `solution.md` Main Theorem (limit) states \(L=\lim\alpha_n=1/2\) **only after** 1–2, and `evidence/share/denseness_path_package.md` § Caveats is empty or records them as closed. Soft-close (sandwich + denseness without E(1)) fails.

## Verification plan
1. **gating:** From `/home/nick/quadratic-minmax-limit` with `PYTHONPATH=src`, run the main four-unit dump **twice**. Both runs: global QVAR True, principal R1 True, residual (ii) True, multi-level Type I True, Lemma D True, `four_e1_units_closed=True`, \(L\) CLOSED, `gsum_disj_lb` False, Aut-Schur/pairing False. Save `{SCRATCH}/predicate_dump_{1,2}.txt`.
2. **gating:** Pytest of 15.272 / 15.270 / 15.249 / 15.167 / 15.236 / 15.237 / 15.170 plus new caveat-unit tests (`pytest -n W`, \(W\) from compute-budget). Capture `{SCRATCH}/pytest_e1.txt`. Exit 0. New units must go False if the identity or the \(k\ge4p\) / multi-level / \(\lambda_{\min}\) argument is wrong.
3. **gating:** Independent cold read of `evidence/share/denseness_path_package.md` **only** (`scientific-critique` + referee). Q1 (the two-level residual-i hinge) must be proved; Q2 (E(1) / \(L=1/2\)) must remain blocked until the other three units close. A missing map, false identity, or soft-close is a fail.
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
Attack only the three open units (spectral floor, residual (ii), then the Type-I remainder unless strong R1 closes it). Prefer one Max+-free algebraic unit per caveat. Wire flags by import after the unit is live and failing-when-wrong. Update the package as a stand-alone argument. Set `solution.md` \(L=1/2\) only after criteria 1–2. Use available compute only behind a predeclared mathematical gate. Referee / scientific-critique when a unit is claimed closed.

## Task checklist
- [ ] Prove GLOBAL mixed-`k` QVAR for all primes \(p\ge5\); do not split into the false per-stratum `k>=7` claim.
- [ ] Prove principal R1 uniformly (the `p=11` census is finite; scalar and broad-channel cones are insufficient).
- [ ] Close residual (ii) for even \(k\ge4p\)
- [ ] Close Type I bad case when Max− is multi-level
- [x] Write Lemma D: every triple + 2-plane (15.276; `A3_PROOF.md`)
- [ ] Import those units so `four_e1_units_closed()` is True for the right reasons
- [ ] Package caveats cleared; `solution.md` states \(L=1/2\) only then
- [ ] Two 15.170 dumps + honesty pytest + cold read of the package

## Risks
- Residual (i)/(ii) leftovers are real. A structure increment that leaves any of the four open fails the goal.
- `e1_closed_general` currently ANDs bi-tight only at \(p=5\) and treats 15.236/237 as “full” residual (ii). That wiring is not acceptance.
- Do not `/goal clear` to “start fresh.”
