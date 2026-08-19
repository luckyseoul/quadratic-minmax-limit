# Fable 5 brief — settle E(1) for MO 413935

Ready-to-use prompt + config for handing this repo's open problem to Claude Fable 5.
Everything below is verified against the repo at `1896c60` unless marked OPEN.

Written deliberately **un-prescriptive**: Fable 5 degrades when given step-by-step
scaffolding written for earlier models. It gets the goal, the verified facts, the
dead ends, and the acceptance bar — it plans the attack itself. Do not re-add a
mandated attack order.

---

## Config

```python
client.beta.messages.create(
    model="claude-fable-5",
    max_tokens=128000,                        # start 64000; raise if truncating
    output_config={"effort": "xhigh"},        # long-horizon => high/xhigh
    betas=["server-side-fallback-2026-06-01"],
    fallbacks=[{"model": "claude-opus-4-8"}],
    messages=[{"role": "user", "content": PROMPT}],
)
```

Hard rules for this model (each returns 400 otherwise):

- **Omit `thinking` entirely.** Thinking is always on. `{"type":"disabled"}` and
  `{"type":"enabled","budget_tokens":N}` both error.
- **No `temperature` / `top_p` / `top_k`.**
- **No assistant prefill** (no trailing `role: "assistant"` message).
- 1M context, 128K max output. **Stream** — single requests on hard tasks run
  many minutes; plan timeouts and progress UX accordingly.
- Raw chain of thought is never returned. `display: "summarized"` gives a summary;
  the default is empty thinking text.

`fallbacks` is near-irrelevant here (the classifiers target cyber/bio; this is
combinatorics) but costs nothing — drop it if you prefer.

Give the whole spec in **one first turn**. Do not build it up interactively —
that measurably reduces both quality and token efficiency on this model.

---

## The prompt

> You are attacking an open problem in combinatorics. Repo:
> `luckyseoul/quadratic-minmax-limit`, branch `prop15586-maxplus-gram-reduction`.
> Read `HANDOFF.md`, `LONG_HORIZON_GOAL.md`, and
> `evidence/SESSION_HANDOFF_2026-08-19_gram_reduction.md` first.
>
> **The problem.** MathOverflow 413935:
> `alpha_n = n^{-3/2} min_{a_ij = +-1} max_{x = +-1} |sum_{i<j} a_ij x_i x_j|`.
> Does `L = lim alpha_n` exist, and what is it? Proved so far: the sandwich
> `1/pi <= liminf <= limsup <= 1/2`; denseness (the limit exists globally iff it
> exists along Paley orders `n = p^2+1`); and `rho = 1` on that Paley family.
> What remains is **E(1)**: the Paley conference matrix `C` is a `Phi`-minimizer
> for every prime `p >= 5`. E(1) reduces to exactly three lemmas, each currently a
> hardcoded `False`:
>
> | # | Predicate | Statement |
> |---|---|---|
> | 1 | `phi_F_ge_6_proved_general` | `lambda_min(Phi) >= 6` on `Z` |
> | 2 | `multilevel_ND_k_ge_4p_proved` | multi-level Max- ND at even `k >= 4p` |
> | 3 | `type_I_aut_e_3AB_positive_general` | `3A+B > 0` on every Aut_e far class |
>
> Everything around all three is proved: leftover 2's theorems A-L all return
> True, and 39/39 lemmas in leftover 3 return True at p = 5,7,11,13,17,19,23.
>
> **Your goal is to close one of the three, or to kill its route by a general
> counter-mechanism.** A new proposition that does not flip a flag or kill a path
> is not progress — the repo has ~500 of those already and the flags never moved.
>
> **Verified structure you can build on** (general `p` unless noted;
> `n = p^2+1`, `P = (I+C/p)/2`, `Max+ = {y in {+-1}^n : Cy = py}`, `N = |Max+|`,
> `Z = {B sym : CB = pB, diag B = 0}`, `Phi(B) = E_y[(y^T B y)^2]`):
>
> - Every `y` in Max+ satisfies `Py = y` and `y_i^2 = 1`. Hence for `pi_i = P e_i`,
>   `<y y^T, pi_i pi_i^T> = y_i^2 = 1` for all `i, y`, so `proj_Zperp(y y^T) = R`
>   is the **same** element for every `y`, with `||R||^2 = 2n`.
> - Therefore `spec(Phi) = nonzero spec(Ghat / N)` where
>   `Ghat_ab = <y_a, y_b>^2 - 2n`. `Ghat` is an **integer** matrix and a genuine
>   Gram matrix, PSD by construction — the package's "not a Gram" objection to
>   `G_{u,disj}` does not apply to it.
> - `dim Z = n(n-6)/8`; `tr Phi = n(n-2)`; `tr K = -4n` for `K = 8I - Phi`;
>   `dim span{1, y_i y_j} = n(n-6)/8 + 1`.
> - `E[y y^T] = I + C/p = 2P`, so the Wick value is exactly `8||B||^2` and
>   **leftover 1 is exactly `lambda_max(K) <= 2`**.
> - Exact spectra: p=5 `{80/13 (mult 26), 144/13 (26), 176/13 (13)}`;
>   p=7 `{3072/409 (50), 3360/409 (100), 3648/409 (50), 4032/409 (50),
>   4320/409 (25)}`. Bottom multiplicity is exactly `n`, top exactly `n/2`.
> - `lambda_max(K) = 48/n` **exactly** at p=5, and `200/409 < 48/50` at p=7. So
>   the floor binds only at the smallest prime — p=5 is a finite check and by p=7
>   there is 4x slack. **A crude bound closes leftover 1 for p >= 7.**
> - Leftover 3's constant is `mu = max_{|kappa(S)|=1} |E_y[y_i y_j y_k y_l]|` with
>   `kappa(S) = C_ij C_kl + C_ik C_jl + C_il C_jk` in `{-3,-1,1,3}`: `3/65` at
>   p=5, `109/2863` at p=7. Targets `L = (p-2)/(2p^2)` and
>   `T = (p-2)/(p(2p-1))`, verified p = 5..43, with `T > L` always (which is why
>   `|mu| <= |T|` cannot close it). Margin `mu/L` = 0.7692, 0.7462 — again ~25%
>   slack, so a crude bound suffices here too.
> - `max|m4|` over **all** four-sets: `21/65` (p=5), `327/2863` (p=7).
> - Denominators throughout are `p*D = N/4`: leftovers 1 and 3 are moments of one
>   tensor.
>
> **The wall.** All three are statistics of Max+, and Max+ is enumerable only for
> `p <= 7`. At p=11 the nullity is 61, so an exhaustive sweep is `2^61 ~ 2.3e18`;
> the structured families give `n_1d = 2772` and `n_{k=3} = 24200` while `n_full`
> is exactly the unclassified family. Getting Max+ moments at general `p` is
> plausibly the single underlying problem. The reduction above narrows leftover
> 1's share of it to 2-point data (`Ghat`); leftover 3 still needs genuine 4-point
> moments on `|kappa|=1` four-sets.
>
> **Do not re-run these — they are settled negatives:**
> - Max+ inner-product classes are **not** an association scheme (max within-class
>   variance of `A_i A_j` is 113 at p=5, 5369 at p=7), so `Ghat`'s spectrum is not
>   determined by the inner-product distribution alone.
> - The 15.237 C 0-1 pair-span classification **survives** an exhaustive test. A
>   0-1 function is in the pair-span iff `q_B(y) = y^T B y` takes <= 2 values on
>   Max+. Over all 5,668,650 support-3 sets at p=5: 3575 two-valued = 2600
>   triangles + 975 extras, every extra spanning 6 vertices (three disjoint edges
>   with degree-6 monomial identically +1) — yet every observed mass is already in
>   `classified_01_pairspan_masses`. No non-triangle extras at p=7 in 598,510
>   samples.
> - `evidence/E1_FAILURE_GRAPH.md` lists failure modes F1-F22 and a long
>   residual-(i) dead-mechanism list. Treat it as evidence about the terrain, not
>   as instructions about method.
> - Chasing an exact Gauss/Jacobi `p`-formula for the ensemble `Q(r)` consumed
>   ~200 commits without moving a flag. `D = |H_+|/(2p)` is 13 and 409 at p=5,7 —
>   not polynomial in `p`. Note the eigenvalues of `Ghat` **are** plain integers
>   (1600/2880/3520 and 86016/.../120960), so if you want a `p`-formula, seek one
>   for `eig(Ghat)`, not for `lambda`.
>
> **Acceptance.** A leftover is closed only when its predicate returns True via a
> real import from a unit that goes False if the argument is wrong — never a
> handwritten `return True`, never a census at `p <= 7` standing in for general
> `p`. Ship fail-when-wrong tests alongside (see `tests/test_prop15586.py` for the
> pattern: perturb the constant, the closed form, and the coefficient, and assert
> each perturbation breaks the match). Soft-close — sandwich + denseness + rho=1
> without E(1) — is not a proof and is explicitly forbidden. If you kill a route
> instead, say so plainly and name the replacement leftover.
>
> **Working notes.** Write what you learn to `evidence/` as you go: one lesson per
> file, a one-line summary at the top, corrections and confirmed approaches alike,
> with the reason each mattered. Update an existing note rather than duplicating
> it; delete notes that turn out to be wrong. Consult them in later sessions.
>
> Delegate independent subtasks to sub-agents and keep working while they run;
> intervene if one goes off track or is missing context. The machine has 88 cores
> and a Tesla V100 — use them. **Note CUDA 13 dropped sm_70, so GPU JIT is broken
> for that V100:** `numba.cuda` fails (`libnvvm: -arch=compute_70 is an
> unsupported option`) and any CuPy path needing NVRTC (axis reductions,
> `cp.unique`, `.sum()` on a bool array, RawKernel) fails compiling
> `cuda_fp4.hpp`. Precompiled cuBLAS, cuSOLVER, elementwise ops and boolean
> masking work; substitute a matvec against a ones-vector for axis reductions.
>
> Before reporting progress, audit each claim against a tool result from this
> session. Only report work you can point to evidence for; if something is not yet
> verified, say so. If tests fail, say so with the output; if a step was skipped,
> say that.
>
> You are operating autonomously and cannot get answers mid-task, so do not ask
> permission for reversible work that follows from this brief — proceed. Before
> ending a turn, check your last paragraph: if it is a plan, a question, or a
> promise about work you have not done, do that work now instead. End only when
> the task is complete or you are blocked on something only a human can supply.
>
> Do not publish, post, or send anything anywhere. Produce the result and the
> writeup; the human decides what goes public.
>
> Lead your final summary with the outcome — what happened, or what you found —
> in one sentence, before any supporting detail. Write it for someone who did not
> watch the run: complete sentences, terms spelled out, no arrow chains or
> shorthand you coined along the way.

---

## Notes for the operator

- **Do not show it a remaining-token countdown.** Surfacing context-budget counts
  causes premature wrap-up.
- **Do not re-prescribe.** If output quality looks off, the fix is usually to
  remove scaffolding, not add it.
- Effort is a dial worth sweeping. `xhigh` is the recommended start for this;
  `high` is a reasonable step down if runs are long without gaining ground.
- If it flips a flag, verify independently before believing it: run the flipped
  predicate, run its fail-when-wrong tests, and confirm the unit goes False under
  perturbation. This repo's entire discipline is that a True flag must be earned.
