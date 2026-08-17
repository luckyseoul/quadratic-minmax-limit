# Session handoff (2026-08-17 afternoon, leftover 1)

**Repo:** `/home/nick/quadratic-minmax-limit` · `luckyseoul/quadratic-minmax-limit` `main` @ `e0c7284`  
**Binding:** `GOAL.md`. Do **not** `/goal clear`. Do **not** unflip `e1` / `L`. Do **not** flip Aut-Schur / Gsum / pairing. House **15.x**.

Scratch leftover graph (not in git): `/tmp/grok-goal-ea1b2cdbd197/implementer/leftover_graph.md`

## Honest status

Leftover 1 (\(Q_\tau\) / \(\lambda_{\min}(\Phi)\ge6\)) is still **OPEN**. Residual (ii) \(k\ge4p\) and Type I multi-level stay **OPEN**. Lemma D flags True. Live `e1` is the old AND. \(L=1/2\) is **not** settled. \(\phi_F\) was **not** imported.

## Last ~4 hours (15.459–15.465)

| unit | commit | what it is | leftover 1? |
|------|--------|------------|-------------|
| 15.459 | `2882abc` | \(k\le3\) type-count is \(3/2\) at \(p=13\); no \((\max k,n_0)\) slice in \([551,617]\) | still open |
| 15.460 | `650afc0` | \(k(R_{\mathrm{aff}})=(p^2-1)/24\); CRRR count \(C(r,3)u\) = 100,1176,70980; \(c_{\mathrm{CRRR}}=1,4,35\) | 35 misses window |
| 15.461 | `57e9c68` | p=7 remainder shapes \(J_6\times J_2^{r-1}\), \(J_4^2\times J_2^{r-2}\), aff+J2\(\times J_2^{r-1}\) do not extend | \(C(r,3)+15=50\) misses |
| 15.462 | `531ec1a` | no single mixed occupancy + \(C(r,3)\) in the window (nearest 455 / 665) | dead family |
| 15.463 | `f72bfcb` | no two-occupancy type-support pair recovers \(c=1,19\) and the \(p=13\) window | \(\{J_6\}+\{J_2,J_4,R\}\) is \(7001/2\) |
| 15.464 | `03396f1` / `d851311` | ns \(E[N]=\mu_-\) in the **uniform** \(\mu\)-subset model; live Paley-\(N^*\) is \(1496/409\), **not** S0 \(Q_n=1440/409\); mean-field \(L\) misses \(Q_{++}\) | OpenAI BLOCK on missing average, fixed |
| 15.465 | `9f9f3ab` / `e0c7284` | ratio Paley-\(N^*/Q_n=1+4p(p-5)/[(p+1)(p-2)(p-1)^2]\) at \(p=5,7\); drop \((p-5)\) \(\to77/72\neq1\) | **2-point form**, not a general \(Q_\tau\) |

Live Paley \(Q\) on 15.290 types: \(Q_{++}=48/13,1544/409\); Paley-\(N^*=32/13,1496/409\); \(Q_{--N1}=1376/409\) at \(p=7\). S0 dual \(Q_n=32/13,1440/409\).

## Referees (Claude first, then OpenAI)

- 15.464 A: Claude PASS-WITH-NOTE; OpenAI **BLOCK** (fixed: uniform model + \(\{0,1\}\) witness).
- 15.465: Claude PASS-WITH-NOTE (fail is literal drop of \(p-5\); Jacobi \(4\) is naming); OpenAI **PASS**.
- `suggest_direction` (stuck leftover 1): name a **character-sum / type-indicator Fourier / Hasse–Davenport** formula for \(Q_\tau\), not another occupancy census or 2-point fit. Do not import \(\phi_F\) from 15.465.

## Do not do

- Occupancy / type-count census (15.458–463 exhausted).
- Treat 15.465 as a general \(Q_\tau\) law or import \(\phi_F\) from it.
- Kill `netcount_opt` PID **712954** (tmux `netcount`, 84 OMP, ckpt `/tmp/grok-goal-ea1b2cdbd197/implementer/ckpt_p11/`). Kill only that exact PID, never `pkill netcount`.
- Flip True flags. Unflip `e1`. Soft-close \(L\). `/goal clear`.
- Public “prize” / “Paata” / “ping”. `pbss-goal-verifier` / `perry-beurling-rh-closeout`.

## Live compute

- **712954** `./netcount_opt -p 11 -j 84 -no-mrv -ckpt …/ckpt_p11` — leave it. Historically 0 shards.
- Budget while it owns the box: `recommended_workers=1`. V100 idle. GPU raw kernels OK (CuPy CUB/`count_nonzero` is broken on this V100; use `RawKernel`).
- Caches: `/tmp/maxplus_p{5,7}.npy`.

## Next

1. Structural Gauss/Jacobi (or type-indicator Fourier) name of \(Q_{++}\) on 15.290 types, **forced by character sums**, not a \(p=5,7\) rational. Fail-when-wrong must break it.
2. Import `phi_F` only after that pin is real, then \(\langle\delta,\psi\rangle\le2\).
3. Residual (ii) \(k\ge4p\) and Type I after the floor.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `use-available-compute` · `claude-referee` (primary, opus) · `openai-referee` (secondary, **after** Claude, same slot) · `handoff` · `session-handoff-packager` · `verification-before-completion` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview`

## Pointers

- Units: `src/e1_gmin_m4_prop15453.py` … `15465.py` + `tests/test_prop1545*.py` `test_prop1546*.py`
- Morning handoff: `evidence/SESSION_HANDOFF_2026-08-17_leftover.md`
- Plan: session `.../019f9af7-3128-71c1-984e-2a7102bec72d/goal/plan.md`
