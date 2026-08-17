# Session handoff (2026-08-17, leftover campaign)

**Repo:** `/home/nick/quadratic-minmax-limit` · `luckyseoul/quadratic-minmax-limit` `main`  
**Binding goal:** `GOAL.md`. Do **not** `/goal clear`. Do **not** unflip `e1` / `L`. Do **not** flip Aut-Schur / Gsum / pairing. Do **not** run `pbss-goal-verifier` or `perry-beurling-rh-closeout`. House numbering: stay **15.x** unless the user calls the campaign failed (16.x) or \(L=1/2\) is honest (17.x).

Scratch leftover graph: `/tmp/grok-goal-ea1b2cdbd197/implementer/leftover_graph.md` (deleted when the goal ends; not in git).

## Goal

Prove \(L=\lim\alpha_n=1/2\) by closing four leftovers so E(1) is Max+-free for all primes \(p\ge5\). Live `e1_closed_general` is still True by the **old AND**. That wiring is not acceptance.

## Live flags (do not change unless a unit actually imports)

| flag | value | note |
|------|-------|------|
| `phi_F_ge_6_proved_general` | False | floor / \(Q_\tau\) open |
| `e1_closed_general` | True | dishonest old AND |
| Aut-Schur / Gsum / pairing | False | stay False |
| `residual_ii_k_eq_4p_empty` | False | interior 4-level + 5+ open |
| Type I multi-level | False | |

## Shipped this continue (15.375–15.382)

- **15.375–15.378:** \(Q_{--}^{\mathrm{form}}=8(A-k)/D\); 1D is the \(Q_{--}=0\) slice; all eight p=7 atoms constructed; interpolants die at p=19. Do **not** import \(\phi_F\) from \(24/p\), \(Q_{\mathrm{form}}\), or the p=7 dictionary.
- **15.379:** 4-level first-moment \(E[S^2]=16+8(a+e)\); lattice min \(20+12/p\) with \(a\) locked to \((p+1)/(2p)\).
- **15.380:** Paley \(E[f_e]=1/p\); adj \(\{\pm1/p\}\) mean 0; disjoint mean \(1/(n-3)\); no regular 4p-graph (\(d=8p/(p^2+1)\notin\mathbb Z\)).
- **15.381:** Max− \(E[f_e]=-1/p\). \(a=\mathrm{thr}\) iff minus-slice. Both 4-level J-corners empty for every \(p\ge5\) (high corner unclassified mass; low corner unclassified except p=7 3-equal, which dies by \(e\in T\) linear independence + certified both-signs \(\max|E[f_e|U]|=259/409<1\)). OpenAI **PASS** after two BLOCKs that were written into the unit. Claude **PASS-WITH-NOTE** (do not branch).
- **15.382:** Interior 4-level needs both \(S=-4\) and \(S=-6\) on \(\{f_e=+1\}\) (distinct mod 4). Pure-pair Aut_e double-star is a single residue class, hence not interior. Minus-slice first moment forces \(|G\cap\boxtimes|=p+1\).

Tests: `tests/test_prop15381.py` 8/8, `tests/test_prop15382.py` 6/6.

## What is still open (the four leftovers)

1. **Name \(Q_\tau\) in \(p\)** (Gauss/Jacobi or equivalent) on 15.290 types matching \(48/13\) and \(1544/409\). Best named form: \(Q_{++}=8A/D\) with \(D_{\mathrm{form}}=2A-3(p+1)/2+C(p,(p-1)/2)\) equal to live \(D\) at p=5,7 only. Do **not** import \(\phi_F\) from the form interval or a p=11 census.
2. **Residual (ii) even \(k\ge4p\).** Leftover is multi-level Max−, \(\max=-2\), \(f_e\equiv-1\) on \(U_{-2}\). Corners empty. Interior 4-level with **singles** or **far Aut_e edges**, and 5+ even levels, stay open.
3. **Type I** multi-level (integer tails exist; maj not a bound at p=7).
4. **Lemma D** flags True; do not unflip / cascade.

## Live compute (do not duplicate)

- **p=11 6-net** `cpu_D11_countnet.py`: 86 Numba workers, genuine backtrack, ~2h+, p=5=130 and p=7=5726 certified. **Zero p=11 shards finished** (`ex.map` only prints per shard of ~5.4 masks). \(N_{\mathrm{pred}}=244332\). Do **not** launch another 86-worker job. Kill only if the owner wants the cores back.
- **GPU** V100 idle except leftover hunters. `gpu_dstar_singles.py`: p=5 sample **0/4000** double-stars with \(S\equiv-2\) on \(U\) (no interior hit). p=7 sampling in flight.
- Budget while 6-net owns the box: `recommended_workers=1`. GPU-first for Grams / Max± scores.

## Next concrete steps

1. **Do not wait on \(D_{\mathrm{form}}=D\) at p=11.** Attack a leftover that can flip a live flag:
   - Interior 4-level with singles: prove \(\mathrm{Var}(S\mid U)>0\) (p=5 sample already never has \(S\equiv-2\) on \(U\)), **or**
   - Far Aut_e leftover, **or**
   - 5+ even levels, **or**
   - Name \(Q_\tau\) independently (not a 3-point interpolant).
2. Ship **15.383** with fail-when-wrong. Do not flip `residual_ii_k_eq_4p_empty` / `phi_F` / `e1` / `L` until a real unit imports.
3. Only after all four leftovers close: import units, clear package caveats, set `solution.md` \(L=1/2\). Verification: two 15.170 dumps, honesty pytest `-n W`, both referees on the package.

## Do not reopen

Invented \(N_+\); \(\phi_F\) from p=5/7 mixtures; \(V_{2\mathrm{orb}}\) unproved; dens interpolants; Lasserre=5; Poincaré=3.85; dead \(c(p)\) split; QR0⊙NC orbit \(Q\); p=7 24-table; 15.341 \(S^{(k)}\); 7260-orbit \(Q_{++}=416/165\); open NC sample; p=5/7 \(n_4\) type table; triangle orbit; p=13 SA sample; T⊙NC; p=7 type dictionary as a general name.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `use-available-compute` · `openai-referee` · `claude-referee` (only if OpenAI is out of tokens, unless the user asks for both) · `handoff` · `session-handoff-packager` · `verification-before-completion` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview`

Subagents: OK for **minor ≤5-minute** independent facts (one owner per leftover). Spawn contract: one fact owed, scratch path, do not flip flags, compute-budget one-liner. GPU if dense fit.

## Pointers

- Goal / plan: `GOAL.md`, session `.../019f9af7-3128-71c1-984e-2a7102bec72d/goal/plan.md`
- Units: `src/e1_gmin_m4_prop15375.py` … `15382.py` + matching `tests/test_prop153xx.py`
- Evidence: `evidence/e1_gmin_m4_prop15375.json` … `15382.json`
- Caches: `/tmp/maxplus_p{5,7}.npy`, `/tmp/maxminus_p{5,7}.npy`
- Package caveats: `evidence/share/denseness_path_package.md`
- Prior session: `evidence/SESSION_HANDOFF_2026-08-16_leftover.md`
