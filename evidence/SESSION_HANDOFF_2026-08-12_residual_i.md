# Session handoff — residual (i) 2026-08-12

## Goal and status

**Goal:** Prove Max+-free residual-(i) hinge for all primes \(p\ge5\) (prefer \(\max_{|\kappa|=1}|\mu_4|\le1/(2p)\); else m4₂ / K₄ / free-e / Path C), then wire `residual_i` / `type_I` via real imports.

**Status: OPEN — goal NOT met.** Incomplete math. Soft-close forbidden. Predicates all **False**. AI-test skipped.

| Flag | Value |
|------|-------|
| `gsum_disj_lb_proved_general` | False |
| `residual_i_dual_eq_empty_proved_general` | False |
| `type_I_k_3p_minus_2_closed_general` | False |
| `e1_closed_general` | False |

## Shipped this arc (structure only — no predicate flip)

| Prop | File | Proved Max+-free | Open |
|------|------|------------------|------|
| 15.218 | `src/e1_gmin_m4_prop15218.py` | Max+⊥Max−; Gsum=2m₁μ₄; Tr(Gsum²); row+ones dual if \(G_{\min}\ge-1/p\) | \|μ\|≤1/(2p) |
| 15.219 | `...15219.py` | R-matrix / Tr(RC) structure | φ-contraction general |
| 15.220 | `...15220.py` | ‖m₄‖₂²=(K₄−9n²+10n)/24 ⇔ K₄≤n(15n−22) | K₄ thr |
| 15.221 | `...15221.py` | (4pI−T)χ=D; wedge ff=0; ‖D‖₂² closed | δ-room |
| 15.222 | `...15222.py` | ⟨D,TD⟩=2p(p−1)(p+1)(p²+1)(p⁴−11p²+40); CS ‖δ‖₂²≤k₂/\|Max+\| | m0, \|Max+\| |
| 15.223 | `...15223.py` | gap_wick closed; Path C ⇒ Wick ⇒ residual (i); gap>room_hyp/24 for p≥5 | Path C / CS general |
| 15.224 | `...15224.py` | T²φ=4(n+2)(φ−2κ); μ_part solves master always | δ∈E_{±4p} for \|μ\| |

Tests: `tests/test_prop15218.py` … `test_prop15224.py`. Evidence JSON alongside.

Scratch (ephemeral): `/tmp/grok-goal-bd1830373544/implementer/e1_l_block.md`, predicate dumps.

## Settled decisions

- Soft-close forbidden; census \(p\le7\) does not flip general predicates.
- Preferred hinge still \(|\mu|\le1/(2p)\) on \(|\kappa|=1\) (⇔ Gsum≥−1/p ⇒ 15.176 Farkas).
- Viable alt targets: Path C \(\delta^2\le\mathrm{room}_{hyp}/24\); \(k_2/|Max+|\le\mathrm{gap}_{wick}\); \(K_4\le\mathrm{Wick}_{hi}\); free-e on true ker.
- \(\mu_{\mathrm{part}}\) majorant ≤1/(2p) proved; \(|\mu|\le|\mu_{\mathrm{part}}|\) **false** pointwise at p=7; viable census target \(|\mu|\le2/n\).
- f4 envelope \(|f4|\le2/n\) proved; \(|\mu|\le|f4|\) fails at p=7.

## Census anchors (not general)

- p=5: max\|μ\|=3/65; δ²=room_hyp/24=1536/65; k₂=6912, \|Max+\|=260, k₂/N≤gap_wick.
- p=7: max\|μ\|=109/2863<2/n; \|Max+\|=11452; CS and Path C hold numerically.

## Open hinge (any one closes residual i)

1. Control \(\delta\in E_{\pm4p}\): \(|\mu|\le1/(2p)\) or \(\le2/n\) on \(|\kappa|=1\) for all \(p\ge5\)
2. Path C: \(\delta^2\le\mathrm{room}_{hyp}/24\) for all \(p\ge5\)
3. Closed \(m_0=\|\chi_{\mathrm{part}}\|_2^2\) and/or \(|Max+|\) so \(k_2/|Max+|\le\mathrm{gap}_{wick}\)
4. \(K_4\le\mathrm{Wick}_{hi}\) Max+-free
5. G₊≻0 (ker=sc) + free-e_sc / cost_D≤8/p

## Next concrete steps

1. User supplies a proof idea for one hinge above, **or** accept residual-(i) open.
2. If proved: wire `gsum_disj_lb_proved_general` / `residual_i_dual_eq_empty_proved_general` / `type_I_k_3p_minus_2_closed_general` only via **real imports** from hinge modules; flip STATUS/package; AI-test ≥2 cold.
3. Residual (ii) full (15.193) still separate open after residual (i).

## Dead / weak paths (do not re-thrash)

- Soft-close; census-only general claims
- \(|\mu|\le|\mu_{\mathrm{part}}|\) pointwise; f4 as general majorant
- \(G_{\min}\ge-2\) (trivial \|\mu\|≤1) dual cost ≫2−α
- D supported only on \(\pm2\sqrt{p^2+15}\); pairwise LP / 4×4 Loewner alone
- R_ke⇒R (convex combo fatal); Max+ classical association scheme

## Suggested skills

**Always (every residual-i / denseness session):**
- **agent-cost-optimization** — short hinge state; no Max± dumps / prop thrash
- **graph-engineered-completion** — nodes = open hinges (δ, |μ|, K₄, m0, ker=sc); edges = 15.218–224; attack one node

**Verification / honesty:**
- **goal-verifier** — goal restatement + criteria + evidence before any “done”
- **verification-before-completion** — fresh pytest / predicate dump before CLOSE claims
- **handoff** — compact next-agent doc (must list always-on skills first)
- **session-handoff-packager** — package session summary to disk

**Math quality / critique:**
- **scientific-critique** — cold 9-dim critique of claims / soft-close risk
- **grill-me** — harden a proof idea one question at a time before formalizing
- **self-refine-loop** — polish a proof draft only when explicitly invoked (cap 5)

**Research / literature:**
- **research** — hybrid research router
- **arxiv** — Paley / conference / association-scheme papers
- **litreview** — formal lit orientation if needed

**Compute / infra:**
- **use-available-compute** — ProcessPool/GPU for p=7 Max+ / sparse T; no serial thrash

**MCP (not a skill file; use when stuck or verifying):**
- **openai-referee** — `health`, `suggest_direction`, `math_review`, `falsify`, `verify`, `deep_review`

**Do not soft-close.** Domain math only until a hinge is proved.

## Pointers

- Plan: session `goal/plan.md` (prove residual-i hinge then wire)
- STATUS: `STATUS.md` · package: `evidence/share/denseness_path_package.md`
- Prior: `evidence/SESSION_HANDOFF_2026-08-06_residual_i.md`
- Props: `src/e1_gmin_m4_prop15218.py`–`15224.py`
