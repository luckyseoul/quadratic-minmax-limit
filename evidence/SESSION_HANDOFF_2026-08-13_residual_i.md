# Session handoff — residual (i) 2026-08-13

## Goal and status

**Goal:** Settle MO 413935 (\(L=\lim\alpha_n\)). Live leftover after this session: residual **(i)** only.

**Status: OPEN — goal NOT met.** Residual (ii) ND is closed. Residual (i) / E1 / L still OPEN. Soft-close forbidden. Do not flip `e1_closed_general` on census \(p\le7\).

| Flag | Value | Source |
|------|-------|--------|
| `residual_ii_a_ND_closed` | **True** | 15.237 |
| `residual_ii_b_ND_closed` | **True** | 15.236 |
| `residual_ii_full_closed` | **True** | 15.179 ∧ 15.237 ∧ 15.236 |
| `deep_s2_freeness_fail_k_ge_3p_ND_closed` | **True** | wires to 15.193 full |
| `gsum_disj_lb_proved_general` | False | 15.170 |
| `residual_i_dual_eq_empty_proved_general` | False | 15.216 |
| `type_I_k_3p_minus_2_closed_general` | False | 15.170 |
| `e1_closed_general` | False | type_I ∧ deep_s2 ∧ bi-tight |

`/goal` hook is the live session plan (`goal/plan.md`), not a repo `.md`.

## Shipped this arc

| Prop | File | Status |
|------|------|--------|
| 15.232–235 | `src/e1_gmin_m4_prop15232.py` … `15235.py` | R̄₄ intersection layers; unsigned dead; residual (i) still OPEN |
| 15.236 | `src/e1_gmin_m4_prop15236.py` | **(ii-b) ND True** (Max− dichotomy + dual-bad empty) |
| 15.237 | `src/e1_gmin_m4_prop15237.py` | **(ii-a) ND True** (L²=L pair-span classification) |
| wiring | `src/e1_gmin_m4_prop15193.py`, `15171.py` | `residual_ii_full` = affine ∧ ii-a ∧ ii-b (not exhaustiveness) |

Tests: `tests/test_prop15236.py`, `tests/test_prop15237.py`, and 15169/15170/15171/15193 updated so deep_s2/ii-a/full expect True; e1 still False.

## Residual-(i) attack state (do not treat as a close)

Preferred hinge: \(|\mu|\le 1/(2p)\) on \(|\kappa|=1\) (⇔ Gsum\(\ge-1/p\) ⇒ 15.176 Farkas). Sufficient: \(|\mu|\le 2/n\); envelope \(|\mu|\le\max(|\mu_{\mathrm{part}}|,|f_4|)\); dual-eq empty; Gsum≽0+master+Fréchet \(\Rightarrow|\mu|\le 2/n\).

**Census (not general):** p=5 max\(\lvert\mu\rvert=3/65<1/10\); p=7 max\(\lvert\mu\rvert=109/2863<2/n\). Dual-eq box-infeasible at p=5 (L∞ \(8/7\)) and p=7 (L∞ \(157/113\)); linearly feasible over \(\mathbb R\). \(\nu=0\) on all \(|\kappa|=1\) at p=5,7.

**Still viable:**
1. Envelope / \(\lvert\mu\rvert\le\max(\lvert\mu_{\mathrm{part}}\rvert,\lvert f_4\rvert)\) (holds at p=5,7; \(j\)-split at extreme \(\varphi\)).
2. Gsum≽0 + master + Fréchet \(\Rightarrow\lvert\mu\rvert\le 2/n\) (p=5 Aut-SOS saturates near \(2/n\); no general dual).
3. Dual-eq Farkas for \(x\ge 0\) (rigid gaps; HiGHS \(y\) not a low-degree function of \(f_e\)).

**Dead / do not re-thrash:** Fréchet alone; \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise (false p=7); \(\lvert\delta\rvert\le\mathrm{room}_\delta\) (false p=5); Gsum PSD near \(\mu_{\mathrm{part}}\) (stays on the wall); 15.237 pair-span list as a Type I kill (Type I *is* a pair-slice; \(\ker F_+\) large); halfspace \(G\)-orbit subframe (tight, 4-point up to 0.6); \(d_H\ge p+1\) / \(N_3=(p^2-5)/4\) (true, too weak for \(K_4\) or L∞); unsigned Per/\(\sum\lvert\mathrm{per}\rvert\); Aut-line dim\(\le 1\).

Scratch (ephemeral): `/tmp/grok-goal-cc538b97808e/implementer/` (`resi_*.py`, `maxplus_p7.npy`). Caches: `/tmp/maxplus_p5.npy`, `/tmp/maxminus_p5.npy`, `/tmp/maxplus_p7.npy`.

## Next concrete steps

1. Prove one Max+-free general-\(p\) residual-(i) hinge (envelope, Aut-SOS dual, or dual-eq Farkas).
2. Wire `gsum_disj_lb_proved_general` / `type_I` / `e1_closed_general` via real imports; set Main Theorem \(L=\tfrac12\).
3. Run plan verification: predicate dump, hinge `main` twice, hinge pytest.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee`
