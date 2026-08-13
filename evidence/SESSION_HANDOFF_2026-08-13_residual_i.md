# Session handoff — residual (i) 2026-08-13 (evening)

## Goal and status

**Goal:** Settle MO 413935 (\(L=\lim\alpha_n\)). Binding leftover: residual **(i)** only.

**Status: OPEN — goal NOT met.** Residual (ii) ND is closed. Residual (i) / E(1) / \(L\) still OPEN. Soft-close forbidden. Do not flip `e1_closed_general` on census \(p\le7\). Do not ship a 15.xxx whose predicates stay False.

HEAD at last code close: `4d89353` (15.232–237). This document records the **post-4d89353 residual-(i) attack** (no new prop).

| Flag | Value | Source |
|------|-------|--------|
| `residual_ii_a_ND_closed` | **True** | 15.237 |
| `residual_ii_b_ND_closed` | **True** | 15.236 |
| `residual_ii_full_closed` | **True** | 15.179 ∧ 15.237 ∧ 15.236 |
| `deep_s2_freeness_fail_k_ge_3p_ND_closed` | **True** | wires to 15.193 full |
| `gsum_disj_lb_proved_general` | False | 15.170 (hardcoded) |
| `residual_i_dual_eq_empty_proved_general` | False | 15.216 |
| `type_I_k_3p_minus_2_closed_general` | False | 15.170 |
| `e1_closed_general` | False | type_I ∧ deep_s2 ∧ bi-tight |

`/goal` hook is the live session plan (`goal/plan.md` under the session dir), **not** a repo `.md`.

## Shipped this arc (already on `origin/main`)

| Prop | File | Status |
|------|------|--------|
| 15.232–235 | `src/e1_gmin_m4_prop15232.py` … `15235.py` | R̄₄ intersection layers; unsigned dead; residual (i) still OPEN |
| 15.236 | `src/e1_gmin_m4_prop15236.py` | **(ii-b) ND True** (Max− dichotomy + dual-bad empty) |
| 15.237 | `src/e1_gmin_m4_prop15237.py` | **(ii-a) ND True** (L²=L pair-span classification) |
| wiring | `src/e1_gmin_m4_prop15193.py`, `15171.py` | `residual_ii_full` = affine ∧ ii-a ∧ ii-b (not exhaustiveness) |

Tests: `tests/test_prop15236.py`, `tests/test_prop15237.py`, and 15169/15170/15171/15193 expect deep_s2/ii-a/full True; e1 still False.

## Residual-(i) hinge (unchanged)

Preferred: \(|\mu|\le 1/(2p)\) on \(|\kappa|=1\) (⇔ Gsum\(\ge-1/p\) ⇒ 15.176 Farkas).  
Sufficient: \(|\mu|\le 2/n\); envelope \(|\mu|\le\max(|\mu_{\mathrm{part}}|,|f_4|)\); dual-eq empty; \(G_+\succ0\) on \(\mathcal W_{++}^0\) (ker=sc) then wire explicit dual; \(K_4\le n(15n-22)\).

**Census (not general — do not flip predicates):**
- \(p=5\): \(|\mathrm{Max}_\pm|=260\), max\(|\mu|=3/65=f_4<1/10\), \(\nu=0\) on all \(|\kappa|=1\).
- \(p=7\): \(|\mathrm{Max}_\pm|=11452\), max\(|\mu|=109/2863<2/n=1/25\), \(\nu=0\) on all \(|\kappa|=1\).
- Interval envelope \(\mu\in\mathrm{conv}\{0,\mu_{\mathrm{part}},f_4\}\) holds on **all** 11700 (\(p=5\)) and 176400 (\(p=7\)) \(|\kappa|=1\) sets; `n_gt_both=0`.
- Dual-eq box-infeasible at \(p=5\) (L∞ \(8/7\)) and \(p=7\) (L∞ \(157/113\)); linearly feasible over \(\mathbb R\).
- L²: \(\|\mu\|_2^2\approx110.6<156\) (\(p=5\)), \(\approx341.6<600\) (\(p=7\)); \(\|\delta\|_2^2\approx11.8,5.2\) vs room \(\approx36.4,227\).

## Proved / certified this session (not shipped as 15.xxx)

These are Max+-free or census facts. **Do not wrap as a new prop** unless one of them flips residual-(i).

1. **Switching.** \(z_\infty=-y_\infty\), \(z_x=y_{nx}\) with \(n\) a nonsquare maps Max+\(\leftrightarrow\)Max− (algebra + census). So \(m_4^-(S)=\sigma(S)\,m_4^+(nS)\). Combined with measured \(\nu=0\), \(m_4^+\) is invariant under nonsquare scaling. Does **not** bound \(|\mu|\).
2. **Far-sum.** \(\sum_{\mathrm{far}}C_{ki}C_{kj}C_{ij}=C_{k0}\) (uses \(\chi(-1)=1\) on \(\mathbb F_{p^2}\), \(T=-1\)). Hence \(\mathrm{diag}(S^+)_k=C_{k0}(n-2)/(p^2(n-3))\).
3. **Comm-repair dual.** Repair \(\Delta=\mathrm{Comm}(\mathrm{diag}((\gamma/\beta)(e_\infty-v)))\) zeros the diagonal. After \(t(J-I)\) the dual has regular degree, Comm, \(W\ge0\), cost\(<2-\alpha\) at \(p=5..23\) (cheaper than \(D(C)\)). **Still needs ker=sc** to close residual (i).
4. **\(\mu\) is a PGL cross-ratio class function** (certified \(p=7\)): CR 9 vs 10 split the extreme-\(\varphi\) class (\(109/2863\) vs \(61/2863\)); CR 7/8 and 16/17 share \(\mu\) within \((\kappa,\varphi)\). No Max+-free \(F(\lambda)\) yet.
5. **Aut \(G\)** for Aut-SOS is squares+Frob only (nonsquare muls flip \(C\)). \(|G|=p^2(p^2-1)\), \(\infty\) fixed.

## Dead / do not re-thrash (this session + prior)

| Attempt | Why dead |
|---------|----------|
| Fréchet alone | only \(\lvert\mu\rvert\le1-2/p\) |
| Joint Fréchet (two 16-atom laws, opposite pairwise, same \(\mu\)) | same weak bound \(1-2/p\) |
| \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise | **false** at \(p=7\) (\(n_{\mathrm{gt}\,f_4}=117600\)) |
| \(\lvert\delta\rvert\le\mathrm{room}_\delta\) pointwise | **false** at \(p=5\) |
| Gsum PSD near \(\mu_{\mathrm{part}}\) | stays on the wall |
| 15.237 pair-span as a Type I kill | Type I *is* a pair-slice; \(\ker F_+\) large |
| Halfspace \(G\)-orbit subframe as full \(\mu\) | differs \(\sim0.6\); \(N=60/168\) vs \(260/11452\) |
| \(d_H\ge p+1\) / \(N_3=(p^2-5)/4\) | true, too weak for \(K_4\) or L∞ |
| Unsigned Per / \(\sum\lvert\mathrm{per}\rvert\) | 15.231 already |
| Aut-line dim\(\le1\) | false |
| W0 (scheme dual) \(\in\mathrm{im}(\mathrm{Gsum})\) | lstsq resid \(0.54\) at \(p=5\); \([S_0,C]\ne0\) |
| Type I + switching \(G=\pi(G)\) | odd rank \(132<600\) at \(p=7\) |
| Low-degree Farkas \(y\) in \(f_e\) | 66 values at \(p=5\); corr \(f_e\sim-0.20\) |
| 3-wise interpolant | \(\mu\) depends on \(\varphi\); \(\varphi\) not a function of the 6 \(K_4\) edges |
| \(\kappa_B\perp E_{\pm4p}\) | rel resid \(0.997\) |
| Matching-PSD / one-edge \(\sum\mu\) | no pointwise |
| max-\(\lvert\mu\rvert\) at extreme \(\lvert\varphi\rvert\) | **false** at \(p=5\) (max at \(\lvert\varphi\rvert=2\)) |
| \(\chi(\)never-zero quadratic\()\) in Max+ | \(0/300\) at \(p=5\) |
| Other Max+ Fourier type | support size \(13=(q+1)/2\), not halfspaces |

## Still viable (next agent: pick one unused)

1. **Prove the interval envelope generally** (\(\mu\in\mathrm{conv}\{0,\mu_{\mathrm{part}},f_4\}\) on \(\lvert\kappa\rvert=1\)). Holds on every certified 4-set at \(p=5,7\).
2. **Prove \(G_+\succ0\) / ker=sc**, then wire the Comm-repair dual (cost already \(<2-\alpha\) at \(p=5..23\)). Rayleigh target \(\lambda_*=8(n-6)/n\).
3. **Prove \(\lvert\mu\rvert\le2/n\)** (weaker than \(1/(2p)\); census-sharp at \(p=5\)).
4. **Prove \(K_4\le n(15n-22)\)** (⇔ residual (i) via 15.197/217).
5. **Find a Max+-free \(F(\lambda)\)** for the CR-class function \(\mu\) (elliptic / Jacobi sum; CR 9 vs 10 is the split to explain).
6. **Dual-eq Farkas for \(x\ge0\)** that is **not** low-degree in \(f_e\).

Un-run at handoff: `{SCRATCH}/resi_gram4.py` (4×4 Gram of \(\{1,f_e,f_f,\chi_S\}\) on Max±). If the joint PSD bound is still only \(1-2/p\), drop it onto the dead list.

## Scratch / caches (ephemeral)

Scratch: `/tmp/grok-goal-51d7ac45c1de/implementer/` (`resi_*.py`, `failed_mechanisms.txt`, `predicate_dump.txt`, `enum_maxplus_p7.py`).  
Prior scratch: `/tmp/grok-goal-cc538b97808e/implementer/`.  
Caches: `/tmp/maxplus_p5.npy` (260×26), `/tmp/maxminus_p5.npy`, `/tmp/maxplus_p7.npy` (11452×50), `/tmp/maxminus_p7.npy`, `/tmp/e1_p7/maxplus.npy`.  
ProcessPool from `python -` stdin **fails** (forkserver cannot import `<stdin>`) — write real scripts.

Compute: 88 cores, idle V100; `recommended_workers` 85. Use ProcessPool/GPU for Max+ contractions.

## Next concrete steps

1. Attempt a **new** Max+-free general-\(p\) residual-(i) bound **not** on the dead list above.
2. Only if that bound is actually proved: wire `gsum_disj_lb_proved_general` / `residual_i_dual_eq_empty_proved_general` / `type_I` / `e1_closed_general` via **real imports**; then set Main Theorem \(L=\tfrac12\).
3. Plan verification (only after a True hinge): predicate dump, hinge `main` twice, hinge pytest → `{SCRATCH}`.

**Do not:** ship a False-predicate 15.xxx; flip e1 on census; soft-close; reopen residual (ii) exhaustiveness; re-thrash the dead list.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee`
