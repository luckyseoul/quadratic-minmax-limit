# Session handoff — residual (i) 2026-08-13 (late)

## Goal and status

**Goal:** Settle MO 413935 (\(L=\lim\alpha_n\)). Binding leftover: residual **(i)** only.

**Status: OPEN — goal NOT met.** Residual (ii) ND is closed. Residual (i) / E(1) / \(L\) still OPEN. Soft-close forbidden. Do not flip `e1_closed_general` on census \(p\le7\). Do not ship a 15.xxx whose predicates stay False.

HEAD: `0443d55` (docs handoff). Last code close: `4d89353` (15.232–237). No new prop this continue.

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

## Residual-(i) hinge (unchanged)

Preferred: \(|\mu|\le 1/(2p)\) on \(|\kappa|=1\) (⇔ Gsum\(\ge-1/p\) ⇒ 15.176 Farkas).  
Sufficient: \(|\mu|\le 2/n\); envelope \(\mu\in\mathrm{conv}\{0,\mu_{\mathrm{part}},f_4\}\); dual-eq empty; \(G_+\succ0\) on \(\mathcal W_{++}^0\) (ker=sc) then Comm-repair dual; \(K_4\le n(15n-22)\).

## Exact 4-point sums (census; GPU + int16)

\(\nu\equiv 0\) and \(\mu\) constant on each CR class at \(p=5,7\). Integer sums \(\sum_{y\in\mathrm{Max}+}\prod_{i\in S}y_i\):

| \(p\) | \(\lvert\mathrm{Max}+\rvert\) | non-split sum | \(\mu\) |
|------|-------------------------------|---------------|--------|
| 3 | 12 | \(4\kappa\) | \(=\mu_{\mathrm{part}}\) (\(\varphi=-2\kappa\)) |
| 5 | 260 | \(8\kappa-2\varphi\) | \(=f_4=(4\kappa-\varphi)/(pn)\) |
| 7 | 11452 | \(240\kappa-10\varphi+96\,\kappa\,\varepsilon\,1_{\lvert\varphi\rvert=10}\) | CR 9 vs 10 split at extreme \(\varphi\) |

At \(p=7\): \(\varepsilon=+1\) on CR 9 (\(\mu=\pm109/2863\)), \(\varepsilon=-1\) on CR 10 (\(\pm61/2863\)). Mean of the split is the linear part \(240\kappa-10\varphi\). Interval envelope holds on **all** \(\lvert\kappa\rvert=1\) sets at \(p=3,5,7\).

Do **not** flip predicates on these formulas. Need a Max+-free \(F(\lambda)\) or a general envelope proof.

## Dead this continue (do not re-run)

| Attempt | Why dead |
|---------|----------|
| 4×4 Gram \(\{1,f_e,f_f,\chi_S\}\) | only \(\lvert m_4\rvert\le1-2/p\) |
| Affine-quadratic level sets on AG(2,5) | only the 60 linear halfspaces; 200 Max+ are not \(Q^{-1}(S)\) |
| CR-class master+diamond LP | 10 / 18 classes; HiGHS \(\max\lvert\mu\rvert=1\) (\(E_{\pm4p}\) survives) |
| Single IP-valency \(K_4\) | regular at \(p=5\) (valencies \(1,13,20,36,60\) + negs); **not** regular at \(p=7\) |
| Global mix \(\mu=(1-t)\mu_{\mathrm{part}}+t f_4\) | \(t\) depends on \((\kappa,\varphi)\) (\(p=7\): \(\approx0.20\) on \(\varphi=\pm2\), \(\approx0.04\) on \(\varphi=\pm6\)) |

## Earlier dead (still dead)

Fréchet / joint Fréchet; \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise; \(\lvert\delta\rvert\le\mathrm{room}_\delta\) (false at \(p=5\) on \((\kappa,\varphi)=(-1,-6)\): \(\lvert\delta\rvert\approx0.039>\mathrm{room}=0.028\) while \(\lvert\mu\rvert=\lvert f_4\rvert\) still safe); W0\(\in\mathrm{im}(\mathrm{Gsum})\); Type I+switching \(G=\pi(G)\); low-degree Farkas \(y\); 3-wise interpolant; \(\kappa_B\perp E_{\pm4p}\); matching-PSD / one-edge \(\sum\mu\); max-\(\lvert\mu\rvert\) at extreme \(\lvert\varphi\rvert\); halfspace orbit as full \(\mu\); \(\chi(\)never-zero quadratic\()\); unsigned Per; Aut-line dim\(\le1\); 15.237 pair-span as Type I kill.

## Still viable

1. **General interval envelope** \(\mu\in\mathrm{conv}\{0,\mu_{\mathrm{part}},f_4\}\) on \(\lvert\kappa\rvert=1\). Holds exactly at \(p=3,5,7\). Cleanest close.
2. **Max+-free \(F(\lambda)\)** matching the three exact sums above (need \(\lvert\mathrm{Max}+\rvert\) and the extreme-\(\varphi\) split character \(\varepsilon\)).
3. **\(G_+\succ0\) / ker=sc**, then wire Comm-repair dual (cost already \(<2-\alpha\) at \(p=5..23\)).
4. **\(\lvert\mu\rvert\le2/n\)** or **\(K_4\le n(15n-22)\)**.
5. Dual-eq Farkas **not** low-degree in \(f_e\).

## Geometry / IP (census)

- \(p=5\) Max+: 60 linear halfspaces + 200 other; 4 AG(2,5) line-intersection types. Other type is not affine-quadratic.
- \(p=5\) IP-regular; \(p=7\) not (row histograms differ).
- \(T\) definition (15.68, Max+-free): \((Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r})\).

## Scratch / caches

Scratch: `/tmp/grok-goal-51d7ac45c1de/implementer/` (`resi_cr_formula.py`, `resi_cr_lp.py`, `resi_mu_closedform.py`, `resi_fit_quad.py`, `resi_ip_int.py`, `resi_gram4.py`). May be deleted when the goal ends.  
Caches: `/tmp/maxplus_p5.npy` (260×26), `/tmp/maxminus_p5.npy`, `/tmp/maxplus_p7.npy` (11452×50), `/tmp/maxminus_p7.npy`.  
ProcessPool from `python -` stdin **fails** — write real scripts. 88 cores + idle V100; `W=86`.

## Next concrete steps

1. Prove a **new** Max+-free general-\(p\) residual-(i) bound **not** on the dead list (preferred: envelope or \(F(\lambda)\) from the exact sums).
2. Only if proved: wire `gsum_disj_lb` / `type_I` / `e1_closed_general` via **real imports**; then Main Theorem \(L=\tfrac12\).
3. Plan verification only after a True hinge.

**Do not:** ship a False-predicate 15.xxx; flip e1 on census; soft-close; reopen residual (ii) exhaustiveness; re-thrash the dead list.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee`
