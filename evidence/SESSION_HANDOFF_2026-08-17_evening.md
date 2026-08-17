# Session handoff — 2026-08-17 evening

**Three of four leftovers remain.** Lemma D is already True. Do not unflip it.

| Leftover | Flag | Status |
|----------|------|--------|
| \(\lambda_{\min}(\Phi)\ge6\) / \(Q_\tau\) | `phi_F_ge_6=False` | OPEN. Live pairings \(\le2\) at p=5,7 (24/13, 200/409). 15.494 names \(A_r=0\) on \(\chi(r+1)=-1\), not \(Q_\tau\). |
| Residual (ii) even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | OPEN. 4-level first-moment exists. leftover+\(s_+=2\) ND at p=3. |
| Type I multi-level | `type_I_multilevel_bad_case_ND_closed=False` | OPEN. Far \(\lvert\kappa\rvert=1\) pairing \(\equiv0\) on Max±; \(3A+B\) still \(G>T\). |
| Lemma D | True | Do not cascade. |

**Latest units:** 15.488–15.494 on `main` (`81ed01e`). Next attack: name \(L_\tau\) (then \(Q\) is the 15.298 half-Gauss image) or \(A_4\) on 15.290 types. Fail-when-wrong. Do not import \(\phi_F\).

**Do not:** soft-close \(L\); flip `e1` / Aut-Schur / Gsum / pairing; p=11 6-net DFS; 2-point interpolants of live \(Q\); jellyfin CPU unless soulkiller is saturated.

**Machines:** soulkiller 88 thr + V100 owns CPU/CUDA. Jellyfin A380 is GPU-only overflow (`clWaitForEvents` busy-waits — use event poll). Orin sm_87 overflow.
