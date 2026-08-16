# Session handoff (2026-08-16, leftover campaign)

**Repo:** `/home/nick/quadratic-minmax-limit` · `luckyseoul/quadratic-minmax-limit` `main`  
**Binding goal:** `GOAL.md`. Do **not** `/goal clear`. Do **not** unflip `e1` / `L`. Do **not** flip Aut-Schur / Gsum / pairing. Do **not** run `pbss-goal-verifier` or `perry-beurling-rh-closeout`.

## Goal

Prove \(L=\lim\alpha_n=1/2\) by closing four leftovers so E(1) is Max+-free for all primes \(p\ge5\). Live `e1_closed_general` is still True by the **old AND** (bi-tight only at \(p=5\); 15.236/237 treated as full residual (ii)). That wiring is not acceptance.

## Live flags (do not change unless a sign actually closes)

| flag | value | note |
|------|-------|------|
| `phi_F_ge_6_proved_general` | False | floor open |
| `phi_F_structure_proved` | True | 15.279 A–**Y** |
| `e1_closed_general` | True | dishonest old AND |
| Aut-Schur / Gsum / pairing | False | stay False |
| `residual_ii_k_ge_4p_ND_closed` | False | multi-level leftover |
| `type_I_aut_e_3AB_positive_general` | False | 3-point unsigned |
| `type_I_multilevel_bad_case_ND_closed` | False | |

## Shipped this wave (structure, not the floor)

- **15.273–15.277:** residual-(ii) / Type I / Aut_e identities. Inequality flags stay False.
- **15.276:** Lemma D written (`evidence/share/A3_PROOF.md` + live `Cy=py`).
- **15.278:** Aut·F=Z certified at \(p=5,7\) (not a general floor).
- **15.279 A–X:** character-pair / regular-set dictionary. **X:** \(L(1/r)=L(r)\); S3 \(r\mapsto-1-r\) is **not** a symmetry; leftover \(n_{\mathrm{aut}}-2\).
- **15.279 Y:** 1D \(\hat N\) closed form. \(N=\lvert D\rvert\) on \(\ker L\), \(p\,n_M(L(\delta))\) off. \((p-1)\mid k\Rightarrow\lvert\hat N\rvert^2=[p(p^2-1)/4]^2\); else \(p^3\lvert\sum_{x,y\in M}\eta(x-y)\rvert^2\). Paley plus-set vanishes on even ntriv \(k\) when \(p\equiv3\pmod4\). Fail-when-wrong: inverted \((p-1)\mid k\) slot. **Does not exhaust Max+ and does not sign ensemble \(\hat L\).**

Tests: `tests/test_prop15279.py::test_1d_Nhat_closed_identity` and `test_structure_proved_floor_and_schur_untouched` exit 0.

## Floor (I1) — still OPEN

Equivalent: \(\hat L(\psi)\ge3q/16\) \(\Leftrightarrow\) \(\mathrm{Term}_4\ge0\) \(\Leftrightarrow\) \(\lambda_{\min}(\Phi\rvert_F)\ge6\) for even \(\psi\notin\{1,\chi\}\).

**Must be an Aut_∞ mixture**, not per-orbit positivity (false: many orbits have \(\hat N=0\)).

- Aut_∞ = \(\langle\)trans, \(\times\square\), Frob\(\rangle\), \(\lvert G\rvert=p^2(p^2-1)\). Same orbits as 15.140.
- \(p=5\): orbits \(30\) (all 1D) + \(100\) (hs\(\odot\)nc). Worst \(k=2\): \((30/130)\times500=1500/13>225/2\). Nonlinear contributes 0.
- \(p=7\): \(56+84+294+588+4\times1176\). Worst \(k=8\): 1D supplies \(\approx60\) of \(441\). Nonlinear required.
- \(N_+(p)\) unknown (\(6,130,5726\) at \(p=3,5,7\)). 15.144 type-list dead for \(p\ge11\).
- Nonlinear regular sets are \(D_A\Delta T\), \(T\) a Hoffman coclique of \(C_{y_A}\). Norm-circle \(T\) empty at \(p=13\).

**Dead (do not reopen):** Paley+\(\omega\) BM (\(p=11\)); Aut box+PSD+simplex (min \(\hat L=0\)); S3 merge (X); pointwise \(\lvert M\rvert^2\) SOS; typewise \(\hat g\mu\ge0\); Paley-Ising DFT as Gauss product; Weil on \(R_{\mathrm{rest}}\) (not unit-conductor; Deligne recycles the floor); 1D-only mixture (\(p=7\) \(k=8\): \(60.38<441\)); Wick-as-equality (\(\lambda=8\)).

Wick of Paley \(C\) is the only 4-point **above** the floor: \(\mathbb E\lvert\hat N\rvert^2_{\mathrm{Wick}}=q(q-1)/4\), room \(q(q-1)/16\). Remainder \(R\) is the leftover Aut-orbit Fourier of \(Q\).

## Type I — still OPEN

Hoffman \(\lvert D\cap C_\boxtimes\rvert=(p-1)/2\) pathwise. Leftover is the \(I_\square\) mixture after \(k_\boxtimes\) is locked.

- \(\beta\in\mathbb F_p\): \(I_\square\) is within-line cyclic auto; unbalance helps.
- \(\beta\notin\mathbb F_p\): \(I_\square\) is between-line transport; balance helps.
- Any \(k\)-only type-weight is dead (opposite \(\mathrm{Cov}(V,\mathrm{slack})\)).
- Both classes have positive slack at \(p=5,7\) (tightest: \(p=7\) \(\beta=4\in\mathbb F_p\), \(92/409\)).
- \(\mathrm{Wick}_\square+\mathrm{Wick}_\boxtimes=-4/p^3\), so \(A\le0\Leftrightarrow\kappa3_\square+\kappa3_\boxtimes\le4/p^3\) (unsigned).

## Residual (ii) \(k\ge4p\) — still OPEN

1D-both 0-1 **exists**. The 28 rank-3 leftovers at \(p=5\) are Type+\(\odot\) finite norm-circle flips (ham \(p+1=6\)). Type+\(\cup\) circle-partners = all Max+ at \(p=5\). Hoffman-line partners are a different family (\(S\ge2\) under 1D-both).

Cut: \(S(y\odot\varepsilon_T)=S(y)-2\Delta\). Both leftover+ \(\Leftrightarrow\mathrm{Sint}\ge2\). On tight Type+ (\(S=2\)), leftover+ on \(z\Leftrightarrow\Delta\le0\). 1D-both has \(\Delta\in\{1,2,3\}\) on the 28. Fractional Type+-both + leftover+ on the 14 is feasible. Two-sided on cylinders is the extra Boolean.

## Live children (if this session is still up)

- I1 \(\hat N(D_A\Delta T)\): `01a007ee-215d-7762-9eeb-b4b05e6dc49f`
- Type I \(\mathbb F_p\) vs not: `01a007e8-ff09-7ab2-96f1-04b8a47dd0d3`
- Res-ii Sint\(\le1\): `01a007f0-6e26-76d1-9982-b877b2275f44`

Join those before launching a same-hypothesis duplicate. Scratch graph: `/tmp/grok-goal-e664043f8ee2/implementer/leftover_graph.md` (scratch only; not in git).

## Next concrete steps

1. Join the three children. Wire a unit only if it is Max+-free and failing-when-wrong. Do not flip `phi_F_ge_6` / Type I / residual-(ii) / `e1` / `L` unless the inequality or emptiness actually closes for all \(p\ge5\).
2. I1: evaluate \(\hat N(D_A\Delta T)\) from Hoffman parameters of \(T\), plus a formula for \(N_+(p)\). Mixture \(\ge3q(q-1)/16\).
3. Type I: sign slack separately on \(\beta\in\mathbb F_p\) and \(\beta\notin\mathbb F_p\) (not a \(k\)-only weight).
4. Res-ii: force \(\mathrm{Sint}\le1\) on Type+\(\odot\) staying-Max+ circle partners for every leftover+ two-sided Type+ \(G\) of size \(4p\).
5. Only then import the four units, clear package caveats, set `solution.md` \(L=1/2\). Verification: two 15.170 dumps, honesty pytest `-n W`, cold Q1/Q2.

## Pointers

- Goal: `GOAL.md`
- Package caveats: `evidence/share/denseness_path_package.md`
- Lemma D writeup: `evidence/share/A3_PROOF.md`
- Floor module: `src/e1_gmin_m4_prop15279.py`, tests `tests/test_prop15279.py`, evidence `evidence/e1_gmin_m4_prop15279.json`
- Aut·F=Z: `src/e1_gmin_m4_prop15278.py`
- Type I multi-level: `src/e1_gmin_m4_prop15275.py`
- Residual (ii) leftover: `src/e1_gmin_m4_prop15274.py`
- Scratch attacks: `/tmp/grok-goal-e664043f8ee2/implementer/attack_term4_*.txt`, `attack_typeI_*.txt`, `attack_resii_*.txt`
- Caches: `/tmp/maxplus_p5.npy`, `/tmp/e1_p7/maxplus.npy`, `/tmp/maxminus_p{5,7}.npy`
- Prior: `evidence/SESSION_HANDOFF_2026-08-15_aitest.md`

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `use-available-compute` · `scientific-critique` · `grill-me` · `self-refine-loop` · `handoff` · `session-handoff-packager` · `verification-before-completion` · `research` · `arxiv` · `litreview` · `openai-referee`

Do **not** load `pbss-goal-verifier` or `perry-beurling-rh-closeout` on this repo.
