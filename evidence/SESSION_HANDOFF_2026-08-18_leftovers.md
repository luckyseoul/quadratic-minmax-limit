# Session handoff — leftover closeout (2026-08-18)

Fresh agent: start here, then `HANDOFF.md`, `STATUS.md`, `GOAL.md`, `src/e1_gmin_m4_prop15562.py`.

## Current goal

Prove \(L=\lim\alpha_n=1/2\) by making E(1) Max+-free for all primes \(p\ge5\).
**Not done.** Live `e1_closed_general` is still the old AND (True). Public \(L\) is OPEN.
Do not soft-close. Do not `/goal clear`. Do not flip leftover flags without import + fail-when-wrong.

Repo: `/home/nick/quadratic-minmax-limit` · `main` · `https://github.com/luckyseoul/quadratic-minmax-limit`

## Four leftovers

| Leftover | Flag | Status |
|---|---|---|
| \(\lambda_{\min}(\Phi)\ge6\) \(\Leftrightarrow\langle\delta,\psi\rangle\le2\) | `phi_F_ge_6=False` | **OPEN.** Reduced to naming \(A_{\mathrm{full}}\) for \(m>3\) (or 16pA / \(Q_\tau\)). Naming \(A_{\mathrm{full}}\) names \(Q_3\), **not** automatically the pairing bound. |
| Residual (ii) even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | **OPEN.** Census only. |
| Type I multi-level | `type_I_* = False` | **OPEN.** Aut_e **DEAD** as a name of \(A_{\mathrm{full}}\) (15.559). |
| Lemma D | True | **CLOSED** (15.276 / 15.272). Do not unflip. |

Aut-Schur / Gsum / pairing stay **False**. House **15.x** unless the user calls the campaign failed (16.x) or \(L=1/2\) is honest (17.x).

## Referee (this wrap)

- Claude `deep_review`: **PASS-WITH-NOTE**, `do_not_branch`.
  - Hinge is \(A_{\mathrm{full}}\) **modulo** \(Q_{3,02}=-4N(2p^2+1)/p\) (from \(F=-2(3p^2+2)\)), certified \(p=5..23\), not proved for all \(p\).
  - Naming \(A_{\mathrm{full}}\) is the last unit for **\(Q_3\)**, not shown to be the last unit for \(\langle\delta,\psi\rangle\le2\).
- OpenAI `deep_review` (same slot, after Claude): **PASS**, `do_not_branch`. Same reading.

## Leftover-1 named chain (do not re-derive)

| Unit | Identity | Still open |
|---|---|---|
| 15.550 | \(S(\lambda)=\mathrm{Kl}(1,\lambda^2/4)\) every odd \(q\) | — |
| 15.553 | Term0 of \(K_{\mathrm{all}}\) | Ω-bulk |
| 15.555–7 | \(G_3=\Phi\circ N=G_{3,1d}+\Phi_{\mathrm{free}}\); pair-sum of \(\Phi_{\mathrm{free}}\) named | pointwise \(\Phi_{\mathrm{free}}\) |
| 15.558 | \(J_{\mathrm{all}}=(1/8)\sum\chi_\Omega^{\|\varepsilon\|}S_\varepsilon\); \(J_T\) via \(S_\Omega\); \(Q_{3,02}=-4N(2p^2+1)/p\) | \(Q_{3,T}\) live coeffs |
| 15.561 | \(n_{k=3}=C(m,3)(p-1)q\); \(Q_{3,\mathrm{generic}}=(N-n_{k=3})A_{1d}\) | — |
| 15.562 | \(\sum\omega^{2t}/(\omega^t-1)^3=-(p^2-1)/24\); \(A_{k=3,\mathrm{dbl}}=0\); \(A_{k=3,n3}=-16p^3/((p-1)(p-3))\); \(n_{k=3}A=-p^5(p^2-1)/3\) (= \(Q_{3,n3}\) at \(p=5\)) | \(A_{\mathrm{full}}\) for \(m>3\) |
| 15.559 | Aut_e inversion mixes full-Ω with nsupp=18 at \(p=7\) (24/2520) | Aut_e DEAD for \(A_{\mathrm{full}}\) |
| 15.563 | p=7 full-Ω two energy types \((1,1,1,3)\) \(n=36q\) and \((1,1,2,2)\) \(n=54q\); \(\hat z\) not Fejer | \(A_{\mathrm{full}}\) still not a p-law |

At \(p=5\) (\(m=3\)): \(Q_{3,\mathrm{dbl}}=n_{1d}A_{1d}\), \(Q_{3,n3}=-p^5(p^2-1)/3\). Done.
At \(p=7\):
\[
Q_{3,\mathrm{dbl}}=n_{1d}A_{1d}+n_{\mathrm{full}}A_{\mathrm{full,dbl}},\quad
Q_{3,n3}=-p^5(p^2-1)/3+n_{\mathrm{full}}A_{\mathrm{full,n3}}.
\]
Ensemble \(A_{\mathrm{full,dbl}}=A_{\mathrm{full,n3}}=-4p^3/15\) is a **mean**, not a per-vector constant. **Do not interpolate \((p-5)/15\).**

## Scratch (not shipped as a p-law)

`/tmp/grok-goal-a558c5f11751/implementer/cpu_full_*.py` and `cpu_full_AB.out`.

- Fibre-sum on a line: \(G=\mathrm{ifft}(\hat z|_L)\cdot p\).
- k=3: \(G=p(2N-p)\) full sawtooth.
- A **subset** of full-Ω line restrictions: \(G(t)=p\bigl(1+2\chi(t-\alpha)-2\chi(t-\beta)\bigr)\) (exact at p=5 on all lines; ~61% of (vector,line) at p=7).
- Every p=7 full-Ω vector has **2 or 3** of its 4 Ω-lines of that 2χ form (never 0,1,4).
  - Type A, 3-of-4: \(n=1764=C(m,2)(p-1)q\). \(A_{\mathrm{dbl}}=+4p^3/9\), \(A_{\mathrm{gen}}=-4p^3/3\), \(A_{n3}=-4p^3/9\).
  - Type B, 2-of-4: \(n=2646\). \(A_{\mathrm{dbl}}=-20p^3/27\), \(A_{\mathrm{gen}}=-4p^3/9\), \(A_{n3}=-4p^3/27\).
  - Mix recovers \(-4p^3/15\). These denominators 9,27 are p=7 rationals until a construction exists.
- Affine-quadratic / \(N(x-c)\) do not name \(z\). Uniform \(\alpha\neq\beta\) Gauss average of the 2χ DFT product is 0, not \(-4p^3/15\) (the live \((\alpha,\beta)\) law is not uniform).
- Unmatched (non-2χ) fibre \(G\) lstsq-fits a different \(\chi\)-linear form (includes \(\pm p^2\)); not named.

## Residual (ii) — do not MIP-clone

Official leftover+\(s_+\) is \(s_+\ge2\), not equality.
- p=5 \(k=20\): leftover+splus empty **all** nF (15.528). leftover-only nF=8 exists.
- p=5 \(k=22\): empty nF=0,3–9,11–14 (15.547+15.552). nF=10 TLE.
- p=5 \(k=26,28,30\): nF=0 empty (15.560). leftover-only empty at k=32.
- k=24 TLE. No p≥7. No structural ND.

## Type I

15.536–15.551 named I and Galois line-unions. 15.546 p=7 mix is census. 15.559 kills Aut_e as \(A_{\mathrm{full}}\) name. Next Type I path is **not** Aut_e.

## Dead (do not reopen)

Occupancy / Aut-involution / n_bar0 / 57+57 interpolants / half-net census / 1-term GJ catalogs / Paley-χ of \(G_3\) / feature-lstsq of \(\Phi_{\mathrm{free}}\) / \(4p\chi\) collinear / quadratic-form seeds / Aut_e for \(A_{\mathrm{full}}\) / \((p-5)/15\).

Do **not** commit: `prop15496`, `prop15530`, `prop15493/495` catalogs, leftover+splus equality JSONs.

## Next session (outcome-changing)

1. **Construct** the p=7 type-A / type-B full-Ω vectors (2χ fibre on 3 or 2 lines). If the construction is Max+-free in \(p\), average the Gauss/DFT product to name \(A_{\mathrm{full},T}\) without the cache value \(-4p^3/15\).
2. Independently: prove \(F=-2(3p^2+2)\) for all odd \(p\) (or find the true p-law). Claude: do not silently promote the p≤23 certificate.
3. After \(Q_3\) is fully named, write the map \(Q_3\to\langle\delta,\psi\rangle\le2\) **before** importing `phi_F_ge_6`.
4. Residual-ii: structural ND only; no new nF=10 / k=24 MIP; no JF MIP clones.
5. Type I: not Aut_e.

Caches: `/tmp/maxplus_p5.npy` (260=2·130), `/tmp/maxplus_p7.npy` (11452=2·5726). No p=11 \(H_+\).

## Compute

Soulkiller: 88 threads, V100. Inventory `~/.grok/skills/use-available-compute/scripts/compute-budget.sh`. GPU if dense fit, else ProcessPool `full_workers≈86`.
Jellyfin `192.168.1.191`: CPU only when SK occupied. Orin `100.67.236.54`: W=4 field jobs. Creds as previously given — **never write them in the repo**.

Scratch only: `/tmp/grok-goal-a558c5f11751/implementer`.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `use-available-compute` · `claude-referee` (primary) · `openai-referee` (after Claude) · `handoff` · `session-handoff-packager` · `verification-before-completion` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview`

Do **not** run `perry-beurling-rh-closeout` or `pbss-goal-verifier`.

## Category board (if asked)

Already-proved path 8×100 · leftover-1 ~45 · \(L=1/2\) 0%.
