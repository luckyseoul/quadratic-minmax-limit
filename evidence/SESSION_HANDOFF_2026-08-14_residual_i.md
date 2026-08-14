# Session handoff — residual (i) 2026-08-14

**Campaign NOT met.** \(L\) OPEN. Residual (i) / Type I / Gsum / E(1) still False.
Do not flip predicates. Soft-close forbidden. Do not re-search the MO/X thread.
After `/goal clear` start a **new** `/goal` (cannot resume a cleared goal).
Token budget is not a constraint this week — no `--budget` flag.

Scratch (this continue, do **not** ship as 15.xxx):
`/tmp/grok-goal-838009eda84d/implementer/`
Source-of-truth writeup there: `PHI_KAPPA_B.md`.

## Live gates (checked 2026-08-14)

| Flag | Value |
|------|-------|
| `residual_ii_full_closed` | **True** |
| `bitight_from_majorization(5)["bitight_empty"]` | **True** |
| `gsum_disj_lb_proved_general` | False |
| `residual_i_dual_eq_empty_proved_general` | False |
| `residual_i_closed_via_249` | False |
| `type_I_k_3p_minus_2_closed_general` | False |
| `e1_closed_general` | False |
| `solution.md` Main Theorem | **\(L\) OPEN** |

E(1) is `type_I ∧ residual_ii ∧ bi-tight(p=5)`. Residual (ii) and bi-tight are already True.
The only official blocker is residual (i).

HEAD before this commit: `208c6ea` (15.238–268, \(\nu=0\)). This commit adds 15.269 (Fourier support / Wick / \(\kappa_3\) **criterion**; \(\kappa_3\) bound still OPEN) plus this handoff. 15.269 does **not** close residual (i).

## Binding leftover (one sentence)

Prove \(N(\varphi)\neq0\) (equivalently \(G_+N_++G_-N_-\neq0\), equivalently \(M(\Psi)\neq0\)) for every even character \(\varphi\) of \(\mathbb F_{p^2}^*\) that is nontrivial on \(\mathbb F_p^*\), for **all primes \(p\ge7\)**. Then \(G_+\succ0\) on \(\mathcal W_{++}^0\), \(\ker=\mathrm{sc}\), 15.249 closes residual (i). Wire predicates by **real import**. Set `solution.md` \(L=1/2\) **only after** live gates are True.

\(p=5\) is already \(G_+\succ0\) by the halfspace Veronese (orbit = all 260 Max+; rank 65/65). Inversion-T has a kernel at \(p=5\) (2 of 3 K-characters, \(N=0\)); do not use inv-T as the \(p=5\) hinge.

## Why this is the hinge

15.207: \(\ker(\mathrm{Gsum})=\mathrm{sc}\) \(\Leftrightarrow G_+\succ0\) on \(\mathcal W_{++}^0\).
15.249: Comm-repair dual \(D_{\mathrm{alg}}\) already gives free-\(e\) over \(\mathrm{sc}\) with \(\mathrm{cost}_D<2-\alpha\) (Weil \(|Q|\le2p\)). So residual (i) reduces to \(\ker=\mathrm{sc}\).

Signed PSL on \(S^2(V_+)\) (\(V_+\) even Weil). No cuspidal constituents of \(\mathcal W_{++}^0\) (inner product with every discrete series \(\theta\) of degree \(q-1\) cancels). Every irrep is principal series or Weil, hence has \(U\)-invariants. G-span of the translation-invariant circulants \(F\) is all of \(\mathcal W_{++}^0\) (certified rank = dim at \(p=5,7\)). Schur: \(G_{+,\mathrm{hs}}\succ0\) on \(\mathcal W_{++}^0\) iff PD on \(F\).

\(F=F_{\mathrm{aff}}\oplus K\), \(\dim F_{\mathrm{aff}}=(p+1)/2-1\), \(\dim K=(p+1)(p-3)/4\). Affine (halfspace) disks span \(F_{\mathrm{aff}}\) (\(\hat\rho(k)=|\hat s(k)|^2\neq0\) all \(k\in\mathbb F_p\)). Affine + \(T\)-orbit of the **inversion disk** spans \(F\) at every prime \(7\le p\le23\) (certified). That span on \(K\) is equivalent to the even-K Mellin of \(|\hat z_{\mathrm{inv}}|^2\) never vanishing.

## Closed form (proved Max+-free, verified \(10^{-14}\) at \(p=5,7,11,13\))

Field encoding: `ia=0`, \(\omega^2=\mathrm{ib}\) nonsquare, \(\mathrm{Tr}(\omega)=0\), \(\chi_q(x)=\chi_p(N(x))\). Halfspace \(s=+1\) on \(\{0,\dots,(p-1)/2\}\).
Inversion disk: \(z(0)=1\), \(z(x)=\chi(x)\,s(L(x^{-1}))\) for \(x\neq0\). \(z\equiv+1\) on \(\mathbb F_p\). \(\sum z=p\), \(R(0)=p^2\), \(\sum_d R(d)=p^2\).

\(\tau(\chi)=p\,\chi(\omega)\). \(\Omega=\{\xi:\chi(\xi)=\chi(\omega)\}=\omega\cdot(\mathbb F_q^*)^2\).
\(\hat z(\xi)=0\) on \(\mathbb F_q^*\setminus\Omega\).

Write \(t=x+y\omega\neq0\), \(\xi=\omega t^2\in\Omega\):
\[
\hat z(\omega t^2)
= p\,1_{x=0}+p\,1_{y=0}
+ G_p\bigl[\chi_p(-2)\,J(2x^2)+\chi_p(-2\,\mathrm{ib})\,J(2\,\mathrm{ib}\,y^2)\bigr],
\]
where \(G_p=\sum_{r\in\mathbb F_p}e^{2\pi i r^2/p}\) (\(|G_p|^2=p\)) and
\[
J(\mu)=\sum_{k\neq0}s(k)\chi_p(k)\,e^{2\pi i \mu k^{-1}/p}.
\]
Script: `zhat_separable.py` (the \(-1\) constant is **wrong**; omit it).

Salié: \(\sum_{u\neq0}\chi(u)\psi(Au+B/u)=0\) if \(AB\) nonsquare (\(A,B\neq0\));
\(=\tau\chi(A)(\psi(2s)+\psi(-2s))\) if \(s^2=AB\neq0\).

Functional equation already proved: \(M=\tau(\xi)N(\xi^{-1})\) (magnitudes), so \(M=0\Leftrightarrow N=0\).
\(N(\varphi)=\sum_{d\neq0}R(d)\varphi(d)\). Even K-characters: \(j\) even and \(j\not\equiv0\pmod{p-1}\).

## Certified, not a proof

2-D additive FFT + multiplicative FFT (`N_fast_sweep.py`): \(N\neq0\) on every even K-bin for every prime \(7\le p\le79\). \(N=0\) only at \(p=5\) (4 of 6 FFT bins = 2 of 3 characters).

| \(p\) | \(\lvert N\rvert_{\min}\) | \(\min/p\) |
|------|---------------------------|-----------|
| 7 | 56 | 8.00 |
| 11 | 83.9 | 7.63 |
| 13 | 28.8 | 2.22 |
| 17 | 12.8 | 0.76 |
| 19 | 222 | 11.7 |
| 37 | 33.5 | 0.91 |
| 79 | 966 | 12.2 |

Smallest relative values at \(p=17,37\). No drift to 0. **Census is not a general proof. Do not flip.**

Even Mellin of \(\lvert J(2x^2)\rvert^2\) on \(\mathbb F_p^*\) vanishes for some even \(\lambda\) when \(p\equiv1\pmod4\) (e.g. \(p=5,13\)). Then \(N\) is carried by the \(B(y)\) piece and/or the cross term \(A(x)\overline{B(y)}\). At \(p=5\) those cancel; at \(p=13\) they do not.

## Wiring when (and only when) \(N\neq0\) is proved

1. New 15.270 (or next free number): \(G_+\succ0\) on \(\mathcal W_{++}^0\) via affine + inv-T, \(p=5\) by Veronese, \(p\ge7\) by \(N\neq0\). Closed-flag True **only if** the bound is actually proved.
2. Import into `residual_i_closed_via_249` (today hardcoded `return False`).
3. 15.216 `residual_i_dual_eq_empty_proved_general` or 15.170 `gsum` — real import.
4. `type_I_k_3p_minus_2_closed_general` becomes True via dual-eq or gsum.
5. `e1_closed_general` becomes True (ii and bi-tight already True).
6. `solution.md` Main Theorem \(L=1/2\), not OPEN.
7. Two 15.170 gate launches + honesty tests under `{SCRATCH}/`.

Do **not** ship another identity whose only `proved: True` items leave residual (i) open.

## Dead this continue (do not re-run)

- Envelope \(\lvert\mu\rvert\le2/n\) (false at \(p=11\)).
- Comm \(S_p^{\mathrm{off}}\) after dropping non-constant diag (does not commute with \(C\)).
- 3-param Aut\(_\infty\) templates (miss \(\lambda_*\) space).
- Half-Radon / inversion+squares not injective at all \(p\).
- Unsigned Aut inversion (does not preserve \(C\)); need signed PSL.
- Halfspace unsigned orbit rank 29/65 at \(p=5\).
- Affine-only disks / inversion-T DFT vanishing at \(p=5,13\) is **not** a kernel of affine+inversion **together**.
- Stickelberger digit-sum \(s(k)=s(k+(q-1)/2)\) only separates \(\sim1/3\) of K-characters.
- Two-term \(M=\tfrac12 G_+N_++\tfrac12 G_-N_-\) implemented as if \(\Omega=\)squares: falsely reported combo \(=0\) at \(p=13,17\) (\(p\equiv1\pmod4\), \(\Omega=\)nonsquares). Use the unified \(M=\tau N\) form.
- Classwise \(A\le0\), \(\lvert\mu_4\rvert\le1/(2p)\) as a class bound (false at \(p=17\)), \(\lvert\delta\rvert\le\mathrm{room}\) (false at \(p=5\)), Aut-line \(\dim\le1\), 4-design, CS on \(\rho\).

## Still viable

1. **Prove \(N\neq0\)** from the separable \(\hat z\) (Jacobi/Gauss expansion of \(\lvert A(x)+B(y)\rvert^2\) against \(\eta(x+y\omega)\); Stickelberger on the resulting Gauss sums; or show \(\hat s^*T_\varphi\hat s\neq0\)).
2. Another listed residual-(i) hinge if this stalls: Aut-SOS / \(P_\pm\) / \(G_+=B^*B\); \(K_4\le\mathrm{Wick}_{hi}\); envelope/reflection (criteria proved, hyp open).
3. 15.269 \(\lvert\kappa_3\rvert\le(p^2-2p-4)/(2p^3)\) is an **alternate** residual-(i) criterion, not the active path.

## Key scratch scripts

`field_fp2.py`, `hatR_inv_formula.py`, `zhat_separable.py`, `N_fast_sweep.py`, `M_is_tau_N.py`, `circ_fq_psl_rank.py`, `gspan_circulants.py`, `inv_vs_affine_ker.py`, `wpp0_character.py`.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee` (optional; use when stuck on \(N\neq0\), not every turn)

## Compute

88 cores, ~60 GiB RAM, idle V100. Prefer ProcessPool \(W\approx86\) or one CUDA context for dense FFTs. Never serial multi-minute loops. Machine was idle at last budget snapshot.

## Do not

- Flip `gsum` / `type_I` / `e1` / `residual_i_closed_via_249` until \(N\neq0\) is proved for all \(p\ge7\).
- Set `solution.md` \(L=1/2\) while any of those is False.
- Re-search MO 413935 / the X thread after compression.
- `/goal resume` after `/goal clear`.
- Ship a 15.xxx that only records the closed form and leaves the bound open.
- Soft-close via sandwich + denseness + \(\rho=1\).
