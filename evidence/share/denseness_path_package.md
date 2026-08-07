# Denseness-path package: lim α_n (honest status)

**Audience:** independent AI / human check (Paata AI-test).  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-07 (revised after cold AI-test BLOCK on residual-i)

## Verdict

**\(L=\lim\alpha_n=\tfrac12\) is NOT proved.**  
Residual **(i)** is **OPEN** (15.216 K₄ thr path has a fatal Rayleigh gap on \(\delta\in E_{4p}\)). Residual (ii) affine is closed; **full residual (ii)** is open. Both block E(1)/L.

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** |
| Denseness (ratio-dense family ⇒ global liminf/limsup) | **Proved** (Prop 6.1) |
| Bi-tight empty via majorization (15.167) | **Candidate** (algebra OK if mult≥d−1 and λ_min≥6) |
| Residual (ii) affine dual two-level (15.179) | **CLOSED** — freeze-to-tight, no Gsum |
| Residual (ii) full / exhaustiveness (15.193) | **OPEN** — multi-level / non-affine freeness-fail |
| Type I freeness-fail residual (i) (15.170/216) | **OPEN** — fatal gap: \(R_{ke}\le2p\not\Rightarrow R(m_4)\le2p\) |
| **E(1) / L = ½** | **OPEN** |

## Setup

\[
m_n=\min_{a_{ij}=\pm1}\max_{x=\pm1}\Bigl|\sum_{i<j}a_{ij}x_ix_j\Bigr|,\qquad
\alpha_n=m_n/n^{3/2}.
\]

## Path (if E(1) holds)

1. **Sandwich:** \(\limsup\alpha_n\le\tfrac12\); \(\liminf\alpha_n\ge1/\pi\).  
2. **Paley family** \(n_k=p_k^2+1\): \(\Phi(C)=\tfrac12 n\sqrt{n-1}\), \(\rho=1\).  
3. **E(1):** \(m_n\ge\Phi(C)-2\) on that family ⇒ \(\alpha_{n_k}\to\tfrac12\).  
4. **Denseness:** \(n_{k+1}/n_k\to1\) ⇒ \(\lim\alpha_n=\tfrac12\).

## What is solid

- **Majorization** (15.167): if mult\((\lambda_{\max})\ge d-1\) and \(\lambda_{\min}\ge6\), then  
  \(L_*=(p^4+24p^2-1)/(2(p^2-1))<2d\) for primes \(p\ge5\) ⇒ bi-tight empty (no residual/16N).  
  Checkable: `L_star_closed`, `two_d_minus_L_star` in `src/e1_gmin_m4_prop15167.py`.  
- **Farkas poly** (15.170): under dual equality need \(6/p-4\) and candidate LB \(-12k/(pn)\),  
  poly \(4p^3-6p^2-32p+18>0\) for \(p\ge5\). Algebra OK **if** the LB holds.
- **15.216 reductions (proved, not a close):** \(R_{ke}=128p/(3p^2+17)\le2p\);  
  \(\|\kappa\|_2^2/p^4\le\) thr-η-budget with gap \((p^2+1)(p^2+9)/24\);  
  under \(Q_e\le10\), only two-level \(k\in\{0,1\}\) and mass-min blocks dual-eq.

## Residual (i) status (15.216–217)

**OPEN.** 15.217 proves \(\Phi=n(n-1)(n-2)/8\) on Max+ and reduces residual-(i) dual-eq empty to
\(\|m_4\|_2^2 \le n(n-2)/4\) (equivalently \(\|\delta\|_2^2 \le R\)-room). That bound is unproved.
Census: holds at \(p=5\). Soft-close forbidden.

## Fatal gaps (named)

**Two denseness-path gaps remain.**

- **Residual (i) OPEN (15.216–217):** Old claim \(R_{ke}\le2p\Rightarrow R\le2p\) was false (AI-test BLOCK). **15.217 repair:** exact \(\Phi=\langle m_4,\kappa\rangle=n(n-1)(n-2)/8\) on Max+ gives the correct criterion \(R\le2p\Leftrightarrow\|m_4\|_2^2\le n(n-2)/4\). That \(m_4\) bound is still unproved Max+-free (census OK at \(p=5\)). Predicates residual_i/type_I/gsum all **False**.

- **Residual (ii) affine branch CLOSED (15.179):** freeness-fail affine \(f_e=3-S\) on \(S\in\{2,4\}\) freezes \(S_H\equiv3\Rightarrow k=3p-1\); impossible for \(k\ge3p\); fail-eq empty under bi-tight. No Gsum LB.
- **Residual (ii) full OPEN (15.193):** freeness-fail does **not** force \(S\in\{2,4\}\) and \(f_e=3-S\). Multi-level (ii-a) and non-affine two-level (ii-b) remain. Predicate `deep_s2_freeness_fail_k_ge_3p_ND_closed()=False`. **Blocks E1/L.**
- Shipped: residual (i) open; residual (ii) affine closed / full open; E1/L **OPEN**. Soft-close forbidden.

## Modules

| Role | Path |
|------|------|
| Bi-tight majorization | `src/e1_gmin_m4_prop15167.py` |
| Gsum hinge + residual (i) | `src/e1_gmin_m4_prop15170.py` |
| Residual (ii) structure + freeze | `src/e1_gmin_m4_prop15171.py`, `prop15179.py` |
| Residual (ii) exhaustiveness audit | `src/e1_gmin_m4_prop15193.py` |
| Residual (i) K₄ thr path (OPEN) | `src/e1_gmin_m4_prop15215.py`, `prop15216.py` |
| Residual (i) Φ identity + m4₂ reduction | `src/e1_gmin_m4_prop15217.py` |
| Residual (i) min-distance + Q(K₄) | `src/e1_gmin_m4_prop15197.py` |
| Residual (i) Wick_hi reductions | `src/e1_gmin_m4_prop15198.py` |
| Main wire | `src/e1_main_chain_status.py` |

## What would close L

1. Residual (i): prove \(R(m_4)\le2p\) (control \(\delta\)), or \(K_4\le\mathrm{Wick}_{hi}\), or \(|\mu|\le2/n\) / ker-box / free-e_sc.  
2. Residual (ii) full: exhaustiveness lemma (freeness-fail \(\Rightarrow S\in\{2,4\}\) and \(f_e=3-S\)), or separate ND for multi-level / non-affine freeness-fail.

Until both: do **not** claim lim α_n = ½.

## AI-test questions (use this file only)

**Q1 (residual i):**  
> Is residual (i) (Type I freeness-fail dual-eq empty for all primes \(p\ge5\)) essentially proved by this package?  
> **Scored 2026-08-07:** NOT_PROVED (2/2 cold BLOCK; fatal Rayleigh gap).  

**Q2 (L — honesty):**  
> Is lim α_n = ½ proved by this package?  
> Answer: **not proved** (residual i open + full residual ii open).
