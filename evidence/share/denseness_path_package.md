# Denseness-path package: lim α_n (honest status)

**Audience:** independent AI / human check (Paata AI-test).  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-13 evening (residual (ii) ND closed; residual (i) still OPEN; post-4d89353 residual-(i) attacks did not close)

## Verdict

**\(L=\lim\alpha_n=\tfrac12\) is NOT proved.**  
Residual **(ii)** ND is **CLOSED** (affine 15.179 + (ii-b) 15.236 + (ii-a) 15.237). Residual **(i)** is **OPEN** (no Max+-free \(|\mu|\le 1/(2p)\) / dual-eq hinge). Only residual (i) blocks E(1)/L.

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** |
| Denseness (ratio-dense family ⇒ global liminf/limsup) | **Proved** (Prop 6.1) |
| Bi-tight empty via majorization (15.167) | **Candidate** (algebra OK if mult≥d−1 and λ_min≥6) |
| Residual (ii) affine dual two-level (15.179) | **CLOSED** — freeze-to-tight, no Gsum |
| Residual (ii-b) ND (15.236) | **CLOSED** |
| Residual (ii-a) ND (15.237) | **CLOSED** |
| Residual (ii) full | **CLOSED** (`residual_ii_full_closed=True`; exhaustiveness still False, not required) |
| Type I freeness-fail residual (i) (15.170/216) | **OPEN** — \(|\mu|\le 1/(2p)\) / dual-eq empty unproved |
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

**15.232–235 structure (not a close):** intersection split of \(\bar R_4\); k=3 linear, k=2 bilinear (per2), k=1 Laplace/trilinear, k=0 S₄ cycle-type + 4-cycle inverse pairing. Unsigned every layer dead except k=3 for \(p\ge89\). Predicates stay False.

**Post-4d89353 residual-(i) (not a close):** switching Max+\(\leftrightarrow\)Max−; far-sum; Comm-repair dual needs ker=sc; envelope holds on all \(|\kappa|=1\) at \(p=5,7\); \(\mu\) is a CR-class function. Dead list and viable hinges: `evidence/SESSION_HANDOFF_2026-08-13_residual_i.md`. No predicate flip.

## Fatal gaps (named)

**One denseness-path gap remains.**

- **Residual (i) OPEN (15.216–217):** Old claim \(R_{ke}\le2p\Rightarrow R\le2p\) was false (AI-test BLOCK). **15.217 repair:** exact \(\Phi=\langle m_4,\kappa\rangle=n(n-1)(n-2)/8\) on Max+ gives the correct criterion \(R\le2p\Leftrightarrow\|m_4\|_2^2\le n(n-2)/4\). That \(m_4\) bound (equivalently \(|\mu|\le 1/(2p)\)) is still unproved Max+-free (census OK at \(p=5,7\)). Predicates residual_i/type_I/gsum all **False**.

- **Residual (ii) CLOSED (15.179+15.236+15.237):** affine freeze + (ii-b) Max− dichotomy + (ii-a) L²=L pair-span. Exhaustiveness (15.193) remains False and is not required. `residual_ii_full_closed=True`, `deep_s2=True`.
- Shipped: residual (i) open; residual (ii) closed; E1/L **OPEN**. Soft-close forbidden.

## Modules

| Role | Path |
|------|------|
| Bi-tight majorization | `src/e1_gmin_m4_prop15167.py` |
| Gsum hinge + residual (i) | `src/e1_gmin_m4_prop15170.py` |
| Residual (ii) structure + freeze | `src/e1_gmin_m4_prop15171.py`, `prop15179.py` |
| Residual (ii) exhaustiveness audit | `src/e1_gmin_m4_prop15193.py` |
| Residual (ii-b) ND | `src/e1_gmin_m4_prop15236.py` |
| Residual (ii-a) ND | `src/e1_gmin_m4_prop15237.py` |
| Residual (i) K₄ thr path (OPEN) | `src/e1_gmin_m4_prop15215.py`, `prop15216.py` |
| Residual (i) Φ identity + m4₂ reduction | `src/e1_gmin_m4_prop15217.py` |
| Residual (i) min-distance + Q(K₄) | `src/e1_gmin_m4_prop15197.py` |
| Residual (i) Wick_hi reductions | `src/e1_gmin_m4_prop15198.py` |
| Main wire | `src/e1_main_chain_status.py` |

## What would close L

1. Residual (i) only: prove \(R(m_4)\le2p\) (control \(\delta\)), or \(K_4\le\mathrm{Wick}_{hi}\), or \(|\mu|\le 1/(2p)\) / \(|\mu|\le2/n\) / dual-eq empty.

Until residual (i) is closed: do **not** claim lim α_n = ½. Residual (ii) is already closed.

## AI-test questions (use this file only)

**Q1 (residual i):**  
> Is residual (i) (Type I freeness-fail dual-eq empty for all primes \(p\ge5\)) essentially proved by this package?  
> **Scored 2026-08-07:** NOT_PROVED (2/2 cold BLOCK; fatal Rayleigh gap).  

**Q2 (L — honesty):**  
> Is lim α_n = ½ proved by this package?  
> Answer: **not proved** (residual i open; residual ii closed).
