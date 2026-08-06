# Denseness-path package: lim α_n (honest status)

**Audience:** independent AI / human check (Paata AI-test).  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-06

## Verdict

**\(L=\lim\alpha_n=\tfrac12\) is NOT proved.**  
Candidate path exists; one **fatal hinge** blocks residual (i)/(ii) and therefore E(1).

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** |
| Denseness (ratio-dense family ⇒ global liminf/limsup) | **Proved** (Prop 6.1) |
| Bi-tight empty via majorization (15.167) | **Candidate** (algebra OK if mult≥d−1 and λ_min≥6) |
| Type I / deep freeness-fail ND (15.170–171) | **OPEN** — needs proved Gsum disj LB |
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

## Fatal gap (named)

**Disjoint Gsum lower bound \(\mathrm{Gsum}_{ab}\ge -12/(pn)\) is not proved for general \(p\ge5\).**

- Residual (i)/(ii) Farkas needs that LB for all primes.  
- Adjacent edges: Gsum = 0 is proved.  
- Disjoint: LB was justified as an association-scheme minimum; **Prop 15.158 proves Max+ is not an IP association scheme.**  
- Census: LB matches at \(p=5\) (\(-6/65=-12/(5\cdot26)\)). Census ≠ general proof.  
- Shipped predicate: `gsum_disj_lb_proved_general() → False`  
  (`src/e1_gmin_m4_prop15170.py`). Consequently  
  `type_I_k_3p_minus_2_closed_general()`,  
  `deep_s2_freeness_fail_k_ge_3p_ND_closed()`,  
  `e1_closed_general()`,  
  `L_status` all report **OPEN**.

## Modules

| Role | Path |
|------|------|
| Bi-tight majorization | `src/e1_gmin_m4_prop15167.py` |
| Gsum hinge + residual (i) | `src/e1_gmin_m4_prop15170.py` |
| Residual (ii) | `src/e1_gmin_m4_prop15171.py` |
| Max+ not scheme | `src/e1_gmin_m4_prop15158.py` |
| Main wire | `src/e1_main_chain_status.py` |

## What would close L

A Max+-free proof that for all primes \(p\ge5\) and disjoint edges \(a\neq b\),  
\(\mathrm{Gsum}_{ab}\ge -12/(p(p^2+1))\) (or any LB strong enough for Farkas),  
**or** a different proof of residual (i)/(ii) that avoids that LB.

Until then: do **not** claim lim α_n = ½.

## AI-test question (use this file only)

> Is lim α_n = ½ proved by this package?  
> Answer: essentially proved / not proved, and name any fatal gap.
