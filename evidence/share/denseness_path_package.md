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

**Pointwise disj Gsum lower bound strong enough for dual-equality Farkas is not proved for general \(p\ge5\).**

- **Correct Farkas threshold (Prop 15.172):** dual equality is impossible if every off-diagonal Gsum entry satisfies \(\mu>-2/p\).  
  (Borderline \(\mu=-2/p\) saturates and needs a separate non-saturation argument.)
- Adjacent edges: Gsum \(=0\) **proved**. Avg disj entry \(=2/(n-3)\) **proved**. Triangular-scheme matrix \(G_0\) with those constraints is PSD **proved**.
- Candidate scheme LB \(\mathrm{Gsum}_{ab}\ge-12/(pn)\): **false at \(p=3\)** (census min \(=-2/p<-12/(pn)\)); at \(p=5\) it holds and is tight. **Prop 15.158:** Max+ is not an IP-scheme, so scheme-min is invalid as a general proof.
- Shipped: `gsum_disj_lb_proved_general() → False`  
  (`src/e1_gmin_m4_prop15170.py`, structure in `prop15172.py`). Consequently  
  residual (i)/(ii), `e1_closed_general()`, `L_status` all **OPEN**.

## Modules

| Role | Path |
|------|------|
| Bi-tight majorization | `src/e1_gmin_m4_prop15167.py` |
| Gsum hinge + residual (i) | `src/e1_gmin_m4_prop15170.py` |
| Residual (ii) | `src/e1_gmin_m4_prop15171.py` |
| Gsum structure / Farkas threshold | `src/e1_gmin_m4_prop15172.py` |
| Gsum vector structure (ξ, stars, K) | `src/e1_gmin_m4_prop15173.py` |
| Max+ not scheme | `src/e1_gmin_m4_prop15158.py` |
| Main wire | `src/e1_main_chain_status.py` |

## What would close L

A Max+-free proof that for all primes \(p\ge5\) and disjoint edges,  
\(\mathrm{Gsum}_{ab}>-2/p\) (or \(\ge-12/(pn)\), or any \(\mu>-2/p\)),  
**or** a different proof of residual (i)/(ii) that avoids dual-equality Farkas.

Until then: do **not** claim lim α_n = ½.

## AI-test question (use this file only)

> Is lim α_n = ½ proved by this package?  
> Answer: essentially proved / not proved, and name any fatal gap.
