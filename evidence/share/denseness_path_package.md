# Denseness-path package: lim α_n (honest status)

**Audience:** independent AI / human check (Paata AI-test).  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-06

## Verdict

**\(L=\lim\alpha_n=\tfrac12\) is NOT proved.**  
Candidate path exists; residual (ii) affine branch is closed, but **full residual (ii)** is not audit-complete; residual (i) Gsum LB still blocks E(1).

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** |
| Denseness (ratio-dense family ⇒ global liminf/limsup) | **Proved** (Prop 6.1) |
| Bi-tight empty via majorization (15.167) | **Candidate** (algebra OK if mult≥d−1 and λ_min≥6) |
| Residual (ii) affine dual two-level (15.179) | **CLOSED** — freeze-to-tight, no Gsum |
| Residual (ii) full / exhaustiveness (15.193) | **OPEN** — multi-level / non-affine freeness-fail |
| Type I freeness-fail residual (i) (15.170) | **OPEN** — needs proved Gsum disj LB / \(\|\mu\|\le2/n\) |
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

## Fatal gaps (named)

**Two denseness-path gaps remain.**

- **Residual (ii) affine branch CLOSED (15.179):** freeness-fail affine \(f_e=3-S\) on \(S\in\{2,4\}\) freezes \(S_H\equiv3\Rightarrow k=3p-1\); impossible for \(k\ge3p\); fail-eq empty under bi-tight. No Gsum LB.
- **Residual (ii) full OPEN (15.193):** freeness-fail does **not** force \(S\in\{2,4\}\) and \(f_e=3-S\). Multi-level (ii-a) and non-affine two-level (ii-b) remain. Predicate `deep_s2_freeness_fail_k_ge_3p_ND_closed()=False`.
- **Residual (i) OPEN:** freeness edge \(e\notin G\), \(|G|=k=3p-2\). Dual equality forces
  \((\mathrm{Gsum}\,x)_e=6/p-4\). Threshold \(\mu_*=(6/p-4)/k\); sufficient \(\mathrm{Gsum}\ge-1/p\) or \(|\mu|\le2/n\) (15.176–192).
- **15.194–196 row/mass/spectral:** pure \(N_e\) fails Paley. Mass-corrected criterion proved; census blocks p=5,7. Spectral: \(Q_e\le 2(n+\lambda_2(n-2))/(n-1)-4\); **\(Q_e\le10\) would close residual (i)** (checked p=5..47). Census \(Q\approx8.17,6.69\). Preferred: Max+-free \(Q_e\le10\) or \(\lambda_2\le6+5/(n-2)\).
- **15.197–201 K₄ / ker path:** \(d_H\ge p+1\); Wick_hi≤thr; frame \(K_4\); **C−2/n∈ker** (15.200). **Proved:** free-e max κ_e≤α(n−2) ⇒ dual-eq empty for p≥5 (15.201). Census free-e max≤α(n−2) at p=5,7; p=7 max|μ|=109/2863<2/n. **OPEN:** free-e max≤α(n−2) general, or |μ|≤2/n, or K₄≤Wick_hi.
- Adjacent Gsum \(=0\), avg disj \(=2/(n-3)\), \(G_0\) PSD **proved**. Candidate \(-12/(pn)\) not general (15.158).
- Shipped: `gsum_disj_lb_proved_general()→False`; residual (ii) affine closed / full open; E1/L **OPEN**.
- Support (15.186–188): \(|\varphi|\le2(p-2)\); \(\mu_{\mathrm{part}}\) majorant \(\le1/(2p)\); target \(|\mu|\le2/n\) sufficient for \(p\ge5\); p=5,7 Max± census beat \(1/(2p)\). **Not** a general Max+-free proof — do not treat as closed.

## Modules

| Role | Path |
|------|------|
| Bi-tight majorization | `src/e1_gmin_m4_prop15167.py` |
| Gsum hinge + residual (i) | `src/e1_gmin_m4_prop15170.py` |
| Residual (ii) structure + freeze | `src/e1_gmin_m4_prop15171.py`, `prop15179.py` |
| Residual (ii) exhaustiveness audit | `src/e1_gmin_m4_prop15193.py` |
| Residual (i) row negative-mass | `src/e1_gmin_m4_prop15194.py` |
| Residual (i) mass-corrected dual-eq | `src/e1_gmin_m4_prop15195.py` |
| Residual (i) spectral Q / λ₂ | `src/e1_gmin_m4_prop15196.py` |
| Residual (i) min-distance + Q(K₄) | `src/e1_gmin_m4_prop15197.py` |
| Residual (i) Wick_hi reductions | `src/e1_gmin_m4_prop15198.py` |
| Residual (i) frame K₄ identity | `src/e1_gmin_m4_prop15199.py` |
| Residual (i) C−2/n ker + free-e box | `src/e1_gmin_m4_prop15200.py` |
| Residual (i) α(n−2) free-e sufficient | `src/e1_gmin_m4_prop15201.py` |
| Residual (i) support constraints | `src/e1_gmin_m4_prop15180.py` |
| Gsum structure / Farkas threshold | `src/e1_gmin_m4_prop15172.py` |
| Gsum vector structure (ξ, stars, K) | `src/e1_gmin_m4_prop15173.py` |
| Hinge dictionary U/cos | `src/e1_gmin_m4_prop15174.py` |
| Gμ reformulation | `src/e1_gmin_m4_prop15175.py` |
| **Correct e∉G Farkas threshold μ_*** | `src/e1_gmin_m4_prop15176.py` |
| Max+ not scheme | `src/e1_gmin_m4_prop15158.py` |
| Main wire | `src/e1_main_chain_status.py` |

## What would close L

1. Residual (i): Max+-free \(K_4\le\mathrm{Wick}_{hi}=12n^2+48n\) (15.198: that \(\Rightarrow Q<10\) and thr for \(p\ge5\)), or \(K_4\le16n^2-14n\), or \(|\mu|\le2/n\) / ker-box / \(\mathrm{Gsum}\ge-1/p\). Pure \(N_e<4-6/p\) fails Paley census.  
2. Residual (ii) full: exhaustiveness lemma (freeness-fail \(\Rightarrow S\in\{2,4\}\) and \(f_e=3-S\)), or separate ND for multi-level / non-affine freeness-fail.

Until both: do **not** claim lim α_n = ½.

## AI-test question (use this file only)

> Is lim α_n = ½ proved by this package?  
> Answer: essentially proved / not proved, and name any fatal gap.
