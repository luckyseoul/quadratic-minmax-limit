# Matching spike criterion via pairing / Γ-form

**Date:** 2026-07-27  
**Status:** Reformulation **proved**; exact criterion certified on the stored \(n=26\) Max-cover and 20 random perfect matchings.  
**Existence of \(\lim\alpha_n\) remains OPEN** (no forall-\(M\) proof for \(p\ge5\); no matching dichotomy; no \(k_\star\) bound).

## Setup

Paley conference \(C\) of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), perfect matching
\(M=\{e_a=(i_a,j_a)\}_{a=1}^{m}\) with \(m=n/2\). Write

\[
S_M(x)=\sum_{a=1}^{m} C_{e_a}\,x_{i_a}x_{j_a},\qquad
Q_C(x)=\sum_{i<j}C_{ij}x_ix_j=S_M(x)+R_M(x).
\]

**Spike criterion (Prop 15.30a).** If \(\max\{Q_C:S_M=-p\}\ge\Phi-2p\), then
\(\Phi(C\oplus M)\ge\Phi(C)\). Equivalently

\[
\max\bigl\{R_M(x):S_M(x)=-p\bigr\}\;\ge\;\Phi-p=p(m-1).
\]

## Pairing coordinates

For each matching edge \(e_a=(i_a,j_a)\) set free absolute sign \(u_a=x_{i_a}\) and
relative sign \(z_a=x_{i_a}x_{j_a}\). Then \(S_M=\sum_a c_a z_a\) with \(c_a=C_{e_a}\),
independent of \(u\), and

\[
R_M=\sum_{1\le a<b\le m}\gamma_{ab}(z_a,z_b)\,u_a u_b,
\]

where, writing \(e_a=(i,j)\), \(e_b=(k,l)\),

\[
\gamma_{ab}(z_a,z_b)
=
C_{ik}+C_{il}z_b+C_{jk}z_a+C_{jl}z_a z_b
\in\{-4,-2,0,2,4\}.
\]

(The four \(K_{2,2}\) patterns on Paley \(n=10\): \((-2,2,2,2)\), \((0,0,0,4)\),
\((-2,-2,-2,2)\), \((-4,0,0,0)\) and sign flips.)

Thus for fixed \(z\), maximising \(R_M\) over free signs is the pure degree-2 problem

\[
\max_{u\in\{\pm1\}^m}\tfrac12\,u^\top\Gamma(z)\,u
\qquad\bigl(\Gamma_{ab}(z)=\gamma_{ab}(z_a,z_b),\;\Gamma_{aa}=0\bigr).
\]

**Criterion reformulation.**

\[
\max_{\substack{z\in\{\pm1\}^m\\ \langle c,z\rangle=-p}}
\max_{u\in\{\pm1\}^m}\tfrac12\,u^\top\Gamma(z)\,u
\;\ge\;p(m-1).
\]

## Certified facts

### \(p=3\) (\(n=10\), \(m=5\), need \(\max R\ge12\))

| Family | \(\max R\) on \(S=-3\) | Spike criterion |
|--------|------------------------:|:----------------|
| 144 undercutters | **10** | False |
| 801 non-undercutters | \(\ge12\) | True |

On undercutters at the maximising \(z\): \(\lambda_{\max}(\Gamma)\approx5.87\), continuous
bound \(\tfrac m2\lambda_{\max}\approx14.68>12\), but discrete max is only \(10\).
Integrality gap kills the continuous bound at \(p=3\).

Also: undercutters **never** meet the second \(Q\)-level \(\{Q=\Phi-2p=9\}\) on
\(\{S=-p\}\); their \(S\)-values on that level lie in \(\{\pm1,\pm3\}\) only.
Clique-flip arithmetic for \(|F|=p=3\) requires admissible \(S_M(y)\equiv3\pmod4\), which
undercutters lack (Prop 15.31).

### \(p=5\) (\(n=26\), \(m=13\), need \(\max R\ge60\))

| Family | \(\max R\) on \(S=-5\) | Notes |
|--------|------------------------:|-------|
| Stored Max-cover (`e1_n26_matching_cover_example.json`) | **60** (tight) | \(\max Q=55=\Phi-2p\); MITM \(\Phi(C\oplus M)=\Phi(C)\) |
| 20 random perfect matchings | \(\in\{60,70\}\), all \(\ge60\) | criterion True |

At the cover optimum: \(\lambda_{\max}(\Gamma)\approx10.55\), continuous
\(\tfrac m2\lambda_{\max}\approx68.6>60\), discrete hits **exactly** the threshold.

Identity \(\mathbb E[R_M\mid S_M]=0\) holds on the cube (matching vs off-matching
Walsh characters are orthogonal).

## Design constants for clique-flip route (unchanged)

For \(p=5\): 390 Seidel-consistent \(p\)-sets; exactly 60 Max\(_{+}\) extensions each;
exactly **90** \(W\)-cliques of size \(p\) per Max\(_{+}\) vector (uniform). Average is
not enough for a soft existence claim: need a transversal clique with
\(\sum\chi_i=3\) on some \(y\) with \(S_M(y)=1\).

## Mod-4 dichotomy (Prop 15.32b–c, proved)

On \(\mathrm{Max}_{+}\), \(\prod_v y_v=(-1)^{p(p-1)/2}\) is constant, hence for any perfect
matching \(S_M\bmod 4\) is constant on \(\mathrm{Max}_{+}\).

| \(p\bmod 4\) | \(\pi=\prod y\) | \(-p\bmod 4\) | Natural maximiser level |
|-------------:|----------------:|--------------:|:------------------------|
| 1 (e.g. 5) | \(+1\) | 3 | class \(r=3\) |
| 3 (e.g. 3,7) | \(-1\) | 1 | class \(r=1\) |

At \(p=5\), census of 500 random matchings: every class-\(r=3\) matching attains
\(S_M=-5\); every class-\(r=1\) matching attains \(S_M=-3\) (1-bit route). Max-covers
are class \(r=1\) with \(\min S=1\) (1-bit blocked; clique-flip required).

## Census at \(p=5\) (`e1_gamma_forall_census.json`)

| Object | Result |
|--------|--------|
| need \(\max R\) | 60 |
| 2 Max-covers (known + SA) | \(\max R=60\) tight; MITM \(\Phi=\Phi(C)\) |
| 6 further SA Max-covers (400 seeds) | all \(\max R=60\), MITM \(\Phi=\Phi(C)\) (`e1_gamma_cover_batch.json`) |
| 80 random PMs | \(\max R\in\{60,70\}\), min 60 |
| SA minimise \(\max R\) (12 seeds) | lowest 60 |
| Counterexample to criterion | **none found** |

Scripts: `src/e1_gamma_criterion_census.py`, `src/e1_gamma_cover_batch.py`, `maxR_matching_level` in `minmax_quadratic.py`.

## Correction (Prop 15.33): criterion is not forall-\(M\)

A certified perfect matching at \(p=5\) has \(\max R=54<60\) (criterion **fails**) but is a
**non-cover** (\(\min_{\mathrm{Max}_{+}}S=-1\)) with MITM \(\Phi=75>\Phi(C)\). See
`e1_criterion_fail_no_undercut.json`.

By Prop 15.33, **every non-cover raises \(\Phi\) by \(\ge2\)**. Only Max-covering matchings can
undercut. The spike criterion remains a useful sufficient test on covers; it is not necessary
for non-undercut, and it is not true for every matching.

## What this does **not** prove

- Every Max-covering perfect matching has \(\Phi(C\oplus M)\ge\Phi(C)\) when \(p\ge5\)
- Clique-flip for every Max-cover
- Matching dichotomy / control of non-matching undercutters
- \(k_\star=O(n^{3/2})\)
- Existence of \(\lim\alpha_n\)

**Do not mark Main Theorem settled.**
