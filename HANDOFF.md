# Research handoff: min-max ±1 quadratic form limit

**Status date:** 2026-07-29 (Prop 15.81 moduli-line GD at p=5 complete; general M_cand still open; L still OPEN)  
**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**Resume file:** `evidence/SESSION_HANDOFF_2026-07-30.md`  
**Problem source:** [MathOverflow 413935](https://mathoverflow.net/questions/413935) / [X prize post](https://x.com/PI010101/status/2081070728422752329)

**Compute (F17, non-negotiable):** 88 cores. Full pytest only via `./scripts/pytest_full.sh` (`W=nproc-2`, default 86). Research scripts: `src/workers.py` `require_workers()` + file-based `ProcessPoolExecutor`, `OMP_NUM_THREADS=1` per worker. **Prefer GPU (CuPy/V100)** for dense m4/Gram batches (`e1_gmin_m4_gpu.py`: one CUDA context, device argmax). **Atomic I/O when the job allows** (`src/io_atomic.py`, Wieferich-style): mmap Max+; evidence via tmp+fsync+`os.replace`; record `gpu.used` + `io` in evidence. Never full-suite / multi-minute residual work on 1 core. Never `python -c` / stdin ProcessPool.

---

## 0. One-line status

**Existence of the limit \(L=\lim\alpha_n\) remains OPEN.**  
Proved sandwich: \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\).  
Dense Paley subsequence with \(\rho=1\) (orders \(n=p^2+1\)) is proved; asymptotic optimality of \(m_n\) vs Paley \(\Phi\) is not.

**P0 start:** prove \(\lambda_2(P\odot P)\le d/(2N)\) for all primes \(p\ge5\) (Prop 15.56.6 / 15.58.3 Veronese form \(\|T(x)\|_F^2\le nN\|x\|^2\) on \(\mathbf1^\perp\)); then deep non-tight ND; only then Main Theorem. Resume: `evidence/SESSION_HANDOFF_2026-07-30.md`. **F16:** do not pin modulus \(c\) by max PSD \(g_{\min}\).

**Props 15.45–15.51 (load-bearing, do not re-derive):**
- Stars never bi-tight (wedge \(G^++G^-=0\)).
- Bi-tight level-2 empty if \(g_{\min}>-(p-2)/(p(2p-1))\) (Prop 15.47). **Certified \(p=5,7\)**.
- 1-bit spike (15.46). Deep at \(p=5\): small-\(k\) infeasible; large-\(k\) **spike** above \(\Phi(C)\).
- **Prop 15.48–15.49:** edge algebra; full CR classification of \(g_{\min}\) (cert p=3,5,7 in `e1_gmin_cr_classify.json`); **uniform LB candidate** \(g_{\min}\ge-(p-2)/(2p^2)\) algebraically beats the bi-tight threshold for all odd \(p>2\) and holds at certified \(p=5,7\) (`e1_gmin_uniform_lb.json`).
- **Prop 15.50 (proved):** Max+ conditional means given two coordinates equal the Gaussian frame interpolant \(\Sigma_{*S}\Sigma_{SS}^{-1}y_S\); conditional \(\mathbb E[y_ky_l\mid y_i,y_j]=\alpha+\delta y_iy_j\); disj average \(G=1/(p^2-2)\). Fréchet on cond cov **too weak** for \(L(p)\) (do not reopen). Evidence: `e1_gmin_cond_mean.json`.
- **Prop 15.51 (proved equiv form):** \(a=(1+pG)/(p+1)\); \(g_{\min}\ge T(p)\Leftrightarrow\min a(e')\ge1/(2p-1)\); disj sum deterministic on slice; residual Loewner cert p=5,7. Evidence: `e1_gmin_structure.json`.
- **Prop 15.52 (proved sum + moduli sketch):** \(\mathbf1^\top y=(p+1)y_\infty\) on Max+; m4 evec system sketch.
- **Prop 15.53 (proved pairing + cert moduli pin):** \(g_{\min}=-\max|m_4|\) on \(|\kappa|=1\); refined classes (37 at \(p=5\)); nullity-1; \(\mathrm{Tr}(G^2)\) pin recovers \(-3/65\); \(\mathbb E[\mathrm{dot}^2]\) from 2-design **proved**. `src/e1_gmin_moduli.py`.
- **Prop 15.54 (moduli \(c\)-calculus + wedge + a-slice + deep refresh):** wedge \(G=\pm1/p\); \(g_{\min}(c)\) pin; a_min p=5,7; deep covers spike at p=5.
- **Prop 15.55–15.56 (tight obstruction + star/cycle/Schur reduction):** If \(\lambda_{\max}(G)=n/2\) simple ⇒ no Max+-tight size-\(2p\) cover. **Proved:** stars map to \(\mathbf1\); \(\lambda_{\max}=\max(n/2,\lambda_{\mathrm{cycle}})\); \(\lambda_{\mathrm{cycle}}=2N\lambda_2(P\odot P)\); avg cycle eig \(<n/2\) for \(p\ge5\). **Certified gap** \(\lambda_2(P\odot P)\le d/(2N)\) at \(p=5,7\) (fails \(p=3\)). Evidence: `e1_gmin_spectral.json`.
- **Prop 15.57 (triangle bound):** \(|\mathrm{ft}|\le2p\), \(\mathbb E[\|By\|^2]\le4\), crude \(\lambda_{\max}(G)\le n\) (factor-2 short of \(n/2\)).
- **Prop 15.58 (Veronese residual form; 2026-07-30):** **Proved:** Max+\(\subset V_+\) (\(Cy=py\)); \((P\odot P)\mathbf1=\alpha\mathbf1\), \(\lambda_{\max}(P\odot P)=\alpha\) (Perron); gap \(\Leftrightarrow\|T(x)\|_F^2\le nN\|x\|^2\) for \(x\perp\mathbf1\) with \(T=\sum x_a y_ay_a^\top\); equiv. \(\lambda_2(W)\le\lambda_1(W)/2\) for \(W_{ab}=(u_a\cdot u_b)^2\). **Certified maximiser:** \(\mathrm{ft}=2p\), \(By\in V_+\), \(\|B_+\|_F^2=2\). Evidence: `e1_gmin_gap_probe.json`.
- **Prop 15.59 (center / rank / two-moment; 2026-07-30):** **Proved:** Max+ centrally symmetric \(\Rightarrow\sum y_a=0\Rightarrow P\mathbf1=0\). **Certified:** \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\) at \(p=3,5\); two-moment on \(W\) forces gap at \(p=7\) only. Evidence: `e1_gmin_veronese.json`.
- **Prop 15.60 (antipodal → projective ENTF; 2026-07-30):** **Proved:** \(T(x)=T(s)\); gap \(\Leftrightarrow\lambda_2(W^{(\mathrm{proj})})\le m/(2d)\); 2×sphere \(\Rightarrow\) gap for \(d\ge6\). Evidence: `e1_gmin_projective.json`.
- **Prop 15.61 (16N \(\Rightarrow\lambda_{\mathrm{cycle}}\le8\Rightarrow\) gap for \(p\ge5\); 2026-07-30):** **Proved equivalences:** \(Q(B):=\sum_y(y^\top By)^2\le16N\|B\|_F^2\) (Tr0 on \(V_+\)) \(\Leftrightarrow\lambda_2(W)\le4N/d^2\Leftrightarrow\lambda_{\mathrm{cycle}}\le8\). **Proved algebra:** \(\lambda_{\mathrm{cycle}}\le8\) and \(p\ge5\) (\(d\ge13>8\)) \(\Rightarrow\lambda_{\max}(G)=n/2\) simple \(\Rightarrow\) bi-tight empty (Prop 15.55). **Certified:** equality \(Q=16N\), \(\lambda_{\mathrm{cycle}}=8\) at \(p=3\); strict at \(p=5\) (\(Q/(16N)=11/13\)), \(p=7\) (\(\approx0.66\)). Evidence: `e1_gmin_16n.json`.
- **Prop 15.62 (typeA+wedge \(=6N\); \(Q=6N\|B\|_F^2+Q_4\); 2026-07-30):** **Proved:** for zero-diag \(B=P_+BP_+\), same-edge contrib \(=2N\|B\|_F^2\), wedge contrib \(=4N\|B\|_F^2\) (via \(\Sigma=2P_+\) on \(V_+\)), hence \(Q=6N\|B\|_F^2+Q_4\) with \(Q_4=Be^\top\mathrm{Gu}_{\mathrm{disj}}Be\) only. **Equivalence:** 16N \(\Leftrightarrow Q_4\le10N\|B\|_F^2\). **Certified** multi-seed p=3,5,7 (`e1_gmin_typeA_wedge.json`); p=3 saturates \(Q\equiv16N\) on whole space. Unrestricted \(\|\mathrm{Gu}_{\mathrm{disj}}\|_{\mathrm{op}}\) too large.
- **Prop 15.63 (H\(\Rightarrow\)16N; 2026-07-30):** **Proved algebra:** \(H(p):=(p+2)^2/d\le5\) for primes \(p\ge3\) (eq only \(p=3\)). **Proved:** \(\mathrm{ray}:=Q_4/(2N)\le H(p)\) \(\Rightarrow Q_4\le10N\) \(\Rightarrow\) 16N \(\Rightarrow\) bi-tight empty for \(p\ge5\). **Certified:** H holds p=3,5,7 with equality at p=3 and p=5 maximiser (\(\mathrm{ray}=49/13\)). Evidence: `e1_gmin_q4_ub.json`.
- **Prop 15.64 (dual \(\Phi\) form; residual reformulation; 2026-07-30):** **Proved:** \(\max Q=N\lambda_{\max}(\Phi|_{\mathcal Z})\); Wick baseline 8; H \(\Leftrightarrow\) residual\(\le(p+1)(p+7)/d\). **Certified** gen.eig p=3,5,7.
- **Prop 15.65 (κ spectrum; clean \(\lambda_2(P\odot P)\le4/N\); boolean essential; 2026-07-30):** **Proved:** \(16N\Leftrightarrow\lambda_2(P\odot P)\le4/N\). **Certified** κ spectrum + random-projector counterexample. Evidence: `e1_gmin_cumulant.json`.
- **Prop 15.66 (zero-diag freeness; pairing residual criterion; 2026-07-30):** **Proved:** \(\Phi\)-maximiser on \(\mathrm{Sym}_0\) automatically has ambient zero diag (\(\lambda_{\max}(\Phi|_{\mathcal Z})=\lambda_{\max}(\Phi|_{\mathrm{Sym}_0})\)). **Proved algebra:** for \(p\ge5\), \(|m_4-\kappa/p^2|\le(p-4)/(2p^2)\) on all \(|\kappa|=1\) 4-sets \(\Rightarrow g_{\min}\ge L(p)>T(p)\). **Certified:** \(g_{\min}\ge L\) at p=5,7; residual triangle holds p=7, **fails p=5** (large resid on small-\(|m_4|\) classes). Exact \(g_{\min}=-3/65,-109/2863\). Evidence: `e1_gmin_m4_residual.json`. **OPEN residual:** prove \(|m_4|\le(p-2)/(2p^2)\) on \(|\kappa|=1\) (or \(\lambda_2(P\odot P)\le4/N\)) for all primes \(p\ge5\). Bi-tight **not** closed for general \(p\).
- **Prop 15.67 (master m4 identity; same-sign Ext; census p=5,7; 2026-07-30):** **Proved:** \(\sigma_{sum}=4\kappa\) (K4 exhaustion); \(m_4=\kappa/p^2+\mathrm{Ext}/(4p)\) with \(\mathrm{Ext}=Tm_4\); same-sign \(|\mathrm{Ext}|\le2(p-4)/p\Rightarrow|m_4|\le L_{abs}\). **Certified W=86:** full \(|\kappa|=1\) census p=5,7 both satisfy direct \(|m_4|\le L\) and Ext criterion. Evidence: `e1_gmin_m4_proof.json`, `src/e1_gmin_m4_proof.py`, `src/workers.py`. **OPEN residual:** same as 15.66.4 for general primes \(p\ge5\).
- **Prop 15.68 (Tκ calculus; residual source on |κ|=3; resolvent reduction; 2026-07-30):** **Proved (conference C² + K4 exhaustion):** Tκ=0 on every |κ|=1 4-set; Tκ∈{±24} on |κ|=3; residual RHS Tκ/p² supported only on |κ|=3. **Proved algebra:** |m4|≤L ⇔ |ρ|≤(p-4)/(2p²) (same-sign) ⇔ resolvent gain |κ|=3→|κ|=1 ≤(p-4)/48; candidate |m4|≤(p-2)/(p(2p+3))≤L (sharp p=5). **Certified:** Paley degrees d3=p²-5, d1=3p²-7 at p=3,5,7. Evidence: `e1_gmin_m4_tkappa.json`. **OPEN:** resolvent gain / candidate for all primes p≥5.
- **Prop 15.71 (κ-stratum counts for any conference matrix; 2026-07-29):** **Proved:** ∑_S κ(S)² = n(n-1)(n-2)(n-5)/8 via C² wedge sum Σ=-n(n-1)(n-2) + K4 ratio-8 exhaustion ⇒ n1=n(n-1)(n-2)²/32, n3=n(n-1)(n-2)(n-6)/96 (Max+-free). **Certified W=86:** full κ census p=3,5,7,11 matches; d1,d3 constancy samples p=3,5,7,11. Evidence: `e1_gmin_m4_stratum.json`, `src/e1_gmin_m4_stratum.py`. Pins resolvent source size (# of |κ|=3 coords). **OPEN residual:** still |m4|≤M_mid/L_abs for all p≥5; lim α_n OPEN.
- **Prop 15.72 (resolvent-gain calculus; reverse degrees; type6 pin; 2026-07-29):** **Proved:** gain≤(p-4)/48 ⇔ same-sign |ρ|≤(p-4)/(2p²) ⇔ |m4|≤L on |κ|=1; Tκ/κ∈{±8} on |κ|=3; reverse degrees d1^(3)=3(p²-1), d3^(3)=p²-9 from n1,n3,d3 handshaking. **Certified W=86:** separate κ-weighted vanishing sk1=sk3=0 on all |κ|=1 at p=3,5,7; reverse deg constancy; type6 Max+-free resolvent |m|≤L and gain≪budget at p=5,7; true Max+ gains 0.0064, 0.036 < budget. Evidence: `e1_gmin_m4_resolvent_gain.json`. **OPEN:** prove gain≤(p-4)/48 for all primes p≥5; lim α_n OPEN.
- **Prop 15.73 (Max+ e4 identity; Paley ∑κ; multi-prime type6; 2026-07-29):** **Proved:** e4=-p(p-1)(p+1)(p+4)/12 from Prop 15.52 + boolean expansion ⇒ ∑_S m4(S)=e4. **Certified:** ∑κ=p²(p²-1)/4 at p=3,5,7,11,13; type6 Max+-free resolvent le_L and le_mid at p=5,7,11,13 with gain≪budget (type6 not exact upper bound — p=7 true m4>type6). Evidence: `e1_gmin_m4_e4_gain.json`. **OPEN residual:** true Max+ |m4|≤M_mid/L for all p≥5; lim α_n OPEN.
- **Prop 15.74 (candidate bound; true Max+ census; signed residual; 2026-07-29):** **Proved:** M_cand algebra; gain_cand=(p²-4p-3)/(24(2p+3)); 4p r=κ(Tρ) on |κ|=1. **Certified W=86 true Max+ (mmap, not type6):** max|m4|≤M_cand at p=5 (sharp 3/65) and p=7 (109/2863<5/119); same-sign gains 1/156 and ≈0.036. Evidence: `e1_gmin_m4_kernel.json`. **OPEN:** prove max|m4|≤M_cand for all primes p≥5; lim α_n OPEN.
- **Prop 15.75 (one-center σ_a=2·star_a; K4 Gram spectrum; GPU cand; 2026-07-29):** **Proved (any conference):** σ_a=∑_r C_ar κ(S_{a→r})=2·star_a on |κ|=1 via C²+Tκ=0; one-center residual form; K4 local G eigenvalues in {1±m4}∪{1±m4±2/p}∪{1±m4±2√2/p} ⇒ weak |m4|≤(p-2)/p. **Certified:** σ census p=3,5,7,11; GPU CuPy/V100 full |κ|=1 m4+cand (mmap Max+, device reduce, atomic JSON) p=5,7 le_cand in ~0.3s. Evidence: `e1_gmin_m4_onecenter.json`, `e1_gmin_m4_gpu.json`. **OPEN residual:** bound ∑C_ar ρ ⇒ |m4|≤M_cand all p≥5; lim α_n OPEN.
- **Prop 15.76 (one-center degrees d1^(1),d3^(1); residual split; 2026-07-29):** **Proved:** (3p²-7)/4 and (p²-5)/4 integral for odd p; 4 d1^(1)=d1. **Certified W=86:** one-center degrees constant on every |κ|=1 set/centre at p=3,5,7,11. **Proved:** abs bootstrap 4pρ≤d1 ρ+d3 R3 fails (4p−d1<0). **GPU residual:** R1,R3 on κ1/κ3 strata p=5,7; still le_cand. Evidence: `e1_gmin_m4_onecenter_deg.json`. **OPEN:** signed S1+S3 bound for M_cand all p≥5; lim α_n OPEN.
- **Prop 15.77 (star·S1≤0 structure; joint cand criterion; 2026-07-29):** **Proved:** star·(S1+S3)=pρ star−2/p²; star=+1 & S1≤0 ⇒ pρ≤2/p²+S3; same-sign ρ=(2/p²+S1+S3)/p at star=+1; cand ⇔ max joint ≤ p ρ_cand−2/p² (negative at p=5 ⇒ need strongly negative S1). **Certified GPU mmap+atomic:** star·S1≤0 on **all** |κ|=1 centres p=5,7 (strict); joint ⇒ ρ≤ρ_cand (sharp p=5); synthetic non-Max+ bump violates star·S1. Evidence: `e1_gmin_m4_S1_star.json`, `src/e1_gmin_m4_S1_star.py`. **OPEN:** prove star·S1≤0 + joint/S3 bound for all primes p≥5; lim α_n OPEN.
- **Prop 15.78 (moment form/GD; 4-set constancy; p=5 exact spectrum; 2026-07-29):** **Proved:** star·S1=E[φ]−E_Wick[φ] (φ=star f0 U1); GD ⇔ star·S1≤0; joint rewrite S1+S3=p m4−1/p−2/p² at κ=1,star=+1. **Certified W=86 pure-C:** star·τ1 constant on every |κ|=1 set p=3,5,7,11. **Certified GPU mmap+atomic:** star·S1 constant on every |κ|=1 Max+ set p=5,7; GD holds; star·S1≤0; **p=5 exact** star·S1∈{−2/65,−42/325} both <0 (⇒ star·S1≤0 **proved at p=5** by full census). Evidence: `e1_gmin_m4_S1_const.json`, `src/e1_gmin_m4_S1_const.py`. **OPEN:** prove GD for all p≥5; joint/S3 ⇒ M_cand; lim α_n OPEN.
- **Prop 15.79 (Aut-constancy of star·τ1 and star·S1; modular τ1; 2026-07-29):** **Proved:** PGL(2,p²) setwise stabilizer of any 4-set contains V4 acting regularly ⇒ transitive on the 4 points; star·τ1 Aut-equivariant ⇒ **constant on centres** (Max+-free); star·S1 Aut-equivariant (Max+ Aut-invariant) ⇒ **constant on centres**; τ1=2A−d1^(1) with d1 odd ⇒ star·τ1 odd. **Certified W=86:** constancy + t1=2A−d1 + star·τ1≡5 (mod 6) + #(values)=(p−1)/2 at p=3,5,7,11. Evidence: `e1_gmin_m4_S1_aut.json`, `src/e1_gmin_m4_S1_aut.py`. **OPEN:** prove star·S1≤0 (GD) all p≥5; joint/S3⇒M_cand; lim α_n OPEN.
- **Prop 15.80 (linear-form Wick; GD ⇔ E[ZU1]≤E_Wick[ZU1]; U1-special; 2026-07-29):** **Proved:** E[L²]=E_Wick[L²] for every linear L on Max+ (pairwise Σ only); GD ⇔ Cov(Z,U1)≤Cov_Wick with Z=star f0; GD fails for generic L (U1 κ1-support essential). **Certified GPU mmap+atomic p=5,7:** GD holds all |κ|=1 sets; U1/Wick ratio≡1; sum star·S1=−1128 (p5), −15271200/2863 (p7); sum star·τ1=ε n1 (ε=±1); generic L violation rate ~45–55%. Evidence: `e1_gmin_m4_S1_gd.json`, `src/e1_gmin_m4_S1_gd.py`. **OPEN:** prove E[ZU1]≤E_Wick all p≥5; joint/S3⇒M_cand; lim α_n OPEN.
- **Prop 15.81 (moduli-line GD criterion; p=5 complete; 2026-07-29):** **Proved form:** on nullity-1 line m=m★+c n, if max star·S1=α+βc with β>0 then GD⇔c≤c_GD=−α/β. **Certified p=5:** 37 constant-m4 classes, nullity 1, max star·S1 exact affine (β>0), c*≈−0.424<c_GD≈−0.296, max|m4|(c*)=M_cand=3/65, GD holds — **cand+GD at p=5 via moduli**. **Certified p=7:** pointwise evec holds; coarse classes 69/82 constant m4 (incomplete for nullity-1 line). Evidence: `e1_gmin_m4_S1_moduli.json`, `src/e1_gmin_m4_S1_moduli.py`. **OPEN:** refine classes all p≥5; prove c*≤c_GD generally; lim α_n OPEN.

- **Residual to \(m_n\ge\Phi-2\Rightarrow L=\tfrac12\):** (1) prove \(g_{\min}\ge L(p)\) or \(\lambda_2(P\odot P)\le4/N\) for all \(p\ge5\); (2) deep non-tight ND. **F13/F16–F18** intact.

**g_min attack status (Prop 15.49–15.50):** \(g_{\min}=-\alpha_\star\) on constant-\(m_4\) \(|\kappa|=1\) CR classes. Values \(-\tfrac13,-\tfrac3{65},-\tfrac{109}{2863}\). Matching form \(\mathbf1_M^\top G\mathbf1_M\ge9.96>4\) at \(p=5\) (no tight matching). **Dead:** \(-3/\Phi\) general LB; 4-point LP; Chebyshev; Wick-as-LB; bare \(C\)-types; affine halfspace orbit; pure deg pigeon; min-norm \(V_+\) interpolation alone; Bose–Mesner alone (span too large); plain CLT residual (overestimates m4); **plain Fréchet on conditional cov** (Prop 15.50.4). **Need:** prove uniform LB for all \(p\ge5\) (character sum / bound \(\delta\) or \(m_4\) on min CR class) + deep residual.

**New (n=10 structure):** exact optima first appear at Hamming distance **5** from Paley \(C_{10}\), and the only 5-edge undercutters are **144 perfect matchings** (of 945). Absolute gap \(\Phi-m_{10}=2\) is consistent with E(1). See `evidence/N10_STRUCTURE.md`.

**New (n=10 classification, N10-C):** those 144 matchings are exactly one \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit, equivalently the matchings that drop every Paley maximizer to \(|Q|\le13\) (six \(+\) maximizers certify). See `evidence/N10_MATCHING_CLASSIFY.md`.

**New (E(2) interval formula):** for Paley \(q\equiv1\bmod4\), interval signing has exact
\(x^\top Cx=2-8\sum_{d\le(q-1)/2}d\chi(d)\) (`evidence/E2_INTERVAL_FORMULA.md`). Constructive
\(\rho_{\mathrm{int}}\) up to \(\approx0.978\) at \(n=1622\). **Not** a full proof of \(\rho\to1\).

**New (E(1) numerics):** exact-Φ SA at Paley \(n=14,18\): no undercut of \(\Phi(C)\) found
(`evidence/e1_paley_gap.json`); \(n=10\) recovers gap \(2\); \(n=14\) exact \(k\le4\) flips stay at \(\Phi=21\);
**\(n=26\) intensive** (172×10k SA + exact rescore + k-flips): **no undercut of \(\Phi=65\)**
(`evidence/e1_n26_intensive.json`, `E1_STATUS.md`). Maximizers of \(\rho=1\) Paley satisfy
\(\mathbb E[yy^\top]=I\) (2-design; \(p=3,5\)). **Not** a proof of E(1).

**Settlement path (shortest):** E(1) alone on the dense \(\rho=1\) family \(n=p^2+1\) \(\Rightarrow L=\tfrac12\) by Prop 6.2
(E(2) not needed for that subfamily). **E(1) is the blocking open problem.**

**E(1) reduction (2026-07-27):** Prop 15.20b edge-counting Lipschitz \(\Phi(A)\ge\Phi(C)-2k\) sharpens the sparse regime (vs Frobenius \(n\sqrt k\)). **E(1) \(\Leftrightarrow\) \(k_\star=o(n^{3/2})\)** for \(\Phi\)-minimisers after switching to Paley on \(n=p^2+1\). At \(n=10\), \(k_\star=5\) (matchings), though some other \(m_{10}\)-matrices have best-\(k\ge15\) (campaign). Rigidity \(k_\star=o(n^{3/2})\) still unproved in general. See `evidence/E1_EDGE_LIPSCHITZ.md`, `evidence/E1_RIGIDITY_ATTACK.md`.

**N10-C6 (2026-07-27):** exhaustive scan of \(\binom{45}{6}\) — all **360** Hamming-6 undercutters of Paley \(C_{10}\) are single **6-cycles** with \(\Phi=13=m_{10}\). Together with N10-S (144 matching undercutters at \(k=5\)), every undercutter at the two smallest cardinalities is a path/cycle graph (\(\Delta\le2\), \(k\le n\)). Sufficient for E(1) if generalised: **path-cycle dichotomy** \(k_\star\le n\) \(\Rightarrow L=\tfrac12\). Far optima need not *contain* a sparse undercutter as a subgraph; only existence of a sparse minimiser matters for \(k_\star\). See `evidence/N10_CYCLE_UNDERCUTTERS.md`, `src/n10_cycle_undercutters.py`.

**Prop 15.24 (2026-07-27):** maximizers of any ρ=1 conference are exactly the boolean ±p-eigenvectors. Paley counts of +p boolean evecs: 12, 260, 11452 at p=3,5,7 (`evidence/BOOLEAN_EVECS_MAX.md`). Ratio #/n^{3/2} increases, so k_⋆≤|Max| is not o(n^{3/2}).

**Prop 15.25 (2026-07-27):** recursive formula \(m_n=\min_B\max_x(|Q_B(x)|+|\sum x_i|)\) over Seidel \(B\) of order \(n-1\) (star-reduction). Certified for \(n\le11\). Does not by itself close E(1) (boost \(\max_L|s|\) often \(O(1)\)).

**Prop 15.26 (2026-07-27):** matching flips of any \(\rho=1\) conference **preserve local maximality** of boolean \(+p\)-eigenvectors (\(y_i(Ay)_i\ge p-2>0\)). At \(n=10\), for all 945 perfect matchings the stronger identity \(\Phi(C\oplus M)=\max_{\mathrm{Max}}|Q_{C\oplus M}|\) holds (cube max attained on boolean eigenset). Star / high-\(\Delta\) flips destroy this. Star-reduction probe: \(f(C')=\Phi(C)\), \(\Phi(C')=m_9\), but matching undercutters give \(d_H(B,C')=12\) (sparsity not preserved). See `evidence/E1_STAR_REDUCTION_PROBE.md`. Route to matching dichotomy still open.

**Prop 15.27 (2026-07-27):** for Paley \(\rho=1\), \(\mathrm{Max}_{+}\) is a tight frame: \(\mathbb E[yy^\top]=I+C/p=2P_{+}\). Consequences: (i) fractional Max-cover LP has value exactly \(p\); (ii) **Max-Lipschitz** \(\Phi(A)\ge\Phi(C)-2k/p\) after best switch (edge lip improved by factor \(p=\sqrt{n-1}\)). E(1) demand weakens from \(k_\star=o(n^{3/2})\) to \(k_\star=o(n^2)\) — still not free (worst-case \(k\sim n^2/4\)). Existence remains **OPEN**.

**Prop 15.28 (2026-07-27):** size-\(p\) Max-covers are tight (\(S_F\equiv1\) on \(\mathrm{Max}_{+}\)). **Proved for Paley:** \(\infty\)-stars whose leaves are nonsquare-direction affine lines of \(\mathrm{AG}(2,p)\) are tight Max-covers with \(\mathbb E_{-}[S^2]=2p-1\), hence \(\max_{\mathrm{Max}_{-}}S\ge1\) and \(\Phi\ge\Phi+2\) (cannot undercut). Square directions are non-covers. Certified: all covering \(p\)-stars (any centre) spike to \(\Phi+2p\) at \(p=3,5\); all 405 size-\(p\) covers at \(n=10\) have \(\Phi\ge17\). See `evidence/E1_SIZE_P_MAXCOVER.md`. **Does not close E(1)** (larger covers / \(k_\star\) still open).

**Prop 15.29 (2026-07-27):** (i) \(n/2\) odd \(\Rightarrow\) matching scores \(S_M\) always odd; non-covering PMs raise \(\Phi\) by \(\ge2\). (ii) Undercutters are two-sided Max-covers. (iii) **Correction:** Max-covering PMs **exist** at \(n=26\); the three SA-found covers are two-sided with Max\(\pm\) height \(63=\Phi-2\) but exact MITM \(\Phi=65=\Phi(C)\) (spike). 0 undercuts. See `evidence/E1_MATCHING_COVER_SPIKE.md`. **Does not close E(1).**

**Prop 15.30 (2026-07-27):** Matching spike criterion **proved**: \(S_M=-p\) \& \(Q_C\ge\Phi-2p\) \(\Rightarrow\Phi(C\oplus M)\ge\Phi(C)\). At \(n=10\), criterion holds iff matching does not undercut (complete over 945 PMs). At \(n=26\), all tested PMs satisfy criterion. Open: prove criterion for every PM when \(p\ge5\). See `evidence/E1_MATCHING_SPIKE_CRITERION.md`. **Does not close E(1).**

**Prop 15.31 (2026-07-27):** Clique-flip sufficiency for matching Max-covers **proved**. Only \(|F|=p\) works for covers; \(p=3\) undercutters block the arithmetic; \(p=5\) forces \(S_M=1\). Design: 390 Seidel \(p\)-sets, 60 Max\(_{+}\) extensions each, \(\ge236\) transversal per matching. All SA covers clique-flip to \(\Phi(C)\). Open: \((y,F)\) for every Max-cover matching, \(p\ge5\). See `evidence/E1_CLIQUE_FLIP.md`. **Does not close E(1).**

**Prop 15.32 (2026-07-27):** Γ-pairing reformulation of the spike criterion **proved**; coordinate product \(\pi\) constant on Max\(_{+}\) **proved**; \(S_M\bmod 4\) constant on Max\(_{+}\) for every matching **proved**. Case split: maximiser / 1-bit / clique-flip according to residue class. Census at \(p=5\): \(\max R\in\{60,70\}\) on all tested PMs (covers tight at 60); 0 counterexamples among random/SA-min-\(\max R\). Open: attainment of \(-p\) or \(-p+2\) in every residue class, and clique-flip on every cover. See `evidence/E1_GAMMA_PAIRING.md`, `e1_gamma_forall_census.json`. **Does not close E(1).**

**Prop 15.33 (2026-07-27):** **Non-covers cannot undercut** (proved): if \(\min_{\mathrm{Max}_{+}}S_F\le-1\) then \(\Phi(C\oplus F)\ge\Phi+2\). Matching non-undercut therefore reduces **only** to Max-covering matchings. Spike criterion is **not necessary**: certified matching with \(\max R=54<60\) but non-cover and MITM \(\Phi=75>\Phi\) (`e1_criterion_fail_no_undercut.json`). All SA Max-covers at \(p=5\) still satisfy criterion and \(\Phi=\Phi(C)\). Open: forall Max-cover matchings, \(p\ge5\). **Does not close E(1).**

**Prop 15.34 (2026-07-27):** Matching flip algebra \(D^2=I\), \(A=C-2D\), \(A^2=(n+3)I-2(CD+DC)\) **proved**. At \(p=5\), every tested Max-cover has \(\|A\|_{\mathrm{op}}=\sqrt{41}\) exactly; random matchings larger. Open: use spectral control of Max-cover flips to force \(\Phi(A)\ge\Phi(C)\). **Does not close E(1).**

**Prop 15.35 (2026-07-27):** At \(p=5\), every Max-cover PM **must** attain \(S_M=1\) with residue \(1\bmod4\) (proved by expectation). Census of **11** SA Max-covers: all two-sided, min+max as covers, \(\mathrm{op}=\sqrt{41}\), clique-flip, \(\Phi=\Phi(C)\). 0/20k random PMs are covers. Open: forall Max-cover \(M\) when \(p\ge5\). See `evidence/e1_maxcover_full_census.json`. **Does not close E(1).**

**Prop 15.36 (2026-07-27):** Matching flip **block algebra** proved: \(B=CD+DC\) always commutes with \(C,D\); \(B|_{V_+}=2p D_{++}\), \(B|_{V_-}=-2p D_{--}\); \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\). At \(p=5\), every SA Max-cover has \(\lambda_{\min}(B)=-6\) and \(\mathrm{op}=\sqrt{41}\), with **≥2** distinct \(D_{++}\) spectral types (both non-undercut, clique-flip, \(\max R=60\)). Random non-covers have more negative \(\lambda_{\min}(B)\). Open: prove \(\lambda_{\min}(B)=-6\Rightarrow\Phi\ge\Phi(C)\) for all Max-covers \(p\ge5\). See `evidence/E1_MAXCOVER_SPECTRUM.md`. **Does not close E(1).**

**Prop 15.37 (2026-07-27):** Continuous Γ-bound pattern: all 11 Max-covers have \(\min_{S=-p}\lambda_{\max}(\Gamma)\ge9.38758>120/13\), so cont bound holds on the **entire** spike level; discrete \(\max R=60\), clique-flip, \(\Phi=\Phi(C)\). GW insufficient for discrete (SDP\(\cdot\alpha<60\)). Residue-1 random PMs also clique-flip in samples. Open: prove cont bound / clique-flip forall Max-covers. See `evidence/E1_MAXCOVER_CONTINUOUS_BOUND.md`. **Does not close E(1).**

**Prop 15.38 (2026-07-28):** **Proved** by exhaustive \(\binom{45}{5}\): among two-sided Max-covers of size 5 on Paley \(C_{10}\), undercutters are **exactly** the 144 perfect matchings (\(\Delta=1\), \(\Phi=13\)); every two-sided cover with \(\Delta\ge2\) has \(\Phi\ge15\). Total 17154 two-sided \(k=5\) covers. Parallel cert: 80 workers. See `evidence/E1_N10_TWOSIDED_K5.md`. Supports low-\(\Delta\) undercutter pattern; **does not close E(1).**

**Prop 15.39 (2026-07-28):** On stored Max-cover PMs at \(p=5\), clique-flip pair count \(N_{\mathrm{flip}}\ge24\) (observed \(24,48,120\)). Open: prove \(N_{\mathrm{flip}}\ge1\) forall Max-cover PMs. See `evidence/E1_CLIQUE_FLIP_COUNT.md`. **Does not close E(1).**

**Prop 15.40 (2026-07-28):** **Proved:** edge-minimal undercutters satisfy \(\Phi(C\oplus F)\ge\Phi(C)-2\) (edge lip + minimality). Open step to E(1): upgrade to \(m_n\ge\Phi(C)-2\) for all Seidel \(A\) (no deeper far undercut). If that holds, gap \(O(1)\) \(\Rightarrow\) E(1) \(\Rightarrow L=\tfrac12\). **Does not close E(1) yet.**

**Prop 15.41 (2026-07-28):** **First-hit lemma proved** (first undercutting prefix on any add-edge chain has gap \(\le2\)). **Dangerous-edge criterion proved** (descent below \(\Phi-2\) requires rigid \(\sigma\)-alignment on all maximisers). **No-descent lemma OPEN** in general; if proved on the \(\rho=1\) family then \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Certified n=10:** all 144 PM undercutters have **0** dangerous edges; every add-1 returns \(\Phi\ge15\); multi-edge samples stay \(\ge13\). Evidence: `evidence/E1_NODESCENT.md`, `e1_n10_nodescent.json`. **Does not close E(1).** **F13:** do not claim \(m_n\ge\Phi-2\) from 15.40 alone.

**Prop 15.42 (2026-07-28):** **Max± dichotomy proved:** \(\Phi(A)<\Phi-2\) only possible for deep two-sided covers (\(s_+\ge2\), \(s_-\le-2\)); every \(s_+\le1\) or \(s_-\ge-1\) matrix has \(\Phi\ge\Phi-2\). **Counting freeness proved.** **Tight \(S\equiv1\) and \(S\equiv2\) no-descent proved** (frame mean \(1/p\)). **Type I freeness when \(N_1>N(p+1)/(2p)\) or \(k\le2p-2\) proved.** **Equivalence:** \(m_n\ge\Phi-2\) iff no-descent on all gap-2 undercutters (given 15.40+parity). **n=10 complete** (PM Type I strict freeness; C6 tight \(S\equiv2\)). Open: no-descent for Type I with large \(k\) and small \(N_1\), and deep non-tight \(k>2p\). **Does not close E(1).**

**Prop 15.43 (2026-07-28):** No-descent **proved** for Type I with freeness (strong) and tight deep \(S\equiv2\) (weak). Type I freeness-failure **isolated** to equality cases reducing to tight size-\(2p\) covers; at \(p=3\) 1-bit spike gives \(\Phi\ge\Phi-2\) for all-even-degree tight covers. **n=10 \(m_{10}=\Phi-2\) closed.** Residual for general \(p\): tight \(S\equiv2\Rightarrow\Phi\ge\Phi-2\) beyond \(p=3\); deep non-tight; \(k=3p-2\) boundary. **Does not close E(1) / \(L\).**

**Bi-tight / Prop 15.45 (2026-07-28):** Stars never bi-tight (all \(p>2\)). Integral bi-tight levels \(2,3,4\) **MILP-infeasible at \(p=5\)**; non-star size-\(p\) tight infeasible; deep two-sided \(k=10,12\) infeasible; size-\(2p\) Max+ covers have \(\max S_{-}\ge2\). Avg degree of level-2 bi-tight is \(4p/(p^2+1)<1\) for \(p\ge5\). Master lemma (Prop 15.44): tight Max+ covers either \(\max_{\mathrm{Max}_{-}}S\ge0\) or bi-tight. See `evidence/E1_BITIGHT.md`, `e1_star_bitight_obstruction.json`. Residual: lift to all \(p\ge5\); deep non-tight \(k>2p\).

**n=26 exact sparse MITM (2026-07-27):** shipped `phi_mitm` (exact \(\Phi\) for even \(n\le28\)). Random matchings, cycles \(C_4\)–\(C_{26}\), stars, and random \(k\le20\) flips of Paley \(C_{26}\): **0 undercuts of \(\Phi=65\)** (min observed 67 on single edges). Matching undercut of \(n=10\) does not lift. Consistent with \(k_\star=0\) at \(n=26\), not a proof. See `evidence/E1_N26_SPARSE_EXACT.md`.

**MO thread re-audit (2026-07-27):** full re-read of [MO 413935](https://mathoverflow.net/questions/413935) + comments via SE API — **0 answers; no hidden proof**. Author already knew liminf \(\ge2^{-5/2}\) and that such bounds do **not** settle existence; not interested in numerics. Starr’s Bowlin/multipartite pointer is already dead for existence (§9–§10). **Time savings:** stop BH re-derivation, stop multipartite reopen, stop pure-SA as deliverable; only load-bearing E(1) (or permanent relative gap / non-existence) counts. Details: `evidence/MO_THREAD_REAUDIT.md`.

**Existence of \(L\) remains OPEN** — sandwich only; E(1) not proved (reduced to path-cycle/\(k_\star=o(n^{3/2})\) dichotomy, not closed).
---

## 1. Exact quantity (do not restate incorrectly)

For \(n\ge2\),

\[
m_n
=
\min_{a_{ij}=\pm1}
\max_{x\in\{\pm1\}^n}
\Biggl|
\sum_{1\le i<j\le n}a_{ij}\,x_i x_j
\Biggr|,
\qquad
\alpha_n=\frac{m_n}{n^{3/2}},
\qquad
L\stackrel{?}{=}\lim_{n\to\infty}\alpha_n.
\]

**Existence of \(L\) is OPEN** — neither proved nor disproved. This handoff does **not** claim \(L\) exists or fails to exist.

### Equivalent matrix form (`solution.md` §1)

Associate the symmetric zero-diagonal Seidel matrix \(A\) with \(A_{ij}=a_{ij}\). Then

\[
\sum_{i<j}a_{ij}x_ix_j=\tfrac12 x^\top Ax,
\quad
\Phi(A)=\max_x\bigl|\tfrac12 x^\top Ax\bigr|,
\quad
m_n=\min_A\Phi(A).
\]

Same \(\alpha_n\). Tests: `test_form_Q_matches_half_xAx`, `test_equivalence_m_vs_half_phi_on_optimum_n6`.

### Factorization (`solution.md` Prop 15.6)

\[
\Phi(A)=\tfrac12\,n\,\|A\|_{\mathrm{op}}\,\rho(A),
\quad
\rho(A)=\frac{\max_{x\in\{\pm1\}^n}|x^\top Ax|}{n\,\|A\|_{\mathrm{op}}},
\quad
r(A)=\rho(A)\cdot\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}.
\]

Conference matrices uniquely minimize \(\|A\|_{\mathrm{op}}\) (floor \(\sqrt{n-1}\)). Asymptotic optimality of \(m_n\) along conference orders \(\Leftrightarrow\) \(\min_A r(A)=\rho(C)+o(1)\).

---

## 2. Proved results (with pointers)

| Item | Statement | Where |
|------|-----------|--------|
| **Sandwich** | \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) | `solution.md` Main Theorem, Prop 5.2, Thm upper via Paley/conference |
| **Dual-Gaussian LB** | For every Seidel \(A\): \(\Phi(A)\ge n\sqrt{n-1}/\pi\); hence \(m_n\ge n\sqrt{n-1}/\pi\) | `solution.md` Prop 5.2; `dual_gaussian_lower_bound` in `src/minmax_quadratic.py` |
| **BH floor (weaker)** | \(\liminf\alpha_n\ge2^{-5/2}\) | `solution.md` Prop 5.1 |
| **Cut-code identity** | \(m_n=\binom n2-2\rho(D_n)\) with \(D_n=\{\pm(x_ix_j)_{i<j}\}\) | `solution.md` Prop 1.2 / cut-code section |
| **Monotonicity / steps** | \(m_n\le m_{n+1}\le m_n+n\); consecutive \(\alpha\) gaps \(O(n^{-1/2})\) | `solution.md` §3 |
| **Denseness / Paley reduction** | Paley orders \(n_k=q_k+1\) (\(q\equiv1\bmod4\)) satisfy \(n_{k+1}/n_k\to1\); existence of \(\lim\alpha_n\) \(\Leftrightarrow\) convergence along Paley orders alone | `solution.md` Prop 6.1–6.2 |
| **Conference spectral UB** | If \(C^2=(n-1)I\), then \(\Phi(C)\le\tfrac12 n\sqrt{n-1}\) | `solution.md` §2; `spherical_half_bound` |
| **Limsup via Paley \(\rho\)** | \(\limsup\alpha_n\le\tfrac12\limsup_k\rho(C_k)\le1/2\) | `solution.md` Prop 15.8 |
| **Unique min-op / tr\(A^4\) / \(\mathbb E[Q^4]\)** | Conference uniquely minimizes op-norm, \(\mathrm{tr}(A^4)\), and \(\mathbb E[Q^4]\); \(L^2\) mass of \(Q\) universal | `solution.md` Props 15.5, 15.10–15.13 |
| **Exact optimality \(n=6\)** | \(m_6=\Phi(C)=5\) | `solution.md` Cor 15.15; `test_exact_optimality_n6_via_Q4_gap` |
| **\(\rho=1\) for \(n=p^2+1\)** | Paley over \(\mathbb F_{p^2}\): halfspace boolean evec \(Cx=px\), \(\rho(C)=1\), \(\Phi(C)=\tfrac12 n\sqrt{n-1}\) | `evidence/PROOF_rho_eq_1.md`; `solution.md` Main Theorem corollary; shipped `paley_conference_prime_power` + `halfspace_boolean_vector` |
| **Limsup \(\rho=1\) along dense family** | On \(n_k=p_k^2+1\), \(\rho=1\) for all \(k\), \(n_{k+1}/n_k\to1\) | Corollary of \(\rho=1\) theorem |

Soft multipartite / Hadamard / annealed / rank-one blow-up **cannot** force \(\liminf=\limsup\) (`solution.md` §9–§10). Dead: global Q4 path for large \(n\) (Props 15.16–15.19).

---

## 3. Open blockers for settling existence

| ID | Blocker | Notes |
|----|---------|--------|
| **E(1)** | Asymptotic optimality: \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) along Paley (or \(\min r=\rho(C)+o(1)\)) | **Exact** optimality fails at \(n=10\): \(m_{10}=13<\Phi_{\mathrm{Paley}}=15\). Must be asymptotic. Local edge-opt of Paley at small orders is not a proof. |
| **E(2)** | \(\rho(C_k)\to\rho_*\) for all Paley (or full conference) orders | On \(n=p^2+1\), \(\rho\equiv1\) already. For general Paley, \(\rho\) is increasing on small exact orders; interval constructions give constructive \(\rho_{\mathrm{int}}\gtrsim0.99\) at large \(n\) (heuristic/constructive, not full \(\rho\to1\)). |
| **Non-existence path** | Two subsequences with unequal \(\lim\alpha\) | Denseness forces any oscillation to appear along Paley orders too. No certified construction. |
| **Thm E (conditional)** | E(1)+E(2) along general Paley \(\Rightarrow L=\rho_*/2\) | Open. **Shortcut:** E(1) on the dense \(\rho=1\) family alone \(\Rightarrow L=\tfrac12\) (Prop 6.1–6.2: \(\alpha_{n_k}\to\tfrac12\Rightarrow\lim\alpha_n=\tfrac12\)). |
| **Thm F (Stolz)** | If \(\delta_n/\sqrt n\to\ell\) then \(\alpha_n\to\tfrac23\ell\) | Extension-cost regularity open. |

**Do not claim existence from sandwich + denseness alone.**

---

## 4. Numerics inventory

### 4.1 Exact \(m_n\) (load-bearing)

| \(n\) | \(m_n\) | \(\alpha_n\) | Certification |
|------:|--------:|-------------:|:--------------|
| 2 | 1 | 0.354 | live `exact_m` |
| 3 | 3 | 0.577 | live `exact_m` |
| 4 | 4 | 0.500 | live `exact_m` |
| 5 | 4 | 0.358 | live `exact_m` |
| 6 | 5 | 0.340 | live `exact_m`; = Paley \(\Phi\) |
| 7 | 9 | 0.486 | live `exact_m` |
| 8 | 10 | 0.442 | live `exact_m` |
| 9 | 12 | 0.444 | recorded multi-worker Gray (constants in tests + `evidence/exact_m_table.json`) |
| 10 | **13** | 0.411 | recorded multi-worker Gray; **\(m_{10}<\Phi_{\mathrm{Paley}}=15\)** |
| 11 | \(\le17\) (claimed \(=17\)) | \(\le0.466\) | UB witness in prior work; external cut-code package claims exact 17 (Esmaeili–Zaghian counterexample). Treat exact equality as **external claim** unless re-verified. |

**Evidence paths:**

- Live \(n=2..8\): `evidence/exact_m_table.json` (also under `{SCRATCH}/evidence/exact_m_table.json`), regenerated by `regen_exact_table.py` / pytest session fixture.
- Recorded \(m_9,m_{10}\): `tests/test_minmax.py` (`M9_RECORDED`, `M10_RECORDED`); `evidence/exact_m_table.json`.
- Dual-Gaussian vs exact: pytest `test_dual_gaussian_lower_bound_holds_for_exact_m`.

### 4.2 \(\rho=1\) family (certified)

| \(p\) | \(n=p^2+1\) | \(\Phi=\frac12 np\) | \(\rho\) | Check |
|------:|------------:|--------------------:|---------:|:------|
| 3 | 10 | 15 | 1 | halfspace evec |
| 5 | 26 | 65 | 1 | halfspace evec |
| 7 | 50 | 175 | 1 | halfspace evec |

**Evidence:** `evidence/rho1_verify.json`, `evidence/PROOF_rho_eq_1.md`, pytest `test_rho_eq_1_paley_prime_power`, `{SCRATCH}/evidence/rho1_verify.log`.

### 4.3 Heuristic / non-certified (do not treat as exact \(m_n\))

| Claim | Status | Trap |
|-------|--------|------|
| SA “best_UB” from `phi_local` | **Lower** bound on \(\Phi(A)\) only | Minimizing `phi_local` does **not** upper-bound \(m_n\) unless finalists are re-scored with **exact** \(\Phi\) |
| Interval \(\rho_{\mathrm{int}}\ge0.99\) at large \(n\) | Constructive lower bound on \(\rho\) for specific sign patterns / intervals | Supports limsup \(\rho\) near 1; not full E(2); not exact \(\Phi\) |
| Local-search \(\rho\) “dips” | Uncertified | Not non-existence |
| Prior SA at \(n=26\): no certified beater of Paley \(\Phi=65\) | Suggestive local min | Hamming distance of true optimizers may be large (at \(n=10\), 1–2 edge flips from Paley do **not** reach \(m_{10}=13\)) |
| **n=10 structure (certified)** | Hamming-5 threshold; 144 perfect-matching optima; 1-edge local opt | `evidence/N10_STRUCTURE.md`, `n10_structure.json`, `n10_matching_optima.json`; tests `test_n10_*` |
| **n=10 classification N10-C** | Maximizer-drop \(\Leftrightarrow\) \(\Phi=13\); single \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) orbit of size 144 | `evidence/N10_MATCHING_CLASSIFY.md`, `n10_matching_classify.json`; `src/n10_matching_classify.py`; `test_n10_matching_classify_*` |
| n=26 matching probe | 86 random perfect-matching flips of Paley: \(\Phi\ge73>65\) | Matchings do **not** undercut at \(n=26\) (`n26_matching_probe.json`) |
| n=26 SA + exact rescore | 86×5k SA finalists, exact \(\Phi\): none \(<65\); certified \(m_{26}\le65\) | `n26_sa_exact_rescore.json` — `phi_local` alone is **not** a UB |
### 4.4 Independent external artifacts (context only)

- Sol / Codex writeup: sandwich \(1/\pi\)…\(1/2\) + cut-code ([Robby955/mo-413935-ai-attempt](https://github.com/Robby955/mo-413935-ai-attempt)) — aligns with our Prop 5.2.
- Curtis cut-code \(M_{11}=17\): [antipodal-cut-code-k11](https://github.com/CurtisAccelerate/antipodal-cut-code-k11) — finite; does not settle the limit (author Paata confirmed).
- X thread: mostly AI hallucinations (Parisi/SK, graphon uniqueness, Wick holomorphy). See prior triage notes if present; ignore as proofs.

---

## 5. Resume playbook

### 5.1 Next attacks (ranked)

1. **E(1) on \(\rho=1\) family** — Prove \(m_n=\frac12 n\sqrt{n-1}-o(n^{3/2})\) for \(n=p^2+1\), or exhibit a permanent relative gap. At \(n=10\) the absolute gap is only \(2\) (rel. \(\approx0.063\)); at \(n=26\) SA+exact-rescore found **no** undercutter of \(\Phi=65\) (certified \(m_{26}\le65\)). Need a general argument, not local edge-flip (optima sit at Hamming \(\ge5\)). **MO re-audit:** author already discarded one-sided bounds as settlement; pure numerics are out of scope for the prize.  
2. **Structural gap from \(n=10\)** — **Mostly closed (N10-S + N10-C).** Hamming-5 threshold; only undercutters at \(k=5\) are 144 perfect matchings; those 144 form one \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit and equal the maximizer-drop set (`evidence/N10_STRUCTURE.md`, `evidence/N10_MATCHING_CLASSIFY.md`). SA also finds Hamming-11–16 optima in the switching metric (same \(r=13/15\)). Remaining: whether a matching-type construction lifts (random matchings at \(n=26\) **raise** \(\Phi\) to \(\ge73\)); classify non-matching distant optima.  
3. **E(2) analytic** — Prove \(\rho(C_n)\to1\) for all large Paley. **Partial:** exact \(\,x^\top Cx=2-8\Sigma_q\); asymptotic \(\rho_{\mathrm{int}}=(8/\pi^2)L_{\mathrm{odd}}+o(1)\) with \(\limsup\rho_{\mathrm{int}}=1\) (`E2_RHO_INT_ASYMPTOTICS.md`). Full pointwise \(\rho\to1\) still open; \(\rho=1\) on \(n=p^2+1\) already gives limsup.  
4. **Non-existence** — Only if two dense subsequences with **proved** unequal \(\alpha\) limits appear; denseness (Prop 6.2) is mandatory.
### 5.2 Traps to avoid

| Trap | Why |
|------|-----|
| Soft multipartite / Hadamard “\(c_k\to0\)” existence | Cross-block error is \(\Theta(N^{3/2})\); abstract sequences can oscillate inside the sandwich |
| Q4 / fourth-moment shell for large \(n\) | Shell vacuous when \(\Delta_*/3>\max\delta\) (`solution.md` Props 15.16–15.19) |
| Treating `phi_local` min as \(m_n\) UB | Local max underestimates \(\Phi\) |
| Uncertified \(\rho\) dips as non-existence | Need certified \(\Phi\), not SA LB |
| Exact conference optimality for all \(n\) | **False** at \(n=10\) |
| Claiming \(L=1/2\) from \(\rho=1\) alone | Needs **E(1)** on that family (Prop 6.1–6.2 then force full limit); dual-Gauss alone only gives liminf \(\ge1/\pi\) |
| ProcessPool via stdin / single-core pegging | Use script files; `W=ncpu-2`, `OMP_NUM_THREADS=1` per worker |
| Parisi / graphon / Wick “existence proofs” from X | Hallucinations |

### 5.3 Compute policy

- Workers: `full_workers = nproc - 2` when machine idle.  
- BLAS: `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` per process.  
- Scripts as files (never `python -c` + ProcessPool from stdin).  
- Scratch for ephemeral runs: goal `{SCRATCH}` only; ship durable evidence under `evidence/` and tests.

### 5.4 Key code entry points

| Function | File | Role |
|----------|------|------|
| `exact_m(n)` | `src/minmax_quadratic.py` | Exact \(m_n\) for \(n\le9\) (n=9 heavy) |
| `paley_conference_matrix(q)` | same | Prime \(q\equiv1\bmod4\) |
| `paley_conference_prime_power(p)` | same | Order \(p^2+1\) |
| `halfspace_boolean_vector(p)` | same | \(\rho=1\) evec |
| `dual_gaussian_lower_bound(n)` | same | Prop 5.2 floor |
| `phi` / `phi_local` | same | Exact / local \(\Phi\) |

### 5.5 Files map

| Path | Role |
|------|------|
| `HANDOFF.md` | **This document** — resume entry point |
| `solution.md` | Full writeup + obstruction analysis |
| `README.md` | Short status |
| `src/minmax_quadratic.py` | Shipped library |
| `tests/test_minmax.py` | Load-bearing tests |
| `evidence/` | Durable numerics + \(\rho=1\) proof note |
| `evidence/PROOF_rho_eq_1.md` | Full \(\rho=1\) proof |
| `evidence/N10_STRUCTURE.md` | Theorem N10-S: matching undercutters of Paley-\(10\) |
| `evidence/N10_MATCHING_CLASSIFY.md` | Theorem N10-C: maximizer criterion + \(\mathrm{P}\Gamma\mathrm{L}\) orbit |
| `src/n10_structure.py` | Maximizer balance + \(k\)-flip + SA structure campaign |
| `src/n10_matching_optima.py` | Perfect-matching census (144/945) |
| `src/n10_matching_classify.py` | Classification: criterion + \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) orbit |
| `src/n26_matching_probe.py` / `n26_sa_exact_rescore.py` | \(\rho=1\) family probes at \(n=26\) |
---

## 6. Conditional landscape (if E(1)+E(2) close)

- Along \(\rho=1\) family: \(\Phi(C_n)/n^{3/2}\to1/2\).  
- If E(1): \(\alpha_{n_k}\to\tfrac12\) along that dense family \(\Rightarrow\lim\alpha_n=\tfrac12\) by Prop 6.1–6.2 (both liminf and limsup).  
- If E(2) gives \(\rho_*<1\) and E(1) on general Paley, conditional Thm E says \(L=\rho_*/2\).

None of this is established. **Leave Main Theorem as sandwich + OPEN.**

---

## 7. Verification commands

```bash
cd /home/nick/quadratic-minmax-limit
# FULL SUITE — never single-core (F17). Default: pytest.ini -n 86 / scripts/pytest_full.sh
./scripts/pytest_full.sh
# targeted single test OK without -n; full suite must use W=nproc-2
python3 -c "from src.minmax_quadratic import dual_gaussian_lower_bound; print(dual_gaussian_lower_bound(10))"
# expect ~ 9.5493 = 10*3/pi
```

Expected: all tests pass; existence not claimed settled in `HANDOFF.md` or `solution.md` Main Theorem.
