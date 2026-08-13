# P0 Engineering Graph (persistent — update every turn)

**Machine:** 88 CPU + Tesla V100-SXM2-16GB (CuPy).  
**Policy:** F17 (no single-core thrash) · F19 (no moduli class-refine loops) · F20 (no GPU theater) · F21 (no False-predicate 15.xxx) · F22 (no e1 flip on census).  
**L = lim α_n:** **OPEN** (2026-08-13). The 2026-08-05 “CLOSED” claim was retracted (disj Gsum LB / 15.158). Residual **(ii) CLOSED**; residual **(i) OPEN**. Trust `HANDOFF.md` / `STATUS.md`, not the DAG below (stale Path-C framing).  
**Optional open:** Path-C residual / 16N (not required for denseness L).  
**Last update:** 2026-08-13 evening (residual-(i) attack log; leftover still \(\lvert\mu\rvert\le1/(2p)\)).

---

## 0. How to use this file

1. Before any heavy job: read **§1 dependency**, **§3 open nodes**, **§4 banned**.
2. Every material change: update node status + add an edge note in **§6 log**.
3. Prefer **proof edges** over census. Census only certifies a node already claimed open-by-census.
4. GPU when the unit is dense batch (m4, dots, Gram, large matmul). ProcessPool for independent irregular shards. **Never** claim GPU if wall is class-key CPU.

---

## 1. Dependency DAG (what closes L)

```
L = lim α_n = 1/2                              [CLOSED 2026-08-05]
 └─ denseness + sandwich                       [PROVED]
     └─ E(1) on ρ=1  OR  Path C + deep ND
         │
         ├─ Path C
         │    ├─ N_ED4_SUF: δ²≤room_hyp/24 (⇔ orth≤room_hyp)  [OPEN ∀p≥5; cert 3,5,7]
         │    │     ρ_min² form sufficient for p≥7 (15.110/15.117)
         │    ├─ mult≥d−1 (15.98) + orth≤room_hyp ⇒ 16N  [form CLOSED]
         │    ├─ bi-tight empty                              [follows 16N+15.61]
         │    └─ N_DEEP non-tight ND                         [OPEN]
         ├─ Path A/B                                         [OPEN; secondary]
```

**Shortest honest settlement:** close **N_ED4_SUF** for all p≥5, then N_DEEP, then Main Theorem; **or** E(1).

---

## 2. Node table (typed)

| ID | Type | Statement | Status | Compute | Evidence / code |
|----|------|-----------|--------|---------|-----------------|
| **G_SAND** | Goal | sandwich 1/π ≤ liminf ≤ limsup ≤ 1/2 | **closed** | algebra | solution.md |
| **G_DENSE** | Goal | dense ρ=1 Paley n=p²+1 | **closed** | algebra | PROOF_rho_eq_1 |
| **G_L** | Goal | lim α_n exists (=1/2 ideal) | **open** | — | HANDOFF |
| **N_E1** | Theorem | m_n = Φ − o(n^{3/2}) on ρ=1 family | **open** | multi-W SA/exact; not pure GPU | e1_n26_* |
| **N_BITIGHT** | Theorem | bi-tight empty if g_min > T(p) | **closed** (form) | algebra | Prop 15.47 |
| **N_MCAND_ALG** | Algebra | M_cand ≤ M_mid ≤ L_abs < T_abs | **closed** | algebra | Prop 15.70/74 |
| **N_MCAND_P5** | Census | max\|m4\| = M_cand = 3/65 at p=5 | **closed** | GPU m4 + mmap | e1_gmin_m4_gpu/kernel |
| **N_MCAND_P7** | Census | max\|m4\| ≤ M_cand at p=7 | **closed** | GPU m4 + mmap | same |
| **N_MCAND_ALL** | Theorem | max\|m4\| ≤ M_cand ∀ primes p≥5 | **open** | proof (not more p=5/7 census) | residual 15.74–15.82 |
| **N_GD_P5** | Census | star·S1 ≤ 0 all \|κ\|=1 centres p=5 | **closed** | GPU+W | S1_* evidence |
| **N_GD_P7** | Census | star·S1 ≤ 0 all \|κ\|=1 centres p=7 | **closed** | GPU+W | same |
| **N_GD_ALL** | Theorem | GD / E[ZU1]≤E_Wick ∀ p≥5 | **open** | **proof** (char sums / Aut) | Prop 15.80 residual |
| **N_MOD_P5** | Pin | type6+CR or coarse: nul=1, c* safe-side c_GD | **closed** | CPU evec (after GPU m4) | 15.81–15.82 |
| **N_MOD_P7_CONST** | Struct | type6+CR m4 constant (48 cls) | **closed** | CPU ProcessPool keys + GPU m4 | refine.json |
| **N_MOD_P7_NUL1** | Pin | nullity ≤1 on constant strat | **blocked** | — | nul=2 even at 130 cls; e4 dependent |
| **N_MOD_P7_MULTI** | Pin | multi-param + TrG2 + 2nd pin ⇒ cand+GD | **open** | algebra + GPU TrG2 moments | pin_extra: e4 no help |
| **N_SPEC_GAP** | Theorem | λ2(P⊙P) ≤ 4/N ∀ p≥5 | **open** | GPU eig / structure | spectral.json; cert p=5,7 |
| **N_ED4_SUF** | Theorem | δ²≤room_hyp/24 (⇔orth≤room_hyp) ∀p≥5 | **open** (closed Max+-free p=5,7) | algebra / design energy | 15.110–15.151; m_e cov formula; need Cov(k,t_e)/t_e(k) bound |
| **N_DEEP** | Theorem | deep non-tight ND | **open** | multi-W | spike at p=5 only |
| **C_MAXP** | Asset | Max+ arrays on disk | **p=5,7 only** | generation cost | /tmp/maxplus*; no p≥11 |
| **C_GPU** | Asset | CuPy V100 available | **yes** | — | gpu_budget |
| **C_GRAPH** | Process | this file maintained | **active** | — | update §6 |

### Status legend
- **closed** — proved or fully certified; do not re-derive  
- **open** — needs new edge  
- **blocked** — known obstruction under current formulation  
- **active** — being worked  

---

## 3. Critical path (ranked — do this order)

| Rank | Node | Why | Do **not** |
|------|------|-----|------------|
| **P0-1** | **N_ED4_SUF** | primary Path C residual orth≤room_hyp (δ²≤room_hyp/24) | class_key thrash (F19); full E_4p energy; p=5/7 census-only |
| **P0-2** | **N_SPEC_GAP** / 16N | alt residual form λ2≤4/N | re-run without new idea |
| **P0-3** | **N_DEEP** | after bi-tight for general p | — |
| **P0-4** | **N_E1** | ultimate; SA only as evidence | claim E(1) from n=26 SA |
| **legacy** | N_GD_ALL / N_MCAND_ALL | older residual tree | reopen as primary |

**Current active node:** **N_ED4_SUF** — residual Max+-free closed at p=5,7. **15.151:** regular-set t-identities; **E[t_e]=n_e(p+3−2e)/(8p)** Max+-free; **m_e=n/2+8p Cov(k,t_e)/(n_e(p+3−2e))**; residual⇔Cov(k,C(k,3)) budget; exact m_e at p=5; CS/Popoviciu too weak. **Next:** closed t_e(k) or character-sum Cov(k,t_e). Dead: CS on Cov, plain ULC, free-orbit type lists p≥11, class_key thrash, soft-close.

---

## 4. Banned edges (do not traverse)

| ID | Ban | Why |
|----|-----|-----|
| F3 | Soft-close L | explicit |
| F14 | Ignore failure graph | session thrash |
| F16 | Pin c by max PSD g_min | wrong root |
| F17 | Single-core multi-minute | user rebukes |
| **F19** | **Moduli class-invariant thrash** | type6+CR/coarse+CR already constant at p=5,7; more keys did not drop p=7 nullity; e4 dependent. **No new proof edge.** |
| **F20** | **GPU theater** | write `use_gpu=True` after 0.1–0.4s CuPy m4, then burn wall on CPU class-build/evec; user never sees GPU util. Either (a) GPU is the bulk of wall or (b) say “CPU multi-W, GPU unused this step”. |
| F15/F18 | Fréchet LB; incomplete Max+ orbits | dead |

---

## 5. Compute routing (honest)

| Work unit | Right backend | Notes |
|-----------|---------------|-------|
| All-quad m4 given Max+ Y | **1× GPU** CuPy batch | one H2D; parent only |
| Pairwise dots / Tr(G²) / Gram | **GPU** if N,n fit 16GB | chunk if needed |
| Class keys (type6, CR, hist) | **ProcessPool W=86** | irregular; not GPU |
| Averaged evec / SVD nullity | **CPU** small | classes ≪ 200 |
| Boolean Max+ generation | **multi-W CPU** (search) | no p≥11 cache yet |
| Character-sum / Aut proof | **algebra** | GPU only to check closed forms |
| Full pytest | `./scripts/pytest_full.sh` W=86 | never bare pytest |

**Audit (this session):**  
- Real GPU wall in refine: ~0.13s (p=5) + ~0.38s (p=7).  
- Subsequent refine_moduli / multi / pin_extra: **GPU util 0%**, pure CPU ProcessPool.  
- That is **F20**. Graph now forbids it as a “GPU run”.

---

## 6. Change log (append only)

| Date | Edge | Note |
|------|------|------|
| 2026-08-01 | **N_ED4_SUF pin + Prop 15.113** | Graph critical path = ED4≤ED4_suf. Proved ⟨f_y,Tκ⟩=2p(p⁴−1), ED4=3n²+4E[W²], Q_δ≤ρ_min²⇒residual. Cert p=3,5,7. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU algebra, no class_key. `prop15113.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.114** | Proved Tf_y=(4p−2γ)⊙f_y; δ⊥γ⊙f; ∑γ=(6/p)C(n,4); ∑γ²=6C(n,4)+n(n−1)(n−2)/4 (adj edges cancel); ‖Tf_y‖² closed. Cert p=3,5. p=7: **3 ED4(y) types** (not 2-pt homogeneous), all ≤ED4_suf; pair-avg ED4≤ED4_suf. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, no class_key. `prop15114.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.115** | Proved E[γ⊙f]=2κ/p (2-design)⇒Aρ=b; κ⊥E_{4p}⇒δ=P_{E_4p}m₄=E_y P f_y; spectral m1=4p−12/p, Var=24(p²−3)(p²−4)/(p²(p²−2)). Cert resolvent+moments p=3,5; ED4≤suf p=3,5,7. Full E_4p energy of f_y too crude. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+. `prop15115.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.116** | Proved e₄ poly; ∑m₄²=ED4/24+const(n); ⟨κ,ρ_min⟩=n(n−1)(n−2)(n−6)/(2p²(p²−5)); flat Pythagoras ID; coherent mass=δ²; Aut-line criterion. Min-dist ED4 envelope **too weak** p≥5. Cert δ²≤ρ_min² p=3,5,7. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+. `prop15116.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.117** | Path C primary=δ²≤room_hyp/24 closed form; ρ_min²≷hyp by p; slack ID; ⟨b,f⟩=⟨ρ_min,b⟩=2(p⁴−1)/p; ⟨ρ_min,m₄⟩ closed; E⟨b,γ⊙f⟩=0; if pointwise ⇒⟨ρ_min,f⟩=4(p⁴−1)/(3(p²−5)). Cert bgf=0 p=3; hyp residual cert when Max+ present. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+. `prop15117.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.118** | Pointwise ⟨b,γ⊙f_y⟩=0 on Max+; ⟨κ,γf⟩=p(p²+1)(p²−1)(p²−4)/4; ⟨T²κ,m₄⟩=⟨T²κ,f_y⟩=8p²(p⁴−1); ⟨ρ_min,f_y⟩=4(p⁴−1)/(3(p²−5)). Cert full Max+ p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+. `prop15118.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.119** | Closed ED4/EW2/m4f flat+bud rationals; residual ⇔ ED4≤bud ⇔ E[W²]≤EW2_bud (⇔ m4f≤m4f_bud if Q_δ const); WE structure (parity, D_max); crude 2n D_max² too weak. Cert WE+saturation p=3,5; orth·N=147456 at p=5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15119.py` + tests 5 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.120** | Pointwise E[W²]_y=EW2_flat+6Q_δ(y); ∑m₄²=F+δ² (F=m4f_flat); maj ED4≤2n³ too weak; CS-γ/LP dead. Cert factorization p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15120.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.121** | ED4=4n²+4N⁻²‖FFT|_{1^⊥}‖_F²; E[W²]=‖FFT‖_F²/N²; EW2_flat=(n²+T²/m)/4; residual ⇔ Frob FFT budget ⇔ Φ variance ≤room_hyp; 16N=op-norm; H/16N·Tr too weak. Cert p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15121.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.122** | Max+ dots: u=(y−z)/2∈V₊, uᵀCu=pk; λ_max(T) threshold; Aut-line δ∈E_4p^{Aut}; B(p)=EW2_bud−d²; LP/PGL still weak. Cert identity p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15122.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.123** | Switch→srg(n,p(p−1)/2,μ−1,μ); Max+↔V₊∩{0,1}ⁿ; regular sets; A'₄=∑m₄²; residual⇔A'₄≤m4f_bud; two-valued form in V₊∩1^⊥. Cert p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15123.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.124** | Proved E[k],E[k²],E[k³] Max+-free; exact_≤3 partition of E[k⁴]; ED4=ed4_from_exact3+16 R₄; residual⇔R₄≤R₄_bud; Hamming Delsarte weak p≥5; W_{p+1}=d. Cert saturation p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+scipy LP, GPU unused. `prop15124.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.125** | W_k=# τ-equitable bipartitions; Max+=spherical 2-design, residual=4-design defect; ed4_from_exact3 factored; R₄_bud closed; Delsarte+moments j≤3 still weak p≥5; antipodal dual. Cert p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+scipy LP, GPU unused. `prop15125.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.126** | Subfield F_p∪{∞} Hoffman seed (cert 3,5,7); 1-design algebra; simplex b≤d (eq p=3); W_{p+1}=d only cert p=3,5. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15126.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.127** | **Closed W_{p+1}** via χ₄ formula; inversive plane S(3,p+1,p²+1); **W=d false** at p=7 (W=11≠25); census p=3,5,7,11. Full W_k still open. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU ProcessPool enum W=86, GPU unused. `prop15127.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.128** | Full W_k census p=3,5,7 (p=7 via 2²⁵ Max+ W=86, ~12s); ED4(p=7)=12835984/409<bud; W_10=0 structural zero; moments j≤3 OK. **N_ED4_SUF still OPEN** (no general closed W_k). L OPEN. F19/F20: CPU multi-W, GPU unused. `prop15128.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.129** | Jensen δ²≤E‖P f‖²; residual dictionary cert p=3,5,7; r̄=W_{p+1}(p+1)/n; **p=7 Hoffman not 1-design** (r̄=44/25); p=5 Hoffman pair geometry. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+, GPU unused. `prop15129.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.130** | P m₄=δ; Jensen refined; **gap room−ρ_min²>0 for all primes p≥7** (closed form); census δ²≤ρ_min² at p=3,5,7; Aut-line program. **N_ED4_SUF still OPEN** (need general δ²≤ρ_min²). L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15130.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.131** | Pair-avg vs basepoint ED4; p=7 **3-type Q_δ** (exact); **true δ²=19180800/1840091** (W-value was Q_hs); max Q≤ρ_min² slack 812048/220451; Var(Q)>0. Residual OK at p=7. **N_ED4_SUF still OPEN** (general max Q / δ²). L OPEN. F19/F20: CPU ProcessPool W=86, GPU unused. `prop15131.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.132** | Max+-free residual ⇔ ∑m₄²≤m4f_bud; Aut-invariance of δ / Q orbit-constancy; γ even + formal 4p-fiber; γ=0 mass pattern p=5/7; **moment-LP + pole+D_max dead** for p=5..19. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+W=86, GPU unused. `prop15132.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.133** | class_key T-spectrum dim E_4p^{ck}=**0,1,0** at p=3,5,7; **F19 quantitative** (p=7 nul=0, δ²>0); PGL CR **dead** as Aut(C); Aut-line OK p=5 only. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU ProcessPool W=86, GPU unused. `prop15133.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.134** | Strict Aut G=\|p²(p²−1)\|; orbits **9/42/128**; dim E_4p^G=**0,2,7**; residual proj matches δ² at p=5,7; G carries p=7 residual; Aut-line dim≤1 **fails** on G. **N_ED4_SUF still OPEN** (need Gauss m₄). L OPEN. F19/F20: CPU W=86, GPU unused. `prop15134.py` + tests 3 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.135** | Spectral form δ²=∑c_j² on E_4p^G; **hs F_p char** Max+-free (∑f_hs=e₄); G·hs **dead** (δ²_Ghs≫room); P1=Pκ=0 (moments don't pin c_j). **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15135.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.136** | **Max+-free flat** on G-orbits (resolvent) for general p; m₄=flat+δ; geom type **insufficient** for m₄; residual = free c_j only; partial p=5 ∞-m₄. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15136.py` + tests 3 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.137** | Q_j G-equivariant; c_j=∑w_t Q_j(y_t); **p=5** c_j=(3/13)Q_j(hs)+(10/13)Q_j(y_*); Q_j(hs) Max+-free; census r=1/2/5 at p=3/5/7. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+Max+ cert, GPU unused. `prop15137.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.138** | Max+-free y_* via hs-switch+norm circle; **p=5 residual Max+-free** (∑c_j²=room, e2e test); p=7 covers 2/5 types; p=13 ansatz empty. **N_ED4_SUF still OPEN for all p≥5.** L OPEN. F19/F20: CPU, GPU unused. `prop15138.py` + tests 5 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.139** | Affine halfspaces Max+ (p=5,7); **p=7 AP dichotomy** 84 vs 56; **double Seidel–norm-circle** covers 588/1176/294; **all five H_+ size classes Max+-free at p=7**. Residual still needs Q_j+weights Max+-free. **N_ED4_SUF OPEN.** L OPEN. F19/F20: CPU W=86, GPU unused. `prop15139.py` + tests 5 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.140** | Weights \|G\|/\|Stab\|; character-sum Q_j; **p=7 ∑c_j²=δ²_pair≤room** (match true); **7/8** Max+-free types; one size-1176 OPEN. **N_ED4_SUF OPEN.** L OPEN. F19/F20: CPU vectorized Q, GPU unused. `prop15140.py` + tests 5 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.141** | Size-12 Seidel partner closes last 1176 type; **all 8 H_+ Max+-free**; **p=7 residual Max+-free**; bi-tight path p=5,7. **N_ED4_SUF still OPEN for general p.** L OPEN. F19/F20: CPU, GPU unused. `prop15141.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.142** | Affine all-S Max+ at p=5,7,11; k-AP split multi-way at p=11; **fourths partners p=5-only** (dead as general size-12); p=7 size-12 fibre 84 T; p=11 type samples. **Residual general OPEN.** L OPEN. F19/F20: CPU W=86, GPU unused. `prop15142.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF structure + Prop 15.143** | p=11 affine **6-orbit census** complete Max+-free (132×1,330×2,660×3); dbl-switch LB \|H₊\|≥28182; nc density class-dependent. **Full H₊ types / residual p≥11 OPEN.** L OPEN. F19/F20: CPU, GPU unused. `prop15143.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.144** | Free H₊ orbit \|O\|=\|G\|=14520; size-2420 chain; **type-enum residual DEAD p≥11**; redirect type-free δ²≤ρ_min². Residual only p=5,7 closed. **L OPEN.** F19/F20: CPU, GPU unused. `prop15144.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.145** | Type-free residual package: δ²≤ρ_min² ⇔ ‖ρ‖²≤2ρ_min² ⇔ ‖m₄‖²≤m4f_suf ⇔ ED4≤ED4_suf; closed m4f_suf/ED4_suf; ρ_min²/room→5/8 (monotone p=7..97); census δ²≤ρ_min² p=5,7 only. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15145.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.146** | Type-free R₄/μ₄ channel: ED4_from_exact3=−n⁴+28n²−40n; δ²=(2/3)(R₄−R4_flat); R4_suf=R4_flat+(3/2)ρ_min²; μ₂=n/2, μ₃=0; spectral w* Max+-free but ≫ budgets (Jensen dead). Cert R₄ p=5,7. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15146.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.147** | Inclusion densities d₁=1/2,d₂=(p²+1)/(4p²),d₃=(p²+3)/(8p²); R₄=E[k^{underline 4}]; residual⇔d₄≤d4_suf; **U=d₃²/d₂ < d4_suf** all p≥5 (P(p²)>0) so ULC would finish residual; **ULC fails** at p=5,7 (~3.6%,1.9%) while residual holds. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15147.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.148** | Relaxed-ULC calculus: C_act=C_flat+κδ²; C_max=Q(p²)/((p²−5)(p²−2)(p²+3)²); residual⇔C_act≤C_max; **uniform target** d₄≤U(1+C₇/p²) for p≥7 with C₇=79923/87373; C_max→1⁻; census window nonempty. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15148.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.149** | Size-bias form: residual⇔E_μ[k]≤k_suf (μ∝k^{underline 3}dW); k_suf/k_flat closed; d₂−1/4, d₃−1/8, d4_flat>1/16; **triple-covered Gauss program** (λ_τ on srg triples). Cert size-bias p=5,7. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15149.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.150** | srg n_e closed; **λ_e=N(p+3−2e)/(8p)** (Aut+affine+moments j≤3; cert p=5,7); **π_e Max+-free**; residual⇔∑π_e m_e≤k_suf. Finite-type m_e attack ready. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15150.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.151** | t-identities; **E[t_e]** Max+-free; **m_e=n/2+8p Cov(k,t_e)/…**; residual as Cov(k,C(k,3)) budget; exact m_e p=5; **CS dead**. Need t_e(k) or char-sum Cov. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15151.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.152** | Free-param **t₃**; p=5 closed t₃(α); **p=7 multi-type** mid-weights (corrects 15.151.B); residual ∑πm ≡ **R₄** channel; Cov(t_e)→Cov(t₃). Pure t_e(k) **dead**. Need R₄≤R4_suf / char-sum m_e. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction + W=86 census, GPU unused. `prop15152.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.153** | Switched **μ₄** dictionary: e1=e3=0,e2=1/(n−1); d4=(1+6/(n−1)+μ4)/16; **μ4_flat/μ4_suf closed**; residual⇔μ4≤μ4_suf; m_e p=7 exact. Char-sum surface ready. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15153.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.154** | **avg(χκ)=3/(n−3)** proved Max+-free (C₂ algebra); μ₄=κ_main+η; **κ_main=3/(p²(p²−2))** under budget; residual⇔**η≤η_suf**. Need Weil on Ext correlation. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15154.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.155** | Aut-line closed: **e₄(s)** poly; **Tχ=χ(4p−2σ_z)**; **Q(s)** on C₂-evecs; **η=c₁R₄+c₀** explicit. Residual = E[s⁴] of one Paley coordinate. Crude 2n³ dead. Need Weil/3-design. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15155.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.156** | **κ₄=E[s⁴]−12n²** residual; closed **κ4_flat/κ4_suf**; bridge κ₄=(n)_4η−16n; 4-design is LB only; crude+moment-LP dead. Need Weil/3-design-UB. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+LP, GPU unused. `prop15156.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.157** | **Gegenbauer defect** μ_G4=E[Q₄(s/n)]; E[s⁴]=n⁴(a₀+a₄μ); residual⇔μ_G4≤μ_G4_suf (closed in p). Census p=5,7. UB 1/h₄ false; d/h₄ too weak. Need Weil/Aut on μ_G4. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction, GPU unused. `prop15157.py` + tests 4 passed. |
| 2026-08-01 | **N_ED4_SUF + Prop 15.158** | **Q₄ closed**; Max+ **1-hom** Gram 2N×d; **IP-scheme blocked** (p=5); pole bound N≥N_* holds p=5,7 but +Q₄ support blocks; Chebyshev dead. Need Weil/coherent config/W_k. **N_ED4_SUF still OPEN.** L OPEN. F19/F20: CPU Fraction+numpy, GPU unused. `prop15158.py` + tests 4 passed. |

| When | Edge | Change |
|------|------|--------|
| 2026-07-29 | init | Graph created after rebuke. Closed: N_MOD_P5, N_MOD_P7_CONST, N_MCAND_P5/P7, N_GD_P5/P7. Open: N_MCAND_ALL, N_GD_ALL, N_SPEC_GAP, N_MOD_P7_MULTI, N_E1. Banned F19/F20. |
| 2026-07-29 | F19 | pin_extra: e4 + dense evec do not cut p=7 nullity (stays 2). multi-strat coarse+CR also nul=2. **Stop refine loop.** |
| 2026-07-29 | F20 | Confirmed: V100 idle during multi/pin_extra; only refine m4 used CuPy briefly. |
| 2026-07-29 | F20+ | Bad recovery attempt: `gd_gpu_moments` added **busywork 4096 GEMMs** to inflate util; CUDA unknown error at p=7. Script **deleted**. Real GPU = necessary dense kernels only, with util logged — not synthetic load. |
| 2026-07-29 | process | Graph + F19/F20 filed in E1_FAILURE_GRAPH + HANDOFF. Active node remains **N_GD_ALL** (algebra), not census/refine. |
| 2026-07-29 | **C_GPU** | **DOWN** after F20+ crash: `nvidia-smi` → “No devices were found” / handle Unknown Error on 0000:03:00.0. Needs host GPU reset (sudo nvidia reload or reboot). Do **not** schedule GPU jobs until C_GPU = yes. |
| 2026-07-30 | **N_MCAND_ALL / N_GD_ALL structure** | **Prop 15.83** (Max+-free algebra): proved \(\mathrm{gain}_L-\mathrm{gain}_{\mathrm{cand}}=3(p-2)/(48(2p+3))>0\) and cascade \(M_{\mathrm{cand}}<M_{\mathrm{mid}}\le L<T\) for all primes \(p\ge5\). Ranks residual targets; does **not** prove gain bound or close N_MCAND_ALL. CPU Fraction only (F20 GPU unused). Evidence: `e1_gmin_m4_prop1583.json`, `src/e1_gmin_m4_prop1583.py`, solution Prop 15.83. **G_L still OPEN.** |
| 2026-07-30 | **C_GPU** | Host reports V100 present again (`nvidia-smi` ok in compute-budget). Prefer real dense GPU only when wall is batch-dominated. |
| 2026-07-30 | **N_GD_ALL structure** | **Prop 15.84**: \(B_{\mathrm{cand}}=(p^3-4p^2-7p-6)/(p^2(2p+3))\); sign \(B_5<0\), \(B_{p\ge7}>0\); GD+\(S_3\le B_{\mathrm{cand}}\)⇒cand; proved \(4p-d_1<0\) (abs bootstrap dead). CPU Fraction (F20 unused). Evidence: `e1_gmin_m4_prop1584.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP structure** | **Prop 15.85**: \(Q_4\) mean/fluctuation split; \(S_1=0\); unit ray \(=-\mu-(\mu/2)S_w+\tfrac12 Be^\top\widehat G Be\); \(\mu\) harmless vs \(H\); Ĝ load-bearing. Links spectral H to signed disj \(m_4\). CPU Fraction. Evidence: `e1_gmin_m4_prop1585.json`. **G_L OPEN.** |
| 2026-07-30 | **N_GD_ALL structure (ε)** | **Prop 15.86**: \(\sum\mathrm{star}\cdot\tau_1=\varepsilon(p)\,n_1\) with \(\varepsilon=(-1)^{(p-1)/2}\) (closes 15.80.4); mean Wick \(=\varepsilon/p^2\); GD⇒mean \(\mathrm{star}\cdot S_1\le0\); \(\tau_1\) AP value set size \((p-1)/2\); \(B_{\mathrm{cand}}\to1/2\), \(B_{\mathrm{cand}}/d_3=\Theta(1/p^2)\). Certified W=86 pure C at p=3,5,7,11. Evidence: `e1_gmin_m4_prop1586.json`. **Pointwise GD still OPEN. G_L OPEN.** |
| 2026-07-30 | **N_GD_ALL structure (K4/g)** | **Prop 15.87**: K4 star theorem \(\sum\mathrm{star}=0\) on \(|\kappa|=1\) (48/64 labelings); \(S_1(a)=g\cdot\mathrm{star}_a\); GD\(\Leftrightarrow g\le0\); CS too weak (\(\Theta(p)\) vs Wick \(O(1)\)); \(E[U_1^2]\approx d_1\) (cert p=3,5,7). Evidence: `e1_gmin_m4_prop1587.json`. **Pointwise \(g\le0\) still OPEN. G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP structure (H-gap)** | **Prop 15.88**: pairwise \(\sum_{i<j}y_iy_j=p\) on Max+; \(n/2-(3+H)=(p^4-8p^2-16p-21)/(2(p^2+1))>0\) for \(p\ge5\) (H\(\Rightarrow\) bi-tight empty); \(g=p\rho\mathrm{star}-2/p^2-\mathrm{star}S_3\). Cert p=5: \(g=g(\tau_1)\). Evidence: `e1_gmin_m4_prop1588.json`. **H / pointwise g still OPEN. G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (Wick split)** | **Prop 15.89**: \(Q_4/N=2+4/p^2+8\sum\rho\kappa_B\) (Wick from \(\sum\kappa_C\kappa_B=(n+1)/4\|B\|^2\), cert p=3,5,7,11); H \(\Leftrightarrow\sum\rho\kappa_B\le(H-1-2/p^2)/4\). **H residual = bound \(\sum\rho\kappa_B\). G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (residual≡H)** | **Prop 15.90**: residual bound algebraically ≡ H; pointwise κ_B identity cert p=3,5; orth form ≡ H; holds p=3,5,7. **No new foothold** — prove ray≤H independently (orth energy or 4th-moment op). Evidence: `e1_gmin_m4_prop1590.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (independent dual H)** | **Prop 15.91**: dim Z=d(d−3)/2 proved; orth/Φ/κ/harm ≡ H; sphere<Wick<6+2H≤16; 2×sphere⇒16N for p≥5. **H still OPEN.** Preferred targets: orth LB, λ_max(κ|Z)≤(p+1)(p+7)/d, harm budget, or 2×sphere. Evidence: `e1_gmin_m4_prop1591.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (pairing + λ₂ reductions)** | **Prop 15.92**: ∑m4 κ_C=n(n−1)(n−2)/8 proved (constant on Max+); 16N⇔λ₂(P⊙P)≤4/N; H⇔λ₂(P⊙P)≤(3+H)/(2N); W saturates H at p=3,5. **Bound still OPEN.** Evidence: `e1_gmin_m4_prop1592.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (Gu/FFT)** | **Prop 15.93**: FFT1=Nd1 proved; 16N⇔λ_max(FFT\|1⊥)≤8N; H⇔≤N(3+H); Gu eigs cert p=3,5 (non-Nd =(N/2)spec Φ\|Z). **Bound OPEN.** Evidence: `e1_gmin_m4_prop1593.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (annihilator + gap mult)** | **Prop 15.94**: P⊙P kills range(P) (central sym); gap criterion mult(λ₂)≥d + ∑M²≤4d²(d+4). Cert p=5 mult=d criterion holds for gap. **OPEN mult≥d ∀p≥5.** Evidence: `e1_gmin_m4_prop1594.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (Wick≤thr + mult structure)** | **Prop 15.95**: Wick_hi≤thr_gap for all primes p≥5 (Fraction); mult≥d+∑M²≤Wick⇒gap∀p≥5; C_diag=4n(11n−14)/p². Cert mult=d and gap_by_mult p=5,7. **OPEN mult≥d and ∑M²≤Wick ∀p≥5.** Evidence: `e1_gmin_m4_prop1595.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (Wick–κ calculus)** | **Prop 15.96**: yᵀΣy=2n; ⟨Wick,M⟩=12n²; ∑M²=12n²−48n+‖κ‖²; ∑M²≤Wick⇔‖κ‖²≤96n; C_diag≤96n. Cert p=3,5. Backend: CPU algebra, GPU unused (F20). **OPEN ‖κ‖²≤96n and mult≥d ∀p≥5.** Evidence: `e1_gmin_m4_prop1596.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (Veronese/Ky Fan mult)** | **Prop 15.97**: mult(λ₂)=mult(λ_max(Γ)); Ky Fan: mult≥d ⇔ d orthonormal maximizers of Var(cᵀBc); Aut-Schur. Cert mult=d p=3,5. Backend: CPU linalg, GPU unused (F20). **OPEN: construct maximizers ∀p≥5.** Evidence: `e1_gmin_m4_prop1597.json`. **G_L OPEN.** |
| 2026-07-30 | **N_SPEC_GAP (mult≥d−1 PSL)** | **Prop 15.98**: mult(λ₂)≥d−1 for Paley Max+ (PSL(2,p²) min irrep); mult≥d−1+‖κ‖²≤96n⇒gap ∀p≥5 (Fraction). Cert p=3,5. Backend: CPU algebra, GPU unused (F20). **OPEN: only ‖κ‖²≤96n.** Evidence: `e1_gmin_m4_prop1598.json`. **G_L OPEN.** |
| 2026-07-31 | **N_SPEC_GAP (Fickus orth reduction)** | **Prop 15.101**: PopP Fickus rank; orth=16N²(∑λ²−S²/m); PSL level count; N-free ε-criterion. Cert bulk p=3,5. Backend: CPU Fraction+Max+, GPU unused (F20). **OPEN: general orth.** Evidence: `e1_gmin_m4_prop15101.json`. **G_L OPEN.** |
| 2026-07-31 | **N_SPEC_GAP (resolvent δ)** | **Prop 15.102**: ρ_min closed form; orth=24‖δ‖²; target⇔δ²≤room_hyp/24. Cert p=3,5. **OPEN: δ ∀p≥5.** Evidence: `e1_gmin_m4_prop15102.json`. **G_L OPEN.** |
| 2026-07-31 | **N_SPEC_GAP (δ cert 3,5,7)** | **Prop 15.103**: δ²≤room_hyp/24 at p=3,5,7 (eq 3,5; ratio≈0.187 at 7); 16N at p=5,7; bi-tight at p=5,7. Backend: p=7 Max+ W=86 + PopP power, GPU unused (F20). **OPEN: general p≥5.** Evidence: `e1_gmin_m4_prop15103.json`. **G_L OPEN.** |

---


- 2026-07-31: Prop 15.99 κ structure (min-distance, budget algebra, closed forms, master source); cert p=3,5; κ bound OPEN for general p; L OPEN; backend CPU Fraction, GPU unused (F20).


- 2026-07-31: Prop 15.100 dual-frame/flat≤Wick/κ_hyp; orth OPEN; L OPEN; CPU Fraction (F20).
- 2026-07-31: Prop 15.101 Fickus Gram residual / PopP bulk variance; N-free ε-criterion; cert p=3,5 orth=room_hyp; general orth OPEN; L OPEN; CPU Fraction+Max+ (F20 GPU unused).
- 2026-07-31: Prop 15.102 resolvent δ-calculus: ρ_min closed form, κ_min=proj, orth=24‖δ‖², target⇔δ²≤room_hyp/24; cert p=3,5; general δ OPEN; L OPEN; CPU Fraction+Max+ (F20).
- 2026-07-31: Prop 15.103 δ-bound cert p=3,5,7 (eq 3,5; ratio≈0.187 at 7); 16N at p=5,7; bi-tight at p=5,7; general p≥5 OPEN; L OPEN; p=7 Max+ ProcessPool + PopP power (F20 GPU unused).

## 7. Completion gate (Path C residual)

Do **not** claim Path C residual closed unless:

- [ ] N_MCAND_ALL **or** N_SPEC_GAP closed for all primes p≥5 (proof, not p∈{5,7} census)
- [ ] N_DEEP addressed or explicitly deferred with reason
- [ ] G_L still marked OPEN unless E(1)+denseness full chain is written and verified
- [ ] Graph §6 log updated; tests green under `./scripts/pytest_full.sh`

---

## 8. Next action (single, from graph)

**Just landed:** Prop 15.103 (δ-bound at p=3,5,7; 16N at p=5,7; bi-tight closed for those primes; general p≥5 OPEN).  

**Active next (closes bi-tight for all Paley p≥5):**  
1. **Primary:** prove $\|\delta\|_2^2 \le \mathrm{room}_{\mathrm{hyp}}/24$ for all primes $p\ge5$, or $\lambda_2(P\odot P)\le4/N$ (16N).  
2. Then N_DEEP + Main Theorem. **G_L stays OPEN.**  

Do **not** replace full Max+ by halfspace sub-orbit (exceeds Wick at p=5).  
Then **N_DEEP**. **G_L stays OPEN until E(1) or full Path C+deep.**

**Not active:** moduli class refine (F19); GPU theater (F20).

### 2026-07-31 — Prop 15.105 Norton/Fickus variance identity
- **Proved:** ∑(λ_α(Φ|Z)−avg)² = ‖κ_orth‖²; orth=0⇒Φ scalar; mult≥d + κ²≤96n ⇒ λ_max(Φ)≤16 (exact Fraction).
- **Cert:** p=3,5,7 variance+16N; mult(top)=d at 3,5,7.
- **OPEN residual:** ‖κ‖²≤96n / δ²≤room_hyp/24 / direct λ_max≤16 for all p≥5.
- **L OPEN** (no soft-close). Backend: CPU Fraction + Max+ cert; GPU unused (F20).
- Evidence: `src/e1_gmin_m4_prop15105.py`, `evidence/e1_gmin_m4_prop15105.json`, `tests/test_prop15105.py` (9 passed).

### 2026-07-31 — Prop 15.106 rest-average-8 + kurtosis form
- **Proved:** mult≥d ⇒ L≤16 ⇔ rest-mean≥8; at L=16 rest-mean=8 exactly; κ≤96n ⇔ κ₄(D)≤3+12/n.
- **Cert:** p=3,5,7 kurtosis+16N+mult=d.
- **OPEN:** κ₄≤3+12/n or λ_max≤16 for all p≥5. **L OPEN** (no soft-close).
- Evidence: `src/e1_gmin_m4_prop15106.py`, `evidence/e1_gmin_m4_prop15106.json`, tests 7 passed. CPU Fraction, GPU unused (F20).

### 2026-07-31 — Prop 15.107 Theorem A (d−1 + room_hyp ⇒ 16N)
- **Proved:** mult≥d−1 + orth≤room_hyp ⇒ λ_max(Φ)≤16 for all primes p≥3 (Fraction majorization).
- **Proved:** Gegenbauer α=−6/(d+4), β=3/((d+2)(d+4)); 4-design ED4≤wick_hi.
- **Cert:** p=3,5,7 orth≤room_hyp + 16N.
- **OPEN residual:** orth≤room_hyp (⇔δ²≤room_hyp/24) for general p≥5. **L OPEN**.
- Evidence: `src/e1_gmin_m4_prop15107.py`, `evidence/e1_gmin_m4_prop15107.json`, tests 7 passed. CPU Fraction, GPU unused (F20).

### 2026-07-31 — Prop 15.108 residual-Gram/Schur dual; algebraic Thm A*; Parseval T_ρ
- **Proved:** Thm A gap closed form 128(p−3)(p+3)(p⁴−12p²−5)/[(p²−5)²(p²+1)²]≥0 all primes p≥3; λ_max(Φ)=4N λ₂(PopP); orth=16N²‖R‖_F²; orth≤room_hyp ⇔ ∑ρ²≤T_ρ(p)=ρ_min²+room_hyp/24 (Max+-free); m4 Parseval expansion.
- **Cert:** p=3,5,7 ∑ρ²≤T_ρ (eq 3,5; ratio≈0.639 at p=7); 16N via PopP.
- **OPEN residual:** ∑ρ²≤T_ρ for general p≥5 (⇔ orth≤room_hyp). **L OPEN** (no soft-close).
- Evidence: `src/e1_gmin_m4_prop15108.py`, `evidence/e1_gmin_m4_prop15108.json`, `tests/test_prop15108.py` (9 passed). CPU Fraction + Max+ cert; GPU unused (F20).

### 2026-07-31 — Prop 15.109 Φ–m4 identity; Aut δ; PF+rank
- **Proved:** E[(yᵀBy)²]=6‖B‖²+8⟨m4,κ_B⟩; ∑κ_B² formula; δ∈E_{4p}^{Aut}; λ₂(PopP)<d/N for p≥5; residual ⇔ c²≤room_hyp/24 when dim Aut-ker≤1.
- **Cert p=5:** nullity Aut E_{4p}=1, c²=room_hyp/24 equality.
- **OPEN:** c²≤room_hyp/24 (⇔∑ρ²≤T_ρ) for all p≥5. **L OPEN**.
- Evidence: `src/e1_gmin_m4_prop15109.py`, `evidence/e1_gmin_m4_prop15109.json`, tests 6 passed. CPU, GPU unused (F20).

### 2026-07-31 — Prop 15.110 closed Max+ identities; ρ_min<budget p≥7
- **Proved:** ∑κ∏ and e₄ closed on Max+; ρ_min²<room_hyp/24 for all primes p≥7; residual for p≥7 ⇐ δ²≤ρ_min²; c=Q_0(halfspace).
- **Cert:** δ²≤ρ_min² at p=3,5,7.
- **OPEN:** δ²≤ρ_min² (or c²≤budget) for general p≥5. **L OPEN**.
- Evidence: `src/e1_gmin_m4_prop15110.py`, tests 7 passed. CPU Fraction, GPU unused (F20).

### 2026-07-31 — Prop 15.111 pair Schur; closed α_κ,α_ρ; Φ residual = 8⟨δ,κ_B⟩
- **Proved:** zero-diag ∑κ_C κ_B identity; α_κ=(p²+2)/(4p²) on Z; pair=(p²+11)/(4(p²−5)); α_ρ=(7p²+5)/(2p²(p²−5)); channel recon β_b=6/p, β_T=6(3p²+5)/p²; E[(yᵀBy)²]=μ̄‖B‖²+8⟨δ,κ_B⟩ so Φ excess is pure δ-channel; 16N⇔max⟨δ,κ_B⟩≤(n−10)/(n−6).
- **Cert:** Schur scalarity of κ/p², b, Tb, ρ_min on Z at p=3,5,7 with closed forms.
- **OPEN:** δ²≤ρ_min² (or max⟨δ,κ_B⟩≤(n−10)/(n−6)) for general p≥5. **L OPEN** (no soft-close).
- Evidence: `src/e1_gmin_m4_prop15111.py`, `evidence/e1_gmin_m4_prop15111.json`, `tests/test_prop15111.py` (8 passed). CPU Fraction + Z cert; GPU unused (F20).

### 2026-07-31 — Prop 15.112 design moments; conf ‖κ‖²; ED4 residual dictionary
- **Proved:** conf ‖κ‖₂² closed; Max+ antipodal; ED4=ED4_flat+24δ² ⇒ δ²≤ρ_min²⇔ED4≤ED4_suf; ED4_suf<ED4_bud for p≥7.
- **Cert:** E[yyᵀ]=2P₊, E[D²]=2n, δ²≤ρ_min² / ED4≤ED4_suf at p=3,5,7.
- **F19 note:** class_key not m₄-equitable at p=7 — do not thrash moduli as closing argument.
- **OPEN:** ED4≤ED4_suf for general p≥5. **L OPEN** (no soft-close).
- Evidence: `src/e1_gmin_m4_prop15112.py`, `evidence/e1_gmin_m4_prop15112.json`, `tests/test_prop15112.py` (6 passed). CPU Fraction + Max+ cert; GPU unused (F20).
