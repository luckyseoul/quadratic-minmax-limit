# P0 Engineering Graph (persistent — update every turn)

**Machine:** 88 CPU + Tesla V100-SXM2-16GB (CuPy).  
**Policy:** F17 (no single-core thrash) · F19 (no moduli class-refine loops) · F20 (no GPU theater).  
**L = lim α_n:** **OPEN** (never soft-close).  
**Last update:** 2026-07-30 (Prop 15.95; session handoff `SESSION_HANDOFF_2026-07-30_prop1595.md` — load skills §0; mult≥d OPEN; L OPEN).

---

## 0. How to use this file

1. Before any heavy job: read **§1 dependency**, **§3 open nodes**, **§4 banned**.
2. Every material change: update node status + add an edge note in **§6 log**.
3. Prefer **proof edges** over census. Census only certifies a node already claimed open-by-census.
4. GPU when the unit is dense batch (m4, dots, Gram, large matmul). ProcessPool for independent irregular shards. **Never** claim GPU if wall is class-key CPU.

---

## 1. Dependency DAG (what closes L)

```
L = lim α_n                                    [OPEN]
 └─ Prop 6.1–6.2 denseness + sandwich          [PROVED]
     └─ E(1) on ρ=1 family n=p²+1              [OPEN]  ← ultimate blocker
         │
         ├─ Path C (bi-tight / Type I empty for p≥5) + deep non-tight ND
         │    ├─ g_min > T(p)  ⇐  max|m4| ≤ L_abs  or λ2(P⊙P)≤4/N
         │    │    ├─ max|m4| ≤ M_cand ≤ M_mid ≤ L_abs     [OPEN ∀p≥5]
         │    │    │    ├─ GD: star·S1 ≤ 0                  [OPEN ∀p≥5; cert p=5,7]
         │    │    │    ├─ joint/S3 + cand criterion        [OPEN ∀p≥5; cert p=5,7]
         │    │    │    └─ moduli pin c* safe-side c_GD     [p=5 DONE; p=7 nul=2 OPEN]
         │    │    └─ λ2(P⊙P) ≤ 4/N  (16N / H / residual)  [OPEN ∀p≥5; cert p=5,7]
         │    └─ deep non-tight ND / gap-2                  [OPEN]
         │
         ├─ Path A: k⋆ = o(n^{3/2}) structure               [OPEN]
         └─ Path B: matching non-undercut + dichotomy      [OPEN; many F-bans]
```

**Shortest honest settlement:** prove E(1) on n=p²+1, **or** close Path C residual for all p≥5 **and** deep ND.

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
| **P0-1** | **N_GD_ALL** or **N_MCAND_ALL** | closes bi-tight residual for Path C | re-census p=5,7; refine more C-keys for fun |
| **P0-2** | **N_SPEC_GAP** | independent Path C residual | re-run gap_probe only on {3,5,7} without new idea |
| **P0-3** | **N_MOD_P7_MULTI** | only if it yields a **general** method | more strat names (coarse+cr+κ+…) |
| **P0-4** | **N_DEEP** | after bi-tight for general p | — |
| **P0-5** | **N_E1** | ultimate; SA only as evidence | claim E(1) from n=26 SA |

**Current active node:** P0-1 — character-sum / boolean structure for GD (or cand), **not** class refine.

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

---

## 7. Completion gate (Path C residual)

Do **not** claim Path C residual closed unless:

- [ ] N_MCAND_ALL **or** N_SPEC_GAP closed for all primes p≥5 (proof, not p∈{5,7} census)
- [ ] N_DEEP addressed or explicitly deferred with reason
- [ ] G_L still marked OPEN unless E(1)+denseness full chain is written and verified
- [ ] Graph §6 log updated; tests green under `./scripts/pytest_full.sh`

---

## 8. Next action (single, from graph)

**Just landed:** Prop 15.95 (Wick_hi≤thr_gap ∀p≥5; strengthened mult+Wick⇒gap; C_diag; mult=d cert p=3,5,7).  

**Active next (closes bi-tight residual) — pick one:**  
1. **Gap path (narrowed):** prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) for all primes \(p\ge5\), and/or \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}=12n^2+48n\). Then Prop 15.95 ⇒ gap ⇒ bi-tight empty. (∑M²≤thr is free once ∑M²≤Wick, for p≥5.)  
2. **16N path:** prove \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\).  
3. **H path:** prove \(\lambda_{\max}(FF^\top|_{1^\perp})\le N(3+H)\).  

Structure: Wick≤thr algebra closed; P⊙P kills range(P); mult=d at p=3,5,7.  
Do **not** re-attack ∑ρ κ_B (dead after 15.90).

Then **N_DEEP**. **G_L stays OPEN until E(1) or full Path C+deep.**

**Not active:** moduli class refine (F19); GPU theater (F20); residual-as-separate-foothold (dead after 15.90).
