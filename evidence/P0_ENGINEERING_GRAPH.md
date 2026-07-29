# P0 Engineering Graph (persistent — update every turn)

**Machine:** 88 CPU + Tesla V100-SXM2-16GB (CuPy).  
**Policy:** F17 (no single-core thrash) · F19 (no moduli class-refine loops) · F20 (no GPU theater).  
**L = lim α_n:** **OPEN** (never soft-close).  
**Last update:** 2026-07-29 (after user rebuke: no graph, GPU theater, moduli thrash).

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

---

## 7. Completion gate (Path C residual)

Do **not** claim Path C residual closed unless:

- [ ] N_MCAND_ALL **or** N_SPEC_GAP closed for all primes p≥5 (proof, not p∈{5,7} census)
- [ ] N_DEEP addressed or explicitly deferred with reason
- [ ] G_L still marked OPEN unless E(1)+denseness full chain is written and verified
- [ ] Graph §6 log updated; tests green under `./scripts/pytest_full.sh`

---

## 8. Next action (single, from graph)

**Active:** P0-1 **N_GD_ALL** — algebraic attack: express \(\mathrm{star}\cdot S_1\) / \(\mathbb E[ZU_1]-\mathbb E_{\mathrm{Wick}}[ZU_1]\) via Paley character sums on boolean \(+p\)-evecs; aim for sign ≤0 independent of census.

**Not active:** any new `e1_gmin_m4_refine*` strat; any ProcessPool-only job labeled “GPU”.

If a verification needs hardware: one CuPy job that is **≥5s GPU-bound** (e.g. multi-prime type6 is not Max+; for Max+ structure use dense moments on p=5,7 with large batch + report `nvidia-smi` util in evidence).
