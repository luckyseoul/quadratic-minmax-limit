# Session handoff — residual (i) attack (2026-08-06)

**Goal:** Close denseness residual (i) → E(1)/L (Max+-free; no soft-close).  
**Outcome:** **NOT CLOSED.** Residual (ii) already CLOSED; residual (i)/E1/L remain OPEN.  
**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**Scratch (ephemeral):** `/tmp/grok-goal-03c774661dc8/implementer/` (p=7 Max± cache, pytest logs, e1_l_block.md)

---

## Predicates (must stay honest until hinge proved)

| Flag | Value |
|------|--------|
| `gsum_disj_lb_proved_general()` | **False** (`src/e1_gmin_m4_prop15170.py`) |
| residual (ii) | **True** (15.179 freeze-to-tight) |
| type_I / residual (i) | **False** |
| `e1_closed_general()` | **False** |
| L = lim α_n = 1/2 | **OPEN** |

Do **not** flip predicates from census (p=5,7). Soft-close forbidden.

---

## What closes residual (i)

Any **one** Max+-free proof for all primes \(p\ge5\):

1. **Moment:** \(|\mu_4|\le 2/n\) or \(\le 1/(2p)\) on every \(|\kappa|=1\) 4-set  
   ⇒ Gsum ≥ −1/p (or −4/n) ⇒ dual-eq Farkas (15.176)  
2. **Ker-box:**  
   \(\max\{\kappa_e : \mathrm{Gsum}\,\kappa=0,\ 1^\top\kappa=0,\ -\alpha\le\kappa_g\le 1-\alpha\ (g\neq e)\} < 2-\alpha\)  
   Sufficient if \(\max\kappa_e \le \tfrac32\cdot 3(n-2)/(pn)\) (15.192: \(\tfrac32\cdot\)scheme < need for \(p\ge5\)).

Then: flip `gsum_disj_lb_proved_general` → type_I → e1 → L with real tests; STATUS/package CLOSED; AI-test ≥2 cold.

---

## Shipped this campaign (15.174–192)

### Residual structure / Farkas (earlier)
- **15.176:** Correct μ_* = (6/p−4)/k; −1/p suffices; −2/p does **not**.
- **15.179:** Residual **(ii) CLOSED** (freeze \(S_H\equiv 3\) ⇒ k=3p−1).
- **15.182:** Dual-eq normal form \(x=\alpha\mathbf{1}-2e_*+\kappa\), \(\kappa_e=2-\alpha\).

### μ₄ / T² / φ path
- **15.184–187:** \(T^2\kappa=-24\varphi+48\kappa\); Paley \(|\varphi|\le 2(p-2)\); \(\mu_{\mathrm{part}}\) majorant ≤1/(2p).
- **15.188:** \(|\mu|\not\le|\mu_{\mathrm{part}}|\) at p=7; viable target **|\mu|≤2/n**.
- **15.189:** π=C/p Max+-free; G+ PSD ⇒ |\mu|≤1−2/p (**too weak**).
- **15.191:** Derang_perm=1, star_sum=0 on |\kappa|=1; Cy size1+size2=−2φ; **|\mu|≤|f₄| fails at p=7**.

### Ker-box path
- **15.190:** scheme-ker ⊆ ker(Gsum); scheme max \(\kappa_e=3(n-2)/(pn)<2-\alpha\) all p≥3.
- **15.192:** Gsum diag=2, row sum=n, avg disj=2/(n−3); **Aut_e averaging** (dual-eq ⇔ Aut_e-inv dual-eq); (3/2)·scheme < need for p≥5.

### Census (not a proof)
| p | max \|\mu\| on \|\kappa\|=1 | 2/n | 1/(2p) | dual-eq ker-box |
|---|---------------------------|-----|--------|-----------------|
| 3 | 1/3 | 1/5 | 1/6 | **FEASIBLE** (max κ_e=14/5) |
| 5 | 3/65 | 1/13 | 1/10 | **empty** (369/455) |
| 7 | 109/2863 | 1/25 | 1/14 | **empty** (11736/19775) |

True C-preserving Aut_e (end of session): |Aut_e|=8,24,48 for p=3,5,7; max κ_e same fractions; C constant on edge orbits.

### Extra ker fact (session notes)
- \(K_{ij}=b(C_{ij}-2/n)\) lies in ker(Gsum) Max+-free (yᵀ(K⊙C)y=0 on Max±).  
- Under box, pure C-dir max is only \(\alpha(n-2)/(n+2)\) (small); does **not** alone upper-bound full max κ_e at p=3.

### Dead / weak for close
- G+ PSD only → 1−2/p  
- |\mu|≤|μ_part| or |\mu|≤|f₄| (false at p=7)  
- Boolean 4-cube LP / 10×10 Gram → still 1−2/p  
- Worst-case e-row with only |\mu|≤1−2/p → does **not** block dual-eq  
- Soft-close / census-only predicate flip  

### Live next attacks
1. Max+-free \(\max\kappa_e\le\tfrac32\cdot\)scheme-max (or any <2−α) via Aut_e orbit algebra + closed-form Gsum between orbits.  
2. Max+-free |\mu|≤2/n: master + Aut-invariant E_{±4p} control, or Cy-expansion Q₄ exact formula.  
3. E-row-only relaxation max κ_e (exact Gsum): ~1.19 (p=5), ~1.17 (p=7) both < need — find Max+-free upper bound on this LP still <2−α.

---

## Files to load first

| Path | Role |
|------|------|
| `STATUS.md` | One-page status |
| `src/e1_gmin_m4_prop15170.py` | gsum / e1 predicates |
| `src/e1_gmin_m4_prop15176.py` | Farkas threshold |
| `src/e1_gmin_m4_prop15179.py` | residual ii CLOSED |
| `src/e1_gmin_m4_prop15182.py` | dual-eq normal form |
| `src/e1_gmin_m4_prop15186–192.py` | μ / ker latest |
| `tests/test_prop15186.py` … `test_prop15192.py` | load-bearing |
| `evidence/share/denseness_path_package.md` | short package (OPEN) |

**PYTHONPATH=src.** Optional p=7 cache: `/tmp/grok-goal-…/implementer/e1_p7/{maxplus,maxminus}.npy` (session scratch; may vanish).

**Tests (example):**
```bash
cd /home/nick/quadratic-minmax-limit
PYTHONPATH=src python3 -m pytest tests/test_prop15179.py tests/test_prop15186.py \
  tests/test_prop15189.py tests/test_prop15190.py tests/test_prop15191.py \
  tests/test_prop15192.py -q
```

---

## Do not

- Soft-close residual (i), E1, or L  
- Flip `gsum_disj_lb_proved_general` without Max+-free hinge  
- Re-open Path-C / 16N thrash unless asked  
- Assert lim α_n=1/2 in package until predicates true  

**Current claim:** residual (ii) closed; residual (i)/E1/L **OPEN**.
