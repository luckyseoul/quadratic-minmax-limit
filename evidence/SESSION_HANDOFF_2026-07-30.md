# Session handoff — 2026-07-30 (after unexpected shutdown / 1% quota)

**Repo:** `/home/nick/quadratic-minmax-limit`  
**HEAD:** `265fccc` — Prop 15.52 (Max+ coord sum + m4 moduli sketch)  
**Branch:** `main` ahead of `origin/main` by **78** commits (do **not** push unless asked)  
**Working tree:** clean at handoff write  
**L = lim α_n:** still **OPEN** — do **not** soft-close Main Theorem  

---

## 0. One-line resume (read this first)

**P0 in flight:** prove \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) (or any LB \(>T(p)\)) for **all primes \(p\ge5\)** via refined-\(C\) m4 moduli + Max+-free pin (\(\mathrm{Tr}(G^2)\) / \(G\)-spectrum). Then P1 deep non-tight; **only if both close**, Main Theorem \(L=\tfrac12\).

**Last interrupted job (not finished):** generalize nullity-1 + \(\mathrm{Tr}(G^2)\) pin beyond \(p=5\); durable proof + tests. Mid-session finding committed below: **PSD-max does not pin true \(c\)** (F16).

---

## 1. Settlement chain (load-bearing; do not re-derive)

\[
m_n\ge\Phi-2\text{ on dense }\rho=1\text{ family }n=p^2+1
\;\Rightarrow\;\mathrm{E}(1)\;\Rightarrow\;L=\tfrac12
\]
via denseness (Prop 6.2). Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) intact.

**F13:** never claim \(m_n\ge\Phi-2\) from Prop 15.40 alone.  
**F3:** never soft-close on sandwich + denseness + ρ=1.

**Path C residual:**
1. Bi-tight empty if \(g_{\min}>T(p)=-(p-2)/(p(2p-1))\) (Prop 15.47) — certified \(p=5,7\).
2. Uniform LB candidate \(L(p)=-(p-2)/(2p^2)\) with \(L>T\) for odd \(p>2\); holds at certified \(g_{\min}(5),g_{\min}(7)\); correctly fails at \(p=3\).
3. Deep non-tight gap-2 ND / always \(\Phi\ge\Phi(C)\) for \(p\ge5\).

---

## 2. What is proved / certified (do not re-derive)

| Item | Where |
|------|--------|
| Props 15.45–15.52 writeup | `solution.md` |
| Bi-tight / 1-bit spike | Props 15.46–15.47; `E1_BITIGHT.md` |
| CR classification of \(g_{\min}\) | Prop 15.49; `e1_gmin_cr_classify.json` |
| Cond means = Gaussian interpolant; Fréchet too weak | Prop 15.50; **F15**; `e1_gmin_cond_mean.json` |
| \(a\leftrightarrow T\) equiv form | Prop 15.51; `e1_gmin_structure.json` |
| \(\mathbf1^\top y=(p+1)y_\infty\); m4 moduli sketch | Prop 15.52; `E1_GMIN_MODULI.md` |
| Certified \(g_{\min}\) | \(p=3\colon-1/3\); \(p=5\colon-3/65\); \(p=7\colon-109/2863\) |
| \(a_{\min}\) | \(p=5\colon5/39>1/9\); \(p=7\colon75/818>1/13\) |
| Uniform LB algebra + p=5,7 check | `e1_gmin_uniform_lb.json` |
| Nullity 1 of refined evec system | cert \(p=5,7\); 14 / 24 classes |
| \(\mathrm{Tr}(G^2)\) pin recovers exact \(g_{\min}\) at \(p=5\) | `E1_GMIN_MODULI.md`, commit `fd4c0ee` |
| G spectrum rank \(\binom{d}{2}-d+1\) | cert p=3,5,7 |
| **80 pytest green** (2026-07-29) | `tests/test_minmax.py` + `tests/test_gmin_residual.py` (13 gmin tests) |
| Failure graph F1–**F16** | `evidence/E1_FAILURE_GRAPH.md` |

Scratch verification (ephemeral): `/tmp/grok-goal-c8bc6d526c68/implementer/pytest_full.txt` — 80 passed in ~94s.

---

## 3. Best attack state (moduli path)

### 3.1 Setup (works)

- Stratify 4-sets by pure \(C\)-invariants \((\mathrm{CR},\kappa,\triangle\text{-type})\).
- Averaged evec identities → \((pI-M)\mathbf m=\mathbf b\) combinatorial.
- Nullity 1: \(\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n\).
- Sum-of-disj-\(G\) already in system (\(s_{\mathrm{null}}=0\)).
- Pin: \(\mathrm{Tr}(G^2)=E+2n_{\mathrm{wedge}}/p^2+6\sum n_A m_A^2\) → quadratic in \(c\).

### 3.2 Critical late finding (was mid-shutdown; now in `E1_GMIN_MODULI.md`)

At \(p=5\) on the line \(G(c)\):
- \(\mathrm{rank}=66\) + \(\lambda_{\min}\approx0\) hold on a **continuum** of \(c\).
- **Max \(g_{\min}\) under PSD** ≈ \(-0.040\) at \(c\approx-0.42\) — **not** the true point (\(c\approx-0.291\), \(g_{\min}=-3/65\approx-0.046\)).
- True pin: \(\mathrm{Tr}(G^2)\) (or full spectrum). **F16** bans PSD-max as selection rule.
- \(\mathrm{Tr}(G^2)\) needs **\(2\times\)** off-diagonals (bug already fixed once).

### 3.3 Open for all \(p\ge5\)

1. Prove nullity 1 for all primes \(p\ge5\).
2. Max+-free closed form for \(\mathrm{Tr}(G^2)\) or \(G\)-eigenvalues as function of \(p\).
3. Select correct quadratic root without Max+ samples.
4. Prove \(g_{\min}\ge L(p)\) (or any LB \(>T\)).

**p=7 note:** not distance-homogeneous (≥2 Max+ types); need full transitions for stable quadratic, not sampled.

---

## 4. Priority order tomorrow

### P0 — \(g_{\min}\) LB (BLOCKER)

**Goal:** theorem in `solution.md` + evidence + test: \(g_{\min}\ge L(p)\) (or \(>T\)) ∀ primes \(p\ge5\).

**Promising (in order):**
1. Finish Max+-free \(\mathrm{Tr}(G^2)\) / spectrum formula (association scheme / frame identities).
2. Character-sum closed form for \(m_4\) on minimizing \(\|\kappa\|=1\) CR class over \(\mathbb F_{p^2}\).
3. Optional numeric: exact \(g_{\min}\) at \(p=11\) (not a proof).

**Do not reopen (dead + F-graph):**
4-point LP, Chebyshev, Wick-as-LB, \(-3/\Phi\) general, bare \(C\)-types, affine/PGL halfspace orbit as full Max+ (60/260 at p=5), min-norm \(V_+\), plain Fréchet (**F15**), PSD-max pin (**F16**), incomplete Tr(G²) (missing 2× off-diag).

**Optional durable code (not yet a file):** implement `src/e1_gmin_moduli.py` (build refined classes, particular+null, Tr(G²) quadratic, cert p=5 vs \(-3/65\)); keep ProcessPool **from script files only**.

### P1 — Deep non-tight (after or independent of gmin)

For \(s_{\min}=2\), \(\max S\ge4\), \(k>2p\): ND or always \(\Phi\ge\Phi(C)\) for \(p\ge5\).  
Scripts/evidence: `src/e1_deep_*.py`, `E1_RESIDUAL.md`.

### P2 — Only if P0+P1 close

1. Main Theorem \(L=\tfrac12\) in `solution.md` with full chain.  
2. `HANDOFF.md` → SETTLED.  
3. Evidence + tests for every new claim; full pytest green.  
4. **No** F3 / F13 soft-close.

---

## 5. Caches / compute

```text
Machine: ~88 cores. full_workers ≈ nproc-2.
BLAS: OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1 per process.
ProcessPool ONLY from script files (stdin fails).
Caches (rebuild after reboot if missing):
  /tmp/maxplus_p5.npy          # present as of last session
  /tmp/e1_p7/maxplus.npy       # + Binv, C, coords, nb
```

```bash
cd /home/nick/quadratic-minmax-limit
python3 src/e1_gmin_cr_classify.py
python3 src/e1_gmin_uniform_lb.py
python3 src/e1_gmin_structure.py
OMP_NUM_THREADS=1 python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py -q
# expect: 80 passed (or more if new tests added)
```

---

## 6. Key file map

| Path | Role |
|------|------|
| **This file** | Primary resume for tomorrow |
| `HANDOFF.md` | Long research handoff (Props through 15.52) |
| `solution.md` | Settled writeup; Main Theorem still OPEN sandwich |
| `evidence/E1_GMIN_MODULI.md` | **Active attack surface** — moduli + pin + PSD continuum |
| `evidence/E1_FAILURE_GRAPH.md` | F1–F16 banned loops |
| `evidence/E1_RESIDUAL.md` | Residual narrative |
| `evidence/e1_gmin_*.json` | Certs (cr_classify, structure, uniform_lb, cond_mean, …) |
| `src/e1_gmin_structure.py` | Prop 15.51 structure |
| `src/e1_gmin_cond_mean.py` | Prop 15.50 |
| `src/e1_gmin_cr_classify.py` | CR classifier |
| `src/e1_gmin_uniform_lb.py` | LB candidate checks |
| `src/e1_gmin_linear_system.py` / `formula.py` / `char_sum.py` | Prior formula hunts |
| `tests/test_gmin_residual.py` | Residual unit tests |
| `minmax_quadratic.py` / `src/minmax_quadratic.py` | Core library |

---

## 7. Explicit non-goals

- Soft-close Main Theorem without full chain  
- Claim \(m_n\ge\Phi-2\) from edge-minimal gap (F13)  
- Sequential MILP / serial SA (F1, F2)  
- Exhaust all PMs (F8)  
- Re-derive sandwich / denseness as deliverable (F9)  
- Treat cont Γ as discrete LB (F4)  
- Push 78 local commits without user request  
- Reopen F15 / F16  

---

## 8. Suggested skills

- `use-available-compute` — before Max+ rebuild or p=11 cert  
- `verification-before-completion` — before any “settled” claim  
- `goal-verifier` — if user asks whether L is done  
- `handoff` — only if packing again under quota pressure  

---

## 9. Acceptance criteria

| Milestone | Criteria |
|-----------|----------|
| g_min closed | Proof in `solution.md`: \(g_{\min}\ge L(p)\) (or better) ∀ primes \(p\ge5\), evidence+test |
| Deep closed | Proof ND / non-undercut for deep non-tight on family |
| **Project settled** | Main Theorem \(L=\tfrac12\); HANDOFF SETTLED; full pytest green; no F3/F13 |

Until then: status line remains **L OPEN**.

---

## 10. What this session did / did not do (honesty)

**Did:** Consolidated post-shutdown state; recorded PSD continuum anti-thrash (**F16**); refreshed moduli evidence; wrote this durable handoff for tomorrow.

**Did not (quota):** Finish the general-\(p\) proof of \(g_{\min}\ge L(p)\); deep non-tight; Main Theorem. Those remain the live research residual.
