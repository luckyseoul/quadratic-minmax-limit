# Session handoff — 2026-07-29 (resume after shutdown / 1% quota cut)

**Repo:** `/home/nick/quadratic-minmax-limit`  
**HEAD:** `6add963` (Prop 15.50 conditional means; branch `main` ahead of origin ~72)  
**L = lim α_n:** still **OPEN** — do **not** soft-close Main Theorem  

**New (Prop 15.50, proved):** Max+ conditional means given two coordinates equal the Gaussian
frame interpolant; \(\mathbb E[y_ky_l\mid y_i,y_j]=\alpha+\delta y_iy_j\); disj avg \(G=1/(p^2-2)\).
Fréchet on cond cov **too weak** for \(L(p)\) (**F15**). Evidence: `e1_gmin_cond_mean.json`.

---

## 0. One-line resume

**Prove** \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) for all primes \(p\ge5\) (character sum / CR class), **then** deep non-tight ND/spike; **only if both close**, update Main Theorem \(L=\tfrac12\).

---

## 1. Settlement chain (load-bearing; do not re-derive)

\[
m_n\ge\Phi-2\text{ on dense }\rho=1\text{ family }n=p^2+1
\;\Rightarrow\;\mathrm{E}(1)\;\Rightarrow\;L=\tfrac12
\]
via denseness (Prop 6.2). Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) intact. **F13:** never claim \(m_n\ge\Phi-2\) from Prop 15.40 alone.

**Path C residual (Props 15.45–15.49):**
1. **Bi-tight empty** for all \(p\ge5\) if \(g_{\min}>T(p)=-(p-2)/(p(2p-1))\) (Prop 15.47). Certified \(p=5,7\).
2. **Uniform LB candidate** \(L(p)=-(p-2)/(2p^2)\): algebra \(L(p)>T(p)\) for all odd \(p>2\); holds at certified \(g_{\min}(5),g_{\min}(7)\); correctly **fails** at \(p=3\).
3. **Deep non-tight** gap-2 ND / always \(\Phi\ge\Phi(C)\) (p=5 deep covers spike above Φ).

---

## 2. What is already proved / certified (do not re-derive)

| Item | Where |
|------|--------|
| Props 15.45–15.49 writeup | `solution.md` (~Prop 15.49) |
| Bi-tight Gsum floor | Prop 15.47; `E1_BITIGHT.md` |
| 1-bit spike | Prop 15.46 |
| Edge algebra + CR structure | Prop 15.48; `src/e1_gmin_*.py` |
| Full CR classification of \(g_{\min}\) | Prop 15.49; `e1_gmin_cr_classify.json` |
| \(g_{\min}=-\alpha_\star\) on const-\(m_4\), \(\|\kappa\|=1\) | same |
| Certified values | \(p=3\colon-1/3\); \(p=5\colon-3/65\); \(p=7\colon-109/2863\) |
| Uniform LB algebra + p=5,7 check | `e1_gmin_uniform_lb.json`, `src/e1_gmin_uniform_lb.py` |
| G spectrum rank \(\binom{d}{2}-d+1\) | cert p=3,5,7 |
| Matching form margin p=5 | \(\mathbf1_M^\top G\mathbf1_M\ge9.96>4\) |
| 74 pytest green (last known) | `tests/test_minmax.py`, `tests/test_gmin_residual.py` |
| Failure graph F1–F14 | `evidence/E1_FAILURE_GRAPH.md` |

---

## 3. Residual attack — priority order

### P0 — Prove \(g_{\min}\ge L(p)\) for all primes \(p\ge5\)

**Facts to use:**
- \(G_{ee'}=E_+[f_e f_{e'}]\) on oriented edges; \(g_{\min}=\min_{\text{disj}}G\).
- \(g_{\min}=-\alpha_\star\) where \(\alpha_\star\) is max of \(m_4\)-type correlation on **minimizing** CR classes with \(\kappa=\pm1\).
- At gmin: pattern like \(\{-3,3,3\}/65\) at p=5 (CR ∈ {5,6}).

**Promising routes (in order):**
1. **Character-sum closed form** for \(m_4\) (or \(\alpha\)) on the min \(\|\kappa\|=1\) PGL(2)-CR class over \(\mathbb F_{p^2}\). Goal: exact fraction or bound \(\ge L(p)\).
2. **Association-scheme / Bose–Mesner** eigenmatrix of the Max+ edge design — prior attempt: G reconstructs; algebra span was **large**, not a small closed scheme that pins gmin alone. Reopen only with a new invariant that collapses rank.
3. **Joint residual CLT + explicit error** that forces \(\alpha_\star\le(p-2)/(2p^2)\). Prior CLT residual **overestimated** m4 — need rigorous error, not heuristic.
4. **SOS / moment hierarchy** on 4-point marginal **with** boolean-evec extension constraints (plain 4-point LP was too weak ~−0.6).

**Do not reopen (dead):**
- 4-point extendable LP alone  
- Chebyshev / Wick-as-LB  
- \(-3/\Phi\) as general LB (fails p=7)  
- Bare \(C\)-types (G not constant)  
- Affine halfspace orbit as proxy for full Max+ (undersamples; wrong gmin)  
- Pure degree pigeon  
- Min-norm \(V_+\) interpolation alone (underestimates)  
- Star tautology / incomplete Max+ samples for gmin  
- Ising maxent  

**Numeric next certs (if proof stalls):** compute exact \(g_{\min}\) at p=11 (needs Max+ enum or smart orbit; caches wipe on reboot). Confirm \(g_{\min}\ge L(11)\). File under `evidence/e1_gmin_p11.json` if obtained. **Not a proof.**

### P1 — Deep non-tight residual

After (or in parallel if independent of gmin):
- For \(s_{\min}=2\), \(\max S\ge4\), \(k>2p\): prove ND or always \(\Phi\ge\Phi(C)\) for \(p\ge5\).
- Lever: Prop 15.46 1-bit spike; p=5 large-k deep covers already spike \(\Phi>\Phi(C)\).
- Scripts: `src/e1_deep_*.py`, evidence under `e1_deep*`, `E1_RESIDUAL.md`.

### P2 — Only if P0+P1 close

1. Update `solution.md` Main Theorem: \(L=\tfrac12\) with full chain citation.  
2. Update `HANDOFF.md` status to SETTLED.  
3. Evidence + tests for every new claim.  
4. Full pytest green; record in scratch.  
5. **Do not** soft-close on sandwich+denseness (F3) or m_n from 15.40 (F13).

---

## 4. Caches / compute

```text
Machine: ~88 cores. full_workers ≈ nproc-2.
BLAS: OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1 per process.
ProcessPool ONLY from script files (stdin fails).
Caches (rebuild after reboot):
  /tmp/maxplus_p5.npy
  /tmp/maxminus_p5.npy   # if needed
  /tmp/e1_p7/maxplus.npy (+ Binv, C, coords, nb)
```

```bash
cd /home/nick/quadratic-minmax-limit
# sanity
python3 src/e1_gmin_cr_classify.py
python3 src/e1_gmin_uniform_lb.py
OMP_NUM_THREADS=1 python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py -q
# boolean Max+ rebuild if cache missing — use existing e1_* scripts, multi-core
```

---

## 5. Key file map

| Path | Role |
|------|------|
| `HANDOFF.md` | Long research handoff (Props through 15.49) |
| `solution.md` | Settled writeup; Main Theorem still OPEN sandwich |
| `evidence/E1_FAILURE_GRAPH.md` | F1–F14 banned loops |
| `evidence/E1_RESIDUAL.md` | Residual narrative through 15.47 |
| `evidence/e1_gmin_cr_classify.json` | CR classification cert |
| `evidence/e1_gmin_uniform_lb.json` | LB candidate vs threshold |
| `evidence/e1_gmin_closed_form_attack.json` | Dead closed-form attempts |
| `src/e1_gmin_cr_classify.py` | CR classifier |
| `src/e1_gmin_uniform_lb.py` | Algebra + known gmin check |
| `src/e1_gmin_char_sum.py` | Character-sum experiments |
| `src/e1_gmin_formula.py` / `linear_system.py` | Formula hunts |
| `tests/test_gmin_residual.py` | Residual unit tests (~80 lines) |
| `minmax_quadratic.py` | Core library |

---

## 6. Explicit non-goals tomorrow

- Soft-close Main Theorem without full chain  
- Claim \(m_n\ge\Phi-2\) from edge-minimal gap (F13)  
- Sequential MILP / serial SA cover enum (F1, F2)  
- Exhaust all PMs (F8)  
- Re-derive sandwich / denseness as deliverable (F9)  
- Treat cont Γ as discrete LB (F4)  
- Push 70 local commits without user request  

---

## 7. Suggested skills (next agent)

- `handoff` / `session-handoff-packager` — only if packing again  
- `use-available-compute` — before any Max+ rebuild or p=11 cert  
- `verification-before-completion` — before any “settled” claim  
- `goal-verifier` — if user asks whether L is done  

---

## 8. Acceptance criteria for “done”

| Milestone | Criteria |
|-----------|----------|
| g_min closed | Proof or theorem in `solution.md`: \(g_{\min}\ge L(p)\) (or better) ∀ primes \(p\ge5\), with evidence+test |
| Deep closed | Proof ND / non-undercut for deep non-tight on family |
| **Project settled** | Main Theorem \(L=\tfrac12\); HANDOFF SETTLED; full pytest green; no F3/F13 |

Until then: status line remains **L OPEN**.
