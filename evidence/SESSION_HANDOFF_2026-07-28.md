# Session handoff — 2026-07-28 (shutdown / ~1% quota)

**For:** next agent / tomorrow morning  
**Repo:** `/home/nick/quadratic-minmax-limit`  
**Branch:** `main` (ahead of origin by ~65+ commits; uncommitted residual scripts were staged this session)  
**L = lim α_n:** still **OPEN** — do **not** soft-close Main Theorem.

---

## 0. One-line status

Settlement chain still open:

\[
m_n\ge\Phi-2 \text{ on }n=p^2+1
\;\Rightarrow\; \mathrm{E}(1)
\;\Rightarrow\; L=\tfrac12
\]

via denseness (Prop 6.2). Equivalence: \(m_n\ge\Phi-2\) iff no-descent on all gap-2 undercutters (Props 15.40–15.42 + parity).

**Locked (do not re-prove):** Props 15.40–15.47, sandwich, ρ=1, denseness, F1–F14 failure graph.  
**Certified only at p=5,7:** \(g_{\min}\) above bi-tight threshold ⇒ Type I freeness-failure + deep-tight dead.  
**Residual for all p≥5:** (i) closed-form \(g_{\min}\) bound; (ii) deep **non-tight** gap-2 ND / spike.

---

## 1. What was interrupted (this residual attack)

After commit `4c0d5e4` (handoff package), the session was attacking residual E(1) with four new scripts (now under `src/` + `evidence/`):

| Script | Evidence | Result |
|--------|----------|--------|
| `src/e1_gmin_formula.py` | `evidence/e1_gmin_formula.json` | Identities: wedge \(G=\pm1/p\); \(C^3\) ⇒ sum of wedge \(G=0\); avg disj \(G>0\); cert \(g_{\min}>\mathrm{thresh}\) at p=5,7. **Closed form still open.** Affine orbit of halfspace **undersamples** Max+ (unreliable alone for \(g_{\min}\)). |
| `src/e1_gmin_linear_system.py` | `evidence/e1_gmin_linear_system.json` | At p=5: bare \(C\)-isomorphism 4-set types do **not** make \(G\) constant (only 1/64 edgepair types constant). Need **finer** scheme relations (not just edge \(C\)-sign patterns). \(g_{\min}=-3/65\) still. |
| `src/e1_gmin_char_sum.py` | (may need re-run) | Affine+Frob orbit of halfspace vs full Max+. Status string: only valid if `orbit_complete` and \(g_{\min}\) matches full — expect **incomplete** (same F as formula script). |
| `src/e1_deep_spike_theory.py` | `evidence/e1_deep_spike_theory.json` | Prop 15.46 thresholds tabulated. Pure degree pigeon at \(k=2p\) **insufficient** (avg deg \(<1\) for \(p\ge5\)). p=5 cover census empty this run (`missing Max± cache` at write time; caches later rebuilt under `/tmp/`). |

**Caches (volatile — wiped on reboot):**
- p=5 Max±: `/tmp/maxplus_p5.npy`, `/tmp/maxminus_p5.npy`
- p=5 deep: `/tmp/e1_deep_p5/{chip,chim}.npy`
- p=7 Max±: `/tmp/e1_p7/maxplus.npy`, `maxminus.npy` (rebuild via `src/e1_gmin_p7.py` / bitight scripts)

**Prior session facts still valid:**
- p=5: \(g_{\min}=-3/65\), \(h_{\min}=-6/65=2g_{\min}\); discrete disj \(G\in\{\pm1,\pm3,9,13,21\}/65\)
- p=7: \(g_{\min}\approx-0.03807\); bi-tight blocked
- 20 tight L2 Max+ covers at p=5: all \(\max S_-\ge4\)
- Deep large-\(k\) p=5: \(\Phi\in\{75..83\}>\Phi(C)=65\) (spike)

---

## 2. Dead ends (do not re-run)

From F-graph + this session + prior residual attack:

1. Soft-close Main Theorem / \(m_n\ge\Phi-2\) from 15.40 alone (**F13**)
2. 4-point boolean LP with pairwise \(E[y_ay_b]=C_{ab}/p\) → ~−0.6 (too weak)
3. Same + Paley patterns + \(C^2\) + projector → still ~−0.6
4. Chebyshev / 1-outlier from disj mean+var → ≪−1
5. Wick \(-1/p^2\) as **lower** bound (false: actual more negative)
6. Claiming \(G\) constant on bare \(C\)-isomorphism 4-set types alone (false at p=5)
7. Affine / Frob halfspace orbit alone for \(g_{\min}\) (undersamples Max+)
8. Pure degree pigeon for deep spike at \(k=2p\) (avg deg <1)
9. Serial MILP enum of all undercutters (F1); ProcessPool from stdin (fails)

---

## 3. Promising next edges (priority)

### A. Closed-form \(g_{\min}\) (preferred — closes Type I + deep-tight for all p≥5)

1. **Association scheme on edge-pairs** with relations finer than \(C\)-sign 6-tuples: Johnson geometry on 2-subsets + Paley/2-graph relations (or strongly regular parameters of the conference graph, lifted to edges). At p=5 Gram spectrum of \(G\) was partially seen: eigenvalues \(n/2, 88/13, 72/13, 40/13, 0\) — extract **min disjoint class** from eigenmatrix / intersection numbers.
2. Hypothesis to check: \(g_{\min}\ge -3/\Phi = -3/(p n/2)\)? At p=5 equality \(-3/65\); at p=7 check against \(-436/11452\). If true and \(-3/\Phi > -(p-2)/(p(2p-1))\) for all \(p\ge5\), Prop 15.47 applies universally.
3. Full character sum over **all** Max+ (not halfspace orbit): use known 2-design \(\mathbb E[yy^\top]=I+C/p\) + higher design if any; or enumerative formula from AG(2,p) boolean evecs structure (counts: 12, 260, 11452 at p=3,5,7).

### B. Deep non-tight gap-2 (needed even after g_min)

1. Force Max− 1-bit: when \(s_-=-2\), existence of \(z\in\mathrm{Max}_{-}\) with some \(\tau_v\le-(p+1)/2\) (Prop 15.46) ⇒ \(\Phi\ge\Phi-2\).
2. Do **not** rely on degree pigeon at minimal \(k=2p\). Use level-set geometry of \(S_F\equiv -2\) on Max−, or matching-cover style spike that already appears in large-\(k\) census.
3. Re-run `e1_deep_spike_theory.py` after Max± rebuild to fill `p5_covers` with actual min-τ stats.

### C. Only after A+B close the chain

Update Main Theorem in `solution.md` / top of `HANDOFF.md` to \(L=1/2\); add evidence + tests; full pytest. **Not before.**

---

## 4. Concrete morning checklist

```bash
cd /home/nick/quadratic-minmax-limit
git status -sb && git log --oneline -5

# Rebuild caches if missing
ls /tmp/maxplus_p5.npy /tmp/maxminus_p5.npy /tmp/e1_p7/ 2>/dev/null

# Sanity
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -m pytest tests/test_minmax.py -q

# Resume residual (file scripts only; W≈nproc-2)
python3 src/e1_gmin_formula.py          # refresh identities if needed
python3 src/e1_gmin_linear_system.py    # p=5 types already known weak
# Next real attack: association scheme extraction OR deep tau forcing
# Consult evidence/E1_FAILURE_GRAPH.md before multi-minute jobs
```

Read order: `HANDOFF.md` §0 → this file → `evidence/E1_RESIDUAL.md` → `evidence/E1_FAILURE_GRAPH.md` → Props 15.45–15.47 in `solution.md`.

---

## 5. Load-bearing file map

| Path | Role |
|------|------|
| `HANDOFF.md` | Master resume entry (status + prop log) |
| `solution.md` | Full writeup; Main Theorem still sandwich+OPEN |
| `evidence/E1_FAILURE_GRAPH.md` | F1–F14 — **read before thrashing** |
| `evidence/E1_RESIDUAL.md` | Residual narrative through 15.47 |
| `evidence/e1_bitight_gsum_obstruction.json` | p=5 bi-tight Gsum cert |
| `evidence/e1_gmin_p7.json` | p=7 g_min cert |
| `evidence/e1_deep_sweep_p5.json`, `e1_deep_cover_phi.json` | deep census / spike |
| `tests/test_minmax.py` | Props 15.40–15.47 tests (last: 68 green on that arc) |
| `src/minmax_quadratic.py` | Library: Paley, halfspace, phi, frames |

---

## 6. Acceptance criteria (unchanged)

- [ ] **ac1:** Settle \(L=\lim\alpha_n\) with load-bearing proof (or honest OPEN)
- [ ] **ac2:** \(g_{\min}\) closed form all \(p\ge5\) **or** deep non-tight ND → \(m_n\ge\Phi-2\)
- [ ] **ac3:** Update Main Theorem/HANDOFF **only if settled**; keep sandwich/ρ=1
- [ ] **ac4:** Evidence+tests for new claims; full pytest; verification logs

**Current:** ac1–ac2 open; ac3 correctly left unsettled; ac4 partial (Props through 15.47 tested; residual scripts untested as load-bearing).

---

## 7. Suggested skills / discipline

- `use-available-compute` before heavy jobs  
- Project `E1_FAILURE_GRAPH.md` before multi-minute compute  
- `verification-before-completion` if claiming settle  
- **Never** soft-close Main Theorem  

## 8. Do not

- Soft-close \(L\) or \(m_n\ge\Phi-2\) from 15.40 alone  
- Re-run p=5 bi-tight MILP / known dead LP/Chebyshev/Wick-LB  
- Exhaust matchings / serial SA as deliverable  
- ProcessPool via stdin; peg one core when dozens free  

---

*Handoff written under ~1% quota after unexpected shutdown. Mid-flight residual scripts + JSON evidence committed or staged with this session so tomorrow does not re-derive dead ends.*
