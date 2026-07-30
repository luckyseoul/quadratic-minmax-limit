# Session handoff — 2026-07-30 (Prop 15.95 / Path C residual)

**Repo:** `/home/nick/quadratic-minmax-limit`  
**Branch:** `main` (local often **far ahead** of `origin/main` — do **not** push unless asked)  
**L = lim α_n:** **OPEN** — never soft-close (F3)  
**P0 graph (mandatory every turn):** `evidence/P0_ENGINEERING_GRAPH.md`  
**Master handoff:** `HANDOFF.md`  
**Failure graph:** `evidence/E1_FAILURE_GRAPH.md` (F1–**F20**)  
**Problem:** [MO 413935](https://mathoverflow.net/questions/413935)

---

## 0. Skills to load at session start (do this first)

Load these **before** free-text thrashing or large file dumps:

| Priority | Skill path | Why |
|----------|------------|-----|
| **1** | `~/.grok/skills/agent-cost-optimization/SKILL.md` | **Token efficiency** — compact ops, no HANDOFF dumps, ponytail minimalism, graph+targeted greps only |
| **2** | `~/.grok/skills/use-available-compute/SKILL.md` | 88 cores + V100; ProcessPool `W≈nproc-2`; never serial thrash (F17/F20) |
| **3** | `~/.grok/skills/graph-engineered-completion/SKILL.md` | Update `P0_ENGINEERING_GRAPH.md` every material step; single next action |
| **4** | `~/.grok/skills/goal-verifier/SKILL.md` | Before any “closed / done” claim — evidence paths, not vibes |
| **5** | `~/.grok/skills/handoff/SKILL.md` | End of session / context risk — compact handoff + skills list |
| Optional | `~/.grok/skills/check-work/SKILL.md` | After shipping prop: tests + diff review |
| Optional | `~/.grok-installed-plugins/superpowers-.../verification-before-completion` | Same as goal-verifier if available |

**Compact workflow (from agent-cost-optimization):**
- No full `HANDOFF.md` / `solution.md` dumps — graph §8 + targeted greps only  
- Ship lemma + `evidence/*.json` + `tests/test_prop*.py` only  
- Scratch for this goal (if still active): use harness scratch only; **never** claim permanent `/tmp/maxplus_*.npy` as repo state (regenerate if missing)  
- Shell: **no `rg`** on this machine — use the grep tool / `python`

---

## 1. One-line resume

**P0:** Close Path C bi-tight residual for **all primes \(p\ge5\)** (proof, not census), then deep ND, then Main Theorem.  
**Just landed:** Prop **15.95** — Wick_hi ≤ thr_gap ∀p≥5; mult≥d + ∑M²≤Wick ⇒ gap; C_diag; mult=d cert p=3,5,7.  
**Active residual (pick one):**
1. Prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) ∀p≥5, and/or \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}=12n^2+48n\)  
2. Prove \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\) (16N)  
3. Prove \(\lambda_{\max}(FF^\top|_{1^\perp})\le N(3+H)\) (H)  

**Do not** re-attack ∑ρ κ_B (dead ≡ H after Prop 15.90).  
**L OPEN. H_proved=false. gap_proved_for_all_p=false.**

---

## 2. Settlement chain (load-bearing)

\[
m_n\ge\Phi-2\text{ on dense }\rho=1\text{ family }n=p^2+1
\;\Rightarrow\;\mathrm{E}(1)\;\Rightarrow\;L=\tfrac12
\]
Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) intact. Denseness alone ≠ closed limit (F3).

**Path C:** bi-tight empty (gap / 16N / H / g_min) → deep non-tight ND → Main Theorem.

---

## 3. What this session shipped (do not re-derive)

| Prop | Content | Evidence / code |
|------|---------|-----------------|
| 15.90 | Residual bound ≡ H; not independent foothold | `src/e1_gmin_m4_prop1590.py` |
| 15.91 | dim Z; orth/Φ/κ/harm ≡ H; 2×sphere⇒16N | `...prop1591.py` |
| 15.92 | ∑m4 κ_C constant; 16N/H ⇔ λ₂(P⊙P) bounds | `...prop1592.py` |
| 15.93 | FFT1=Nd1; 16N ⇔ λ_max(FFT\|1⊥)≤8N | `...prop1593.py` |
| 15.94 | P⊙P annihilates range(P); gap criterion mult≥d + ∑M²≤4d²(d+4) | `...prop1594.py` |
| **15.95** | **Wick≤thr ∀p≥5; mult+Wick⇒gap; C_diag; mult=d @ 3,5,7** | **`...prop1595.py`**, `evidence/e1_gmin_m4_prop1595.json`, `tests/test_prop1595.py` (6 green) |

**15.95 certs:**

| p | mult(λ₂) | ∑M²≤thr | gap_by_mult | 16N (numerical) |
|---|----------|---------|-------------|-----------------|
| 3 | =d=5 | no | no (correct) | yes |
| 5 | =d=13 | yes | **yes** | yes |
| 7 | =d=25 | yes | **yes** | yes |

**Algebra (15.95.1):** `thr−Wick = n(n²−16n−96)/2 ≥ 0` for `n=p²+1≥26` (p≥5).

**Identity:** `∑M² = 12n²−48n + C_diag + 24∑ρ²`, `C_diag=4n(11n−14)/p²`.

Writeups: `solution.md` Props 15.90–15.95; `HANDOFF.md` bullets; graph §6 log + §8 next.

---

## 4. Dead ends / bans

- **F3:** no soft-close L on sandwich + denseness  
- **F16:** do not pin modulus c by max PSD g_min  
- **F19:** no moduli class-key thrash  
- **F20:** label CPU multi-W / algebra if GPU unused; no GPU theater  
- **∑ρ κ_B:** dead as separate foothold (Prop 15.90)  
- CS / hypercube interlacing / unrestricted Gu op / two-moment **without** mult: too weak for 16N at p=5  
- Rank+two-moment **without** mult: fails gap at p=5; works at p=7 only  

---

## 5. Max+ arrays (ephemeral — regenerate if missing)

| p | Path | N |
|---|------|---|
| 3 | `load_maxplus(3)` in `e1_gmin_cr_classify` | 12 |
| 5 | `/tmp/maxplus_p5.npy` | 260 |
| 7 | `/tmp/e1_p7/maxplus.npy` | 11452 |

p=7 full N×N eig too large → matvec P⊙P + subspace iteration (see `certify_*` in prop1594/1595).

---

## 6. Next concrete steps (single focus)

1. **Primary:** Prove `mult(λ₂(P⊙P))≥d` for all primes p≥5  
   - Hint: Aut ≥ PSL(2,p²), q=p², has irrep of degree (q+1)/2 = d; PopP Aut-equivariant; λ₂-space is even Veronese (quadratic forms on V₊), mult=d at p=3,5,7  
   - Or explicit d-dim test space with Rayleigh ≡ λ₂  
2. **And/or:** Prove `∑M² ≤ Wick_hi` (boolean 4th moment ≤ Gaussian Wick on Max+)  
3. **Or:** 16N / H via λ_max(FFT\|1⊥)  
4. Ship as **Prop 15.96+**, evidence JSON, tests; keep **L OPEN** until full chain  
5. Only after bi-tight ∀p≥5: **N_DEEP**, then Main Theorem → G_L  

---

## 7. Suggested first commands for next agent

```bash
# After loading skills (§0):
cd /home/nick/quadratic-minmax-limit
sed -n '1,20p;160,180p' evidence/P0_ENGINEERING_GRAPH.md
python3 -c "import json;print(json.load(open('evidence/e1_gmin_m4_prop1595.json'))['settlement_note'])"
# Targeted tests only:
python3 -m pytest tests/test_prop1594.py tests/test_prop1595.py -q --tb=line
```

Do **not** start with full `HANDOFF.md` or full `solution.md` read.

---

## 8. Goal harness notes (if still active)

- Close bi-tight residual for **general** p≥5 with proof + evidence + tests  
- Certification only at {5,7} is **not** closed  
- G_L stays OPEN until E(1) or full Path C + deep  
- Scratch (if any): harness implementer dir only — not permanent repo state  

---

*Handoff written for next agent. Skills table is mandatory load at start.*
