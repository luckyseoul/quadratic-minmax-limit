# Research handoff: min-max ±1 quadratic form limit

**Status date:** 2026-07-26 (n=10 structure campaign)  
**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**Problem source:** [MathOverflow 413935](https://mathoverflow.net/questions/413935) / [X prize post](https://x.com/PI010101/status/2081070728422752329)

---

## 0. One-line status

**Existence of the limit \(L=\lim\alpha_n\) remains OPEN.**  
Proved sandwich: \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\).  
Dense Paley subsequence with \(\rho=1\) (orders \(n=p^2+1\)) is proved; asymptotic optimality of \(m_n\) vs Paley \(\Phi\) is not.

**New (n=10 structure):** exact optima first appear at Hamming distance **5** from Paley \(C_{10}\), and the only 5-edge undercutters are **144 perfect matchings** (of 945). Absolute gap \(\Phi-m_{10}=2\) is consistent with E(1). See `evidence/N10_STRUCTURE.md`.

**New (n=10 classification, N10-C):** those 144 matchings are exactly one \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit, equivalently the matchings that drop every Paley maximizer to \(|Q|\le13\) (six \(+\) maximizers certify). See `evidence/N10_MATCHING_CLASSIFY.md`. Existence of \(L\) remains **OPEN**.
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
| **Thm E (conditional)** | E(1)+E(2) \(\Rightarrow\lim\alpha_n=\rho_*/2\) | Hypotheses open. Along \(\rho=1\) family, E(1) alone would give \(\limsup\alpha=1/2\); still need matching liminf for \(L=1/2\). |
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

1. **E(1) on \(\rho=1\) family** — Prove \(m_n=\frac12 n\sqrt{n-1}-o(n^{3/2})\) for \(n=p^2+1\), or exhibit a permanent relative gap. At \(n=10\) the absolute gap is only \(2\) (rel. \(\approx0.063\)); at \(n=26\) SA+exact-rescore found **no** undercutter of \(\Phi=65\) (certified \(m_{26}\le65\)). Need a general argument, not local edge-flip (optima sit at Hamming \(\ge5\)).  
2. **Structural gap from \(n=10\)** — **Mostly closed (N10-S + N10-C).** Hamming-5 threshold; only undercutters at \(k=5\) are 144 perfect matchings; those 144 form one \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit and equal the maximizer-drop set (`evidence/N10_STRUCTURE.md`, `evidence/N10_MATCHING_CLASSIFY.md`). SA also finds Hamming-11–16 optima in the switching metric (same \(r=13/15\)). Remaining: whether a matching-type construction lifts (random matchings at \(n=26\) **raise** \(\Phi\) to \(\ge73\)); classify non-matching distant optima.  
3. **E(2) analytic** — Prove \(\rho(C_n)=1-O(n^{-1/2})\) (or \(\to1\)) for all large Paley, not only \(p^2+1\). Interval constructions are evidence, not a proof.  
4. **Non-existence** — Only if two dense subsequences with **proved** unequal \(\alpha\) limits appear; denseness (Prop 6.2) is mandatory.
### 5.2 Traps to avoid

| Trap | Why |
|------|-----|
| Soft multipartite / Hadamard “\(c_k\to0\)” existence | Cross-block error is \(\Theta(N^{3/2})\); abstract sequences can oscillate inside the sandwich |
| Q4 / fourth-moment shell for large \(n\) | Shell vacuous when \(\Delta_*/3>\max\delta\) (`solution.md` Props 15.16–15.19) |
| Treating `phi_local` min as \(m_n\) UB | Local max underestimates \(\Phi\) |
| Uncertified \(\rho\) dips as non-existence | Need certified \(\Phi\), not SA LB |
| Exact conference optimality for all \(n\) | **False** at \(n=10\) |
| Claiming \(L=1/2\) from \(\rho=1\) alone | Needs E(1) **and** liminf \(\to1/2\) (dual-Gauss only gives \(1/\pi\)) |
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
- If E(1): \(\alpha_{n_k}\to1/2\) along that dense family \(\Rightarrow\limsup\alpha_n=1/2\).  
- Full \(L=1/2\) still needs \(\liminf\ge1/2\) (stronger than dual-Gauss \(1/\pi\)).  
- If E(2) gives \(\rho_*<1\) and E(1), conditional Thm E says \(L=\rho_*/2\).

None of this is established. **Leave Main Theorem as sandwich + OPEN.**

---

## 7. Verification commands

```bash
cd /home/nick/quadratic-minmax-limit
python3 -m pytest tests/test_minmax.py -v
python3 -c "from src.minmax_quadratic import dual_gaussian_lower_bound; print(dual_gaussian_lower_bound(10))"
# expect ~ 9.5493 = 10*3/pi
python3 /tmp/grok-goal-*/implementer/verify_rho1.py  # or: pytest -k rho_eq_1
```

Expected: all tests pass; existence not claimed settled in `HANDOFF.md` or `solution.md` Main Theorem.
