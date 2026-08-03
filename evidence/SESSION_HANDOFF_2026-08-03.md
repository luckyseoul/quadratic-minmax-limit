# Session handoff — 2026-08-03

**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**MO:** [413935](https://mathoverflow.net/questions/413935)  
**L = lim α_n:** **OPEN** (do not soft-close)

---

## Current goal / status

Prize-level existence of \(\lim\alpha_n\). Sandwich + dense \(\rho=1\) Paley family are in. Settlement needs **Path C residual** or **E(1)** on \(n=p^2+1\), then Main.

**P0 — single open inequality (general primes \(p\ge5\)):**
\[
\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24
\;\Longleftrightarrow\;
\mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}
\;\Longleftrightarrow\;
\|\kappa\|_F^2\le\kappa_{\mathrm{hyp}}
\]
(equiv. forms: ED4≤ED4_bud, μ_G4≤μ_G4_suf).  
Then: **16N** (15.107+15.98) → bi-tight empty → deep ND → Main.

**Certified residual:** p=5 equality, p=7 strict. **Not proved for general p.**

Props shipped through **15.158** (code + tests + evidence JSON). Last arc: literature scan + residual attacks; no general-p close.

---

## Settled (do not re-derive)

| Piece | Where |
|--------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) | `solution.md` Main |
| ρ=1 halfspace on \(n=p^2+1\) | `PROOF_rho_eq_1.md`, `minmax_quadratic` |
| mult(λ₂)≥d−1 (PSL) | Prop 15.98 |
| mult + orth≤room_hyp ⇒ 16N | Prop 15.107 |
| δ-calculus, room_hyp closed forms | 15.102–15.117 |
| Q₄ closed; Max+ 1-hom, not IP-scheme | 15.157–15.158 |
| W census p=3,5,7; residual holds there | 15.128 |
| Bi-tight level-2 empty if \(g_{\min}>T(p)\); cert p=5,7 | 15.47 |
| n=10 undercutters = 144 PM + C6; m₁₀=Φ−2 | N10-S/C, 15.40–15.43 |

---

## Dead (F-graph + session)

Do **not** reopen: soft-close sandwich alone (**F3**); class_key thrash at p=7 (**F19**); GPU theater (**F20**); serial multi-minute (**F17**); IP association scheme of Max+; Chebyshev/μ≤1/moment-LP only; affine halfspace orbit as full Max+ (**F18**); pure t_e(k) freeness for p≥11; CS on δ·κ_B channel (too weak at p=5).

Full list: `evidence/E1_FAILURE_GRAPH.md`.  
Session blocker log: `evidence/RESIDUAL_BLOCKER_2026-08-03.md`.  
Lit scan: `evidence/LITERATURE_SCAN_2026-08.md` (MO still 0 answers; no arXiv close).

---

## Open residual surface (honest)

Prefer **one** of:

1. **General-p** δ²≤room_hyp/24 (or μ_G4≤μ_G4_suf / orth / κ²≤κ_hyp) — Weil/Gauss on Max+ or G-orbits (15.134), not prop re-encoding.  
2. **E(1):** no-descent for remaining gap-2 undercutters on ρ=1, or k_⋆=O(n^{3/2}); matching spike for all M at p≥5 still open.  
3. **g_min > T(p)** for all p≥5 → bi-tight empty without full 16N (still need deep non-tight ND).

**Do not** start Prop 15.159 thrash that only renames the residual.

---

## Next concrete steps

1. Resume from **this file** + `HANDOFF.md` §0 + **`evidence/STRATEGY_REFRAME_2026-08-03.md`** — not from chat.  
2. **Do not** start Prop 15.159 thrash. Prefer **structure identification** (what forces residual) over bounding δ.  
3. Highest continuity opens: Weil/Gauss **m₄ on G-orbits** (15.134); SDP dual for 16N/orth; G-irrep spectrum of E_{4p}^G; edge BM (not Max+ IP-scheme).  
4. Alternate: E(1) no-descent / k_⋆ / g_min>T(p) for all p≥5.  
5. Compute: `W=nproc-2`; Max+ caches `/tmp/maxplus_p5.npy`, `/tmp/maxplus_p7.npy`.

---

## Suggested skills

`graph-engineered-completion` · `use-available-compute` · `agent-cost-optimization` · `goal-verifier` · `handoff`

---

## Pointers

| Path | Role |
|------|------|
| `HANDOFF.md` | Master status |
| `solution.md` | Full writeup; Main still OPEN |
| `evidence/P0_ENGINEERING_GRAPH.md` | P0 graph |
| `evidence/E1_FAILURE_GRAPH.md` | F1–F20 |
| `evidence/LITERATURE_SCAN_2026-08.md` | External lit |
| `evidence/RESIDUAL_BLOCKER_2026-08-03.md` | Why not closed this session |
| `src/e1_gmin_m4_prop15102.py` … `15158.py` | Residual stack |
| `src/e1_gmin_m4_prop1598.py` / `15107.py` | mult / 16N chain |
| `tests/test_prop15*.py` | Prop tests |
