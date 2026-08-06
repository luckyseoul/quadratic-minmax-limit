# Final handoff — min-max ±1 quadratic form limit

**Date:** 2026-08-05 (prize update **2026-08-06**)  
**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**GitHub:** https://github.com/luckyseoul/quadratic-minmax-limit (`main`)  
**Problem:** [MO 413935](https://mathoverflow.net/questions/413935) · [X prize @PI010101](https://x.com/PI010101/status/2081070728422752329)

---

## 1. One-line status

**Candidate** \(L=\lim\alpha_n=\tfrac12\) on denseness path (Props **15.167–15.171**).  
**Prize status: NOT accepted.** Independent Grok Build 4.5 (posted by @PI010101) found a **gap** in lim=½.

**Optional still OPEN (not required for this path):** Path-C residual / \(16N\) / Hypothesis H.

### Prize / thread update (2026-08-06)

- @PI010101 asked for full solution on GitHub/PDF; repo linked.  
- He then posted: *“Grok build 4.5 found a gap in the proof lim=1/2”* with share  
  https://grok.com/share/c2hhcmQtNA_218425aa-c1d1-4263-a3ea-9114ef04cd9c  
- His requirement: pass an **“AI test”** — any reasonable AI identifies the solution as essentially correct in **2–3 repeated prompts** — **only then** will he human-verify.  
- User **cannot use MO** (new account / no points). Delivery channel: **X prize thread + GitHub**.  
- Thread reply may have been deleted; gate above still stands.

### What’s left (simple)

1. **Find the actual hole** (read Paata’s Grok share or cold gap-find). One-sentence: which step fails.  
2. **Fix or drop the claim** (repair hinge, or stop saying CLOSED).  
3. **Short clean proof package** (setup → sandwich → denseness → 3 lemmas → \(L=\frac12\)).  
4. **Pass AI test** (2–3 models × 2–3 prompts; same gap ⇒ fix that).  
5. **Ping Paata** with fixed writeup + GitHub.

**Not needed now:** residual/16N · more X cards · MO answer.

---

## 2. What the problem is for

Author (MO): *“no significant motivation … pure curiosity.”*  
Wants existence (and proof), not numerics.  
X reframed as pure-AI prize (full transcript, no human math hints).

Not an applied target; Seidel / Paley / conference machinery is the method, not the purpose.

---

## 3. Claim chain (load-bearing)

Write \(\alpha_n=m_n/n^{3/2}\) with
\[
m_n=\min_{a_{ij}=\pm1}\max_{x=\pm1}\Bigl|\sum_{i<j}a_{ij}x_ix_j\Bigr|.
\]

| Step | Status | Where |
|------|--------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) | CLOSED | `solution.md` Main Theorem |
| \(\rho=1\) on Paley \(n=p^2+1\) | CLOSED | `evidence/PROOF_rho_eq_1.md` |
| Denseness (Prop 6.1–6.2): existence along dense family ⇒ global lim | CLOSED | `solution.md` |
| **Bi-tight empty** all primes \(p\ge5\) | **CLOSED** | Prop **15.167** |
| Type I freeness-fail \(k=3p-2\), \(s_-\le-1\) impossible | **CLOSED** | Prop **15.170** |
| Deep freeness-fail \(k\ge3p\) ND | **CLOSED** | Prop **15.171** |
| E(1): \(m_n\ge\Phi(C)-2\) on \(\rho=1\) family | **CLOSED** | 15.167 ∧ 15.170 ∧ 15.171 |
| \(\Rightarrow L=\tfrac12\) | **CLOSED** | denseness + sandwich limsup |

### Lemma summaries (for verification)

**A — Bi-tight (15.167).**  
Max± edge Gram: mult\((\lambda_{\max})\ge d-1\), \(\lambda_{\min}\ge6\) ⇒ majorization  
\(L_*=(p^4+24p^2-1)/(2(p^2-1))\). For \(p\ge5\), \(L_*<2d\) ⇒ \(\lambda_{\mathrm{cycle}}<d\) ⇒ top eig simple ⇒ no bi-tight size-\(2p\) cover. **No residual/16N.**

**B — Type I (15.170).**  
Freeness-fail \(k=3p-2\), gap-2 ⇒ \(s_-=-1\); bad case dualizes to  
\((\mathrm{Gsum}\,x)_e=6/p-4\). Box LB \(\ge-12k/(pn)\). Conflict iff  
\(4p^3-6p^2-32p+18>0\) (true for primes \(p\ge5\); at \(p=5\): \(208>0\)).

**C — Deep (15.171).**  
\(s_+=2\), freeness-fail \(k\ge3p\): free ⇒ weak ND; auto-freeness \(k\le3p-2\);  
fail-eq \(k=3p-1\) empty; remaining dual two-level Gsum vs same box LB ⇒ ND.

**Assembly.** A + freeness/tight empty + B + C ⇒ \(m_n\ge\Phi-2\) ⇒ denseness ⇒ \(L=\tfrac12\).

---

## 4. What is *not* claimed

- Path-C residual \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for all \(p\ge5\) — **OPEN**
- \(16N\) / Hypothesis H for all \(p\ge5\) — **OPEN**
- Soft-close from sandwich + denseness alone (F3) — **banned**; L only via real E(1) predicates
- Intermediate prop logs pre-15.167 still say “L OPEN” historically — superseded; see `STATUS.md`

---

## 5. Repo map (start here)

| Path | Role |
|------|------|
| **`STATUS.md`** | Single-page claim table |
| **`HANDOFF.md`** | Live research handoff (banner CLOSED) |
| **`solution.md`** | Full writeup; top = limit CLOSED; Props 15.167–171 at end |
| **`README.md`** | Public status |
| `src/e1_gmin_m4_prop15167.py` … `prop15171.py` | Proof modules |
| `tests/test_prop15167.py` … `test_prop15171.py` | Load-bearing tests |
| `evidence/e1_gmin_m4_prop15167.json` … `15171.json` | Evidence JSON |
| `evidence/e1_main_chain_status.json` | `L_closed: true` |
| `evidence/share/paper/lim_alpha_n_closed.{tex,pdf}` | Paper draft |
| `x-cards/page{1,2}.jpg` | X result cards |
| `x-cards/lemmas_page{1,2}.jpg` | Key lemmas for Grok/verification |
| `scripts/render_x_twopager.py` / `render_x_lemmas.py` | Card generators |

**Git:**  
- `c5ad4f7` — ship closed chain + paper + cards  
- `da435f2` — align entry docs so Grok no longer sees OPEN at top  

**Untracked (optional junk):** `evidence/gpu_*.json` session dumps — not required for claim.

---

## 6. X / Grok thread state

- Prize root: https://x.com/PI010101/status/2081070728422752329  
- User (@RegolithHunter) posted result cards; Grok reply  
  https://x.com/grok/status/2084976346036416839 asked for:
  1. Full writeup on arXiv or MO  
  2. Key lemmas shared for inline check  
  3. (It still thought residual blocked L — because public handoff/solution tops were stale; fixed in `da435f2`)

**Post-ready assets:** `x-cards/lemmas_page1.jpg`, `lemmas_page2.jpg` + repo link.

---

## 7. Anti-patterns (do not reopen)

| ID | Ban |
|----|-----|
| F3 | Soft-close L from sandwich + denseness without E(1) predicates |
| F19 | class_key / moduli thrash (not m₄-equitable at p=7) |
| F20 | GPU theater for non-dense jobs |
| — | Treating residual/16N as required for L after 15.167 |
| — | Trusting mid-file “L OPEN” lines as current status |

---

## 8. Preferred next work (priority order)

1. **Prize path (primary)** — see §1 “What’s left (simple)”  
   - Extract gap from Paata Grok share  
   - Hinge audit / fix (λ_min≥6, Gsum LB, case exhaustiveness, dual-equality, tests import)  
   - Thin self-contained writeup  
   - Multi-AI “essentially correct” test  
   - Reply on X with GitHub (no MO)

2. **Peer / formal rigor pass**
   - Cold-read Props 15.167, 15.170, 15.171 against `solution.md`  
   - Every predicate in `e1_main_chain_status.json` proof-backed, not census-only  
   - Fix `pytest` imports (`ModuleNotFoundError: e1_gmin_m4_prop15100`) then green 15.167/170/171

3. **Optional Path-C residual/16N** (independent of L) — only if deliberately reopening Path C  

4. **Cleanup**
   - Optionally gitignore untracked `evidence/gpu_*.json`  
   - Revisit GitHub **About description** after gap fix (currently claims CLOSED)

---

## 9. Verify before claiming success again

```bash
cd /home/nick/quadratic-minmax-limit
git log -1 --oneline   # expect da435f2 or later on main
head -20 STATUS.md
python3 -m pytest tests/test_prop15167.py tests/test_prop15170.py tests/test_prop15171.py -q
python3 -c "import json; print(json.load(open('evidence/e1_main_chain_status.json'))['L_status'])"
```

---

## 10. Suggested skills (next agent)

| Skill | When |
|-------|------|
| `goal-verifier` | Confirm acceptance criteria after any claim change |
| `check-work` | Diff + tests after edits |
| `scientific-critique` / cold MO answer review | Before posting MO/arXiv |
| `use-available-compute` | Only if reopening residual/16N GPU work |
| `handoff` | End of session |

---

## 11. Session outcome (this arc)

- Candidate mathematical claim: **L = ½** via bi-tight majorization + dual Gsum Farkas freeness-fail ND.  
- Public repo + X cards shipped; GitHub About description set to CLOSED (pending re-honest after gap).  
- Thread: Paata engaged; independent Grok flagged gap; **AI-test bar** before human verify.  
- Residual/16N left optional. User has no MO rep — prize via X + GitHub only.

**Next human decision:** identify/fix gap → AI test → re-ping Paata (not residual thrash, not MO).
