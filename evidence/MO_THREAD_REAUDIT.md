# MathOverflow 413935 — re-audit for hints and time savings

**Date:** 2026-07-27  
**Source:** [MO 413935](https://mathoverflow.net/questions/413935) (API + HTML; 0 answers, 11 comments)  
**X prize:** [PI010101/status/2081070728422752329](https://x.com/PI010101/status/2081070728422752329)

This note is a second pass over the **original** thread only (not X hallucinations). Goal: avoid re-work the author already framed as insufficient, and harvest any unused pointers.

---

## 1. What the question actually asks

> Does \(\lim_{n\to\infty} n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\bigl|\sum_{i<j}a_{ij}x_ix_j\bigr|\) exist?

Author constraints (verbatim intent):

- **No significant motivation** — pure curiosity.
- **Not interested** in numerical computations.
- **Not interested** in upper/lower bounds **unless** they give significant evidence the limit **does not** exist.
- If the limit exists, **what is the proof?**

**Time saving:** further pure SA / exact-\(m_n\) tables do not move the prize needle. Use numerics only as structural probes (e.g. N10 matching class), not as a substitute for E(1).

---

## 2. Author’s own comments (load-bearing)

| Comment | Content | Status in this repo |
|---------|---------|---------------------|
| Liminf | Author can show \(\liminf\alpha_n\ge C>0\); calculus gives \(C=2^{-5/2}\). Better constants possible; **these arguments do not prove/disprove existence**. | Prop 5.1 (BH / DMP); superseded floor Prop 5.2 dual-Gaussian \(1/\pi\). **Do not re-derive BH.** |
| Upper | Agrees limsup is finite (random \(A\) / Talagrand). | Conference / Paley \(\Rightarrow\) limsup \(\le1/2\). |
| Tags | Wants eigenvalues, spin-glasses, maxCut, Littlewood 4/3, polynomials — tag limit <6. | Seidel / \(\Phi\) / \(\rho\cdot\mathrm{op}\) programme matches “eigenvalues + hypercube + Fourier”. |
| DMP pointer | Points to Defant–Mastyło–Pérez via [blog 2019-12-06](https://extremal010101.wordpress.com/2019/12/06/) for a lower bound on 2-homogeneous forms. | Already Prop 5.1; blog is SK/Parisi + BH exposition, **not** a solution of the min-max limit. |

**Key meta-hint from author:** one-sided bounds (even with improved constants) were already known to him as **insufficient for existence**. Our sandwich is progress; it is still not settlement without E(1) (or a certified non-existence construction).

---

## 3. Other comments (ranked by usefulness)

### 3.1 domotorp — reformulation (keep)

Absolute value removable by symmetry:
\[
\tfrac12\lim\min_A\max_x n^{-3/2}\,x^\top A x
\]
with \(A\) symmetric \(\pm1\) off-diagonal (Seidel). **Already** `solution.md` §1.

### 3.2 Shannon Starr — multipartite / Bowlin (mostly dead for existence)

1. Initially suggested limit might be 0; retracted after author’s liminf.
2. Limsup finite by Talagrand (random \(A\)); min \(\le\) typical max.
3. Suggested recursive bipartition of \(K_n\) and focus on cross edges; cited Garry Bowlin, *Maximum Frustration in Bipartite Signed Graphs* (Electron. J. Combin. 19 (2012), #P10). Starr himself: “probably does not have good enough results… potentially interesting connection.”
4. Named possible audiences: Cris Moore, Lenka Zdeborová, Aukosh Jagannath (spin glass / large deviations).

**Time saving:** Bowlin is about **maximum frustration** of signed \(K_{\ell,r}\) (bipartite). Recursive bipartition / soft multipartite is already analyzed and **killed for existence** in `solution.md` §9–§10 (cross-block \(\Theta(N^{3/2})\) error). Do **not** reopen multipartite comparison as an existence proof. Optional only if a **certified** two-density construction for non-existence appears (none known).

### 3.3 Andrei Z. — OEIS

Submit \(m_n\) sequence. Harmless; does not settle the limit. Exact table \(n\le10\) already in `evidence/exact_m_table.json`.

### 3.4 Rafael Hipólito — antisymmetric confusion

Misread upper-triangle encoding; not useful.

---

## 4. X prize post (meta only)

- Author spent **~5 hours** without a solution → no cheap trick in the bare statement.
- Wants a **purely AI** transcript; will check the proof himself.
- Monetary prize for a **correct** solution (existence with proof, or non-existence with proof).

**Time saving:** do not “close OPEN” with sandwich + denseness + \(\rho=1\) alone (forbidden by plan and by author’s own “bounds don’t settle existence” stance).

---

## 5. Actionable checklist (do / don’t)

| Do | Don’t |
|----|-------|
| Attack **E(1)** on dense \(\rho=1\) family \(n=p^2+1\) (shortest path to \(L=1/2\)) | Re-derive BH \(2^{-5/2}\) or re-prove dual-Gaussian for its own sake |
| Or build a **permanent relative gap** \(\Phi-m_n\ge\varepsilon n^{3/2}\) on that family (disproves E(1); evidence against lim = 1/2) | Heavy numerics as the main deliverable (author not interested) |
| Use N10 structure only as **structural** constraint on proofs (matchings, Hamming-5, local opt) | Soft multipartite / Bowlin recursion for existence |
| Keep writeup **OPEN** until a load-bearing E(1) (or non-existence) argument lands | Claim limit from \(\rho=1\) alone |
| Ignore Parisi / graphon / Wick “proofs” from X | Treat random-SK constant \(\approx0.76\) as the min-max limit |

---

## 6. Mapping to current blockers

| Thread signal | Our status |
|---------------|------------|
| Existence is the question | **Still OPEN** (`HANDOFF.md`) |
| Author liminf \(2^{-5/2}\) | Improved to \(1/\pi\) (Prop 5.2) |
| Author: bounds ≠ existence | Sandwich + denseness + \(\rho=1\) still need **E(1)** |
| Starr multipartite | Dead path §9–§10 |
| DMP / Fourier tags | Used; not the E(1) engine |
| Eigenvalue framing | Active: \(\Phi=\tfrac12 n\rho\,\mathrm{op}\), product \(r\), rigidity programme Prop 15.20–15.21 |

**Bottom line after re-audit:** the original thread contains **no hidden proof sketch** of existence or non-existence. The only unused soft pointer (Bowlin / recursive bipartition) is already ruled out for settlement. Time is best spent on a **general E(1)** argument (rigidity + Lipschitz, or uniform \(\Phi\ge\Phi(C)-o(n^{3/2})\) on \(n=p^2+1\)), not on further one-sided bounds or multipartite comparisons.
