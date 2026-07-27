# E(1) status on the \(\rho=1\) family (blocking for \(L=\tfrac12\))

**Date:** 2026-07-27  
**Claim E(1):** For \(n=p^2+1\) (Paley over \(\mathbb F_{p^2}\), \(\rho=1\), \(\Phi=\tfrac12 np\)),
\[
m_n=\Phi+o(n^{3/2})=\tfrac12 n\sqrt{n-1}-o(n^{3/2}).
\]
**Consequence:** Prop 6.1–6.2 \(\Rightarrow\lim\alpha_n=\tfrac12\).

**Status: NOT PROVED.** Existence of \(\lim\alpha_n\) remains **OPEN**.

---

## Evidence (certified upper bounds on \(m_n\), exact \(\Phi\))

| \(n\) | \(p\) | \(\Phi\) | certified \(m_n\le\) | gap \(\Phi-\mathrm{UB}\) | method |
|------:|------:|---------:|--------------------:|-------------------------:|:-------|
| 6 | — (q=5) | 5 | 5 | 0 | exact `exact_m` |
| 10 | 3 | 15 | 13 | 2 | exact Gray |
| 14 | — (q=13) | 21 | 21 | 0 | SA+exact; \(k\le4\) flips exact min \(\Phi=21\) |
| 18 | — (q=17) | 33 | 33 | 0 | SA+exact |
| 26 | 5 | 65 | 65 | 0 | intensive SA×10k (172 workers) + exact rescore; k-flips k=1..40 all \(\Phi\ge67\) |

Sources: `e1_paley_gap.json`, `e1_n26_intensive.json`, `n14_kflip.json`, `exact_m_table.json`.

**Trap:** `phi_local` at \(n=26\) can read 61 while exact \(\Phi=65+\). Only exact \(\Phi\) is valid for UBs.

## Structural facts useful for a future proof

1. **Edge-counting Lipschitz (Prop 15.20b, new):** \(k\) edge flips change each \(Q(\cdot)\) by at most \(2k\), hence \(\Phi(A)\ge\Phi(C)-2k\). **E(1) holds if \(k_\star=o(n^{3/2})\)** for \(\Phi\)-minimisers after switching (weaker than the old Frobenius demand \(k=o(n)\)). See `E1_EDGE_LIPSCHITZ.md`.
2. **Degree Lipschitz (Prop 15.20c):** max-degree \(D\) disagreement gives gap \(\le Dn\). Matchings: \(D=1\), gap \(\le n\).
3. **At \(n=10\):** optimizers at Hamming **5** from \(C\) (perfect matchings); gap \(2=o(n^{3/2})\). So the only known failure of *exact* optimality still obeys E(1). All 144 undercutting matchings have \(D=1\).
4. **Maximizers of \(\rho=1\) Paley:** \(\mathbb E[yy^\top]=I\) (exact 2-design; verified \(p=3,5\)). Adjacent-edge 4th moments vanish; disjoint-edge moments are not constant — not a full 4-design. Moment method alone does **not** force E(1).
5. **Local optimality:** single-edge flips raise \(\Phi\) under maximizer balance (Prop 15.21; holds for Paley \(n\le18\) and \(n=26\) k-flip samples).

## What would finish E(1)

- **Rigidity (edge form):** every \(\Phi\)-minimiser on \(n=p^2+1\) is switching-close to Paley/conference at Hamming scale \(k_\star=o(n^{3/2})\); then Prop 15.20b \(\Rightarrow\) E(1) \(\Rightarrow\lim\alpha_n=\tfrac12\).
- **Or** uniform bound \(\Phi(A)\ge\Phi(C)-o(n^{3/2})\) for all Seidel \(A\) on these orders (stronger).
- **Or** permanent relative gap \(\Phi-m_n\ge\varepsilon n^{3/2}\) for infinitely many such \(n\) (disproves E(1); then limsup \(\alpha\le(1-\varepsilon)/2\) along the family).

None of these is established. **Do not mark \(\lim\alpha_n\) settled.**
