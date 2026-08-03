# Residual / lim α_n — honest blocker (2026-08-03)

## Goal acceptance vs reality

| Criterion | Required | Status |
|-----------|----------|--------|
| Residual δ²≤room_hyp/24 ∀p≥5 | general-p proof + tests | **OPEN** (census p=5,7 only) |
| Bi-tight via 16N chain | needs residual or ‖κ‖²≤96n | **blocked** |
| E(1) / deep ND | m_n≥Φ−o(n^{3/2}) on ρ=1 | **OPEN** |
| L closed in docs | no "L OPEN" | **still OPEN** (correct) |

## Attacks executed this goal arc

1. **Literature scan** (`evidence/LITERATURE_SCAN_2026-08.md`): MO 413935 unsolved; author prize wants existence; no arXiv close.
2. **Algebra inventory**: 15.98 mult≥d−1 proved; 15.107 mult+orth⇒16N proved; missing only orth≤room_hyp / δ²≤budget.
3. **Aut-line / halfspace c**: at p=5, ⟨δ,f_hs⟩=‖δ‖²=room_hyp/24 (dim Aut-ker effective 1). p=7 multi-orbit / non-scheme blocks class_key (F19). No Gauss-sum closed form for c found.
4. **Formula fit** for δ² from census (3,5,7): no simple rational candidate proved; ratios 0, 1, 0.327.
5. **Orbit structure**: Max+ ≠ single affine-Frob orbit of hs at p=5 (60/260); |PΓL|/N not integer-constant at p=7.
6. **E(1)**: nodescent OPEN; k_⋆=o(n^{3/2}) OPEN; dual-Gaussian gap is Θ(n^{3/2}) not o.

## Why not soft-close (F3)

Sandwich 1/π≤liminf≤limsup≤1/2 + ρ=1 denseness is already shipped and explicitly **insufficient** for existence (Prop 6.2 needs E(1) or residual path). Shipping L closed without residual/E(1) is F3.

## Why not another Prop 15.159

Re-encoding residual (μ_G4, R4, η, …) without a general-p inequality is the thrash pattern already run through 15.158. Dead ends catalogued in HANDOFF + E1_FAILURE_GRAPH F1–F20.

## What would unblock

- Closed form / Weil bound: μ_G4≤μ_G4_suf or δ²≤room_hyp/24 for all primes p≥5
- OR nodescent Lemma B on ρ=1 family (⇒ gap O(1) ⇒ E(1) ⇒ L=1/2)
- OR k_⋆=O(n^{3/2}) rigidity for Φ-minimisers on n=p²+1
- OR named external theorem importing one of the above

## Verdict

**GOAL NOT MET.** Prize-level open problem; no honest totality close available from current algebra + literature.

## Additional attacks (same session)

7. **CS on δ-channel** (15.111): ‖δ‖·max‖κ_B‖ exceeds (n−10)/(n−6) at p=5 — too weak.
8. **Fickus arXiv:2605.28738**: Singer–Zauner ETF gap — not a Max+ residual close.
9. **Strict Aut G Bose–Mesner** (15.134 already): dim E_{4p}^G = 0,2,7 at p=3,5,7; Aut-line dim≤1 fails for G; Gauss m4 on orbits still needs Max+ or closed form — affine-orbit m4 ≠ true m4 at p=5.
10. **Matching spike** (15.30a): implication proved; all-M criterion for p≥5 open; non-matching undercutters open.

**Still blocked:** general-p δ²≤room_hyp/24 (or equivalent). No soft-close.

## Dual-gap preferred target (Props 15.159–15.160)

Named structural objective (not δ thrash):

- Dual gap: \(G=(d/32)(16I-\Phi)\succeq I\) ⇔ \(\lambda_{\max}\le 16(d-2)/d\) ⇔ 16N.
- Hypothesis H: \(\mathrm{ray}_{\max}\le H(p)=(p+2)^2/d\).
- **Proved (15.160):** \(\mathrm{thr\_ray}-H=(3p+7)(p-5)/(2d)\), so for all primes \(p\ge5\),
  \(H\le\mathrm{thr\_ray}\) (eq only \(p=5\)). Thus **H ⇒ G≽I ⇒ 16N**.
- **Certified:** p=5 ray=H=thr_ray; p=7 ray<H<thray; G≽I at both.
- **Still OPEN:** H (or any residual equivalent) for general primes \(p\ge5\).

G-Schurian coherent config on Max+ at p=5 has **124** pair orbits (≫ 3 Φ levels),
so Schurian BM from strict Aut alone does not pin the Φ spectrum.

**GOAL NOT MET** for lim α_n totality; L remains OPEN.
