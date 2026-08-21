# Route kills (2026-08-21) — no leftover flag flipped

One-line: leftover 3 is L∞ over Aut_e classes (L2 conversion fails);
a uniform Paley bound E_-[S²]<20+12/p is false; QVAR Wick already
beats the exceptional threshold but does not close leftover 1.

## Leftover 3 is not an L2 problem

3A+B>0 is required on every Aut_e far class. Y+3 is linear in the
cover weights μ_O≥0, so the worst case is a single class. An L2 bound
on μ (including tr(Φ²)=4‖M‖_F²−3n²+2n²(n−1)/p²) cannot replace
|μ|≤L. Claude-referee suggested that conversion; it is rejected.

## Uniform Paley majorant E_-[S²]<20+12/p is false

15.379 C hoped E_-[S²]<20+12/p would empty 4-level leftover (lattice
min E[S²]=20+12/p). The explicit p=5 leftover double-star (15.406 D,
20 edges, support {−12,−8,−4,−2}) has

    E_-[S²] = 25.169… > 22.4 = 20+12/p.

It is leftover-only (min_+=−4), not leftover+splus. It still kills any
majorant claimed for every G of size 4p. A 4-level-only majorant is
circular. Do not reopen 15.379 C as a uniform bound.

## Delsarte 2-design + min-distance does not give principal room

Max+ is an antipodal spherical 2-design in V_+ with d_H≥p+1
(15.197). The LP max of E[(y·z)⁴] under only those constraints
exceeds the leftover-1 principal room: 12688>9256 at p=5,
116187>31400 at p=7. Paley structure beyond 2-design+distance is
required. Do not reopen unprojected Delsarte for ‖δ‖².

## QVAR Wick already exceeds the exceptional threshold

A_ψ=P K_ψ P/4 with K_ψ(a,b)=ψ(a−b) on finite coordinates
(15.589 H). Certified ‖A‖_HS²=q(q−1)/32. Wick 8‖A‖²=q(q−1)/4
strictly exceeds 3q(q−1)/16. Live E|yᵀ A y|² equals the census
E|Z_ψ|² (3300/13 at p=5, 317520/409 at p=7) and sits above Wick
(253.8>150, 776>588). A nonnegative 4-harmonic excess would close
the exceptional scalar only. Leftover 1 still needs principal room.
λ_exc is the top Φ-eigenvalue at p=5,7,11; that ordering is not a
theorem (15.589 C).

## Leftover 2 Walsh does not empty the witness family

Interior 4-level needs both residues mod 4 on Uᶜ (15.406). The p=5
leftover witness is same-residue on Uᶜ ({−12,−8,−4}≡0 mod 4) and is
Walsh-compatible. Walsh, even if proved for all p, does not empty
5+/same-residue leftover, and must not flip residual_ii alone.

## Status

All three leftover predicates stay False. Aut-Schur / Gsum / pairing
stay False. L OPEN.
