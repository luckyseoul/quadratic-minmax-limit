# No-descent candidate toward m_n ≥ Φ−2

**Date:** 2026-07-28
**Status:** Structural lemmas partial; n=10 cert for matching undercutters;
**forall / general n OPEN. Existence of lim α_n remains OPEN.**

## Goal (graph-allowed path)

Upgrade Prop 15.40 → global \(m_n\ge\Phi(C)-2\) **without** F13 shortcut.

## Lemmas

### Lemma A (first-hit; proved)
Along any edge-adding chain, the first undercutting prefix F has
Φ(C⊕F) ≥ Φ(C)−2 (Prop 15.20b with k=1).

### Lemma B (no-descent; OPEN in general)
If Φ(C⊕F)=Φ(C)−2, then for every edge e∉F, Φ(C⊕(F∪{e}))≥Φ(C)−2.

If A+B hold on the ρ=1 family, induction on Hamming distance yields
m_n≥Φ(C)−2, gap O(1)⇒E(1)⇒L=1/2 by denseness.

**F13 reminder:** abstract 2-Lipschitz functions need not satisfy B.
The claim is special to Φ of Seidel matrices.

## Certified at n=10 (matching undercutters)

| Check | Result |
|-------|--------|
| # PM undercutters | **144** (all Φ=13) |
| add-1 violations Φ<13 | **0** |
| min Φ over add-1 | 15.0 |
| multi-edge sample violations | **0** |
| multi-edge sample min Φ | 13.0 |
| workers | 86 |

JSON: `evidence/e1_n10_nodescent.json`. Script: `src/e1_n10_nodescent.py`.

## Dangerous-edge obstruction (partial structure)

Write A=C⊕F, Φ(A)=Φ−2, edge e=(u,v), B=A⊕e.
Then Q_B(x)=Q_A(x)−2σ_x with σ_x=A_uv x_u x_v ∈ {±1}.

If some maximiser x with Q_A(x)=Φ−2 has σ_x=−1, then Q_B(x)=Φ,
so Φ(B)≥Φ≥Φ−2. Similarly for negative maximisers with σ=+1.

**Descent to Φ−4 requires** rigid alignment σ≡+1 on all + maximisers
and σ≡−1 on all − maximisers. Open: prove this alignment is impossible
for undercutters of ρ=1 conferences (or non-maximiser spikes save Φ−2).

## What this does **not** prove

- Lemma B for general n (only n=10 PM cert)
- Clique-flip / matching non-undercut for p≥5
- k_⋆ / E(1) / existence of L

**Do not mark Main Theorem settled.**
