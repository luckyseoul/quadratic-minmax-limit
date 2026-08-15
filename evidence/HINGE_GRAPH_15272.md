# Residual-(i) hinge graph (2026-08-15)

## Nodes

| ID | Type | Status | Notes |
|----|------|--------|-------|
| L | claim | resolved | solution.md \(L=1/2\) |
| E1 | claim | resolved | type_I ∧ ii ∧ bitight |
| II | lemma | resolved | 15.179+236+237 |
| BITIGHT | lemma | resolved | 15.167 |
| RI | lemma | resolved | residual (i) via 15.249 ← 15.272 |
| KERSC | lemma | resolved | 15.207 ← gplus |
| GPLUS | lemma | resolved | 15.272 k=1∪k=3 span, not Aut-Schur |
| PDF | lemma | resolved | Singer PD on F, p≥7 |
| D249 | lemma | resolved | cost_D<2−α on sc |
| NU0 | lemma | resolved | 15.268 |
| P5VER | lemma | resolved | Veronese rank 65/65 |
| BADMU | lemma | resolved | 15.272 B + every triple (M3) |
| K1HYP | lemma | resolved | Johnson 1-swap + live DFT rank m−1 |
| COMP | lemma | resolved | complementary K_{m−1} mixed |
| L0 | lemma | resolved | through-L0 Fejer last dim |
| ISO | lemma | resolved | isotypic sum = dim F^⊥ |
| KINV | lemma | resolved | K=cot invertible (unused for close) |
| PAIR | lemma | open | 1ᵀK⁻¹v ≠0; unused |
| AUTSCH | lemma | failed | Jacquet ⇏ span; p=5 k=3 is 61/65 |
| GSUM | lemma | open | disj LB unused |

## Edges (depends-on)

K1HYP → GPLUS
BADMU → GPLUS
COMP → GPLUS
L0 → GPLUS
ISO → GPLUS
PDF → GPLUS
P5VER → GPLUS
GPLUS → KERSC → D249 → RI → E1 → L
II → E1
BITIGHT → E1
AUTSCH ↛ GPLUS
PAIR ↛ GPLUS
GSUM ↛ RI

## Open issues (not on the close)

- PAIR: cotangent pairing leftover; do not census further.
- GSUM disj LB unused.
- Path-C / 16N optional.
- Prize AI-test after writeup.
