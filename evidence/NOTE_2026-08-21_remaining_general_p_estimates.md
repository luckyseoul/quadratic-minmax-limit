# Isolated remaining general-p estimates (no flag flipped)

Date: 2026-08-21. No leftover predicate is True. `phi_F_ge_6` imports
`leftover1_qvar_and_principal_proved` (`src/e1_gmin_leftover1_qvar_principal.py`);
that AND is False. p=13 k=7 orbit exhaustion is not a leftover close. Aut-Schur / Gsum disj LB / pairing
stay False. L=1/2 stays OPEN.

The three leftovers reduce to the following Max+-free estimates. Census
at p≤7 / p=11 Max+ / p=13 orbits cannot be the only reason a leftover
is True.

## Leftover 1 (`phi_F_ge_6_proved_general`)

Two independent blocks, both open:

1. QVAR on every k≥7 stratum, all primes p≥13:
   `E|Z_ψ|^2 ≥ 3q(q-1)/16`. Closed on k=1..6 for every prime (15.589).
   Pointwise and orbitwise QVAR are false. Coarse profile constraints
   cannot prove it (15.589 I).
2. After QVAR, principal room
   `||δ||^2 ≤ n(n+10)^2 / [6(n-14)(n-6)]`, equivalently
   `E[(y·z)^4] ≤ 4n(3n^2-37n+2)/(n-14)`.
   Crude `E[s^4]≤2n^3` is Θ(n^3) vs a Θ(n^2) budget (15.120 / 15.156).
   fable.md "crude bound closes leftover 1 for p≥7" is slack folklore,
   not a theorem.

## Leftover 2 (`residual_ii_k_eq_4p_empty` / `multilevel_ND_k_ge_4p_proved`)

Official residual-(ii) at k=4p is leftover+splus (Max− dual-bad 3+ even
scores AND min_+=2). Leftover-only with min_+<2 exists (p=5 L8) and is
not residual (ii). Two-value, three-level minus-slice, and 4-level
J-corners are empty.

Remaining estimate, minus-slice interior 4-level (and 5+ / unclassified):

- Paley majorant `E_-[S^2] < 20+12/p` is **false** (p=5 leftover witness;
  15.379 C killed 2026-08-21). Do not reopen.
- Walsh: affine F_2-span of a Max− pair-slice equals the xor-hyperplane
  section of affine_span(Max−). Certified p=3,5,7 rref (15.406 C) and
  p=11 full ensemble (15.596). **15.598 proved:** square-direction
  affine lines cut Max− by \(\sum_S y=0\), so U is the xor cut of a
  named affine-geometry space. 15.602: unique 1-dim G_aff^□-invariant
  subspace of H0 is ⟨1⟩. 15.604: 1_QR ∈ H0 iff p≡1 (mod 4);
  ker(D−I)∩H0 dim 2. 15.605: H0=⟨1⟩⊕W with W=span of extra
  translates, Paley A²=A over F2. 15.606: W=⊕ nsq W^H, M transits;
  G_aff-irred if 2 is a primitive root mod p. 15.607: F_p^× mixes
  Φ_p-factors, so W is G_aff-irreducible for every odd p and
  dir(affine_span(Max−))=H0. 15.608: 1∈dir(U) by antipodes; two
  PSL-orbits of F_p-sublines. 15.609: I(H0)=H0. Walsh spanning of
  V/⟨1⟩ still OPEN. 15.610: Aut({0,∞}) uniqueness DEAD.
  15.611: W ≅ F2[X]/(X^N+1); ker2 dim 2 is a p-law.
  15.612: Walsh ⇔ W1 ∧ W2; CLASS p-law, W1 not.
  15.613: named z∈U; ε(y+Dy) constant on U.
  15.614: W1 for p≡3; named D-spans miss g-orbits at p=11.
  Spanning of the xor-slice (Walsh ∀p) still open. That empties
  interior 4-level only; leftover-only / 5+ remain.

## Leftover 3 (`type_I_multilevel_bad_case_ND_closed`)

On |κ|=1, ν=0 (15.268). Then 3A+B>0 iff G>T. Sufficient:

```
|μ| ≤ |L| = (p-2)/(2p^2)
```

on every |κ|=1 four-set. |μ|≤|T| is not a close (G=T ⇒ 3A+B=0).
|μ|≤maj is false at p=7 and p=11. Census |μ|≤|L| at p=5,7,11 with slack
0.769 → 0.746 → 0.307.

p=5 is already a finite from-C theorem (15.275 L). For every prime
p≥7, the 15.191 target `|μ|≤2/n` is strictly stronger than `|μ|≤L`
(`2/n < L` ⇔ p≥7) and would close Type I. Census: p=7 109/2863 vs
2/50 (~5% slack, tight); p=11 17827/1560713 vs 2/122. Do **not**
prove `|μ|≤|f4|` (false at p=7). Triangle `|R̄₄|+2|φ|` is too weak
for `|μ|≤L` even at p=5. L2 conversion of 3A+B is rejected. Not
imported.

Equivalent size-4 remainder (15.229 F):

```
(p^4 − 1)μ + 2φ = R̄₄,     |R̄₄| ≤ |L|(p^4−1) + 4(p−2).
```

Crude ∑|per(C[S,T])| is ~47500 vs budget ~49 at p=5; unsigned |ν|≤1
on d3 neighbours is likewise too weak. Residual-(i) |μ|≤1/(2p) is
strictly weaker than |μ|≤|L| and does not close Type I.

## 2026-08-21 census on leftover 3 (not a close)

`G_pairmoment_p11.npy` plus Max+ at p=5,7: on |κ|=1, μ is a function of
(κ,φ) only at p=5 (`μ=(4κ−φ)/(p n)`, 4 classes, all constant). At p=7
and p=11 it splits inside (κ,φ) with `star=0` (15.247 B), so a (κ,φ)
linear formula is a p=5 accident. |μ|≤|L| still holds:

| p | max\|μ\| | L | μ/L |
|---|---|---|---|
| 5 | 3/65 | 3/50 | 0.769 |
| 7 | 109/2863 | 5/98 | 0.746 |
| 11 | 17827/1560713 | 9/242 | 0.307 |

JSON: `evidence/mu_kappa_phi_census.json`, `evidence/mu_star_census.json`.
This does not import leftover 3.

## Git (so the next session does not re-read empty branches)

origin has exactly four heads. Local `main` is `origin/main` = merge of
PR #3 (`1dc9bb1`). The other three are ancestors of that tip:

- `codex/leftover-moment-attack` `81ae118` (then this commit)
- `maxplus-p11-enumeration` `5f0ac34`
- `prop15586-maxplus-gram-reduction` `1fa0301`

Packed `.npy` orbits stay on
`/mnt/storage/e1work/maxplus_p13/orbit_attack_2026-08-20/` (not git).
IP backups under `/mnt/storage/cluster-backups/` are NODE.json only.

## What was not done

No leftover flag was imported. No handwritten `return True`. A unit
that does not flip a leftover is not progress on this goal.
