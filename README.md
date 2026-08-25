# Min-max quadratic form of ±1 coefficients

MathOverflow [413935](https://mathoverflow.net/questions/413935) /
[X challenge](https://x.com/PI010101/status/2081070728422752329):

```
m_n = min_{a_ij = ±1}  max_{x_j = ±1}  | Σ_{1≤i<j≤n} a_ij · x_i · x_j |

α_n = m_n / n^(3/2)
```

## About

Machine-assisted attack on a 2022 MathOverflow problem: the limiting constant
of the min-max ±1 quadratic form. The limit **L is OPEN**. This repo is a
fully-audited proof ledger — every claim is a Python predicate that returns
`True`/`False`, ~600 propositions, no prose-only results, and soft-closing is
banned by test (`tests/test_main_chain_docs.py`).

## Status

**Goal:** settle the limit (see **`LONG_HORIZON_GOAL.md`**). Not done until L is proved or disproved.

**Main claim:** L = lim_n α_n is **OPEN** (2026-08-25).

Sandwich and Paley ρ=1 are proved. E(1) on n=p²+1 is **not**. The live
`four_e1_units_closed()` ledger is:

| GOAL unit | live predicate | status |
|---|---|---|
| spectral floor | `phi_F_ge_6` | **OPEN** — needs global QVAR and principal R1 |
| residual (ii), even `k≥4p` | `residual_ii_k_ge_4p` | **OPEN** — Walsh slice closed; non-Walsh multi-level part open |
| Type I, multi-level Max− | `type_I_multilevel` | **OPEN** — `3A+B>0` remains unproved in general |
| Lemma D | `lemma_D` | **TRUE** — construction and two-plane amplitudes checked |

Thus three top-level predicates are false, but the unfinished mathematics is
organized into two fronts: the spectral/QVAR–R1 front and the non-Walsh
multi-level Max− front. Soft-close is forbidden. The acceptance package is
**`evidence/share/denseness_path_package.md`**.

**Proved (sandwich):**
```
1/π  ≤  liminf_n α_n  ≤  limsup_n α_n  ≤  1/2
```

**Also proved:** ρ=1 for Paley conference matrices of order n=p²+1.

See **`STATUS.md`**, `HANDOFF.md`, denseness package, `solution.md`.

---

## Discovery map — what has moved

The problem reduces to **E(1) on Paley conference matrices of order n = p²+1**.
The map below records the current dependency structure. Prop. 15.628 closed
Walsh, W1, and W2; it did not close the non-Walsh multi-level cases. The
principal R1 inequality also remains open, and the current floor wiring
requires the separate global-QVAR estimate:

```mermaid
graph TD
    L["L = lim α_n<br/><b>OPEN</b>"] --> E1["E(1) on n = p²+1"]
    E1 --> D["lemma_D<br/><b>TRUE</b>"]
    E1 --> FLOOR["spectral floor<br/>φ_F ≥ 6<br/><b>OPEN</b>"]
    FLOOR --> QVAR["global mixed-k QVAR<br/><b>OPEN</b>"]
    FLOOR --> R1["principal <b>R1</b><br/>‖δ‖² ≤ n(λ̄−6)²/48<br/><b>OPEN</b>"]
    E1 --> TYPEI["Type I multi-level<br/>3A+B > 0<br/><b>OPEN</b>"]
    R1 -. sufficient .-> TYPEI
    E1 --> RES["residual (ii), even k≥4p<br/><b>OPEN</b>"]
    RES --> WALSH["Walsh / W1 / W2<br/><b>CLOSED</b> (15.628)"]
    RES --> MULTI["non-Walsh multi-level<br/><b>OPEN</b>"]
    style L fill:#ffe6e6
    style D fill:#e6ffe6
    style FLOOR fill:#fff4e6
    style QVAR fill:#fff4e6
    style R1 fill:#fff4e6
    style TYPEI fill:#fff4e6
    style RES fill:#fff4e6
    style WALSH fill:#e6ffe6
    style MULTI fill:#fff4e6
```

The older “two roots, R1 and R2” shorthand now needs two qualifications.
First, the live spectral-floor predicate is `global QVAR ∧ principal R1`, not
R1 alone. Second, only the Walsh component of R2 is closed. A proof of the
strong `n/12` R1 bound would also imply the weaker Type-I `3A+B` estimate,
but no such bound has been proved.

### The R1 collapse (props 15.590–15.597)

A chain of exact identities, each verified as rationals, not numerics:

| step | identity | status |
|---|---|---|
| ν on the ‖κ‖=3 locus | `Σ_S ν(S)² = ½‖m₄⁺‖² − n(n−2)/16` | exact |
| Es4 | `Es4 = 4n² + tr(Φ²)`, Φ = the 15.589 Gram operator | exact |
| design floor | `Es4 ≥ 12n² + 16n + 128n/(n−6)`, equality iff Φ scalar | **proved** |
| particular part | **`Φ_part = λ̄·I`** — the explicit half is spectrally flat | **proved ∀p** |
| residual | `V := ‖Φ − λ̄I‖²_F = 24‖δ‖²` | exact |

so the principal spectral floor and the Type-I sufficient estimate are
bounds on the same scalar `δ`, the master-equation residual tracked since
15.217:

| implication | needs ‖δ‖² ≤ | limit |
|---|---|---|
| principal part of the spectral floor | n(λ̄−6)²/48 | **n/12** ← binding |
| Type-I `3A+B>0` sufficient bound | c₃(p)·n/24 | ~2.9n |
| residual-(i) | `delta_room_for_R` (15.217) | ~n²/8 |

This hierarchy does **not** prove global QVAR, and it does not import any
of the three false GOAL predicates.

### Measured vs. required

`‖δ‖²/n` against the binding threshold ≈ 0.083 — fails at the two census
primes, clears at p=11 with 4.3× margin:

```
p= 5  ██████████████████████████████████████████████  0.9089   (census)
p= 7  ██████████▌                                     0.2085   (census)
p=11  █                                               0.0194   ✓ 4.3× margin
      └─ threshold ≈ 0.083
```

A rigorous **data-free lower bound** on the same scalar, over ten primes, is
flat and converging — no computable quantity threatens the requirement:

| p | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 37 | 47 |
|---|---|---|---|---|---|---|---|---|---|---|
| LB·p⁴/p | 10.00 | 8.91 | 8.34 | 8.24 | 8.14 | 8.11 | 8.08 | 8.05 | 8.03 | **8.02** |

### The Φ spectrum (what leftover 1 actually asks)

`Φ_part = λ̄I` is proved, so **all** spectral deviation comes from δ:

| p | λ_min(Φ) | λ̄ = 8(n−2)/(n−6) | target | margin |
|---|---|---|---|---|
| 5 | 6.1538 | 9.600 | 6 | +0.15 |
| 7 | 7.5110 | 8.727 | 6 | +1.51 |
| 11 | 8.0544 | 8.276 | 6 | +2.05 |

Unconditionally proved: `0 ≤ λ_min(Φ) ≤ λ̄` (lower since Φ is a Gram operator,
upper since `tr Φ_δ = 0`). **The entire open content of leftover 1 is the
window [0, 6)** — no argument short of a genuine δ bound reaches it.

### R2 close (leftover 2 Walsh slice, props 15.598–15.628)

Independent root. Square-direction affine lines cut Max−, so `U` is the
xor-hyperplane of `affine_span(Max−)`; `rank(S) = n/2` is now a **theorem for
every odd prime** (15.600).  Prop. 15.628 proves that edge-eligible
nonsquare GQR circles span the target code and constructs every such circle
as an actual `U`-difference using arbitrary affine halfspaces.  Therefore
**Walsh spanning, W1, and W2 are proved for every odd prime**.  The p=11
37,457,112-point scan remains an independent holdout; the explicit p=19
affine witness supersedes the earlier generic-solver timeout.

### Exact Paley-lattice structure (props 15.629–15.631)

The post-Walsh attack exposed a precise lattice behind R1. Let
`L = ker_Z(C−pI)`, let `P=(I+C/p)/2`, and let `A` be generated by the
square-direction affine-circle words.

| proposition | proved result | boundary |
|---|---|---|
| 15.629 | the profile glue gives `[L:A]=p^((m−1)(m−2)/2)`, `det(L)=2p^(m²)`, `L*=P Z^n`, discriminant `Z/2 ⊕ (Z/p)^(m²)`, and level `4p` | identifies the exact lattice; no R1 bound |
| 15.630 | `min(L*)=1/2`; the complete minimum shell is `{±Pe_i}` with kissing number `2(p²+1)`; every other nonzero dual vector has norm at least `(p−1)/p` | ordinary dual shell, not the odd Max+ coset shell |
| 15.631 | the Max+ coset phase is radial: `<u,y₀> ≡ 2p‖u‖² (mod 2)`; the first transformed degree-four harmonic shell has a positive exact coefficient | higher dual-shell harmonic sums remain uncontrolled |

These are general theorems for every odd prime, including the standard
Paley `(25,50)` adjacent-ETF minimum-shell case. They convert R1 into a
level-`4p` norm-parity-twisted harmonic theta problem with its first shell
known exactly. They are a substantial structural advance, but they do not
prove R1, global QVAR, E(1), or the limit.

### Route kills — do not re-tread

Recorded with counterexamples so they are not reopened:

| killed | why |
|---|---|
| level-4 moment/SDP relaxation | feasible points beat both thresholds (p=5, 7) |
| Delsarte 2-design + min distance | LP min far below the target |
| degree escalation of the contraction kernel | K₄ grows; degree 6 adds nothing at p=7 |
| any `(12+ε)n²` majorant for Es4 | structurally insufficient — 12 is forced |
| uniform `M ≤ C/p⁴` | **falsified** at p=17: true scaling is `M ≳ 8/p³` |
| L² δ-bound for leftover 2 | error/signal ≈ p/11 → ∞, crosses 1 at p=11 |
| linear 4-point and 6-point LPs | feasible-but-negative while true pairing is positive |
| Γ_δ quantization | p=5 integrality was a single-orbit artifact; dies at p=7 |

The older class-function plan
(`evidence/PLAN_2026-08-22_class_function_route.md`) remains a detailed
record of the PSL/Hecke compression and its killed shortcuts. The current
R1 structure is sharper: Props. 15.629–15.631 identify the integral glue,
the exact dual minimum shell, and the radial Poisson phase. The missing step
is now explicitly the signed degree-four contribution of the higher dual
shells, not the first shell or an unidentified glue-class phase.

### What is left

1. **Spectral floor:** prove global mixed-`k` QVAR and the principal R1
   bound `‖δ‖² ≤ n(λ̄−6)²/48` (the simpler `n/12` bound is sufficient).
   The lattice/shadow theorems identify the exact theta object but do not
   control its higher dual shells.
2. **Non-Walsh multi-level Max−:** close residual (ii) for even `k≥4p`.
   Walsh/W1/W2 are done; the 5+-level branch is not. The related Type-I
   `3A+B>0` gate also remains false, although strong R1 would imply it.

Lemma D is complete and is no longer on the work list.

---

## Files

| Path | Role |
|------|------|
| `HANDOFF.md` | Research handoff / resume entry point |
| `evidence/HISTORY_AND_REFERENCES.md` | MO/X/Paata education and pre-internet sources (not a close) |
| `solution.md` | Full mathematical writeup |
| `src/e1_gmin_m4_prop15167.py` … `prop15171.py` | Bi-tight + E(1) residual ND modules |
| `src/e1_gmin_m4_prop15590.py` … `prop15597.py` | R1 collapse: ν → Es4 → Φ → δ; principal/Type-I bound hierarchy |
| `src/e1_gmin_m4_prop15598.py` … `prop15601.py` | R2: square-direction lines, rank(S)=n/2, Walsh |
| `src/e1_gmin_m4_prop15628.py`, `scripts/w2_affine_circle_close.py` | R2 close: eligible GQR circle span + explicit affine completions |
| `src/e1_gmin_m4_prop15629.py` | Profile-glued integral Paley eigenspace lattice |
| `src/e1_gmin_m4_prop15630.py` | Exact dual minimum shell and kissing number |
| `src/e1_gmin_m4_prop15631.py` | Radial dual-shadow transform of the Max+ odd coset |
| `evidence/NOTE_2026-08-24_r1_profile_glue_lattice.md` | Proof note for the lattice quotient, determinant, dual, and level |
| `evidence/NOTE_2026-08-25_dual_minimum_shell.md` | MDS/Newton proof of the exact dual shell |
| `evidence/NOTE_2026-08-25_radial_dual_shadow.md` | Poisson phase, dual gap, and first harmonic shell |
| `evidence/PLAN_2026-08-22_class_function_route.md` | PSL/Hecke route ledger and killed shortcuts |
| `scripts/frame_line_system.py` | Data-free frame-line solver (any p, no Max± ensemble) |
| `src/minmax_quadratic.py` | Exact `m_n`, Paley, Φ, bounds, ρ=1 evec |
| `tests/test_prop15167.py` … `test_prop15171.py` | Load-bearing E(1)/L tests |
| `x-cards/` | X summary + key-lemmas JPEGs |
| `evidence/share/` | Paper PDF/TeX + share assets |
| `evidence/` | Verification JSON and session notes |

## Quick check

```bash
python3 -m pytest tests/test_minmax.py -v
python3 -m pytest tests/test_prop15628.py tests/test_prop15629.py tests/test_prop15630.py tests/test_prop15631.py -q
python3 -c "from src.minmax_quadratic import exact_m; print([exact_m(n) for n in range(2,9)])"
```

## Exact small values

| n | m_n | α_n (approx) |
|---|-----|--------------|
| 2 | 1 | 0.354 |
| 3 | 3 | 0.577 |
| 4 | 4 | 0.500 |
| 5 | 4 | 0.358 |
| 6 | 5 | 0.340 |
| 7 | 9 | 0.486 |
| 8 | 10 | 0.442 |
| 9 | 12 | 0.444 |
| 10 | 13 | 0.411 |

At n=10, Paley (order p²+1, p=3) has Φ=15 > m_10: conference is not exactly optimal.
Exact optima first appear at Hamming distance 5 from Paley, and the only 5-edge undercutters are 144 perfect matchings — see `evidence/N10_STRUCTURE.md`. Those 144 form one PΓL(2,9)-orbit (maximizer-drop criterion) — see `evidence/N10_MATCHING_CLASSIFY.md`.
