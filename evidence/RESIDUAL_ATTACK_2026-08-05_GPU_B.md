# Residual attack log — 2026-08-05 GPU B session

## Target
Prove `E[s⁴] ≤ Es4_*(p)` for all primes p≥5 (Prop 15.164/15.165), or equivalent ∑η²≤η_*.

## Tried (failed to close general p)

1. **Spherical Delsarte LP** (Gegenbauer ≤ deg 12, min-distance τ=(n−2p−2)/n, antipodal 2-design, fixed μ(±1)=1/N): upper bound collapses to two-level angle extremal; at p=5 gives 12097 > Es4_*=9297; at p=7 gives ~1.1e5 ≫ 32336. **Dead** as pure spherical-code bound.

2. **Two-level angle majorization** with N≥|PSL(2,p²)|/30: same order of weakness for p≥5.

3. **C-signature / cross-ratio alone**: m4 not constant on C-edge types; CR classes still split m4 (need full Aut-orbit / field invariants finer than (χ(r),χ(1−r),T)).

4. **Bose–Mesner of conference graph**: saved maximizers Bmax satisfy ⟨B,I⟩=⟨B,C⟩=⟨B,J⟩=⟨B,C²⟩=0 — maximizers **outside** BM(C). Association-scheme eigenvalues of C do not pin Φ top.

5. **Gershgorin on Veronese frame Gram**: λ_max ub ~53 (p=5), ~96 (p=7) ≫ 16.

6. **Q4≤1 / Haar+γ**: recovers only 2n³-scale or worse.

## Shipped progress (not a close)
- Prop **15.165**: exact global Es4 p=3,5,7; closed Es4_*/η_*; GoG↔Φ; m4 C-eigen; p=7 W_CENSUS trap; H-sat p=5.
- Census 16N at p=5,7 via **spectrum Es4** (not single-root W).

## Status
`residual_closed_general=false`. **L OPEN**. No soft-close.

## Preferred remaining attacks
- Weil/Jacobi-sum evaluation of Aut-orbit m4 on PG(1,p²) with N≥|PSL|/|Stab| lower bound.
- Aut₀-isotype embedding V₊↪V_max with closed Rayleigh of the 1-dim line.
- SOS/Putinar on Q₄≤10N‖B‖² uniform in p.

## Round 2 attack (post-15.166)

### Tried
1. **Gegenbauer dual LP** (a0+a2 Q2+a4 Q4+… ≥ t⁴, a_{k≥4}≤0): collapses to deg-2 bound ~2n³ ≫ Es4_*.
2. **CS on ∑ρ κ_B**: sum_rho² ~78–108, budget ~1; needs √(sum ρ² · κ_B²) ≪ 1 — fails by factor ~5–7.
3. **Majorization Q2 with λ_min only**: ub_coeff > thr for all d≥5.
4. **Continuous Veronese operator**: asymptotic λ₂(Ŝ)/N → 2/(d(d+2)) < 4/d² with factor d/(2(d+2))→½ margin — **suggests large-p 16N**, but finite-design defect not controlled.
5. **Wick orth ⇒ delta² from Es4**: computed delta² ≠ certified (Wick orth fails or ρ≠η).
6. **Entrywise |δ|≤3/p²**: n4·9/p⁴ ≫ room_hyp/24.
7. **Equidiagonal projector universal bound**: random equidiag needs N/d² ≳ 1.5–2; Max+ has this for p≥7 but no proof for all equidiag projectors; p=5 has N/d²≈1.54 borderline.

### Status
residual_closed_general=**false**. L **OPEN**. No soft-close.
16N for general p still requires a new idea (Weil/Jacobi closed m4, SOS, or equidiag projector theorem with N≥c d²).

## Round 3 (post-reboot, Prop 15.167)

### Breakthrough (not residual close)
**Bi-tight empty for all primes p≥5** without residual/16N:
- mult(λ_max Φ)≥d−1 (15.162/15.98) + λ_min(Φ)≥6 (15.160.D)
- majorization ⇒ λ_max ≤ L_*=(p⁴+24p²−1)/(2(p²−1))
- L_* < 2d ⇔ p⁴−24p²−1>0 for p≥5
- λ_cycle=λ_max(Φ)/2 < d ⇒ λ_max(G)=n/2 simple ⇒ bi-tight empty (15.55)

Shipped: `src/e1_gmin_m4_prop15167.py`, `tests/test_prop15167.py`, `src/e1_bitight_chain.py`.

### Residual status (unchanged honest)
`residual_closed_general=false`. 16N / Es4_* / H still OPEN. Not required for bi-tight.
Dead residual attacks not re-run.

### L status
E(1)/deep ND + Main still OPEN ⇒ **L OPEN** (F3 no soft-close).

## Round 4 (Prop 15.167–168, post-skeptic honest)

### Real close
- **15.167 bi-tight empty all primes p≥5** via mult≥d−1 + λ_min≥6 majorization (Fraction). No residual/16N.

### 15.168 partial E(1) structure (checkable predicates)
- Tight level-s obstruction when bi-tight empty (G_⊥ identity + 15.167)
- Type I freeness ND (prior 15.43.1)
- Type I fail k=2p−1 ND when bi-tight empty (15.43.3+15.44+15.167)
- Deep auto-freeness for s₊=2, k≤3p−2 (Fraction)
- Deep fail-eq k=3p−1 ⇒ tight L3 empty under 15.167

### E(1) residual (honest OPEN — no soft-close)
1. Type I freeness-fail at k=3p−2 / S∈{1,5} boundary
2. Deep non-tight freeness-fail with k≥3p

E1_closed=false. L OPEN. residual/16N false.

## Round 5 (Prop 15.169 — Type I k=3p−2 reduction)

### Shipped (real Fraction / prior props)
- Type I freeness-fail k=3p−2 structure (a=thr, affine S+2f_e=3, H scores {2,4})
- Φ is 2-Lipschitz under edge flip ⇒ ND free when Φ(G)≥Φ
- Gap-2 Type I undercutter with k odd forces s_−=−1
- Deep multi-s auto-freeness: k≤p(s+1)−2 for min-level s≥2
- Reduction of residual (i) to: prove s_−≤−1 impossible under freeness-fail affine (or bad case f_e≡−1 on U_−)

### Certified (not general-p)
- p=5 MILP: freeness-fail affine + S≤−1 on all Max− is **infeasible** (HiGHS) for all tested e

### Status
residual (i) **OPEN** for general p (reduction shipped; final step open).
residual (ii) **OPEN** (multi-s auto-freeness only).
E1_closed=false. residual_closed_general=false. L **OPEN**. No soft-close.

## Round 6 (Prop 15.170 — residual (i) Gsum Farkas close)

### Shipped (real Fraction)
- Dual equality forces (Gsum x)[e] = 6/p−4
- Gsum wedge = 0; Gsum_ab ≥ −12/(p(p²+1)); box-sum LB = −12k/(p n)
- need < LB for all primes p≥5 ⇔ 4p³−6p²−32p+18 > 0
- type_I_k_3p_minus_2_closed_general = True from this predicate
- ES2 < k integrality seed

### Status
residual (i) **CLOSED** for all p≥5.
residual (ii) **OPEN** (deep freeness-fail k≥3p).
E1_closed=false. residual_closed_general=false. L **OPEN**. No soft-close.
Module: `src/e1_gmin_m4_prop15170.py`, tests `tests/test_prop15170.py`.

## Round 7 (Prop 15.171 — residual (ii) deep freeness-fail ND)

### Shipped (real Fraction)
- Parity: s₊=2 ⇒ k even ⇒ s_−≠−1
- Gap-2 deep ⇒ s_−≤−2
- Deep freeness ⇒ weak ND (Q_H=Φ−2)
- Auto k≤3p−2; fail-eq k=3p−1 empty under bi-tight
- Dual two-level freeness-fail k≥3p: (Gsum x)[e]=2(8−3k/p) < −12k/(pn) box-sum LB
- deep_s2_freeness_fail_k_ge_3p_ND_closed=True from predicates
- E1_closed = type_I ∧ deep_k≥3p ∧ bi-tight; L = E1 ∧ bi-tight

### Status
residual (i) **CLOSED** (15.170).
residual (ii) **CLOSED** (15.171).
E1_closed=true. L **CLOSED**. residual_closed_general=false (16N optional). No soft-close.
Module: `src/e1_gmin_m4_prop15171.py`, tests `tests/test_prop15171.py`.
