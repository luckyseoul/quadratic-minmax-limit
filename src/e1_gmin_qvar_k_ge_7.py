#!/usr/bin/env python3
"""
QVAR on k≥7 — ensemble floor; pointwise/orbitwise already false.

Does **not** number a leftover-1 hinge.  Does **not** close leftover 1
(principal δ-room still open).  p=13 orbits are not a general close.

FLOORS (Max+-free Fraction identities)
  A. λ_exc = 32 E|Z_ψ|²/[q(q-1)]  (15.589 E).  QVAR
        E|Z_ψ|² ≥ 3q(q-1)/16  ⇔  λ_exc ≥ 6.
     Fail: drop 16; replace 32 by 16.
  B. For p≡3 (mod 4): a_L ∈ 2pℤ, b_L=a_L/(2p), T=(p²−1)/8,
        E|∑ ψ(g_L) b_L|² ≥ 3T/8
     is the same inequality (multiply by 4p²).  Fail: T↦(p²−1)/4; drop 8.
  C. Weil/Fourier (15.589 M): the k-stratum is empty whenever k≥4 and
        p > 4k².
     Hence for every k≥7 and every prime p>4k², QVAR holds vacuously on
     that (k,p).  Fail: 4k² ↦ k².  This does **not** cover the top
     stratum k=m=(p+1)/2 (never p>4m²).

COUNTEREXAMPLES (not a close)
  D. Pointwise QVAR is false at the first live pair: p=13 k=7 CP-SAT
     witness |Z_ψ|² = 2548 < 10647/2 = 3q(q-1)/16.  Fail: claim 2548
     meets the floor.
  E. A single signed-PSL orbit can still average above the floor
     (806468/85 at that p=13 k=7 lift).  That is not a p-law and is
     **not imported**.

SIEVE LINEAR ALGEBRA (proved; not a QVAR close)
  G. Kernel ladder.  For k distinct points of P¹(F_p) the degree-d
     homogeneous evaluation matrix is (d+1)×k of rank min(d+1,k).
     The leading-coefficient kernel in F_p^k therefore has dimension
        k−d−1     whenever k≥d+1.
     Fail: drop the −1 (claim dim k−d).  At top degree d=k−2 the kernel
     is 1-dimensional with full support: a vanishing coordinate would
     leave a square Vandermonde of k−1 distinct points, hence only the
     zero vector.  Translation spans the 2-dimensional degree-(k−3)
     kernel.  This is the recursive coupled-sieve normal form, not QVAR.

  H. Zero top-scalar energy.  If that 1-dimensional top scalar vanishes,
     every active profile has reduced degree ≤k−3.  Weil’s bound at
     r=k−3 empties the k-stratum once
        (k−3)² p > 4k²(k−4)²
     (for k=7: p≥113).  Fail: use r=k−2 (the 15.589 M cutoff p>4k²).
     Nonzero top scalar, and every top stratum k=m, remain open.

OPEN
  F. Ensemble QVAR on every k≥7 for every prime p≥13 with p≤4k²
     (includes the top stratum for all p≥13, and nonzero-λ classes
     below the zero-scalar cutoff).  After V–AB, F is E|Q|^2 on V_{p,k}
     with Q=Q_<r+2λ C+λ²α (p≡3 top: α=0).  Equivalently E[⟨y,F⟩²]≥3pS.
     k=7 is energy-empty at tabled p≥53 (AE); live k=7 is p=13..47 plus
     any later p where 7 b_min≤T.  p=41 k=7 is nonempty with E|Z|²=0
     (AF), so per-stratum QVAR is false at that pair.  qvar stays False.

ATTACK GRAPH (graph-engineered-completion; not a close)
  FloorA --iff--> FloorB
  WeilC --vacuous_QVAR--> (k,p) with p>4k²
  WeilC --does_not_cover--> top k=m
  PointwiseD --refutes--> pointwise/single-vector QVAR
  KernelG --enables--> coupled sieve; --does_not_prove--> F
  ZeroScalarH --empties--> λ=0 classes in a p-range; --does_not_prove--> F
  TopExactSpherical --killed--> p=11 k=6 moment ≠ V_sph
  Bochner --too_weak--> E|Z|²≥0
  UniqueHarm4 --killed--> GAP dim Harm_4^G ∈ {2,3,6} at p=5,7,11
  ApsiRank1 --killed--> tr A_ψ=0 and ||A||_HS>0 ⇒ rank≥2
  UniformLambda --killed--> p=11 full-support masses (15.589 K)
  CSHarm4 --killed--> p=5 |E−V_sph| > gap so ||μ||||f|| ≰ gap
  PairingId --global_only--> E=V_sph+⟨μ4,f4⟩ ⇔ ⟨μ4,f4⟩≥−gap; k-stratum is not a 2-design
  LambdaReduce --partial--> λ=0 empty ⇒ need only λ≠0; top k=m never Weil-empties λ=0
  TopVarP --iff_not_a_bound--> top p=3 mod 4: E|Z|²=4 E[(T₊−S/2)²]; Var bound OPEN
  EqualEnergyQ --not_a_lower_bound--> E_eq≥floor for k≤m−1 but actual can be < E_eq (p=11 k=4)
  TwoTypeR --killed--> p=11 opposite-pair E[a_i a_j] splits (Walsh-style uniqueness fail)
  LambdaNonzeroS --killed--> p=13 k=7 witness has top_scalar=7≠0 and |Z|²=2548<floor
  TwoDesignT --killed--> E[Q_L²] has a 4-distinct remainder; 2-design collision 2p(p−1)(3p−4) is not E[Q²]
  SecondMomentU --reduces_not_proves--> E|Z|²=wᵀKw on k-stratum; QVAR iff Rayleigh of K on ψ; λ_min bound OPEN
  SingerV --iff_not_a_bound--> C_m regular ⇒ K circulant every k; p≡3 mod 4 sign multiplicity one ⇒ w eigenvector; QVAR iff Nyquist ≥ floor/m. Bound OPEN. p≡1 mod 4: no real sign (m odd)
  TopKernelW --reduces_not_proves--> Singer-gauge top kernel is (−1)^k (p≡3 mod 4: only rational characters 1 and sign); leading energy sign-blind. Bound OPEN
  TopLeadX --reduces_not_proves--> 15.588 C + W: top F_p-leading ∈ span{w}; p=7 all saturate; p=11 has λ=0. Remainder unbound
  PhiModeY --iff_not_a_bound--> E|Z|²=E|∑_Ω ψ ρ|²/(16q); QVAR iff ψ-mode of Φ ≥ 3q²(q-1); Φ constant ⇒ EZ=0 false. Bound OPEN
  Phi3levelZ --killed--> k-stratum CS Φ(1)≥4q² m/k; 3-level S_□ is not a lower bound; Wick 8q² is global-only (p=7 top Φ(1)=266168/15<19208)
  P41K7AF --stratum_qvar_false--> unique 7×30=T; 6720 Boolean; E|Z_ψ|²=0 < 529515; not a k≥7 close
  LineDftAA --iff_not_a_bound--> ŷ(λα_j)=σ̂_j(λ); ξ_j=∑_{ℓ_j}ψ ρ; η=ψ|F_p^* is 1 (p≡3) or Legendre (p≡1); K_ξ circulant, QVAR is the trivial mode. Bound OPEN
  JacobiGramAB --kills_unrestricted--> W=c^*Jc on monomials; η=1 J=pμ_{d+e}−μ_dμ_e PSD (discrete Chebyshev); η=χ indefinite at p=13 deg≤5 so unrestricted Weil is not a certificate; affine restriction recovers L(2,χ); λ=0 is ker J_{00}. Coupled LDL on V_{p,k} OPEN
  CoupledAC --leading_pivot--> Q=Q_<r + 2λ C(c_low) + λ² α on V_{p,k}; ψ(η)=i^{p-1}; p≡3 top α=0 (CD last pivot vanishes); p≡1 ψ constant on Singer lines α=m ψ_0 J_rr≠0; drop 2λ killed. Bound OPEN
  RidgeAD --iff_not_a_bound--> p≡3: Z=⟨y,F⟩/2, ||F||²=pS (off-diag fibers vanish); QVAR iff E[⟨y,F⟩²]≥3pS. Pythagoras E[Q²]=E[Q_low²]+4E[λ²C²]+4E[λ Q_low C]. Actual E[Z²] not ∝λ². Bound OPEN
  EnergyAE --vacuous_k7--> 7 b_min>T empties k=7 at p=53,59,61,67,71 (exact quintic lifts). Live p=13..47. Not a top / k≥8 close
  F --required_by--> qvar_k_ge_7_proved_general
  qvar ∧ principal --required_by--> leftover1 / phi_F
  leftover1 ∧ leftover2/3 ∧ lemmaD --required_by--> L
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15589 import (  # noqa: E402
    lambda_exc_from_quartic_variance,
    normalized_profile_energy_total,
    normalized_quartic_variance_threshold,
    q_of,
    quartic_variance_floor_threshold,
    spherical_quartic_variance,
    theorem_E_exceptional_quartic_variance,
    theorem_U_k6_QVAR_all_primes,
    weil_activity_barrier_excludes,
)

# p=13 k=7 CP-SAT witness (evidence/k7_p13_cpsat_witness.json).
P13_K7_POINTWISE_ABS_ZPSI_SQ = 2548
P13_K7_ORBIT_MEAN_ABS_ZPSI_SQ = Fraction(806468, 85)


def T_of(p: int) -> int:
    """T=(p²−1)/8.  Integral for odd p."""
    return (p * p - 1) // 8


def T_wrong_drop8(p: int) -> int:
    """Fail-eq: /8 ↦ /4."""
    return (p * p - 1) // 4


def integer_qvar_threshold(p: int) -> Fraction:
    """3T/8."""
    return Fraction(3 * T_of(p), 8)


def integer_qvar_threshold_wrong_drop8(p: int) -> Fraction:
    return Fraction(3 * T_of(p), 4)


def qvar_floors_equivalent(p: int) -> bool:
    """4p² · (3T/8) = 3q(q−1)/16."""
    q = q_of(p)
    ez = quartic_variance_floor_threshold(p)
    return Fraction(4 * p * p) * integer_qvar_threshold(p) == ez == Fraction(
        3 * q * (q - 1), 16
    )


def p13_k7_qvar_threshold() -> Fraction:
    return quartic_variance_floor_threshold(13)


def pointwise_qvar_false_p13_k7() -> bool:
    """The CP-SAT witness is strictly below the ensemble floor."""
    return P13_K7_POINTWISE_ABS_ZPSI_SQ < p13_k7_qvar_threshold()


def pointwise_qvar_wrong_claim_meets_floor() -> bool:
    """Fail-eq: claim 2548 ≥ 10647/2."""
    return P13_K7_POINTWISE_ABS_ZPSI_SQ >= p13_k7_qvar_threshold()


def single_orbit_mean_clears_but_not_a_close() -> bool:
    """Orbit average > floor, not imported as a p-law."""
    return P13_K7_ORBIT_MEAN_ABS_ZPSI_SQ >= p13_k7_qvar_threshold()


def lambda_exc_from_pointwise_p13_k7() -> Fraction:
    return lambda_exc_from_quartic_variance(
        13, Fraction(P13_K7_POINTWISE_ABS_ZPSI_SQ)
    )


def weil_vacuous_qvar_k_ge_7(p: int, k: int) -> bool:
    """True iff the (k,p) stratum is empty by 15.589 M, hence QVAR holds
    vacuously there.  Not a general k≥7 close (top stratum never empty)."""
    if k < 7:
        return False
    return weil_activity_barrier_excludes(p, k)


def weil_barrier_wrong_drop4(p: int, k: int) -> bool:
    """Fail-eq: p>4k² ↦ p>k²."""
    return k >= 4 and p > k * k


def kernel_dim_leading(k: int, degree: int) -> int:
    """dim ker = k−d−1 for k distinct P¹ points, k≥d+1."""
    return k - degree - 1


def kernel_dim_wrong_drop1(k: int, degree: int) -> int:
    """Fail-eq: k−d−1 ↦ k−d."""
    return k - degree


def _pow_mod(base: int, exp: int, p: int) -> int:
    if exp == 0:
        return 1
    return pow(base % p, exp, p)


def homogeneous_eval_matrix(
    points: list[tuple[int, int]], degree: int, p: int
) -> list[list[int]]:
    """(d+1)×k matrix, entry (e,j) = x_j^e y_j^{d−e}."""
    return [
        [_pow_mod(x, e, p) * _pow_mod(y, degree - e, p) % p for x, y in points]
        for e in range(degree + 1)
    ]


def matrix_rank_modp(matrix: list[list[int]], p: int) -> int:
    """Row-rank over F_p."""
    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    n_r, n_c = len(rows), len(rows[0])
    rank = 0
    col = 0
    for r in range(n_r):
        pivot = None
        while col < n_c:
            for i in range(r, n_r):
                if rows[i][col] % p:
                    pivot = i
                    break
            if pivot is not None:
                break
            col += 1
        if col >= n_c:
            break
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = pow(rows[r][col] % p, p - 2, p)
        rows[r] = [(v * inv) % p for v in rows[r]]
        for i in range(n_r):
            if i == r:
                continue
            factor = rows[i][col] % p
            if factor:
                rows[i] = [
                    (rows[i][c] - factor * rows[r][c]) % p for c in range(n_c)
                ]
        rank += 1
        col += 1
    return rank


def distinct_p1_points(k: int, p: int) -> list[tuple[int, int]]:
    """k distinct points of P¹(F_p): (1,0),… then (0,1)."""
    if k > p + 1:
        raise ValueError("P¹(F_p) has only p+1 points")
    pts = [(1, j) for j in range(min(k, p))]
    if k == p + 1:
        pts.append((0, 1))
    return pts


def kernel_dim_of_points(
    points: list[tuple[int, int]], degree: int, p: int
) -> int:
    k = len(points)
    rank = matrix_rank_modp(homogeneous_eval_matrix(points, degree, p), p)
    return k - rank


def top_kernel_full_support(k: int, p: int) -> bool:
    """The unique top kernel vector (d=k−2) has no zero coordinate."""
    d = k - 2
    points = distinct_p1_points(k, p)
    M = homogeneous_eval_matrix(points, d, p)
    # M is (k-1)×k of rank k-1; kernel spanned by the signed maximal minors.
    # A zero coordinate ⇔ remaining (k-1)×(k-1) block is singular.
    for drop in range(k):
        minor = [row[:drop] + row[drop + 1 :] for row in M]
        if matrix_rank_modp(minor, p) < k - 1:
            return False
    return True


def zero_top_scalar_weil_excludes(p: int, k: int) -> bool:
    """λ=0 (degree ≤k−3) empty by Weil.  Integer form of
    (k−3)² p > 4 k² (k−4)².  Requires k≥7 so k−4>0."""
    if k < 7 or p < 13:
        return False
    return (k - 3) * (k - 3) * p > 4 * k * k * (k - 4) * (k - 4)


def theorem_A_floors(primes=(5, 7, 11, 13, 17, 19, 23, 31, 43, 47)) -> dict:
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        ez = quartic_variance_floor_threshold(p)
        lam6 = lambda_exc_from_quartic_variance(p, ez)
        row_ok = (
            qvar_floors_equivalent(p)
            and lam6 == 6
            and T_wrong_drop8(p) != T_of(p)
            and integer_qvar_threshold_wrong_drop8(p) != integer_qvar_threshold(p)
        )
        if p % 4 == 3:
            row_ok = (
                row_ok
                and T_of(p) == normalized_profile_energy_total(p)
                and integer_qvar_threshold(p)
                == normalized_quartic_variance_threshold(p)
            )
        rows[str(p)] = {
            "T": T_of(p),
            "E_Z_floor": str(ez),
            "integer_3T_8": str(integer_qvar_threshold(p)),
            "equivalent": qvar_floors_equivalent(p),
            "lambda_exc_at_floor": str(lam6),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "theorem": (
            "E|Z_ψ|²≥3q(q-1)/16 ⇔ λ_exc≥6 ⇔ (p≡3 mod 4) E|B|²≥3T/8 "
            "with T=(p²-1)/8.  Fail: drop 16; T/8↦T/4."
        ),
        "by_p": rows,
    }


def theorem_D_pointwise_orbit_counterexample() -> dict:
    thr = p13_k7_qvar_threshold()
    lam = lambda_exc_from_pointwise_p13_k7()
    return {
        "proved": bool(
            pointwise_qvar_false_p13_k7()
            and not pointwise_qvar_wrong_claim_meets_floor()
            and lam < 6
            and thr == Fraction(10647, 2)
            and T_of(13) == 21
            and integer_qvar_threshold(13) == Fraction(63, 8)
            and single_orbit_mean_clears_but_not_a_close()
            and P13_K7_ORBIT_MEAN_ABS_ZPSI_SQ != Fraction(
                P13_K7_POINTWISE_ABS_ZPSI_SQ
            )
        ),
        "pointwise_abs_Zpsi_sq": P13_K7_POINTWISE_ABS_ZPSI_SQ,
        "floor": str(thr),
        "pointwise_below_floor": pointwise_qvar_false_p13_k7(),
        "lambda_exc_pointwise": str(lam),
        "orbit_mean": str(P13_K7_ORBIT_MEAN_ABS_ZPSI_SQ),
        "orbit_mean_clears": single_orbit_mean_clears_but_not_a_close(),
        "orbit_mean_imported": False,
        "theorem": (
            "p=13 k=7 CP-SAT witness |Z_ψ|²=2548<10647/2, so pointwise "
            "and single-vector QVAR are false.  A free signed-PSL orbit "
            "averages to 806468/85>floor; that is not imported."
        ),
    }


def theorem_G_kernel_ladder(
    samples=((13, 7), (17, 7), (17, 8), (19, 7), (19, 9)),
) -> dict:
    """Vandermonde kernel ladder.  Enables the sieve; does not prove QVAR."""
    ok = True
    rows = {}
    for p, k in samples:
        if k > p + 1 or k < 7:
            continue
        pts = distinct_p1_points(k, p)
        row_ok = True
        dims = {}
        for d in range(1, k - 1):
            got = kernel_dim_of_points(pts, d, p)
            expect = kernel_dim_leading(k, d)
            dims[str(d)] = got
            if got != expect:
                row_ok = False
            if kernel_dim_wrong_drop1(k, d) == got:
                row_ok = False
        top_ok = top_kernel_full_support(k, p)
        row_ok = row_ok and top_ok and kernel_dim_leading(k, k - 2) == 1
        rows[f"p{p}_k{k}"] = {
            "kernel_dims": dims,
            "top_full_support": top_ok,
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "theorem": (
            "For k distinct P¹(F_p) points, deg-d evaluation has rank "
            "min(d+1,k), so leading-coefficient ker dim = k−d−1 (k≥d+1). "
            "Fail: k−d.  Top d=k−2 is 1-dimensional with full support. "
            "Sieve normal form, not QVAR."
        ),
        "samples": rows,
    }


def theorem_H_zero_top_scalar_weil() -> dict:
    """λ=0 classes empty in a p-range; nonzero λ and top stratum remain."""
    # k=7: p≥113 by (4)^2 p > 4*49*9 = 1764 ⇒ p>110.25
    k7_cutoff = zero_top_scalar_weil_excludes(113, 7)
    k7_below = not zero_top_scalar_weil_excludes(109, 7)
    # fail-eq: full Weil p>4k² is strictly weaker (196 vs 113)
    fail_ok = (
        k7_cutoff
        and k7_below
        and not weil_vacuous_qvar_k_ge_7(113, 7)
        and weil_vacuous_qvar_k_ge_7(197, 7)
    )
    top_never = all(
        not zero_top_scalar_weil_excludes(p, (p + 1) // 2)
        for p in (13, 17, 19, 23, 47, 101, 113, 197, 1009)
        if (p + 1) // 2 >= 7
    )
    return {
        "proved": bool(fail_ok and top_never),
        "covers_general_k_ge_7": False,
        "k7_empty_lambda0_from": 113,
        "top_stratum_never_excluded": top_never,
        "theorem": (
            "If the top kernel scalar vanishes, r≤k−3 and Weil empties "
            "the class once (k−3)²p > 4k²(k−4)² (k=7: p≥113).  Fail: "
            "the r=k−2 cutoff p>4k² (p=113,k=7 is λ=0-empty but not "
            "Weil-vacuous).  Does not cover nonzero λ or k=m."
        ),
    }


def theorem_I_top_stratum_not_exactly_spherical() -> dict:
    """p=11 k=6 is the top stratum m=(p+1)/2. Its QVAR moment (15.589 U)
    is strictly above V_sph, so 'top stratum = spherical 4-design value'
    is false.  Not a k≥7 close; kills an exactness shortcut."""
    U = theorem_U_k6_QVAR_all_primes()
    p = 11
    m = (p + 1) // 2
    eb2 = Fraction(U["p11_complete_census_E_B2"])
    ez = Fraction(4 * p * p) * eb2
    vs = spherical_quartic_variance(p)
    thr = quartic_variance_floor_threshold(p)
    ok = (
        m == 6
        and bool(U["proved"])
        and ez >= thr
        and ez != vs
        and vs > thr
        and eb2 != 0
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p": p,
        "k": m,
        "E_abs_Zpsi_sq": str(ez),
        "V_sph": str(vs),
        "QVAR_threshold": str(thr),
        "equals_spherical": ez == vs,
        "theorem": (
            "Top stratum k=m at p=11 (k=6, 15.589 U) has "
            "E|Z_ψ|²=4p²·114771/14903 ≠ V_sph, so the top stratum is "
            "not exactly the spherical 4-design value.  Fail: claim "
            "equality with V_sph.  Does not prove k≥7 QVAR."
        ),
    }


def theorem_J_bochner_too_weak_for_qvar(
    primes=(5, 7, 11, 13, 17, 19, 23),
) -> dict:
    """Autocorrelation of |D̂|² is positive-definite on F_q^*, so the
    ψ-Fourier coefficient (hence E|Z_ψ|²) is ≥0.  The QVAR floor is
    strictly positive, so Bochner does not prove QVAR."""
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        thr = quartic_variance_floor_threshold(p)
        row_ok = thr > 0 and Fraction(0) < thr
        rows[str(p)] = {"floor": str(thr), "floor_positive": thr > 0, "ok": row_ok}
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "bochner_gives": "E|Z_ψ|² ≥ 0",
        "meets_qvar_floor": False,
        "theorem": (
            "u=|1̂_D|², Q(r)=E[u(ξ)u(rξ)] is an averaged autocorrelation, "
            "hence positive-definite on F_q^*.  Its ψ-mode is a nonnegative "
            "multiple of E|Z_ψ|², so E|Z_ψ|²≥0.  Fail: claim 0 meets "
            "3q(q-1)/16 (floor is positive for every prime p≥5).  "
            "Does not prove k≥7 QVAR."
        ),
        "by_p": rows,
    }


# GAP CharacterTable("L2(p^2)") + Sym^4 character (independent audit).
# dim Harm_4(W_e)^G = dim(Sym^4 W_e)^G − dim(Sym^2 W_e)^G, and dim(Sym^2)^G=1.
GAP_HARM4_G = {
    5: 2,   # L2(25), |G|=7800
    7: 3,   # L2(49), |G|=58800
    11: 6,  # L2(121), |G|=885720
}


def theorem_K_harm4_not_one_dimensional() -> dict:
    """The G-invariant 4-harmonics on W_e are not 1-dimensional, so
    E|Z_ψ|² − V_sph is not a multiple of the Es4 / 4-design defect.
    Fail: claim dim Harm_4^G = 1.  Not a QVAR close."""
    ok = True
    rows = {}
    for p, dim in GAP_HARM4_G.items():
        row_ok = dim > 1 and dim != 1
        rows[str(p)] = {"Harm4_G": dim, "is_one": dim == 1, "ok": row_ok}
        ok = ok and row_ok
    # fail-eq visible: unique-harmonic would force all three equal to 1
    unique_killed = not (GAP_HARM4_G[5] == GAP_HARM4_G[7] == GAP_HARM4_G[11] == 1)
    return {
        "proved": bool(ok and unique_killed),
        "covers_general_k_ge_7": False,
        "unique_invariant_4_harmonic": False,
        "theorem": (
            "GAP L2(p²) Weil degree (q+1)/2: dim Harm_4(W_e)^G = 2,3,6 "
            "at p=5,7,11 (Sym^4^G − 1).  Fail: claim dimension 1.  "
            "Hence E|Z_ψ|²−V_sph is not determined by Es4−ED4, so leftover 3 "
            "does not imply QVAR.  Does not prove k≥7 QVAR."
        ),
        "by_p": rows,
    }


def theorem_L_A_psi_not_rank_one() -> dict:
    """tr A_ψ=0 and ||A_ψ||_HS²=q(q-1)/32>0, so A_ψ cannot be rank-1
    (a real rank-1 matrix has trace equal to its sole eigenvalue).
    Fail: claim rank 1.  QVAR is not a single halfspace fourth moment."""
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        if not is_prime(p):
            continue
        q = q_of(p)
        hs = Fraction(q * (q - 1), 32)
        row_ok = hs > 0
        rows[str(p)] = {"HS_sq": str(hs), "ok": row_ok}
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "rank_ge_2": True,
        "theorem": (
            "tr A_ψ=0 (15.589 H) and ||A_ψ||_HS²=q(q-1)/32>0, hence "
            "rank(A_ψ)≥2 over R.  Fail: claim rank 1.  Therefore "
            "yᵀA_ψ y is not c(y·v)² and QVAR is not E[(y·v)⁴].  "
            "Does not prove k≥7 QVAR."
        ),
        "by_p": rows,
    }


def theorem_N_lambda_mass_not_uniform() -> dict:
    """On a 1-dim top kernel, λ∈F_p is not equally massive: at p=11
    full support (15.589 K) the zero class is strictly lighter than
    each nonzero class.  Fail: claim mass(λ=0)=1/p.  Mixture identity
    E=p0 E0+p≠ E≠ is tautological and does not prove QVAR."""
    from e1_gmin_m4_prop15589 import theorem_K_full_support_top_degree_mixing

    K = theorem_K_full_support_top_degree_mixing()
    p11 = K["p11"]
    n0 = p11["top_zero_count"]
    n1 = p11["vectors_per_nonzero_class"]
    n_nz = p11["n_nonzero_scalar_classes"]
    n_tot = p11["full_support_count"]
    p0 = Fraction(n0, n_tot)
    uniform = Fraction(1, 11)
    ok = (
        bool(K["proved_counterexample"])
        and n0 + n_nz * n1 == n_tot
        and p0 != uniform
        and n0 < n1
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p11_p0": str(p0),
        "uniform_1_over_p": str(uniform),
        "uniform_lambda_mass": False,
        "theorem": (
            "λ-mixture E_k=p0 E_{λ=0}+p≠ E_{λ≠0} is an identity.  "
            "Uniform λ∈F_p would give p0=1/p, false at p=11 full support "
            "(15.589 K: 2090880 zero vs 3397438 per nonzero class).  "
            "Fail: claim p0=1/11.  Does not prove k≥7 QVAR."
        ),
    }


def theorem_O_cs_harm4_pairing_too_weak() -> dict:
    """Cauchy–Schwarz on Harm_4^G cannot prove QVAR: at p=5 the pairing
    excess |E|Z_ψ|² − V_sph| already exceeds the spherical gap, so
    ||μ_4|| ||f_4|| > gap and excess ≥ −||μ||||f|| is weaker than
    excess ≥ −gap.  Fail: claim |excess| ≤ gap at p=5.  Census from
    15.589 E; not a k≥7 close."""
    from e1_gmin_m4_prop15589 import (
        spherical_QVAR_gap,
        theorem_E_exceptional_quartic_variance,
    )

    E = theorem_E_exceptional_quartic_variance()
    p = 5
    ez = Fraction(E["by_p"][str(p)]["E_abs_Zpsi_sq"])
    vs = spherical_quartic_variance(p)
    thr = quartic_variance_floor_threshold(p)
    gap = spherical_QVAR_gap(p)
    excess = ez - vs
    ok = (
        bool(E["proved_census"])
        and gap > 0
        and abs(excess) > gap
        and ez >= thr
        and vs > thr
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p": p,
        "excess": str(excess),
        "gap": str(gap),
        "cs_hypothesis_excess_le_gap": abs(excess) <= gap,
        "theorem": (
            "On the 2-design, E|Z_ψ|² = V_sph + ⟨μ_4,f_4⟩.  CS gives "
            "⟨μ_4,f_4⟩ ≥ −||μ_4|| ||f_4||, and ||μ||||f|| ≥ |excess|.  "
            "At p=5, |3300/13 − V_sph| > q(q-1)(q-11)/16(q+5), so the "
            "CS hypothesis ||μ||||f|| ≤ gap is false.  Fail: claim "
            "|excess|≤gap.  Does not prove k≥7 QVAR."
        ),
    }


def theorem_C_weil_vacuous_range() -> dict:
    """Vacuous QVAR on k≥7 when p>4k².  Top stratum never included."""
    ok = True
    rows = {}
    samples = ((7, 13, False), (7, 197, True), (8, 257, True), (7, 193, False))
    for k, p, expect in samples:
        got = weil_vacuous_qvar_k_ge_7(p, k)
        row_ok = got is expect
        # fail-eq visible: at (k,p)=(7,197), p>k² already (197>49), so
        # drop-4 still True; use (7,50): 50>49=k² but 50<196=4k².
        rows[f"k{k}_p{p}"] = {
            "vacuous": got,
            "expect": expect,
            "ok": row_ok,
        }
        ok = ok and row_ok
    # dedicated fail-eq pair
    fail_ok = (
        not weil_vacuous_qvar_k_ge_7(50, 7)
        and weil_barrier_wrong_drop4(50, 7)
    )
    m = lambda p: (p + 1) // 2
    top_never = all(
        not weil_vacuous_qvar_k_ge_7(p, m(p))
        for p in (13, 17, 19, 23, 47, 101, 1009)
        if m(p) >= 7
    )
    return {
        "proved": bool(ok and fail_ok and top_never),
        "covers_general_k_ge_7": False,
        "top_stratum_never_vacuous": top_never,
        "theorem": (
            "15.589 M: k-stratum empty for k≥4 and p>4k², hence QVAR "
            "vacuous on those (k,p) for k≥7.  Fail: 4k²↦k² (p=50,k=7).  "
            "Does not cover k=m=(p+1)/2."
        ),
        "samples": rows,
    }


def profile_energy_total_S(p: int) -> int:
    """S=p(p²−1)/4.  Pointwise sum of a_L for p=3 mod 4."""
    return p * (p * p - 1) // 4


def equal_energy_prediction(p: int, k: int) -> Fraction:
    """S² (m−k)/(k(m−1)) for sum w=0 and equal a_L=S/k on a k-subset."""
    if p % 4 != 3:
        raise ValueError("p=3 mod 4")
    m = (p + 1) // 2
    S = profile_energy_total_S(p)
    if k < 1 or k > m:
        raise ValueError("k in 1..m")
    if k == m:
        return Fraction(0)
    return Fraction(S * S * (m - k), k * (m - 1))


def equal_energy_meets_floor_k_le_m_minus_1(p: int, k: int) -> bool:
    """Algebra: E_eq ≥ QVAR floor for every 1≤k≤m−1, p=3 mod 4, p≥7."""
    return equal_energy_prediction(p, k) >= quartic_variance_floor_threshold(p)


def top_plus_energy_qvar_threshold(p: int) -> Fraction:
    """QVAR on the top ⇔ E[(T₊ − S/2)²] ≥ floor/4."""
    return quartic_variance_floor_threshold(p) / 4


def top_plus_energy_qvar_threshold_wrong_drop4(p: int) -> Fraction:
    """Fail-eq: drop the 4."""
    return quartic_variance_floor_threshold(p)


def theorem_P_top_stratum_plus_energy_identity() -> dict:
    """On the top stratum for p=3 mod 4, ENERGY + Theorem F give
    |Z_ψ|² = (2 T₊ − S)² pointwise, hence
        E|Z_ψ|² = 4 E[(T₊ − S/2)²].
    QVAR is E[(T₊ − S/2)²] ≥ floor/4.  Certified exact at p=7 (the
    whole Max+ is k≤m=4; the k=4 slice is the top).  Does **not**
    bound the plus-energy moment; not a k≥7 close.
    Fail: drop the 4 in the threshold."""
    import numpy as np
    from e1_gmin_m4_prop15588 import maxplus, profiles_of

    p = 7
    S = profile_energy_total_S(p)
    Y = maxplus(p)
    P = profiles_of(p, Y)
    m = P.shape[1]
    eps = Y[:, 0]
    a = np.empty((len(Y), m), dtype=np.int64)
    for j in range(m):
        a[:, j] = ((P[:, j, :].astype(np.int64) - eps[:, None]) ** 2).sum(axis=1) // 4
    # p=7 square-direction psi signs, same order as directions()/profiles_of:
    # gens 1,7,9,12 → (+,+,−,−).  Independent of this module's psi table.
    w = np.array([1, 1, -1, -1], dtype=np.int64)
    B = a @ w
    Tplus = a[:, w > 0].sum(axis=1)
    active = (P != eps[:, None, None]).any(axis=2)
    k = active.sum(1)
    top = k == m
    Bt, Tt = B[top], Tplus[top]
    n_top = int(top.sum())
    # exact histogram: 3234 at ±28, 2352 at 0
    hist_ok = (
        int((Bt == 28).sum()) == 3234
        and int((Bt == -28).sum()) == 3234
        and int((Bt == 0).sum()) == 2352
        and n_top == 8820
        and np.all(Tt + a[top][:, w < 0].sum(axis=1) == S)
        and np.all(Bt == 2 * Tt - S)
    )
    EB2 = Fraction(int((Bt.astype(object) ** 2).sum()), n_top)
    plus_mom = Fraction(int(((Tt - S // 2) ** 2).sum()), n_top)
    floor = quartic_variance_floor_threshold(p)
    ok = (
        hist_ok
        and EB2 == 4 * plus_mom == Fraction(8624, 15)
        and plus_mom >= top_plus_energy_qvar_threshold(p)
        and top_plus_energy_qvar_threshold_wrong_drop4(p) != top_plus_energy_qvar_threshold(p)
        and EB2 >= floor
        and int((Bt ** 2).min()) == 0  # pointwise QVAR false on the top
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p": p,
        "n_top": n_top,
        "E_abs_Zpsi_sq": str(EB2),
        "plus_energy_second_moment": str(plus_mom),
        "four_times_plus_moment": str(4 * plus_mom),
        "pointwise_min_on_top": 0,
        "qvar_iff_plus_moment_ge_floor_over_4": True,
        "plus_moment_bound_proved": False,
        "theorem": (
            "Top stratum, p=3 mod 4: ENERGY and Theorem F give "
            "|Z_ψ|²=(2T₊−S)², so E|Z_ψ|²=4 E[(T₊−S/2)²].  QVAR iff "
            "that plus-energy moment ≥ floor/4.  Fail: drop the 4.  "
            "p=7 k=4: exact 8624/15=4·2156/15, pointwise min 0.  "
            "The variance bound is open; not a k≥7 close."
        ),
    }


def theorem_Q_equal_energy_not_a_lower_bound() -> dict:
    """If signs sum to 0 and every active profile has energy S/k, then
        E_eq = S² (m−k)/(k(m−1)).
    This reproduces the exact k=1 and k=3 formulas, and E_eq ≥ floor
    for every k≤m−1 (p=3 mod 4, p≥7).  It is **not** a lower bound:
    p=11 k=4 has E|Z_ψ|²=9438 < 10890=E_eq (15.589 J).  Fail: claim
    9438 ≥ 10890, or claim E_eq < floor at k=m−1.  Not a k≥7 close."""
    p11_k4_actual = Fraction(9438)
    p11_k4_eq = equal_energy_prediction(11, 4)
    p7_k1 = equal_energy_prediction(7, 1)
    p7_k3 = equal_energy_prediction(7, 3)
    from e1_gmin_m4_prop15589 import k1_quartic_variance, k3_quartic_variance_p3mod4

    k_le_m1_ok = all(
        equal_energy_meets_floor_k_le_m_minus_1(p, k)
        for p in (7, 11, 19, 23, 31, 43, 47)
        for k in range(1, (p + 1) // 2)
        if k != 2
    )
    top_zero = all(
        equal_energy_prediction(p, (p + 1) // 2) == 0
        for p in (7, 11, 19, 23)
    )
    # k=m-1 still meets the floor
    m1_ok = all(
        equal_energy_meets_floor_k_le_m_minus_1(p, (p + 1) // 2 - 1)
        for p in (7, 11, 19, 23, 31, 47)
    )
    ok = (
        p7_k1 == k1_quartic_variance(7)
        and p7_k3 == k3_quartic_variance_p3mod4(7)
        and p11_k4_eq == 10890
        and p11_k4_actual < p11_k4_eq
        and p11_k4_actual >= quartic_variance_floor_threshold(11)
        and k_le_m1_ok
        and top_zero
        and m1_ok
        and not (p11_k4_actual >= p11_k4_eq)
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "equal_energy_is_lower_bound": False,
        "p11_k4_actual": str(p11_k4_actual),
        "p11_k4_equal_energy": str(p11_k4_eq),
        "k_le_m_minus_1_E_eq_meets_floor": k_le_m1_ok,
        "top_E_eq_is_zero": top_zero,
        "theorem": (
            "Equal-energy signed prediction E_eq=S²(m−k)/(k(m−1)) matches "
            "k=1,3 exactly and meets the QVAR floor for every k≤m−1, but "
            "is not a lower bound: p=11 k=4 has 9438<10890=E_eq (15.589 J).  "
            "Fail: claim actual ≥ E_eq.  Does not prove k≥7 QVAR."
        ),
    }


# p=11 k=6 opposite-pair gram entries from k6_p11_full.npy (37,925,570 rows).
# Same-sign cluster is tight; opposite-sign splits into two classes.
P11_K6_OPP_GRAM_LOW = (
    2917.806242,
    2918.18948,
    2918.306289,
    2918.400229,
    2918.829256,
    2919.306179,
)
P11_K6_OPP_GRAM_HIGH = (2927.091322, 2927.795188, 2928.117552)
P11_K6_SAME_GRAM = (
    2923.676126,
    2923.732151,
    2923.892286,
    2924.023669,
    2924.103137,
    2924.289779,
)


def p13_k7_witness_record() -> dict:
    """Live CP-SAT witness; source of the pointwise/λ≠0 negatives."""
    path = ROOT / "evidence" / "k7_p13_cpsat_witness.json"
    return json.loads(path.read_text())


def theorem_S_nonzero_top_scalar_not_pointwise_qvar() -> dict:
    """The first live pair is λ≠0 and still below the floor, so a
    λ≠0 restriction does not restore pointwise (or single-vector)
    QVAR.  p=13≡1 (mod 4), so the p=3 (mod 4) T₊ identity does not
    apply to this witness either.  Fail: claim top_scalar=0.
    Ensemble λ≠0 remains open; not a k≥7 close."""
    rec = p13_k7_witness_record()
    top = int(rec["top_scalar"])
    abs2 = int(rec["min_abs_Zpsi_sq_seen"])
    thr = p13_k7_qvar_threshold()
    leading = rec["leading"]
    ok = (
        rec["p"] == 13
        and rec["k"] == 7
        and top == 7
        and top != 0
        and abs2 == P13_K7_POINTWISE_ABS_ZPSI_SQ == 2548
        and abs2 < thr
        and 13 % 4 == 1
        and all(int(x) != 0 for x in leading)
        and len(leading) == 7
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "top_scalar": top,
        "pointwise_abs_Zpsi_sq": abs2,
        "p_mod_4": 13 % 4,
        "lambda_nonzero_pointwise_qvar": False,
        "theorem": (
            "p=13 k=7 CP-SAT witness has top_scalar=7≠0 and full-support "
            "leading (4,7,1,3,3,1,7), with |Z_ψ|²=2548<10647/2.  Fail: "
            "claim top_scalar=0.  λ≠0 does not give pointwise QVAR; "
            "p≡1 (mod 4) so the T₊ identity is the wrong congruence.  "
            "Does not prove k≥7 QVAR."
        ),
    }


def theorem_R_two_type_opposite_pair_not_unique() -> dict:
    """A proof that treats E[a_i a_j] as a single number for every
    opposite-sign pair of top-stratum directions is false at p=11:
    the k=6 gram splits into a low class (~2918) and a high class (~2928).
    Same-sign pairs stay in a 0.6-wide cluster.  Analog of 15.596:
    uniqueness/full-rank of the pair-type is not the mechanism.
    Fail: claim all nine opposite entries equal.  Not a k≥7 close."""
    low = P11_K6_OPP_GRAM_LOW
    high = P11_K6_OPP_GRAM_HIGH
    same = P11_K6_SAME_GRAM
    spread_opp = max(high) - min(low)
    spread_same = max(same) - min(same)
    all_opp_equal = len({round(x, 3) for x in low + high}) == 1
    ok = (
        spread_opp > 5
        and spread_same < 1
        and min(high) - max(low) > 5
        and not all_opp_equal
        and len(low) == 6
        and len(high) == 3
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "unique_opposite_pair_moment": False,
        "p11_k6_opp_spread": spread_opp,
        "p11_k6_same_spread": spread_same,
        "theorem": (
            "p=11 top (k=6) directional-energy gram: opposite-sign pairs "
            "split 2917.8–2919.3 vs 2927.1–2928.1; same-sign stay inside "
            "2923.7–2924.3.  Fail: claim a unique opposite-pair moment.  "
            "Do not prove QVAR by two-type uniqueness.  Not a k≥7 close."
        ),
    }


def two_design_collision_EQ2(p: int) -> int:
    """Sum of 2-design E[y_e y_f] over directed same-fiber pair-pairs
    with fewer than 4 distinct indices.  On a square-direction pencil
    every collinear pair has C=1, so E[y_i y_j]=1/p, and the count of
    uniq=2 and uniq=3 pair-pairs is combinatorial:
        n_two = 2 p² (p−1),   n_three = 4 p² (p−1)(p−2),
        C_coll = n_two · 1 + n_three · (1/p) = 2p(p−1)(3p−4).
    Fail: drop the (3p−4) or the 2p."""
    return 2 * p * (p - 1) * (3 * p - 4)


def two_design_collision_EQ2_wrong_drop3p4(p: int) -> int:
    return 2 * p * (p - 1)


def theorem_T_two_design_does_not_determine_profile_energy_variance() -> dict:
    """Var(a_L) — and therefore the top-stratum T₊ identity — is not a
    2-design quantity.  Q_L = ∑_{x≠z, t_L(x)=t_L(z)} y_x y_z satisfies
    a_L = (p(p−1)+Q_L)/4 on Max+, and E[Q_L²] = C_coll + R_4 with
    C_coll = 2p(p−1)(3p−4) the collision (uniq<4) 2-design evaluation
    and R_4 the 4-distinct remainder.  At p=7, C_coll=1428 but
    E[Q²] = 1575252/409 on the eps=+1 half (equivalently from a).
    Fail: claim E[Q²]=C_coll (R_4=0).  Not a k≥7 close."""
    import numpy as np
    from e1_gmin_m4_prop15588 import directions, maxplus, profiles_of

    p = 7
    C_coll = two_design_collision_EQ2(p)
    Y = maxplus(p)
    Yp = Y[Y[:, 0] == 1]
    P = profiles_of(p, Yp)
    eps = Yp[:, 0]
    # first square direction
    a0 = ((P[:, 0, :].astype(np.int64) - eps[:, None]) ** 2).sum(axis=1) // 4
    Q0 = 4 * a0 - p * (p - 1)
    EQ2 = Fraction(int((Q0.astype(object) ** 2).sum()), len(Q0))
    ok = (
        C_coll == 1428
        and C_coll == 2 * p * (p - 1) * (3 * p - 4)
        and two_design_collision_EQ2_wrong_drop3p4(p) != C_coll
        and EQ2 != C_coll
        and EQ2 == Fraction(1575252, 409)
        and EQ2 - C_coll != 0
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "C_coll": C_coll,
        "E_Q2_p7_eps_plus": str(EQ2),
        "R4_nonzero": EQ2 != C_coll,
        "two_design_determines_Var_a": False,
        "theorem": (
            "On a square-direction pencil, 2-design collisions give "
            "C_coll=2p(p−1)(3p−4) (n_two=2p²(p−1), n_three=4p²(p−1)(p−2), "
            "E[y_i y_j]=1/p).  Fail: drop (3p−4).  At p=7, E[Q²]=1575252/409 "
            "≠ 1428, so the 4-distinct remainder is the main term.  "
            "2-design does not determine Var(a_L) or top-stratum QVAR."
        ),
    }


def theorem_U_second_moment_rayleigh_not_a_close() -> dict:
    """On a fixed k-stratum, p=3 mod 4, ENERGY+Theorem F give
        E|Z_ψ|² = wᵀ K w,   K_{L,M}=E[a_L a_M]
    (second-moment matrix, not centered covariance).  QVAR iff
        wᵀ K w ≥ 3q(q−1)/16.
    A sufficient condition is λ_min(K restricted to the sign-isotypic
    subspace) ≥ floor/‖w‖², but wᵀKw/‖w‖² is only a Rayleigh quotient
    unless that isotypic has multiplicity one.  At p=7 top, one
    opposite-sign orbital recovers the T₊ identity (8624/15).  At p=11
    opposite pairs split, so a single β_opp is false.  Fail: claim
    wᵀKw equals the all-ones Rayleigh S²/m.  Bound OPEN; not a k≥7 close.

    Referee: openai gpt-5.6-sol suggest_direction (orbital/Terwilliger
    blocks of K) + math_review PASS-WITH-NOTE on this reduction."""
    import numpy as np
    from e1_gmin_m4_prop15588 import maxplus, profiles_of

    p = 7
    S = profile_energy_total_S(p)
    floor = quartic_variance_floor_threshold(p)
    Y = maxplus(p)
    P = profiles_of(p, Y)
    m = P.shape[1]
    eps = Y[:, 0]
    a = np.empty((len(Y), m), dtype=np.int64)
    for j in range(m):
        a[:, j] = ((P[:, j, :].astype(np.int64) - eps[:, None]) ** 2).sum(axis=1) // 4
    w = np.array([1, 1, -1, -1], dtype=np.int64)
    active = (P != eps[:, None, None]).any(axis=2)
    top = active.sum(1) == m
    at = a[top]
    n_top = int(top.sum())
    Bt = at @ w
    wKw = Fraction(int(sum(int(x) * int(x) for x in Bt.tolist())), n_top)
    S2_over_m = Fraction(S * S, m)
    ok = (
        n_top == 8820
        and wKw == Fraction(8624, 15)
        and wKw >= floor
        and wKw != S2_over_m
        and int(w.sum()) == 0
        and S2_over_m == Fraction(1764)
        and floor == 441
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p": p,
        "wKw": str(wKw),
        "all_ones_rayleigh_S2_over_m": str(S2_over_m),
        "equals_all_ones_rayleigh": wKw == S2_over_m,
        "lambda_min_bound_proved": False,
        "rayleigh_not_automatically_eigenvalue": True,
        "theorem": (
            "k-stratum second-moment matrix K=E[aaᵀ]: E|Z_ψ|²=wᵀKw.  "
            "QVAR iff that Rayleigh ≥ floor.  Fail: claim wᵀKw=S²/m "
            f"(p=7 top: 8624/15 ≠ 1764).  λ_min of the sign-isotypic "
            "block is a sufficient unproved bound (multiplicity-one "
            "not assumed).  Not a k≥7 close."
        ),
    }


def singer_cycle_order(p: int) -> int:
    """|C_m| = |(F_q^{*2})/F_p^*| = (p+1)/2."""
    return (p + 1) // 2


def singer_cycle_order_wrong_p1(p: int) -> int:
    """Fail-eq: claim |C_m|=p+1 (that is |P¹(F_p)|)."""
    return p + 1


def singer_has_real_sign_character(p: int) -> bool:
    """Unique homomorphism C_m→{±1} besides 1 iff m even iff p≡3 (mod 4)."""
    return p % 4 == 3


def cyclic_sign_multiplicity_in_regular(m: int) -> int:
    """⟨sign, regular⟩ = 1 when m is even.  Cyclic ⇒ every irrep once."""
    if m % 2:
        return 0
    return 1


def cyclic_sign_multiplicity_wrong_claim_two(m: int) -> int:
    """Fail-eq: claim the sign appears twice in the regular representation."""
    return 2


def _field_primitive_root(p: int):
    from e1_gmin_m4_prop15588 import field_ctx

    q, mul, chi, tr = field_ctx(p)
    n = q - 1
    fac = set()
    t = n
    d = 2
    while d * d <= t:
        while t % d == 0:
            fac.add(d)
            t //= d
        d += 1
    if t > 1:
        fac.add(t)

    def powm(u, e):
        r, b = 1, u
        while e:
            if e & 1:
                r = mul(r, b)
            b = mul(b, b)
            e >>= 1
        return r

    g = next(
        x
        for x in range(2, q)
        if all(powm(x, n // f) != 1 for f in fac)
    )
    return q, mul, chi, tr, g


def singer_square_reps(p: int) -> list[int]:
    """η^k, k=0..m-1, η=γ^{p-1}: the m square F_p-lines."""
    q, mul, chi, tr, g = _field_primitive_root(p)
    eta = 1
    for _ in range(p - 1):
        eta = mul(eta, g)
    m = singer_cycle_order(p)
    reps = []
    x = 1
    for _k in range(m):
        reps.append(x)
        x = mul(x, eta)
    return reps


def _in_fp_star(r: int, p: int) -> bool:
    return r != 0 and (r // p) == 0


def singer_action_free_transitive(p: int) -> bool:
    """C_m acts regularly on the square directions: m distinct lines."""
    from e1_gmin_m4_prop15588 import field_ctx

    q, mul, chi, tr = field_ctx(p)
    reps = singer_square_reps(p)
    m = singer_cycle_order(p)
    if len(reps) != m or len(set(reps)) != m:
        return False
    if any(chi(r) != 1 for r in reps):
        return False

    def same_line(a, b):
        invb = next(x for x in range(1, q) if mul(b, x) == 1)
        return _in_fp_star(mul(a, invb), p)

    for i in range(m):
        for j in range(i + 1, m):
            if same_line(reps[i], reps[j]):
                return False
    return True


def singer_direction_perm(p: int) -> list[int]:
    """Permutation sending directions() order to Singer order."""
    from e1_gmin_m4_prop15588 import directions, field_ctx

    q, mul, chi, tr = field_ctx(p)
    reps = singer_square_reps(p)
    dir_reps = []
    seen = set()
    for g in range(1, q):
        if g in seen:
            continue
        line = [mul(t, g) for t in range(1, p)]
        seen.update(line)
        if chi(g) == 1:
            dir_reps.append(g)

    def same_line(a, b):
        invb = next(x for x in range(1, q) if mul(b, x) == 1)
        return _in_fp_star(mul(a, invb), p)

    perm = []
    for r in reps:
        perm.append(next(i for i, g in enumerate(dir_reps) if same_line(g, r)))
    if len(set(perm)) != singer_cycle_order(p):
        raise RuntimeError("Singer perm is not a bijection")
    return perm


def singer_sign_vector(p: int) -> list[int]:
    """w_k = (-1)^k in Singer order.  Only a character when m is even."""
    m = singer_cycle_order(p)
    return [1 if k % 2 == 0 else -1 for k in range(m)]


def _p7_singer_stratum_K() -> dict:
    """Exact K=E[aaᵀ] in Singer order on every nonempty p=7 stratum."""
    import numpy as np
    from e1_gmin_m4_prop15588 import maxplus, profiles_of

    p = 7
    m = 4
    S = profile_energy_total_S(p)
    perm = singer_direction_perm(p)
    w = np.array(singer_sign_vector(p), dtype=np.int64)
    Y = maxplus(p)
    P = profiles_of(p, Y)
    eps = Y[:, 0]
    a = np.empty((len(Y), m), dtype=np.int64)
    for j in range(m):
        a[:, j] = ((P[:, j, :].astype(np.int64) - eps[:, None]) ** 2).sum(axis=1) // 4
    a = a[:, perm]
    act = (P != eps[:, None, None]).any(axis=2)[:, perm]
    kvec = act.sum(1)
    out = {}
    for kv in sorted(set(int(x) for x in kvec.tolist())):
        mask = kvec == kv
        ak = a[mask]
        n = int(mask.sum())
        gram = ak.T.astype(object) @ ak.astype(object)
        K = [
            [Fraction(int(gram[i, j]), n) for j in range(m)]
            for i in range(m)
        ]
        c = K[0]
        circ = all(K[i][j] == c[(j - i) % m] for i in range(m) for j in range(m))
        Kw = [sum(K[i][j] * int(w[j]) for j in range(m)) for i in range(m)]
        lams = [Kw[i] * int(w[i]) for i in range(m)]
        evec = len(set(lams)) == 1
        lam = lams[0]
        EZ = Fraction(int(sum(int(x) * int(x) for x in (ak @ w).tolist())), n)
        nyq_c = c[0] - c[1] + c[2] - c[3]
        off_zero = all(c[d] == 0 for d in range(1, m))
        out[int(kv)] = {
            "n": n,
            "circulant": circ,
            "sign_eigenvector": evec,
            "lambda_nyquist": lam,
            "E_Z2": EZ,
            "c": c,
            "K_multiple_of_I": off_zero,
            "nyquist_identity": bool(
                circ and evec and lam * m == EZ and nyq_c == lam
            ),
        }
    return {"S": S, "perm": perm, "by_k": out}


def theorem_V_singer_circulant_nyquist() -> dict:
    """Singer C_m makes K circulant; p≡3 mod 4 the sign is an eigenvector.

    PROVED (algebra, every odd p)
      V1. C_m=(F_q^{*2})/F_p^* is cyclic of order m=(p+1)/2 and acts
          regularly on the m square directions.  Fail: |C_m|=p+1.
      V2. Square multiplications lie in Aut (15.588 II) and permute the
          m profile coordinates as the regular representation, preserving
          every k-stratum.  Hence K^{(k)}=E[aaᵀ|activity=k] is circulant
          in Singer order for every k, not only the top.

    PROVED (algebra, p≡3 mod 4)
      V3. m even ⇒ unique sign character of C_m, multiplicity one in the
          regular representation.  So w_j=(-1)^j is an eigenvector of
          every K^{(k)}.  QVAR iff the Nyquist eigenvalue
          λ_{m/2} ≥ floor/m.  Fail: claim a real sign character at p=13
          (m=7 odd).  Fail: multiplicity 2.

    CERTIFIED (p=7 Max+, k=1,3,4)
      V4. K is circulant, w is an eigenvector, λ_{m/2}·m = E|Z|² exactly.
          On the top, K is not a multiple of I and Nyquist ≠ S²/m.
          Fail: claim K=(S²/m)I at k=4; claim λ_Nyquist=S²/m.

    OPEN
      V5. λ_{m/2} ≥ floor/m on every k≥7, all p≥13.  For p≡1 mod 4 the
          real sign/Nyquist form is unavailable (V3); use 15.473.
          qvar_k_ge_7 stays False.
    """
    alg_ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        m = singer_cycle_order(p)
        has = singer_has_real_sign_character(p)
        mult = cyclic_sign_multiplicity_in_regular(m)
        reg = singer_action_free_transitive(p)
        row_ok = (
            m == (p + 1) // 2
            and singer_cycle_order_wrong_p1(p) != m
            and has is (p % 4 == 3)
            and (mult == 1 if has else mult == 0)
            and cyclic_sign_multiplicity_wrong_claim_two(m) != 1
            and reg
        )
        if p == 13:
            row_ok = row_ok and not has and m == 7 and m % 2 == 1
        rows[str(p)] = {
            "m": m,
            "sign_character": has,
            "sign_multiplicity": mult,
            "regular_action": reg,
            "ok": row_ok,
        }
        alg_ok = alg_ok and row_ok
    P7 = _p7_singer_stratum_K()
    S = P7["S"]
    by = P7["by_k"]
    floor = quartic_variance_floor_threshold(7)
    k1 = by[1]
    k3 = by[3]
    k4 = by[4]
    p7_ok = (
        k1["circulant"]
        and k3["circulant"]
        and k4["circulant"]
        and k1["sign_eigenvector"]
        and k3["sign_eigenvector"]
        and k4["sign_eigenvector"]
        and k1["nyquist_identity"]
        and k3["nyquist_identity"]
        and k4["nyquist_identity"]
        and k1["lambda_nyquist"] == Fraction(S * S, 4)
        and k1["E_Z2"] == Fraction(S * S)
        and k4["E_Z2"] == Fraction(8624, 15)
        and k4["lambda_nyquist"] * 4 == Fraction(8624, 15)
        and k1["K_multiple_of_I"]
        and not k4["K_multiple_of_I"]
        and k4["lambda_nyquist"] != Fraction(S * S, 4)
        and k4["E_Z2"] >= floor
    )
    return {
        "proved": bool(alg_ok and p7_ok),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "sign_is_eigenvector_p3mod4": True,
        "nyquist_iff_qvar_p3mod4": True,
        "p13_has_real_sign": singer_has_real_sign_character(13),
        "p7_top_K_multiple_of_I": k4["K_multiple_of_I"],
        "p7_top_lambda": str(k4["lambda_nyquist"]),
        "p7_top_EZ": str(k4["E_Z2"]),
        "p7_S2_over_m": str(Fraction(S * S, 4)),
        "by_p": rows,
        "p7_by_k": {
            str(k): {
                "n": rec["n"],
                "circulant": rec["circulant"],
                "sign_eigenvector": rec["sign_eigenvector"],
                "nyquist_identity": rec["nyquist_identity"],
                "lambda_nyquist": str(rec["lambda_nyquist"]),
                "E_Z2": str(rec["E_Z2"]),
                "K_multiple_of_I": rec["K_multiple_of_I"],
            }
            for k, rec in by.items()
        },
        "theorem": (
            "C_m acts regularly on square directions, so K^{(k)} is "
            "circulant every k.  p≡3 mod 4: unique sign character, "
            "multiplicity one, w eigenvector, QVAR iff λ_{m/2}≥floor/m.  "
            "Fail: |C_m|=p+1; sign at p=13; multiplicity 2; "
            "K=(S²/m)I on p=7 top (8624/15 ≠ 1764).  Bound OPEN.  "
            "p≡1 mod 4 has no real sign (m odd); use 15.473."
        ),
    }


def gcd_m_pminus1(p: int) -> int:
    """gcd((p+1)/2, p−1).  Equals 2 for every p≡3 (mod 4)."""
    from math import gcd

    return gcd((p + 1) // 2, p - 1)


def singer_linear_forms(p: int) -> list[tuple[int, int]]:
    """t_k(x)=Tr(α η^{-k} x) with α=γ^m, η=γ^{p-1}."""
    q, mul, chi, tr, g = _field_primitive_root(p)
    m = singer_cycle_order(p)
    alpha = 1
    for _ in range(m):
        alpha = mul(alpha, g)
    forms = []
    for r in singer_square_reps(p):
        inv = next(x for x in range(1, q) if mul(r, x) == 1)
        ck = mul(alpha, inv)
        forms.append((tr(mul(ck, 1)) % p, tr(mul(ck, p)) % p))
    return forms


def top_kernel_in_singer_gauge(p: int) -> list[int]:
    """Unique kernel of ∑ a_j t_j^{m-2} in Singer gauge, as a vector in F_p^m."""
    from math import comb

    m = singer_cycle_order(p)
    deg = m - 2
    forms = singer_linear_forms(p)
    M = [
        [
            (comb(deg, i) * pow(a, deg - i, p) * pow(b, i, p)) % p
            for a, b in forms
        ]
        for i in range(deg + 1)
    ]
    # right-nullspace via the existing F_p rank routine: try standard basis
    # differences against a particular solution from signed minors.
    rank = matrix_rank_modp(M, p)
    if rank != m - 1:
        raise RuntimeError(f"expected rank m-1, got {rank} at p={p}")
    A = [row[:] for row in M]
    n_r, n_c = len(A), m
    pivcol = [-1] * n_r
    r = 0
    used = []
    for c in range(n_c):
        pivot = None
        for i in range(r, n_r):
            if A[i][c] % p:
                pivot = i
                break
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        inv = pow(A[r][c] % p, p - 2, p)
        A[r] = [(v * inv) % p for v in A[r]]
        for i in range(n_r):
            if i == r:
                continue
            fac = A[i][c] % p
            if fac:
                A[i] = [(A[i][cc] - fac * A[r][cc]) % p for cc in range(n_c)]
        pivcol[r] = c
        used.append(c)
        r += 1
        if r == n_r:
            break
    free = [c for c in range(m) if c not in set(used)]
    if len(free) != 1:
        raise RuntimeError(f"expected 1-dim kernel, free={free}")
    f = free[0]
    v = [0] * m
    v[f] = 1
    for row, pc in enumerate(pivcol):
        if pc < 0:
            break
        v[pc] = (-A[row][f]) % p
    return v


def scale_to_sign(v: list[int], p: int) -> list[int]:
    m = len(v)
    w = singer_sign_vector(p)
    w_p = [1 if x == 1 else p - 1 for x in w]
    nz = next(i for i in range(m) if v[i] % p)
    scale = (w_p[nz] * pow(v[nz] % p, p - 2, p)) % p
    return [(x * scale) % p for x in v]


def theorem_W_top_kernel_is_alternating(
    primes_3=(7, 11, 19, 23, 31),
    primes_1=(5, 13, 17),
) -> dict:
    """In Singer gauge t_k=Tr(α η^{-k} x), the unique top kernel (degree
    m−2 among m square forms) is the alternating vector (−1)^k.

    PROVED (p≡3 mod 4, all such p)
      gcd(m, p−1)=2, so the only F_p-rational characters of C_m are
      trivial and sign.  Kernel is 1-dimensional with full support
      (theorem G).  All-ones is not in the kernel (fail: claim it is).
      Hence the kernel is the sign.  On the top stratum the leading
      coefficients are λ w, so the degree-(m−2) energy is equal in
      every direction and drops from Z=∑ w a (∑ w=0).  Top QVAR is
      a lower-degree remainder.  Not a bound.

    CERTIFIED also at p≡1 mod 4 (p=5,13,17): the same vector aligns,
    even though it is not a C_m-character (m odd).  Unused for Nyquist.

    Fail: gcd=p−1; kernel=all-ones at p=7.
    """
    gcd_ok = all(gcd_m_pminus1(p) == 2 for p in primes_3)
    gcd_fail = any(gcd_m_pminus1(p) == p - 1 for p in primes_3)
    rows = {}
    ok = gcd_ok and not gcd_fail
    for p in list(primes_3) + list(primes_1):
        m = singer_cycle_order(p)
        kap = top_kernel_in_singer_gauge(p)
        w_p = [1 if x == 1 else p - 1 for x in singer_sign_vector(p)]
        ones = [1] * m
        aligned = scale_to_sign(kap, p) == w_p
        not_ones = scale_to_sign(kap, p) != ones
        full = all(x % p for x in kap)
        row_ok = aligned and not_ones and full and len(kap) == m
        rows[str(p)] = {
            "aligned": aligned,
            "not_all_ones": not_ones,
            "full_support": full,
            "gcd": gcd_m_pminus1(p),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p3mod4_gcd_always_2": gcd_ok,
        "leading_energy_sign_blind_p3mod4": True,
        "qvar_bound_proved": False,
        "by_p": rows,
        "theorem": (
            "Singer-gauge top kernel is (−1)^k.  p≡3 mod 4: gcd(m,p−1)=2 "
            "so the only F_p-rational C_m-characters are 1 and sign; "
            "all-ones is not the kernel, hence sign.  Fail: all-ones at "
            "p=7.  Leading energy is sign-blind and drops from Z.  "
            "Top QVAR is a lower-degree remainder.  Bound OPEN."
        ),
    }


def _interp_fp_coeff(vals, p: int) -> list[int]:
    """Coefficients of the unique F_p interpolant of vals[s], s=0..p-1."""
    A = [[pow(s, d, p) for d in range(p)] for s in range(p)]
    b = [int(v) % p for v in vals]
    n = p
    piv = [-1] * n
    r = 0
    used = []
    for c in range(n):
        pivot = None
        for i in range(r, n):
            if A[i][c] % p:
                pivot = i
                break
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        b[r], b[pivot] = b[pivot], b[r]
        inv = pow(A[r][c] % p, p - 2, p)
        A[r] = [(v * inv) % p for v in A[r]]
        b[r] = (b[r] * inv) % p
        for i in range(n):
            if i == r:
                continue
            fac = A[i][c] % p
            if fac:
                A[i] = [(A[i][cc] - fac * A[r][cc]) % p for cc in range(n)]
                b[i] = (b[i] - fac * b[r]) % p
        piv[r] = c
        used.append(c)
        r += 1
    coeff = [0] * n
    for row, pc in enumerate(piv):
        if pc < 0:
            break
        coeff[pc] = b[row] % p
    return coeff


def p7_top_singer_leading_census() -> dict:
    """Complete p=7 top: F_p-leading in Singer gauge is λw, λ≠0."""
    import numpy as np
    from e1_gmin_m4_prop15588 import field_ctx, maxplus

    p = 7
    q, mul, chi, tr = field_ctx(p)
    m = 4
    forms = singer_linear_forms(p)
    w_p = [1 if x == 1 else p - 1 for x in singer_sign_vector(p)]
    T = np.zeros((m, q), dtype=np.int64)
    for j, (aa, bb) in enumerate(forms):
        for x in range(q):
            T[j, x] = (aa * (x % p) + bb * (x // p)) % p
    Y = maxplus(p)
    Yf = Y[:, 1:].astype(np.int64)
    eps = Y[:, 0].astype(np.int64)
    P = np.zeros((len(Y), m, p), dtype=np.int64)
    for j in range(m):
        for s in range(p):
            P[:, j, s] = Yf[:, T[j] == s].sum(axis=1)
    h = (P - eps[:, None, None]) // 2
    act = (P != eps[:, None, None]).any(axis=2)
    top = act.sum(1) == m
    ht = h[top]
    n_top = int(top.sum())
    n_align = n_zero = n_other = 0
    for n in range(n_top):
        A = [_interp_fp_coeff(ht[n, j], p)[2] for j in range(m)]
        if all(x == 0 for x in A):
            n_zero += 1
            continue
        if scale_to_sign(A, p) == w_p:
            n_align += 1
        else:
            n_other += 1
    return {
        "n_top": n_top,
        "n_align": n_align,
        "n_zero": n_zero,
        "n_other": n_other,
    }


def theorem_X_top_leading_is_sign_isotypic() -> dict:
    """On the top stratum, p≡3 mod 4, Singer-gauge F_p-leading coefficients
    lie in span{w}.  15.588 C: level-(m−2) vector is in the unique kernel.
    Theorem W: that kernel is the sign.  Fail: claim the kernel is all-ones;
    claim every top vector at p=7 has λ=0.

    p=7 complete Max+: all 8820 top vectors saturate (λ≠0) and align.
    p=11 top has a λ=0 class (15.589 K), so saturation is not the general
    top.  Does not bound the integer-energy remainder; not a k≥7 close.
    """
    W = theorem_W_top_kernel_is_alternating()
    C7 = p7_top_singer_leading_census()
    from e1_gmin_m4_prop15589 import theorem_K_full_support_top_degree_mixing

    K = theorem_K_full_support_top_degree_mixing()
    p11_zero = int(K["p11"]["top_zero_count"])
    p11_tot = int(K["p11"]["full_support_count"])
    ok = (
        bool(W["proved"])
        and C7["n_top"] == 8820
        and C7["n_align"] == 8820
        and C7["n_zero"] == 0
        and C7["n_other"] == 0
        and p11_zero > 0
        and p11_zero < p11_tot
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "p7_all_saturate": C7["n_zero"] == 0,
        "p7_all_align": C7["n_align"] == C7["n_top"],
        "p11_has_lambda0": p11_zero > 0,
        "p7_census": C7,
        "qvar_bound_proved": False,
        "theorem": (
            "Top, p≡3 mod 4, Singer gauge: F_p-leading ∈ span{w} "
            "(15.588 C + W).  Fail: all-ones kernel; fail: p=7 top has "
            "λ=0 (0 of 8820).  p=11 top has a λ=0 class (15.589 K), so "
            "saturation is not general.  Integer-energy remainder unbound.  "
            "Not a k≥7 close."
        ),
    }


def spectral_EZ_prefactor(p: int) -> int:
    """E|Z_ψ|² = E|∑_{Ω} ψ ρ|² / (16q).  16q = 16p²."""
    return 16 * q_of(p)


def spectral_EZ_prefactor_wrong_16p(p: int) -> int:
    """Fail-eq: 16q ↦ 16p."""
    return 16 * p


def phi_sum_on_squares(p: int) -> int:
    """∑_{r∈F_q^{*2}} Φ(r) = 2 q² (q−1).

    Plancherel: ∑_Ω ρ = q(q−1) pointwise.  Square-multiplication Aut
    makes Φ(u,v)=E[ρ(u)ρ(v)] a class function of r=v/u on the squares.
    |Ω| pairs per square r, so |Ω| ∑ Φ = E[(∑ ρ)²] = q²(q−1)².
    """
    q = q_of(p)
    return 2 * q * q * (q - 1)


def spectral_qvar_floor_on_psi_mode(p: int) -> int:
    """QVAR iff E|∑_Ω ψ ρ|² ≥ 3 q² (q−1)."""
    q = q_of(p)
    return 3 * q * q * (q - 1)


def spectral_qvar_floor_wrong_drop_q(p: int) -> int:
    """Fail-eq: 3q²(q−1) ↦ 3q(q−1)."""
    q = q_of(p)
    return 3 * q * (q - 1)


def theorem_Y_psi_mode_of_Phi_not_a_close() -> dict:
    """QVAR is the ψ-Fourier mode of Φ(r)=E[|ŷ(u)|² |ŷ(ru)|²].

    PROVED (Max+-free, every odd p)
      Y1. ŷ(0)²=p² and ∑_Ω |ŷ|²=q(q−1) pointwise (15.588 A + Plancherel).
      Y2. Square multiplications ∈ Aut ⇒ Φ(u,v) depends only on r=v/u
          among squares.  Then ∑_{r square} Φ(r)=2q²(q−1).
      Y3. Autocorrelation: Φ(r)≤Φ(1).  Bochner: the ψ-mode is ≥0, i.e.
          E|Z_ψ|²≥0, strictly weaker than QVAR (theorem J).
      Y4. 15.473 + zhat=2Dhat on Ω: E|Z_ψ|²=E|∑_{u∈Ω} ψ(u) ρ(u)|²/(16q)
          with ρ=|ŷ|².  QVAR iff E|∑ ψ ρ|² ≥ 3q²(q−1).
          Fail: 16q↦16p; drop a q in 3q²(q−1).
      Y5. If Φ were constant on squares then |P|=|N|=(q−1)/4 ⇒ the
          ψ-mode vanishes ⇒ E|Z|²=0, false at p=5,7 (15.589 E).
          Fail: claim Φ constant at p=7 (then EZ=0).

    OPEN: a lower bound on the ψ-mode (equivalently Nyquist of K).
    Does not prove k≥7 QVAR.
    """
    E = theorem_E_exceptional_quartic_variance()
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 47):
        q = q_of(p)
        pref = spectral_EZ_prefactor(p)
        floor_z = quartic_variance_floor_threshold(p)
        mode_floor = spectral_qvar_floor_on_psi_mode(p)
        row_ok = (
            pref == 16 * q == 16 * p * p
            and spectral_EZ_prefactor_wrong_16p(p) != pref
            and phi_sum_on_squares(p) == 2 * q * q * (q - 1)
            and Fraction(mode_floor, pref) == floor_z
            and spectral_qvar_floor_wrong_drop_q(p) != mode_floor
            and (q - 1) % 4 == 0
        )
        rows[str(p)] = {
            "prefactor_16q": pref,
            "phi_sum": phi_sum_on_squares(p),
            "mode_floor": mode_floor,
            "ok": row_ok,
        }
        ok = ok and row_ok
    ez5 = Fraction(E["by_p"]["5"]["E_abs_Zpsi_sq"])
    ez7 = Fraction(E["by_p"]["7"]["E_abs_Zpsi_sq"])
    const_killed = (
        ez5 > 0
        and ez7 > 0
        and ez7 == Fraction(317520, 409)
        and ez7 >= quartic_variance_floor_threshold(7)
    )
    return {
        "proved": bool(ok and const_killed and bool(E["proved_census"])),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "constant_Phi_implies_EZ_zero": True,
        "p7_EZ": str(ez7),
        "by_p": rows,
        "theorem": (
            "E|Z_ψ|²=E|∑_Ω ψ ρ|²/(16q).  QVAR iff the ψ-mode of Φ is "
            "at least 3q²(q−1).  Fail: 16q↦16p; drop q in the mode floor; "
            "Φ constant (would force EZ=0, false at p=7: 317520/409).  "
            "∑_{squares} Φ=2q²(q−1); Φ(r)≤Φ(1); Bochner too weak.  "
            "Bound OPEN.  Not a k≥7 close."
        ),
    }


def phi1_cs_kstratum(p: int, k: int) -> Fraction:
    """CS: Φ(1)=E[ρ(u)²] ≥ 4q² m/k on a C_m-invariant k-stratum.
    Inactive lines have ρ=0; P(u active)=k/m; E[ρ|active]=2qm/k."""
    q = q_of(p)
    m = singer_cycle_order(p)
    return Fraction(4 * q * q * m, k)


def three_level_S_box(p: int, phi1: Fraction) -> Fraction:
    """S_□ if Φ is 3-level (Φ(±1)=φ1, constant off).  Not a lower bound."""
    q = q_of(p)
    return Fraction(2 * (q - 1) * (phi1 - 4 * q * q), q - 5)


def three_level_phi1_for_qvar(p: int) -> Fraction:
    """Φ(1) that would make 3-level S_□=6q².  q²(7q-19)/(q-1)."""
    q = q_of(p)
    return Fraction(q * q * (7 * q - 19), q - 1)


def p7_top_phi1() -> Fraction:
    """Φ(1) on the p=7 top, one Ω frequency, complete Max+."""
    import numpy as np
    from e1_gmin_m4_prop15588 import field_ctx, maxplus, profiles_of

    p = 7
    q, mul, chi, tr = field_ctx(p)
    Y = maxplus(p)
    P = profiles_of(p, Y)
    eps = Y[:, 0]
    act = (P != eps[:, None, None]).any(axis=2)
    top = act.sum(1) == (p + 1) // 2
    Yf = Y[top, 1:].astype(np.int64)
    u = next(c for c in range(1, q) if chi(c) == 1)
    trx = np.array([tr(mul(u, x)) for x in range(q)], dtype=np.int8)
    omega = np.exp(2j * np.pi * np.arange(p) / p)
    rho = np.abs((Yf * omega[trx]).sum(axis=1)) ** 2
    n = int(top.sum())
    s2 = Fraction(int(np.round((rho ** 2).sum())), n)
    return s2


def theorem_Z_three_level_not_a_lower_bound() -> dict:
    """k-stratum CS and the 3-level formula do not prove QVAR.

    PROVED
      Z1. On a C_m-invariant k-stratum, Φ(1)≥4q² m/k (CS; inactive
          lines carry ρ=0).  Fail: drop m/k (on k=1 this is 4q² m).
      Z2. Square Aut includes −1, so Φ(1)=Φ(−1).
      Z3. If Φ were 3-level, S_□=2(q−1)(Φ(1)−4q²)/(q−5) and QVAR
          would be Φ(1)≥q²(7q−19)/(q−1).  Fail: claim this 3-level
          S_□ lower-bounds the true ψ-mode.
      Z4. p=7 top: Φ(1)=266168/15=17744.5̅ < 8q²=19208, so the
          global Wick Q(1)=8q² of 15.279 D is not per-stratum.
          CS=4q²=9604 < Φ(1).  3-level would already clear 6q²,
          but actual Φ has ≥4 value clusters (Y probe), and the
          true mode is larger than 3-level — 3-level is not worst
          case (adversary puts mass on {ψ=−1}).  Bound OPEN.

    Does not prove k≥7 QVAR.  Do not use 15.279's global 3-level Wick
    as a k-stratum identity.
    """
    p = 7
    q = 49
    phi1 = p7_top_phi1()
    cs_top = phi1_cs_kstratum(p, 4)
    eight = 8 * q * q
    three = three_level_S_box(p, phi1)
    need = three_level_phi1_for_qvar(p)
    s_box_floor = 6 * q * q
    ez = Fraction(8624, 15)
    # true S_□ from EZ: E|M|²=16q EZ, S_□=E|M|²/|Ω|, |Ω|=(q-1)/2
    true_S = Fraction(16 * q * ez * 2, q - 1)
    ok = (
        cs_top == Fraction(4 * q * q)
        and phi1_cs_kstratum(p, 1) == Fraction(4 * q * q * 4)
        and phi1 > cs_top
        and phi1 < eight
        and phi1 > need
        and three > s_box_floor
        and true_S > three
        and true_S / (q * q) >= 6
        and phi1 == Fraction(266168, 15)
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "three_level_is_lower_bound": False,
        "p7_top_phi1": str(phi1),
        "p7_top_8q2": eight,
        "p7_top_cs": str(cs_top),
        "p7_three_level_S": str(three),
        "p7_true_S": str(true_S),
        "wick_8q2_not_per_stratum": phi1 != eight,
        "theorem": (
            "k-stratum CS Φ(1)≥4q² m/k (fail: drop m/k).  Φ(1)=Φ(−1).  "
            "3-level S_□=2(q−1)(Φ(1)−4q²)/(q−5) is not a lower bound "
            "(p=7 top: true S_□>3-level; Wick 8q² is global-only, "
            "Φ(1)=266168/15<19208).  Bound OPEN.  Not a k≥7 close."
        ),
    }


def quartic_psi_table(p: int) -> list[complex]:
    """ψ(γ^k)=i^k on F_q^*, the unique (up to conj) order-4 character
    with ψ²=χ."""
    q, mul, _chi, _tr, g = _field_primitive_root(p)
    psi = [0j] * q
    x = 1
    for k in range(q - 1):
        psi[x] = (1j) ** k
        x = mul(x, g)
    return psi


def psi_on_fp_generator(p: int) -> complex:
    """ψ(γ^{p+1})=i^{p+1}=−i^{p-1}.  F_p^*=⟨γ^{p+1}⟩."""
    return (1j) ** (p + 1)


def psi_on_fp_generator_wrong_drop_i2(p: int) -> complex:
    """Fail-eq: i^{p+1} ↦ i^{p-1} (drop the extra i²=−1)."""
    return (1j) ** (p - 1)


def eta_is_trivial(p: int) -> bool:
    """ψ|F_p^*=1  ⇔  i^{p+1}=1  ⇔  p≡3 (mod 4)."""
    return psi_on_fp_generator(p) == 1


def eta_is_legendre(p: int) -> bool:
    """ψ|F_p^*=χ_p  ⇔  i^{p+1}=−1  ⇔  p≡1 (mod 4)."""
    return psi_on_fp_generator(p) == -1


def plancherel_off_zero_energy(sigma: list[int]) -> int:
    """∑_{λ≠0}|σ̂(λ)|² = p∑σ²−(∑σ)².  Integer, no DFT."""
    p = len(sigma)
    s1 = sum(sigma)
    s2 = sum(x * x for x in sigma)
    return p * s2 - s1 * s1


def four_p_times_a(sigma: list[int], eps: int) -> int:
    """4p a, a=(1/4)∑(σ−eps)²."""
    p = len(sigma)
    a = sum((x - eps) * (x - eps) for x in sigma) // 4
    return 4 * p * a


def four_p_times_a_wrong_drop_p(sigma: list[int], eps: int) -> int:
    """Fail-eq: 4p a ↦ 4 a."""
    return four_p_times_a(sigma, eps) // len(sigma)


def _k1_two_valued_profile(p: int, eps: int = 1) -> list[int]:
    """k=1 profile: n_+=(p+eps)/2 entries +p, rest −p, sum p eps."""
    n_plus = (p + eps) // 2
    return [p] * n_plus + [-p] * (p - n_plus)


def _line_dft_census(p: int) -> dict:
    """Complete Max+ check of ŷ=σ̂ and the ξ-circulant, p=5 or 7."""
    import numpy as np
    from e1_gmin_m4_prop15588 import directions, field_ctx, maxplus, profiles_of

    q, mul, chi, tr = field_ctx(p)
    psi = np.array(quartic_psi_table(p), dtype=np.complex128)
    Y = maxplus(p)
    N = len(Y)
    Yf = Y[:, 1:].astype(np.int64)
    P = profiles_of(p, Y)
    m = P.shape[1]
    eps = Y[:, 0]
    sq, _ = directions(p)
    alphas = []
    for t_of, _form in sq:
        a = next(
            c
            for c in range(1, q)
            if all(tr(mul(c, x)) == int(t_of[x]) for x in range(q))
        )
        alphas.append(a)
    omega = np.exp(2j * np.pi * np.arange(p) / p)
    s = np.arange(p)
    lam = np.arange(p)
    phase = omega[(lam[:, None] * s[None, :]) % p]
    sigma_hat = np.einsum(
        "njs,ls->njl", P.astype(np.float64), phase, optimize=True
    )
    trx = np.empty((q, q), dtype=np.int8)
    for c in range(q):
        for x in range(q):
            trx[c, x] = tr(mul(c, x))
    yhat = np.empty((N, q), dtype=np.complex128)
    for c in range(q):
        yhat[:, c] = (Yf * omega[trx[c]]).sum(axis=1)
    dft_err = 0.0
    for j, a in enumerate(alphas):
        for lv in range(1, p):
            u = mul(lv, a)
            dft_err = max(
                dft_err,
                float(np.abs(yhat[:, u] - sigma_hat[:, j, lv]).max()),
            )
    rho_hat = np.abs(sigma_hat) ** 2
    xi = np.zeros((N, m), dtype=np.complex128)
    for j, a in enumerate(alphas):
        for lv in range(1, p):
            u = mul(lv, a)
            xi[:, j] += psi[u] * rho_hat[:, j, lv]
    spec = np.zeros(N, dtype=np.complex128)
    chi_t = np.array([chi(c) for c in range(q)])
    mean_pow = np.abs(yhat).mean(axis=0)
    sq_idx = np.where(chi_t == 1)[0]
    nsq_idx = np.where(chi_t == -1)[0]
    e_sq = float((np.abs(yhat[:, sq_idx]) ** 2).mean())
    e_nsq = float((np.abs(yhat[:, nsq_idx]) ** 2).mean())
    Omega = sq_idx if e_sq > e_nsq else nsq_idx
    rho = np.abs(yhat) ** 2
    spec = (psi[Omega][None, :] * rho[:, Omega]).sum(axis=1)
    xi_err = float(np.max(np.abs(xi.sum(1) - spec)))
    act = (P != eps[:, None, None]).any(axis=2)
    kvec = act.sum(1)
    perm = singer_direction_perm(p)
    pref = spectral_EZ_prefactor(p)
    by_k = {}
    for kv in sorted(set(int(x) for x in kvec.tolist())):
        mask = kvec == kv
        n_k = int(mask.sum())
        xk = xi[mask][:, perm]
        gram = (xk.conj().T @ xk) / n_k
        circ_err = 0.0
        c0 = gram[0]
        for i in range(m):
            for j in range(m):
                circ_err = max(
                    circ_err, abs(gram[i, j] - c0[(j - i) % m])
                )
        sabs = np.abs(xk.sum(1)) ** 2
        Esum = Fraction(int(np.round(sabs.sum())), n_k)
        by_k[int(kv)] = {
            "n": n_k,
            "circ_err": circ_err,
            "E_abs_sum_xi_sq": Esum,
            "E_Z2": Esum / pref,
            "clears": Esum / pref >= quartic_variance_floor_threshold(p),
        }
    a_prof = ((P.astype(np.int64) - eps[:, None, None]) ** 2).sum(axis=2) // 4
    W_off = rho_hat[:, :, 1:].sum(axis=2)
    w_vs_a = float(np.max(np.abs(W_off - 4 * p * a_prof)))
    return {
        "N": N,
        "m": m,
        "dft_err": dft_err,
        "xi_err": xi_err,
        "w_vs_4pa": w_vs_a,
        "by_k": by_k,
        "Omega_chi": int(chi(int(Omega[0]))),
    }


def theorem_AA_line_dft_both_congruences() -> dict:
    """Line-DFT form of the ψ-mode, both residue classes.

    PROVED (Max+-free algebra, every odd p)
      AA3. ψ(γ^{p+1})=i^{p+1}.  p≡3 (mod 4) ⇒ this is 1, so
           η=ψ|F_p^*=1.  p≡1 (mod 4) ⇒ this is −1, and
           ⟨γ^{p+1}⟩=F_p^* so η=Legendre.  Fail: i^{p+1}↦i^{p-1}
           (then p=7 would look like Legendre).  Fail: claim η
           trivial at p=13.
      AA4. Plancherel on each profile: ∑_{λ≠0}|σ̂|²=p∑σ²−(∑σ)².
           If ∑σ=p eps then this equals 4p a, a=(1/4)∑(σ−eps)².
           Fail: drop the p (4p a ↦ 4a).

    PROVED (15.588 B + Poisson on dual lines; certified p=5,7)
      AA1. ŷ(λ α_j)=σ̂_j(λ) for λ∈F_p^*.  Inactive profiles have
           σ̂(λ≠0)=0, so ρ lives only on active dual lines.
      AA2. ξ_j:=∑_{u∈ℓ_j} ψ(u) ρ(u)=∑_{λ≠0} ψ(λ α_j)|σ̂_j(λ)|²
           and ∑_Ω ψ ρ=∑_j ξ_j.  For p≡3, ξ_j=ψ(α_j)·4p a_j
           (15.589 F).  For p≡1, ξ_j is the Legendre-twisted
           1D energy, not 4p a_j.
      AA5. p≡1: W_j=G(η)∑_{d≠0} η(d) A_σ(d) with G(η)=√p and
           A_σ(d)=∑_s σ(s)σ(s−d).  Certified p=5 (G=√5).
      AA6. K_ξ=E[ξξ^*|k] is circulant in Singer order (p=5 k=1,3;
           p=7 k=1,3,4).  QVAR iff m λ_0(K_ξ)=E|∑ξ|²≥3q²(q−1),
           the trivial mode.  p≡3: this is Nyquist of K_a.
           Fail: claim a real sign eigenvector of K_a at p=5
           (m=3 odd).  Certified EZ: p=5 k=1,3 equal 500 and 180;
           p=7 k=1,3,4 equal S², S²(m−3)/(3(m−1)), 8624/15.

    OPEN: a lower bound on λ_0(K_ξ) on every live k≥7, including
    the top, for every p≥13.  The unrestricted Paley/Weil operator
    (no profile constraint) is the wrong form.  Not a k≥7 close.
    """
    from e1_gmin_m4_prop15589 import (
        k1_quartic_variance,
        k3_quartic_variance_p3mod4,
    )

    alg_ok = True
    eta_rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        g = psi_on_fp_generator(p)
        g_wrong = psi_on_fp_generator_wrong_drop_i2(p)
        row_ok = (
            g == (1j) ** (p + 1)
            and g_wrong != g
            and eta_is_trivial(p) == (p % 4 == 3)
            and eta_is_legendre(p) == (p % 4 == 1)
            and (p % 4 != 1 or not eta_is_trivial(p))
            and (p != 13 or eta_is_legendre(p))
            and singer_has_real_sign_character(p) == (p % 4 == 3)
        )
        eta_rows[str(p)] = {"psi_Fp_gen": str(g), "ok": row_ok}
        alg_ok = alg_ok and row_ok
    sig = _k1_two_valued_profile(7, 1)
    pl = plancherel_off_zero_energy(sig)
    four = four_p_times_a(sig, 1)
    pl_ok = (
        sum(sig) == 7
        and pl == four
        and four == 4 * 7 * (sum((x - 1) * (x - 1) for x in sig) // 4)
        and four_p_times_a_wrong_drop_p(sig, 1) != four
    )
    C5 = _line_dft_census(5)
    C7 = _line_dft_census(7)
    k1_5 = C5["by_k"][1]["E_Z2"]
    k3_5 = C5["by_k"][3]["E_Z2"]
    k1_7 = C7["by_k"][1]["E_Z2"]
    k3_7 = C7["by_k"][3]["E_Z2"]
    k4_7 = C7["by_k"][4]["E_Z2"]
    cert_ok = (
        C5["dft_err"] < 1e-10
        and C7["dft_err"] < 1e-10
        and C5["xi_err"] < 1e-8
        and C7["xi_err"] < 1e-8
        and C7["w_vs_4pa"] < 1e-6
        and C5["by_k"][1]["circ_err"] < 1e-8
        and C5["by_k"][3]["circ_err"] < 1e-8
        and C7["by_k"][4]["circ_err"] < 1e-8
        and k1_5 == k1_quartic_variance(5) == Fraction(500)
        and k3_5 == Fraction(180)
        and k1_7 == k1_quartic_variance(7) == Fraction(84 * 84)
        and k3_7 == k3_quartic_variance_p3mod4(7)
        and k4_7 == Fraction(8624, 15)
        and C5["by_k"][1]["clears"]
        and C5["by_k"][3]["clears"]
        and C7["by_k"][4]["clears"]
        and C5["by_k"][1]["E_abs_sum_xi_sq"]
        == k1_5 * spectral_EZ_prefactor(5)
        and not singer_has_real_sign_character(5)
    )
    return {
        "proved": bool(alg_ok and pl_ok and cert_ok),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "eta_trivial_iff_p3mod4": True,
        "p13_eta_legendre": eta_is_legendre(13),
        "p5_k1_EZ": str(k1_5),
        "p5_k3_EZ": str(k3_5),
        "p7_k4_EZ": str(k4_7),
        "p5_census": {str(k): {"n": rec["n"], "E_Z2": str(rec["E_Z2"])} for k, rec in C5["by_k"].items()},
        "p7_census": {str(k): {"n": rec["n"], "E_Z2": str(rec["E_Z2"])} for k, rec in C7["by_k"].items()},
        "eta_rows": eta_rows,
        "theorem": (
            "ŷ(λ α_j)=σ̂_j(λ); ∑_Ω ψ ρ=∑_j ξ_j with "
            "ξ_j=∑_{ℓ_j} ψ ρ.  η=ψ|F_p^* is 1 (p≡3) or Legendre "
            "(p≡1); fail: i^{p+1}↦i^{p-1}; fail: η trivial at p=13.  "
            "Plancherel W=4p a (fail: drop p).  K_ξ circulant, QVAR "
            "is the trivial mode (p=5 k=1,3: 500, 180; p=7 top: "
            "8624/15).  Unrestricted Weil is the wrong operator.  "
            "Bound OPEN.  Not a k≥7 close."
        ),
    }


def monomial_power_sum(p: int, d: int) -> int:
    """μ_d = ∑_{s=0}^{p-1} s^d."""
    return sum(s**d for s in range(p))


def plancherel_J_eta1(p: int, d: int, e: int) -> int:
    """η=1: J_{de}=p μ_{d+e}−μ_d μ_e."""
    return p * monomial_power_sum(p, d + e) - monomial_power_sum(
        p, d
    ) * monomial_power_sum(p, e)


def plancherel_J_eta1_wrong_uncentered(p: int, d: int, e: int) -> int:
    """Fail-eq: drop the rank-1 centering μ_d μ_e."""
    return p * monomial_power_sum(p, d + e)


def plancherel_J11_closed(p: int) -> int:
    """J_{11}(η=1)=p²(p²−1)/12."""
    return p * p * (p * p - 1) // 12


def coupled_profile_dim(k: int) -> int:
    """dim V_{p,k} (p>k−2): k constants plus ker dims k−d−1 for d=1..k−2.

    Includes the 1-dimensional leading slot (λ=0 allowed).  Independent
    of p once reduced polynomials of degree <p exist.
    """
    if k < 1:
        return 0
    return k + (k - 2) * (k - 1) // 2 if k >= 2 else k


def coupled_profile_dim_wrong_drop_lambda0(k: int) -> int:
    """Fail-eq: exclude vanishing leading (subtract dim ker_{d=k-2}=1)."""
    return coupled_profile_dim(k) - 1


def coupled_profile_dim_wrong_drop_constants(k: int) -> int:
    """Fail-eq: drop the unconstrained constants (subtract k)."""
    return coupled_profile_dim(k) - k


def _monomial_dft(p: int, deg: int):
    """G[d,λ]=∑_s s^d ω^{λ s}.  Tiny (deg+1)×p."""
    import numpy as np

    omega = np.exp(2j * np.pi * np.arange(p) / p)
    s = np.arange(p, dtype=np.float64)
    G = np.empty((deg + 1, p), dtype=np.complex128)
    for d in range(deg + 1):
        sd = s**d
        G[d] = np.array(
            [(sd * omega[(lam * np.arange(p)) % p]).sum() for lam in range(p)]
        )
    return omega, G


def _jacobi_gram_eta(p: int, G, eta) -> "object":
    import numpy as np

    Gd = G[:, 1:]
    w = eta[1:]
    return (Gd * w[None, :]) @ Gd.conj().T


def _fp_legendre_vec(p: int):
    import numpy as np

    eta = np.zeros(p, dtype=np.float64)
    for x in range(1, p):
        eta[x] = 1.0 if pow(x, (p - 1) // 2, p) == 1 else -1.0
    return eta


def jacobi_legendre_inertia(p: int, deg: int) -> dict:
    """Signed inertia of J_χ on monomials of degree ≤deg.  No Max+ census."""
    import numpy as np

    _omega, G = _monomial_dft(p, deg)
    J = _jacobi_gram_eta(p, G, _fp_legendre_vec(p))
    ev = np.linalg.eigvalsh(0.5 * (J + J.conj().T).real)
    return {
        "eigs": [float(x) for x in ev],
        "n_neg": int((ev < -1e-6).sum()),
        "n_pos": int((ev > 1e-6).sum()),
        "min_eig": float(ev.min()),
        "J11": float(J[1, 1].real),
        "J00": float(J[0, 0].real),
    }


def G1_closed_form_err(p: int) -> float:
    """max |G_1(λ) − p/(ω^λ−1)| on λ≠0."""
    import numpy as np

    omega, G = _monomial_dft(p, 1)
    closed = p / (omega[1:] - 1)
    return float(np.max(np.abs(G[1, 1:] - closed)))


def G1_closed_form_wrong_drop_p_err(p: int) -> float:
    """Fail-eq: p/(ω−1) ↦ 1/(ω−1)."""
    import numpy as np

    omega, G = _monomial_dft(p, 1)
    wrong = 1.0 / (omega[1:] - 1)
    return float(np.max(np.abs(G[1, 1:] - wrong)))


def S_p_from_G1(p: int) -> float:
    """∑_{λ≠0} χ(λ)/|ω^λ−1|² = J_{11}(χ)/p²."""
    import numpy as np

    omega, G = _monomial_dft(p, 1)
    eta = _fp_legendre_vec(p)
    g1 = G[1, 1:]
    return float((eta[1:] * (np.abs(g1) ** 2) / (p * p)).sum().real)


def theorem_AB_jacobi_gram_not_unrestricted_weil() -> dict:
    """Jacobi Gram of W on polynomial profiles; unrestricted Weil killed.

    PROVED (Max+-free algebra / tiny DFT, no Max+ census, no p=13 orbits)
      AB1. G_0(λ≠0)=0 and G_1(λ)=p/(ω^λ−1).  Fail: drop p in G_1.
      AB2. J_{de}(η)=∑_{λ≠0} η(λ) G_d(λ) conj(G_e).  For a real
           polynomial σ=∑ c_d s^d, W(σ)=c^T J c.  Fail: include λ=0
           in J (then J_{00}=p² and W picks up ŷ(0)²).
      AB3. η=1: J_{de}=p μ_{d+e}−μ_d μ_e with μ_d=∑_{s=0}^{p-1} s^d.
           J_{00}=0 (constant mode = λ=0 of the 1D DFT).  J_{11}=
           p²(p²−1)/12.  The Gram is PSD: discrete Chebyshev /
           Christoffel–Darboux on {0,…,p−1}.  Fail: drop μ_d μ_e;
           fail: claim a negative pivot of J_{η=1}.
      AB4. p≡1 (η=χ): unrestricted J on deg≤5 is indefinite
           (p=13: ≥1 negative eigenvalue).  Affine restriction
           deg≤1 is PSD at p=5, with J_{11}=p² S_p and
           S_p=∑ χ(λ)/|ω^λ−1|² recovering 15.589 G / L(2,χ).
           Fail: claim J_χ ⪰ 0 at p=13 on deg≤5 (the live top
           degree).  Unrestricted Weil is not a QVAR certificate.
      AB5. Coupled coefficient space V_{p,k} has dimension
           k+(k−2)(k−1)/2 (constants free; d=1..k−2 kernel
           dims k−d−1).  Includes the leading 1-dim slot
           (λ=0 allowed).  Fail: drop λ=0 (dim−1); fail: drop
           constants (dim−k).  k=7: 22.  k=4: 7.
      AB6. The negative directions of J_χ at p=13 live in
           degree ≥2, so keeping λ=0 (constants / vanishing
           leading) does not restore positivity of the
           *unrestricted* Gram.  Coupling across the k
           directions is required.

    OPEN: LDL*/Christoffel–Darboux of the coupled form
        Q(c)=∑_j ψ(α_j) c^{(j)T} J c^{(j)}
    on V_{p,k}, equivalently λ_min(H|_{S_{p,k}}, G_profile).
    Named quantity, not a bound.  Not a k≥7 close.
    """
    alg = True
    j11_rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        j00 = plancherel_J_eta1(p, 0, 0)
        j11 = plancherel_J_eta1(p, 1, 1)
        row_ok = (
            j00 == 0
            and j11 == plancherel_J11_closed(p)
            and plancherel_J_eta1_wrong_uncentered(p, 1, 1) != j11
            and G1_closed_form_err(p) < 1e-8
            and G1_closed_form_wrong_drop_p_err(p) > 1e-3
        )
        j11_rows[str(p)] = {"J11_eta1": j11, "ok": row_ok}
        alg = alg and row_ok
    dim7 = coupled_profile_dim(7)
    dim4 = coupled_profile_dim(4)
    dim_ok = (
        dim7 == 22
        and dim4 == 7
        and coupled_profile_dim(3) == 4
        and coupled_profile_dim_wrong_drop_lambda0(7) == 21
        and coupled_profile_dim_wrong_drop_constants(7) == 15
        and coupled_profile_dim_wrong_drop_lambda0(7) != dim7
        and coupled_profile_dim_wrong_drop_constants(7) != dim7
    )
    I13 = jacobi_legendre_inertia(13, 5)
    I5 = jacobi_legendre_inertia(5, 1)
    Sp5 = S_p_from_G1(5)
    # η=1 PSD: J11>0 and J00=0 already; check a 3×3 leading principal
    p7_psd = all(
        plancherel_J_eta1(7, d, d) >= 0 for d in range(6)
    ) and plancherel_J_eta1(7, 0, 0) == 0
    cert = (
        I13["n_neg"] >= 1
        and I5["n_neg"] == 0
        and I5["n_pos"] == 1
        and I5["J11"] > 0
        and abs(I5["J11"] - 25 * Sp5) < 1e-6
        and p7_psd
        and eta_is_legendre(13)
        and not eta_is_trivial(13)
    )
    return {
        "proved": bool(alg and dim_ok and cert),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "unrestricted_weil_is_certificate": False,
        "p13_deg5_n_neg": I13["n_neg"],
        "p13_deg5_min_eig": I13["min_eig"],
        "p5_affine_n_neg": I5["n_neg"],
        "p5_affine_J11": I5["J11"],
        "p5_S_p": Sp5,
        "coupled_dim_k7": dim7,
        "coupled_dim_k4": dim4,
        "J11_eta1": j11_rows,
        "theorem": (
            "W=c^T J c on monomials, J_{de}(η)=∑_{λ≠0} η(λ)G_d conj(G_e).  "
            "η=1: J=pμ_{d+e}−μ_dμ_e PSD with J_{00}=0 (fail: drop centering; "
            "fail: drop p in G_1=p/(ω−1); fail: J_{η=1} has a negative pivot).  "
            "η=χ at p=13 deg≤5 is indefinite (fail: claim J_χ⪰0) so unrestricted "
            "Weil is not a certificate.  Affine p=5 recovers p² S_p.  "
            "dim V_{p,k}=k+(k−2)(k−1)/2 includes λ=0 (fail: drop leading; "
            "fail: drop constants).  Coupled LDL OPEN.  Not a k≥7 close."
        ),
    }


def psi_of_singer_eta(p: int) -> complex:
    """ψ(η) for η=γ^{p-1}, the Singer generator of C_m.  i^{p-1}."""
    return (1j) ** (p - 1)


def psi_of_singer_eta_wrong_as_fp_gen(p: int) -> complex:
    """Fail-eq: confuse with ψ(γ^{p+1})=i^{p+1}=−i^{p-1}."""
    return psi_on_fp_generator(p)


def singer_psi_constant_on_lines(p: int) -> bool:
    """ψ(ck) independent of Singer index  ⇔  ψ(η)=1  ⇔  p≡1 (mod 4)."""
    return psi_of_singer_eta(p) == 1


def top_leading_self_energy_vanishes(p: int) -> bool:
    """Top: κ=sign, κ²=1.  α=J_rr ∑ψ_j =0  ⇔  ∑ψ=0  ⇔  p≡3 (m even)."""
    return p % 4 == 3


def singer_vandermonde(p: int, k: int, d: int) -> list[list[int]]:
    """(d+1)×k matrix of (aX+bY)^d on the first k Singer forms."""
    from math import comb

    forms = singer_linear_forms(p)[:k]
    return [
        [
            (comb(d, i) * pow(a, d - i, p) * pow(b, i, p)) % p
            for a, b in forms
        ]
        for i in range(d + 1)
    ]


def degree_d_kernel_dim(p: int, k: int, d: int) -> int:
    """dim ker ∑ c_j t_j^d  (d≥1).  Equals k−d−1 when k≥d+1 and p>d."""
    M = singer_vandermonde(p, k, d)
    return k - matrix_rank_modp(M, p)


def coupled_Q_direct_p7(lam: int, clow: tuple) -> int:
    """∑_j w_j c^{(j)T} J c^{(j)} on p=7 top, c^{(j)}=(c0,c1,λ w_j)."""
    p, r = 7, 2
    J = [
        [plancherel_J_eta1(p, d, e) for e in range(r + 1)]
        for d in range(r + 1)
    ]
    w = singer_sign_vector(p)
    acc = 0
    for j in range(4):
        vec = [clow[j][0], clow[j][1], lam * w[j]]
        acc += w[j] * sum(
            vec[d] * J[d][e] * vec[e] for d in range(3) for e in range(3)
        )
    return acc


def coupled_Q_expansion_p7(lam: int, clow: tuple) -> int:
    """Q_low + 2λ ∑_j (J_{2,*}·c_low^{(j)})  (α=0, ψκ=1 on p=7 top)."""
    p = 7
    J20 = plancherel_J_eta1(p, 2, 0)
    J21 = plancherel_J_eta1(p, 2, 1)
    J22 = plancherel_J_eta1(p, 2, 2)
    w = singer_sign_vector(p)
    Qlow = coupled_Q_direct_p7(0, clow)
    cross = sum(J20 * clow[j][0] + J21 * clow[j][1] for j in range(4))
    alpha = J22 * sum(w)  # ∑ ψ κ² = ∑ w = 0
    return Qlow + 2 * lam * cross + lam * lam * alpha


def coupled_Q_expansion_wrong_drop_cross_p7(lam: int, clow: tuple) -> int:
    """Fail-eq: drop the 2λ interference."""
    return coupled_Q_direct_p7(0, clow)


def coupled_Q_expansion_wrong_unsigned_lead_p7(lam: int, clow: tuple) -> int:
    """Fail-eq: λ² m J_rr instead of α=0."""
    p = 7
    J22 = plancherel_J_eta1(p, 2, 2)
    return coupled_Q_direct_p7(0, clow) + lam * lam * 4 * J22


def theorem_AC_coupled_leading_pivot_not_a_cd_close() -> dict:
    """Coupled Q on V_{p,k}: leading CD pivot vanishes on p≡3 tops.

    PROVED (Max+-free algebra; integer expansion at p=7; no p=13 orbits)
      AC1. ψ(η)=i^{p-1} for η=γ^{p-1}.  Equals 1 iff p≡1 (mod 4),
           −1 iff p≡3.  Hence ψ(c_k)=ψ(α)ψ(η)^{-k} is constant on
           Singer lines iff p≡1, and equals the sign character iff
           p≡3.  Fail: confuse with i^{p+1}=ψ(γ^{p+1}); fail: claim
           ψ constant on lines at p=7.
      AC2. Degree-d kernel on k Singer forms has dim k−d−1 (d≥1,
           k≥d+1, p>d).  Summing d=1..k−2 plus k constants recovers
           dim V.  Certified p=7 k=4 and p=13 k=7.  Fail: dim k−d.
      AC3. Split c^{(j)}=(c_low^{(j)}, λ κ_j) with κ the degree-(k−2)
           kernel (top: κ=sign, λ=0 allowed).  Then
             Q=Q_<r + 2λ C(c_low) + λ² α,
           α=J_{rr}∑_j ψ_j κ_j².  On top p≡3: ψ=κ=sign, m even,
           ∑ψ=0 so α=0.  The last Christoffel–Darboux / Schur pivot
           of the degree filtration VANISHES.  Fail: claim α=m J_rr
           at p=7 (unsigned leading).  Fail: drop the 2λ cross
           (then Q would ignore λ; false for the p=7 integer sample).
      AC4. On top p≡1: ψ constant, κ²=1, α=m ψ_0 J_{rr} with m odd
           so α≠0 (p=13).  Fail: claim α=0 at p=13 (the p≡3 identity).
      AC5. p≡3 top therefore has Q=Q_low+2λ C: QVAR is the second
           moment of the leading–lower interference, which is why
           λ=0 cannot be split off and why equal-energy (C=0, α=0)
           gives Q=0.  Named quantity still E|Q|^2 on V_{p,k}.

    OPEN: a lower bound on E|Q_low+2λ C+λ² α|².  The vanishing
    last pivot kills a uniformly-positive-Schur-pivot CD certificate
    on every p≡3 top (the infinite family).  Not a k≥7 close.
    """
    eta_ok = True
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        eta_ok = eta_ok and (
            psi_of_singer_eta(p) == (1j) ** (p - 1)
            and singer_psi_constant_on_lines(p) == (p % 4 == 1)
            and top_leading_self_energy_vanishes(p) == (p % 4 == 3)
            and psi_of_singer_eta_wrong_as_fp_gen(p) != psi_of_singer_eta(p)
            and (p != 13 or not top_leading_self_energy_vanishes(p))
            and (p != 7 or top_leading_self_energy_vanishes(p))
            and (p != 7 or not singer_psi_constant_on_lines(p))
        )
    dim_ok = True
    for p, k in ((7, 4), (13, 7), (11, 4), (5, 3)):
        s = k
        for d in range(1, k - 1):
            dim_d = degree_d_kernel_dim(p, k, d)
            expect = k - d - 1
            s += dim_d
            dim_ok = dim_ok and dim_d == expect
            dim_ok = dim_ok and dim_d != kernel_dim_wrong_drop1(k, d)
        dim_ok = dim_ok and s == coupled_profile_dim(k)
    clow = ((1, 0), (0, 1), (2, -1), (-1, 3))
    exp_ok = True
    cross_seen = False
    for lam in (0, 1, -2, 3):
        direct = coupled_Q_direct_p7(lam, clow)
        expn = coupled_Q_expansion_p7(lam, clow)
        drop = coupled_Q_expansion_wrong_drop_cross_p7(lam, clow)
        unsigned = coupled_Q_expansion_wrong_unsigned_lead_p7(lam, clow)
        exp_ok = exp_ok and direct == expn
        if lam != 0:
            exp_ok = exp_ok and drop != direct
            exp_ok = exp_ok and unsigned != direct
            cross_seen = cross_seen or drop != direct
    exp_ok = exp_ok and cross_seen
    # α=0 at p=7: expansion's λ² coeff is ∑w=0
    w = singer_sign_vector(7)
    alpha0 = sum(w) == 0 and len(w) == 4
    return {
        "proved": bool(eta_ok and dim_ok and exp_ok and alpha0),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "cd_last_pivot_vanishes_p3_top": True,
        "p13_top_alpha_zero": False,
        "p7_expansion_ok": exp_ok,
        "coupled_dim_p13_k7": coupled_profile_dim(7),
        "theorem": (
            "Q=Q_<r+2λ C(c_low)+λ²α on V_{p,k}, λ=0 allowed.  "
            "ψ(η)=i^{p-1} (fail: i^{p+1}; fail: ψ constant at p=7).  "
            "Top p≡3: α=0, last CD pivot vanishes (fail: α=m J_rr).  "
            "Top p≡1: α=m ψ_0 J_rr≠0 (fail: α=0 at p=13).  "
            "Fail: drop 2λ (p=7 integer expansion).  Equal-energy top "
            "has Q=0 because α=C=0.  E|Q|² bound OPEN.  Not a k≥7 close."
        ),
    }


def ridge_qvar_floor_inner(p: int) -> int:
    """3 p S.  QVAR iff E[⟨y_fin,F⟩²] ≥ this, since ⟨y,F⟩=2Z and
    4·(3q(q−1)/16)=3pS."""
    return 3 * p * profile_energy_total_S(p)


def ridge_qvar_floor_inner_wrong_one(p: int) -> int:
    """Fail-eq: 3pS ↦ pS (spherical one-mode)."""
    return p * profile_energy_total_S(p)


def ridge_F_norm_sq(p: int) -> int:
    """||F||²=pS on every Max+ (p≡3)."""
    return p * profile_energy_total_S(p)


def ridge_F_norm_sq_wrong_drop_p(p: int) -> int:
    """Fail-eq: pS ↦ S."""
    return profile_energy_total_S(p)


def pythagoras_Q_sq(lam: int, clow: tuple) -> tuple[int, int, int, int]:
    """Per-config identity Q² = Q_low² + 4λ² C² + 4λ Q_low C on p=7 top model."""
    p = 7
    Q = coupled_Q_direct_p7(lam, clow)
    Qlow = coupled_Q_direct_p7(0, clow)
    J21 = plancherel_J_eta1(p, 2, 1)
    # J20=0; C = ∑_j J21 c1_j
    C = J21 * sum(clow[j][1] for j in range(4))
    return Q * Q, Qlow * Qlow, (lam * C) * (lam * C), lam * Qlow * C


def _p7_ridge_and_lambda_census() -> dict:
    """p=7 top: ridge pairing for actual Z, interpolant Pythagoras, E[Z²] vs |λ|."""
    import numpy as np
    from e1_gmin_m4_prop15588 import maxplus, profiles_of

    p = 7
    q = p * p
    S = profile_energy_total_S(p)
    Y = maxplus(p)
    P = profiles_of(p, Y)
    eps = Y[:, 0]
    act = (P != eps[:, None, None]).any(axis=2)
    top = act.sum(1) == 4
    Pt = P[top].astype(np.int64)
    ep = eps[top].astype(np.int64)
    Yf = Y[top, 1:].astype(np.int64)
    n = int(top.sum())
    perm = singer_direction_perm(p)
    w = np.array(singer_sign_vector(p), dtype=np.int64)
    from e1_gmin_m4_prop15588 import directions as _dirs

    sq, _nsq = _dirs(p)
    tarr = np.stack([np.asarray(sq[perm[j]][0], dtype=np.int64) for j in range(4)])
    h = (Pt[:, perm, :] - ep[:, None, None]) // 2
    a = (h.astype(np.int64) ** 2).sum(axis=2)
    Z = a @ w
    # F(x)=∑_j w_j h_j(t_j(x))
    F = np.zeros((n, q), dtype=np.int64)
    for j in range(4):
        F += w[j] * h[:, j, tarr[j]]
    yfin = Yf
    inner = (yfin * F).sum(axis=1)
    F2 = (F.astype(np.int64) ** 2).sum(axis=1)
    # interpolant leading λ (F_p, centered)
    rho = Pt[:, perm, :] % p
    lam = np.zeros(n, dtype=np.int64)
    s = np.arange(p)
    V = np.stack([s * 0 + 1, s, s * s], axis=1) % p
    for i in range(n):
        # interpolate first direction; λ = c2 / w_0, w_0=1
        # solve V c = rho[i,0] over F_p with 3 points s=0,1,2
        A = [[1, 0, 0], [1, 1, 1], [1, 2, 4]]
        b = [int(rho[i, 0, 0]) % p, int(rho[i, 0, 1]) % p, int(rho[i, 0, 2]) % p]
        # 3×3
        M = [A[r][:] + [b[r]] for r in range(3)]
        for c in range(3):
            piv = next(r for r in range(c, 3) if M[r][c] % p)
            M[c], M[piv] = M[piv], M[c]
            inv = pow(M[c][c] % p, p - 2, p)
            M[c] = [(v * inv) % p for v in M[c]]
            for r in range(3):
                if r == c:
                    continue
                fac = M[r][c] % p
                if fac:
                    M[r] = [(M[r][cc] - fac * M[c][cc]) % p for cc in range(4)]
        c2 = M[2][3] % p
        if c2 > p // 2:
            c2 -= p
        lam[i] = c2
    by_abs = {}
    for L in (1, 2, 3):
        msk = np.abs(lam) == L
        nn = int(msk.sum())
        eZ = Fraction(int((Z[msk].astype(object) ** 2).sum()), nn) if nn else Fraction(0)
        by_abs[L] = {"n": nn, "E_Z2": eZ}
    EZ = Fraction(int((Z.astype(object) ** 2).sum()), n)
    Einner = Fraction(int((inner.astype(object) ** 2).sum()), n)
    return {
        "n": n,
        "ridge_inner_eq_2Z": bool(np.all(inner == 2 * Z)),
        "F2_eq_pS": bool(np.all(F2 == p * S)),
        "E_Z2": EZ,
        "E_inner2": Einner,
        "four_EZ": 4 * EZ,
        "three_pS": 3 * p * S,
        "by_abs_lam": by_abs,
        "lam_abs_independent": all(by_abs[L]["E_Z2"] == EZ for L in (1, 2, 3)),
        "n_lam_nonzero": int((lam != 0).sum()),
    }


def theorem_AD_ridge_pythagoras_not_a_bound() -> dict:
    """Signed-ridge pairing and Pythagoras of Q_low+2λC.  Bound OPEN.

    PROVED (Max+-free algebra, p≡3; certified p=7 top)
      AD1. h_j=(σ_j−eps)/2, ∑_s h_j=0.  Distinct directions: t_j,t_ℓ
           independent so ∑_x h_j(t_j)h_ℓ(t_ℓ)=(∑h_j)(∑h_ℓ)=0.
           F=∑_j w_j h_j∘t_j has ||F||²=pS (ENERGY).  Fail: drop p.
      AD2. ε_H=∑_j h_j∘t_j=(p y_fin−eps)/2, ∑F=0, ⟨ε_H,F⟩=pZ, hence
           Z=⟨y_fin,F⟩/2.  Fail: drop the 2.
      AD3. QVAR ⇔ E[⟨y,F⟩²]≥3pS.  (4·floor_Z=3pS.)  Fail: 3pS↦pS
           (one spherical mode).  Certified p=7: inner=2Z, ||F||²=pS
           on all 8820 top configs; 4 E Z²=34496/15>1764=3pS.
      AD4. On V_{p,k} the expansion is an identity
           Q²=Q_low²+4λ²C²+4λ Q_low C.  Fail: drop a 4.
      AD5. Actual E[Z²] is independent of |λ| at p=7 top (8624/15
           on each |λ|∈{1,2,3}).  Fail: claim E[Z²]∝λ²
           (would be 1:4:9).  The unbounded polynomial model is
           not actual Z (wrap).  λ=0 cannot be split: Pythagoras
           remainder 4E[λ²C²] is the interference, not a p-law census.

    OPEN: E[⟨y,F⟩²]≥3pS (equivalently 4E[λ²C²]+E[Q_low²] after
    matching actual Z).  Not a k≥7 close.
    """
    clow = ((1, 0), (0, 1), (2, -1), (-1, 3))
    py_ok = True
    for lam in (0, 1, -2, 3):
        Q2, Ql2, l2C2, lQlC = pythagoras_Q_sq(lam, clow)
        py_ok = py_ok and Q2 == Ql2 + 4 * l2C2 + 4 * lQlC
        drop4 = Ql2 + l2C2 + 4 * lQlC
        if lam != 0:
            py_ok = py_ok and drop4 != Q2
    iff_ok = True
    for p in (7, 11, 19, 23, 31):
        S = profile_energy_total_S(p)
        four_floor = 4 * quartic_variance_floor_threshold(p)
        iff_ok = (
            iff_ok
            and four_floor == ridge_qvar_floor_inner(p)
            and ridge_qvar_floor_inner(p) == 3 * p * S
            and ridge_qvar_floor_inner_wrong_one(p) != 3 * p * S
            and ridge_F_norm_sq(p) == p * S
            and ridge_F_norm_sq_wrong_drop_p(p) != p * S
        )
    C7 = _p7_ridge_and_lambda_census()
    floor = quartic_variance_floor_threshold(7)
    cert = (
        C7["ridge_inner_eq_2Z"]
        and C7["F2_eq_pS"]
        and C7["E_inner2"] == 4 * C7["E_Z2"]
        and C7["four_EZ"] >= C7["three_pS"]
        and C7["E_Z2"] >= floor
        and C7["lam_abs_independent"]
        and C7["by_abs_lam"][1]["E_Z2"] != 9 * C7["by_abs_lam"][3]["E_Z2"]
        and C7["n_lam_nonzero"] == C7["n"]
        and py_ok
        and iff_ok
    )
    return {
        "proved": bool(cert),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "p7_inner_eq_2Z": C7["ridge_inner_eq_2Z"],
        "p7_F2_eq_pS": C7["F2_eq_pS"],
        "p7_E_Z2": str(C7["E_Z2"]),
        "p7_lam_independent": C7["lam_abs_independent"],
        "ridge_floor_3pS": True,
        "theorem": (
            "p≡3: Z=⟨y,F⟩/2 with ||F||²=pS (fail: drop 2; fail: drop p).  "
            "QVAR iff E[⟨y,F⟩²]≥3pS (fail: 3↦1).  Pythagoras "
            "Q²=Q_low²+4λ²C²+4λ Q_low C (fail: drop a 4).  Actual "
            "E[Z²] independent of |λ| at p=7 (fail: ∝λ²).  "
            "Bound OPEN.  Not a k≥7 close."
        ),
    }


# Exact depressed-quintic lift minima (normalized b=a/(2p)).  Same
# enumerator as evidence/k7_quintic_profile_probe.py; p=13,17,19 match
# that probe (1,3,4).  7*min>T empties k=7 (vacuous QVAR).
K7_QUINTIC_MIN_B = {
    13: 1,
    17: 3,
    19: 4,
    23: 6,
    29: 12,
    31: 10,
    37: 21,
    41: 30,
    43: 31,
    47: 38,
    53: 56,
    59: 68,
    61: 78,
    67: 92,
    71: 108,
    73: 116,
    79: 139,
    83: 155,
    89: 178,
}


def k7_seven_min_exceeds_T(p: int) -> bool:
    """True iff 7 * min_quintic_b > T, so k=7 is empty by energy."""
    if p not in K7_QUINTIC_MIN_B:
        return False
    return 7 * K7_QUINTIC_MIN_B[p] > T_of(p)


def k7_seven_min_exceeds_T_wrong_drop7to6(p: int) -> bool:
    """Fail-eq: 7 min ↦ 6 min (k=6 test)."""
    if p not in K7_QUINTIC_MIN_B:
        return False
    return 6 * K7_QUINTIC_MIN_B[p] > T_of(p)


def k7_energy_empty_primes() -> list[int]:
    return sorted(p for p in K7_QUINTIC_MIN_B if k7_seven_min_exceeds_T(p))


def k7_live_primes_in_table() -> list[int]:
    return sorted(p for p in K7_QUINTIC_MIN_B if not k7_seven_min_exceeds_T(p))


def theorem_AE_k7_energy_empty_pge53() -> dict:
    """k=7 is empty by exact quintic minima for every tabled p≥53.

    PROVED (exact lift enumeration; Max+-free arithmetic)
      AE1. A k=7 profile has reduced degree ≤5.  Nonzero top scalar
           plus translation (15.588 C+) gives the depressed quintic
           a s^5+c s^3+d s^2+e s+f, a≠0.  Fail: drop the depression
           (keep s^4) — not used as a QVAR bound.
      AE2. Zero-sum odd lifts of those reduced polynomials have a
           minimum normalized energy b_min(p).  If 7 b_min > T then
           the seven profiles cannot sum to T, so k=7 is empty and
           QVAR holds vacuously.  Fail: 7↦6 (at p=53, 6*56=336<351
           so the k=6 test would not empty).  Fail: claim p=13 empty
           (b_min=1, 7<T=21).  Fail: 7 b_min ≥ T at p=41 (7*30=T,
           equality is not emptiness).
      AE3. Exact b_min: p=13,17,19 are 1,3,4 (matches the existing
           quintic probe).  First emptiness is p=53 (56, 7*56=392>351).
           Tabled empty: 53,59,61,67,71,73,79,83,89.  Live in the table:
           13,17,19,23,29,31,37,41,43,47.  Weil still empties p>4*7^2
           =196.  Does not empty the top k=m for these p (m>7 except
           p=13).  Not a general k≥7 close.

    OPEN: k=7 at the ten live primes (ensemble, not pointwise);
    k≥8 including every top.  qvar stays False.
    """
    rows = {}
    ok = True
    for p, minimum in K7_QUINTIC_MIN_B.items():
        T = T_of(p)
        empty = 7 * minimum > T
        row_ok = (
            empty == k7_seven_min_exceeds_T(p)
            and (p != 13 or minimum == 1 and not empty)
            and (p != 41 or 7 * minimum == T and not empty)
            and (p != 53 or empty and 6 * minimum <= T)
        )
        rows[str(p)] = {
            "minimum_quintic_b": minimum,
            "T": T,
            "seven_min": 7 * minimum,
            "empty_by_energy": empty,
            "ok": row_ok,
        }
        ok = ok and row_ok
    empty_ps = k7_energy_empty_primes()
    live_ps = k7_live_primes_in_table()
    table_ok = (
        empty_ps == [53, 59, 61, 67, 71, 73, 79, 83, 89]
        and live_ps == [13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        and K7_QUINTIC_MIN_B[13] == 1
        and K7_QUINTIC_MIN_B[17] == 3
        and K7_QUINTIC_MIN_B[19] == 4
        and not weil_vacuous_qvar_k_ge_7(53, 7)
        and weil_vacuous_qvar_k_ge_7(197, 7)
        and 6 * K7_QUINTIC_MIN_B[53] <= T_of(53) < 7 * K7_QUINTIC_MIN_B[53]
    )
    return {
        "proved": bool(ok and table_ok),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "k7_empty_p": empty_ps,
        "k7_live_p_in_table": live_ps,
        "by_p": rows,
        "theorem": (
            "k=7 depressed-quintic b_min: 7 b_min>T empties the stratum "
            "(fail: 7↦6 at p=53; fail: empty at p=13; fail: ≥ instead of "
            "> at p=41 where 7*30=T).  Exact empty "
            "p=53,59,61,67,71,73,79,83,89.  Live p=13..47 in the table; "
            "Weil p>196.  Not a top close and not a k≥7 close."
        ),
    }


def theorem_AF_p41_k7_stratum_fails_qvar() -> dict:
    """p=41 k=7 census: nonempty, a_L∈2pℤ, E|Z_ψ|² below the floor.

    Computed by evidence/k7_p41_coefficient_sieve.scan_p41_k7 (coupled
    cubic/quadratic/linear sieve, Boolean endpoint, quartic-kernel Z_ψ
    audited against ∑_{d≠0} ψ(d)N(d)).  p=41≡1 (mod 4), so the p≡3
    signed-sum ∑ψ(g_L)a_L is not Z_ψ.
    """
    sys.path.insert(0, str(ROOT / "evidence"))
    from k7_p41_coefficient_sieve import scan_p41_k7

    report = scan_p41_k7()
    q = 41 * 41
    threshold = 3 * q * (q - 1) // 16
    energy = report["E_abs_Zpsi_sq"]
    n_bool = report["boolean_representatives_mod_translation"]
    empty_wrong = 7 * K7_QUINTIC_MIN_B[41] > T_of(41)
    ok = (
        report["p"] == 41
        and report["normalized_total_T"] == T_of(41)
        and 7 * K7_QUINTIC_MIN_B[41] == T_of(41)
        and not empty_wrong
        and report["k7_empty"] is False
        and n_bool > 0
        and energy < threshold
        and threshold == report["QVAR_threshold"]
        and threshold != 3 * q * (q - 1) // 8
        and report["official_Z_matches_kernel"] is True
        and report["boolean_a_L_all_match"] is True
        and report["a_L_expected"] == 2 * 41 * 30
        and report["p_mod_4"] == 1
        and report["integer_signed_sum_is_Z"] is False
        and report["stratum_qvar"] is False
        and report.get("maxplus_Cy_eq_py") is True
    )
    return {
        "proved": bool(ok),
        "covers_general_k_ge_7": False,
        "lambda_min_bound_proved": False,
        "k7_empty": bool(report["k7_empty"]),
        "stratum_qvar": bool(report["stratum_qvar"]),
        "boolean_mod_translation": int(n_bool),
        "E_abs_Zpsi_sq": energy,
        "QVAR_threshold": int(threshold),
        "official_Z_matches_kernel": bool(report["official_Z_matches_kernel"]),
        "a_L_in_2pZ": bool(report["boolean_a_L_all_match"]),
        "maxplus_Cy_eq_py": bool(report.get("maxplus_Cy_eq_py")),
        "theorem": (
            "p=41 k=7 unique partition 7×30=T (fail: >).  Coupled sieve "
            "+ Boolean: nonempty, a_L=2p·30∈2pℤ, official Z_ψ matches "
            "the quartic kernel, E|Z_ψ|² < 3q(q−1)/16 (fail: empty; "
            "fail: E≥floor; fail: drop 16).  p≡1 so ∑ψ(g_L)a_L is not "
            "Z_ψ.  Per-stratum QVAR false at this pair.  Not a k≥7 close."
        ),
    }


def qvar_k_ge_7_proved_general() -> bool:
    """Universal k≥7 QVAR, imported from the (41,7) census and AE cover.

    False when the computed (41,7) stratum misses the floor, and false
    while AE does not cover tops / remaining live k=7.  No handwritten
    True or False.
    """
    AF = theorem_AF_p41_k7_stratum_fails_qvar()
    AE = theorem_AE_k7_energy_empty_pge53()
    return bool(AF["stratum_qvar"]) and bool(AE["covers_general_k_ge_7"])


def live_L_status() -> str:
    """Global limit status is independent of this optional Paley route."""
    from original_mo_status import original_mo_status

    return original_mo_status()["limit_status"]


def main() -> dict:
    from io_atomic import write_json_atomic

    A = theorem_A_floors()
    C = theorem_C_weil_vacuous_range()
    D = theorem_D_pointwise_orbit_counterexample()
    G = theorem_G_kernel_ladder()
    H = theorem_H_zero_top_scalar_weil()
    I = theorem_I_top_stratum_not_exactly_spherical()
    J = theorem_J_bochner_too_weak_for_qvar()
    K = theorem_K_harm4_not_one_dimensional()
    Lrk = theorem_L_A_psi_not_rank_one()
    N = theorem_N_lambda_mass_not_uniform()
    O = theorem_O_cs_harm4_pairing_too_weak()
    P = theorem_P_top_stratum_plus_energy_identity()
    Qeq = theorem_Q_equal_energy_not_a_lower_bound()
    R = theorem_R_two_type_opposite_pair_not_unique()
    S = theorem_S_nonzero_top_scalar_not_pointwise_qvar()
    T = theorem_T_two_design_does_not_determine_profile_energy_variance()
    U = theorem_U_second_moment_rayleigh_not_a_close()
    V = theorem_V_singer_circulant_nyquist()
    W = theorem_W_top_kernel_is_alternating()
    X = theorem_X_top_leading_is_sign_isotypic()
    Y = theorem_Y_psi_mode_of_Phi_not_a_close()
    Z3 = theorem_Z_three_level_not_a_lower_bound()
    AA = theorem_AA_line_dft_both_congruences()
    AB = theorem_AB_jacobi_gram_not_unrestricted_weil()
    AC = theorem_AC_coupled_leading_pivot_not_a_cd_close()
    AD = theorem_AD_ridge_pythagoras_not_a_bound()
    AE = theorem_AE_k7_energy_empty_pge53()
    AF = theorem_AF_p41_k7_stratum_fails_qvar()
    out = {
        "title": "QVAR k>=7 (not a close)",
        "numbered": False,
        "A_floors": A,
        "C_weil_vacuous": C,
        "D_p13_k7_counterexample": D,
        "G_kernel_ladder": G,
        "H_zero_top_scalar_weil": H,
        "I_top_not_exactly_spherical": I,
        "J_bochner_too_weak": J,
        "K_harm4_not_1dim": K,
        "L_A_psi_not_rank1": Lrk,
        "N_lambda_mass_not_uniform": N,
        "O_cs_harm4_too_weak": O,
        "P_top_plus_energy_identity": P,
        "Q_equal_energy_not_lower_bound": Qeq,
        "R_two_type_not_unique": R,
        "S_lambda_nonzero_not_pointwise": S,
        "T_two_design_not_Var_a": T,
        "U_second_moment_rayleigh": U,
        "V_singer_circulant_nyquist": V,
        "W_top_kernel_alternating": W,
        "X_top_leading_sign_isotypic": X,
        "Y_psi_mode_of_Phi": Y,
        "Z_three_level_not_lb": Z3,
        "AA_line_dft_both_congruences": AA,
        "AB_jacobi_gram": AB,
        "AC_coupled_leading_pivot": AC,
        "AD_ridge_pythagoras": AD,
        "AE_k7_energy_empty": AE,
        "AF_p41_k7_stratum_fails_qvar": AF,
        "qvar_k_ge_7_proved_general": qvar_k_ge_7_proved_general(),
        "p13_orbits_not_a_close": True,
        "L_status": live_L_status(),
    }
    path = ROOT / "evidence" / "e1_gmin_qvar_k_ge_7.json"
    write_json_atomic(path, out)
    print("QVAR k>=7 (unnumbered, not a close)", flush=True)
    print(f"  A floors: {A['proved']}", flush=True)
    print(f"  C Weil vacuous range: {C['proved']}", flush=True)
    print(f"  D p=13 k=7 pointwise false: {D['proved']}", flush=True)
    print(f"  G kernel ladder: {G['proved']}", flush=True)
    print(f"  H zero-top-scalar Weil: {H['proved']}", flush=True)
    print(f"  I top ≠ V_sph (p=11 k=6): {I['proved']}", flush=True)
    print(f"  J Bochner too weak: {J['proved']}", flush=True)
    print(f"  K Harm4^G not 1-dim: {K['proved']}", flush=True)
    print(f"  L A_ψ not rank-1: {Lrk['proved']}", flush=True)
    print(f"  N λ-mass not uniform: {N['proved']}", flush=True)
    print(f"  O CS Harm4 too weak: {O['proved']}", flush=True)
    print(f"  P top plus-energy identity: {P['proved']}", flush=True)
    print(f"  Q equal-energy not a lower bound: {Qeq['proved']}", flush=True)
    print(f"  R two-type uniqueness killed: {R['proved']}", flush=True)
    print(f"  S λ≠0 not pointwise QVAR: {S['proved']}", flush=True)
    print(f"  T 2-design does not determine Var(a): {T['proved']}", flush=True)
    print(f"  U second-moment Rayleigh (not a close): {U['proved']}", flush=True)
    print(f"  V Singer circulant / Nyquist iff (bound open): {V['proved']}", flush=True)
    print(f"  W top kernel is alternating (bound open): {W['proved']}", flush=True)
    print(f"  X top leading ∈ span{{w}} (bound open): {X['proved']}", flush=True)
    print(f"  Y ψ-mode of Φ (bound open): {Y['proved']}", flush=True)
    print(f"  Z 3-level not a lower bound: {Z3['proved']}", flush=True)
    print(f"  AA line-DFT both congruences (bound open): {AA['proved']}", flush=True)
    print(f"  AB Jacobi Gram (unrestricted Weil killed): {AB['proved']}", flush=True)
    print(f"  AC coupled leading pivot (CD last pivot vanishes p≡3 top): {AC['proved']}", flush=True)
    print(f"  AD ridge Pythagoras (bound open): {AD['proved']}", flush=True)
    print(f"  AE k=7 energy-empty p>=53 (not a k>=7 close): {AE['proved']}", flush=True)
    print(f"  AF p=41 k=7 stratum E|Z|^2=0 (QVAR fails here): {AF['proved']}", flush=True)
    print(f"  qvar_k_ge_7_proved_general: {qvar_k_ge_7_proved_general()}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
