#!/usr/bin/env python3
"""
R1 principal L² room — ‖P_{E_{4p}} m₄⁺‖² ≤ n(λ̄−6)²/48  (→ n/12).

Unnumbered.  Does **not** claim the bound.  Does **not** flip leftover
flags, QVAR k≥7, phi_F, L, Aut-Schur, Gsum, or pairing.

WHAT IS PROVED Max+-free
  A. The 15.595 leftover-1 threshold is the exact Fraction
        R1(p) = n(λ̄−6)²/48,  λ̄=8(n−2)/(n−6),  n=p²+1
     and R1(p)/n → 1/12.  Fail: replace 48 by 24; claim R1 = n/12 exactly
     at p=5.
  B. Measured ‖δ‖² from 15.594 (exact, Max+ enum):
        p=5:  1536/65
        p=7:  19180800/1840091
        p=11: 1382747375360/583792784981
     vs R1(p): EXCEEDS at p=5 and p=7; HOLDS at p=11.  Fail: claim the
     p=5 (resp. p=7) measurement is ≤ R1.  p=11 is census, not a p-law.
  C. 2-design + sharp min-distance cannot prove R1.  Max+-free:
        E s² = 2n (15.593),  |s| ≤ A = (p−1)²−2 for y≠±z (15.197, sharp
        at p=5,7).  The moment-problem majorant E s⁴ ≤ 2n A² is Θ(n³) and
        strictly exceeds the leftover-1 Es4 budget (design floor +
        (n/2)(λ̄−6)²) at every prime p≥5.  Fail: claim 2n A² ≤ 12 n²
        at p=11.

HECKE / PSL (15.589 A–C, 15.593 D)
  D. Primitive idempotents of End_G(Z) are the projectors P_c onto
        W_e (dim n/2) and the (q−9)/8 principal-series ρ(α) (dim n).
     Φ = ∑_c λ_c P_c,  V = ∑_c dim(c)(λ_c−λ̄)² = 24‖δ‖², and
        D := n(λ̄−6)² − 48‖δ‖²
          = n(λ̄−6)² − 2V
          = 2 ∑_c dim(c) (τ − (λ_c−λ̄)²)
     with τ = 16(n+10)²/(n−6)³.  Equivalently
        D = n(λ̄−6)² − 2 ∑_c ‖P_c(Φ−λ̄ I)‖_F².
     The coefficient of every isotypic deviation norm is −2, not ≥0.
     Fail: replace −2 by +2 (would force D>0 even at p=5); claim those
     Hecke coefficients are nonnegative.  Certified vs 15.586 spectra
     at p=5,7 and vs measured ‖δ‖² at p=5,7,11.
  E. A dual SOS in End_G(Z) from traces + multiplicity-freeness alone
     cannot prove D≥0: on the hyperplane ∑ dim(c)(λ_c−λ̄)=0 the form
     D = const − 2∑ dim(c) μ_c² is negative definite in the μ's.
     Forbidden-configuration slacks (y_i²−1, (C−pI)y) vanish on Max+
     and do not cancel the −2.  Hence
        hecke_dual_nonneg_coeffs_p_ge_11() is False,
     and r1_l2_bound_for_p_ge_11() stays False.  p=11 D>0 is census.

KNOWN MAJORANTS (proved too weak for n/12 at every p≥5, including p≥11)
  F. 15.100/15.102  ‖δ‖² ≤ room/24 = 4n(p²−9)/(3(p²−5))  (⇔ ‖κ‖²≤96n)
     and the tighter κ_hyp form room/24·((d−1)/d)² both exceed n/12
     at every prime p≥5.  Fail: claim room/24 ≤ n/12 at p=11.
     Combined with C, no shipped Max+-free majorant proves R1 for p≥11.

DEAD INTERPOLANT (killed as a retained δ-bound; not imported)
  G. B(p) := κ_hyp_δ(p)·4/(p−3)² still has Max+-free algebra B≤n/12 for
     p≥11 (cubic h), and B(5)=measured.  That does **not** make B a bound
     on δ.  Kill:
       (i)  the *equality* law δ² = κ_hyp_δ·4/(p−3)² is false at p=7 and
            p=11 (exact Fractions).  Fail: claim equality at p=7.
       (ii) the Aut-dim interpolant κ_hyp_δ·(ν_G(5)/ν_G(p))² with
            ν_G(5)=2, ν_G(7)=7 (15.134 C) is false at p=7: measured
            19180800/1840091 > 12288/2695.  Fail: claim that p=7
            measurement is ≤ the ν_G-ratio majorant.
       (iii) no Max+-free operator identity produces 4/(p−3)²
            (Hecke dual −2; 2-design+distance Θ(n³); room/24 and κ_hyp
            are Θ(n) > B for p>5).
     So 4/(p−3)² is **not retained as a bound on ‖δ‖²**.  r1_l2 stays
     False.  p=5,7 still exceed n/12 (census).

principal_delta_room_moment_proved imports r1_l2 (so it is True only via
this unit).  qvar_k_ge_7 stays False.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15100 import d_of, room_orth  # noqa: E402
from e1_gmin_m4_prop15589 import dim_Z, n_principal_constituents  # noqa: E402
from e1_gmin_m4_prop15593 import design_floor, lambda_bar, threshold_leftover1  # noqa: E402
from e1_gmin_m4_prop15595 import req_leftover1  # noqa: E402


# Exact ‖δ‖² = V/24 from 15.594 (p=5,7 Max+ enum; p=11 from 15.593 V).
MEASURED_DELTA_SQ = {
    5: Fraction(1536, 65),
    7: Fraction(19180800, 1840091),
    11: Fraction(1382747375360, 583792784981),
}


def r1_threshold(p: int) -> Fraction:
    """‖δ‖² ≤ n(λ̄−6)²/48.  Same as 15.595 req_leftover1 and V/24 of 15.593 D."""
    return req_leftover1(p)


def r1_threshold_wrong_drop48(p: int) -> Fraction:
    """Fail-eq: 48 ↦ 24."""
    n = n_of(p)
    lbar = Fraction(8 * (n - 2), n - 6)
    return Fraction(n, 24) * (lbar - 6) ** 2


def r1_n_over_12_not_exact(p: int) -> Fraction:
    """Fail-eq: claim the finite-p leftover-1 threshold is exactly n/12."""
    return Fraction(n_of(p), 12)


def leftover1_minus_n_over_12_per_n(p: int) -> Fraction:
    """Exact gap: leftover-1 threshold/n − 1/12 = 8(n+2)/(3(n-6)²) > 0.

    So n/12 is strictly stronger than leftover 1 at every finite n>6.
    Fail: drop the 8; drop the 3 in the denominator.
    """
    n = n_of(p)
    return Fraction(8 * (n + 2), 3 * (n - 6) ** 2)


def leftover1_minus_n_over_12_per_n_wrong_drop8(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(n + 2, 3 * (n - 6) ** 2)


def leftover1_minus_n_over_12_per_n_wrong_drop3(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n + 2), (n - 6) ** 2)


def min_distance_A(p: int) -> int:
    """(p−1)² − 2 = n − 2(p+1).  Sharp at p=5,7 (15.197 D)."""
    return (p - 1) ** 2 - 2


def moment_problem_Es4_majorant(p: int) -> Fraction:
    """E s⁴ ≤ 2n A² from E s²=2n and |s|≤A (ignoring ±n atoms, which only
    increase a true upper bound).  Θ(n³)."""
    n = n_of(p)
    A = min_distance_A(p)
    return Fraction(2 * n * A * A)


def leftover1_Es4_budget(p: int) -> Fraction:
    """Design floor + leftover-1 variance room: Es4 ≤ floor + (n/2)(λ̄−6)²."""
    return design_floor(p) + threshold_leftover1(p)


def measured_meets_r1(p: int) -> bool:
    if p not in MEASURED_DELTA_SQ:
        raise KeyError(p)
    return MEASURED_DELTA_SQ[p] <= r1_threshold(p)


def r1_l2_bound_for_p_ge_11() -> bool:
    """True only if a proved δ-bound is imported. The 4/(p−3)² interpolant
    is killed as a retained bound (theorem H)."""
    return bool(interpolant_retained_as_delta_bound())


# 15.586 E exact Φ-spectra (multiplicity-free PSL isotypes).
SPECTRUM_P5: dict[Fraction, int] = {
    Fraction(80, 13): 26,
    Fraction(144, 13): 26,
    Fraction(176, 13): 13,
}
SPECTRUM_P7: dict[Fraction, int] = {
    Fraction(3072, 409): 50,
    Fraction(3360, 409): 100,
    Fraction(3648, 409): 50,
    Fraction(4032, 409): 50,
    Fraction(4320, 409): 25,
}


def hecke_tau(p: int) -> Fraction:
    """Var threshold τ with D = 2 ∑ dim(c)(τ − (λ_c−λ̄)²)."""
    n = n_of(p)
    return Fraction(16 * (n + 10) ** 2, (n - 6) ** 3)


def hecke_tau_wrong_drop16(p: int) -> Fraction:
    """Fail-eq: 16 ↦ 8."""
    n = n_of(p)
    return Fraction(8 * (n + 10) ** 2, (n - 6) ** 3)


def hecke_deviation_coeff() -> int:
    """Coefficient of each ‖P_c(Φ−λ̄ I)‖_F² in D.  Exact: −2."""
    return -2


def hecke_deviation_coeff_wrong_sign() -> int:
    """Fail-eq: −2 ↦ +2, which would make D = n(λ̄−6)² + 2V > 0 at p=5."""
    return +2


def V_from_spectrum(evals: dict[Fraction, int]) -> Fraction:
    tr = sum(lam * m for lam, m in evals.items())
    dim = sum(evals.values())
    lbar = Fraction(tr, dim)
    return sum(m * (lam - lbar) ** 2 for lam, m in evals.items())


def D_from_delta(p: int) -> Fraction:
    """n(λ̄−6)² − 48‖δ‖² from the measured δ of 15.594."""
    return Fraction(n_of(p)) * (lambda_bar(p) - 6) ** 2 - 48 * MEASURED_DELTA_SQ[p]


def D_from_spectrum(evals: dict[Fraction, int], p: int) -> Fraction:
    """n(λ̄−6)² − 2V from an exact Φ-spectrum."""
    V = V_from_spectrum(evals)
    return Fraction(n_of(p)) * (lambda_bar(p) - 6) ** 2 - 2 * V


def D_hecke_expansion(evals: dict[Fraction, int], p: int) -> Fraction:
    """2 ∑ dim(c)(τ − (λ_c−λ̄)²)."""
    tr = sum(lam * m for lam, m in evals.items())
    dim = sum(evals.values())
    lbar = Fraction(tr, dim)
    tau = hecke_tau(p)
    return 2 * sum(m * (tau - (lam - lbar) ** 2) for lam, m in evals.items())


def D_wrong_sign_expansion(evals: dict[Fraction, int], p: int) -> Fraction:
    """Fail-eq +2 on V: n(λ̄−6)² + 2V."""
    V = V_from_spectrum(evals)
    return Fraction(n_of(p)) * (lambda_bar(p) - 6) ** 2 + 2 * V


def hecke_dual_nonneg_coeffs_p_ge_11() -> bool:
    """A dual SOS in End_G(Z) with nonnegative isotypic coeffs does not exist:
    the identity D = n(λ̄−6)² − 2∑‖P_c(Φ−λ̄ I)‖_F² has coeff −2."""
    return False


def room_orth_delta_sq(p: int) -> Fraction:
    """15.102: ‖κ‖²≤96n ⇔ ‖δ‖² ≤ room/24 = 4n(p²−9)/(3(p²−5))."""
    return room_orth(p) / 24


def kappa_hyp_delta_sq(p: int) -> Fraction:
    """15.102 κ_hyp form: room/24 · ((d−1)/d)² with d=n/2."""
    d = d_of(p)
    return room_orth_delta_sq(p) * Fraction((d - 1) ** 2, d * d)


def room_orth_delta_sq_wrong_drop24(p: int) -> Fraction:
    """Fail-eq: /24 ↦ /12, which would look like n/12-scale at small p."""
    return room_orth(p) / 12


def interpolating_majorant(p: int) -> Fraction:
    """B(p) = κ_hyp_δ(p) · 4/(p−3)².  p>3."""
    if p <= 3:
        raise ValueError(p)
    return kappa_hyp_delta_sq(p) * 4 / (p - 3) ** 2


def interpolating_majorant_wrong_drop4(p: int) -> Fraction:
    """Fail-eq: 4 ↦ 2."""
    return kappa_hyp_delta_sq(p) * 2 / (p - 3) ** 2


def interpolating_majorant_wrong_p_minus_1(p: int) -> Fraction:
    """Fail-eq: (p−3) ↦ (p−1)."""
    return kappa_hyp_delta_sq(p) * 4 / (p - 1) ** 2


def n12_minus_B_sufficient_cubic(p: int) -> int:
    """h(p)=p³−3p²−69p−177.  h≥0 ⇒ B≤n/12 (p>3)."""
    return p ** 3 - 3 * p * p - 69 * p - 177


def n12_minus_B_sufficient_cubic_wrong_drop69(p: int) -> int:
    return p ** 3 - 3 * p * p - p - 177


def n12_minus_B_sufficient_cubic_wrong_drop3(p: int) -> int:
    return p ** 3 - 69 * p - 177


def interpolating_majorant_le_n_over_12_for_p_ge_11() -> bool:
    """Max+-free: h(11)=32>0 and h'(p)=3(p²−2p−23)>0 for p≥7, so h>0
    for every prime p≥11, and h≥0 is sufficient for B≤n/12."""
    h11 = n12_minus_B_sufficient_cubic(11)
    if h11 != 32:
        return False
    # h'(p)=3(p²−2p−23)=3((p−1)²−24); (p−1)²≥36 for p≥7 ⇒ h'>0
    if 3 * ((7 - 1) ** 2 - 24) <= 0:
        return False
    for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        if not is_prime(p):
            continue
        if n12_minus_B_sufficient_cubic(p) <= 0:
            return False
        if interpolating_majorant(p) > Fraction(n_of(p), 12):
            return False
        if interpolating_majorant_wrong_drop4(p) == interpolating_majorant(p):
            return False
        if interpolating_majorant_wrong_p_minus_1(p) == interpolating_majorant(p):
            return False
    # p=5: B exceeds n/12 (must not claim the algebra for p=5)
    if interpolating_majorant(5) <= Fraction(n_of(5), 12):
        return False
    return True


def delta_sq_le_interpolating_majorant() -> bool:
    """Not a proved bound. Equality law already false at p=7,11."""
    return False


# 15.134 C: dim E_{4p}^{G} on strict Aut orbit-space.
NU_G = {3: 0, 5: 2, 7: 7}


def nu_G_ratio_sq_majorant(p: int) -> Fraction:
    """κ_hyp_δ · (ν_G(5)/ν_G(p))².  Interpolates p=5 equality; dies at p=7."""
    if p not in NU_G or NU_G[p] == 0:
        raise KeyError(p)
    return kappa_hyp_delta_sq(p) * Fraction(NU_G[5] ** 2, NU_G[p] ** 2)


def nu_G_ratio_sq_majorant_wrong_nu5(p: int) -> Fraction:
    """Fail-eq: ν_G(5)↦1."""
    return kappa_hyp_delta_sq(p) * Fraction(1, NU_G[p] ** 2)


def interpolant_equality_law_fails_p7_p11() -> bool:
    """Tight law δ² = κ_hyp_δ·4/(p−3)² is false at p=7 and p=11."""
    for p in (7, 11):
        lhs = MEASURED_DELTA_SQ[p] * (p - 3) ** 2
        rhs = 4 * kappa_hyp_delta_sq(p)
        if lhs == rhs:
            return False
    return True


def interpolant_equality_law_wrong_claim_p7() -> bool:
    """Fail-eq: claim the p=7 measurement saturates B."""
    p = 7
    return MEASURED_DELTA_SQ[p] * (p - 3) ** 2 == 4 * kappa_hyp_delta_sq(p)


def nu_G_ratio_sq_killed_at_p7() -> bool:
    """measured_p7 > κ_hyp_δ(7)·(2/7)²."""
    return MEASURED_DELTA_SQ[7] > nu_G_ratio_sq_majorant(7)


def interpolant_retained_as_delta_bound() -> bool:
    """Killed: not retained as a bound on ‖δ‖²."""
    return False


def theorem_A_threshold_formula(primes=(5, 7, 11, 13, 17, 19, 23, 47)) -> dict:
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        r = r1_threshold(p)
        n = n_of(p)
        gap = leftover1_minus_n_over_12_per_n(p)
        row_ok = (
            r == req_leftover1(p)
            and r1_threshold_wrong_drop48(p) != r
            and r1_n_over_12_not_exact(p) != r
            and r / n - Fraction(1, 12) == gap
            and gap > 0
            and leftover1_minus_n_over_12_per_n_wrong_drop8(p) != gap
            and leftover1_minus_n_over_12_per_n_wrong_drop3(p) != gap
        )
        rows[str(p)] = {
            "R1": str(r),
            "R1_per_n": str(r / n),
            "n_over_12": str(Fraction(n, 12)),
            "gap_R1_minus_n12_per_n": str(gap),
            "ok": row_ok,
        }
        ok = ok and row_ok
    p5_exact_n12 = r1_threshold(5) == Fraction(n_of(5), 12)
    return {
        "proved": bool(ok and not p5_exact_n12),
        "theorem": (
            "Leftover-1 threshold R1(p)=n(λ̄−6)²/48, λ̄=8(n−2)/(n−6).  "
            "R1(p)/n − 1/12 = 8(n+2)/(3(n-6)²)>0, so n/12 is strictly "
            "stronger at every finite n>6 (limit 1/12).  Fail: 48↦24; "
            "claim R1=n/12 exactly at p=5; drop 8 or 3 in the gap."
        ),
        "by_p": rows,
        "limit_1_over_12": True,
        "exact_n_over_12_at_p5": p5_exact_n12,
        "n_over_12_strictly_below_leftover1": True,
    }


def measured_delta_sq_per_n(p: int) -> Fraction:
    return MEASURED_DELTA_SQ[p] / n_of(p)


def theorem_B_measured_vs_r1() -> dict:
    rows = {}
    for p, dsq in MEASURED_DELTA_SQ.items():
        n = n_of(p)
        thr = r1_threshold(p)
        per = dsq / n
        n12 = Fraction(1, 12)
        meets = dsq <= thr
        exceeds_n12 = per > n12
        rows[str(p)] = {
            "delta_sq": str(dsq),
            "R1": str(thr),
            "delta_sq_per_n": str(per),
            "R1_per_n": str(thr / n),
            "n_over_12": str(n12),
            "meets_R1": meets,
            "exceeds_n_over_12": exceeds_n12,
        }
    # fail-eq: claiming p=5 or p=7 meets R1, or claiming n/12 at p=5,7
    fake5 = MEASURED_DELTA_SQ[5] <= r1_threshold(5)
    fake7 = MEASURED_DELTA_SQ[7] <= r1_threshold(7)
    fake5_n12 = measured_delta_sq_per_n(5) <= Fraction(1, 12)
    fake7_n12 = measured_delta_sq_per_n(7) <= Fraction(1, 12)
    return {
        "proved": bool(
            rows["5"]["meets_R1"] is False
            and rows["7"]["meets_R1"] is False
            and rows["11"]["meets_R1"] is True
            and rows["5"]["exceeds_n_over_12"] is True
            and rows["7"]["exceeds_n_over_12"] is True
            and rows["11"]["exceeds_n_over_12"] is False
            and fake5 is False
            and fake7 is False
            and fake5_n12 is False
            and fake7_n12 is False
        ),
        "theorem": (
            "Exact ‖δ‖²/n exceeds both R1/n and 1/12 at p=5 and p=7; "
            "holds vs both at p=11.  Fail: claim the p=5 or p=7 measurement "
            "is ≤ R1 or ≤ n/12.  p=11 is census, not a proof for p≥11.  "
            "A uniform n/12 claim for all p≥5 is false."
        ),
        "by_p": rows,
        "p11_is_census_not_plaw": True,
        "n_over_12_not_claimed_at_p5_p7": True,
    }


def theorem_C_moment_problem_too_weak(
    primes=(5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47),
) -> dict:
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        maj = moment_problem_Es4_majorant(p)
        bud = leftover1_Es4_budget(p)
        n = n_of(p)
        too_weak = maj > bud
        not_12n2 = maj != Fraction(12 * n * n)
        row_ok = too_weak and not_12n2
        rows[str(p)] = {
            "majorant_2n_A2": str(maj),
            "leftover1_Es4_budget": str(bud),
            "too_weak": too_weak,
            "ok": row_ok,
        }
        ok = ok and row_ok
    p11_fake = moment_problem_Es4_majorant(11) <= Fraction(12 * n_of(11) ** 2)
    return {
        "proved": bool(ok and not p11_fake),
        "theorem": (
            "E s²=2n and |s|≤A=(p−1)²−2 (sharp at p=5,7) give "
            "E s⁴ ≤ 2n A² = Θ(n³), which exceeds the leftover-1 Es4 "
            "budget at every prime p≥5.  Fail: claim 2n A² ≤ 12 n² at p=11."
        ),
        "by_p": rows,
        "A_sharp_at_p5_p7": True,
    }


def theorem_D_hecke_gap_expansion() -> dict:
    """PSL primitive-idempotent expansion of D; dual SOS coeffs are −2."""
    n5, n7 = n_of(5), n_of(7)
    ok_dim = (
        dim_Z(5) == 65
        and dim_Z(7) == 275
        and n_principal_constituents(5) == 2
        and n_principal_constituents(7) == 5
        and n_principal_constituents(11) == 14
        and sum(SPECTRUM_P5.values()) == dim_Z(5)
        and sum(SPECTRUM_P7.values()) == dim_Z(7)
        and n_of(5) // 2 in SPECTRUM_P5.values()
        and n_of(7) // 2 in SPECTRUM_P7.values()
    )
    D5_delta = D_from_delta(5)
    D7_delta = D_from_delta(7)
    D11_delta = D_from_delta(11)
    D5_spec = D_from_spectrum(SPECTRUM_P5, 5)
    D7_spec = D_from_spectrum(SPECTRUM_P7, 7)
    D5_hecke = D_hecke_expansion(SPECTRUM_P5, 5)
    D7_hecke = D_hecke_expansion(SPECTRUM_P7, 7)
    ident5 = D5_delta == D5_spec == D5_hecke
    ident7 = D7_delta == D7_spec == D7_hecke
    # 2 τ dim Z = n(λ̄−6)²
    tau_id = True
    for p in (5, 7, 11, 13, 17, 23, 47):
        n = n_of(p)
        lhs = 2 * hecke_tau(p) * dim_Z(p)
        rhs = n * (lambda_bar(p) - 6) ** 2
        tau_id = tau_id and lhs == rhs and hecke_tau_wrong_drop16(p) != hecke_tau(p)
    # fail-eq: +2 on V makes D>0 at p=5
    fake_pos5 = D_wrong_sign_expansion(SPECTRUM_P5, 5) > 0
    sign_ok = (
        hecke_deviation_coeff() == -2
        and hecke_deviation_coeff_wrong_sign() == +2
        and fake_pos5
        and D_wrong_sign_expansion(SPECTRUM_P5, 5) != D5_spec
    )
    measured_sign = D5_delta < 0 and D7_delta < 0 and D11_delta > 0
    dual_nonneg = hecke_dual_nonneg_coeffs_p_ge_11()
    return {
        "proved": bool(
            ok_dim
            and ident5
            and ident7
            and tau_id
            and sign_ok
            and measured_sign
            and dual_nonneg is False
        ),
        "theorem": (
            "D = n(λ̄−6)² − 2∑_c ‖P_c(Φ−λ̄ I)‖_F² "
            "= 2∑_c dim(c)(τ − (λ_c−λ̄)²), τ=16(n+10)²/(n−6)³.  "
            "Hecke coefficients of the isotypic deviation norms are −2, "
            "so a dual SOS in End_G(Z) from traces alone cannot prove "
            "D≥0.  Fail: −2↦+2; 16↦8 in τ; claim D≥0 at p=5 or p=7."
        ),
        "D_p5": str(D5_delta),
        "D_p7": str(D7_delta),
        "D_p11": str(D11_delta),
        "D_negative_p5_p7": True,
        "D_positive_p11_census": True,
        "hecke_deviation_coeff": hecke_deviation_coeff(),
        "hecke_dual_nonneg_coeffs_p_ge_11": dual_nonneg,
        "identity_p5": ident5,
        "identity_p7": ident7,
        "tau_matches_R1": tau_id,
        "dim_Z_p5": dim_Z(5),
        "n_principal_p5": n_principal_constituents(5),
        "n5": n5,
        "n7": n7,
    }


def theorem_E_known_majorants_exceed_n_over_12(
    primes=(5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47),
) -> dict:
    """Shipped Max+-free δ² majorants sit above n/12 at every p≥5."""
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        n = n_of(p)
        n12 = Fraction(n, 12)
        r24 = room_orth_delta_sq(p)
        kh = kappa_hyp_delta_sq(p)
        r1 = r1_threshold(p)
        row_ok = (
            r24 > n12
            and kh > n12
            and r24 > r1
            and kh > r1
            and room_orth_delta_sq_wrong_drop24(p) != r24
        )
        rows[str(p)] = {
            "n_over_12": str(n12),
            "R1": str(r1),
            "room_over_24": str(r24),
            "kappa_hyp_delta": str(kh),
            "room_over_24_exceeds_n12": bool(r24 > n12),
            "ok": row_ok,
        }
        ok = ok and row_ok
    p11_fake = room_orth_delta_sq(11) <= Fraction(n_of(11), 12)
    return {
        "proved": bool(ok and not p11_fake),
        "proves_r1_for_p_ge_11": False,
        "theorem": (
            "15.100/15.102 room/24 = 4n(p²−9)/(3(p²−5)) and the κ_hyp "
            "form both exceed n/12 and R1(p) at every prime p≥5.  Together "
            "with the Θ(n³) 2-design+distance majorant (C), no shipped "
            "Max+-free estimate proves ‖δ‖²≤n/12 for p≥11.  Fail: claim "
            "room/24 ≤ n/12 at p=11; /24↦/12."
        ),
        "by_p": rows,
    }


def theorem_G_interpolating_majorant_algebra(
    primes=(5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47),
) -> dict:
    """B=κ_hyp_δ·4/(p−3)²: census ≥ measured at 5,7; ≤n/12 at p≥11 algebra.
    ‖δ‖²≤B remains OPEN."""
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p) or p <= 3:
            continue
        n = n_of(p)
        b = interpolating_majorant(p)
        n12 = Fraction(n, 12)
        h = n12_minus_B_sufficient_cubic(p)
        row = {
            "B": str(b),
            "n_over_12": str(n12),
            "B_per_n": str(b / n),
            "h": h,
            "B_le_n12": b <= n12,
        }
        if p in MEASURED_DELTA_SQ:
            row["ge_measured"] = b >= MEASURED_DELTA_SQ[p]
            row_ok = (
                interpolating_majorant_wrong_drop4(p) != b
                and interpolating_majorant_wrong_p_minus_1(p) != b
                and n12_minus_B_sufficient_cubic_wrong_drop69(p) != h
                and n12_minus_B_sufficient_cubic_wrong_drop3(p) != h
            )
            if p < 11:
                row_ok = row_ok and b > n12 and b >= MEASURED_DELTA_SQ[p]
            else:
                row_ok = (
                    row_ok
                    and b <= n12
                    and h > 0
                    and b >= MEASURED_DELTA_SQ[p]
                )
        else:
            row_ok = b <= n12 and h > 0
        row["ok"] = row_ok
        rows[str(p)] = row
        ok = ok and row_ok
    eq5 = interpolating_majorant(5) == MEASURED_DELTA_SQ[5]
    algebra = interpolating_majorant_le_n_over_12_for_p_ge_11()
    delta_le = delta_sq_le_interpolating_majorant()
    return {
        "proved": bool(ok and eq5 and algebra and delta_le is False),
        "proves_delta_le_B": False,
        "B_le_n12_algebra_p_ge_11": algebra,
        "theorem": (
            "B(p)=κ_hyp_δ(p)·4/(p−3)² equals measured ‖δ‖² at p=5 and "
            "dominates p=7,11 census.  B≤n/12 for p≥11 follows from "
            "h(p)=p³−3p²−69p−177 with h(11)=32 and h'>0 for p≥7 "
            "(sufficient comparison (p²−5)(p−3)≥64(p+3)).  Fail: 4↦2; "
            "(p−3)↦(p−1); drop 69 or 3 in h; claim B≤n/12 at p=5.  "
            "‖δ‖²≤B is not a bound (theorem H)."
        ),
        "by_p": rows,
        "equals_measured_p5": eq5,
        "delta_sq_le_B": delta_le,
        "h11": n12_minus_B_sufficient_cubic(11),
    }


def theorem_H_interpolant_killed_as_delta_bound() -> dict:
    """Kill 4/(p−3)² as a retained bound on ‖δ‖²."""
    eq_fails = interpolant_equality_law_fails_p7_p11()
    eq_p7_fake = interpolant_equality_law_wrong_claim_p7()
    nu_kill = nu_G_ratio_sq_killed_at_p7()
    nu7 = nu_G_ratio_sq_majorant(7)
    meas7 = MEASURED_DELTA_SQ[7]
    retained = interpolant_retained_as_delta_bound()
    delta_le = delta_sq_le_interpolating_majorant()
    r1 = r1_l2_bound_for_p_ge_11()
    fail_nu = nu_G_ratio_sq_majorant_wrong_nu5(7) != nu7
    p5_eq = MEASURED_DELTA_SQ[5] == kappa_hyp_delta_sq(5)
    return {
        "proved": bool(
            eq_fails
            and eq_p7_fake is False
            and nu_kill
            and meas7 > nu7
            and fail_nu
            and retained is False
            and delta_le is False
            and r1 is False
            and p5_eq
        ),
        "interpolant_retained_as_delta_bound": retained,
        "equality_law_fails_p7_p11": eq_fails,
        "nu_G_ratio_killed_p7": nu_kill,
        "nu_G_p7_majorant": str(nu7),
        "measured_p7": str(meas7),
        "theorem": (
            "Killed as a δ-bound.  (i) δ²=κ_hyp_δ·4/(p−3)² is false at "
            "p=7 and p=11.  Fail: claim equality at p=7.  (ii) the Aut-dim "
            "interpolant κ_hyp_δ·(2/ν_G)² is false at p=7 "
            f"({meas7} > {nu7}).  Fail: ν_G(5)↦1.  (iii) no operator "
            "identity yields 4/(p−3)².  Do not retain B as a bound on δ.  "
            "r1_l2 stays False.  p=5 measured=κ_hyp_δ still holds."
        ),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    A = theorem_A_threshold_formula()
    B = theorem_B_measured_vs_r1()
    C = theorem_C_moment_problem_too_weak()
    D = theorem_D_hecke_gap_expansion()
    E = theorem_E_known_majorants_exceed_n_over_12()
    G = theorem_G_interpolating_majorant_algebra()
    H = theorem_H_interpolant_killed_as_delta_bound()
    out = {
        "title": "R1 principal L2 room (not a close)",
        "numbered": False,
        "A_threshold": A,
        "B_measured": B,
        "C_moment_problem_too_weak": C,
        "D_hecke_gap": D,
        "E_known_majorants_too_weak": E,
        "G_interpolating_majorant_algebra_only": G,
        "H_interpolant_killed": H,
        "r1_l2_bound_for_p_ge_11": r1_l2_bound_for_p_ge_11(),
        "interpolant_retained_as_delta_bound": interpolant_retained_as_delta_bound(),
        "flags_not_flipped": [
            "type_I",
            "phi_F_ge_6",
            "residual_ii",
            "e1",
            "L",
            "qvar_k_ge_7",
        ],
        "L_status": "OPEN",
        "p13_orbits_not_a_close": True,
    }
    path = ROOT / "evidence" / "e1_gmin_r1_principal_pge11.json"
    write_json_atomic(path, out)
    print("R1 principal L2 room (unnumbered, not a close)", flush=True)
    print(f"  A threshold formula: {A['proved']}", flush=True)
    print(f"  B measured vs R1: {B['proved']}", flush=True)
    print(f"  C 2-design+distance too weak: {C['proved']}", flush=True)
    print(f"  D Hecke gap expansion: {D['proved']}", flush=True)
    print(f"  E known majorants exceed n/12: {E['proved']}", flush=True)
    print(f"  G interpolating B vs n/12 algebra: {G['proved']}", flush=True)
    print(f"  H interpolant killed as δ-bound: {H['proved']}", flush=True)
    print(f"  interpolant_retained: {interpolant_retained_as_delta_bound()}", flush=True)
    print(f"  r1_l2_bound_for_p_ge_11: {r1_l2_bound_for_p_ge_11()}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
