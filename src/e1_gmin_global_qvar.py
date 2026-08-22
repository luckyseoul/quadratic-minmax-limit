#!/usr/bin/env python3
"""
Global QVAR — leftover-1 quartic-variance conjunct, all k mixed.

Does **not** number leftover 1 (principal δ-room is a separate import).
Does **not** require QVAR on each k-stratum.  p=13 orbits / mesh k=6
are not a general close.  Aut-Schur / Gsum / pairing stay their own
units.  L follows the four-leftover AND.

CLAIM (Max+-free)
  Over the full Max+ ensemble (every activity k, λ=0 not split),
      E|Z_ψ|² ≥ 3q(q-1)/16
  equivalently λ_exc ≥ 6, equivalently S_□ ≥ 6q², equivalently the
  Gauss 4-distinct pairing of m₄ is ≥0, equivalently F̂(ψ)≥0
  (15.589 E, 15.279 L/M).  Fail: drop 16; replace 32 by 16;
  demand the bound on each k-stratum.

PER-STRATUM IS NOT GLOBAL (recorded, not p-laws)
  D. p=13 k=7 pointwise |Z_ψ|²=2548 < 10647/2.  Fail: claim 2548
     meets the floor.
  AF. p=41 k=7 nonempty Max+ (Cy=py, a_L=2p·30∈2pℤ) has ensemble
     E|Z_ψ|²=0 < 529515.  Fail: empty; fail: E≥floor; fail: use as
     a p-law.  Universal k≥7 QVAR is false; global mixing may still
     hold.

CLOSED INPUTS (not a global close)
  k=1..6 QVAR for every prime (15.589 FG, O, R, U).  For p≤11 one
  has m=(p+1)/2≤6, so those primes have no k≥7 stratum and global
  QVAR follows.  First open primes are p≥13.  Fail: claim k=1..6
  covers p=13 (m=7).  Do not split the λ=0 class: p=11 k=6 λ=0 has
  E B²=137/36<45/8 and the mixture still clears (15.589 K).

OPEN
  Gauss 4-dist pairing of m₄ / F̂(ψ)≥0 on the full Max+ fiber
  (15.279 L inequality).  Census exceed-floor at p=5,7 is not a
  p-law.  Wick Q(1)=8q² is global-only and does not sign F.
  Aut-inv master equation + |m₄|≤1 cannot prove the pairing
  (`e1_gmin_qvar_box_master`: LP min=−285/4 at p=5).
  For p≡3 the square-direction Gram is the Singer circulant and
  ψ is Nyquist (theorem H); the Nyquist eigenvalue is unnamed.
  Theorem I names the Ω-orbit mass S_j=p∑n_s²−k² from OA occupancy
  (Goryainov–Lin / Fourier of a line) and reduces p≡3 QVAR to the
  Nyquist of occupancy energy; the inequality is still OPEN.

Until the pairing is proved, global_qvar_proved_general stays False
by importing 15.279 L's inequality_proved / G / H (no handwritten
True/False).
"""
from __future__ import annotations

import inspect
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    EN_on_nonsquares,
    EN_on_squares,
    L_on_nonsquares,
    theorem_bool_coll_orthogonality,
    theorem_floor_iff_gauss_4dist,
)
from e1_gmin_m4_prop15589 import (  # noqa: E402
    lambda_exc_from_quartic_variance,
    n_of,
    q_of,
    quartic_variance_floor_threshold,
    spherical_QVAR_gap,
    spherical_quartic_variance,
    theorem_E_exceptional_quartic_variance,
    theorem_FG_profile_energy_and_low_strata,
    theorem_K_full_support_top_degree_mixing,
    theorem_O_k4_QVAR_all_primes,
    theorem_R_k5_QVAR_all_primes,
    theorem_U_k6_QVAR_all_primes,
)
from e1_gmin_qvar_k_ge_7 import (  # noqa: E402
    P13_K7_POINTWISE_ABS_ZPSI_SQ,
    T_of,
    pointwise_qvar_false_p13_k7,
    pointwise_qvar_wrong_claim_meets_floor,
    qvar_k_ge_7_proved_general,
    theorem_AF_p41_k7_stratum_fails_qvar,
    theorem_D_pointwise_orbit_counterexample,
)


def qvar_threshold_wrong_drop16(p: int) -> Fraction:
    q = q_of(p)
    return Fraction(3 * q * (q - 1), 8)


def top_activity(p: int) -> int:
    """m=(p+1)/2, the top profile stratum."""
    return (p + 1) // 2


def mean_Phi(p: int) -> Fraction:
    """Trace mean of Φ on Z: μ=8(n−2)/(n−6)."""
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def mean_Phi_wrong_drop_n6(p: int) -> Fraction:
    """Fail-eq: drop n−6."""
    n = n_of(p)
    return Fraction(8 * (n - 2), 1)


def squares_group_order(p: int) -> int:
    """|T|=(q−1)/2, T=F_q^×²."""
    return (q_of(p) - 1) // 2


def off_pm1_count(p: int) -> int:
    """|T \\ {±1}| = (q−5)/2."""
    return (q_of(p) - 5) // 2


def Q_pm1(p: int) -> int:
    """Boolean/Wick Q(±1)=8q² (15.279 D, 15.290)."""
    q = q_of(p)
    return 8 * q * q


def wick_off_value(p: int) -> int:
    """Wick Q(r)=4q² on squares r≠±1."""
    q = q_of(p)
    return 4 * q * q


def sum_T_Q(p: int) -> int:
    """∑_{r∈T} Q(r)=2q²(q−1).  Parseval: ∑_Ω u=q(q−1), E u=2q."""
    q = q_of(p)
    return 2 * q * q * (q - 1)


def total_wick_deficit(p: int) -> int:
    """A+B = ∑_{T\\{±1}} (4q²−Q(r)) = 8q².  Max+-free Parseval+Q(±1)."""
    return wick_off_value(p) * off_pm1_count(p) - (
        sum_T_Q(p) - 2 * Q_pm1(p)
    )


def total_wick_deficit_wrong_drop8(p: int) -> int:
    """Fail-eq: 8q² ↦ 4q²."""
    q = q_of(p)
    return 4 * q * q


def qvar_iff_B_ge_3q2(B: int, p: int) -> bool:
    """λ=2B/q² ≥ 6 ⇔ B≥3q²."""
    q = q_of(p)
    return B * 2 >= 6 * q * q


def qvar_iff_B_wrong_2q2(B: int, p: int) -> bool:
    """Fail-eq: B≥3q² ↦ B≥2q² (λ≥4, the 2-pairing Wick)."""
    q = q_of(p)
    return B * 2 >= 4 * q * q


def singer_theta_psi(p: int) -> int:
    """ψ(π^{p-1}) = i^{p-1} for a primitive quartic with ψ(π)=i.

    p≡3 (mod 4) ⇒ i^{p-1}=−1 (Nyquist on the square-direction cycle).
    p≡1 (mod 4) ⇒ i^{p-1}=+1 (the cycle is not alternating).
    """
    return -1 if p % 4 == 3 else 1


def line_intersection_var(p: int) -> Fraction:
    """Var(|D ∩ ℓ|) for an affine F_p-line ℓ in a square direction.

    Max+-free: theorem V gives E[N|χ=1]=p(p−1)/4, Aut 2-transitive
    on finite points, and every pair on a square line has χ=1.
    Fail-eq: hypergeometric (p−1)/4 of a random |D|-set.
    """
    return Fraction(p - 1, 2)


def line_intersection_var_wrong_hypergeometric(p: int) -> Fraction:
    """Fail-eq: unconstrained |D|-subset line variance."""
    return Fraction(p - 1, 4)


def profile_energy_mean(p: int) -> int:
    """E[a_L]=S/m=p(p−1)/2 on each square direction."""
    return p * (p - 1) // 2


def nyquist_eig_threshold(p: int) -> Fraction:
    """QVAR ⇔ m λ_η ≥ 3q(q−1)/16 ⇔ λ_η ≥ 3p²(p−1)/8 for p≡3."""
    q = q_of(p)
    m = top_activity(p)
    return Fraction(3 * q * (q - 1), 16 * m)


def lpp_pairing_threshold(p: int) -> Fraction:
    """For p≡3: E|Z_ψ|²=((q−1)/2)∑_{□}ψ L_{++}, so QVAR ⇔ ∑ψ L_{++}≥3q/8."""
    return Fraction(3 * q_of(p), 8)


def lpp_pairing_threshold_wrong_drop2(p: int) -> Fraction:
    """Fail-eq: 3q/8 ↦ 3q/4."""
    return Fraction(3 * q_of(p), 4)


def named_sum_N2_off0(p: int) -> Fraction:
    """15.317 B: E[∑_{δ≠0} N(δ)²]=k⁴/q+|Ω|q/2−k²."""
    q = q_of(p)
    k = p * (p - 1) // 2
    omega = (q - 1) // 2
    return Fraction(k**4, q) + Fraction(omega * q, 2) - k * k


def regular_set_size(p: int) -> int:
    return p * (p - 1) // 2


def square_N_sum_pointwise(p: int) -> int:
    """∑_{χ(δ)=1} N(δ) = [k(k−1)+k(p+1)/2]/2 on every Max+ (15.318 A)."""
    k = regular_set_size(p)
    return (k * (k - 1) + k * (p + 1) // 2) // 2


def chi_N_sum_pointwise(p: int) -> int:
    """∑_{δ≠0} χ(δ) N(δ)=k(p+1)/2, Max+-free."""
    k = regular_set_size(p)
    return k * (p + 1) // 2


def chi_N_sum_wrong_pk(p: int) -> int:
    """Fail-eq: ∑ χ N = p k."""
    return p * regular_set_size(p)


def plancherel_hatD_mass(p: int) -> int:
    """U_++U_-=k(q-k)=∑_{α≠0}|1̂_D(α)|², pointwise (15.474 B)."""
    q = q_of(p)
    k = regular_set_size(p)
    return k * (q - k)


def delta_imbalance_qvar_threshold(p: int) -> Fraction:
    """E Δ² ≥ 3q²(q-1)/16  ⇔  E|Z_ψ|² ≥ 3q(q-1)/16, via |Z|²=Δ²/q."""
    q = q_of(p)
    return Fraction(3 * q * q * (q - 1), 16)


def Uplus_Uminus_qvar_ceiling(p: int) -> Fraction:
    """QVAR ⇔ E[U_+ U_-] ≤ (K² − T)/4, K=k(q−k), T=3q²(q−1)/16.

    Pointwise Δ² = K² − 4 U_+ U_- (15.474 conservation).  Z=0 vectors
    saturate U_+=U_-=K/2 and are compatible with this *upper* bound
    on the product, unlike a pointwise SOS of Δ²−T.  Fail: drop 4.
    """
    K = plancherel_hatD_mass(p)
    T = delta_imbalance_qvar_threshold(p)
    return (K * K - T) / 4


def Uplus_Uminus_qvar_ceiling_wrong_drop4(p: int) -> Fraction:
    """Fail-eq: drop the 4."""
    K = plancherel_hatD_mass(p)
    T = delta_imbalance_qvar_threshold(p)
    return K * K - T


def delta_sq_from_wick_B(B, p: int) -> Fraction:
    """E Δ² = (q−1)B/16.  QVAR B≥3q² ⇔ E Δ² ≥ T.  Fail: drop 16."""
    return Fraction((q_of(p) - 1) * B, 16)


def delta_sq_from_wick_B_wrong_drop16(B, p: int) -> Fraction:
    return Fraction((q_of(p) - 1) * B, 1)


def sigma_hd(p: int) -> int:
    """Hasse–Davenport G(χ_{p²})=σ p, σ=−(−1)^{(p−1)/2} (15.474 A)."""
    return -((-1) ** ((p - 1) // 2))


def chi_on_trace_zero(p: int) -> int:
    """χ(α)=σ for every α∈ker Tr \\ {0}.

    α^p=−α ⇒ α^{p−1}=−1 ⇒ χ(α)=α^{(p−1)(p+1)/2}=(−1)^{(p+1)/2}=σ.
    Fail: χ=−σ (puts square duals off Ω).
    """
    return (-1) ** ((p + 1) // 2)


def chi_on_trace_zero_wrong_neg_sigma(p: int) -> int:
    """Fail-eq: χ|ker Tr = −σ."""
    return -sigma_hd(p)


def omega_size(p: int) -> int:
    return (q_of(p) - 1) // 2


def dual_orbit_size(p: int) -> int:
    """|{α≠0 : Tr(α v)=0}| = p−1 for any direction v."""
    return p - 1


def square_duals_fill_omega(p: int) -> bool:
    """m·(p−1)=|Ω|.  Square-direction duals partition Ω."""
    return top_activity(p) * dual_orbit_size(p) == omega_size(p)


def occupancy_energy_mean(p: int) -> int:
    """E[∑_s n_s²]=p(μ²+Var)=p(p²−1)/4 on a square class.

    μ=Var=(p−1)/2 (15.305 C mean + theorem H line variance).
    """
    return p * (p * p - 1) // 4


def occupancy_energy_sum_pointwise(p: int) -> int:
    """∑_{square dirs j} E_j = (k(q−k)+m k²)/p = p(p+1)(p²−1)/8.

    Pointwise from Plancherel (15.474 B) + S_j=p E_j−k².  Fail: m E[E_j]
    only in expectation.
    """
    return p * (p + 1) * (p * p - 1) // 8


def occupancy_energy_sum_wrong_drop_p(p: int) -> int:
    """Fail-eq: forget the /p in (k(q−k)+m k²)/p."""
    k = regular_set_size(p)
    return plancherel_hatD_mass(p) + top_activity(p) * k * k


def orbit_mass_from_energy(energy: int, p: int) -> int:
    """S_j=p E_j−k².  Fail: drop p."""
    k = regular_set_size(p)
    return p * energy - k * k


def orbit_mass_from_energy_wrong_drop_p(energy: int, p: int) -> int:
    k = regular_set_size(p)
    return energy - k * k


def nyquist_occupancy_energy_threshold(p: int) -> Fraction:
    """p≡3: Δ=p ∑ψ(L_j) E_j, QVAR ⇔ E[Nyquist E]² ≥ 3p²(p²−1)/16."""
    return Fraction(3 * p * p * (p * p - 1), 16)


def singer_cycle_even(p: int) -> bool:
    """m=(p+1)/2 is even ⇔ p≡3 (mod 4).  Odd cycle is not 2-colorable."""
    return (p + 1) // 2 % 2 == 0


def theorem_A_global_floor_iff(
    primes=(5, 7, 11, 13, 17, 19, 23, 41),
) -> dict:
    """λ_exc≥6 ⇔ E|Z_ψ|²≥3q(q-1)/16 ⇔ S_□≥6q² ⇔ Gauss 4-dist ≥0."""
    E = theorem_E_exceptional_quartic_variance()
    L = theorem_floor_iff_gauss_4dist(list(primes))
    M = theorem_bool_coll_orthogonality(list(primes))
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        thr = quartic_variance_floor_threshold(p)
        lam6 = lambda_exc_from_quartic_variance(p, thr)
        row_ok = lam6 == 6 and qvar_threshold_wrong_drop16(p) != thr
        if p in (5, 7):
            rec = E["by_p"][str(p)]
            var = Fraction(rec["E_abs_Zpsi_sq"])
            row_ok = (
                row_ok
                and lambda_exc_from_quartic_variance(p, var)
                == Fraction(rec["lambda_exc"])
                and var >= thr
            )
        rows[str(p)] = {
            "thr": str(thr),
            "lambda_exc_at_floor": str(lam6),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(
            ok
            and E["proved_reduction"]
            and L["proved"]
            and M["proved"]
        ),
        "inequality_proved": bool(L["inequality_proved"]),
        "per_stratum_equivalent": False,
        "theorem": (
            "On the full Max+ ensemble, λ_exc=32 E|Z_ψ|²/[q(q-1)] "
            "⇔ E|Z_ψ|²≥3q(q-1)/16 ⇔ S_□≥6q² ⇔ the Gauss 4-distinct "
            "pairing of m₄ is ≥0 ⇔ F̂(ψ)≥0 (15.589 E, 15.279 L/M).  "
            "Fail: drop 16.  The k-stratum restriction is not a "
            "2-design and is not equivalent to this identity."
        ),
        "by_p": rows,
        "census_formula": E["proved_census"],
    }


def theorem_B_per_stratum_is_not_global() -> dict:
    """(13,7) pointwise and (41,7) ensemble misses are not p-laws."""
    D = theorem_D_pointwise_orbit_counterexample()
    AF = theorem_AF_p41_k7_stratum_fails_qvar()
    q41 = 41 * 41
    thr41 = 3 * q41 * (q41 - 1) // 16
    ok = (
        D["proved"]
        and pointwise_qvar_false_p13_k7()
        and not pointwise_qvar_wrong_claim_meets_floor()
        and D["pointwise_abs_Zpsi_sq"] == P13_K7_POINTWISE_ABS_ZPSI_SQ
        and AF["proved"]
        and AF["stratum_qvar"] is False
        and AF["E_abs_Zpsi_sq"] == 0
        and AF["maxplus_Cy_eq_py"] is True
        and AF["a_L_in_2pZ"] is True
        and AF["E_abs_Zpsi_sq"] < AF["QVAR_threshold"]
        and AF["QVAR_threshold"] == thr41
        and qvar_k_ge_7_proved_general() is False
        and top_activity(13) == 7
        and top_activity(41) == 21
    )
    return {
        "proved": bool(ok),
        "p13_k7_pointwise_abs_Zpsi_sq": P13_K7_POINTWISE_ABS_ZPSI_SQ,
        "p13_k7_pointwise_below_floor": True,
        "p41_k7_E_abs_Zpsi_sq": 0,
        "p41_k7_Cy_eq_py": True,
        "p41_k7_a_L_in_2pZ": True,
        "p41_k7_stratum_qvar": False,
        "imported_as_p_law": False,
        "qvar_k_ge_7_proved_general": False,
        "theorem": (
            "p=13 k=7 has a Max+ witness |Z_ψ|²=2548<10647/2 (fail: "
            "meets floor).  p=41 k=7 is a nonempty Max+ stratum "
            "(Cy=py, a_L=2p·30) with ensemble E|Z_ψ|²=0<529515 "
            "(fail: empty; fail: E≥floor).  Neither is a p-law and "
            "neither is required of the global mixture.  Per-stratum "
            "k≥7 QVAR is false, so leftover-1 must not import "
            "qvar_k_ge_7_proved_general."
        ),
    }


def theorem_C_k1_through_k6_closed_not_a_global_cover() -> dict:
    """k=1..6 QVAR all primes; covers p≤11 only, not p≥13."""
    FG = theorem_FG_profile_energy_and_low_strata()
    O = theorem_O_k4_QVAR_all_primes()
    R = theorem_R_k5_QVAR_all_primes()
    U = theorem_U_k6_QVAR_all_primes()
    small = all(top_activity(p) <= 6 for p in (5, 7, 11))
    p13_has_k7 = top_activity(13) >= 7
    ok = (
        FG["proved_k1_k3_QVAR_all_primes"]
        and O["proved"]
        and R["proved"]
        and U["proved"]
        and small
        and p13_has_k7
        and top_activity(5) == 3
        and top_activity(7) == 4
        and top_activity(11) == 6
        and top_activity(13) == 7
    )
    return {
        "proved": bool(ok),
        "covers_p_le_11": True,
        "covers_p_ge_13": False,
        "k1_k3": FG["proved_k1_k3_QVAR_all_primes"],
        "k4": O["proved"],
        "k5": R["proved"],
        "k6": U["proved"],
        "theorem": (
            "QVAR holds on every k=1..6 stratum for every prime "
            "(15.589 FG/O/R/U).  For p≤11 one has m≤6, so the global "
            "ensemble is a mixture of closed strata.  Fail: claim this "
            "covers p=13 (m=7) or every p≥13."
        ),
    }


def theorem_D_do_not_split_lambda_zero() -> dict:
    """p=11 k=6 λ=0 misses QVAR; the unsplit mixture clears."""
    K = theorem_K_full_support_top_degree_mixing()
    return {
        "proved": bool(K["proved_counterexample"]),
        "lambda_zero_alone_can_miss": True,
        "mixture_clears_p11": True,
        "theorem": (
            "The λ=0 class of a single stratum can miss QVAR "
            "(p=11 k=6: E B²=137/36<45/8) while the unsplit mixture "
            "clears (15.589 K).  Global QVAR does not split λ=0."
        ),
        "p11_k6_lambda0_E_B2": "137/36",
        "p11_k6_mixture_E_B2": "114771/14903",
    }


def theorem_E_census_exceeds_floor_not_a_close() -> dict:
    """p=5,7 exact ensembles exceed the floor; not a p≥13 close."""
    E = theorem_E_exceptional_quartic_variance()
    ok = bool(E["proved_census"] and not E["proved_general_inequality"])
    rows = {}
    for p in (5, 7):
        rec = E["by_p"][str(p)]
        var = Fraction(rec["E_abs_Zpsi_sq"])
        thr = quartic_variance_floor_threshold(p)
        vs = spherical_quartic_variance(p)
        gap = spherical_QVAR_gap(p)
        rows[str(p)] = {
            "E_abs_Zpsi_sq": str(var),
            "thr": str(thr),
            "V_sph": str(vs),
            "clears_floor": var >= thr,
            "clears_V_sph": var >= vs,
            "V_sph_gap": str(gap),
        }
        ok = ok and var >= thr and vs - thr == gap and gap > 0
    return {
        "proved": bool(ok),
        "covers_general": False,
        "by_p": rows,
        "theorem": (
            "Exact Max+ censuses at p=5,7 have E|Z_ψ|² above both the "
            "QVAR floor and V_sph.  Fail: interpolate those two values "
            "as a p-law.  V_sph>floor is a 2-design spherical average, "
            "not a 4-design, so it does not prove the floor at p≥13."
        ),
    }


def theorem_F_mean_above_floor_ordering_open(
    primes=(5, 7, 11, 13, 17, 19, 23, 41),
) -> dict:
    """μ>6 always.  λ_exc=λ_max would close QVAR; ordering is not a p-law."""
    E = theorem_E_exceptional_quartic_variance()
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        mu = mean_Phi(p)
        row_ok = mu > 6 and mean_Phi_wrong_drop_n6(p) != mu
        if p in (5, 7):
            lam = Fraction(E["by_p"][str(p)]["lambda_exc"])
            row_ok = row_ok and lam > mu
        rows[str(p)] = {
            "mu": str(mu),
            "mu_gt_6": mu > 6,
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "ordering_proved_general": False,
        "suffices_if_exc_is_max": True,
        "theorem": (
            "tr(Φ)/dim Z = μ=8(n−2)/(n−6)>6 for every prime p≥5 "
            "(fail: drop n−6).  If the exceptional scalar were the "
            "maximal Φ-eigenvalue then λ_exc≥μ>6 and global QVAR "
            "would follow.  That ordering holds at p=5,7 (census "
            "λ_exc>μ) and is reported at p=11 (15.589 C) but is not "
            "a p-law: the exceptional block may drop below the "
            "principal series.  Fail: infer λ_exc>6 from μ>6 alone."
        ),
        "by_p": rows,
    }


def theorem_G_nyquist_deficit_split(
    primes=(5, 7, 11, 13, 17, 19, 23, 41),
) -> dict:
    """Exceptional ψ is the Nyquist character of T.  A+B=8q²; QVAR iff B≥3q².

    T=F_q^×² cyclic of even order (q−1)/2.  The order-4 character with
    ψ²=χ restricts to the unique order-2 character η of T (constant
    phase-step π in dlog).  Q(r)=E[|ẑ(ξ)|² |ẑ(rξ)|²] is the lag-r
    autocorrelation on the T-torsor Ω, hence positive definite on T
    (PSL square-multiplication stationarity).  Bochner gives only
    S_□=∑_T η Q ≥ 0, i.e. λ≥0, too weak.

    Parseval: ẑ(0)=p pointwise, ∑_Ω u = q(q−1), E u=2q, so
    ∑_T Q = 2q²(q−1).  Q(±1)=8q².  The Wick deficit
        A = ∑_{H\\{±1}} (4q²−Q),   B = ∑_K (4q²−Q)
    on fourth-powers-off-pm1 vs the complementary squares therefore
    satisfies A+B=8q² (fail: 8↦4).  Then S_□=2B and
        λ_exc = 2B/q²  ⇔  QVAR  B≥3q²
    (fail: B≥2q², which is λ≥4).  Equal density on |K|=(q−1)/4 vs
    |H\\{±1}|=(q−9)/4 would give B=4q²(q−1)/(q−5)>3q² at every
    p≥5, but that comparison is not proved.  Sign of A−B OPEN.
    """
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        q = q_of(p)
        ab = total_wick_deficit(p)
        off = off_pm1_count(p)
        nT = squares_group_order(p)
        nH_off = (q - 9) // 4
        nK = (q - 1) // 4
        # equal-density B would suffice
        eqB = Fraction(4 * q * q * (q - 1), q - 5)
        row_ok = (
            nT == nH_off + nK + 2
            and off == nH_off + nK
            and ab == 8 * q * q
            and ab == total_wick_deficit(p)
            and ab != total_wick_deficit_wrong_drop8(p)
            and Q_pm1(p) == 8 * q * q
            and wick_off_value(p) == 4 * q * q
            and sum_T_Q(p) == 2 * q * q * (q - 1)
            and qvar_iff_B_ge_3q2(3 * q * q, p)
            and not qvar_iff_B_ge_3q2(3 * q * q - 1, p)
            and qvar_iff_B_wrong_2q2(2 * q * q, p)
            and eqB > 3 * q * q
        )
        rows[str(p)] = {
            "A_plus_B": ab,
            "eight_q2": 8 * q * q,
            "off": off,
            "n_H_off": nH_off,
            "n_K": nK,
            "equal_density_B": str(eqB),
            "equal_density_suffices": eqB > 3 * q * q,
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "qvar_iff_B_ge_3q2": True,
        "bochner_only_lambda_ge_0": True,
        "equal_density_would_suffice": True,
        "theorem": (
            "On T, exceptional ψ is Nyquist (η).  ∑_T Q=2q²(q−1) and "
            "Q(±1)=8q² force the Wick deficit A+B=8q² (fail: 8↦4).  "
            "λ_exc=2B/q², so QVAR ⇔ B≥3q² (fail: B≥2q²).  PD/Bochner "
            "only gives B≥0.  Equal density B=4q²(q−1)/(q−5)>3q² at "
            "every p≥5, unproved.  A−B sign OPEN."
        ),
        "by_p": rows,
    }


def theorem_H_singer_circulant_p_eq_3_mod_4(
    primes=(7, 11, 19, 23, 31, 43),
) -> dict:
    """p≡3: square-direction Gram is the Singer circulant; ψ is Nyquist.

    Max+-free structure, inequality OPEN.

    Multiplication by a primitive (p+1)-st root θ=π^{p-1} cycles the
    p+1 F_p-lines in two orbits of size m=(p+1)/2 (orbit-stabilizer:
    |⟨θ⟩|=p+1, stabilizer {±1} on lines).  χ_q(θ)=1, so the square
    directions are one orbit.  A primitive quartic has ψ(π)=i, hence
    ψ(θ)=i^{p-1}=−1 precisely when p≡3 (mod 4) (fail: p≡1 gives +1).
    Profile energies (a_L) are therefore a palindromic circulant of
    m/2=(p+1)/4 pair-classes, and the quartic sign vector is the
    Nyquist mode of that cycle.  Paley Aut (F_q^×-multiplication)
    forces the Gram of the Max+ ensemble to lie in this Bose–Mesner
    algebra.  QVAR becomes
        m λ_η ≥ 3q(q−1)/16
    with λ_η the Nyquist eigenvalue of E[a a^T].

    Ridge occupancy: for a square line pencil, every pair has χ=1,
    so theorem V plus 2-transitivity give Var(|D∩ℓ|)=(p−1)/2
    (fail: hypergeometric (p−1)/4).  Expanding a_L = ∑_{δ∈L^×} N(δ)
    + p(p−1)(3−p)/4, the lags are F_p^×-periodizations of L(r)=
    E[N(1)N(r)].  Theorem W names L on nonsquares; L on squares
    (hence the lags) stays OPEN.  Coarse {t,t+1} cyclic b-orbits
    are not this Aut-circulant and miss the floor (NOTE 2026-08-20
    §5).  Census lags at p=7,11 clear the floor and are not p-laws.
    """
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p) or p % 4 != 3:
            continue
        q = q_of(p)
        m = top_activity(p)
        n_classes = m // 2
        mu = EN_on_squares(p)
        nD = p * (p - 1) // 2
        # Var(|D∩ℓ|) from 2-point of D on a square line
        mu_pt = Fraction(p - 1, 2 * p)
        var1 = mu_pt * (1 - mu_pt)
        cov_sq = Fraction(p - 1, 4 * p * p)
        var_k = p * var1 + p * (p - 1) * cov_sq
        row_ok = (
            singer_theta_psi(p) == -1
            and singer_theta_psi(5) == 1
            and m == (p + 1) // 2
            and n_classes == (p + 1) // 4
            and var_k == line_intersection_var(p)
            and var_k != line_intersection_var_wrong_hypergeometric(p)
            and mu == Fraction(nD, 2)
            and profile_energy_mean(p) == p * (p - 1) // 2
            and nyquist_eig_threshold(p)
            == Fraction(3 * q * (q - 1), 16 * m)
            and L_on_nonsquares(p) == EN_on_squares(p) * EN_on_nonsquares(p)
            and lpp_pairing_threshold(p) * (q - 1) / 2
            == Fraction(3 * q * (q - 1), 16)
            and lpp_pairing_threshold(p) != lpp_pairing_threshold_wrong_drop2(p)
            and named_sum_N2_off0(p) / (q - 1) * (q - 1)
            == named_sum_N2_off0(p)
            and square_N_sum_pointwise(p)
            == (
                regular_set_size(p) * (regular_set_size(p) - 1)
                + chi_N_sum_pointwise(p)
            )
            // 2
            and chi_N_sum_pointwise(p) != chi_N_sum_wrong_pk(p)
            and square_N_sum_pointwise(p) == (q - 1) * p * (p - 1) // 8
            and plancherel_hatD_mass(p)
            == regular_set_size(p) * (q - regular_set_size(p))
            and delta_imbalance_qvar_threshold(p)
            == lpp_pairing_threshold(p) * (q - 1) / 2 * q
            and Uplus_Uminus_qvar_ceiling(p) * 4
            == plancherel_hatD_mass(p) ** 2
            - delta_imbalance_qvar_threshold(p)
            and Uplus_Uminus_qvar_ceiling(p)
            != Uplus_Uminus_qvar_ceiling_wrong_drop4(p)
            and delta_sq_from_wick_B(3 * q * q, p)
            == delta_imbalance_qvar_threshold(p)
            and delta_sq_from_wick_B(3 * q * q, p)
            != delta_sq_from_wick_B_wrong_drop16(3 * q * q, p)
        )
        rows[str(p)] = {
            "m": m,
            "n_pair_classes": n_classes,
            "psi_theta": singer_theta_psi(p),
            "Var_k": str(var_k),
            "mean_a": profile_energy_mean(p),
            "lambda_eta_threshold": str(nyquist_eig_threshold(p)),
            "ok": row_ok,
        }
        ok = ok and row_ok
    # recorded census, not a p-law
    census = {
        "7": {
            "lags_E_aa": {
                "0": "233730/409",
                "1": "160524/409",
                "2": "166698/409",
            },
            "E_abs_Zpsi_sq": "317520/409",
            "lambda_exc": "4320/409",
            "circulant": True,
            "weights_nyquist_alt": True,
            "nunique_Lpp_on_squares": 6,
            "Lpp_two_level_on_squares": False,
            "N_S_even_dft_ratio_to_nyquist": {
                "2": "1/9",
                "4": "14/15",
                "6": "38/45",
                "8": "32/105",
                "12": "1",
            },
        },
        "11": {
            "lags_E_aa": {
                "0": "1517186330/425649",
                "1": "413328740/141883",
                "2": "177502160/60807",
                "3": "414446780/141883",
            },
            "E_abs_Zpsi_sq": "557807580/141883",
            "lambda_exc": "1229328/141883",
            "circulant": True,
            "weights_nyquist_alt": True,
        },
    }
    # p=7,11 census clears; not imported as a law
    c7 = Fraction(census["7"]["E_abs_Zpsi_sq"])
    c11 = Fraction(census["11"]["E_abs_Zpsi_sq"])
    ok = (
        ok
        and c7 >= Fraction(3 * 49 * 48, 16)
        and c11 >= Fraction(3 * 121 * 120, 16)
        and singer_theta_psi(47) == -1
        and singer_theta_psi(29) == 1
    )
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "covers_p_eq_1_mod_4": False,
        "imported_as_p_law": False,
        "qvar_iff_nyquist_eig": True,
        "qvar_iff_Lpp_pairing": True,
        "qvar_iff_delta_imbalance": True,
        "qvar_iff_Uplus_Uminus_ceiling": True,
        "eta_block_1dim_RR_tautological": True,
        "r_star_r_constructed": False,
        "EN2_two_level_by_chi": True,
        "Lpp_two_level_on_squares": False,
        "nonsquare_L_named": True,
        "census_p7_p11": census,
        "theorem": (
            "For p≡3 (mod 4) the square F_p-directions form one "
            "Singer cycle of length m=(p+1)/2 under θ=π^{p-1}, and "
            "ψ(θ)=i^{p-1}=−1 (fail: p≡1).  The Max+ Gram of a_L is "
            "the palindromic circulant of this cycle, so QVAR is the "
            "Nyquist eigenvalue λ_η ≥ 3p²(p−1)/8.  Line occupancy "
            "has Var(|D∩ℓ|)=(p−1)/2 (fail: hypergeometric (p−1)/4).  "
            "Z_ψ=∑_{□}ψ(δ)N(δ), so E|Z|²=((q−1)/2)∑_{□}ψ L_{++} and "
            "QVAR ⇔ ∑ψ L_{++}≥3q/8 (fail: 3q/4).  Aut F_q^{×2} is "
            "transitive on squares, so E[N(δ)²] is 2-level by χ "
            "(values unnamed; 15.317 B names only the mixed sum).  "
            "∑_{□}N is k-pointwise constant (15.318 A; fail: ∑χN=pk).  "
            "15.474: U_++U_-=k(q−k) pointwise, Δ=U_+−U_-, "
            "Δ²=K²−4U_+U_- so QVAR ⇔ E[U_+U_-]≤(K²−T)/4 "
            "(fail: drop 4; Z=0 saturates the product and is allowed).  "
            "L_{++} is not 2-level on squares (p=7: 6 values).  "
            "Lags are F_p-periodizations of L_{++}; theorem W names L "
            "on nonsquares, L_{++} on squares OPEN.  p=7,11 census "
            "lags are not a p-law.  Coarse {t,t+1} cyclic b-orbits miss "
            "the floor and are not this Aut-circulant."
        ),
        "by_p": rows,
    }


def theorem_I_oa_occupancy_orbit_mass(
    primes=(5, 7, 11, 13, 17, 19, 23, 31, 43),
) -> dict:
    """Ω-orbit mass S_j is p times square-direction occupancy energy.

    Max+-free structure (Goryainov–Lin OA + Fourier of a line).
    Inequality OPEN.  Not a p-law from p=5,7 energy lags.

    Paley P(p²) is the block graph of the orthogonal array
    OA(m,p), m=(p+1)/2 quadratic slopes (Goryainov–Lin 2021,
    arXiv:2104.08839; Godsil–Royle).  Canonical cliques are the
    affine F_p-lines in square directions.  Balanced clique
    indicators span the (p−1)/2-eigenspace, which is exactly the
    Ω-space of 15.474 A (hat 1_D supported on {χ=σ}).

    Dual of direction v: {α≠0: Tr(α v)=0} is an F_p^×-orbit of
    size p−1.  χ|F_p^×=1 so χ is constant on the orbit.  The
    subfield direction (v=1, a square) has dual ker Tr \\ {0}, and
    χ=σ there (fail: χ=−σ).  Counting m(p−1)=|Ω| puts every
    square-direction dual in Ω and partitions Ω.

    e_α is constant on each parallel line in direction v iff
    Tr(α v)=0, so hat 1_D(α)=∑_s n_s e_α(b_s).  Parseval on the
    F_p-quotient: ∑_{Tr(α v)=0} |hat|² = p ∑_s n_s².  Subtract
    α=0 (value k²):
        S_j := ∑_{a∈F_p^×} |hat 1_D(a ξ_j)|² = p E_j − k²
    (fail: drop p).  ∑_j S_j=k(q−k) pointwise, so
        ∑_j E_j = p(p+1)(p²−1)/8
    pointwise (fail: only in expectation).  Nonsquare occupancy
    is already locked at (p−1)/2 (15.305 C).  Square classes have
    E[E_j]=p(p²−1)/4 from Var=(p−1)/2.

    p≡3: m even, ψ|F_p^×=1, ψ(θ)=−1, so ψ is Nyquist-constant on
    each dual orbit and Δ=∑_j ψ(L_j) S_j = p ∑_j ψ(L_j) E_j.
    QVAR ⇔ E[Nyquist E]² ≥ 3p²(p²−1)/16.  p≡1: m odd (Singer not
    2-colorable); F_p has an order-4 character and ψ splits each
    dual orbit, so this pairing does not cover p≡1.
    """
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        q = q_of(p)
        k = regular_set_size(p)
        m = top_activity(p)
        mu = Fraction(p - 1, 2)
        mean_E = occupancy_energy_mean(p)
        sum_E = occupancy_energy_sum_pointwise(p)
        row_ok = (
            sigma_hd(p) == -((-1) ** ((p - 1) // 2))
            and chi_on_trace_zero(p) == sigma_hd(p)
            and chi_on_trace_zero(p) != chi_on_trace_zero_wrong_neg_sigma(p)
            and square_duals_fill_omega(p)
            and m * (p - 1) == (q - 1) // 2
            and mean_E == p * (p * p - 1) // 4
            and mean_E == p * (mu * mu + Fraction(p - 1, 2))
            and sum_E == p * (p + 1) * (p * p - 1) // 8
            and sum_E * p
            == plancherel_hatD_mass(p) + m * k * k
            and sum_E != occupancy_energy_sum_wrong_drop_p(p)
            and orbit_mass_from_energy(mean_E, p)
            != orbit_mass_from_energy_wrong_drop_p(mean_E, p)
            and singer_cycle_even(p) == (p % 4 == 3)
            and singer_cycle_even(p) != (p % 4 == 1)
            and nyquist_occupancy_energy_threshold(p)
            == Fraction(3 * p * p * (p * p - 1), 16)
            and (
                delta_imbalance_qvar_threshold(p)
                == nyquist_occupancy_energy_threshold(p) * p * p
            )
        )
        # Parseval check: m * named mean S = k(q-k)
        mean_S = orbit_mass_from_energy(mean_E, p)
        row_ok = row_ok and m * mean_S == plancherel_hatD_mass(p)
        # Cauchy min E_j = k²/p, S_j ≥ 0
        cauchy = k * k // p
        row_ok = (
            row_ok
            and orbit_mass_from_energy(cauchy, p) == 0
            and mean_E > cauchy
        )
        rows[str(p)] = {
            "m": m,
            "singer_even": singer_cycle_even(p),
            "mean_E": mean_E,
            "sum_E": sum_E,
            "cauchy_min_E": cauchy,
            "chi_kerTr": chi_on_trace_zero(p),
            "sigma": sigma_hd(p),
            "ok": row_ok,
        }
        ok = ok and row_ok
    # recorded census, not a p-law
    census = {
        "5": {
            "E_energy_lags": {"0": "12300/13", "1": "11400/13", "2": "11400/13"},
            "unique_E": [20, 30, 50],
            "cauchy_min_attained": True,
            "psi_constant_on_Fp_orbit": False,
            "m_odd": True,
        },
        "7": {
            "E_energy_lags": {
                "0": "2939265/409",
                "1": "2866059/409",
                "2": "2872233/409",
                "3": "2866059/409",
            },
            "unique_E": [63, 77, 91, 105, 147],
            "cauchy_min_attained": True,
            "psi_constant_on_Fp_orbit": True,
            "m_even": True,
        },
    }
    c5_0 = Fraction(census["5"]["E_energy_lags"]["0"])
    c7_0 = Fraction(census["7"]["E_energy_lags"]["0"])
    ok = (
        ok
        and c5_0 == Fraction(12300, 13)
        and c7_0 == Fraction(2939265, 409)
        and occupancy_energy_mean(5) == 30
        and occupancy_energy_mean(7) == 84
        and occupancy_energy_sum_pointwise(5) == 90
        and occupancy_energy_sum_pointwise(7) == 336
        and not singer_cycle_even(5)
        and singer_cycle_even(7)
    )
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "covers_p_eq_1_mod_4": False,
        "imported_as_p_law": False,
        "qvar_iff_nyquist_occupancy_energy": True,
        "S_j_eq_pE_minus_k2": True,
        "square_duals_partition_Omega": True,
        "psi_constant_on_Fp_orbit_p_eq_3": True,
        "eta_block_1dim_RR_tautological": False,
        "goryainov_lin_oa_basis": True,
        "census_p5_p7_energy_lags": census,
        "theorem": (
            "Paley P(p²) is the block graph of OA(m,p) with "
            "m=(p+1)/2 quadratic slopes (Goryainov–Lin; fail: "
            "nonsquare slopes).  Square-direction duals partition "
            "Ω: χ|ker Tr=σ (fail: −σ) and m(p−1)=|Ω|.  Fourier of "
            "a line gives S_j=p∑_s n_{j,s}²−k² (fail: drop p).  "
            "∑_j E_j=p(p+1)(p²−1)/8 pointwise (fail: expectation "
            "only).  Nonsquare n=(p−1)/2 locked (15.305 C).  For "
            "p≡3, ψ is Nyquist-constant on each dual orbit so "
            "Δ=p∑ψ(L_j) E_j and QVAR ⇔ E[Nyquist E]²≥3p²(p²−1)/16.  "
            "p≡1: m odd, ψ splits F_p-orbits, this pairing does not "
            "apply.  Occupancy-energy lags at p=5,7 are not a p-law.  "
            "Cauchy min E_j=k²/p (S_j=0) is attained on some Max+ "
            "rows and is compatible with Z=0; the ensemble Nyquist "
            "is OPEN."
        ),
        "by_p": rows,
    }


def theorem_P_pairing_positivity() -> dict:
    """F̂(ψ)≥0 / Gauss 4-dist pairing ≥0.  Identity proved; sign OPEN."""
    L = theorem_floor_iff_gauss_4dist()
    M = theorem_bool_coll_orthogonality()
    return {
        "proved": bool(L["proved"] and M["proved"]),
        "inequality_proved": bool(
            L["inequality_proved"] and M["inequality_proved"]
        ),
        "theorem": (
            "S_□−6q² equals the Gauss 4-distinct pairing of m₄ and "
            "equals F̂(ψ) for F=Q−Coll (15.279 L/M).  Wick already "
            "gives S_□,Wick=8q².  Positivity of the Boolean remainder "
            "is OPEN: no SoS for F, Paley+ω is not a scheme at p≥11, "
            "CS on ⟨ρ,K⟩ overshoots, 1D lifts do not exhaust Max+ and "
            "some Aut-orbits have S_□=0."
        ),
    }


def global_qvar_proved_general() -> bool:
    """True only by importing a pairing / Nyquist-deficit inequality.

    No handwritten True or False.  False while 15.279 L/M leave
    F̂(ψ)≥0 open, while B≥3q² is unproved, and while the p≡3
    Singer-circulant Nyquist eigenvalue is unnamed.
    """
    A = theorem_A_global_floor_iff()
    P = theorem_P_pairing_positivity()
    G = theorem_G_nyquist_deficit_split()
    H = theorem_H_singer_circulant_p_eq_3_mod_4()
    I = theorem_I_oa_occupancy_orbit_mass()
    return bool(
        A["proved"]
        and (
            P["inequality_proved"]
            or G["inequality_proved"]
            or H["inequality_proved"]
            or I["inequality_proved"]
        )
    )


def live_L_status() -> str:
    """L follows the four-leftover AND.  Not baked OPEN."""
    from e1_main_chain_status import four_e1_units_closed

    return "CLOSED" if four_e1_units_closed().get("closed") else "OPEN"


def main() -> dict:
    from io_atomic import write_json_atomic

    A = theorem_A_global_floor_iff()
    B = theorem_B_per_stratum_is_not_global()
    C = theorem_C_k1_through_k6_closed_not_a_global_cover()
    D = theorem_D_do_not_split_lambda_zero()
    E = theorem_E_census_exceeds_floor_not_a_close()
    F = theorem_F_mean_above_floor_ordering_open()
    Gdef = theorem_G_nyquist_deficit_split()
    Hcirc = theorem_H_singer_circulant_p_eq_3_mod_4()
    Iocc = theorem_I_oa_occupancy_orbit_mass()
    P = theorem_P_pairing_positivity()
    src = inspect.getsource(global_qvar_proved_general)
    out = {
        "title": "Global QVAR (not a close)",
        "numbered": False,
        "A_global_floor_iff": A,
        "B_per_stratum_is_not_global": B,
        "C_k1_through_k6_not_a_cover": C,
        "D_do_not_split_lambda_zero": D,
        "E_census_not_a_close": E,
        "F_mean_above_floor_ordering_open": F,
        "G_nyquist_deficit_split": Gdef,
        "H_singer_circulant_p_eq_3_mod_4": Hcirc,
        "I_oa_occupancy_orbit_mass": Iocc,
        "P_pairing_positivity": P,
        "global_qvar_proved_general": global_qvar_proved_general(),
        "qvar_k_ge_7_proved_general": qvar_k_ge_7_proved_general(),
        "p13_orbits_not_a_close": True,
        "L_status": live_L_status(),
        "no_handwritten_true": "return True" not in src,
    }
    path = ROOT / "evidence" / "e1_gmin_global_qvar.json"
    write_json_atomic(path, out)
    print("Global QVAR (unnumbered, not a close)", flush=True)
    print(f"  A floor iff: {A['proved']} inequality={A['inequality_proved']}", flush=True)
    print(f"  B per-stratum ≠ global: {B['proved']}", flush=True)
    print(f"  C k=1..6 not a p≥13 cover: {C['proved']}", flush=True)
    print(f"  D do not split λ=0: {D['proved']}", flush=True)
    print(f"  E census not a close: {E['proved']}", flush=True)
    print(
        f"  F μ>6 ordering open: {F['proved']} "
        f"ordering={F['ordering_proved_general']}",
        flush=True,
    )
    print(
        f"  G Nyquist deficit A+B=8q2: {Gdef['proved']} "
        f"inequality={Gdef['inequality_proved']}",
        flush=True,
    )
    print(
        f"  H Singer circulant p≡3: {Hcirc['proved']} "
        f"inequality={Hcirc['inequality_proved']}",
        flush=True,
    )
    print(
        f"  I OA occupancy orbit mass: {Iocc['proved']} "
        f"inequality={Iocc['inequality_proved']}",
        flush=True,
    )
    print(
        f"  P pairing positivity: {P['proved']} "
        f"inequality={P['inequality_proved']}",
        flush=True,
    )
    print(f"  global_qvar_proved_general: {global_qvar_proved_general()}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
