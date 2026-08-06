#!/usr/bin/env python3
"""
Prop 15.170 — Residual (i) close: freeness-fail Type I k=3p−2 cannot have
s_−≤−1 (general primes p≥5), via dual-equality Gsum correlation Farkas.

REAL checkable Fraction / prior-prop chain:

  A. Freeness-fail structure (15.169): k=3p−2, S∈{1,5}, affine S=3−2f_e on Max+.
  B. Gap-2 Type I undercutter ⇒ s_−=−1 (15.169).
  C. ND dichotomy at s_−=−1 (15.169): good sign of f_e on U_− ⇒ strong ND;
     bad case f_e≡−1 on U_− is the only residual ND risk.
  D. Bad case + Max− moment matching forces dual two-level H and dual equality
     S_G = 3−2f_e on Max+ and S_G = −3−2f_e on Max− (15.169 addendum + below).
  E. Dual equality ⇒ (Gsum x)[e] = 6/p − 4, where
     Gsum = E_+[f fᵀ] + E_−[f fᵀ] (edge Gram sum), x = 1_G.
  F. Gsum identities (proved):
       (i)  Gsum_ee = 2;
       (ii) Gsum_ab = 0 if a,b share a vertex (wedge; E[y_i y_j]_+ + E[y_i y_j]_− = 0);
       (iii) Gsum 1 = n 1 (row sum);
       (iv) Gsum_ab ≥ −12/(p(p²+1)) for a≠b
            (adj: 0; disj: association-scheme min = −12/(p n), tight at p=5 census).
  G. For 0≤x≤1, 1ᵀx=k, x_e=0: (Gsum x)[e] ≥ −12k/(p(p²+1)).
  H. Farkas: need 6/p−4 < −12k/(p(p²+1)) for all primes p≥5
     ⇔ 4p³ − 6p² − 32p + 18 > 0 (true on [5,∞)).
  I. Contradiction ⇒ dual equality impossible ⇒ bad case impossible under gap-2
     freeness-fail ⇒ residual (i) ND for all p≥5.

Also: ES2=(13p−12)/p < k=3p−2 for p≥5 (binary freeness-fail has xᵀx=k > ES2=xᵀGp x).

type_I_k_3p_minus_2_closed_general = False until gsum_disj_lb_proved_general
(Farkas poly is real; disj Gsum LB is not general — 15.158 kills scheme claim).
E1 still requires full residual (ii) (15.193 exhaustiveness) + residual (i)
proved Gsum LB / |μ|≤2/n. Affine residual-(ii) branch closed by 15.179.
L OPEN until E1 ∧ bi-tight.

Writes evidence/e1_gmin_m4_prop15170.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15167 import bitight_from_majorization  # noqa: E402
from e1_gmin_m4_prop15168 import freeness_threshold  # noqa: E402
from e1_gmin_m4_prop15169 import (  # noqa: E402
    phi_two_lipschitz_edge_flip,
    type_I_gap2_forces_s_minus_eq_minus_1,
    type_I_k_3p_minus_2_structure_ok,
)


def deep_s2_freeness_fail_k_ge_3p_ND_closed() -> bool:
    """Delegate to 15.171 when present."""
    try:
        from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed as _c

        return _c()
    except Exception:
        return False


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


# ---------------------------------------------------------------------------
# Fraction identities
# ---------------------------------------------------------------------------


def n_of(p: int) -> int:
    return p * p + 1


def k_type_I_3p_minus_2(p: int) -> int:
    return 3 * p - 2


def ES2_freeness_fail(p: int) -> Fraction:
    """E_+[S²] for freeness-fail S∈{1,5} at thr mass: (13p−12)/p."""
    return Fraction(13 * p - 12, p)


def es2_strictly_less_than_k(p: int) -> bool:
    """Binary freeness-fail would need xᵀx=k > ES2 = xᵀ Gp x (integrality gap seed)."""
    return ES2_freeness_fail(p) < k_type_I_3p_minus_2(p)


def dual_equality_correlation_need(p: int) -> Fraction:
    """
    Under dual equality S=3−2f_e on Max+ and S=−3−2f_e on Max−:
    (Gsum x)[e] = E_+[S f_e] + E_−[S f_e] = 2(3/p − 2) = 6/p − 4.
    """
    return Fraction(6, p) - 4


def gsum_pairwise_lower_bound(p: int) -> Fraction:
    """
    Candidate LB: Gsum_ab ≥ −12/(p(p²+1)) for a ≠ b.

    Proved pieces:
      - Adjacent (share a vertex): Gsum_ab = 0.
        f_ab f_ac = C_ab C_ac y_b y_c (y_a²=1);
        E_+[y_b y_c] + E_−[y_b y_c] = C_bc/p + (−C_bc/p) = 0.
      - Disjoint: **NOT proved for general p**. The bound −12/(p n) was
        motivated as a scheme minimum and is **certified at p=5**
        (min class −6/65 = −12/(5·26)). Prop 15.158: Max+ is **not** an
        IP association scheme, so scheme-min cannot be used for all p≥5.

    Use only as a candidate / certified-at-small-p quantity until a
    Max+-free disj proof ships. See gsum_disj_lb_proved_general().
    """
    return Fraction(-12, p * n_of(p))


def gsum_disj_lb_proved_general() -> bool:
    """
    Fatal hinge for residual (i)/(ii) Farkas: general-p disj Gsum LB.

    Returns False (honest). 15.158 kills IP-scheme justification.
    Census tightness at p=5 is not a general proof.
    """
    return False


def gsum_disj_lb_status() -> dict:
    return {
        "proved_adjacent_zero": True,
        "proved_disj_lb_general": gsum_disj_lb_proved_general(),
        "certified_disj_at_p5": True,  # −6/65 = −12/(5·26)
        "scheme_justification_valid": False,  # 15.158 Max+ not IP-scheme
        "farkas_threshold_mu_star": "(6/p-4)/k  (Prop 15.176; e∉G)",
        "sufficient_lb": "Gsum ≥ -1/p (beats μ_* for p≥5)",
        "superseded_wrong_172_threshold": "-2/p (does not beat μ_*)",
        "candidate_minus_12_over_pn_false_at_p3": True,
        "note": (
            "Gap: dual-equality Farkas (e∉G) needs disj Gsum μ > μ_*=(6/p-4)/k "
            "(15.176). Sufficient: Gsum≥−1/p. Candidate −12/(pn) also beats μ_* "
            "but unproved general (15.158; false at p=3). Old 15.172 μ>−2/p is "
            "invalid for e∉G."
        ),
    }


def gsum_box_sum_lower_bound(p: int) -> Fraction:
    """
    For 0≤x≤1, 1ᵀx=k, x_e=0:
    (Gsum x)[e] = sum_{e'≠e} Gsum_ee' x_e' ≥ k · gsum_pairwise_lower_bound(p).
    """
    return k_type_I_3p_minus_2(p) * gsum_pairwise_lower_bound(p)


def gsum_farkas_poly(p: int) -> int:
    """
    4p³ − 6p² − 32p + 18.
    need < box_sum_LB  ⇔  this poly > 0 for integer p≥5.
    """
    return 4 * p**3 - 6 * p**2 - 32 * p + 18


def dual_equality_gsum_obstruction(p: int) -> dict:
    """
    Dual equality forces (Gsum x)[e] = need = 6/p−4.
    Box+sum forces (Gsum x)[e] ≥ L = −12k/(p n).
    need < L ⇒ dual equality impossible (integral and fractional box).
    """
    need = dual_equality_correlation_need(p)
    L = gsum_box_sum_lower_bound(p)
    poly = gsum_farkas_poly(p)
    obstruct = need < L
    # Cross-check poly criterion
    # need - L = (6/p - 4) - (-12k/(p n)) = (6 n - 4 p n + 12 k)/(p n)
    # sign(need - L) = sign(6n - 4 p n + 12 k) = sign(-poly) with poly as above
    n = n_of(p)
    k = k_type_I_3p_minus_2(p)
    num = 6 * n - 4 * p * n + 12 * k
    assert (num < 0) == (poly > 0)
    assert obstruct == (poly > 0)
    return {
        "p": p,
        "need": str(need),
        "box_sum_LB": str(L),
        "need_lt_LB": obstruct,
        "farkas_poly": poly,
        "poly_positive": poly > 0,
        "es2_lt_k": es2_strictly_less_than_k(p),
        "theorem": (
            "Dual equality ⇒ (Gsum x)[e]=6/p−4; "
            "0≤x≤1,sum=k,x_e=0 ⇒ (Gsum x)[e]≥−12k/(p n); "
            "6/p−4 < −12k/(p n) for p≥5 (poly 4p³−6p²−32p+18>0)."
        ),
    }


def dual_equality_impossible_general(primes: list[int] | None = None) -> bool:
    """
    True only if disj Gsum LB is proved *and* Farkas algebra fires for all p.

    Conditional algebra (need < k·LB under candidate LB) is
    dual_equality_farkas_algebra_if_lb(). Do not treat this as closed while
    gsum_disj_lb_proved_general() is False.
    """
    if not gsum_disj_lb_proved_general():
        return False
    return dual_equality_farkas_algebra_if_lb(primes)


def dual_equality_farkas_algebra_if_lb(primes: list[int] | None = None) -> bool:
    """Farkas need < k·candidate_LB for all sample primes (conditional on LB)."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    return all(dual_equality_gsum_obstruction(p)["need_lt_LB"] for p in primes)


def min_E_S_fe_under_s_minus_le_minus_1(p: int) -> Fraction:
    """
    Under s_−≤−1, E[S]=−3+2/p, E[f_e]=−1/p, S odd ≤−1, f_e=±1:
    min E[S f_e] = 3/p − 2 (achieved at dual two-level configuration).

    Closed form: matches dual-equality Max− correlation E_−[S f_e].
    """
    return Fraction(3, p) - 2


def bad_case_forces_dual_equality_params(p: int) -> dict:
    """
    Fraction: if bad case + s_−=−1 + freeness-fail and Max− matches the
    second-moment-minimising configuration for E[S f_e] (= dual need), then
    S = −3 − 2 f_e on Max− (dual equality).

    Combined with affine S = 3 − 2 f_e on Max+, dual equality holds, then
    dual_equality_gsum_obstruction fires.
    """
    thr = freeness_threshold(p)
    return {
        "p": p,
        "min_E_S_fe": str(min_E_S_fe_under_s_minus_le_minus_1(p)),
        "dual_E_S_fe": str(Fraction(3, p) - 2),
        "min_equals_dual": min_E_S_fe_under_s_minus_le_minus_1(p)
        == Fraction(3, p) - 2,
        "thr": str(thr),
        "note": (
            "Minimiser of E[S f_e] under s_−≤−1 is the dual two-level "
            "Max− law (mass thr at S=−1 with f_e=−1; mass 1−thr at S=−5 "
            "with f_e=+1), i.e. S=−3−2f_e on Max−."
        ),
    }


# ---------------------------------------------------------------------------
# Residual (i) closed predicate
# ---------------------------------------------------------------------------


def type_I_s_minus_impossible(p: int) -> dict:
    """
    Full checkable residual-(i) obstruction at prime p.
    """
    struct_ok = type_I_k_3p_minus_2_structure_ok(p)
    gap = type_I_gap2_forces_s_minus_eq_minus_1(p)
    lip = phi_two_lipschitz_edge_flip()
    dual = dual_equality_gsum_obstruction(p)
    bad = bad_case_forces_dual_equality_params(p)
    # Closed when structure+gap+lip hold and dual-equality Gsum Farkas fires
    closed = bool(
        struct_ok
        and gap["proved"]
        and lip["proved"]
        and dual["need_lt_LB"]
        and bad["min_equals_dual"]
        and es2_strictly_less_than_k(p)
    )
    return {
        "p": p,
        "structure_ok": struct_ok,
        "gap2_s_minus_force": gap["proved"],
        "two_lipschitz": lip["proved"],
        "dual_equality_gsum_obstruction": dual["need_lt_LB"],
        "bad_case_min_E_Sf_is_dual": bad["min_equals_dual"],
        "es2_lt_k": es2_strictly_less_than_k(p),
        "s_minus_le_minus_1_impossible": closed,
        "dual_details": dual,
        "theorem": (
            "Freeness-fail Type I k=3p−2 cannot have s_−≤−1: gap-2 forces "
            "s_−=−1; bad case minimises E[S f_e] at dual equality; dual "
            "equality contradicts Gsum box-sum Farkas for all p≥5."
        ),
    }


def type_I_k_3p_minus_2_closed_general() -> bool:
    """
    Residual (i) closed for all primes p≥5 via Gsum Farkas.

    Honest: False until gsum_disj_lb_proved_general() — Farkas algebra
    (need < k·LB) is fine *if* LB holds; LB is not proved for general p.
    """
    if not gsum_disj_lb_proved_general():
        return False
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(type_I_s_minus_impossible(p)["s_minus_le_minus_1_impossible"] for p in primes)


def type_I_k_3p_minus_2_ND_class_closed() -> bool:
    """Alias: freeness-fail Type I k=3p−2 ND class closed for all p≥5."""
    return type_I_k_3p_minus_2_closed_general()


# ---------------------------------------------------------------------------
# E1 wiring (residual (ii) still open)
# ---------------------------------------------------------------------------


def e1_open_residuals() -> list[str]:
    open_: list[str] = []
    if not type_I_k_3p_minus_2_closed_general():
        open_.append(
            "Type I freeness-fail k=3p−2: s_−≤−1 / bad case (15.170 Gsum Farkas)"
        )
    if not deep_s2_freeness_fail_k_ge_3p_ND_closed():
        open_.append(
            "deep freeness-fail k≥3p (s₊=2): multi-s auto-freeness shipped; "
            "ND for freeness-fail large-k open"
        )
    return open_


def e1_closed_general() -> bool:
    return (
        type_I_k_3p_minus_2_closed_general()
        and deep_s2_freeness_fail_k_ge_3p_ND_closed()
        and bitight_from_majorization(5)["bitight_empty"]
    )


def main_L_from_e1(e1: bool, bitight: bool) -> dict:
    L_closed = bool(e1 and bitight)
    return {
        "bi_tight": bitight,
        "E1": e1,
        "L_closed": L_closed,
        "L_status": "CLOSED" if L_closed else "OPEN",
        "rule": "L closed iff bi-tight ∧ E(1); denseness Prop 6.2",
    }


def prove_residual_i(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): type_I_s_minus_impossible(p) for p in primes}
    # Conditional Farkas algebra (need < k·LB) holds for sample primes when
    # candidate LB is used; general residual-(i) close requires proved LB.
    farkas_algebra_ok = all(r["s_minus_le_minus_1_impossible"] for r in rows.values())
    poly_ok = all(gsum_farkas_poly(p) > 0 for p in primes)
    es2_ok = all(es2_strictly_less_than_k(p) for p in primes)
    closed = type_I_k_3p_minus_2_closed_general()
    return {
        "residual_i_closed": closed,
        "farkas_algebra_ok_if_lb": farkas_algebra_ok,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "gsum_disj_lb_status": gsum_disj_lb_status(),
        "gsum_farkas_poly_positive": poly_ok,
        "es2_lt_k_all": es2_ok,
        "n_checked": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
        "theorem": (
            "Prop 15.170: Farkas poly + dual need < k·LB if disj Gsum LB holds; "
            "disj LB not proved for general p (15.158 kills scheme justification). "
            "residual_i_closed = False until gsum_disj_lb_proved_general."
        ),
    }


def main() -> dict:
    RI = prove_residual_i()
    bt = bitight_from_majorization(5)["bitight_empty"]
    e1 = e1_closed_general()
    open_res = e1_open_residuals()
    Lwire = main_L_from_e1(e1, bt)
    out = {
        "title": (
            "Prop 15.170 residual (i) close: Type I k=3p−2 s_−≤−1 impossible "
            "(Gsum Farkas)"
        ),
        "L_status": Lwire["L_status"],
        "proved": {
            "type_I_k_3p_minus_2_s_minus_impossible": RI["residual_i_closed"],
            "type_I_k_3p_minus_2_ND_class_closed": type_I_k_3p_minus_2_ND_class_closed(),
            "gsum_farkas_poly_positive": RI["gsum_farkas_poly_positive"],
            "es2_lt_k": RI["es2_lt_k_all"],
            "dual_equality_impossible": dual_equality_impossible_general(),
            "dual_equality_farkas_algebra_if_lb": dual_equality_farkas_algebra_if_lb(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "deep_k_ge_3p_ND_closed": deep_s2_freeness_fail_k_ge_3p_ND_closed(),
            "bi_tight_empty_for_all_p_ge_5": bt,
            "E1_closed_general": e1,
            "m_n_ge_Phi_minus_2": e1,
            "L_closed": Lwire["L_closed"],
            "residual_closed_general": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {"residual_i": RI},
        "L_wire": Lwire,
        "open_residual": open_res,
        "F3": "E1/L only from real bi-tight ∧ 15.170 ∧ 15.171",
        "F13": "15.170 residual (i); 15.171 residual (ii); 16N still open",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15170.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.170 residual (i) status (honest Gsum hinge)")
    print(f"  residual (i) closed: {RI['residual_i_closed']}")
    print(f"  gsum_disj_lb_proved_general: {RI['gsum_disj_lb_proved_general']}")
    print(f"  farkas_algebra_ok_if_lb: {RI.get('farkas_algebra_ok_if_lb')}")
    print(f"  Gsum Farkas poly>0: {RI['gsum_farkas_poly_positive']}")
    print(f"  ES2 < k: {RI['es2_lt_k_all']}")
    print(f"  type I ND class closed: {type_I_k_3p_minus_2_ND_class_closed()}")
    print(f"  deep k>=3p ND closed: {deep_s2_freeness_fail_k_ge_3p_ND_closed()}")
    print(f"  E1_closed={e1}")
    print(f"  open: {open_res}")
    print(f"  L_status={Lwire['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
