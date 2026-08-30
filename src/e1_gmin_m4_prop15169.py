#!/usr/bin/env python3
"""
Prop 15.169 — Type I freeness-fail k=3p−2 ND reduction + deep multi-s auto-freeness
(honest; full E1 still OPEN until residual (i) final step and residual (ii) close).

REAL (checkable Fraction / prior props):

TYPE I k=3p−2 freeness-fail (residual (i) of 15.168):
  A. Structure at freeness equality: k=3p−2, S∈{1,5}, N1/N=thr=(p+1)/(2p),
     affine S+2f_e=3 on Max+ (15.43.3 parallel boundary). Fraction identities.
  B. H:=G∪e has |H|=3p−1, S_H∈{2,4} on Max+, min S_H=2.
  C. Phi is 2-Lipschitz under single-edge flip: each Q_y changes by ±2, so
     |Phi(A⊕f)−Phi(A)|≤2. Hence if Phi(C⊕G)≥Phi then Phi(C⊕G⊕e)≥Phi−2 (ND).
  D. Type I ⇒ Phi≥Phi−2 (15.42.1). Gap-2 undercutter (Phi=Phi−2) requires
     s_−≤−1. With k odd and E_−[S]=−k/p=−3+2/p>−3, odd scores force s_−=−1.
  E. When s_−=−1: U_−={S_G=−1} nonempty. For z∈U_−, Q after adding e:
     f_e(z)=+1 ⇒ |Q_H(z)|=Phi (strong ND); f_e≡−1 on U_− ⇒ s_−^H≤−2
     (dichotomy only gives Phi(H)≥Phi−4).
  F. (OPEN hinge, Prop 15.170/172) Bad case / s_−≤−1 Farkas is conditional on
     proved disj Gsum LB (μ > −2/p). Candidate −12/(pn) algebra is real but
     gsum_disj_lb_proved_general()=False (15.158; candidate false at p=3).

  Consequence for residual (i): structure+gap+lip shipped; ND class OPEN
  until gsum_disj_lb_proved_general.

DEEP multi-s auto-freeness (partial residual (ii)):
  G. For s₊=s≥2 with scores ≥s same parity, N_s/N lb = (s+2 − k/p)/2
     (from E[S]≥s·a+(s+2)(1−a)). Auto-freeness when lb>thr:
     k ≤ p(s+1)−2. Recovers s=2 ⇒ k≤3p−2 (15.168.F).
  H. Deep freeness-fail only possible for k ≥ p(s+1)−1 at that s.
  I. Deep freeness-fail k≥3p for s₊=2: affine dual two-level CLOSED (15.179);
     full residual (ii) OPEN until exhaustiveness (15.193) — multi-level /
     non-affine freeness-fail not reduced to freeze.

E1_closed_general = False until (i) F and (ii) full close (exhaustiveness or
separate ND). residual_closed_general = False. L OPEN.
Writes evidence/e1_gmin_m4_prop15169.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15720 import (  # noqa: E402
    required_bitight_levels_empty_all_primes,
)
from e1_gmin_m4_prop15168 import (  # noqa: E402
    freeness_threshold,
    tight_cover_obstruction_applicable,
    type_I_fail_k_2p_minus_1_ND,
    deep_fail_k_3p_minus_1_impossible,
    deep_s2_auto_freeness,
    k_max_auto_freeness_s2,
    e1_closed_general as e1_closed_15168,
)


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
# Type I freeness-fail k=3p−2 structure (Fraction)
# ---------------------------------------------------------------------------


def type_I_k_3p_minus_2_params(p: int) -> dict:
    """Fraction structure of freeness-fail equality k=3p−2 / S∈{1,5}."""
    thr = freeness_threshold(p)
    k = 3 * p - 2
    # If S only in {1,5}: a + 5(1-a) = k/p ⇒ a = (5 - k/p)/4
    a = (Fraction(5) - Fraction(k, p)) / 4
    # Affine S+2f_e=3: on S=1 ⇒ f_e=1; on S=5 ⇒ f_e=-1
    # H scores: S_H = S+f_e ∈ {2,4}
    kH = k + 1  # 3p-1
    ES = Fraction(k, p)
    ES_H = Fraction(kH, p)
    return {
        "p": p,
        "k": k,
        "k_H": kH,
        "freeness_thr": str(thr),
        "N1_N_at_equality": str(a),
        "matches_thr": a == thr,
        "E_S_Max_plus": str(ES),
        "E_S_H_Max_plus": str(ES_H),
        "affine": "S+2f_e=3 on Max+",
        "S_support": [1, 5],
        "S_H_support": [2, 4],
        "s_plus_H": 2,
        "mean_Max_minus_S": str(-ES),
        "mean_Max_minus_S_H": str(-ES_H),
        "k_odd": k % 2 == 1,
        "k_H_even": kH % 2 == 0,
    }


def type_I_k_3p_minus_2_structure_ok(p: int) -> bool:
    r = type_I_k_3p_minus_2_params(p)
    return bool(r["matches_thr"] and r["k_odd"] and r["k_H_even"] and r["s_plus_H"] == 2)


def phi_two_lipschitz_edge_flip() -> dict:
    """
    Phi is 2-Lipschitz under single edge flip (proved).

    For any zero-diag Seidel A and edge e: Q_y(A⊕e) = Q_y(A) − 2 f_e^A(y)
    with f_e=±1, so each |Q_y| changes by at most 2 in the value of Q_y,
    and therefore max_y |Q_y| changes by at most 2.
    """
    return {
        "proved": True,
        "lipschitz_constant": 2,
        "theorem": (
            "Q_y(A⊕e)=Q_y(A)−2f_e(y) with f_e=±1 ⇒ |Phi(A⊕e)−Phi(A)|≤2."
        ),
        "corollary_ND_if_Phi_G_ge_Phi": (
            "If Phi(C⊕G)≥Phi then Phi(C⊕G⊕e)≥Phi−2 (weak ND for every e)."
        ),
    }


def type_I_gap2_forces_s_minus_eq_minus_1(p: int) -> dict:
    """
    Type I gap-2 undercutter with k=3p−2 (odd) forces s_−=−1 (proved Fraction).
    """
    k = 3 * p - 2
    mean_m = -Fraction(k, p)  # E_-[S]
    # mean > -3: -3+2/p > -3
    mean_gt_minus_3 = mean_m > -3
    # odd scores (k odd): values in {...,-5,-3,-1,1,...}
    # if max <= -3: mean <= -3, contradict mean_gt_minus_3
    # if max >= 0: dichotomy Phi >= Phi, not gap-2 undercutter
    # gap-2 Type I undercutter: Phi=Phi-2 requires s_- <= -1 (else Phi>=Phi)
    # combined with max >= -1 from mean: s_- = -1
    return {
        "p": p,
        "k": k,
        "mean_Max_minus": str(mean_m),
        "mean_gt_minus_3": mean_gt_minus_3,
        "k_odd": k % 2 == 1,
        "s_minus_forced_if_gap2_undercutter": -1,
        "proved": bool(mean_gt_minus_3 and k % 2 == 1),
        "theorem": (
            "k odd ⇒ odd scores; E_-[S]=−3+2/p>−3 ⇒ s_−≥−1; "
            "gap-2 Type I undercutter ⇒ s_−≤−1; hence s_−=−1."
        ),
    }




def bad_case_dual_two_level_H_params(p: int) -> dict:
    """
    Fraction identities IF bad case holds: s_-=-1, Max- two-level {-1,-3},
    f_e≡-1 on U_- (mass 1/p), freeness-fail affine on Max+.

    Then H=G∪e has dual two-level:
      Max+: S_H ∈ {2,4} mass thr at 2
      Max-: S_H ∈ {-2,-4} mass thr at -2
    |H|=3p-1, E[S_H^2]_+ = E[S_H^2]_- = 10 - 6/p (checkable).
    """
    thr = freeness_threshold(p)
    k = 3 * p - 2
    kH = 3 * p - 1
    # Max- two-level for G: mass b=1/p at S=-1
    b = Fraction(1, p)
    # after bad f_e: mass at S_H=-2 equals thr
    mass_m2 = b + (1 - b) / 2  # U_- all go to -2; half of rest with f=+1 go to -2
    # Simplify: b + (1-b)/2 = (1+b)/2 = (p+1)/(2p) = thr
    mass_m2_exact = (1 + b) / 2
    ES2_H = 4 * thr + 16 * (1 - thr)  # same formula as Max+ for {2,4}
    return {
        "p": p,
        "mass_at_S_H_minus_2": str(mass_m2_exact),
        "equals_thr": mass_m2_exact == thr,
        "E_S2_H": str(ES2_H),
        "E_S2_H_simplified": str(Fraction(10) - Fraction(6, p)),
        "ES2_match": ES2_H == Fraction(10) - Fraction(6, p),
        "k_H": kH,
        "bi_tight_level2_size": 2 * p,
        "size_defect": kH - 2 * p,  # p-1
        "theorem": (
            "Bad case + two-level Max- ⇒ dual two-level H with mass thr at ±2, "
            "E[S_H^2]=10-6/p, |H|=3p-1=2p+(p-1). Bi-tight level-2 would need |H|=2p "
            "and E[S^2]=4; defect p-1 edges and score leakage to ±4."
        ),
    }


def sum_pairs_Gplus_from_ES2(p: int) -> Fraction:
    """sum_{e'<e'' in G} G+_{e'e''} = (E_+[S^2] - k)/2 (exact for freeness-fail)."""
    k = 3 * p - 2
    ES2 = Fraction(13 * p - 12, p)
    return (ES2 - k) / 2


def sum_pairs_Gplus_negative(p: int) -> bool:
    """For p≥5, sum of G+ over pairs in G is negative (checkable)."""
    # (ES2 - k)/2 = -3(p-1)(p-4)/(2p) for p≥5
    return sum_pairs_Gplus_from_ES2(p) < 0


def type_I_k_3p_minus_2_ND_reduction(p: int) -> dict:
    """
    Residual-(i) reduction + conditional 15.170 Gsum Farkas.

    ND_for_this_class / bad_case_impossible_general_p are True only when
    structure+gap+lip hold **and** gsum_disj_lb_proved_general() (15.170 hinge)
    **and** the dual-equality algebra fires. Candidate-LB algebra alone is not enough.
    """
    struct = type_I_k_3p_minus_2_params(p)
    gap = type_I_gap2_forces_s_minus_eq_minus_1(p)
    lip = phi_two_lipschitz_edge_flip()
    bad_tl = bad_case_dual_two_level_H_params(p)
    sum_gp = sum_pairs_Gplus_from_ES2(p)
    # Prop 15.170: dual-equality Gsum Farkas (lazy import avoids cycles)
    try:
        from e1_gmin_m4_prop15170 import (
            dual_equality_gsum_obstruction,
            dual_equality_farkas_algebra_if_lb,
            gsum_disj_lb_proved_general,
        )

        dual = dual_equality_gsum_obstruction(p)
        lb_proved = bool(gsum_disj_lb_proved_general())
        farkas_algebra = bool(dual["need_lt_LB"])
        # Conditional algebra under candidate LB (not a close)
        dual["farkas_algebra_if_lb"] = dual_equality_farkas_algebra_if_lb([p])
        dual["gsum_disj_lb_proved_general"] = lb_proved
    except Exception:
        dual = {
            "need_lt_LB": False,
            "farkas_algebra_if_lb": False,
            "gsum_disj_lb_proved_general": False,
            "error": "prop15170 import failed",
        }
        lb_proved = False
        farkas_algebra = False
    # Honest: bad case impossible only with proved LB + algebra
    bad_case_impossible_general = bool(lb_proved and farkas_algebra)
    ND_class_closed = bool(
        struct["matches_thr"]
        and gap["proved"]
        and lip["proved"]
        and bad_case_impossible_general
    )
    return {
        "p": p,
        "structure_ok": type_I_k_3p_minus_2_structure_ok(p),
        "structure": struct,
        "gap2_forces_s_minus": gap,
        "two_lipschitz": lip,
        "bad_case_dual_two_level": bad_tl,
        "sum_pairs_Gplus": str(sum_gp),
        "sum_pairs_Gplus_negative": sum_gp < 0,
        "farkas_algebra_if_lb": farkas_algebra,
        "gsum_disj_lb_proved_general": lb_proved,
        "bad_case_impossible_general_p": bad_case_impossible_general,
        "dual_equality_gsum_obstruction": dual,
        "ND_for_this_class": ND_class_closed,
        "open_step": (
            ""
            if ND_class_closed
            else (
                "Prove disj Gsum μ > −2/p Max+-free (15.172) so dual-equality "
                "Farkas closes residual (i); gsum_disj_lb_proved_general=False."
            )
        ),
        "certified_p5_milp_s_minus_le_m1_infeasible": True,
        "note_status": (
            "Prop 15.170/172: Farkas algebra under candidate LB is real; "
            "residual (i) ND class OPEN until gsum_disj_lb_proved_general."
        ),
    }


def type_I_k_3p_minus_2_closed_general() -> bool:
    """True when Prop 15.170 Gsum Farkas closes residual (i) for all p≥5."""
    try:
        from e1_gmin_m4_prop15170 import type_I_k_3p_minus_2_closed_general as _c170

        return _c170()
    except Exception:
        primes = [p for p in range(5, 40) if is_prime(p)]
        return all(
            type_I_k_3p_minus_2_ND_reduction(p)["ND_for_this_class"] for p in primes
        )


# ---------------------------------------------------------------------------
# Deep multi-s auto-freeness (extends 15.168.F)
# ---------------------------------------------------------------------------


def deep_s_freeness_lb(p: int, s: int, k: int) -> Fraction:
    """
    For scores ≥ s with step-2 (same parity), a=N_s/N:
    E[S] ≥ s·a + (s+2)(1−a) = s+2−2a ⇒ a ≥ (s+2 − k/p)/2.
    """
    return (Fraction(s + 2) - Fraction(k, p)) / 2


def k_max_auto_freeness_s(p: int, s: int) -> int:
    """Largest k with auto-freeness at min level s: p(s+1)−2."""
    return p * (s + 1) - 2


def deep_s_auto_freeness(p: int, s: int, k: int) -> dict:
    thr = freeness_threshold(p)
    lb = deep_s_freeness_lb(p, s, k)
    return {
        "p": p,
        "s": s,
        "k": k,
        "N_s_N_lb": str(lb),
        "freeness_thr": str(thr),
        "auto_freeness": lb > thr,
        "k_max_auto": k_max_auto_freeness_s(p, s),
    }


def deep_multi_s_auto_freeness_table(p: int) -> dict:
    """Auto-freeness k-bounds for s=2,3,4,5."""
    rows = {}
    for s in (2, 3, 4, 5):
        kmax = k_max_auto_freeness_s(p, s)
        at = deep_s_auto_freeness(p, s, kmax)
        past = deep_s_auto_freeness(p, s, kmax + 1)
        rows[str(s)] = {
            "k_max_auto": kmax,
            "auto_at_kmax": at["auto_freeness"],
            "not_auto_at_kmax_plus_1": not past["auto_freeness"],
            "ok": at["auto_freeness"] and not past["auto_freeness"],
        }
    return {
        "p": p,
        "by_s": rows,
        "all_boundaries_ok": all(r["ok"] for r in rows.values()),
        "theorem": (
            "For min-level s≥2, scores ≥s step 2: auto-freeness for k≤p(s+1)−2 "
            "(Fraction lb>thr). s=2 recovers k≤3p−2."
        ),
    }


def deep_freeness_fail_k_ge_3p_open() -> bool:
    """Residual (ii): open only if 15.171 has not closed it."""
    try:
        from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed as _c171

        return not _c171()
    except Exception:
        return True


def deep_s2_freeness_fail_k_ge_3p_ND_closed() -> bool:
    """True when Prop 15.171 dual two-level Gsum Farkas + freeness chain closes residual (ii)."""
    try:
        from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed as _c171

        return _c171()
    except Exception:
        return False


# ---------------------------------------------------------------------------
# E1 status after 15.169 (honest)
# ---------------------------------------------------------------------------


def e1_open_residuals() -> list[str]:
    open_ = []
    if not type_I_k_3p_minus_2_closed_general():
        open_.append(
            "Type I freeness-fail k=3p−2: s_−≤−1 / bad case "
            "(structure+2-Lipschitz+gap2 shipped; 15.170 Gsum Farkas open/fail)"
        )
    if deep_freeness_fail_k_ge_3p_open():
        open_.append(
            "deep freeness-fail k≥3p (s₊=2): multi-s auto-freeness shipped; "
            "ND for freeness-fail large-k open"
        )
    return open_


def e1_closed_general() -> bool:
    """Both residuals must close — no soft-close."""
    return (
        type_I_k_3p_minus_2_closed_general()
        and deep_s2_freeness_fail_k_ge_3p_ND_closed()
        and required_bitight_levels_empty_all_primes()
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


def prove_type_I_3p_minus_2(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {str(p): type_I_k_3p_minus_2_ND_reduction(p) for p in primes}
    struct_ok = all(r["structure_ok"] for r in rows.values())
    gap_ok = all(r["gap2_forces_s_minus"]["proved"] for r in rows.values())
    return {
        "structure_identities_ok": struct_ok,
        "gap2_s_minus_force_ok": gap_ok,
        "two_lipschitz_ok": True,
        "ND_class_closed_general": type_I_k_3p_minus_2_closed_general(),
        "n_checked": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
        "open_step": rows[str(primes[0])]["open_step"],
    }


def prove_deep_multi_s(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {str(p): deep_multi_s_auto_freeness_table(p) for p in primes}
    ok = all(r["all_boundaries_ok"] for r in rows.values())
    # s=2 matches 15.168
    s2_ok = all(
        deep_multi_s_auto_freeness_table(p)["by_s"]["2"]["k_max_auto"] == 3 * p - 2
        for p in primes
    )
    return {
        "multi_s_auto_freeness_ok": ok,
        "s2_matches_15168": s2_ok,
        "deep_k_ge_3p_ND_closed": deep_s2_freeness_fail_k_ge_3p_ND_closed(),
        "n_checked": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:3]},
    }


def main() -> dict:
    TI = prove_type_I_3p_minus_2()
    DP = prove_deep_multi_s()
    # Keep prior 15.168 predicates true
    prior_k2pm1 = all(
        type_I_fail_k_2p_minus_1_ND(p)["ND_for_this_class"]
        for p in (5, 7, 11, 13)
    )
    prior_k3pm1 = all(
        deep_fail_k_3p_minus_1_impossible(p)["impossible_when_bitight_empty"]
        for p in (5, 7, 11, 13)
    )
    bt = required_bitight_levels_empty_all_primes()
    e1 = e1_closed_general()
    open_res = e1_open_residuals()
    Lwire = main_L_from_e1(e1, bt)
    ti_closed = type_I_k_3p_minus_2_closed_general()
    out = {
        "title": (
            "Prop 15.169 Type I k=3p−2 reduction + deep multi-s auto-freeness; "
            "residuals via 15.170–171"
        ),
        "L_status": Lwire["L_status"],
        "proved": {
            "type_I_k_3p_minus_2_structure": TI["structure_identities_ok"],
            "type_I_gap2_s_minus_force": TI["gap2_s_minus_force_ok"],
            "phi_two_lipschitz": TI["two_lipschitz_ok"],
            "type_I_k_3p_minus_2_ND_class_closed": ti_closed,
            "deep_multi_s_auto_freeness": DP["multi_s_auto_freeness_ok"],
            "deep_s2_matches_15168": DP["s2_matches_15168"],
            "deep_k_ge_3p_ND_closed": DP["deep_k_ge_3p_ND_closed"],
            "prior_type_I_k_2p_minus_1_ND": prior_k2pm1,
            "prior_deep_k_3p_minus_1_impossible": prior_k3pm1,
            "bi_tight_empty_for_all_p_ge_5": bt,
            "E1_closed_general": e1,
            "m_n_ge_Phi_minus_2": e1,
            "L_closed": Lwire["L_closed"],
            "residual_closed_general": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {"type_I_3pm2": TI, "deep_multi_s": DP},
        "L_wire": Lwire,
        "open_residual": open_res,
        "F3": "E1/L only from real bi-tight ∧ 15.170 ∧ 15.171",
        "F13": "15.169–171 E1 ND; residual/16N still open",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15169.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.169 (15.170–171 residual chain)")
    print(f"  type I structure: {TI['structure_identities_ok']}")
    print(f"  gap2 s_- force: {TI['gap2_s_minus_force_ok']}")
    print(f"  2-Lipschitz: {TI['two_lipschitz_ok']}")
    print(f"  type I k=3p-2 ND class closed: {ti_closed}")
    print(f"  deep multi-s auto: {DP['multi_s_auto_freeness_ok']}")
    print(f"  deep k>=3p ND closed: {DP['deep_k_ge_3p_ND_closed']}")
    print(f"  E1_closed={e1}")
    print(f"  open: {open_res}")
    print(f"  L_status={Lwire['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
