#!/usr/bin/env python3
"""Prop. 15.773: residual (ii) at k=5p-1 is empty for all primes p>=29.

The fresh flat branch is closed by the universal common-row identity, not
by classifying mean-2p quadratics. Previous LOCAL classifications carry to
the new edge count; their former all-layer conclusions are not assumed.
See evidence/NOTE_2026-09-04_JOINT_5P_MINUS_ONE_CLOSE.md for the proof.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import baseline_coefficient_rules, residual_even_floor_table
from e1_gmin_m4_prop15751 import (
    density_profile_certificate, height_at_least_two_certificate,
    height_one_junta_certificate,
)
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion
from e1_gmin_m4_prop15768 import p_plus_fifteen_local_exclusion
from e1_gmin_m4_prop15769 import (
    hard_family_catalog as p3_hard_family_catalog,
    p_plus_thirteen_local_exclusion as p3_p_plus_thirteen_local_exclusion,
)
from e1_gmin_m4_prop15770 import (
    p1_p_plus_thirteen_local_exclusion, p1_sharp_family_catalog,
    p3_next_residue_ledger, p_minus_one_local_exclusion,
)
from e1_gmin_m4_prop15772 import (
    hard_family_catalog as p1_hard_family_catalog,
    p1_p_plus_eleven_local_exclusion, p1_third_residue_ledger,
)
from io_atomic import write_json_atomic

ROOT = Path(__file__).resolve().parents[1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int) -> None:
    if (not isinstance(p, int) or isinstance(p, bool) or p < 29
            or not is_prime(p)):
        raise ValueError("need prime p>=29")


def common_row_identity(p: int) -> dict[str, object]:
    """Exact first/second moments establish the identity for either sign."""
    _check_prime(p)
    m = (p + 1) // 2
    one = 2 * Fraction(m, p) - 1
    two = 4 * Fraction(m * (m - 1), p * (p - 1)) - 4 * Fraction(m, p) + 1
    _require(one == Fraction(1, p) and two == Fraction(-1, p),
             "middle-slice moments failed")
    _require(p * p + 1 - 10 * p > 0, "isolated chart unavailable")
    return {"p": p, "H_edge_count": 5 * p,
            "guaranteed_isolated_vertices": p * p - 10 * p + 1,
            "I": 0, "E_z_i": str(one), "E_z_i_z_j": str(two),
            "hard_sign_relative_to_transported_c_H": 1 if p % 4 == 1 else -1,
            "transported_c_H_must_be_recomputed": True,
            "hard_parallel_edge_signed_coefficient": 1,
            "off_fibre_coefficient_sum": "hT-P_L",
            "hard_mean_identity": "a_L=(p+1)P_L-hT-3p",
            "opposite_mean_identity": "a_opp=(p+1)Q+hT-3p",
            "valid_for_both_h_signs": True,
            "equality_classification_or_offset_assumed": False, "proved": True}


def joint_residue_ledger(p: int) -> dict[str, object]:
    """Reuse low-cell theorems, but recompute all new quotient counts."""
    _check_prime(p)
    q, m = (p - 1) // 2, (p + 1) // 2
    prior = p1_third_residue_ledger(p) if p % 4 == 1 else p3_next_residue_ledger(p)
    _require(prior["proved"] and len(prior["rows"]) == m,
             "carried local classification failed")
    minus_one = p_minus_one_local_exclusion(p) if p % 4 == 3 else None
    _require(minus_one is None or minus_one["proved"], "p3 mass p-1 exclusion failed")
    rows = []
    for u in range(q):
        old = prior["rows"][u]
        _require(old["u"] == u and old["forced_low_quotient"] == 1
                 and old["forced_low_mean"] == 2 * u + p + 1,
                 "carried low-cell mean changed")
        live = old["live_rows"]
        if p % 4 == 3 and u == q - 1:
            _require(live == [{"b": 2, "classification": "p_minus_one"},
                              {"b": p - 1, "classification": "p_minus_one"}],
                     "p3 mass p-1 branches were silently discarded")
            live = []
        rows.append({"u": u, "quotient_sum": p - u,
                     "all_quotients_positive": True,
                     "forced_low_quotient": 1, "forced_low_mean": 2 * u + p + 1,
                     "forced_low_count_at_least": u + 1,
                     "candidate_rows": old["candidate_rows"], "live_rows": live})
    zero = prior["rows"][q]
    _require(zero["forced_low_quotient"] == 0 and zero["forced_low_mean"] == p - 1,
             "quotient-zero local classification changed")
    rows.append({"u": q, "quotient_sum": m,
                 "exhaustive_alternatives": ["some_quotient_zero", "all_quotients_one"],
                 "zero_case_low_mean": p - 1, "zero_case_live_rows": zero["live_rows"],
                 "flat_case_low_mean": 2 * p, "flat_case_common_parallel_count": True,
                 "arbitrarily_high_quotients_retained_in_zero_case": True})
    expected = [0, q - 3, q - 2, q - 1] if p % 4 == 1 else [q - 2]
    _require([r["u"] for r in rows[:-1] if r["live_rows"]] == expected,
             "joint residue partition is not exhaustive")
    return {"p": p, "q": q, "m": m, "t": q, "k": 5 * p - 1,
            "surviving_positive_quotient_residues": expected, "rows": rows,
            "quotient_zero_is_not_forced_at_u_equals_q": True, "proved": True}


def normalized_parallel_profile(p: int, u: int, low_quotient: int,
                                low_parallel: int, quotients: list[int]) -> list[int]:
    """Apply the identity to any admissible quotient profile, without a cap."""
    _check_prime(p)
    m = (p + 1) // 2
    if (not 0 <= u < m or low_quotient not in (0, 1)
            or not isinstance(low_parallel, int) or low_parallel < 0
            or len(quotients) != m or any(type(k) is not int or k < 0 for k in quotients)
            or sum(quotients) != p - u or low_quotient not in quotients):
        raise ValueError("invalid common-row quotient profile")
    parallel = [low_parallel + k - low_quotient for k in quotients]
    _require(min(parallel) >= 0 and sum(parallel) <= 5 * p,
             "parallel profile exceeds edge budget")
    hT = (p + 1) * low_parallel - 3 * p - 2 * u - (p + 1) * low_quotient
    _require(all((p + 1) * P - hT - 3 * p == 2 * u + (p + 1) * k
                 for P, k in zip(parallel, quotients))
             and 2 * sum(parallel) == 5 * p + hT, "common-row profile mismatch")
    return parallel


def opposite_mass_exclusion(p: int, shift: int, local: dict[str, object]) -> dict[str, object]:
    """Boundary-zero local theorem plus genuine parity minima on other rows."""
    _check_prime(p)
    floor_certificate = residual_even_floor_table(p)
    _require(floor_certificate["proved"], "phase-zero floor dependency failed")
    floors = floor_certificate["phase_zero_floors"]
    sharp = sharp_integral_quadratic_lift_floor(p)
    mass = p + shift
    nonzero = [[b, floor, mass - floor] for b, floor in floors.items() if b and floor <= mass]
    _require(sharp["proved"] and local["proved"]
             and [b for b, _, _ in nonzero] == [2, p - 1]
             and all(0 < excess < sharp["sharp_scaled_floor"] for _, _, excess in nonzero),
             "opposite all-boundary local exclusion failed")
    return {"scaled_mass": mass, "nonzero_boundary_floor_excess_rows": nonzero,
            "b2_pointwise_baseline": "(x_i-x_j)^2",
            "last_pointwise_baseline": "1-x_j" if p % 4 == 1 else "x_j",
            "half_differences_genuinely_nonnegative_integral": True,
            "b0_is_twice_an_integral_nonnegative_quadratic": True,
            "local_dependency_proved": True, "proved": True}


def _branch_ledger(p: int, u: int, low: int, offset: int,
                   forbidden_Q: int, small: int, shift: int,
                   local: dict[str, object]) -> dict[str, object]:
    q, m, edges = (p - 1) // 2, (p + 1) // 2, 5 * p
    quotient_sum = p - u
    upper = (edges - quotient_sum + m * low) // m
    candidates = [P for P in range(upper + 1) if (P - offset) % q == 0]
    _require(candidates == [offset] and upper <= 9 < q,
             "common low-row coefficient normalization failed")
    hard = m * (offset - low) + quotient_sum
    hT = (p + 1) * offset - 3 * p - 2 * u - (p + 1) * low
    opposite = edges - hard
    next_Q = forbidden_Q + 1
    mean = lambda Q: (p + 1) * Q + hT - 3 * p
    floor_certificate = residual_even_floor_table(p)
    _require(floor_certificate["proved"], "phase-zero floor dependency failed")
    floors = floor_certificate["phase_zero_floors"]
    surplus = opposite - m * next_Q
    checks = {
        "signed_total_agreement": 2 * hard == edges + hT,
        "smaller_Q_negative": mean(forbidden_Q - 1) < 0,
        "forbidden_mass_identity": mean(forbidden_Q) == small,
        "small_mass_below_both_floors": 0 < small < min(p - 3, min(f for b, f in floors.items() if b)),
        "next_mass_identity": mean(next_Q) == p + shift,
        "positive_forcing_count": 0 <= surplus < m,
    }
    _require(all(checks.values()), "carried opposite-row ledger failed")
    exclusion = opposite_mass_exclusion(p, shift, local)
    _require(exclusion["proved"], "carried opposite exclusion dependency failed")
    return {"u": u, "low_quotient": low, "low_scaled_mean": 2 * u + (p + 1) * low,
            "coefficient_offset": offset, "common_low_parallel_candidates": candidates,
            "common_low_parallel_upper_bound": upper,
            "every_hard_parallel_formula": f"P_L={offset-low}+k_L",
            "quotient_sum": quotient_sum, "quotient_excess_above_one": quotient_sum - m,
            "hard_edge_count": hard, "hard_sign_times_global_T": hT,
            "opposite_edge_count": opposite, "forbidden_Q": forbidden_Q,
            "forbidden_scaled_mean": small, "forced_next_Q": next_Q,
            "forced_next_scaled_mean": p + shift, "surplus": surplus,
            "forced_next_row_count_at_least": m - surplus,
            "opposite_local_exclusion": exclusion, "checks": checks, "proved": True}


def flat_mean_2p_exclusion(p: int) -> dict[str, object]:
    """No offset congruence or mean-2p equality classification is used."""
    _check_prime(p)
    m = (p + 1) // 2
    identity = common_row_identity(p)
    local = opposite_mass_exclusion(p, 9, p_plus_nine_local_exclusion(p))
    _require(identity["proved"] and local["proved"], "flat branch dependency failed")
    _require(5 * p // m == 9 and 7 - p < 0 and 8 < p - 3,
             "flat R=P+Q threshold failed")
    rows = []
    for P in range(10):
        hard, hT, next_Q = m * P, (p + 1) * P - 5 * p, 9 - P
        opposite = 5 * p - hard
        surplus = opposite - m * next_Q
        _require(2 * hard == 5 * p + hT and next_Q >= 0
                 and surplus == m - 5 and 0 <= surplus < m
                 and (p + 1) * next_Q + hT - 3 * p == p + 9,
                 "flat opposite-row pigeonhole failed")
        rows.append({"P": P, "hard_edge_count": hard, "hard_sign_times_global_T": hT,
                     "opposite_edge_count": opposite, "formal_forbidden_Q": 8 - P,
                     "forbidden_Q_is_in_domain": P <= 8,
                     "forced_next_Q": next_Q, "forced_next_scaled_mean": p + 9,
                     "surplus": surplus, "forced_next_row_count_at_least": 5, "proved": True})
    return {"p": p, "hard_mean": 2 * p, "all_quotients": 1,
            "common_parallel_range": [0, 9], "common_row_identity": identity,
            "mean_in_R_equals_P_plus_Q": "(p+1)R-8p", "minimum_R": 9,
            "R_at_most_seven_negative": True, "R_eight_mass_eight_excluded": True,
            "P9_uses_Q_nonnegative_not_a_negative_index_row": True,
            "new_equality_classification_used": False, "coefficient_offset_assumed": False,
            "parallel_cases": rows, "opposite_local_exclusion": local, "proved": True}


def joint_layer_exclusion(p: int) -> dict[str, object]:
    _check_prime(p)
    q, m = (p - 1) // 2, (p + 1) // 2
    identity, residues = common_row_identity(p), joint_residue_ledger(p)
    _require(identity["proved"] and residues["proved"], "joint identity or residue dependency failed")
    baselines = baseline_coefficient_rules(p)
    _require(baselines["proved"], "baseline coefficient theorem failed")
    p9 = p_plus_nine_local_exclusion(p)
    if p % 4 == 1:
        new_hard, sharp = p1_hard_family_catalog(p), p1_sharp_family_catalog(p)
        _require(new_hard["proved"] and new_hard["coefficient_offsets"] == [4, 6]
                 and sharp["proved"]
                 and [r["coefficient_offset"] for r in sharp["families"]] == [3, 5],
                 "carried p1 hard offsets failed")
        dependencies = {7: {"proved": all(r["proved"] for r in (
            height_at_least_two_certificate(p), height_one_junta_certificate(p),
            density_profile_certificate(p)))}, 9: p9,
            11: p1_p_plus_eleven_local_exclusion(p),
            13: p1_p_plus_thirteen_local_exclusion(p), 15: p_plus_fifteen_local_exclusion(p)}
        specs = [("complement_literal", 0, 1, 5, 2, 6, 7),
                 ("complement_triple", q - 3, 1, 2, 6, 14, 15),
                 ("XNOR_sharp_P3", q - 2, 1, 3, 5, 12, 13),
                 ("XNOR_sharp_P5", q - 2, 1, 5, 3, 12, 13),
                 ("gap_four_or_literal_P4", q - 1, 1, 4, 4, 10, 11),
                 ("literal_sharp_P6", q - 1, 1, 6, 2, 10, 11),
                 ("zero_quotient_XNOR", q, 0, 4, 3, 8, 9)]
        counts = [4, 8, 7, 7, 6, 6, 5]
    else:
        sharp = p3_hard_family_catalog(p)
        _require(sharp["proved"] and sorted(r["coefficient_offset"] for r in sharp["families"])
                 == [2, 3, 4, 5], "carried p3 sharp offsets failed")
        dependencies = {9: p9, 13: p3_p_plus_thirteen_local_exclusion(p)}
        specs = [(f"sharp_P{P}", q - 2, 1, P, 8 - P, 12, 13) for P in (2, 3, 4, 5)]
        specs += [("zero_quotient_XNOR", q, 0, 4, 3, 8, 9),
                  ("zero_quotient_literal", q, 0, 3, 4, 8, 9)]
        counts = [7, 7, 7, 7, 5, 5]
    branches = {name: _branch_ledger(p, u, low, P, Q, small, shift, dependencies[shift])
                for name, u, low, P, Q, small, shift in specs}
    _require([r["forced_next_row_count_at_least"] for r in branches.values()] == counts,
             "carried forcing counts changed")
    flat = flat_mean_2p_exclusion(p)
    _require(flat["proved"] and all(row["proved"] for row in branches.values()),
             "joint branch dependency failed")
    return {"p": p, "p_mod_4": p % 4, "m": m, "layer_index_t": q,
            "original_k": 5 * p - 1, "H_edge_count": 5 * p,
            "common_row_identity": identity, "residue_ledger": residues,
            "carried_branch_exclusions": branches, "flat_branch_exclusion": flat,
            "all_boundary_sizes_excluded": True, "residual_ii_layer_excluded": True,
            "finite_prime_graph_or_slice_census_used": False, "proved": True}


def proposition_15773() -> dict[str, object]:
    records = [joint_layer_exclusion(p) for p in (29, 31, 37, 43)]
    return {"prop": "15.773", "status": "PROVED_INFINITE_FAMILY",
            "statement": "For every prime p>=29, residual(ii) at k=5p-1 is empty",
            "proof_note": "evidence/NOTE_2026-09-04_JOINT_5P_MINUS_ONE_CLOSE.md",
            "symbolic_range": {"minimum_prime": 29, "p_mod_4": [1, 3],
                               "t": "(p-1)/2", "k": "5p-1"},
            "new_generic_frontier": "p>=29,t>=(p+1)/2,k>=5p+1",
            "records_are_identity_replays_not_exhaustive_prime_evidence": True,
            "new_mean_2p_equality_classification_used": False, "records": records,
            "proved": all(r["proved"] for r in records),
            "residual_ii_closed_general": False, "e1_closed_general": False,
            "original_MO_limit_closed": False}


if __name__ == "__main__":
    result = proposition_15773()
    write_json_atomic(ROOT / "evidence/e1_gmin_m4_prop15773.json", result)
    print(json.dumps({"prop": "15.773", "proved": result["proved"]}, sort_keys=True))
