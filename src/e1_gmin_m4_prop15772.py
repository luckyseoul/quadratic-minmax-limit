#!/usr/bin/env python3
"""Proposition 15.772: p=1 mod4 third post-band layer, k=5p-3.

The proof is uniform for primes p>=29.  The new ingredient is the
punctured complement-triple gap theorem: gap two is impossible and gap
four has precisely three pair-plus-complement-literal forms, all of
coefficient offset four.  Existing sharp Boolean lifts and the new
mass-(p-1), mass-(p+11) exclusions finish the common-row ledger.

Numerical records below replay exact identities, not a prime census or a
substitute for the proof in NOTE_2026-09-04_P1_THIRD_POST_BAND_CLOSE.md.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_complement_triple_gap import (
    complement_triple_gap_certificate,
    p1_p_minus_one_local_exclusion,
    p1_p_plus_eleven_local_exclusion,
)
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import residual_even_floor_table
from e1_gmin_m4_prop15751 import (
    density_profile_certificate,
    height_at_least_two_certificate,
    height_one_junta_certificate,
)
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion
from e1_gmin_m4_prop15768 import p_plus_fifteen_local_exclusion
from e1_gmin_m4_prop15770 import (
    p1_p_plus_thirteen_local_exclusion,
    p1_sharp_family_catalog,
    sharp_p_minus_three_classification_all_odd,
)
from io_atomic import write_json_atomic

ROOT = Path(__file__).resolve().parents[1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int) -> None:
    if (not isinstance(p, int) or isinstance(p, bool) or p < 29
            or p % 4 != 1 or not is_prime(p)):
        raise ValueError("need prime p>=29 with p=1 mod4")


def p1_third_residue_ledger(p: int) -> dict[str, object]:
    """All residues; only genuine parity-minimum differences use 15.688."""
    _check_prime(p)
    q, m = (p - 1) // 2, (p + 1) // 2
    floors = residual_even_floor_table(p)["phase_one_floors"]
    gap = complement_triple_gap_certificate(p)
    _require(
        gap["proved"] and gap["excess_two_excluded"]
        and not gap["excess_four_excluded"]
        and gap["allowed_excesses_in_zero_to_four"] == [0, 4]
        and gap["excess_four_coefficient_offsets"] == [4],
        "punctured complement-triple gap or equality classification failed",
    )
    minus_one = p1_p_minus_one_local_exclusion(p)
    expected = {
        0: [(p - 1, "exact_complement_literal")],
        q - 3: [(p - 3, "exact_complement_triple")],
        q - 2: [(2, "XNOR_sharp_lift")],
        q - 1: [(p - 3, "triple_gap_four"), (p - 1, "literal_sharp_lift")],
        q: [(2, "exact_XNOR")],
    }
    rows = []
    for u in range(m):
        quotient_sum = 2 * q - u
        low_k = 1 if u < q else 0
        low_mean = 2 * u + (p + 1) * low_k
        low_count = 2 * m - quotient_sum if low_k else m - quotient_sum
        candidates, live = [], []
        for boundary, floor in floors.items():
            if floor > low_mean:
                continue
            excess = low_mean - floor
            if boundary == p - 3:
                # Unlike XNOR/literal, (r-2)^2 is not parity-minimal at r=0.
                _require(excess in (0, 2, 4) and gap["proved"],
                         "punctured complement-triple branch lacks its theorem")
                classification = {0: "exact_complement_triple",
                                  2: "excluded_punctured_gap_two",
                                  4: "triple_gap_four"}[excess]
            elif boundary in (2, p - 1):
                if excess == 0:
                    classification = "exact_XNOR" if boundary == 2 else "exact_complement_literal"
                elif 0 < excess < p - 3:
                    classification = "excluded_genuine_subsharp_lift"
                elif excess == p - 3:
                    classification = "XNOR_sharp_lift" if boundary == 2 else "literal_sharp_lift"
                elif boundary == 2 and excess == p - 1:
                    _require(minus_one["proved"], "mass p-1 exclusion missing")
                    classification = "excluded_mass_p_minus_one"
                else:
                    raise ArithmeticError("unclassified parity-minimum lift")
            else:
                raise ArithmeticError("unexpected middle boundary below 2p")
            candidates.append({"b": boundary, "floor": floor, "excess": excess,
                               "classification": classification})
            if not classification.startswith("excluded_"):
                live.append((boundary, classification))
        _require(live == expected.get(u, []) and low_count > 0,
                 f"residue u={u} has an unclassified low cell")
        rows.append({"u": u, "quotient_sum": quotient_sum,
                     "forced_low_quotient": low_k, "forced_low_mean": low_mean,
                     "forced_low_count_at_least": low_count,
                     "candidate_rows": candidates,
                     "live_rows": [{"b": b, "classification": c} for b, c in live]})
    return {"p": p, "q": q, "m": m, "t": q - 1, "k": 5 * p - 3,
            "H_edge_count": 5 * p - 2,
            "guaranteed_isolated_vertices": p * p - 10 * p + 5,
            "all_positive_quotients_when_u_less_than_q": True,
            "quotient_zero_forced_when_u_equals_q": True,
            "mean_2p_high_rows_need_no_separate_classification": True,
            "surviving_residues": sorted(expected), "rows": rows,
            "proved": True}


def hard_family_catalog(p: int) -> dict[str, object]:
    """The new all-low mean 2p-2 families, allowing overlaps of supports."""
    _check_prime(p)
    gap = complement_triple_gap_certificate(p)
    sharp = sharp_p_minus_three_classification_all_odd(p)
    excluded = p1_p_minus_one_local_exclusion(p)
    families = [
        {"b": p - 3, "family": "pair_plus_complement_literal",
         "signed_target": "5+z_i*z_j-z_k", "coefficient_offset": 4},
        {"b": p - 1, "family": "literal_plus_omitted_pair",
         "signed_target": "(4+z_a)+(1-z_i-z_j+z_i*z_j)",
         "coefficient_offset": 4},
        {"b": p - 1, "family": "literal_plus_all_equal_triple",
         "signed_target": "(4+z_a)+(1+z_i*z_j+z_i*z_k+z_j*z_k)",
         "coefficient_offset": 6},
    ]
    _require(gap["proved"] and not gap["excess_four_excluded"]
             and gap["excess_four_coefficient_offsets"] == [4]
             and len(gap["excess_four_forms"]) == 3
             and sharp["proved"] and excluded["proved"],
             "new hard family dependency failed")
    return {"p": p, "scaled_mean": 2 * p - 2, "families": families,
            "coefficient_offsets": [4, 6],
            "b2_mass_p_minus_one_excluded": True,
            "arbitrary_support_overlap_allowed": True,
            "proved": True}


def _branch_ledger(p: int, u: int, low_k: int, offset: int,
                   forbidden_Q: int, forbidden_mass: int,
                   mass_shift: int, local: dict[str, object]) -> dict[str, object]:
    """Normalize before pigeonhole; no all-low assumption on elevated rows."""
    q, m, edge_count = (p - 1) // 2, (p + 1) // 2, 5 * p - 2
    quotient_sum = 2 * q - u
    low_mean = 2 * u + (p + 1) * low_k
    # Common row identity gives P_L=P_low+k_L-low_k for EVERY hard row.
    # Sum first, then use coefficient congruence to identify P_low.
    upper_P = (edge_count - quotient_sum + m * low_k) // m
    candidates = [P for P in range(upper_P + 1) if (P - offset) % q == 0]
    _require(candidates == [offset] and upper_P <= 9 < q,
             "coefficient congruence no longer fixes the common low parallel count")
    hard_edges = m * (offset - low_k) + quotient_sum
    hT = (p + 1) * offset - 3 * p - low_mean
    opposite_edges = edge_count - hard_edges
    next_Q, next_mass = forbidden_Q + 1, p + mass_shift
    opposite_mean = lambda Q: (p + 1) * Q + hT - 3 * p
    floors = residual_even_floor_table(p)["phase_zero_floors"]
    sharp_floor = sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"]
    nonzero_rows = [(b, f, next_mass - f) for b, f in floors.items()
                    if b != 0 and f <= next_mass]
    surplus = opposite_edges - m * next_Q
    checks = {
        "common_row_sum_equals_signed_total": 2 * hard_edges == edge_count + hT,
        "lower_Q_means_negative": opposite_mean(forbidden_Q - 1) < 0,
        "forbidden_Q_mean": opposite_mean(forbidden_Q) == forbidden_mass,
        "forbidden_positive_below_both_floors": 0 < forbidden_mass < min(
            sharp_floor, min(f for b, f in floors.items() if b)),
        "next_Q_mean": opposite_mean(next_Q) == next_mass,
        "positive_number_of_next_rows": 0 <= surplus < m,
        "only_nonzero_boundary_candidates_are_pointwise_baselines":
            [b for b, _, _ in nonzero_rows] == [2, p - 1],
        "every_nonzero_boundary_leaves_genuine_subsharp_lift":
            all(0 < excess < sharp_floor for _, _, excess in nonzero_rows),
        "boundary_zero_local_exclusion": bool(local["proved"]),
    }
    _require(all(checks.values()), "common opposite-row contradiction failed")
    return {"u": u, "low_quotient": low_k, "low_scaled_mean": low_mean,
            "coefficient_offset": offset, "common_low_parallel_candidates": candidates,
            "common_low_parallel_upper_bound": upper_P,
            "every_hard_parallel_formula": f"P_L={offset-low_k}+k_L",
            "hard_edge_count": hard_edges, "hard_sign_times_global_T": hT,
            "opposite_edge_count": opposite_edges,
            "forbidden_Q": forbidden_Q, "forbidden_scaled_mean": forbidden_mass,
            "forced_next_Q": next_Q, "forced_next_scaled_mean": next_mass,
            "surplus": surplus, "forced_next_row_count_at_least": m - surplus,
            "nonzero_boundary_floor_excess_rows": [list(row) for row in nonzero_rows],
            "checks": checks, "proved": True}


def p1_third_layer_exclusion(p: int) -> dict[str, object]:
    _check_prime(p)
    q, m = (p - 1) // 2, (p + 1) // 2
    residue_ledger = p1_third_residue_ledger(p)
    hard_catalog = hard_family_catalog(p)
    p7 = {"proved": all(row["proved"] for row in (
        height_at_least_two_certificate(p), height_one_junta_certificate(p),
        density_profile_certificate(p)))}
    dependencies = {
        7: p7, 9: p_plus_nine_local_exclusion(p),
        11: p1_p_plus_eleven_local_exclusion(p),
        13: p1_p_plus_thirteen_local_exclusion(p),
        15: p_plus_fifteen_local_exclusion(p),
    }
    specifications = [
        ("old_complement_literal", 0, 1, 5, 2, 6, 7),
        ("carried_complement_triple", q - 3, 1, 2, 6, 14, 15),
        ("carried_XNOR_omitted_pair", q - 2, 1, 3, 5, 12, 13),
        ("carried_XNOR_all_equal_triple", q - 2, 1, 5, 3, 12, 13),
        ("new_offset_four", q - 1, 1, 4, 4, 10, 11),
        ("new_offset_six", q - 1, 1, 6, 2, 10, 11),
        ("quotient_zero_XNOR", q, 0, 4, 3, 8, 9),
    ]
    branches = {name: _branch_ledger(p, u, low, offset, Q, small, shift, dependencies[shift])
                for name, u, low, offset, Q, small, shift in specifications}
    _require([row["forced_next_row_count_at_least"] for row in branches.values()]
             == [5, 9, 8, 8, 7, 7, 6], "uniform pigeonhole counts changed")
    sharp_XNOR = p1_sharp_family_catalog(p)
    _require(sharp_XNOR["proved"] and
             [f["coefficient_offset"] for f in sharp_XNOR["families"]] == [3, 5],
             "carried sharp XNOR classification failed")
    _require(p * p - 10 * p + 5 > 0 and all(row["proved"] for row in branches.values()),
             "third layer proof failed")
    return {"p": p, "p_mod_4": 1, "m": m, "layer_index_t": q - 1,
            "original_k": 5 * p - 3, "H_edge_count": 5 * p - 2,
            "residue_ledger": residue_ledger, "new_hard_family_catalog": hard_catalog,
            "branch_exclusions": branches,
            "all_boundary_sizes_excluded": True,
            "residual_ii_layer_excluded": True,
            "finite_prime_graph_or_slice_census_used": False, "proved": True}


def proposition_15772() -> dict[str, object]:
    records = [p1_third_layer_exclusion(p) for p in (29, 37, 41, 53)]
    return {"prop": "15.772", "status": "PROVED_INFINITE_FAMILY",
            "statement": "For every prime p=1 mod4, p>=29, residual(ii) at k=5p-3 is empty",
            "proof_note": "evidence/NOTE_2026-09-04_P1_THIRD_POST_BAND_CLOSE.md",
            "new_local_theorem_note": "evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md",
            "symbolic_range": {"minimum_prime": 29, "p_mod_4": 1,
                               "t": "(p-3)/2", "k": "5p-3"},
            "new_p1_frontier": "p=1 mod4,p>=29,t>=(p-1)/2",
            "records_are_identity_replays_not_exhaustive_prime_evidence": True,
            "records": records, "proved": all(row["proved"] for row in records),
            "residual_ii_closed_general": False, "e1_closed_general": False,
            "original_MO_limit_closed": False}


if __name__ == "__main__":
    result = proposition_15772()
    write_json_atomic(ROOT / "evidence/e1_gmin_m4_prop15772.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
