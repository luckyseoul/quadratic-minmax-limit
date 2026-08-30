#!/usr/bin/env python3
"""Prop. 15.708 -- exclude all 151 p=17 pair-slack-twenty-four profiles.

The 151 profiles split into 142 rows with ``(u_0,u_1)=(0,8)`` and nine
rows with ``(u_0,u_1)=(8,8)``.  Quotient accounting retains a rigid
phase-zero floor of weight zero in every row of the first block and
at least eight rigid phase-one weight-two floors in every row.  The global
Paley-sign identity of Proposition 15.706 therefore excludes all 142 rows.

For each of the remaining nine rows, comparing a rigid phase-zero weight-16
floor with a rigid phase-one weight-two floor forces the infinity degree to
be four.  The cell identities incident with the unique even fibre of that
weight-16 direction then give a negative value for a nonnegative edge count.
Thus all 151 rows are impossible without a solver or a new arc census.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15707 import (
    _rigid_b2_lower_bound,
    p17_slack_twenty_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]
P = 17
PAIR_SLACK = 24


def p17_slack_twenty_four_arithmetic_reduction() -> dict[str, object]:
    """Exclude 45 profiles and reduce the other nine to infinity degree four."""
    census = p17_second_boundary_profile_census()
    previous = p17_slack_twenty_exclusion()
    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    profiles = [
        (index, census["profiles"][index])
        for index in sorted(previous_indices)
        if int(census["profiles"][index]["pair_slack"]) == PAIR_SLACK
    ]
    residue_split = Counter(
        (int(row["u0"]), int(row["u1"])) for _index, row in profiles
    )
    if len(profiles) != 151 or residue_split != {(0, 8): 142, (8, 8): 9}:
        raise ArithmeticError("p=17 slack-twenty-four residue ledger changed")
    undetermined_histogram = Counter(
        sum(
            int(row["phase_profiles_b"][phase].get(16, 0))
            for phase in ("0", "1")
        )
        for _index, row in profiles
    )
    if undetermined_histogram != {0: 40, 1: 47, 2: 45, 3: 19}:
        raise ArithmeticError("p=17 slack-twenty-four direction histogram changed")

    easy_rows = []
    hard_rows = []
    rigid_low_histogram: Counter[int] = Counter()
    hard_undetermined_histogram: Counter[int] = Counter()
    for local_index, (census_index, profile) in enumerate(profiles):
        phase_zero = _rigid_b2_lower_bound(profile, 0)
        phase_one = _rigid_b2_lower_bound(profile, 1)
        if int(phase_one["rigid_b2_lower_bound"]) < 8:
            raise ArithmeticError("slack-24 phase one lost its rigid b=2 core")
        phase_zero_profile = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"]["0"].items()
        }
        row = {
            "profile_index_within_slack_twenty_four": local_index,
            "census_index": census_index,
            "u0": int(profile["u0"]),
            "u1": int(profile["u1"]),
            "phase_profiles_b": profile["phase_profiles_b"],
            "phase_zero_quotient_certificate": phase_zero,
            "phase_one_quotient_certificate": phase_one,
        }
        if int(profile["u0"]) == 0:
            # A b=0 direction is rigid exactly at its minimum quotient zero.
            # Spending one free quotient increment can spoil at most one such
            # direction, so this lower bound is independent of the b=2 rows.
            low_at_minimum = phase_zero_profile.get(0, 0)
            retained = max(
                0,
                low_at_minimum - int(phase_zero["free_quotient_increments"]),
            )
            if retained < 1:
                raise ArithmeticError("slack-24 easy row lost every rigid low floor")
            row["phase_zero_rigid_b0_lower_bound"] = retained
            rigid_low_histogram[retained] += 1
            easy_rows.append(row)
        else:
            undetermined = sum(
                int(profile["phase_profiles_b"][phase].get(16, 0))
                for phase in ("0", "1")
            )
            row["undetermined_direction_count"] = undetermined
            hard_undetermined_histogram[undetermined] += 1
            hard_rows.append(row)

    if (
        len(easy_rows) != 142
        or dict(sorted(rigid_low_histogram.items()))
        != {2: 3, 3: 40, 4: 59, 5: 40}
        or len(hard_rows) != 9
        or dict(sorted(hard_undetermined_histogram.items())) != {2: 3, 3: 6}
    ):
        raise ArithmeticError("p=17 slack-24 rigid-direction split changed")

    modulus = 72
    inverse = pow(P, -1, modulus)
    easy_candidates = [
        value for value in range(70) if value % modulus == 4 * inverse % modulus
    ]
    hard_candidates = [
        value for value in range(70) if value % modulus == -4 * inverse % modulus
    ]
    easy_gauge_sum = (P * easy_candidates[0] - 4) // modulus
    hard_gauge_sum = (P * hard_candidates[0] + 4) // modulus
    easy_b0_minimum_gauge = (easy_candidates[0] - 3 + 7) // 8
    easy_b2_minimum_gauge = (easy_candidates[0] - 4 + 7) // 8
    hard_b16_minimum_gauge = (hard_candidates[0] - 3 + 7) // 8
    hard_b2_minimum_gauge = (hard_candidates[0] - 4 + 7) // 8
    hard_b16_parallel_count = 3 + 8 * hard_b16_minimum_gauge - hard_candidates[0]
    hard_b2_parallel_count = 4 + 8 * hard_b2_minimum_gauge - hard_candidates[0]
    if (
        inverse != 17
        or easy_candidates != [68]
        or hard_candidates != [4]
        or easy_gauge_sum != 16
        or (easy_b0_minimum_gauge, easy_b2_minimum_gauge) != (9, 8)
        or hard_gauge_sum != 1
        or hard_b16_minimum_gauge != 1
        or hard_b2_minimum_gauge != 0
        or (hard_b16_parallel_count, hard_b2_parallel_count) != (7, 0)
    ):
        raise ArithmeticError("p=17 slack-24 global-sign arithmetic changed")

    histogram = dict(previous["remaining_pair_slack_histogram"])
    if int(previous["profile_count_after"]) != 1020 or histogram.get(PAIR_SLACK) != 151:
        raise ArithmeticError("pre-15.708 p=17 ledger changed")
    excluded_indices = {int(row["census_index"]) for row in easy_rows}
    remaining_indices = sorted(previous_indices - excluded_indices)
    histogram = dict(
        sorted(
            Counter(
                int(census["profiles"][index]["pair_slack"])
                for index in remaining_indices
            ).items()
        )
    )
    after = len(remaining_indices)
    if after != 878 or histogram.get(PAIR_SLACK) != 9 or sum(histogram.values()) != after:
        raise ArithmeticError("post-arithmetic p=17 ledger changed")

    return {
        "proposition": "15.708-arithmetic",
        "p": P,
        "boundary_size": 16,
        "pair_slack_treated": PAIR_SLACK,
        "profile_count_before": int(previous["profile_count_after"]),
        "slack_twenty_four_profiles_before": len(profiles),
        "undetermined_direction_histogram_before": dict(
            sorted(undetermined_histogram.items())
        ),
        "profiles_excluded_by_global_sign_identity": len(easy_rows),
        "slack_twenty_four_profiles_after_arithmetic": len(hard_rows),
        "profile_count_after_arithmetic": after,
        "excluded_profile_indices_by_global_sign": sorted(excluded_indices),
        "remaining_profile_indices_after_arithmetic": remaining_indices,
        "remaining_pair_slack_histogram_after_arithmetic": histogram,
        "residue_split": {"u0=0,u1=8": 142, "u0=8,u1=8": 9},
        "phase_one_rigid_b2_lower_bound": 8,
        "easy_rigid_b0_lower_bound_histogram": dict(
            sorted(rigid_low_histogram.items())
        ),
        "easy_global_sign_identity": "17*I=4+72*(g_0+g_1)",
        "easy_infinity_degree_candidate": easy_candidates[0],
        "easy_parallel_nonnegativity_contradiction": (
            "g_0>=9 and g_1>=8, but g_0+g_1=16"
        ),
        "hard_undetermined_direction_histogram": dict(
            sorted(hard_undetermined_histogram.items())
        ),
        "hard_global_sign_identity": "17*I=-4+72*(g_16+g_2)",
        "hard_infinity_degree_candidate": hard_candidates[0],
        "hard_forced_gauges_and_parallel_counts": {
            "phase_zero_b16": {
                "gauge": hard_b16_minimum_gauge,
                "parallel_count": hard_b16_parallel_count,
            },
            "phase_one_b2": {
                "gauge": hard_b2_minimum_gauge,
                "parallel_count": hard_b2_parallel_count,
            },
        },
        "easy_rows": easy_rows,
        "hard_rows": hard_rows,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_analytically": True,
    }


def p17_slack_twenty_four_exclusion() -> dict[str, object]:
    """Exclude the nine I=4 rows by one rigid b=16 cross-cell identity."""
    reduction = p17_slack_twenty_four_arithmetic_reduction()
    census = p17_second_boundary_profile_census()
    hard_rows = list(reduction["hard_rows"])
    rigid_b16_counts = []
    for row in hard_rows:
        b16_count = int(row["phase_profiles_b"]["0"].get(16, 0))
        free_increments = int(
            row["phase_zero_quotient_certificate"]["free_quotient_increments"]
        )
        rigid_b16_counts.append(max(0, b16_count - free_increments))
    if rigid_b16_counts != [3, 3, 3, 3, 2, 2, 2, 2, 2]:
        raise ArithmeticError("p=17 slack-24 rigid b=16 core changed")

    # For a rigid phase-zero b=16 direction, let j be its unique even
    # fibre.  The canonical floor is 1-x_j.  With z_s infinity-star counts,
    # g=1, and L_st the eps-weighted finite cross-cell sum, coefficient
    # comparison gives L_jt=-z_j-z_t.  There are 65 finite edges: 64 have
    # phase zero and contribute +1 to L, while the unique phase-one edge
    # contributes -1.  Hence N_j-delta_j=-15*z_j-I, where N_j>=0 and
    # delta_j is zero or one.  I=4 makes N_j<=-3.
    infinity_degree = int(reduction["hard_infinity_degree_candidate"])
    finite_edges = 4 * P + 1 - infinity_degree
    normalized_finite_sign_sum = 5 - P * infinity_degree
    phase_one_finite_edges = (finite_edges + normalized_finite_sign_sum) // 2
    phase_zero_finite_edges = finite_edges - phase_one_finite_edges
    delta_maximum = phase_one_finite_edges
    incident_count_candidates = [
        delta - 15 * z_j - infinity_degree
        for z_j in range(infinity_degree + 1)
        for delta in range(delta_maximum + 1)
    ]
    impossible_upper_bound = max(incident_count_candidates)
    if (
        infinity_degree != 4
        or finite_edges != 65
        or normalized_finite_sign_sum != -63
        or (phase_zero_finite_edges, phase_one_finite_edges) != (64, 1)
        or impossible_upper_bound != -3
        or any(value >= 0 for value in incident_count_candidates)
    ):
        raise ArithmeticError("p=17 slack-24 exceptional-edge count changed")

    arithmetic_indices = set(
        int(index) for index in reduction["remaining_profile_indices_after_arithmetic"]
    )
    hard_indices = {int(row["census_index"]) for row in hard_rows}
    if not hard_indices <= arithmetic_indices or len(hard_indices) != 9:
        raise ArithmeticError("p=17 slack-24 hard index ledger changed")
    remaining_indices = sorted(arithmetic_indices - hard_indices)
    remaining_histogram = dict(
        sorted(
            Counter(
                int(census["profiles"][index]["pair_slack"])
                for index in remaining_indices
            ).items()
        )
    )
    if PAIR_SLACK in remaining_histogram:
        raise ArithmeticError("p=17 slack-24 hard ledger changed")
    before = int(reduction["profile_count_before"])
    excluded = int(reduction["slack_twenty_four_profiles_before"])
    after = len(remaining_indices)
    if (
        before != 1020
        or excluded != 151
        or after != 869
        or sum(remaining_histogram.values()) != after
    ):
        raise ArithmeticError("post-15.708 p=17 ledger changed")

    return {
        "proposition": "15.708",
        "p": P,
        "boundary_size": 16,
        "pair_slack_treated": PAIR_SLACK,
        "profile_count_before": before,
        "slack_twenty_four_profiles_before": excluded,
        "profiles_excluded_by_global_sign_identity": int(
            reduction["profiles_excluded_by_global_sign_identity"]
        ),
        "profiles_excluded_by_unique_even_fibre_identity": len(hard_rows),
        "profiles_excluded_here": excluded,
        "slack_twenty_four_profiles_after": 0,
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(
            set(int(index) for index in reduction["excluded_profile_indices_by_global_sign"])
            | hard_indices
        ),
        "remaining_profile_indices": remaining_indices,
        "remaining_pair_slack_histogram": remaining_histogram,
        "arithmetic_reduction": reduction,
        "hard_rigid_phase_zero_b16_lower_bounds": rigid_b16_counts,
        "unique_even_fibre_certificate": {
            "rigid_floor": "A_d=1-x_j",
            "cell_identity": "L_st=g_d-z_s-z_t-1_{j in {s,t}}",
            "forced_gauge": 1,
            "incident_cell_identity": "L_jt=-z_j-z_t",
            "phase_zero_finite_edges": phase_zero_finite_edges,
            "phase_one_finite_edges": phase_one_finite_edges,
            "incident_signed_count": "N_j-delta_j=-15*z_j-I",
            "constraints": "N_j>=0, delta_j in {0,1}, z_j>=0, I=4",
            "nonnegative_count_upper_bound": impossible_upper_bound,
            "all_integral_incident_count_candidates": incident_count_candidates,
            "contradiction": True,
        },
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "uses_new_arc_classification": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_slack_twenty_four_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15708.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.708: p=17 slack-24 profiles "
        f"{theorem['slack_twenty_four_profiles_before']} -> "
        f"{theorem['slack_twenty_four_profiles_after']}"
    )


if __name__ == "__main__":
    main()
