#!/usr/bin/env python3
"""Prop. 15.707 -- exclude all 193 p=17 pair-slack-twenty profiles.

Proposition 15.706's contradiction needs only one rigid b=2 direction of
each quadratic type, not pair slack zero. Every slack-twenty row has at least
eight rigid phase-one b=2 directions. In all 184 rows with u_0=0, quotient
accounting retains at least three rigid phase-zero directions with b=0 or 2.
Both b values give the same global-sign identity, so all 184 are impossible.
The nine u_0=8 rows all have at least two undetermined directions; repair
plus the already-audited complete 13-/14-arc data excludes those as well.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15701 import p17_fifteen_arc_classification
from e1_gmin_m4_prop15702 import complete_fourteen_arc_secant_index_certificate
from e1_gmin_m4_prop15703 import (
    complete_fourteen_minus_one_certificate,
    complete_thirteen_arc_certificate,
)
from e1_gmin_m4_prop15706 import (
    p17_slack_zero_global_sign_certificate,
    p17_slack_zero_profile_exclusion,
)
from e1_gmin_m4_prop15723 import floor_excess_admissible


ROOT = Path(__file__).resolve().parents[1]
P = 17
M = 9
PERIOD = 18
PAIR_SLACK = 20


def _minimum_quotient(phase: int, residue: int, b: int) -> int:
    """Least admissible quotient in the exact 15.700 ledger."""
    floor = full_symbolic_floor(P, b, phase)
    for quotient in range(M - residue + 1):
        excess = 2 * residue + PERIOD * quotient - floor
        if floor_excess_admissible(P, b, phase, excess):
            return quotient
    raise ArithmeticError("profile direction has no admissible quotient")


def _rigid_b2_lower_bound(row: dict[str, object], phase: int) -> dict[str, int]:
    residue = int(row[f"u{phase}"])
    profile = {
        int(b): int(count)
        for b, count in row["phase_profiles_b"][str(phase)].items()
    }
    minimum_sum = sum(
        count * _minimum_quotient(phase, residue, b)
        for b, count in profile.items()
    )
    quotient_sum = M - residue
    free_increments = quotient_sum - minimum_sum
    if free_increments < 0:
        raise ArithmeticError("minimum quotients exceed their exact sum")
    b2_count = profile.get(2, 0)
    b2_minimum = _minimum_quotient(phase, residue, 2)
    b2_floor = full_symbolic_floor(P, 2, phase)
    b2_at_minimum_is_rigid = 2 * residue + PERIOD * b2_minimum == b2_floor
    rigid_lower_bound = (
        max(0, b2_count - free_increments)
        if b2_at_minimum_is_rigid
        else 0
    )
    return {
        "residue": residue,
        "quotient_sum": quotient_sum,
        "minimum_quotient_sum": minimum_sum,
        "free_quotient_increments": free_increments,
        "b2_count": b2_count,
        "b2_minimum_quotient": b2_minimum,
        "b2_at_minimum_is_rigid": b2_at_minimum_is_rigid,
        "rigid_b2_lower_bound": rigid_lower_bound,
    }


def p17_slack_twenty_two_direction_geometric_certificate() -> dict[str, object]:
    """Profile-blind exclusion for slack twenty with two free directions."""
    fifteen = p17_fifteen_arc_classification()
    fourteen = complete_fourteen_arc_secant_index_certificate()
    thirteen = complete_thirteen_arc_certificate()
    fourteen_minus_one = complete_fourteen_minus_one_certificate()
    if int(fifteen["pgl_class_count_of_15_arcs"]) != 1:
        raise ArithmeticError("p17 fifteen-arc classification changed")
    if int(fourteen["minimum_outside_secant_index"]) != 2:
        raise ArithmeticError("p17 complete-fourteen secant floor changed")
    if max(int(x) for x in thirteen["index_one_point_counts"]) != 3:
        raise ArithmeticError("p17 complete-thirteen c1 maximum changed")
    if max(int(x) for x in fourteen_minus_one["index_one_count_histogram"]) != 4:
        raise ArithmeticError("p17 complete-fourteen-minus-one c1 maximum changed")
    conic_floors = {h: 4 * h * (7 - h) for h in range(1, 4)}
    if conic_floors != {1: 24, 2: 40, 3: 48}:
        raise ArithmeticError("two-direction conic slack floor changed")
    return {
        "repair_depth_range": [1, 5],
        "two_undetermined_points_adjoined": True,
        "repair_depth_at_most_three": {
            "classified_arc_size_at_least": 15,
            "conic_off_point_slack_floors": conic_floors,
            "minimum_positive_slack": 24,
        },
        "repair_depth_four": {
            "complete_fourteen_minimum_outside_secant_index": 2,
            "four_deleted_point_slack_floor": 32,
        },
        "repair_depth_five": {
            "five_deleted_points_force_secant_index_one": True,
            "complete_thirteen_maximum_index_one_points": 3,
            "complete_fourteen_minus_one_maximum_index_one_points": 4,
            "required_index_one_points": 5,
        },
        "proved_conditional_on_previously_audited_arc_classifications": True,
    }
def p17_slack_twenty_exclusion() -> dict[str, object]:
    previous = p17_slack_zero_profile_exclusion()
    global_sign = p17_slack_zero_global_sign_certificate()
    census = p17_second_boundary_profile_census()
    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    profiles = [
        (index, census["profiles"][index])
        for index in sorted(previous_indices)
        if int(census["profiles"][index]["pair_slack"]) == PAIR_SLACK
    ]
    if len(profiles) != 193 or Counter(
        (row["u0"], row["u1"]) for _index, row in profiles
    ) != {
        (0, 8): 184,
        (8, 8): 9,
    }:
        raise ArithmeticError("p=17 slack-twenty residue ledger changed")
    undetermined_histogram = Counter(
        sum(
            int(row["phase_profiles_b"][phase].get(16, 0))
            for phase in ("0", "1")
        )
        for _index, row in profiles
    )
    if undetermined_histogram != {0: 59, 1: 74, 2: 50, 3: 10}:
        raise ArithmeticError("p=17 slack-twenty direction histogram changed")

    rows = []
    global_sign_excluded = []
    geometric_excluded = []
    rigid_low_histogram: Counter[int] = Counter()
    for local_index, (census_index, profile) in enumerate(profiles):
        phase_zero = _rigid_b2_lower_bound(profile, 0)
        phase_one = _rigid_b2_lower_bound(profile, 1)
        if phase_one["rigid_b2_lower_bound"] < 8:
            raise ArithmeticError("slack-twenty phase one lost its rigid b=2 core")
        row = {
            "profile_index_within_slack_twenty": local_index,
            "census_index": census_index,
            "u0": int(profile["u0"]),
            "u1": int(profile["u1"]),
            "phase_profiles_b": profile["phase_profiles_b"],
            "phase_zero_quotient_certificate": phase_zero,
            "phase_one_quotient_certificate": phase_one,
        }
        phase_zero_profile = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"]["0"].items()
        }
        low_at_minimum = (
            phase_zero_profile.get(0, 0) + phase_zero_profile.get(2, 0)
            if int(profile["u0"]) == 0
            else 0
        )
        rigid_low_lower_bound = max(
            0,
            low_at_minimum - phase_zero["free_quotient_increments"],
        )
        row["phase_zero_rigid_b0_or_b2_lower_bound"] = rigid_low_lower_bound
        rows.append(row)
        if rigid_low_lower_bound >= 1:
            global_sign_excluded.append(row)
            rigid_low_histogram[rigid_low_lower_bound] += 1
        else:
            undetermined = sum(
                int(profile["phase_profiles_b"][phase].get(16, 0))
                for phase in ("0", "1")
            )
            row["undetermined_direction_count"] = undetermined
            if undetermined >= 2:
                geometric_excluded.append(row)

    if (
        len(global_sign_excluded) != 184
        or dict(sorted(rigid_low_histogram.items()))
        != {3: 11, 4: 42, 5: 88, 6: 38, 7: 5}
        or Counter(row["undetermined_direction_count"] for row in geometric_excluded)
        != {2: 5, 3: 4}
    ):
        raise ArithmeticError("p=17 rigid-low-direction profile split changed")
    if not global_sign["proved"] or global_sign["infinity_degree_candidates_in_range"] != [68]:
        raise ArithmeticError("15.706 global-sign contradiction changed")

    histogram = dict(previous["remaining_pair_slack_histogram"])
    if int(previous["profile_count_after"]) != 1213 or histogram.get(PAIR_SLACK) != 193:
        raise ArithmeticError("pre-15.707 p=17 ledger changed")
    geometry = p17_slack_twenty_two_direction_geometric_certificate()
    excluded_indices = {
        int(row["census_index"])
        for row in global_sign_excluded + geometric_excluded
    }
    excluded_count = len(excluded_indices)
    if excluded_count != len(global_sign_excluded) + len(geometric_excluded):
        raise ArithmeticError("15.707 exclusion branches overlap")
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
    if after != 1020 or PAIR_SLACK in histogram or sum(histogram.values()) != after:
        raise ArithmeticError("post-15.707 p=17 ledger changed")

    return {
        "proposition": "15.707",
        "p": P,
        "boundary_size": 16,
        "pair_slack_treated": PAIR_SLACK,
        "profile_count_before": int(previous["profile_count_after"]),
        "slack_twenty_profiles_before": len(profiles),
        "undetermined_direction_histogram_before": dict(
            sorted(undetermined_histogram.items())
        ),
        "profiles_excluded_here": excluded_count,
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(excluded_indices),
        "remaining_profile_indices": remaining_indices,
        "slack_twenty_profiles_after": 0,
        "remaining_pair_slack_histogram": histogram,
        "phase_one_rigid_b2_lower_bound": 8,
        "phase_zero_split": {
            "u0_zero_profiles_forced_to_retain_rigid_b0_or_b2": len(global_sign_excluded),
            "rigid_b0_or_b2_lower_bound_histogram": dict(
                sorted(rigid_low_histogram.items())
            ),
            "u0_eight_profiles": 9,
        },
        "global_sign_excluded_profiles": global_sign_excluded,
        "two_direction_geometric_excluded_profiles": geometric_excluded,
        "phase_zero_rigid_floor_identities": {
            "b0": {"mean": 0, "pair_target_sum": 0},
            "b2": {"mean": 18, "pair_target_sum": -1},
            "common_global_constant": 3,
        },
        "reused_global_sign_identity": global_sign["opposite_type_comparison"],
        "forced_infinity_degree": 68,
        "impossible_affine_boundary_sizes": [66, 68, 70],
        "two_direction_geometric_certificate": geometry,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_analytically": True,
        "proved_conditional_on_previously_audited_arc_classifications": True,
    }


def main() -> None:
    theorem = p17_slack_twenty_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15707.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.707: p=17 slack-twenty profiles "
        f"{theorem['slack_twenty_profiles_before']} -> "
        f"{theorem['slack_twenty_profiles_after']}; "
        f"{theorem['profile_count_after']} total profiles remain"
    )


if __name__ == "__main__":
    main()
