#!/usr/bin/env python3
"""Prop. 15.711 -- exclude the five p=17 residue-zero profiles.

If any phase-zero weight-zero direction attains its floor, Proposition
15.710 already gives a contradiction.  Avoiding that anchor in one of the
five residue-zero profiles spends every free quotient increment, forcing
directional mean 18 in both phases.

The global directional means then leave four possible infinity degrees.
Every finite edge is forced into phase one, making each phase-one weight-16
cross-cell coefficient a nonnegative edge count.  The resulting fibre
capacity bound contradicts all four infinity degrees.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15707 import _rigid_b2_lower_bound
from e1_gmin_m4_prop15710 import p17_phase_one_b16_global_sign_reduction


ROOT = Path(__file__).resolve().parents[1]
P = 17
EDGE_COUNT = 4 * P + 1


def p17_residue_zero_uniform_mean_exclusion() -> dict[str, object]:
    """Exclude all five `(u_0,u_1)=(0,0)` rows without a solver."""
    previous = p17_phase_one_b16_global_sign_reduction()
    profiles = [
        row
        for row in previous["surviving_profiles"]
        if (int(row["u0"]), int(row["u1"])) == (0, 0)
    ]
    expected_profiles = [
        (96, {"0": {0: 9}, "1": {16: 9}}),
        (100, {"0": {0: 7, 2: 2}, "1": {16: 9}}),
        (104, {"0": {0: 5, 2: 4}, "1": {16: 9}}),
        (108, {"0": {0: 3, 2: 6}, "1": {16: 9}}),
        (112, {"0": {0: 1, 2: 8}, "1": {16: 9}}),
    ]
    if [
        (int(row["pair_slack"]), row["phase_profiles_b"]) for row in profiles
    ] != expected_profiles:
        raise ArithmeticError("p=17 residue-zero five-profile block changed")

    allocation_rows = []
    for row in profiles:
        phase_zero = _rigid_b2_lower_bound(
            {
                "u0": row["u0"],
                "u1": row["u1"],
                "phase_profiles_b": row["phase_profiles_b"],
            },
            0,
        )
        b0_count = int(row["phase_profiles_b"]["0"].get(0, 0))
        b2_count = int(row["phase_profiles_b"]["0"].get(2, 0))
        free_increments = int(phase_zero["free_quotient_increments"])
        minimum_quotient_sum = int(phase_zero["minimum_quotient_sum"])
        if (
            b0_count + b2_count != 9
            or minimum_quotient_sum != b2_count
            or free_increments != b0_count
        ):
            raise ArithmeticError("p=17 residue-zero allocation saturation changed")
        allocation_rows.append(
            {
                "pair_slack": int(row["pair_slack"]),
                "phase_zero_b0_count": b0_count,
                "phase_zero_b2_count": b2_count,
                "minimum_quotient_sum": minimum_quotient_sum,
                "free_quotient_increments": free_increments,
                "avoiding_rigid_b0_forces_every_phase_zero_quotient": 1,
                "all_directional_means": 18,
            }
        )

    # With common parallel counts P_+,P_- and normalized finite sign sum Sbar,
    # mean 18 gives 18P_±=69-I±Sbar. Hence 69-I is divisible by nine.
    # Infinity is not in the 16-point odd boundary, so I is even.
    infinity_candidates = [
        infinity_degree
        for infinity_degree in range(EDGE_COUNT + 1)
        if infinity_degree % 2 == 0 and (EDGE_COUNT - infinity_degree) % 9 == 0
    ]
    if infinity_candidates != [6, 24, 42, 60]:
        raise ArithmeticError("uniform-mean infinity candidates changed")

    candidate_rows = []
    for infinity_degree in infinity_candidates:
        k = (infinity_degree - 6) // 18
        parallel_sum = (EDGE_COUNT - infinity_degree) // 9
        # A rigid phase-one b=16 direction has P_+=5+8g-I. Its residue
        # modulo eight equals the upper bound P_++P_-, so P_-=0.
        phase_one_parallel = parallel_sum
        phase_zero_parallel = 0
        gauge = (phase_one_parallel - 5 + infinity_degree) // 8
        if (
            parallel_sum != 7 - 2 * k
            or phase_one_parallel % 8 != (5 - infinity_degree) % 8
            or gauge != 1 + 2 * k
        ):
            raise ArithmeticError("uniform-mean parallel arithmetic changed")

        # For a special fibre j and ordinary-fibre maximum m, positivity of
        # L_st=g-z_s-z_t+1_{j in {s,t}} gives
        # I <= g+1+15*floor(g/2).
        cell_upper_bound = gauge + 1 + 15 * (gauge // 2)
        if not infinity_degree > cell_upper_bound:
            raise ArithmeticError("phase-one b16 fibre bound lost contradiction")
        candidate_rows.append(
            {
                "k": k,
                "infinity_degree": infinity_degree,
                "parallel_count_sum": parallel_sum,
                "phase_one_parallel_count": phase_one_parallel,
                "phase_zero_parallel_count": phase_zero_parallel,
                "phase_one_b16_gauge": gauge,
                "cell_upper_bound_on_infinity_degree": cell_upper_bound,
                "contradiction_margin": infinity_degree - cell_upper_bound,
            }
        )

    histogram = {
        int(slack): int(count)
        for slack, count in previous["remaining_pair_slack_histogram"].items()
    }
    for row in profiles:
        slack = int(row["pair_slack"])
        histogram[slack] -= 1
        if histogram[slack] == 0:
            del histogram[slack]
    after = int(previous["profile_count_after"]) - len(profiles)
    if (
        after != 14
        or sum(histogram.values()) != after
        or histogram != {96: 2, 100: 3, 104: 3, 108: 2, 112: 2, 116: 1, 128: 1}
    ):
        raise ArithmeticError("post-15.711 p=17 ledger changed")
    survivor_residues = Counter(
        (int(row["u0"]), int(row["u1"]))
        for row in previous["surviving_profiles"]
        if (int(row["u0"]), int(row["u1"])) != (0, 0)
    )
    if survivor_residues != {(7, 0): 9, (8, 0): 5}:
        raise ArithmeticError("post-15.711 residue ledger changed")
    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    excluded_indices = {int(row["census_index"]) for row in profiles}
    remaining_indices = sorted(previous_indices - excluded_indices)
    if (
        len(previous_indices) != 19
        or len(excluded_indices) != 5
        or len(remaining_indices) != 14
    ):
        raise ArithmeticError("post-15.711 profile-index ledger changed")

    return {
        "proposition": "15.711",
        "p": P,
        "boundary_size": 16,
        "profile_count_before": int(previous["profile_count_after"]),
        "profiles_excluded_here": len(profiles),
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(excluded_indices),
        "remaining_profile_indices": remaining_indices,
        "remaining_pair_slack_histogram": histogram,
        "remaining_residue_pair_histogram": {"u0=7,u1=0": 9, "u0=8,u1=0": 5},
        "allocation_saturation_rows": allocation_rows,
        "uniform_directional_mean": 18,
        "directional_mean_identities": {
            "phase_one": "18*P_+=69-I+Sbar",
            "phase_zero": "18*P_-=69-I-Sbar",
            "sum": "69-I=9*(P_++P_-)",
        },
        "infinity_degree_even_because_infinity_not_in_odd_boundary": True,
        "candidate_rows": candidate_rows,
        "phase_one_b16_cell_identity": (
            "L_st=g-z_s-z_t+1_{j in {s,t}}"
        ),
        "all_finite_edges_forced_to_phase_one": True,
        "fibre_capacity_bound": "I<=g+1+15*floor(g/2)",
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "uses_new_arc_classification": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_residue_zero_uniform_mean_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15711.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.711: p=17 profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )


if __name__ == "__main__":
    main()
