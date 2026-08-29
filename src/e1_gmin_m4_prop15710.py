#!/usr/bin/env python3
"""Prop. 15.710 -- reduce the p=17 endpoint from 227 profiles to nineteen.

Every profile left by Proposition 15.709 has nine rigid phase-one weight-16
directions.  Of the 227 rows, 176 retain a rigid phase-zero weight-zero
direction.  Comparing those two rigid floors forces infinity degree 60 and
gauge sum 14, while parallel nonnegativity requires gauge sum at least 15.

Thirty-two more rows retain a rigid phase-zero weight-16 direction. Comparing
weight-16 floors in the two phases forces infinity degree 68 and gauge sum
16, while parallel nonnegativity requires gauge sum at least 17. Exactly
nineteen profiles survive.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15707 import _minimum_quotient, _rigid_b2_lower_bound
from e1_gmin_m4_prop15709 import p17_phase_one_residue_eight_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 17


def _rigid_floor_lower_bound(profile: dict[str, object], phase: int, b: int) -> int:
    """Guaranteed directions attaining the actual symbolic floor."""
    quotient = _rigid_b2_lower_bound(profile, phase)
    residue = int(profile[f"u{phase}"])
    minimum_quotient = _minimum_quotient(phase, residue, b)
    if (
        2 * residue + 18 * minimum_quotient
        != full_symbolic_floor(P, b, phase)
    ):
        return 0
    count = int(profile["phase_profiles_b"][str(phase)].get(b, 0))
    return max(0, count - int(quotient["free_quotient_increments"]))


def p17_phase_one_b16_global_sign_reduction() -> dict[str, object]:
    """Apply complementary b0/b16 global-sign comparisons to all 227 rows."""
    previous = p17_phase_one_residue_eight_exclusion()
    previous_histogram = {
        int(slack): int(count)
        for slack, count in previous["remaining_pair_slack_histogram"].items()
    }
    census = p17_second_boundary_profile_census()
    profiles = [
        (census_index, row)
        for census_index, row in enumerate(census["profiles"])
        if int(row["pair_slack"]) in previous_histogram and int(row["u1"]) == 0
    ]
    if len(profiles) != 227 or Counter(
        int(row["pair_slack"]) for _index, row in profiles
    ) != Counter(previous_histogram):
        raise ArithmeticError("pre-15.710 p=17 profile ledger changed")

    b0_rows = []
    b16_rows = []
    survivors = []
    b0_anchor_histogram: Counter[int] = Counter()
    b16_anchor_histogram: Counter[int] = Counter()
    for census_index, profile in profiles:
        rigid_b0 = _rigid_floor_lower_bound(profile, 0, 0)
        rigid_phase_zero_b16 = _rigid_floor_lower_bound(profile, 0, 16)
        rigid_phase_one_b16 = _rigid_floor_lower_bound(profile, 1, 16)
        record = {
            "census_index": census_index,
            "pair_slack": int(profile["pair_slack"]),
            "u0": int(profile["u0"]),
            "u1": int(profile["u1"]),
            "phase_profiles_b": profile["phase_profiles_b"],
            "phase_zero_rigid_b0_lower_bound": rigid_b0,
            "phase_zero_rigid_b16_lower_bound": rigid_phase_zero_b16,
            "phase_one_rigid_b16_lower_bound": rigid_phase_one_b16,
        }
        if rigid_phase_one_b16 != 9:
            raise ArithmeticError("phase-one rigid b=16 core changed")
        if rigid_b0:
            record["exclusion"] = "phase-zero-b0-versus-phase-one-b16"
            b0_rows.append(record)
            b0_anchor_histogram[rigid_b0] += 1
        elif rigid_phase_zero_b16:
            record["exclusion"] = "opposite-phase-b16-comparison"
            b16_rows.append(record)
            b16_anchor_histogram[rigid_phase_zero_b16] += 1
        else:
            survivors.append(record)

    survivor_slack_histogram = Counter(
        int(row["pair_slack"]) for row in survivors
    )
    survivor_residue_histogram = Counter(
        (int(row["u0"]), int(row["u1"])) for row in survivors
    )
    if (
        len(b0_rows) != 176
        or b0_anchor_histogram
        != {1: 8, 2: 26, 3: 32, 4: 56, 5: 38, 6: 16}
        or Counter((row["u0"], row["u1"]) for row in b0_rows)
        != {(0, 0): 176}
        or len(b16_rows) != 32
        or b16_anchor_histogram
        != {1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 4}
        or Counter((row["u0"], row["u1"]) for row in b16_rows)
        != {(8, 0): 32}
        or len(survivors) != 19
        or survivor_slack_histogram
        != {96: 3, 100: 4, 104: 4, 108: 3, 112: 3, 116: 1, 128: 1}
        or survivor_residue_histogram
        != {(0, 0): 5, (7, 0): 9, (8, 0): 5}
    ):
        raise ArithmeticError("p=17 complementary rigid-anchor sweep changed")

    modulus = 72
    inverse = pow(P, -1, modulus)
    b0_candidate = 12 * inverse % modulus
    b16_candidate = 4 * inverse % modulus
    b0_gauge_sum = (P * b0_candidate - 12) // modulus
    b16_gauge_sum = (P * b16_candidate - 4) // modulus
    b0_minimum_gauge = (b0_candidate - 3 + 7) // 8
    phase_one_b16_at_60_minimum_gauge = (b0_candidate - 5 + 7) // 8
    phase_zero_b16_at_68_minimum_gauge = (b16_candidate - 3 + 7) // 8
    phase_one_b16_at_68_minimum_gauge = (b16_candidate - 5 + 7) // 8
    if (
        inverse != 17
        or b0_candidate != 60
        or b0_gauge_sum != 14
        or (b0_minimum_gauge, phase_one_b16_at_60_minimum_gauge) != (8, 7)
        or b16_candidate != 68
        or b16_gauge_sum != 16
        or (
            phase_zero_b16_at_68_minimum_gauge,
            phase_one_b16_at_68_minimum_gauge,
        )
        != (9, 8)
    ):
        raise ArithmeticError("p=17 complementary global-sign arithmetic changed")

    survivor_histogram = survivor_slack_histogram
    excluded_histogram = Counter(previous_histogram)
    excluded_histogram.subtract(survivor_histogram)
    if any(count < 0 for count in excluded_histogram.values()):
        raise ArithmeticError("p=17 excluded-slack accounting changed")
    excluded_histogram = +excluded_histogram

    return {
        "proposition": "15.710",
        "p": P,
        "boundary_size": 16,
        "profile_count_before": len(profiles),
        "profiles_excluded_here": len(b0_rows) + len(b16_rows),
        "profiles_excluded_by_b0_b16_comparison": len(b0_rows),
        "profiles_excluded_by_b16_b16_comparison": len(b16_rows),
        "profile_count_after": len(survivors),
        "remaining_pair_slack_histogram": dict(sorted(survivor_histogram.items())),
        "excluded_pair_slack_histogram": dict(sorted(excluded_histogram.items())),
        "remaining_residue_pair_histogram": {
            "u0=0,u1=0": 5,
            "u0=7,u1=0": 9,
            "u0=8,u1=0": 5,
        },
        "rigid_phase_one_b16_lower_bound": 9,
        "rigid_phase_zero_b0_lower_bound_histogram": dict(
            sorted(b0_anchor_histogram.items())
        ),
        "rigid_phase_zero_b16_lower_bound_histogram": dict(
            sorted(b16_anchor_histogram.items())
        ),
        "b0_b16_global_identity": "17*I=12+72*(g_0+g_16)",
        "b0_b16_contradiction": {
            "infinity_degree": b0_candidate,
            "forced_gauge_sum": b0_gauge_sum,
            "minimum_gauges": [
                b0_minimum_gauge,
                phase_one_b16_at_60_minimum_gauge,
            ],
            "minimum_gauge_sum": (
                b0_minimum_gauge + phase_one_b16_at_60_minimum_gauge
            ),
        },
        "b16_b16_global_identity": "17*I=4+72*(g_16^-+g_16^+)",
        "b16_b16_contradiction": {
            "infinity_degree": b16_candidate,
            "forced_gauge_sum": b16_gauge_sum,
            "minimum_gauges": [
                phase_zero_b16_at_68_minimum_gauge,
                phase_one_b16_at_68_minimum_gauge,
            ],
            "minimum_gauge_sum": (
                phase_zero_b16_at_68_minimum_gauge
                + phase_one_b16_at_68_minimum_gauge
            ),
        },
        "excluded_b0_rows": b0_rows,
        "excluded_b16_rows": b16_rows,
        "surviving_profiles": survivors,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "uses_new_arc_classification": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_phase_one_b16_global_sign_reduction()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15710.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.710: p=17 profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )


if __name__ == "__main__":
    main()
