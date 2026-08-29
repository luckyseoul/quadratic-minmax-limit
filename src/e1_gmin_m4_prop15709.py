#!/usr/bin/env python3
"""Prop. 15.709 -- exclude every remaining p=17 profile with u_1=8.

After Proposition 15.708, 507 exact profiles remain.  Exactly 280 have
phase-one residue eight, and every one retains at least eight rigid
phase-one weight-two directions.  The 66 rows with phase-zero residue zero
retain a rigid weight-zero direction and fail the global gauge comparison.
The 214 rows with phase-zero residue eight retain a rigid weight-16
direction and fail Proposition 15.708's unique-even-fibre cell identity.

The exact remainder is therefore 227 profiles, all with phase-one residue
zero and pair slack at least 96.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15707 import _rigid_b2_lower_bound
from e1_gmin_m4_prop15708 import p17_slack_twenty_four_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 17


def p17_phase_one_residue_eight_exclusion() -> dict[str, object]:
    """Apply the two rigid-anchor contradictions to the full live ledger."""
    previous = p17_slack_twenty_four_exclusion()
    previous_histogram = {
        int(slack): int(count)
        for slack, count in previous["remaining_pair_slack_histogram"].items()
    }
    census = p17_second_boundary_profile_census()
    profiles = [
        (census_index, row)
        for census_index, row in enumerate(census["profiles"])
        if int(row["pair_slack"]) in previous_histogram
    ]
    if len(profiles) != 507 or Counter(
        int(row["pair_slack"]) for _index, row in profiles
    ) != Counter(previous_histogram):
        raise ArithmeticError("pre-15.709 p=17 profile ledger changed")

    b0_rows = []
    b16_rows = []
    survivors = []
    b0_anchor_histogram: Counter[int] = Counter()
    b16_anchor_histogram: Counter[int] = Counter()
    excluded_slack_histogram: Counter[int] = Counter()
    survivor_slack_histogram: Counter[int] = Counter()
    for census_index, profile in profiles:
        phase_zero = _rigid_b2_lower_bound(profile, 0)
        phase_one = _rigid_b2_lower_bound(profile, 1)
        phase_zero_profile = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"]["0"].items()
        }
        free_increments = int(phase_zero["free_quotient_increments"])
        rigid_b0 = max(0, phase_zero_profile.get(0, 0) - free_increments)
        rigid_b16 = max(0, phase_zero_profile.get(16, 0) - free_increments)
        rigid_phase_one_b2 = int(phase_one["rigid_b2_lower_bound"])
        record = {
            "census_index": census_index,
            "pair_slack": int(profile["pair_slack"]),
            "u0": int(profile["u0"]),
            "u1": int(profile["u1"]),
            "phase_zero_rigid_b0_lower_bound": rigid_b0,
            "phase_zero_rigid_b16_lower_bound": rigid_b16,
            "phase_one_rigid_b2_lower_bound": rigid_phase_one_b2,
        }

        if int(profile["u1"]) == 8:
            if rigid_phase_one_b2 < 8:
                raise ArithmeticError("u1=8 row lost its rigid phase-one core")
            if int(profile["u0"]) == 0 and rigid_b0:
                record["exclusion"] = "global-gauge-b0-versus-b2"
                b0_rows.append(record)
                b0_anchor_histogram[rigid_b0] += 1
            elif int(profile["u0"]) == 8 and rigid_b16:
                record["exclusion"] = "unique-even-fibre-b16-versus-b2"
                b16_rows.append(record)
                b16_anchor_histogram[rigid_b16] += 1
            else:
                raise ArithmeticError("u1=8 row has no rigid phase-zero anchor")
            excluded_slack_histogram[int(profile["pair_slack"])] += 1
        else:
            survivors.append(record)
            survivor_slack_histogram[int(profile["pair_slack"])] += 1

    expected_excluded_slack = {
        28: 35,
        32: 26,
        36: 20,
        40: 15,
        44: 13,
        48: 11,
        52: 9,
        56: 9,
        60: 9,
        64: 9,
        68: 9,
        72: 9,
        76: 9,
        80: 9,
        84: 9,
        88: 9,
        92: 9,
        96: 9,
        100: 9,
        104: 9,
        108: 9,
        112: 9,
        116: 7,
        120: 5,
        124: 3,
        128: 1,
    }
    expected_survivor_slack = {
        96: 3,
        100: 5,
        104: 8,
        108: 13,
        112: 24,
        116: 28,
        120: 30,
        124: 26,
        128: 22,
        132: 16,
        136: 11,
        140: 7,
        144: 5,
        148: 3,
        152: 2,
        156: 2,
        160: 2,
        164: 1,
        168: 1,
        172: 1,
        176: 1,
        180: 1,
        184: 1,
        188: 1,
        192: 1,
        196: 1,
        200: 1,
        204: 1,
        208: 1,
        212: 1,
        216: 1,
        220: 1,
        224: 1,
        228: 1,
        232: 1,
        236: 1,
        240: 1,
    }
    if (
        len(b0_rows) != 66
        or b0_anchor_histogram != {3: 10, 4: 27, 5: 29}
        or len(b16_rows) != 214
        or b16_anchor_histogram != {2: 4, 3: 30, 4: 36, 5: 36, 6: 36, 7: 36, 8: 36}
        or dict(sorted(excluded_slack_histogram.items())) != expected_excluded_slack
        or len(survivors) != 227
        or dict(sorted(survivor_slack_histogram.items())) != expected_survivor_slack
        or Counter((row["u0"], row["u1"]) for row in survivors)
        != {(0, 0): 181, (8, 0): 37, (7, 0): 9}
    ):
        raise ArithmeticError("p=17 full rigid-anchor sweep changed")

    hard_certificate = previous["unique_even_fibre_certificate"]
    if (
        hard_certificate["contradiction"] is not True
        or int(hard_certificate["nonnegative_count_upper_bound"]) != -3
    ):
        raise ArithmeticError("15.708 unique-even-fibre lemma changed")
    after = len(survivors)
    if sum(expected_survivor_slack.values()) != after or min(expected_survivor_slack) != 96:
        raise ArithmeticError("post-15.709 p=17 accounting changed")

    return {
        "proposition": "15.709",
        "p": P,
        "boundary_size": 16,
        "profile_count_before": len(profiles),
        "profiles_excluded_here": len(b0_rows) + len(b16_rows),
        "profiles_excluded_by_global_gauge_identity": len(b0_rows),
        "profiles_excluded_by_unique_even_fibre_identity": len(b16_rows),
        "profile_count_after": after,
        "minimum_remaining_pair_slack": min(expected_survivor_slack),
        "remaining_pair_slack_histogram": expected_survivor_slack,
        "remaining_residue_pair_histogram": {
            "u0=0,u1=0": 181,
            "u0=7,u1=0": 9,
            "u0=8,u1=0": 37,
        },
        "excluded_pair_slack_histogram": expected_excluded_slack,
        "rigid_phase_zero_b0_lower_bound_histogram": dict(
            sorted(b0_anchor_histogram.items())
        ),
        "rigid_phase_zero_b16_lower_bound_histogram": dict(
            sorted(b16_anchor_histogram.items())
        ),
        "rigid_phase_one_b2_lower_bound": 8,
        "global_gauge_identity": "17*I=4+72*(g_0+g_1)",
        "global_gauge_contradiction": "g_0+g_1=16 but g_0>=9,g_1>=8",
        "unique_even_fibre_identity": hard_certificate,
        "excluded_b0_rows": b0_rows,
        "excluded_b16_rows": b16_rows,
        "surviving_rows": survivors,
        "all_phase_one_residue_eight_profiles_excluded": True,
        "all_survivors_have_phase_one_residue_zero": True,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "uses_new_arc_classification": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_phase_one_residue_eight_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15709.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.709: p=17 profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )


if __name__ == "__main__":
    main()
