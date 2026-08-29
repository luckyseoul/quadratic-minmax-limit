#!/usr/bin/env python3
"""Prop. 15.701 -- p=17 low-positive-slack conic-core reduction.

Proposition 15.700 leaves 1,330 phase-labelled profiles at the second
all-finite boundary ``p=17, s=16``: two of slack zero and every profile of
positive pair slack.  The exact positive-slack counts begin with 227 rows of
slack four, 195 of slack eight, and 155 of slack twelve.

Pair slack ``4r`` permits deleting at most ``r`` boundary points to obtain an
arc.  Sticker's exhaustive classification has one PGL class of 15-arcs in
``PG(2,17)``.  Since a conic with three points removed is a 15-arc, every
15-arc is conic-contained.

At slack four the repaired arc already has size at least 15.  At slack eight,
one undetermined direction completes a worst-case 14-arc to a 15-arc.  At
slack twelve, two undetermined directions complete a worst-case 13-arc to a
15-arc.  The exact profile ledger therefore puts 227, 128, and 43 rows,
respectively, on a conic core.

If ``h`` of the original sixteen affine points lie off that conic, then
``1<=h<=r<=3``.  Each off-conic point has at least eight full conic secants,
while the retained conic set omits ``2+h`` points.  Hence it lies on at least
``6-h`` retained conic secants.  A line counted for ``a`` off-conic points
has occupancy at least ``2+a`` and charges at least ``4a`` pair slack.  Thus

    pair slack >= 4*h*(6-h) >= 20,

contradicting slack four, eight, or twelve.  This excludes 398 rows and leaves
932.  It is a strict reduction, not endpoint closure: two tangent-conic
slack-zero rows and 930 positive-slack rows remain.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15678 import p17_arc_classification_ledger
from e1_gmin_m4_prop15700 import (
    p17_second_boundary_profile_census,
    p17_slack_zero_conic_reduction,
)


ROOT = Path(__file__).resolve().parents[1]
P = 17
BOUNDARY_SIZE = 16


def line_pair_slack(occupancy: int) -> int:
    """Contribution of one affine line to the global pair slack."""
    if occupancy < 0:
        raise ValueError("occupancy must be nonnegative")
    return 2 * (math.comb(occupancy, 2) - occupancy // 2)


def p17_fifteen_arc_classification() -> dict[str, object]:
    """Extract the unique 15-arc consequence from Sticker's census."""
    classification = p17_arc_classification_ledger()
    classes = classification["pgl_classes_in_pg2_17"]
    if int(classes["15"]) != 1:
        raise ArithmeticError("p=17 15-arc class count changed")
    return {
        "external_dependency": True,
        "source": classification["source"],
        "source_url": classification["source_url"],
        "location": classification["location"],
        "pgl_class_count_of_15_arcs": int(classes["15"]),
        "known_representative": "a nondegenerate conic with three points deleted",
        "consequence": "every arc of size at least 15 in PG(2,17) is conic-contained",
        "extension_reason": (
            "a 15-subarc lies on the unique conic-derived class; an additional "
            "off-conic point sees at least 8-3=5 retained conic secants and "
            "would create a 3-secant"
        ),
        "proved_conditional_on_external_classification": True,
    }


def p17_conic_core_slack_lemma() -> dict[str, object]:
    """Record the repair charge and the off-conic secant inequality."""
    line_rows = [
        {
            "occupancy": n,
            "pair_slack": line_pair_slack(n),
            "deletions_to_occupancy_two": n - 2,
        }
        for n in range(3, BOUNDARY_SIZE + 1)
    ]
    if not all(
        int(row["deletions_to_occupancy_two"])
        <= int(row["pair_slack"]) // 4
        for row in line_rows
    ):
        raise ArithmeticError("p=17 bad-line deletion charge changed")

    off_conic = {
        h: {
            "off_conic_points": h,
            "retained_conic_points": BOUNDARY_SIZE - h,
            "omitted_conic_points": (P + 1) - (BOUNDARY_SIZE - h),
            "minimum_retained_secants_per_off_conic_point": 6 - h,
            "minimum_secant_incidence_count": h * (6 - h),
            "pair_slack_floor": 4 * h * (6 - h),
        }
        for h in range(1, 4)
    }
    if [off_conic[h]["pair_slack_floor"] for h in range(1, 4)] != [20, 32, 36]:
        raise ArithmeticError("p=17 off-conic slack floor changed")
    return {
        "line_charge_rows": line_rows,
        "arc_repair_bound": (
            "pair slack 4r permits deleting at most r points to obtain an arc"
        ),
        "full_conic_secants_through_off_conic_point_at_least": 8,
        "off_conic_count_rows": off_conic,
        "counting_argument": (
            "if h<=3 of the 16 original points are off a conic, each sees at "
            "least 8-(18-(16-h))=6-h secants whose two conic points remain; "
            "a line counted for a off-conic points has occupancy at least "
            "2+a and contributes at least 4a pair slack"
        ),
        "positive_slack_below_twenty_impossible_after_conic_core": True,
        "proved": True,
    }


def _undetermined_directions(row: dict[str, object]) -> int:
    profiles = row["phase_profiles_b"]
    return sum(int(profiles[phase].get(BOUNDARY_SIZE, 0)) for phase in ("0", "1"))


def p17_low_positive_slack_profile_ledger() -> dict[str, object]:
    """Audit the exact low-slack rows that reach a classified 15-arc."""
    census = p17_second_boundary_profile_census()
    rules = {
        4: {"delete_at_most": 1, "required_undetermined": 0, "adjoin": 0},
        8: {"delete_at_most": 2, "required_undetermined": 1, "adjoin": 1},
        12: {"delete_at_most": 3, "required_undetermined": 2, "adjoin": 2},
    }
    rows = []
    for slack, rule in rules.items():
        profiles = [
            row for row in census["profiles"] if int(row["pair_slack"]) == slack
        ]
        t_histogram = dict(sorted(Counter(_undetermined_directions(row) for row in profiles).items()))
        qualifying = [
            row
            for row in profiles
            if _undetermined_directions(row) >= int(rule["required_undetermined"])
        ]
        classified_size_floor = (
            BOUNDARY_SIZE - int(rule["delete_at_most"]) + int(rule["adjoin"])
        )
        if classified_size_floor < 15:
            raise ArithmeticError("p=17 repair failed to reach a 15-arc")
        rows.append(
            {
                "pair_slack": slack,
                "profile_count": len(profiles),
                "undetermined_direction_histogram": t_histogram,
                "delete_at_most": int(rule["delete_at_most"]),
                "required_undetermined_directions": int(
                    rule["required_undetermined"]
                ),
                "adjoined_infinity_points_at_most": int(rule["adjoin"]),
                "classified_arc_size_floor": classified_size_floor,
                "excluded_profile_count": len(qualifying),
                "remaining_profile_count": len(profiles) - len(qualifying),
            }
        )
    observed = {
        int(row["pair_slack"]): (
            int(row["profile_count"]),
            dict(row["undetermined_direction_histogram"]),
            int(row["excluded_profile_count"]),
            int(row["remaining_profile_count"]),
        )
        for row in rows
    }
    expected = {
        4: (227, {0: 113, 1: 102, 2: 12}, 227, 0),
        8: (195, {0: 67, 1: 104, 2: 24}, 128, 67),
        12: (155, {0: 33, 1: 79, 2: 43}, 43, 112),
    }
    if observed != expected:
        raise ArithmeticError("p=17 low-slack profile ledger changed")
    return {
        "rows": rows,
        "excluded_profile_count": sum(int(row["excluded_profile_count"]) for row in rows),
        "remaining_profile_count_in_handled_slacks": sum(
            int(row["remaining_profile_count"]) for row in rows
        ),
        "proved": True,
    }


def p17_low_positive_slack_conic_reduction() -> dict[str, object]:
    """Proposition 15.701."""
    previous = p17_slack_zero_conic_reduction()
    census = p17_second_boundary_profile_census()
    classification = p17_fifteen_arc_classification()
    lemma = p17_conic_core_slack_lemma()
    ledger = p17_low_positive_slack_profile_ledger()

    before = int(previous["profile_count_after"])
    excluded = int(ledger["excluded_profile_count"])
    after = before - excluded
    histogram = dict(census["pair_slack_histogram"])
    histogram[0] = int(previous["slack_zero_profile_count_after"])
    histogram.pop(4)
    histogram[8] = 67
    histogram[12] = 112
    histogram = dict(sorted(histogram.items()))
    if excluded != 398 or after != 932 or sum(histogram.values()) != after:
        raise ArithmeticError("p=17 post-conic reduction accounting changed")
    return {
        "proposition": "15.701",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 2,
        "remaining_positive_slack_profiles": after - 2,
        "profile_ledger": ledger,
        "repair_and_off_conic_lemma": lemma,
        "external_classification": classification,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p17_low_positive_slack_conic_reduction()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15701.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.701: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
