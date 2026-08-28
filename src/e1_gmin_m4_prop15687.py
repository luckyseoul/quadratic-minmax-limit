#!/usr/bin/env python3
"""Prop. 15.687 -- exclude all 68 p=23 pair-slack-twenty profiles.

Every slack-20 profile has at least two undetermined directions. Repair
deletes at most five points.

With at least three undetermined directions, use two overlapping pairs of
their infinity points. If repair uses at most four deletions, each pair gives
an arc of size at least 18 and hence a conic; the two conics share the
repaired arc and coincide, forcing three collinear infinity points onto one
conic. If all five deletions are needed, each pair gives a 17-arc. A complete
one would force five multiplicity-one outside points, impossible by the
five-class certificate. Therefore both pair arcs extend to conics, and the
same common-conic contradiction applies.

The only hard arithmetic rows have exactly two undetermined directions.
Fewer than five deletions again reach an 18-arc. With five deletions, the
repaired 15-arc plus the two infinity points is a 17-arc K. It must be
complete, or an extension to size 18 gives the conic contradiction. The
undetermined directions ensure that every K-secant through a deleted point
is already a secant of the repaired arc. Slack equality then requires five
outside points of K with secant multiplicity one. Proposition 15.685's
exhaustive five-class certificate proves that the maximum is one.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15684 import p23_residue_zero_profile_census
from e1_gmin_m4_prop15685 import complete_17_arc_classification_certificate
from e1_gmin_m4_prop15686 import p23_slack_sixteen_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 23
BOUNDARY_SIZE = 20
PAIR_SLACK = 20
REPAIR_DELETION_BOUND = 5


def five_point_conic_core_certificate() -> dict[str, object]:
    """Extend 15.684's off-conic slack count through h=5."""
    rows = {
        h: {
            "off_conic_points": h,
            "retained_conic_points": BOUNDARY_SIZE - h,
            "omitted_conic_points": (P + 1) - (BOUNDARY_SIZE - h),
            "retained_secants_per_off_conic_point": 7 - h,
            "pair_slack_floor": 4 * h * (7 - h),
        }
        for h in range(1, 6)
    }
    minimum = min(int(row["pair_slack_floor"]) for row in rows.values())
    if minimum != 24 or int(rows[5]["pair_slack_floor"]) != 40:
        raise ArithmeticError("five-point conic-core floor changed")
    return {
        "off_conic_count_rows": rows,
        "minimum_positive_pair_slack": minimum,
        "slack_twenty_impossible_after_conic_core": PAIR_SLACK < minimum,
        "proof": (
            "an off-conic point has at least 11 full conic secants; omitting "
            "4+h conic points leaves at least 7-h retained secants, each "
            "charging four slack per off-conic incidence"
        ),
        "proved": True,
    }


def slack_twenty_profile_certificate() -> dict[str, object]:
    """Record and split all 68 exact arithmetic rows by t0."""
    census = p23_residue_zero_profile_census()
    profile_count = int(census["pair_slack_histogram"][PAIR_SLACK])
    t0_histogram = dict(
        census["undetermined_direction_histogram_by_slack"][PAIR_SLACK]
    )
    if profile_count != 68 or t0_histogram != {2: 2, 3: 36, 4: 30}:
        raise ArithmeticError("slack-twenty profile block changed")
    return {
        "pair_slack": PAIR_SLACK,
        "profile_count": profile_count,
        "repair_deletion_bound": REPAIR_DELETION_BOUND,
        "undetermined_direction_histogram": t0_histogram,
        "profiles_with_at_least_three_undetermined_directions": (
            t0_histogram[3] + t0_histogram[4]
        ),
        "hard_two_direction_profiles": t0_histogram[2],
        "proved": True,
    }


def p23_slack_twenty_exclusion() -> dict[str, object]:
    """Proposition 15.687."""
    profiles = slack_twenty_profile_certificate()
    conic = five_point_conic_core_certificate()
    classification = complete_17_arc_classification_certificate()
    required = REPAIR_DELETION_BOUND
    available = int(classification["maximum_one_secant_point_count"])
    if not available < required:
        raise ArithmeticError("complete-17-arc obstruction changed")

    previous = p23_slack_sixteen_exclusion()
    before = int(previous["p23_profile_count_after"])
    removed = int(profiles["profile_count"])
    after = before - removed
    remaining_histogram = dict(previous["remaining_pair_slack_histogram"])
    remaining_histogram.pop(PAIR_SLACK)
    if before != 201 or removed != 68 or after != 133:
        raise ArithmeticError("post-15.687 profile count changed")
    if sum(remaining_histogram.values()) != after or min(remaining_histogram) != 24:
        raise ArithmeticError("post-15.687 histogram changed")

    return {
        "proposition": "15.687",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "excluded_pair_slack": PAIR_SLACK,
        "profile_count_excluded": removed,
        "at_least_three_direction_branch": {
            "profile_count": int(
                profiles["profiles_with_at_least_three_undetermined_directions"]
            ),
            "pairwise_extension": (
                "choose U1,U2,U3 and compare A+{U1,U2} with A+{U1,U3}; "
                "three infinity points are never adjoined simultaneously"
            ),
            "at_most_four_deletions": (
                "both pair arcs have size at least 18, so their conics "
                "coincide on A and contain three collinear infinity points"
            ),
            "five_deletions": (
                "each pair arc has size 17; completeness would require five "
                "multiplicity-one outside points, so both are incomplete "
                "and extend to the same impossible conic"
            ),
            "required_one_secant_points_if_pair_arc_complete": required,
            "maximum_available_in_any_complete_17_arc_class": available,
            "excluded": True,
        },
        "two_direction_branch": {
            "profile_count": int(profiles["hard_two_direction_profiles"]),
            "five_deletion_repaired_arc_size": 15,
            "adjoined_infinity_points": 2,
            "resulting_arc_size": 17,
            "resulting_arc_must_be_complete": True,
            "required_one_secant_points": required,
            "maximum_available_in_any_complete_17_arc_class": available,
            "excluded": True,
        },
        "p23_profile_count_before": before,
        "p23_profile_count_after": after,
        "remaining_pair_slack_histogram": remaining_histogram,
        "all_remaining_profiles_have_pair_slack_at_least": min(
            remaining_histogram
        ),
        "p23_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "profile_certificate": profiles,
        "conic_core_certificate": conic,
        "classification_certificate": classification,
        "external_dependencies": previous["external_dependencies"],
        "proved": True,
    }


def main() -> None:
    theorem = p23_slack_twenty_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15687.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.687: all 68 p=23 slack-20 profiles excluded; "
        f"exact remainder {theorem['p23_profile_count_before']} -> "
        f"{theorem['p23_profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
