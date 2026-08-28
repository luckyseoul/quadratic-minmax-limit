#!/usr/bin/env python3
"""Prop. 15.686 -- exclude the unique p=23 pair-slack-sixteen profile.

The unique slack-16 profile left after Propositions 15.684--15.685 has one
undetermined direction. Repairing a hypothetical 20-set takes at most four
deletions. Fewer than four deletions, followed by adjoining the undetermined
infinity point, gives an arc of size at least 18 and the conic-core
contradiction.

In the hard branch, write S=A union D with A a 16-arc and |D|=4, and adjoin
the undetermined infinity point U. Then K=A union {U} is a 17-arc. It must
be complete, or it extends to an 18-arc and gives the same conic-core
contradiction.

For d in D, the line Ud contains no second point of S, because U is
undetermined for S. Thus every secant of K through d is already a secant of
A. Completeness gives mu_A(d)=mu_K(d)>=1. The line-slack inequality gives

    16 = slack(S) >= 4 sum_{d in D} mu_A(d),

so all four points have secant multiplicity one outside K. Proposition
15.685 exhausts the five complete-17-arc classes and proves that no class has
more than one such outside point. This contradiction removes the profile.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15684 import (
    conic_core_repair_lemma,
    line_pair_slack,
    p23_reduction_theorem,
)
from e1_gmin_m4_prop15685 import (
    complete_17_arc_classification_certificate,
    p23_slack_twelve_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]
P = 23
BOUNDARY_SIZE = 20
PAIR_SLACK = 16


def slack_sixteen_repair_certificate() -> dict[str, object]:
    """Reduce the unique row to four 1-covered points of a complete 17-arc."""
    base = p23_reduction_theorem()
    rows = [
        row
        for row in base["exceptional_low_slack_profiles"]
        if int(row["pair_slack"]) == PAIR_SLACK
    ]
    if len(rows) != 1:
        raise ArithmeticError("unique slack-sixteen profile changed")
    profile = rows[0]
    if int(profile["undetermined_directions"]) != 1:
        raise ArithmeticError("slack-sixteen undetermined direction changed")

    repair = conic_core_repair_lemma()
    deletion_bound = PAIR_SLACK // 4
    if (
        deletion_bound != 4
        or int(repair["classification_threshold"]) != 18
        or min(
            int(row["pair_slack_floor"])
            for row in repair["off_conic_count_rows"].values()
        )
        != 24
    ):
        raise ArithmeticError("conic-core input changed")

    secant_line_rows = []
    for deleted_on_line in (1, 2, 3, 4):
        occupancy = 2 + deleted_on_line
        slack = line_pair_slack(occupancy)
        secant_line_rows.append(
            {
                "deleted_points_on_arc_secant": deleted_on_line,
                "final_line_occupancy": occupancy,
                "line_pair_slack": slack,
                "four_per_incidence_floor": 4 * deleted_on_line,
                "bound_holds": slack >= 4 * deleted_on_line,
            }
        )
    if not all(row["bound_holds"] for row in secant_line_rows):
        raise ArithmeticError("secant incidence floor changed")

    return {
        "unique_profile": profile,
        "repair_deletion_bound": deletion_bound,
        "undetermined_direction_count": 1,
        "at_most_three_deletions": (
            "adjoin the undetermined infinity point to obtain an arc of size "
            "at least 18; conic containment contradicts positive slack 16<24"
        ),
        "four_deletion_branch": {
            "repaired_arc_size": 16,
            "adjoined_undetermined_infinity_points": 1,
            "resulting_arc_size": 17,
            "if_resulting_arc_incomplete": (
                "extend to an 18-arc and apply the conic-core contradiction"
            ),
            "therefore_resulting_arc_complete": True,
        },
        "undetermined_direction_consequence": (
            "for each deleted d, line Ud contains no point of A; hence every "
            "secant of K=A union {U} through d is a secant of A"
        ),
        "secant_line_slack_rows": secant_line_rows,
        "global_incidence_bound": "slack(S)>=4*sum_{d in D} mu_A(d)",
        "completeness_floor": "mu_A(d)=mu_K(d)>=1",
        "equality_forced_by_slack_sixteen": (
            "four deleted points and slack 16 force mu_K(d)=1 for all four"
        ),
        "required_one_secant_points": 4,
        "proved": True,
    }


def p23_slack_sixteen_exclusion() -> dict[str, object]:
    """Proposition 15.686."""
    repair = slack_sixteen_repair_certificate()
    classification = complete_17_arc_classification_certificate()
    required = int(repair["required_one_secant_points"])
    available = int(classification["maximum_one_secant_point_count"])
    if not available < required:
        raise ArithmeticError("classification no longer excludes slack sixteen")

    previous = p23_slack_twelve_exclusion()
    before = int(previous["p23_profile_count_after"])
    after = before - 1
    remaining_histogram = dict(previous["remaining_pair_slack_histogram"])
    remaining_histogram.pop(PAIR_SLACK)
    if before != 202 or after != 201 or sum(remaining_histogram.values()) != after:
        raise ArithmeticError("post-15.686 profile accounting changed")

    return {
        "proposition": "15.686",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "excluded_pair_slack": PAIR_SLACK,
        "required_one_secant_points": required,
        "maximum_available_in_any_complete_17_arc_class": available,
        "slack_sixteen_profile_excluded": True,
        "p23_profile_count_before": before,
        "p23_profiles_excluded_here": 1,
        "p23_profile_count_after": after,
        "remaining_pair_slack_histogram": remaining_histogram,
        "all_remaining_profiles_have_pair_slack_at_least": min(
            remaining_histogram
        ),
        "p23_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "repair_certificate": repair,
        "classification_certificate": classification,
        "external_dependencies": previous["external_dependencies"],
        "proved": True,
    }


def main() -> None:
    theorem = p23_slack_sixteen_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15686.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.686: p=23 slack-16 profile excluded; "
        f"exact remainder {theorem['p23_profile_count_before']} -> "
        f"{theorem['p23_profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
