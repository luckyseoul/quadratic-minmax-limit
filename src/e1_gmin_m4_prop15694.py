#!/usr/bin/env python3
"""Prop. 15.694 -- exact equality normal form for p=19 slack twenty.

Proposition 15.693 forces every surviving slack-twenty boundary to use all
five repair deletions.  If ``S=A disjoint_union D`` is such a repair, then
``|A|=11`` and ``|D|=5``.  Restoring a deleted point that lies on no secant
of ``A`` would give a forbidden four-deletion repair.  Hence every deleted
point has positive secant multiplicity.  The line-slack inequality

    slack(S) >= 4 sum_{x in D} mu_A(x)

and total slack twenty force ``mu_A(x)=1`` for all five points and equality
on every affine line.

Writing ``a=|A cap l|`` and ``d=|D cap l|``, the repaired arc gives
``a<=2``.  Equality in the line bound permits exactly

    (a,d) in {(0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(2,2)}.

Thus ``D`` is also an affine arc, no line has more than four boundary
points, and a line through two deleted points contains zero or two core
points.  The five charged deleted/secant incidences split into either five
3-lines, one 4-line plus three 3-lines, or two 4-lines plus one 3-line.

For any two undetermined infinity points, adjoining them to ``A`` gives a
13-arc ``K``.  The five deleted points have secant index one relative to
``K``, as do the other ``t-2`` undetermined infinity points.  Therefore
``c1(K)>=5+(t-2)``, namely seven for ``t=4`` and eight for ``t=5``.  The
published all-13-arc maximum is nine, so this is a strict class filter but
not yet a contradiction.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15684 import line_pair_slack
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles
from e1_gmin_m4_prop15693 import p19_repair_depth_reduction


ROOT = Path(__file__).resolve().parents[1]
P = 19
BOUNDARY_SIZE = 16
PAIR_SLACK = 20
CORE_SIZE = 11
DELETED_SIZE = 5


def p19_thirteen_arc_secant_index_classification() -> dict[str, object]:
    """Published all-13-arc class count and c1 range."""
    return {
        "external_dependency": True,
        "source": (
            "E. B. Al-Zangana, The Geometry of the Plane of Order Nineteen "
            "and its Application to Error-Correcting Codes, PhD thesis, "
            "University of Sussex, 2011, Chapter 4, Section 4.21, pages 103-104"
        ),
        "projective_thirteen_arc_classes": 2733,
        "incomplete_thirteen_arc_classes": 501,
        "complete_thirteen_arc_classes": 2232,
        "minimum_c1_over_all_thirteen_arcs": 0,
        "maximum_c1_over_all_thirteen_arcs": 9,
        "complete_class_count_cross_check": (
            "H. Sticker, Classification of Arcs in Small Desarguesian "
            "Projective Planes, PhD thesis, Ghent University, 2012"
        ),
        "proved_conditional_on_external_classification": True,
    }


def p19_slack_twenty_profile_ledger() -> dict[str, object]:
    profiles = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) == PAIR_SLACK
    ]
    t_histogram = dict(
        sorted(Counter(int(row["undetermined_directions"]) for row in profiles).items())
    )
    if len(profiles) != 4 or t_histogram != {4: 2, 5: 2}:
        raise ArithmeticError("p=19 slack-twenty profile block changed")
    depth = p19_repair_depth_reduction()["rows"][PAIR_SLACK]
    if int(depth["minimum_repair_deletions_after_classification"]) != 5:
        raise ArithmeticError("slack-twenty repair depth changed")
    return {
        "pair_slack": PAIR_SLACK,
        "profile_count": len(profiles),
        "undetermined_direction_histogram": t_histogram,
        "repair_deletion_upper_bound": PAIR_SLACK // 4,
        "repair_deletion_lower_bound_from_15_693": 5,
        "therefore_exact_repair_deletions": 5,
        "core_size": CORE_SIZE,
        "deleted_size": DELETED_SIZE,
        "profiles": profiles,
        "proved": True,
    }


def p19_slack_twenty_line_equality_ledger() -> dict[str, object]:
    rows = []
    allowed = []
    for core_occupancy in range(3):
        for deleted_occupancy in range(DELETED_SIZE + 1):
            occupancy = core_occupancy + deleted_occupancy
            exact_slack = line_pair_slack(occupancy)
            charged_floor = (
                4 * deleted_occupancy if core_occupancy == 2 else 0
            )
            equality = exact_slack == charged_floor
            row = {
                "core_occupancy": core_occupancy,
                "deleted_occupancy": deleted_occupancy,
                "boundary_occupancy": occupancy,
                "exact_line_slack": exact_slack,
                "charged_secant_floor": charged_floor,
                "equality": equality,
            }
            rows.append(row)
            if equality:
                allowed.append((core_occupancy, deleted_occupancy))
    expected = [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (2, 0),
        (2, 1),
        (2, 2),
    ]
    if allowed != expected:
        raise ArithmeticError("slack-twenty line equality types changed")
    return {
        "line_rows": rows,
        "allowed_core_deleted_occupancies": [list(row) for row in allowed],
        "deleted_set_is_arc": True,
        "maximum_boundary_line_occupancy": 4,
        "two_deleted_points_force_core_occupancy_zero_or_two": True,
        "proved": True,
    }


def p19_slack_twenty_bad_line_patterns() -> list[dict[str, int]]:
    """Solutions r+2q=5 for 3-lines and 4-lines charged by D."""
    rows = []
    for four_lines in range(3):
        three_lines = DELETED_SIZE - 2 * four_lines
        rows.append(
            {
                "three_point_lines_core2_deleted1": three_lines,
                "four_point_lines_core2_deleted2": four_lines,
                "charged_deleted_secant_incidences": (
                    three_lines + 2 * four_lines
                ),
                "line_slack": 4 * three_lines + 8 * four_lines,
            }
        )
    if not all(
        row["charged_deleted_secant_incidences"] == DELETED_SIZE
        and row["line_slack"] == PAIR_SLACK
        for row in rows
    ):
        raise ArithmeticError("slack-twenty bad-line patterns changed")
    return rows


def p19_slack_twenty_equality_normal_form() -> dict[str, object]:
    """Proposition 15.694."""
    profiles = p19_slack_twenty_profile_ledger()
    line_types = p19_slack_twenty_line_equality_ledger()
    patterns = p19_slack_twenty_bad_line_patterns()
    classification = p19_thirteen_arc_secant_index_classification()

    c1_floors = {
        t: DELETED_SIZE + (t - 2)
        for t in profiles["undetermined_direction_histogram"]
    }
    if c1_floors != {4: 7, 5: 8}:
        raise ArithmeticError("repaired 13-arc c1 floors changed")
    if max(c1_floors.values()) > int(
        classification["maximum_c1_over_all_thirteen_arcs"]
    ):
        raise ArithmeticError("classification would now close a profile")

    return {
        "proposition": "15.694",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "pair_slack": PAIR_SLACK,
        "profile_ledger": profiles,
        "repair_normal_form": {
            "core_size": CORE_SIZE,
            "deleted_size": DELETED_SIZE,
            "core_is_arc": True,
            "deleted_is_arc": True,
            "deleted_core_secant_multiplicities": [1] * DELETED_SIZE,
            "global_slack_equality": (
                "slack(S)=4*sum_{x in D} mu_A(x)=20"
            ),
        },
        "line_equality": line_types,
        "bad_line_patterns": patterns,
        "adjoin_any_two_undetermined_infinity_points": {
            "resulting_arc_size": 13,
            "deleted_index_one_points": DELETED_SIZE,
            "unused_undetermined_index_one_points": "t-2",
            "c1_floors_by_t": c1_floors,
            "classified_maximum_c1": classification[
                "maximum_c1_over_all_thirteen_arcs"
            ],
            "strict_class_filter_but_not_contradiction": True,
        },
        "classification": classification,
        "bounded_sat_and_cpsat_trials": (
            "UNKNOWN results are diagnostics only and are not proposition evidence"
        ),
        "profile_count_before": 7,
        "profile_count_after": 7,
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p19_slack_twenty_equality_normal_form()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15694.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.694: p=19 slack-twenty equality normal form; "
        "four profiles remain in this block"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
