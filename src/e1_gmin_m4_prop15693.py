#!/usr/bin/env python3
"""Prop. 15.693 -- exclude all seven p=19 slack-sixteen profiles.

Every surviving slack-sixteen boundary has three or four undetermined
directions.  Repair deletes at most four points.  With at most three
deletions, adjoining two undetermined infinity points gives an arc of size
at least fifteen and hence a conic, contradicting a third undetermined
direction.

In the four-deletion branch, let ``A`` be the repaired 12-arc and adjoin two
undetermined infinity points to obtain a 14-arc ``K``.  If ``K`` is
incomplete, it extends through the classified gap to the 20-point conic.
The third undetermined infinity point has exactly one ``K``-secant (the
line at infinity), whereas conic secant counting leaves at least three, a
contradiction.  Thus ``K`` is complete.

Each deleted point lies on at least one secant of ``K``.  Undeterminedness
prevents those secants from using either adjoined infinity point, so every
such secant is already a secant of ``A`` and charges four units of boundary
slack.  Total slack sixteen forces all four deleted points to have secant
multiplicity one.  Every unused undetermined infinity point also has
secant multiplicity exactly one with respect to ``K``: its sole secant is
the line at infinity through the two adjoined points.  Hence ``K`` has at
least ``4+(t-2)>=5`` outside points of secant index one.

Al-Zangana's exhaustive PG(2,19) classification gives 83 projective classes
of 14-arcs, 70 complete, and at most four index-one outside points for any
14-arc.  This contradiction excludes all seven profiles.  The same count
forces every slack-twenty survivor to use all five allowed repair deletions,
and gives smaller deletion-depth reductions at slack 24,28,32.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15688 import p19_residue_zero_profiles
from e1_gmin_m4_prop15689 import p19_low_slack_geometric_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 19
BOUNDARY_SIZE = 16
PAIR_SLACK = 16


def p19_fourteen_arc_secant_index_classification() -> dict[str, object]:
    """The exhaustive external class bound used in the contradiction."""
    return {
        "external_dependency": True,
        "source": (
            "E. B. Al-Zangana, The Geometry of the Plane of Order Nineteen "
            "and its Application to Error-Correcting Codes, PhD thesis, "
            "University of Sussex, 2011, Chapter 4, Section 4.22, page 105"
        ),
        "projective_fourteen_arc_classes": 83,
        "incomplete_fourteen_arc_classes": 13,
        "complete_fourteen_arc_classes": 70,
        "outside_secant_index_notation": (
            "c_i is the number of outside points on exactly i arc secants"
        ),
        "maximum_c1_over_all_fourteen_arcs": 4,
        "therefore_maximum_c1_over_complete_fourteen_arcs": 4,
        "class_count_cross_check": (
            "H. Sticker, Classification of Arcs in Small Desarguesian "
            "Projective Planes, PhD thesis, Ghent University, 2012"
        ),
        "proved_conditional_on_external_classification": True,
    }


def p19_slack_sixteen_profile_ledger() -> dict[str, object]:
    profiles = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) == PAIR_SLACK
    ]
    t_histogram = dict(
        sorted(Counter(int(row["undetermined_directions"]) for row in profiles).items())
    )
    if len(profiles) != 7 or t_histogram != {3: 1, 4: 6}:
        raise ArithmeticError("p=19 slack-sixteen block changed")
    return {
        "pair_slack": PAIR_SLACK,
        "profile_count": len(profiles),
        "undetermined_direction_histogram": t_histogram,
        "repair_deletion_bound": PAIR_SLACK // 4,
        "proved": True,
    }


def _c1_floor_for_four_deletions(slack: int, undetermined: int) -> int:
    """Minimum index-one points on a complete repaired 14-arc.

    Four deleted points have positive integral secant multiplicities with
    total at most slack/4.  At least ``8-slack/4`` of them equal one.  The
    unused undetermined infinity points add ``undetermined-2`` more.
    """
    if slack % 4 or slack < 16 or undetermined < 3:
        raise ValueError("invalid four-deletion row")
    deleted_index_one = max(0, 8 - slack // 4)
    return deleted_index_one + undetermined - 2


def p19_repair_depth_reduction() -> dict[str, object]:
    """Record what the c1<=4 classification forces at every live slack."""
    census = p19_residue_zero_profiles()
    rows = {}
    expected = {
        16: ({3: 1, 4: 6}, None),
        20: ({4: 2, 5: 2}, 5),
        24: ({5: 1}, 5),
        28: ({5: 1}, 4),
        32: ({5: 1}, 4),
    }
    for slack, (expected_t, minimum_deletions) in expected.items():
        profiles = [
            row
            for row in census["profiles"]
            if int(row["pair_slack"]) == slack
        ]
        t_histogram = dict(
            sorted(
                Counter(
                    int(row["undetermined_directions"]) for row in profiles
                ).items()
            )
        )
        if t_histogram != expected_t:
            raise ArithmeticError("p=19 repair-depth profile changed")
        c1_floors = {
            t: _c1_floor_for_four_deletions(slack, t)
            for t in t_histogram
        }
        rows[slack] = {
            "profile_count": len(profiles),
            "undetermined_direction_histogram": t_histogram,
            "four_deletion_c1_floors": c1_floors,
            "four_deletion_excluded_for_every_profile": all(
                floor > 4 for floor in c1_floors.values()
            ),
            "minimum_repair_deletions_after_classification": minimum_deletions,
        }
    if not rows[20]["four_deletion_excluded_for_every_profile"]:
        raise ArithmeticError("slack-twenty depth reduction changed")
    if not rows[24]["four_deletion_excluded_for_every_profile"]:
        raise ArithmeticError("slack-twenty-four depth reduction changed")
    return {"rows": rows, "proved": True}


def p19_slack_sixteen_exclusion() -> dict[str, object]:
    """Proposition 15.693."""
    previous = p19_low_slack_geometric_exclusion()
    classification = p19_fourteen_arc_secant_index_classification()
    ledger = p19_slack_sixteen_profile_ledger()
    depths = p19_repair_depth_reduction()

    before = int(previous["profile_count_after"])
    removed = int(ledger["profile_count"])
    after = before - removed
    remaining_histogram = dict(previous["remaining_pair_slack_histogram"])
    remaining_histogram.pop(PAIR_SLACK)
    if before != 14 or removed != 7 or after != 7:
        raise ArithmeticError("post-15.693 profile count changed")
    if remaining_histogram != {20: 4, 24: 1, 28: 1, 32: 1}:
        raise ArithmeticError("post-15.693 histogram changed")

    minimum_c1 = min(
        4 + (t - 2) for t in ledger["undetermined_direction_histogram"]
    )
    maximum_c1 = int(
        classification["maximum_c1_over_all_fourteen_arcs"]
    )
    if minimum_c1 <= maximum_c1:
        raise ArithmeticError("14-arc secant-index contradiction changed")

    return {
        "proposition": "15.693",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "excluded_pair_slack": PAIR_SLACK,
        "profile_count_before": before,
        "profile_count_excluded": removed,
        "profile_count_after": after,
        "remaining_pair_slack_histogram": remaining_histogram,
        "at_most_three_deletions": (
            "two undetermined infinity points give an arc of size at least "
            "15 and hence a conic; a third undetermined infinity point has "
            "only the infinity-line secant but retained-conic counting forces more"
        ),
        "four_deletion_branch": {
            "repaired_arc_size": 12,
            "adjoined_infinity_points": 2,
            "resulting_arc_size": 14,
            "resulting_arc_must_be_complete": True,
            "deleted_index_one_points": 4,
            "unused_undetermined_index_one_points": "t-2",
            "minimum_total_index_one_points": minimum_c1,
            "classified_maximum": maximum_c1,
            "excluded": True,
        },
        "slack_twenty_profiles_now_force_exactly_five_repair_deletions": True,
        "repair_depth_reduction": depths,
        "profile_ledger": ledger,
        "classification": classification,
        "p19_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p19_slack_sixteen_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15693.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.693: p=19 slack-sixteen profiles excluded; "
        f"exact remainder {theorem['profile_count_before']} -> "
        f"{theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
