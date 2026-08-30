#!/usr/bin/env python3
"""Prop. 15.704 -- exclude p=17 slack-sixteen rows with a free direction.

Proposition 15.703 leaves 227 pair-slack-sixteen profiles.  Their numbers of
undetermined directions are ``{0:87,1:88,2:47,3:5}``.  Repairing a boundary
by at most four deletions gives an arc ``A``.

For two undetermined directions, adjoining their infinity points gives a
14-arc in the four-deletion branch.  If complete, each deleted point lies on
at least two core secants and forces slack at least 32; if incomplete, it
extends to the unique conic class and the off-conic floor is at least 20.
Three directions are handled by two overlapping pairs: both pair arcs must
extend to the same conic, which cannot contain three collinear infinity
points.

With one undetermined point ``U``, fewer than four deletions reach either a
conic or a complete 14-arc whose three deleted points force slack at least
24.  At depth four, ``K=A+U`` is a 13-arc.  If complete, slack sixteen would
require four outside secant-index-one points, but the eight complete classes
have at most three.  If incomplete, extend by ``j``.  The conic branch again
fails.  In the unique complete-14 branch, ``j`` in the deleted set forces
slack at least 24.  Otherwise all four deleted points must have index one
outside ``K``.  An exact audit of all 14 deletions, candidate quadruples, and
choices of ``U`` finds eight placements where ``U`` is genuinely
undetermined; every one has reconstructed slack 32.

Thus all 140 rows with at least one undetermined direction are impossible.
The p17 remainder drops from 1,368 to 1,228: two slack-zero rows, 87
slack-sixteen rows with no undetermined direction, and 1,139 rows of slack at
least twenty.  The endpoint remains open.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path

from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15702 import (
    COMPLETE_14_ARC,
    complete_fourteen_arc_secant_index_certificate,
)
from e1_gmin_m4_prop15703 import (
    _arc_data,
    _set_pair_slack,
    complete_thirteen_arc_certificate,
    p17_slack_twelve_exclusion,
)


ROOT = Path(__file__).resolve().parents[1]
P = 17
BOUNDARY_SIZE = 16
PAIR_SLACK = 16
Point = tuple[int, int, int]


def _undetermined_directions(row: dict[str, object]) -> int:
    profiles = row["phase_profiles_b"]
    return sum(int(profiles[phase].get(BOUNDARY_SIZE, 0)) for phase in ("0", "1"))


def _determinant(a: Point, b: Point, c: Point) -> int:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % P


def _is_undetermined_for(boundary: set[Point], infinity: Point) -> bool:
    """No chord of the boundary has projective direction ``infinity``."""
    return all(
        _determinant(infinity, first, second)
        for first, second in combinations(boundary, 2)
    )


def _line_pair_slack(occupancy: int) -> int:
    return 2 * (math.comb(occupancy, 2) - occupancy // 2)


def slack_sixteen_profile_ledger() -> dict[str, object]:
    """Split the exact 227-row block by undetermined-direction count."""
    profiles = p17_second_boundary_profile_census()["profiles"]
    indexed_profiles = [
        (index, row)
        for index, row in enumerate(profiles)
        if int(row["pair_slack"]) == PAIR_SLACK
    ]
    histogram = dict(
        sorted(
            Counter(
                _undetermined_directions(row) for _index, row in indexed_profiles
            ).items()
        )
    )
    excluded_indices = [
        index
        for index, row in indexed_profiles
        if _undetermined_directions(row) >= 1
    ]
    remaining_indices = [
        index
        for index, row in indexed_profiles
        if _undetermined_directions(row) == 0
    ]
    if len(indexed_profiles) != 227 or histogram != {0: 87, 1: 88, 2: 47, 3: 5}:
        raise ArithmeticError("p17 slack-sixteen profile ledger changed")
    return {
        "profile_count": len(indexed_profiles),
        "undetermined_direction_histogram": histogram,
        "profiles_with_at_least_one_undetermined_direction": len(indexed_profiles)
        - histogram[0],
        "remaining_zero_direction_profiles": histogram[0],
        "excluded_profile_indices": excluded_indices,
        "remaining_zero_direction_profile_indices": remaining_indices,
        "proved": True,
    }


def slack_sixteen_repair_lemma() -> dict[str, object]:
    """Audit the incidence and conic floors used in every branch."""
    line_rows = []
    for deleted_on_secant in range(1, 5):
        occupancy = 2 + deleted_on_secant
        slack = _line_pair_slack(occupancy)
        line_rows.append(
            {
                "deleted_points_on_core_secant": deleted_on_secant,
                "line_occupancy": occupancy,
                "line_pair_slack": slack,
                "four_per_incidence_floor": 4 * deleted_on_secant,
                "bound_holds": slack >= 4 * deleted_on_secant,
            }
        )
    conic_rows = {
        h: {
            "off_conic_points": h,
            "retained_conic_secants_per_off_conic_point": 6 - h,
            "pair_slack_floor": 4 * h * (6 - h),
        }
        for h in range(1, 5)
    }
    if not all(row["bound_holds"] for row in line_rows):
        raise ArithmeticError("p17 secant incidence floor changed")
    if [conic_rows[h]["pair_slack_floor"] for h in range(1, 5)] != [20, 32, 36, 32]:
        raise ArithmeticError("p17 four-point conic floor changed")
    return {
        "repair_deletion_bound": PAIR_SLACK // 4,
        "line_secant_incidence_rows": line_rows,
        "global_incidence_bound": "slack(S) >= 4 sum_{d in D} mu_A(d)",
        "off_conic_rows": conic_rows,
        "minimum_positive_conic_core_slack": min(
            int(row["pair_slack_floor"]) for row in conic_rows.values()
        ),
        "proved": True,
    }


def one_direction_complete_arc_certificate() -> dict[str, object]:
    """Audit complete-13 and complete-14-minus-one hard branches."""
    complete13 = complete_thirteen_arc_certificate()
    c1_counts = list(complete13["index_one_point_counts"])
    if c1_counts != [0, 0, 0, 0, 0, 0, 2, 3]:
        raise ArithmeticError("p17 complete-13 index-one census changed")

    complete14_certificate = complete_fourteen_arc_secant_index_certificate()
    if int(complete14_certificate["minimum_outside_secant_index"]) != 2:
        raise ArithmeticError("p17 complete-14 secant floor changed")

    complete14 = set(COMPLETE_14_ARC)
    c1_count_histogram: Counter[int] = Counter()
    raw_candidate_quadruples = 0
    raw_infinity_placements: Counter[int] = Counter()
    undetermined_placements: Counter[int] = Counter()
    for extension in COMPLETE_14_ARC:
        core = complete14 - {extension}
        data = _arc_data(core)
        index_one = [
            point
            for point, index in data["outside_indices"].items()
            if index == 1
        ]
        c1_count_histogram[len(index_one)] += 1
        for deleted in combinations(index_one, 4):
            raw_candidate_quadruples += 1
            for infinity in core:
                boundary = (core - {infinity}) | set(deleted)
                slack = _set_pair_slack(boundary)
                raw_infinity_placements[slack] += 1
                if _is_undetermined_for(boundary, infinity):
                    undetermined_placements[slack] += 1
    if (
        dict(sorted(c1_count_histogram.items())) != {0: 4, 1: 8, 4: 2}
        or raw_candidate_quadruples != 2
        or dict(sorted(raw_infinity_placements.items())) != {16: 2, 28: 16, 32: 8}
        or dict(sorted(undetermined_placements.items())) != {32: 8}
    ):
        raise ArithmeticError("p17 one-direction complete-14 census changed")
    return {
        "complete_13_class_count": int(complete13["class_count"]),
        "complete_13_index_one_point_counts": c1_counts,
        "maximum_complete_13_index_one_point_count": max(c1_counts),
        "complete_14_minus_one_index_one_count_histogram": dict(
            sorted(c1_count_histogram.items())
        ),
        "raw_candidate_index_one_quadruples": raw_candidate_quadruples,
        "raw_candidate_infinity_placement_slack_histogram": dict(
            sorted(raw_infinity_placements.items())
        ),
        "genuinely_undetermined_infinity_placement_slack_histogram": dict(
            sorted(undetermined_placements.items())
        ),
        "genuinely_undetermined_slack_sixteen_placements": int(
            undetermined_placements[PAIR_SLACK]
        ),
        "proved_conditional_on_published_complete_arc_class_counts": True,
    }


def p17_slack_sixteen_free_direction_exclusion() -> dict[str, object]:
    """Proposition 15.704."""
    previous = p17_slack_twelve_exclusion()
    census = p17_second_boundary_profile_census()
    profiles = slack_sixteen_profile_ledger()
    repair = slack_sixteen_repair_lemma()
    arcs = one_direction_complete_arc_certificate()
    if (
        int(repair["minimum_positive_conic_core_slack"]) <= PAIR_SLACK
        or int(arcs["maximum_complete_13_index_one_point_count"]) >= 4
        or int(arcs["genuinely_undetermined_slack_sixteen_placements"]) != 0
    ):
        raise ArithmeticError("p17 slack-sixteen obstruction changed")

    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    excluded_indices = set(int(index) for index in profiles["excluded_profile_indices"])
    if not excluded_indices <= previous_indices:
        raise ArithmeticError("15.704 tried to exclude an already absent profile")
    remaining_indices = sorted(previous_indices - excluded_indices)
    before = len(previous_indices)
    excluded = len(excluded_indices)
    after = len(remaining_indices)
    histogram = dict(
        sorted(
            Counter(
                int(census["profiles"][index]["pair_slack"])
                for index in remaining_indices
            ).items()
        )
    )
    if (
        before != 1368
        or excluded != 140
        or after != 1228
        or histogram.get(PAIR_SLACK) != 87
        or sum(histogram.values()) != after
    ):
        raise ArithmeticError("p17 post-slack-sixteen accounting changed")
    return {
        "proposition": "15.704",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(excluded_indices),
        "remaining_profile_indices": remaining_indices,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 2,
        "remaining_slack_sixteen_profiles": 87,
        "remaining_profiles_of_slack_at_least_twenty": after - 15,
        "two_direction_branch": (
            "the four-deletion pair arc has size 14; if complete its four "
            "deleted points force slack at least 32, and if incomplete its "
            "conic extension forces positive slack at least 20"
        ),
        "three_direction_branch": (
            "two overlapping infinity-point pairs must extend to the same "
            "conic, which cannot contain three collinear infinity points"
        ),
        "one_direction_branch": (
            "a complete 13-arc has at most three index-one outside points; "
            "all eight valid complete-14-minus-one placements have slack 32"
        ),
        "profile_ledger": profiles,
        "repair_lemma": repair,
        "one_direction_complete_arc_certificate": arcs,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p17_slack_sixteen_free_direction_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15704.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.704: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
