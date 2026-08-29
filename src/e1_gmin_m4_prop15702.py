#!/usr/bin/env python3
"""Prop. 15.702 -- complete-14-arc exclusion at p=17.

Proposition 15.701 leaves 932 profiles at the second all-finite boundary
``p=17,s=16``.  Its unresolved low-slack block consists of 67 slack-eight
profiles with no undetermined direction and 112 slack-twelve profiles, 79 of
which have one undetermined direction.

Sticker's exhaustive complete-arc classification records one complete
14-arc class in ``PG(2,17)``.  An explicit representative is audited here.
Its outside-point secant-index histogram is

    {2:4, 3:4, 4:76, 5:128, 6:75, 7:6};

in particular it has no outside point of secant index zero or one.

For slack eight, repair deletes at most two points.  A repair with fewer than
two deletions reaches a conic core and is already excluded by 15.701.  A
14-arc repair is either incomplete, hence extends to a conic-contained
15-arc, or is the unique complete class.  In the latter case each of the two
deleted points lies on at least two core secants and forces slack at least
``4*(2+2)=16``.  Thus all remaining slack-eight profiles are impossible.

For slack twelve, any surviving repair must delete exactly three points.  If
the profile has an undetermined direction, adjoining its point at infinity
turns the repaired 13-arc into a 14-arc.  An incomplete 14-arc again reaches
a conic.  If it is complete, the infinity direction guarantees that no line
from a deleted affine point to the adjoined infinity point contains a core
point.  Consequently all complete-14-arc secants through each deleted point
are already core secants, forcing slack at least ``4*3*2=24``.

This excludes 67+79=146 additional profiles and leaves 786: two at slack
zero, 33 at slack twelve, and 751 at slack at least sixteen.  The endpoint
remains open.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15701 import p17_low_positive_slack_conic_reduction


ROOT = Path(__file__).resolve().parents[1]
P = 17
BOUNDARY_SIZE = 16
Point = tuple[int, int, int]


COMPLETE_14_ARC: tuple[Point, ...] = (
    (1, 0, 0),
    (1, 1, 1),
    (1, 2, 7),
    (1, 4, 5),
    (1, 5, 16),
    (1, 7, 10),
    (1, 8, 2),
    (1, 9, 12),
    (1, 10, 8),
    (1, 12, 11),
    (1, 13, 14),
    (1, 15, 13),
    (0, 1, 0),
    (0, 0, 1),
)


def _projective_points_or_lines() -> tuple[Point, ...]:
    return tuple(
        [(1, y, z) for y in range(P) for z in range(P)]
        + [(0, 1, z) for z in range(P)]
        + [(0, 0, 1)]
    )


def _incident(point: Point, line: Point) -> bool:
    return sum(a * b for a, b in zip(point, line)) % P == 0


def p17_complete_fourteen_arc_classification() -> dict[str, object]:
    """Record the published unique complete-14-arc class."""
    return {
        "external_dependency": True,
        "source": (
            "H. Sticker, Classification of Arcs in Small Desarguesian "
            "Projective Planes, PhD thesis, Ghent University, 2012"
        ),
        "source_url": (
            "https://cage.ugent.be/geometry/Theses/57/PhDHeideSticker.pdf"
        ),
        "location": "Section 5.1, printed page 102 (PDF page 111)",
        "complete_14_arc_class_count": 1,
        "reported_automorphism_group": "D8",
        "classification_scope": "PGL-inequivalent complete (k,2)-arcs",
        "proved_conditional_on_external_classification": True,
    }


def complete_fourteen_arc_secant_index_certificate() -> dict[str, object]:
    """Audit one complete 14-arc and every outside point's secant index."""
    points = _projective_points_or_lines()
    lines = points
    arc = set(COMPLETE_14_ARC)
    if len(points) != 307 or len(arc) != 14 or not arc <= set(points):
        raise ArithmeticError("p=17 complete-14-arc point ledger changed")

    occupancies = [sum(point in arc for point in points if _incident(point, line)) for line in lines]
    occupancy_histogram = dict(sorted(Counter(occupancies).items()))
    if occupancy_histogram != {0: 146, 1: 70, 2: 91}:
        raise ArithmeticError("p=17 complete-14-arc line census changed")
    secants = [line for line, occupancy in zip(lines, occupancies) if occupancy == 2]
    outside_indices = {
        point: sum(_incident(point, line) for line in secants)
        for point in points
        if point not in arc
    }
    index_histogram = dict(sorted(Counter(outside_indices.values()).items()))
    expected = {2: 4, 3: 4, 4: 76, 5: 128, 6: 75, 7: 6}
    if index_histogram != expected or min(outside_indices.values()) != 2:
        raise ArithmeticError("p=17 complete-14-arc secant indices changed")
    return {
        "representative": [list(point) for point in COMPLETE_14_ARC],
        "projective_point_and_line_count": len(points),
        "line_occupancy_histogram": occupancy_histogram,
        "outside_point_count": len(outside_indices),
        "outside_secant_index_histogram": index_histogram,
        "minimum_outside_secant_index": min(outside_indices.values()),
        "complete": 0 not in index_histogram,
        "no_index_one_outside_points": 1 not in index_histogram,
        "proved": True,
    }


def _undetermined_directions(row: dict[str, object]) -> int:
    profiles = row["phase_profiles_b"]
    return sum(int(profiles[phase].get(BOUNDARY_SIZE, 0)) for phase in ("0", "1"))


def p17_complete_fourteen_arc_profile_ledger() -> dict[str, object]:
    """Extract the 67+79 rows newly excluded by the certificate."""
    profiles = p17_second_boundary_profile_census()["profiles"]
    slack_eight = [row for row in profiles if int(row["pair_slack"]) == 8]
    slack_twelve = [row for row in profiles if int(row["pair_slack"]) == 12]
    t8 = dict(sorted(Counter(_undetermined_directions(row) for row in slack_eight).items()))
    t12 = dict(sorted(Counter(_undetermined_directions(row) for row in slack_twelve).items()))
    newly_eight = [row for row in slack_eight if _undetermined_directions(row) == 0]
    newly_twelve = [row for row in slack_twelve if _undetermined_directions(row) == 1]
    if t8 != {0: 67, 1: 104, 2: 24} or t12 != {0: 33, 1: 79, 2: 43}:
        raise ArithmeticError("p=17 low-slack undetermined ledger changed")
    if len(newly_eight) != 67 or len(newly_twelve) != 79:
        raise ArithmeticError("p=17 complete-14-arc exclusion count changed")
    return {
        "slack_eight": {
            "all_profile_count": len(slack_eight),
            "undetermined_direction_histogram": t8,
            "previously_excluded": len(slack_eight) - len(newly_eight),
            "newly_excluded": len(newly_eight),
            "remaining": 0,
        },
        "slack_twelve": {
            "all_profile_count": len(slack_twelve),
            "undetermined_direction_histogram": t12,
            "previously_excluded": t12[2],
            "newly_excluded": len(newly_twelve),
            "remaining_without_undetermined_direction": t12[0],
        },
        "newly_excluded_profile_count": len(newly_eight) + len(newly_twelve),
        "proved": True,
    }


def p17_complete_fourteen_arc_exclusion() -> dict[str, object]:
    """Proposition 15.702."""
    previous = p17_low_positive_slack_conic_reduction()
    census = p17_second_boundary_profile_census()
    classification = p17_complete_fourteen_arc_classification()
    certificate = complete_fourteen_arc_secant_index_certificate()
    ledger = p17_complete_fourteen_arc_profile_ledger()
    if int(classification["complete_14_arc_class_count"]) != 1:
        raise ArithmeticError("p=17 complete-14-arc class count changed")
    if int(certificate["minimum_outside_secant_index"]) != 2:
        raise ArithmeticError("p=17 complete-14-arc minimum index changed")

    before = int(previous["profile_count_after"])
    excluded = int(ledger["newly_excluded_profile_count"])
    after = before - excluded
    histogram = dict(previous["remaining_pair_slack_histogram"])
    histogram.pop(8)
    histogram[12] = 33
    histogram = dict(sorted(histogram.items()))
    high_slack_count = sum(count for slack, count in histogram.items() if int(slack) >= 16)
    if (
        before != 932
        or excluded != 146
        or after != 786
        or high_slack_count != 751
        or sum(histogram.values()) != after
    ):
        raise ArithmeticError("p=17 post-complete-14-arc accounting changed")
    return {
        "proposition": "15.702",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 2,
        "remaining_slack_twelve_profiles": 33,
        "remaining_profiles_of_slack_at_least_sixteen": high_slack_count,
        "slack_eight_exclusion": (
            "repair depth two gives a complete 14-arc whose two deleted "
            "points each have secant index at least two, forcing slack >=16"
        ),
        "slack_twelve_one_direction_exclusion": (
            "three-point repair plus the undetermined infinity point gives a "
            "complete 14-arc; each deleted point retains secant index at least "
            "two, forcing slack >=24"
        ),
        "profile_ledger": ledger,
        "complete_14_arc_certificate": certificate,
        "external_classification": classification,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p17_complete_fourteen_arc_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15702.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.702: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
