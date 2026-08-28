#!/usr/bin/env python3
"""Prop. 15.685 -- exclude the unique p=23 pair-slack-twelve profile.

Proposition 15.684 leaves one arithmetic profile of pair slack 12 at the
second all-finite endpoint p=23,s=20. Its exclusion reduces the exact
remainder from 203 profiles to 202.

The repair lemma of 15.684 deletes at most three points from a hypothetical
boundary S to obtain an arc. Two deletions would already give an 18-arc and
hence a conic core, contradicting slack below 24. Thus the only case to
consider is S=A union D with A a 17-arc and |D|=3. If A were incomplete,
adjoining one point would again produce an 18-point conic core, so A must
be complete.

For an outside point x, let mu_A(x) be the number of secants of A through x.
Every secant containing r points of D contributes at least 4r pair slack.
Completeness gives mu_A(x)>=1 for every x outside A. Since the total slack
is 12, all three deleted points would have to satisfy mu_A(x)=1.

Coolsaet--Sticker classify exactly five projective classes of complete
17-arcs in PG(2,23). The explicit representatives below are checked directly
to be complete arcs. Their full outside-point secant-multiplicity histograms
are pairwise distinct, so they represent five inequivalent classes and hence
exhaust the classification. They contain respectively 0, 0, 1, 0, and 0
points with mu_A=1. No class supplies the required three points.

This is an exact classification-assisted proof, not a heuristic search.
The coordinates, arc/completeness checks, all 536 outside-point
multiplicities per representative, and the classification count are
recorded in the generated evidence file.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

from e1_gmin_m4_prop15684 import (
    conic_core_repair_lemma,
    line_pair_slack,
    p23_complete_arc_classification,
    p23_reduction_theorem,
)


ROOT = Path(__file__).resolve().parents[1]
P = 23
ARC_SIZE = 17
BOUNDARY_SIZE = 20
PAIR_SLACK = 12

Point = tuple[int, int, int]


# Five explicit complete 17-arcs. Coordinates are homogeneous over F_23.
COMPLETE_17_ARC_REPRESENTATIVES: tuple[tuple[Point, ...], ...] = (
    (
        (19, 19, 1), (3, 10, 1), (19, 4, 1), (8, 18, 1), (6, 2, 1),
        (1, 1, 1), (22, 10, 1), (1, 22, 0), (17, 12, 1), (5, 12, 1),
        (3, 7, 1), (20, 7, 1), (1, 17, 0), (9, 5, 1), (20, 1, 1),
        (5, 2, 1), (16, 19, 1),
    ),
    (
        (20, 20, 1), (15, 6, 1), (9, 11, 1), (17, 13, 1), (17, 1, 1),
        (14, 0, 1), (2, 10, 1), (12, 21, 1), (4, 11, 1), (11, 6, 1),
        (14, 7, 1), (13, 10, 1), (10, 13, 1), (11, 19, 1), (0, 18, 1),
        (13, 7, 1), (0, 1, 1),
    ),
    (
        (16, 18, 1), (22, 1, 1), (15, 10, 1), (8, 12, 1), (11, 9, 1),
        (10, 1, 1), (6, 9, 1), (10, 4, 1), (12, 10, 1), (15, 21, 1),
        (22, 13, 1), (19, 16, 1), (9, 6, 1), (20, 5, 1), (9, 5, 1),
        (18, 13, 1), (19, 12, 1),
    ),
    (
        (22, 19, 1), (21, 1, 1), (7, 19, 1), (4, 11, 1), (10, 20, 1),
        (19, 4, 1), (2, 22, 1), (2, 3, 1), (12, 1, 1), (22, 8, 1),
        (9, 2, 1), (13, 11, 1), (8, 15, 1), (4, 2, 1), (1, 12, 1),
        (1, 10, 1), (9, 13, 1),
    ),
    (
        (16, 15, 1), (8, 8, 1), (18, 13, 1), (15, 11, 1), (16, 2, 1),
        (14, 21, 1), (6, 15, 1), (3, 10, 1), (21, 14, 1), (15, 3, 1),
        (11, 10, 1), (11, 17, 1), (3, 21, 1), (21, 19, 1), (6, 11, 1),
        (9, 17, 1), (13, 9, 1),
    ),
)


EXPECTED_SECANT_MULTIPLICITY_HISTOGRAMS: tuple[dict[int, int], ...] = (
    {2: 2, 3: 6, 4: 68, 5: 172, 6: 190, 7: 86, 8: 12},
    {2: 1, 3: 15, 4: 59, 5: 159, 6: 208, 7: 86, 8: 8},
    {1: 1, 3: 6, 4: 69, 5: 171, 6: 196, 7: 78, 8: 15},
    {3: 14, 4: 58, 5: 170, 6: 206, 7: 72, 8: 16},
    {2: 1, 3: 8, 4: 63, 5: 185, 6: 176, 7: 91, 8: 12},
)


def normalize_projective(vector: Iterable[int], p: int = P) -> Point:
    """Normalize a nonzero homogeneous vector by its first nonzero entry."""
    row = tuple(int(value) % p for value in vector)
    for value in row:
        if value:
            inverse = pow(value, -1, p)
            return tuple((entry * inverse) % p for entry in row)  # type: ignore[return-value]
    raise ValueError("the zero vector is not a projective point")


def projective_points(p: int = P) -> tuple[Point, ...]:
    """All points of PG(2,p), once each in first-nonzero-one form."""
    rows = (
        [(1, y, z) for y in range(p) for z in range(p)]
        + [(0, 1, z) for z in range(p)]
        + [(0, 0, 1)]
    )
    if len(rows) != p * p + p + 1 or len(set(rows)) != len(rows):
        raise ArithmeticError("projective point enumeration changed")
    return tuple(rows)


def line_through(first: Point, second: Point, p: int = P) -> Point:
    """Normalized homogeneous line through two distinct projective points."""
    cross = (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )
    return normalize_projective(cross, p)


def incident(point: Point, line: Point, p: int = P) -> bool:
    return sum(point[index] * line[index] for index in range(3)) % p == 0


def complete_arc_certificate(points: Iterable[Point], p: int = P) -> dict[str, object]:
    """Verify an arc, completeness, and every outside secant multiplicity."""
    arc = tuple(sorted({normalize_projective(point, p) for point in points}))
    if len(arc) != ARC_SIZE:
        raise ArithmeticError("representative does not have 17 distinct points")

    secants = Counter(
        line_through(arc[first], arc[second], p)
        for first in range(len(arc))
        for second in range(first)
    )
    is_arc = len(secants) == math.comb(ARC_SIZE, 2) and max(secants.values()) == 1
    if not is_arc:
        raise ArithmeticError("representative contains three collinear points")

    all_points = projective_points(p)
    arc_set = set(arc)
    multiplicities = {
        point: sum(incident(point, line, p) for line in secants)
        for point in all_points
        if point not in arc_set
    }
    histogram = dict(sorted(Counter(multiplicities.values()).items()))
    complete = min(multiplicities.values()) >= 1
    if (
        len(secants) != 136
        or len(multiplicities) != p * p + p + 1 - ARC_SIZE
        or sum(histogram.values()) != 536
        or sum(key * value for key, value in histogram.items())
        != math.comb(ARC_SIZE, 2) * (p - 1)
    ):
        raise ArithmeticError("secant incidence accounting changed")

    return {
        "coordinates": [list(point) for point in arc],
        "point_count": len(arc),
        "secant_line_count": len(secants),
        "outside_point_count": len(multiplicities),
        "secant_multiplicity_histogram": histogram,
        "one_secant_point_count": histogram.get(1, 0),
        "minimum_outside_secant_multiplicity": min(multiplicities.values()),
        "maximum_outside_secant_multiplicity": max(multiplicities.values()),
        "secant_outside_incidence_count": sum(multiplicities.values()),
        "is_arc": is_arc,
        "is_complete": complete,
    }


def complete_17_arc_classification_certificate() -> dict[str, object]:
    """Exhaust the five classified complete-17-arc classes."""
    external = p23_complete_arc_classification()
    if external["complete_arc_counts"][17] != 5:
        raise ArithmeticError("external complete-17-arc class count changed")

    representatives = [
        complete_arc_certificate(points)
        for points in COMPLETE_17_ARC_REPRESENTATIVES
    ]
    observed = [
        row["secant_multiplicity_histogram"] for row in representatives
    ]
    if observed != list(EXPECTED_SECANT_MULTIPLICITY_HISTOGRAMS):
        raise ArithmeticError("complete-17-arc multiplicity ledger changed")
    fingerprints = {
        tuple(sorted(histogram.items()))
        for histogram in EXPECTED_SECANT_MULTIPLICITY_HISTOGRAMS
    }
    if len(fingerprints) != 5 or not all(
        row["is_arc"] and row["is_complete"] for row in representatives
    ):
        raise ArithmeticError("representatives do not certify five classes")

    one_secant_counts = [
        int(row["one_secant_point_count"]) for row in representatives
    ]
    return {
        "classification_source": external,
        "classified_projective_class_count": 5,
        "verified_representative_count": len(representatives),
        "inequivalence_invariant": (
            "the full outside-point secant-multiplicity histogram is preserved "
            "by projective equivalence"
        ),
        "pairwise_distinct_invariants": len(fingerprints) == 5,
        "therefore_exhaustive": len(representatives)
        == int(external["complete_arc_counts"][17]),
        "one_secant_point_counts_by_class": one_secant_counts,
        "maximum_one_secant_point_count": max(one_secant_counts),
        "representatives": representatives,
        "proved": True,
    }


def slack_twelve_repair_certificate() -> dict[str, object]:
    """Reduce a slack-12 set to three 1-covered points of a complete 17-arc."""
    theorem = p23_reduction_theorem()
    rows = [
        row
        for row in theorem["exceptional_low_slack_profiles"]
        if int(row["pair_slack"]) == PAIR_SLACK
    ]
    if len(rows) != 1:
        raise ArithmeticError("unique slack-twelve profile changed")
    profile = rows[0]
    repair = conic_core_repair_lemma()
    if (
        int(repair["classification_threshold"]) != 18
        or min(
            int(row["pair_slack_floor"])
            for row in repair["off_conic_count_rows"].values()
        )
        != 24
    ):
        raise ArithmeticError("conic-core repair input changed")

    secant_line_rows = []
    for deleted_on_line in (1, 2, 3):
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
        raise ArithmeticError("secant slack incidence bound changed")

    return {
        "unique_profile": profile,
        "repair_deletion_bound": PAIR_SLACK // 4,
        "two_or_fewer_deletions": (
            "the repaired arc has size at least 18, hence is conic-contained; "
            "positive slack 12 contradicts the conic-core floor 24"
        ),
        "three_deletion_branch": {
            "arc_size": ARC_SIZE,
            "deleted_point_count": 3,
            "if_arc_incomplete": (
                "adjoin one point to make an 18-arc, hence obtain the same "
                "conic-core contradiction"
            ),
            "therefore_arc_complete": True,
        },
        "secant_line_slack_rows": secant_line_rows,
        "global_incidence_bound": "slack(S)>=4*sum_{x in D} mu_A(x)",
        "completeness_floor": "mu_A(x)>=1 for every x outside A",
        "equality_forced_by_slack_twelve": (
            "three deleted points and slack 12 force mu_A(x)=1 for all three"
        ),
        "required_one_secant_points": 3,
        "proved": True,
    }


def p23_slack_twelve_exclusion() -> dict[str, object]:
    """Proposition 15.685."""
    repair = slack_twelve_repair_certificate()
    classification = complete_17_arc_classification_certificate()
    required = int(repair["required_one_secant_points"])
    available = int(classification["maximum_one_secant_point_count"])
    if not available < required:
        raise ArithmeticError("complete-arc classification no longer excludes profile")

    previous = p23_reduction_theorem()
    before = int(previous["residue_zero_profile_count_after"])
    after = before - 1
    remaining_histogram = dict(previous["remaining_pair_slack_histogram"])
    remaining_histogram.pop(PAIR_SLACK)
    if before != 203 or after != 202 or sum(remaining_histogram.values()) != after:
        raise ArithmeticError("post-15.685 profile accounting changed")

    return {
        "proposition": "15.685",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "excluded_pair_slack": PAIR_SLACK,
        "required_one_secant_points": required,
        "maximum_available_in_any_complete_17_arc_class": available,
        "slack_twelve_profile_excluded": True,
        "p23_profile_count_before": before,
        "p23_profiles_excluded_here": 1,
        "p23_profile_count_after": after,
        "remaining_pair_slack_histogram": remaining_histogram,
        "p23_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "repair_certificate": repair,
        "classification_certificate": classification,
        "external_dependencies": [
            {
                "authors": "K. Coolsaet and H. Sticker",
                "title": (
                    "A full classification of the complete k-arcs of "
                    "PG(2,23) and PG(2,25)"
                ),
                "journal": "Journal of Combinatorial Designs 17 (2009), 459--477",
                "doi": "10.1002/jcd.20211",
                "input_used": (
                    "there are exactly five projective classes of complete "
                    "17-arcs in PG(2,23)"
                ),
            }
        ],
        "proved": True,
    }


def main() -> None:
    theorem = p23_slack_twelve_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15685.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.685: p=23 slack-12 profile excluded; "
        f"exact remainder {theorem['p23_profile_count_before']} -> "
        f"{theorem['p23_profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
