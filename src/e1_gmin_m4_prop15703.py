#!/usr/bin/env python3
"""Prop. 15.703 -- close the p=17 slack-twelve block.

Proposition 15.702 leaves 113 slack-twelve profiles, all without an
undetermined direction.  Any realization repairs by exactly three deletions
to a 13-arc ``A``; shallower repair was excluded by the complete-14/conic
argument.

Sticker's exhaustive table has eight complete 13-arc classes in
``PG(2,17)``.  A normalized PGL generator fixes a quadrangle, finds a complete
arc, blocks every normalized image of its projective orbit, and repeats.  The
eight generated representatives have stabilizer orders
``1,2,2,2,2,3,4,6``, exactly the published fingerprint, and eight distinct
outside secant-index histograms.  Six classes have no index-one point, one
has two, and one has three.  The sole possible index-one triple in the last
class produces pair slack sixteen, not twelve.

If the repaired 13-arc is incomplete, extend it to a 14-arc.  An incomplete
14-arc extends to a conic-contained 15-arc and was excluded by 15.701.  The
unique complete 14-arc is handled by deleting each of its fourteen points.
The resulting 13-arcs have index-one counts ``{0:4,1:8,4:2}``; their eight
candidate triples all produce slack twenty.  If the deleted set contains the
extension point itself, the other two points each have complete-14-arc
secant index at least two and already force slack at least sixteen.

Thus no slack-twelve boundary exists.  All 113 rows are excluded and the p17
remainder falls from 1,481 to 1,368: two slack-zero tangent-conic rows and
1,366 rows of slack at least sixteen.  The endpoint remains open.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path

from e1_gmin_m4_prop15700 import p17_second_boundary_profile_census
from e1_gmin_m4_prop15702 import COMPLETE_14_ARC, p17_complete_fourteen_arc_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 17
BOUNDARY_SIZE = 16
Point = tuple[int, int, int]
CLASS_EVIDENCE = ROOT / "evidence" / "e1_gmin_m4_prop15703_complete13_classes.json"


def _projective_points_or_lines() -> tuple[Point, ...]:
    return tuple(
        [(1, y, z) for y in range(P) for z in range(P)]
        + [(0, 1, z) for z in range(P)]
        + [(0, 0, 1)]
    )


def _incident(point: Point, line: Point) -> bool:
    return sum(a * b for a, b in zip(point, line)) % P == 0


def _line_pair_slack(occupancy: int) -> int:
    return 2 * (math.comb(occupancy, 2) - occupancy // 2)


def _arc_data(arc: set[Point]) -> dict[str, object]:
    points = _projective_points_or_lines()
    occupancies = [
        sum(point in arc for point in points if _incident(point, line))
        for line in points
    ]
    if max(occupancies) > 2:
        raise ArithmeticError("purported p17 arc has a 3-secant")
    secants = [line for line, count in zip(points, occupancies) if count == 2]
    indices = {
        point: sum(_incident(point, line) for line in secants)
        for point in points
        if point not in arc
    }
    return {
        "line_occupancy_histogram": dict(sorted(Counter(occupancies).items())),
        "outside_indices": indices,
        "outside_secant_index_histogram": dict(sorted(Counter(indices.values()).items())),
    }


def _set_pair_slack(points_set: set[Point]) -> int:
    points = _projective_points_or_lines()
    return sum(
        _line_pair_slack(
            sum(point in points_set for point in points if _incident(point, line))
        )
        for line in points
    )


def p17_complete_thirteen_arc_classification() -> dict[str, object]:
    """Published class count and stabilizer-order fingerprint."""
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
        "complete_13_arc_class_count": 8,
        "stabilizer_group_orders": [1, 2, 2, 2, 2, 3, 4, 6],
        "reported_groups": ["1", "C2", "C2", "C2", "C2", "C3", "C4", "S3"],
        "proved_conditional_on_external_classification": True,
    }


def complete_thirteen_arc_certificate() -> dict[str, object]:
    """Audit all eight representatives and their index-one triples."""
    payload = json.loads(CLASS_EVIDENCE.read_text())
    classification = p17_complete_thirteen_arc_classification()
    representatives = payload["representatives"]
    if (
        int(payload["p"]) != P
        or int(payload["arc_size"]) != 13
        or int(payload["class_count"]) != 8
        or len(representatives) != 8
        or payload["final_status"] != "PUBLISHED_CLASS_COUNT_REACHED"
    ):
        raise ArithmeticError("p17 complete-13 class evidence changed")

    rows = []
    histograms = []
    candidate_slacks = []
    for class_index, representative in enumerate(representatives):
        arc = {tuple(point) for point in representative["representative"]}
        if len(arc) != 13:
            raise ArithmeticError("p17 complete-13 representative size changed")
        data = _arc_data(arc)
        histogram = data["outside_secant_index_histogram"]
        expected = {int(key): int(value) for key, value in representative["outside_secant_index_histogram"].items()}
        if (
            data["line_occupancy_histogram"] != {0: 151, 1: 78, 2: 78}
            or histogram != expected
            or 0 in histogram
        ):
            raise ArithmeticError("p17 complete-13 representative audit changed")
        histograms.append(tuple(sorted(histogram.items())))
        c1 = [point for point, index in data["outside_indices"].items() if index == 1]
        triple_slacks = [
            _set_pair_slack(arc | set(triple)) for triple in combinations(c1, 3)
        ]
        candidate_slacks.extend(triple_slacks)
        rows.append(
            {
                "class_index": class_index,
                "inferred_pgl_stabilizer_order": int(
                    representative["inferred_pgl_stabilizer_order"]
                ),
                "outside_secant_index_histogram": histogram,
                "index_one_point_count": len(c1),
                "candidate_index_one_triple_slacks": triple_slacks,
            }
        )
    observed_orders = sorted(row["inferred_pgl_stabilizer_order"] for row in rows)
    c1_counts = sorted(row["index_one_point_count"] for row in rows)
    if (
        len(set(histograms)) != 8
        or observed_orders != classification["stabilizer_group_orders"]
        or c1_counts != [0, 0, 0, 0, 0, 0, 2, 3]
        or candidate_slacks != [16]
    ):
        raise ArithmeticError("p17 complete-13 class fingerprint changed")
    return {
        "class_count": len(rows),
        "pairwise_invariant_distinct_secant_index_histograms": True,
        "observed_stabilizer_group_orders": observed_orders,
        "index_one_point_counts": c1_counts,
        "all_candidate_index_one_triple_slacks": candidate_slacks,
        "class_rows": rows,
        "generator_evidence": str(CLASS_EVIDENCE.relative_to(ROOT)),
        "generator_method": payload["classification_basis"],
        "proved_conditional_on_published_class_count": True,
    }


def complete_fourteen_minus_one_certificate() -> dict[str, object]:
    """Audit every 13-subarc of the unique complete 14-arc."""
    complete = set(COMPLETE_14_ARC)
    c1_count_histogram = Counter()
    candidate_slacks = []
    for removed in COMPLETE_14_ARC:
        arc = complete - {removed}
        data = _arc_data(arc)
        c1 = [point for point, index in data["outside_indices"].items() if index == 1]
        c1_count_histogram[len(c1)] += 1
        candidate_slacks.extend(
            _set_pair_slack(arc | set(triple)) for triple in combinations(c1, 3)
        )
    if dict(sorted(c1_count_histogram.items())) != {0: 4, 1: 8, 4: 2}:
        raise ArithmeticError("p17 complete-14-minus-one c1 census changed")
    if candidate_slacks != [20] * 8:
        raise ArithmeticError("p17 complete-14-minus-one triple census changed")
    return {
        "deleted_point_cases": len(COMPLETE_14_ARC),
        "index_one_count_histogram": dict(sorted(c1_count_histogram.items())),
        "candidate_index_one_triple_count": len(candidate_slacks),
        "candidate_index_one_triple_slack_histogram": dict(
            sorted(Counter(candidate_slacks).items())
        ),
        "extension_point_in_deleted_set_case": (
            "the other two deleted points are outside the complete 14-arc and "
            "each have secant index at least two, forcing slack at least 16"
        ),
        "proved": True,
    }


def p17_slack_twelve_exclusion() -> dict[str, object]:
    """Proposition 15.703."""
    previous = p17_complete_fourteen_arc_exclusion()
    census = p17_second_boundary_profile_census()
    classification = p17_complete_thirteen_arc_classification()
    complete_certificate = complete_thirteen_arc_certificate()
    incomplete_certificate = complete_fourteen_minus_one_certificate()
    profiles = [row for row in census["profiles"] if int(row["pair_slack"]) == 12]
    t_histogram = Counter(
        sum(
            int(row["phase_profiles_b"][phase].get(BOUNDARY_SIZE, 0))
            for phase in ("0", "1")
        )
        for row in profiles
    )
    if dict(sorted(t_histogram.items())) != {0: 113, 1: 111, 2: 43}:
        raise ArithmeticError("p17 slack-twelve profile ledger changed")

    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    excluded_indices = {
        index
        for index in previous_indices
        if int(census["profiles"][index]["pair_slack"]) == 12
    }
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
        before != 1481
        or excluded != 113
        or after != 1368
        or 12 in histogram
        or sum(histogram.values()) != after
    ):
        raise ArithmeticError("p17 post-slack-twelve accounting changed")
    return {
        "proposition": "15.703",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(excluded_indices),
        "remaining_profile_indices": remaining_indices,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 2,
        "remaining_profiles_of_slack_at_least_sixteen": after - 2,
        "slack_twelve_profile_undetermined_histogram_before_15701_to_15703": dict(
            sorted(t_histogram.items())
        ),
        "repair_normal_form": (
            "exactly three deletions to a 13-arc; a complete core is one of "
            "eight classes, while an incomplete core extends either to a "
            "conic branch or to a subarc of the unique complete 14-arc"
        ),
        "complete_13_arc_certificate": complete_certificate,
        "complete_14_minus_one_certificate": incomplete_certificate,
        "external_classification": classification,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p17_slack_twelve_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15703.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.703: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
