#!/usr/bin/env python3
r"""Prop. 15.727 -- rigidity and the first four closes at the slack endpoint.

Continue Proposition 15.726 with a ``p+1``-point affine set ``D`` and

    R = floor((p-1)/3).

Choose a *minimum-cardinality* set ``T`` such that ``A=D\T`` is an arc,
and write ``t=|T|``.  The deletion lemma gives ``t<=R``.  If ``t<R``, the
Ball--Lavrauw tangent envelope used in Proposition 15.726 still applies.
For ``I=sum_(z in T) s_A(z)`` it gives

    I >= t(p-1-3t)/2.

This concave quadratic is greater than ``R`` at both ends of
``1<=t<=R-1``, while the linewise slack inequality gives ``I<=R``.
Therefore ``t=R``.  Minimum-cardinality implies inclusion-minimality, so
every deleted point has positive integral secant index.  Hence

    I=R, and s_A(z)=1 for every z in T.

Equality in the linewise slack comparison is rigid.  Every line of ``D``
with at least three points is either a trisecant containing two points of
``A`` and one of ``T``, or a 4-secant containing two points of each.  No two
such rich lines share a point of ``D``: otherwise their deletion demands can
share that point, producing an arc after at most ``R-1`` deletions.  If
``x`` and ``y`` count the trisecants and 4-secants, respectively, then

    x+2y=R.

In particular ``T`` consists of ``R`` distinct outside points of secant
index one for the arc ``A``.  Thus ``c_1(A)>=R``.

Existing exhaustive planar-arc classifications now exclude the endpoint
for the first four in-scope primes.

* ``p=17``: ``A`` is a 13-arc and would need ``c_1>=5``.  A complete
  13-arc has at most three index-one outside points.  An incomplete 13-arc
  extends either into the unique complete 14-arc, whose 13-subarcs have at
  most four, or into a 15-arc and hence a conic, where it has none.
* ``p=19``: every classified 14-arc has ``c_1<=4``, versus the required six.
* ``p=23``: a complete 17-arc has ``c_1<=1``, versus the required seven.
  An incomplete 17-arc extends to an 18-arc and hence to the 24-point conic,
  where it has no index-one outside point.
* ``p=29``: the two complete 21-arcs both have ``c_1=0``.  An incomplete
  21-arc extends either to the unique complete 24-arc (the Klein quartic),
  whose outside secant index is at least six, or to the conic.  Deleting
  three Klein points or nine conic points still leaves no index-one point.

Consequently the first unexcluded slack is at least 6, 7, 8, and 10 at
``p=17,19,23,29``.  From ``p=31`` onward the endpoint is reduced to the
disjoint 3/4-secant normal form; it is not excluded here.  The ``p+1`` shell,
residual (ii), Type I, and the quadratic-minmax limit remain open.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

from e1_gmin_m4_prop15684 import p23_complete_arc_classification
from e1_gmin_m4_prop15685 import (
    Point,
    complete_17_arc_classification_certificate,
    incident,
    line_through,
    normalize_projective,
    projective_points,
)
from e1_gmin_m4_prop15693 import p19_fourteen_arc_secant_index_classification
from e1_gmin_m4_prop15701 import p17_fifteen_arc_classification
from e1_gmin_m4_prop15702 import p17_complete_fourteen_arc_classification
from e1_gmin_m4_prop15703 import (
    complete_fourteen_minus_one_certificate,
    complete_thirteen_arc_certificate,
)
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15722 import occupancy_slack_term
from e1_gmin_m4_prop15726 import tangent_envelope_dependency


ROOT = Path(__file__).resolve().parents[1]
CLASSIFIED_ENDPOINT_PRIMES = (17, 19, 23, 29)


# Sticker, thesis, Section 6.14, equations (6.9) and (6.10).  Coordinates
# are homogeneous over F_29; negative entries are normalized by the audit.
P29_COMPLETE_21_ARC_REPRESENTATIVES: tuple[tuple[Point, ...], ...] = (
    (
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 5, 10), (1, 10, 5), (5, 1, 10), (10, 1, 5),
        (5, 10, 1), (10, 5, 1),
        (1, 4, 9), (1, 9, 4), (4, 1, 9), (9, 1, 4),
        (4, 9, 1), (9, 4, 1),
        (1, -3, -2), (1, -2, -3), (-3, 1, -2), (-2, 1, -3),
        (-3, -2, 1), (-2, -3, 1),
    ),
    (
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 2, 8), (1, 8, 2), (2, 1, 8), (8, 1, 2),
        (2, 8, 1), (8, 2, 1),
        (1, 5, 13), (1, 13, 5), (5, 1, 13), (13, 1, 5),
        (5, 13, 1), (13, 5, 1),
        (1, -3, -5), (1, -5, -3), (-3, 1, -5), (-5, 1, -3),
        (-3, -5, 1), (-5, -3, 1),
    ),
)

P29_EXPECTED_COMPLETE_21_HISTOGRAMS: tuple[dict[int, int], ...] = (
    {4: 18, 5: 75, 6: 190, 7: 312, 8: 189, 9: 63, 10: 3},
    {3: 3, 4: 21, 5: 66, 6: 187, 7: 294, 8: 243, 9: 27, 10: 9},
)
P29_EXPECTED_KLEIN_24_HISTOGRAM = {
    6: 28,
    8: 126,
    9: 504,
    10: 84,
    11: 84,
    12: 21,
}


def _check_prime_parameter(p: int) -> None:
    if p < 17 or not is_prime(p):
        raise ValueError("need an odd prime parameter p>=17")


def first_unexcluded_endpoint(p: int) -> int:
    """The first integer left by Proposition 15.726."""
    _check_prime_parameter(p)
    return (p - 1) // 3


def endpoint_residue_data(p: int) -> dict[str, int]:
    """Write ``p=3R+c`` with ``c`` equal to one or two."""
    R = first_unexcluded_endpoint(p)
    c = p - 3 * R
    if c not in (1, 2):
        raise ArithmeticError("prime endpoint has an impossible mod-three residue")
    return {"p": p, "R": R, "c": c}


def proper_deletion_size_contradiction(p: int, t: int) -> dict[str, object]:
    """Audit one ``t<R`` case at the endpoint.

    The tangent-envelope lower bound is the same as in Proposition 15.726;
    the only new point is that its size hypothesis continues to hold through
    ``t=R-1``.
    """
    data = endpoint_residue_data(p)
    R = data["R"]
    if not 1 <= t < R:
        raise ValueError("need 1<=t<R at the first unexcluded endpoint")
    arc_size = p + 1 - t
    tau = t + 1
    size_threshold = 2 * tau + 2
    incidence_lower_twice = t * (p - 1 - 3 * t)
    contradicts_I_upper = incidence_lower_twice > 2 * R
    return {
        "p": p,
        "R": R,
        "t": t,
        "arc_size": arc_size,
        "tangent_deficiency_tau": tau,
        "envelope_degree": 2 * tau,
        "envelope_size_threshold": size_threshold,
        "size_hypothesis_met": arc_size >= size_threshold,
        "twice_incidence_lower_bound": incidence_lower_twice,
        "twice_slack_upper_bound": 2 * R,
        "incidence_lower_bound_exceeds_R": contradicts_I_upper,
        "contradiction": arc_size >= size_threshold and contradicts_I_upper,
        "proved": arc_size >= size_threshold and contradicts_I_upper,
    }


def endpoint_concavity_ledger(p: int) -> dict[str, object]:
    """Record the symbolic endpoint margins for every proper deletion size."""
    data = endpoint_residue_data(p)
    R = data["R"]
    c = data["c"]
    rows = [proper_deletion_size_contradiction(p, t) for t in range(1, R)]
    f1_twice_minus_2R = p - 4 - 2 * R
    f_last_twice_minus_2R = (R - 1) * (p - 1 - 3 * (R - 1)) - 2 * R
    expected = {
        1: (R - 3, R - 3),
        2: (R - 2, 2 * (R - 2)),
    }[c]
    if (f1_twice_minus_2R, f_last_twice_minus_2R) != expected:
        raise ArithmeticError("endpoint concavity margins changed")
    proved = bool(
        R >= 5
        and f1_twice_minus_2R > 0
        and f_last_twice_minus_2R > 0
        and all(row["proved"] for row in rows)
    )
    return {
        **data,
        "proper_deletion_interval": [1, R - 1],
        "envelope_size_margin_at_t_R_minus_one": c,
        "twice_F_1_minus_R": f1_twice_minus_2R,
        "twice_F_R_minus_one_minus_R": f_last_twice_minus_2R,
        "concavity": (
            "F(t)=t(p-1-3t)/2 is concave, so its minimum on "
            "1<=t<=R-1 is attained at an endpoint"
        ),
        "all_proper_deletion_sizes_excluded": all(row["proved"] for row in rows),
        "rows": rows,
        "proved": proved,
    }


def linewise_endpoint_equality(a: int, u: int) -> dict[str, object]:
    """Classify equality in the linewise comparison ``I<=R``."""
    if a not in (0, 1, 2) or u < 0:
        raise ValueError("need 0<=a<=2 and u>=0")
    occupancy = a + u
    slack = occupancy_slack_term(occupancy)
    incidence = u if a == 2 else 0
    equality = slack == incidence
    rich = occupancy >= 3
    allowed_rich_equality = a == 2 and u in (1, 2)
    return {
        "arc_points_a": a,
        "deleted_points_u": u,
        "line_occupancy": occupancy,
        "slack_contribution": slack,
        "secant_incidence_contribution": incidence,
        "equality": equality,
        "rich": rich,
        "rich_equality_case": equality and rich,
        "allowed_rich_equality_case": allowed_rich_equality,
        "classification_correct": (equality and rich) == allowed_rich_equality,
    }


def endpoint_block_row(p: int, four_secants: int) -> dict[str, object]:
    """One possible disjoint 3/4-secant count at endpoint equality."""
    data = endpoint_residue_data(p)
    R = data["R"]
    c = data["c"]
    y = four_secants
    if y < 0 or 2 * y > R:
        raise ValueError("invalid number of disjoint 4-secants")
    x = R - 2 * y
    points_on_rich_lines = 3 * x + 4 * y
    singleton_points = p + 1 - points_on_rich_lines
    maximum_arc_size = p + 1 - R
    line_counts = {
        0: p * (p - 1) // 2 - R - y,
        1: p + 1 + 3 * R + 2 * y,
        2: p * (p + 1) // 2 - 3 * R,
        3: x,
        4: y,
    }
    projective_line_count = p * p + p + 1
    point_line_incidence_count = (p + 1) * (p + 1)
    point_pair_count = p * (p + 1) // 2
    line_moments_match = bool(
        sum(line_counts.values()) == projective_line_count
        and sum(size * count for size, count in line_counts.items())
        == point_line_incidence_count
        and sum(
            size * (size - 1) // 2 * count
            for size, count in line_counts.items()
        )
        == point_pair_count
    )

    # Remove every singleton and one point from every 4-secant.  What remains
    # is a regular semiarc whose trisecants partition its points.
    core_size = 3 * (R - y)
    core_tangents_per_point = c + 3 + 3 * y
    core_identity = core_size + core_tangents_per_point == p + 3
    return {
        **data,
        "trisecants_x": x,
        "four_secants_y": y,
        "slack_check_x_plus_2y": x + 2 * y,
        "points_on_rich_lines": points_on_rich_lines,
        "singleton_points": singleton_points,
        "singleton_formula": c + 1 + 2 * y,
        "minimum_deletions": x + 2 * y,
        "maximum_arc_size": maximum_arc_size,
        "maximum_arc_choice_count": 3**x * 6**y,
        "rich_lines_pairwise_D_disjoint": True,
        "projective_line_occupancy_counts": line_counts,
        "line_moments_match": line_moments_match,
        "mod_three_parity": (
            "R is even" if c == 1 else "R is odd, so at least one trisecant remains"
        ),
        "regular_trisecant_core": {
            "construction": (
                "remove every singleton point and one point from each 4-secant"
            ),
            "point_count": core_size,
            "point_count_formula": "3(R-y)=p-c-3y",
            "trisecant_count": R - y,
            "every_point_on_exactly_one_trisecant": True,
            "tangents_per_point": core_tangents_per_point,
            "tangent_formula": "c+3+3y",
            "size_plus_tangents": p + 3,
        },
        "proved": bool(
            singleton_points == c + 1 + 2 * y
            and x + 2 * y == R
            and line_moments_match
            and core_size == p - c - 3 * y
            and core_identity
            and (c != 2 or x >= 1)
        ),
    }


def universal_endpoint_rigidity(p: int) -> dict[str, object]:
    """Package the all-prime endpoint equality reduction."""
    data = endpoint_residue_data(p)
    R = data["R"]
    concavity = endpoint_concavity_ledger(p)
    block_rows = [endpoint_block_row(p, y) for y in range(R // 2 + 1)]
    equality_rows = [
        linewise_endpoint_equality(a, u)
        for a in range(3)
        for u in range(0, 8)
    ]
    proved = bool(
        concavity["proved"]
        and all(row["classification_correct"] for row in equality_rows)
        and all(row["proved"] for row in block_rows)
    )
    return {
        **data,
        "minimum_cardinality_arc_deletion_size": R,
        "incidence_sum_I": R,
        "every_deleted_point_secant_index": 1,
        "required_c1_of_arc_at_least": R,
        "arc_size": p + 1 - R,
        "rich_line_types": [
            "trisecant: two arc points plus one deleted point",
            "4-secant: two arc points plus two deleted points",
        ],
        "rich_lines_pairwise_D_disjoint": True,
        "disjointness_reason": (
            "if two rich lines shared a D-point, choose that point in both "
            "line-deletion demands and repair every rich line with at most R-1 "
            "distinct deletions, contradicting minimum cardinality R"
        ),
        "block_equation": "x+2y=R",
        "block_rows": block_rows,
        "proper_deletion_contradiction": concavity,
        "linewise_equality_audit": equality_rows,
        "tangent_envelope_dependency": tangent_envelope_dependency(),
        "finite_search_used_for_universal_reduction": False,
        "endpoint_excluded_by_universal_reduction_alone": False,
        "result_status": "proved structural reduction",
        "proved": proved,
    }


def _conic_c1_exclusion(p: int, arc_size: int) -> dict[str, object]:
    """Show that a conic-contained arc has no index-one outside point."""
    omitted = p + 1 - arc_size
    full_secants_min = (p - 1) // 2
    retained_secants_min = full_secants_min - omitted
    proved = retained_secants_min >= 2
    return {
        "p": p,
        "arc_size": arc_size,
        "conic_size": p + 1,
        "omitted_conic_points": omitted,
        "off_conic_full_secants_at_least": full_secants_min,
        "off_conic_retained_secants_at_least": retained_secants_min,
        "missing_conic_point_secant_index": 0,
        "conclusion": "a conic-contained arc has c1=0",
        "proved": proved,
    }


def _arc_secant_index_certificate(
    points: Iterable[Point], p: int, expected_size: int
) -> dict[str, object]:
    """Verify an arc and its full outside secant-index histogram."""
    arc = tuple(sorted({normalize_projective(point, p) for point in points}))
    if len(arc) != expected_size:
        raise ArithmeticError("arc representative has the wrong point count")
    secants = Counter(
        line_through(arc[first], arc[second], p)
        for first in range(len(arc))
        for second in range(first)
    )
    is_arc = bool(
        len(secants) == math.comb(expected_size, 2)
        and max(secants.values(), default=0) == 1
    )
    if not is_arc:
        raise ArithmeticError("arc representative contains three collinear points")
    arc_set = set(arc)
    outside_indices = {
        point: sum(incident(point, line, p) for line in secants)
        for point in projective_points(p)
        if point not in arc_set
    }
    histogram = dict(sorted(Counter(outside_indices.values()).items()))
    expected_outside = p * p + p + 1 - expected_size
    expected_incidence = math.comb(expected_size, 2) * (p - 1)
    if (
        len(outside_indices) != expected_outside
        or sum(histogram.values()) != expected_outside
        or sum(index * count for index, count in histogram.items())
        != expected_incidence
    ):
        raise ArithmeticError("outside secant incidence accounting changed")
    return {
        "p": p,
        "coordinates": [list(point) for point in arc],
        "point_count": expected_size,
        "secant_line_count": len(secants),
        "outside_point_count": len(outside_indices),
        "outside_secant_index_histogram": histogram,
        "index_one_point_count": histogram.get(1, 0),
        "minimum_outside_secant_index": min(outside_indices.values()),
        "maximum_outside_secant_index": max(outside_indices.values()),
        "outside_secant_incidence_count": sum(outside_indices.values()),
        "is_arc": is_arc,
        "is_complete": min(outside_indices.values()) >= 1,
        "proved": True,
    }


def p29_complete_arc_spectrum() -> dict[str, object]:
    """The published complete-arc spectrum needed above size 21."""
    return {
        "p": 29,
        "source": (
            "K. Coolsaet and H. Sticker, The complete k-arcs of PG(2,27) "
            "and PG(2,29), J. Combin. Des. 19 (2011), 111--130"
        ),
        "doi": "10.1002/jcd.20261",
        "locations": [
            "Section 5 and Table 3, article pages 15--16 (complete spectrum)",
            "Sections 5.6--5.7, article page 20 (representatives and Klein quartic)",
        ],
        "complete_arc_counts_used": {21: 2, 24: 1, 30: 1},
        "no_complete_arc_sizes": [22, 23, 25, 26, 27, 28, 29],
        "maximum_arc_size": 30,
        "maximum_size_reason": (
            "Segre's odd-order q+1 arc bound; the classified size-30 arc "
            "attains it"
        ),
        "size_24_identification": "Klein quartic x^3 y+y^3 z+z^3 x=0",
        "size_30_identification": "nondegenerate conic",
        "classification_external_dependency": True,
        "proved_conditional_on_external_classification": True,
    }


def p29_complete_twenty_one_arc_certificate() -> dict[str, object]:
    """Audit both classified complete 21-arcs in ``PG(2,29)``."""
    spectrum = p29_complete_arc_spectrum()
    rows = [
        _arc_secant_index_certificate(points, 29, 21)
        for points in P29_COMPLETE_21_ARC_REPRESENTATIVES
    ]
    observed = [row["outside_secant_index_histogram"] for row in rows]
    expected = list(P29_EXPECTED_COMPLETE_21_HISTOGRAMS)
    distinct_invariants = len(
        {tuple(sorted(histogram.items())) for histogram in observed}
    ) == len(rows)
    proved = bool(
        int(spectrum["complete_arc_counts_used"][21]) == 2
        and len(rows) == 2
        and observed == expected
        and distinct_invariants
        and all(row["is_arc"] and row["is_complete"] for row in rows)
        and all(int(row["index_one_point_count"]) == 0 for row in rows)
    )
    return {
        "p": 29,
        "arc_size": 21,
        "classified_projective_class_count": 2,
        "verified_representative_count": len(rows),
        "pairwise_distinct_secant_index_histograms": distinct_invariants,
        "therefore_exhaustive": len(rows)
        == int(spectrum["complete_arc_counts_used"][21]),
        "index_one_point_counts_by_class": [
            int(row["index_one_point_count"]) for row in rows
        ],
        "maximum_index_one_point_count": max(
            int(row["index_one_point_count"]) for row in rows
        ),
        "representatives": rows,
        "classification_source": spectrum,
        "proved": proved,
    }


def p29_klein_twenty_four_arc_certificate() -> dict[str, object]:
    """Audit the unique complete 24-arc and its outside secant indices."""
    points = tuple(
        point
        for point in projective_points(29)
        if (
            point[0] ** 3 * point[1]
            + point[1] ** 3 * point[2]
            + point[2] ** 3 * point[0]
        )
        % 29
        == 0
    )
    row = _arc_secant_index_certificate(points, 29, 24)
    histogram = row["outside_secant_index_histogram"]
    proved = bool(
        histogram == P29_EXPECTED_KLEIN_24_HISTOGRAM
        and int(row["minimum_outside_secant_index"]) == 6
        and row["is_arc"]
        and row["is_complete"]
    )
    return {
        **row,
        "equation": "x^3 y+y^3 z+z^3 x=0 over F_29",
        "classified_projective_class_count": 1,
        "classification_source": p29_complete_arc_spectrum(),
        "proved": proved,
    }


def p17_endpoint_exclusion() -> dict[str, object]:
    """Exclude ``R=5`` by the audited 13/14/15-arc classifications."""
    rigidity = universal_endpoint_rigidity(17)
    complete_13 = complete_thirteen_arc_certificate()
    complete_14 = p17_complete_fourteen_arc_classification()
    fourteen_minus_one = complete_fourteen_minus_one_certificate()
    fifteen = p17_fifteen_arc_classification()
    conic = _conic_c1_exclusion(17, 13)
    required = int(rigidity["required_c1_of_arc_at_least"])
    complete_13_max = max(int(value) for value in complete_13["index_one_point_counts"])
    complete_14_subarc_max = max(
        int(value) for value in fourteen_minus_one["index_one_count_histogram"]
    )
    available = max(complete_13_max, complete_14_subarc_max, 0)
    proved = bool(
        rigidity["proved"]
        and required == 5
        and complete_13_max == 3
        and complete_13["proved_conditional_on_published_class_count"]
        and int(complete_14["complete_14_arc_class_count"]) == 1
        and complete_14["proved_conditional_on_external_classification"]
        and complete_14_subarc_max == 4
        and fourteen_minus_one["proved"]
        and int(fifteen["pgl_class_count_of_15_arcs"]) == 1
        and fifteen["proved_conditional_on_external_classification"]
        and conic["proved"]
        and available < required
    )
    return {
        "p": 17,
        "R": 5,
        "arc_size": 13,
        "required_c1": required,
        "complete_13_arc_c1_maximum": complete_13_max,
        "complete_14_minus_one_c1_maximum": complete_14_subarc_max,
        "conic_contained_c1": 0,
        "maximum_available_c1": available,
        "incomplete_branch": (
            "extend to a 14-arc; a complete extension is the unique audited "
            "class, while an incomplete extension reaches the unique conic-"
            "contained 15-arc class"
        ),
        "complete_13_certificate": complete_13,
        "complete_14_classification": complete_14,
        "complete_14_minus_one_certificate": fourteen_minus_one,
        "fifteen_arc_classification": fifteen,
        "conic_branch": conic,
        "excluded": proved,
        "proved": proved,
    }


def p19_endpoint_exclusion() -> dict[str, object]:
    """Exclude ``R=6`` with the exhaustive 14-arc secant-index table."""
    rigidity = universal_endpoint_rigidity(19)
    classification = p19_fourteen_arc_secant_index_classification()
    required = int(rigidity["required_c1_of_arc_at_least"])
    available = int(classification["maximum_c1_over_all_fourteen_arcs"])
    proved = bool(
        rigidity["proved"]
        and classification["proved_conditional_on_external_classification"]
        and required == 6
        and available == 4
        and available < required
    )
    return {
        "p": 19,
        "R": 6,
        "arc_size": 14,
        "required_c1": required,
        "maximum_available_c1": available,
        "classification": classification,
        "excluded": proved,
        "proved": proved,
    }


def p23_endpoint_exclusion() -> dict[str, object]:
    """Exclude ``R=7`` by the complete-17 classification and extension gap."""
    rigidity = universal_endpoint_rigidity(23)
    complete_17 = complete_17_arc_classification_certificate()
    spectrum = p23_complete_arc_classification()
    conic = _conic_c1_exclusion(23, 17)
    required = int(rigidity["required_c1_of_arc_at_least"])
    available = int(complete_17["maximum_one_secant_point_count"])
    proved = bool(
        rigidity["proved"]
        and required == 7
        and available == 1
        and complete_17["proved"]
        and spectrum["proved_conditional_on_external_classification"]
        and spectrum["no_complete_arc_sizes"] == [18, 19, 20, 21, 22, 23]
        and int(spectrum["complete_arc_counts"][24]) == 1
        and conic["proved"]
        and available < required
    )
    return {
        "p": 23,
        "R": 7,
        "arc_size": 17,
        "required_c1": required,
        "complete_17_arc_c1_maximum": available,
        "conic_contained_c1": 0,
        "incomplete_branch": (
            "extend to size 18; the complete-size gap 18..23 forces extension "
            "to the unique 24-point conic"
        ),
        "complete_17_certificate": complete_17,
        "complete_arc_spectrum": spectrum,
        "conic_branch": conic,
        "excluded": proved,
        "proved": proved,
    }


def p29_endpoint_exclusion() -> dict[str, object]:
    """Exclude ``R=9`` using the exhaustive complete-arc spectrum."""
    rigidity = universal_endpoint_rigidity(29)
    spectrum = p29_complete_arc_spectrum()
    complete_21 = p29_complete_twenty_one_arc_certificate()
    klein_24 = p29_klein_twenty_four_arc_certificate()
    conic = _conic_c1_exclusion(29, 21)
    required = int(rigidity["required_c1_of_arc_at_least"])
    klein_retained_minimum = int(klein_24["minimum_outside_secant_index"]) - 3
    proved = bool(
        rigidity["proved"]
        and required == 9
        and complete_21["proved"]
        and spectrum["proved_conditional_on_external_classification"]
        and int(complete_21["maximum_index_one_point_count"]) == 0
        and spectrum["no_complete_arc_sizes"]
        == [22, 23, 25, 26, 27, 28, 29]
        and int(spectrum["complete_arc_counts_used"][24]) == 1
        and int(spectrum["complete_arc_counts_used"][30]) == 1
        and int(spectrum["maximum_arc_size"]) == 30
        and klein_24["proved"]
        and klein_retained_minimum == 3
        and conic["proved"]
    )
    return {
        "p": 29,
        "R": 9,
        "arc_size": 21,
        "required_c1": required,
        "complete_21_arc_c1_maximum": int(
            complete_21["maximum_index_one_point_count"]
        ),
        "klein_24_minus_three_outside_secant_index_minimum": (
            klein_retained_minimum
        ),
        "conic_contained_c1": 0,
        "incomplete_branch": (
            "extend stepwise: the complete-size gaps force either the unique "
            "24-point Klein arc or the unique 30-point conic; a 21-subarc "
            "has outside secant index zero on omitted arc points and at least "
            "three or five, respectively, off the containing arc"
        ),
        "complete_21_certificate": complete_21,
        "klein_24_certificate": klein_24,
        "complete_arc_spectrum": spectrum,
        "conic_branch": conic,
        "excluded": proved,
        "proved": proved,
    }


def classified_endpoint_exclusions() -> dict[int, dict[str, object]]:
    """Return all four classification-assisted endpoint closes."""
    rows = {
        17: p17_endpoint_exclusion(),
        19: p19_endpoint_exclusion(),
        23: p23_endpoint_exclusion(),
        29: p29_endpoint_exclusion(),
    }
    if not all(row["proved"] for row in rows.values()):
        raise ArithmeticError("a classified endpoint exclusion failed")
    return rows


def active_first_possible_slack(p: int) -> int:
    """Current first possible positive slack after Propositions 15.726--.727."""
    R = first_unexcluded_endpoint(p)
    return R + 1 if p in CLASSIFIED_ENDPOINT_PRIMES else R


def proposition_15727() -> dict[str, object]:
    """Package the universal reduction and its four finite-prime closes."""
    classified = classified_endpoint_exclusions()
    sample_primes = (17, 19, 23, 29, 31, 37, 41)
    universal = {str(p): universal_endpoint_rigidity(p) for p in sample_primes}
    proved = bool(
        all(row["proved"] for row in universal.values())
        and all(row["proved"] for row in classified.values())
    )
    return {
        "prop": "15.727",
        "statement": (
            "at R=floor((p-1)/3), every minimum arc repair has size R, "
            "all deleted points have secant index one, and the rich lines "
            "are point-disjoint trisecants/4-secants; the endpoint is "
            "excluded at p=17,19,23,29"
        ),
        "universal_sample_ledgers": universal,
        "classified_endpoint_exclusions": {
            str(p): row for p, row in classified.items()
        },
        "first_possible_positive_slack_after": {
            str(p): active_first_possible_slack(p) for p in sample_primes
        },
        "first_prime_not_endpoint_excluded_here": 31,
        "new_long_solver_run_used": False,
        "external_classification_assisted": True,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "the disjoint 3/4-secant endpoint normal form from p=31 onward, "
            "larger outside slack, the rest of residual (ii), Type I, and L"
        ),
        "result_status": "proved theorem with published classification inputs",
        "proved": proved,
    }


def write_evidence() -> Path:
    output = ROOT / "evidence" / "e1_gmin_m4_prop15727.json"
    output.write_text(json.dumps(proposition_15727(), indent=2, sort_keys=True) + "\n")
    return output


def main() -> None:
    result = proposition_15727()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.727 endpoint audit failed")
    path = write_evidence()
    print("Prop 15.727 endpoint rigidity: p=17,19,23,29 excluded")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
