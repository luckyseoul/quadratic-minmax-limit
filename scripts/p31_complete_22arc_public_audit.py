#!/usr/bin/env python3
"""Audit the eleven publicly listed complete 22-arcs in PG(2,31).

Gerzson Keri's public ``31q3x19.txt`` A-list predates Coolsaet's corrected
classification: it contains eleven representatives, while the corrected
projective class count is twelve.  This script exhaustively checks the
eleven rows that are actually public.  It deliberately does not synthesize,
assume, or certify the missing twelfth representative.

For a prime field, Keri's compact A-list records a 2 by 18 matrix.  Following
his published instructions, its two rows ``a_i,b_i`` expand to the 22 points

    e_1, e_2, e_3, (1,1,1), (1,a_i,b_i)  (1 <= i <= 18).

The audit verifies the arc and completeness conditions, computes every
outside secant index, and computes the largest matching among distinct
secants supporting index-one outside points.  Output is deterministic JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Iterable


P = 31
ARC_SIZE = 22
PUBLIC_CLASS_COUNT = 11
CORRECTED_CLASS_COUNT = 12
ENDPOINT_REQUIRED_C1 = 10
EXPECTED_C1_SEQUENCE = (0, 2, 0, 0, 0, 0, 0, 2, 2, 0, 0)
EXPECTED_MATCHING_SEQUENCE = (0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "evidence" / "data" / "31q3x19.txt"
DEFAULT_OUTPUT = ROOT / "evidence" / "p31_complete_22arc_public_11_audit.json"

Point = tuple[int, int, int]
Line = tuple[int, int, int]


def normalize_projective(vector: Iterable[int]) -> Point:
    """Normalize a nonzero vector over F_31 by its first nonzero entry."""
    values = tuple(value % P for value in vector)
    for value in values:
        if value:
            inverse = pow(value, -1, P)
            return tuple(entry * inverse % P for entry in values)  # type: ignore[return-value]
    raise ValueError("the zero vector has no projective normalization")


def projective_points() -> tuple[Point, ...]:
    """The 993 normalized points of PG(2,31)."""
    return tuple(
        [(1, y, z) for y in range(P) for z in range(P)]
        + [(0, 1, z) for z in range(P)]
        + [(0, 0, 1)]
    )


def line_through(first: Point, second: Point) -> Line:
    return normalize_projective(
        (
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        )
    )


def incident(point: Point, line: Line) -> bool:
    return sum(a * b for a, b in zip(point, line)) % P == 0


def parse_public_a_list(path: Path = DEFAULT_INPUT) -> tuple[tuple[Point, ...], ...]:
    """Parse Keri's eleven compact 2-by-18 A-list representatives."""
    text = path.read_text(encoding="ascii")
    if not text.startswith("A-list:"):
        raise ValueError("expected Keri A-list header")
    compact_rows: list[tuple[int, ...]] = []
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if len(fields) == 18 and all(field.isdigit() for field in fields):
            row = tuple(int(field) for field in fields)
            if any(not 0 <= value < P for value in row):
                raise ValueError("A-list entry is not a residue modulo 31")
            compact_rows.append(row)
    if len(compact_rows) != 2 * PUBLIC_CLASS_COUNT:
        raise ValueError(
            f"expected {2 * PUBLIC_CLASS_COUNT} compact rows, got {len(compact_rows)}"
        )

    arcs: list[tuple[Point, ...]] = []
    frame: tuple[Point, ...] = (
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 1),
    )
    for offset in range(0, len(compact_rows), 2):
        top, bottom = compact_rows[offset : offset + 2]
        arc = frame + tuple((1, a, b) for a, b in zip(top, bottom))
        arcs.append(tuple(normalize_projective(point) for point in arc))
    return tuple(arcs)


def maximum_matching(
    vertices: tuple[Point, ...], edges: Iterable[tuple[Point, Point]]
) -> tuple[int, list[list[list[int]]]]:
    """Return an exact maximum matching and one deterministic witness."""
    index = {point: position for position, point in enumerate(vertices)}
    adjacency = [0] * len(vertices)
    for first, second in edges:
        i, j = index[first], index[second]
        adjacency[i] |= 1 << j
        adjacency[j] |= 1 << i

    @lru_cache(maxsize=None)
    def solve(mask: int) -> tuple[int, tuple[tuple[int, int], ...]]:
        if not mask:
            return 0, ()
        first_bit = mask & -mask
        i = first_bit.bit_length() - 1
        best_size, best_matching = solve(mask ^ first_bit)
        choices = adjacency[i] & mask
        while choices:
            second_bit = choices & -choices
            j = second_bit.bit_length() - 1
            size, matching = solve(mask ^ first_bit ^ second_bit)
            candidate = tuple(sorted(matching + ((i, j),)))
            if size + 1 > best_size or (
                size + 1 == best_size and candidate < best_matching
            ):
                best_size = size + 1
                best_matching = candidate
            choices ^= second_bit
        return best_size, best_matching

    size, matching = solve((1 << len(vertices)) - 1)
    witness = [
        [list(vertices[first]), list(vertices[second])]
        for first, second in matching
    ]
    return size, witness


def _string_keyed(counter: Counter[int]) -> dict[str, int]:
    return {str(key): counter[key] for key in sorted(counter)}


def audit_representative(arc: tuple[Point, ...], class_index: int) -> dict[str, object]:
    """Compute an exact outside-secant certificate for one representative."""
    if len(arc) != ARC_SIZE or len(set(arc)) != ARC_SIZE:
        raise ArithmeticError("representative does not contain 22 distinct points")

    secant_pairs: dict[Line, tuple[Point, Point]] = {}
    for first in range(len(arc)):
        for second in range(first):
            line = line_through(arc[first], arc[second])
            if line in secant_pairs:
                raise ArithmeticError("representative contains three collinear points")
            secant_pairs[line] = tuple(sorted((arc[first], arc[second])))  # type: ignore[assignment]
    if len(secant_pairs) != math.comb(ARC_SIZE, 2):
        raise ArithmeticError("secant count changed")

    arc_set = set(arc)
    outside_indices: dict[Point, int] = {}
    unique_line_by_point: dict[Point, Line] = {}
    for point in projective_points():
        if point in arc_set:
            continue
        incident_secants = tuple(
            line for line in secant_pairs if incident(point, line)
        )
        outside_indices[point] = len(incident_secants)
        if len(incident_secants) == 1:
            unique_line_by_point[point] = incident_secants[0]

    histogram = Counter(outside_indices.values())
    expected_outside = P * P + P + 1 - ARC_SIZE
    expected_incidence = math.comb(ARC_SIZE, 2) * (P - 1)
    if len(outside_indices) != expected_outside:
        raise ArithmeticError("outside point count changed")
    if sum(histogram.values()) != expected_outside:
        raise ArithmeticError("outside histogram count changed")
    if sum(index * count for index, count in histogram.items()) != expected_incidence:
        raise ArithmeticError("outside secant incidence moment changed")
    if min(outside_indices.values()) < 1:
        raise ArithmeticError("the publicly listed representative is not complete")

    points_by_unique_line: dict[Line, list[Point]] = defaultdict(list)
    for point, line in unique_line_by_point.items():
        points_by_unique_line[line].append(point)
    matching_size, matching_witness = maximum_matching(
        tuple(sorted(arc)),
        (secant_pairs[line] for line in sorted(points_by_unique_line)),
    )
    degree = Counter({point: 0 for point in arc})
    for line in points_by_unique_line:
        first, second = secant_pairs[line]
        degree[first] += 1
        degree[second] += 1

    unique_secants = [
        {
            "line": list(line),
            "arc_pair": [list(point) for point in secant_pairs[line]],
            "index_one_outside_points": [
                list(point) for point in sorted(points_by_unique_line[line])
            ],
        }
        for line in sorted(points_by_unique_line)
    ]
    c1 = len(unique_line_by_point)
    return {
        "public_class_index": class_index,
        "coordinates": [list(point) for point in arc],
        "point_count": len(arc),
        "secant_line_count": len(secant_pairs),
        "outside_point_count": len(outside_indices),
        "outside_secant_incidence_count": sum(outside_indices.values()),
        "outside_secant_index_histogram": _string_keyed(histogram),
        "minimum_outside_secant_index": min(outside_indices.values()),
        "maximum_outside_secant_index": max(outside_indices.values()),
        "index_one_point_count": c1,
        "unique_secant_line_count": len(points_by_unique_line),
        "index_one_points_per_unique_secant_line_histogram": _string_keyed(
            Counter(len(points) for points in points_by_unique_line.values())
        ),
        "unique_secant_graph_degree_histogram": _string_keyed(
            Counter(degree.values())
        ),
        "maximum_disjoint_unique_secant_matching": matching_size,
        "one_maximum_matching_arc_pairs": matching_witness,
        "index_one_points_share_one_secant": (
            len(points_by_unique_line) == 1 if c1 else None
        ),
        "unique_secants": unique_secants,
        "is_arc": True,
        "is_complete": True,
    }


def build_certificate(input_path: Path = DEFAULT_INPUT) -> dict[str, object]:
    """Build the finite certificate for exactly the eleven public classes."""
    representatives = [
        audit_representative(arc, class_index)
        for class_index, arc in enumerate(parse_public_a_list(input_path), start=1)
    ]
    c1_sequence = tuple(
        int(row["index_one_point_count"]) for row in representatives
    )
    matching_sequence = tuple(
        int(row["maximum_disjoint_unique_secant_matching"])
        for row in representatives
    )
    if c1_sequence != EXPECTED_C1_SEQUENCE:
        raise ArithmeticError("public-class c1 sequence changed")
    if matching_sequence != EXPECTED_MATCHING_SEQUENCE:
        raise ArithmeticError("public-class unique-secant matching sequence changed")
    positive_c1_classes = [
        int(row["public_class_index"])
        for row in representatives
        if int(row["index_one_point_count"]) > 0
    ]
    same_secant_classes = [
        int(row["public_class_index"])
        for row in representatives
        if row["index_one_points_share_one_secant"] is True
    ]
    return {
        "experiment": "p31_complete_22arc_public_audit",
        "result_status": "exhaustive finite certificate",
        "scope": (
            "all eleven complete 22-arc representatives in Keri's public "
            "31q3x19 A-list; the corrected twelfth projective class is absent"
        ),
        "sources": {
            "public_coordinate_list": (
                "https://old.sztaki.hu/~keri/n-arcs/31q3x19.txt"
            ),
            "coordinate_format_instructions": (
                "https://old.sztaki.hu/~keri/n-arcs/instructions.pdf"
            ),
            "original_eleven_class_table": (
                "https://old.sztaki.hu/~keri/n-arcs/n-arcs_in_PG(2,31).pdf"
            ),
            "corrected_twelve_class_classification": (
                "https://doi.org/10.1002/jcd.21410"
            ),
            "corrected_classification_repository_record": (
                "https://biblio.ugent.be/publication/7076091"
            ),
        },
        "input": {
            "path": str(input_path.resolve().relative_to(ROOT)),
            "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
            "format": "Keri A-list: eleven compact 2-by-18 matrices over F_31",
        },
        "field_order": P,
        "arc_size": ARC_SIZE,
        "secants_per_arc": math.comb(ARC_SIZE, 2),
        "outside_points_per_arc": P * P + P + 1 - ARC_SIZE,
        "outside_secant_incidence_moment": math.comb(ARC_SIZE, 2) * (P - 1),
        "corrected_projective_class_count": CORRECTED_CLASS_COUNT,
        "public_representative_count": len(representatives),
        "missing_representative_count": CORRECTED_CLASS_COUNT - len(representatives),
        "endpoint_required_c1": ENDPOINT_REQUIRED_C1,
        "index_one_point_counts_by_public_class": list(c1_sequence),
        "maximum_c1_over_public_classes": max(c1_sequence),
        "maximum_disjoint_unique_secant_matching_by_public_class": list(
            matching_sequence
        ),
        "classes_with_positive_c1": positive_c1_classes,
        "positive_c1_classes_with_all_index_one_points_on_one_secant": (
            same_secant_classes
        ),
        "all_public_classes_excluded_by_c1_requirement": all(
            c1 < ENDPOINT_REQUIRED_C1 for c1 in c1_sequence
        ),
        "all_public_classes_exactly_audited": len(representatives)
        == PUBLIC_CLASS_COUNT,
        "twelfth_class_representative_available": False,
        "twelfth_class_audited": False,
        "all_twelve_classes_excluded": False,
        "p31_endpoint_closed": False,
        "representatives": representatives,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="write JSON to stdout instead of the output path",
    )
    args = parser.parse_args()
    certificate = build_certificate(args.input)
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.stdout:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
