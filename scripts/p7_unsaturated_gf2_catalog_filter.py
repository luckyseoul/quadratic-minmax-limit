#!/usr/bin/env python3
"""GF(2) dependency filter for fixed p=7 unsaturated catalog choices.

Every exact score cardinality, degree parity, total edge count, distinguished
edge, and Paley-product condition has a mod-two shadow.  The coefficient
matrix is independent of the boundary and catalog values.  We row-reduce it
once with dependency witnesses, then test a complete direction catalog by
bitset syndrome evaluation.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_no_infinity_unsaturated_cpsat import (  # noqa: E402
    atomic_write,
    direction_target_options,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def equation_matrix() -> tuple[list[int], list[dict], int]:
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    rows = []
    metadata = []

    def add(indices, kind: str, **extra) -> None:
        mask = 0
        for index in indices:
            mask |= 1 << int(index)
        rows.append(mask)
        metadata.append({"kind": kind, **extra})

    add(range(len(edges)), "edge_count")
    add([edges.index((0, 1))], "distinguished_edge")
    for vertex in range(50):
        add(
            [index for index, edge in enumerate(edges) if vertex in edge],
            "degree",
            vertex=vertex,
        )
    add(
        [index for index, sign in enumerate(data["edge_signs"]) if int(sign) == -1],
        "paley_product",
    )
    for direction_index, direction in enumerate(projective_directions(7)):
        eps, labels = field_direction_data(7, direction)
        for point_index, X in enumerate(itertools.combinations(range(7), 4)):
            X_set = set(X)
            bad = []
            for edge_index, (a, endpoint) in enumerate(edges):
                y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                y_b = 1 if labels[endpoint - 1] in X_set else -1
                if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                    bad.append(edge_index)
            add(
                bad,
                "score_bad_count",
                direction_index=direction_index,
                point_index=point_index,
            )
    if len(rows) != 333:
        raise AssertionError(f"expected 333 parity equations, got {len(rows)}")
    return rows, metadata, len(edges)


def row_dependencies(rows: list[int]) -> tuple[int, tuple[int, ...]]:
    pivots: dict[int, tuple[int, int]] = {}
    dependencies = []
    for row_index, original in enumerate(rows):
        mask = original
        witness = 1 << row_index
        while mask:
            pivot = mask.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (mask, witness)
                break
            basis_mask, basis_witness = pivots[pivot]
            mask ^= basis_mask
            witness ^= basis_witness
        if mask == 0:
            dependencies.append(witness)
    if len(dependencies) != len(rows) - len(pivots):
        raise AssertionError("GF(2) dependency dimension mismatch")
    return len(pivots), tuple(dependencies)


def direction_scope(
    c_h: int,
    boundary: tuple[int, ...],
    elevated: tuple[int, ...],
) -> tuple[list[dict], dict[int, int]]:
    rows = []
    type_floors = {-1: 0, 1: 0}
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, len(B), phase)
        type_floors[eps] += floor
        rows.append(
            {
                "direction": direction,
                "eps": eps,
                "labels": labels,
                "B": B,
                "phase": phase,
                "floor": floor,
            }
        )
    if any(value not in (24, 32) for value in type_floors.values()):
        raise ValueError(f"boundary is outside the surviving scope: {type_floors}")
    elevated_set = set(elevated)
    for eps in (-1, 1):
        expected = 1 if type_floors[eps] == 24 else 0
        observed = sum(
            index in elevated_set and int(row["eps"]) == eps
            for index, row in enumerate(rows)
        )
        if observed != expected:
            raise ValueError(f"type {eps} needs {expected} elevated directions")
    return rows, type_floors


def scan_catalog(
    c_h: int,
    boundary: tuple[int, ...],
    elevated: tuple[int, ...],
    scan_direction: int,
    fixed_indices: dict[int, int],
) -> dict:
    started = time.time()
    matrix, metadata, n_variables = equation_matrix()
    rank, dependencies = row_dependencies(matrix)
    direction_rows, type_floors = direction_scope(c_h, boundary, elevated)
    if not 0 <= scan_direction < 8 or scan_direction not in elevated:
        raise ValueError("scan direction must be an elevated direction")

    option_tables = []
    for index, row in enumerate(direction_rows):
        options = direction_target_options(
            len(row["B"]),
            int(row["phase"]),
            set(row["B"]),
            type_floors[int(row["eps"])],
            index in set(elevated),
        )
        if index != scan_direction and len(options) > 1:
            if index not in fixed_indices:
                raise ValueError(
                    f"direction {index} has {len(options)} choices and needs --fixed-index"
                )
            fixed_index = int(fixed_indices[index])
            if not 0 <= fixed_index < len(options):
                raise ValueError(f"fixed index outside direction {index}")
        option_tables.append(options)
    scan_options = option_tables[scan_direction]

    base_rhs = 0

    def set_rhs(row_index: int, value: int) -> None:
        nonlocal base_rhs
        if value & 1:
            base_rhs |= 1 << row_index

    set_rhs(0, 29)
    set_rhs(1, 1)
    for vertex in range(50):
        set_rhs(2 + vertex, int(vertex in set(boundary)))
    set_rhs(52, int(c_h == -1))

    pair_order = tuple(itertools.combinations(range(7), 2))

    def score_rhs_bits(direction_index: int, option: tuple[int, ...]) -> int:
        bits = 0
        constant = int(option[1])
        pairs = tuple(int(value) for value in option[2:])
        for point_index, X in enumerate(itertools.combinations(range(7), 4)):
            normalized_score = constant + sum(
                pairs[pair_index]
                * (1 if ((s in X) == (t in X)) else -1)
                for pair_index, (s, t) in enumerate(pair_order)
            )
            if normalized_score < 3 or normalized_score % 2 == 0:
                raise AssertionError("catalog row has an invalid normalized score")
            bad_count = (29 - normalized_score) // 2
            if bad_count & 1:
                bits |= 1 << (53 + 35 * direction_index + point_index)
        return bits

    for direction_index, options in enumerate(option_tables):
        if direction_index == scan_direction:
            continue
        option_index = int(fixed_indices.get(direction_index, 0))
        base_rhs ^= score_rhs_bits(direction_index, options[option_index])

    fixed_syndrome = tuple(
        (dependency & base_rhs).bit_count() & 1 for dependency in dependencies
    )
    consistent = []
    syndrome_histogram = Counter()
    for option_index, option in enumerate(scan_options):
        contribution = score_rhs_bits(scan_direction, option)
        syndrome = tuple(
            fixed_syndrome[index]
            ^ ((dependency & contribution).bit_count() & 1)
            for index, dependency in enumerate(dependencies)
        )
        packed = sum(bit << index for index, bit in enumerate(syndrome))
        syndrome_histogram[packed] += 1
        if packed == 0:
            consistent.append(option_index)

    return {
        "experiment": "p7_unsaturated_gf2_catalog_filter",
        "status": "complete_exact_mod_two_catalog_scan",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "fixed_elevated_directions": list(elevated),
        "scan_direction": scan_direction,
        "catalog_total": len(scan_options),
        "gf2_variables": n_variables,
        "gf2_equations": len(matrix),
        "gf2_rank": rank,
        "gf2_dependency_dimension": len(dependencies),
        "consistent_catalog_rows": len(consistent),
        "inconsistent_catalog_rows": len(scan_options) - len(consistent),
        "consistent_catalog_indices": (
            "ALL" if len(consistent) == len(scan_options) else consistent
        ),
        "syndrome_count": len(syndrome_histogram),
        "syndrome_histogram": [
            {"syndrome": str(value), "count": count}
            for value, count in sorted(syndrome_histogram.items())
        ],
        "fixed_indices": {
            str(direction): index for direction, index in sorted(fixed_indices.items())
        },
        "elapsed_seconds": time.time() - started,
        "metadata_kinds": dict(Counter(row["kind"] for row in metadata)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--elevated-directions", type=int, nargs="+", required=True)
    parser.add_argument("--scan-direction", type=int, required=True)
    parser.add_argument(
        "--fixed-index",
        type=int,
        nargs=2,
        action="append",
        metavar=("DIRECTION", "INDEX"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixed_indices = {}
    for direction, index in args.fixed_index or []:
        if direction in fixed_indices:
            raise ValueError(f"duplicate fixed index for direction {direction}")
        fixed_indices[direction] = index
    out = scan_catalog(
        args.c_h,
        tuple(sorted(args.fixed_boundary)),
        tuple(sorted(args.elevated_directions)),
        args.scan_direction,
        fixed_indices,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
