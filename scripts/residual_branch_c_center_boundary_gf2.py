#!/usr/bin/env python3
"""Exact GF(2) degree-boundary filter for p=31 Möbius centre fibres.

For a candidate common graph, the source-graph degree boundary must equal the
XOR of the sixteen hard affine centre lines.  A clean collision fixes two
half centers.  The remaining relaxation has ``14*30`` center columns and 15
fixed-antipodal-edge columns.  It includes one parity equation for every
mutable half and one for the fixed-edge choice.

This is only a necessary GF(2) relaxation: a consistent system need not have
an exactly-one solution over the Boolean variables, and inconsistency of a
bounded input catalog is not global closure of residual (ii).
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Edge,
    Functional,
    _functional_value,
    _negative_edge,
    projective_functionals,
)
from e1_gmin_m4_p31_direct_mobius_parallel_design import (  # noqa: E402
    P,
    _oriented_orbit_coefficients,
    _spatial_direction_index,
)
from io_atomic import write_json_atomic  # noqa: E402


POINT_COUNT = P * P
CENTER_OPTIONS = P - 1


def _functional(value: Sequence[int]) -> Functional:
    if len(value) != 2:
        raise ValueError("a functional needs two coordinates")
    out = int(value[0]) % P, int(value[1]) % P
    if out == (0, 0):
        raise ValueError("the zero functional is invalid")
    return out


def _parse_halves(value: object) -> tuple[tuple[Functional, Functional], ...]:
    if not isinstance(value, list) or len(value) != 16:
        raise ValueError("a p31 top design needs sixteen halves")
    out = []
    for row in value:
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError("each half must be [target, auxiliary]")
        out.append((_functional(row[0]), _functional(row[1])))
    return tuple(out)


def _point_index(point: tuple[int, int]) -> int:
    return point[0] * P + point[1]


def _boundary_mask(edges: Iterable[Edge]) -> int:
    mask = 0
    for edge in edges:
        mask ^= 1 << _point_index(edge[0])
        mask ^= 1 << _point_index(edge[1])
    return mask


def _line_mask(target: Functional, center: int) -> int:
    mask = 0
    for point in product(range(P), repeat=2):
        if _functional_value(P, target, point) == center % P:
            mask |= 1 << _point_index(point)
    if mask.bit_count() != P:
        raise ArithmeticError("an affine center line changed size")
    return mask


def _half_edges(
    target: Functional, auxiliary: Functional, center: int
) -> tuple[Edge, ...]:
    orbit_map = _oriented_orbit_coefficients(target, auxiliary, center)
    edges = tuple(
        orbit if coefficient == 1 else _negative_edge(P, orbit)
        for orbit, coefficient in orbit_map.items()
    )
    if len(edges) != P - 1 or len(set(edges)) != P - 1:
        raise ArithmeticError("a physical half lost an edge")
    return edges


def _fixed_edges(direction_index: int) -> tuple[Edge, ...]:
    direction = projective_functionals(P)[direction_index]
    out = set()
    for point in product(range(P), repeat=2):
        if point == (0, 0) or _functional_value(P, direction, point):
            continue
        negative = -point[0] % P, -point[1] % P
        out.add(tuple(sorted((point, negative))))
    edges = tuple(sorted(out))
    if len(edges) != (P - 1) // 2 or any(
        _spatial_direction_index(edge) != direction_index for edge in edges
    ):
        raise ArithmeticError("the fixed antipodal fibre changed")
    return edges


def _rank(rows: Iterable[int], column_count: int) -> int:
    basis: dict[int, int] = {}
    coefficient_mask = (1 << column_count) - 1
    for row in rows:
        value = row & coefficient_mask
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def _rank_and_consistency(
    rows: Sequence[int], column_count: int
) -> tuple[int, int, bool]:
    coefficient_rank = _rank(rows, column_count)
    augmented_rank = _rank(rows, column_count + 1)
    return coefficient_rank, augmented_rank, coefficient_rank == augmented_rank


def _row_contradiction_witness(
    rows: Sequence[int], column_count: int
) -> int | None:
    """Return a bitmask of input equations XORing to ``0 = 1``."""
    coefficient_mask = (1 << column_count) - 1
    basis: dict[int, tuple[int, int, int]] = {}
    for row_index, row in enumerate(rows):
        value = row & coefficient_mask
        right = (row >> column_count) & 1
        provenance = 1 << row_index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value, right, provenance
                break
            other, other_right, other_provenance = basis[pivot]
            value ^= other
            right ^= other_right
            provenance ^= other_provenance
        if not value and right:
            return provenance
    return None


class BoundaryColumnCache:
    """Cache the 30 exact boundary-plus-line masks for each physical half."""

    def __init__(self) -> None:
        self.halves: dict[
            tuple[Functional, Functional], tuple[int, ...]
        ] = {}
        self.fixed: dict[int, tuple[int, ...]] = {}

    def half_masks(
        self, target: Functional, auxiliary: Functional
    ) -> tuple[int, ...]:
        key = target, auxiliary
        if key not in self.halves:
            self.halves[key] = tuple(
                _boundary_mask(_half_edges(target, auxiliary, center))
                ^ _line_mask(target, center)
                for center in range(1, P)
            )
        return self.halves[key]

    def fixed_masks(self, direction_index: int) -> tuple[int, ...]:
        if direction_index not in self.fixed:
            self.fixed[direction_index] = tuple(
                _boundary_mask((edge,))
                for edge in _fixed_edges(direction_index)
            )
        return self.fixed[direction_index]


def _vector_rank_and_consistency(
    columns: Sequence[int], right_hand_side: int
) -> tuple[int, int, bool]:
    basis: dict[int, int] = {}

    def reduce(value: int, *, insert: bool) -> int:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                if insert:
                    basis[pivot] = value
                return value
            value ^= basis[pivot]
        return 0

    for column in columns:
        reduce(column, insert=True)
    coefficient_rank = len(basis)
    residual = reduce(right_hand_side, insert=False)
    return coefficient_rank, coefficient_rank + bool(residual), not residual


def boundary_system_columns(
    halves: Sequence[tuple[Functional, Functional]],
    collision_seed: dict[str, object],
    fixed_direction_index: int,
    cache: BoundaryColumnCache,
) -> dict[str, object]:
    """Build the same system by columns in equation space for fast batches."""
    collision_halves = tuple(
        int(value) for value in collision_seed["half_indices"]
    )
    collision_centers = tuple(
        int(value) for value in collision_seed["centers"]
    )
    if (
        len(halves) != 16
        or len(collision_halves) != 2
        or len(set(collision_halves)) != 2
        or len(collision_centers) != 2
    ):
        raise ValueError("the column system needs 16 halves and one fixed pair")
    frozen = dict(zip(collision_halves, collision_centers, strict=True))
    mutable = tuple(index for index in range(len(halves)) if index not in frozen)
    fixed_masks = cache.fixed_masks(fixed_direction_index)
    variable_count = len(mutable) * CENTER_OPTIONS + len(fixed_masks)
    if variable_count != 435:
        raise ArithmeticError("the column system lost 435 variables")

    collision_orbit = tuple(
        tuple(int(coordinate) for coordinate in point)
        for point in collision_seed["orbit"]
    )
    stated_coefficients = tuple(
        int(value) for value in collision_seed["coefficients"]
    )
    replayed_coefficients = []
    collision_edges = []
    right_hand_side = 0
    for half_index, center in zip(
        collision_halves, collision_centers, strict=True
    ):
        target, auxiliary = halves[half_index]
        orbit_map = _oriented_orbit_coefficients(target, auxiliary, center)
        if collision_orbit not in orbit_map:
            raise ArithmeticError("the frozen half lost its collision orbit")
        coefficient = orbit_map[collision_orbit]
        replayed_coefficients.append(coefficient)
        collision_edges.append(
            collision_orbit
            if coefficient == 1
            else _negative_edge(P, collision_orbit)
        )
        right_hand_side ^= cache.half_masks(target, auxiliary)[center - 1]
    if (
        tuple(replayed_coefficients) != stated_coefficients
        or sum(replayed_coefficients) != 0
    ):
        raise ArithmeticError("the collision coefficients failed replay")
    right_hand_side ^= _boundary_mask(collision_edges)

    columns = []
    for mutable_position, half_index in enumerate(mutable):
        target, auxiliary = halves[half_index]
        one_hot_bit = 1 << (POINT_COUNT + mutable_position)
        columns.extend(
            mask | one_hot_bit
            for mask in cache.half_masks(target, auxiliary)
        )
    fixed_hot_bit = 1 << (POINT_COUNT + len(mutable))
    columns.extend(mask | fixed_hot_bit for mask in fixed_masks)
    right_hand_side |= ((1 << (len(mutable) + 1)) - 1) << POINT_COUNT
    coefficient_rank, augmented_rank, consistent = (
        _vector_rank_and_consistency(columns, right_hand_side)
    )
    byte_count = (POINT_COUNT + len(mutable) + 1 + 7) // 8
    payload = b"".join(
        value.to_bytes(byte_count, "little")
        for value in (*columns, right_hand_side)
    )
    return {
        "variable_count": variable_count,
        "mutable_half_count": len(mutable),
        "center_variable_count": len(mutable) * CENTER_OPTIONS,
        "fixed_edge_variable_count": len(fixed_masks),
        "vertex_equation_count": POINT_COUNT,
        "one_hot_parity_equation_count": len(mutable) + 1,
        "equation_count": POINT_COUNT + len(mutable) + 1,
        "coefficient_rank": coefficient_rank,
        "augmented_rank": augmented_rank,
        "consistent": consistent,
        "collision_half_indices": collision_halves,
        "collision_centers": collision_centers,
        "collision_orbit": collision_orbit,
        "collision_coefficients": tuple(replayed_coefficients),
        "fixed_direction_index": fixed_direction_index,
        "matrix_augmented_columns_sha256": hashlib.sha256(payload).hexdigest(),
    }


def boundary_system(
    halves: Sequence[tuple[Functional, Functional]],
    collision_seed: dict[str, object],
    fixed_direction_index: int,
) -> dict[str, object]:
    if len(halves) != 16:
        raise ValueError("the system needs sixteen halves")
    collision_halves = tuple(
        int(value) for value in collision_seed["half_indices"]
    )
    collision_centers = tuple(
        int(value) for value in collision_seed["centers"]
    )
    if (
        len(collision_halves) != 2
        or len(set(collision_halves)) != 2
        or len(collision_centers) != 2
    ):
        raise ValueError("the clean collision must fix exactly two halves")
    frozen = dict(zip(collision_halves, collision_centers, strict=True))
    mutable = tuple(index for index in range(len(halves)) if index not in frozen)
    fixed_edges = _fixed_edges(fixed_direction_index)
    column_count = len(mutable) * CENTER_OPTIONS + len(fixed_edges)
    if column_count != 435:
        raise ArithmeticError("the corrected boundary system lost 435 columns")

    vertex_rows = [0] * POINT_COUNT
    constant = 0
    collision_edges = []
    collision_orbit = tuple(
        tuple(int(coordinate) for coordinate in point)
        for point in collision_seed["orbit"]
    )
    stated_coefficients = tuple(
        int(value) for value in collision_seed["coefficients"]
    )
    replayed_coefficients = []
    for half_index, center in zip(
        collision_halves, collision_centers, strict=True
    ):
        target, auxiliary = halves[half_index]
        orbit_map = _oriented_orbit_coefficients(target, auxiliary, center)
        if collision_orbit not in orbit_map:
            raise ArithmeticError("the frozen half lost the collision orbit")
        coefficient = orbit_map[collision_orbit]
        replayed_coefficients.append(coefficient)
        physical = (
            collision_orbit
            if coefficient == 1
            else _negative_edge(P, collision_orbit)
        )
        collision_edges.append(physical)
        constant ^= _boundary_mask(_half_edges(target, auxiliary, center))
        constant ^= _line_mask(target, center)
    if tuple(replayed_coefficients) != stated_coefficients or sum(
        replayed_coefficients
    ) != 0:
        raise ArithmeticError("the stated collision orientation failed replay")
    # Removing the two colliding physical edges is addition over GF(2).
    constant ^= _boundary_mask(collision_edges)

    column = 0
    mutable_column_ranges = []
    for half_index in mutable:
        start = column
        target, auxiliary = halves[half_index]
        for center in range(1, P):
            mask = _boundary_mask(_half_edges(target, auxiliary, center))
            mask ^= _line_mask(target, center)
            remaining = mask
            while remaining:
                low = remaining & -remaining
                vertex_rows[low.bit_length() - 1] |= 1 << column
                remaining ^= low
            column += 1
        mutable_column_ranges.append((half_index, start, column))
    fixed_start = column
    for edge in fixed_edges:
        mask = _boundary_mask((edge,))
        remaining = mask
        while remaining:
            low = remaining & -remaining
            vertex_rows[low.bit_length() - 1] |= 1 << column
            remaining ^= low
        column += 1
    if column != column_count:
        raise ArithmeticError("column assembly changed size")

    augmented_bit = 1 << column_count
    remaining = constant
    while remaining:
        low = remaining & -remaining
        vertex_rows[low.bit_length() - 1] |= augmented_bit
        remaining ^= low
    rows = list(vertex_rows)
    for _half_index, start, stop in mutable_column_ranges:
        rows.append(((1 << stop) - (1 << start)) | augmented_bit)
    rows.append(
        ((1 << column_count) - (1 << fixed_start)) | augmented_bit
    )

    coefficient_rank, augmented_rank, consistent = _rank_and_consistency(
        rows, column_count
    )
    contradiction = _row_contradiction_witness(rows, column_count)
    if consistent != (contradiction is None):
        raise ArithmeticError("rank and tracked contradiction disagree")
    contradiction_vertices = tuple(
        index
        for index in range(POINT_COUNT)
        if contradiction is not None and (contradiction >> index) & 1
    )
    contradiction_one_hot_rows = tuple(
        index - POINT_COUNT
        for index in range(POINT_COUNT, len(rows))
        if contradiction is not None and (contradiction >> index) & 1
    )
    row_payload = b"".join(
        int(row).to_bytes((column_count + 8) // 8, "little") for row in rows
    )
    return {
        "variable_count": column_count,
        "mutable_half_count": len(mutable),
        "center_variable_count": len(mutable) * CENTER_OPTIONS,
        "fixed_edge_variable_count": len(fixed_edges),
        "vertex_equation_count": POINT_COUNT,
        "one_hot_parity_equation_count": len(mutable) + 1,
        "equation_count": len(rows),
        "coefficient_rank": coefficient_rank,
        "augmented_rank": augmented_rank,
        "consistent": consistent,
        "contradiction_vertex_indices": contradiction_vertices,
        "contradiction_one_hot_row_indices": contradiction_one_hot_rows,
        "contradiction_equation_count": (
            0 if contradiction is None else contradiction.bit_count()
        ),
        "collision_half_indices": collision_halves,
        "collision_centers": collision_centers,
        "collision_orbit": collision_orbit,
        "collision_coefficients": tuple(replayed_coefficients),
        "fixed_direction_index": fixed_direction_index,
        "matrix_augmented_rows_sha256": hashlib.sha256(row_payload).hexdigest(),
    }


def _records(payload: dict[str, object]) -> tuple[dict[str, object], ...]:
    if isinstance(payload.get("designs"), list):
        source = payload["designs"]
    elif isinstance(payload.get("searchable_designs"), list):
        source = payload["searchable_designs"]
    else:
        raise ValueError("catalog has no supported design list")
    out = []
    for offset, row in enumerate(source):
        if not isinstance(row, dict):
            raise ValueError("a design record must be an object")
        seeds = row.get("clean_collision_seeds")
        if seeds is None and row.get("canonical_collision_seed") is not None:
            seeds = [row["canonical_collision_seed"]]
        if not isinstance(seeds, list) or not seeds:
            continue
        normalized = dict(row)
        normalized["design_index"] = int(
            row.get("design_index", row.get("component_id", offset))
        )
        normalized["halves"] = _parse_halves(row["halves"])
        normalized["collision_seeds"] = tuple(seeds)
        out.append(normalized)
    return tuple(out)


def _canonical_seed(record: dict[str, object]) -> dict[str, object]:
    seeds = record["collision_seeds"]
    assert isinstance(seeds, tuple)
    seed = min(
        seeds,
        key=lambda row: (
            tuple(int(value) for value in row["half_indices"]),
            tuple(int(value) for value in row["centers"]),
            tuple(tuple(int(x) for x in point) for point in row["orbit"]),
        ),
    )
    assert isinstance(seed, dict)
    return seed


def filter_catalog(
    path: Path,
    *,
    parity: str,
    minimum_depth: int | None,
) -> dict[str, object]:
    raw = path.read_bytes()
    payload = json.loads(raw)
    records = _records(payload)
    if parity != "all":
        residue = 1 if parity == "odd" else 0
        records = tuple(
            row for row in records if int(row["design_index"]) % 2 == residue
        )
    if minimum_depth is not None:
        records = tuple(
            row for row in records if int(row.get("depth", 0)) >= minimum_depth
        )
    results = []
    cache = BoundaryColumnCache()
    for ordinal, record in enumerate(records):
        seed = _canonical_seed(record)
        system = boundary_system_columns(
            record["halves"],
            seed,
            int(record.get("fixed_direction_index", seed["spatial_direction_index"])),
            cache,
        )
        results.append(
            {
                "design_index": int(record["design_index"]),
                "component_id": record.get("component_id"),
                "depth": record.get("depth"),
                "design_sha256": record.get("design_sha256"),
                "collision_seed": seed,
                "system": system,
            }
        )
        if ordinal % 100 == 0:
            print(
                json.dumps(
                    {
                        "filtered": ordinal + 1,
                        "last_design_index": record["design_index"],
                        "consistent_so_far": sum(
                            row["system"]["consistent"] for row in results
                        ),
                    }
                ),
                flush=True,
            )
    rank_histogram = Counter(
        (row["system"]["coefficient_rank"], row["system"]["augmented_rank"])
        for row in results
    )
    return {
        "schema": "residual_branch_c_center_boundary_gf2_v1",
        "classification": "exact necessary GF2 degree-boundary filter",
        "catalog": str(path),
        "catalog_sha256": hashlib.sha256(raw).hexdigest(),
        "catalog_schema": payload.get("schema"),
        "parity": parity,
        "minimum_depth": minimum_depth,
        "design_count": len(results),
        "consistent_design_count": sum(
            row["system"]["consistent"] for row in results
        ),
        "inconsistent_design_count": sum(
            not row["system"]["consistent"] for row in results
        ),
        "rank_histogram": {
            f"{coefficient}/{augmented}": count
            for (coefficient, augmented), count in sorted(rank_histogram.items())
        },
        "cached_distinct_half_count": len(cache.halves),
        "results": results,
        "catalog_exhausted": bool(payload.get("component_exhausted", True)),
        "residual_ii_closed": False,
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--parity", choices=("all", "odd", "even"), default="all")
    parser.add_argument("--minimum-depth", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.minimum_depth is not None and args.minimum_depth < 0:
        parser.error("minimum depth must be nonnegative")
    return args


def main(argv: Iterable[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = filter_catalog(
        args.catalog, parity=args.parity, minimum_depth=args.minimum_depth
    )
    if args.output:
        write_json_atomic(args.output, result)
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "schema",
                    "catalog_sha256",
                    "design_count",
                    "consistent_design_count",
                    "inconsistent_design_count",
                    "rank_histogram",
                    "catalog_exhausted",
                    "residual_ii_closed",
                )
            },
            indent=2,
            sort_keys=True,
        )
    )
    return result


if __name__ == "__main__":
    main()
