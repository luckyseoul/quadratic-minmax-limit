#!/usr/bin/env python3
"""Exact multi-prime modular sieve for ordinary p=7 six-finite orbits.

The floor quotient has two strata.  This script handles the ordinary stratum,
where each quadratic type has floor sum 24 or 32.  A sum-24 type must elevate
exactly one of its four directions by eight scaled-mean units; a sum-32 type
stays at its floor.  Consequently each type contributes at most one
non-singleton complete Johnson-slice catalog.  The 135 left dependencies of
the common 282-by-1225 edge-count system can therefore be checked by a
one-table lookup or a two-table hash join, without a Cartesian expansion.

Every arithmetic operation in the exclusion is exact over the requested
prime fields.  A catalog tuple must use the same row indices in every field,
so the joined signatures are stronger than separate per-prime counts.  Cases
that survive this necessary condition are retained for later exact
constraints; this script does not claim that a survivor is feasible.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_unsaturated_mod7_batch import mapped_slack_catalog  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direction_scope(c_h: int, boundary: tuple[int, ...]) -> tuple[list[dict], dict[int, int]]:
    """Rebuild all eight exact direction rows for a six-finite boundary."""
    if c_h not in (-1, 1) or len(boundary) != 6 or len(set(boundary)) != 6:
        raise ValueError("need c_H=+/-1 and six distinct finite vertices")
    rows = []
    type_floors = {-1: 0, 1: 0}
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            if not 1 <= vertex <= 49:
                raise ValueError("finite vertices must lie in 1..49")
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        floor = scaled_direction_floor(7, len(B), phase)
        type_floors[int(eps)] += int(floor)
        rows.append(
            {
                "eps": int(eps),
                "B": B,
                "phase": int(phase),
                "floor": int(floor),
            }
        )
    return rows, type_floors


def elevation_cases(rows: list[dict], type_floors: dict[int, int]) -> tuple[tuple[int, ...], ...]:
    """All exact +8 allocations for a floor-(24 or 32) ordinary orbit."""
    choices: list[tuple[int | None, ...]] = []
    for eps in (-1, 1):
        directions = tuple(index for index, row in enumerate(rows) if row["eps"] == eps)
        if len(directions) != 4:
            raise AssertionError("each quadratic type must have four directions")
        if type_floors[eps] == 24:
            choices.append(directions)
        elif type_floors[eps] == 32:
            choices.append((None,))
        else:
            raise ValueError(f"non-ordinary type floor sum {type_floors[eps]}")
    return tuple(
        tuple(sorted(value for value in choice if value is not None))
        for choice in itertools.product(*choices)
    )


def contribution_matrices(
    dependencies: dict[int, np.ndarray],
    direction_index: int,
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
    B: set[int],
) -> tuple[np.ndarray, ...]:
    """Return one exact catalog-syndrome matrix for every requested prime."""
    values = mapped_slack_catalog(odd_fibres, phase, scaled_mean, B).astype(np.int64)
    bad_counts = 13 - values
    if np.any(bad_counts < 0):
        raise AssertionError("catalog score exceeds the 29-edge range")
    block = slice(2 + 35 * direction_index, 2 + 35 * (direction_index + 1))
    return tuple(
        (
            dependency[:, block] @ (bad_counts.T % modulus) % modulus
        ).astype(np.uint8)
        for modulus, dependency in dependencies.items()
    )


def count_zero_sum_tuples_multi(
    moduli: tuple[int, ...],
    base: tuple[np.ndarray, ...],
    contributions: list[tuple[np.ndarray, ...]],
) -> int:
    """Count catalog tuples satisfying every prime-modulus dependency."""
    if not contributions:
        return int(all(not np.any(row % modulus) for modulus, row in zip(moduli, base)))
    if len(contributions) == 1:
        passing = np.ones(contributions[0][0].shape[1], dtype=bool)
        for modulus, base_row, matrix in zip(moduli, base, contributions[0]):
            target = (-base_row.astype(np.int16)) % modulus
            passing &= np.all(matrix == target[:, None], axis=0)
        return int(np.count_nonzero(passing))
    if len(contributions) != 2:
        raise AssertionError("ordinary cases must have at most two variable catalogs")
    first, second = contributions
    first_rows = np.concatenate([matrix.T for matrix in first], axis=1)
    first_counts = Counter(bytes(row) for row in first_rows)
    needed_blocks = []
    for modulus, base_row, matrix in zip(moduli, base, second):
        needed_blocks.append(
            ((-base_row[:, None].astype(np.int16) - matrix.astype(np.int16)) % modulus).T.astype(np.uint8)
        )
    needed_rows = np.concatenate(needed_blocks, axis=1)
    return sum(first_counts.get(bytes(row), 0) for row in needed_rows)


def run(
    source: Path,
    c_h: int,
    moduli: tuple[int, ...],
    shard_index: int,
    shard_count: int,
    max_orbits: int | None,
) -> dict:
    started = time.time()
    payload = json.loads(source.read_text())
    if (
        int(payload.get("p", 0)) != 7
        or int(payload.get("c_H", 0)) != c_h
        or int(payload.get("boundary_size", 0)) != 6
        or int(payload.get("infinity_value", -1)) != 0
        or payload.get("status") != "complete_boundary_only_enumeration"
    ):
        raise ValueError("source is not the complete p=7 six-finite orbit quotient")

    if not moduli or len(set(moduli)) != len(moduli):
        raise ValueError("need distinct prime moduli")
    matrix = equation_matrix()
    dependencies = {}
    linear_rows = []
    manufactured = np.arange(matrix.shape[1], dtype=np.int64) % 2
    for modulus in moduli:
        rank, dependency = left_dependencies(matrix, modulus)
        if np.any(dependency @ (matrix % modulus) % modulus):
            raise AssertionError(f"left-null dependency audit failed modulo {modulus}")
        if np.any(dependency @ ((matrix @ manufactured) % modulus) % modulus):
            raise AssertionError(f"manufactured right side rejected modulo {modulus}")
        dependencies[modulus] = dependency
        linear_rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "left_dependency_dimension": int(len(dependency)),
                "left_null_audit": True,
                "manufactured_rhs_calibration": True,
            }
        )

    base_edge = tuple(
        (dependency[:, :2] @ np.asarray([29, 1], dtype=np.int64) % modulus).astype(np.uint8)
        for modulus, dependency in dependencies.items()
    )
    cache: dict[tuple, tuple[np.ndarray, ...]] = {}
    ordinary = []
    deep_orbits = 0
    for orbit_index, orbit in enumerate(payload["orbits"]):
        costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
        if all(value in (24, 32) for value in costs.values()):
            if orbit_index % shard_count == shard_index:
                ordinary.append((orbit_index, orbit))
        else:
            deep_orbits += 1
    total_ordinary = sum(
        all(int(value) in (24, 32) for value in orbit["type_costs"].values())
        for orbit in payload["orbits"]
    )
    if max_orbits is not None:
        ordinary = ordinary[:max_orbits]

    pattern_counts: Counter[tuple[int, ...]] = Counter()
    floor_pair_counts: Counter[tuple[int, int]] = Counter()
    elevation_case_count = 0
    rejected_case_count = 0
    surviving = []
    for orbit_index, orbit in ordinary:
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        rows, type_floors = direction_scope(c_h, boundary)
        recorded = {int(key): int(value) for key, value in orbit["type_costs"].items()}
        if type_floors != recorded:
            raise AssertionError("rebuilt type-floor sums disagree with orbit source")
        floor_pair_counts[(type_floors[-1], type_floors[1])] += 1
        for elevated in elevation_cases(rows, type_floors):
            elevation_case_count += 1
            elevated_set = set(elevated)
            base = tuple(row.astype(np.int16).copy() for row in base_edge)
            variable = []
            variable_metadata = []
            for direction_index, row in enumerate(rows):
                scaled_mean = row["floor"] + (8 if direction_index in elevated_set else 0)
                key = (
                    direction_index,
                    len(row["B"]),
                    row["phase"],
                    scaled_mean,
                    tuple(sorted(row["B"])),
                )
                if key not in cache:
                    cache[key] = contribution_matrices(
                        dependencies,
                        direction_index,
                        len(row["B"]),
                        row["phase"],
                        scaled_mean,
                        set(row["B"]),
                    )
                contribution = cache[key]
                if contribution[0].shape[1] == 1:
                    base = tuple(
                        (base_row + matrix_row[:, 0]) % modulus
                        for modulus, base_row, matrix_row in zip(moduli, base, contribution)
                    )
                else:
                    variable.append(contribution)
                    variable_metadata.append(
                        {
                            "direction": direction_index,
                            "eps": row["eps"],
                            "b": len(row["B"]),
                            "phase": row["phase"],
                            "scaled_mean": scaled_mean,
                            "catalog_rows": int(contribution[0].shape[1]),
                        }
                    )
            pattern = tuple(sorted((item[0].shape[1] for item in variable), reverse=True))
            if len(variable) > 2:
                raise AssertionError("ordinary case produced more than two variable catalogs")
            pattern_counts[pattern] += 1
            consistent = count_zero_sum_tuples_multi(moduli, base, variable)
            if consistent == 0:
                rejected_case_count += 1
            else:
                surviving.append(
                    {
                        "orbit_index": orbit_index,
                        "orbit_size": int(orbit["size"]),
                        "representative_vertices": list(boundary),
                        "type_floor_sums": {str(key): value for key, value in type_floors.items()},
                        "elevated_directions": list(elevated),
                        "catalog_pattern": list(pattern),
                        "variable_catalogs": variable_metadata,
                        "mod7_consistent_catalog_tuples": int(consistent),
                    }
                )

    return {
        "experiment": "p7_size6_finite_ordinary_modular",
        "status": "complete_exact_multi_prime_ordinary_orbit_sieve",
        "p": 7,
        "c_H": c_h,
        "source": str(source),
        "source_sha256": source_hash(source),
        "source_orbits": int(payload["orbit_count"]),
        "ordinary_orbits_in_source": total_ordinary,
        "deep_deficit_orbits_in_source": deep_orbits,
        "processed_ordinary_orbits": len(ordinary),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "max_orbits": max_orbits,
        "floor_pair_counts": {f"{a},{b}": count for (a, b), count in sorted(floor_pair_counts.items())},
        "elevation_cases": elevation_case_count,
        "catalog_pattern_counts": {
            "x".join(map(str, pattern)) if pattern else "singleton": count
            for pattern, count in sorted(pattern_counts.items())
        },
        "syndrome_contribution_cache_entries": len(cache),
        "modular_infeasible_cases": rejected_case_count,
        "surviving_cases": len(surviving),
        "all_processed_cases_modularly_infeasible": not surviving,
        "linear_system": {
            "equations": int(matrix.shape[0]),
            "edge_variables": int(matrix.shape[1]),
            "moduli": linear_rows,
        },
        "survivors": surviving,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--moduli", type=int, nargs="+", default=(7,))
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-orbits", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run(
        args.source,
        args.c_h,
        tuple(args.moduli),
        args.shard_index,
        args.shard_count,
        args.max_orbits,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "survivors"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
