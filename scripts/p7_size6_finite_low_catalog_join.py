#!/usr/bin/env python3
"""Exact modular hash joins for unresolved low-catalog deep allocations."""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_size6_finite_deep_modular_cpsat import load_source, source_hash  # noqa: E402
from p7_size6_finite_ordinary_mod7 import (  # noqa: E402
    contribution_matrices,
    count_zero_sum_tuples_multi,
    direction_scope,
)
from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies  # noqa: E402


def count_up_to_three(
    moduli: tuple[int, ...],
    base: tuple[np.ndarray, ...],
    contributions: list[tuple[np.ndarray, ...]],
) -> int:
    if len(contributions) <= 2:
        return count_zero_sum_tuples_multi(moduli, base, contributions)
    if len(contributions) != 3:
        raise ValueError("low-catalog join expects at most three variable catalogs")
    contributions = sorted(contributions, key=lambda item: item[0].shape[1])
    first, second, third = contributions
    pair_counts: Counter[bytes] = Counter()
    for second_index in range(second[0].shape[1]):
        blocks = [
            ((a.astype(np.int16) + b[:, second_index, None].astype(np.int16)) % modulus).T.astype(np.uint8)
            for modulus, a, b in zip(moduli, first, second)
        ]
        pair_counts.update(bytes(row) for row in np.concatenate(blocks, axis=1))
    needed_blocks = [
        ((-base_row[:, None].astype(np.int16) - matrix.astype(np.int16)) % modulus).T.astype(np.uint8)
        for modulus, base_row, matrix in zip(moduli, base, third)
    ]
    return sum(pair_counts.get(bytes(row), 0) for row in np.concatenate(needed_blocks, axis=1))


def run(
    source: Path,
    allocation_shards: tuple[Path, ...],
    moduli: tuple[int, ...],
    shard_index: int,
    shard_count: int,
) -> dict:
    started = time.time()
    payload = load_source(source)
    digest = source_hash(source)
    unknown = []
    seen = set()
    for path in allocation_shards:
        recording = json.loads(path.read_text())
        if recording["source_sha256"] != digest:
            raise ValueError("allocation shard source hash mismatch")
        for orbit in recording["rows"]:
            for leaf in orbit["leaves"]:
                if leaf["solver_status"] != "UNKNOWN":
                    continue
                key = (
                    int(leaf["orbit_index"]),
                    tuple(sorted((int(k), int(v)) for k, v in leaf["fixed_scaled_means"].items())),
                )
                if key in seen:
                    raise ValueError("duplicate unknown allocation leaf")
                seen.add(key)
                unknown.append((key, leaf))
    unknown.sort(key=lambda item: item[0])
    selected = [item for index, item in enumerate(unknown) if index % shard_count == shard_index]

    matrix = equation_matrix()
    dependencies = {}
    linear_rows = []
    for modulus in moduli:
        rank, dependency = left_dependencies(matrix, modulus)
        if np.any(dependency @ (matrix % modulus) % modulus):
            raise AssertionError("left-null dependency audit failed")
        dependencies[modulus] = dependency
        linear_rows.append(
            {"modulus": modulus, "rank": rank, "left_dependency_dimension": int(len(dependency))}
        )
    base_edge = tuple(
        (dependency[:, :2] @ np.asarray([29, 1], dtype=np.int64) % modulus).astype(np.uint8)
        for modulus, dependency in dependencies.items()
    )
    cache: dict[tuple, tuple[np.ndarray, ...]] = {}
    rows = []
    for (orbit_index, fixed_items), leaf in selected:
        fixed_means = dict(fixed_items)
        orbit = payload["orbits"][orbit_index]
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        direction_rows, type_floors = direction_scope(int(payload["c_H"]), boundary)
        if any(
            fixed_means[index] - int(row["floor"]) not in (0, 8)
            for index, row in enumerate(direction_rows)
        ):
            raise ValueError("unknown leaf is not a floor-or-plus-eight catalog case")
        if any(
            sum(fixed_means[index] for index, row in enumerate(direction_rows) if int(row["eps"]) == eps) != 32
            for eps in (-1, 1)
        ):
            raise ValueError("fixed means do not satisfy the exact type sums")
        base = tuple(row.astype(np.int16).copy() for row in base_edge)
        variable = []
        metadata = []
        for direction_index, row in enumerate(direction_rows):
            mean = fixed_means[direction_index]
            key = (
                direction_index,
                len(row["B"]),
                row["phase"],
                mean,
                tuple(sorted(row["B"])),
            )
            if key not in cache:
                cache[key] = contribution_matrices(
                    dependencies,
                    direction_index,
                    len(row["B"]),
                    row["phase"],
                    mean,
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
                metadata.append(
                    {
                        "direction": direction_index,
                        "b": len(row["B"]),
                        "phase": row["phase"],
                        "scaled_mean": mean,
                        "catalog_rows": int(contribution[0].shape[1]),
                    }
                )
        consistent = count_up_to_three(moduli, base, variable)
        rows.append(
            {
                "orbit_index": orbit_index,
                "fixed_scaled_means": {str(key): value for key, value in fixed_items},
                "type_floor_sums": {str(key): value for key, value in type_floors.items()},
                "variable_catalogs": metadata,
                "catalog_pattern": sorted((item[0].shape[1] for item in variable), reverse=True),
                "modular_consistent_catalog_tuples": int(consistent),
                "modularly_infeasible": consistent == 0,
            }
        )
    return {
        "experiment": "p7_size6_finite_low_catalog_join",
        "status": "complete_exact_multi_prime_up_to_three_catalog_join_shard",
        "source": str(source),
        "source_sha256": digest,
        "allocation_shards": [str(path) for path in allocation_shards],
        "moduli": list(moduli),
        "linear_system": linear_rows,
        "unknown_leaves_in_source": len(unknown),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "processed_leaves": len(rows),
        "modularly_infeasible_leaves": sum(row["modularly_infeasible"] for row in rows),
        "surviving_leaves": sum(not row["modularly_infeasible"] for row in rows),
        "syndrome_contribution_cache_entries": len(cache),
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--allocation-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 7))
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run(
        args.source,
        tuple(args.allocation_shards),
        tuple(args.moduli),
        args.shard_index,
        args.shard_count,
    )
    atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
