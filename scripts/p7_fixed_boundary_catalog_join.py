#!/usr/bin/env python3
"""Exact multi-prime slack-catalog join for one fixed p=7 mean allocation."""
from __future__ import annotations

import argparse
from collections import Counter
import functools
import hashlib
import itertools
import json
import math
import sys
import time
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
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_unsaturated_slack_catalog import exact_slack_catalog_values  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))
POINT_INDEX = {point: index for index, point in enumerate(POINTS)}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@functools.lru_cache(maxsize=4)
def cached_zero_catalog(catalog_path_text: str, summary_path_text: str) -> np.ndarray:
    catalog_path = Path(catalog_path_text)
    summary_path = Path(summary_path_text)
    summary = json.loads(summary_path.read_text())
    digest = file_sha256(catalog_path)
    if (
        summary.get("experiment") != "p7_slack_catalog_first_lift_shards"
        or summary.get("status") != "complete_exact_disjoint_shard_enumeration"
        or int(summary.get("scaled_mean", -1)) != 16
        or int(summary.get("solution_count", -1)) != 575_407
        or summary.get("merged_sha256") != digest
    ):
        raise ValueError("mean-16 catalog cache failed its exact shard summary")
    values = np.load(catalog_path, allow_pickle=False)
    if values.shape != (575_407, 35) or np.any(
        values.sum(axis=1, dtype=np.int64) != 40
    ):
        raise AssertionError("mean-16 catalog cache failed its shape/mean audit")
    return values.astype(np.int16, copy=False)


def canonical_catalog(
    b: int,
    phase: int,
    mean: int,
    cache: tuple[Path, Path] | None,
) -> np.ndarray:
    if b in (1, 2, 5, 6):
        minimum_mean = 8 if phase == 0 else 6
        if mean >= minimum_mean:
            B = set(range(b))
            parity = np.asarray(
                [
                    (sum(value in B for value in point) + phase) & 1
                    for point in POINTS
                ],
                dtype=np.int16,
            )
            excess_mean = mean - minimum_mean
            excess = canonical_catalog(0, 0, excess_mean, cache)
            return excess + parity[None, :]
    if b == 0 and phase == 0 and mean == 16 and cache is not None:
        return cached_zero_catalog(str(cache[0]), str(cache[1]))
    return np.asarray(exact_slack_catalog_values(b, phase, mean), dtype=np.int16)


def mapped_catalog(
    b: int,
    phase: int,
    mean: int,
    B: set[int],
    cache: tuple[Path, Path] | None,
) -> np.ndarray:
    canonical = canonical_catalog(b, phase, mean, cache)
    actual_B = sorted(B)
    actual_complement = sorted(set(range(7)) - B)
    permutation = dict(zip(range(b), actual_B)) | dict(
        zip(range(b, 7), actual_complement)
    )
    inverse = {target: source for source, target in permutation.items()}
    source_columns = [
        POINT_INDEX[tuple(sorted(inverse[value] for value in point))]
        for point in POINTS
    ]
    mapped = canonical[:, source_columns]
    mapped = mapped[np.all(mapped <= 13, axis=1)]
    if len({row.tobytes() for row in mapped}) != len(mapped):
        raise AssertionError("mapped catalog contains duplicate rows")
    return mapped


def direction_rows(c_h: int, boundary: tuple[int, ...]) -> list[dict]:
    rows = []
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        counts = [0] * 7
        for vertex in boundary:
            counts[labels[vertex - 1]] += 1
        B = {index for index, count in enumerate(counts) if count & 1}
        phase = int(-eps * c_h == -1)
        rows.append(
            {
                "direction": list(direction),
                "eps": int(eps),
                "B": B,
                "b": len(B),
                "phase": phase,
                "floor": int(scaled_direction_floor(7, len(B), phase)),
            }
        )
    return rows


def signature_bytes(blocks: tuple[np.ndarray, ...]) -> bytes:
    return b"".join(np.ascontiguousarray(block, dtype=np.uint8).tobytes() for block in blocks)


def signature_sums(
    contributions: tuple[tuple[np.ndarray, ...], ...],
    moduli: tuple[int, ...],
):
    zero = tuple(np.zeros(matrix.shape[0], dtype=np.uint8) for matrix in contributions[0])

    def recurse(index: int, current: tuple[np.ndarray, ...]):
        if index == len(contributions):
            yield current
            return
        matrices = contributions[index]
        for column in range(matrices[0].shape[1]):
            updated = tuple(
                ((row.astype(np.int16) + matrix[:, column].astype(np.int16)) % modulus).astype(np.uint8)
                for modulus, row, matrix in zip(moduli, current, matrices)
            )
            yield from recurse(index + 1, updated)

    yield from recurse(0, zero)


def choose_partition(sizes: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if len(sizes) < 2:
        raise ValueError("partition requires at least two catalogs")
    best = None
    indices = tuple(range(len(sizes)))
    for count in range(1, len(sizes)):
        for first in itertools.combinations(indices, count):
            if 0 not in first:
                continue
            second = tuple(index for index in indices if index not in first)
            products = (
                math.prod(sizes[index] for index in first),
                math.prod(sizes[index] for index in second),
            )
            score = (max(products), min(products), first)
            if best is None or score < best[0]:
                best = (score, first, second)
    if best is None:
        raise AssertionError("catalog partition search failed")
    return best[1], best[2]


def count_join(
    moduli: tuple[int, ...],
    base: tuple[np.ndarray, ...],
    contributions: tuple[tuple[np.ndarray, ...], ...],
) -> tuple[int, dict]:
    if not contributions:
        consistent = all(not np.any(row % modulus) for modulus, row in zip(moduli, base))
        return int(consistent), {"partition": [], "enumerated_signatures": [1, 1]}
    if len(contributions) == 1:
        count = 0
        matrix = contributions[0]
        for column in range(matrix[0].shape[1]):
            if all(
                np.all((base_row.astype(np.int16) + block[:, column]) % modulus == 0)
                for modulus, base_row, block in zip(moduli, base, matrix)
            ):
                count += 1
        return count, {
            "partition": [[0], []],
            "enumerated_signatures": [matrix[0].shape[1], 1],
        }

    sizes = tuple(item[0].shape[1] for item in contributions)
    first_indices, second_indices = choose_partition(sizes)
    first = tuple(contributions[index] for index in first_indices)
    second = tuple(contributions[index] for index in second_indices)
    first_counts = Counter(signature_bytes(row) for row in signature_sums(first, moduli))
    total = 0
    second_count = 0
    if len(first) == 1 and len(second) == 2:
        inner, outer = second
        for outer_column in range(outer[0].shape[1]):
            needed_blocks = [
                np.mod(
                    -base_row[:, None].astype(np.int16)
                    - inner_matrix.astype(np.int16)
                    - outer_matrix[:, outer_column, None].astype(np.int16),
                    modulus,
                ).T.astype(np.uint8)
                for modulus, base_row, inner_matrix, outer_matrix in zip(
                    moduli, base, inner, outer
                )
            ]
            packed = np.concatenate(needed_blocks, axis=1)
            total += sum(first_counts.get(row.tobytes(), 0) for row in packed)
            second_count += len(packed)
        return total, {
            "partition": [list(first_indices), list(second_indices)],
            "catalog_sizes": list(sizes),
            "enumerated_signatures": [sum(first_counts.values()), second_count],
            "distinct_first_signatures": len(first_counts),
            "vectorized_two_catalog_probe": True,
        }
    for row in signature_sums(second, moduli):
        needed = tuple(
            ((-base_row.astype(np.int16) - block.astype(np.int16)) % modulus).astype(np.uint8)
            for modulus, base_row, block in zip(moduli, base, row)
        )
        total += first_counts.get(signature_bytes(needed), 0)
        second_count += 1
    return total, {
        "partition": [list(first_indices), list(second_indices)],
        "catalog_sizes": list(sizes),
        "enumerated_signatures": [sum(first_counts.values()), second_count],
        "distinct_first_signatures": len(first_counts),
    }


def run(
    c_h: int,
    boundary: tuple[int, ...],
    means: tuple[int, ...],
    moduli: tuple[int, ...],
    cache: tuple[Path, Path] | None = None,
) -> dict:
    started = time.time()
    boundary = tuple(sorted(int(value) for value in boundary))
    if c_h not in (-1, 1) or len(means) != 8:
        raise ValueError("need c_H=+/-1 and eight fixed means")
    if not boundary or len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("boundary must be a distinct nonempty even set")
    if not all(1 <= vertex <= 49 for vertex in boundary):
        raise ValueError("this join accepts finite vertices 1..49 only")
    rows = direction_rows(c_h, boundary)
    if any(means[index] < row["floor"] for index, row in enumerate(rows)):
        raise ValueError("fixed mean lies below its parity floor")
    if any(
        sum(means[index] for index, row in enumerate(rows) if row["eps"] == eps)
        != 32
        for eps in (-1, 1)
    ):
        raise ValueError("fixed means violate an exact type sum")
    if any(
        len(
            {
                means[index] % 8
                for index, row in enumerate(rows)
                if row["eps"] == eps
            }
        )
        != 1
        for eps in (-1, 1)
    ):
        raise ValueError("fixed means violate a common type residue")

    matrix, dependencies, linear_rows = linear_data(moduli)
    base = tuple(
        (
            dependency[:, :2]
            @ np.asarray([29, 1], dtype=np.int64)
            % modulus
        ).astype(np.uint8)
        for modulus, dependency in dependencies.items()
    )
    variable = []
    metadata = []
    for direction_index, (row, mean) in enumerate(zip(rows, means)):
        values = mapped_catalog(
            row["b"], row["phase"], mean, set(row["B"]), cache
        )
        bad = 13 - values.astype(np.int64)
        contribution = tuple(
            (
                dependency[
                    :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
                ]
                @ (bad.T % modulus)
                % modulus
            ).astype(np.uint8)
            for modulus, dependency in dependencies.items()
        )
        metadata.append(
            {
                "direction_index": direction_index,
                "eps": row["eps"],
                "b": row["b"],
                "phase": row["phase"],
                "scaled_mean": mean,
                "catalog_rows_with_pointwise_bound": int(len(values)),
            }
        )
        if len(values) == 1:
            base = tuple(
                ((base_row.astype(np.int16) + block[:, 0]) % modulus).astype(np.uint8)
                for modulus, base_row, block in zip(moduli, base, contribution)
            )
        else:
            variable.append(contribution)

    consistent, join = count_join(moduli, base, tuple(variable))
    return {
        "experiment": "p7_fixed_boundary_catalog_join",
        "status": "complete_exact_multimodular_catalog_join",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "fixed_scaled_means": list(means),
        "moduli": list(moduli),
        "direction_catalogs": metadata,
        "variable_catalog_count": len(variable),
        "join": join,
        "consistent_catalog_tuples": int(consistent),
        "modularly_infeasible": consistent == 0,
        "finite_mean_allocation_exclusion": consistent == 0,
        "linear_system": linear_rows,
        "edge_equation_matrix_shape": list(matrix.shape),
        "catalog_cache": (
            {
                "path": str(cache[0]),
                "summary": str(cache[1]),
                "sha256": file_sha256(cache[0]),
            }
            if cache is not None
            else None
        ),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs="+", required=True)
    parser.add_argument("--means", type=int, nargs=8, required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 7))
    parser.add_argument("--catalog-cache", type=Path)
    parser.add_argument("--catalog-cache-summary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.catalog_cache is None) != (args.catalog_cache_summary is None):
        raise ValueError("catalog cache and summary must be supplied together")
    cache = (
        (args.catalog_cache, args.catalog_cache_summary)
        if args.catalog_cache is not None
        else None
    )
    out = run(
        args.c_h,
        tuple(args.fixed_boundary),
        tuple(args.means),
        tuple(args.moduli),
        cache,
    )
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
