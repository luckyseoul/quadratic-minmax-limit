#!/usr/bin/env python3
"""Build exact omission tables for every post-15.664 raised support.

The remaining p=7 size-eight leaves raise between one and five direction
catalogs.  For each of the 59 possible raised supports, this constructs the
subspace of modular dependencies vanishing on every raised direction and
materializes all remaining floor-catalog contributions.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import scaled_direction_floor  # noqa: E402
from p7_exceptional_omit_high_catalogs import (  # noqa: E402
    IncrementalRowBasis,
    modular_rank,
    modular_right_nullspace,
)
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402


PROJECTION_ROWS = 40
MAX_FLOOR_CATALOG_ROWS = 36


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def atomic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.stem}.{os.getpid()}.{time.time_ns()}.tmp.npz"
    )
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def valid_support_masks(epsilons: tuple[int, ...]) -> tuple[int, ...]:
    negative = tuple(index for index, eps in enumerate(epsilons) if eps == -1)
    positive = tuple(index for index, eps in enumerate(epsilons) if eps == 1)
    if len(negative) != 4 or len(positive) != 4:
        raise AssertionError("quadratic direction types changed")
    supports = set()
    for size in range(1, 5):
        for subset in itertools.combinations(positive, size):
            supports.add(sum(1 << index for index in subset))
    for negative_index in negative:
        for size in (1, 2, 4):
            for subset in itertools.combinations(positive, size):
                supports.add((1 << negative_index) | sum(1 << index for index in subset))
    output = tuple(sorted(supports))
    if len(output) != 59:
        raise AssertionError(f"expected 59 raised supports, found {len(output)}")
    return output


def choose_balanced_projection_rows(
    conditioned: np.ndarray, support_mask: int, modulus: int
) -> tuple[np.ndarray, list[int]]:
    """Greedily maximize marginal rank on every retained direction block."""
    blocks = [
        conditioned[:, 2 + 35 * direction : 2 + 35 * (direction + 1)]
        for direction in range(8)
    ]
    retained = [
        direction for direction in range(8) if not support_mask & (1 << direction)
    ]
    bases = [IncrementalRowBasis(modulus) for _ in range(8)]
    selected: list[int] = []
    remaining = set(range(len(conditioned)))
    while len(selected) < PROJECTION_ROWS:
        best = None
        for index in remaining:
            gains = {
                direction: int(bases[direction].gains_rank(blocks[direction][index]))
                for direction in retained
            }
            support = sum(
                int(np.count_nonzero(blocks[direction][index]))
                for direction in retained
            )
            score = (
                sum(gains.values()),
                sum((35 - len(bases[direction].rows)) * gains[direction] for direction in retained),
                support,
                -index,
            )
            if best is None or score > best[0]:
                best = (score, index)
        if best is None:
            raise AssertionError("balanced projection selection failed")
        index = best[1]
        selected.append(index)
        remaining.remove(index)
        for direction in retained:
            bases[direction].add(blocks[direction][index])
    selected_array = np.asarray(selected, dtype=np.int64)
    if modular_rank(conditioned[selected_array], modulus) != PROJECTION_ROWS:
        raise AssertionError("balanced selected rows lost global rank")
    return selected_array, [len(basis.rows) for basis in bases]


def floor_for(eps: int, b: int) -> int:
    return int(scaled_direction_floor(7, b, int(eps == -1)))


def contribution(
    projected: np.ndarray,
    direction_index: int,
    mask: int,
    eps: int,
    modulus: int,
) -> np.ndarray:
    b = mask.bit_count()
    phase = int(eps == -1)
    mean = floor_for(eps, b)
    odd_fibres = {value for value in range(7) if mask & (1 << value)}
    values = mapped_catalog(b, phase, mean, odd_fibres, None).astype(np.int64)
    bad = 13 - values
    block = projected[
        :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
    ].astype(np.int64)
    return (block @ (bad.T % modulus) % modulus).astype(np.uint8)


def run(output_path: Path, summary_path: Path, modulus: int) -> dict:
    started = time.time()
    matrix, dependencies, linear_rows = linear_data((modulus,))
    dependency = dependencies[modulus].astype(np.uint8)
    dependency_dimension = dependency.shape[0]
    if matrix.shape != (282, 1225) or dependency.shape[1] != 282:
        raise AssertionError("common score system changed")
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    support_masks = valid_support_masks(epsilons)
    support_index = np.full(256, -1, dtype=np.int16)
    for index, mask in enumerate(support_masks):
        support_index[mask] = index

    support_count = len(support_masks)
    selected_coefficients = np.zeros(
        (support_count, PROJECTION_ROWS, dependency_dimension), dtype=np.uint8
    )
    projected_dependencies = np.zeros(
        (support_count, PROJECTION_ROWS, 282), dtype=np.uint8
    )
    conditioned_dimensions = np.zeros(support_count, dtype=np.uint8)
    projection_records = []
    for support_index_value, support_mask in enumerate(support_masks):
        directions = tuple(index for index in range(8) if support_mask & (1 << index))
        columns = np.asarray(
            [
                column
                for direction in directions
                for column in range(2 + 35 * direction, 2 + 35 * (direction + 1))
            ],
            dtype=np.int64,
        )
        coefficients, block_rank = modular_right_nullspace(
            dependency[:, columns].astype(np.int64).T, modulus
        )
        conditioned = coefficients @ dependency.astype(np.int64) % modulus
        dimension = dependency_dimension - block_rank
        if (
            coefficients.shape != (dimension, dependency_dimension)
            or conditioned.shape != (dimension, 282)
            or dimension < PROJECTION_ROWS
            or np.any(conditioned[:, columns])
            or modular_rank(coefficients, modulus) != dimension
            or modular_rank(conditioned, modulus) != dimension
        ):
            raise AssertionError(f"conditioned space failed for support {support_mask}")
        selected, selected_block_ranks = choose_balanced_projection_rows(
            conditioned, support_mask, modulus
        )
        selected_coefficients[support_index_value] = coefficients[selected] % modulus
        projected_dependencies[support_index_value] = conditioned[selected] % modulus
        conditioned_dimensions[support_index_value] = dimension
        if (
            modular_rank(selected_coefficients[support_index_value], modulus)
            != PROJECTION_ROWS
            or modular_rank(projected_dependencies[support_index_value], modulus)
            != PROJECTION_ROWS
            or np.any(projected_dependencies[support_index_value, :, columns])
        ):
            raise AssertionError(f"selected projection failed for support {support_mask}")
        projection_records.append(
            {
                "support_mask": support_mask,
                "directions": list(directions),
                "type_support": [
                    sum(epsilons[index] == -1 for index in directions),
                    sum(epsilons[index] == 1 for index in directions),
                ],
                "conditioned_dimension": dimension,
                "selected_rank": PROJECTION_ROWS,
                "selected_direction_block_ranks": selected_block_ranks,
            }
        )

    base = np.zeros((support_count, PROJECTION_ROWS), dtype=np.uint8)
    singleton = np.full(
        (support_count, 8, 128, PROJECTION_ROWS), 255, dtype=np.uint8
    )
    variable = np.zeros(
        (
            support_count,
            8,
            128,
            MAX_FLOOR_CATALOG_ROWS,
            PROJECTION_ROWS,
        ),
        dtype=np.uint8,
    )
    variable_count = np.zeros((support_count, 8, 128), dtype=np.uint8)
    catalog_histogram: Counter[int] = Counter()
    for support_index_value, support_mask in enumerate(support_masks):
        projected = projected_dependencies[support_index_value]
        base[support_index_value] = (
            projected[:, :2].astype(np.int64)
            @ np.asarray([29, 1], dtype=np.int64)
            % modulus
        ).astype(np.uint8)
        for direction, eps in enumerate(epsilons):
            for mask in range(128):
                if mask.bit_count() not in (0, 2, 4, 6):
                    continue
                if support_mask & (1 << direction):
                    singleton[support_index_value, direction, mask] = 0
                    continue
                values = contribution(projected, direction, mask, eps, modulus)
                rows = int(values.shape[1])
                catalog_histogram[rows] += 1
                if rows == 1:
                    singleton[support_index_value, direction, mask] = values[:, 0]
                elif rows == 36 and eps == -1 and mask.bit_count() == 4:
                    variable[support_index_value, direction, mask] = values.T
                    variable_count[support_index_value, direction, mask] = rows
                else:
                    raise AssertionError(
                        f"unexpected floor catalog support={support_mask} "
                        f"direction={direction} mask={mask} rows={rows}"
                    )

    arrays = {
        "labels": labels.astype(np.int8),
        "epsilons": epsilon_array.astype(np.int8),
        "dependency": dependency,
        "support_masks": np.asarray(support_masks, dtype=np.uint8),
        "support_index": support_index,
        "conditioned_dimensions": conditioned_dimensions,
        "selected_coefficients": selected_coefficients,
        "projected_dependencies": projected_dependencies,
        "base": base,
        "singleton": singleton,
        "variable": variable,
        "variable_count": variable_count,
    }
    atomic_npz(output_path, arrays)
    out = {
        "experiment": "p7_size8_multi_elevation_tables",
        "status": "complete_exact_post_15664_omission_tables",
        "p": 7,
        "c_H": -1,
        "modulus": modulus,
        "projection_rows": PROJECTION_ROWS,
        "valid_support_count": support_count,
        "minimum_conditioned_dimension": int(conditioned_dimensions.min()),
        "maximum_conditioned_dimension": int(conditioned_dimensions.max()),
        "linear_system": linear_rows,
        "projection_records": projection_records,
        "catalog_count_histogram": dict(sorted(catalog_histogram.items())),
        "array_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "array_sha256": {key: array_sha256(value) for key, value in arrays.items()},
        "output": str(output_path),
        "output_sha256": sha256(output_path),
        "output_bytes": output_path.stat().st_size,
        "elapsed_seconds": time.time() - started,
    }
    atomic_json(summary_path, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--modulus", type=int, choices=(3, 7), default=7)
    args = parser.parse_args()
    print(json.dumps(run(args.output, args.summary, args.modulus), indent=2), flush=True)


if __name__ == "__main__":
    main()
