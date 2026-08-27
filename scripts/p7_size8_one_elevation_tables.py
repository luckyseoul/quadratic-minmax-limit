#!/usr/bin/env python3
"""Build exact omission tables for the p=7 size-eight four-allocation stratum.

For each possible elevated direction, take the 112-dimensional subspace of
the common mod-seven left dependencies that is identically zero on that
direction's 35 score columns.  Select 22 deterministic dependencies and
materialize every remaining floor-catalog contribution.  The resulting
cache lets a CUDA scan test a mean allocation without enumerating its
1,764- or 2,233-row elevated catalog.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import scaled_direction_floor  # noqa: E402
from p7_exceptional_omit_high_catalogs import (  # noqa: E402
    choose_projection_rows,
    modular_rank,
    modular_right_nullspace,
)
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402


MODULUS = 7
PROJECTION_ROWS = 22
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


def floor_for(eps: int, b: int) -> int:
    return int(scaled_direction_floor(7, b, int(eps == -1)))


def contribution(
    projected: np.ndarray,
    direction_index: int,
    mask: int,
    eps: int,
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
    return (block @ (bad.T % MODULUS) % MODULUS).astype(np.uint8)


def run(output: Path, summary: Path) -> dict:
    started = time.time()
    matrix, dependencies, linear_rows = linear_data((MODULUS,))
    dependency = dependencies[MODULUS].astype(np.uint8)
    if (
        matrix.shape != (282, 1225)
        or dependency.shape != (135, 282)
        or linear_rows
        != [
            {
                "modulus": 7,
                "rank": 147,
                "left_dependency_dimension": 135,
                "left_null_audit": True,
            }
        ]
    ):
        raise AssertionError("common mod-seven score system changed")
    if np.any(dependency.astype(np.int64) @ (matrix.astype(np.int64) % 7) % 7):
        raise AssertionError("full dependency basis is not left-null")

    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    coefficients = np.zeros((8, PROJECTION_ROWS, 135), dtype=np.uint8)
    projected = np.zeros((8, PROJECTION_ROWS, 282), dtype=np.uint8)
    projection_rows = []
    for omitted_direction in range(8):
        columns = np.arange(
            2 + 35 * omitted_direction,
            2 + 35 * (omitted_direction + 1),
        )
        conditioned_coefficients, block_rank = modular_right_nullspace(
            dependency[:, columns].astype(np.int64).T
        )
        conditioned = (
            conditioned_coefficients @ dependency.astype(np.int64) % MODULUS
        )
        if (
            block_rank != 23
            or conditioned_coefficients.shape != (112, 135)
            or conditioned.shape != (112, 282)
            or np.any(conditioned[:, columns])
            or modular_rank(conditioned_coefficients) != 112
            or modular_rank(conditioned) != 112
        ):
            raise AssertionError(
                f"conditioned dependency space changed for direction {omitted_direction}"
            )
        selected_full = choose_projection_rows(conditioned, omitted_direction)
        selected_indices = np.asarray(selected_full[:PROJECTION_ROWS], dtype=np.int64)
        selected_coefficients = conditioned_coefficients[selected_indices] % MODULUS
        selected_rows = conditioned[selected_indices] % MODULUS
        if (
            modular_rank(selected_coefficients) != PROJECTION_ROWS
            or modular_rank(selected_rows) != PROJECTION_ROWS
            or np.any(selected_rows[:, columns])
        ):
            raise AssertionError("selected omission projection lost rank")
        coefficients[omitted_direction] = selected_coefficients.astype(np.uint8)
        projected[omitted_direction] = selected_rows.astype(np.uint8)
        block_ranks = [
            modular_rank(
                selected_rows[
                    :, 2 + 35 * direction : 2 + 35 * (direction + 1)
                ]
            )
            for direction in range(8)
        ]
        projection_rows.append(
            {
                "omitted_direction": omitted_direction,
                "omitted_block_rank": block_rank,
                "conditioned_dependency_dimension": 112,
                "selected_conditioned_basis_rows": [
                    int(value) for value in selected_indices
                ],
                "selected_rank_mod7": modular_rank(selected_rows),
                "direction_block_ranks_mod7": block_ranks,
                "selected_coefficients_sha256": array_sha256(
                    selected_coefficients.astype(np.uint8)
                ),
                "selected_dependencies_sha256": array_sha256(
                    selected_rows.astype(np.uint8)
                ),
            }
        )

    base = np.zeros((8, PROJECTION_ROWS), dtype=np.uint8)
    singleton = np.full(
        (8, 8, 128, PROJECTION_ROWS), 255, dtype=np.uint8
    )
    variable = np.zeros(
        (8, 8, 128, MAX_FLOOR_CATALOG_ROWS, PROJECTION_ROWS), dtype=np.uint8
    )
    variable_count = np.zeros((8, 8, 128), dtype=np.uint8)
    catalog_histogram: Counter[int] = Counter()
    for omitted_direction in range(8):
        base[omitted_direction] = (
            projected[omitted_direction, :, :2].astype(np.int64)
            @ np.asarray([29, 1], dtype=np.int64)
            % MODULUS
        ).astype(np.uint8)
        for direction_index, eps in enumerate(epsilons):
            for mask in range(128):
                if mask.bit_count() not in (0, 2, 4, 6):
                    continue
                values = contribution(
                    projected[omitted_direction], direction_index, mask, eps
                )
                rows = int(values.shape[1])
                catalog_histogram[rows] += 1
                if rows == 1:
                    singleton[omitted_direction, direction_index, mask] = values[:, 0]
                elif (
                    rows == MAX_FLOOR_CATALOG_ROWS
                    and eps == -1
                    and mask.bit_count() == 4
                ):
                    variable[omitted_direction, direction_index, mask] = values.T
                    variable_count[omitted_direction, direction_index, mask] = rows
                else:
                    raise AssertionError(
                        "unexpected floor catalog in omission table: "
                        f"omit={omitted_direction} d={direction_index} "
                        f"mask={mask} rows={rows}"
                    )

    arrays = {
        "labels": labels.astype(np.int8),
        "epsilons": epsilon_array.astype(np.int8),
        "dependency": dependency,
        "selected_coefficients": coefficients,
        "projected_dependencies": projected,
        "base": base,
        "singleton": singleton,
        "variable": variable,
        "variable_count": variable_count,
    }
    atomic_npz(output, arrays)
    out = {
        "experiment": "p7_size8_one_elevation_tables",
        "status": "complete_exact_elevated_direction_omission_tables",
        "p": 7,
        "c_H": -1,
        "modulus": MODULUS,
        "projection_rows": PROJECTION_ROWS,
        "linear_system": linear_rows,
        "dependency_shape": list(dependency.shape),
        "dependency_sha256": array_sha256(dependency),
        "projection_rows_by_omitted_direction": projection_rows,
        "catalog_count_histogram": dict(sorted(catalog_histogram.items())),
        "array_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "array_sha256": {key: array_sha256(value) for key, value in arrays.items()},
        "output": str(output),
        "output_sha256": sha256(output),
        "output_bytes": output.stat().st_size,
        "elapsed_seconds": time.time() - started,
    }
    atomic_json(summary, out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.output, args.summary)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
