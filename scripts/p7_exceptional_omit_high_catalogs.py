#!/usr/bin/env python3
"""Exact mod-seven projections that eliminate one high-mean direction.

Every exceptional high-mean leaf has exactly one ``b=0, phase=0`` direction
with scaled mean 20, 24, or 32.  Enumerating that catalog is unnecessary.
Starting from the 135 left dependencies of the common edge system, take the
left kernel of the selected direction block.  For p=7 this gives 112 exact
dependencies whose entire 35-column selected-direction block is zero.  The
large catalog therefore contributes the zero syndrome, whatever its values.

This script deterministically chooses 66 conditioned dependencies, splits
them into three disjoint 22-coordinate projections, and materializes only
the seven small directional catalogs needed by the corresponding leaves.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_exceptional_projected_catalogs import atomic_json, atomic_npz, sha256  # noqa: E402
from p7_exceptional_tail22_catalogs import store22  # noqa: E402
from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402


OMITTED_DIRECTIONS = (0, 2, 5, 7)
MODULUS = 7
PROFILE_SIZE = 22
PROFILE_COUNT = 3


def matrix_sha256(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(matrix).tobytes()).hexdigest()


def modular_rank(matrix: np.ndarray, modulus: int = MODULUS) -> int:
    work = np.asarray(matrix, dtype=np.int64).copy() % modulus
    if work.ndim != 2:
        raise ValueError("rank input must be a matrix")
    row = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = (
            work[row] * pow(int(work[row, column]), -1, modulus)
        ) % modulus
        for target in range(row + 1, work.shape[0]):
            if work[target, column]:
                work[target] = (
                    work[target] - work[target, column] * work[row]
                ) % modulus
        row += 1
        if row == work.shape[0]:
            break
    return row


def modular_right_nullspace(
    matrix: np.ndarray, modulus: int = MODULUS
) -> tuple[np.ndarray, int]:
    """Return a row basis of ``{x: matrix @ x = 0}`` over F_modulus."""
    work = np.asarray(matrix, dtype=np.int64).copy() % modulus
    rows, columns = work.shape
    pivots: list[int] = []
    row = 0
    for column in range(columns):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = (
            work[row] * pow(int(work[row, column]), -1, modulus)
        ) % modulus
        for target in range(rows):
            if target != row and work[target, column]:
                work[target] = (
                    work[target] - work[target, column] * work[row]
                ) % modulus
        pivots.append(column)
        row += 1
        if row == rows:
            break
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = np.zeros(columns, dtype=np.int64)
        vector[free_column] = 1
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -work[pivot_row, free_column] % modulus
        basis.append(vector)
    return np.asarray(basis, dtype=np.int64), row


class IncrementalRowBasis:
    def __init__(self, modulus: int = MODULUS):
        self.modulus = modulus
        self.rows: list[tuple[int, np.ndarray]] = []

    def reduce(self, vector: np.ndarray) -> np.ndarray:
        out = np.asarray(vector, dtype=np.int64).copy() % self.modulus
        for pivot, row in self.rows:
            if out[pivot]:
                out = (out - out[pivot] * row) % self.modulus
        return out

    def gains_rank(self, vector: np.ndarray) -> bool:
        return bool(np.any(self.reduce(vector)))

    def add(self, vector: np.ndarray) -> bool:
        out = self.reduce(vector)
        if not np.any(out):
            return False
        pivot = int(np.flatnonzero(out)[0])
        out = out * pow(int(out[pivot]), -1, self.modulus) % self.modulus
        updated = []
        for old_pivot, old_row in self.rows:
            if old_row[pivot]:
                old_row = (old_row - old_row[pivot] * out) % self.modulus
            updated.append((old_pivot, old_row))
        updated.append((pivot, out))
        updated.sort(key=lambda item: item[0])
        self.rows = updated
        return True


def conditioned_dependencies(omitted_direction: int) -> dict:
    if omitted_direction not in OMITTED_DIRECTIONS:
        raise ValueError(f"omitted direction must be one of {OMITTED_DIRECTIONS}")
    _matrix, dependencies, _rows = linear_data((MODULUS,))
    dependency = dependencies[MODULUS].astype(np.int64)
    omitted_columns = np.arange(
        2 + 35 * omitted_direction,
        2 + 35 * (omitted_direction + 1),
    )
    coefficients, block_rank = modular_right_nullspace(
        dependency[:, omitted_columns].T
    )
    conditioned = coefficients @ dependency % MODULUS
    if block_rank != 23 or coefficients.shape != (112, 135):
        raise AssertionError("conditioned dependency dimension changed")
    if np.any(conditioned[:, omitted_columns]):
        raise AssertionError("conditioned dependencies touch omitted direction")
    if modular_rank(coefficients) != 112 or modular_rank(conditioned) != 112:
        raise AssertionError("conditioned dependency basis lost rank")
    return {
        "dependency": dependency,
        "coefficients": coefficients,
        "conditioned": conditioned,
        "omitted_block_rank": block_rank,
    }


def choose_projection_rows(conditioned: np.ndarray, omitted_direction: int) -> list[int]:
    """Greedily balance marginal block rank over all seven retained directions."""
    blocks = [
        conditioned[:, 2 + 35 * index : 2 + 35 * (index + 1)]
        for index in range(8)
    ]
    bases = [IncrementalRowBasis() for _ in range(8)]
    selected: list[int] = []
    remaining = set(range(len(conditioned)))
    while len(selected) < PROFILE_SIZE * PROFILE_COUNT:
        best = None
        for index in remaining:
            gains = [
                0
                if direction == omitted_direction
                else int(bases[direction].gains_rank(blocks[direction][index]))
                for direction in range(8)
            ]
            support = sum(
                int(np.count_nonzero(blocks[direction][index]))
                for direction in range(8)
                if direction != omitted_direction
            )
            score = (
                sum(gains),
                sum(
                    (21 - len(bases[direction].rows)) * gains[direction]
                    for direction in range(8)
                ),
                support,
                -index,
            )
            if best is None or score > best[0]:
                best = (score, index)
        if best is None:
            raise AssertionError("projection-row selection failed")
        index = best[1]
        selected.append(index)
        remaining.remove(index)
        for direction in range(8):
            if direction != omitted_direction:
                bases[direction].add(blocks[direction][index])
    expected_ranks = [0 if index == omitted_direction else 14 for index in range(8)]
    actual_ranks = [len(basis.rows) for basis in bases]
    if actual_ranks != expected_ranks:
        raise AssertionError(
            f"balanced conditioned projection ranks changed: {actual_ranks}"
        )
    return selected


def high_leaves(batch: dict, omitted_direction: int) -> list[dict]:
    rows = batch["direction_rows"]
    selected = []
    for leaf in batch["leaves"]:
        if leaf.get("solver_status") == "INFEASIBLE":
            continue
        hits = [
            index
            for index, (row, mean) in enumerate(
                zip(rows, leaf["scaled_means_direction_order"])
            )
            if int(row["b"]) == 0
            and int(row["phase"]) == 0
            and int(mean) > 16
        ]
        if hits:
            if len(hits) != 1 or hits[0] not in OMITTED_DIRECTIONS:
                raise AssertionError("unexpected exceptional high-mean pattern")
            if hits[0] == omitted_direction:
                selected.append(leaf)
    return selected


def build_one(
    orbit_index: int,
    boundary: tuple[int, ...],
    mean_path_text: str,
    cache_path_text: str,
    cache_summary_text: str,
    output_dir_text: str,
    omitted_direction: int,
    profile: int,
    selected_full: tuple[int, ...],
    conditioned_text: bytes,
) -> dict:
    started = time.time()
    mean_path = Path(mean_path_text)
    batch = json.loads(mean_path.read_text())
    if (
        batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
        or tuple(batch.get("fixed_boundary", [])) != boundary
        or int(batch.get("allocation_count", -1)) != 180
    ):
        raise ValueError("exceptional mean batch has incompatible metadata")
    leaves = high_leaves(batch, omitted_direction)
    if not leaves:
        raise AssertionError("orbit has no leaves for requested omitted direction")

    conditioned = np.frombuffer(conditioned_text, dtype=np.uint8).reshape(112, 282)
    selected = selected_full[
        profile * PROFILE_SIZE : (profile + 1) * PROFILE_SIZE
    ]
    projected = conditioned[np.asarray(selected, dtype=np.int64)].astype(np.int64)
    arrays: dict[str, np.ndarray] = {}
    base_digits = projected[:, :2] @ np.asarray([29, 1], dtype=np.int64) % MODULUS
    store22(arrays, "base_p7", base_digits.reshape(PROFILE_SIZE, 1))

    cache = (Path(cache_path_text), Path(cache_summary_text))
    rows = direction_rows(-1, boundary)
    metadata = []
    for direction_index, row in enumerate(rows):
        means = sorted(
            {
                int(leaf["scaled_means_direction_order"][direction_index])
                for leaf in leaves
            }
        )
        block = projected[
            :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
        ]
        for mean in means:
            if direction_index == omitted_direction:
                if mean not in (20, 24, 32) or np.any(block):
                    raise AssertionError("invalid omitted high-mean catalog")
                digits = np.zeros((PROFILE_SIZE, 1), dtype=np.int64)
                catalog_rows = "eliminated_exactly_by_zero_dependency_block"
                source_rows = 1
            else:
                values = mapped_catalog(
                    int(row["b"]), int(row["phase"]), mean, set(row["B"]), cache
                )
                bad = 13 - values.astype(np.int64)
                digits = block @ (bad.T % MODULUS) % MODULUS
                catalog_rows = int(len(values))
                source_rows = int(len(values))
            stem = f"d{direction_index}_m{mean}_p7"
            store22(arrays, stem, digits)
            signatures = set(
                zip(
                    arrays[f"{stem}_lo"].tolist(),
                    arrays[f"{stem}_hi"].tolist(),
                )
            )
            metadata.append(
                {
                    "direction_index": direction_index,
                    "scaled_mean": mean,
                    "catalog_rows": catalog_rows,
                    "stored_aligned_rows": source_rows,
                    "distinct_projected_signatures": len(signatures),
                }
            )

    output = Path(output_dir_text) / (
        f"cminus_exceptional_omitd{omitted_direction}_p{profile}_"
        f"orbit{orbit_index:02d}.npz"
    )
    atomic_npz(output, arrays)
    return {
        "orbit_index": orbit_index,
        "fixed_boundary": list(boundary),
        "mean_source": str(mean_path),
        "high_leaf_indices": [int(leaf["leaf_index"]) for leaf in leaves],
        "high_leaf_count": len(leaves),
        "output": str(output),
        "sha256": sha256(output),
        "bytes": output.stat().st_size,
        "catalogs": metadata,
        "elapsed_seconds": time.time() - started,
    }


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    source = json.loads(args.source.read_text())
    orbits = exceptional_orbits(source)
    conditioned_data = conditioned_dependencies(args.omitted_direction)
    conditioned = conditioned_data["conditioned"].astype(np.uint8)
    coefficients = conditioned_data["coefficients"].astype(np.uint8)
    selected_full = choose_projection_rows(conditioned, args.omitted_direction)
    selected = selected_full[
        args.profile * PROFILE_SIZE : (args.profile + 1) * PROFILE_SIZE
    ]
    selected_rows = conditioned[np.asarray(selected)]
    selected_coefficients = coefficients[np.asarray(selected)]
    if modular_rank(selected_rows) != PROFILE_SIZE:
        raise AssertionError("selected conditioned rows are not independent")

    coverage = []
    profile_ranks = []
    full_ranks = []
    for direction_index in range(8):
        columns = np.arange(
            2 + 35 * direction_index,
            2 + 35 * (direction_index + 1),
        )
        block = selected_rows[:, columns]
        full_block = conditioned[np.asarray(selected_full)][:, columns]
        coverage.append(int(np.sum(np.any(block, axis=1))))
        profile_ranks.append(modular_rank(block))
        full_ranks.append(modular_rank(full_block))
    expected_full = [
        0 if index == args.omitted_direction else 14 for index in range(8)
    ]
    if full_ranks != expected_full:
        raise AssertionError("full 66-row block-rank audit failed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    task_rows = [
        (
            orbit_index,
            boundary,
            str(
                args.mean_dir
                / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"
            ),
            str(args.catalog_cache),
            str(args.catalog_cache_summary),
            str(args.output_dir),
            args.omitted_direction,
            args.profile,
            tuple(selected_full),
            conditioned.tobytes(),
        )
        for orbit_index, boundary in orbits
    ]
    rows = []

    def record(row: dict) -> None:
        rows.append(row)
        print(
            json.dumps({key: value for key, value in row.items() if key != "catalogs"}),
            flush=True,
        )

    if args.workers == 1:
        for task in task_rows:
            record(build_one(*task))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(build_one, *task) for task in task_rows]
            for future in as_completed(futures):
                record(future.result())
    rows.sort(key=lambda row: row["orbit_index"])
    family_rows = conditioned[np.asarray(selected_full)]
    return {
        "experiment": "p7_exceptional_omit_high_catalogs",
        "status": "complete_exact_high_direction_eliminating_projection",
        "p": 7,
        "c_H": -1,
        "modulus": MODULUS,
        "omitted_direction": args.omitted_direction,
        "omitted_direction_block_rank": conditioned_data["omitted_block_rank"],
        "conditioned_dependency_dimension": len(conditioned),
        "conditioned_dependency_block_is_zero": True,
        "conditioned_dependency_basis_rank": modular_rank(conditioned),
        "selection_algorithm": "greedy_balanced_marginal_direction_block_rank_v1",
        "selected_full_conditioned_basis_rows": selected_full,
        "selected_full_dependency_sha256": matrix_sha256(family_rows),
        "selected_full_block_ranks": full_ranks,
        "projection_profile": args.profile,
        "selected_conditioned_basis_rows": selected,
        "selected_dependency_coefficients_mod7": selected_coefficients.tolist(),
        "selected_dependency_rows_sha256": matrix_sha256(selected_rows),
        "selected_row_coverage_by_direction": coverage,
        "selected_row_block_ranks": profile_ranks,
        "projection_group_order": 7**PROFILE_SIZE,
        "projection_encoding": "two_packed_nibble_words_to_injective_uint64_base7",
        "large_catalog_treatment": "exactly_eliminated_by_zero_dependency_block",
        "projection_is_necessary_not_sufficient": True,
        "orbit_count": len(rows),
        "total_high_leaves": sum(row["high_leaf_count"] for row in rows),
        "elapsed_seconds": time.time() - started,
        "orbits": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mean-dir", type=Path, required=True)
    parser.add_argument("--catalog-cache", type=Path, required=True)
    parser.add_argument("--catalog-cache-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--omitted-direction", type=int, choices=OMITTED_DIRECTIONS, required=True
    )
    parser.add_argument("--profile", type=int, choices=range(PROFILE_COUNT), required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(args)
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(
        json.dumps({key: value for key, value in out.items() if key != "orbits"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
