#!/usr/bin/env python3
"""Exact CUDA joins under full-rank linear maps F7^21 -> F7^k.

This is a linear-projection companion to
``p7_infinity7_positive_z7_semigroup_case_join_gpu.py`` for the 51 positive
``p=7, z=7`` grade-three cases.  The established CPU/GPU module remains the
source of the reconstructed cases, Hilbert basis, exact torsion quotient,
direct catalog calibration, CUDA collision-resolved set convolution, cap
semantics, and manufactured audits.

Every projection retains all six F3 quotient coordinates verbatim and maps
the 21 F7 quotient coordinates by one explicit full-row-rank matrix
``M in F7^(k x 21)``.  The identical rowwise map is applied to Hilbert
generators, direct-catalog rows, and the projected target.  Thus mod-3 and
mod-7 values always come from the same source row before deduplication.

The deterministic matrix catalog includes seeded dense maps, block sums,
Vandermonde and bivariate evaluation maps, a pointed-system geometry map,
and a coordinate-selector baseline.  Every matrix is emitted explicitly,
hash-pinned, and accompanied by an exact nonzero k-minor rank certificate.
Cases and matrices can be independently sharded.  A completed convolution
with a missing target is a rigorous rejection; target presence is necessary
only.  State, memory, allocation, and hash-cap failures are explicit skips.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join_gpu as coordinate_gpu  # noqa: E402


cpu_join = coordinate_gpu.cpu_join
cuda_join = coordinate_gpu.cuda_join

EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_join_linear_gpu"
SUPPORTED_K = coordinate_gpu.SUPPORTED_K
EXPECTED_TARGET_CASES = coordinate_gpu.EXPECTED_TARGET_CASES
EXACT_DIRECTION_SUPPORT_CAP = coordinate_gpu.EXACT_DIRECTION_SUPPORT_CAP
DEFAULT_PAIR_CHUNK_CAP = coordinate_gpu.DEFAULT_PAIR_CHUNK_CAP
DEFAULT_GPU_MEMORY_CAP_MIB = coordinate_gpu.DEFAULT_GPU_MEMORY_CAP_MIB
ALL_MOD7_COORDINATES = tuple(range(21))
SELECTOR_CROSSCHECK_COORDINATES = coordinate_gpu.DEFAULT_CROSSCHECK_PROJECTION
DEFAULT_MATRIX_FAMILIES = (
    "seeded_dense",
    "geometry",
    "evaluation",
    "vandermonde",
    "block_sums",
    "selector",
)
DEFAULT_DENSE_SEEDS = (
    "qml-p7-z7-linear-v1-dense-0",
    "qml-p7-z7-linear-v1-dense-1",
    "qml-p7-z7-linear-v1-dense-2",
)
MATRIX_FAMILY_ORDER = {name: index for index, name in enumerate(DEFAULT_MATRIX_FAMILIES)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def exact_rank_pivots(matrix: np.ndarray, modulus: int = 7) -> tuple[int, tuple[int, ...]]:
    """Return exact row rank and pivot columns by modular Gauss-Jordan."""
    source = np.ascontiguousarray(matrix, dtype=np.int64) % modulus
    require(source.ndim == 2, "rank input is not a matrix")
    reduced = source.copy()
    row = 0
    pivots = []
    for column in range(reduced.shape[1]):
        candidates = np.flatnonzero(reduced[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            reduced[[row, pivot]] = reduced[[pivot, row]]
        inverse = pow(int(reduced[row, column]), -1, modulus)
        reduced[row] = reduced[row] * inverse % modulus
        for other in range(reduced.shape[0]):
            if other == row or reduced[other, column] == 0:
                continue
            reduced[other] = (
                reduced[other] - reduced[other, column] * reduced[row]
            ) % modulus
        pivots.append(column)
        row += 1
        if row == reduced.shape[0]:
            break
    return row, tuple(pivots)


def determinant_mod_prime(matrix: np.ndarray, modulus: int = 7) -> int:
    """Compute a square determinant exactly in the prime field."""
    source = np.ascontiguousarray(matrix, dtype=np.int64) % modulus
    require(
        source.ndim == 2 and source.shape[0] == source.shape[1],
        "determinant input is not square",
    )
    work = source.copy()
    determinant = 1
    for column in range(len(work)):
        candidates = np.flatnonzero(work[column:, column])
        if not len(candidates):
            return 0
        pivot = column + int(candidates[0])
        if pivot != column:
            work[[column, pivot]] = work[[pivot, column]]
            determinant = -determinant
        pivot_value = int(work[column, column])
        determinant = determinant * pivot_value % modulus
        inverse = pow(pivot_value, -1, modulus)
        for row in range(column + 1, len(work)):
            factor = int(work[row, column]) * inverse % modulus
            work[row, column:] = (
                work[row, column:] - factor * work[column, column:]
            ) % modulus
    return determinant % modulus


def matrix_record(
    *,
    name: str,
    family: str,
    matrix: np.ndarray,
    construction: dict,
) -> dict:
    """Canonicalize and certify one explicit F7 projection matrix."""
    source = np.ascontiguousarray(matrix, dtype=np.int64) % 7
    require(
        source.ndim == 2 and source.shape[1] == 21,
        "linear projection matrix must have 21 columns",
    )
    k = source.shape[0]
    require(k in SUPPORTED_K or k == len(SELECTOR_CROSSCHECK_COORDINATES), "unsupported matrix row count")
    rank, pivots = exact_rank_pivots(source, 7)
    imported_rank = cpu_join.affine.modular_rank(source, 7)
    require(rank == imported_rank == k, f"matrix {name} is not full row rank")
    require(len(pivots) == k, "full-rank matrix has wrong pivot census")
    minor = source[:, list(pivots)]
    determinant = determinant_mod_prime(minor, 7)
    require(determinant != 0, "rank-certificate pivot minor is singular")
    matrix_u8 = np.ascontiguousarray(source, dtype=np.uint8)
    explicit = matrix_u8.tolist()
    record = {
        "name": name,
        "family": family,
        "domain": "F7^21",
        "codomain": f"F7^{k}",
        "shape": [k, 21],
        "matrix_mod7": explicit,
        "matrix_sha256_uint8": cpu_join.array_sha256(matrix_u8),
        "matrix_canonical_json_sha256": json_sha256(explicit),
        "rank_mod7": rank,
        "pivot_columns_zero_based": list(pivots),
        "pivot_minor_determinant_mod7": determinant,
        "exact_full_row_rank_certificate": True,
        "construction": construction,
    }
    record["record_sha256"] = json_sha256(record)
    return record


def selector_matrix(k: int, coordinates: tuple[int, ...] | None = None) -> np.ndarray:
    selected = tuple(range(k)) if coordinates is None else tuple(coordinates)
    require(len(selected) == k, "selector coordinate count changed")
    require(
        len(set(selected)) == k and all(0 <= value < 21 for value in selected),
        "selector coordinates are invalid",
    )
    matrix = np.zeros((k, 21), dtype=np.uint8)
    matrix[np.arange(k), np.asarray(selected, dtype=np.int64)] = 1
    return matrix


def deterministic_field_stream(label: str, count: int) -> np.ndarray:
    """Version-stable SHA-256 stream with unbiased residues modulo seven."""
    values = []
    counter = 0
    while len(values) < count:
        block = hashlib.sha256(f"{label}|{counter}".encode("utf-8")).digest()
        counter += 1
        for value in block:
            if value < 252:
                values.append(value % 7)
                if len(values) == count:
                    break
    return np.ascontiguousarray(values, dtype=np.uint8)


def seeded_dense_matrix(k: int, seed: str) -> tuple[np.ndarray, dict]:
    """Generate the first full-rank SHA-stream matrix for a named seed."""
    for attempt in range(256):
        label = f"{seed}|k={k}|attempt={attempt}"
        matrix = deterministic_field_stream(label, k * 21).reshape(k, 21)
        rank, _pivots = exact_rank_pivots(matrix, 7)
        if rank == k:
            return matrix, {
                "strategy": "seeded_dense_SHA256_rejection_stream",
                "seed": seed,
                "attempt": attempt,
                "byte_acceptance_rule": "accept byte<252 then reduce modulo7",
                "counter_encoding": "UTF8(seed|k=<k>|attempt=<a>)|decimal_counter",
            }
    raise AssertionError(f"seeded dense generator did not find rank {k}: {seed}")


def contiguous_block_sum_matrix(k: int) -> tuple[np.ndarray, dict]:
    matrix = np.zeros((k, 21), dtype=np.uint8)
    base, remainder = divmod(21, k)
    blocks = []
    start = 0
    for row in range(k):
        stop = start + base + (1 if row < remainder else 0)
        matrix[row, start:stop] = 1
        blocks.append(list(range(start, stop)))
        start = stop
    require(start == 21, "contiguous blocks do not cover the domain")
    return matrix, {
        "strategy": "disjoint_contiguous_block_sums",
        "blocks": blocks,
        "blocks_form_partition_of_0_through_20": True,
    }


def interleaved_block_sum_matrix(k: int) -> tuple[np.ndarray, dict]:
    matrix = np.zeros((k, 21), dtype=np.uint8)
    blocks = []
    for row in range(k):
        block = list(range(row, 21, k))
        matrix[row, block] = 1
        blocks.append(block)
    return matrix, {
        "strategy": "disjoint_interleaved_block_sums",
        "blocks": blocks,
        "block_rule": "coordinate_mod_k_equals_row",
        "blocks_form_partition_of_0_through_20": True,
    }


def vandermonde_matrix(k: int) -> tuple[np.ndarray, dict]:
    evaluation_points = [column % 7 for column in range(21)]
    matrix = np.asarray(
        [
            [pow(value, exponent, 7) for value in evaluation_points]
            for exponent in range(k)
        ],
        dtype=np.uint8,
    )
    return matrix, {
        "strategy": "univariate_evaluation_Vandermonde",
        "column_evaluation_points_mod7": evaluation_points,
        "row_exponents": list(range(k)),
        "first_seven_columns_form_full_rank_Vandermonde_for_k_at_most_7": True,
    }


def evaluation_grid_matrix(k: int) -> tuple[np.ndarray, dict]:
    """Evaluate independent low-degree monomials on F7 x {0,1,2}."""
    points = [(column % 7, column // 7) for column in range(21)]
    candidates = sorted(
        ((x_power, y_power) for x_power in range(7) for y_power in range(3)),
        key=lambda row: (sum(row), row[1], row[0]),
    )
    selected_exponents = []
    selected_rows = []
    current_rank = 0
    for x_power, y_power in candidates:
        row = np.asarray(
            [
                pow(x, x_power, 7) * pow(y, y_power, 7) % 7
                for x, y in points
            ],
            dtype=np.uint8,
        )
        trial = np.stack((*selected_rows, row)) if selected_rows else row[None, :]
        rank, _pivots = exact_rank_pivots(trial, 7)
        if rank > current_rank:
            selected_rows.append(row)
            selected_exponents.append((x_power, y_power))
            current_rank = rank
        if current_rank == k:
            break
    require(current_rank == k, "evaluation-grid monomials did not reach full rank")
    return np.ascontiguousarray(np.stack(selected_rows), dtype=np.uint8), {
        "strategy": "bivariate_monomial_evaluation_grid",
        "column_points_x_mod7_y_in_0_1_2": [list(row) for row in points],
        "selected_monomial_exponents_x_y": [list(row) for row in selected_exponents],
        "selection_rule": "graded order; retain exactly when modular row rank increases",
    }


def greedy_independent_rows(
    candidates: Iterable[tuple[str, np.ndarray]], k: int
) -> tuple[np.ndarray, list[str]]:
    selected = []
    labels = []
    rank = 0
    for label, candidate in candidates:
        row = np.ascontiguousarray(candidate, dtype=np.int64) % 7
        require(row.shape == (21,), "geometry candidate has wrong width")
        trial = np.stack((*selected, row)) if selected else row[None, :]
        trial_rank, _pivots = exact_rank_pivots(trial, 7)
        if trial_rank > rank:
            selected.append(np.ascontiguousarray(row, dtype=np.uint8))
            labels.append(label)
            rank = trial_rank
        if rank == k:
            break
    require(rank == k, "geometry candidates did not span the requested codomain")
    return np.ascontiguousarray(np.stack(selected), dtype=np.uint8), labels


def geometry_profile_matrix(context: dict[str, Any], k: int) -> tuple[np.ndarray, dict]:
    """Select independent F7 quotient profiles of Johnson-geometry features."""
    complement = np.ascontiguousarray(
        context["quotient_data"][7]["complement"], dtype=np.int64
    ) % 7
    require(complement.shape[0] == 21, "F7 quotient complement dimension changed")
    directional_equation_end = 1 + 8 * cpu_join.AMBIENT_DIMENSION
    require(
        complement.shape[1] == directional_equation_end + 1,
        "pointed-system equation-block width changed",
    )
    points = tuple(tuple(int(value) for value in point) for point in cpu_join.affine.POINTS)
    require(len(points) == cpu_join.AMBIENT_DIMENSION == 35, "Johnson point census changed")

    candidates: list[tuple[str, np.ndarray]] = []

    def profile(columns: list[int]) -> np.ndarray:
        return np.ascontiguousarray(complement[:, columns].sum(axis=1) % 7)

    for direction in range(8):
        block_start = 1 + cpu_join.AMBIENT_DIMENSION * direction
        all_columns = [block_start + index for index in range(35)]
        candidates.append((f"direction_{direction}_all_points", profile(all_columns)))

    for direction in range(8):
        block_start = 1 + cpu_join.AMBIENT_DIMENSION * direction
        for vertex in range(7):
            columns = [
                block_start + index
                for index, point in enumerate(points)
                if vertex in point
            ]
            candidates.append(
                (f"direction_{direction}_vertex_{vertex}_star", profile(columns))
            )

    for direction in range(8):
        block_start = 1 + cpu_join.AMBIENT_DIMENSION * direction
        for left in range(7):
            for right in range(left + 1, 7):
                columns = [
                    block_start + index
                    for index, point in enumerate(points)
                    if left in point and right in point
                ]
                candidates.append(
                    (
                        f"direction_{direction}_pair_{left}_{right}_star",
                        profile(columns),
                    )
                )

    # Column profiles must span F7^21 because the 21 complement rows are
    # independent.  They provide a deterministic geometry-labelled fallback.
    candidates.append(("constant_equation_profile", complement[:, 0]))
    for direction in range(8):
        block_start = 1 + cpu_join.AMBIENT_DIMENSION * direction
        for point_index, point in enumerate(points):
            candidates.append(
                (
                    f"direction_{direction}_point_{'_'.join(map(str, point))}",
                    complement[:, block_start + point_index],
                )
            )
    for extra_index, column in enumerate(
        range(directional_equation_end, complement.shape[1])
    ):
        candidates.append(
            (f"pointed_normalization_equation_{extra_index}_profile", complement[:, column])
        )

    matrix, labels = greedy_independent_rows(candidates, k)
    return matrix, {
        "strategy": "pointed_system_Johnson_geometry_quotient_profiles",
        "selected_feature_labels": labels,
        "candidate_order": (
            "direction totals; direction/vertex stars; direction/pair stars; "
            "constant, individual direction/4-subset, and pointed-normalization "
            "equation profiles"
        ),
        "F7_complement_sha256_int64": cpu_join.array_sha256(complement),
        "individual_equation_profile_fallback_is_complete_because_complement_has_row_rank_21": True,
    }


def parse_matrix_families(value: str) -> tuple[str, ...]:
    requested = tuple(item.strip() for item in value.split(",") if item.strip())
    require(requested, "matrix-family specification is empty")
    require(len(requested) == len(set(requested)), "matrix family repeated")
    unknown = set(requested) - set(DEFAULT_MATRIX_FAMILIES)
    require(not unknown, f"unknown matrix families: {sorted(unknown)}")
    return tuple(sorted(requested, key=MATRIX_FAMILY_ORDER.__getitem__))


def build_matrix_catalog(
    context: dict[str, Any], k: int, families: tuple[str, ...]
) -> tuple[list[dict], dict]:
    """Construct, rank-certify, explicitly emit, and deduplicate all maps."""
    records = []

    if "seeded_dense" in families:
        for seed_index, seed in enumerate(DEFAULT_DENSE_SEEDS):
            matrix, construction = seeded_dense_matrix(k, seed)
            records.append(
                matrix_record(
                    name=f"dense_sha256_{seed_index}",
                    family="seeded_dense",
                    matrix=matrix,
                    construction=construction,
                )
            )

    if "geometry" in families:
        matrix, construction = geometry_profile_matrix(context, k)
        records.append(
            matrix_record(
                name="geometry_star_profiles",
                family="geometry",
                matrix=matrix,
                construction=construction,
            )
        )

    if "evaluation" in families:
        matrix, construction = evaluation_grid_matrix(k)
        records.append(
            matrix_record(
                name="evaluation_grid_low_degree",
                family="evaluation",
                matrix=matrix,
                construction=construction,
            )
        )

    if "vandermonde" in families:
        matrix, construction = vandermonde_matrix(k)
        records.append(
            matrix_record(
                name="vandermonde_residue_classes",
                family="vandermonde",
                matrix=matrix,
                construction=construction,
            )
        )

    if "block_sums" in families:
        for name, factory in (
            ("block_sums_interleaved", interleaved_block_sum_matrix),
            ("block_sums_contiguous", contiguous_block_sum_matrix),
        ):
            matrix, construction = factory(k)
            records.append(
                matrix_record(
                    name=name,
                    family="block_sums",
                    matrix=matrix,
                    construction=construction,
                )
            )

    if "selector" in families:
        matrix = selector_matrix(k)
        records.append(
            matrix_record(
                name="selector_prefix_baseline",
                family="selector",
                matrix=matrix,
                construction={
                    "strategy": "coordinate_selector_baseline",
                    "selected_coordinates": list(range(k)),
                },
            )
        )

    require(records, "requested matrix families produced no matrices")
    hashes = [row["matrix_sha256_uint8"] for row in records]
    require(len(hashes) == len(set(hashes)), "matrix catalog contains duplicate maps")
    require(all(row["rank_mod7"] == k for row in records), "matrix rank audit failed")
    return records, {
        "requested_families": list(families),
        "matrix_count": len(records),
        "matrix_names_in_audited_order": [row["name"] for row in records],
        "matrix_hashes_in_audited_order": hashes,
        "matrix_catalog_sha256": json_sha256(records),
        "all_matrices_explicit_hash_pinned_and_exact_full_row_rank": True,
        "all_rank_certificates_use_nonzero_k_by_k_minors_mod7": True,
    }


def shard_indexed_rows(
    rows: list[Any], shard_index: int, shard_count: int, label: str
) -> tuple[list[tuple[int, Any]], dict]:
    require(shard_count > 0, f"{label} shard count must be positive")
    require(0 <= shard_index < shard_count, f"{label} shard index is outside count")
    indexed = list(enumerate(rows))
    selected = [row for row in indexed if row[0] % shard_count == shard_index]
    require(selected, f"selected {label} shard is empty")
    all_shards = [
        [index for index in range(shard, len(rows), shard_count)]
        for shard in range(shard_count)
    ]
    flattened = [index for shard in all_shards for index in shard]
    require(
        sorted(flattened) == list(range(len(rows))) and len(flattened) == len(set(flattened)),
        f"{label} shards do not form a disjoint cover",
    )
    return selected, {
        "label": label,
        "rule": f"audited_{label}_index_mod_shard_count_equals_shard_index",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "full_count": len(rows),
        "selected_count": len(selected),
        "selected_indices": [index for index, _row in selected],
        "all_shards_form_a_disjoint_cover": True,
    }


def matrix_from_record(record: dict) -> np.ndarray:
    matrix = np.ascontiguousarray(record["matrix_mod7"], dtype=np.uint8)
    require(matrix.shape == tuple(record["shape"]), "explicit matrix shape changed")
    require(
        cpu_join.array_sha256(matrix) == record["matrix_sha256_uint8"],
        "explicit matrix failed its pinned byte hash",
    )
    rank, pivots = exact_rank_pivots(matrix, 7)
    require(rank == record["rank_mod7"] == len(matrix), "explicit matrix rank changed")
    require(list(pivots) == record["pivot_columns_zero_based"], "matrix pivots changed")
    require(
        determinant_mod_prime(matrix[:, list(pivots)], 7)
        == record["pivot_minor_determinant_mod7"]
        != 0,
        "matrix minor certificate changed",
    )
    return matrix


def apply_linear_projection(full_rows: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Retain F3^6 and apply M to F7^21, row by identical source row."""
    source = np.ascontiguousarray(full_rows, dtype=np.uint8)
    projection = np.ascontiguousarray(matrix, dtype=np.uint8)
    require(
        source.ndim == 2 and source.shape[1] == 27,
        "full torsion rows must have width 6+21",
    )
    require(
        projection.ndim == 2 and projection.shape[1] == 21,
        "linear F7 matrix width changed",
    )
    require(np.all(source[:, :6] < 3), "F3 source digit escaped its field")
    require(np.all(source[:, 6:] < 7), "F7 source digit escaped its field")
    mapped_f7 = (
        source[:, 6:].astype(np.int64) @ projection.astype(np.int64).T
    ) % 7
    result = np.ascontiguousarray(
        np.concatenate((source[:, :6], mapped_f7.astype(np.uint8)), axis=1),
        dtype=np.uint8,
    )
    require(np.array_equal(result[:, :6], source[:, :6]), "F3 coordinates changed")
    return result


def project_equation_vector_linear(
    vector: np.ndarray, quotient_data: dict[int, dict], matrix: np.ndarray
) -> np.ndarray:
    full = cpu_join.project_equation_vector(
        vector, quotient_data, ALL_MOD7_COORDINATES
    )
    return apply_linear_projection(full[None, :], matrix)[0]


def build_linear_direction_supports(
    context: dict[str, Any],
    matrix_record_row: dict,
    pair_chunk_cap: int,
) -> tuple[dict[int, dict[int, np.ndarray]], dict]:
    """Construct all eight exact, direct-catalog-calibrated linear supports."""
    matrix = matrix_from_record(matrix_record_row)
    k = len(matrix)
    moduli = (3,) * 6 + (7,) * k
    codec = cpu_join.semigroup.MixedRadixCodec(moduli)
    effective_pair_chunk_cap = max(
        1, min(pair_chunk_cap, EXACT_DIRECTION_SUPPORT_CAP)
    )
    support_table: dict[int, dict[int, np.ndarray]] = {}
    direction_audits = []

    for direction in range(8):
        full_generators = np.ascontiguousarray(
            context["full_generator_rows"][direction], dtype=np.uint8
        )
        joint_generators = apply_linear_projection(full_generators, matrix)
        generator_codes = {
            degree: codec.unique_codes(
                joint_generators[context["degrees"] == degree]
            )
            for degree in (1, 2, 3)
        }
        raw_supports, recurrence_records = cpu_join.semigroup.support_recurrence(
            generator_codes,
            codec,
            EXACT_DIRECTION_SUPPORT_CAP,
            effective_pair_chunk_cap,
            max_grade=3,
        )

        mask = int(context["orbit"]["masks"][direction])
        translated: dict[int, np.ndarray] = {}
        grade_audits = []
        for grade in range(4):
            recurrence_record = recurrence_records[grade]
            if grade not in raw_supports:
                require(
                    not recurrence_record["completed"],
                    "missing linear support was labeled complete",
                )
                grade_audits.append(
                    {
                        "excess_grade": grade,
                        "decision_status": recurrence_record["decision_status"],
                        "completed": False,
                        "skipped": True,
                        "partial_support_used": False,
                        "rigorous_rejection_allowed": False,
                    }
                )
                continue

            _grade, floor_mean, floor = cpu_join.semigroup.excess_grade(
                mask, 0 if mask.bit_count() == 7 else 8
            )
            mean = floor_mean + 4 * grade
            anchor = context["rebuilt"]["anchors"].get(mask, mean)
            full_offset = cpu_join.semigroup.project_rows(
                (anchor - floor)[None, :],
                direction,
                context["quotient_data"],
                ALL_MOD7_COORDINATES,
            )
            offset = apply_linear_projection(full_offset, matrix)[0]
            observed = codec.translate(raw_supports[grade], offset)

            # Independently enumerate the complete direct catalog and apply
            # exactly the same rowwise linear map before deduplication.
            catalog = cpu_join.affine.mapped_catalog(mask, mean).astype(np.int64)
            require(
                len(catalog) == cpu_join.EXPECTED_GRADE_CATALOG_ROWS[grade],
                f"grade-{grade} complete catalog census changed",
            )
            delta = np.ascontiguousarray(anchor[None, :] - catalog, dtype=np.int64)
            require(not np.any(delta.sum(axis=1)), "direct catalog delta changed mean")
            require(
                not np.any(context["rebuilt"]["anchors"].kernel_rows @ delta.T),
                "direct catalog left degree two",
            )
            require(not np.any(delta % 2), "direct catalog delta changed parity")
            full_direct_rows = cpu_join.semigroup.project_rows(
                delta,
                direction,
                context["quotient_data"],
                ALL_MOD7_COORDINATES,
            )
            direct_rows = apply_linear_projection(full_direct_rows, matrix)
            direct = codec.unique_codes(direct_rows)
            require(
                np.array_equal(observed, direct),
                f"direction {direction} grade-{grade} linear semigroup/direct mismatch",
            )
            translated[grade] = observed
            grade_audits.append(
                {
                    "excess_grade": grade,
                    "scaled_mean": mean,
                    "complete_catalog_rows": len(catalog),
                    "projected_unique_states": len(observed),
                    "support_sha256_uint64": cpu_join.array_sha256(
                        observed.astype("<u8", copy=False)
                    ),
                    "full_projected_group_saturated": len(observed)
                    == codec.group_size,
                    "semigroup_recurrence_equals_complete_direct_catalog_projection": True,
                    "box_exact": True,
                    "same_source_row_used_mod3_mod7_before_linear_map_and_deduplication": True,
                }
            )

        support_table[direction] = translated
        direction_audits.append(
            {
                "direction": direction,
                "mask": mask,
                "linear_matrix_name": matrix_record_row["name"],
                "linear_matrix_sha256_uint8": matrix_record_row[
                    "matrix_sha256_uint8"
                ],
                "projected_generator_unique_states_by_degree": {
                    str(degree): len(generator_codes[degree])
                    for degree in (1, 2, 3)
                },
                "full_F3_6_F7_21_generator_rows_sha256_uint8": cpu_join.array_sha256(
                    full_generators
                ),
                "linearly_projected_joint_generator_sha256_uint8": cpu_join.array_sha256(
                    joint_generators
                ),
                "recurrence_records_sha256": cpu_join.json_sha256(
                    recurrence_records
                ),
                "grades": grade_audits,
                "all_grades_zero_through_three_completed": set(translated)
                == set(range(4)),
            }
        )

    audit = {
        "linear_matrix": matrix_record_row,
        "projected_group": f"F3^6 x F7^{k}",
        "projected_group_size": codec.group_size,
        "state_cap": EXACT_DIRECTION_SUPPORT_CAP,
        "pair_chunk_cap": effective_pair_chunk_cap,
        "directions": direction_audits,
        "all_direction_grades_zero_through_three_completed": all(
            row["all_grades_zero_through_three_completed"]
            for row in direction_audits
        ),
        "all_completed_grade_supports_match_complete_direct_catalogs": True,
        "all_six_F3_coordinates_retained_verbatim": True,
        "identical_explicit_F7_matrix_applied_to_generators_and_direct_catalogs": True,
    }
    require(
        audit["all_direction_grades_zero_through_three_completed"],
        "an exact grade-zero-through-three linear support did not complete",
    )
    require(
        all(set(support_table[direction]) == {0, 1, 2, 3} for direction in range(8)),
        "linear direction support table is incomplete",
    )
    return support_table, audit


def exact_linear_case_inputs(
    context: dict[str, Any],
    target_row: dict,
    matrix_record_row: dict,
    support_table: dict[int, dict[int, np.ndarray]],
) -> tuple[
    list[np.ndarray],
    np.ndarray,
    tuple[int, ...],
    tuple[tuple[int, np.ndarray], ...],
    dict,
]:
    """Map one exact target and its eight factors with the same matrix."""
    matrix = matrix_from_record(matrix_record_row)
    k = len(matrix)
    rebuilt = context["rebuilt"]
    orbit, leaf, system, _factory = cpu_join.old_join.validate_parent_survivor(
        target_row, rebuilt
    )
    grades = tuple(
        cpu_join.leaf_grade(orbit, leaf, direction) for direction in range(8)
    )
    require(max(grades) == 3, "selected case is not a grade-three case")
    require(sum(grades) == 14, "selected case total directional grade changed")

    moduli = (3,) * 6 + (7,) * k
    codec = cpu_join.semigroup.MixedRadixCodec(moduli)
    factors = tuple(
        (direction, support_table[direction][grade])
        for direction, grade in enumerate(grades)
    )
    require(all(len(codes) for _direction, codes in factors), "empty linear support")
    ordered = tuple(
        sorted(
            factors,
            key=lambda row: (
                0 if len(row[1]) == codec.group_size else 1,
                len(row[1]),
                row[0],
            ),
        )
    )

    anchor_rhs, _raw_syndromes = cpu_join.affine.anchor_rhs_and_raw_syndromes(
        orbit, leaf, system, rebuilt["anchors"]
    )
    full_base_digits = cpu_join.project_equation_vector(
        anchor_rhs, context["quotient_data"], ALL_MOD7_COORDINATES
    )
    base_digits = apply_linear_projection(full_base_digits[None, :], matrix)[0]
    target_digits = np.ascontiguousarray(
        (-base_digits.astype(np.int16)) % np.asarray(moduli, dtype=np.int16),
        dtype=np.uint8,
    )

    exact_rhs = anchor_rhs.copy()
    for direction, grade in enumerate(grades):
        mask = int(orbit["masks"][direction])
        mean = int(leaf["scaled_means"][direction])
        catalog_row = cpu_join.affine.mapped_catalog(mask, mean)[0].astype(np.int64)
        anchor = rebuilt["anchors"].get(mask, mean)
        block = slice(
            1 + cpu_join.AMBIENT_DIMENSION * direction,
            1 + cpu_join.AMBIENT_DIMENSION * (direction + 1),
        )
        exact_rhs[block] += anchor - catalog_row
    require(
        not np.any(context["common"] @ exact_rhs),
        "common exact dependency syndrome did not vanish",
    )

    support_rows = []
    interface_records = []
    for direction, codes in ordered:
        rows = codec.decode(codes)
        recoded = cuda_join.unique_support_codes(rows, moduli)
        require(
            np.array_equal(recoded, codes),
            "linear CPU support changed at the CUDA interface",
        )
        support_rows.append(rows)
        interface_records.append(
            {
                "direction": direction,
                "excess_grade": grades[direction],
                "exact_state_count": len(codes),
                "support_sha256_uint64": cpu_join.array_sha256(
                    codes.astype("<u8", copy=False)
                ),
            }
        )

    audit = {
        "case_key": str(target_row["case_key"]),
        "catalog_pattern": context["current_by_key"][str(target_row["case_key"])][
            "catalog_pattern"
        ],
        "linear_matrix_name": matrix_record_row["name"],
        "linear_matrix_sha256_uint8": matrix_record_row["matrix_sha256_uint8"],
        "linear_matrix_rank_mod7": matrix_record_row["rank_mod7"],
        "directional_excess_grades": list(grades),
        "convolution_direction_order": [row[0] for row in ordered],
        "convolution_order_rule": "full_group_first_else_support_size_then_direction",
        "direction_support_records": interface_records,
        "full_F3_6_F7_21_base_digits_sha256_uint8": cpu_join.array_sha256(
            full_base_digits
        ),
        "linearly_mapped_base_digits_sha256_uint8": cpu_join.array_sha256(
            base_digits
        ),
        "target_digits": target_digits.tolist(),
        "target_digits_sha256_uint8": cpu_join.array_sha256(target_digits),
        "common_exact_dependency_syndrome_checked_with_direct_catalog_rows": True,
        "CPU_codes_decode_and_reencode_identically_in_CUDA_engine_codec": True,
        "all_six_F3_coordinates_retained_verbatim": True,
        "identical_matrix_maps_generator_catalog_and_target_F7_coordinates": True,
        "same_direction_rows_couple_mod3_mod7_before_linear_map_and_convolution": True,
    }
    return support_rows, target_digits, moduli, ordered, audit


def manufactured_linear_projection_audit() -> dict:
    """Audit row identity, selector equivalence, and mixed-field linearity."""
    selector = selector_matrix(2, SELECTOR_CROSSCHECK_COORDINATES)
    selector_record = matrix_record(
        name="manufactured_selector_0_1",
        family="selector",
        matrix=selector,
        construction={
            "strategy": "coordinate_selector_crosscheck",
            "selected_coordinates": list(SELECTOR_CROSSCHECK_COORDINATES),
        },
    )
    source = np.asarray(
        [
            [*(index % 3 for index in range(6)), *(index % 7 for index in range(21))],
            [*((2 * index + 1) % 3 for index in range(6)), *((3 * index + 2) % 7 for index in range(21))],
            [*((index + 2) % 3 for index in range(6)), *((5 * index + 4) % 7 for index in range(21))],
        ],
        dtype=np.uint8,
    )
    observed = apply_linear_projection(source, selector)
    expected = np.ascontiguousarray(
        np.concatenate(
            (
                source[:, :6],
                source[:, 6 + np.asarray(SELECTOR_CROSSCHECK_COORDINATES)],
            ),
            axis=1,
        ),
        dtype=np.uint8,
    )
    require(np.array_equal(observed, expected), "selector matrix changed coordinate projection")

    dense, dense_construction = seeded_dense_matrix(5, DEFAULT_DENSE_SEEDS[0])
    left = source[:2]
    right = source[1:]
    mixed_sum = np.empty_like(left)
    mixed_sum[:, :6] = (left[:, :6] + right[:, :6]) % 3
    mixed_sum[:, 6:] = (left[:, 6:] + right[:, 6:]) % 7
    mapped_left = apply_linear_projection(left, dense)
    mapped_right = apply_linear_projection(right, dense)
    expected_sum = np.empty_like(mapped_left)
    expected_sum[:, :6] = (mapped_left[:, :6] + mapped_right[:, :6]) % 3
    expected_sum[:, 6:] = (mapped_left[:, 6:] + mapped_right[:, 6:]) % 7
    require(
        np.array_equal(apply_linear_projection(mixed_sum, dense), expected_sum),
        "mixed-field linearity audit failed",
    )
    return {
        "status": "passed",
        "selector_matrix": selector_record,
        "manufactured_source_rows_sha256_uint8": cpu_join.array_sha256(source),
        "selector_matrix_equals_coordinate_selection_elementwise": True,
        "all_six_F3_coordinates_retained_elementwise": True,
        "seeded_dense_mixed_F3_F7_linearity_checked_elementwise": True,
        "seeded_dense_construction": dense_construction,
        "same_input_row_supplies_retained_F3_and_mapped_F7_output": True,
    }


def small_real_selector_crosscheck(
    cp: Any,
    kernel: Any,
    gpu_memory_cap_bytes: int,
    context: dict[str, Any],
    pair_chunk_cap: int,
) -> dict:
    """Compare the linear selector path with the established coordinate path."""
    projection = SELECTOR_CROSSCHECK_COORDINATES
    selector = selector_matrix(len(projection), projection)
    selector_record = matrix_record(
        name="real_selector_0_1",
        family="selector",
        matrix=selector,
        construction={
            "strategy": "coordinate_selector_crosscheck",
            "selected_coordinates": list(projection),
        },
    )

    coordinate_table, coordinate_direction_audit = coordinate_gpu.build_direction_supports(
        context, projection, pair_chunk_cap
    )
    linear_table, linear_direction_audit = build_linear_direction_supports(
        context, selector_record, pair_chunk_cap
    )
    support_comparisons = []
    for direction in range(8):
        for grade in range(4):
            coordinate_codes = coordinate_table[direction][grade]
            linear_codes = linear_table[direction][grade]
            require(
                np.array_equal(coordinate_codes, linear_codes),
                f"selector support differs at direction {direction}, grade {grade}",
            )
            support_comparisons.append(
                {
                    "direction": direction,
                    "excess_grade": grade,
                    "state_count": len(linear_codes),
                    "support_sha256_uint64": cpu_join.array_sha256(
                        linear_codes.astype("<u8", copy=False)
                    ),
                }
            )

    target_row = context["targets"][0]
    linear_rows, linear_target, linear_moduli, linear_ordered, linear_case_audit = (
        exact_linear_case_inputs(
            context, target_row, selector_record, linear_table
        )
    )
    (
        coordinate_rows,
        coordinate_target,
        coordinate_moduli,
        coordinate_ordered,
        coordinate_case_audit,
    ) = coordinate_gpu.exact_case_inputs(
        context, target_row, projection, coordinate_table
    )
    require(linear_moduli == coordinate_moduli, "selector moduli differ")
    require(np.array_equal(linear_target, coordinate_target), "selector targets differ")
    require(
        [direction for direction, _codes in linear_ordered]
        == [direction for direction, _codes in coordinate_ordered],
        "selector convolution orders differ",
    )
    require(
        all(
            np.array_equal(linear_codes, coordinate_codes)
            for (_left_direction, linear_codes), (_right_direction, coordinate_codes)
            in zip(linear_ordered, coordinate_ordered)
        ),
        "selector case factors differ",
    )
    require(
        all(np.array_equal(left, right) for left, right in zip(linear_rows, coordinate_rows)),
        "selector CUDA interface rows differ",
    )

    codec = cpu_join.semigroup.MixedRadixCodec(linear_moduli)
    cpu_codes, cpu_convolution = cpu_join.convolve_support_sequence(
        linear_ordered,
        codec,
        state_cap=codec.group_size,
        pair_chunk_cap=min(pair_chunk_cap, codec.group_size),
    )
    require(cpu_codes is not None and cpu_convolution["completed"], "CPU crosscheck capped")
    gpu_decision, gpu_codes = cuda_join.gpu_exact_support_convolution(
        cp,
        kernel,
        linear_rows,
        linear_target,
        linear_moduli,
        state_cap=codec.group_size,
        pair_chunk_cap=pair_chunk_cap,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    coordinate_gpu.validate_cuda_decision(gpu_decision)
    require(gpu_codes is not None, "CUDA selector crosscheck skipped")
    require(np.array_equal(cpu_codes, gpu_codes), "linear selector CPU/CUDA sets differ")
    target_code = int(codec.encode(linear_target[None, :])[0])
    target_position = int(np.searchsorted(cpu_codes, np.uint64(target_code)))
    cpu_present = bool(
        target_position < len(cpu_codes)
        and int(cpu_codes[target_position]) == target_code
    )
    require(
        cpu_present == gpu_decision["projected_target_present"],
        "linear selector CPU/CUDA target decisions differ",
    )
    return {
        "status": "passed",
        "selector_matrix": selector_record,
        "coordinate_projection_mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "coordinate_direction_support_audit": coordinate_direction_audit,
        "linear_direction_support_audit": linear_direction_audit,
        "support_comparison_records": support_comparisons,
        "support_comparison_records_sha256": json_sha256(support_comparisons),
        "linear_case_input_audit": linear_case_audit,
        "coordinate_case_input_audit": coordinate_case_audit,
        "CPU_convolution_audit": cpu_convolution,
        "CUDA_convolution_audit": gpu_decision,
        "complete_support_state_count": len(cpu_codes),
        "complete_support_sha256_uint64": cpu_join.array_sha256(
            cpu_codes.astype("<u8", copy=False)
        ),
        "target_present": cpu_present,
        "all_32_linear_selector_and_coordinate_supports_equal_elementwise": True,
        "linear_selector_and_coordinate_target_equal_elementwise": True,
        "linear_selector_and_coordinate_case_factors_equal_elementwise": True,
        "linear_selector_CPU_and_CUDA_complete_support_sets_equal_elementwise": True,
        "linear_selector_CPU_and_CUDA_target_decisions_equal": True,
    }


def run_production(
    *,
    cp: Any,
    kernel: Any,
    gpu_audit: dict,
    context: dict[str, Any],
    context_audit: dict,
    manufactured_cross_engine: dict,
    manufactured_linear: dict,
    selector_crosscheck: dict | None,
    matrix_catalog: list[dict],
    matrix_catalog_audit: dict,
    selected_matrices: list[tuple[int, dict]],
    matrix_shard_audit: dict,
    selected_cases: list[tuple[int, dict]],
    case_shard_audit: dict,
    state_cap: int,
    pair_chunk_cap: int,
    gpu_memory_cap_bytes: int,
    output_path: Path,
) -> dict:
    started = time.time()
    support_tables: dict[str, dict[int, dict[int, np.ndarray]]] = {}
    support_audits = []
    for matrix_index, matrix_row in selected_matrices:
        table, audit = build_linear_direction_supports(
            context, matrix_row, pair_chunk_cap
        )
        matrix_hash = str(matrix_row["matrix_sha256_uint8"])
        require(matrix_hash not in support_tables, "selected matrix hash repeated")
        support_tables[matrix_hash] = table
        support_audits.append(
            {
                "matrix_index_in_audited_catalog": matrix_index,
                **audit,
            }
        )

    projection_totals: Counter[str] = Counter()
    case_results = []
    for target_index, target_row in selected_cases:
        projection_results = []
        for matrix_index, matrix_row in selected_matrices:
            matrix_hash = str(matrix_row["matrix_sha256_uint8"])
            support_rows, target, moduli, _ordered, case_input_audit = (
                exact_linear_case_inputs(
                    context,
                    target_row,
                    matrix_row,
                    support_tables[matrix_hash],
                )
            )
            decision, _final_codes = cuda_join.gpu_exact_support_convolution(
                cp,
                kernel,
                support_rows,
                target,
                moduli,
                state_cap=state_cap,
                pair_chunk_cap=pair_chunk_cap,
                gpu_memory_cap_bytes=gpu_memory_cap_bytes,
            )
            coordinate_gpu.validate_cuda_decision(decision)
            projection_totals[str(decision["decision_status"])] += 1
            projection_results.append(
                {
                    "matrix_index_in_audited_catalog": matrix_index,
                    "matrix_name": matrix_row["name"],
                    "matrix_family": matrix_row["family"],
                    "matrix_sha256_uint8": matrix_hash,
                    "matrix_rank_mod7": matrix_row["rank_mod7"],
                    "retained_all_six_mod3_coordinates_verbatim": True,
                    "case_input_audit": case_input_audit,
                    **decision,
                }
            )
            cp.get_default_memory_pool().free_all_blocks()

        rejected = any(row["rigorous_rejection"] for row in projection_results)
        skipped = not rejected and any(
            str(row["decision_status"]).startswith("skipped_")
            for row in projection_results
        )
        necessary = not rejected and not skipped
        require(sum((rejected, skipped, necessary)) == 1, "case decision is ambiguous")
        case_key = str(target_row["case_key"])
        row = {
            "target_index_in_audited_51_case_order": target_index,
            "case_key": case_key,
            "catalog_pattern": context["current_by_key"][case_key]["catalog_pattern"],
            "prior_global_join_decision": context["current_by_key"][case_key][
                "decision_status"
            ],
            "linear_projection_results": projection_results,
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
            "decision_status": (
                "rigorous_exact_CUDA_linear_semigroup_projection_rejection"
                if rejected
                else "explicit_cap_skip_without_negative_decision"
                if skipped
                else "necessary_only_survivor_of_all_completed_linear_CUDA_projections"
            ),
        }
        row["decision_certificate_sha256"] = json_sha256(row)
        case_results.append(row)

    counts = {
        "selected": len(case_results),
        "rejected": sum(row["rigorously_rejected"] for row in case_results),
        "surviving": sum(row["necessary_only_survivor"] for row in case_results),
        "skipped": sum(row["skipped"] for row in case_results),
    }
    require(
        counts["rejected"] + counts["surviving"] + counts["skipped"]
        == counts["selected"],
        "case result census is not a partition",
    )
    k = int(matrix_catalog[0]["shape"][0])
    script_path = Path(__file__).resolve()
    return {
        "experiment": EXPERIMENT,
        "status": "complete_sharded_exact_CUDA_linear_semigroup_case_join",
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "51 orbit0/branch-A grade-three-only representatives",
        "cuda_invoked": True,
        "solver_invoked": False,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "coordinate_CUDA_companion_path": str(Path(coordinate_gpu.__file__).resolve()),
            "coordinate_CUDA_companion_sha256": file_sha256(
                Path(coordinate_gpu.__file__)
            ),
            "CPU_semigroup_case_join_path": str(Path(cpu_join.__file__).resolve()),
            "CPU_semigroup_case_join_sha256": file_sha256(Path(cpu_join.__file__)),
            "CUDA_torsion_engine_path": str(Path(cuda_join.__file__).resolve()),
            "CUDA_torsion_engine_sha256": file_sha256(Path(cuda_join.__file__)),
        },
        "configuration": {
            "k": k,
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "gpu_memory_cap_bytes": gpu_memory_cap_bytes,
            "exact_direction_support_cap": EXACT_DIRECTION_SUPPORT_CAP,
            "selected_matrix_indices": [index for index, _row in selected_matrices],
            "selected_matrix_hashes": [
                row["matrix_sha256_uint8"] for _index, row in selected_matrices
            ],
        },
        "expected_memory_per_linear_projection": coordinate_gpu.expected_engine_memory(
            k, state_cap, pair_chunk_cap
        ),
        "gpu_engine_audit": gpu_audit,
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_CPU_CUDA_audit": manufactured_cross_engine,
        "manufactured_linear_projection_audit": manufactured_linear,
        "small_real_selector_matrix_crosscheck": selector_crosscheck,
        "small_real_selector_crosscheck_skipped_by_explicit_flag": selector_crosscheck
        is None,
        "full_linear_matrix_catalog_audit": matrix_catalog_audit,
        "full_linear_matrix_catalog": matrix_catalog,
        "matrix_sharding_audit": matrix_shard_audit,
        "case_sharding_audit": case_shard_audit,
        "linear_direction_support_audits": support_audits,
        "result_counts": counts,
        "projection_decision_census": dict(sorted(projection_totals.items())),
        "case_results_sha256": cpu_join.old_join.canonical_case_digest(case_results),
        "case_results": case_results,
        "logical_semantics": {
            "all_51_cases_reconstructed_before_case_sharding": True,
            "full_matrix_catalog_rank_audited_before_matrix_sharding": True,
            "all_six_F3_coordinates_retained_verbatim": True,
            "every_F7_map_is_an_explicit_hash_pinned_full_rank_k_by_21_matrix": True,
            "same_Hilbert_or_catalog_row_supplies_mod3_and_mod7_before_linear_map": True,
            "identical_linear_matrix_maps_generators_catalogs_and_target": True,
            "grade_zero_through_three_direction_supports_are_exact_and_direct_catalog_calibrated": True,
            "all_eight_direction_supports_convolved_in_every_completed_projection": True,
            "CUDA_hash_collisions_resolved_by_exact_key_comparison": True,
            "exact_full_group_saturation_only_proves_target_presence": True,
            "missing_target_after_completed_convolution_is_a_rigorous_rejection": True,
            "target_presence_is_necessary_only": True,
            "all_state_memory_allocation_or_hash_capacity_caps_are_explicit_skips": True,
            "partial_support_after_any_cap_is_never_used": True,
            "binary_edge_feasibility_claimed": False,
            "positive_z7_closure_claimed": False,
        },
        "positive_z7_excluded": False,
        "full_theorem_claimed": False,
        "output_path": str(output_path.resolve()),
        "elapsed_seconds": time.time() - started,
    }


def catalog_only_result(
    *,
    context_audit: dict,
    matrix_catalog: list[dict],
    matrix_catalog_audit: dict,
    manufactured_linear: dict,
    output_path: Path,
) -> dict:
    script_path = Path(__file__).resolve()
    return {
        "experiment": f"{EXPERIMENT}_matrix_catalog",
        "status": "complete_deterministic_full_rank_linear_matrix_catalog",
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "coordinate_CUDA_companion_path": str(Path(coordinate_gpu.__file__).resolve()),
            "coordinate_CUDA_companion_sha256": file_sha256(
                Path(coordinate_gpu.__file__)
            ),
        },
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_linear_projection_audit": manufactured_linear,
        "full_linear_matrix_catalog_audit": matrix_catalog_audit,
        "full_linear_matrix_catalog": matrix_catalog,
        "cuda_loaded": False,
        "production_cases_processed": 0,
        "output_path": str(output_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--parent-input", type=Path, default=cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument("--k", type=int, choices=SUPPORTED_K, default=5)
    parser.add_argument(
        "--matrix-families",
        default=",".join(DEFAULT_MATRIX_FAMILIES),
        help=(
            "comma-separated subset of seeded_dense,geometry,evaluation,"
            "vandermonde,block_sums,selector"
        ),
    )
    parser.add_argument("--matrix-shard-index", type=int, default=0)
    parser.add_argument("--matrix-shard-count", type=int, default=1)
    parser.add_argument("--case-shard-index", type=int, default=0)
    parser.add_argument("--case-shard-count", type=int, default=1)
    parser.add_argument(
        "--state-cap",
        type=int,
        default=0,
        help="global support cap; 0 means the complete projected-group size",
    )
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument(
        "--gpu-memory-cap-mib", type=int, default=DEFAULT_GPU_MEMORY_CAP_MIB
    )
    parser.add_argument("--skip-small-real-crosscheck", action="store_true")
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--matrix-catalog-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    require(args.gpu_memory_cap_mib > 0, "--gpu-memory-cap-mib must be positive")
    require(args.state_cap >= 0, "--state-cap cannot be negative")
    require(
        not (args.self_audit_only and args.matrix_catalog_only),
        "self-audit-only and matrix-catalog-only are mutually exclusive",
    )
    families = parse_matrix_families(args.matrix_families)

    context, context_audit = coordinate_gpu.build_cpu_context(
        args.parent_input, args.current_join, args.hilbert_basis
    )
    require(len(context["targets"]) == EXPECTED_TARGET_CASES, "target census changed")
    matrix_catalog, matrix_catalog_audit = build_matrix_catalog(
        context, args.k, families
    )
    manufactured_linear = manufactured_linear_projection_audit()

    if args.matrix_catalog_only:
        result = catalog_only_result(
            context_audit=context_audit,
            matrix_catalog=matrix_catalog,
            matrix_catalog_audit=matrix_catalog_audit,
            manufactured_linear=manufactured_linear,
            output_path=args.output,
        )
        cuda_join.atomic_write(args.output, result)
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "output": str(args.output),
                    "k": args.k,
                    "matrix_count": len(matrix_catalog),
                    "matrix_names": [row["name"] for row in matrix_catalog],
                    "matrix_hashes": [
                        row["matrix_sha256_uint8"] for row in matrix_catalog
                    ],
                },
                indent=2,
                sort_keys=True,
            ),
            flush=True,
        )
        return

    group_size = cuda_join.group_size_for((3,) * 6 + (7,) * args.k)
    state_cap = group_size if args.state_cap == 0 else min(args.state_cap, group_size)
    gpu_memory_cap_bytes = args.gpu_memory_cap_mib * (1 << 20)

    cp, kernel, gpu_audit = cuda_join.load_cupy(args.device)
    with cp.cuda.Device(args.device):
        manufactured_cross_engine = coordinate_gpu.manufactured_cross_engine_audit(
            cp, kernel, gpu_memory_cap_bytes
        )
        selector_crosscheck = (
            None
            if args.skip_small_real_crosscheck and not args.self_audit_only
            else small_real_selector_crosscheck(
                cp,
                kernel,
                gpu_memory_cap_bytes,
                context,
                args.pair_chunk_cap,
            )
        )

        if args.self_audit_only:
            result = {
                "experiment": f"{EXPERIMENT}_self_audit",
                "status": "manufactured_and_real_selector_CPU_CUDA_self_audits_passed",
                "source_provenance": {
                    "this_script_path": str(Path(__file__).resolve()),
                    "this_script_sha256": file_sha256(Path(__file__).resolve()),
                    "coordinate_CUDA_companion_sha256": file_sha256(
                        Path(coordinate_gpu.__file__)
                    ),
                },
                "gpu_engine_audit": gpu_audit,
                "CPU_construction_and_input_audits": context_audit,
                "manufactured_CPU_CUDA_audit": manufactured_cross_engine,
                "manufactured_linear_projection_audit": manufactured_linear,
                "small_real_selector_matrix_crosscheck": selector_crosscheck,
                "full_linear_matrix_catalog_audit": matrix_catalog_audit,
                "full_linear_matrix_catalog": matrix_catalog,
                "production_cases_processed": 0,
                "output_path": str(args.output.resolve()),
            }
        else:
            selected_matrices, matrix_shard_audit = shard_indexed_rows(
                matrix_catalog,
                args.matrix_shard_index,
                args.matrix_shard_count,
                "matrix",
            )
            matrix_shard_audit.update(
                {
                    "selected_matrix_names": [
                        row["name"] for _index, row in selected_matrices
                    ],
                    "selected_matrix_hashes": [
                        row["matrix_sha256_uint8"]
                        for _index, row in selected_matrices
                    ],
                }
            )
            selected_cases, case_shard_audit = coordinate_gpu.shard_cases(
                context["targets"], args.case_shard_index, args.case_shard_count
            )
            case_shard_audit["label"] = "case"
            result = run_production(
                cp=cp,
                kernel=kernel,
                gpu_audit=gpu_audit,
                context=context,
                context_audit=context_audit,
                manufactured_cross_engine=manufactured_cross_engine,
                manufactured_linear=manufactured_linear,
                selector_crosscheck=selector_crosscheck,
                matrix_catalog=matrix_catalog,
                matrix_catalog_audit=matrix_catalog_audit,
                selected_matrices=selected_matrices,
                matrix_shard_audit=matrix_shard_audit,
                selected_cases=selected_cases,
                case_shard_audit=case_shard_audit,
                state_cap=state_cap,
                pair_chunk_cap=args.pair_chunk_cap,
                gpu_memory_cap_bytes=gpu_memory_cap_bytes,
                output_path=args.output,
            )

    cuda_join.atomic_write(args.output, result)
    summary = {
        "status": result["status"],
        "output": str(args.output),
        "gpu": gpu_audit["device_name"],
        "k": args.k,
        "manufactured_CPU_CUDA_audit": manufactured_cross_engine["status"],
        "small_real_selector_matrix_crosscheck": (
            None if selector_crosscheck is None else selector_crosscheck["status"]
        ),
        "processed_cases": result.get("result_counts", {}).get("selected", 0),
        "rigorously_rejected_cases": result.get("result_counts", {}).get(
            "rejected", 0
        ),
        "necessary_only_survivors": result.get("result_counts", {}).get(
            "surviving", 0
        ),
        "skipped_cases": result.get("result_counts", {}).get("skipped", 0),
        "selected_matrices": result.get("matrix_sharding_audit", {}).get(
            "selected_matrix_names", []
        ),
        "expected_memory": coordinate_gpu.expected_engine_memory(
            args.k, state_cap, args.pair_chunk_cap
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
