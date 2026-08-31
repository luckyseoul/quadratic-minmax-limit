#!/usr/bin/env python3
"""Exact F3 and coordinate-product audit for the 51 positive p=7, z=7 cases.

The imported coordinate artifacts contain exact supports

    S_k <= F3^6 x F7^k.

This script independently reconstructs their common mathematical context and
the complete F3-only directional and global supports.  It then applies the
finite-fiber cardinality lemma: if ``S <= T x B`` and
``|S| = |T| |B|``, every fiber has its maximum size and hence
``S = T x B``.  Thus a completed coordinate artifact with
``|S_k| = |T| 7^k`` certifies that its F7 coordinates impose no restriction.

Coordinate artifacts are hostile input.  Their source hashes, reconstruction
fingerprints, case certificates, shard metadata, exact completion, prefix
coordinates, and disjoint coverage are checked before any product claim.
Partial disjoint coverage is reported as partial and never promoted to a
51-case statement.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join as cpu_join  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_join_gpu as case_gpu  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_coordinate_product_audit"
EXPECTED_CASES = 51
EXPECTED_DIRECTIONS = 8
F3_DIMENSION = 6
SUPPORTED_K = (5, 6)
LEFT_DIRECTIONS = (1, 3, 4, 6)
RIGHT_DIRECTIONS = (0, 2, 5, 7)
ENGINE_SPECS = {
    "p7_infinity7_positive_z7_semigroup_case_join_gpu": {
        "status": "complete_sharded_exact_CUDA_semigroup_case_join",
        "script": ROOT / "scripts/p7_infinity7_positive_z7_semigroup_case_join_gpu.py",
        "collision_flag": "cuda_hash_set_is_exact_not_probabilistic",
    },
    "p7_infinity7_positive_z7_semigroup_case_join_opencl": {
        "status": "complete_sharded_exact_OpenCL_semigroup_case_join",
        "script": ROOT / "scripts/p7_infinity7_positive_z7_semigroup_case_join_opencl.py",
        "collision_flag": "OpenCL_uint64_hash_collisions_are_resolved_by_exact_key_comparison",
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def require_true(value: Any, message: str) -> None:
    require(value is True, message)


def canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def json_sha256(value: object) -> str:
    return hashlib.sha256(canonical_json(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def atomic_write_json(path: Path, payload: dict) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def gf_rank(matrix: np.ndarray, prime: int = 3) -> tuple[int, list[int]]:
    rows = np.ascontiguousarray(matrix, dtype=np.int64) % prime
    require(rows.ndim == 2, "GF rank input is not a matrix")
    work = rows.copy()
    rank = 0
    pivots: list[int] = []
    for column in range(work.shape[1]):
        pivot = next(
            (row for row in range(rank, work.shape[0]) if int(work[row, column]) != 0),
            None,
        )
        if pivot is None:
            continue
        work[[rank, pivot]] = work[[pivot, rank]]
        inverse = pow(int(work[rank, column]), -1, prime)
        work[rank] = work[rank] * inverse % prime
        for row in range(work.shape[0]):
            if row != rank and int(work[row, column]):
                work[row] = (work[row] - work[row, column] * work[rank]) % prime
        pivots.append(column)
        rank += 1
        if rank == work.shape[0]:
            break
    return rank, pivots


def normalized_line_vector(rows: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(rows, dtype=np.uint8) % 3
    nonzero = rows[np.any(rows, axis=1)]
    require(len(nonzero) > 0, "direction has zero F3 image")
    normalized = []
    for row in nonzero:
        pivot = int(np.flatnonzero(row)[0])
        normalized.append((row.astype(np.int64) * pow(int(row[pivot]), -1, 3)) % 3)
    distinct = np.unique(np.asarray(normalized, dtype=np.uint8), axis=0)
    require(len(distinct) == 1, "rank-one direction has inconsistent normalized vectors")
    return np.ascontiguousarray(distinct[0], dtype=np.uint8)


def scalar_on_line(row: np.ndarray, vector: np.ndarray) -> int:
    row = np.ascontiguousarray(row, dtype=np.uint8) % 3
    vector = np.ascontiguousarray(vector, dtype=np.uint8) % 3
    pivot = int(np.flatnonzero(vector)[0])
    scalar = int(row[pivot])
    require(
        np.array_equal(row, (scalar * vector.astype(np.int64) % 3).astype(np.uint8)),
        "F3 support row escaped its audited directional line",
    )
    return scalar


def selected_fingerprint(audit: dict) -> dict:
    inputs = audit["input_and_current_decision_audit"]
    basis = audit["normaliz_Hilbert_basis"]
    quotient = audit["varying_torsion_quotient_audit"]
    generators = audit["full_generator_projection_audit"]
    grades = audit["target_grade_audit"]
    return {
        "parent_file_sha256": inputs["parent_input"]["file_sha256"],
        "parent_case_results_sha256": inputs["parent_input"]["all_case_results_sha256"],
        "parent_survivor_case_keys_sha256": inputs["parent_input"][
            "survivor_case_keys_sha256"
        ],
        "current_join_file_sha256": inputs["current_global_join"]["file_sha256"],
        "current_join_case_results_sha256": inputs["current_global_join"][
            "case_results_sha256"
        ],
        "target_case_count": inputs["target_case_count"],
        "target_case_keys_sha256": inputs["target_case_keys_sha256"],
        "target_pattern_counts": inputs["target_pattern_counts"],
        "target_grade_vectors_sha256": grades["directional_grade_vectors_sha256"],
        "target_grade_records_sha256": grades["case_grade_records_sha256"],
        "Hilbert_file_sha256": basis["file_sha256"],
        "Hilbert_basis_sha256_int64": basis["basis_sha256_int64"],
        "Hilbert_degrees_sha256_int64": basis["degrees_sha256_int64"],
        "quotient_dimensions": quotient["actual_effective_dimensions"],
        "pointed_matrix_sha256_int64": quotient["integer_pointed_matrix_sha256_int64"],
        "full_generator_records_sha256": generators["records_sha256"],
    }


def build_context(
    parent_path: Path, current_join_path: Path, hilbert_basis_path: Path
) -> tuple[dict[str, Any], dict]:
    context, audit = case_gpu.build_cpu_context(
        parent_path, current_join_path, hilbert_basis_path
    )
    require(len(context["targets"]) == EXPECTED_CASES, "target census changed")
    require_true(
        audit["all_51_grade_three_cases_reconstructed_before_sharding"],
        "context was not reconstructed before sharding",
    )
    return context, audit


def build_f3_supports(context: dict[str, Any]) -> tuple[dict, dict]:
    table, imported_audit = case_gpu.build_direction_supports(
        context, (), case_gpu.EXACT_DIRECTION_SUPPORT_CAP
    )
    codec = cpu_join.semigroup.MixedRadixCodec((3,) * F3_DIMENSION)
    direction_records = []
    line_vectors = []
    for direction in range(EXPECTED_DIRECTIONS):
        generator_rows = np.ascontiguousarray(
            context["full_generator_rows"][direction][:, :F3_DIMENSION],
            dtype=np.uint8,
        )
        rank, pivots = gf_rank(generator_rows, 3)
        require(rank == 1, f"direction {direction} F3 image rank is not one")
        vector = normalized_line_vector(generator_rows)
        line_vectors.append(vector)
        grade_records = []
        for grade in range(4):
            codes = np.ascontiguousarray(table[direction][grade], dtype=np.uint64)
            require(np.array_equal(codes, np.unique(codes)), "direction support is not sorted unique")
            rows = codec.decode(codes)
            scalars = sorted({scalar_on_line(row, vector) for row in rows})
            rebuilt_rows = np.asarray(
                [(scalar * vector.astype(np.int64) % 3).tolist() for scalar in scalars],
                dtype=np.uint8,
            )
            require(
                np.array_equal(np.unique(rows, axis=0), np.unique(rebuilt_rows, axis=0)),
                "scalar directional support does not reconstruct exact F3 rows",
            )
            grade_records.append(
                {
                    "grade": grade,
                    "state_count": len(codes),
                    "support_codes_sha256_uint64": array_sha256(codes.astype("<u8", copy=False)),
                    "support_rows": rows.tolist(),
                    "scalar_support": scalars,
                    "exact_rows_equal_scalar_multiples_of_direction_vector": True,
                }
            )
        direction_records.append(
            {
                "direction": direction,
                "F3_generator_image_rank": rank,
                "rref_pivot_columns": pivots,
                "normalized_line_vector": vector.tolist(),
                "generator_rows_sha256_uint8": array_sha256(generator_rows),
                "grades": grade_records,
            }
        )

    scalar_matrix = np.ascontiguousarray(np.stack(line_vectors, axis=1), dtype=np.uint8)
    total_rank, total_pivots = gf_rank(scalar_matrix, 3)
    left_rank, left_pivots = gf_rank(scalar_matrix[:, LEFT_DIRECTIONS], 3)
    right_rank, right_pivots = gf_rank(scalar_matrix[:, RIGHT_DIRECTIONS], 3)
    intersection_dimension = left_rank + right_rank - total_rank
    require(total_rank == 6, "total F3 directional image rank changed")
    require(left_rank == right_rank == 3, "declared F3 split is not 3+3")
    require(intersection_dimension == 0, "declared F3 split is not a direct sum")
    kernel_dimension = EXPECTED_DIRECTIONS - total_rank
    require(kernel_dimension == 2, "F3 scalar-map kernel dimension changed")
    require(3**kernel_dimension == 9, "F3 target bucket bound changed")
    return table, {
        "group": "F3^6",
        "group_size": codec.group_size,
        "imported_exact_direction_support_audit_sha256": json_sha256(imported_audit),
        "directions": direction_records,
        "scalar_map_matrix_6_by_8_columns_are_direction_vectors": scalar_matrix.tolist(),
        "scalar_map_matrix_sha256_uint8": array_sha256(scalar_matrix),
        "direction_image_ranks": [row["F3_generator_image_rank"] for row in direction_records],
        "total_direction_image_rank": total_rank,
        "total_rref_pivot_columns": total_pivots,
        "direct_sum_split": {
            "left_directions": list(LEFT_DIRECTIONS),
            "right_directions": list(RIGHT_DIRECTIONS),
            "left_rank": left_rank,
            "right_rank": right_rank,
            "sum_rank": total_rank,
            "intersection_dimension": intersection_dimension,
            "left_rref_pivot_columns_within_split": left_pivots,
            "right_rref_pivot_columns_within_split": right_pivots,
            "certified_3_plus_3_direct_sum": True,
        },
        "scalar_map_domain_dimension": EXPECTED_DIRECTIONS,
        "scalar_map_codomain_dimension": F3_DIMENSION,
        "scalar_map_rank": total_rank,
        "scalar_map_kernel_dimension": kernel_dimension,
        "maximum_target_bucket_patterns": 3**kernel_dimension,
        "all_eight_direction_images_have_rank_one": True,
    }


def enumerate_target_buckets(matrix: np.ndarray, target: np.ndarray) -> list[list[int]]:
    rows = []
    for scalars in itertools.product(range(3), repeat=EXPECTED_DIRECTIONS):
        vector = np.asarray(scalars, dtype=np.int64)
        if np.array_equal((matrix.astype(np.int64) @ vector) % 3, target):
            rows.append(list(scalars))
    require(len(rows) <= 9, "target has more than nine F3 scalar buckets")
    return rows


def exact_f3_cases(
    context: dict[str, Any], support_table: dict, support_audit: dict
) -> tuple[list[dict], dict[str, dict]]:
    codec = cpu_join.semigroup.MixedRadixCodec((3,) * F3_DIMENSION)
    matrix = np.asarray(
        support_audit["scalar_map_matrix_6_by_8_columns_are_direction_vectors"],
        dtype=np.uint8,
    )
    line_vectors = [matrix[:, direction] for direction in range(EXPECTED_DIRECTIONS)]
    records = []
    by_key: dict[str, dict] = {}
    for target_index, target_row in enumerate(context["targets"]):
        support_rows, target, moduli, ordered, input_audit = case_gpu.exact_case_inputs(
            context, target_row, (), support_table
        )
        require(moduli == (3,) * F3_DIMENSION, "F3 case input retained a non-F3 coordinate")
        global_codes, convolution_audit = cpu_join.convolve_support_sequence(
            ordered, codec, codec.group_size, codec.group_size
        )
        require(global_codes is not None, "F3 global support unexpectedly capped")
        require_true(convolution_audit["completed"], "F3 global support is incomplete")
        target_code = int(codec.encode(target[None, :])[0])
        target_position = int(np.searchsorted(global_codes, np.uint64(target_code)))
        target_present = bool(
            target_position < len(global_codes)
            and int(global_codes[target_position]) == target_code
        )

        grades = tuple(int(value) for value in input_audit["directional_excess_grades"])
        allowed_scalars = []
        for direction, grade in enumerate(grades):
            decoded = codec.decode(support_table[direction][grade])
            allowed_scalars.append(
                sorted({scalar_on_line(row, line_vectors[direction]) for row in decoded})
            )
        target_buckets = enumerate_target_buckets(matrix, target.astype(np.int64))
        require(len(target_buckets) == 9, "surjective F3 scalar map did not give nine target buckets")
        compatible = [
            row
            for row in target_buckets
            if all(row[direction] in allowed_scalars[direction] for direction in range(8))
        ]
        scalar_image_rows = {
            tuple((matrix.astype(np.int64) @ np.asarray(row, dtype=np.int64) % 3).tolist())
            for row in itertools.product(*allowed_scalars)
        }
        scalar_codes = codec.unique_codes(np.asarray(sorted(scalar_image_rows), dtype=np.uint8))
        require(
            np.array_equal(global_codes, scalar_codes),
            "scalar bucket image differs from exact F3 support convolution",
        )
        require(target_present is bool(compatible), "target presence and scalar buckets disagree")
        require(target_present, "a grade-three case target is absent from its exact F3 support")

        key = str(target_row["case_key"])
        record = {
            "target_index_in_audited_51_case_order": target_index,
            "case_key": key,
            "catalog_pattern": context["current_by_key"][key]["catalog_pattern"],
            "directional_excess_grades": list(grades),
            "directional_scalar_supports": allowed_scalars,
            "target_digits_F3_6": target.tolist(),
            "target_code": target_code,
            "target_digits_sha256_uint8": array_sha256(target),
            "global_support_state_count": len(global_codes),
            "global_support_sha256_uint64": array_sha256(
                global_codes.astype("<u8", copy=False)
            ),
            "global_support_equals_scalar_bucket_image": True,
            "scalar_target_bucket_count_before_directional_restrictions": len(target_buckets),
            "compatible_scalar_target_bucket_count": len(compatible),
            "compatible_scalar_target_buckets": compatible,
            "maximum_target_bucket_patterns": 9,
            "target_present": target_present,
            "exact_F3_convolution_completed": True,
        }
        record["record_sha256"] = json_sha256(record)
        records.append(record)
        by_key[key] = record
    require(len(records) == EXPECTED_CASES, "F3 case census changed")
    require(len(by_key) == EXPECTED_CASES, "F3 case keys are not unique")
    return records, by_key


def product_cardinality_holds(base_size: int, fiber_size: int, joint_size: int) -> bool:
    require(base_size > 0 and fiber_size > 0 and joint_size >= 0, "invalid finite-fiber sizes")
    return joint_size == base_size * fiber_size


def manufactured_fiber_audit() -> dict:
    base = (0, 2)
    fiber = tuple(range(7))
    proper_product = {(a, b) for a in base for b in fiber}
    require(len(base) < 3, "manufactured base is not proper")
    require(
        product_cardinality_holds(len(base), len(fiber), len(proper_product)),
        "proper-product trap was not certified",
    )
    non_product = proper_product - {(0, 0)}
    non_product_projection = {a for a, _b in non_product}
    require(non_product_projection == set(base), "non-product trap changed its base projection")
    require(
        not product_cardinality_holds(len(base), len(fiber), len(non_product)),
        "non-product fiber trap passed the cardinality lemma",
    )
    fiber_sizes = {a: sum(x == a for x, _b in non_product) for a in base}
    require(sorted(fiber_sizes.values()) == [6, 7], "non-product trap fiber sizes changed")
    return {
        "passed": True,
        "proper_base": list(base),
        "ambient_base_size": 3,
        "fiber": "F7",
        "proper_product_size": len(proper_product),
        "proper_nonfull_base_product_certified": True,
        "non_product_size": len(non_product),
        "non_product_projection_unchanged": True,
        "non_product_fiber_sizes": {str(key): value for key, value in fiber_sizes.items()},
        "non_product_rejected_by_cardinality": True,
    }


def verify_case_certificate(row: dict) -> None:
    stored = row.get("decision_certificate_sha256")
    require(isinstance(stored, str) and len(stored) == 64, "missing case decision certificate")
    unsigned = {key: value for key, value in row.items() if key != "decision_certificate_sha256"}
    require(json_sha256(unsigned) == stored, "case decision certificate failed")


def verify_artifact(
    path: Path,
    payload: dict,
    context: dict[str, Any],
    context_audit: dict,
    f3_by_key: dict[str, dict],
) -> dict:
    experiment = payload.get("experiment")
    require(experiment in ENGINE_SPECS, "coordinate artifact has an unsupported experiment")
    spec = ENGINE_SPECS[experiment]
    require(payload.get("status") == spec["status"], "coordinate artifact is not complete")
    require(payload.get("p") == 7 and payload.get("z") == 7, "coordinate artifact problem changed")
    require(payload.get("phase") == 0 and payload.get("c_H") == 1, "coordinate artifact phase changed")
    require_true(payload.get("solver_invoked") is False, "coordinate artifact solver flag changed")

    source = payload["source_provenance"]
    require(
        source["CPU_semigroup_case_join_sha256"] == file_sha256(Path(cpu_join.__file__)),
        "coordinate artifact used another CPU construction script",
    )
    require(
        source["this_script_sha256"] == file_sha256(spec["script"]),
        "coordinate artifact engine source hash differs from this checkout",
    )
    if experiment.endswith("_opencl"):
        require(
            source["CUDA_semigroup_companion_sha256"] == file_sha256(Path(case_gpu.__file__)),
            "OpenCL artifact used another CPU companion construction",
        )

    artifact_audit = payload["CPU_construction_and_input_audits"]
    require(
        selected_fingerprint(artifact_audit) == selected_fingerprint(context_audit),
        "coordinate artifact reconstruction fingerprint differs from current context",
    )
    require_true(
        artifact_audit["all_51_grade_three_cases_reconstructed_before_sharding"],
        "artifact did not reconstruct all 51 cases before sharding",
    )

    configuration = payload["configuration"]
    k = configuration["k"]
    require(k in SUPPORTED_K, "coordinate artifact k is unsupported")
    prefix = list(range(k))
    require(
        configuration["explicit_mod7_projection_subsets"] == [prefix],
        "coordinate artifact is not the single width-k coordinate prefix",
    )
    require(
        configuration["exact_direction_support_cap"] == case_gpu.EXACT_DIRECTION_SUPPORT_CAP,
        "coordinate artifact direction-support cap changed",
    )
    expected_group_size = 3**F3_DIMENSION * 7**k
    require(configuration["state_cap"] == expected_group_size, "artifact did not permit the full group")

    shard = payload["case_sharding_audit"]
    require_true(shard["all_shards_form_a_disjoint_cover_of_the_51_cases"], "shard rule audit failed")
    require(shard["full_target_case_count"] == EXPECTED_CASES, "shard target census changed")
    shard_index = shard["shard_index"]
    shard_count = shard["shard_count"]
    require(isinstance(shard_count, int) and shard_count > 0, "invalid shard count")
    require(isinstance(shard_index, int) and 0 <= shard_index < shard_count, "invalid shard index")
    expected_indices = list(range(shard_index, EXPECTED_CASES, shard_count))
    require(shard["selected_target_indices"] == expected_indices, "artifact shard indices violate rule")
    expected_keys = [str(context["targets"][index]["case_key"]) for index in expected_indices]
    require(shard["selected_case_keys"] == expected_keys, "artifact shard keys violate audited order")
    require(shard["selected_case_count"] == len(expected_indices), "artifact shard count changed")
    require(json_sha256(expected_keys) == shard["selected_case_keys_sha256"], "shard key hash failed")

    rows = payload["case_results"]
    require(isinstance(rows, list) and len(rows) == len(expected_indices), "artifact result census changed")
    require(
        cpu_join.old_join.canonical_case_digest(rows) == payload["case_results_sha256"],
        "artifact case-result digest failed",
    )
    counts = payload["result_counts"]
    require(counts == {"selected": len(rows), "rejected": 0, "surviving": len(rows), "skipped": 0},
            "artifact is not an exact no-skip survivor shard")

    certified_rows = []
    seen_indices = set()
    for expected_index, row in zip(expected_indices, rows, strict=True):
        verify_case_certificate(row)
        index = row["target_index_in_audited_51_case_order"]
        key = str(row["case_key"])
        require(index == expected_index, "artifact case order differs from shard order")
        require(index not in seen_indices, "artifact repeats a case index")
        seen_indices.add(index)
        require(key == str(context["targets"][index]["case_key"]), "artifact case key/index mismatch")
        require_true(row["necessary_only_survivor"], "artifact row is not a survivor")
        require(row["rigorously_rejected"] is False and row["skipped"] is False,
                "artifact row is rejected or skipped")
        projections = row["projection_results"]
        require(isinstance(projections, list) and len(projections) == 1,
                "artifact row does not contain exactly one prefix projection")
        projection = projections[0]
        require(projection["mod7_coordinates"] == prefix, "artifact row prefix changed")
        require_true(projection["retained_all_six_mod3_coordinates"], "artifact dropped an F3 coordinate")
        require_true(projection["completed_exact_convolution"], "artifact convolution is incomplete")
        require_true(projection["projected_target_present"], "artifact target is absent")
        require(projection["rigorous_rejection"] is False, "artifact projection rejected its target")
        require_true(projection[spec["collision_flag"]], "artifact exact-set collision audit failed")
        require(projection["finite_group_size"] == expected_group_size, "artifact group size changed")
        require(projection["state_cap"] == expected_group_size, "projection state cap changed")
        require(projection["completed_catalog_factors"] == EXPECTED_DIRECTIONS,
                "artifact did not convolve all directions")
        require(len(projection["completed_state_sizes"]) == EXPECTED_DIRECTIONS,
                "artifact completion trace has the wrong length")
        final_size = projection["final_state_count"]
        require(final_size == projection["completed_state_sizes"][-1],
                "artifact final state count differs from completion trace")
        f3 = f3_by_key[key]
        fiber_size = 7**k
        require(
            product_cardinality_holds(f3["global_support_state_count"], fiber_size, final_size),
            "coordinate support failed the finite-fiber product cardinality",
        )
        require_true(f3["target_present"], "artifact target is not present in exact F3 support")
        require(
            projection["case_input_audit"]["target_digits"][:F3_DIMENSION]
            == f3["target_digits_F3_6"],
            "artifact and reconstructed F3 targets differ",
        )
        require(
            projection["case_input_audit"]["directional_excess_grades"]
            == f3["directional_excess_grades"],
            "artifact and reconstructed case grades differ",
        )
        certified_rows.append(
            {
                "target_index": index,
                "case_key": key,
                "F3_projection_state_count": f3["global_support_state_count"],
                "F7_fiber_size": fiber_size,
                "joint_support_state_count": final_size,
                "joint_equals_F3_count_times_F7_fiber_count": True,
                "certified_joint_support": f"pi_F3(S_{k}) x F7^{k}",
                "target_presence_depends_only_on_F3": True,
            }
        )
    require(seen_indices == set(expected_indices), "artifact result indices changed")
    return {
        "path": str(path.resolve()),
        "file_bytes": path.stat().st_size,
        "file_sha256": file_sha256(path),
        "experiment": experiment,
        "engine_source_sha256": source["this_script_sha256"],
        "k": k,
        "prefix_mod7_coordinates": prefix,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "selected_indices": expected_indices,
        "selected_case_keys_sha256": json_sha256(expected_keys),
        "case_count": len(certified_rows),
        "exact_completed_cases": len(certified_rows),
        "skipped_cases": 0,
        "all_cases_certify_coordinate_product": True,
        "case_product_records_sha256": json_sha256(certified_rows),
        "case_product_records": certified_rows,
    }


def audit_artifact_groups(records: list[dict]) -> dict:
    require(records, "at least one --coordinate-artifact is required")
    grouped: dict[int, list[dict]] = defaultdict(list)
    for record in records:
        grouped[record["k"]].append(record)
    summaries = []
    for k in sorted(grouped):
        rows = grouped[k]
        prefixes = {tuple(row["prefix_mod7_coordinates"]) for row in rows}
        require(prefixes == {tuple(range(k))}, f"k={k} artifacts do not share the exact prefix")
        all_indices = [index for row in rows for index in row["selected_indices"]]
        require(len(all_indices) == len(set(all_indices)), f"k={k} artifact shards overlap")
        require(all(0 <= index < EXPECTED_CASES for index in all_indices), "artifact index out of range")
        covered = sorted(all_indices)
        complete = covered == list(range(EXPECTED_CASES))
        missing = sorted(set(range(EXPECTED_CASES)) - set(covered))
        summaries.append(
            {
                "k": k,
                "prefix_mod7_coordinates": list(range(k)),
                "artifact_count": len(rows),
                "covered_case_count": len(covered),
                "covered_indices": covered,
                "missing_indices": missing,
                "shards_are_disjoint": True,
                "coverage_status": "complete_disjoint_51_case_cover" if complete else "partial_disjoint_cover",
                "full_51_case_coordinate_product_claim": complete,
                "partial_input_promoted_to_full_claim": False,
                "all_covered_cases_certify_coordinate_product": True,
            }
        )
    return {
        "widths_present": sorted(grouped),
        "groups": summaries,
        "k5_full_claim": next(
            (row["full_51_case_coordinate_product_claim"] for row in summaries if row["k"] == 5),
            False,
        ),
        "k6_full_claim": next(
            (row["full_51_case_coordinate_product_claim"] for row in summaries if row["k"] == 6),
            False,
        ),
    }


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    paths = [path.resolve() for path in args.coordinate_artifact]
    require(len(paths) == len(set(paths)), "the same coordinate artifact path was supplied twice")
    require(all(path.is_file() for path in paths), "a coordinate artifact path is not a file")
    context, context_audit = build_context(
        args.parent_input, args.current_join, args.hilbert_basis
    )
    support_table, f3_support_audit = build_f3_supports(context)
    f3_cases, f3_by_key = exact_f3_cases(context, support_table, f3_support_audit)
    expected_fingerprint = selected_fingerprint(context_audit)

    artifact_records = []
    for path in paths:
        raw = path.read_bytes()
        try:
            payload = json.loads(raw)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise AssertionError(f"coordinate artifact is not valid JSON: {path}") from error
        require(isinstance(payload, dict), "coordinate artifact root is not an object")
        artifact_records.append(
            verify_artifact(path, payload, context, context_audit, f3_by_key)
        )
    group_audit = audit_artifact_groups(artifact_records)
    script_path = Path(__file__).resolve()
    result = {
        "experiment": EXPERIMENT,
        "status": "complete_exact_F3_coordinate_product_audit",
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "51 orbit0/branch-A positive z=7 grade-three cases",
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "CPU_semigroup_case_join_sha256": file_sha256(Path(cpu_join.__file__)),
            "GPU_semigroup_companion_sha256": file_sha256(Path(case_gpu.__file__)),
            "reconstructed_context_fingerprint": expected_fingerprint,
            "reconstructed_context_fingerprint_sha256": json_sha256(expected_fingerprint),
        },
        "configuration": {
            "parent_input": str(args.parent_input.resolve()),
            "current_join": str(args.current_join.resolve()),
            "hilbert_basis": str(args.hilbert_basis.resolve()),
            "coordinate_artifacts": [str(path) for path in paths],
            "output": str(args.output.resolve()),
        },
        "manufactured_finite_fiber_audit": manufactured_fiber_audit(),
        "finite_fiber_cardinality_lemma": {
            "statement": (
                "For finite sets S subset T x B, every fiber S_t has size at most |B|. "
                "Therefore |S| <= |T||B|, and equality forces every fiber to equal B; "
                "hence S = T x B and pi_T(S) = T."
            ),
            "application": (
                "Exact same-row construction gives S_k subset T x F7^k, where T is the "
                "independently reconstructed exact F3 global support.  The audited equality "
                "|S_k|=|T|7^k therefore certifies S_k=T x F7^k."
            ),
            "target_consequence": (
                "After product certification, (t3,t7) lies in S_k exactly when t3 lies in T; "
                "the retained F7 prefix cannot affect target presence."
            ),
        },
        "exact_F3_direction_and_structure_audit": f3_support_audit,
        "exact_F3_case_count": len(f3_cases),
        "all_51_F3_targets_present": all(row["target_present"] for row in f3_cases),
        "F3_case_records_sha256": json_sha256(f3_cases),
        "F3_case_records": f3_cases,
        "coordinate_artifact_count": len(artifact_records),
        "coordinate_artifact_records_sha256": json_sha256(artifact_records),
        "coordinate_artifact_records": artifact_records,
        "coordinate_width_cover_audit": group_audit,
        "logical_semantics": {
            "F3_directional_supports_are_exact_complete_direct_catalog_projections": True,
            "F3_global_supports_are_exact_eight_direction_convolutions": True,
            "each_direction_F3_image_rank_is_one": True,
            "total_F3_image_rank_is_six": True,
            "declared_direction_split_is_a_3_plus_3_direct_sum": True,
            "F3_scalar_map_kernel_dimension_is_two": True,
            "each_target_has_at_most_nine_scalar_bucket_patterns": True,
            "artifact_cardinalities_are_used_only_after_provenance_and_completion_audits": True,
            "coordinate_product_is_claimed_only_for_cases_present_in_disjoint_exact_artifacts": True,
            "full_width_k_claim_requires_disjoint_cover_of_all_51_cases": True,
            "partial_k6_input_is_labeled_partial": True,
            "target_presence_in_a_certified_product_depends_only_on_F3": True,
            "binary_edge_feasibility_claimed": False,
            "positive_z7_closure_claimed": False,
            "full_theorem_claimed": False,
        },
        "elapsed_seconds": time.time() - started,
    }
    result["certificate_payload_sha256_before_this_field"] = json_sha256(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-input", type=Path, default=cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument(
        "--coordinate-artifact",
        type=Path,
        action="append",
        required=True,
        help="completed coordinate-prefix k=5 or k=6 shard JSON; repeat for each shard",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args)
    atomic_write_json(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output.resolve()),
                "script_sha256": result["source_provenance"]["this_script_sha256"],
                "F3_cases": result["exact_F3_case_count"],
                "all_F3_targets_present": result["all_51_F3_targets_present"],
                "coordinate_width_cover_audit": result["coordinate_width_cover_audit"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
