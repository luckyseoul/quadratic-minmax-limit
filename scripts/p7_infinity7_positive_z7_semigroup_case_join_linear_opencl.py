#!/usr/bin/env python3
"""Exact OpenCL joins under full-rank linear maps F7^21 -> F7^k.

This hybrid preserves the audited arbitrary-map construction, explicit matrix
catalog, modular-rank certificates, and same-source-row projection path from
``p7_infinity7_positive_z7_semigroup_case_join_linear_gpu.py``.  It sends the
resulting F3^6 x F7^k supports to the exact collision-resolved OpenCL engine in
``p7_infinity7_positive_z7_semigroup_case_join_opencl.py``.

All six F3 coordinates are retained verbatim.  The same explicit full-row-rank
matrix is applied to the 21 F7 coordinates of every Hilbert generator, direct
catalog row, and target, before deduplication.  Matrix and case sharding are
independent deterministic covers.  A missing target after a completed full
map convolution is a rigorous rejection; target presence is necessary only.
State, memory, allocation, and hash-cap failures are explicit nondecisions.

The device audit inherits exact uint64 compare-and-exchange with linear-probe
collision resolution, forced collision/present/absent/cross-prime/cap traps,
and k=5/k=6 manufactured CPU comparisons.  The hybrid audit additionally
proves selector-matrix equivalence to the coordinate OpenCL path on all 32
real direction/grade supports and compares one complete small real CPU/OpenCL
support set elementwise.  JSON output is replaced atomically.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

import sys

sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join_linear_gpu as linear  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_join_opencl as opencl  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_join_linear_opencl"
SUPPORTED_K = linear.SUPPORTED_K
DEFAULT_MATRIX_FAMILIES = linear.DEFAULT_MATRIX_FAMILIES
DEFAULT_PAIR_CHUNK_CAP = opencl.DEFAULT_PAIR_CHUNK_CAP
DEFAULT_OPENCL_MEMORY_CAP_MIB = opencl.DEFAULT_OPENCL_MEMORY_CAP_MIB
EXPECTED_TARGET_CASES = linear.EXPECTED_TARGET_CASES


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    return opencl.file_sha256(path)


def json_sha256(value: object) -> str:
    return opencl.json_sha256(value)


def selector_coordinate_opencl_crosscheck(
    runtime: dict[str, Any],
    opencl_memory_cap_bytes: int,
    context: dict[str, Any],
    pair_chunk_cap: int,
) -> dict:
    """Prove selector/coordinate equivalence and one real CPU/OpenCL set equality."""
    projection = linear.SELECTOR_CROSSCHECK_COORDINATES
    selector = linear.selector_matrix(len(projection), projection)
    selector_record = linear.matrix_record(
        name="real_selector_0_1",
        family="selector",
        matrix=selector,
        construction={
            "strategy": "coordinate_selector_crosscheck",
            "selected_coordinates": list(projection),
        },
    )

    coordinate_table, coordinate_direction_audit = opencl.build_direction_supports(
        context, projection, pair_chunk_cap, direction_workers=1
    )
    linear_table, linear_direction_audit = linear.build_linear_direction_supports(
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
                    "support_sha256_uint64": linear.cpu_join.array_sha256(
                        linear_codes.astype("<u8", copy=False)
                    ),
                }
            )

    target_row = context["targets"][0]
    (
        linear_rows,
        linear_target,
        linear_moduli,
        linear_ordered,
        linear_case_audit,
    ) = linear.exact_linear_case_inputs(
        context, target_row, selector_record, linear_table
    )
    (
        coordinate_rows,
        coordinate_target,
        coordinate_moduli,
        coordinate_ordered,
        coordinate_case_audit,
    ) = linear.coordinate_gpu.exact_case_inputs(
        context, target_row, projection, coordinate_table
    )
    require(linear_moduli == coordinate_moduli, "selector moduli differ")
    require(
        np.array_equal(linear_target, coordinate_target),
        "selector targets differ",
    )
    require(
        [direction for direction, _codes in linear_ordered]
        == [direction for direction, _codes in coordinate_ordered],
        "selector convolution orders differ",
    )
    require(
        all(
            np.array_equal(left_codes, right_codes)
            for (_left_direction, left_codes), (_right_direction, right_codes)
            in zip(linear_ordered, coordinate_ordered)
        ),
        "selector case factors differ",
    )
    require(
        len(linear_rows) == len(coordinate_rows)
        and all(
            np.array_equal(left_rows, right_rows)
            for left_rows, right_rows in zip(linear_rows, coordinate_rows)
        ),
        "selector OpenCL interface rows differ",
    )

    codec = linear.cpu_join.semigroup.MixedRadixCodec(linear_moduli)
    cpu_codes, cpu_convolution = linear.cpu_join.convolve_support_sequence(
        linear_ordered,
        codec,
        state_cap=codec.group_size,
        pair_chunk_cap=min(pair_chunk_cap, codec.group_size),
    )
    require(
        cpu_codes is not None and cpu_convolution["completed"],
        "small real CPU crosscheck capped",
    )
    decision, opencl_codes = opencl.opencl_exact_support_convolution(
        runtime,
        linear_rows,
        linear_target,
        linear_moduli,
        state_cap=codec.group_size,
        pair_chunk_cap=pair_chunk_cap,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    opencl.validate_opencl_decision(decision)
    require(opencl_codes is not None, "small real OpenCL crosscheck skipped")
    require(
        np.array_equal(cpu_codes, opencl_codes),
        "linear-selector CPU/OpenCL support sets differ",
    )
    target_code = int(codec.encode(linear_target[None, :])[0])
    target_position = int(np.searchsorted(cpu_codes, np.uint64(target_code)))
    cpu_present = bool(
        target_position < len(cpu_codes)
        and int(cpu_codes[target_position]) == target_code
    )
    require(
        cpu_present == decision["projected_target_present"],
        "linear-selector CPU/OpenCL target decisions differ",
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
        "OpenCL_convolution_audit": decision,
        "complete_support_state_count": len(cpu_codes),
        "complete_support_sha256_uint64": linear.cpu_join.array_sha256(
            cpu_codes.astype("<u8", copy=False)
        ),
        "target_present": cpu_present,
        "all_32_selector_matrix_and_coordinate_supports_equal_elementwise": True,
        "selector_matrix_and_coordinate_target_equal_elementwise": True,
        "selector_matrix_and_coordinate_case_factors_equal_elementwise": True,
        "selector_matrix_and_coordinate_OpenCL_inputs_equal_elementwise": True,
        "small_real_CPU_and_OpenCL_complete_support_sets_equal_elementwise": True,
        "small_real_CPU_and_OpenCL_target_decisions_equal": True,
    }


def run_production(
    *,
    runtime: dict[str, Any],
    context: dict[str, Any],
    context_audit: dict,
    manufactured_opencl: dict,
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
    opencl_memory_cap_bytes: int,
    output_path: Path,
) -> dict:
    started = time.time()
    support_tables: dict[str, dict[int, dict[int, np.ndarray]]] = {}
    support_audits = []
    for matrix_index, matrix_row in selected_matrices:
        table, audit = linear.build_linear_direction_supports(
            context, matrix_row, pair_chunk_cap
        )
        matrix_hash = str(matrix_row["matrix_sha256_uint8"])
        require(matrix_hash not in support_tables, "selected matrix hash repeated")
        support_tables[matrix_hash] = table
        support_audits.append(
            {"matrix_index_in_audited_catalog": matrix_index, **audit}
        )

    projection_totals: Counter[str] = Counter()
    case_results = []
    for target_index, target_row in selected_cases:
        projection_results = []
        for matrix_index, matrix_row in selected_matrices:
            matrix_hash = str(matrix_row["matrix_sha256_uint8"])
            support_rows, target, moduli, _ordered, case_input_audit = (
                linear.exact_linear_case_inputs(
                    context,
                    target_row,
                    matrix_row,
                    support_tables[matrix_hash],
                )
            )
            decision, _final_codes = opencl.opencl_exact_support_convolution(
                runtime,
                support_rows,
                target,
                moduli,
                state_cap=state_cap,
                pair_chunk_cap=pair_chunk_cap,
                opencl_memory_cap_bytes=opencl_memory_cap_bytes,
            )
            opencl.validate_opencl_decision(decision)
            projection_totals[str(decision["decision_status"])] += 1
            projection_results.append(
                {
                    "matrix_index_in_audited_catalog": matrix_index,
                    "matrix_name": matrix_row["name"],
                    "matrix_family": matrix_row["family"],
                    "matrix_sha256_uint8": matrix_hash,
                    "matrix_rank_mod7": matrix_row["rank_mod7"],
                    "retained_all_six_F3_coordinates_verbatim": True,
                    "case_input_audit": case_input_audit,
                    **decision,
                }
            )

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
                "rigorous_exact_OpenCL_linear_semigroup_projection_rejection"
                if rejected
                else "explicit_cap_skip_without_negative_decision"
                if skipped
                else "necessary_only_survivor_of_all_completed_linear_OpenCL_projections"
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
        sum(counts[key] for key in ("rejected", "surviving", "skipped"))
        == counts["selected"],
        "case result census is not a partition",
    )
    k = int(matrix_catalog[0]["shape"][0])
    script_path = Path(__file__).resolve()
    return {
        "experiment": EXPERIMENT,
        "status": (
            "complete_sharded_exact_OpenCL_linear_semigroup_case_join"
            if counts["skipped"] == 0
            else "exact_OpenCL_linear_semigroup_case_join_with_explicit_skips"
        ),
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "51 orbit0/branch-A grade-three-only representatives",
        "opencl_invoked": True,
        "solver_invoked": False,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "linear_CUDA_companion_path": str(Path(linear.__file__).resolve()),
            "linear_CUDA_companion_sha256": file_sha256(Path(linear.__file__)),
            "coordinate_OpenCL_companion_path": str(Path(opencl.__file__).resolve()),
            "coordinate_OpenCL_companion_sha256": file_sha256(Path(opencl.__file__)),
            "CPU_semigroup_case_join_path": str(Path(linear.cpu_join.__file__).resolve()),
            "CPU_semigroup_case_join_sha256": file_sha256(
                Path(linear.cpu_join.__file__)
            ),
        },
        "configuration": {
            "k": k,
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "opencl_memory_cap_bytes": opencl_memory_cap_bytes,
            "exact_direction_support_cap": linear.EXACT_DIRECTION_SUPPORT_CAP,
            "selected_matrix_indices": [index for index, _row in selected_matrices],
            "selected_matrix_hashes": [
                row["matrix_sha256_uint8"] for _index, row in selected_matrices
            ],
        },
        "expected_memory_per_linear_projection": opencl.expected_engine_memory(
            k, state_cap, pair_chunk_cap
        ),
        "opencl_engine_audit": runtime["audit"],
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_CPU_OpenCL_audit": manufactured_opencl,
        "manufactured_linear_projection_audit": manufactured_linear,
        "selector_matrix_coordinate_OpenCL_crosscheck": selector_crosscheck,
        "selector_crosscheck_skipped_by_explicit_flag": selector_crosscheck is None,
        "full_linear_matrix_catalog_audit": matrix_catalog_audit,
        "full_linear_matrix_catalog": matrix_catalog,
        "matrix_sharding_audit": matrix_shard_audit,
        "case_sharding_audit": case_shard_audit,
        "linear_direction_support_audits": support_audits,
        "result_counts": counts,
        "projection_decision_census": dict(sorted(projection_totals.items())),
        "case_results_sha256": linear.cpu_join.old_join.canonical_case_digest(
            case_results
        ),
        "case_results": case_results,
        "logical_semantics": {
            "all_51_cases_reconstructed_before_case_sharding": True,
            "full_matrix_catalog_rank_audited_before_matrix_sharding": True,
            "all_six_F3_coordinates_retained_verbatim": True,
            "every_F7_map_is_an_explicit_hash_pinned_full_rank_k_by_21_matrix": True,
            "same_Hilbert_or_catalog_row_supplies_F3_and_F7_before_linear_map": True,
            "identical_linear_matrix_maps_generators_catalogs_and_target": True,
            "grade_zero_through_three_supports_are_exact_and_direct_catalog_calibrated": True,
            "all_eight_direction_supports_convolved_in_every_completed_map": True,
            "OpenCL_uint64_hash_collisions_resolved_by_exact_key_comparison": True,
            "cl_khr_int64_base_atomics_required_and_audited": True,
            "exact_full_group_saturation_only_proves_target_presence": True,
            "missing_target_after_completed_full_map_convolution_is_rigorous": True,
            "target_presence_is_necessary_only": True,
            "all_state_memory_allocation_or_hash_caps_are_explicit_skips": True,
            "partial_support_after_any_cap_is_never_used": True,
            "probabilistic_membership_used": False,
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
            "linear_CUDA_companion_path": str(Path(linear.__file__).resolve()),
            "linear_CUDA_companion_sha256": file_sha256(Path(linear.__file__)),
            "coordinate_OpenCL_companion_path": str(Path(opencl.__file__).resolve()),
            "coordinate_OpenCL_companion_sha256": file_sha256(Path(opencl.__file__)),
        },
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_linear_projection_audit": manufactured_linear,
        "full_linear_matrix_catalog_audit": matrix_catalog_audit,
        "full_linear_matrix_catalog": matrix_catalog,
        "opencl_loaded": False,
        "production_cases_processed": 0,
        "output_path": str(output_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform-index", type=int, default=0)
    parser.add_argument("--device-index", type=int, default=0)
    parser.add_argument("--device-name-contains")
    parser.add_argument(
        "--parent-input", type=Path, default=linear.cpu_join.DEFAULT_PARENT_INPUT
    )
    parser.add_argument(
        "--current-join", type=Path, default=linear.cpu_join.DEFAULT_CURRENT_JOIN
    )
    parser.add_argument(
        "--hilbert-basis", type=Path, default=linear.cpu_join.DEFAULT_HILBERT_BASIS
    )
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
        "--opencl-memory-cap-mib", type=int, default=DEFAULT_OPENCL_MEMORY_CAP_MIB
    )
    parser.add_argument("--skip-small-real-crosscheck", action="store_true")
    parser.add_argument("--manufactured-only", action="store_true")
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--matrix-catalog-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    require(args.opencl_memory_cap_mib > 0, "--opencl-memory-cap-mib must be positive")
    require(args.state_cap >= 0, "--state-cap cannot be negative")
    require(
        sum((args.manufactured_only, args.self_audit_only, args.matrix_catalog_only))
        <= 1,
        "manufactured-only, self-audit-only, and matrix-catalog-only are mutually exclusive",
    )
    families = linear.parse_matrix_families(args.matrix_families)
    manufactured_linear = linear.manufactured_linear_projection_audit()

    context = None
    context_audit = None
    matrix_catalog = None
    matrix_catalog_audit = None
    if not args.manufactured_only:
        context, context_audit = linear.coordinate_gpu.build_cpu_context(
            args.parent_input, args.current_join, args.hilbert_basis
        )
        require(len(context["targets"]) == EXPECTED_TARGET_CASES, "target census changed")
        matrix_catalog, matrix_catalog_audit = linear.build_matrix_catalog(
            context, args.k, families
        )

    if args.matrix_catalog_only:
        require(
            context_audit is not None
            and matrix_catalog is not None
            and matrix_catalog_audit is not None,
            "matrix catalog context missing",
        )
        result = catalog_only_result(
            context_audit=context_audit,
            matrix_catalog=matrix_catalog,
            matrix_catalog_audit=matrix_catalog_audit,
            manufactured_linear=manufactured_linear,
            output_path=args.output,
        )
        opencl.cuda_join.atomic_write(args.output, result)
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

    group_size = opencl.cuda_join.group_size_for((3,) * 6 + (7,) * args.k)
    state_cap = group_size if args.state_cap == 0 else min(args.state_cap, group_size)
    opencl_memory_cap_bytes = args.opencl_memory_cap_mib * (1 << 20)
    runtime = opencl.load_opencl(
        args.platform_index, args.device_index, args.device_name_contains
    )
    manufactured_opencl = opencl.manufactured_cross_engine_audit(
        runtime, opencl_memory_cap_bytes
    )

    if args.manufactured_only:
        result = {
            "experiment": f"{EXPERIMENT}_manufactured_self_audit",
            "status": "manufactured_linear_and_CPU_OpenCL_self_audits_passed",
            "source_provenance": {
                "this_script_path": str(Path(__file__).resolve()),
                "this_script_sha256": file_sha256(Path(__file__).resolve()),
                "linear_CUDA_companion_sha256": file_sha256(Path(linear.__file__)),
                "coordinate_OpenCL_companion_sha256": file_sha256(Path(opencl.__file__)),
            },
            "opencl_engine_audit": runtime["audit"],
            "manufactured_CPU_OpenCL_audit": manufactured_opencl,
            "manufactured_linear_projection_audit": manufactured_linear,
            "selector_real_crosscheck_run": False,
            "production_cases_processed": 0,
            "output_path": str(args.output.resolve()),
        }
    else:
        require(
            context is not None
            and context_audit is not None
            and matrix_catalog is not None
            and matrix_catalog_audit is not None,
            "production context or matrix catalog missing",
        )
        selector_crosscheck = (
            None
            if args.skip_small_real_crosscheck and not args.self_audit_only
            else selector_coordinate_opencl_crosscheck(
                runtime,
                opencl_memory_cap_bytes,
                context,
                args.pair_chunk_cap,
            )
        )
        if args.self_audit_only:
            require(selector_crosscheck is not None, "self-audit omitted real crosscheck")
            result = {
                "experiment": f"{EXPERIMENT}_self_audit",
                "status": "manufactured_and_real_linear_CPU_OpenCL_self_audits_passed",
                "source_provenance": {
                    "this_script_path": str(Path(__file__).resolve()),
                    "this_script_sha256": file_sha256(Path(__file__).resolve()),
                    "linear_CUDA_companion_sha256": file_sha256(Path(linear.__file__)),
                    "coordinate_OpenCL_companion_sha256": file_sha256(Path(opencl.__file__)),
                },
                "opencl_engine_audit": runtime["audit"],
                "CPU_construction_and_input_audits": context_audit,
                "manufactured_CPU_OpenCL_audit": manufactured_opencl,
                "manufactured_linear_projection_audit": manufactured_linear,
                "selector_matrix_coordinate_OpenCL_crosscheck": selector_crosscheck,
                "full_linear_matrix_catalog_audit": matrix_catalog_audit,
                "full_linear_matrix_catalog": matrix_catalog,
                "production_cases_processed": 0,
                "output_path": str(args.output.resolve()),
            }
        else:
            selected_matrices, matrix_shard_audit = linear.shard_indexed_rows(
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
            selected_cases, case_shard_audit = linear.coordinate_gpu.shard_cases(
                context["targets"], args.case_shard_index, args.case_shard_count
            )
            case_shard_audit["label"] = "case"
            result = run_production(
                runtime=runtime,
                context=context,
                context_audit=context_audit,
                manufactured_opencl=manufactured_opencl,
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
                opencl_memory_cap_bytes=opencl_memory_cap_bytes,
                output_path=args.output,
            )

    opencl.cuda_join.atomic_write(args.output, result)
    summary = {
        "status": result["status"],
        "output": str(args.output),
        "opencl_device": runtime["audit"]["device_name"],
        "k": args.k,
        "manufactured_CPU_OpenCL_audit": manufactured_opencl["status"],
        "selector_matrix_coordinate_OpenCL_crosscheck": (
            None
            if args.manufactured_only or result.get("selector_matrix_coordinate_OpenCL_crosscheck") is None
            else result["selector_matrix_coordinate_OpenCL_crosscheck"]["status"]
        ),
        "processed_cases": result.get("result_counts", {}).get("selected", 0),
        "rigorously_rejected_cases": result.get("result_counts", {}).get("rejected", 0),
        "necessary_only_survivors": result.get("result_counts", {}).get("surviving", 0),
        "skipped_cases": result.get("result_counts", {}).get("skipped", 0),
        "selected_matrices": result.get("matrix_sharding_audit", {}).get(
            "selected_matrix_names", []
        ),
        "expected_memory": opencl.expected_engine_memory(
            args.k, state_cap, args.pair_chunk_cap
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
