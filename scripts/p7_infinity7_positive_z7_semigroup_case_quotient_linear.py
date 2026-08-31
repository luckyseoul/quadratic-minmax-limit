#!/usr/bin/env python3
"""Exact CPU stabilizer-quotient joins under linear maps F7^21 -> F7^k.

This standalone driver composes the audited full-rank matrix catalog from
``p7_infinity7_positive_z7_semigroup_case_join_linear_gpu`` with the exact
translation-stabilizer and quotient engine from
``p7_infinity7_positive_z7_semigroup_case_quotient``.  Every map retains all
six F3 coordinates and applies one explicit full-row-rank matrix to the 21 F7
coordinates.  The same rowwise map is used for Hilbert generators, complete
direct catalogs, and targets before any deduplication.

For each selected matrix and case, all eight factor supports are retained.
Only the sum of exact H3 translation stabilizers is quotiented out.  A missing
target after a completed exact quotient join is a rigorous rejection; target
presence is necessary only.  State-cap hits are explicit skips and partial
supports are never used.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join_linear_gpu as linear  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_quotient as quotient  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_quotient_linear"
SUPPORTED_K = (5, 6)
DEFAULT_STATE_CAP = 2_000_000
DEFAULT_PAIR_CHUNK_CAP = 200_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def transformed_quotient_data(
    source: dict[int, dict], matrix_record: dict
) -> tuple[dict[int, dict], dict]:
    """Replace the F7 complement by M times it, retaining the F3 complement."""
    matrix = linear.matrix_from_record(matrix_record).astype(np.int64)
    original = np.ascontiguousarray(source[7]["complement"], dtype=np.int64)
    require(original.shape[0] == 21, "source F7 quotient does not have dimension 21")
    mapped = np.ascontiguousarray(matrix @ original % 7, dtype=np.int64)
    transformed = {prime: dict(data) for prime, data in source.items()}
    transformed[3]["complement"] = np.ascontiguousarray(
        source[3]["complement"], dtype=np.int64
    )
    transformed[7]["complement"] = mapped
    require(
        len(transformed[3]["complement"]) == 6 and len(mapped) == len(matrix),
        "transformed quotient dimensions changed",
    )
    return transformed, {
        "matrix_name": matrix_record["name"],
        "matrix_family": matrix_record["family"],
        "matrix_sha256_uint8": matrix_record["matrix_sha256_uint8"],
        "matrix_rank_mod7": matrix_record["rank_mod7"],
        "rank_certificate": {
            "pivot_columns_zero_based": matrix_record["pivot_columns_zero_based"],
            "pivot_minor_determinant_mod7": matrix_record[
                "pivot_minor_determinant_mod7"
            ],
            "exact_full_row_rank_certificate": matrix_record[
                "exact_full_row_rank_certificate"
            ],
        },
        "source_F7_dimension": 21,
        "mapped_F7_dimension": len(matrix),
        "transformed_F7_complement_sha256_int64": quotient.array_sha256(mapped),
        "all_six_F3_complement_rows_retained_elementwise": bool(
            np.array_equal(
                transformed[3]["complement"], source[3]["complement"]
            )
        ),
        "construction": "mapped_F7_complement_equals_M_times_full_F7_complement_mod7",
    }


def matrix_mapping_audit(
    context: dict[str, Any], matrix_record: dict
) -> tuple[dict[int, dict], dict[int, np.ndarray], dict]:
    """Prove complement composition equals direct rowwise linear projection."""
    mapped_data, transform_audit = transformed_quotient_data(
        context["quotient_data"], matrix_record
    )
    mapped_rows, generator_audit = quotient.naive_join.full_generator_projections(
        context["basis"], context["orbit"], mapped_data
    )
    matrix = linear.matrix_from_record(matrix_record)
    direction_rows = []
    for direction in range(8):
        expected = linear.apply_linear_projection(
            context["full_generator_rows"][direction], matrix
        )
        require(
            np.array_equal(mapped_rows[direction], expected),
            f"direction {direction} complement composition differs from linear path",
        )
        direction_rows.append(
            {
                "direction": direction,
                "row_count": len(expected),
                "mapped_rows_sha256_uint8": quotient.array_sha256(expected),
            }
        )

    target = context["targets"][0]
    orbit, leaf, system, _factory = quotient.naive_join.old_join.validate_parent_survivor(
        target, context["rebuilt"]
    )
    anchor_rhs, _raw = quotient.naive_join.affine.anchor_rhs_and_raw_syndromes(
        orbit, leaf, system, context["rebuilt"]["anchors"]
    )
    composed_target = quotient.naive_join.project_equation_vector(
        anchor_rhs, mapped_data, tuple(range(len(matrix)))
    )
    direct_target = linear.project_equation_vector_linear(
        anchor_rhs, context["quotient_data"], matrix
    )
    require(
        np.array_equal(composed_target, direct_target),
        "composed and direct linear target maps differ",
    )
    return mapped_data, mapped_rows, {
        **transform_audit,
        "generator_projection_audit": generator_audit,
        "direction_generator_equivalence": direction_rows,
        "first_case_target_sha256_uint8": quotient.array_sha256(direct_target),
        "all_eight_generator_tables_equal_direct_rowwise_linear_projection": True,
        "first_case_target_equals_direct_rowwise_linear_projection": True,
        "same_composed_map_is_used_by_imported_support_and_target_paths": True,
    }


def build_stabilizer_payload(
    *,
    context: dict[str, Any],
    matrix_record: dict,
    state_cap: int,
    pair_chunk_cap: int,
) -> tuple[dict | None, dict]:
    mapped_data, mapped_rows, mapping_audit = matrix_mapping_audit(
        context, matrix_record
    )
    k = int(matrix_record["shape"][0])
    payload, support_audit = quotient.build_projection_supports(
        projection=tuple(range(k)),
        full_generator_rows=mapped_rows,
        degrees=context["degrees"],
        rebuilt=context["rebuilt"],
        quotient_data=mapped_data,
        state_cap=state_cap,
        pair_chunk_cap=pair_chunk_cap,
    )
    if payload is not None:
        payload["quotient_data"] = mapped_data
    return payload, {
        "linear_matrix": matrix_record,
        "linear_mapping_equivalence_audit": mapping_audit,
        "stabilizer_support_audit": support_audit,
    }


def real_selector_equivalence_audit(
    context: dict[str, Any], k: int, pair_chunk_cap: int
) -> dict:
    """Compare all 32 exact supports and one target with coordinate selection."""
    projection = tuple(range(k))
    selector_record = linear.matrix_record(
        name=f"selector_prefix_{k}_CPU_equivalence",
        family="selector",
        matrix=linear.selector_matrix(k, projection),
        construction={
            "strategy": "coordinate_selector_CPU_crosscheck",
            "selected_coordinates": list(projection),
        },
    )
    coordinate_table, coordinate_audit = (
        linear.coordinate_gpu.build_direction_supports(
            context, projection, pair_chunk_cap
        )
    )
    linear_table, linear_audit = linear.build_linear_direction_supports(
        context, selector_record, pair_chunk_cap
    )
    support_rows = []
    for direction in range(8):
        for grade in range(4):
            left = coordinate_table[direction][grade]
            right = linear_table[direction][grade]
            require(
                np.array_equal(left, right),
                f"selector/coordinate support differs at direction {direction}, grade {grade}",
            )
            support_rows.append(
                {
                    "direction": direction,
                    "grade": grade,
                    "states": len(left),
                    "sha256_uint64": quotient.array_sha256(
                        left.astype("<u8", copy=False)
                    ),
                }
            )

    target = context["targets"][0]
    linear_inputs = linear.exact_linear_case_inputs(
        context, target, selector_record, linear_table
    )
    coordinate_inputs = linear.coordinate_gpu.exact_case_inputs(
        context, target, projection, coordinate_table
    )
    linear_rows, linear_target, linear_moduli, linear_order, _linear_case = linear_inputs
    coordinate_rows, coordinate_target, coordinate_moduli, coordinate_order, _coord_case = (
        coordinate_inputs
    )
    require(linear_moduli == coordinate_moduli, "selector moduli differ")
    require(np.array_equal(linear_target, coordinate_target), "selector target differs")
    require(
        [row[0] for row in linear_order] == [row[0] for row in coordinate_order],
        "selector factor order differs",
    )
    require(
        all(np.array_equal(left, right) for left, right in zip(linear_rows, coordinate_rows)),
        "selector factor rows differ",
    )
    return {
        "status": "passed",
        "k": k,
        "selector_matrix": selector_record,
        "coordinate_projection": list(projection),
        "coordinate_direction_support_audit": coordinate_audit,
        "linear_direction_support_audit": linear_audit,
        "support_comparisons": support_rows,
        "support_comparisons_sha256": quotient.json_sha256(support_rows),
        "first_case_target_sha256_uint8": quotient.array_sha256(linear_target),
        "all_32_selector_and_coordinate_supports_equal_elementwise": True,
        "selector_and_coordinate_target_equal_elementwise": True,
        "selector_and_coordinate_factor_rows_equal_elementwise": True,
        "all_six_F3_coordinates_preserved": True,
    }


def selected_cases(
    targets: list[dict],
    *,
    all_cases: bool,
    case_count: int,
    shard_index: int,
    shard_count: int,
) -> tuple[list[tuple[int, dict]], dict]:
    require(shard_count > 0, "--case-shard-count must be positive")
    require(0 <= shard_index < shard_count, "--case-shard-index is outside count")
    require(case_count > 0, "--case-count must be positive")
    universe = list(enumerate(targets if all_cases else targets[:case_count]))
    selected = [row for row in universe if row[0] % shard_count == shard_index]
    require(selected, "selected case shard is empty")
    flattened = [
        index
        for shard in range(shard_count)
        for index, _row in universe
        if index % shard_count == shard
    ]
    require(
        sorted(flattened) == [index for index, _row in universe]
        and len(flattened) == len(set(flattened)),
        "case shards do not form a disjoint cover",
    )
    keys = [str(row[1]["case_key"]) for row in selected]
    return selected, {
        "selection_mode": "all_51_cases" if all_cases else "audited_order_prefix",
        "unsharded_case_count": len(universe),
        "case_shard_index": shard_index,
        "case_shard_count": shard_count,
        "selected_case_count": len(selected),
        "selected_target_indices": [row[0] for row in selected],
        "selected_case_keys": keys,
        "selected_case_keys_sha256": quotient.json_sha256(keys),
        "all_case_shards_form_a_disjoint_cover": True,
    }


def run(
    *,
    context: dict[str, Any],
    context_audit: dict,
    matrix_catalog: list[dict],
    matrix_catalog_audit: dict,
    selected_matrices: list[tuple[int, dict]],
    matrix_shard_audit: dict,
    cases: list[tuple[int, dict]],
    case_shard_audit: dict,
    state_cap: int,
    pair_chunk_cap: int,
    selector_audit: dict,
    manufactured_audits: dict,
    output_path: Path,
) -> dict:
    started = time.time()
    usage_started = resource.getrusage(resource.RUSAGE_SELF)
    payloads: dict[str, dict | None] = {}
    matrix_support_audits = []
    for matrix_index, matrix_record in selected_matrices:
        payload, audit = build_stabilizer_payload(
            context=context,
            matrix_record=matrix_record,
            state_cap=state_cap,
            pair_chunk_cap=pair_chunk_cap,
        )
        matrix_hash = str(matrix_record["matrix_sha256_uint8"])
        payloads[matrix_hash] = payload
        matrix_support_audits.append(
            {"matrix_index_in_audited_catalog": matrix_index, **audit}
        )

    decision_census: Counter[str] = Counter()
    case_results = []
    for target_index, target in cases:
        key = str(target["case_key"])
        projection_results = []
        for matrix_index, matrix_record in selected_matrices:
            matrix_hash = str(matrix_record["matrix_sha256_uint8"])
            payload = payloads[matrix_hash]
            if payload is None:
                decision = {
                    "decision_status": "skipped_linear_support_or_stabilizer_state_cap",
                    "rigorously_rejected": False,
                    "necessary_only": False,
                    "skipped": True,
                    "exact_quotient_membership_completed": False,
                    "partial_support_used": False,
                }
            else:
                k = int(matrix_record["shape"][0])
                decision = quotient.evaluate_quotient_case(
                    target_row=target,
                    current_row=context["current_by_key"][key],
                    rebuilt=context["rebuilt"],
                    common=context["common"],
                    quotient_data=payload["quotient_data"],
                    projection=tuple(range(k)),
                    projection_data=payload,
                    state_cap=state_cap,
                    pair_chunk_cap=pair_chunk_cap,
                )
            decision_census[str(decision["decision_status"])] += 1
            projection_results.append(
                {
                    "matrix_index_in_audited_catalog": matrix_index,
                    "matrix_name": matrix_record["name"],
                    "matrix_family": matrix_record["family"],
                    "matrix_sha256_uint8": matrix_hash,
                    "matrix_rank_mod7": matrix_record["rank_mod7"],
                    "matrix_rank_certificate": {
                        "pivot_columns_zero_based": matrix_record[
                            "pivot_columns_zero_based"
                        ],
                        "pivot_minor_determinant_mod7": matrix_record[
                            "pivot_minor_determinant_mod7"
                        ],
                        "exact_full_row_rank_certificate": True,
                    },
                    "retained_all_six_F3_coordinates_verbatim": True,
                    **decision,
                }
            )

        rejected = any(row["rigorously_rejected"] for row in projection_results)
        skipped = not rejected and any(row["skipped"] for row in projection_results)
        necessary = not rejected and not skipped
        require(sum((rejected, skipped, necessary)) == 1, "case decision is ambiguous")
        row = {
            "target_index_in_audited_51_case_order": target_index,
            "case_key": key,
            "catalog_pattern": context["current_by_key"][key]["catalog_pattern"],
            "prior_global_join_decision": context["current_by_key"][key][
                "decision_status"
            ],
            "linear_stabilizer_quotient_results": projection_results,
            "decision_status": (
                "rigorous_linear_stabilizer_quotient_rejection"
                if rejected
                else "explicit_cap_skip_without_negative_decision"
                if skipped
                else "necessary_only_survivor_of_all_completed_linear_quotients"
            ),
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
        }
        row["decision_certificate_sha256"] = quotient.json_sha256(row)
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
    usage = resource.getrusage(resource.RUSAGE_SELF)
    script_path = Path(__file__).resolve()
    return {
        "experiment": EXPERIMENT,
        "status": (
            "complete_sharded_exact_CPU_linear_stabilizer_quotient_join"
            if counts["skipped"] == 0
            else "exact_CPU_linear_stabilizer_quotient_join_with_explicit_skips"
        ),
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": quotient.file_sha256(script_path),
            "quotient_engine_path": str(Path(quotient.__file__).resolve()),
            "quotient_engine_sha256": quotient.file_sha256(Path(quotient.__file__)),
            "linear_matrix_catalog_path": str(Path(linear.__file__).resolve()),
            "linear_matrix_catalog_sha256": quotient.file_sha256(Path(linear.__file__)),
        },
        "configuration": {
            "k": int(matrix_catalog[0]["shape"][0]),
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "selected_matrix_indices": [row[0] for row in selected_matrices],
            "selected_matrix_hashes": [
                row[1]["matrix_sha256_uint8"] for row in selected_matrices
            ],
        },
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_audits": manufactured_audits,
        "real_selector_matrix_equivalence_audit": selector_audit,
        "full_linear_matrix_catalog_audit": matrix_catalog_audit,
        "full_linear_matrix_catalog": matrix_catalog,
        "matrix_sharding_audit": matrix_shard_audit,
        "case_sharding_audit": case_shard_audit,
        "linear_matrix_support_and_stabilizer_audits": matrix_support_audits,
        "projection_decision_census": dict(sorted(decision_census.items())),
        "result_counts": counts,
        "case_results_sha256": quotient.naive_join.old_join.canonical_case_digest(
            case_results
        ),
        "case_results": case_results,
        "logical_semantics": {
            "all_six_F3_coordinates_retained_verbatim": True,
            "every_F7_map_explicit_hash_pinned_and_full_row_rank": True,
            "same_linear_map_used_for_generators_complete_catalogs_and_targets": True,
            "same_source_row_supplies_mod3_and_mod7_before_deduplication": True,
            "every_H3_translation_stabilizer_computed_exactly": True,
            "only_summed_exact_H3_translation_stabilizers_quotiented_out": True,
            "all_eight_factor_images_retained": True,
            "missing_target_after_completed_quotient_join_is_rigorous_rejection": True,
            "target_presence_is_necessary_only": True,
            "state_cap_hit_is_explicit_skip": True,
            "partial_support_after_cap_is_discarded": True,
            "different_matrix_witnesses_are_not_assumed_compatible": True,
            "positive_z7_closure_claimed": False,
        },
        "positive_z7_excluded": False,
        "full_theorem_claimed": False,
        "resource_profile": {
            "elapsed_seconds": time.time() - started,
            "user_cpu_seconds": usage.ru_utime - usage_started.ru_utime,
            "system_cpu_seconds": usage.ru_stime - usage_started.ru_stime,
            "maximum_resident_set_kib": usage.ru_maxrss,
            "GPU_used": False,
        },
        "output_path": str(output_path.resolve()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-input", type=Path, default=linear.cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=linear.cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=linear.cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument("--k", type=int, choices=SUPPORTED_K, default=5)
    parser.add_argument(
        "--matrix-families",
        default=",".join(linear.DEFAULT_MATRIX_FAMILIES),
        help=(
            "comma-separated subset of seeded_dense,geometry,evaluation,"
            "vandermonde,block_sums,selector"
        ),
    )
    parser.add_argument("--matrix-shard-index", type=int, default=0)
    parser.add_argument("--matrix-shard-count", type=int, default=1)
    parser.add_argument("--case-shard-index", type=int, default=0)
    parser.add_argument("--case-shard-count", type=int, default=1)
    cases = parser.add_mutually_exclusive_group()
    cases.add_argument("--all-cases", action="store_true")
    cases.add_argument("--case-count", type=int, default=1)
    parser.add_argument(
        "--state-cap",
        type=int,
        default=DEFAULT_STATE_CAP,
        help="exact support cap; 0 means the complete projected-group size",
    )
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    require(args.state_cap >= 0, "--state-cap cannot be negative")
    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    families = linear.parse_matrix_families(args.matrix_families)
    context, context_audit = linear.coordinate_gpu.build_cpu_context(
        args.parent_input, args.current_join, args.hilbert_basis
    )
    matrix_catalog, matrix_catalog_audit = linear.build_matrix_catalog(
        context, args.k, families
    )
    manufactured_audits = {
        "linear_projection": linear.manufactured_linear_projection_audit(),
        "stabilizer_quotient_present_absent_same_row_stabilizer_and_cap": (
            quotient.manufactured_self_audit()
        ),
    }
    selector_audit = real_selector_equivalence_audit(
        context, args.k, args.pair_chunk_cap
    )

    if args.self_audit_only:
        result = {
            "experiment": f"{EXPERIMENT}_self_audit",
            "status": "bounded_CPU_self_audits_passed",
            "k": args.k,
            "source_provenance": {
                "this_script_path": str(Path(__file__).resolve()),
                "this_script_sha256": quotient.file_sha256(Path(__file__).resolve()),
            },
            "manufactured_audits": manufactured_audits,
            "real_selector_matrix_equivalence_audit": selector_audit,
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
        selected, case_shard_audit = selected_cases(
            context["targets"],
            all_cases=args.all_cases,
            case_count=args.case_count,
            shard_index=args.case_shard_index,
            shard_count=args.case_shard_count,
        )
        group_size = 3**6 * 7**args.k
        state_cap = group_size if args.state_cap == 0 else min(args.state_cap, group_size)
        result = run(
            context=context,
            context_audit=context_audit,
            matrix_catalog=matrix_catalog,
            matrix_catalog_audit=matrix_catalog_audit,
            selected_matrices=selected_matrices,
            matrix_shard_audit=matrix_shard_audit,
            cases=selected,
            case_shard_audit=case_shard_audit,
            state_cap=state_cap,
            pair_chunk_cap=args.pair_chunk_cap,
            selector_audit=selector_audit,
            manufactured_audits=manufactured_audits,
            output_path=args.output,
        )

    linear.cuda_join.atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "k": args.k,
                "processed_cases": result.get("result_counts", {}).get("selected", 0),
                "rejected": result.get("result_counts", {}).get("rejected", 0),
                "surviving": result.get("result_counts", {}).get("surviving", 0),
                "skipped": result.get("result_counts", {}).get("skipped", 0),
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
