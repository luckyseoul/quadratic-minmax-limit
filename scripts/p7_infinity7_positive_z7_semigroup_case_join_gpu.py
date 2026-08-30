#!/usr/bin/env python3
"""Exact CUDA companion for the 51 positive p=7, z=7 grade-three cases.

The CPU semigroup case-join module remains the source of every mathematical
object: the 51-case census, pointed-system reconstruction, torsion quotient,
Normaliz Hilbert basis, joint same-row mod-3/mod-7 directional supports, and
exact projected targets.  This companion only replaces the final eight-factor
global support convolution with the repository's collision-resolved CUDA
open-addressing engine.

Production projections retain all six characteristic-three coordinates and
an explicit width-five or width-six subset of the 21 characteristic-seven
coordinates.  Cases can be deterministically sharded by their audited order.
A missing target after a completed convolution is a rigorous rejection;
target presence is necessary only.  State, memory, allocation, or hash-table
caps are explicit skips and never become negative decisions.

The self-audit combines the CPU and CUDA manufactured same-row/cap traps and
then compares the complete exact CPU and CUDA support sets for one real
grade-three case in F3^6 x F7^2.
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
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join as cpu_join  # noqa: E402
import p7_infinity7_positive_z7_torsion_support_gpu as cuda_join  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_join_gpu"
SUPPORTED_K = (5, 6)
EXPECTED_TARGET_CASES = 51
EXACT_DIRECTION_SUPPORT_CAP = cpu_join.EXPECTED_GRADE_CATALOG_ROWS[3]
DEFAULT_PAIR_CHUNK_CAP = 1_000_000
DEFAULT_GPU_MEMORY_CAP_MIB = 12_000
DEFAULT_CROSSCHECK_PROJECTION = (0, 1)


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
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_cpu_context(
    parent_path: Path,
    current_join_path: Path,
    hilbert_basis_path: Path,
) -> tuple[dict[str, Any], dict]:
    """Run the CPU module's complete construction and provenance audits."""
    parent, representatives, current_by_key, targets, input_audit = (
        cpu_join.load_current_problem(parent_path, current_join_path)
    )
    require(len(targets) == EXPECTED_TARGET_CASES, "grade-three case census changed")
    rebuilt, reconstruction_audit = cpu_join.reconstruct_and_validate(
        parent, representatives
    )
    target_grade_audit = cpu_join.audit_target_grades(
        targets, current_by_key, rebuilt
    )

    system = rebuilt["systems"][0]["A"]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    kernel_rows = np.ascontiguousarray(rebuilt["kernel_rows"], dtype=np.int64)
    basis, degrees, basis_audit = cpu_join.semigroup.load_hilbert_basis(
        hilbert_basis_path, kernel_rows
    )
    common, common_audit = cpu_join.torsion.exact_common_dependency_basis(
        matrix, kernel_rows
    )
    quotient_data, quotient_audit = cpu_join.semigroup.derive_torsion_quotients(
        matrix, common
    )
    require(
        len(quotient_data[3]["complement"]) == 6
        and len(quotient_data[7]["complement"]) == 21,
        "effective torsion dimensions changed",
    )
    orbit = rebuilt["orbits"][0]
    full_generator_rows, generator_audit = cpu_join.full_generator_projections(
        basis, orbit, quotient_data
    )
    context: dict[str, Any] = {
        "representatives": representatives,
        "current_by_key": current_by_key,
        "targets": targets,
        "rebuilt": rebuilt,
        "system": system,
        "basis": basis,
        "degrees": degrees,
        "common": common,
        "quotient_data": quotient_data,
        "orbit": orbit,
        "full_generator_rows": full_generator_rows,
    }
    audit = {
        "input_and_current_decision_audit": input_audit,
        "representative_reconstruction_audit": reconstruction_audit,
        "target_grade_audit": target_grade_audit,
        "normaliz_Hilbert_basis": basis_audit,
        "exact_common_rational_dependency_audit": common_audit,
        "varying_torsion_quotient_audit": quotient_audit,
        "full_generator_projection_audit": generator_audit,
        "all_51_grade_three_cases_reconstructed_before_sharding": True,
    }
    return context, audit


def parse_projection_subsets(value: str | None, k: int) -> tuple[tuple[int, ...], ...]:
    if value is None:
        projections = (tuple(range(k)),)
    else:
        projections = cpu_join.parse_projection_spec(value)
    require(
        1 <= len(projections) <= cuda_join.MAX_PROJECTIONS,
        "projection count lies outside the CUDA engine bound",
    )
    require(
        all(len(row) == k for row in projections),
        f"every explicit mod-7 projection must contain exactly k={k} coordinates",
    )
    require(len(projections) == len(set(projections)), "projection subset repeated")
    return tuple(sorted(projections))


def shard_cases(
    targets: list[dict], shard_index: int, shard_count: int
) -> tuple[list[tuple[int, dict]], dict]:
    require(shard_count > 0, "--shard-count must be positive")
    require(0 <= shard_index < shard_count, "--shard-index lies outside shard count")
    indexed = list(enumerate(targets))
    selected = [row for row in indexed if row[0] % shard_count == shard_index]
    require(selected, "selected shard contains no grade-three cases")
    selected_keys = [str(row[1]["case_key"]) for row in selected]
    all_shards = [
        [str(targets[index]["case_key"]) for index in range(shard, len(targets), shard_count)]
        for shard in range(shard_count)
    ]
    flattened = [key for shard in all_shards for key in shard]
    require(
        len(flattened) == len(targets) and len(set(flattened)) == len(targets),
        "case sharding is not a disjoint cover",
    )
    return selected, {
        "rule": "audited_target_index_mod_shard_count_equals_shard_index",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "full_target_case_count": len(targets),
        "selected_case_count": len(selected),
        "selected_target_indices": [row[0] for row in selected],
        "selected_case_keys": selected_keys,
        "selected_case_keys_sha256": json_sha256(selected_keys),
        "all_shards_form_a_disjoint_cover_of_the_51_cases": True,
    }


def build_direction_supports(
    context: dict[str, Any],
    projection: tuple[int, ...],
    pair_chunk_cap: int,
) -> tuple[dict[int, dict[int, np.ndarray]], dict]:
    """Delegate exact same-row support construction and catalog audits to CPU."""
    table, audit = cpu_join.projected_direction_supports(
        projection=projection,
        full_generator_rows=context["full_generator_rows"],
        basis=context["basis"],
        degrees=context["degrees"],
        orbit=context["orbit"],
        anchors=context["rebuilt"]["anchors"],
        quotient_data=context["quotient_data"],
        state_cap=EXACT_DIRECTION_SUPPORT_CAP,
        pair_chunk_cap=max(1, min(pair_chunk_cap, EXACT_DIRECTION_SUPPORT_CAP)),
    )
    require(
        audit["all_direction_grades_zero_through_three_completed"],
        "an exact grade-zero-through-three direction support did not complete",
    )
    require(
        audit["all_completed_grade_supports_match_complete_direct_catalogs"],
        "a semigroup direction support failed its direct-catalog calibration",
    )
    require(
        all(set(table[direction]) == {0, 1, 2, 3} for direction in range(8)),
        "direction support table is incomplete",
    )
    return table, audit


def exact_case_inputs(
    context: dict[str, Any],
    target_row: dict,
    projection: tuple[int, ...],
    support_table: dict[int, dict[int, np.ndarray]],
) -> tuple[list[np.ndarray], np.ndarray, tuple[int, ...], tuple[tuple[int, np.ndarray], ...], dict]:
    """Obtain one exact target and its eight exact directional support factors."""
    rebuilt = context["rebuilt"]
    orbit, leaf, system, _factory = cpu_join.old_join.validate_parent_survivor(
        target_row, rebuilt
    )
    grades = tuple(
        cpu_join.leaf_grade(orbit, leaf, direction) for direction in range(8)
    )
    require(max(grades) == 3, "selected case is not a grade-three case")
    require(sum(grades) == 14, "selected case total directional grade changed")

    moduli = (3,) * 6 + (7,) * len(projection)
    codec = cpu_join.semigroup.MixedRadixCodec(moduli)
    factors = tuple(
        (direction, support_table[direction][grade])
        for direction, grade in enumerate(grades)
    )
    require(all(len(codes) for _direction, codes in factors), "empty direction support")
    # A full factor is placed first so the imported CUDA engine's exact
    # saturation shortcut fires immediately; otherwise small-first minimizes
    # intermediate supports.  The group operation is commutative.
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
    base_digits = cpu_join.project_equation_vector(
        anchor_rhs, context["quotient_data"], projection
    )
    target_digits = np.ascontiguousarray(
        (-base_digits.astype(np.int16)) % np.asarray(moduli, dtype=np.int16),
        dtype=np.uint8,
    )

    # Repeat the CPU module's exact-dependency audit for this precise target.
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
            "CPU mixed-radix support changed at the CUDA interface",
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
        "directional_excess_grades": list(grades),
        "convolution_direction_order": [row[0] for row in ordered],
        "convolution_order_rule": "full_group_first_else_support_size_then_direction",
        "direction_support_records": interface_records,
        "base_digits_sha256_uint8": cpu_join.array_sha256(base_digits),
        "target_digits": target_digits.tolist(),
        "target_digits_sha256_uint8": cpu_join.array_sha256(target_digits),
        "common_exact_dependency_syndrome_checked_with_direct_catalog_rows": True,
        "CPU_codes_decode_and_reencode_identically_in_CUDA_engine_codec": True,
        "same_direction_rows_couple_mod3_mod7_before_global_convolution": True,
    }
    return support_rows, target_digits, moduli, ordered, audit


def validate_cuda_decision(decision: dict) -> None:
    status = str(decision["decision_status"])
    if status.startswith("skipped_"):
        require(not decision["completed_exact_convolution"], "skip labeled complete")
        require(not decision["rigorous_rejection"], "cap skip became a rejection")
        require(decision["projected_target_present"] is None, "skip guessed target status")
        require(decision["skip_is_explicit_not_approximation"], "skip semantics weakened")
    else:
        require(decision["completed_exact_convolution"], "non-skip did not complete")
        require(
            decision["rigorous_rejection"] is (not decision["projected_target_present"]),
            "completed CUDA decision is inconsistent",
        )
        require(decision["cuda_hash_set_is_exact_not_probabilistic"], "CUDA set weakened")


def manufactured_cross_engine_audit(
    cp: Any, kernel: Any, gpu_memory_cap_bytes: int
) -> dict:
    cpu_audit = cpu_join.manufactured_case_join_audit()
    gpu_audit = cuda_join.manufactured_cpu_gpu_self_audit(
        cp, kernel, gpu_memory_cap_bytes
    )
    require(cpu_audit["passed"], "CPU manufactured case-join audit failed")
    require(
        cpu_audit["same_row_cross_prime_false_positive_trap_rejected"]
        and gpu_audit["same_catalog_row_cross_prime_trap_rejected"],
        "a CPU/GPU same-row trap failed",
    )
    require(
        cpu_audit["state_cap_is_explicit_skip_and_partial_support_is_discarded"]
        and gpu_audit["state_cap_produces_explicit_skip"],
        "a CPU/GPU cap trap failed",
    )
    return {
        "status": "passed",
        "CPU_manufactured_case_join": cpu_audit,
        "CUDA_manufactured_exact_set": gpu_audit,
        "same_row_cross_prime_traps_passed_on_CPU_and_GPU": True,
        "state_cap_traps_are_explicit_nondecision_skips_on_CPU_and_GPU": True,
        "CUDA_pair_chunk_invariance_and_collision_resolution_audited": True,
    }


def real_case_cpu_gpu_crosscheck(
    cp: Any,
    kernel: Any,
    gpu_memory_cap_bytes: int,
    context: dict[str, Any],
    pair_chunk_cap: int,
) -> dict:
    """Compare complete support sets for one real case in a small quotient."""
    projection = DEFAULT_CROSSCHECK_PROJECTION
    table, direction_audit = build_direction_supports(
        context, projection, pair_chunk_cap
    )
    target_row = context["targets"][0]
    support_rows, target, moduli, ordered, case_audit = exact_case_inputs(
        context, target_row, projection, table
    )
    codec = cpu_join.semigroup.MixedRadixCodec(moduli)
    cpu_codes, cpu_convolution = cpu_join.convolve_support_sequence(
        ordered,
        codec,
        state_cap=codec.group_size,
        pair_chunk_cap=min(pair_chunk_cap, codec.group_size),
    )
    require(cpu_codes is not None, "real-case CPU crosscheck unexpectedly capped")
    require(cpu_convolution["completed"], "real-case CPU crosscheck is incomplete")

    gpu_decision, gpu_codes = cuda_join.gpu_exact_support_convolution(
        cp,
        kernel,
        support_rows,
        target,
        moduli,
        state_cap=codec.group_size,
        pair_chunk_cap=pair_chunk_cap,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    validate_cuda_decision(gpu_decision)
    require(gpu_codes is not None, "real-case CUDA crosscheck unexpectedly skipped")
    require(np.array_equal(cpu_codes, gpu_codes), "real-case CPU/CUDA support sets differ")
    target_code = int(codec.encode(target[None, :])[0])
    target_index = int(np.searchsorted(cpu_codes, np.uint64(target_code)))
    cpu_present = bool(
        target_index < len(cpu_codes) and int(cpu_codes[target_index]) == target_code
    )
    require(
        cpu_present == gpu_decision["projected_target_present"],
        "real-case CPU/CUDA target decisions differ",
    )
    return {
        "status": "passed",
        "projection_mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "case_input_audit": case_audit,
        "direction_support_audit": direction_audit,
        "CPU_convolution_audit": cpu_convolution,
        "CUDA_convolution_audit": gpu_decision,
        "complete_support_state_count": len(cpu_codes),
        "complete_support_sha256_uint64": cpu_join.array_sha256(
            cpu_codes.astype("<u8", copy=False)
        ),
        "target_present": cpu_present,
        "CPU_and_CUDA_complete_support_sets_equal_elementwise": True,
        "CPU_and_CUDA_target_decisions_equal": True,
    }


def expected_engine_memory(
    k: int, state_cap: int, pair_chunk_cap: int
) -> dict:
    """Conservative bound using the imported CUDA engine's sizing formula."""
    moduli = (3,) * 6 + (7,) * k
    group_size = cuda_join.group_size_for(moduli)
    effective_state_cap = min(state_cap, group_size)
    maximum_insertions = min(
        group_size, effective_state_cap + pair_chunk_cap + 1
    )
    table_size = cuda_join.next_power_of_two(max(8, 2 * maximum_insertions))
    estimated_bytes = (
        table_size * 24
        + (effective_state_cap + EXACT_DIRECTION_SUPPORT_CAP) * 8
        + maximum_insertions * 16
        + len(moduli) * 4
    )
    return {
        "projected_group_size": group_size,
        "effective_state_cap": effective_state_cap,
        "maximum_direction_support_rows": EXACT_DIRECTION_SUPPORT_CAP,
        "hash_table_entries_upper_bound": table_size,
        "conservative_engine_peak_bytes_upper_bound": estimated_bytes,
        "conservative_engine_peak_mib_upper_bound": estimated_bytes / (1 << 20),
        "formula_matches_imported_CUDA_engine_peak_estimator": True,
        "host_final_code_array_bytes_if_full_group": group_size * 8,
    }


def run_production(
    cp: Any,
    kernel: Any,
    gpu_audit: dict,
    context: dict[str, Any],
    context_audit: dict,
    manufactured_audit: dict,
    real_crosscheck: dict | None,
    projections: tuple[tuple[int, ...], ...],
    selected_cases: list[tuple[int, dict]],
    shard_audit: dict,
    state_cap: int,
    pair_chunk_cap: int,
    gpu_memory_cap_bytes: int,
    output_path: Path,
) -> dict:
    started = time.time()
    projection_tables: dict[tuple[int, ...], dict[int, dict[int, np.ndarray]]] = {}
    support_audits = []
    for projection in projections:
        table, audit = build_direction_supports(context, projection, pair_chunk_cap)
        projection_tables[projection] = table
        support_audits.append(audit)

    projection_totals: Counter[str] = Counter()
    case_results = []
    for target_index, target_row in selected_cases:
        projection_results = []
        for projection in projections:
            support_rows, target, moduli, _ordered, case_input_audit = exact_case_inputs(
                context, target_row, projection, projection_tables[projection]
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
            validate_cuda_decision(decision)
            projection_totals[str(decision["decision_status"])] += 1
            projection_results.append(
                {
                    "mod7_coordinates": list(projection),
                    "retained_all_six_mod3_coordinates": True,
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
            "projection_results": projection_results,
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
            "decision_status": (
                "rigorous_exact_CUDA_semigroup_projection_rejection"
                if rejected
                else "explicit_cap_skip_without_negative_decision"
                if skipped
                else "necessary_only_survivor_of_all_completed_CUDA_projections"
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
    k = len(projections[0])
    script_path = Path(__file__).resolve()
    return {
        "experiment": EXPERIMENT,
        "status": "complete_sharded_exact_CUDA_semigroup_case_join",
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
            "CPU_semigroup_case_join_path": str(Path(cpu_join.__file__).resolve()),
            "CPU_semigroup_case_join_sha256": file_sha256(Path(cpu_join.__file__)),
            "CUDA_torsion_engine_path": str(Path(cuda_join.__file__).resolve()),
            "CUDA_torsion_engine_sha256": file_sha256(Path(cuda_join.__file__)),
        },
        "configuration": {
            "k": k,
            "explicit_mod7_projection_subsets": [list(row) for row in projections],
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "gpu_memory_cap_bytes": gpu_memory_cap_bytes,
            "exact_direction_support_cap": EXACT_DIRECTION_SUPPORT_CAP,
        },
        "expected_memory": expected_engine_memory(k, state_cap, pair_chunk_cap),
        "gpu_engine_audit": gpu_audit,
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_CPU_CUDA_audit": manufactured_audit,
        "small_real_case_CPU_CUDA_crosscheck": real_crosscheck,
        "small_real_case_crosscheck_skipped_by_explicit_flag": real_crosscheck is None,
        "case_sharding_audit": shard_audit,
        "projected_direction_support_audits": support_audits,
        "result_counts": counts,
        "projection_decision_census": dict(sorted(projection_totals.items())),
        "case_results_sha256": cpu_join.old_join.canonical_case_digest(case_results),
        "case_results": case_results,
        "logical_semantics": {
            "all_51_cases_were_reconstructed_before_case_sharding": True,
            "all_six_mod3_coordinates_are_retained": True,
            "mod7_coordinates_are_explicitly_selected_from_the_derived_21D_quotient": True,
            "same_Hilbert_or_catalog_row_supplies_mod3_and_mod7_before_deduplication": True,
            "grade_zero_through_three_direction_supports_are_exact_and_direct_catalog_calibrated": True,
            "all_eight_direction_supports_are_convolved_in_every_completed_projection": True,
            "CUDA_hash_collisions_are_resolved_by_exact_key_comparison": True,
            "exact_full_group_saturation_shortcut_only_proves_target_presence": True,
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--parent-input", type=Path, default=cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument("--k", type=int, choices=SUPPORTED_K, default=5)
    parser.add_argument(
        "--mod7-projections",
        help=(
            "explicit semicolon-separated width-k subsets, e.g. "
            "'0,1,2,3,4;5,6,7,8,9'; default is the width-k prefix"
        ),
    )
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    require(args.gpu_memory_cap_mib > 0, "--gpu-memory-cap-mib must be positive")
    projections = parse_projection_subsets(args.mod7_projections, args.k)
    group_size = cuda_join.group_size_for((3,) * 6 + (7,) * args.k)
    require(args.state_cap >= 0, "--state-cap cannot be negative")
    state_cap = group_size if args.state_cap == 0 else min(args.state_cap, group_size)
    gpu_memory_cap_bytes = args.gpu_memory_cap_mib * (1 << 20)

    cp, kernel, gpu_audit = cuda_join.load_cupy(args.device)
    with cp.cuda.Device(args.device):
        manufactured = manufactured_cross_engine_audit(
            cp, kernel, gpu_memory_cap_bytes
        )
        context, context_audit = build_cpu_context(
            args.parent_input, args.current_join, args.hilbert_basis
        )
        real_crosscheck = (
            None
            if args.skip_small_real_crosscheck and not args.self_audit_only
            else real_case_cpu_gpu_crosscheck(
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
                "status": "manufactured_and_real_case_CPU_CUDA_self_audits_passed",
                "gpu_engine_audit": gpu_audit,
                "CPU_construction_and_input_audits": context_audit,
                "manufactured_CPU_CUDA_audit": manufactured,
                "small_real_case_CPU_CUDA_crosscheck": real_crosscheck,
                "production_cases_processed": 0,
                "output_path": str(args.output.resolve()),
            }
        else:
            selected_cases, shard_audit = shard_cases(
                context["targets"], args.shard_index, args.shard_count
            )
            result = run_production(
                cp,
                kernel,
                gpu_audit,
                context,
                context_audit,
                manufactured,
                real_crosscheck,
                projections,
                selected_cases,
                shard_audit,
                state_cap,
                args.pair_chunk_cap,
                gpu_memory_cap_bytes,
                args.output,
            )

    cuda_join.atomic_write(args.output, result)
    summary = {
        "status": result["status"],
        "output": str(args.output),
        "gpu": gpu_audit["device_name"],
        "manufactured_CPU_CUDA_audit": manufactured["status"],
        "small_real_case_CPU_CUDA_crosscheck": (
            None if real_crosscheck is None else real_crosscheck["status"]
        ),
        "processed_cases": result.get("result_counts", {}).get("selected", 0),
        "rigorously_rejected_cases": result.get("result_counts", {}).get("rejected", 0),
        "necessary_only_survivors": result.get("result_counts", {}).get("surviving", 0),
        "skipped_cases": result.get("result_counts", {}).get("skipped", 0),
        "expected_memory": expected_engine_memory(
            args.k, state_cap, args.pair_chunk_cap
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
