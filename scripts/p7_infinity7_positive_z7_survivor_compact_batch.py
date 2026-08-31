#!/usr/bin/env python3
"""Resume-safe compact CP-SAT batch for the 1,296 z=7 hull survivors.

The input is the complete pointed mod-3/mod-7 affine-hull evidence.  Its
full 4,320-case certificate and ordered 1,296 survivor-key certificate are
recomputed and anchored before exact orbit leaves are reconstructed.  Each
selected survivor gets the compact pointed relaxation with all eight means
fixed to that leaf.  INFEASIBLE rigorously rejects the leaf; FEASIBLE is a
directly audited necessary-only relaxation witness; UNKNOWN stays unknown.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine_hull  # noqa: E402
import p7_infinity7_positive_z7_pointed_compact_cpsat as compact  # noqa: E402


EXPECTED_EXPERIMENT = "p7_infinity7_positive_z7_pointed_affine_hull_multimod"
EXPECTED_STATUS = "complete_rigorous_pointed_mod3_mod7_affine_hull_necessary_sieve_with_survivors"
EXPECTED_SOURCE_LEAVES = 2_160
EXPECTED_POINTED_CASES = 4_320
EXPECTED_REJECTED = 3_024
EXPECTED_SURVIVORS = 1_296
EXPECTED_CASE_CERTIFICATE = "3f65e57cf2f09bc4c674711e3bda3a46503c49169ca54aee4728de5619976aaf"
EXPECTED_SURVIVOR_CERTIFICATE = "f756e0128e12c78ad7a17f85dd621e4e9ff00f0be80c06ac60a71f74045fc784"
EXPECTED_PROJECTION_CENSUS = {"positive_hit": 11_304, "tested": 17_824, "zero_hit": 6_520}
EXPERIMENT = "p7_infinity7_positive_z7_survivor_compact_batch"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ordered_rows_sha256(rows: list[dict]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def parse_positive(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def load_and_audit_evidence(
    path: Path,
    orbits: list[dict],
    leaves_by_orbit: list[list[dict]],
    leaf_audit: dict,
    systems_by_orbit: list[dict[str, dict]],
    compact_input_audit: dict,
) -> tuple[list[dict], dict]:
    """Load the large JSON once and independently certify its full contents."""
    raw = path.read_bytes()
    file_sha256 = hashlib.sha256(raw).hexdigest()
    evidence = json.loads(raw)
    require(evidence["experiment"] == EXPECTED_EXPERIMENT, "wrong affine-hull experiment")
    require(evidence["status"] == EXPECTED_STATUS, "affine-hull evidence is not the canonical full survivor run")
    require(
        evidence["p"] == 7
        and evidence["c_H"] == 1
        and evidence["z"] == 7
        and evidence["phase"] == 0
        and evidence["moduli"] == [3, 7],
        "affine-hull evidence scope changed",
    )
    require(evidence["full_run"] is True and evidence["smoke_test"] is False, "input is not a full run")
    require(evidence["full_source_exact_mean_leaves"] == EXPECTED_SOURCE_LEAVES, "full leaf target changed")
    require(evidence["processed_source_exact_mean_leaves"] == EXPECTED_SOURCE_LEAVES, "source leaf coverage incomplete")
    require(evidence["full_pointed_branch_cases"] == EXPECTED_POINTED_CASES, "full pointed target changed")
    require(evidence["processed_pointed_branch_cases"] == EXPECTED_POINTED_CASES, "pointed coverage incomplete")
    require(evidence["rejected_pointed_branch_cases"] == EXPECTED_REJECTED, "rejection census changed")
    require(evidence["surviving_pointed_branch_cases"] == EXPECTED_SURVIVORS, "survivor census changed")
    require(EXPECTED_REJECTED + EXPECTED_SURVIVORS == EXPECTED_POINTED_CASES, "canonical census inconsistent")
    require(evidence["projection_test_census"] == EXPECTED_PROJECTION_CENSUS, "projection census changed")
    require(evidence["z7_branch_excluded"] is False, "survivor evidence incorrectly claims closure")
    require(evidence["orbit_source"] == compact_input_audit["orbit_source"], "orbit source changed")
    require(evidence["mean_leaf_coverage"] == leaf_audit, "mean-leaf coverage certificate changed")
    require(
        evidence["translation_equivariant_linear_system"]
        == compact_input_audit["translation_system_audit"],
        "translation-system audit changed",
    )
    require(
        evidence["pointed_linear_systems"] == compact_input_audit["pointed_system_audits"],
        "pointed-system audits changed",
    )

    semantics = evidence["logical_semantics"]
    for key in (
        "all_2160_exact_mean_leaves_reconstructed_before_selection",
        "both_pointed_branches_are_disjoint_and_exhaustive_for_nonempty_stars",
        "all_planned_retained_pairs_are_tested_even_after_an_earlier_zero_hit",
        "same_retained_catalog_indices_are_used_modulo_3_and_7",
        "zero_hit_affine_hull_projection_is_rigorous_exclusion",
    ):
        require(semantics[key] is True, f"required affine-hull semantic audit failed: {key}")
    for key in (
        "affine_hull_passing_is_exact_catalog_tuple_feasibility",
        "all_pair_tests_share_one_global_catalog_assignment",
        "modular_passing_is_binary_edge_feasibility",
        "exact_edge_lift_claimed",
        "smoke_run_can_claim_z7_exclusion",
    ):
        require(semantics[key] is False, f"affine-hull survivor semantics changed: {key}")

    case_results = evidence["case_results"]
    require(len(case_results) == EXPECTED_POINTED_CASES, "case-result list is incomplete")
    observed_case_certificate = ordered_rows_sha256(case_results)
    require(evidence["all_case_results_sha256"] == EXPECTED_CASE_CERTIFICATE, "stored case certificate changed")
    require(observed_case_certificate == EXPECTED_CASE_CERTIFICATE, "case-result certificate recomputation failed")

    expected_keys = []
    reconstructed = {}
    for orbit_index, (orbit, leaves) in enumerate(zip(orbits, leaves_by_orbit)):
        require(len(leaves) == 1_080, "per-orbit reconstructed leaf count changed")
        for leaf_index, leaf in enumerate(leaves):
            for branch in ("A", "B"):
                key = f"orbit{orbit_index}_leaf{leaf_index}_branch{branch}"
                expected_keys.append(key)
                reconstructed[key] = (orbit_index, leaf_index, branch, orbit, leaf)
    require(len(expected_keys) == EXPECTED_POINTED_CASES, "reconstructed pointed key count changed")
    require([row["case_key"] for row in case_results] == expected_keys, "case-result order or coverage changed")

    filtered_survivors = []
    for row in case_results:
        key = row["case_key"]
        orbit_index, leaf_index, branch, orbit, leaf = reconstructed[key]
        require(row["branch_orbit_index"] == orbit_index, "case orbit index changed")
        require(row["orbit_leaf_index"] == leaf_index, "case leaf index changed")
        require(row["pointed_star_branch"] == branch, "case branch changed")
        require(row["source_orbit_index"] == orbit["source_orbit_index"], "source orbit changed")
        require(row["representative_finite_field"] == list(orbit["representative"]), "representative changed")
        require(row["fixed_edge_rows"] == systems_by_orbit[orbit_index][branch]["fixed_edge_rows"], "fixed rows changed")
        require(row["scaled_means"] == list(leaf["scaled_means"]), "case means disagree with reconstructed leaf")
        require(row["q_values"] == list(leaf["q_values"]), "case q-values changed")
        require(row["residue_pair_minus_plus"] == leaf["residue_pair"], "case residue pair changed")
        require(row["catalog_levels"] == list(leaf["catalog_levels"]), "case levels changed")
        require(row["catalog_classes"] == list(leaf["catalog_classes"]), "case classes changed")
        require(row["catalog_pattern_H_S_M"] == list(leaf["pattern"]), "case pattern changed")
        require(row["all_planned_retained_pairs_tested"] is True, "case skipped a planned projection")
        require(row["tested_projection_count"] == row["planned_projection_count"], "projection count mismatch")
        passing = bool(row["passes_all_necessary_projections"])
        require(passing == (not row["failed_projection_indices"]), "passing/failure indices disagree")
        if passing:
            filtered_survivors.append(key)

    survivor_keys = evidence["survivor_case_keys"]
    require(len(survivor_keys) == len(set(survivor_keys)) == EXPECTED_SURVIVORS, "survivor keys are not 1,296 unique cases")
    require(survivor_keys == filtered_survivors, "survivor key list disagrees with case decisions")
    observed_survivor_certificate = affine_hull.json_sha256(survivor_keys)
    require(evidence["survivor_case_keys_sha256"] == EXPECTED_SURVIVOR_CERTIFICATE, "stored survivor certificate changed")
    require(observed_survivor_certificate == EXPECTED_SURVIVOR_CERTIFICATE, "survivor certificate recomputation failed")
    survivor_histogram = Counter(
        (reconstructed[key][0], reconstructed[key][2]) for key in survivor_keys
    )
    require(
        survivor_histogram
        == Counter({(0, "A"): 324, (0, "B"): 324, (1, "A"): 324, (1, "B"): 324}),
        "survivor orbit/branch census changed",
    )

    survivors = []
    for ordinal, key in enumerate(survivor_keys):
        orbit_index, leaf_index, branch, _orbit, leaf = reconstructed[key]
        survivors.append(
            {
                "survivor_ordinal": ordinal,
                "case_key": key,
                "branch_orbit_index": orbit_index,
                "orbit_leaf_index": leaf_index,
                "pointed_star_branch": branch,
                "scaled_means": list(leaf["scaled_means"]),
            }
        )
    return survivors, {
        "path": str(path),
        "file_size_bytes": len(raw),
        "file_sha256": file_sha256,
        "experiment": evidence["experiment"],
        "status": evidence["status"],
        "full_source_exact_mean_leaves": EXPECTED_SOURCE_LEAVES,
        "full_pointed_branch_cases": EXPECTED_POINTED_CASES,
        "rejected_pointed_branch_cases": EXPECTED_REJECTED,
        "surviving_pointed_branch_cases": EXPECTED_SURVIVORS,
        "all_case_results_sha256": observed_case_certificate,
        "survivor_case_keys_sha256": observed_survivor_certificate,
        "projection_test_census": EXPECTED_PROJECTION_CENSUS,
        "survivor_orbit_branch_histogram": {
            f"orbit{orbit}_{branch}": count
            for (orbit, branch), count in sorted(survivor_histogram.items())
        },
        "all_4320_ordered_case_rows_cross_audited_against_reconstructed_leaves": True,
        "all_1296_survivor_keys_unique_reconstructed_and_certified": True,
    }


def process_survivor(
    survivor: dict,
    case_template: dict,
    system: dict,
    kernel: np.ndarray,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool,
) -> dict:
    """Build one compact model, fix all means, and optionally solve it."""
    from ortools.sat.python import cp_model

    started = time.time()
    model, variables, construction = compact.build_model(case_template, system, kernel)
    fixed_means = tuple(int(value) for value in survivor["scaled_means"])
    require(len(fixed_means) == 8, "survivor mean vector changed length")
    for direction, value in enumerate(fixed_means):
        model.add(variables["means"][direction] == value)
    validation_error = model.validate()
    require(not validation_error, f"fixed-mean compact model is invalid: {validation_error}")
    proto = model.Proto()
    require(
        len(proto.variables) == compact.EXPECTED_VARIABLES,
        "fixed-mean model variable count changed",
    )
    require(
        len(proto.constraints) == compact.EXPECTED_CONSTRAINTS + 8,
        "fixed-mean equality count changed",
    )
    base = {
        "survivor_ordinal": int(survivor["survivor_ordinal"]),
        "case_key": survivor["case_key"],
        "branch_orbit_index": int(survivor["branch_orbit_index"]),
        "orbit_leaf_index": int(survivor["orbit_leaf_index"]),
        "pointed_star_branch": survivor["pointed_star_branch"],
        "fixed_scaled_means": list(fixed_means),
        "seed": seed,
        "timeout_seconds": timeout,
        "workers": workers,
        "model_audit": {
            "variables": len(proto.variables),
            "constraints_before_fixed_means": int(construction["total_model_constraints"]),
            "fixed_mean_equalities": 8,
            "constraints_after_fixed_means": len(proto.constraints),
            "slack_intvars_0_to_13": int(construction["slack_intvars_0_to_13"]),
            "mod3_dependency_constraints": int(construction["modular_dependency_constraints"]["3"]),
            "mod7_dependency_constraints": int(construction["modular_dependency_constraints"]["7"]),
            "model_validation": "passed",
            "all_eight_means_fixed_exactly": True,
        },
    }
    if construction_only:
        return {
            **base,
            "solver_status": "NOT_RUN_CONSTRUCTION_ONLY",
            "rigorous_status": "fixed_mean_model_constructed_no_solver_decision",
            "feasible": None,
            "necessary_only_relaxation_survivor": None,
            "finite_infeasibility_certificate": False,
            "case_decided": False,
            "conflicts": None,
            "branches": None,
            "solver_wall_time_seconds": 0.0,
            "elapsed_seconds": time.time() - started,
        }

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(timeout)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 2
    status = solver.solve(model)
    require(status != cp_model.MODEL_INVALID, "validated fixed-mean model became invalid")
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    infeasible = status == cp_model.INFEASIBLE
    if infeasible:
        rigorous_status = "rigorously_excluded_fixed_mean_affine_hull_survivor"
    elif feasible:
        rigorous_status = "feasible_fixed_mean_compact_relaxation_necessary_only"
    else:
        rigorous_status = "unknown_fixed_mean_compact_relaxation"
    result = {
        **base,
        "solver_status": solver.status_name(status),
        "rigorous_status": rigorous_status,
        "feasible": feasible,
        "necessary_only_relaxation_survivor": feasible,
        "finite_infeasibility_certificate": infeasible,
        "case_decided": feasible or infeasible,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "solver_wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        witness, audit = compact.direct_witness_audit(
            solver,
            case_template,
            system,
            kernel,
            variables,
        )
        fixed_mean_valid = witness["scaled_means_direction_order"] == list(fixed_means)
        require(fixed_mean_valid and audit["valid"], "fixed-mean witness audit failed")
        audit["fixed_scaled_means_match_reconstructed_leaf"] = True
        witness["fixed_scaled_means"] = list(fixed_means)
        witness["necessary_only_not_an_edge_witness"] = True
        result["relaxation_witness"] = witness
    return result


def result_certificate(results: list[dict]) -> str:
    return ordered_rows_sha256(results)


def validate_resume_results(
    results: list[dict],
    target: list[dict],
    construction_only: bool,
) -> None:
    require(len(results) <= len(target), "checkpoint contains more cases than current target")
    require(
        [row["case_key"] for row in results]
        == [row["case_key"] for row in target[: len(results)]],
        "checkpoint results are not the canonical target prefix",
    )
    for result, survivor in zip(results, target):
        require(result["survivor_ordinal"] == survivor["survivor_ordinal"], "checkpoint ordinal changed")
        require(result["fixed_scaled_means"] == survivor["scaled_means"], "checkpoint fixed means changed")
        require(result["pointed_star_branch"] == survivor["pointed_star_branch"], "checkpoint branch changed")
        if construction_only:
            require(result["solver_status"] == "NOT_RUN_CONSTRUCTION_ONLY", "construction checkpoint contains solve")
        elif result["feasible"] is True:
            require(
                result["relaxation_witness"]["direct_witness_audit"]["valid"] is True
                and result["relaxation_witness"]["direct_witness_audit"][
                    "fixed_scaled_means_match_reconstructed_leaf"
                ]
                is True,
                "checkpoint feasible witness lacks direct fixed-mean audit",
            )


def checkpoint_payload(
    *,
    evidence_audit: dict,
    full_survivors: list[dict],
    selected_unlimited: list[dict],
    target: list[dict],
    results: list[dict],
    orbit: str,
    branch: str,
    case_limit: int | None,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool,
    kernel_audit: dict,
    reconstruction_audit: dict,
    elapsed_seconds: float,
) -> dict:
    infeasible = sum(bool(row["finite_infeasibility_certificate"]) for row in results)
    feasible = sum(row["feasible"] is True for row in results)
    not_run = sum(row["solver_status"] == "NOT_RUN_CONSTRUCTION_ONLY" for row in results)
    unknown = len(results) - infeasible - feasible - not_run
    remaining = len(target) - len(results)
    full_scope = (
        orbit == "all"
        and branch == "ALL"
        and case_limit is None
        and len(target) == EXPECTED_SURVIVORS
    )
    closure = (
        not construction_only
        and full_scope
        and len(results) == EXPECTED_SURVIVORS
        and infeasible == EXPECTED_SURVIVORS
    )
    if remaining:
        status = "running_resume_safe_checkpoint"
        conclusion = "selected fixed-mean batch remains incomplete"
    elif construction_only:
        status = "complete_fixed_mean_construction_only"
        conclusion = "all targeted fixed-mean models constructed; no solver decisions requested"
    elif closure:
        status = "complete_rigorous_positive_z7_survivor_compact_exclusion"
        conclusion = "all 1,296 certified affine-hull survivors are compact-infeasible; positive z=7 is closed"
    elif infeasible == len(results):
        status = "complete_rigorous_selected_survivor_compact_exclusion"
        conclusion = "every selected survivor leaf is compact-infeasible; unselected survivors remain"
    elif feasible:
        status = "complete_selected_batch_with_necessary_only_survivors"
        conclusion = "at least one fixed-mean compact relaxation is feasible; no edge witness is claimed"
    else:
        status = "complete_selected_batch_with_unknown_cases"
        conclusion = "at least one fixed-mean solve is unknown; no exclusion follows for that leaf"

    full_keys = [row["case_key"] for row in full_survivors]
    selected_keys = [row["case_key"] for row in selected_unlimited]
    target_keys = [row["case_key"] for row in target]
    return {
        "experiment": EXPERIMENT,
        "status": status,
        "rigorous_conclusion": conclusion,
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "model_semantics": {
            "input_survivor": "passed every necessary affine-hull projection",
            "fixed_model": "compact pointed edge relaxation with all eight exact means fixed",
            "infeasible": "rigorously excludes that affine-hull survivor leaf",
            "feasible": "necessary-only compact relaxation witness, directly audited; not an edge witness",
            "unknown": "honest unresolved case",
        },
        "input_evidence_audit": evidence_audit,
        "kernel_audit": kernel_audit,
        "reconstruction_and_pointed_system_audit": reconstruction_audit,
        "survivor_universe": {
            "count": len(full_keys),
            "case_keys_sha256": json_sha256(full_keys),
            "canonical_case_keys": full_keys,
            "all_keys_reconstructed_from_exact_orbit_leaves": True,
        },
        "selection": {
            "orbit": orbit,
            "branch": branch,
            "case_limit": case_limit,
            "matching_survivors_before_limit": len(selected_keys),
            "matching_survivor_keys_sha256": json_sha256(selected_keys),
            "target_cases_after_limit": len(target_keys),
            "target_case_keys_sha256": json_sha256(target_keys),
            "full_1296_case_scope": full_scope,
        },
        "configuration": {
            "timeout_seconds_per_case": timeout,
            "workers_per_case": workers,
            "seed_base": seed,
            "construction_only": construction_only,
        },
        "resume_audit": {
            "atomic_checkpoint_after_every_case": True,
            "canonical_prefix_processing": True,
            "case_seed_formula": "seed_base + survivor_ordinal",
            "processed_case_results_sha256": result_certificate(results),
        },
        "targeted_survivor_cases": len(target),
        "processed_survivor_cases": len(results),
        "remaining_survivor_cases": remaining,
        "infeasible_survivor_cases": infeasible,
        "necessary_only_feasible_survivor_cases": feasible,
        "unknown_survivor_cases": unknown,
        "construction_only_cases": not_run,
        "all_targeted_cases_processed": remaining == 0,
        "all_targeted_cases_infeasible": (
            not construction_only and remaining == 0 and infeasible == len(target)
        ),
        "closure_requires_all_1296_canonical_survivors_infeasible": True,
        "all_1296_survivors_infeasible": closure,
        "z7_branch_excluded": closure,
        "case_results": results,
        "elapsed_seconds": elapsed_seconds,
    }


def load_resume_checkpoint(
    output: Path,
    evidence_audit: dict,
    selected_unlimited: list[dict],
    target: list[dict],
    orbit: str,
    branch: str,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool,
) -> tuple[list[dict], float]:
    if not output.exists():
        return [], 0.0
    checkpoint = json.loads(output.read_text(encoding="utf-8"))
    require(checkpoint["experiment"] == EXPERIMENT, "output is not a survivor compact checkpoint")
    old_input = checkpoint["input_evidence_audit"]
    require(old_input["file_sha256"] == evidence_audit["file_sha256"], "checkpoint input file changed")
    require(
        old_input["all_case_results_sha256"] == EXPECTED_CASE_CERTIFICATE
        and old_input["survivor_case_keys_sha256"] == EXPECTED_SURVIVOR_CERTIFICATE,
        "checkpoint input certificates changed",
    )
    selection = checkpoint["selection"]
    require(selection["orbit"] == orbit and selection["branch"] == branch, "checkpoint selectors changed")
    require(
        selection["matching_survivor_keys_sha256"]
        == json_sha256([row["case_key"] for row in selected_unlimited]),
        "checkpoint selected survivor universe changed",
    )
    configuration = checkpoint["configuration"]
    require(float(configuration["timeout_seconds_per_case"]) == timeout, "checkpoint timeout changed")
    require(int(configuration["workers_per_case"]) == workers, "checkpoint worker count changed")
    require(int(configuration["seed_base"]) == seed, "checkpoint seed changed")
    require(bool(configuration["construction_only"]) == construction_only, "checkpoint construction mode changed")
    results = list(checkpoint["case_results"])
    require(
        checkpoint["resume_audit"]["processed_case_results_sha256"]
        == result_certificate(results),
        "checkpoint case-result certificate failed",
    )
    validate_resume_results(results, target, construction_only)
    return results, float(checkpoint.get("elapsed_seconds", 0.0))


def run(
    input_path: Path,
    output_path: Path,
    orbit: str,
    branch: str,
    case_limit: int | None,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool,
) -> dict:
    started = time.time()
    require(input_path.is_file(), "input affine-hull evidence does not exist")
    require(input_path.resolve() != output_path.resolve(), "input and output paths must differ")
    require(timeout > 0, "timeout must be positive")
    require(workers > 0, "workers must be positive")

    kernel, kernel_audit = compact.kernel_data()
    case_templates, systems_by_orbit, compact_input_audit = compact.audited_inputs()
    orbits, orbit_source = affine_hull.parent.load_z7_orbits()
    leaves_by_orbit, leaf_audit = affine_hull.parent.exact_mean_leaves(orbits)
    require(orbit_source == compact_input_audit["orbit_source"], "reconstructed orbit source changed")
    require(sum(len(rows) for rows in leaves_by_orbit) == EXPECTED_SOURCE_LEAVES, "leaf reconstruction incomplete")
    survivors, evidence_audit = load_and_audit_evidence(
        input_path,
        orbits,
        leaves_by_orbit,
        leaf_audit,
        systems_by_orbit,
        compact_input_audit,
    )
    template_lookup = {
        (int(row["branch_orbit_index"]), str(row["pointed_star_branch"])): row
        for row in case_templates
    }
    require(len(template_lookup) == 4, "compact pointed template count changed")

    selected_unlimited = [
        row
        for row in survivors
        if (orbit == "all" or int(orbit) == int(row["branch_orbit_index"]))
        and (branch == "ALL" or branch == row["pointed_star_branch"])
    ]
    require(selected_unlimited, "selectors chose no affine-hull survivors")
    target = selected_unlimited if case_limit is None else selected_unlimited[:case_limit]
    require(target, "case limit produced an empty target")
    results, prior_elapsed = load_resume_checkpoint(
        output_path,
        evidence_audit,
        selected_unlimited,
        target,
        orbit,
        branch,
        timeout,
        workers,
        seed,
        construction_only,
    )
    reconstruction_audit = {
        "exact_mean_leaf_coverage": leaf_audit,
        "orbit_source": orbit_source,
        "compact_pointed_inputs": compact_input_audit,
        "all_2160_exact_mean_leaves_reconstructed": True,
        "all_four_compact_pointed_systems_audited_mod3_and_mod7": True,
        "mean_variables_fixed_from_reconstructed_leaf_not_input_payload": True,
    }

    for survivor in target[len(results) :]:
        orbit_index = int(survivor["branch_orbit_index"])
        pointed_branch = str(survivor["pointed_star_branch"])
        result = process_survivor(
            survivor,
            template_lookup[(orbit_index, pointed_branch)],
            systems_by_orbit[orbit_index][pointed_branch],
            kernel,
            timeout,
            workers,
            seed + int(survivor["survivor_ordinal"]),
            construction_only,
        )
        results.append(result)
        payload = checkpoint_payload(
            evidence_audit=evidence_audit,
            full_survivors=survivors,
            selected_unlimited=selected_unlimited,
            target=target,
            results=results,
            orbit=orbit,
            branch=branch,
            case_limit=case_limit,
            timeout=timeout,
            workers=workers,
            seed=seed,
            construction_only=construction_only,
            kernel_audit=kernel_audit,
            reconstruction_audit=reconstruction_audit,
            elapsed_seconds=prior_elapsed + time.time() - started,
        )
        compact.pointed.atomic_write(output_path, payload)

    payload = checkpoint_payload(
        evidence_audit=evidence_audit,
        full_survivors=survivors,
        selected_unlimited=selected_unlimited,
        target=target,
        results=results,
        orbit=orbit,
        branch=branch,
        case_limit=case_limit,
        timeout=timeout,
        workers=workers,
        seed=seed,
        construction_only=construction_only,
        kernel_audit=kernel_audit,
        reconstruction_audit=reconstruction_audit,
        elapsed_seconds=prior_elapsed + time.time() - started,
    )
    compact.pointed.atomic_write(output_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--orbit",
        choices=("0", "1", "all"),
        default="all",
        help="select one z=7 line orbit, or both",
    )
    parser.add_argument(
        "--branch",
        type=str.upper,
        choices=("A", "B", "ALL"),
        default="ALL",
        help="select pointed branch A or B, or both",
    )
    parser.add_argument(
        "--case-limit",
        type=parse_positive,
        help="process only the first N selected canonical survivor cases",
    )
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--workers", type=parse_positive, default=16)
    parser.add_argument("--seed", type=int, default=15721001)
    parser.add_argument(
        "--construction-only",
        action="store_true",
        help="construct one fixed-mean model per target without solving",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    result = run(
        args.input,
        args.output,
        args.orbit,
        args.branch,
        args.case_limit,
        args.timeout,
        args.workers,
        args.seed,
        args.construction_only,
    )
    if not args.quiet:
        summary = {
            key: result[key]
            for key in (
                "status",
                "rigorous_conclusion",
                "targeted_survivor_cases",
                "processed_survivor_cases",
                "remaining_survivor_cases",
                "infeasible_survivor_cases",
                "necessary_only_feasible_survivor_cases",
                "unknown_survivor_cases",
                "construction_only_cases",
                "z7_branch_excluded",
                "elapsed_seconds",
            )
        }
        summary["output"] = str(args.output)
        print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
