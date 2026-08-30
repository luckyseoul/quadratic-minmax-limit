#!/usr/bin/env python3
"""Bounded exact four-prime global catalog join for positive p=7, z=7.

This is a strict continuation of
``p7_infinity7_positive_z7_global_catalog_join.py``.  Its input universe is
the 159 orbit-0/branch-A representatives that survived the completed exact
same-index mod-3/mod-7 global join.  The prior 87 rigorous rejections and 78
budget skips are imported verbatim as decisions and are never recomputed.

For each of the 159 candidates, the complete pointed left dependencies are
used modulo 3, 5, 7, and 11.  Every high direction is conditioned out only by
the reduction of the exact zero-mean integral Johnson degree-two lattice.
Every enumerable S/M direction chooses one exact catalog row, with one row
index shared across all four modular components.  The raw residue signatures
are concatenated before exact deduplication and a balanced meet-in-the-middle
join.

The characteristic-five subtlety is handled explicitly.  The naive modular
kernel of the displayed primitive Johnson equations has one spurious mod-5
direction.  A Smith-normal-form audit proves that the integer degree-two
evaluation lattice is saturated; its exact-zero-mean sublattice has rank 20.
The 55 exact catalog differences have rank 20 modulo both 5 and 11, and hence
span the complete liftable constrained space at both primes.

A zero join rigorously rejects the parent candidate.  A positive join remains
necessary only: high catalogs are still replaced by their complete affine
hulls, and modular right-hand-side consistency does not provide a binary edge
lift.  Smoke mode executes a deterministic bounded subset and never transfers
its decisions to the full four-case symmetry classes.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p7_infinity7_positive_z7_global_catalog_join as global37  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine  # noqa: E402
import p7_size_four_slack_classify as johnson  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_global_catalog_join_fourmod"
EXPECTED_GLOBAL_EXPERIMENT = "p7_infinity7_positive_z7_global_catalog_join"
EXPECTED_REPRESENTATIVES = 324
EXPECTED_PRIOR_REJECTIONS = 87
EXPECTED_PRIOR_SURVIVORS = 159
EXPECTED_PRIOR_SKIPS = 78
EXPECTED_NO_ENUMERABLE = 10
DEFAULT_MAX_SIDE_STATES = 6_000_000
DEFAULT_SMOKE_MAX_SIDE_STATES = 200_000
DEFAULT_CHUNK_STATES = 20_000
DEFAULT_CHECKPOINT_EVERY = 10
EXPECTED_MAX_PRIOR_SURVIVOR_RAW_SIDE = 5_531_904
MODULI = (3, 5, 7, 11)
NEW_MODULI = (5, 11)
PRIOR_MODULI = (3, 7)

PRIOR_SURVIVOR_STATUS = "necessary_only_global_catalog_join_survivor"
PRIOR_REJECTION_STATUS = "rigorous_global_catalog_join_rejection"
PRIOR_SKIP_STATUS = "skipped_side_state_budget"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def json_sha256(value: object) -> str:
    return global37.json_sha256(value)


def file_sha256(path: Path) -> str:
    return global37.file_sha256(path)


def canonical_case_digest(rows: list[dict]) -> str:
    return global37.canonical_case_digest(rows)


def catalog_size(catalog_class: str) -> int:
    return global37.catalog_size(catalog_class)


def load_completed_global_input(
    path: Path,
    parent_payload: dict,
    parent_representatives: list[dict],
    parent_provenance: dict,
) -> tuple[dict, list[dict], list[dict], dict, dict]:
    """Validate every imported representative decision and its provenance."""
    raw = path.read_bytes()
    payload = json.loads(raw)
    require(payload["experiment"] == EXPECTED_GLOBAL_EXPERIMENT, "global input experiment changed")
    require(
        payload["status"]
        == "complete_bounded_global_catalog_join_sieve_with_skips_or_survivors",
        "global input is not a completed bounded run",
    )
    require(payload["p"] == 7 and payload["z"] == 7 and payload["phase"] == 0, "global scope changed")
    require(payload["c_H"] == 1 and payload["moduli"] == [3, 7], "global sign/moduli changed")
    require(payload["smoke_test"] is False, "global input is a smoke run")
    require(payload["complete_selected_universe_run"] is True, "global representative run is incomplete")
    require(payload["selected_case_count"] == EXPECTED_REPRESENTATIVES, "global selected census changed")
    require(payload["processed_exact_global_join_cases"] == 246, "global processed census changed")
    require(
        payload["rigorously_rejected_parent_survivors"] == EXPECTED_PRIOR_REJECTIONS,
        "global rejection census changed",
    )
    require(
        payload["necessary_only_global_join_survivors"] == EXPECTED_PRIOR_SURVIVORS,
        "global survivor census changed",
    )
    require(payload["skipped_parent_survivors"] == EXPECTED_PRIOR_SKIPS, "global skip census changed")
    representative_audit = payload["representative_selection"]
    require(representative_audit["enabled"] is True, "global input did not select representatives")
    require(representative_audit["four_case_transfer_claimed"] is True, "global transfer was not certified")
    require(
        payload["symmetry_certificate_audit"]["all_required_transfer_claims_validated"] is True,
        "global input's symmetry validation failed",
    )

    imported_parent = payload["input_provenance"]
    for key in (
        "file_sha256",
        "all_case_results_sha256",
        "survivor_case_keys_sha256",
        "processed_pointed_branch_cases",
        "surviving_pointed_branch_cases",
    ):
        require(imported_parent[key] == parent_provenance[key], f"global parent provenance changed: {key}")
    require(
        parent_payload["all_case_results_sha256"] == imported_parent["all_case_results_sha256"],
        "parent case certificate disagrees with global input",
    )

    rows = payload["case_results"]
    require(len(rows) == EXPECTED_REPRESENTATIVES, "global case-result length changed")
    require(payload["case_results_sha256"] == canonical_case_digest(rows), "global case digest failed")
    expected_keys = [str(row["case_key"]) for row in parent_representatives]
    observed_keys = [str(row["case_key"]) for row in rows]
    require(observed_keys == expected_keys, "global representative order changed")
    require(payload["selected_case_keys_sha256"] == json_sha256(observed_keys), "global selected-key hash failed")

    parent_by_key = {str(row["case_key"]): row for row in parent_representatives}
    decision_counts = Counter()
    raw_side_histogram = Counter()
    candidates: list[dict] = []
    no_enumerable = []
    for row in rows:
        certificate = str(row["decision_certificate_sha256"])
        body = {key: value for key, value in row.items() if key != "decision_certificate_sha256"}
        require(certificate == json_sha256(body), f"global decision certificate failed: {row['case_key']}")
        parent = parent_by_key[str(row["case_key"])]
        require(int(row["branch_orbit_index"]) == 0, "global row is not orbit 0")
        require(str(row["pointed_star_branch"]) == "A", "global row is not branch A")
        require(tuple(row["scaled_means"]) == tuple(parent["scaled_means"]), "global means changed")
        require(tuple(row["catalog_classes"]) == tuple(parent["catalog_classes"]), "global classes changed")
        require(
            tuple(row["high_directions_relaxed_to_exact_affine_hulls"])
            == tuple(parent["high_directions_relaxed_to_exact_affine_hulls"]),
            "global high directions changed",
        )
        require(
            tuple(row["enumerated_directions_joined_globally"])
            == tuple(parent["enumerated_directions"]),
            "global enumerable directions changed",
        )
        status = str(row["decision_status"])
        decision_counts[status] += 1
        if status == PRIOR_REJECTION_STATUS:
            require(row["rigorously_rejected"] is True, "imported rejection flag changed")
            require(row["exact_global_join_processed"] is True, "imported rejection was not processed")
            require(row["join"]["matching_unique_signature_pairs"] == 0, "imported rejection has a join hit")
            require(row["recovered_same_index_catalog_witness"] is None, "imported rejection has a witness")
        elif status == PRIOR_SKIP_STATUS:
            require(row["skipped"] is True, "imported skip flag changed")
            require(row["exact_global_join_processed"] is False, "imported skip was processed")
        elif status == PRIOR_SURVIVOR_STATUS:
            require(row["necessary_only_survivor"] is True, "imported survivor flag changed")
            require(row["exact_global_join_processed"] is True, "imported survivor was not processed")
            require(row["join"]["matching_unique_signature_pairs"] > 0, "imported survivor has no join hit")
            require(row["recovered_same_index_catalog_witness"] is not None, "imported survivor lacks witness")
            raw_side = int(
                row["complete_catalog_balanced_partition_for_comparison"][
                    "maximum_projected_side_product"
                ]
            )
            raw_side_histogram[raw_side] += 1
            require(raw_side <= DEFAULT_MAX_SIDE_STATES, "prior survivor raw side exceeds six million")
            candidates.append(row)
            if not row["enumerated_directions_joined_globally"]:
                no_enumerable.append(str(row["case_key"]))
        else:
            raise AssertionError(f"unexpected global decision status: {status}")

    require(decision_counts[PRIOR_REJECTION_STATUS] == EXPECTED_PRIOR_REJECTIONS, "prior rejection count failed")
    require(decision_counts[PRIOR_SURVIVOR_STATUS] == EXPECTED_PRIOR_SURVIVORS, "prior survivor count failed")
    require(decision_counts[PRIOR_SKIP_STATUS] == EXPECTED_PRIOR_SKIPS, "prior skip count failed")
    require(len(no_enumerable) == EXPECTED_NO_ENUMERABLE, "no-enumerable survivor census changed")
    require(max(raw_side_histogram) == EXPECTED_MAX_PRIOR_SURVIVOR_RAW_SIDE, "maximum raw side changed")

    imported_script_hash = payload["source_provenance"]["this_script_sha256"]
    local_helper_hash = file_sha256(
        ROOT / "scripts" / "p7_infinity7_positive_z7_global_catalog_join.py"
    )
    provenance = {
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "experiment": payload["experiment"],
        "status": payload["status"],
        "case_results_sha256": payload["case_results_sha256"],
        "run_identity_sha256": payload["run_identity_sha256"],
        "source_script_sha256": imported_script_hash,
        "current_local_helper_sha256": local_helper_hash,
        "current_local_helper_byte_identical_to_producing_script": (
            local_helper_hash == imported_script_hash
        ),
    }
    audit = {
        "all_324_decision_certificates_recomputed": True,
        "representative_case_order_matches_parent": True,
        "prior_decision_histogram": dict(sorted(decision_counts.items())),
        "candidate_survivors_selected_for_four_prime_refinement": len(candidates),
        "no_enumerable_affine_only_candidates_included": len(no_enumerable),
        "no_enumerable_case_keys": no_enumerable,
        "prior_survivor_complete_catalog_balanced_side_histogram": {
            str(key): value for key, value in sorted(raw_side_histogram.items())
        },
        "largest_prior_survivor_undeduplicated_balanced_side": max(raw_side_histogram),
        "all_prior_survivor_undeduplicated_balanced_sides_at_most_six_million": True,
        "prior_rejections_are_preserved_without_recomputation": True,
        "prior_budget_skips_are_preserved_without_recomputation": True,
        "completed_run_is_bound_to_its_own_producing_script_hash": True,
        "current_helper_outputs_are_independently_reconstructed_and_audited": True,
    }
    return payload, rows, candidates, provenance, audit


def build_exact_zero_mean_hulls(rebuilt: dict) -> tuple[dict[int, np.ndarray], dict]:
    """Prove the complete liftable integral hull, including characteristic 5."""
    from sympy.matrices.normalforms import smith_normal_form
    from sympy.polys.domains import ZZ

    points, monomials, evaluation, _left_kernel = johnson.johnson_space()
    require(tuple(tuple(int(v) for v in point) for point in points) == affine.POINTS, "Johnson point order changed")
    feature_matrix = np.asarray(evaluation.tolist(), dtype=np.int64)
    require(feature_matrix.shape == (35, 29), "Johnson evaluation shape changed")
    require(int(evaluation.rank()) == 21, "exact Johnson feature rank changed")
    smith = smith_normal_form(evaluation, domain=ZZ)
    invariants = [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index] != 0
    ]
    require(len(invariants) == 21 and invariants == [1] * 21, "Johnson evaluation lattice is not saturated")

    column_sums = np.asarray(feature_matrix.sum(axis=0), dtype=np.int64)
    nonzero_sum_gcd = math.gcd(*(abs(int(value)) for value in column_sums if value))
    require(nonzero_sum_gcd == 5, "degree-two exact mean lattice content changed")
    divided_mean_gcd = math.gcd(*(abs(int(value // 5)) for value in column_sums if value))
    require(divided_mean_gcd == 1, "divided exact mean functional is not primitive")

    source_catalog = affine.canonical_catalog(7, 4).astype(np.int64)
    differences = np.ascontiguousarray(source_catalog[1:] - source_catalog[0], dtype=np.int64)
    kernel_rows = np.asarray(rebuilt["kernel_rows"], dtype=np.int64)
    from sympy import Matrix

    exact_kernel_rank = int(Matrix(kernel_rows.tolist()).rank())
    require(exact_kernel_rank == 14, "exact Johnson left-kernel rank changed")
    require(not np.any(kernel_rows @ feature_matrix), "exact Johnson kernel does not annihilate features")
    require(source_catalog.shape == (56, 35), "source hull catalog changed")
    require(differences.shape == (55, 35), "source hull differences changed")
    require(not np.any(differences.sum(axis=1)), "hull differences lost exact zero mean")
    require(not np.any(differences % 2), "hull differences lost exact parity")
    require(not np.any(kernel_rows @ differences.T), "hull differences left exact degree two")

    constraints = np.vstack((np.ones((1, 35), dtype=np.int64), kernel_rows))
    bases: dict[int, np.ndarray] = {}
    prime_audits = []
    for modulus in MODULI:
        difference_rank = affine.modular_rank(differences, modulus)
        require(difference_rank == 20, f"mod-{modulus} exact catalog-difference rank changed")
        basis, basis_indices = affine.independent_exact_rows(differences, modulus)
        require(affine.modular_rank(basis, modulus) == 20, f"mod-{modulus} selected hull basis lost rank")
        require(not np.any(basis.sum(axis=1)), "selected hull basis lost exact zero mean")
        require(not np.any(kernel_rows @ basis.T), "selected hull basis left exact degree two")
        naive_constraint_rank = affine.modular_rank(constraints, modulus)
        naive_nullspace, nullspace_rank = affine.modular_right_nullspace(constraints, modulus)
        require(naive_constraint_rank == nullspace_rank, "naive constraint rank audit disagrees")
        feature_rank = affine.modular_rank(feature_matrix, modulus)
        require(feature_rank == 21, f"mod-{modulus} saturated feature reduction lost rank")
        bases[modulus] = basis
        prime_audits.append(
            {
                "modulus": modulus,
                "exact_catalog_difference_rank": difference_rank,
                "selected_exact_basis_rank": affine.modular_rank(basis, modulus),
                "selected_source_difference_indices": list(basis_indices),
                "selected_exact_basis_sha256_int64": affine.matrix_sha256(basis),
                "saturated_degree_two_lattice_reduction_rank": feature_rank,
                "exact_zero_mean_sublattice_reduction_dimension": 20,
                "catalog_differences_span_full_liftable_exact_zero_mean_degree_two_space": True,
                "naive_primitive_kernel_plus_mean_rank": naive_constraint_rank,
                "naive_modular_constraint_nullity": len(naive_nullspace),
                "naive_modular_kernel_equals_liftable_exact_space": len(naive_nullspace) == 20,
                "characteristic_five_spurious_naive_kernel_dimension": (
                    len(naive_nullspace) - 20 if modulus == 5 else 0
                ),
            }
        )

    require(
        next(row for row in prime_audits if row["modulus"] == 5)[
            "characteristic_five_spurious_naive_kernel_dimension"
        ]
        == 1,
        "expected one non-liftable naive mod-5 direction",
    )
    for modulus in PRIOR_MODULI:
        require(
            np.array_equal(bases[modulus], rebuilt["hull_bases"][modulus]),
            f"mod-{modulus} exact hull basis changed from parent reconstruction",
        )

    audit = {
        "johnson_points": 35,
        "degree_at_most_two_feature_columns": len(monomials),
        "exact_rational_feature_rank": 21,
        "exact_primitive_left_kernel_rank": exact_kernel_rank,
        "exact_primitive_left_kernel_annihilates_feature_matrix": True,
        "feature_matrix_sha256_int64": affine.matrix_sha256(feature_matrix),
        "smith_nonzero_invariant_factors": invariants,
        "all_smith_nonzero_invariant_factors_are_one": True,
        "integer_degree_two_evaluation_lattice_is_saturated": True,
        "exact_coordinate_sum_image_gcd": nonzero_sum_gcd,
        "primitive_divided_mean_functional_gcd": divided_mean_gcd,
        "exact_zero_mean_sublattice_rank": 20,
        "exact_zero_mean_sublattice_is_saturated": True,
        "source_catalog_rows": len(source_catalog),
        "source_catalog_sha256_int16": affine.matrix_sha256(source_catalog.astype(np.int16)),
        "all_difference_rows": len(differences),
        "all_differences_sha256_int64": affine.matrix_sha256(differences),
        "prime_audits": prime_audits,
        "mod5_full_space_semantics": (
            "full constrained space means the reduction of the liftable exact integral "
            "zero-mean degree-two lattice, not the one-dimension-larger naive kernel of "
            "the rank-dropping primitive equations modulo 5"
        ),
        "proof": (
            "Smith invariant factors all equal one, so the integer feature-column image is "
            "the saturated lattice of all integral vectors in the rational Johnson degree-two "
            "space. Its exact-zero-mean kernel is saturated of rank 20. The exact catalog "
            "differences lie in that kernel and have modular rank 20, hence span its complete "
            "reduction modulo each audited prime."
        ),
    }
    return bases, audit


def augment_complete_pointed_dependencies(
    rebuilt: dict,
    symmetry_payload: dict,
) -> dict:
    """Construct direct mod-5/mod-11 left-null bases on all four exact systems."""
    certificate_rows = {
        (int(row["branch_orbit_index"]), str(row["pointed_star_branch"])): row
        for row in symmetry_payload["canonical_compact_dependency_audit"]["four_case_audits"]
    }
    require(set(certificate_rows) == {(o, b) for o in range(2) for b in ("A", "B")}, "certificate system rows changed")
    system_rows = []
    for orbit_index in range(2):
        for branch in ("A", "B"):
            system = rebuilt["systems"][orbit_index][branch]
            matrix = np.asarray(system["matrix"], dtype=np.int64)
            certificate = certificate_rows[(orbit_index, branch)]
            require(certificate["case_key"] == f"orbit{orbit_index}_{branch}", "certificate case key changed")
            require(certificate["pointed_matrix_shape"] == list(matrix.shape), "certificate matrix shape changed")
            require(
                certificate["pointed_matrix_sha256_int16"]
                == affine.matrix_sha256(matrix.astype(np.int16)),
                "certificate pointed matrix hash changed",
            )
            require(
                certificate["base_rhs_sha256_int64"]
                == affine.matrix_sha256(np.asarray(system["base_rhs"], dtype=np.int64)[None, :]),
                "certificate base RHS hash changed",
            )
            expected_by_prime = {
                int(row["modulus"]): row for row in certificate["prime_audits"]
            }
            require(set(expected_by_prime) == set(MODULI), "certificate prime audit set changed")
            prime_rows = []
            calibration = (17 * np.arange(matrix.shape[1], dtype=np.int64) + 3) % 2
            for modulus in MODULI:
                if modulus in NEW_MODULI:
                    dependencies, rank_from_nullspace = affine.modular_right_nullspace(
                        matrix.T, modulus
                    )
                    dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
                    independently_recomputed_rank = affine.modular_rank(matrix, modulus)
                    require(
                        independently_recomputed_rank == rank_from_nullspace,
                        f"mod-{modulus} pointed rank/nullspace disagreement",
                    )
                    system["moduli"][modulus] = {
                        "dependencies": dependencies,
                        "rank": independently_recomputed_rank,
                    }
                    source = "direct_new_complete_left_nullspace_on_exact_pointed_matrix"
                else:
                    dependencies = np.asarray(
                        system["moduli"][modulus]["dependencies"], dtype=np.int64
                    )
                    independently_recomputed_rank = int(system["moduli"][modulus]["rank"])
                    rank_from_nullspace = matrix.shape[0] - len(dependencies)
                    source = "parent_reconstructed_complete_left_nullspace_reaudited"

                expected = expected_by_prime[modulus]
                require(
                    independently_recomputed_rank == int(expected["edge_matrix_rank"]),
                    f"mod-{modulus} pointed rank disagrees with symmetry certificate",
                )
                dependency_dimension = matrix.shape[0] - independently_recomputed_rank
                require(
                    dependencies.shape == (dependency_dimension, matrix.shape[0]),
                    f"mod-{modulus} complete left dependency shape changed",
                )
                require(
                    dependency_dimension == int(expected["complete_left_dependency_dimension"]),
                    f"mod-{modulus} dependency dimension disagrees with certificate",
                )
                require(
                    affine.modular_rank(dependencies, modulus) == dependency_dimension,
                    f"mod-{modulus} complete left dependencies lost row rank",
                )
                require(
                    not np.any(dependencies @ (matrix % modulus) % modulus),
                    f"mod-{modulus} complete left-null identity failed",
                )
                manufactured_rhs = matrix @ calibration % modulus
                manufactured_syndrome = dependencies @ manufactured_rhs % modulus
                require(not np.any(manufactured_syndrome), "manufactured consistent RHS rejected")
                first_dependency = dependencies[0] % modulus
                nonzero_coordinates = np.flatnonzero(first_dependency)
                require(len(nonzero_coordinates) > 0, "dependency basis contains zero row")
                corrupted_coordinate = int(nonzero_coordinates[0])
                inconsistent_rhs = manufactured_rhs.copy()
                inconsistent_rhs[corrupted_coordinate] += 1
                inconsistent_syndrome = dependencies @ inconsistent_rhs % modulus
                require(np.any(inconsistent_syndrome), "manufactured inconsistent RHS was accepted")
                prime_rows.append(
                    {
                        "modulus": modulus,
                        "construction_source": source,
                        "pointed_matrix_rank": independently_recomputed_rank,
                        "rank_returned_by_complete_left_nullspace": rank_from_nullspace,
                        "complete_left_dependency_dimension": dependency_dimension,
                        "complete_left_dependency_basis_rank": affine.modular_rank(
                            dependencies, modulus
                        ),
                        "dependency_sha256_uint8": affine.matrix_sha256(
                            dependencies.astype(np.uint8)
                        ),
                        "complete_left_null_identity": True,
                        "rank_nullity_proves_dependency_completeness": True,
                        "manufactured_binary_edge_vector_sha256_uint8": affine.matrix_sha256(
                            calibration[None, :].astype(np.uint8)
                        ),
                        "manufactured_consistent_rhs_sha256_uint8": affine.matrix_sha256(
                            manufactured_rhs[None, :].astype(np.uint8)
                        ),
                        "manufactured_consistent_rhs_syndrome_zero": True,
                        "manufactured_inconsistent_rhs_corrupted_coordinate": corrupted_coordinate,
                        "manufactured_inconsistent_rhs_detected": True,
                        "matches_explicit_symmetry_certificate_rank_and_dimension": True,
                    }
                )
            system_rows.append(
                {
                    "branch_orbit_index": orbit_index,
                    "pointed_star_branch": branch,
                    "matrix_shape": list(matrix.shape),
                    "matrix_sha256_int16": affine.matrix_sha256(matrix.astype(np.int16)),
                    "prime_audits": prime_rows,
                }
            )
    return {
        "exact_pointed_systems_audited": len(system_rows),
        "new_moduli_directly_constructed": list(NEW_MODULI),
        "all_four_orbit_branch_systems_constructed_at_new_moduli": True,
        "all_complete_left_nullspaces_rank_nullity_and_manufactured_rhs_audited": True,
        "systems": system_rows,
        "all_system_audits_sha256": json_sha256(system_rows),
    }


def validate_symmetry_certificate(
    path: Path,
    payload: dict,
    parent_survivors: list[dict],
    parent_representatives: list[dict],
    parent_provenance: dict,
    hull_audit: dict,
) -> tuple[dict, dict[str, tuple[str, ...]]]:
    raw = path.read_bytes()
    certificate = json.loads(raw)
    require(
        certificate["experiment"] == "p7_infinity7_positive_z7_compact_symmetry_audit",
        "symmetry certificate experiment changed",
    )
    require(
        certificate["status"]
        == "complete_exact_four_case_compact_symmetry_and_survivor_partition_audit",
        "symmetry certificate is incomplete",
    )
    require(certificate["p"] == 7 and certificate["z"] == 7, "symmetry scope changed")
    require(certificate["phase"] == 0 and certificate["c_H"] == 1, "symmetry sign/phase changed")
    require(certificate["moduli"] == list(MODULI), "symmetry certificate is not four-prime")
    require(certificate["solver_invoked"] is False, "symmetry certificate unexpectedly invoked a solver")
    require(certificate["all_required_structural_audits_passed"] is True, "symmetry structural audit failed")

    partition = certificate["optional_affine_survivor_partition_audit"]
    require(partition["performed"] is True and partition["survivor_partition_claimed"] is True, "survivor partition missing")
    require(partition["survivor_count"] == len(parent_survivors) == 1_296, "symmetry survivor census changed")
    evidence = partition["input_evidence_audit"]
    require(evidence["file_sha256"] == parent_provenance["file_sha256"], "symmetry parent file hash changed")
    require(
        evidence["all_case_results_sha256"] == parent_provenance["all_case_results_sha256"],
        "symmetry parent case hash changed",
    )
    require(
        evidence["survivor_case_keys_sha256"] == parent_provenance["survivor_case_keys_sha256"],
        "symmetry parent survivor hash changed",
    )

    compact = certificate["canonical_compact_dependency_audit"]
    require(compact["all_four_augmented_compact_row_spaces_identical"] is True, "compact row spaces differ")
    common = {int(row["modulus"]): row for row in compact["common_prime_audits"]}
    require(set(common) == set(MODULI), "compact common-prime audit changed")
    require(
        all(row["all_four_canonical_row_spaces_identical"] is True for row in common.values()),
        "a four-prime compact row-space identity failed",
    )
    affine_maps = certificate["affine_compact_isomorphism_audit"]
    require(affine_maps["all_four_compact_cases_are_exactly_isomorphic"] is True, "compact isomorphism failed")
    for mapping_row in affine_maps["maps"]:
        prime_rows = {int(row["modulus"]): row for row in mapping_row["augmented_modular_row_space_transport"]}
        require(set(prime_rows) == set(MODULI), "affine map omitted a prime")
        require(all(row["equals_target_augmented_row_space"] is True for row in prime_rows.values()), "affine modular transport failed")
        require(mapping_row["all_14_exact_primitive_kernel_equations_per_direction_transport"] is True, "exact kernel transport failed")
        require(mapping_row["directional_sum_and_exact_mean_identity_transport"] is True, "exact mean transport failed")

    transfer = certificate["global_catalog_join_transfer"]
    require(transfer["scope"] == "positive_p7_z7_1296_affine_hull_survivors", "transfer scope changed")
    require(transfer["equivalence_class_count"] == 324, "symmetry class count changed")
    require(transfer["equivalence_class_size"] == 4, "symmetry class size changed")
    require(transfer["all_1296_survivors_partitioned_into_324_disjoint_complete_classes"] is True, "class partition failed")
    require(transfer["transfer_valid_for_global_same_index_catalog_join"] is True, "global transfer claim failed")

    catalog_transport = transfer["catalog_row_transport_audit"]
    require(catalog_transport["all_complete_S_M_catalogs_transported_as_exact_row_sets"] is True, "exact catalog transport failed")
    require(
        all(
            row["exact_row_set_transport_bijective"] is True
            and row["anchor_base_and_catalog_contribution_shifts_cancel_exactly"] is True
            for row in catalog_transport["domain_audits"]
        ),
        "a characteristic-independent exact catalog transport failed",
    )
    high_transport = transfer["high_affine_hull_and_anchor_target_audit"]
    require(high_transport["all_mapped_anchors_have_exact_target_parity_mean_and_degree_two"] is True, "anchor transport failed")
    require(high_transport["mapped_high_anchor_plus_full_hull_equals_target_anchor_plus_full_hull"] is True, "high affine transport failed")
    require(hull_audit["integer_degree_two_evaluation_lattice_is_saturated"] is True, "new hull proof missing")
    require(
        all(
            row["catalog_differences_span_full_liftable_exact_zero_mean_degree_two_space"]
            for row in hull_audit["prime_audits"]
            if int(row["modulus"]) in NEW_MODULI
        ),
        "new-prime full hull proof failed",
    )

    expected_members = {str(row["case_key"]) for row in parent_survivors}
    expected_representatives = {str(row["case_key"]) for row in parent_representatives}
    observed_members: set[str] = set()
    observed_representatives: set[str] = set()
    mapping: dict[str, tuple[str, ...]] = {}
    canonical_classes = []
    for row in transfer["equivalence_classes"]:
        representative = str(row["representative_case_key"])
        members = tuple(str(value) for value in row["member_case_keys"])
        require(row["global_same_index_join_equivalence_proved"] is True, "class equivalence failed")
        require(len(members) == len(set(members)) == 4, "class member count changed")
        require(representative in members, "class omits its representative")
        require(representative.startswith("orbit0_leaf") and representative.endswith("_branchA"), "class representative changed")
        require(not observed_members.intersection(members), "symmetry classes overlap")
        observed_members.update(members)
        observed_representatives.add(representative)
        mapping[representative] = members
        canonical_classes.append(
            {"representative_case_key": representative, "member_case_keys": list(members)}
        )
    require(observed_members == expected_members, "symmetry classes do not cover parent survivors exactly")
    require(observed_representatives == expected_representatives, "symmetry representatives changed")
    require(
        transfer["parent_survivor_case_keys_sha256"] == parent_provenance["survivor_case_keys_sha256"],
        "transfer survivor hash changed",
    )
    require(transfer["equivalence_classes_sha256"] == json_sha256(canonical_classes), "equivalence class hash failed")

    audit = {
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "experiment": certificate["experiment"],
        "status": certificate["status"],
        "moduli": list(MODULI),
        "equivalence_class_count": len(mapping),
        "class_size": 4,
        "covered_parent_survivors": len(observed_members),
        "equivalence_classes_sha256": transfer["equivalence_classes_sha256"],
        "all_four_compact_augmented_row_spaces_identical": True,
        "exact_catalog_row_bijections_are_characteristic_independent": True,
        "new_mod5_mod11_high_hull_transport_follows_from_exact_coordinate_transport_and_full_lattice_hull_proof": True,
        "transfer_valid_for_this_four_prime_global_join": True,
        "transfer_claimed_only_after_a_complete_non_smoke_representative_run": True,
    }
    return audit, mapping


class PrimeProjection:
    """Condition one complete pointed dependency space by exact high hulls."""

    def __init__(self, system: dict, modulus: int, hull_basis: np.ndarray):
        self.system = system
        self.modulus = int(modulus)
        self.dependencies = np.asarray(system["moduli"][modulus]["dependencies"], dtype=np.int64)
        self.hull_basis = np.asarray(hull_basis, dtype=np.int64)
        self.conditioners: dict[tuple[int, ...], tuple[np.ndarray, dict]] = {}

    def direction_block(self, direction: int) -> np.ndarray:
        require(0 <= direction < 8, "direction index out of range")
        return self.dependencies[:, 1 + 35 * direction : 1 + 35 * (direction + 1)]

    def conditioner(self, high: tuple[int, ...]) -> tuple[np.ndarray, dict]:
        high = tuple(sorted(int(value) for value in high))
        require(len(high) == len(set(high)), "high direction repeated")
        if high not in self.conditioners:
            dependency_dimension = len(self.dependencies)
            if high:
                images = [
                    self.direction_block(direction)
                    @ (self.hull_basis.T % self.modulus)
                    % self.modulus
                    for direction in high
                ]
                hull_image = np.ascontiguousarray(np.concatenate(images, axis=1), dtype=np.int64)
                coefficients, image_rank = affine.modular_right_nullspace(
                    hull_image.T, self.modulus
                )
            else:
                hull_image = np.empty((dependency_dimension, 0), dtype=np.int64)
                coefficients = np.eye(dependency_dimension, dtype=np.int64)
                image_rank = 0
            coefficients = np.ascontiguousarray(coefficients, dtype=np.int64)
            conditioned = np.ascontiguousarray(
                coefficients @ self.dependencies % self.modulus, dtype=np.int64
            )
            dimension = len(coefficients)
            require(dimension == dependency_dimension - image_rank, "conditioner rank-nullity failed")
            require(affine.modular_rank(coefficients, self.modulus) == dimension, "conditioner lost rank")
            require(not np.any(coefficients @ hull_image % self.modulus), "conditioner missed a high hull")
            require(affine.modular_rank(conditioned, self.modulus) == dimension, "conditioned dependencies lost rank")
            require(
                not np.any(
                    conditioned
                    @ (np.asarray(self.system["matrix"], dtype=np.int64) % self.modulus)
                    % self.modulus
                ),
                "conditioned dependency left the pointed nullspace",
            )
            metadata = {
                "modulus": self.modulus,
                "high_directions": list(high),
                "high_direction_count": len(high),
                "complete_dependency_dimension": dependency_dimension,
                "exact_hull_basis_rows_per_high_direction": len(self.hull_basis),
                "concatenated_hull_image_columns": int(hull_image.shape[1]),
                "concatenated_hull_image_rank": int(image_rank),
                "conditioned_dependency_dimension": dimension,
                "coefficient_sha256_int64": affine.matrix_sha256(coefficients),
                "hull_image_sha256_int64": affine.matrix_sha256(hull_image),
                "conditioned_dependencies_sha256_uint8": affine.matrix_sha256(
                    conditioned.astype(np.uint8)
                ),
                "rank_nullity_audited": True,
                "full_exact_high_hulls_annihilated": True,
                "conditioned_complete_left_nullspace_audited": True,
            }
            self.conditioners[high] = coefficients, metadata
        return self.conditioners[high]

    def conditioned_direction_block(self, high: tuple[int, ...], direction: int) -> np.ndarray:
        coefficients, _metadata = self.conditioner(high)
        return np.ascontiguousarray(
            coefficients @ self.direction_block(direction) % self.modulus,
            dtype=np.int64,
        )

    def audit(self) -> dict:
        rows = [self.conditioners[key][1] for key in sorted(self.conditioners)]
        histogram = Counter(
            (row["high_direction_count"], row["conditioned_dependency_dimension"])
            for row in rows
        )
        return {
            "modulus": self.modulus,
            "complete_dependency_dimension": len(self.dependencies),
            "conditioner_count": len(rows),
            "conditioned_dimension_histogram": {
                f"H{high}_D{dimension}": count
                for (high, dimension), count in sorted(histogram.items())
            },
            "all_conditioners_exact_rank_nullity_hull_annihilation_and_left_null_audited": True,
            "conditioners": rows,
        }


class FourPrimeProjectionFactory:
    def __init__(self, system: dict, hull_bases: dict[int, np.ndarray], anchors: object):
        self.system = system
        self.anchors = anchors
        self.primes = {
            modulus: PrimeProjection(system, modulus, hull_bases[modulus])
            for modulus in MODULI
        }
        self.signature_records: dict[tuple, dict] = {}

    def anchor_rhs_and_raw_syndromes(
        self, orbit: dict, leaf: dict
    ) -> tuple[np.ndarray, dict[int, np.ndarray]]:
        rhs = np.asarray(self.system["base_rhs"], dtype=np.int64).copy()
        for direction, (mask, mean) in enumerate(
            zip(orbit["masks"], leaf["scaled_means"])
        ):
            anchor = self.anchors.get(int(mask), int(mean))
            block = slice(1 + 35 * direction, 1 + 35 * (direction + 1))
            rhs[block] = 13 - anchor
        raw = {
            modulus: np.ascontiguousarray(
                self.primes[modulus].dependencies @ rhs % modulus, dtype=np.int64
            )
            for modulus in MODULI
        }
        return rhs, raw

    def projected_bases(
        self, raw_syndromes: dict[int, np.ndarray], high: tuple[int, ...]
    ) -> dict[int, np.ndarray]:
        bases = {}
        for modulus in MODULI:
            coefficients, _metadata = self.primes[modulus].conditioner(high)
            bases[modulus] = np.ascontiguousarray(
                coefficients @ raw_syndromes[modulus] % modulus,
                dtype=np.uint8,
            )
        return bases

    def catalog_signature(
        self,
        high: tuple[int, ...],
        direction: int,
        mask: int,
        mean: int,
    ) -> tuple[np.ndarray, dict]:
        key = (tuple(high), int(direction), int(mask), int(mean))
        catalog = affine.mapped_catalog(int(mask), int(mean)).astype(np.int64)
        anchor = self.anchors.get(int(mask), int(mean))
        delta = np.ascontiguousarray(anchor[None, :] - catalog, dtype=np.int64)
        require(not np.any(delta.sum(axis=1)), "catalog delta lost exact zero mean")
        require(not np.any(delta % 2), "catalog delta lost exact parity")
        require(not np.any(self.anchors.kernel_rows @ delta.T), "catalog delta left exact degree two")
        components = []
        dimensions = {}
        hashes = {}
        nonzero_entries = {}
        offsets = {}
        offset = 0
        for modulus in MODULI:
            block = self.primes[modulus].conditioned_direction_block(high, direction)
            signatures = block @ (delta.T % modulus) % modulus
            signatures = np.ascontiguousarray(signatures.T, dtype=np.uint8)
            require(len(signatures) == len(catalog), "component signature lost a catalog row")
            components.append(signatures)
            dimensions[modulus] = signatures.shape[1]
            hashes[modulus] = affine.matrix_sha256(signatures)
            nonzero_entries[modulus] = int(np.count_nonzero(signatures))
            offsets[modulus] = (offset, offset + signatures.shape[1])
            offset += signatures.shape[1]
        combined = np.ascontiguousarray(np.concatenate(components, axis=1), dtype=np.uint8)
        legacy37 = np.ascontiguousarray(
            np.concatenate((components[0], components[2]), axis=1), dtype=np.uint8
        )
        identifier = (
            f"omit-{','.join(map(str, high)) or 'none'}_d{direction}_mask{mask}_mean{mean}"
        )
        record = {
            "identifier": identifier,
            "high_directions": list(high),
            "direction": int(direction),
            "mask": int(mask),
            "scaled_mean": int(mean),
            "catalog_rows": len(catalog),
            "exact_catalog_sha256_int16": affine.matrix_sha256(catalog.astype(np.int16)),
            "anchor_sha256_int64": affine.matrix_sha256(anchor[None, :]),
            "delta_sha256_int64": affine.matrix_sha256(delta),
            "component_dimensions": {str(p): dimensions[p] for p in MODULI},
            "component_offsets": {str(p): list(offsets[p]) for p in MODULI},
            "component_signature_sha256_uint8": {str(p): hashes[p] for p in MODULI},
            "component_nonzero_entry_count": {
                str(p): nonzero_entries[p] for p in MODULI
            },
            "component_signature_is_identically_zero": {
                str(p): nonzero_entries[p] == 0 for p in MODULI
            },
            "prior_mod3_mod7_concatenated_signature_sha256_uint8": affine.matrix_sha256(
                legacy37
            ),
            "four_prime_concatenated_signature_sha256_uint8": affine.matrix_sha256(combined),
            "same_exact_catalog_matrix_and_row_index_used_for_all_four_primes": True,
        }
        prior = self.signature_records.setdefault(key, record)
        require(prior == record, "recomputed four-prime signature changed")
        return combined, record

    def audit(self) -> dict:
        records = [self.signature_records[key] for key in sorted(self.signature_records)]
        return {
            "pointed_star_branch": self.system["branch"],
            "modulus_conditioners": [self.primes[p].audit() for p in MODULI],
            "unique_complete_catalog_signatures_computed": len(records),
            "all_signature_records_sha256": json_sha256(records),
            "all_signatures_use_one_shared_exact_catalog_index_across_3_5_7_11": True,
        }


class FourPrimeContributionStore:
    """Deduplicate complete four-prime contribution signatures exactly."""

    def __init__(self, factory: FourPrimeProjectionFactory):
        self.factory = factory
        self.cache: dict[tuple, dict] = {}
        self.records: dict[tuple, dict] = {}

    def get(
        self,
        high: tuple[int, ...],
        direction: int,
        mask: int,
        mean: int,
        catalog_class: str,
    ) -> dict:
        key = (tuple(high), int(direction), int(mask), int(mean))
        if key not in self.cache:
            complete, source = self.factory.catalog_signature(
                high, direction, mask, mean
            )
            expected_size = catalog_size(catalog_class)
            require(len(complete) == expected_size, "complete catalog size changed")
            unique, first_indices = global37.deduplicate_contribution_rows(complete)
            require(np.all(first_indices < expected_size), "catalog representative index escaped")
            metadata = {
                **source,
                "catalog_class": catalog_class,
                "complete_catalog_rows": expected_size,
                "four_prime_projected_unique_contribution_states": len(unique),
                "discarded_four_prime_projected_duplicate_rows": expected_size - len(unique),
                "sorted_unique_four_prime_signature_sha256_uint8": affine.matrix_sha256(unique),
                "representative_exact_catalog_indices_sha256_uint64": affine.matrix_sha256(
                    first_indices[None, :]
                ),
                "exact_raw_byte_deduplication_before_partition_and_budget": True,
                "multiplicity_retained": False,
            }
            self.cache[key] = {
                "rows": unique,
                "first_catalog_indices": first_indices,
                "metadata": metadata,
            }
            self.records[key] = metadata
        else:
            require(
                self.records[key]["catalog_class"] == catalog_class,
                "one contribution key acquired two catalog classes",
            )
        return self.cache[key]

    def audit(self) -> dict:
        rows = [self.records[key] for key in sorted(self.records)]
        compression = Counter(
            (
                len(row["high_directions"]),
                row["catalog_class"],
                row["complete_catalog_rows"],
                row["four_prime_projected_unique_contribution_states"],
            )
            for row in rows
        )
        return {
            "unique_projected_contributions": len(rows),
            "complete_catalog_rows_total": sum(row["complete_catalog_rows"] for row in rows),
            "four_prime_unique_states_total": sum(
                row["four_prime_projected_unique_contribution_states"] for row in rows
            ),
            "compression_histogram": {
                f"H{h}_{c}_{n}_to_{u}": count
                for (h, c, n, u), count in sorted(compression.items())
            },
            "all_complete_signatures_and_representative_indices_sha256": json_sha256(rows),
            "nonzero_complete_component_signature_records_by_modulus": {
                str(modulus): sum(
                    not row["component_signature_is_identically_zero"][str(modulus)]
                    for row in rows
                )
                for modulus in MODULI
            },
            "all_new_prime_complete_contribution_components_identically_zero": all(
                row["component_signature_is_identically_zero"][str(modulus)]
                for row in rows
                for modulus in NEW_MODULI
            ),
            "deduplicated_before_partition_and_budget": True,
            "multiplicity_is_irrelevant_to_existence": True,
            "records": rows,
        }


def make_candidate_plans(
    candidate_rows: list[dict],
    parent_by_key: dict[str, dict],
    rebuilt: dict,
    factory: FourPrimeProjectionFactory,
    contributions: FourPrimeContributionStore,
    max_side_states: int,
) -> tuple[list[dict], dict]:
    plans = []
    raw_histogram = Counter()
    four_histogram = Counter()
    compression_pairs = Counter()
    projected_base_nonzero_cases = Counter()
    new_prime_active_case_keys: list[str] = []
    for prior_row in candidate_rows:
        parent_row = parent_by_key[str(prior_row["case_key"])]
        orbit, leaf, system, _old_factory = global37.validate_parent_survivor(
            parent_row, rebuilt
        )
        require(system is factory.system, "candidate left the representative pointed system")
        high = tuple(int(value) for value in leaf["high_directions"])
        enumerated = tuple(int(value) for value in leaf["enumerated_directions"])
        complete_sizes = tuple(
            catalog_size(str(leaf["catalog_classes"][direction]))
            for direction in enumerated
        )
        projected = {
            direction: contributions.get(
                high,
                direction,
                int(orbit["masks"][direction]),
                int(leaf["scaled_means"][direction]),
                str(leaf["catalog_classes"][direction]),
            )
            for direction in enumerated
        }
        four_sizes = tuple(len(projected[direction]["rows"]) for direction in enumerated)
        prior_sizes = tuple(int(value) for value in prior_row["projected_unique_contribution_sizes"])
        require(len(prior_sizes) == len(four_sizes), "prior/four-prime contribution arity changed")
        require(
            all(old <= new <= complete for old, new, complete in zip(prior_sizes, four_sizes, complete_sizes)),
            "adding mod-5/mod-11 either merged a mod-3/mod-7 class or exceeded the catalog",
        )
        prior_audits = {
            int(row["direction"]): row for row in prior_row["retained_catalog_signature_audits"]
        }
        require(set(prior_audits) == set(enumerated), "prior retained signature audit coverage changed")
        for position, direction in enumerate(enumerated):
            metadata = projected[direction]["metadata"]
            imported = prior_audits[direction]
            require(
                metadata["prior_mod3_mod7_concatenated_signature_sha256_uint8"]
                == imported["complete_signature_sha256_uint8"],
                "recomputed same-index mod-3/mod-7 contribution differs from completed run",
            )
            require(
                metadata["exact_catalog_sha256_int16"] == imported["exact_catalog_sha256_int16"],
                "exact catalog differs from completed run",
            )
            require(
                prior_sizes[position] == imported["projected_unique_contribution_states"],
                "prior contribution compression fields disagree",
            )
            compression_pairs[(prior_sizes[position], four_sizes[position])] += 1

        raw_partition = global37.balanced_partition(enumerated, complete_sizes)
        require(
            raw_partition == prior_row["complete_catalog_balanced_partition_for_comparison"],
            "recomputed undeduplicated balanced partition differs from completed run",
        )
        require(
            int(raw_partition["maximum_projected_side_product"])
            <= DEFAULT_MAX_SIDE_STATES,
            "a prior survivor lacks the six-million undeduplicated-side bound",
        )
        partition = global37.balanced_partition(enumerated, four_sizes)
        require(
            int(partition["maximum_projected_side_product"])
            <= int(raw_partition["maximum_projected_side_product"]),
            "four-prime deduplication exceeded the undeduplicated balanced bound",
        )
        conditioners = {
            modulus: factory.primes[modulus].conditioner(high)[1]
            for modulus in MODULI
        }
        dimensions = {
            modulus: int(conditioners[modulus]["conditioned_dependency_dimension"])
            for modulus in MODULI
        }
        joint_dimension = sum(dimensions.values())
        anchor_rhs, raw_syndromes = factory.anchor_rhs_and_raw_syndromes(orbit, leaf)
        anchor_hash = affine.matrix_sha256(anchor_rhs[None, :])
        require(anchor_hash == parent_row["anchor_rhs_sha256_int64"], "preflight parent anchor changed")
        require(anchor_hash == prior_row["anchor_rhs_sha256_int64"], "preflight prior anchor changed")
        projected_bases = factory.projected_bases(raw_syndromes, high)
        base_nonzero = {
            modulus: bool(np.any(projected_bases[modulus] % modulus))
            for modulus in MODULI
        }
        for modulus, nonzero in base_nonzero.items():
            projected_base_nonzero_cases[modulus] += int(nonzero)
        new_component_nonzero = any(
            not projected[direction]["metadata"][
                "component_signature_is_identically_zero"
            ][str(modulus)]
            for direction in enumerated
            for modulus in NEW_MODULI
        )
        if any(base_nonzero[modulus] for modulus in NEW_MODULI) or new_component_nonzero:
            new_prime_active_case_keys.append(str(prior_row["case_key"]))
        raw_histogram[int(raw_partition["maximum_projected_side_product"])] += 1
        four_histogram[int(partition["maximum_projected_side_product"])] += 1
        eligible = int(partition["maximum_projected_side_product"]) <= max_side_states
        h, s, m = (int(value) for value in leaf["pattern"])
        plans.append(
            {
                "case_key": str(prior_row["case_key"]),
                "prior_row": prior_row,
                "parent_row": parent_row,
                "orbit": orbit,
                "leaf": leaf,
                "system": system,
                "factory": factory,
                "branch_orbit_index": 0,
                "source_orbit_index": int(parent_row["source_orbit_index"]),
                "orbit_leaf_index": int(parent_row["orbit_leaf_index"]),
                "pointed_star_branch": "A",
                "catalog_pattern": f"H{h}_S{s}_M{m}",
                "high_directions": high,
                "enumerated_directions": enumerated,
                "complete_catalog_sizes": complete_sizes,
                "prior_mod3_mod7_unique_contribution_sizes": prior_sizes,
                "four_prime_unique_contribution_sizes": four_sizes,
                "projected_contributions": projected,
                "partition": partition,
                "complete_catalog_balanced_partition": raw_partition,
                "conditioned_dimensions": dimensions,
                "joint_signature_dimension": joint_dimension,
                "preflight_projected_base_sha256_by_modulus": {
                    str(p): affine.matrix_sha256(projected_bases[p][None, :])
                    for p in MODULI
                },
                "preflight_projected_base_nonzero_by_modulus": {
                    str(p): base_nonzero[p] for p in MODULI
                },
                "estimated_four_prime_signature_states": int(
                    partition["left_projected_state_product"]
                    + partition["right_projected_state_product"]
                ),
                "estimated_four_prime_signature_bytes": int(
                    (
                        partition["left_projected_state_product"]
                        + partition["right_projected_state_product"]
                    )
                    * joint_dimension
                ),
                "budget_eligible": eligible,
                "budget_skip_reason": (
                    None
                    if eligible
                    else "balanced_four_prime_unique_side_product_exceeds_max_side_states"
                ),
            }
        )
    require(len(plans) == EXPECTED_PRIOR_SURVIVORS, "four-prime candidate plan census changed")
    require(len({plan["case_key"] for plan in plans}) == len(plans), "candidate plan repeated")
    require(sum(not plan["enumerated_directions"] for plan in plans) == EXPECTED_NO_ENUMERABLE, "no-enumerable plan count changed")
    if max_side_states == DEFAULT_MAX_SIDE_STATES:
        require(all(plan["budget_eligible"] for plan in plans), "default six-million budget acquired a new skip")
    audit = {
        "candidate_cases": len(plans),
        "no_enumerable_cases": sum(not plan["enumerated_directions"] for plan in plans),
        "default_max_side_states": DEFAULT_MAX_SIDE_STATES,
        "configured_max_side_states": max_side_states,
        "prior_undeduplicated_balanced_side_histogram": {
            str(key): value for key, value in sorted(raw_histogram.items())
        },
        "four_prime_deduplicated_balanced_side_histogram": {
            str(key): value for key, value in sorted(four_histogram.items())
        },
        "largest_prior_undeduplicated_balanced_side": max(raw_histogram),
        "largest_four_prime_deduplicated_balanced_side": max(four_histogram),
        "budget_eligible_candidates": sum(plan["budget_eligible"] for plan in plans),
        "new_budget_skips": sum(not plan["budget_eligible"] for plan in plans),
        "all_four_prime_balanced_sides_bounded_by_prior_undeduplicated_balanced_sides": True,
        "all_candidates_fit_default_six_million_even_if_four_prime_deduplication_vanishes": True,
        "mod3_mod7_to_four_prime_contribution_size_pairs": {
            f"{old}_to_{new}": count
            for (old, new), count in sorted(compression_pairs.items())
        },
        "all_four_prime_contribution_counts_are_between_prior_unique_and_complete_counts": True,
        "all_four_prime_unique_contribution_counts_equal_prior_mod3_mod7_counts": all(
            old == new for old, new in compression_pairs
        ),
        "projected_base_nonzero_candidate_counts_by_modulus": {
            str(modulus): projected_base_nonzero_cases[modulus]
            for modulus in MODULI
        },
        "new_prime_active_candidate_count": len(new_prime_active_case_keys),
        "new_prime_active_case_keys": new_prime_active_case_keys,
        "new_prime_active_case_keys_sha256": json_sha256(new_prime_active_case_keys),
        "all_mod5_mod11_projected_bases_and_complete_catalog_contributions_zero_on_candidate_universe": (
            not new_prime_active_case_keys
        ),
        "four_prime_quotient_is_strictly_stronger_than_mod3_mod7_on_this_candidate_universe": bool(
            new_prime_active_case_keys
        ),
    }
    return plans, audit


def component_layout(dimensions: dict[int, int]) -> dict[int, tuple[int, int]]:
    require(set(dimensions) == set(MODULI), "component dimension set changed")
    result = {}
    offset = 0
    for modulus in MODULI:
        dimension = int(dimensions[modulus])
        require(dimension > 0, f"mod-{modulus} conditioned dimension vanished")
        result[modulus] = (offset, offset + dimension)
        offset += dimension
    return result


def enumerate_signature_set(
    matrices: tuple[np.ndarray, ...],
    component_dimensions: dict[int, int],
    state_limit: int,
    chunk_states: int,
) -> dict:
    """Enumerate and exactly deduplicate one four-prime MITM side."""
    layout = component_layout(component_dimensions)
    width = sum(component_dimensions.values())
    for matrix in matrices:
        require(matrix.ndim == 2 and matrix.shape[1] == width, "side signature width changed")
        require(
            len(np.unique(global37.raw_row_keys(matrix))) == len(matrix),
            "side contribution was not deduplicated before enumeration",
        )
        for modulus, (start, stop) in layout.items():
            require(np.all(matrix[:, start:stop] < modulus), f"mod-{modulus} residue byte escaped")
    sizes = tuple(len(matrix) for matrix in matrices)
    raw_count = math.prod(sizes) if sizes else 1
    require(raw_count <= state_limit, "side enumeration exceeded audited state limit")
    require(chunk_states > 0, "chunk state count must be positive")
    strides = global37.mixed_radix_strides(sizes)
    states = np.empty((raw_count, width), dtype=np.uint8)
    for start_index in range(0, raw_count, chunk_states):
        stop_index = min(raw_count, start_index + chunk_states)
        flat = np.arange(start_index, stop_index, dtype=np.int64)
        chunk = np.zeros((stop_index - start_index, width), dtype=np.uint8)
        for matrix, size, stride in zip(matrices, sizes, strides):
            indices = (flat // stride) % size
            for modulus, (start, stop) in layout.items():
                np.add(chunk[:, start:stop], matrix[indices, start:stop], out=chunk[:, start:stop])
                np.remainder(chunk[:, start:stop], modulus, out=chunk[:, start:stop])
        states[start_index:stop_index] = chunk
    raw_hash = affine.matrix_sha256(states)
    keys = global37.raw_row_keys(states)
    unique_keys, first_indices = np.unique(keys, return_index=True)
    unique_keys = np.ascontiguousarray(unique_keys)
    unique_rows = unique_keys.view(np.uint8).reshape(-1, width)
    first_indices = np.ascontiguousarray(first_indices, dtype=np.uint64)
    require(np.array_equal(global37.raw_row_keys(unique_rows), unique_keys), "signature key round trip failed")
    return {
        "rows": unique_rows,
        "keys": unique_keys,
        "first_flat_indices": first_indices,
        "contribution_unique_sizes": sizes,
        "raw_cartesian_states": raw_count,
        "unique_signature_states": len(unique_rows),
        "discarded_duplicate_states": raw_count - len(unique_rows),
        "signature_width_bytes": width,
        "component_dimensions": {str(p): component_dimensions[p] for p in MODULI},
        "raw_cartesian_states_sha256_uint8": raw_hash,
        "sorted_unique_signatures_sha256_uint8": affine.matrix_sha256(unique_rows),
        "representative_flat_indices_sha256_uint64": affine.matrix_sha256(
            first_indices[None, :]
        ),
        "exact_raw_byte_deduplication": True,
        "multiplicity_retained": False,
    }


def public_signature_set(state_set: dict) -> dict:
    return {
        key: value
        for key, value in state_set.items()
        if key not in {"rows", "keys", "first_flat_indices"}
    }


def meet_signature_sets(
    bases: dict[int, np.ndarray],
    left: dict,
    right: dict,
    component_dimensions: dict[int, int],
    chunk_states: int,
) -> dict:
    layout = component_layout(component_dimensions)
    width = sum(component_dimensions.values())
    left_rows = left["rows"]
    right_rows = right["rows"]
    left_keys = left["keys"]
    require(left_rows.shape[1] == right_rows.shape[1] == width, "join width changed")
    require(np.array_equal(global37.raw_row_keys(left_rows), left_keys), "left key order changed")
    match_count = 0
    first_match: tuple[int, int] | None = None
    certificate = hashlib.sha256()
    for start_index in range(0, len(right_rows), chunk_states):
        stop_index = min(len(right_rows), start_index + chunk_states)
        rows = right_rows[start_index:stop_index]
        needed_components = []
        for modulus, (start, stop) in layout.items():
            needed_components.append(
                (-bases[modulus][None, :].astype(np.int16) - rows[:, start:stop].astype(np.int16))
                % modulus
            )
        needed = np.ascontiguousarray(np.concatenate(needed_components, axis=1), dtype=np.uint8)
        needed_keys = global37.raw_row_keys(needed)
        positions = np.searchsorted(left_keys, needed_keys)
        candidates = np.flatnonzero(positions < len(left_keys))
        if not len(candidates):
            continue
        equal = left_keys[positions[candidates]] == needed_keys[candidates]
        hits = candidates[equal]
        if not len(hits):
            continue
        left_indices = positions[hits].astype(np.uint64)
        right_indices = (start_index + hits).astype(np.uint64)
        match_count += len(hits)
        certificate.update(left_indices.astype("<u8", copy=False).tobytes())
        certificate.update(right_indices.astype("<u8", copy=False).tobytes())
        if first_match is None:
            first_match = (int(left_indices[0]), int(right_indices[0]))
    base_row = np.concatenate(tuple(bases[p] for p in MODULI)).astype(np.uint8)[None, :]
    return {
        "matching_unique_signature_pairs": match_count,
        "first_matching_unique_indices": list(first_match) if first_match is not None else None,
        "matching_pair_index_certificate_sha256": certificate.hexdigest(),
        "concatenated_base_sha256_uint8": affine.matrix_sha256(base_row),
        "component_dimensions": {str(p): component_dimensions[p] for p in MODULI},
        "joint_signature_dimension": width,
        "exact_raw_byte_intersection": True,
        "hash_collision_assumption_used": False,
    }


def verify_recovered_witness(
    plan: dict,
    bases: dict[int, np.ndarray],
    left: dict,
    right: dict,
    join: dict,
    kernel_rows: np.ndarray,
) -> dict | None:
    matched = join["first_matching_unique_indices"]
    if matched is None:
        return None
    left_flat = int(left["first_flat_indices"][matched[0]])
    right_flat = int(right["first_flat_indices"][matched[1]])
    partition = plan["partition"]
    left_directions = tuple(partition["left_directions"])
    right_directions = tuple(partition["right_directions"])
    left_sizes = tuple(partition["left_contribution_unique_sizes"])
    right_sizes = tuple(partition["right_contribution_unique_sizes"])
    selected = dict(
        zip(left_directions, global37.decode_flat_index(left_flat, left_sizes))
    )
    selected.update(
        zip(right_directions, global37.decode_flat_index(right_flat, right_sizes))
    )
    require(set(selected) == set(plan["enumerated_directions"]), "witness missed a direction")
    layout = component_layout(plan["conditioned_dimensions"])
    syndromes = {p: bases[p].astype(np.int64).copy() for p in MODULI}
    exact_rows = []
    direction_rows = []
    for direction in plan["enumerated_directions"]:
        contribution = plan["projected_contributions"][direction]
        unique_index = int(selected[direction])
        require(0 <= unique_index < len(contribution["rows"]), "witness unique index escaped")
        signature = contribution["rows"][unique_index]
        catalog_index = int(contribution["first_catalog_indices"][unique_index])
        for modulus, (start, stop) in layout.items():
            syndromes[modulus] += signature[start:stop]
        mask = int(plan["orbit"]["masks"][direction])
        mean = int(plan["leaf"]["scaled_means"][direction])
        catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
        require(0 <= catalog_index < len(catalog), "witness exact catalog index escaped")
        exact_row = catalog[catalog_index]
        require(2 * int(exact_row.sum()) == 5 * mean, "witness exact mean changed")
        require(not np.any(kernel_rows @ exact_row), "witness row left exact degree two")
        require(np.array_equal(exact_row % 2, affine.parity_for_mask(mask)), "witness parity changed")
        exact_rows.append(exact_row)
        direction_rows.append(
            {
                "direction": direction,
                "catalog_class": plan["leaf"]["catalog_classes"][direction],
                "catalog_size": len(catalog),
                "four_prime_unique_contribution_index": unique_index,
                "representative_exact_catalog_row_index": catalog_index,
                "mask": mask,
                "scaled_mean": mean,
                "selected_exact_row_sha256_int64": affine.matrix_sha256(exact_row[None, :]),
            }
        )
    for modulus in MODULI:
        require(not np.any(syndromes[modulus] % modulus), f"recovered witness fails mod {modulus}")
    exact_matrix = np.stack(exact_rows) if exact_rows else np.empty((0, 35), dtype=np.int64)
    witness = {
        "left_representative_flat_index": left_flat,
        "right_representative_flat_index": right_flat,
        "catalog_rows_by_direction": direction_rows,
        "selected_exact_catalog_rows_sha256_int64": affine.matrix_sha256(exact_matrix),
        "same_exact_catalog_row_index_used_for_each_direction_across_3_5_7_11": True,
        "all_selected_rows_have_exact_leaf_means_parity_and_degree_two": True,
        "recovered_joint_syndrome_is_zero_at_all_four_primes": True,
    }
    witness["witness_certificate_sha256"] = json_sha256(witness)
    return witness


def execute_case(
    plan: dict,
    max_side_states: int,
    execution_side_cap: int,
    chunk_states: int,
    kernel_rows: np.ndarray,
    smoke_test: bool,
) -> dict:
    base_record = {
        "case_key": plan["case_key"],
        "branch_orbit_index": 0,
        "source_orbit_index": plan["source_orbit_index"],
        "orbit_leaf_index": plan["orbit_leaf_index"],
        "pointed_star_branch": "A",
        "catalog_pattern": plan["catalog_pattern"],
        "prior_global_join_decision_status": PRIOR_SURVIVOR_STATUS,
        "prior_global_join_decision_certificate_sha256": plan["prior_row"][
            "decision_certificate_sha256"
        ],
        "scaled_means": list(plan["leaf"]["scaled_means"]),
        "catalog_classes": list(plan["leaf"]["catalog_classes"]),
        "high_directions_relaxed_to_exact_zero_mean_hulls": list(plan["high_directions"]),
        "enumerated_directions_joined_globally": list(plan["enumerated_directions"]),
        "complete_catalog_sizes": list(plan["complete_catalog_sizes"]),
        "prior_mod3_mod7_unique_contribution_sizes": list(
            plan["prior_mod3_mod7_unique_contribution_sizes"]
        ),
        "four_prime_unique_contribution_sizes": list(
            plan["four_prime_unique_contribution_sizes"]
        ),
        "four_prime_balanced_partition": plan["partition"],
        "complete_catalog_balanced_partition": plan[
            "complete_catalog_balanced_partition"
        ],
        "conditioned_dimensions": {
            str(p): plan["conditioned_dimensions"][p] for p in MODULI
        },
        "joint_signature_dimension": plan["joint_signature_dimension"],
        "max_side_states": max_side_states,
    }
    if not plan["budget_eligible"]:
        result = {
            **base_record,
            "decision_status": "skipped_four_prime_side_state_budget",
            "exact_four_prime_global_join_processed": False,
            "rigorously_rejected": False,
            "necessary_only_survivor": False,
            "skipped": True,
            "skip_reason": plan["budget_skip_reason"],
        }
        result["decision_certificate_sha256"] = json_sha256(result)
        return result
    if plan["partition"]["maximum_projected_side_product"] > execution_side_cap:
        require(smoke_test, "complete run tried to apply a smoke execution cap")
        result = {
            **base_record,
            "decision_status": "skipped_four_prime_smoke_execution_cap",
            "exact_four_prime_global_join_processed": False,
            "rigorously_rejected": False,
            "necessary_only_survivor": False,
            "skipped": True,
            "skip_reason": "balanced_four_prime_unique_side_product_exceeds_smoke_execution_cap",
            "smoke_execution_side_cap": execution_side_cap,
        }
        result["decision_certificate_sha256"] = json_sha256(result)
        return result

    anchor_rhs, raw_syndromes = plan["factory"].anchor_rhs_and_raw_syndromes(
        plan["orbit"], plan["leaf"]
    )
    anchor_hash = affine.matrix_sha256(anchor_rhs[None, :])
    require(anchor_hash == plan["parent_row"]["anchor_rhs_sha256_int64"], "parent anchor RHS changed")
    require(anchor_hash == plan["prior_row"]["anchor_rhs_sha256_int64"], "global input anchor RHS changed")
    bases = plan["factory"].projected_bases(raw_syndromes, plan["high_directions"])
    require(
        {p: len(bases[p]) for p in MODULI} == plan["conditioned_dimensions"],
        "conditioned base dimensions changed",
    )
    signature_audits = []
    for direction in plan["enumerated_directions"]:
        contribution = plan["projected_contributions"][direction]
        metadata = contribution["metadata"]
        mask = int(plan["orbit"]["masks"][direction])
        mean = int(plan["leaf"]["scaled_means"][direction])
        catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
        require(np.all(2 * catalog.sum(axis=1) == 5 * mean), "complete catalog exact mean changed")
        require(not np.any(kernel_rows @ catalog.T), "complete catalog left exact degree two")
        signature_audits.append(
            {
                "direction": direction,
                "catalog_class": plan["leaf"]["catalog_classes"][direction],
                "complete_catalog_rows": len(catalog),
                "four_prime_projected_unique_contribution_states": len(
                    contribution["rows"]
                ),
                "identifier": metadata["identifier"],
                "exact_catalog_sha256_int16": metadata["exact_catalog_sha256_int16"],
                "component_signature_sha256_uint8": metadata[
                    "component_signature_sha256_uint8"
                ],
                "four_prime_complete_signature_sha256_uint8": metadata[
                    "four_prime_concatenated_signature_sha256_uint8"
                ],
                "four_prime_sorted_unique_signature_sha256_uint8": metadata[
                    "sorted_unique_four_prime_signature_sha256_uint8"
                ],
                "same_exact_catalog_row_order_used_across_all_four_primes": True,
                "all_catalog_rows_have_exact_target_mean_parity_and_degree_two": True,
            }
        )

    left_matrices = tuple(
        plan["projected_contributions"][direction]["rows"]
        for direction in plan["partition"]["left_directions"]
    )
    right_matrices = tuple(
        plan["projected_contributions"][direction]["rows"]
        for direction in plan["partition"]["right_directions"]
    )
    left = enumerate_signature_set(
        left_matrices,
        plan["conditioned_dimensions"],
        execution_side_cap,
        chunk_states,
    )
    right = enumerate_signature_set(
        right_matrices,
        plan["conditioned_dimensions"],
        execution_side_cap,
        chunk_states,
    )
    require(
        left["raw_cartesian_states"] == plan["partition"]["left_projected_state_product"]
        and right["raw_cartesian_states"] == plan["partition"]["right_projected_state_product"],
        "four-prime side state count changed",
    )
    join = meet_signature_sets(
        bases,
        left,
        right,
        plan["conditioned_dimensions"],
        chunk_states,
    )
    rejected = join["matching_unique_signature_pairs"] == 0
    witness = verify_recovered_witness(
        plan,
        bases,
        left,
        right,
        join,
        kernel_rows,
    )
    require((witness is None) == rejected, "four-prime join/witness status mismatch")
    result = {
        **base_record,
        "decision_status": (
            "rigorous_four_prime_global_catalog_join_rejection"
            if rejected
            else "necessary_only_four_prime_global_catalog_join_survivor"
        ),
        "exact_four_prime_global_join_processed": True,
        "rigorously_rejected": rejected,
        "necessary_only_survivor": not rejected,
        "skipped": False,
        "skip_reason": None,
        "anchor_rhs_sha256_int64": anchor_hash,
        "retained_catalog_signature_audits": signature_audits,
        "left_signature_set": public_signature_set(left),
        "right_signature_set": public_signature_set(right),
        "join": join,
        "recovered_same_index_catalog_witness": witness,
        "high_directions_are_still_affine_hull_relaxations": True,
        "passing_is_binary_edge_feasibility": False,
    }
    result["decision_certificate_sha256"] = json_sha256(result)
    return result


def four_prime_join_self_audit() -> dict:
    dimensions = {p: 1 for p in MODULI}
    zeros = np.zeros((1, 4), dtype=np.uint8)
    choices = np.asarray([[0, 0, 0, 0], [1, 1, 1, 1]], dtype=np.uint8)
    bases = {p: np.asarray([0], dtype=np.uint8) for p in MODULI}
    left = enumerate_signature_set((choices,), dimensions, 10, 10)
    right = enumerate_signature_set((zeros,), dimensions, 10, 10)
    positive = meet_signature_sets(bases, left, right, dimensions, 10)
    require(positive["matching_unique_signature_pairs"] == 1, "manufactured four-prime join failed")

    impossible_bases = {p: np.asarray([1], dtype=np.uint8) for p in MODULI}
    empty_left = enumerate_signature_set((), dimensions, 10, 10)
    empty_right = enumerate_signature_set((), dimensions, 10, 10)
    impossible = meet_signature_sets(
        impossible_bases, empty_left, empty_right, dimensions, 10
    )
    require(impossible["matching_unique_signature_pairs"] == 0, "impossible join passed")

    trap = np.asarray(
        [
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ],
        dtype=np.uint8,
    )
    require(all(np.any(trap[:, index] == 0) for index in range(4)), "same-index trap lacks marginal hits")
    trap_left = enumerate_signature_set((trap,), dimensions, 10, 10)
    trap_join = meet_signature_sets(bases, trap_left, empty_right, dimensions, 10)
    require(trap_join["matching_unique_signature_pairs"] == 0, "same-index trap false-positive")

    duplicates = np.asarray([[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1]], dtype=np.uint8)
    unique, first = global37.deduplicate_contribution_rows(duplicates)
    require(len(unique) == 2 and list(first) == [0, 2], "contribution dedup audit failed")
    return {
        "passed": True,
        "manufactured_positive_four_prime_join": positive,
        "manufactured_impossible_four_prime_join": impossible,
        "same_index_four_prime_false_positive_trap": trap_join,
        "each_prime_has_a_marginal_trap_hit_but_no_shared_row_hits_all_primes": True,
        "contribution_exact_deduplication_test": {
            "input_rows": len(duplicates),
            "unique_rows": len(unique),
            "first_indices": [int(value) for value in first],
        },
    }


def smoke_selection(plans: list[dict]) -> tuple[list[dict], dict]:
    by_pattern: dict[str, list[dict]] = {}
    for plan in plans:
        by_pattern.setdefault(plan["catalog_pattern"], []).append(plan)
    for rows in by_pattern.values():
        rows.sort(key=lambda row: row["case_key"])
    requested = (
        "H2_S0_M0",
        "H3_S0_M1",
        "H2_S3_M1",
        "H2_S6_M0",
        "H2_S2_M3",
    )
    selected = [by_pattern[pattern][0] for pattern in requested]
    require(len({row["case_key"] for row in selected}) == len(selected), "smoke selection repeated")
    require(all(row["prior_row"]["decision_status"] == PRIOR_SURVIVOR_STATUS for row in selected), "smoke selected a non-survivor")
    return selected, {
        "selected_case_count": len(selected),
        "requested_catalog_patterns": list(requested),
        "selected_case_keys": [row["case_key"] for row in selected],
        "selected_case_keys_sha256": json_sha256([row["case_key"] for row in selected]),
        "includes_no_enumerable_affine_only_case": any(
            not row["enumerated_directions"] for row in selected
        ),
        "all_selected_cases_are_prior_exact_mod3_mod7_global_join_survivors": True,
        "includes_deliberate_smoke_execution_cap_case": any(
            row["partition"]["maximum_projected_side_product"]
            > DEFAULT_SMOKE_MAX_SIDE_STATES
            for row in selected
        ),
        "full_159_candidate_coverage_claimed": False,
    }


def compact_inherited_result(row: dict) -> dict:
    status = str(row["decision_status"])
    if status == PRIOR_REJECTION_STATUS:
        decision = "preserved_prior_rigorous_mod3_mod7_global_join_rejection"
        rigorous = True
        skipped = False
    elif status == PRIOR_SKIP_STATUS:
        decision = "preserved_prior_mod3_mod7_side_state_budget_skip"
        rigorous = False
        skipped = True
    else:
        raise AssertionError("attempted to inherit a prior survivor without four-prime processing")
    result = {
        "case_key": row["case_key"],
        "catalog_pattern": row["catalog_pattern"],
        "prior_global_join_decision_status": status,
        "prior_global_join_decision_certificate_sha256": row[
            "decision_certificate_sha256"
        ],
        "decision_status": decision,
        "four_prime_processing_attempted": False,
        "rigorously_rejected": rigorous,
        "necessary_only_survivor": False,
        "skipped": skipped,
        "semantics": (
            "A rigorous mod-3/mod-7 rejection remains rigorous after adding moduli."
            if rigorous
            else "This representative remains unprocessed because the completed predecessor skipped it."
        ),
    }
    result["decision_certificate_sha256"] = json_sha256(result)
    return result


def unselected_smoke_result(row: dict) -> dict:
    result = {
        "case_key": row["case_key"],
        "catalog_pattern": row["catalog_pattern"],
        "prior_global_join_decision_status": PRIOR_SURVIVOR_STATUS,
        "prior_global_join_decision_certificate_sha256": row[
            "decision_certificate_sha256"
        ],
        "decision_status": "not_selected_by_bounded_four_prime_smoke",
        "four_prime_processing_attempted": False,
        "rigorously_rejected": False,
        "necessary_only_survivor": False,
        "skipped": True,
        "semantics": "Bounded smoke mode makes no four-prime decision for this candidate.",
    }
    result["decision_certificate_sha256"] = json_sha256(result)
    return result


def compose_all_representative_results(
    prior_rows: list[dict],
    candidate_results: list[dict],
    smoke_test: bool,
) -> list[dict]:
    new_by_key = {str(row["case_key"]): row for row in candidate_results}
    require(len(new_by_key) == len(candidate_results), "new candidate result repeated")
    rows = []
    for prior_row in prior_rows:
        status = str(prior_row["decision_status"])
        key = str(prior_row["case_key"])
        if status in (PRIOR_REJECTION_STATUS, PRIOR_SKIP_STATUS):
            rows.append(compact_inherited_result(prior_row))
        elif key in new_by_key:
            rows.append(new_by_key[key])
        else:
            require(smoke_test, "complete run omitted a prior survivor candidate")
            rows.append(unselected_smoke_result(prior_row))
    require(len(rows) == EXPECTED_REPRESENTATIVES, "composed representative census changed")
    return rows


def summarize_results(rows: list[dict]) -> dict:
    counts = Counter(str(row["decision_status"]) for row in rows)
    by_pattern: dict[str, Counter[str]] = {}
    for row in rows:
        by_pattern.setdefault(str(row["catalog_pattern"]), Counter())[str(row["decision_status"])] += 1
    return {
        "decision_status_histogram": dict(sorted(counts.items())),
        "coverage_by_catalog_pattern": {
            pattern: dict(sorted(counter.items()))
            for pattern, counter in sorted(by_pattern.items())
        },
        "preserved_prior_rigorous_rejections": counts[
            "preserved_prior_rigorous_mod3_mod7_global_join_rejection"
        ],
        "preserved_prior_budget_skips": counts[
            "preserved_prior_mod3_mod7_side_state_budget_skip"
        ],
        "new_four_prime_rigorous_rejections": counts[
            "rigorous_four_prime_global_catalog_join_rejection"
        ],
        "four_prime_necessary_only_survivors": counts[
            "necessary_only_four_prime_global_catalog_join_survivor"
        ],
        "new_four_prime_budget_skips": counts[
            "skipped_four_prime_side_state_budget"
        ],
        "smoke_execution_cap_skips": counts[
            "skipped_four_prime_smoke_execution_cap"
        ],
        "not_selected_by_smoke": counts[
            "not_selected_by_bounded_four_prime_smoke"
        ],
    }


def checkpoint_payload(
    run_identity: str,
    input_provenance: dict,
    configuration: dict,
    selected_keys: list[str],
    rows: list[dict],
) -> dict:
    return {
        "experiment": EXPERIMENT,
        "status": "in_progress_atomic_checkpoint",
        "run_identity_sha256": run_identity,
        "input_provenance": input_provenance,
        "configuration": configuration,
        "selected_candidate_case_keys": selected_keys,
        "selected_candidate_case_keys_sha256": json_sha256(selected_keys),
        "selected_candidate_case_count": len(selected_keys),
        "completed_candidate_case_count": len(rows),
        "next_candidate_case_index": len(rows),
        "completed_candidate_case_keys_sha256": json_sha256(
            [row["case_key"] for row in rows]
        ),
        "candidate_case_results_sha256": canonical_case_digest(rows),
        "candidate_case_results": rows,
        "checkpoint_is_resume_safe_and_atomically_replaced": True,
    }


def load_resume_rows(
    output: Path,
    run_identity: str,
    selected_keys: list[str],
) -> tuple[list[dict], dict | None]:
    if not output.exists():
        return [], None
    payload = json.loads(output.read_text(encoding="utf-8"))
    require(payload["experiment"] == EXPERIMENT, "resume output belongs to another experiment")
    require(payload["run_identity_sha256"] == run_identity, "resume inputs/configuration changed")
    if payload["status"] != "in_progress_atomic_checkpoint":
        require(payload["case_results_sha256"] == canonical_case_digest(payload["case_results"]), "completed output digest failed")
        return [], payload
    rows = payload["candidate_case_results"]
    require(payload["candidate_case_results_sha256"] == canonical_case_digest(rows), "checkpoint digest failed")
    require([row["case_key"] for row in rows] == selected_keys[: len(rows)], "checkpoint is not a selected prefix")
    require(payload["next_candidate_case_index"] == len(rows), "checkpoint next index changed")
    for row in rows:
        certificate = row["decision_certificate_sha256"]
        body = {key: value for key, value in row.items() if key != "decision_certificate_sha256"}
        require(certificate == json_sha256(body), "checkpoint decision certificate failed")
    return rows, None


def run(
    parent_input: Path,
    global_input: Path,
    symmetry_certificate: Path,
    output: Path,
    max_side_states: int = DEFAULT_MAX_SIDE_STATES,
    chunk_states: int = DEFAULT_CHUNK_STATES,
    checkpoint_every: int = DEFAULT_CHECKPOINT_EVERY,
    resume: bool = False,
    smoke_test: bool = False,
    smoke_max_side_states: int = DEFAULT_SMOKE_MAX_SIDE_STATES,
    quiet: bool = False,
) -> dict:
    started = time.time()
    require(max_side_states > 0, "max side states must be positive")
    require(chunk_states > 0, "chunk states must be positive")
    require(checkpoint_every > 0, "checkpoint interval must be positive")
    require(smoke_max_side_states > 0, "smoke side cap must be positive")

    parent_payload, parent_survivors, parent_provenance = global37.load_parent_input(
        parent_input
    )
    parent_representatives = global37.representative_survivors(parent_survivors)
    global_payload, prior_rows, candidate_rows, global_provenance, global_audit = (
        load_completed_global_input(
            global_input,
            parent_payload,
            parent_representatives,
            parent_provenance,
        )
    )
    symmetry_raw = symmetry_certificate.read_bytes()
    symmetry_payload = json.loads(symmetry_raw)

    if not quiet:
        print("reconstructing exact parent geometry and four pointed systems", flush=True)
    rebuilt = global37.reconstruction()
    reconstruction_audit = global37.validate_reconstruction_against_parent(
        parent_payload, rebuilt
    )
    hull_bases, hull_audit = build_exact_zero_mean_hulls(rebuilt)
    pointed_dependency_audit = augment_complete_pointed_dependencies(
        rebuilt, symmetry_payload
    )
    symmetry_audit, symmetry_mapping = validate_symmetry_certificate(
        symmetry_certificate,
        symmetry_payload,
        parent_survivors,
        parent_representatives,
        parent_provenance,
        hull_audit,
    )

    representative_system = rebuilt["systems"][0]["A"]
    factory = FourPrimeProjectionFactory(
        representative_system, hull_bases, rebuilt["anchors"]
    )
    contribution_store = FourPrimeContributionStore(factory)
    parent_by_key = {str(row["case_key"]): row for row in parent_representatives}
    if not quiet:
        print("building actual four-prime contribution quotients for 159 candidates", flush=True)
    plans, budget_audit = make_candidate_plans(
        candidate_rows,
        parent_by_key,
        rebuilt,
        factory,
        contribution_store,
        max_side_states,
    )
    plans_by_key = {plan["case_key"]: plan for plan in plans}
    require(
        [plan["case_key"] for plan in plans]
        == [str(row["case_key"]) for row in candidate_rows],
        "candidate plan order changed",
    )
    if smoke_test:
        selected_plans, smoke_audit = smoke_selection(plans)
        execution_side_cap = min(max_side_states, smoke_max_side_states)
    else:
        selected_plans = plans
        smoke_audit = None
        execution_side_cap = max_side_states
    selected_keys = [plan["case_key"] for plan in selected_plans]

    input_provenance = {
        "parent_affine_hull": parent_provenance,
        "completed_mod3_mod7_global_join": global_provenance,
        "explicit_four_prime_symmetry_certificate": {
            "path": str(symmetry_certificate.resolve()),
            "file_bytes": len(symmetry_raw),
            "file_sha256": hashlib.sha256(symmetry_raw).hexdigest(),
            "experiment": symmetry_payload.get("experiment"),
            "status": symmetry_payload.get("status"),
        },
    }
    configuration = {
        "max_side_states": max_side_states,
        "execution_side_cap": execution_side_cap,
        "chunk_states": chunk_states,
        "checkpoint_every": checkpoint_every,
        "smoke_test": smoke_test,
        "smoke_max_side_states": smoke_max_side_states if smoke_test else None,
        "candidate_universe_count": len(plans),
        "selected_candidate_count": len(selected_plans),
    }
    run_identity = json_sha256(
        {
            "experiment": EXPERIMENT,
            "input_file_hashes": {
                key: value["file_sha256"] for key, value in input_provenance.items()
            },
            "configuration": configuration,
            "selected_keys": selected_keys,
            "script_sha256": file_sha256(Path(__file__).resolve()),
        }
    )
    completed_rows: list[dict] = []
    resumed_candidate_count = 0
    completed_output = None
    if resume:
        completed_rows, completed_output = load_resume_rows(
            output, run_identity, selected_keys
        )
        resumed_candidate_count = len(completed_rows)
        if completed_output is not None:
            return completed_output
    elif output.exists():
        # The final atomic write is authorized, but silently consuming stale
        # checkpoints without --resume would make provenance ambiguous.
        existing = json.loads(output.read_text(encoding="utf-8"))
        require(
            existing.get("status") != "in_progress_atomic_checkpoint",
            "output contains a checkpoint; pass --resume or choose another output",
        )

    self_audit = four_prime_join_self_audit()
    for index, plan in enumerate(selected_plans[len(completed_rows) :], start=len(completed_rows)):
        if not quiet:
            print(
                f"[{index + 1}/{len(selected_plans)}] {plan['case_key']} "
                f"{plan['catalog_pattern']} side={plan['partition']['maximum_projected_side_product']}",
                flush=True,
            )
        result = execute_case(
            plan,
            max_side_states,
            execution_side_cap,
            chunk_states,
            rebuilt["kernel_rows"],
            smoke_test,
        )
        completed_rows.append(result)
        if len(completed_rows) % checkpoint_every == 0 and len(completed_rows) < len(selected_plans):
            affine.pointed.atomic_write(
                output,
                checkpoint_payload(
                    run_identity,
                    input_provenance,
                    configuration,
                    selected_keys,
                    completed_rows,
                ),
            )

    require(len(completed_rows) == len(selected_plans), "selected candidate execution incomplete")
    all_results = compose_all_representative_results(
        prior_rows, completed_rows, smoke_test
    )
    summary = summarize_results(all_results)
    require(
        summary["preserved_prior_rigorous_rejections"] == EXPECTED_PRIOR_REJECTIONS,
        "prior rejections were not preserved",
    )
    require(
        summary["preserved_prior_budget_skips"] == EXPECTED_PRIOR_SKIPS,
        "prior skips were not preserved",
    )
    complete_candidate_run = not smoke_test and len(completed_rows) == EXPECTED_PRIOR_SURVIVORS
    symmetry_transfer_claimed = complete_candidate_run
    transferred_counts = None
    if symmetry_transfer_claimed:
        require(len(symmetry_mapping) == EXPECTED_REPRESENTATIVES, "symmetry mapping census changed")
        transferred_counts = {
            "full_parent_survivor_cases": 1_296,
            "rigorous_rejections": 4
            * (
                summary["preserved_prior_rigorous_rejections"]
                + summary["new_four_prime_rigorous_rejections"]
            ),
            "necessary_only_four_prime_survivors": 4
            * summary["four_prime_necessary_only_survivors"],
            "preserved_prior_budget_skips": 4 * summary["preserved_prior_budget_skips"],
            "new_four_prime_budget_skips": 4 * summary["new_four_prime_budget_skips"],
        }

    output_payload = {
        "experiment": EXPERIMENT,
        "status": (
            "bounded_smoke_test_only"
            if smoke_test
            else "complete_bounded_four_prime_global_catalog_join_sieve_with_preserved_prior_decisions"
        ),
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "moduli": list(MODULI),
        "smoke_test": smoke_test,
        "run_identity_sha256": run_identity,
        "configuration": configuration,
        "input_provenance": input_provenance,
        "source_provenance": {
            "this_script_path": str(Path(__file__).resolve()),
            "this_script_sha256": file_sha256(Path(__file__).resolve()),
            "trusted_global_helper_path": str(
                (ROOT / "scripts" / "p7_infinity7_positive_z7_global_catalog_join.py").resolve()
            ),
            "trusted_global_helper_sha256": file_sha256(
                ROOT / "scripts" / "p7_infinity7_positive_z7_global_catalog_join.py"
            ),
            "trusted_affine_helper_sha256": file_sha256(
                ROOT / "scripts" / "p7_infinity7_positive_z7_pointed_affine_hull_multimod.py"
            ),
            "reused_safe_generic_helpers": [
                "load_parent_input",
                "representative_survivors",
                "reconstruction",
                "validate_reconstruction_against_parent",
                "validate_parent_survivor",
                "balanced_partition",
                "deduplicate_contribution_rows",
                "mixed_radix_strides",
                "decode_flat_index",
                "raw_row_keys",
                "atomic_write",
            ],
        },
        "completed_global_input_audit": global_audit,
        "reconstruction_audit": reconstruction_audit,
        "exact_zero_mean_degree_two_lattice_hull_audit": hull_audit,
        "complete_pointed_left_dependency_audit": pointed_dependency_audit,
        "symmetry_certificate_audit": symmetry_audit,
        "symmetry_transfer_claimed": symmetry_transfer_claimed,
        "transferred_full_1296_counts": transferred_counts,
        "four_prime_projection_factory_audit": factory.audit(),
        "four_prime_contribution_deduplication_audit": contribution_store.audit(),
        "preflight_budget_audit": budget_audit,
        "four_prime_global_join_self_audit": self_audit,
        "smoke_selection_audit": smoke_audit,
        "candidate_universe_count": len(plans),
        "selected_candidate_count": len(selected_plans),
        "completed_candidate_count": len(completed_rows),
        "candidate_case_results_sha256": canonical_case_digest(completed_rows),
        "candidate_case_results": completed_rows,
        "representative_case_count": len(all_results),
        "case_results_sha256": canonical_case_digest(all_results),
        "case_results": all_results,
        "result_summary": summary,
        "logical_semantics": {
            "zero_four_prime_join_is_a_rigorous_rejection": True,
            "positive_four_prime_join_is_necessary_only": True,
            "positive_join_is_not_binary_edge_feasibility": True,
            "all_enumerable_directions_share_one_exact_catalog_row_index_across_all_four_primes": True,
            "high_directions_are_relaxed_only_to_complete_liftable_exact_zero_mean_degree_two_hulls": True,
            "prior_rigorous_mod3_mod7_rejections_remain_rigorous": True,
            "prior_budget_skips_are_not_silently_reclassified": True,
            "smoke_mode_transfers_no_decisions": True,
        },
        "checkpoint": {
            "resume_requested": resume,
            "resumed_candidate_count": resumed_candidate_count,
            "atomic_output_path": str(output.resolve()),
            "checkpoint_every": checkpoint_every,
            "final_output_is_atomic": True,
        },
        "elapsed_seconds": time.time() - started,
    }
    return output_payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-input", type=Path, required=True, help="33 MB full affine-hull JSON")
    parser.add_argument("--global-input", type=Path, required=True, help="completed 324-representative mod-3/mod-7 global join JSON")
    parser.add_argument("--symmetry-certificate", type=Path, required=True, help="explicit four-prime compact symmetry certificate")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-side-states", type=int, default=DEFAULT_MAX_SIDE_STATES)
    parser.add_argument("--chunk-states", type=int, default=DEFAULT_CHUNK_STATES)
    parser.add_argument("--checkpoint-every", type=int, default=DEFAULT_CHECKPOINT_EVERY)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="execute only a deterministic bounded candidate subset")
    parser.add_argument("--smoke-max-side-states", type=int, default=DEFAULT_SMOKE_MAX_SIDE_STATES)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = run(
        parent_input=args.parent_input,
        global_input=args.global_input,
        symmetry_certificate=args.symmetry_certificate,
        output=args.output,
        max_side_states=args.max_side_states,
        chunk_states=args.chunk_states,
        checkpoint_every=args.checkpoint_every,
        resume=args.resume,
        smoke_test=args.smoke,
        smoke_max_side_states=args.smoke_max_side_states,
        quiet=args.quiet,
    )
    affine.pointed.atomic_write(args.output, result)
    if not args.quiet:
        print(json.dumps({
            "status": result["status"],
            "selected_candidate_count": result["selected_candidate_count"],
            "completed_candidate_count": result["completed_candidate_count"],
            "result_summary": result["result_summary"],
            "elapsed_seconds": result["elapsed_seconds"],
            "output": str(args.output.resolve()),
        }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
