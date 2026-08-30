#!/usr/bin/env python3
"""Bounded exact global catalog join for positive p=7, z=7 survivors.

The parent affine-hull sieve tests every retained pair separately, so distinct
pair tests may be witnessed by mutually incompatible catalog rows.  This
script removes that relaxation for every enumerable direction at once.

For each parent survivor, all high directions are still relaxed to their
audited exact-zero-mean degree-two affine hulls.  In the resulting quotient,
one complete catalog row is chosen simultaneously for every S/M direction.
Each row's mod-3 and mod-7 components come from the same retained catalog
index and are concatenated before any sum is formed.  A balanced exact
meet-in-the-middle join then asks whether the global sum equals the exact
anchor target.  Each direction's projected catalog is first deduplicated
exactly; signature multiplicity is irrelevant to existence.  The side-state
budget and balanced partition therefore use products of projected unique
contribution states, not products of raw catalog cardinalities.  Side sums
are deduplicated again.  Raw-byte signatures themselves are retained
throughout, so there is no probabilistic hashing or collision assumption.

The join is deliberately bounded.  A case is processed exactly when both
projected unique-state products in its optimal balanced partition fit
``--max-side-states``.  Every other case is emitted with an explicit
deterministic skip reason.  ``--representatives-only`` selects the 324
orbit-0/branch-A compact representatives.  Results are not transferred to
the other 972 cases unless an explicit four-case symmetry certificate is
supplied and validated.  A zero intersection rigorously rejects the selected
parent survivor.  A nonzero intersection is still necessary only, because
high catalogs remain affine-hull relaxations and modular right-hand-side
consistency is not a binary edge lift.
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

import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_global_catalog_join"
PARENT_EXPERIMENT = "p7_infinity7_positive_z7_pointed_affine_hull_multimod"
EXPECTED_PARENT_SURVIVORS = 1_296
EXPECTED_PARENT_CASES = 4_320
EXPECTED_REPRESENTATIVES = 324
DEFAULT_MAX_SIDE_STATES = 6_000_000
DEFAULT_SMOKE_MAX_SIDE_STATES = 200_000
DEFAULT_CHUNK_STATES = 20_000
DEFAULT_CHECKPOINT_EVERY = 10
EXPECTED_DEFAULT_ELIGIBLE = 984
EXPECTED_DEFAULT_SKIPPED = 312
EXPECTED_DEFAULT_ENUMERABLE_ELIGIBLE_REPRESENTATIVES = 236
EXPECTED_DEFAULT_SKIPPED_REPRESENTATIVES = 78
EXPECTED_NO_ENUMERABLE_REPRESENTATIVES = 10
MODULI = (3, 7)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def raw_row_keys(rows: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(rows, dtype=np.uint8)
    require(rows.ndim == 2 and rows.shape[1] > 0, "signature rows have bad shape")
    return rows.view(np.dtype((np.void, rows.shape[1]))).reshape(-1)


def canonical_case_digest(rows: list[dict]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def load_parent_input(path: Path) -> tuple[dict, list[dict], dict]:
    raw = path.read_bytes()
    payload = json.loads(raw)
    require(payload["experiment"] == PARENT_EXPERIMENT, "input is not the affine-hull parent output")
    require(payload["full_run"] is True and payload["smoke_test"] is False, "input is not a full parent run")
    require(payload["p"] == 7 and payload["z"] == 7 and payload["phase"] == 0, "input scope changed")
    require(payload["c_H"] == 1 and payload["moduli"] == [3, 7], "input sign/moduli changed")
    require(
        payload["processed_pointed_branch_cases"] == EXPECTED_PARENT_CASES,
        "parent did not cover all pointed cases",
    )
    require(payload["z7_branch_excluded"] is False, "parent unexpectedly already excluded z=7")

    case_rows = payload["case_results"]
    require(len(case_rows) == EXPECTED_PARENT_CASES, "parent case-result census changed")
    require(
        canonical_case_digest(case_rows) == payload["all_case_results_sha256"],
        "parent case-result certificate failed",
    )
    case_keys = [str(row["case_key"]) for row in case_rows]
    require(len(case_keys) == len(set(case_keys)), "parent case keys are not unique")
    survivors = [row for row in case_rows if row["passes_all_necessary_projections"]]
    survivor_keys = [str(row["case_key"]) for row in survivors]
    require(len(survivors) == payload["surviving_pointed_branch_cases"], "parent survivor count disagrees")
    require(len(survivors) == EXPECTED_PARENT_SURVIVORS, "expected 1,296 parent survivors")
    require(survivor_keys == payload["survivor_case_keys"], "parent survivor order/key list changed")
    require(
        affine.json_sha256(survivor_keys) == payload["survivor_case_keys_sha256"],
        "parent survivor-key certificate failed",
    )
    require(
        all(
            projection["matching_joint_projected_catalog_tuples"] > 0
            for row in survivors
            for projection in row["projections"]
        ),
        "a purported parent survivor contains a failed pair projection",
    )
    return payload, survivors, {
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "experiment": payload["experiment"],
        "status": payload["status"],
        "all_case_results_sha256": payload["all_case_results_sha256"],
        "survivor_case_keys_sha256": payload["survivor_case_keys_sha256"],
        "processed_pointed_branch_cases": payload["processed_pointed_branch_cases"],
        "surviving_pointed_branch_cases": len(survivors),
        "full_parent_case_and_survivor_certificates_recomputed": True,
    }


def representative_survivors(survivors: list[dict]) -> list[dict]:
    rows = [
        row
        for row in survivors
        if int(row["branch_orbit_index"]) == 0
        and str(row["pointed_star_branch"]) == "A"
    ]
    require(len(rows) == EXPECTED_REPRESENTATIVES, "orbit0/branchA representative census changed")
    require(len({int(row["orbit_leaf_index"]) for row in rows}) == len(rows), "representative leaf repeated")
    return rows


def symmetry_certificate_contract(parent_survivor_keys_sha256: str) -> dict:
    return {
        "required_scope": "positive_p7_z7_1296_affine_hull_survivors",
        "required_parent_survivor_case_keys_sha256": parent_survivor_keys_sha256,
        "required_equivalence_class_count": EXPECTED_REPRESENTATIVES,
        "required_class_size": 4,
        "required_class_row_fields": [
            "representative_case_key",
            "member_case_keys",
        ],
        "required_top_level_or_transfer_block_booleans": [
            "all_class_maps_preserve_exact_catalog_row_identity_mod3_mod7",
            "all_class_maps_preserve_high_affine_hulls_and_anchor_targets",
            "transfer_valid_for_global_same_index_catalog_join",
        ],
        "semantics": (
            "The 324 classes must partition the exact 1,296 parent survivor keys; "
            "each representative must be the orbit0/branchA member."
        ),
    }


def validate_compact_symmetry_audit_certificate(
    path: Path,
    raw: bytes,
    payload: dict,
    survivors: list[dict],
    representatives: list[dict],
    input_provenance: dict,
    contract: dict,
) -> tuple[dict, dict[str, tuple[str, ...]]]:
    require(
        payload["status"]
        == "complete_exact_four_case_compact_symmetry_and_survivor_partition_audit",
        "compact symmetry certificate did not audit the survivor input",
    )
    require(
        payload["p"] == 7
        and payload["c_H"] == 1
        and payload["z"] == 7
        and payload["phase"] == 0,
        "compact symmetry certificate scope changed",
    )
    require(payload["all_required_structural_audits_passed"] is True, "structural symmetry audit failed")
    row_space = payload["canonical_compact_dependency_audit"]
    require(
        row_space["all_four_augmented_compact_row_spaces_identical"] is True,
        "four compact dependency row spaces are not identical",
    )
    require(
        all(
            row["all_four_canonical_row_spaces_identical"] is True
            for row in row_space["common_prime_audits"]
            if int(row["modulus"]) in MODULI
        ),
        "mod-3/mod-7 compact row-space identity failed",
    )
    isomorphism = payload["affine_compact_isomorphism_audit"]
    for field in (
        "identity_maps_A_to_B_within_each_orbit_after_compact_dependency_projection",
        "branch_A_map_connects_orbit0_A_to_orbit1_A",
        "branch_B_map_connects_orbit0_B_to_orbit1_B",
        "all_four_compact_cases_are_exactly_isomorphic",
    ):
        require(isomorphism[field] is True, f"compact symmetry field failed: {field}")
    for map_row in isomorphism["maps"]:
        for field in (
            "maps_branch_fixed_edge_normalization",
            "all_boundary_masks_and_280_parity_domains_transport",
            "all_14_exact_primitive_kernel_equations_per_direction_transport",
            "directional_sum_and_exact_mean_identity_transport",
            "all_compact_constraint_families_transport",
        ):
            require(map_row[field] is True, f"affine compact map failed: {field}")
        require(
            all(
                row["equals_target_augmented_row_space"] is True
                for row in map_row["augmented_modular_row_space_transport"]
                if int(row["modulus"]) in MODULI
            ),
            "affine map failed mod-3/mod-7 augmented row-space transport",
        )

    leaf_transport = payload["leaf_transport_audit"]
    require(leaf_transport["leaf_transport_is_bijective"] is True, "leaf transport is not bijective")
    require(
        leaf_transport["all_scaled_means_q_values_levels_and_classes_transport"] is True,
        "leaf catalog metadata does not transport",
    )
    require(leaf_transport["all_4320_pointed_leaf_cases_partitioned_once"] is True, "leaf classes are incomplete")
    require(
        leaf_transport["leaf_permutation_sha256_int64"]
        == leaf_transport["expected_leaf_permutation_sha256_int64"],
        "leaf permutation certificate changed",
    )
    survivor_audit = payload["optional_affine_survivor_partition_audit"]
    require(survivor_audit["performed"] is True, "survivor partition audit was not performed")
    require(
        survivor_audit[
            "all_1296_survivors_partition_into_exactly_324_complete_four_case_classes"
        ]
        is True,
        "survivors do not form 324 complete four-case classes",
    )
    require(
        survivor_audit["complete_four_case_classes"] == EXPECTED_REPRESENTATIVES
        and survivor_audit["partial_four_case_classes"] == 0,
        "survivor class census changed",
    )
    evidence = survivor_audit["input_evidence_audit"]
    require(evidence["file_sha256"] == input_provenance["file_sha256"], "symmetry audit used another input file")
    require(
        evidence["all_case_results_sha256"]
        == input_provenance["all_case_results_sha256"],
        "symmetry audit parent case certificate changed",
    )
    require(
        evidence["survivor_case_keys_sha256"]
        == input_provenance["survivor_case_keys_sha256"],
        "symmetry audit survivor certificate changed",
    )
    transfer_scope = payload["exact_transfer_scope"]
    require(
        transfer_scope["compact_infeasibility_transfers_across_each_four_case_class"]
        is True
        and transfer_scope["compact_feasibility_transfers_across_each_four_case_class"]
        is True,
        "compact decisions do not transfer",
    )

    sigma = tuple(int(value) for value in leaf_transport["direction_permutation_source_to_target"])
    require(sorted(sigma) == list(range(8)), "certificate direction transport is not bijective")
    target_lookup = {
        tuple(int(value) for value in row["scaled_means"]): int(row["orbit_leaf_index"])
        for row in survivors
        if int(row["branch_orbit_index"]) == 1
        and str(row["pointed_star_branch"]) == "A"
    }
    require(len(target_lookup) == EXPECTED_REPRESENTATIVES, "target representative means are not unique")
    survivor_keys = {str(row["case_key"]) for row in survivors}
    mapping: dict[str, tuple[str, ...]] = {}
    canonical_classes = []
    for row in representatives:
        source_leaf = int(row["orbit_leaf_index"])
        transported = [0] * 8
        for source_direction, target_direction in enumerate(sigma):
            transported[target_direction] = int(row["scaled_means"][source_direction])
        target_leaf = target_lookup.get(tuple(transported))
        require(target_leaf is not None, "representative has no certified orbit1 survivor image")
        representative = str(row["case_key"])
        members = (
            f"orbit0_leaf{source_leaf}_branchA",
            f"orbit0_leaf{source_leaf}_branchB",
            f"orbit1_leaf{target_leaf}_branchA",
            f"orbit1_leaf{target_leaf}_branchB",
        )
        require(set(members) <= survivor_keys, "derived symmetry class is not a complete survivor class")
        mapping[representative] = members
        canonical_classes.append(
            {
                "representative_case_key": representative,
                "member_case_keys": list(members),
            }
        )
    require(
        {key for members in mapping.values() for key in members} == survivor_keys,
        "derived symmetry classes do not cover the survivor set",
    )
    canonical_hash = json_sha256(canonical_classes)
    transfer = payload.get("global_catalog_join_transfer")
    require(isinstance(transfer, dict), "compact symmetry certificate lacks global-join transfer audit")
    require(transfer.get("scope") == contract["required_scope"], "global-join transfer scope changed")
    require(
        transfer.get("parent_survivor_case_keys_sha256")
        == input_provenance["survivor_case_keys_sha256"],
        "global-join transfer is for another survivor set",
    )
    for field in contract["required_top_level_or_transfer_block_booleans"]:
        require(transfer.get(field) is True, f"global-join transfer did not prove {field}")
    declared_classes = transfer.get("equivalence_classes")
    require(
        isinstance(declared_classes, list)
        and len(declared_classes) == EXPECTED_REPRESENTATIVES,
        "global-join transfer equivalence-class census changed",
    )
    declared_canonical = [
        {
            "representative_case_key": str(row["representative_case_key"]),
            "member_case_keys": [str(value) for value in row["member_case_keys"]],
        }
        for row in declared_classes
    ]
    require(
        declared_canonical == canonical_classes,
        "global-join transfer classes disagree with independently derived compact classes",
    )
    require(
        transfer.get("equivalence_classes_sha256") == canonical_hash,
        "global-join transfer equivalence-class certificate changed",
    )
    return {
        "provided": True,
        "validated": True,
        "transfer_claimed": True,
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "experiment": payload["experiment"],
        "status": payload["status"],
        "scope": contract["required_scope"],
        "equivalence_class_count": len(mapping),
        "class_size": 4,
        "covered_parent_survivors": len(survivor_keys),
        "leaf_permutation_sha256_int64": leaf_transport[
            "leaf_permutation_sha256_int64"
        ],
        "canonical_representative_member_classes_sha256": canonical_hash,
        "source_declared_equivalence_classes_sha256": transfer[
            "equivalence_classes_sha256"
        ],
        "global_join_transfer_proof_composition": transfer["proof_composition"],
        "global_join_transfer_inferred_from_exact_certified_coordinate_bijections": False,
        "global_join_transfer_explicitly_audited_by_source_certificate": True,
        "all_required_transfer_claims_validated": True,
        "required_contract": contract,
    }, mapping


def load_symmetry_certificate(
    path: Path | None,
    survivors: list[dict],
    representatives: list[dict],
    input_provenance: dict,
) -> tuple[dict, dict[str, tuple[str, ...]] | None]:
    parent_survivor_keys_sha256 = input_provenance["survivor_case_keys_sha256"]
    contract = symmetry_certificate_contract(parent_survivor_keys_sha256)
    if path is None:
        return {
            "provided": False,
            "validated": False,
            "transfer_claimed": False,
            "path": None,
            "required_contract": contract,
            "reason": "no symmetry certificate supplied",
        }, None

    raw = path.read_bytes()
    payload = json.loads(raw)
    if payload.get("experiment") == "p7_infinity7_positive_z7_compact_symmetry_audit":
        return validate_compact_symmetry_audit_certificate(
            path,
            raw,
            payload,
            survivors,
            representatives,
            input_provenance,
            contract,
        )
    transfer = payload.get("global_catalog_join_transfer", payload)
    require("complete" in str(payload.get("status", "")).lower(), "symmetry certificate is not complete")
    require(
        str(transfer.get("scope")) == contract["required_scope"],
        "symmetry certificate scope changed",
    )
    require(
        transfer.get("parent_survivor_case_keys_sha256")
        == parent_survivor_keys_sha256,
        "symmetry certificate is for another parent survivor set",
    )
    for field in contract["required_top_level_or_transfer_block_booleans"]:
        require(transfer.get(field) is True, f"symmetry certificate did not prove {field}")
    classes = transfer.get("equivalence_classes")
    require(isinstance(classes, list) and len(classes) == EXPECTED_REPRESENTATIVES, "symmetry class census changed")

    survivor_keys = {str(row["case_key"]) for row in survivors}
    representative_keys = {str(row["case_key"]) for row in representatives}
    observed_members: set[str] = set()
    mapping: dict[str, tuple[str, ...]] = {}
    canonical_classes = []
    for row in classes:
        representative = str(row["representative_case_key"])
        members = tuple(str(value) for value in row["member_case_keys"])
        require(representative in representative_keys, "symmetry representative is not orbit0/branchA")
        require(len(members) == 4 and len(set(members)) == 4, "symmetry class is not four distinct cases")
        require(representative in members, "symmetry class omits its representative")
        require(set(members) <= survivor_keys, "symmetry class contains a non-survivor key")
        require(not (set(members) & observed_members), "symmetry classes overlap")
        require(row.get("global_same_index_join_equivalence_proved", True) is True, "class lacks global-join equivalence")
        observed_members.update(members)
        mapping[representative] = members
        canonical_classes.append(
            {
                "representative_case_key": representative,
                "member_case_keys": list(members),
            }
        )
    require(observed_members == survivor_keys, "symmetry classes do not partition all 1,296 survivors")
    require(set(mapping) == representative_keys, "symmetry classes miss a representative")
    canonical_class_hash = json_sha256(canonical_classes)
    return {
        "provided": True,
        "validated": True,
        "transfer_claimed": True,
        "path": str(path.resolve()),
        "file_bytes": len(raw),
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "experiment": payload.get("experiment"),
        "status": payload.get("status"),
        "scope": transfer["scope"],
        "equivalence_class_count": len(classes),
        "class_size": 4,
        "covered_parent_survivors": len(observed_members),
        "source_declared_equivalence_classes_sha256": transfer.get(
            "equivalence_classes_sha256"
        ),
        "canonical_representative_member_classes_sha256": canonical_class_hash,
        "all_required_transfer_claims_validated": True,
        "required_contract": contract,
    }, mapping


def balanced_partition(
    directions: tuple[int, ...], sizes: tuple[int, ...]
) -> dict:
    """Minimize the larger projected unique-state product."""
    require(len(directions) == len(sizes), "direction/catalog-size arity changed")
    require(len(set(directions)) == len(directions), "enumerated direction repeated")
    require(all(size > 0 for size in sizes), "catalog size must be positive")
    if not directions:
        return {
            "left_positions": [],
            "right_positions": [],
            "left_directions": [],
            "right_directions": [],
            "left_contribution_unique_sizes": [],
            "right_contribution_unique_sizes": [],
            "left_projected_state_product": 1,
            "right_projected_state_product": 1,
            "maximum_projected_side_product": 1,
            "partitions_considered_after_side_symmetry": 1,
        }

    best = None
    considered = 0
    positions = tuple(range(len(directions)))
    for count in range(1, len(directions) + 1):
        for left in itertools.combinations(positions, count):
            if 0 not in left:
                continue
            considered += 1
            left_set = set(left)
            right = tuple(position for position in positions if position not in left_set)
            left_product = math.prod(sizes[position] for position in left)
            right_product = math.prod(sizes[position] for position in right)
            score = (
                max(left_product, right_product),
                left_product + right_product,
                abs(left_product - right_product),
                abs(len(left) - len(right)),
                tuple(directions[position] for position in left),
            )
            if best is None or score < best[0]:
                best = (score, left, right, left_product, right_product)
    require(best is not None, "balanced partition search found no partition")
    _score, left, right, left_product, right_product = best
    require(considered == 2 ** (len(directions) - 1), "side-symmetry partition count changed")
    return {
        "left_positions": list(left),
        "right_positions": list(right),
        "left_directions": [directions[position] for position in left],
        "right_directions": [directions[position] for position in right],
        "left_contribution_unique_sizes": [sizes[position] for position in left],
        "right_contribution_unique_sizes": [sizes[position] for position in right],
        "left_projected_state_product": left_product,
        "right_projected_state_product": right_product,
        "maximum_projected_side_product": max(left_product, right_product),
        "partitions_considered_after_side_symmetry": considered,
    }


def reconstruction() -> dict:
    orbits, orbit_source = affine.parent.load_z7_orbits()
    leaves_by_orbit, leaf_audit = affine.parent.exact_mean_leaves(orbits)
    require(sum(len(rows) for rows in leaves_by_orbit) == 2_160, "leaf reconstruction changed")
    kernel_rows, hull_bases, hull_audit = affine.build_hull_audit()
    source_catalog = affine.canonical_catalog(7, 4)
    anchors = affine.AnchorFactory(kernel_rows, source_catalog[1:] - source_catalog[0])
    normalizations, star_split_audit = affine.pointed.stabilizer_and_star_audit(
        orbits, leaves_by_orbit
    )
    _translation, systems, translation_audit, system_audits = affine.build_pointed_systems(
        orbits, normalizations
    )
    factories = [
        {
            branch: affine.JointProjectionFactory(
                systems[orbit_index][branch], hull_bases, anchors
            )
            for branch in ("A", "B")
        }
        for orbit_index in range(2)
    ]
    return {
        "orbits": orbits,
        "orbit_source": orbit_source,
        "leaves_by_orbit": leaves_by_orbit,
        "leaf_audit": leaf_audit,
        "kernel_rows": kernel_rows,
        "hull_bases": hull_bases,
        "hull_audit": hull_audit,
        "anchors": anchors,
        "star_split_audit": star_split_audit,
        "systems": systems,
        "translation_audit": translation_audit,
        "system_audits": system_audits,
        "factories": factories,
    }


def validate_reconstruction_against_parent(parent_payload: dict, rebuilt: dict) -> dict:
    comparisons = {
        "degree_two_zero_mean_hull_audit": rebuilt["hull_audit"],
        "translation_equivariant_linear_system": rebuilt["translation_audit"],
        "pointed_linear_systems": rebuilt["system_audits"],
        "pointed_star_split_audit": rebuilt["star_split_audit"],
        "mean_leaf_coverage": rebuilt["leaf_audit"],
        "orbit_source": rebuilt["orbit_source"],
    }
    rows = []
    for key, recomputed in comparisons.items():
        stored = parent_payload[key]
        require(recomputed == stored, f"reconstructed parent audit changed: {key}")
        rows.append(
            {
                "field": key,
                "stored_sha256": json_sha256(stored),
                "recomputed_sha256": json_sha256(recomputed),
                "exact_match": True,
            }
        )
    return {
        "parent_audit_fields_recomputed": len(rows),
        "all_recomputed_fields_match_parent_exactly": True,
        "fields": rows,
    }


def catalog_size(catalog_class: str) -> int:
    if catalog_class == "S":
        return 56
    if catalog_class == "M":
        return 1_764
    raise AssertionError(f"non-enumerable catalog class {catalog_class}")


def deduplicate_contribution_rows(signatures: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    signatures = np.ascontiguousarray(signatures, dtype=np.uint8)
    require(signatures.ndim == 2 and signatures.shape[1] > 0, "contribution shape changed")
    keys = raw_row_keys(signatures)
    unique_keys, first_indices = np.unique(keys, return_index=True)
    unique_keys = np.ascontiguousarray(unique_keys)
    unique_rows = unique_keys.view(np.uint8).reshape(-1, signatures.shape[1])
    first_indices = np.ascontiguousarray(first_indices, dtype=np.uint64)
    require(len(unique_rows) == len(first_indices), "contribution representative count changed")
    require(
        np.array_equal(raw_row_keys(unique_rows), unique_keys),
        "contribution raw-key round trip failed",
    )
    return unique_rows, first_indices


class ProjectedContributionStore:
    """Exactly deduplicate each direction before partitioning or budgeting."""

    def __init__(self) -> None:
        self.cache: dict[tuple, dict] = {}
        self.records: dict[tuple, dict] = {}

    def get(
        self,
        orbit_index: int,
        branch: str,
        orbit: dict,
        leaf: dict,
        factory: object,
        high: tuple[int, ...],
        direction: int,
    ) -> dict:
        mask = int(orbit["masks"][direction])
        mean = int(leaf["scaled_means"][direction])
        catalog_class = str(leaf["catalog_classes"][direction])
        key = (
            int(orbit_index),
            str(branch),
            tuple(high),
            int(direction),
            mask,
            mean,
        )
        if key not in self.cache:
            signatures, identifier = factory.catalog_signature(
                high, direction, mask, mean
            )
            signatures = np.ascontiguousarray(signatures, dtype=np.uint8)
            complete_size = catalog_size(catalog_class)
            require(len(signatures) == complete_size, "complete contribution size changed")
            unique_rows, first_indices = deduplicate_contribution_rows(signatures)
            require(
                np.all(first_indices < complete_size),
                "contribution representative catalog index out of range",
            )
            audit_key = (tuple(high), int(direction), mask, mean)
            source_audit = factory.signature_audit_records[audit_key]
            require(source_audit["identifier"] == identifier, "source contribution identifier changed")
            record = {
                "orbit_index": int(orbit_index),
                "pointed_star_branch": str(branch),
                "omitted_high_directions": list(high),
                "omitted_high_direction_count": len(high),
                "direction": int(direction),
                "catalog_class": catalog_class,
                "mask": mask,
                "scaled_mean": mean,
                "identifier": identifier,
                "complete_catalog_rows": complete_size,
                "projected_unique_contribution_states": len(unique_rows),
                "discarded_projected_duplicate_rows": complete_size - len(unique_rows),
                "compression_ratio": len(unique_rows) / complete_size,
                "complete_signature_sha256_uint8": affine.matrix_sha256(signatures),
                "sorted_unique_signature_sha256_uint8": affine.matrix_sha256(unique_rows),
                "representative_catalog_indices_sha256_uint64": affine.matrix_sha256(
                    first_indices[None, :]
                ),
                "source_exact_catalog_sha256_int16": source_audit[
                    "exact_catalog_sha256_int16"
                ],
                "same_exact_catalog_row_order_used_mod3_and_mod7": True,
                "exact_raw_byte_deduplication": True,
                "multiplicity_retained": False,
            }
            self.cache[key] = {
                "rows": unique_rows,
                "first_catalog_indices": first_indices,
                "metadata": record,
            }
            self.records[key] = record
        return self.cache[key]

    def audit(self) -> dict:
        records = [self.records[key] for key in sorted(self.records)]
        compression = Counter(
            (
                row["omitted_high_direction_count"],
                row["catalog_class"],
                row["complete_catalog_rows"],
                row["projected_unique_contribution_states"],
            )
            for row in records
        )
        expected_orbit0_a = {
            (1, "M"): 1_716,
            (2, "M"): 927,
            (3, "M"): 136,
            (2, "S"): 56,
            (3, "S"): 44,
        }
        observed_orbit0_a: dict[tuple[int, str], set[int]] = {}
        for row in records:
            if row["orbit_index"] == 0 and row["pointed_star_branch"] == "A":
                observed_orbit0_a.setdefault(
                    (row["omitted_high_direction_count"], row["catalog_class"]),
                    set(),
                ).add(row["projected_unique_contribution_states"])
        for key, expected in expected_orbit0_a.items():
            require(
                observed_orbit0_a.get(key) == {expected},
                f"measured orbit0/A contribution compression changed for {key}",
            )
        return {
            "unique_projected_contributions": len(records),
            "complete_catalog_rows_total_over_unique_contributions": sum(
                row["complete_catalog_rows"] for row in records
            ),
            "projected_unique_states_total_over_unique_contributions": sum(
                row["projected_unique_contribution_states"] for row in records
            ),
            "compression_histogram": {
                f"H{high}_{catalog_class}_{complete}_to_{unique}": count
                for (high, catalog_class, complete, unique), count in sorted(compression.items())
            },
            "measured_orbit0_branchA_compression": {
                f"H{high}_{catalog_class}": sorted(values)
                for (high, catalog_class), values in sorted(observed_orbit0_a.items())
            },
            "required_measured_orbit0_branchA_compression": {
                f"H{high}_{catalog_class}": value
                for (high, catalog_class), value in expected_orbit0_a.items()
            },
            "all_required_measured_compressions_match": True,
            "all_contribution_records_sha256": json_sha256(records),
            "multiplicity_discarded_before_partition_and_budget": True,
            "records": records,
        }


def validate_parent_survivor(row: dict, rebuilt: dict) -> tuple[dict, dict, dict, object]:
    orbit_index = int(row["branch_orbit_index"])
    leaf_index = int(row["orbit_leaf_index"])
    branch = str(row["pointed_star_branch"])
    require(orbit_index in (0, 1) and branch in ("A", "B"), "survivor branch key changed")
    orbit = rebuilt["orbits"][orbit_index]
    leaf = rebuilt["leaves_by_orbit"][orbit_index][leaf_index]
    system = rebuilt["systems"][orbit_index][branch]
    factory = rebuilt["factories"][orbit_index][branch]
    expected_key = f"orbit{orbit_index}_leaf{leaf_index}_branch{branch}"
    require(row["case_key"] == expected_key, "survivor case key disagrees with indices")
    require(tuple(row["representative_finite_field"]) == tuple(orbit["representative"]), "orbit representative changed")
    require(tuple(row["q_values"]) == tuple(leaf["q_values"]), "survivor q-vector changed")
    require(tuple(row["scaled_means"]) == tuple(leaf["scaled_means"]), "survivor means changed")
    require(tuple(row["catalog_classes"]) == tuple(leaf["catalog_classes"]), "survivor classes changed")
    require(tuple(row["catalog_levels"]) == tuple(leaf["catalog_levels"]), "survivor levels changed")
    require(tuple(row["high_directions_relaxed_to_exact_affine_hulls"]) == tuple(leaf["high_directions"]), "high directions changed")
    require(tuple(row["enumerated_directions"]) == tuple(leaf["enumerated_directions"]), "enumerated directions changed")
    require(row["fixed_edge_rows"] == system["fixed_edge_rows"], "pointed fixed-edge rows changed")
    require(row["passes_all_necessary_projections"] is True, "input row is not a parent survivor")
    return orbit, leaf, system, factory


def make_case_plans(
    survivor_rows: list[dict],
    rebuilt: dict,
    max_side_states: int,
    contributions: ProjectedContributionStore,
) -> list[dict]:
    plans = []
    for row in survivor_rows:
        orbit, leaf, system, factory = validate_parent_survivor(row, rebuilt)
        orbit_index = int(row["branch_orbit_index"])
        branch = str(row["pointed_star_branch"])
        enumerated = tuple(int(value) for value in leaf["enumerated_directions"])
        high = tuple(int(value) for value in leaf["high_directions"])
        complete_sizes = tuple(
            catalog_size(leaf["catalog_classes"][direction])
            for direction in enumerated
        )
        projected = {
            direction: contributions.get(
                orbit_index,
                branch,
                orbit,
                leaf,
                factory,
                high,
                direction,
            )
            for direction in enumerated
        }
        unique_sizes = tuple(len(projected[direction]["rows"]) for direction in enumerated)
        partition = balanced_partition(enumerated, unique_sizes)
        complete_catalog_partition = balanced_partition(enumerated, complete_sizes)
        conditioned = {
            modulus: factory.primes[modulus].conditioner(high)[1]
            for modulus in MODULI
        }
        dimensions = {
            modulus: int(conditioned[modulus]["conditioned_dependency_dimension"])
            for modulus in MODULI
        }
        joint_dimension = dimensions[3] + dimensions[7]
        projected_states = (
            int(partition["left_projected_state_product"])
            + int(partition["right_projected_state_product"])
        )
        eligible = int(partition["maximum_projected_side_product"]) <= max_side_states
        plans.append(
            {
                "case_key": row["case_key"],
                "parent_row": row,
                "orbit": orbit,
                "leaf": leaf,
                "system": system,
                "factory": factory,
                "branch_orbit_index": int(row["branch_orbit_index"]),
                "source_orbit_index": int(row["source_orbit_index"]),
                "orbit_leaf_index": int(row["orbit_leaf_index"]),
                "pointed_star_branch": row["pointed_star_branch"],
                "pattern": tuple(int(value) for value in leaf["pattern"]),
                "high_directions": high,
                "enumerated_directions": enumerated,
                "complete_catalog_sizes": complete_sizes,
                "projected_unique_contribution_sizes": unique_sizes,
                "projected_contributions": projected,
                "partition": partition,
                "complete_catalog_balanced_partition": complete_catalog_partition,
                "conditioned_dimensions": dimensions,
                "joint_signature_dimension": joint_dimension,
                "estimated_projected_signature_states": projected_states,
                "estimated_projected_signature_bytes": projected_states * joint_dimension,
                "budget_eligible": eligible,
                "budget_skip_reason": None
                if eligible
                else "balanced_projected_unique_side_product_exceeds_max_side_states",
            }
        )
    require(len(plans) == len(survivor_rows), "case-plan census changed")
    require(len({plan["case_key"] for plan in plans}) == len(plans), "case plan repeated")
    return plans


def public_plan(plan: dict) -> dict:
    h, s, m = plan["pattern"]
    return {
        "case_key": plan["case_key"],
        "branch_orbit_index": plan["branch_orbit_index"],
        "source_orbit_index": plan["source_orbit_index"],
        "orbit_leaf_index": plan["orbit_leaf_index"],
        "pointed_star_branch": plan["pointed_star_branch"],
        "catalog_pattern": f"H{h}_S{s}_M{m}",
        "high_directions": list(plan["high_directions"]),
        "enumerated_directions": list(plan["enumerated_directions"]),
        "enumerated_catalog_classes": [
            plan["leaf"]["catalog_classes"][direction]
            for direction in plan["enumerated_directions"]
        ],
        "complete_catalog_sizes": list(plan["complete_catalog_sizes"]),
        "projected_unique_contribution_sizes": list(
            plan["projected_unique_contribution_sizes"]
        ),
        "complete_catalog_balanced_partition_for_comparison": plan[
            "complete_catalog_balanced_partition"
        ],
        "projected_unique_state_balanced_partition": plan["partition"],
        "conditioned_dimensions": {
            str(modulus): dimension
            for modulus, dimension in plan["conditioned_dimensions"].items()
        },
        "joint_signature_dimension": plan["joint_signature_dimension"],
        "estimated_projected_signature_states": plan[
            "estimated_projected_signature_states"
        ],
        "estimated_projected_signature_bytes": plan[
            "estimated_projected_signature_bytes"
        ],
        "budget_eligible": plan["budget_eligible"],
        "budget_skip_reason": plan["budget_skip_reason"],
    }


def preflight_census(plans: list[dict], max_side_states: int) -> dict:
    eligible = [plan for plan in plans if plan["budget_eligible"]]
    skipped = [plan for plan in plans if not plan["budget_eligible"]]
    no_enumerable = [plan for plan in plans if not plan["enumerated_directions"]]
    enumerable_eligible = [
        plan for plan in eligible if plan["enumerated_directions"]
    ]
    if max_side_states == DEFAULT_MAX_SIDE_STATES:
        if len(plans) == EXPECTED_PARENT_SURVIVORS:
            require(len(eligible) == EXPECTED_DEFAULT_ELIGIBLE, "default eligible census changed")
            require(len(skipped) == EXPECTED_DEFAULT_SKIPPED, "default budget-skip census changed")
            require(len(no_enumerable) == 4 * EXPECTED_NO_ENUMERABLE_REPRESENTATIVES, "full no-enumerable census changed")
        elif len(plans) == EXPECTED_REPRESENTATIVES:
            require(
                len(enumerable_eligible)
                == EXPECTED_DEFAULT_ENUMERABLE_ELIGIBLE_REPRESENTATIVES,
                "default tractable representative census changed",
            )
            require(
                len(skipped) == EXPECTED_DEFAULT_SKIPPED_REPRESENTATIVES,
                "default large-product representative census changed",
            )
            require(
                len(no_enumerable) == EXPECTED_NO_ENUMERABLE_REPRESENTATIVES,
                "default no-enumerable representative census changed",
            )

    by_pattern: dict[str, Counter[str]] = {}
    maximum_histogram: Counter[int] = Counter()
    for plan in plans:
        h, s, m = plan["pattern"]
        key = f"H{h}_S{s}_M{m}"
        rows = by_pattern.setdefault(key, Counter())
        rows["parent_survivors"] += 1
        rows["budget_eligible" if plan["budget_eligible"] else "budget_skipped"] += 1
        rows["estimated_projected_signature_states"] += plan[
            "estimated_projected_signature_states"
        ]
        rows["estimated_projected_signature_bytes"] += plan[
            "estimated_projected_signature_bytes"
        ]
        maximum_histogram[
            int(plan["partition"]["maximum_projected_side_product"])
        ] += 1

    public_plans = [public_plan(plan) for plan in plans]
    over_budget = [
        {
            "case_key": plan["case_key"],
            "catalog_pattern": public_plan(plan)["catalog_pattern"],
            "left_projected_state_product": plan["partition"][
                "left_projected_state_product"
            ],
            "right_projected_state_product": plan["partition"][
                "right_projected_state_product"
            ],
            "maximum_projected_side_product": plan["partition"][
                "maximum_projected_side_product"
            ],
            "max_side_states": max_side_states,
            "reason": plan["budget_skip_reason"],
        }
        for plan in skipped
    ]
    return {
        "selected_parent_survivor_universe": len(plans),
        "budget_eligible_cases": len(eligible),
        "budget_eligible_enumerable_cases": len(enumerable_eligible),
        "no_enumerable_direction_cases": len(no_enumerable),
        "budget_skipped_cases": len(skipped),
        "max_side_states": max_side_states,
        "largest_eligible_projected_side_product": max(
            (
                plan["partition"]["maximum_projected_side_product"]
                for plan in eligible
            ),
            default=0,
        ),
        "largest_eligible_estimated_projected_signature_bytes": max(
            (plan["estimated_projected_signature_bytes"] for plan in eligible),
            default=0,
        ),
        "maximum_projected_side_product_histogram": {
            str(product): count for product, count in sorted(maximum_histogram.items())
        },
        "coverage_by_catalog_pattern": {
            key: dict(sorted(rows.items())) for key, rows in sorted(by_pattern.items())
        },
        "all_case_budget_plans_sha256": json_sha256(public_plans),
        "over_budget_case_keys_sha256": json_sha256([row["case_key"] for row in over_budget]),
        "over_budget_cases": over_budget,
        "case_budget_plans": public_plans,
    }


def smoke_selection(plans: list[dict], smoke_max_side_states: int) -> tuple[list[dict], dict]:
    selected: list[dict] = []
    reasons: dict[str, list[str]] = {}

    def add(plan: dict | None, reason: str) -> None:
        if plan is None:
            return
        if plan not in selected:
            selected.append(plan)
        reasons.setdefault(plan["case_key"], []).append(reason)

    eligible_smoke = [
        plan
        for plan in plans
        if plan["budget_eligible"]
        and plan["partition"]["maximum_projected_side_product"]
        <= smoke_max_side_states
    ]
    add(next((plan for plan in eligible_smoke if len(plan["enumerated_directions"]) == 0), None), "zero_catalog_base_case")
    add(next((plan for plan in eligible_smoke if len(plan["enumerated_directions"]) == 1), None), "one_catalog_case")
    add(
        next(
            (
                plan
                for plan in eligible_smoke
                if len(plan["enumerated_directions"]) >= 3
                and {plan["leaf"]["catalog_classes"][d] for d in plan["enumerated_directions"]}
                == {"S", "M"}
            ),
            None,
        ),
        "mixed_small_medium_global_join",
    )
    add(
        max(
            (
                plan
                for plan in eligible_smoke
                if len(plan["enumerated_directions"]) >= 2
                and all(plan["leaf"]["catalog_classes"][d] == "S" for d in plan["enumerated_directions"])
            ),
            key=lambda plan: plan["partition"]["maximum_projected_side_product"],
            default=None,
        ),
        "largest_all_small_join_within_smoke_cap",
    )
    add(next((plan for plan in plans if not plan["budget_eligible"]), None), "explicit_global_budget_skip")
    require(selected, "smoke selection is empty")
    require(len(selected) <= 5, "smoke selection bound changed")
    return selected, {
        "smoke_max_side_states": smoke_max_side_states,
        "selected_case_count": len(selected),
        "selection_reasons": [
            {"case_key": plan["case_key"], "reasons": reasons[plan["case_key"]]}
            for plan in selected
        ],
        "selected_case_keys_sha256": json_sha256([plan["case_key"] for plan in selected]),
        "full_parent_survivor_coverage_claimed": False,
    }


def mixed_radix_strides(sizes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(math.prod(sizes[index + 1 :]) for index in range(len(sizes)))


def decode_flat_index(flat_index: int, sizes: tuple[int, ...]) -> tuple[int, ...]:
    strides = mixed_radix_strides(sizes)
    return tuple((int(flat_index) // stride) % size for size, stride in zip(sizes, strides))


def enumerate_signature_set(
    matrices: tuple[np.ndarray, ...],
    component_dimensions: dict[int, int],
    state_limit: int,
    chunk_states: int,
) -> dict:
    """Enumerate one side exactly, then discard signature multiplicity."""
    d3 = int(component_dimensions[3])
    d7 = int(component_dimensions[7])
    width = d3 + d7
    require(d3 > 0 and d7 > 0, "signature component dimension vanished")
    for matrix in matrices:
        require(matrix.ndim == 2 and matrix.shape[1] == width, "side signature width changed")
        require(np.all(matrix[:, :d3] < 3), "mod-3 signature byte out of range")
        require(np.all(matrix[:, d3:] < 7), "mod-7 signature byte out of range")
        require(
            len(np.unique(raw_row_keys(matrix))) == len(matrix),
            "side received a contribution that was not deduplicated first",
        )
    sizes = tuple(len(matrix) for matrix in matrices)
    raw_count = math.prod(sizes) if sizes else 1
    require(raw_count <= state_limit, "side enumeration exceeded its audited state limit")
    require(chunk_states > 0, "chunk size must be positive")
    strides = mixed_radix_strides(sizes)
    states = np.empty((raw_count, width), dtype=np.uint8)
    for start in range(0, raw_count, chunk_states):
        stop = min(raw_count, start + chunk_states)
        flat = np.arange(start, stop, dtype=np.int64)
        mod3 = np.zeros((stop - start, d3), dtype=np.uint8)
        mod7 = np.zeros((stop - start, d7), dtype=np.uint8)
        for matrix, size, stride in zip(matrices, sizes, strides):
            indices = (flat // stride) % size
            np.add(mod3, matrix[indices, :d3], out=mod3)
            np.remainder(mod3, 3, out=mod3)
            np.add(mod7, matrix[indices, d3:], out=mod7)
            np.remainder(mod7, 7, out=mod7)
        states[start:stop, :d3] = mod3
        states[start:stop, d3:] = mod7

    raw_hash = affine.matrix_sha256(states)
    keys = raw_row_keys(states)
    unique_keys, first_indices = np.unique(keys, return_index=True)
    unique_keys = np.ascontiguousarray(unique_keys)
    unique_rows = unique_keys.view(np.uint8).reshape(-1, width)
    first_indices = np.ascontiguousarray(first_indices, dtype=np.uint64)
    require(len(unique_rows) == len(first_indices), "deduplicated representative index count changed")
    require(np.array_equal(raw_row_keys(unique_rows), unique_keys), "unique signature key round-trip failed")
    return {
        "rows": unique_rows,
        "first_flat_indices": first_indices,
        "keys": unique_keys,
        "contribution_unique_sizes": sizes,
        "raw_cartesian_states": raw_count,
        "unique_signature_states": len(unique_rows),
        "discarded_duplicate_states": raw_count - len(unique_rows),
        "signature_width_bytes": width,
        "raw_cartesian_states_sha256_uint8": raw_hash,
        "sorted_unique_signatures_sha256_uint8": affine.matrix_sha256(unique_rows),
        "representative_flat_indices_sha256_uint64": affine.matrix_sha256(
            first_indices[None, :]
        ),
        "exact_raw_byte_deduplication": True,
        "every_input_contribution_was_already_exactly_deduplicated": True,
        "multiplicity_retained": False,
    }


def public_signature_set(state_set: dict) -> dict:
    return {
        key: value
        for key, value in state_set.items()
        if key not in {"rows", "first_flat_indices", "keys"}
    }


def meet_signature_sets(
    bases: dict[int, np.ndarray],
    left: dict,
    right: dict,
    chunk_states: int,
) -> dict:
    """Intersect exact signature sets and preserve one deterministic witness."""
    d3 = len(bases[3])
    d7 = len(bases[7])
    width = d3 + d7
    left_rows = left["rows"]
    right_rows = right["rows"]
    left_keys = left["keys"]
    require(left_rows.shape[1] == right_rows.shape[1] == width, "join width changed")
    require(np.array_equal(raw_row_keys(left_rows), left_keys), "left key ordering changed")
    match_count = 0
    first_match: tuple[int, int] | None = None
    certificate = hashlib.sha256()
    for start in range(0, len(right_rows), chunk_states):
        stop = min(len(right_rows), start + chunk_states)
        rows = right_rows[start:stop]
        needed3 = (-bases[3][None, :].astype(np.int16) - rows[:, :d3].astype(np.int16)) % 3
        needed7 = (-bases[7][None, :].astype(np.int16) - rows[:, d3:].astype(np.int16)) % 7
        needed = np.ascontiguousarray(np.concatenate((needed3, needed7), axis=1), dtype=np.uint8)
        needed_keys = raw_row_keys(needed)
        positions = np.searchsorted(left_keys, needed_keys)
        candidates = np.flatnonzero(positions < len(left_keys))
        if not len(candidates):
            continue
        equal = left_keys[positions[candidates]] == needed_keys[candidates]
        hits = candidates[equal]
        if not len(hits):
            continue
        left_indices = positions[hits].astype(np.uint64)
        right_indices = (start + hits).astype(np.uint64)
        match_count += len(hits)
        certificate.update(left_indices.astype("<u8", copy=False).tobytes())
        certificate.update(right_indices.astype("<u8", copy=False).tobytes())
        if first_match is None:
            first_match = (int(left_indices[0]), int(right_indices[0]))

    base_row = np.concatenate((bases[3], bases[7])).astype(np.uint8)[None, :]
    return {
        "matching_unique_signature_pairs": match_count,
        "first_matching_unique_indices": list(first_match) if first_match is not None else None,
        "matching_pair_index_certificate_sha256": certificate.hexdigest(),
        "concatenated_base_sha256_uint8": affine.matrix_sha256(base_row),
        "mod3_dimension": d3,
        "mod7_dimension": d7,
        "joint_signature_dimension": width,
        "exact_raw_byte_intersection": True,
        "hash_collision_assumption_used": False,
    }


def verify_recovered_witness(
    plan: dict,
    bases: dict[int, np.ndarray],
    contributions: dict[int, dict],
    left: dict,
    right: dict,
    join: dict,
    kernel_rows: np.ndarray,
) -> dict | None:
    indices = join["first_matching_unique_indices"]
    if indices is None:
        return None
    left_flat = int(left["first_flat_indices"][indices[0]])
    right_flat = int(right["first_flat_indices"][indices[1]])
    partition = plan["partition"]
    left_directions = tuple(partition["left_directions"])
    right_directions = tuple(partition["right_directions"])
    left_sizes = tuple(partition["left_contribution_unique_sizes"])
    right_sizes = tuple(partition["right_contribution_unique_sizes"])
    selected = dict(zip(left_directions, decode_flat_index(left_flat, left_sizes)))
    selected.update(zip(right_directions, decode_flat_index(right_flat, right_sizes)))
    require(set(selected) == set(plan["enumerated_directions"]), "witness missed an enumerable direction")

    d3 = len(bases[3])
    syndrome3 = bases[3].astype(np.int64).copy()
    syndrome7 = bases[7].astype(np.int64).copy()
    exact_rows = []
    direction_rows = []
    for direction in plan["enumerated_directions"]:
        unique_index = int(selected[direction])
        contribution = contributions[direction]
        signature = contribution["rows"]
        require(
            0 <= unique_index < len(signature),
            "decoded projected contribution index out of range",
        )
        catalog_index = int(contribution["first_catalog_indices"][unique_index])
        syndrome3 += signature[unique_index, :d3]
        syndrome7 += signature[unique_index, d3:]
        mask = int(plan["orbit"]["masks"][direction])
        mean = int(plan["leaf"]["scaled_means"][direction])
        catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
        require(0 <= catalog_index < len(catalog), "catalog representative index out of range")
        exact_row = catalog[catalog_index]
        require(2 * int(exact_row.sum()) == 5 * mean, "witness catalog row changed exact mean")
        require(not np.any(kernel_rows @ exact_row), "witness catalog row is not exactly degree two")
        require(
            np.array_equal(exact_row % 2, affine.parity_for_mask(mask)),
            "witness catalog row changed parity",
        )
        exact_rows.append(exact_row)
        direction_rows.append(
            {
                "direction": direction,
                "catalog_class": plan["leaf"]["catalog_classes"][direction],
                "catalog_size": len(catalog),
                "projected_unique_contribution_index": unique_index,
                "representative_exact_catalog_row_index": catalog_index,
                "mask": mask,
                "scaled_mean": mean,
                "selected_exact_row_sha256_int64": affine.matrix_sha256(exact_row[None, :]),
            }
        )
    require(not np.any(syndrome3 % 3), "recovered witness fails mod 3")
    require(not np.any(syndrome7 % 7), "recovered witness fails mod 7")
    exact_matrix = (
        np.stack(exact_rows)
        if exact_rows
        else np.empty((0, 35), dtype=np.int64)
    )
    witness = {
        "left_representative_flat_index": left_flat,
        "right_representative_flat_index": right_flat,
        "catalog_rows_by_direction": direction_rows,
        "selected_exact_catalog_rows_sha256_int64": affine.matrix_sha256(exact_matrix),
        "same_row_index_used_for_each_direction_in_mod3_and_mod7": True,
        "all_selected_rows_have_exact_leaf_means_parity_and_degree_two": True,
        "recovered_joint_syndrome_is_zero_mod3_and_mod7": True,
    }
    witness["witness_certificate_sha256"] = json_sha256(witness)
    return witness


def global_join_self_audit() -> dict:
    """Brute-force a manufactured same-index witness and an impossible target."""
    dimensions = {3: 2, 7: 2}
    matrices = (
        np.asarray([[0, 0, 0, 0], [0, 1, 1, 0], [0, 1, 1, 0]], dtype=np.uint8),
        np.asarray([[0, 0, 0, 0], [0, 2, 6, 0]], dtype=np.uint8),
        np.asarray([[0, 0, 0, 0], [0, 1, 2, 0]], dtype=np.uint8),
    )
    deduplicated = tuple(deduplicate_contribution_rows(matrix) for matrix in matrices)
    unique_matrices = tuple(row[0] for row in deduplicated)
    representative_indices = tuple(row[1] for row in deduplicated)
    require(tuple(len(row) for row in unique_matrices) == (2, 2, 2), "synthetic contribution dedup changed")
    chosen = (1, 1, 1)
    total3 = sum((matrix[index, :2].astype(np.int64) for matrix, index in zip(unique_matrices, chosen)), np.zeros(2, dtype=np.int64))
    total7 = sum((matrix[index, 2:].astype(np.int64) for matrix, index in zip(unique_matrices, chosen)), np.zeros(2, dtype=np.int64))
    bases = {3: (-total3 % 3).astype(np.uint8), 7: (-total7 % 7).astype(np.uint8)}
    partition = balanced_partition((0, 1, 2), tuple(len(matrix) for matrix in unique_matrices))
    left_matrices = tuple(unique_matrices[position] for position in partition["left_positions"])
    right_matrices = tuple(unique_matrices[position] for position in partition["right_positions"])
    left = enumerate_signature_set(left_matrices, dimensions, 100, 10)
    right = enumerate_signature_set(right_matrices, dimensions, 100, 10)
    joined = meet_signature_sets(bases, left, right, 10)

    brute_positive = 0
    for indices in itertools.product(*(range(len(matrix)) for matrix in unique_matrices)):
        syndrome3 = bases[3].astype(np.int64).copy()
        syndrome7 = bases[7].astype(np.int64).copy()
        for matrix, index in zip(unique_matrices, indices):
            syndrome3 += matrix[index, :2]
            syndrome7 += matrix[index, 2:]
        brute_positive += int(not np.any(syndrome3 % 3) and not np.any(syndrome7 % 7))
    require(brute_positive > 0 and joined["matching_unique_signature_pairs"] > 0, "manufactured join witness was lost")
    require(sum(len(raw) - len(unique) for raw, unique in zip(matrices, unique_matrices)) > 0, "self-audit did not deduplicate contributions")

    impossible_bases = {3: np.asarray([1, 0], dtype=np.uint8), 7: np.asarray([0, 0], dtype=np.uint8)}
    impossible = meet_signature_sets(impossible_bases, left, right, 10)
    brute_impossible = 0
    for indices in itertools.product(*(range(len(matrix)) for matrix in unique_matrices)):
        syndrome3 = impossible_bases[3].astype(np.int64).copy()
        syndrome7 = impossible_bases[7].astype(np.int64).copy()
        for matrix, index in zip(unique_matrices, indices):
            syndrome3 += matrix[index, :2]
            syndrome7 += matrix[index, 2:]
        brute_impossible += int(not np.any(syndrome3 % 3) and not np.any(syndrome7 % 7))
    require(brute_impossible == impossible["matching_unique_signature_pairs"] == 0, "impossible target passed")
    same_index_audit = affine.joint_join_self_audit()
    require(same_index_audit["same_index_false_positive_trap"]["joint_same_index_matches"] == 0, "same-index trap changed")
    return {
        "passed": True,
        "manufactured_catalog_row_indices": list(chosen),
        "manufactured_representative_original_row_indices": [
            int(representative_indices[direction][index])
            for direction, index in enumerate(chosen)
        ],
        "manufactured_brute_force_catalog_tuple_matches": brute_positive,
        "manufactured_unique_signature_pair_matches": joined["matching_unique_signature_pairs"],
        "manufactured_witness_found": True,
        "impossible_brute_force_catalog_tuple_matches": brute_impossible,
        "impossible_unique_signature_pair_matches": impossible["matching_unique_signature_pairs"],
        "contribution_level_exact_deduplication_exercised_before_partition": True,
        "complete_contribution_sizes": [len(row) for row in matrices],
        "projected_unique_contribution_sizes": [len(row) for row in unique_matrices],
        "side_sum_exact_deduplication_exercised": True,
        "duplicate_multiplicity_intentionally_discarded": True,
        "same_index_cross_prime_false_positive_trap": same_index_audit[
            "same_index_false_positive_trap"
        ],
        "balanced_partition": partition,
        "left_signature_set": public_signature_set(left),
        "right_signature_set": public_signature_set(right),
        "manufactured_join": joined,
        "impossible_join": impossible,
    }


def execute_case(
    plan: dict,
    max_side_states: int,
    execution_side_cap: int,
    chunk_states: int,
    kernel_rows: np.ndarray,
    smoke_test: bool,
) -> dict:
    public = public_plan(plan)
    base_record = {
        "case_key": plan["case_key"],
        "branch_orbit_index": plan["branch_orbit_index"],
        "source_orbit_index": plan["source_orbit_index"],
        "orbit_leaf_index": plan["orbit_leaf_index"],
        "pointed_star_branch": plan["pointed_star_branch"],
        "catalog_pattern": public["catalog_pattern"],
        "residue_pair_minus_plus": plan["leaf"]["residue_pair"],
        "scaled_means": list(plan["leaf"]["scaled_means"]),
        "catalog_classes": list(plan["leaf"]["catalog_classes"]),
        "high_directions_relaxed_to_exact_affine_hulls": list(plan["high_directions"]),
        "enumerated_directions_joined_globally": list(plan["enumerated_directions"]),
        "complete_catalog_sizes": list(plan["complete_catalog_sizes"]),
        "projected_unique_contribution_sizes": list(
            plan["projected_unique_contribution_sizes"]
        ),
        "projected_unique_state_balanced_partition": plan["partition"],
        "complete_catalog_balanced_partition_for_comparison": plan[
            "complete_catalog_balanced_partition"
        ],
        "joint_signature_dimension": plan["joint_signature_dimension"],
        "max_side_states": max_side_states,
    }
    if not plan["budget_eligible"]:
        result = {
            **base_record,
            "decision_status": "skipped_side_state_budget",
            "exact_global_join_processed": False,
            "rigorously_rejected": False,
            "necessary_only_survivor": False,
            "skipped": True,
            "skip_reason": plan["budget_skip_reason"],
        }
        result["decision_certificate_sha256"] = json_sha256(result)
        return result
    if plan["partition"]["maximum_projected_side_product"] > execution_side_cap:
        require(smoke_test, "full run tried to bypass an eligible case with an execution cap")
        result = {
            **base_record,
            "decision_status": "skipped_smoke_execution_cap",
            "exact_global_join_processed": False,
            "rigorously_rejected": False,
            "necessary_only_survivor": False,
            "skipped": True,
            "skip_reason": "balanced_projected_unique_side_product_exceeds_smoke_execution_cap",
            "smoke_execution_side_cap": execution_side_cap,
        }
        result["decision_certificate_sha256"] = json_sha256(result)
        return result

    anchor_rhs, raw_syndromes = affine.anchor_rhs_and_raw_syndromes(
        plan["orbit"], plan["leaf"], plan["system"], plan["factory"].anchors
    )
    require(
        affine.matrix_sha256(anchor_rhs[None, :])
        == plan["parent_row"]["anchor_rhs_sha256_int64"],
        "recomputed anchor RHS changed from parent",
    )
    bases = plan["factory"].projected_bases(raw_syndromes, plan["high_directions"])
    require(
        {modulus: len(bases[modulus]) for modulus in MODULI}
        == plan["conditioned_dimensions"],
        "global-only-high conditioner dimension changed",
    )

    contributions: dict[int, dict] = {}
    signature_rows = []
    for direction, complete_size, unique_size in zip(
        plan["enumerated_directions"],
        plan["complete_catalog_sizes"],
        plan["projected_unique_contribution_sizes"],
    ):
        mask = int(plan["orbit"]["masks"][direction])
        mean = int(plan["leaf"]["scaled_means"][direction])
        contribution = plan["projected_contributions"][direction]
        signature = contribution["rows"]
        metadata = contribution["metadata"]
        require(len(signature) == unique_size, "projected contribution size changed")
        require(
            signature.shape[1] == plan["joint_signature_dimension"],
            "global signature width changed",
        )
        catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
        require(len(catalog) == complete_size, "complete catalog size changed")
        require(
            np.all(2 * catalog.sum(axis=1, dtype=np.int64) == 5 * mean),
            "complete catalog row changed exact mean",
        )
        require(
            not np.any(plan["factory"].anchors.kernel_rows @ catalog.T),
            "complete catalog row left exact degree two",
        )
        contributions[direction] = contribution
        signature_rows.append(
            {
                "direction": direction,
                "catalog_class": plan["leaf"]["catalog_classes"][direction],
                "complete_catalog_rows": complete_size,
                "projected_unique_contribution_states": unique_size,
                "discarded_projected_duplicate_rows": complete_size - unique_size,
                "mask": mask,
                "scaled_mean": mean,
                "identifier": metadata["identifier"],
                "exact_catalog_sha256_int16": metadata[
                    "source_exact_catalog_sha256_int16"
                ],
                "complete_signature_sha256_uint8": metadata[
                    "complete_signature_sha256_uint8"
                ],
                "sorted_unique_signature_sha256_uint8": metadata[
                    "sorted_unique_signature_sha256_uint8"
                ],
                "same_exact_catalog_row_order_used_mod3_and_mod7": True,
                "all_catalog_rows_have_exact_target_mean_and_degree_two": True,
            }
        )

    left_matrices = tuple(
        contributions[direction]["rows"]
        for direction in plan["partition"]["left_directions"]
    )
    right_matrices = tuple(
        contributions[direction]["rows"]
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
        left["raw_cartesian_states"]
        == plan["partition"]["left_projected_state_product"]
        and right["raw_cartesian_states"]
        == plan["partition"]["right_projected_state_product"],
        "enumerated side-state count changed",
    )
    join = meet_signature_sets(bases, left, right, chunk_states)
    rejected = join["matching_unique_signature_pairs"] == 0
    witness = verify_recovered_witness(
        plan,
        bases,
        contributions,
        left,
        right,
        join,
        kernel_rows,
    )
    require((witness is None) == rejected, "join/witness status mismatch")
    result = {
        **base_record,
        "decision_status": (
            "rigorous_global_catalog_join_rejection"
            if rejected
            else "necessary_only_global_catalog_join_survivor"
        ),
        "exact_global_join_processed": True,
        "rigorously_rejected": rejected,
        "necessary_only_survivor": not rejected,
        "skipped": False,
        "skip_reason": None,
        "anchor_rhs_sha256_int64": affine.matrix_sha256(anchor_rhs[None, :]),
        "retained_catalog_signature_audits": signature_rows,
        "left_signature_set": public_signature_set(left),
        "right_signature_set": public_signature_set(right),
        "join": join,
        "recovered_same_index_catalog_witness": witness,
        "high_directions_are_still_affine_hull_relaxations": True,
        "passing_is_binary_edge_feasibility": False,
    }
    result["decision_certificate_sha256"] = json_sha256(result)
    return result


def summarize_results(rows: list[dict]) -> dict:
    counts = Counter()
    by_pattern: dict[str, Counter[str]] = {}
    for row in rows:
        counts["selected"] += 1
        if row["exact_global_join_processed"]:
            counts["processed"] += 1
        if row["rigorously_rejected"]:
            counts["rejected"] += 1
        if row["necessary_only_survivor"]:
            counts["surviving"] += 1
        if row["skipped"]:
            counts["skipped"] += 1
            counts[row["decision_status"]] += 1
        pattern = str(row["catalog_pattern"])
        local = by_pattern.setdefault(pattern, Counter())
        local["selected"] += 1
        local[
            "rejected"
            if row["rigorously_rejected"]
            else "surviving"
            if row["necessary_only_survivor"]
            else "skipped"
        ] += 1
    require(counts["selected"] == counts["processed"] + counts["skipped"], "result processing census changed")
    require(counts["processed"] == counts["rejected"] + counts["surviving"], "join decision census changed")
    return {
        "counts": dict(sorted(counts.items())),
        "coverage_by_catalog_pattern": {
            key: dict(sorted(value.items())) for key, value in sorted(by_pattern.items())
        },
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
        "selected_case_keys_sha256": json_sha256(selected_keys),
        "selected_case_count": len(selected_keys),
        "completed_case_count": len(rows),
        "next_case_index": len(rows),
        "completed_case_keys_sha256": json_sha256([row["case_key"] for row in rows]),
        "case_results_sha256": canonical_case_digest(rows),
        "case_results": rows,
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
    require(payload["experiment"] == EXPERIMENT, "resume output is from another experiment")
    require(payload["run_identity_sha256"] == run_identity, "resume configuration/provenance changed")
    if payload["status"] != "in_progress_atomic_checkpoint":
        require(payload["case_results_sha256"] == canonical_case_digest(payload["case_results"]), "completed output certificate failed")
        return [], payload
    rows = payload["case_results"]
    require(payload["case_results_sha256"] == canonical_case_digest(rows), "checkpoint case certificate failed")
    require([row["case_key"] for row in rows] == selected_keys[: len(rows)], "checkpoint is not a selected-key prefix")
    require(payload["next_case_index"] == len(rows), "checkpoint next index changed")
    for row in rows:
        certificate = row.pop("decision_certificate_sha256")
        require(certificate == json_sha256(row), "checkpoint decision certificate failed")
        row["decision_certificate_sha256"] = certificate
    return rows, None


def run(
    input_path: Path,
    output_path: Path,
    max_side_states: int = DEFAULT_MAX_SIDE_STATES,
    representatives_only: bool = False,
    symmetry_certificate_path: Path | None = None,
    smoke_test: bool = False,
    smoke_max_side_states: int = DEFAULT_SMOKE_MAX_SIDE_STATES,
    chunk_states: int = DEFAULT_CHUNK_STATES,
    checkpoint_every: int = DEFAULT_CHECKPOINT_EVERY,
    resume: bool = False,
    quiet: bool = False,
) -> dict:
    started = time.time()
    require(max_side_states > 0, "max side-state budget must be positive")
    require(smoke_max_side_states > 0, "smoke side-state cap must be positive")
    require(chunk_states > 0, "chunk state count must be positive")
    require(checkpoint_every >= 0, "checkpoint interval cannot be negative")

    parent_payload, all_survivor_rows, input_provenance = load_parent_input(input_path)
    representative_rows = representative_survivors(all_survivor_rows)
    symmetry_audit, symmetry_mapping = load_symmetry_certificate(
        symmetry_certificate_path,
        all_survivor_rows,
        representative_rows,
        input_provenance,
    )
    require(
        symmetry_certificate_path is None or representatives_only,
        "a symmetry certificate is only meaningful with --representatives-only",
    )
    survivor_rows = representative_rows if representatives_only else all_survivor_rows
    rebuilt = reconstruction()
    reconstruction_audit = validate_reconstruction_against_parent(parent_payload, rebuilt)
    contribution_store = ProjectedContributionStore()
    plans = make_case_plans(
        survivor_rows,
        rebuilt,
        max_side_states,
        contribution_store,
    )
    preflight = preflight_census(plans, max_side_states)
    if smoke_test:
        selected, smoke_audit = smoke_selection(plans, smoke_max_side_states)
        execution_side_cap = min(max_side_states, smoke_max_side_states)
    else:
        selected = plans
        smoke_audit = None
        execution_side_cap = max_side_states
    selected_keys = [plan["case_key"] for plan in selected]
    require(len(selected_keys) == len(set(selected_keys)), "selected global case repeated")

    script_path = Path(__file__).resolve()
    affine_path = Path(affine.__file__).resolve()
    source_provenance = {
        "this_script_path": str(script_path),
        "this_script_sha256": file_sha256(script_path),
        "affine_parent_script_path": str(affine_path),
        "affine_parent_script_sha256": file_sha256(affine_path),
        "trusted_imported_helpers": [
            "load_z7_orbits",
            "exact_mean_leaves",
            "build_hull_audit",
            "AnchorFactory",
            "build_pointed_systems",
            "JointProjectionFactory.catalog_signature",
            "anchor_rhs_and_raw_syndromes",
            "mapped_catalog",
            "modular_right_nullspace",
        ],
    }
    configuration = {
        "max_side_states": max_side_states,
        "representatives_only": representatives_only,
        "symmetry_certificate_path": (
            str(symmetry_certificate_path.resolve())
            if symmetry_certificate_path is not None
            else None
        ),
        "symmetry_transfer_validated": symmetry_audit["validated"],
        "smoke_test": smoke_test,
        "smoke_max_side_states": smoke_max_side_states if smoke_test else None,
        "execution_side_cap": execution_side_cap,
        "chunk_states": chunk_states,
        "checkpoint_every": checkpoint_every,
    }
    run_identity = json_sha256(
        {
            "input_file_sha256": input_provenance["file_sha256"],
            "input_all_case_results_sha256": input_provenance["all_case_results_sha256"],
            "input_survivor_keys_sha256": input_provenance["survivor_case_keys_sha256"],
            "source_provenance": source_provenance,
            "symmetry_certificate_audit": symmetry_audit,
            "configuration_without_checkpoint_cadence": {
                key: value for key, value in configuration.items() if key != "checkpoint_every"
            },
            "selected_case_keys": selected_keys,
        }
    )

    rows: list[dict] = []
    if resume:
        rows, completed = load_resume_rows(output_path, run_identity, selected_keys)
        if completed is not None:
            return completed
    resumed_case_count = len(rows)
    for index, plan in enumerate(selected[resumed_case_count:], start=resumed_case_count):
        row = execute_case(
            plan,
            max_side_states,
            execution_side_cap,
            chunk_states,
            rebuilt["kernel_rows"],
            smoke_test,
        )
        rows.append(row)
        if checkpoint_every and len(rows) % checkpoint_every == 0 and len(rows) < len(selected):
            affine.pointed.atomic_write(
                output_path,
                checkpoint_payload(
                    run_identity,
                    input_provenance,
                    configuration,
                    selected_keys,
                    rows,
                ),
            )
            if not quiet:
                print(
                    f"checkpoint {len(rows)}/{len(selected)} cases -> {output_path}",
                    file=sys.stderr,
                    flush=True,
                )
        require(row["case_key"] == selected_keys[index], "processed case order changed")

    summary = summarize_results(rows)
    counts = summary["counts"]
    complete_selected_universe_run = not smoke_test
    if complete_selected_universe_run:
        expected_selected = (
            EXPECTED_REPRESENTATIVES
            if representatives_only
            else EXPECTED_PARENT_SURVIVORS
        )
        require(len(rows) == expected_selected, "complete run missed a selected parent survivor")
        require(counts["processed"] == preflight["budget_eligible_cases"], "eligible case was not processed")
        require(counts["skipped"] == preflight["budget_skipped_cases"], "over-budget skip census changed")
        require(
            all(
                (row["decision_status"] == "skipped_side_state_budget")
                == (not plan["budget_eligible"])
                for row, plan in zip(rows, plans)
            ),
            "full processing did not exactly follow the preflight budget partition",
        )

    all_selected_rejected = (
        complete_selected_universe_run and counts["rejected"] == len(plans)
    )
    transfer_from_representatives_claimed = bool(
        complete_selected_universe_run
        and representatives_only
        and symmetry_mapping is not None
        and symmetry_audit["validated"]
    )
    if representatives_only and transfer_from_representatives_claimed:
        representative_result_by_key = {row["case_key"]: row for row in rows}
        transferred_case_decisions = []
        require(symmetry_mapping is not None, "validated transfer mapping disappeared")
        for representative in [plan["case_key"] for plan in plans]:
            source_result = representative_result_by_key[representative]
            for member in symmetry_mapping[representative]:
                transferred_case_decisions.append(
                    {
                        "case_key": member,
                        "representative_case_key": representative,
                        "representative_decision_certificate_sha256": source_result[
                            "decision_certificate_sha256"
                        ],
                        "decision_status": source_result["decision_status"],
                        "rigorously_rejected": source_result["rigorously_rejected"],
                        "necessary_only_survivor": source_result[
                            "necessary_only_survivor"
                        ],
                        "skipped": source_result["skipped"],
                        "transfer_justified_by_validated_symmetry_certificate": True,
                    }
                )
        require(
            len(transferred_case_decisions) == EXPECTED_PARENT_SURVIVORS,
            "transferred case-decision census changed",
        )
        require(
            {row["case_key"] for row in transferred_case_decisions}
            == {row["case_key"] for row in all_survivor_rows},
            "transferred decisions do not cover the parent survivor set",
        )
        transferred_counts = {
            "processed": counts["processed"] * 4,
            "rejected": counts["rejected"] * 4,
            "surviving": counts["surviving"] * 4,
            "skipped": counts["skipped"] * 4,
        }
        require(sum(transferred_counts[key] for key in ("rejected", "surviving", "skipped")) == EXPECTED_PARENT_SURVIVORS, "transferred result census changed")
    elif not representatives_only and complete_selected_universe_run:
        transferred_case_decisions = None
        transferred_counts = {
            "processed": counts["processed"],
            "rejected": counts["rejected"],
            "surviving": counts["surviving"],
            "skipped": counts["skipped"],
        }
    else:
        transferred_case_decisions = None
        transferred_counts = None
    all_parent_cases_covered = bool(
        complete_selected_universe_run
        and (not representatives_only or transfer_from_representatives_claimed)
    )
    z7_excluded = bool(all_selected_rejected and all_parent_cases_covered)
    self_audit = global_join_self_audit()
    decision_digest = canonical_case_digest(rows)
    rejected_rows = [row for row in rows if row["rigorously_rejected"]]
    surviving_rows = [row for row in rows if row["necessary_only_survivor"]]
    skipped_rows = [row for row in rows if row["skipped"]]
    conditioner_audits = [
        {
            "branch_orbit_index": orbit_index,
            "source_orbit_index": int(rebuilt["orbits"][orbit_index]["source_orbit_index"]),
            "branches": [
                {
                    "pointed_star_branch": branch,
                    "modulus_conditioners": [
                        rebuilt["factories"][orbit_index][branch].primes[modulus].audit()
                        for modulus in MODULI
                    ],
                }
                for branch in ("A", "B")
            ],
        }
        for orbit_index in range(2)
    ]

    result = {
        "experiment": EXPERIMENT,
        "status": (
            "bounded_smoke_test_only"
            if smoke_test
            else "complete_bounded_global_catalog_join_exclusion"
            if z7_excluded
            else "complete_representative_global_catalog_join_without_transfer"
            if representatives_only and not transfer_from_representatives_claimed
            else "complete_bounded_global_catalog_join_sieve_with_skips_or_survivors"
        ),
        "run_identity_sha256": run_identity,
        "p": 7,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": 7,
        "z": 7,
        "phase": 0,
        "moduli": [3, 7],
        "input_provenance": input_provenance,
        "representative_selection": {
            "enabled": representatives_only,
            "representative_definition": "parent survivor with branch_orbit_index=0 and pointed_star_branch=A",
            "representative_case_count": EXPECTED_REPRESENTATIVES,
            "selected_universe_case_count": len(plans),
            "four_case_transfer_claimed": transfer_from_representatives_claimed,
            "transfer_requires_validated_symmetry_certificate": True,
        },
        "symmetry_certificate_audit": symmetry_audit,
        "source_provenance": source_provenance,
        "configuration": configuration,
        "reconstruction_audit": reconstruction_audit,
        "degree_two_zero_mean_hull_audit": rebuilt["hull_audit"],
        "translation_equivariant_linear_system": rebuilt["translation_audit"],
        "pointed_linear_systems": rebuilt["system_audits"],
        "pointed_star_split_audit": rebuilt["star_split_audit"],
        "conditioner_audits": conditioner_audits,
        "projected_contribution_deduplication_audit": contribution_store.audit(),
        "global_join_self_audit": self_audit,
        "preflight_budget_census": preflight,
        "algorithm": {
            "high_directions": "condition away only their full audited degree-two exact-zero-mean affine hulls",
            "enumerable_directions": "retain every S/M direction and choose one complete exact catalog row simultaneously",
            "exact_means": "each complete catalog is fixed at the leaf's exact scaled mean; every row is re-audited before joining",
            "same_index_rule": "each direction's mod-3 and mod-7 signature components are computed from one shared catalog matrix and row index",
            "contribution_deduplication": "deduplicate every direction's concatenated mod-3/mod-7 projected rows exactly before budgeting; retain one exact catalog index per unique row",
            "partition": "exhaustively minimize the larger projected unique-contribution product, breaking side symmetry and ties deterministically",
            "side_enumeration": "lexicographic mixed-radix exact sums in bounded chunks",
            "deduplication": "injective raw-byte signature keys; multiplicity discarded because only existence matters",
            "join": "exact sorted raw-byte lookup of -base-right in the left signature set",
            "rejection": "empty exact global signature intersection rigorously rejects the parent survivor",
            "passing": "one simultaneous exact enumerable-catalog witness exists in the modular quotient; high catalogs and binary edges remain relaxed",
            "representatives_only": "compute only orbit0/branchA signatures; transfer to the other three compact cases only under a validated external symmetry certificate",
        },
        "logical_semantics": {
            "all_enumerable_directions_share_one_global_catalog_assignment": True,
            "same_catalog_row_index_used_mod3_and_mod7_per_direction": True,
            "signature_hash_collision_assumption_used": False,
            "signature_multiplicity_needed_for_existence": False,
            "budget_uses_complete_catalog_products": False,
            "budget_uses_projected_unique_contribution_products": True,
            "high_catalogs_are_exactly_enumerated": False,
            "high_catalogs_are_relaxed_only_to_their_exact_zero_mean_affine_hulls": True,
            "zero_global_join_is_rigorous_rejection": True,
            "passing_is_exact_high_catalog_feasibility": False,
            "passing_is_binary_edge_feasibility": False,
            "smoke_run_can_claim_full_coverage": False,
            "representative_transfer_without_validated_certificate": False,
        },
        "smoke_test": smoke_test,
        "smoke_selection_audit": smoke_audit,
        "complete_selected_universe_run": complete_selected_universe_run,
        "full_1296_cases_computed_directly": bool(
            complete_selected_universe_run and not representatives_only
        ),
        "selected_case_count": len(selected),
        "resumed_case_count": resumed_case_count,
        "processed_exact_global_join_cases": counts["processed"],
        "rigorously_rejected_parent_survivors": counts["rejected"],
        "necessary_only_global_join_survivors": counts["surviving"],
        "skipped_parent_survivors": counts["skipped"],
        "skipped_side_state_budget": counts.get("skipped_side_state_budget", 0),
        "skipped_smoke_execution_cap": counts.get("skipped_smoke_execution_cap", 0),
        "transferred_full_1296_counts": transferred_counts,
        "transferred_case_decisions_sha256": (
            json_sha256(transferred_case_decisions)
            if transferred_case_decisions is not None
            else None
        ),
        "transferred_case_decisions": transferred_case_decisions,
        "representative_results_transferred_to_equivalent_cases": transfer_from_representatives_claimed,
        "result_coverage_by_catalog_pattern": summary["coverage_by_catalog_pattern"],
        "selected_case_keys_sha256": json_sha256(selected_keys),
        "case_results_sha256": decision_digest,
        "rejection_certificates_sha256": json_sha256(
            [row["decision_certificate_sha256"] for row in rejected_rows]
        ),
        "survivor_witnesses_sha256": json_sha256(
            [row["decision_certificate_sha256"] for row in surviving_rows]
        ),
        "skipped_cases_sha256": json_sha256(
            [row["decision_certificate_sha256"] for row in skipped_rows]
        ),
        "case_results": rows,
        "all_budget_eligible_selected_survivors_processed": complete_selected_universe_run,
        "all_parent_survivors_covered_by_join_or_explicit_budget_skip": all_parent_cases_covered,
        "z7_branch_excluded": z7_excluded,
        "checkpoint": {
            "atomic_output_path": str(output_path.resolve()),
            "checkpoint_every_cases": checkpoint_every,
            "resume_supported": True,
            "checkpoint_identity_sha256": run_identity,
            "final_output_is_atomic": True,
        },
        "elapsed_seconds": time.time() - started,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="33 MB full affine-hull JSON")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--max-side-states",
        type=int,
        default=DEFAULT_MAX_SIDE_STATES,
        help=(
            "maximum projected unique-contribution product on either balanced "
            f"side (default: {DEFAULT_MAX_SIDE_STATES})"
        ),
    )
    parser.add_argument(
        "--representatives-only",
        action="store_true",
        help="process only the 324 orbit0/branchA compact representatives",
    )
    parser.add_argument(
        "--symmetry-certificate",
        type=Path,
        help="optional standalone four-case transfer certificate; requires --representatives-only",
    )
    parser.add_argument("--chunk-states", type=int, default=DEFAULT_CHUNK_STATES)
    parser.add_argument("--checkpoint-every", type=int, default=DEFAULT_CHECKPOINT_EVERY)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="run a deterministic bounded subset only")
    parser.add_argument(
        "--smoke-max-side-states",
        type=int,
        default=DEFAULT_SMOKE_MAX_SIDE_STATES,
        help=f"additional smoke execution cap (default: {DEFAULT_SMOKE_MAX_SIDE_STATES})",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = run(
        input_path=args.input,
        output_path=args.output,
        max_side_states=args.max_side_states,
        representatives_only=args.representatives_only,
        symmetry_certificate_path=args.symmetry_certificate,
        smoke_test=args.smoke,
        smoke_max_side_states=args.smoke_max_side_states,
        chunk_states=args.chunk_states,
        checkpoint_every=args.checkpoint_every,
        resume=args.resume,
        quiet=args.quiet,
    )
    affine.pointed.atomic_write(args.output, result)
    if not args.quiet:
        print(
            json.dumps(
                {
                    "experiment": result["experiment"],
                    "status": result["status"],
                    "selected_case_count": result["selected_case_count"],
                    "processed_exact_global_join_cases": result[
                        "processed_exact_global_join_cases"
                    ],
                    "rigorously_rejected_parent_survivors": result[
                        "rigorously_rejected_parent_survivors"
                    ],
                    "necessary_only_global_join_survivors": result[
                        "necessary_only_global_join_survivors"
                    ],
                    "skipped_parent_survivors": result["skipped_parent_survivors"],
                    "z7_branch_excluded": result["z7_branch_excluded"],
                    "elapsed_seconds": result["elapsed_seconds"],
                    "output": str(args.output),
                },
                indent=2,
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
