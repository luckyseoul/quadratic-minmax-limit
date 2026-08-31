#!/usr/bin/env python3
"""Pointed rigorous mod-seven projection sieve for positive p=7, z=7.

For each of the two audited affine-line boundary representatives, the
line stabilizer gives an exhaustive split of the nonempty infinity star:

* branch A meets the line, so a line edge is moved to ``(infinity, 0)``;
* branch B avoids the line, so all seven line edges are zero and a present
  outside edge is moved to one audited outside representative.

Each normalization is imposed by appending exact edge rows to the
translation-equivariant 281-row system.  Complete small/medium catalogs are
then tested with the same necessary omitted-block projections as the parent
z=7 sieve.  A zero projection is a rigorous rejection; every reported
survivor remains only a necessary mod-seven survivor, never an edge lift.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_mod7_projection as parent  # noqa: E402
from p7_infinity7_positive_z2_mod7_join import (  # noqa: E402
    ContributionFactory,
    DIRECTION_TYPES,
    LABELS,
    matrix_sha256,
    modular_rank,
    modular_right_nullspace,
    translation_equivariant_system,
)
from p7_infinity7_positive_zge2_orbits import (  # noqa: E402
    direction_metadata,
    square_affine_semilinear_group,
)
from p7_unsaturated_modular_catalog_filter import equation_matrix  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


P = 7
Q = P * P
MODULUS = 7
EDGE_COUNT = 29
EXPECTED_LEAVES = 2_160
EXPECTED_POINTED_CASES = 2 * EXPECTED_LEAVES
EXPECTED_BRANCH_SHAPES = {"A": (282, 1_225), "B": (289, 1_225)}
EXPECTED_BRANCH_RANKS = {"A": 147, "B": 154}
EXPECTED_DEPENDENCIES = 135
EXPECTED_LINE_STABILIZER = 84


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def atomic_write(path: Path, payload: dict) -> None:
    """Atomically replace ``path`` with fsynced canonical JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _partition_signature(rows: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(tuple(sorted(row)) for row in rows))


def induced_direction_action(
    permutation: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
    partition_lookup: dict[tuple[tuple[int, ...], ...], int],
) -> tuple[int, ...]:
    """Return source-to-target action on the eight affine directions."""
    action = []
    for source in range(P + 1):
        transformed = _partition_signature(
            tuple(
                tuple(permutation[u] for u in range(Q) if labels[source][u] == fibre)
                for fibre in range(P)
            )
        )
        target = partition_lookup.get(transformed)
        require(target is not None, "stabilizer element did not permute directions")
        action.append(int(target))
    require(sorted(action) == list(range(P + 1)), "induced direction action is not bijective")
    return tuple(action)


def stabilizer_and_star_audit(
    orbits: list[dict], leaves_by_orbit: list[list[dict]]
) -> tuple[list[dict], dict]:
    """Audit the two line stabilizers and the exhaustive pointed-star split."""
    require(len(orbits) == len(leaves_by_orbit) == 2, "z=7 orbit/leaf source changed")
    directions, types, labels = direction_metadata()
    require(tuple(types) == tuple(DIRECTION_TYPES), "group and catalog direction types disagree")
    require(tuple(labels) == tuple(LABELS), "group and catalog fibre labels disagree")
    group, recomputed_group_audit = square_affine_semilinear_group(directions, types, labels)

    source = json.loads(parent.ORBIT_EVIDENCE.read_text(encoding="utf-8"))
    source_group = source["group_audit"]
    for key in (
        "group_size",
        "permutation_sha256",
        "induced_direction_action_sha256",
        "all_preserve_paley_difference_signs",
        "all_permute_projective_directions",
        "all_preserve_direction_types",
    ):
        require(recomputed_group_audit[key] == source_group[key], f"group audit field changed: {key}")

    fibres = tuple(
        tuple(tuple(u for u in range(Q) if labels[direction][u] == fibre) for fibre in range(P))
        for direction in range(P + 1)
    )
    partition_lookup = {
        _partition_signature(fibres[direction]): direction for direction in range(P + 1)
    }
    require(len(partition_lookup) == P + 1, "direction partitions are not distinct")

    normalizations = []
    orbit_rows = []
    for orbit_index, (orbit, leaves) in enumerate(zip(orbits, leaves_by_orbit)):
        line = tuple(int(u) for u in orbit["representative"])
        line_set = set(line)
        outside = tuple(u for u in range(Q) if u not in line_set)
        require(len(line) == P and len(line_set) == P, "z=7 representative is not a seven-point line")
        require(0 in line_set, "z=7 line representative does not contain finite point zero")
        require(len(outside) == Q - P, "line/outside partition has wrong size")

        stabilizer = tuple(
            permutation
            for permutation in group
            if {permutation[u] for u in line_set} == line_set
        )
        require(len(stabilizer) == EXPECTED_LINE_STABILIZER, "line stabilizer size changed")
        require(len(stabilizer) == Q * (Q - 1) // int(orbit["size"]), "orbit-stabilizer failed")
        outside_representative = min(outside)
        line_point_orbit = {permutation[0] for permutation in stabilizer}
        outside_point_orbit = {
            permutation[outside_representative] for permutation in stabilizer
        }
        require(line_point_orbit == line_set, "line stabilizer is not transitive on line points")
        require(outside_point_orbit == set(outside), "line stabilizer is not transitive outside the line")

        actions = {
            induced_direction_action(permutation, labels, partition_lookup)
            for permutation in stabilizer
        }
        require(actions, "line stabilizer induced no direction actions")
        for action in actions:
            require(
                all(types[source_direction] == types[target_direction] for source_direction, target_direction in enumerate(action)),
                "line stabilizer changed a direction type",
            )
            require(
                all(
                    orbit["b_values"][source_direction]
                    == orbit["b_values"][target_direction]
                    for source_direction, target_direction in enumerate(action)
                ),
                "line stabilizer changed the representative's direction b-values",
            )

        mean_vectors = {tuple(int(value) for value in leaf["scaled_means"]) for leaf in leaves}
        require(len(mean_vectors) == len(leaves) == 1_080, "z=7 mean leaves are not unique")
        for action in actions:
            for mean in mean_vectors:
                transformed = [0] * (P + 1)
                for source_direction, target_direction in enumerate(action):
                    transformed[target_direction] = mean[source_direction]
                require(tuple(transformed) in mean_vectors, "mean-leaf set is not stabilizer invariant")

        action_digest = hashlib.sha256()
        for action in sorted(actions):
            action_digest.update(bytes(action))
            action_digest.update(b"\xff")
        normalizations.append(
            {
                "line": line,
                "outside": outside,
                "outside_representative": outside_representative,
                "stabilizer": stabilizer,
                "direction_actions": tuple(sorted(actions)),
            }
        )
        orbit_rows.append(
            {
                "branch_orbit_index": orbit_index,
                "source_orbit_index": int(orbit["source_orbit_index"]),
                "representative_finite_field": list(line),
                "line_stabilizer_size": len(stabilizer),
                "line_point_orbit_size": len(line_point_orbit),
                "outside_point_orbit_size": len(outside_point_orbit),
                "branch_A_representative_finite_point": 0,
                "branch_B_representative_finite_point": outside_representative,
                "induced_direction_action_count": len(actions),
                "induced_direction_action_sha256": action_digest.hexdigest(),
                "all_1080_mean_vectors_closed_under_stabilizer": True,
                "mean_vector_orbit_reduction_applied": False,
            }
        )

    # This is stronger than needed: it partitions every nonempty 49-edge
    # star, while an actual boundary star has odd cardinality.
    all_nonempty_stars = 2**Q - 1
    branch_a_nonempty_stars = (2**P - 1) * 2 ** (Q - P)
    branch_b_nonempty_stars = 2 ** (Q - P) - 1
    require(
        branch_a_nonempty_stars + branch_b_nonempty_stars == all_nonempty_stars,
        "A/B nonempty-star split is not exhaustive",
    )
    all_odd_stars = 2 ** (Q - 1)
    branch_b_odd_stars = 2 ** (Q - P - 1)
    branch_a_odd_stars = all_odd_stars - branch_b_odd_stars
    require(branch_a_odd_stars + branch_b_odd_stars == all_odd_stars, "odd-star split failed")

    return normalizations, {
        "group_source": {
            "path": str(parent.ORBIT_EVIDENCE.relative_to(ROOT)),
            "group_size": source_group["group_size"],
            "permutation_sha256": source_group["permutation_sha256"],
            "recomputed_group_matches_committed_evidence": True,
        },
        "branches": {
            "A": "star meets line; move one present line edge to infinity--finite-point-0",
            "B": "star avoids line; fix all seven line edges to zero and move one present outside edge to the audited representative",
        },
        "nonempty_star_partition": {
            "all": all_nonempty_stars,
            "A_meets_line": branch_a_nonempty_stars,
            "B_avoids_line_but_is_nonempty": branch_b_nonempty_stars,
        },
        "odd_boundary_star_partition": {
            "all": all_odd_stars,
            "A_meets_line": branch_a_odd_stars,
            "B_avoids_line": branch_b_odd_stars,
        },
        "A_and_B_disjoint": True,
        "A_and_B_cover_every_nonempty_star": True,
        "boundary_membership_makes_the_infinity_star_odd_and_nonempty": True,
        "point_normalizations_valid_by_stabilizer_transitivity": True,
        "mean_vector_set_invariant_under_every_induced_stabilizer_action": True,
        "mean_vector_orbit_reduction_applied": False,
        "per_orbit": orbit_rows,
    }


def edge_row(edge_index: int, edge_variables: int) -> np.ndarray:
    row = np.zeros(edge_variables, dtype=np.int16)
    row[edge_index] = 1
    return row


def build_branch_systems(
    translation_matrix: np.ndarray,
    translation_audit: dict,
    orbits: list[dict],
    normalizations: list[dict],
) -> tuple[list[dict[str, dict]], list[dict]]:
    """Append exact branch rows and recompute all left dependencies mod seven."""
    data = geometry(P, "affine")
    edges = tuple(tuple(int(value) for value in edge) for edge in data["edges"])
    edge_to_index = {edge: index for index, edge in enumerate(edges)}
    require(len(edges) == 1_225 and len(edge_to_index) == len(edges), "edge geometry changed")
    require(
        tuple(edges[edge_to_index[(0, u + 1)]] for u in range(Q))
        == tuple((0, u + 1) for u in range(Q)),
        "infinity-star edge indexing changed",
    )

    full_matrix = equation_matrix()
    require(full_matrix.shape == (282, 1_225), "full source matrix shape changed")
    require(
        np.array_equal(translation_matrix, np.concatenate((full_matrix[:1], full_matrix[2:]), axis=0)),
        "281-row matrix is not exactly the full system without its fixed-edge row",
    )
    require(translation_audit["rank"] == 146, "translation-equivariant rank changed")

    systems_by_orbit: list[dict[str, dict]] = []
    audits = []
    for orbit_index, (orbit, normalization) in enumerate(zip(orbits, normalizations)):
        line = tuple(normalization["line"])
        outside_representative = int(normalization["outside_representative"])
        require(tuple(orbit["representative"]) == line, "normalization line changed")

        specs = {
            "A": ((0,), (1,)),
            "B": ((*line, outside_representative), (*(0 for _ in line), 1)),
        }
        orbit_systems: dict[str, dict] = {}
        orbit_audits = []
        for branch in ("A", "B"):
            finite_points, tail_rhs_values = specs[branch]
            require(len(finite_points) == len(tail_rhs_values), "branch row/RHS count mismatch")
            require(len(set(finite_points)) == len(finite_points), "branch repeats a fixed edge")
            appended = np.stack(
                [edge_row(edge_to_index[(0, point + 1)], len(edges)) for point in finite_points]
            )
            matrix = np.ascontiguousarray(np.vstack((translation_matrix, appended)), dtype=np.int16)
            dependencies, rank_from_nullspace = modular_right_nullspace(matrix.T, MODULUS)
            dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
            rank = modular_rank(matrix, MODULUS)

            require(matrix.shape == EXPECTED_BRANCH_SHAPES[branch], f"branch {branch} matrix shape changed")
            require(rank == rank_from_nullspace == EXPECTED_BRANCH_RANKS[branch], f"branch {branch} rank changed")
            require(dependencies.shape == (EXPECTED_DEPENDENCIES, matrix.shape[0]), f"branch {branch} dependencies changed")
            require(
                not np.any(dependencies @ (matrix.astype(np.int64) % MODULUS) % MODULUS),
                f"branch {branch} left-null audit failed",
            )
            require(
                rank - int(translation_audit["rank"]) == len(finite_points),
                f"branch {branch} appended edge rows are not independent",
            )
            if branch == "A":
                require(np.array_equal(appended[0], full_matrix[1]), "branch A is not the distinguished point-zero edge")
            else:
                require(set(finite_points[:-1]) == set(line), "branch B did not fix every line edge")
                require(outside_representative not in line, "branch B representative lies on the line")
                require(tuple(tail_rhs_values[:-1]) == (0,) * P and tail_rhs_values[-1] == 1, "branch B RHS changed")

            calibration = (17 * np.arange(matrix.shape[1], dtype=np.int64) + 3) % 2
            manufactured_rhs = matrix.astype(np.int64) @ calibration % MODULUS
            require(
                not np.any(dependencies @ manufactured_rhs % MODULUS),
                f"branch {branch} rejected a manufactured consistent RHS",
            )

            base_rhs = np.zeros(matrix.shape[0], dtype=np.int64)
            base_rhs[0] = EDGE_COUNT
            base_rhs[len(translation_matrix) :] = np.asarray(tail_rhs_values, dtype=np.int64)
            fixed_edge_rows = [
                {
                    "finite_field_point": int(point),
                    "graph_edge": [0, int(point) + 1],
                    "edge_variable_index": int(edge_to_index[(0, int(point) + 1)]),
                    "rhs": int(rhs),
                }
                for point, rhs in zip(finite_points, tail_rhs_values)
            ]
            system = {
                "branch": branch,
                "matrix": matrix,
                "dependencies": dependencies,
                "base_rhs": base_rhs,
                "fixed_edge_rows": fixed_edge_rows,
            }
            audit = {
                "branch": branch,
                "equations": int(matrix.shape[0]),
                "edge_variables": int(matrix.shape[1]),
                "rank": rank,
                "left_dependency_dimension": len(dependencies),
                "translation_equivariant_prefix_rows": len(translation_matrix),
                "appended_exact_edge_rows": len(appended),
                "rank_increment_from_appended_rows": rank - int(translation_audit["rank"]),
                "direction_block_offset": 1,
                "direction_block_width": 35,
                "edge_count_rhs": EDGE_COUNT,
                "fixed_edge_rows": fixed_edge_rows,
                "left_null_audit": True,
                "manufactured_rhs_calibration": True,
                "matrix_sha256": matrix_sha256(matrix),
                "dependency_sha256": matrix_sha256(dependencies.astype(np.uint8)),
                "base_rhs_sha256": matrix_sha256(base_rhs[None, :].astype(np.uint8)),
            }
            system["audit"] = audit
            orbit_systems[branch] = system
            orbit_audits.append(audit)
        systems_by_orbit.append(orbit_systems)
        audits.append(
            {
                "branch_orbit_index": orbit_index,
                "source_orbit_index": int(orbit["source_orbit_index"]),
                "systems": orbit_audits,
            }
        )

    require(len(systems_by_orbit) == 2, "did not build systems for both z=7 lines")
    require(all(set(rows) == {"A", "B"} for rows in systems_by_orbit), "a pointed branch system is missing")
    return systems_by_orbit, audits


class BranchProjectionFactory:
    """Catalog projections for one exact pointed branch row system."""

    def __init__(self, system: dict):
        self.system = system
        self.dependencies = np.asarray(system["dependencies"], dtype=np.int64)
        require(self.dependencies.shape[0] == EXPECTED_DEPENDENCIES, "branch dependency dimension changed")
        self.contributions = ContributionFactory(self.dependencies)
        self.conditioners: dict[tuple[int, ...], tuple[np.ndarray, dict]] = {}
        self.projected_catalogs: dict[tuple[tuple[int, ...], int, int, int], np.ndarray] = {}

    def conditioner(self, omitted: tuple[int, ...]) -> tuple[np.ndarray, dict]:
        omitted = tuple(sorted(int(direction) for direction in omitted))
        require(len(set(omitted)) == len(omitted), "omitted direction set has duplicates")
        require(all(0 <= direction < P + 1 for direction in omitted), "bad omitted direction")
        if omitted not in self.conditioners:
            if omitted:
                columns = np.concatenate(
                    [
                        np.arange(1 + 35 * direction, 1 + 35 * (direction + 1))
                        for direction in omitted
                    ]
                )
                block = self.dependencies[:, columns]
                coefficients, rank_from_nullspace = modular_right_nullspace(block.T, MODULUS)
                block_rank = modular_rank(block, MODULUS)
                require(block_rank == rank_from_nullspace, "omitted-block rank computations disagree")
            else:
                columns = np.empty(0, dtype=np.int64)
                coefficients = np.eye(EXPECTED_DEPENDENCIES, dtype=np.int64)
                block_rank = 0
            coefficients = np.ascontiguousarray(coefficients, dtype=np.int64)
            conditioned = np.ascontiguousarray(coefficients @ self.dependencies % MODULUS, dtype=np.int64)
            dimension = int(len(coefficients))
            require(dimension == EXPECTED_DEPENDENCIES - block_rank, "conditioner rank-nullity failed")
            require(dimension > 0, "omitted blocks consumed every dependency coordinate")
            require(modular_rank(coefficients, MODULUS) == dimension, "conditioner basis lost rank")
            require(modular_rank(conditioned, MODULUS) == dimension, "conditioned dependencies lost rank")
            require(
                not np.any(conditioned @ (self.system["matrix"].astype(np.int64) % MODULUS) % MODULUS),
                "conditioned row is not a dependency of the pointed matrix",
            )
            if omitted:
                require(not np.any(conditioned[:, columns]), "conditioner failed to annihilate omitted blocks")
            metadata = {
                "omitted_directions": list(omitted),
                "omitted_direction_count": len(omitted),
                "omitted_block_rank": int(block_rank),
                "conditioned_dependency_dimension": dimension,
                "all_omitted_blocks_annihilated": True,
                "conditioned_left_null_audit": True,
                "coefficient_sha256": matrix_sha256(coefficients.astype(np.uint8)),
            }
            self.conditioners[omitted] = (coefficients, metadata)
        return self.conditioners[omitted]

    def catalog(self, omitted: tuple[int, ...], direction: int, mask: int, mean: int) -> np.ndarray:
        key = (tuple(sorted(omitted)), int(direction), int(mask), int(mean))
        if key not in self.projected_catalogs:
            coefficients, _metadata = self.conditioner(key[0])
            complete_catalog = self.contributions.get(direction, mask, mean).astype(np.int64)
            self.projected_catalogs[key] = np.ascontiguousarray(
                coefficients @ complete_catalog % MODULUS,
                dtype=np.uint8,
            )
        return self.projected_catalogs[key]

    def audit(self) -> dict:
        histogram: Counter[tuple[int, int, int]] = Counter()
        digest = hashlib.sha256()
        for omitted, (_coefficients, row) in sorted(self.conditioners.items()):
            key = (
                int(row["omitted_direction_count"]),
                int(row["omitted_block_rank"]),
                int(row["conditioned_dependency_dimension"]),
            )
            histogram[key] += 1
            digest.update(
                json.dumps(
                    {"omitted": omitted, **row},
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("ascii")
            )
            digest.update(b"\n")
        return {
            "branch": self.system["branch"],
            "conditioner_count": len(self.conditioners),
            "projected_catalog_cache_entries": len(self.projected_catalogs),
            "observed_conditioner_histogram": {
                f"omit{count}_rank{rank}_dim{dimension}": value
                for (count, rank, dimension), value in sorted(histogram.items())
            },
            "conditioner_audit_sha256": digest.hexdigest(),
            "every_used_omitted_set_rank_audited": True,
            "every_used_omitted_block_fully_annihilated": True,
            "every_conditioned_row_left_null_for_branch_matrix": True,
            "affine_span_substitution_for_omitted_catalogs_used": False,
        }


def pointed_relaxed_state(
    orbit: dict,
    leaf: dict,
    system: dict,
    factory: BranchProjectionFactory,
) -> dict:
    """Run every necessary retained-pair projection for one pointed leaf."""
    classes = tuple(leaf["catalog_classes"])
    high = tuple(int(direction) for direction in leaf["high_directions"])
    enumerated = tuple(int(direction) for direction in leaf["enumerated_directions"])
    require(len(classes) == P + 1, "catalog class vector changed length")
    require(all(classes[direction] == "H" for direction in high), "high catalog escaped omission")
    require(
        all(classes[direction] in ("S", "M") for direction in enumerated),
        "enumerated catalog is neither small nor medium",
    )

    dependencies = factory.dependencies
    fixed = dependencies @ np.asarray(system["base_rhs"], dtype=np.int64) % MODULUS
    for direction, catalog_class in enumerate(classes):
        if catalog_class != "U":
            continue
        contribution = factory.contributions.get(
            direction,
            int(orbit["masks"][direction]),
            int(leaf["scaled_means"][direction]),
        )
        require(contribution.shape[1] == 1, "fixed direction catalog is not unique")
        fixed = (fixed + contribution[:, 0]) % MODULUS

    plans = []
    for retained in parent.retained_sets(enumerated):
        omitted = tuple(sorted(set(high) | (set(enumerated) - set(retained))))
        require(set(high) <= set(omitted), "a high block was not omitted")
        require(set(retained).isdisjoint(omitted), "retained and omitted directions overlap")
        require(
            set(high) | set(enumerated) == set(omitted) | set(retained),
            "variable direction coverage gap",
        )
        coefficients, conditioning = factory.conditioner(omitted)
        plans.append(
            (
                -int(conditioning["conditioned_dependency_dimension"]),
                tuple(retained),
                omitted,
                coefficients,
                conditioning,
            )
        )
    plans.sort(key=lambda row: (row[0], row[1]))
    require(plans, "leaf generated no necessary projection")

    projection_rows = []
    for _negative_dimension, retained, omitted, coefficients, conditioning in plans:
        base = np.ascontiguousarray(coefficients @ fixed % MODULUS, dtype=np.uint8)
        contributions = tuple(
            factory.catalog(
                omitted,
                direction,
                int(orbit["masks"][direction]),
                int(leaf["scaled_means"][direction]),
            )
            for direction in retained
        )
        expected_sizes = tuple(56 if classes[direction] == "S" else 1_764 for direction in retained)
        require(tuple(row.shape[1] for row in contributions) == expected_sizes, "retained catalog size changed")
        count, join = parent.exact_projection_join(base, contributions)
        projection = {
            "retained_directions": list(retained),
            "retained_catalog_classes": [classes[direction] for direction in retained],
            "omitted_directions": list(omitted),
            "omitted_semantics": "entire 35-column direction blocks annihilated",
            "conditioned_dependency_dimension": conditioning["conditioned_dependency_dimension"],
            **join,
        }
        projection_rows.append(projection)
        if count == 0:
            return {
                "passes_all_projections": False,
                "planned_projection_count": len(plans),
                "tested_projection_count": len(projection_rows),
                "projections": projection_rows,
                "first_failed_projection": projection,
            }

    require(len(projection_rows) == len(plans), "survivor skipped a planned projection")
    return {
        "passes_all_projections": True,
        "planned_projection_count": len(plans),
        "tested_projection_count": len(projection_rows),
        "projections": projection_rows,
        "first_failed_projection": None,
    }


def render_case(
    orbit_index: int,
    leaf_index: int,
    orbit: dict,
    leaf: dict,
    system: dict,
    decision: dict,
) -> dict:
    return {
        "branch_orbit_index": orbit_index,
        "source_orbit_index": int(orbit["source_orbit_index"]),
        "orbit_leaf_index": leaf_index,
        "orbit_size": int(orbit["size"]),
        "representative_finite_field": list(orbit["representative"]),
        "pointed_star_branch": system["branch"],
        "fixed_infinity_star_edges": system["fixed_edge_rows"],
        "residue_pair_minus_plus": leaf["residue_pair"],
        "q_values": list(leaf["q_values"]),
        "scaled_means": list(leaf["scaled_means"]),
        "catalog_levels": list(leaf["catalog_levels"]),
        "catalog_classes": list(leaf["catalog_classes"]),
        "catalog_pattern_H_S_M": list(leaf["pattern"]),
        "high_directions_relaxed_by_full_block_annihilation": list(leaf["high_directions"]),
        "enumerated_directions": list(leaf["enumerated_directions"]),
        "necessary_mod7_survivor_only": bool(decision["passes_all_projections"]),
        "exact_edge_lift_claimed": False,
        "planned_projection_count": int(decision["planned_projection_count"]),
        "tested_projection_count": int(decision["tested_projection_count"]),
        "projections": decision["projections"],
    }


def run(smoke_test: bool = False) -> dict:
    started = time.time()
    require(Counter(DIRECTION_TYPES) == Counter({-1: 4, 1: 4}), "direction type census changed")
    translation_matrix, translation_dependencies, translation_audit = translation_equivariant_system()
    require(translation_matrix.shape == (281, 1_225), "translation-equivariant matrix shape changed")
    require(
        translation_audit["rank"] == 146 and translation_dependencies.shape == (135, 281),
        "translation-equivariant rank/dependency census changed",
    )
    require(
        translation_audit["direction_block_offset"] == 1
        and translation_audit["edge_count_rhs"] == EDGE_COUNT,
        "translation-equivariant RHS layout changed",
    )

    orbits, orbit_source = parent.load_z7_orbits()
    leaves_by_orbit, leaf_audit = parent.exact_mean_leaves(orbits)
    require(sum(len(leaves) for leaves in leaves_by_orbit) == EXPECTED_LEAVES, "corrected leaf census changed")
    catalogs = parent.catalog_audit()
    join_audit = parent.join_self_audit()
    normalizations, star_split_audit = stabilizer_and_star_audit(orbits, leaves_by_orbit)
    systems_by_orbit, pointed_linear_audits = build_branch_systems(
        translation_matrix,
        translation_audit,
        orbits,
        normalizations,
    )
    factories = [
        {
            branch: BranchProjectionFactory(systems_by_orbit[orbit_index][branch])
            for branch in ("A", "B")
        }
        for orbit_index in range(len(orbits))
    ]

    selected = (
        parent.smoke_selection(leaves_by_orbit)
        if smoke_test
        else [
            (orbit_index, leaf_index)
            for orbit_index, leaves in enumerate(leaves_by_orbit)
            for leaf_index in range(len(leaves))
        ]
    )
    require(len(set(selected)) == len(selected), "source leaf selection contains duplicates")
    require(all(0 <= orbit_index < 2 for orbit_index, _leaf_index in selected), "bad selected orbit")

    state_cache: dict[tuple, dict] = {}
    case_keys: set[tuple[int, int, str]] = set()
    survivors = []
    rejection_samples = []
    decision_digest = hashlib.sha256()
    source_residues: Counter[str] = Counter()
    source_patterns: Counter[tuple[int, int, int]] = Counter()
    pointed_residues: Counter[str] = Counter()
    rejected_residues: dict[str, Counter[str]] = {"A": Counter(), "B": Counter()}
    logical_projection_histogram: Counter[tuple[str, int, int, int]] = Counter()
    branch_totals: dict[str, Counter[str]] = {"A": Counter(), "B": Counter()}
    per_orbit_branch: list[dict[str, Counter[str]]] = [
        {"A": Counter(), "B": Counter()} for _orbit in orbits
    ]

    for orbit_index, leaf_index in selected:
        orbit = orbits[orbit_index]
        leaf = leaves_by_orbit[orbit_index][leaf_index]
        source_residues[leaf["residue_pair"]] += 1
        source_patterns[tuple(leaf["pattern"])] += 1
        for branch in ("A", "B"):
            case_key = (orbit_index, leaf_index, branch)
            require(case_key not in case_keys, "pointed branch case processed twice")
            case_keys.add(case_key)
            pointed_residues[leaf["residue_pair"]] += 1
            branch_totals[branch]["processed"] += 1
            per_orbit_branch[orbit_index][branch]["processed"] += 1

            system = systems_by_orbit[orbit_index][branch]
            factory = factories[orbit_index][branch]
            state_key = (
                orbit_index,
                branch,
                tuple(leaf["catalog_classes"]),
                tuple(
                    None if catalog_class == "H" else int(mean)
                    for catalog_class, mean in zip(leaf["catalog_classes"], leaf["scaled_means"])
                ),
            )
            if state_key not in state_cache:
                state_cache[state_key] = pointed_relaxed_state(orbit, leaf, system, factory)
            decision = state_cache[state_key]
            for projection in decision["projections"]:
                logical_projection_histogram[
                    (
                        branch,
                        int(projection["retained_catalog_count"]),
                        len(projection["omitted_directions"]),
                        int(projection["conditioned_dependency_dimension"]),
                    )
                ] += 1

            digest_row = {
                "orbit": orbit_index,
                "leaf": leaf_index,
                "branch": branch,
                "residue": leaf["residue_pair"],
                "q": leaf["q_values"],
                "classes": leaf["catalog_classes"],
                "passing": decision["passes_all_projections"],
                "tested": decision["tested_projection_count"],
                "failure": decision["first_failed_projection"],
            }
            decision_digest.update(
                json.dumps(digest_row, sort_keys=True, separators=(",", ":")).encode("ascii")
            )
            decision_digest.update(b"\n")

            if decision["passes_all_projections"]:
                branch_totals[branch]["surviving"] += 1
                per_orbit_branch[orbit_index][branch]["surviving"] += 1
                survivors.append(render_case(orbit_index, leaf_index, orbit, leaf, system, decision))
            else:
                branch_totals[branch]["rejected"] += 1
                per_orbit_branch[orbit_index][branch]["rejected"] += 1
                rejected_residues[branch][leaf["residue_pair"]] += 1
                if sum(sample["pointed_star_branch"] == branch for sample in rejection_samples) < 16:
                    sample = render_case(orbit_index, leaf_index, orbit, leaf, system, decision)
                    sample["projections"] = [decision["first_failed_projection"]]
                    rejection_samples.append(sample)

    source_processed = len(selected)
    pointed_processed = len(case_keys)
    rejected = pointed_processed - len(survivors)
    require(pointed_processed == 2 * source_processed, "A/B branch coverage is not two per source leaf")
    require(pointed_processed == rejected + len(survivors), "pointed decision census mismatch")
    require(
        all(rows["processed"] == source_processed for rows in branch_totals.values()),
        "one pointed branch missed a selected source leaf",
    )
    require(
        all(rows["processed"] == rows["rejected"] + rows["surviving"] for rows in branch_totals.values()),
        "per-branch decision census mismatch",
    )

    full_run = not smoke_test
    if full_run:
        require(source_processed == EXPECTED_LEAVES, "full run did not process all corrected 2,160 leaves")
        require(pointed_processed == EXPECTED_POINTED_CASES, "full run did not process both pointed branches")
        require(dict(sorted(source_residues.items())) == parent.EXPECTED_RESIDUES, "full residue census changed")
        require(
            dict(sorted(pointed_residues.items()))
            == {key: 2 * value for key, value in parent.EXPECTED_RESIDUES.items()},
            "full pointed residue census changed",
        )
        require(
            all(
                per_orbit_branch[orbit_index][branch]["processed"] == 1_080
                for orbit_index in range(2)
                for branch in ("A", "B")
            ),
            "full run missed an orbit/branch leaf",
        )
        require(len(case_keys) == 2 * sum(len(leaves) for leaves in leaves_by_orbit), "coverage key census failed")

    excluded = full_run and not survivors
    weighted_source_processed = sum(orbits[orbit_index]["size"] for orbit_index, _leaf_index in selected)
    conditioning_audits = [
        {
            "branch_orbit_index": orbit_index,
            "source_orbit_index": int(orbits[orbit_index]["source_orbit_index"]),
            "branches": [factories[orbit_index][branch].audit() for branch in ("A", "B")],
        }
        for orbit_index in range(2)
    ]

    return {
        "experiment": "p7_infinity7_positive_z7_pointed_mod7",
        "status": (
            "complete_rigorous_pointed_mod7_projection_exclusion"
            if excluded
            else "complete_rigorous_pointed_mod7_necessary_sieve_with_survivors"
            if full_run
            else "bounded_smoke_test_only"
        ),
        "p": P,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": P,
        "z": 7,
        "phase": 0,
        "translation_equivariant_linear_system": translation_audit,
        "pointed_linear_systems": pointed_linear_audits,
        "orbit_source": orbit_source,
        "star_split_audit": star_split_audit,
        "mean_leaf_coverage": leaf_audit,
        "catalog_row_counts": catalogs,
        "catalog_source": "complete exact Johnson-slice catalogs imported from the audited z=2 implementation",
        "join_self_audit": join_audit,
        "conditioning_audits": conditioning_audits,
        "projection_input": {
            "scope": "all corrected 2,160 exact mean leaves",
            "parent_projection_survivor_prefilter_used": False,
            "parent_projection_evidence_path": None,
            "mean_vector_stabilizer_reduction_used": False,
        },
        "sieve_method": {
            "pointed_branch_A": "append x_(infinity,finite-point-0)=1",
            "pointed_branch_B": "append seven infinity-to-line zeros and one infinity-to-outside representative one",
            "high_blocks": "annihilate the full 35-column direction block, a superset relaxation of the exact high catalog",
            "more_than_two_enumerated_catalogs": "test every retained pair and annihilate every other variable direction block",
            "at_most_two_enumerated_catalogs": "retain all and exact-join complete catalogs",
            "rejection_semantics": "one zero-hit necessary projection rigorously rejects that pointed leaf",
            "survivor_semantics": "passed every planned necessary projection; not an edge lift",
        },
        "smoke_test": smoke_test,
        "smoke_limitations": (
            "Processes one representative of each orbit/residue/catalog-pattern class in both pointed branches; validates plumbing and assertions but cannot exclude z=7."
            if smoke_test
            else None
        ),
        "full_run": full_run,
        "full_source_exact_mean_leaves": EXPECTED_LEAVES,
        "processed_source_exact_mean_leaves": source_processed,
        "full_pointed_branch_cases": EXPECTED_POINTED_CASES,
        "processed_pointed_branch_cases": pointed_processed,
        "rejected_pointed_branch_cases": rejected,
        "surviving_pointed_branch_cases": len(survivors),
        "processed_source_residue_pair_histogram": dict(sorted(source_residues.items())),
        "processed_pointed_residue_pair_histogram": dict(sorted(pointed_residues.items())),
        "rejected_residue_pair_histogram_by_branch": {
            branch: dict(sorted(rows.items())) for branch, rows in rejected_residues.items()
        },
        "processed_source_pattern_histogram": {
            f"H{h}_S{s}_M{m}": count for (h, s, m), count in sorted(source_patterns.items())
        },
        "processed_weighted_boundary_mean_allocations": weighted_source_processed,
        "full_weighted_boundary_mean_allocations": parent.EXPECTED_WEIGHTED_CASES,
        "computed_relaxed_states": len(state_cache),
        "logical_projection_test_histogram": {
            f"branch{branch}_retain{retain}_omit{omit}_dim{dimension}": count
            for (branch, retain, omit, dimension), count in sorted(logical_projection_histogram.items())
        },
        "all_case_decisions_sha256": decision_digest.hexdigest(),
        "all_selected_leaves_tested_in_both_pointed_branches": True,
        "all_corrected_mean_leaves_covered": full_run,
        "exhaustive_star_split_audited": True,
        "zero_projection_is_rigorous_rejection": True,
        "modular_passing_is_edge_feasibility": False,
        "exact_edge_lift_claimed": False,
        "z7_branch_excluded": excluded,
        "per_branch_summary": [
            {
                "pointed_star_branch": branch,
                "processed_cases": rows["processed"],
                "rejected_cases": rows["rejected"],
                "surviving_cases": rows["surviving"],
            }
            for branch, rows in branch_totals.items()
        ],
        "per_orbit_branch_summary": [
            {
                "branch_orbit_index": orbit_index,
                "source_orbit_index": int(orbits[orbit_index]["source_orbit_index"]),
                "orbit_size": int(orbits[orbit_index]["size"]),
                "branches": [
                    {
                        "pointed_star_branch": branch,
                        "processed_cases": per_orbit_branch[orbit_index][branch]["processed"],
                        "rejected_cases": per_orbit_branch[orbit_index][branch]["rejected"],
                        "surviving_cases": per_orbit_branch[orbit_index][branch]["surviving"],
                    }
                    for branch in ("A", "B")
                ],
            }
            for orbit_index in range(2)
        ],
        "rejection_samples": rejection_samples,
        "survivor_cases": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--smoke",
        "--smoke-test",
        dest="smoke_test",
        action="store_true",
        help="run a bounded pattern-representative subset; never claims z=7 closure",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = run(smoke_test=args.smoke_test)
    atomic_write(args.output, result)
    if not args.quiet:
        omitted = {
            "pointed_linear_systems",
            "conditioning_audits",
            "per_branch_summary",
            "per_orbit_branch_summary",
            "rejection_samples",
            "survivor_cases",
        }
        print(
            json.dumps(
                {key: value for key, value in result.items() if key not in omitted},
                indent=2,
                sort_keys=True,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
