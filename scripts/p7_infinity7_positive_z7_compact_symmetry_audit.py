#!/usr/bin/env python3
"""Exact four-case symmetry audit for the positive p=7, z=7 compact model.

This is a structural certificate, not a solver.  It reconstructs the two
pointed A/B systems, canonicalizes every augmented compact dependency row
space modulo 3, 5, 7, and 11, and audits the nonsquare affine transports
between the two boundary representatives.  With ``--input`` it also audits
the canonical full affine-hull evidence, proves that its 1,296 survivors are
exactly 324 complete four-case symmetry classes, and certifies transfer of
the global mod-3/mod-7 same-index catalog join.

The isomorphism transfers compact infeasibility and compact feasibility.
It does not transfer raw A/B binary edge witnesses, and UNKNOWN transfers
nothing.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
import p7_infinity7_positive_z7_mod7_projection as parent  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine_hull  # noqa: E402
import p7_infinity7_positive_z7_pointed_full_cpsat as pointed_cases  # noqa: E402
import p7_infinity7_positive_z7_pointed_mod7 as pointed  # noqa: E402
import p7_infinity7_positive_zge2_orbits as orbit_helpers  # noqa: E402
import p7_size_four_slack_classify as johnson  # noqa: E402
from p7_infinity7_positive_z2_mod7_join import POINTS  # noqa: E402
from p7_unsaturated_modular_catalog_filter import left_dependencies  # noqa: E402


P = 7
Q = P * P
EDGE_COUNT = 4 * P + 1
MODULI = (3, 5, 7, 11)
BRANCHES = ("A", "B")
SIGMA = (4, 7, 1, 0, 2, 3, 5, 6)
AFFINE_PARAMETERS = {
    "A": {"multiplier": 8, "translation": 0},
    "B": {"multiplier": 32, "translation": 24},
}
EXPECTED_EDGE_RANKS = {
    3: {"A": 162, "B": 169},
    5: {"A": 168, "B": 175},
    7: {"A": 147, "B": 154},
    11: {"A": 168, "B": 175},
}
EXPECTED_COMMON_RANKS = {3: 120, 5: 114, 7: 135, 11: 114}
EXPECTED_CANONICAL_RREF_HASHES = {
    3: "4896adb1116b115941a784a3ae3bf5310f2dc73125a70b4a65f5e9e422197ba4",
    5: "f4ec7f92643c521afe188f77d247b9ae2970d7d398e6e07ae1192a46b3c58810",
    7: "10c0a4c13efbd4573bd75c2f24225a907eddeccbdbed37aead5ebbd3c174d609",
    11: "5c9ed22d7a49521cb7323e0b707ce70dc319d4fceb49b745564b8f5b82d853f8",
}
EXPECTED_LEAF_PERMUTATION_SHA256 = (
    "cc5af4ca13a2d16713a2b84e98f693c9825e20d7dba5844c87c05207e0642bae"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def partition_signature(rows: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(tuple(sorted(row)) for row in rows))


def canonical_rref(rows: np.ndarray, modulus: int) -> tuple[np.ndarray, tuple[int, ...]]:
    """Return the unique reduced row-echelon basis over F_modulus."""
    work = np.ascontiguousarray(rows, dtype=np.int64).copy() % modulus
    require(work.ndim == 2, "RREF input is not a matrix")
    pivot_row = 0
    pivots = []
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[pivot_row:, column])
        if not len(candidates):
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            work[[pivot_row, selected]] = work[[selected, pivot_row]]
        inverse = pow(int(work[pivot_row, column]), -1, modulus)
        work[pivot_row] = work[pivot_row] * inverse % modulus
        factors = work[:, column].copy()
        factors[pivot_row] = 0
        active = np.flatnonzero(factors)
        if len(active):
            work[active] = (
                work[active] - factors[active, None] * work[pivot_row]
            ) % modulus
        pivots.append(column)
        pivot_row += 1
        if pivot_row == work.shape[0]:
            break
    require(not np.any(work[pivot_row:]), "RREF left a nonzero row below its rank")
    canonical = np.ascontiguousarray(work[:pivot_row], dtype=np.uint8)
    return canonical, tuple(pivots)


def construct_pointed_systems() -> dict:
    """Reconstruct the exact four pointed systems through repository helpers."""
    cases, case_orbit_source, normalization_audit = pointed_cases.audited_pointed_cases()
    orbits, orbit_source = parent.load_z7_orbits()
    require(len(cases) == 4 and len(orbits) == 2, "four-case source census changed")
    require(
        json_sha256(case_orbit_source) == json_sha256(orbit_source),
        "pointed case and orbit helper sources disagree",
    )
    case_lookup = {
        (int(case["branch_orbit_index"]), str(case["pointed_star_branch"])): case
        for case in cases
    }
    require(
        set(case_lookup) == {(orbit, branch) for orbit in range(2) for branch in BRANCHES},
        "pointed case keys changed",
    )

    normalizations = []
    for orbit_index, orbit in enumerate(orbits):
        line = tuple(int(value) for value in orbit["representative"])
        branch_b = case_lookup[(orbit_index, "B")]
        outside = [
            int(row["finite_field_point"])
            for row in branch_b["fixed_infinity_star_edges"]
            if int(row["value"]) == 1 and int(row["finite_field_point"]) not in line
        ]
        require(len(outside) == 1, "branch B outside normalization changed")
        normalizations.append(
            {"line": line, "outside_representative": outside[0]}
        )

    translation, systems, translation_audit, system_audits = (
        affine_hull.build_pointed_systems(orbits, normalizations)
    )
    require(translation.shape == (281, 1_225), "translation system shape changed")
    require(len(systems) == 2, "pointed system orbit count changed")
    for orbit_index in range(2):
        for branch in BRANCHES:
            system = systems[orbit_index][branch]
            case = case_lookup[(orbit_index, branch)]
            expected_fixed = [
                (tuple(int(value) for value in row["edge"]), int(row["value"]))
                for row in case["fixed_infinity_star_edges"]
            ]
            observed_fixed = [
                (tuple(int(value) for value in row["graph_edge"]), int(row["rhs"]))
                for row in system["fixed_edge_rows"]
            ]
            require(expected_fixed == observed_fixed, "case/system fixed rows disagree")
            require(int(system["base_rhs"][0]) == EDGE_COUNT, "edge-count RHS changed")
            require(
                not np.any(np.asarray(system["base_rhs"])[1:281]),
                "base RHS polluted a direction block",
            )
    return {
        "cases": cases,
        "case_lookup": case_lookup,
        "orbits": orbits,
        "orbit_source": orbit_source,
        "normalization_audit": normalization_audit,
        "normalizations": normalizations,
        "translation": translation,
        "translation_audit": translation_audit,
        "systems": systems,
        "system_audits": system_audits,
    }


def compact_augmented_rows(system: dict, modulus: int) -> tuple[np.ndarray, dict]:
    """Build [280 slack coefficients | constant] for every left dependency."""
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    base_rhs = np.ascontiguousarray(system["base_rhs"], dtype=np.int64)
    edge_rank, dependencies = left_dependencies(matrix, modulus)
    dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
    branch = str(system["branch"])
    expected_dimension = EXPECTED_COMMON_RANKS[modulus]
    require(
        edge_rank == EXPECTED_EDGE_RANKS[modulus][branch],
        f"mod-{modulus} branch {branch} edge rank changed",
    )
    require(
        dependencies.shape == (expected_dimension, matrix.shape[0]),
        f"mod-{modulus} branch {branch} dependency shape changed",
    )
    dependency_rref, _ = canonical_rref(dependencies, modulus)
    require(
        len(dependency_rref) == expected_dimension,
        f"mod-{modulus} dependency basis is not complete",
    )
    require(
        not np.any(dependencies @ (matrix % modulus) % modulus),
        f"mod-{modulus} dependency basis is not left-null",
    )

    direction_blocks = dependencies[:, 1:281]
    require(direction_blocks.shape == (expected_dimension, 280), "direction slice changed")
    constants = dependencies @ base_rhs + 13 * np.sum(direction_blocks, axis=1)
    zero_slack_rhs = base_rhs.copy()
    zero_slack_rhs[1:281] = 13
    require(
        np.array_equal(constants % modulus, dependencies @ zero_slack_rhs % modulus),
        "fixed-edge/base constant construction changed",
    )
    augmented = np.ascontiguousarray(
        np.column_stack((-direction_blocks, constants)) % modulus,
        dtype=np.int64,
    )
    canonical, pivots = canonical_rref(augmented, modulus)
    require(
        len(canonical) == expected_dimension,
        f"mod-{modulus} augmented compact rank changed",
    )
    return augmented, {
        "edge_matrix_rank": int(edge_rank),
        "complete_left_dependency_dimension": int(len(dependencies)),
        "complete_left_dependency_basis_rank": int(len(dependency_rref)),
        "left_null_audit": True,
        "augmented_columns_280_slacks_plus_constant": int(augmented.shape[1]),
        "augmented_compact_rank": int(len(canonical)),
        "pivot_columns": list(pivots),
        "raw_augmented_basis_sha256_uint8": array_sha256(augmented.astype(np.uint8)),
        "canonical_rref_sha256_uint8": array_sha256(canonical),
        "fixed_edge_rhs_included_in_constant": True,
        "bad_count_formula": "13 - slack",
    }


def audit_compact_row_spaces(context: dict) -> tuple[dict, dict, dict]:
    augmented: dict[tuple[int, str, int], np.ndarray] = {}
    canonical: dict[tuple[int, str, int], np.ndarray] = {}
    case_rows = []
    for orbit_index in range(2):
        for branch in BRANCHES:
            system = context["systems"][orbit_index][branch]
            prime_rows = []
            for modulus in MODULI:
                rows, row_audit = compact_augmented_rows(system, modulus)
                reduced, _ = canonical_rref(rows, modulus)
                expected_hash = EXPECTED_CANONICAL_RREF_HASHES[modulus]
                require(
                    row_audit["canonical_rref_sha256_uint8"] == expected_hash,
                    f"mod-{modulus} canonical compact hash changed",
                )
                augmented[(orbit_index, branch, modulus)] = rows
                canonical[(orbit_index, branch, modulus)] = reduced
                prime_rows.append({"modulus": modulus, **row_audit})
            case_rows.append(
                {
                    "case_key": f"orbit{orbit_index}_{branch}",
                    "branch_orbit_index": orbit_index,
                    "pointed_star_branch": branch,
                    "pointed_matrix_shape": list(system["matrix"].shape),
                    "pointed_matrix_sha256_int16": affine_hull.matrix_sha256(
                        system["matrix"]
                    ),
                    "base_rhs_sha256_int64": affine_hull.matrix_sha256(
                        np.asarray(system["base_rhs"], dtype=np.int64)[None, :]
                    ),
                    "prime_audits": prime_rows,
                }
            )

    common_rows = []
    for modulus in MODULI:
        bases = [canonical[(orbit, branch, modulus)] for orbit in range(2) for branch in BRANCHES]
        require(
            all(np.array_equal(bases[0], basis) for basis in bases[1:]),
            f"mod-{modulus} four-case canonical row spaces differ",
        )
        common_hash = array_sha256(bases[0])
        require(
            len(bases[0]) == EXPECTED_COMMON_RANKS[modulus]
            and common_hash == EXPECTED_CANONICAL_RREF_HASHES[modulus],
            f"mod-{modulus} common compact row-space certificate changed",
        )
        common_rows.append(
            {
                "modulus": modulus,
                "common_rank": len(bases[0]),
                "canonical_rref_sha256_uint8": common_hash,
                "all_four_canonical_row_spaces_identical": True,
            }
        )
    return {
        "augmented_row_convention": (
            "columns 0..279 are coefficients of direction-major Johnson slack variables; "
            "column 280 is dependency*fixed_base_rhs + 13*sum(direction coefficients)"
        ),
        "constants_and_all_branch_fixed_edge_rhs_included": True,
        "expected_common_ranks_by_modulus": {
            str(key): value for key, value in EXPECTED_COMMON_RANKS.items()
        },
        "expected_canonical_hashes_by_modulus": {
            str(key): value for key, value in EXPECTED_CANONICAL_RREF_HASHES.items()
        },
        "common_prime_audits": common_rows,
        "four_case_audits": case_rows,
        "all_four_augmented_compact_row_spaces_identical": True,
    }, augmented, canonical


def transported_columns(rows: np.ndarray, coordinate_map: tuple[int, ...]) -> np.ndarray:
    require(len(coordinate_map) == 280, "compact coordinate map length changed")
    transported = np.zeros_like(rows)
    transported[:, -1] = rows[:, -1]
    for source, target in enumerate(coordinate_map):
        transported[:, target] = rows[:, source]
    return transported


def feature_and_kernel_data() -> tuple[np.ndarray, np.ndarray, dict]:
    points, _monomials, exact_evaluation, _left_kernel = johnson.johnson_space()
    require(tuple(points) == tuple(POINTS), "Johnson coordinate order changed")
    require(exact_evaluation.shape == (35, 29), "Johnson evaluation shape changed")
    require(exact_evaluation.rank() == 21, "exact Johnson feature rank changed")
    evaluation = np.asarray(exact_evaluation.tolist(), dtype=np.int64)
    kernel = np.asarray(johnson._primitive_left_kernel_rows(), dtype=np.int64)  # noqa: SLF001
    require(kernel.shape == (14, 35), "primitive Johnson kernel shape changed")
    require(not np.any(kernel @ evaluation), "primitive kernel is not exactly left-null")
    kernel_rref, _ = canonical_rref(kernel, 101)
    require(len(kernel_rref) == 14, "primitive kernel lost exact independence")
    return evaluation, kernel, {
        "johnson_points": 35,
        "degree_at_most_two_feature_columns": 29,
        "exact_feature_rank": 21,
        "exact_left_kernel_dimension": 14,
        "primitive_kernel_sha256_int64": array_sha256(kernel),
    }


def affine_map_audit(
    branch: str,
    context: dict,
    augmented: dict,
    canonical: dict,
    directions: tuple[tuple[int, int], ...],
    types: tuple[int, ...],
    labels: tuple[tuple[int, ...], ...],
    evaluation: np.ndarray,
    kernel: np.ndarray,
) -> tuple[dict, tuple[int, ...]]:
    """Audit one branch's nonsquare affine compact-model transport."""
    q, multiply, add, character, _frob, _norm, irreducible_a, irreducible_b = field_ctx(P)
    require(q == Q, "field order changed")
    parameters = AFFINE_PARAMETERS[branch]
    multiplier = int(parameters["multiplier"])
    translation = int(parameters["translation"])
    permutation = tuple(add(multiply(multiplier, u), translation) for u in range(Q))
    require(len(set(permutation)) == Q, f"f_{branch} is not affine-bijective")
    require(character(multiplier) == -1, f"f_{branch} multiplier is not nonsquare")

    source_line = set(context["orbits"][0]["representative"])
    target_line = set(context["orbits"][1]["representative"])
    require(
        {permutation[u] for u in source_line} == target_line,
        f"f_{branch} does not carry orbit0 line to orbit1 line",
    )
    source_fixed = {
        int(row["finite_field_point"]): int(row["rhs"])
        for row in context["systems"][0][branch]["fixed_edge_rows"]
    }
    target_fixed = {
        int(row["finite_field_point"]): int(row["rhs"])
        for row in context["systems"][1][branch]["fixed_edge_rows"]
    }
    require(
        {permutation[point]: value for point, value in source_fixed.items()} == target_fixed,
        f"f_{branch} does not carry the pointed fixed-edge normalization",
    )

    fibres = tuple(
        tuple(
            tuple(u for u in range(Q) if labels[direction][u] == fibre)
            for fibre in range(P)
        )
        for direction in range(P + 1)
    )
    lookup = {
        partition_signature(fibres[direction]): direction for direction in range(P + 1)
    }
    induced = []
    fibre_maps = []
    johnson_maps = []
    compact_coordinates = [-1] * 280
    point_lookup = {tuple(point): index for index, point in enumerate(POINTS)}
    for source_direction in range(P + 1):
        image_partition = partition_signature(
            tuple(
                tuple(permutation[u] for u in fibres[source_direction][fibre])
                for fibre in range(P)
            )
        )
        target_direction = lookup.get(image_partition)
        require(target_direction is not None, "affine map did not permute directions")
        induced.append(int(target_direction))
        source_to_target_fibre = []
        for source_fibre in range(P):
            images = {
                labels[target_direction][permutation[u]]
                for u in fibres[source_direction][source_fibre]
            }
            require(len(images) == 1, "one source fibre split under affine transport")
            source_to_target_fibre.append(int(next(iter(images))))
        require(
            sorted(source_to_target_fibre) == list(range(P)),
            "induced fibre map is not bijective",
        )
        require(
            all(
                labels[target_direction][permutation[u]]
                == source_to_target_fibre[labels[source_direction][u]]
                for u in range(Q)
            ),
            "pointwise fibre-label transport failed",
        )
        fibre_maps.append(tuple(source_to_target_fibre))

        local_johnson = []
        for source_index, point in enumerate(POINTS):
            target_point = tuple(sorted(source_to_target_fibre[value] for value in point))
            target_index = point_lookup[target_point]
            local_johnson.append(target_index)
            compact_coordinates[35 * source_direction + source_index] = (
                35 * target_direction + target_index
            )
        require(sorted(local_johnson) == list(range(35)), "Johnson map is not bijective")
        johnson_maps.append(tuple(local_johnson))

    require(tuple(induced) == SIGMA, f"f_{branch} direction permutation changed")
    require(sorted(compact_coordinates) == list(range(280)), "compact slack map is not bijective")
    require(
        all(types[target] == -types[source] for source, target in enumerate(SIGMA)),
        "nonsquare direction-type swap failed",
    )

    source_masks = tuple(int(value) for value in context["orbits"][0]["masks"])
    target_masks = tuple(int(value) for value in context["orbits"][1]["masks"])
    source_parity = np.empty(280, dtype=np.uint8)
    target_parity = np.empty(280, dtype=np.uint8)
    floor_rows = []
    for direction in range(8):
        for point_index, point in enumerate(POINTS):
            source_parity[35 * direction + point_index] = (
                sum((source_masks[direction] >> value) & 1 for value in point) & 1
            )
            target_parity[35 * direction + point_index] = (
                sum((target_masks[direction] >> value) & 1 for value in point) & 1
            )
    for source_direction, target_direction in enumerate(SIGMA):
        fibre_map = fibre_maps[source_direction]
        mapped_mask = sum(
            1 << fibre_map[fibre]
            for fibre in range(P)
            if source_masks[source_direction] & (1 << fibre)
        )
        require(mapped_mask == target_masks[target_direction], "boundary mask transport failed")
        source_b = source_masks[source_direction].bit_count()
        target_b = target_masks[target_direction].bit_count()
        source_floor = int(pointed_cases.scaled_direction_floor(P, source_b, 0))
        target_floor = int(pointed_cases.scaled_direction_floor(P, target_b, 0))
        require(source_b == target_b and source_floor == target_floor, "floor transport failed")
        floor_rows.append(
            {
                "source_direction": source_direction,
                "target_direction": target_direction,
                "b": source_b,
                "phase0_floor": source_floor,
            }
        )
    require(
        all(source_parity[source] == target_parity[target] for source, target in enumerate(compact_coordinates)),
        "all 280 exact parity domains did not transport",
    )

    original_kernel_rref, _ = canonical_rref(kernel, 101)
    kernel_rows = []
    for source_direction, local_map in enumerate(johnson_maps):
        transported_kernel = np.zeros_like(kernel)
        for source, target in enumerate(local_map):
            transported_kernel[:, target] = kernel[:, source]
        require(
            not np.any(transported_kernel @ evaluation),
            "transported primitive kernel is not exactly degree-two left-null",
        )
        transported_rref, _ = canonical_rref(transported_kernel, 101)
        require(
            len(transported_rref) == 14
            and np.array_equal(transported_rref, original_kernel_rref),
            "transported primitive kernel does not span the exact same equations",
        )
        kernel_rows.append(
            {
                "source_direction": source_direction,
                "target_direction": SIGMA[source_direction],
                "johnson_permutation_sha256_int64": array_sha256(
                    np.asarray(local_map, dtype=np.int64)
                ),
                "transported_kernel_sha256_int64": array_sha256(transported_kernel),
                "exact_kernel_row_space_preserved": True,
            }
        )

    modular_rows = []
    coordinate_map = tuple(int(value) for value in compact_coordinates)
    for modulus in MODULI:
        source_rows = augmented[(0, branch, modulus)]
        transported = transported_columns(source_rows, coordinate_map)
        transported_rref, _ = canonical_rref(transported, modulus)
        target_rref = canonical[(1, branch, modulus)]
        require(
            np.array_equal(transported_rref, target_rref),
            f"f_{branch} failed mod-{modulus} augmented-row transport",
        )
        modular_rows.append(
            {
                "modulus": modulus,
                "transported_rank": len(transported_rref),
                "transported_canonical_rref_sha256_uint8": array_sha256(
                    transported_rref
                ),
                "equals_target_augmented_row_space": True,
            }
        )

    return {
        "branch": branch,
        "formula": (
            f"f_{branch}(u)={multiplier}u"
            if translation == 0
            else f"f_{branch}(u)={multiplier}u+{translation}"
        ),
        "field_encoding": "u=c0+7*c1",
        "irreducible_polynomial_parameters": {"a": irreducible_a, "b": irreducible_b},
        "multiplier": multiplier,
        "translation": translation,
        "multiplier_quadratic_character": -1,
        "finite_point_permutation_sha256_uint8": array_sha256(
            np.asarray(permutation, dtype=np.uint8)
        ),
        "maps_orbit0_line_to_orbit1_line": True,
        "maps_branch_fixed_edge_normalization": True,
        "direction_permutation_source_to_target": list(induced),
        "direction_types_are_swapped": True,
        "fibre_maps_source_to_target": [list(row) for row in fibre_maps],
        "fibre_maps_sha256_uint8": array_sha256(np.asarray(fibre_maps, dtype=np.uint8)),
        "compact_slack_coordinate_permutation_sha256_int64": array_sha256(
            np.asarray(coordinate_map, dtype=np.int64)
        ),
        "source_parity_sha256_uint8": array_sha256(source_parity),
        "target_parity_sha256_uint8": array_sha256(target_parity),
        "all_boundary_masks_and_280_parity_domains_transport": True,
        "phase0_floor_transport": floor_rows,
        "johnson_kernel_transport": kernel_rows,
        "all_14_exact_primitive_kernel_equations_per_direction_transport": True,
        "directional_sum_and_exact_mean_identity_transport": True,
        "four_per_type_sum_32_constraints_transport_with_type_swap": True,
        "common_residue_variables_0_or_4_swap_between_types": True,
        "augmented_modular_row_space_transport": modular_rows,
        "all_compact_constraint_families_transport": True,
    }, coordinate_map


def audit_leaf_transport(
    orbits: list[dict], leaves_by_orbit: list[list[dict]]
) -> tuple[dict, tuple[int, ...], tuple[tuple[str, ...], ...]]:
    require([len(rows) for rows in leaves_by_orbit] == [1_080, 1_080], "leaf census changed")
    target_lookup = {
        tuple(int(value) for value in leaf["scaled_means"]): index
        for index, leaf in enumerate(leaves_by_orbit[1])
    }
    require(len(target_lookup) == 1_080, "orbit1 mean leaves are not unique")
    permutation = []
    residue_transitions: Counter[str] = Counter()
    for source_index, source_leaf in enumerate(leaves_by_orbit[0]):
        transported_means = [0] * 8
        for source_direction, target_direction in enumerate(SIGMA):
            transported_means[target_direction] = int(source_leaf["scaled_means"][source_direction])
        target_index = target_lookup.get(tuple(transported_means))
        require(target_index is not None, f"orbit0 leaf {source_index} has no orbit1 image")
        target_leaf = leaves_by_orbit[1][target_index]
        expected_residue = source_leaf["residue_pair"][1] + source_leaf["residue_pair"][0]
        require(
            target_leaf["residue_pair"] == expected_residue,
            "minus/plus common residues did not swap",
        )
        for source_direction, target_direction in enumerate(SIGMA):
            require(
                int(target_leaf["q_values"][target_direction])
                == int(source_leaf["q_values"][source_direction]),
                "leaf q-value transport failed",
            )
            require(
                int(target_leaf["catalog_levels"][target_direction])
                == int(source_leaf["catalog_levels"][source_direction]),
                "leaf catalog-level transport failed",
            )
            require(
                target_leaf["catalog_classes"][target_direction]
                == source_leaf["catalog_classes"][source_direction],
                "leaf catalog-class transport failed",
            )
        residue_transitions[
            f"{source_leaf['residue_pair']}->{target_leaf['residue_pair']}"
        ] += 1
        permutation.append(int(target_index))

    require(sorted(permutation) == list(range(1_080)), "leaf transport is not bijective")
    leaf_permutation = tuple(permutation)
    permutation_hash = array_sha256(np.asarray(leaf_permutation, dtype=np.int64))
    require(
        permutation_hash == EXPECTED_LEAF_PERMUTATION_SHA256,
        "Pauli leaf-permutation certificate changed",
    )
    classes = tuple(
        (
            f"orbit0_leaf{source}_branchA",
            f"orbit0_leaf{source}_branchB",
            f"orbit1_leaf{target}_branchA",
            f"orbit1_leaf{target}_branchB",
        )
        for source, target in enumerate(leaf_permutation)
    )
    universe = {
        f"orbit{orbit}_leaf{leaf}_branch{branch}"
        for orbit in range(2)
        for leaf in range(1_080)
        for branch in BRANCHES
    }
    flattened = [key for row in classes for key in row]
    require(
        len(flattened) == len(set(flattened)) == len(universe) == 4_320
        and set(flattened) == universe,
        "four-case classes do not partition all 4,320 pointed leaves",
    )
    return {
        "source_orbit": 0,
        "target_orbit": 1,
        "source_leaves": 1_080,
        "target_leaves": 1_080,
        "direction_permutation_source_to_target": list(SIGMA),
        "leaf_permutation_storage": "target leaf indices in source-leaf order as int64",
        "leaf_permutation_sha256_int64": permutation_hash,
        "expected_leaf_permutation_sha256_int64": EXPECTED_LEAF_PERMUTATION_SHA256,
        "residue_swap_minus_plus": True,
        "residue_transition_histogram": dict(sorted(residue_transitions.items())),
        "all_scaled_means_q_values_levels_and_classes_transport": True,
        "leaf_transport_is_bijective": True,
        "four_case_class_count": len(classes),
        "four_case_class_size": 4,
        "all_4320_pointed_leaf_cases_partitioned_once": True,
        "ordered_four_case_classes_sha256": json_sha256(classes),
    }, leaf_permutation, classes


def optional_survivor_audit(
    input_path: Path | None,
    context: dict,
    leaves_by_orbit: list[list[dict]],
    leaf_audit: dict,
    classes: tuple[tuple[str, ...], ...],
) -> tuple[dict, list[dict] | None, tuple[int, ...]]:
    if input_path is None:
        return (
            {
                "performed": False,
                "reason": "no optional --input affine-hull survivor evidence was supplied",
                "survivor_partition_claimed": False,
            },
            None,
            (),
        )
    require(input_path.is_file(), f"optional survivor evidence does not exist: {input_path}")
    import p7_infinity7_positive_z7_survivor_compact_batch as survivor_batch

    survivors, evidence_audit = survivor_batch.load_and_audit_evidence(
        input_path,
        context["orbits"],
        leaves_by_orbit,
        leaf_audit,
        context["systems"],
        {
            "orbit_source": context["orbit_source"],
            "translation_system_audit": context["translation_audit"],
            "pointed_system_audits": context["system_audits"],
        },
    )
    survivor_keys = {row["case_key"] for row in survivors}
    require(len(survivors) == len(survivor_keys) == 1_296, "survivor key census changed")
    class_size_histogram: Counter[int] = Counter()
    complete_source_indices = []
    covered = set()
    for source_index, row in enumerate(classes):
        present = tuple(key for key in row if key in survivor_keys)
        class_size_histogram[len(present)] += 1
        if present:
            require(len(present) == 4, "survivor evidence contains a partial symmetry class")
            complete_source_indices.append(source_index)
            covered.update(present)
    require(
        class_size_histogram == Counter({0: 756, 4: 324}),
        "survivor four-case class histogram changed",
    )
    require(
        len(complete_source_indices) == 324 and covered == survivor_keys,
        "complete four-case classes do not exactly cover all survivors",
    )
    return (
        {
            "performed": True,
            "input_evidence_audit": evidence_audit,
            "survivor_count": len(survivor_keys),
            "four_case_class_size_histogram": {
                str(size): count for size, count in sorted(class_size_histogram.items())
            },
            "complete_four_case_classes": len(complete_source_indices),
            "empty_four_case_classes": class_size_histogram[0],
            "partial_four_case_classes": 0,
            "complete_class_source_indices_sha256_int64": array_sha256(
                np.asarray(complete_source_indices, dtype=np.int64)
            ),
            "all_1296_survivors_partition_into_exactly_324_complete_four_case_classes": True,
            "experiment_status_full_coverage_and_case_survivor_certificates_reaudited": True,
            "survivor_partition_claimed": True,
        },
        survivors,
        tuple(complete_source_indices),
    )


def local_johnson_maps(coordinate_map: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    rows = []
    for source_direction, target_direction in enumerate(SIGMA):
        local = tuple(
            int(coordinate_map[35 * source_direction + source]) - 35 * target_direction
            for source in range(35)
        )
        require(sorted(local) == list(range(35)), "local Johnson transport is not bijective")
        rows.append(local)
    return tuple(rows)


def transport_johnson_values(values: np.ndarray, local_map: tuple[int, ...]) -> np.ndarray:
    source = np.ascontiguousarray(values)
    require(source.ndim in (1, 2) and source.shape[-1] == 35, "Johnson transport shape changed")
    target = np.empty_like(source)
    for source_index, target_index in enumerate(local_map):
        target[..., target_index] = source[..., source_index]
    return np.ascontiguousarray(target)


def unavailable_global_catalog_join_transfer() -> dict:
    return {
        "scope": "positive_p7_z7_1296_affine_hull_survivors",
        "parent_survivor_case_keys_sha256": None,
        "equivalence_class_count": 0,
        "equivalence_classes": [],
        "equivalence_classes_sha256": None,
        "all_class_maps_preserve_exact_catalog_row_identity_mod3_mod7": False,
        "all_class_maps_preserve_high_affine_hulls_and_anchor_targets": False,
        "transfer_valid_for_global_same_index_catalog_join": False,
        "reason": "global same-index transfer requires the optional audited full survivor input",
    }


def audit_global_catalog_join_transfer(
    *,
    context: dict,
    leaves_by_orbit: list[list[dict]],
    leaf_permutation: tuple[int, ...],
    classes: tuple[tuple[str, ...], ...],
    survivors: list[dict],
    complete_source_indices: tuple[int, ...],
    survivor_audit: dict,
    coordinate_maps: dict[str, tuple[int, ...]],
    augmented: dict,
    canonical: dict,
) -> dict:
    """Prove transfer for the global mod-3/mod-7 same-index catalog join."""
    require(len(survivors) == 1_296, "global-transfer survivor census changed")
    require(len(complete_source_indices) == 324, "global-transfer class census changed")
    survivor_keys = {str(row["case_key"]) for row in survivors}
    representative_keys = {
        str(row["case_key"])
        for row in survivors
        if int(row["branch_orbit_index"]) == 0
        and str(row["pointed_star_branch"]) == "A"
    }
    require(len(representative_keys) == 324, "global-transfer representative census changed")

    canonical_classes = []
    equivalence_classes = []
    covered = set()
    for source_index in complete_source_indices:
        members = tuple(str(value) for value in classes[source_index])
        representative = f"orbit0_leaf{source_index}_branchA"
        require(members[0] == representative, "four-case representative order changed")
        require(
            len(members) == len(set(members)) == 4
            and set(members) <= survivor_keys
            and not (set(members) & covered),
            "global-transfer class is not one disjoint complete survivor class",
        )
        covered.update(members)
        canonical_row = {
            "representative_case_key": representative,
            "member_case_keys": list(members),
        }
        canonical_classes.append(canonical_row)
        equivalence_classes.append(
            {
                **canonical_row,
                "global_same_index_join_equivalence_proved": True,
            }
        )
    require(covered == survivor_keys, "global-transfer classes do not cover all survivors")
    require(
        {row["representative_case_key"] for row in equivalence_classes}
        == representative_keys,
        "global-transfer classes miss an orbit0/branchA representative",
    )

    kernel_rows, hull_bases, hull_audit = affine_hull.build_hull_audit()
    source_hull_catalog = affine_hull.canonical_catalog(7, 4).astype(np.int64)
    anchors = affine_hull.AnchorFactory(
        kernel_rows, source_hull_catalog[1:] - source_hull_catalog[0]
    )
    require(
        hull_audit["exact_zero_mean_dimension"] == 20,
        "global-transfer hull dimension changed",
    )
    branch_local_maps = {
        branch: local_johnson_maps(coordinate_maps[branch]) for branch in BRANCHES
    }

    hull_transport_rows = []
    for branch in BRANCHES:
        for source_direction, target_direction in enumerate(SIGMA):
            local_map = branch_local_maps[branch][source_direction]
            prime_rows = []
            for modulus in (3, 7):
                basis = np.asarray(hull_bases[modulus], dtype=np.int64)
                require(basis.shape == (20, 35), "selected exact hull basis shape changed")
                transported_basis = transport_johnson_values(basis, local_map)
                require(
                    not np.any(transported_basis.sum(axis=1))
                    and not np.any(kernel_rows @ transported_basis.T),
                    "transported hull basis lost exact zero mean or degree two",
                )
                require(
                    affine_hull.modular_rank(transported_basis, modulus) == 20
                    and affine_hull.modular_rank(
                        np.vstack((basis, transported_basis)), modulus
                    )
                    == 20,
                    f"branch {branch} direction {source_direction} changed the mod-{modulus} full hull",
                )
                prime_rows.append(
                    {
                        "modulus": modulus,
                        "full_hull_rank": 20,
                        "transported_basis_sha256_int64": array_sha256(
                            transported_basis
                        ),
                        "transported_span_equals_target_full_hull": True,
                    }
                )
            hull_transport_rows.append(
                {
                    "branch": branch,
                    "source_direction": source_direction,
                    "target_direction": target_direction,
                    "johnson_permutation_sha256_int64": array_sha256(
                        np.asarray(local_map, dtype=np.int64)
                    ),
                    "prime_audits": prime_rows,
                }
            )

    domain_occurrences: Counter[tuple] = Counter()
    class_occurrence_histogram: Counter[str] = Counter()
    for source_index in complete_source_indices:
        source_leaf = leaves_by_orbit[0][source_index]
        target_index = int(leaf_permutation[source_index])
        target_leaf = leaves_by_orbit[1][target_index]
        require(
            classes[source_index][2]
            == f"orbit1_leaf{target_index}_branchA",
            "class target leaf disagrees with leaf permutation",
        )
        for branch in BRANCHES:
            for source_direction, target_direction in enumerate(SIGMA):
                catalog_class = str(source_leaf["catalog_classes"][source_direction])
                require(
                    catalog_class
                    == str(target_leaf["catalog_classes"][target_direction]),
                    "global-transfer catalog class changed",
                )
                source_mask = int(context["orbits"][0]["masks"][source_direction])
                target_mask = int(context["orbits"][1]["masks"][target_direction])
                mean = int(source_leaf["scaled_means"][source_direction])
                require(
                    mean == int(target_leaf["scaled_means"][target_direction]),
                    "global-transfer exact mean changed",
                )
                key = (
                    branch,
                    source_direction,
                    target_direction,
                    source_mask,
                    target_mask,
                    mean,
                    catalog_class,
                )
                domain_occurrences[key] += 1
                class_occurrence_histogram[catalog_class] += 1
    require(
        sum(domain_occurrences.values()) == 324 * 2 * 8,
        "global-transfer directional occurrence coverage changed",
    )
    require(
        set(class_occurrence_histogram) <= {"U", "S", "M", "H"},
        "global-transfer saw an unknown catalog class",
    )

    anchor_internal = {}
    anchor_records = []
    nonliteral_anchor_transports = 0
    for key in sorted(domain_occurrences):
        (
            branch,
            source_direction,
            target_direction,
            source_mask,
            target_mask,
            mean,
            catalog_class,
        ) = key
        local_map = branch_local_maps[branch][source_direction]
        source_anchor = anchors.get(source_mask, mean)
        mapped_anchor = transport_johnson_values(source_anchor, local_map)
        target_anchor = anchors.get(target_mask, mean)
        anchor_shift = np.ascontiguousarray(mapped_anchor - target_anchor, dtype=np.int64)
        require(
            2 * int(mapped_anchor.sum()) == 5 * mean
            and np.array_equal(
                mapped_anchor % 2, affine_hull.parity_for_mask(target_mask)
            )
            and not np.any(kernel_rows @ mapped_anchor),
            "mapped anchor lost exact target mean, parity, or degree two",
        )
        require(
            not np.any(anchor_shift.sum())
            and not np.any(anchor_shift % 2)
            and not np.any(kernel_rows @ anchor_shift),
            "anchor representatives differ by more than an exact zero-mean hull vector",
        )
        for modulus in (3, 7):
            require(
                affine_hull.modular_rank(
                    np.vstack((hull_bases[modulus], anchor_shift)), modulus
                )
                == 20,
                f"anchor shift escaped the full mod-{modulus} affine hull",
            )
        nonliteral_anchor_transports += int(np.any(anchor_shift))
        anchor_internal[key] = (
            mapped_anchor,
            target_anchor,
            anchor_shift,
            local_map,
        )
        anchor_records.append(
            {
                "branch": branch,
                "source_direction": source_direction,
                "target_direction": target_direction,
                "source_mask": source_mask,
                "target_mask": target_mask,
                "scaled_mean": mean,
                "catalog_class": catalog_class,
                "occurrences": domain_occurrences[key],
                "mapped_anchor_sha256_int64": array_sha256(mapped_anchor),
                "target_anchor_sha256_int64": array_sha256(target_anchor),
                "anchor_shift_sha256_int64": array_sha256(anchor_shift),
                "anchors_literally_equal": not np.any(anchor_shift),
                "anchor_shift_is_exact_even_zero_mean_degree_two": True,
                "anchor_shift_lies_in_full_hull_mod3_mod7": True,
            }
        )

    catalog_records = []
    exact_catalog_occurrences = 0
    required_sm_occurrences = 0
    for key in sorted(domain_occurrences):
        catalog_class = key[-1]
        if catalog_class not in ("U", "S", "M"):
            continue
        (
            branch,
            source_direction,
            target_direction,
            source_mask,
            target_mask,
            mean,
            _catalog_class,
        ) = key
        mapped_anchor, target_anchor, anchor_shift, local_map = anchor_internal[key]
        source_catalog = affine_hull.mapped_catalog(source_mask, mean).astype(np.int64)
        mapped_source_catalog = transport_johnson_values(source_catalog, local_map)
        target_catalog = affine_hull.mapped_catalog(target_mask, mean).astype(np.int64)
        target_lookup = {row.tobytes(): index for index, row in enumerate(target_catalog)}
        require(
            len(target_lookup) == len(target_catalog),
            "target complete catalog unexpectedly repeats an exact row",
        )
        row_map = np.asarray(
            [target_lookup.get(row.tobytes(), -1) for row in mapped_source_catalog],
            dtype=np.int64,
        )
        require(
            np.all(row_map >= 0)
            and len(row_map) == len(target_catalog)
            and len(set(int(value) for value in row_map)) == len(target_catalog),
            "affine transport is not a bijection of complete exact catalog rows",
        )
        matched_target_catalog = target_catalog[row_map]
        require(
            np.array_equal(mapped_source_catalog, matched_target_catalog),
            "exact transformed catalog row set differs from target",
        )
        for modulus in (3, 7):
            require(
                np.array_equal(
                    mapped_source_catalog % modulus,
                    target_catalog[row_map] % modulus,
                ),
                f"shared exact row map failed modulo {modulus}",
            )

        mapped_source_delta = mapped_anchor[None, :] - mapped_source_catalog
        target_delta = target_anchor[None, :] - matched_target_catalog
        require(
            np.array_equal(
                mapped_source_delta - target_delta,
                np.broadcast_to(anchor_shift, mapped_source_delta.shape),
            ),
            "anchor/catalog contribution shift is not row-independent",
        )
        mapped_anchor_base = 13 - mapped_anchor
        target_anchor_base = 13 - target_anchor
        require(
            np.array_equal(
                (mapped_anchor_base - target_anchor_base)[None, :]
                + (mapped_source_delta - target_delta),
                np.zeros_like(mapped_source_delta),
            ),
            "anchor base-target shift did not cancel catalog-contribution shift",
        )
        expected_size = {"U": 1, "S": 56, "M": 1_764}[catalog_class]
        require(len(source_catalog) == expected_size, "catalog class size changed")
        occurrences = int(domain_occurrences[key])
        exact_catalog_occurrences += occurrences
        if catalog_class in ("S", "M"):
            required_sm_occurrences += occurrences
        catalog_records.append(
            {
                "branch": branch,
                "source_direction": source_direction,
                "target_direction": target_direction,
                "source_mask": source_mask,
                "target_mask": target_mask,
                "scaled_mean": mean,
                "catalog_class": catalog_class,
                "occurrences": occurrences,
                "complete_catalog_rows": len(source_catalog),
                "source_exact_catalog_sha256_int16": array_sha256(
                    source_catalog.astype(np.int16)
                ),
                "mapped_source_exact_catalog_sha256_int16": array_sha256(
                    mapped_source_catalog.astype(np.int16)
                ),
                "target_exact_catalog_sha256_int16": array_sha256(
                    target_catalog.astype(np.int16)
                ),
                "shared_mod3_mod7_row_index_permutation_sha256_int64": array_sha256(
                    row_map
                ),
                "exact_row_set_transport_bijective": True,
                "one_row_index_map_shared_by_mod3_and_mod7": True,
                "anchor_base_and_catalog_contribution_shifts_cancel_exactly": True,
            }
        )
    require(required_sm_occurrences > 0, "survivor classes contain no S/M catalog occurrence")
    catalogs_by_branch_identity: dict[tuple, list[dict]] = {}
    for row in catalog_records:
        identity_key = (
            row["source_direction"],
            row["target_direction"],
            row["source_mask"],
            row["target_mask"],
            row["scaled_mean"],
            row["catalog_class"],
        )
        catalogs_by_branch_identity.setdefault(identity_key, []).append(row)
    for rows in catalogs_by_branch_identity.values():
        require(
            {row["branch"] for row in rows} == {"A", "B"}
            and len({row["source_exact_catalog_sha256_int16"] for row in rows}) == 1
            and len({row["target_exact_catalog_sha256_int16"] for row in rows}) == 1,
            "A/B identity did not use the same exact mask/mean catalog domains",
        )

    high_keys = [key for key in domain_occurrences if key[-1] == "H"]
    require(high_keys, "survivor classes contain no high affine-hull occurrence")
    for key in high_keys:
        _mapped_anchor, _target_anchor, anchor_shift, _local_map = anchor_internal[key]
        for modulus in (3, 7):
            require(
                affine_hull.modular_rank(
                    np.vstack((hull_bases[modulus], anchor_shift)), modulus
                )
                == 20,
                "high anchor coset transport failed",
            )

    for branch in BRANCHES:
        for modulus in (3, 7):
            source_rows = augmented[(0, branch, modulus)]
            transported = transported_columns(source_rows, coordinate_maps[branch])
            transported_rref, _ = canonical_rref(transported, modulus)
            require(
                np.array_equal(transported_rref, canonical[(1, branch, modulus)]),
                "global-transfer compact augmented base target changed",
            )
    require(
        all(
            np.array_equal(canonical[(orbit, "A", modulus)], canonical[(orbit, "B", modulus)])
            for orbit in range(2)
            for modulus in (3, 7)
        ),
        "A/B compact identity failed for global transfer",
    )

    equivalence_classes_sha256 = json_sha256(canonical_classes)
    parent_survivor_hash = survivor_audit["input_evidence_audit"][
        "survivor_case_keys_sha256"
    ]
    require(
        isinstance(parent_survivor_hash, str) and len(parent_survivor_hash) == 64,
        "audited parent survivor certificate is missing",
    )
    return {
        "scope": "positive_p7_z7_1296_affine_hull_survivors",
        "parent_survivor_case_keys_sha256": parent_survivor_hash,
        "equivalence_class_count": len(equivalence_classes),
        "equivalence_class_size": 4,
        "equivalence_classes": equivalence_classes,
        "equivalence_classes_sha256": equivalence_classes_sha256,
        "catalog_row_transport_audit": {
            "audited_primes": [3, 7],
            "unique_exact_U_S_M_transport_domains": len(catalog_records),
            "exact_U_S_M_direction_occurrences_in_branch_maps": exact_catalog_occurrences,
            "required_S_M_direction_occurrences_in_branch_maps": required_sm_occurrences,
            "catalog_class_occurrence_histogram": dict(
                sorted(class_occurrence_histogram.items())
            ),
            "domain_audits": catalog_records,
            "A_B_identity_uses_identical_exact_mask_mean_catalog_matrices": True,
            "all_complete_S_M_catalogs_transported_as_exact_row_sets": True,
            "one_exact_row_bijection_per_domain_used_for_both_mod3_and_mod7": True,
        },
        "high_affine_hull_and_anchor_target_audit": {
            "full_exact_zero_mean_degree_two_hull_dimension": 20,
            "hull_transport_audits": hull_transport_rows,
            "unique_anchor_transport_domains": len(anchor_records),
            "anchor_transport_domain_audits_sha256": json_sha256(anchor_records),
            "nonliteral_anchor_transports_explicitly_audited": nonliteral_anchor_transports,
            "high_transport_domains": len(high_keys),
            "high_direction_occurrences_in_branch_maps": sum(
                domain_occurrences[key] for key in high_keys
            ),
            "all_mapped_anchors_have_exact_target_parity_mean_and_degree_two": True,
            "every_anchor_difference_explicitly_audited_as_a_hull_vector_mod3_mod7": True,
            "mapped_high_anchor_plus_full_hull_equals_target_anchor_plus_full_hull": True,
            "retained_catalog_anchor_shifts_cancel_between_base_and_contribution": True,
        },
        "compact_augmented_equation_transport_audit": {
            "audited_primes": [3, 7],
            "all_branch_fixed_edge_rhs_and_constants_included": True,
            "A_B_identity_preserves_augmented_row_spaces": True,
            "branch_affine_maps_preserve_augmented_row_spaces_and_base_targets": True,
        },
        "proof_composition": (
            "A/B identity preserves the compact augmented systems.  Each branch affine map "
            "bijects exact S/M catalog rows with one row permutation shared by mod 3 and mod 7; "
            "it preserves the full high zero-mean degree-two hull.  A transformed anchor may "
            "differ from the target anchor by an explicitly audited hull vector; for retained "
            "catalogs that shift cancels exactly between the anchor base target and the "
            "anchor-minus-catalog contribution, while high conditioning annihilates it."
        ),
        "all_1296_survivors_partitioned_into_324_disjoint_complete_classes": True,
        "no_partial_survivor_classes": True,
        "all_class_maps_preserve_exact_catalog_row_identity_mod3_mod7": True,
        "all_class_maps_preserve_high_affine_hulls_and_anchor_targets": True,
        "transfer_valid_for_global_same_index_catalog_join": True,
    }


def run(input_path: Path | None) -> dict:
    started = time.time()
    context = construct_pointed_systems()
    row_space_audit, augmented, canonical = audit_compact_row_spaces(context)

    directions, types, labels = orbit_helpers.direction_metadata()
    require(
        tuple(directions) == tuple([(1, value) for value in range(P)] + [(0, 1)]),
        "direction order changed",
    )
    require(tuple(types) == tuple(parent.DIRECTION_TYPES), "direction types changed")
    require(tuple(labels) == tuple(parent.LABELS), "fibre labels changed")
    require(tuple(SIGMA) == (4, 7, 1, 0, 2, 3, 5, 6), "Pauli sigma changed")
    evaluation, kernel, kernel_audit = feature_and_kernel_data()
    map_rows = []
    coordinate_maps = {}
    for branch in BRANCHES:
        map_audit, coordinate_map = affine_map_audit(
            branch,
            context,
            augmented,
            canonical,
            directions,
            types,
            labels,
            evaluation,
            kernel,
        )
        map_rows.append(map_audit)
        coordinate_maps[branch] = coordinate_map
    require(
        coordinate_maps["A"] != coordinate_maps["B"],
        "the two branch-specific fibre transports unexpectedly coincided",
    )

    leaves_by_orbit, leaf_coverage = parent.exact_mean_leaves(context["orbits"])
    leaf_transport, leaf_permutation, classes = audit_leaf_transport(
        context["orbits"], leaves_by_orbit
    )
    survivor_audit, survivors, complete_source_indices = optional_survivor_audit(
        input_path, context, leaves_by_orbit, leaf_coverage, classes
    )
    survivor_complete = bool(survivor_audit["performed"])
    if survivor_complete:
        require(survivors is not None, "audited survivor rows disappeared")
        global_catalog_join_transfer = audit_global_catalog_join_transfer(
            context=context,
            leaves_by_orbit=leaves_by_orbit,
            leaf_permutation=leaf_permutation,
            classes=classes,
            survivors=survivors,
            complete_source_indices=complete_source_indices,
            survivor_audit=survivor_audit,
            coordinate_maps=coordinate_maps,
            augmented=augmented,
            canonical=canonical,
        )
    else:
        global_catalog_join_transfer = unavailable_global_catalog_join_transfer()
    status = (
        "complete_exact_four_case_compact_symmetry_and_survivor_partition_audit"
        if survivor_complete
        else "complete_exact_four_case_compact_symmetry_audit_without_optional_survivor_input"
    )
    return {
        "experiment": "p7_infinity7_positive_z7_compact_symmetry_audit",
        "status": status,
        "p": P,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "moduli": list(MODULI),
        "pointed_system_construction": {
            "repository_helpers": [
                "pointed_full_cpsat.audited_pointed_cases",
                "z7_mod7_projection.load_z7_orbits",
                "pointed_affine_hull_multimod.build_pointed_systems",
                "p7_unsaturated_modular_catalog_filter.left_dependencies",
            ],
            "orbit_source": context["orbit_source"],
            "translation_equivariant_matrix_shape": list(context["translation"].shape),
            "translation_equivariant_matrix_sha256_int16": affine_hull.matrix_sha256(
                context["translation"]
            ),
            "pointed_case_count": 4,
            "branches": list(BRANCHES),
            "all_four_exact_pointed_systems_reconstructed_and_cross_audited": True,
        },
        "canonical_compact_dependency_audit": row_space_audit,
        "johnson_degree_two_audit": kernel_audit,
        "affine_compact_isomorphism_audit": {
            "maps": map_rows,
            "direction_permutation_sigma": list(SIGMA),
            "identity_maps_A_to_B_within_each_orbit_after_compact_dependency_projection": True,
            "branch_A_map_connects_orbit0_A_to_orbit1_A": True,
            "branch_B_map_connects_orbit0_B_to_orbit1_B": True,
            "all_four_compact_cases_are_exactly_isomorphic": True,
        },
        "leaf_transport_audit": leaf_transport,
        "mean_leaf_coverage": leaf_coverage,
        "optional_affine_survivor_partition_audit": survivor_audit,
        "global_catalog_join_transfer": global_catalog_join_transfer,
        "exact_transfer_scope": {
            "compact_infeasibility_transfers_across_each_four_case_class": True,
            "compact_feasibility_transfers_across_each_four_case_class": True,
            "global_same_index_catalog_join_transfer_with_audited_input": survivor_complete,
            "compact_feasible_witness_semantics": "necessary-only compact relaxation witness",
            "raw_A_or_B_binary_edge_witness_transfers": False,
            "reason_raw_edge_witness_does_not_transfer": (
                "the certificate identifies projected compact dependency systems, not the raw "
                "branch-specific binary edge-variable systems"
            ),
            "unknown_transfers_any_decision": False,
            "unknown_semantics": "UNKNOWN transfers nothing",
        },
        "solver_invoked": False,
        "all_required_structural_audits_passed": True,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        help="optional canonical full pointed affine-hull survivor JSON",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.input)
    pointed.atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "experiment": result["experiment"],
                "status": result["status"],
                "output": str(args.output),
                "common_ranks": EXPECTED_COMMON_RANKS,
                "leaf_permutation_sha256": EXPECTED_LEAF_PERMUTATION_SHA256,
                "survivor_partition_audited": result[
                    "optional_affine_survivor_partition_audit"
                ]["performed"],
                "global_same_index_catalog_join_transfer": result[
                    "global_catalog_join_transfer"
                ]["transfer_valid_for_global_same_index_catalog_join"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
