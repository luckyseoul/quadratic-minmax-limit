#!/usr/bin/env python3
"""Compact pointed mod-3/5/7/11 CP-SAT relaxation for positive p=7, z=7.

The four audited pointed line cases are reused without enumerating mean
leaves.  For each selected orbit/branch this model has 280 Johnson-slice
slack variables with exact phase-zero parity, all primitive degree-two
kernel equations, exact directional means and type budgets, and common
per-type residues restricted to zero or four modulo eight.  Every complete
left dependency of the exact pointed edge system is imposed modulo 3, 5, 7,
and 11, always against the same 280 slack variables.

This is an edge relaxation.  CP-SAT ``INFEASIBLE`` rigorously excludes the
selected pointed case.  ``FEASIBLE`` supplies only a directly audited
catalog/right-side witness and is necessary-only, never an edge witness.
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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_global_modular_cpsat as global_compact  # noqa: E402
import p7_infinity7_positive_z7_mod7_projection as z7_parent  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as pointed_systems  # noqa: E402
import p7_infinity7_positive_z7_pointed_full_cpsat as pointed_cases  # noqa: E402
import p7_infinity7_positive_z7_pointed_mod7 as pointed  # noqa: E402


P = 7
Q = P * P
EDGE_COUNT = 4 * P + 1
HELPER_MODULI = (3, 7)
EXTRA_MODULI = (5, 11)
MODULI = (3, 5, 7, 11)
EXPECTED_ORBITS = 2
EXPECTED_BRANCHES = ("A", "B")
EXPECTED_CASES = 4
EXPECTED_SLACKS = 8 * 35
EXPECTED_KERNEL_EQUATIONS = 8 * 14
EXPECTED_BRANCH_RANKS = {
    3: {"A": 162, "B": 169},
    5: {"A": 168, "B": 175},
    7: {"A": 147, "B": 154},
    11: {"A": 168, "B": 175},
}
EXPECTED_DEPENDENCY_COUNTS = {3: 120, 5: 114, 7: 135, 11: 114}
EXPECTED_MODULAR_CONSTRAINTS = 483
EXPECTED_VARIABLES = 2 * EXPECTED_SLACKS + 8 + 8 + 2
EXPECTED_CONSTRAINTS = (
    EXPECTED_SLACKS
    + EXPECTED_KERNEL_EQUATIONS
    + 8
    + 2
    + 8
    + EXPECTED_MODULAR_CONSTRAINTS
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def kernel_data() -> tuple[np.ndarray, dict]:
    """Load and independently audit the primitive Johnson kernel rows."""
    kernel = np.ascontiguousarray(
        np.asarray(global_compact._primitive_left_kernel_rows(), dtype=np.int64)  # noqa: SLF001
    )
    require(kernel.shape == (14, 35), "primitive Johnson kernel shape changed")
    row_gcds = [math.gcd(*(abs(int(value)) for value in row)) for row in kernel]
    require(row_gcds == [1] * 14, "Johnson kernel row is not primitive")
    ranks = {
        modulus: pointed_systems.modular_rank(kernel, modulus) for modulus in MODULI
    }
    require(
        ranks == {3: 14, 5: 14, 7: 14, 11: 14},
        "primitive Johnson kernel lost modular rank",
    )
    return kernel, {
        "rows": 14,
        "columns": 35,
        "integer_row_gcds": row_gcds,
        "modular_ranks": {str(modulus): rank for modulus, rank in ranks.items()},
        "matrix_sha256_int64": array_sha256(kernel),
        "all_rows_primitive": True,
    }


def attach_and_audit_four_modulus_dependencies(
    systems_by_orbit: list[dict[str, dict]],
    orbits: list[dict],
) -> list[dict]:
    """Attach direct mod-5/11 bases and audit complete bases at all four primes."""
    require(len(systems_by_orbit) == len(orbits) == EXPECTED_ORBITS, "orbit system count changed")
    audits = []
    for orbit_index, orbit_systems in enumerate(systems_by_orbit):
        require(set(orbit_systems) == set(EXPECTED_BRANCHES), "pointed branch system missing")
        branch_audits = []
        for branch in EXPECTED_BRANCHES:
            system = orbit_systems[branch]
            matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
            calibration = (17 * np.arange(matrix.shape[1], dtype=np.int64) + 3) % 2
            prime_audits = []
            for modulus in MODULI:
                if modulus in EXTRA_MODULI:
                    rank, dependencies = global_compact.left_dependencies(matrix, modulus)
                    dependencies = np.ascontiguousarray(dependencies, dtype=np.int64)
                    system["moduli"][modulus] = {
                        "dependencies": dependencies,
                        "rank": int(rank),
                    }
                    source = "direct_left_dependencies_on_exact_pointed_matrix"
                else:
                    rank = int(system["moduli"][modulus]["rank"])
                    dependencies = np.ascontiguousarray(
                        system["moduli"][modulus]["dependencies"], dtype=np.int64
                    )
                    source = "audited_pointed_affine_helper"

                expected_rank = EXPECTED_BRANCH_RANKS[modulus][branch]
                direct_rank = pointed_systems.modular_rank(matrix, modulus)
                expected_dependencies = EXPECTED_DEPENDENCY_COUNTS[modulus]
                require(
                    rank == direct_rank == expected_rank,
                    f"mod-{modulus} branch {branch} rank changed",
                )
                require(
                    matrix.shape[0] - rank == expected_dependencies,
                    f"mod-{modulus} branch {branch} nullity changed",
                )
                require(
                    dependencies.shape == (expected_dependencies, matrix.shape[0]),
                    f"mod-{modulus} branch {branch} dependency shape changed",
                )
                dependency_basis_rank = pointed_systems.modular_rank(dependencies, modulus)
                require(
                    dependency_basis_rank == expected_dependencies,
                    f"mod-{modulus} branch {branch} dependency basis is not full",
                )
                left_product = dependencies @ (matrix % modulus) % modulus
                require(
                    not np.any(left_product),
                    f"mod-{modulus} branch {branch} left-null audit failed",
                )
                manufactured_rhs = matrix @ calibration % modulus
                manufactured_syndrome = dependencies @ manufactured_rhs % modulus
                require(
                    not np.any(manufactured_syndrome),
                    f"mod-{modulus} branch {branch} manufactured RHS was rejected",
                )
                prime_audits.append(
                    {
                        "modulus": modulus,
                        "dependency_source": source,
                        "rank": rank,
                        "direct_recomputed_matrix_rank": direct_rank,
                        "left_dependency_dimension": len(dependencies),
                        "left_dependency_basis_rank": dependency_basis_rank,
                        "complete_left_dependency_basis": True,
                        "left_null_audit": True,
                        "dependency_sha256_uint8": pointed_systems.matrix_sha256(
                            dependencies.astype(np.uint8)
                        ),
                        "manufactured_binary_edge_vector_sha256_uint8": pointed_systems.matrix_sha256(
                            calibration[None, :].astype(np.uint8)
                        ),
                        "manufactured_rhs_sha256_uint8": pointed_systems.matrix_sha256(
                            manufactured_rhs[None, :].astype(np.uint8)
                        ),
                        "manufactured_witness_syndrome_is_zero": True,
                    }
                )
            branch_audits.append(
                {
                    "pointed_star_branch": branch,
                    "equations": int(matrix.shape[0]),
                    "edge_variables": int(matrix.shape[1]),
                    "fixed_edge_rows": system["fixed_edge_rows"],
                    "matrix_sha256_int16": pointed_systems.matrix_sha256(system["matrix"]),
                    "base_rhs_sha256_int64": pointed_systems.matrix_sha256(
                        np.asarray(system["base_rhs"], dtype=np.int64)[None, :]
                    ),
                    "prime_audits": prime_audits,
                }
            )
        audits.append(
            {
                "branch_orbit_index": orbit_index,
                "source_orbit_index": int(orbits[orbit_index]["source_orbit_index"]),
                "systems": branch_audits,
            }
        )
    return audits


def audited_inputs() -> tuple[list[dict], list[dict[str, dict]], dict]:
    """Rebuild the four audited cases and their complete pointed systems."""
    require(P == pointed.P == pointed_cases.P == pointed_systems.P, "p=7 constants changed")
    require(Q == pointed.Q == pointed_cases.Q == pointed_systems.Q, "GF(49) constants changed")
    require(HELPER_MODULI == pointed_systems.MODULI, "pointed helper moduli changed")
    require(EDGE_COUNT == pointed.EDGE_COUNT == pointed_systems.EDGE_COUNT, "edge count changed")
    require(tuple(global_compact.POINTS) == tuple(pointed_systems.POINTS), "Johnson point order changed")

    cases, orbit_source, normalization_audit = pointed_cases.audited_pointed_cases()
    require(len(cases) == EXPECTED_CASES, "audited pointed case count changed")
    orbits, second_orbit_source = z7_parent.load_z7_orbits()
    require(len(orbits) == EXPECTED_ORBITS, "z=7 orbit count changed")
    require(
        json_sha256(orbit_source) == json_sha256(second_orbit_source),
        "pointed case and system orbit sources disagree",
    )

    normalizations = []
    for orbit_index, orbit in enumerate(orbits):
        line = tuple(int(value) for value in orbit["representative"])
        rows = [row for row in cases if int(row["branch_orbit_index"]) == orbit_index]
        require({row["pointed_star_branch"] for row in rows} == {"A", "B"}, "orbit branch missing")
        branch_b = next(row for row in rows if row["pointed_star_branch"] == "B")
        outside_present = [
            int(row["finite_field_point"])
            for row in branch_b["fixed_infinity_star_edges"]
            if int(row["value"]) == 1 and int(row["finite_field_point"]) not in line
        ]
        require(len(outside_present) == 1, "branch B outside representative changed")
        normalizations.append(
            {"line": line, "outside_representative": outside_present[0]}
        )

    translation_matrix, systems_by_orbit, translation_audit, system_audits = (
        pointed_systems.build_pointed_systems(orbits, normalizations)
    )
    require(translation_matrix.shape == (281, 1_225), "translation system shape changed")
    require(len(systems_by_orbit) == EXPECTED_ORBITS, "pointed system orbit count changed")
    four_modulus_system_audits = attach_and_audit_four_modulus_dependencies(
        systems_by_orbit, orbits
    )

    case_keys = set()
    for case in cases:
        orbit_index = int(case["branch_orbit_index"])
        branch = str(case["pointed_star_branch"])
        key = (orbit_index, branch)
        require(key not in case_keys, "pointed case key repeated")
        case_keys.add(key)
        require(branch in EXPECTED_BRANCHES, "unknown pointed branch")
        system = systems_by_orbit[orbit_index][branch]
        expected_fixed = [
            (tuple(int(value) for value in row["edge"]), int(row["value"]))
            for row in case["fixed_infinity_star_edges"]
        ]
        observed_fixed = [
            (tuple(int(value) for value in row["graph_edge"]), int(row["rhs"]))
            for row in system["fixed_edge_rows"]
        ]
        require(observed_fixed == expected_fixed, "case and pointed system fixed rows disagree")
        require(np.all(system["base_rhs"][1:281] == 0), "pointed base RHS polluted direction blocks")
        require(int(system["base_rhs"][0]) == EDGE_COUNT, "pointed edge-count RHS changed")
        for modulus in MODULI:
            dependencies = np.asarray(system["moduli"][modulus]["dependencies"])
            expected_rank = EXPECTED_BRANCH_RANKS[modulus][branch]
            require(
                int(system["moduli"][modulus]["rank"]) == expected_rank,
                f"mod-{modulus} branch rank changed",
            )
            require(
                dependencies.shape
                == (EXPECTED_DEPENDENCY_COUNTS[modulus], system["matrix"].shape[0]),
                f"mod-{modulus} dependency shape changed",
            )
            require(
                not np.any(
                    dependencies
                    @ (np.asarray(system["matrix"], dtype=np.int64) % modulus)
                    % modulus
                ),
                f"mod-{modulus} pointed left-null audit failed",
            )
    require(case_keys == {(orbit, branch) for orbit in range(2) for branch in EXPECTED_BRANCHES}, "four-case coverage failed")

    return cases, systems_by_orbit, {
        "orbit_source": orbit_source,
        "pointed_normalization_audit": normalization_audit,
        "translation_system_audit": translation_audit,
        "pointed_system_audits": system_audits,
        "four_modulus_pointed_system_audits": four_modulus_system_audits,
        "all_four_case_systems_cross_audited": True,
        "matrix_rank_nullspace_and_manufactured_rhs_audits_reused": True,
        "affine_helper_mod3_mod7_audit_payload_preserved": True,
        "mod5_mod11_dependencies_computed_directly_on_each_exact_pointed_matrix": True,
        "all_four_moduli_left_null_full_basis_and_manufactured_rhs_audited": True,
        "mean_leaf_enumeration_used": False,
    }


def direction_specs(case: dict) -> list[dict]:
    """Build exact fixed-boundary phase-zero parity and floor metadata."""
    directions = tuple(global_compact.projective_directions(P))
    require(len(directions) == P + 1, "projective direction count changed")
    specs = []
    for direction_index, direction in enumerate(directions):
        eps, labels = global_compact.field_direction_data(P, direction)
        mask = int(case["direction_masks"][direction_index])
        b = mask.bit_count()
        require(b == int(case["b_values"][direction_index]), "stored b-value changed")
        fibres = frozenset(fibre for fibre in range(P) if mask & (1 << fibre))
        parity = tuple(
            sum(value in fibres for value in point) & 1
            for point in global_compact.POINTS
        )
        floor = int(pointed_cases.scaled_direction_floor(P, b, 0))
        specs.append(
            {
                "direction_index": direction_index,
                "direction": tuple(int(value) for value in direction),
                "eps": int(eps),
                "mask": mask,
                "b": b,
                "phase": 0,
                "floor": floor,
                "parity": parity,
                "labels_sha256_int16": array_sha256(np.asarray(labels, dtype=np.int16)),
            }
        )
    require(Counter(row["eps"] for row in specs) == Counter({-1: 4, 1: 4}), "direction types changed")
    require(sorted(row["b"] for row in specs) == [1] + [7] * 7, "z=7 b-profile changed")
    require(sorted(row["floor"] for row in specs) == [0] * 7 + [8], "phase-zero floor profile changed")
    require(
        Counter(bit for row in specs for bit in row["parity"]) == Counter({0: 260, 1: 20}),
        "fixed z=7 parity census changed",
    )
    return specs


def rendered_direction_specs(specs: list[dict]) -> list[dict]:
    rows = []
    for spec in specs:
        parity = np.asarray(spec["parity"], dtype=np.uint8)
        rows.append(
            {
                key: value
                for key, value in spec.items()
                if key != "parity"
            }
            | {
                "parity_zero_coordinates": int(np.count_nonzero(parity == 0)),
                "parity_one_coordinates": int(np.count_nonzero(parity == 1)),
                "parity_sha256_uint8": array_sha256(parity),
                "possible_scaled_means_before_common_residue_coupling": list(
                    range(int(spec["floor"]), 33, 4)
                ),
            }
        )
    return rows


def build_model(case: dict, system: dict, kernel: np.ndarray) -> tuple[object, dict, dict]:
    """Construct one compact pointed relaxation without invoking a solver."""
    from ortools.sat.python import cp_model

    orbit_index = int(case["branch_orbit_index"])
    branch = str(case["pointed_star_branch"])
    require(system["orbit_index"] == orbit_index and system["branch"] == branch, "wrong pointed system")
    specs = direction_specs(case)
    model = cp_model.CpModel()

    slacks = []
    lifts = []
    means = []
    quotients = []
    means_by_type: dict[int, list] = {-1: [], 1: []}
    for spec in specs:
        direction_index = int(spec["direction_index"])
        values = []
        direction_lifts = []
        for point_index, parity in enumerate(spec["parity"]):
            lift = model.new_int_var(
                0,
                (13 - int(parity)) // 2,
                f"lift_d{direction_index}_x{point_index}",
            )
            value = model.new_int_var(0, 13, f"slack_d{direction_index}_x{point_index}")
            model.add(value == 2 * lift + int(parity))
            values.append(value)
            direction_lifts.append(lift)
        for kernel_index, row in enumerate(kernel):
            model.add(
                sum(int(row[index]) * values[index] for index in range(35)) == 0
            )

        mean = model.new_int_var(
            int(spec["floor"]),
            32,
            f"scaled_mean_d{direction_index}",
        )
        model.add(2 * sum(values) == 5 * mean)
        slacks.append(values)
        lifts.append(direction_lifts)
        means.append(mean)
        means_by_type[int(spec["eps"])].append((direction_index, mean))

    common_residues = {}
    for eps in (-1, 1):
        typed = means_by_type[eps]
        require(len(typed) == 4, "direction type is not four-by-four")
        model.add(sum(mean for _direction, mean in typed) == 32)
        residue = model.new_int_var_from_domain(
            cp_model.Domain.from_values([0, 4]),
            f"common_residue_{eps}",
        )
        common_residues[eps] = residue
        for direction_index, mean in typed:
            quotient = model.new_int_var(0, 4, f"mean_quotient_d{direction_index}")
            model.add(mean == residue + 8 * quotient)
            quotients.append((direction_index, quotient))

    modular_counts = {}
    base_rhs = np.asarray(system["base_rhs"], dtype=np.int64)
    require(base_rhs.shape == (system["matrix"].shape[0],), "pointed base RHS shape changed")
    for modulus in MODULI:
        dependencies = np.asarray(system["moduli"][modulus]["dependencies"], dtype=np.int64)
        require(dependencies.shape[1] == len(base_rhs), "dependency/base RHS width mismatch")
        for dependency in dependencies:
            constant = int(dependency @ base_rhs)
            terms = []
            for direction_index in range(P + 1):
                block = dependency[
                    1 + 35 * direction_index : 1 + 35 * (direction_index + 1)
                ]
                constant += 13 * int(np.sum(block, dtype=np.int64))
                terms.extend(
                    -int(block[point_index]) * slacks[direction_index][point_index]
                    for point_index in range(35)
                    if int(block[point_index])
                )
            model.add_modulo_equality(0, constant + sum(terms), modulus)
        modular_counts[modulus] = len(dependencies)
    require(
        modular_counts == EXPECTED_DEPENDENCY_COUNTS,
        "complete dependency constraint census changed",
    )
    require(
        sum(modular_counts.values()) == EXPECTED_MODULAR_CONSTRAINTS,
        "total modular dependency constraint census changed",
    )

    validation_error = model.validate()
    require(not validation_error, f"constructed CP-SAT model is invalid: {validation_error}")
    proto = model.Proto()
    require(len(proto.variables) == EXPECTED_VARIABLES, "compact variable count changed")
    require(len(proto.constraints) == EXPECTED_CONSTRAINTS, "compact constraint count changed")
    stats = model.model_stats()
    construction = {
        "case_key": f"orbit{orbit_index}_{branch}",
        "slack_intvars_0_to_13": EXPECTED_SLACKS,
        "parity_lift_intvars": EXPECTED_SLACKS,
        "scaled_mean_intvars": len(means),
        "mean_quotient_intvars": len(quotients),
        "common_residue_intvars_domain_0_or_4": len(common_residues),
        "primitive_kernel_equations": EXPECTED_KERNEL_EQUATIONS,
        "directional_mean_identities": len(means),
        "four_per_type_sum_32_equations": 2,
        "common_residue_equations": len(quotients),
        "modular_dependency_constraints": {
            str(modulus): count for modulus, count in modular_counts.items()
        },
        "modular_dependency_constraints_total": sum(modular_counts.values()),
        "same_280_slack_variables_coupled_across_all_moduli": True,
        "total_model_variables": len(proto.variables),
        "total_model_constraints": len(proto.constraints),
        "model_validation": "passed",
        "model_stats_sha256": hashlib.sha256(stats.encode("utf-8")).hexdigest(),
        "direction_rows": rendered_direction_specs(specs),
        "mean_leaf_enumeration_used": False,
        "every_left_dependency_imposed": True,
        "bad_count_block_formula": "13 - slack",
        "edge_count_and_fixed_edge_values_enter_through_exact_base_rhs": True,
    }
    variables = {
        "slacks": slacks,
        "lifts": lifts,
        "means": means,
        "quotients": dict(quotients),
        "common_residues": common_residues,
        "specs": specs,
    }
    return model, variables, construction


def direct_witness_audit(
    solver: object,
    case: dict,
    system: dict,
    kernel: np.ndarray,
    variables: dict,
) -> tuple[dict, dict]:
    """Recompute every relaxation identity from a feasible assignment."""
    slack_values = np.asarray(
        [
            [int(solver.value(variable)) for variable in direction]
            for direction in variables["slacks"]
        ],
        dtype=np.int64,
    )
    means = [int(solver.value(variable)) for variable in variables["means"]]
    residues = {
        eps: int(solver.value(variable))
        for eps, variable in variables["common_residues"].items()
    }
    quotients = {
        direction: int(solver.value(variable))
        for direction, variable in variables["quotients"].items()
    }
    specs = variables["specs"]

    shape_ok = slack_values.shape == (8, 35)
    bounds_ok = bool(np.all((0 <= slack_values) & (slack_values <= 13)))
    parity_ok = all(
        np.array_equal(
            slack_values[direction] % 2,
            np.asarray(specs[direction]["parity"], dtype=np.int64),
        )
        for direction in range(8)
    )
    kernel_ok = all(
        not np.any(kernel @ slack_values[direction]) for direction in range(8)
    )
    mean_identity_ok = all(
        2 * int(slack_values[direction].sum()) == 5 * means[direction]
        for direction in range(8)
    )
    floor_ok = all(
        int(specs[direction]["floor"]) <= means[direction] <= 32
        for direction in range(8)
    )
    quotient_ok = all(
        means[direction] == residues[int(specs[direction]["eps"])] + 8 * quotients[direction]
        for direction in range(8)
    )

    type_rows = {}
    type_ok = True
    for eps in (-1, 1):
        indices = [index for index, spec in enumerate(specs) if int(spec["eps"]) == eps]
        values = [means[index] for index in indices]
        observed_residues = [value % 8 for value in values]
        valid = bool(
            len(indices) == 4
            and sum(values) == 32
            and residues[eps] in (0, 4)
            and observed_residues == [residues[eps]] * 4
        )
        type_ok &= valid
        type_rows[str(eps)] = {
            "direction_indices": indices,
            "scaled_means": values,
            "sum": sum(values),
            "common_residue": residues[eps],
            "observed_residues_mod_8": observed_residues,
            "valid": valid,
        }

    rhs = np.asarray(system["base_rhs"], dtype=np.int64).copy()
    for direction in range(8):
        rhs[1 + 35 * direction : 1 + 35 * (direction + 1)] = 13 - slack_values[direction]
    base_rhs = np.asarray(system["base_rhs"], dtype=np.int64)
    base_ok = bool(
        rhs[0] == EDGE_COUNT
        and np.array_equal(rhs[281:], base_rhs[281:])
        and np.array_equal(base_rhs[1:281], np.zeros(280, dtype=np.int64))
    )
    syndrome_rows = []
    syndromes_ok = True
    for modulus in MODULI:
        dependencies = np.asarray(system["moduli"][modulus]["dependencies"], dtype=np.int64)
        syndrome = dependencies @ (rhs % modulus) % modulus
        valid = not np.any(syndrome)
        syndromes_ok &= bool(valid)
        syndrome_rows.append(
            {
                "modulus": modulus,
                "dependency_rows": len(dependencies),
                "nonzero_syndrome_coordinates": int(np.count_nonzero(syndrome)),
                "syndrome_sha256_uint8": array_sha256(syndrome.astype(np.uint8)),
                "valid": bool(valid),
            }
        )

    checks = {
        "slack_shape_8_by_35": shape_ok,
        "slack_bounds_0_to_13": bounds_ok,
        "exact_phase_zero_parities": parity_ok,
        "all_14_primitive_kernel_equations_per_direction": kernel_ok,
        "all_exact_directional_mean_identities": mean_identity_ok,
        "all_exact_phase_zero_floors": floor_ok,
        "four_per_type_sum_32_and_common_residue_0_or_4": type_ok,
        "mean_residue_quotient_identities": quotient_ok,
        "edge_count_and_fixed_edge_base_rhs": base_ok,
        "all_mod3_mod5_mod7_mod11_dependency_syndromes_zero": syndromes_ok,
    }
    audit = {
        "valid": all(checks.values()),
        "checks": checks,
        "type_rows": type_rows,
        "modular_syndromes": syndrome_rows,
        "slack_sha256_int64": array_sha256(slack_values),
        "pointed_rhs_sha256_int64": array_sha256(rhs),
        "pointed_rhs_length": len(rhs),
        "fixed_edge_rows": system["fixed_edge_rows"],
        "direct_recomputation_used": True,
    }
    witness = {
        "necessary_only_not_an_edge_witness": True,
        "scaled_means_direction_order": means,
        "common_residues_by_type": {
            str(eps): value for eps, value in residues.items()
        },
        "slack_values": slack_values.tolist(),
        "direct_witness_audit": audit,
    }
    require(audit["valid"], "feasible compact witness failed direct audit")
    return witness, audit


def process_case(
    case: dict,
    system: dict,
    kernel: np.ndarray,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool,
) -> dict:
    """Construct, and unless requested otherwise solve, one pointed case."""
    from ortools.sat.python import cp_model

    started = time.time()
    require(timeout > 0, "timeout must be positive")
    require(workers > 0, "workers must be positive")
    model, variables, construction = build_model(case, system, kernel)
    orbit_index = int(case["branch_orbit_index"])
    branch = str(case["pointed_star_branch"])
    case_key = f"orbit{orbit_index}_{branch}"
    base = {
        "case_key": case_key,
        "branch_orbit_index": orbit_index,
        "source_orbit_index": int(case["source_orbit_index"]),
        "orbit_size": int(case["orbit_size"]),
        "representative_finite_field": list(case["representative_finite_field"]),
        "pointed_star_branch": branch,
        "fixed_infinity_star_edges": list(case["fixed_infinity_star_edges"]),
        "pointed_matrix_shape": list(system["matrix"].shape),
        "pointed_matrix_sha256_int16": pointed_systems.matrix_sha256(system["matrix"]),
        "pointed_base_rhs_sha256_int64": pointed_systems.matrix_sha256(
            np.asarray(system["base_rhs"], dtype=np.int64)[None, :]
        ),
        "construction": construction,
        "timeout_seconds": timeout,
        "workers": workers,
        "seed": seed,
    }
    if construction_only:
        return {
            **base,
            "solver_status": "NOT_RUN_CONSTRUCTION_ONLY",
            "rigorous_status": "model_constructed_and_audited_no_solver_decision",
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
    status_name = solver.status_name(status)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    infeasible = status == cp_model.INFEASIBLE
    require(status != cp_model.MODEL_INVALID, "validated compact model was rejected as invalid")
    if infeasible:
        rigorous_status = "rigorously_infeasible_compact_edge_relaxation"
    elif feasible:
        rigorous_status = "feasible_compact_relaxation_necessary_only"
    else:
        rigorous_status = "solver_unknown_no_rigorous_case_decision"
    result = {
        **base,
        "solver_status": status_name,
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
        witness, audit = direct_witness_audit(
            solver,
            case,
            system,
            kernel,
            variables,
        )
        require(audit["valid"], "direct witness audit failed")
        result["relaxation_witness"] = witness
    return result


def run(
    orbit: str,
    branch: str,
    timeout: float,
    workers: int,
    seed: int,
    construction_only: bool = False,
) -> dict:
    started = time.time()
    kernel, kernel_audit = kernel_data()
    cases, systems_by_orbit, input_audit = audited_inputs()
    selected = pointed_cases.selected_cases(cases, orbit, branch)
    results = []
    for index, case in enumerate(selected):
        orbit_index = int(case["branch_orbit_index"])
        pointed_branch = str(case["pointed_star_branch"])
        results.append(
            process_case(
                case,
                systems_by_orbit[orbit_index][pointed_branch],
                kernel,
                timeout,
                workers,
                seed + index,
                construction_only,
            )
        )

    full_scope = len(selected) == EXPECTED_CASES
    infeasible_count = sum(bool(row["finite_infeasibility_certificate"]) for row in results)
    feasible_count = sum(row["feasible"] is True for row in results)
    unknown_count = sum(
        not row["case_decided"] and row["solver_status"] != "NOT_RUN_CONSTRUCTION_ONLY"
        for row in results
    )
    if construction_only:
        status = "complete_pointed_compact_construction_audit_only"
        conclusion = "all selected compact models constructed and audited; no solver decision requested"
    elif feasible_count:
        status = "pointed_compact_relaxation_with_necessary_survivors"
        conclusion = "at least one selected compact relaxation is feasible; this is not an edge witness"
    elif infeasible_count == len(results) and full_scope:
        status = "complete_rigorous_positive_z7_compact_relaxation_exclusion"
        conclusion = "all four exhaustive pointed edge relaxations are infeasible, excluding positive z=7"
    elif infeasible_count == len(results):
        status = "complete_rigorous_selected_pointed_compact_exclusion"
        conclusion = "every selected pointed edge relaxation is infeasible; unselected cases remain"
    else:
        status = "pointed_compact_relaxation_with_unknown_cases"
        conclusion = "at least one selected solve is unknown; no conclusion follows for that case"

    case_digest_rows = [
        {
            "case_key": row["case_key"],
            "solver_status": row["solver_status"],
            "rigorous_status": row["rigorous_status"],
            "finite_infeasibility_certificate": row["finite_infeasibility_certificate"],
            "feasible": row["feasible"],
            "construction": row["construction"],
            "witness_audit": (
                row.get("relaxation_witness", {}).get("direct_witness_audit")
            ),
        }
        for row in results
    ]
    per_branch = []
    for pointed_branch in EXPECTED_BRANCHES:
        rows = [row for row in results if row["pointed_star_branch"] == pointed_branch]
        if rows:
            per_branch.append(
                {
                    "pointed_star_branch": pointed_branch,
                    "selected_cases": len(rows),
                    "infeasible_cases": sum(bool(row["finite_infeasibility_certificate"]) for row in rows),
                    "necessary_only_feasible_cases": sum(row["feasible"] is True for row in rows),
                    "unknown_or_not_run_cases": sum(not row["case_decided"] for row in rows),
                }
            )

    return {
        "experiment": "p7_infinity7_positive_z7_pointed_compact_cpsat",
        "status": status,
        "rigorous_conclusion": conclusion,
        "p": P,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": P,
        "z": 7,
        "phase": 0,
        "model_semantics": {
            "kind": "compact exact degree-two catalog and pointed edge-right-side relaxation",
            "infeasible": "rigorous exclusion of the selected pointed edge case",
            "feasible": "necessary-only modular catalog/right-side witness; not an edge witness",
            "mean_leaf_enumeration_used": False,
            "moduli": list(MODULI),
            "same_280_slack_variables_coupled_across_all_moduli": True,
            "all_complete_pointed_left_dependency_bases_imposed": True,
        },
        "kernel_audit": kernel_audit,
        **input_audit,
        "selection": {
            "orbit": orbit,
            "branch": branch,
            "selected_case_count": len(selected),
            "full_four_case_scope": full_scope,
        },
        "construction_only": construction_only,
        "timeout_seconds_per_case": timeout,
        "workers_per_case": workers,
        "seed_base": seed,
        "processed_pointed_cases": len(results),
        "infeasible_pointed_cases": infeasible_count,
        "necessary_only_feasible_pointed_cases": feasible_count,
        "unknown_pointed_cases": unknown_count,
        "all_selected_cases_decided": not construction_only and unknown_count == 0,
        "all_selected_cases_infeasible": not construction_only and infeasible_count == len(results),
        "z7_branch_excluded": (
            not construction_only and full_scope and infeasible_count == EXPECTED_CASES
        ),
        "all_case_results_sha256": json_sha256(case_digest_rows),
        "per_branch_summary": per_branch,
        "case_results": results,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orbit",
        choices=("0", "1", "all"),
        default="all",
        help="select one audited z=7 line orbit, or both",
    )
    parser.add_argument(
        "--branch",
        type=str.upper,
        choices=("A", "B", "ALL"),
        default="ALL",
        help="select pointed star branch A or B, or both",
    )
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15719001)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--construction-only",
        action="store_true",
        help="construct and audit selected models without invoking CP-SAT solve",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    result = run(
        args.orbit,
        args.branch,
        args.timeout,
        args.workers,
        args.seed,
        args.construction_only,
    )
    pointed.atomic_write(args.output, result)
    if not args.quiet:
        compact = {
            key: value
            for key, value in result.items()
            if key
            not in {
                "orbit_source",
                "pointed_normalization_audit",
                "translation_system_audit",
                "pointed_system_audits",
                "case_results",
            }
        }
        compact["case_summary"] = [
            {
                "case_key": row["case_key"],
                "solver_status": row["solver_status"],
                "rigorous_status": row["rigorous_status"],
                "model_variables": row["construction"]["total_model_variables"],
                "model_constraints": row["construction"]["total_model_constraints"],
                "elapsed_seconds": row["elapsed_seconds"],
            }
            for row in result["case_results"]
        ]
        print(json.dumps(compact, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
