#!/usr/bin/env python3
"""Exact full-edge CP-SAT attack on the four positive z=7 H0_S0_M7 leaves.

The scope is deliberately narrow: orbit 0, pointed branch A, and exactly the
four reconstructed leaves whose catalog pattern is ``H0_S0_M7``.  The model
is assembled from existing audited constructors:

* ``pointed_compact.build_model`` supplies the integer Johnson-slice slack
  words, exact parity/bounds/kernel/mean constraints, and pointed-system
  provenance;
* the complete exact U/M catalog tables restrict every slack word to one
  actual catalog row at its fixed leaf mean;
* the audited pointed matrix couples those same slack words, over the
  integers, to 1,225 binary edge variables;
* ``pointed_full`` supplies the exact boundary, Paley-product, direction-mean,
  complete eigenshell, and direct witness audits.

Thus SAT is an independently audited exact edge witness, UNSAT rigorously
excludes the selected finite model, and UNKNOWN has no mathematical force.
The only explicit problem-level symmetry normalization is the already-audited
branch-A constraint ``(infinity, 0)=1``.  No additional leaf symmetry is
assumed.  Checkpoints are between cases; CP-SAT itself is not resumable in the
middle of one case.
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

import p7_infinity7_positive_z7_mod7_projection as parent  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine  # noqa: E402
import p7_infinity7_positive_z7_pointed_compact_cpsat as compact  # noqa: E402
import p7_infinity7_positive_z7_pointed_full_cpsat as full  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_h0m7_cpsat"
SCHEMA_VERSION = 1
P = 7
Q = 49
N = 50
EDGE_COUNT = 29
EDGE_VARIABLES = 1_225
TARGET_PATTERN = (0, 0, 7)
EXPECTED_TARGET_LEAF_INDICES = (405, 411, 413, 414)
EXPECTED_TARGET_LEAF_SHA256 = (
    "5359fd6dba9da8729cb444f5b637bea6cb6b48732a5382be5fec8fc98c55973c"
)
EXPECTED_KERNEL_SHA256 = (
    "f90bafe2de158e7a4f08cd32e603631bfa0dab7fdf4fcae2bb013642294d7ccb"
)
EXPECTED_POINTED_MATRIX_SHA256 = (
    "724bb3db3db7972f73d64bffe0e13b352259d2443684ec2374a0091a850d03ba"
)
EXPECTED_POINTED_BASE_RHS_SHA256 = (
    "8660e465871a764a5ee3a45ff7654a1a9e41420c8ba079e01e9dd22029d6da7a"
)
EXPECTED_POINTED_CASE_SHA256 = (
    "a290182cc9aada46ee68c3b33ba3fd64d6ef8df30c1b942c3e3a80d0c04e2678"
)
EXPECTED_ORBIT_SOURCE_SHA256 = (
    "78c29c51313ad434ac26f7a30b60e7e9e2561948cb52df0578cb8143887f25ac"
)
EXPECTED_LEAF_AUDIT_SHA256 = (
    "f2cf33e94777e04dc4fa8efb00dc17d61f164796651f5fb3f2a80233c7c7c6f5"
)

# These are fail-closed provenance pins, not merely values copied to output.
EXPECTED_FILE_SHA256 = {
    "scripts/p7_infinity7_positive_z7_pointed_full_cpsat.py": (
        "2991c7291d1980cef9576eb86fb678563e281673672abae77c4390c36f7859a9"
    ),
    "scripts/p7_infinity7_positive_z7_pointed_compact_cpsat.py": (
        "e4c1de66df2caa5ddaf4ad0e6989d63a6e20a7e20af29a1abc07c3602b10d9ad"
    ),
    "scripts/p7_infinity7_positive_z7_pointed_affine_hull_multimod.py": (
        "0ae46bf3a0ad64975b7e0ac55aea562c9efcdee43828cb205e4e7cabd196d8b8"
    ),
    "scripts/p7_infinity7_positive_z7_mod7_projection.py": (
        "af8ab78d543642a42714e8ba7eff057c93e337fd960c8259230b29d80e0043ac"
    ),
    "scripts/p7_infinity7_positive_z7_pointed_mod7.py": (
        "87deadfedfe4c690275d320190d01759375e6932963916c48f71b4c6cc42c46d"
    ),
    "scripts/p7_infinity7_positive_z2_mod7_join.py": (
        "c5d2968e7ee96d2dbbf938f6ab05d34fbb6e648cc3d797c4511181fbff5ee78f"
    ),
    "scripts/p7_infinity7_global_modular_cpsat.py": (
        "01db07a28ef51f55b8d8d6a0db6a82471566495ad7afc72812cfe73db93e8b07"
    ),
    "scripts/p7_unsaturated_modular_catalog_filter.py": (
        "f9b2781984ab3e2336977d43b657fe337bb09b37baea600b6fdb1f94483d135a"
    ),
    "scripts/residual_boundary_four_lift_cpsat.py": (
        "ea3a561b73c9917a6ec19e3af8821a29630878ace91963ddb96050c9129cfd00"
    ),
    "scripts/p7_size_four_slack_classify.py": (
        "fd262752a3278a7461d4cec3f29a5ee16c66bf81176d99aca1069648b966f217"
    ),
    "evidence/p7_infinity7_positive_zge2_orbits.json": (
        "4cf79ffb19b0e6961f6f187efe513d0fc5503ff2f0c63317603f4b0205d7cd63"
    ),
}
EXPECTED_SHELL_SHA256 = (
    "1d86b18651a4cfdb798bbb743c7686b909804e946a226c1407873c60e0a78427"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def ordered_rows_sha256(rows: list[dict]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(
            json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def positive_float(value: str) -> float:
    parsed = float(value)
    if not np.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive finite number")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def provenance_audit() -> dict:
    rows = []
    for relative, expected in EXPECTED_FILE_SHA256.items():
        path = ROOT / relative
        require(path.is_file(), f"missing pinned provenance file: {relative}")
        observed = file_sha256(path)
        require(observed == expected, f"provenance hash changed: {relative}")
        rows.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": observed,
                "matches_expected": True,
            }
        )

    shell = full.SHELL_PATHS[1]
    require(shell == Path("/tmp/maxplus_p7.npy"), "full-shell cache path changed")
    require(shell.is_file(), f"missing complete full-shell cache: {shell}")
    shell_sha256 = file_sha256(shell)
    require(shell_sha256 == EXPECTED_SHELL_SHA256, "complete full-shell cache hash changed")
    rows.append(
        {
            "path": str(shell),
            "size_bytes": shell.stat().st_size,
            "sha256": shell_sha256,
            "matches_expected": True,
        }
    )
    script_path = Path(__file__).resolve()
    return {
        "fail_closed_hash_pins": True,
        "files": rows,
        "script_path": str(script_path.relative_to(ROOT)),
        "script_sha256": file_sha256(script_path),
        "all_expected_hashes_match": True,
    }


def rendered_leaf(leaf_index: int, leaf: dict) -> dict:
    return {
        "leaf_index": int(leaf_index),
        "residue_pair": str(leaf["residue_pair"]),
        "q_values": [int(value) for value in leaf["q_values"]],
        "scaled_means": [int(value) for value in leaf["scaled_means"]],
        "catalog_levels": [int(value) for value in leaf["catalog_levels"]],
        "catalog_classes": [str(value) for value in leaf["catalog_classes"]],
        "pattern": [int(value) for value in leaf["pattern"]],
    }


def audited_inputs() -> dict:
    """Reconstruct and pin the one pointed case and its four target leaves."""
    provenance = provenance_audit()
    kernel, kernel_audit = compact.kernel_data()
    require(
        kernel_audit["matrix_sha256_int64"] == EXPECTED_KERNEL_SHA256,
        "primitive Johnson kernel certificate changed",
    )

    cases, systems_by_orbit, compact_input_audit = compact.audited_inputs()
    orbits, orbit_source = parent.load_z7_orbits()
    leaves_by_orbit, leaf_audit = parent.exact_mean_leaves(orbits)
    require(len(orbits) == 2, "z=7 orbit census changed")
    require([len(rows) for rows in leaves_by_orbit] == [1_080, 1_080], "leaf census changed")
    require(
        json_sha256(orbit_source) == json_sha256(compact_input_audit["orbit_source"]),
        "compact and leaf orbit provenance disagree",
    )
    require(json_sha256(orbit_source) == EXPECTED_ORBIT_SOURCE_SHA256, "orbit source changed")
    require(json_sha256(leaf_audit) == EXPECTED_LEAF_AUDIT_SHA256, "leaf audit changed")

    target_pairs = [
        (leaf_index, leaf)
        for leaf_index, leaf in enumerate(leaves_by_orbit[0])
        if tuple(int(value) for value in leaf["pattern"]) == TARGET_PATTERN
    ]
    require(
        tuple(index for index, _leaf in target_pairs) == EXPECTED_TARGET_LEAF_INDICES,
        "orbit0 H0_S0_M7 target indices changed",
    )
    require(len(target_pairs) == 4, "target leaf count changed")
    for _leaf_index, leaf in target_pairs:
        require(not leaf["high_directions"], "H0 target unexpectedly has a high direction")
        require(
            Counter(leaf["catalog_classes"]) == Counter({"M": 7, "U": 1}),
            "target is not exactly seven M catalogs and one U catalog",
        )
        require(leaf["residue_pair"] == "00", "H0_S0_M7 residue pair changed")

    leaf_certificate_rows = [rendered_leaf(index, leaf) for index, leaf in target_pairs]
    require(
        json_sha256(leaf_certificate_rows) == EXPECTED_TARGET_LEAF_SHA256,
        "target leaf certificate changed",
    )

    selected_cases = [
        row
        for row in cases
        if int(row["branch_orbit_index"]) == 0
        and str(row["pointed_star_branch"]) == "A"
    ]
    require(len(selected_cases) == 1, "orbit0/branchA pointed case is not unique")
    case = selected_cases[0]
    require(json_sha256(case) == EXPECTED_POINTED_CASE_SHA256, "pointed case changed")
    require(
        case["fixed_infinity_star_edges"]
        == [{"edge": [0, 1], "finite_field_point": 0, "value": 1}],
        "branch-A symmetry normalization changed",
    )
    system = systems_by_orbit[0]["A"]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int16)
    base_rhs = np.ascontiguousarray(system["base_rhs"], dtype=np.int64)
    require(matrix.shape == (282, EDGE_VARIABLES), "orbit0/A pointed matrix shape changed")
    require(base_rhs.shape == (282,), "orbit0/A pointed base RHS shape changed")
    require(
        affine.matrix_sha256(matrix) == EXPECTED_POINTED_MATRIX_SHA256,
        "orbit0/A pointed matrix certificate changed",
    )
    require(
        affine.matrix_sha256(base_rhs[None, :]) == EXPECTED_POINTED_BASE_RHS_SHA256,
        "orbit0/A pointed RHS certificate changed",
    )
    require(np.all(matrix[0] == 1), "pointed edge-count row changed")
    require(int(base_rhs[0]) == EDGE_COUNT, "pointed edge-count RHS changed")
    require(np.all(base_rhs[1:281] == 0), "pointed direction base RHS changed")
    require(int(base_rhs[281]) == 1, "pointed branch-A fixed-edge RHS changed")

    return {
        "provenance": provenance,
        "kernel": kernel,
        "kernel_audit": kernel_audit,
        "case": case,
        "system": system,
        "orbits": orbits,
        "orbit_source": orbit_source,
        "leaf_audit": leaf_audit,
        "target_pairs": target_pairs,
        "leaf_certificate_rows": leaf_certificate_rows,
        "compact_input_audit_sha256": json_sha256(compact_input_audit),
    }


def catalog_domain_audit(case: dict, target_pairs: list[tuple[int, dict]], kernel: np.ndarray):
    """Materialize and audit every exact U/M table used by the four cases."""
    specs = compact.direction_specs(case)
    cache: dict[tuple[int, int], np.ndarray] = {}
    records = []
    uses = Counter()
    for _leaf_index, leaf in target_pairs:
        for direction, (mean, catalog_class) in enumerate(
            zip(leaf["scaled_means"], leaf["catalog_classes"])
        ):
            mask = int(case["direction_masks"][direction])
            key = (mask, int(mean))
            uses[key] += 1
            if key in cache:
                continue
            catalog = np.ascontiguousarray(affine.mapped_catalog(*key), dtype=np.int16)
            expected_rows = {"U": 1, "M": 1_764}[str(catalog_class)]
            require(catalog.shape == (expected_rows, 35), "exact catalog shape changed")
            require(len(np.unique(catalog, axis=0)) == len(catalog), "exact catalog repeats a row")
            require(np.all((0 <= catalog) & (catalog <= 13)), "exact catalog left score bounds")
            parity = np.asarray(specs[direction]["parity"], dtype=np.int16)
            require(
                np.all(catalog % 2 == parity[None, :]),
                "exact catalog parity disagrees with pointed case",
            )
            require(
                not np.any(np.asarray(kernel, dtype=np.int64) @ catalog.T),
                "exact catalog left the primitive degree-two kernel",
            )
            require(
                np.all(2 * catalog.sum(axis=1, dtype=np.int64) == 5 * int(mean)),
                "exact catalog mean changed",
            )
            cache[key] = catalog
            records.append(
                {
                    "direction": direction,
                    "mask": mask,
                    "b": mask.bit_count(),
                    "scaled_mean": int(mean),
                    "catalog_class": str(catalog_class),
                    "rows": len(catalog),
                    "columns": 35,
                    "sha256_int16": array_sha256(catalog),
                    "rows_unique": True,
                    "bounds_parity_kernel_and_mean_audited": True,
                }
            )

    records.sort(key=lambda row: (row["direction"], row["scaled_mean"], row["mask"]))
    return cache, {
        "unique_mask_mean_domains": len(cache),
        "direction_domain_uses_across_four_leaves": sum(uses.values()),
        "expected_direction_domain_uses": 4 * 8,
        "records": records,
        "records_sha256": json_sha256(records),
        "complete_exact_catalog_tables_materialized": True,
        "every_domain_has_exactly_1_U_or_1764_M_rows": True,
    }


def exact_rhs(system: dict, slacks: list[list]) -> list:
    base_rhs = np.asarray(system["base_rhs"], dtype=np.int64)
    rhs: list = [int(value) for value in base_rhs]
    for direction in range(8):
        for point in range(35):
            rhs[1 + 35 * direction + point] = 13 - slacks[direction][point]
    return rhs


def add_full_edge_lift(
    model,
    variables: dict,
    data: dict,
    case: dict,
    system: dict,
) -> tuple[dict, dict]:
    """Couple the exact catalog model to the audited full binary edge model."""
    edges = tuple(tuple(int(value) for value in edge) for edge in data["edges"])
    require(len(edges) == len(set(edges)) == EDGE_VARIABLES, "full edge order changed")
    affine_edges = tuple(
        tuple(int(value) for value in edge) for edge in affine.geometry(P, "affine")["edges"]
    )
    require(edges == affine_edges, "full and pointed edge-variable orders disagree")
    edge_index = {edge: index for index, edge in enumerate(edges)}
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]

    matrix = np.asarray(system["matrix"], dtype=np.int64)
    rhs = exact_rhs(system, variables["slacks"])
    require(matrix.shape == (len(rhs), len(selected)), "pointed edge coupling shape changed")
    for row_index, (row, target) in enumerate(zip(matrix, rhs)):
        support = np.flatnonzero(row)
        require(len(support) > 0, f"empty pointed equation row {row_index}")
        model.add(
            sum(int(row[index]) * selected[int(index)] for index in support) == target
        )

    # Exact full-model constraints are retained explicitly even where a
    # pointed matrix row makes one redundant; this gives propagation and a
    # direct correspondence with the audited reference implementation.
    model.add(sum(selected) == EDGE_COUNT)
    for row in case["fixed_infinity_star_edges"]:
        edge = tuple(int(value) for value in row["edge"])
        require(edge in edge_index, f"fixed branch edge is missing: {edge}")
        model.add(selected[edge_index[edge]] == int(row["value"]))

    boundary = tuple(int(value) for value in case["fixed_boundary_graph_vertices"])
    boundary_set = set(boundary)
    require(len(boundary) == len(boundary_set) == 8, "fixed z=7 boundary changed")
    incident = [[] for _ in range(N)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(selected[index])
        incident[b].append(selected[index])
    require(all(len(row) == N - 1 for row in incident), "complete incidence changed")
    for vertex in range(N):
        if vertex in boundary_set:
            model.add_bool_xor(incident[vertex])
        else:
            model.add_bool_xor([~incident[vertex][0], *incident[vertex][1:]])

    signs = np.asarray(data["edge_signs"], dtype=np.int8)
    negative = [selected[index] for index in np.flatnonzero(signs < 0)]
    require(negative, "Paley graph has no negative edges")
    model.add_bool_xor([~negative[0], *negative[1:]])

    specs = full.direction_specs(data, case)
    compact_specs = variables["specs"]
    require(len(specs) == len(compact_specs) == 8, "direction count changed")
    for direction, (spec, compact_spec, mean) in enumerate(
        zip(specs, compact_specs, variables["means"])
    ):
        require(int(spec["direction_index"]) == direction, "full direction order changed")
        require(int(compact_spec["direction_index"]) == direction, "compact direction order changed")
        require(int(spec["eps"]) == int(compact_spec["eps"]), "direction type mismatch")
        require(int(spec["floor"]) == int(compact_spec["floor"]), "direction floor mismatch")
        coefficients = np.asarray(spec["coefficients"], dtype=np.int64)
        support = np.flatnonzero(coefficients)
        model.add(
            mean
            == sum(
                int(coefficients[index]) * selected[int(index)] for index in support
            )
            - 21
        )

    complete_shell_constraints = {}
    for eps in (-1, 1):
        normalized = eps * np.asarray(data["features"][eps], dtype=np.int8)
        for feature in normalized:
            bad = np.flatnonzero(feature < 0)
            model.add(sum(selected[int(index)] for index in bad) <= 13)
        complete_shell_constraints[str(eps)] = len(normalized)

    saturated_specs = full.saturated_affine_specs(data, specs)
    for row in saturated_specs:
        coefficients = np.asarray(row["coefficients"], dtype=np.int64)
        support = np.flatnonzero(coefficients)
        model.add(
            sum(int(coefficients[index]) * selected[int(index)] for index in support)
            == int(row["rhs"])
        )

    return {
        "selected": selected,
        "edges": edges,
        "full_specs": specs,
        "saturated_specs": saturated_specs,
    }, {
        "binary_edge_variables": len(selected),
        "exact_pointed_integer_edge_catalog_coupling_equations": len(matrix),
        "explicit_exact_edge_count_equations": 1,
        "explicit_pointed_fixed_edge_equations": len(case["fixed_infinity_star_edges"]),
        "native_boundary_xor_constraints": N,
        "native_positive_paley_product_xor_constraints": 1,
        "exact_full_direction_mean_links": len(specs),
        "complete_full_eigenshell_constraints": complete_shell_constraints,
        "complete_full_eigenshell_constraints_total": sum(complete_shell_constraints.values()),
        "saturated_affine_equalities": len(saturated_specs),
        "full_and_pointed_edge_orders_identical": True,
        "pointed_rhs_direction_formula": "bad_count = 13 - exact_catalog_slack",
        "only_explicit_problem_symmetry_break": "audited branch-A normalization edge (infinity,0)=1",
        "additional_unproved_leaf_symmetry_breaks": 0,
    }


def build_case_model(
    context: dict,
    data: dict,
    leaf: dict,
    catalog_cache: dict[tuple[int, int], np.ndarray],
) -> tuple[object, dict, dict]:
    model, variables, compact_construction = compact.build_model(
        context["case"], context["system"], context["kernel"]
    )
    proto_before = model.Proto()
    base_variables = len(proto_before.variables)
    base_constraints = len(proto_before.constraints)

    fixed_means = tuple(int(value) for value in leaf["scaled_means"])
    require(len(fixed_means) == 8, "fixed leaf mean vector changed length")
    for direction, value in enumerate(fixed_means):
        model.add(variables["means"][direction] == value)

    catalog_rows = []
    for direction, (mean, catalog_class) in enumerate(
        zip(fixed_means, leaf["catalog_classes"])
    ):
        mask = int(context["case"]["direction_masks"][direction])
        catalog = catalog_cache[(mask, mean)]
        expected_rows = {"U": 1, "M": 1_764}[str(catalog_class)]
        require(len(catalog) == expected_rows, "leaf catalog table size changed")
        model.add_allowed_assignments(
            variables["slacks"][direction],
            [[int(value) for value in row] for row in catalog],
        )
        catalog_rows.append(expected_rows)

    edge_variables, edge_construction = add_full_edge_lift(
        model,
        variables,
        data,
        context["case"],
        context["system"],
    )
    validation_error = model.validate()
    require(not validation_error, f"exact H0_S0_M7 model is invalid: {validation_error}")
    proto_after = model.Proto()
    construction = {
        "compact_constructor": "p7_infinity7_positive_z7_pointed_compact_cpsat.build_model",
        "compact_base_variables": base_variables,
        "compact_base_constraints": base_constraints,
        "compact_construction": compact_construction,
        "fixed_mean_equalities": 8,
        "fixed_scaled_means": list(fixed_means),
        "complete_exact_catalog_table_constraints": 8,
        "complete_exact_catalog_table_row_counts": catalog_rows,
        "complete_exact_catalog_rows_total_with_direction_multiplicity": sum(catalog_rows),
        "catalog_table_semantics": "one complete mapped exact catalog row per direction",
        "edge_model": edge_construction,
        "total_variables": len(proto_after.variables),
        "total_constraints": len(proto_after.constraints),
        "model_validation": "passed",
        "model_stats_sha256": hashlib.sha256(model.model_stats().encode("utf-8")).hexdigest(),
    }
    variables = variables | edge_variables
    return model, variables, construction


def exact_edge_catalog_audit(
    solver,
    context: dict,
    leaf: dict,
    variables: dict,
    catalog_cache: dict[tuple[int, int], np.ndarray],
) -> dict:
    selected = np.asarray(
        [int(solver.value(variable)) for variable in variables["selected"]],
        dtype=np.int64,
    )
    slacks = np.asarray(
        [
            [int(solver.value(variable)) for variable in direction]
            for direction in variables["slacks"]
        ],
        dtype=np.int64,
    )
    matrix = np.asarray(context["system"]["matrix"], dtype=np.int64)
    numeric_rhs = np.asarray(context["system"]["base_rhs"], dtype=np.int64).copy()
    numeric_rhs[1:281] = (13 - slacks).reshape(-1)
    lhs = matrix @ selected
    equation_ok = np.array_equal(lhs, numeric_rhs)

    memberships = []
    membership_ok = True
    for direction, mean in enumerate(leaf["scaled_means"]):
        mask = int(context["case"]["direction_masks"][direction])
        catalog = catalog_cache[(mask, int(mean))].astype(np.int64)
        matches = np.flatnonzero(np.all(catalog == slacks[direction][None, :], axis=1))
        valid = len(matches) == 1
        membership_ok &= valid
        memberships.append(
            {
                "direction": direction,
                "mask": mask,
                "scaled_mean": int(mean),
                "catalog_class": str(leaf["catalog_classes"][direction]),
                "catalog_rows": len(catalog),
                "matching_exact_catalog_row_indices": [int(value) for value in matches],
                "unique_exact_catalog_membership": valid,
                "selected_slack_sha256_int64": array_sha256(slacks[direction]),
            }
        )

    checks = {
        "selected_vector_binary": bool(np.all((selected == 0) | (selected == 1))),
        "exactly_29_selected_edges": int(selected.sum()) == EDGE_COUNT,
        "all_282_pointed_integer_equations": equation_ok,
        "all_eight_slacks_have_unique_exact_catalog_membership": membership_ok,
        "fixed_leaf_means_match": [
            int(solver.value(variable)) for variable in variables["means"]
        ]
        == [int(value) for value in leaf["scaled_means"]],
    }
    return {
        "valid": all(checks.values()),
        "checks": checks,
        "selected_edge_vector_sha256_uint8": array_sha256(selected.astype(np.uint8)),
        "slack_matrix_sha256_int64": array_sha256(slacks),
        "pointed_lhs_sha256_int64": array_sha256(lhs),
        "pointed_rhs_sha256_int64": array_sha256(numeric_rhs),
        "catalog_memberships": memberships,
    }


def solve_case(
    context: dict,
    data: dict,
    target_ordinal: int,
    leaf_index: int,
    leaf: dict,
    catalog_cache: dict[tuple[int, int], np.ndarray],
    seconds: float,
    workers: int,
    seed: int,
    log_search_progress: bool,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    model, variables, construction = build_case_model(context, data, leaf, catalog_cache)
    build_seconds = time.time() - started
    case_key = f"orbit0_leaf{leaf_index}_branchA"

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(seconds)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 2
    solver.parameters.log_search_progress = bool(log_search_progress)
    status = solver.solve(model)
    require(status != cp_model.MODEL_INVALID, "validated exact model became MODEL_INVALID")
    status_name = solver.status_name(status)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    infeasible = status == cp_model.INFEASIBLE
    if feasible:
        mathematical_status = "SAT"
        rigorous_status = "rigorously_feasible_exact_full_edge_witness_audited"
    elif infeasible:
        mathematical_status = "UNSAT"
        rigorous_status = "rigorously_infeasible_exact_full_edge_catalog_model"
    else:
        mathematical_status = "UNKNOWN"
        rigorous_status = "unknown_or_timeout_no_mathematical_conclusion"

    result = {
        "target_ordinal": target_ordinal,
        "case_key": case_key,
        "branch_orbit_index": 0,
        "orbit_leaf_index": int(leaf_index),
        "pointed_star_branch": "A",
        "catalog_pattern": "H0_S0_M7",
        "scaled_means": [int(value) for value in leaf["scaled_means"]],
        "catalog_classes": [str(value) for value in leaf["catalog_classes"]],
        "solver_status": status_name,
        "mathematical_status": mathematical_status,
        "rigorous_status": rigorous_status,
        "feasible": feasible,
        "finite_infeasibility_certificate": infeasible,
        "case_decided": feasible or infeasible,
        "seconds_per_case": seconds,
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "solver_wall_time_seconds": solver.wall_time,
        "model_build_seconds": build_seconds,
        "elapsed_seconds": time.time() - started,
        "construction": construction,
    }
    if feasible:
        chosen_edges = [
            list(edge)
            for edge, variable in zip(variables["edges"], variables["selected"])
            if solver.value(variable)
        ]
        solver_means = [int(solver.value(variable)) for variable in variables["means"]]
        solver_residues = {
            eps: int(solver.value(variable))
            for eps, variable in variables["common_residues"].items()
        }
        full_audit = full.direct_witness_audit(
            data,
            context["case"],
            chosen_edges,
            variables["full_specs"],
            variables["saturated_specs"],
            solver_means,
            solver_residues,
        )
        compact_witness, compact_audit = compact.direct_witness_audit(
            solver,
            context["case"],
            context["system"],
            context["kernel"],
            variables,
        )
        coupling_audit = exact_edge_catalog_audit(
            solver, context, leaf, variables, catalog_cache
        )
        require(full_audit["valid"], "SAT edge witness failed the full direct audit")
        require(compact_audit["valid"], "SAT catalog witness failed compact direct audit")
        require(coupling_audit["valid"], "SAT witness failed exact edge/catalog coupling audit")
        result["chosen_edges_H"] = chosen_edges
        result["solver_common_residues"] = {
            str(eps): value for eps, value in solver_residues.items()
        }
        result["exact_witness_audit"] = {
            "valid": True,
            "full_edge_audit": full_audit,
            "compact_catalog_audit": compact_audit,
            "exact_edge_catalog_coupling_audit": coupling_audit,
            "compact_witness": compact_witness,
        }
    return result


def run_identity(
    provenance_certificate: str,
    target_keys: list[str],
    seconds: float,
    workers: int,
    seed: int,
    smoke: bool,
) -> dict:
    row = {
        "provenance_certificate_sha256": provenance_certificate,
        "target_case_keys": target_keys,
        "seconds_per_case": seconds,
        "workers_per_case": workers,
        "seed_base": seed,
        "smoke": smoke,
    }
    return row | {"sha256": json_sha256(row)}


def summary_status(results: list[dict], full_target: bool, complete: bool, smoke: bool):
    sat = [row["case_key"] for row in results if row["mathematical_status"] == "SAT"]
    unsat = [row["case_key"] for row in results if row["mathematical_status"] == "UNSAT"]
    unknown = [row["case_key"] for row in results if row["mathematical_status"] == "UNKNOWN"]
    if not complete:
        status = "checkpoint_partial_exact_h0m7_run"
        conclusion = "selected target processing is incomplete"
    elif sat:
        status = "rigorous_exact_h0m7_full_edge_witness_found"
        conclusion = "at least one selected H0_S0_M7 leaf has an audited exact full-edge witness"
    elif unknown:
        status = "exact_h0m7_run_with_unknown_cases"
        conclusion = "at least one selected case is UNKNOWN; no exclusion follows for it"
    elif full_target and len(unsat) == 4:
        status = "complete_rigorous_orbit0_branchA_h0m7_representative_exclusion"
        conclusion = "all four orbit0/branchA H0_S0_M7 representative leaves are excluded"
    else:
        status = "complete_rigorous_selected_h0m7_case_exclusion"
        conclusion = "every selected H0_S0_M7 case is excluded; unselected target cases remain"
    if smoke and not sat:
        conclusion += "; bounded smoke mode makes no claim about unprocessed target cases"
    return status, conclusion, sat, unsat, unknown


def make_payload(
    context: dict,
    geometry_audit: dict,
    catalog_audit: dict,
    identity: dict,
    selected_targets: list[tuple[int, int, dict]],
    results: list[dict],
    smoke: bool,
    started: float,
) -> dict:
    target_keys = [f"orbit0_leaf{leaf_index}_branchA" for _ordinal, leaf_index, _leaf in selected_targets]
    require(
        [row["case_key"] for row in results] == target_keys[: len(results)],
        "case results are not the canonical selected prefix",
    )
    complete = len(results) == len(selected_targets)
    full_target = target_keys == [
        f"orbit0_leaf{leaf_index}_branchA" for leaf_index in EXPECTED_TARGET_LEAF_INDICES
    ]
    status, conclusion, sat, unsat, unknown = summary_status(
        results, full_target, complete, smoke
    )
    return {
        "experiment": EXPERIMENT,
        "schema_version": SCHEMA_VERSION,
        "status": status,
        "rigorous_conclusion": conclusion,
        "p": P,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "exactly orbit0/branchA H0_S0_M7 representative leaves",
        "smoke": smoke,
        "host_lucky_used": False,
        "mesh_hosts_used": [],
        "model_semantics": {
            "binary_edges": "all 1,225 complete-graph edges",
            "selected_edge_count": 29,
            "boundary": "native XOR at all 50 vertices",
            "paley_product": "positive via native parity constraint",
            "catalogs": "all eight fixed to complete exact U/M row tables at the reconstructed means",
            "edge_catalog_coupling": "all 282 orbit0/A pointed equations imposed over the integers",
            "eigenshells": "all unique edge-feature rows from both complete cached shells",
            "SAT": "exact witness, accepted only after three direct audits",
            "UNSAT": "rigorous finite exclusion of that exact selected case",
            "UNKNOWN": "no mathematical conclusion",
        },
        "symmetry": {
            "explicit_problem_level_break": "branch A meets the line, normalized by its audited stabilizer to (infinity,0)=1",
            "normalization_already_proved_by_reused_constructor": True,
            "additional_leaf_symmetry_breaks": 0,
            "cp_sat_internal_symmetry_level": 3,
        },
        "checkpoint_semantics": "atomic checkpoint after each completed case; no within-case CP-SAT resume",
        "run_identity": identity,
        "provenance": context["provenance"],
        "provenance_certificate_sha256": identity["provenance_certificate_sha256"],
        "input_certificates": {
            "target_leaf_indices": list(EXPECTED_TARGET_LEAF_INDICES),
            "target_leaf_rows": context["leaf_certificate_rows"],
            "target_leaf_sha256": EXPECTED_TARGET_LEAF_SHA256,
            "kernel_audit": context["kernel_audit"],
            "pointed_matrix_sha256_int16": EXPECTED_POINTED_MATRIX_SHA256,
            "pointed_base_rhs_sha256_int64": EXPECTED_POINTED_BASE_RHS_SHA256,
            "pointed_case_sha256": EXPECTED_POINTED_CASE_SHA256,
            "orbit_source": context["orbit_source"],
            "mean_leaf_coverage": context["leaf_audit"],
            "compact_input_audit_sha256": context["compact_input_audit_sha256"],
            "full_geometry_audit": geometry_audit,
            "exact_catalog_domain_audit": catalog_audit,
        },
        "selected_target_case_keys": target_keys,
        "selected_target_count": len(selected_targets),
        "full_four_target_scope": full_target,
        "completed_case_count": len(results),
        "run_complete": complete,
        "all_selected_cases_decided": complete and not unknown,
        "all_four_target_leaves_excluded": full_target and complete and len(unsat) == 4,
        "SAT_case_keys": sat,
        "UNSAT_case_keys": unsat,
        "UNKNOWN_case_keys": unknown,
        "case_results_sha256": ordered_rows_sha256(results),
        "case_results": results,
        "elapsed_seconds": time.time() - started,
    }


def validate_checkpoint(payload: dict, identity: dict, target_keys: list[str]) -> list[dict]:
    require(payload.get("experiment") == EXPERIMENT, "checkpoint experiment changed")
    require(payload.get("schema_version") == SCHEMA_VERSION, "checkpoint schema changed")
    require(payload.get("run_identity") == identity, "checkpoint run identity changed")
    require(payload.get("selected_target_case_keys") == target_keys, "checkpoint targets changed")
    results = list(payload.get("case_results", []))
    require(
        payload.get("case_results_sha256") == ordered_rows_sha256(results),
        "checkpoint case-result digest failed",
    )
    require(
        [row.get("case_key") for row in results] == target_keys[: len(results)],
        "checkpoint results are not a canonical target prefix",
    )
    for row in results:
        require(
            row.get("mathematical_status") in {"SAT", "UNSAT", "UNKNOWN"},
            "checkpoint contains an invalid mathematical status",
        )
        if row["mathematical_status"] == "SAT":
            require(
                row.get("exact_witness_audit", {}).get("valid") is True,
                "checkpoint SAT row lacks a valid exact witness audit",
            )
        if row["mathematical_status"] == "UNSAT":
            require(row.get("finite_infeasibility_certificate") is True, "bad UNSAT row")
        if row["mathematical_status"] == "UNKNOWN":
            require(row.get("case_decided") is False, "UNKNOWN checkpoint row claims a decision")
    return results


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    context = audited_inputs()
    data, geometry_audit = full.load_full_geometry()
    catalog_cache, catalog_audit = catalog_domain_audit(
        context["case"], context["target_pairs"], context["kernel"]
    )
    provenance_certificate = json_sha256(
        {
            "provenance": context["provenance"],
            "target_leaf_sha256": EXPECTED_TARGET_LEAF_SHA256,
            "kernel_sha256": EXPECTED_KERNEL_SHA256,
            "pointed_matrix_sha256": EXPECTED_POINTED_MATRIX_SHA256,
            "pointed_base_rhs_sha256": EXPECTED_POINTED_BASE_RHS_SHA256,
            "catalog_records_sha256": catalog_audit["records_sha256"],
            "full_geometry_audit": geometry_audit,
        }
    )

    by_ordinal = [
        (ordinal, leaf_index, leaf)
        for ordinal, (leaf_index, leaf) in enumerate(context["target_pairs"])
    ]
    requested = sorted(set(args.case_index)) if args.case_index else list(range(4))
    if args.smoke and not args.case_index:
        requested = [0]
    if args.smoke:
        require(len(requested) == 1, "smoke mode accepts exactly one target case")
    selected_targets = [by_ordinal[index] for index in requested]
    seconds = min(args.seconds, 2.0) if args.smoke else args.seconds
    workers = min(args.workers, 2) if args.smoke else args.workers
    target_keys = [
        f"orbit0_leaf{leaf_index}_branchA" for _ordinal, leaf_index, _leaf in selected_targets
    ]
    identity = run_identity(
        provenance_certificate,
        target_keys,
        seconds,
        workers,
        args.seed,
        args.smoke,
    )

    results: list[dict] = []
    if args.resume:
        require(args.checkpoint is not None, "--resume requires --checkpoint")
        require(args.checkpoint.is_file(), "resume checkpoint does not exist")
        results = validate_checkpoint(
            json.loads(args.checkpoint.read_text(encoding="utf-8")),
            identity,
            target_keys,
        )

    for position in range(len(results), len(selected_targets)):
        ordinal, leaf_index, leaf = selected_targets[position]
        print(
            f"[{position + 1}/{len(selected_targets)}] {target_keys[position]} "
            f"seconds={seconds:g} workers={workers}",
            flush=True,
        )
        result = solve_case(
            context,
            data,
            ordinal,
            leaf_index,
            leaf,
            catalog_cache,
            seconds,
            workers,
            args.seed + ordinal,
            args.log_search_progress,
        )
        results.append(result)
        print(
            f"{result['case_key']}: {result['mathematical_status']} "
            f"solver={result['solver_status']} wall={result['solver_wall_time_seconds']:.3f}s",
            flush=True,
        )
        if args.checkpoint is not None:
            checkpoint = make_payload(
                context,
                geometry_audit,
                catalog_audit,
                identity,
                selected_targets,
                results,
                args.smoke,
                started,
            )
            full.pointed.atomic_write(args.checkpoint, checkpoint)

    payload = make_payload(
        context,
        geometry_audit,
        catalog_audit,
        identity,
        selected_targets,
        results,
        args.smoke,
        started,
    )
    full.pointed.atomic_write(args.output, payload)
    if args.checkpoint is not None and args.checkpoint != args.output:
        full.pointed.atomic_write(args.checkpoint, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seconds",
        "--timeout",
        dest="seconds",
        type=positive_float,
        default=3_600.0,
        help="CP-SAT wall-time limit for each selected case",
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=16,
        help="CP-SAT search workers for each selected case",
    )
    parser.add_argument("--seed", type=int, default=15718001)
    parser.add_argument(
        "--case-index",
        type=int,
        choices=range(4),
        action="append",
        default=[],
        help="select target ordinal 0..3; repeat to shard, default is all four",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        help="atomically save completed-case prefixes here after every case",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="resume a hash- and configuration-matched between-case checkpoint",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="run one target only, capped at two seconds and two workers",
    )
    parser.add_argument("--log-search-progress", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = run(args)
    print(
        json.dumps(
            {
                "experiment": result["experiment"],
                "status": result["status"],
                "rigorous_conclusion": result["rigorous_conclusion"],
                "selected_target_case_keys": result["selected_target_case_keys"],
                "completed_case_count": result["completed_case_count"],
                "SAT_case_keys": result["SAT_case_keys"],
                "UNSAT_case_keys": result["UNSAT_case_keys"],
                "UNKNOWN_case_keys": result["UNKNOWN_case_keys"],
                "case_results_sha256": result["case_results_sha256"],
                "output": str(args.output),
                "checkpoint": str(args.checkpoint) if args.checkpoint else None,
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
