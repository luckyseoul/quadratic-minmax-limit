#!/usr/bin/env python3
"""Exact pointed full-shell CP-SAT lift for positive p=7, z=7.

The two finite boundary representatives are the audited affine-line orbits
from the exact z>=2 orbit census.  For each line, its stabilizer gives the
following exhaustive normalization of the odd, hence nonempty, infinity
star:

* A: the star meets the line; move a present line edge to ``(infinity, 0)``;
* B: the star avoids the line; set all seven line edges to zero and move a
  present outside edge to the audited outside representative.

Unlike the modular projection sieve, this script builds one exact edge model
per selected orbit/branch.  Direction means remain integer variables, so no
mean-allocation leaves are enumerated.  ``INFEASIBLE`` is a rigorous finite
exclusion of that pointed case; any feasible result is audited directly
against the complete cached eigenshells before it is reported.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
import p7_infinity7_positive_z7_mod7_projection as z7_parent  # noqa: E402
import p7_infinity7_positive_z7_pointed_mod7 as pointed  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


P = 7
Q = P * P
N = Q + 1
EDGE_COUNT = 4 * P + 1
EXPECTED_EDGE_VARIABLES = N * (N - 1) // 2
EXPECTED_ORBITS = 2
EXPECTED_BRANCHES = ("A", "B")
EXPECTED_CASES = EXPECTED_ORBITS * len(EXPECTED_BRANCHES)
SHELL_PATHS = {
    1: Path("/tmp/maxplus_p7.npy"),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def audited_pointed_cases() -> tuple[list[dict], dict, dict]:
    """Recompute the pointed line normalizations without mean leaves.

    This deliberately reuses the group, direction-action, and partition
    routines from the audited pointed mod-seven implementation.  The full
    CP-SAT model does not need, and does not generate, its 2,160 mean leaves.
    """
    require(P == pointed.P and Q == pointed.Q, "pointed p=7 constants changed")
    require(EDGE_COUNT == pointed.EDGE_COUNT, "pointed edge count changed")
    require(
        pointed.EXPECTED_LINE_STABILIZER == 84,
        "pointed line-stabilizer order changed",
    )
    orbits, orbit_source = z7_parent.load_z7_orbits()
    require(len(orbits) == EXPECTED_ORBITS, "z=7 orbit count changed")

    directions, types, labels = pointed.direction_metadata()
    require(tuple(types) == tuple(pointed.DIRECTION_TYPES), "direction types changed")
    require(tuple(labels) == tuple(pointed.LABELS), "direction labels changed")
    group, recomputed_group_audit = pointed.square_affine_semilinear_group(
        directions, types, labels
    )

    source = json.loads(z7_parent.ORBIT_EVIDENCE.read_text(encoding="utf-8"))
    source_group = source["group_audit"]
    group_fields = (
        "group_size",
        "permutation_sha256",
        "induced_direction_action_sha256",
        "all_preserve_paley_difference_signs",
        "all_permute_projective_directions",
        "all_preserve_direction_types",
    )
    for key in group_fields:
        require(
            recomputed_group_audit[key] == source_group[key],
            f"group audit field changed: {key}",
        )

    fibres = tuple(
        tuple(
            tuple(u for u in range(Q) if labels[direction][u] == fibre)
            for fibre in range(P)
        )
        for direction in range(P + 1)
    )
    partition_lookup = {
        pointed._partition_signature(fibres[direction]): direction  # noqa: SLF001
        for direction in range(P + 1)
    }
    require(len(partition_lookup) == P + 1, "affine direction partitions changed")

    cases: list[dict] = []
    orbit_audits = []
    for orbit_index, orbit in enumerate(orbits):
        require(
            int(orbit["branch_orbit_index"]) == orbit_index,
            "z=7 branch-orbit ordering changed",
        )
        line = tuple(int(value) for value in orbit["representative"])
        line_set = set(line)
        outside = tuple(value for value in range(Q) if value not in line_set)
        require(len(line) == len(line_set) == P, "z=7 representative is not a line")
        require(0 in line_set, "z=7 representative no longer contains point zero")
        require(len(outside) == Q - P, "line complement has the wrong size")

        stabilizer = tuple(
            permutation
            for permutation in group
            if {permutation[value] for value in line_set} == line_set
        )
        require(
            len(stabilizer) == pointed.EXPECTED_LINE_STABILIZER,
            "line stabilizer order changed",
        )
        require(
            len(stabilizer) == Q * (Q - 1) // int(orbit["size"]),
            "line orbit-stabilizer identity failed",
        )
        outside_representative = min(outside)
        line_point_orbit = {permutation[0] for permutation in stabilizer}
        outside_point_orbit = {
            permutation[outside_representative] for permutation in stabilizer
        }
        require(line_point_orbit == line_set, "stabilizer is not transitive on the line")
        require(
            outside_point_orbit == set(outside),
            "stabilizer is not transitive off the line",
        )

        actions = {
            pointed.induced_direction_action(permutation, labels, partition_lookup)
            for permutation in stabilizer
        }
        require(actions, "line stabilizer induced no direction actions")
        for action in actions:
            require(
                all(types[source_index] == types[target_index] for source_index, target_index in enumerate(action)),
                "line stabilizer changed a direction type",
            )
            require(
                all(
                    orbit["b_values"][source_index]
                    == orbit["b_values"][target_index]
                    for source_index, target_index in enumerate(action)
                ),
                "line stabilizer changed the boundary fibre-parity profile",
            )

        action_digest = hashlib.sha256()
        for action in sorted(actions):
            action_digest.update(bytes(action))
            action_digest.update(b"\xff")

        boundary = (0, *(value + 1 for value in line))
        require(len(boundary) == len(set(boundary)) == P + 1, "boundary lift changed")
        require(tuple(sorted(boundary)) == boundary, "boundary lift is not canonical")
        branch_specs = {
            "A": ((0, 1, 1),),
            "B": tuple((0, value + 1, 0) for value in line)
            + ((0, outside_representative + 1, 1),),
        }
        require(branch_specs["A"] == ((0, 1, 1),), "branch A representative changed")
        require(
            tuple(row[1] - 1 for row in branch_specs["B"][:-1]) == line,
            "branch B does not fix every line edge absent",
        )
        require(
            branch_specs["B"][-1] == (0, outside_representative + 1, 1),
            "branch B outside representative changed",
        )
        require(
            len({(a, b) for a, b, _value in branch_specs["B"]}) == P + 1,
            "branch B repeats a fixed edge",
        )

        for branch in EXPECTED_BRANCHES:
            cases.append(
                {
                    "branch_orbit_index": orbit_index,
                    "source_orbit_index": int(orbit["source_orbit_index"]),
                    "orbit_size": int(orbit["size"]),
                    "representative_finite_field": list(line),
                    "direction_masks": list(orbit["masks"]),
                    "b_values": list(orbit["b_values"]),
                    "fixed_boundary_graph_vertices": list(boundary),
                    "pointed_star_branch": branch,
                    "fixed_infinity_star_edges": [
                        {
                            "edge": [a, b],
                            "finite_field_point": b - 1,
                            "value": value,
                        }
                        for a, b, value in branch_specs[branch]
                    ],
                }
            )

        orbit_audits.append(
            {
                "branch_orbit_index": orbit_index,
                "source_orbit_index": int(orbit["source_orbit_index"]),
                "orbit_size": int(orbit["size"]),
                "representative_finite_field": list(line),
                "line_stabilizer_size": len(stabilizer),
                "line_point_orbit_size": len(line_point_orbit),
                "outside_point_orbit_size": len(outside_point_orbit),
                "branch_A_representative_finite_point": 0,
                "branch_B_representative_finite_point": outside_representative,
                "induced_direction_action_count": len(actions),
                "induced_direction_action_sha256": action_digest.hexdigest(),
            }
        )

    require(len(cases) == EXPECTED_CASES, "pointed case census changed")
    require(
        len(
            {
                (row["branch_orbit_index"], row["pointed_star_branch"])
                for row in cases
            }
        )
        == EXPECTED_CASES,
        "pointed case key repeated",
    )

    all_nonempty_stars = 2**Q - 1
    branch_a_nonempty = (2**P - 1) * 2 ** (Q - P)
    branch_b_nonempty = 2 ** (Q - P) - 1
    require(
        branch_a_nonempty + branch_b_nonempty == all_nonempty_stars,
        "A/B nonempty-star partition failed",
    )
    all_odd_stars = 2 ** (Q - 1)
    branch_b_odd = 2 ** (Q - P - 1)
    branch_a_odd = all_odd_stars - branch_b_odd
    require(branch_a_odd + branch_b_odd == all_odd_stars, "A/B odd-star partition failed")

    audit = {
        "pointed_implementation": str(
            (ROOT / "scripts" / "p7_infinity7_positive_z7_pointed_mod7.py").relative_to(ROOT)
        ),
        "pointed_implementation_sha256": file_sha256(
            ROOT / "scripts" / "p7_infinity7_positive_z7_pointed_mod7.py"
        ),
        "reused_pointed_helpers": [
            "direction_metadata",
            "square_affine_semilinear_group",
            "_partition_signature",
            "induced_direction_action",
        ],
        "mean_allocation_enumeration_used": False,
        "group_size": len(group),
        "group_matches_exact_orbit_evidence": True,
        "all_group_elements_preserve_paley_difference_signs": True,
        "branches": {
            "A": "star meets line; normalize one present line edge to infinity--finite-point-0",
            "B": "star avoids line; fix all seven line edges absent and normalize one present outside edge",
        },
        "nonempty_star_partition": {
            "all": all_nonempty_stars,
            "A_meets_line": branch_a_nonempty,
            "B_avoids_line_but_is_nonempty": branch_b_nonempty,
        },
        "odd_boundary_star_partition": {
            "all": all_odd_stars,
            "A_meets_line": branch_a_odd,
            "B_avoids_line": branch_b_odd,
        },
        "A_and_B_disjoint": True,
        "A_and_B_cover_every_nonempty_star": True,
        "boundary_membership_makes_infinity_star_odd_and_nonempty": True,
        "normalizations_valid_by_stabilizer_transitivity": True,
        "per_orbit": orbit_audits,
    }
    return cases, orbit_source, audit


def load_full_geometry() -> tuple[dict, dict]:
    missing = [str(path) for path in SHELL_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "missing complete p=7 eigenshell cache(s): " + ", ".join(missing)
        )
    C = np.rint(paley_conference_prime_power(P)).astype(np.int16)
    require(C.shape == (N, N), "Paley conference matrix shape changed")
    require(np.array_equal(C, C.T), "Paley conference matrix is not symmetric")
    require(np.count_nonzero(np.diag(C)) == 0, "Paley conference diagonal changed")
    off_diagonal = C[~np.eye(N, dtype=bool)]
    require(set(np.unique(off_diagonal).tolist()) == {-1, 1}, "Paley signs are not +/-1")
    require(
        np.array_equal(C.astype(np.int32) @ C.astype(np.int32), Q * np.eye(N, dtype=np.int32)),
        "Paley conference identity failed",
    )
    edges = tuple((a, b) for a in range(N) for b in range(a + 1, N))
    require(len(edges) == EXPECTED_EDGE_VARIABLES, "complete graph edge count changed")
    require(len(set(edges)) == EXPECTED_EDGE_VARIABLES, "complete graph repeats an edge")

    plus = np.rint(np.load(SHELL_PATHS[1])).astype(np.int8)
    require(
        plus.ndim == 2 and plus.shape[0] > 0 and plus.shape[1] == N,
        "full plus shell shape changed",
    )
    require(set(np.unique(plus).tolist()) <= {-1, 1}, "full plus shell is not Boolean")
    require(
        np.array_equal(plus.astype(np.int16) @ C, P * plus.astype(np.int16)),
        "canonical plus cache failed its exact eigenshell audit",
    )

    # Multiplication by a nonsquare negates every finite Paley difference.
    # Switching the infinity coordinate then gives a signed permutation S
    # with S C S^T = -C, hence an exact bijection Max+ -> Max-.
    q, multiply, _add, character, *_rest = field_ctx(P)
    nonsquare = next(value for value in range(1, q) if character(value) == -1)
    permutation = (0, *(1 + multiply(nonsquare, value) for value in range(q)))
    require(sorted(permutation) == list(range(N)), "nonsquare action is not a permutation")
    switching = np.ones(N, dtype=np.int8)
    switching[0] = -1
    transformed_C = (
        C[np.ix_(permutation, permutation)]
        * switching[:, None].astype(np.int16)
        * switching[None, :].astype(np.int16)
    )
    require(np.array_equal(transformed_C, -C), "signed nonsquare action is not anti-Paley")
    minus = np.empty_like(plus)
    minus[:, permutation] = plus * switching[None, :]
    require(
        np.array_equal(minus.astype(np.int16) @ C, -P * minus.astype(np.int16)),
        "derived minus shell failed its exact eigenshell audit",
    )
    require(
        len(np.unique(minus, axis=0)) == len(minus),
        "derived minus shell contains duplicate vectors",
    )

    left = np.fromiter((a for a, _b in edges), dtype=np.int32)
    right = np.fromiter((b for _a, b in edges), dtype=np.int32)
    signs = C[left, right].astype(np.int8)
    shells = {1: plus, -1: minus}
    features = {
        eps: np.unique(
            np.ascontiguousarray(
                rows[:, left] * rows[:, right] * signs[None, :],
                dtype=np.int8,
            ),
            axis=0,
        )
        for eps, rows in shells.items()
    }
    data = {
        "C": C,
        "n": N,
        "edges": edges,
        "edge_signs": signs,
        "shells": shells,
        "features": features,
        "shell_mode": "full_signed_transfer",
    }

    shell_rows = {}
    feature_rows = {}
    for eps in (-1, 1):
        shell = np.asarray(data["shells"][eps])
        features = np.asarray(data["features"][eps])
        require(
            shell.ndim == 2 and shell.shape[0] > 0 and shell.shape[1] == N,
            "full shell shape changed",
        )
        require(set(np.unique(shell).tolist()) <= {-1, 1}, "full shell is not Boolean")
        require(
            features.ndim == 2
            and features.shape[0] > 0
            and features.shape[1] == EXPECTED_EDGE_VARIABLES,
            "feature shape changed",
        )
        require(set(np.unique(features).tolist()) <= {-1, 1}, "edge features are not Boolean")
        shell_rows[str(eps)] = int(shell.shape[0])
        feature_rows[str(eps)] = int(features.shape[0])
    return data, {
        "p": P,
        "vertices": N,
        "edge_variables": EXPECTED_EDGE_VARIABLES,
        "paley_conference_identity": "C^2=49I",
        "full_shell_rows": shell_rows,
        "unique_edge_feature_rows": feature_rows,
        "shell_file_sha256": {
            str(eps): file_sha256(path) for eps, path in SHELL_PATHS.items()
        },
        "minus_shell_derivation": {
            "source": "canonical complete plus shell",
            "nonsquare_field_element": nonsquare,
            "finite_action": "multiplication by the recorded nonsquare",
            "switching": "negate the infinity coordinate",
            "signed_permutation_anti_paley_identity": "S C S^T = -C",
            "anti_paley_identity_audited": True,
            "derived_minus_eigenshell_audited": True,
            "derived_minus_rows_unique": True,
        },
        "geometry_full_eigenshell_audit_passed": True,
    }


def direction_specs(data: dict, case: dict) -> list[dict]:
    C = data["C"]
    edges = data["edges"]
    boundary = tuple(int(value) for value in case["fixed_boundary_graph_vertices"])
    infinity_value = int(0 in boundary)
    specs = []
    for direction_index, direction in enumerate(projective_directions(P)):
        eps, labels = field_direction_data(P, direction)
        counts = [0] * P
        for vertex in boundary:
            if vertex:
                counts[labels[vertex - 1]] += 1
        odd_fibres = {index for index, value in enumerate(counts) if value & 1}
        parity_sign = -eps  # c_H = +1
        if infinity_value:
            parity_sign *= eps
        if len(odd_fibres) & 1:
            parity_sign *= -1
        phase = int(parity_sign == -1)
        floor = int(scaled_direction_floor(P, len(odd_fibres), phase))
        coefficients = np.fromiter(
            (
                1
                if a == 0
                else P
                if labels[a - 1] == labels[b - 1]
                else -eps * int(C[a, b])
                for a, b in edges
            ),
            dtype=np.int16,
            count=len(edges),
        )
        specs.append(
            {
                "direction_index": direction_index,
                "direction": tuple(int(value) for value in direction),
                "eps": int(eps),
                "labels": tuple(int(value) for value in labels),
                "B": frozenset(odd_fibres),
                "b": len(odd_fibres),
                "phase": phase,
                "floor": floor,
                "coefficients": coefficients,
            }
        )

    require(len(specs) == P + 1, "projective direction count changed")
    require(
        tuple(row["eps"] for row in specs) == tuple(pointed.DIRECTION_TYPES),
        "full-shell and pointed direction types disagree",
    )
    require(
        tuple(row["b"] for row in specs) == tuple(case["b_values"]),
        "full-graph boundary lift changed the pointed b-values",
    )
    require(
        tuple(sum(1 << value for value in row["B"]) for row in specs)
        == tuple(case["direction_masks"]),
        "full-graph boundary lift changed the pointed direction masks",
    )
    require(all(row["phase"] == 0 for row in specs), "positive z=7 phase changed")
    require(sorted(row["b"] for row in specs) == [1] + [7] * 7, "z=7 b-profile changed")
    return specs


def saturated_affine_specs(data: dict, specs: list[dict]) -> list[dict]:
    """Return every valid floor-saturation equality from the reference model."""
    C = np.asarray(data["C"], dtype=np.int16)
    edges = data["edges"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    rows = []
    for eps in (-1, 1):
        records = [row for row in specs if row["eps"] == eps]
        if sum(int(row["floor"]) for row in records) != 32:
            continue
        if any(len(row["B"]) not in (0, 2) for row in records):
            continue
        for record in records:
            labels = record["labels"]
            B = record["B"]
            phase = int(record["phase"])
            for chosen_fibres in itertools.combinations(range(P), 4):
                chosen_set = set(chosen_fibres)
                y = np.empty(N, dtype=np.int8)
                y[0] = eps
                y[1:] = np.fromiter(
                    (1 if labels[value] in chosen_set else -1 for value in range(Q)),
                    dtype=np.int8,
                    count=Q,
                )
                overlap = len(B & chosen_set)
                if len(B) == 0:
                    slack = phase
                elif phase == 0:
                    slack = overlap * (2 - overlap)
                else:
                    slack = (overlap - 1) ** 2
                coefficients = (
                    eps
                    * y[left].astype(np.int16)
                    * y[right].astype(np.int16)
                    * C[left, right]
                )
                rows.append(
                    {
                        "direction_index": int(record["direction_index"]),
                        "chosen_fibres": tuple(chosen_fibres),
                        "coefficients": coefficients,
                        "rhs": 3 + 2 * slack,
                    }
                )
    return rows


def direct_witness_audit(
    data: dict,
    case: dict,
    chosen_edges: list[list[int]],
    specs: list[dict],
    saturated_specs: list[dict],
    solver_means: list[int],
    solver_residues: dict[int, int],
) -> dict:
    edges = data["edges"]
    C = np.asarray(data["C"], dtype=np.int16)
    chosen = {tuple(int(value) for value in edge) for edge in chosen_edges}
    chosen_rows_were_unique = len(chosen) == len(chosen_edges)
    all_chosen_edges_exist = chosen <= set(edges)
    selected = np.asarray([int(edge in chosen) for edge in edges], dtype=np.int16)

    degrees = [0] * N
    for a, b in chosen:
        degrees[a] += 1
        degrees[b] += 1
    observed_boundary = tuple(index for index, degree in enumerate(degrees) if degree & 1)
    expected_boundary = tuple(int(value) for value in case["fixed_boundary_graph_vertices"])
    product = math.prod(int(C[a, b]) for a, b in chosen)

    normalized_score_support = {}
    shell_extrema = {}
    shell_valid = True
    for eps in (-1, 1):
        scores = np.asarray(data["features"][eps], dtype=np.int16) @ selected
        normalized = eps * scores
        normalized_score_support[str(eps)] = sorted(int(value) for value in np.unique(normalized))
        shell_extrema[str(eps)] = {
            "raw_min": int(scores.min()),
            "raw_max": int(scores.max()),
            "normalized_min": int(normalized.min()),
        }
        shell_valid &= bool(normalized.min() >= 3)

    branch_edge_results = []
    branch_valid = True
    for row in case["fixed_infinity_star_edges"]:
        edge = tuple(int(value) for value in row["edge"])
        observed = int(edge in chosen)
        expected = int(row["value"])
        branch_valid &= observed == expected
        branch_edge_results.append(
            {"edge": list(edge), "expected": expected, "observed": observed}
        )
    line_graph_vertices = {
        int(value) + 1 for value in case["representative_finite_field"]
    }
    present_line_star = {
        b for a, b in chosen if a == 0 and b in line_graph_vertices
    }
    present_outside_star = {
        b for a, b in chosen if a == 0 and b not in line_graph_vertices
    }
    if case["pointed_star_branch"] == "A":
        branch_semantics_valid = 1 in present_line_star
    else:
        branch_semantics_valid = not present_line_star and bool(present_outside_star)
    branch_valid &= branch_semantics_valid

    direction_rows = []
    computed_means = []
    direction_valid = True
    for spec, solver_mean in zip(specs, solver_means):
        computed = int(np.asarray(spec["coefficients"], dtype=np.int32) @ selected) - 21
        computed_means.append(computed)
        row_valid = bool(
            computed == int(solver_mean)
            and int(spec["floor"]) <= computed <= 32
            and computed % 2 == 0
        )
        direction_valid &= row_valid
        direction_rows.append(
            {
                "direction_index": int(spec["direction_index"]),
                "eps": int(spec["eps"]),
                "b": int(spec["b"]),
                "phase": int(spec["phase"]),
                "floor": int(spec["floor"]),
                "computed_scaled_mean": computed,
                "solver_scaled_mean": int(solver_mean),
                "valid": row_valid,
            }
        )

    type_rows = {}
    for eps in (-1, 1):
        values = [
            computed_means[index]
            for index, spec in enumerate(specs)
            if int(spec["eps"]) == eps
        ]
        residues = [value % 8 for value in values]
        solver_residue = int(solver_residues[eps])
        valid = bool(
            len(values) == 4
            and sum(values) == 32
            and len(set(residues)) == 1
            and residues[0] == solver_residue
        )
        direction_valid &= valid
        type_rows[str(eps)] = {
            "scaled_means": values,
            "sum": sum(values),
            "residues_mod_8": residues,
            "solver_common_residue": solver_residue,
            "valid": valid,
        }

    saturated_failures = []
    for index, row in enumerate(saturated_specs):
        lhs = int(np.asarray(row["coefficients"], dtype=np.int32) @ selected)
        if lhs != int(row["rhs"]):
            saturated_failures.append(
                {"constraint_index": index, "lhs": lhs, "rhs": int(row["rhs"])}
            )
    saturated_valid = not saturated_failures

    checks = {
        "chosen_rows_unique": chosen_rows_were_unique,
        "all_chosen_edges_exist": all_chosen_edges_exist,
        "exact_29_edges": len(chosen) == EDGE_COUNT,
        "exact_boundary_xor": observed_boundary == expected_boundary,
        "positive_paley_product": product == 1,
        "complete_full_eigenshell_inequalities": shell_valid,
        "pointed_branch_constraints": branch_valid,
        "direction_mean_identities_floors_budgets_residues": direction_valid,
        "saturated_affine_equalities": saturated_valid,
    }
    return {
        "valid": all(checks.values()),
        "checks": checks,
        "edge_count": len(chosen),
        "boundary": list(observed_boundary),
        "expected_boundary": list(expected_boundary),
        "c_H": product,
        "pointed_star_branch": case["pointed_star_branch"],
        "branch_fixed_edge_audit": branch_edge_results,
        "present_infinity_to_line_graph_vertices": sorted(present_line_star),
        "present_infinity_to_outside_graph_vertices": sorted(present_outside_star),
        "normalized_complete_shell_score_support": normalized_score_support,
        "complete_shell_extrema": shell_extrema,
        "direction_rows": direction_rows,
        "type_mean_audit": type_rows,
        "saturated_affine_equality_count": len(saturated_specs),
        "saturated_affine_equality_failures": saturated_failures,
    }


def solve_case(
    data: dict,
    case: dict,
    timeout: float,
    workers: int,
    seed: int,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    require(timeout > 0, "timeout must be positive")
    require(workers > 0, "workers must be positive")
    edges = data["edges"]
    signs = data["edge_signs"]
    edge_index = {edge: index for index, edge in enumerate(edges)}
    require(len(edge_index) == EXPECTED_EDGE_VARIABLES, "edge index changed")

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == EDGE_COUNT)

    for row in case["fixed_infinity_star_edges"]:
        edge = tuple(int(value) for value in row["edge"])
        require(edge in edge_index, f"pointed fixed edge is absent: {edge}")
        model.add(selected[edge_index[edge]] == int(row["value"]))

    boundary = tuple(int(value) for value in case["fixed_boundary_graph_vertices"])
    boundary_set = set(boundary)
    incident = [[] for _ in range(N)]
    for index, (a, b) in enumerate(edges):
        incident[a].append(selected[index])
        incident[b].append(selected[index])
    require(all(len(row) == N - 1 for row in incident), "complete-graph incidence changed")
    for vertex in range(N):
        if vertex in boundary_set:
            model.add_bool_xor(incident[vertex])
        else:
            model.add_bool_xor([~incident[vertex][0], *incident[vertex][1:]])

    negative = [selected[index] for index, sign in enumerate(signs) if int(sign) == -1]
    require(negative, "Paley graph has no negative edges")
    model.add_bool_xor([~negative[0], *negative[1:]])

    specs = direction_specs(data, case)
    means_by_type: dict[int, list] = {-1: [], 1: []}
    mean_variables = []
    for spec in specs:
        mean = model.new_int_var(
            int(spec["floor"]),
            32,
            f"scaled_mean_{spec['direction_index']}",
        )
        coefficients = spec["coefficients"]
        model.add(
            mean
            == sum(
                int(coefficient) * selected[index]
                for index, coefficient in enumerate(coefficients)
            )
            - 21
        )
        model.add_modulo_equality(0, mean, 2)
        means_by_type[int(spec["eps"])].append(mean)
        mean_variables.append(mean)

    common_residues = {}
    for eps in (-1, 1):
        require(len(means_by_type[eps]) == 4, "direction type is not four-by-four")
        model.add(sum(means_by_type[eps]) == 32)
        residue = model.new_int_var(0, 7, f"common_residue_{eps}")
        for mean in means_by_type[eps]:
            model.add_modulo_equality(residue, mean, 8)
        common_residues[eps] = residue

    complete_shell_constraints = {}
    for eps in (-1, 1):
        normalized = eps * np.asarray(data["features"][eps], dtype=np.int8)
        for row in normalized:
            bad = np.flatnonzero(row < 0).tolist()
            model.add(sum(selected[index] for index in bad) <= 13)
        complete_shell_constraints[str(eps)] = int(len(normalized))

    saturated_specs = saturated_affine_specs(data, specs)
    for row in saturated_specs:
        model.add(
            sum(
                int(coefficient) * selected[index]
                for index, coefficient in enumerate(row["coefficients"])
            )
            == int(row["rhs"])
        )

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(timeout)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.random_seed = int(seed)
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.linearization_level = 2
    status = solver.solve(model)
    status_name = solver.status_name(status)
    require(status != cp_model.MODEL_INVALID, "CP-SAT reported MODEL_INVALID")
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    infeasible = status == cp_model.INFEASIBLE

    if infeasible:
        rigorous_status = "rigorously_infeasible_exact_finite_model"
    elif feasible:
        rigorous_status = "feasible_direct_witness_audit_required"
    else:
        rigorous_status = "unknown_or_timeout_no_mathematical_conclusion"

    out = {
        "case_key": f"orbit{case['branch_orbit_index']}_branch{case['pointed_star_branch']}",
        **case,
        "c_H": 1,
        "edge_variables": len(selected),
        "exact_selected_edge_count": EDGE_COUNT,
        "direction_rows": [
            {
                "direction_index": int(row["direction_index"]),
                "direction": list(row["direction"]),
                "eps": int(row["eps"]),
                "B": sorted(int(value) for value in row["B"]),
                "b": int(row["b"]),
                "phase": int(row["phase"]),
                "floor": int(row["floor"]),
            }
            for row in specs
        ],
        "type_floor_sums": {
            str(eps): sum(int(row["floor"]) for row in specs if row["eps"] == eps)
            for eps in (-1, 1)
        },
        "complete_shell_constraints": complete_shell_constraints,
        "saturated_affine_slack_equalities": len(saturated_specs),
        "solver_status": status_name,
        "rigorous_status": rigorous_status,
        "feasible": feasible,
        "finite_infeasibility_certificate": infeasible,
        "case_decided": feasible or infeasible,
        "timeout_seconds": timeout,
        "workers": workers,
        "seed": seed,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "solver_wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        chosen = [
            list(edge)
            for edge, variable in zip(edges, selected)
            if solver.value(variable)
        ]
        solver_means = [int(solver.value(variable)) for variable in mean_variables]
        solver_residues = {
            eps: int(solver.value(variable)) for eps, variable in common_residues.items()
        }
        witness_audit = direct_witness_audit(
            data,
            case,
            chosen,
            specs,
            saturated_specs,
            solver_means,
            solver_residues,
        )
        require(witness_audit["valid"], "CP-SAT witness failed the direct full-shell audit")
        out["chosen_edges_H"] = chosen
        out["solver_scaled_means"] = solver_means
        out["solver_common_residues"] = {
            str(eps): value for eps, value in solver_residues.items()
        }
        out["witness_audit"] = witness_audit
        out["rigorous_status"] = "rigorously_feasible_direct_witness_audited"
    return out


def selected_cases(cases: list[dict], orbit: str, branch: str) -> list[dict]:
    selected = [
        row
        for row in cases
        if (orbit == "all" or int(orbit) == int(row["branch_orbit_index"]))
        and (branch == "ALL" or branch == row["pointed_star_branch"])
    ]
    require(selected, "orbit/branch selector chose no cases")
    return selected


def run(
    orbit: str,
    branch: str,
    timeout: float,
    workers: int,
    seed: int,
) -> dict:
    started = time.time()
    cases, orbit_source, normalization_audit = audited_pointed_cases()
    chosen_cases = selected_cases(cases, orbit, branch)
    data, geometry_audit = load_full_geometry()
    results = [
        solve_case(data, case, timeout, workers, seed + index)
        for index, case in enumerate(chosen_cases)
    ]

    full_scope = len(chosen_cases) == EXPECTED_CASES
    all_infeasible = all(row["finite_infeasibility_certificate"] for row in results)
    witnesses = [row["case_key"] for row in results if row["feasible"]]
    unknown = [row["case_key"] for row in results if not row["case_decided"]]
    if witnesses:
        status = "rigorous_positive_z7_full_shell_witness_found"
        conclusion = "positive z=7 is feasible in the exact full-shell model"
    elif full_scope and all_infeasible:
        status = "complete_rigorous_positive_z7_pointed_full_shell_exclusion"
        conclusion = "positive z=7 is excluded by all four exhaustive pointed cases"
    elif all_infeasible:
        status = "complete_rigorous_selected_pointed_case_exclusion"
        conclusion = "every selected pointed case is excluded; unselected cases remain"
    else:
        status = "exact_pointed_full_shell_run_with_unknown_cases"
        conclusion = "at least one selected case remains unknown; no exclusion follows"

    return {
        "experiment": "p7_infinity7_positive_z7_pointed_full_cpsat",
        "status": status,
        "rigorous_conclusion": conclusion,
        "p": P,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": P,
        "z": 7,
        "phase": 0,
        "model": {
            "edge_count": "exactly 29",
            "boundary": "native XOR at all 50 vertices",
            "paley_product": "positive; even parity of negative Paley edges",
            "eigenshells": "every unique edge-feature row from both complete cached shells",
            "direction_means": "exact identities, floors, type budgets, and common residues modulo eight",
            "mean_allocation_enumeration_used": False,
            "saturated_equalities": "all floor-saturated affine equalities valid for the selected boundary",
        },
        "orbit_source": orbit_source,
        "pointed_normalization_audit": normalization_audit,
        "full_geometry_audit": geometry_audit,
        "selection": {
            "orbit": orbit,
            "branch": branch,
            "selected_case_count": len(chosen_cases),
            "full_four_case_scope": full_scope,
        },
        "timeout_seconds_per_case": timeout,
        "workers_per_case": workers,
        "seed_base": seed,
        "all_selected_cases_decided": not unknown,
        "all_selected_cases_infeasible": all_infeasible,
        "full_z7_exclusion_proved": full_scope and all_infeasible,
        "direct_witness_case_keys": witnesses,
        "unknown_case_keys": unknown,
        "case_results": results,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orbit",
        choices=("all", "0", "1"),
        default="all",
        help="select one audited z=7 boundary orbit, or all two",
    )
    parser.add_argument(
        "--branch",
        type=str.upper,
        choices=("ALL", "A", "B"),
        default="ALL",
        help="select pointed star branch A or B, or both",
    )
    parser.add_argument(
        "--timeout",
        "--seconds",
        dest="timeout",
        type=float,
        default=600.0,
        help="CP-SAT time limit in seconds for each selected case",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15717001)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    result = run(args.orbit, args.branch, args.timeout, args.workers, args.seed)
    pointed.atomic_write(args.output, result)
    if not args.quiet:
        compact = {
            key: value
            for key, value in result.items()
            if key not in {"pointed_normalization_audit", "full_geometry_audit", "case_results"}
        }
        compact["case_summary"] = [
            {
                "case_key": row["case_key"],
                "solver_status": row["solver_status"],
                "rigorous_status": row["rigorous_status"],
                "solver_wall_time_seconds": row["solver_wall_time_seconds"],
            }
            for row in result["case_results"]
        ]
        print(json.dumps(compact, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
