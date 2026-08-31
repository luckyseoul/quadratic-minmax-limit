#!/usr/bin/env python3
"""Independent HiGHS MILP for one fixed p=5 full-shell boundary.

The selected variables are the 325 edges of ``K_26``.  Integer quotient
variables linearize the degree and Paley-product parities, while one
integer lift per antipodal full-shell row imposes

    bad_count(y) + 2 L(y) = 9 - P(y),    0 <= L(y) <= 4.

The six exact directional mean identities from Proposition 15.632 are
included as redundant strengthening.  A HiGHS ``Infeasible`` result is an
exact finite fixed-boundary exclusion; a feasible edge set is independently
audited with the CP-SAT model's direct graph/shell checker.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from p5_full_shell_fixed_boundary_cpsat import (  # noqa: E402
    audit_witness,
    parity_vector,
    shell_rows,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def solve_case(source_path: Path, orbit_index: int, seconds: float) -> dict:
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    started = time.time()
    source = json.loads(source_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("source must be p=5")
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    if len(boundary) % 2 or len(set(boundary)) != len(boundary):
        raise ValueError("boundary must have distinct vertices and even size")
    boundary_set = set(boundary)
    c_h = int(source["c_H"])

    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    n = int(data["n"])
    n_edges = len(edges)

    degree_offset = n_edges
    sign_quotient = degree_offset + n
    lift_offset = sign_quotient + 1
    n_lifts = 260
    direction_offset = lift_offset + n_lifts
    n_directions = 6
    n_variables = direction_offset + n_directions

    shell_data: list[tuple[int, int, list[int], int]] = []
    shell_metadata: dict[str, dict[str, int]] = {}
    lift_index = 0
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        parity_mass = int(parity.sum())
        shell_metadata[str(eps)] = {
            "rows": 130,
            "parity_mass": parity_mass,
            "lift_mass": (78 - parity_mass) // 2,
        }
        for row_index, row in enumerate(normalized):
            shell_data.append(
                (
                    eps,
                    row_index,
                    np.flatnonzero(row < 0).tolist(),
                    int(parity[row_index]),
                )
            )
            lift_index += 1
    if lift_index != n_lifts:
        raise AssertionError("unexpected full-shell row count")

    direction_data = [
        field_direction_data(5, direction) for direction in projective_directions(5)
    ]
    finite_boundary = tuple(vertex - 1 for vertex in boundary if vertex != 0)
    infinity_value = int(0 in boundary_set)
    direction_rows = []
    direction_coefficients = []
    for eps, labels in direction_data:
        counts = [0] * 5
        for value in finite_boundary:
            counts[labels[value]] += 1
        odd_fibres = sum(count & 1 for count in counts)
        sign = -eps * c_h
        if infinity_value:
            sign *= eps
        if odd_fibres & 1:
            sign *= -1
        phase = int(sign == -1)
        floor = scaled_direction_floor(5, odd_fibres, phase)
        coefficients = []
        for a, b in edges:
            if a == 0:
                coefficient = 1
            else:
                la, lb = labels[a - 1], labels[b - 1]
                coefficient = 5 if la == lb else -eps * int(C[a, b])
            coefficients.append(coefficient)
        direction_coefficients.append(coefficients)
        direction_rows.append(
            {"eps": eps, "odd_fibres": odd_fibres, "phase": phase, "floor": floor}
        )

    # Global count, fixed edge, degree parities, product parity, 260 shell
    # lift rows, two shell lift sums, six directional means, and two type sums.
    n_constraints = 2 + n + 1 + n_lifts + 2 + n_directions + 2
    A = lil_matrix((n_constraints, n_variables), dtype=np.float64)
    lower = np.full(n_constraints, -np.inf, dtype=np.float64)
    upper = np.full(n_constraints, np.inf, dtype=np.float64)
    row = 0

    A[row, :n_edges] = 1
    lower[row] = upper[row] = 21
    row += 1
    A[row, edges.index((0, 1))] = 1
    lower[row] = upper[row] = 1
    row += 1

    for vertex in range(n):
        for edge_index, edge in enumerate(edges):
            if vertex in edge:
                A[row, edge_index] = 1
        A[row, degree_offset + vertex] = -2
        lower[row] = upper[row] = int(vertex in boundary_set)
        row += 1

    for edge_index, (a, b) in enumerate(edges):
        if int(C[a, b]) == -1:
            A[row, edge_index] = 1
    A[row, sign_quotient] = -2
    lower[row] = upper[row] = int(c_h == -1)
    row += 1

    for local_lift, (_eps, _shell_row, bad_indices, parity) in enumerate(shell_data):
        A[row, bad_indices] = 1
        A[row, lift_offset + local_lift] = 2
        lower[row] = upper[row] = 9 - parity
        row += 1

    shell_start = 0
    for eps in (-1, 1):
        A[row, lift_offset + shell_start : lift_offset + shell_start + 130] = 1
        lift_mass = shell_metadata[str(eps)]["lift_mass"]
        lower[row] = upper[row] = lift_mass
        shell_start += 130
        row += 1

    for direction_index, coefficients in enumerate(direction_coefficients):
        A[row, :n_edges] = coefficients
        A[row, direction_offset + direction_index] = -2
        lower[row] = upper[row] = 15
        row += 1

    for eps in (-1, 1):
        for direction_index, record in enumerate(direction_rows):
            if int(record["eps"]) == eps:
                A[row, direction_offset + direction_index] = 1
        lower[row] = upper[row] = 9
        row += 1
    if row != n_constraints:
        raise AssertionError("constraint count mismatch")

    variable_lower = np.zeros(n_variables, dtype=np.float64)
    variable_upper = np.ones(n_variables, dtype=np.float64)
    variable_upper[degree_offset:sign_quotient] = 12
    variable_upper[sign_quotient] = 65
    variable_upper[lift_offset:direction_offset] = 4
    for direction_index, record in enumerate(direction_rows):
        variable_lower[direction_offset + direction_index] = int(record["floor"]) // 2
        variable_upper[direction_offset + direction_index] = 105

    result = milp(
        np.zeros(n_variables, dtype=np.float64),
        integrality=np.ones(n_variables, dtype=np.uint8),
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(A.tocsr(), lower, upper),
        options={
            "time_limit": float(seconds),
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )
    feasible = result.x is not None and int(result.status) == 0
    output = {
        "experiment": "p5_full_shell_fixed_boundary_highs",
        "status": "independent_exact_edge_parity_and_shell_lift_milp",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "c_H": c_h,
        "solver_status": int(result.status),
        "message": str(result.message),
        "finite_infeasibility_certificate": int(result.status) == 2,
        "feasible": feasible,
        "edge_variables": n_edges,
        "integer_variables": n_variables - n_edges,
        "constraints": n_constraints,
        "shells": shell_metadata,
        "direction_rows": direction_rows,
        "mip_node_count": getattr(result, "mip_node_count", None),
        "mip_gap": getattr(result, "mip_gap", None),
        "elapsed_seconds": time.time() - started,
    }
    if result.x is not None:
        chosen_edges = [
            list(edge)
            for edge, value in zip(edges, result.x[:n_edges])
            if value > 0.5
        ]
        output["chosen_edges_H"] = chosen_edges
        output["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if feasible and not output["witness_audit"]["valid"]:
            raise AssertionError("HiGHS witness failed independent audit")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve_case(args.source, args.orbit_index, args.seconds)
    atomic_write(args.output, result)
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "chosen_edges_H"},
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
