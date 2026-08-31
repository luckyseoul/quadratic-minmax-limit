#!/usr/bin/env python3
"""HiGHS MILP cross-check for one fixed residual size-four boundary.

This is the same direct ``H`` model as
``residual_boundary_four_lift_cpsat.py`` with parity linearized by integer
quotient variables.  It is useful as an independent backend when CP-SAT's
XOR/cardinality search stalls.  An `INFEASIBLE` result is a finite branch
certificate; a feasible result is audited directly against the requested
affine or full Max shells.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def solve_case(
    p: int,
    c_h: int,
    fixed_boundary: tuple[int, ...],
    shell_mode: str,
    seconds: float,
) -> dict:
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    if p not in (5, 7) or c_h not in (-1, 1):
        raise ValueError("need p in {5,7} and c_h in {+-1}")
    started = time.time()
    fixed_boundary = tuple(sorted(fixed_boundary))
    data = geometry(p, shell_mode)
    n = int(data["n"])
    edges = data["edges"]
    signs = data["edge_signs"]
    if len(fixed_boundary) != 4 or len(set(fixed_boundary)) != 4:
        raise ValueError("fixed_boundary must have four distinct vertices")
    if not all(0 <= vertex < n for vertex in fixed_boundary):
        raise ValueError("fixed boundary vertex outside graph")

    n_edges = len(edges)
    degree_offset = n_edges
    sign_quotient = n_edges + n
    n_variables = sign_quotient + 1
    score_rows = [
        (eps, row)
        for eps in (-1, 1)
        for row in data["features"][eps]
    ]
    n_constraints = 2 + n + 1 + len(score_rows)
    A = lil_matrix((n_constraints, n_variables), dtype=np.float64)
    lower = np.full(n_constraints, -np.inf, dtype=np.float64)
    upper = np.full(n_constraints, np.inf, dtype=np.float64)
    row_index = 0

    A[row_index, :n_edges] = 1
    lower[row_index] = upper[row_index] = 4 * p + 1
    row_index += 1
    distinguished = edges.index((0, 1))
    A[row_index, distinguished] = 1
    lower[row_index] = upper[row_index] = 1
    row_index += 1

    fixed_set = set(fixed_boundary)
    for vertex in range(n):
        for edge_index, (a, b) in enumerate(edges):
            if vertex in (a, b):
                A[row_index, edge_index] = 1
        A[row_index, degree_offset + vertex] = -2
        lower[row_index] = upper[row_index] = int(vertex in fixed_set)
        row_index += 1

    for edge_index, sign in enumerate(signs):
        if sign == -1:
            A[row_index, edge_index] = 1
    A[row_index, sign_quotient] = -2
    lower[row_index] = upper[row_index] = int(c_h == -1)
    row_index += 1

    score_limit = 2 * p - 1
    for eps, feature in score_rows:
        bad = np.flatnonzero(eps * feature < 0)
        A[row_index, bad] = 1
        upper[row_index] = score_limit
        row_index += 1
    if row_index != n_constraints:
        raise AssertionError("constraint count mismatch")

    variable_lower = np.zeros(n_variables, dtype=np.float64)
    variable_upper = np.ones(n_variables, dtype=np.float64)
    variable_upper[degree_offset:sign_quotient] = (n - 1) // 2
    variable_upper[sign_quotient] = n_edges // 2
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
    out = {
        "experiment": "residual_size_four_fixed_highs",
        "status": "exact_fixed_boundary_edge_model",
        "p": p,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "shell_mode": shell_mode,
        "solver_status": int(result.status),
        "message": str(result.message),
        "feasible": feasible,
        "finite_infeasibility_certificate": int(result.status) == 2,
        "n_variables": n_variables,
        "n_constraints": n_constraints,
        "mip_node_count": getattr(result, "mip_node_count", None),
        "mip_gap": getattr(result, "mip_gap", None),
        "elapsed_seconds": time.time() - started,
    }
    if result.x is not None:
        chosen = [
            list(edge)
            for edge, value in zip(edges, result.x[:n_edges])
            if value > 0.5
        ]
        out["chosen_edges_H"] = chosen
        out["witness_audit"] = verify_witness(
            p,
            c_h,
            chosen,
            int(0 in fixed_boundary),
            fixed_boundary,
            shell_mode,
        )
        if int(result.status) == 0 and not out["witness_audit"]["valid"]:
            raise AssertionError("HiGHS witness failed independent audit")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, choices=(5, 7), required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--shell-mode", choices=("affine", "full"), default="affine")
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = solve_case(
        args.p,
        args.c_h,
        tuple(args.fixed_boundary),
        args.shell_mode,
        args.seconds,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
