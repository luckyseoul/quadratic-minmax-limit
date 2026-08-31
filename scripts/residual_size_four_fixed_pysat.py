#!/usr/bin/env python3
"""Native cardinality-SAT model for one fixed size-four boundary.

MiniCard receives the score bounds and edge cardinality natively.  Boundary
and edge-product parities are encoded by linear-size Tseitin XOR chains.
This backend is independent of CP-SAT, HiGHS, and SCIP.
"""
from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def add_xor_chain(solver, literals: list[int], target: int, top_id: int) -> int:
    """Add XOR(literals)=target and return the last allocated variable."""
    if not literals:
        if target:
            solver.add_clause([])
        return top_id
    current = literals[0]
    for literal in literals[1:]:
        top_id += 1
        result = top_id
        # result <-> current XOR literal.
        solver.add_clause([current, literal, -result])
        solver.add_clause([current, -literal, result])
        solver.add_clause([-current, literal, result])
        solver.add_clause([-current, -literal, -result])
        current = result
    solver.add_clause([current if target else -current])
    return top_id


def solve_case(
    p: int,
    c_h: int,
    fixed_boundary: tuple[int, ...],
    shell_mode: str,
    seconds: float,
) -> dict:
    from pysat.solvers import Solver

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

    edge_literals = list(range(1, len(edges) + 1))
    top_id = len(edges)
    solver = Solver(name="minicard")
    solver.add_atmost(edge_literals, 4 * p + 1)
    solver.add_atmost([-literal for literal in edge_literals], len(edges) - (4 * p + 1))
    solver.add_clause([edge_literals[edges.index((0, 1))]])

    fixed_set = set(fixed_boundary)
    for vertex in range(n):
        incident = [
            edge_literals[j] for j, (a, b) in enumerate(edges) if vertex in (a, b)
        ]
        top_id = add_xor_chain(solver, incident, int(vertex in fixed_set), top_id)
    negative = [
        edge_literals[j] for j, sign in enumerate(signs) if sign == -1
    ]
    top_id = add_xor_chain(solver, negative, int(c_h == -1), top_id)

    score_limit = 2 * p - 1
    n_score_constraints = 0
    for eps in (-1, 1):
        for feature in data["features"][eps]:
            bad = [
                edge_literals[j] for j in np.flatnonzero(eps * feature < 0).tolist()
            ]
            solver.add_atmost(bad, score_limit)
            n_score_constraints += 1

    timer = threading.Timer(float(seconds), solver.interrupt)
    timer.start()
    try:
        satisfiable = solver.solve_limited(expect_interrupt=True)
    finally:
        timer.cancel()
    out = {
        "experiment": "residual_size_four_fixed_pysat",
        "status": "exact_fixed_boundary_native_cardinality_sat",
        "p": p,
        "c_H": c_h,
        "fixed_boundary": list(fixed_boundary),
        "shell_mode": shell_mode,
        "solver_status": (
            "SATISFIABLE" if satisfiable is True else "UNSATISFIABLE"
            if satisfiable is False
            else "UNKNOWN"
        ),
        "feasible": satisfiable is True,
        "finite_infeasibility_certificate": satisfiable is False,
        "n_edge_variables": len(edges),
        "n_total_variables": top_id,
        "n_score_constraints": n_score_constraints,
        "accumulated_stats": solver.accum_stats(),
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        assignment = {literal for literal in solver.get_model() if literal > 0}
        chosen = [
            list(edge)
            for edge, literal in zip(edges, edge_literals)
            if literal in assignment
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
        if not out["witness_audit"]["valid"]:
            raise AssertionError("SAT witness failed independent audit")
    solver.delete()
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
