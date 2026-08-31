#!/usr/bin/env python3
"""Independent exact CNF solver for one p=5 full-shell boundary case.

This is a second encoding of the model in
``p5_full_shell_fixed_boundary_cpsat.py``.  Sequential counters impose the
21-edge equality and every normalized shell score bound.  Explicit XOR
chains impose vertex degrees, the Paley product, and each shell bad-count
parity.  Modern CDCL solvers can therefore attack hard CP-SAT tails without
sharing implementation machinery.  SAT witnesses are independently
audited; UNSAT is an exact finite CNF certificate result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

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


def add_xor_value(clauses: list[list[int]], literals: list[int], value: int, top_id: int) -> int:
    if not literals or value not in (0, 1):
        raise ValueError("XOR needs literals and a Boolean value")
    current = literals[0]
    for literal in literals[1:]:
        top_id += 1
        result = top_id
        clauses.extend(
            (
                [current, literal, -result],
                [current, -literal, result],
                [-current, literal, result],
                [-current, -literal, -result],
            )
        )
        current = result
    clauses.append([current if value else -current])
    return top_id


def solve_case(
    source_path: Path,
    orbit_index: int,
    solver_name: str,
    seconds: float,
) -> dict:
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Solver

    started = time.time()
    source = json.loads(source_path.read_text())
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    boundary_set = set(boundary)
    c_h = int(source["c_H"])
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    edge_literals = list(range(1, len(edges) + 1))
    top_id = len(edge_literals)
    clauses: list[list[int]] = []
    native_cardinality = solver_name == "minicard"
    native_atmosts: list[tuple[list[int], int]] = []

    if native_cardinality:
        native_atmosts.append((edge_literals, 21))
        native_atmosts.append(([-literal for literal in edge_literals], len(edges) - 21))
    else:
        edge_count = CardEnc.equals(
            lits=edge_literals,
            bound=21,
            top_id=top_id,
            encoding=EncType.seqcounter,
        )
        clauses.extend(edge_count.clauses)
        top_id = edge_count.nv
    clauses.append([edge_literals[edges.index((0, 1))]])

    for vertex in range(int(data["n"])):
        incident = [
            edge_literals[index]
            for index, (a, b) in enumerate(edges)
            if vertex in (a, b)
        ]
        top_id = add_xor_value(
            clauses,
            incident,
            int(vertex in boundary_set),
            top_id,
        )
    negative = [
        edge_literals[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    top_id = add_xor_value(clauses, negative, int(c_h == -1), top_id)

    score_constraints = 0
    shell_parity_constraints = 0
    shell_metadata = {}
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        shell_metadata[str(eps)] = {
            "parity_mass": int(parity.sum()),
            "lift_mass": int((78 - int(parity.sum())) // 2),
        }
        for row_index, row in enumerate(normalized):
            bad = [
                edge_literals[index]
                for index in np.flatnonzero(row < 0).tolist()
            ]
            if native_cardinality:
                native_atmosts.append((bad, 9))
            else:
                atmost = CardEnc.atmost(
                    lits=bad,
                    bound=9,
                    top_id=top_id,
                    encoding=EncType.seqcounter,
                )
                clauses.extend(atmost.clauses)
                top_id = atmost.nv
            score_constraints += 1
            # bad = 9 - P - 2L, hence bad parity is 1-P.
            top_id = add_xor_value(
                clauses,
                bad,
                1 - int(parity[row_index]),
                top_id,
            )
            shell_parity_constraints += 1

    build_seconds = time.time() - started
    solver = Solver(name=solver_name, bootstrap_with=clauses)
    for literals, bound in native_atmosts:
        solver.add_atmost(literals, bound)
    timer = threading.Timer(float(seconds), solver.interrupt)
    solve_started = time.time()
    timer.start()
    try:
        result = solver.solve_limited(expect_interrupt=True)
    finally:
        timer.cancel()
    status = "SATISFIABLE" if result is True else "UNSATISFIABLE" if result is False else "UNKNOWN"
    output = {
        "experiment": "p5_full_shell_fixed_boundary_cnf",
        "status": "independent_exact_cnf_edge_score_and_parity_model",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "c_H": c_h,
        "solver": solver_name,
        "solver_status": status,
        "finite_infeasibility_certificate": result is False,
        "feasible": result is True,
        "edge_variables": len(edge_literals),
        "total_variables": top_id,
        "clauses": len(clauses),
        "native_cardinality_constraints": len(native_atmosts),
        "score_constraints": score_constraints,
        "shell_parity_constraints": shell_parity_constraints,
        "shells": shell_metadata,
        "build_seconds": build_seconds,
        "solve_seconds": time.time() - solve_started,
        "elapsed_seconds": time.time() - started,
        "stats": {key: int(value) for key, value in solver.accum_stats().items()},
    }
    if result is True:
        assignment = {literal for literal in solver.get_model() if literal > 0}
        chosen_edges = [
            list(edge)
            for edge, literal in zip(edges, edge_literals)
            if literal in assignment
        ]
        output["chosen_edges_H"] = chosen_edges
        output["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if not output["witness_audit"]["valid"]:
            raise AssertionError("CNF witness failed independent audit")
    solver.delete()
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--solver", default="kissat404")
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve_case(args.source, args.orbit_index, args.solver, args.seconds)
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "chosen_edges_H"}, indent=2))


if __name__ == "__main__":
    main()
