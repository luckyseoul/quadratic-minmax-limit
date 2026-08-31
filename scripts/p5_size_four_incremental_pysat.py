#!/usr/bin/env python3
"""Incremental native-cardinality SAT sweep of p=5 size-four orbits.

One solver instance holds edge cardinality and all requested Max-shell score
cuts.  XOR chains expose each vertex-degree parity and the Paley-product
parity as assumption literals.  Boundary-orbit cases are then solved
incrementally, retaining learned clauses across cases and writing progress
atomically after every result.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def add_xor_output(solver, literals: list[int], top_id: int) -> tuple[int, int]:
    if not literals:
        raise ValueError("XOR output needs at least one literal")
    current = literals[0]
    for literal in literals[1:]:
        top_id += 1
        result = top_id
        solver.add_clause([current, literal, -result])
        solver.add_clause([current, -literal, result])
        solver.add_clause([-current, literal, result])
        solver.add_clause([-current, -literal, -result])
        current = result
    return current, top_id


def run(
    orbit_path: Path,
    output_path: Path,
    shell_mode: str,
    seconds_per_case: float,
    seed_order: int,
) -> dict:
    from pysat.solvers import Solver

    source = json.loads(orbit_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("this incremental sweep is scoped to p=5")
    c_h = int(source["c_H"])
    infinity_value = int(source["infinity_value"])
    data = geometry(5, shell_mode)
    n = int(data["n"])
    edges = data["edges"]
    signs = data["edge_signs"]
    edge_literals = list(range(1, len(edges) + 1))
    top_id = len(edges)
    solver = Solver(name="minicard")
    solver.add_atmost(edge_literals, 21)
    solver.add_atmost([-literal for literal in edge_literals], len(edges) - 21)
    solver.add_clause([edge_literals[edges.index((0, 1))]])

    boundary_outputs = []
    for vertex in range(n):
        incident = [
            edge_literals[j] for j, (a, b) in enumerate(edges) if vertex in (a, b)
        ]
        output, top_id = add_xor_output(solver, incident, top_id)
        boundary_outputs.append(output)
    negative = [edge_literals[j] for j, sign in enumerate(signs) if sign == -1]
    product_output, top_id = add_xor_output(solver, negative, top_id)

    n_score_constraints = 0
    for eps in (-1, 1):
        for feature in data["features"][eps]:
            bad = [
                edge_literals[j] for j in np.flatnonzero(eps * feature < 0).tolist()
            ]
            solver.add_atmost(bad, 9)
            n_score_constraints += 1

    prior_rows = []
    if output_path.exists():
        prior = json.loads(output_path.read_text())
        if (
            int(prior["c_H"]) != c_h
            or int(prior["infinity_value"]) != infinity_value
            or prior["shell_mode"] != shell_mode
        ):
            raise ValueError("existing output belongs to another sweep")
        prior_rows = list(prior["rows"])
    completed = {int(row["orbit_index"]) for row in prior_rows}
    rows = {int(row["orbit_index"]): row for row in prior_rows}
    cases = list(enumerate(source["orbits"]))
    # A deterministic stride changes which geometry feeds early learned clauses.
    if cases:
        offset = seed_order % len(cases)
        cases = cases[offset:] + cases[:offset]
    started = time.time()

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        counts = {
            status: sum(row["solver_status"] == status for row in ordered)
            for status in ("UNSATISFIABLE", "SATISFIABLE", "UNKNOWN")
        }
        return {
            "experiment": "p5_size_four_incremental_pysat",
            "status": "incremental_native_cardinality_orbit_sweep",
            "source": str(orbit_path),
            "p": 5,
            "c_H": c_h,
            "infinity_value": infinity_value,
            "shell_mode": shell_mode,
            "source_orbit_count": len(cases),
            "completed": len(ordered),
            "pending": len(cases) - len(ordered),
            "status_counts": counts,
            "all_unsatisfiable": len(ordered) == len(cases)
            and counts["UNSATISFIABLE"] == len(cases),
            "n_edge_variables": len(edges),
            "n_total_variables": top_id,
            "n_score_constraints": n_score_constraints,
            "seconds_per_case": seconds_per_case,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output_path, snapshot())
    for index, orbit in cases:
        if index in completed:
            continue
        boundary = tuple(int(v) for v in orbit["representative_vertices"])
        boundary_set = set(boundary)
        assumptions = [
            output if vertex in boundary_set else -output
            for vertex, output in enumerate(boundary_outputs)
        ]
        assumptions.append(product_output if c_h == -1 else -product_output)
        before = solver.accum_stats().copy()
        case_started = time.time()
        timer = threading.Timer(float(seconds_per_case), solver.interrupt)
        timer.start()
        try:
            result = solver.solve_limited(
                assumptions=assumptions,
                expect_interrupt=True,
            )
        finally:
            timer.cancel()
        after = solver.accum_stats().copy()
        row = {
            "orbit_index": index,
            "boundary": list(boundary),
            "orbit_size": int(orbit["size"]),
            "solver_status": (
                "SATISFIABLE" if result is True else "UNSATISFIABLE"
                if result is False
                else "UNKNOWN"
            ),
            "finite_infeasibility_certificate": result is False,
            "elapsed_seconds": time.time() - case_started,
            "stats_delta": {
                key: int(after.get(key, 0) - before.get(key, 0)) for key in after
            },
        }
        if result is True:
            assignment = {literal for literal in solver.get_model() if literal > 0}
            chosen = [
                list(edge)
                for edge, literal in zip(edges, edge_literals)
                if literal in assignment
            ]
            row["chosen_edges_H"] = chosen
            row["witness_audit"] = verify_witness(
                5,
                c_h,
                chosen,
                infinity_value,
                boundary,
                shell_mode,
            )
            if not row["witness_audit"]["valid"]:
                raise AssertionError("incremental SAT witness failed audit")
        rows[index] = row
        atomic_write(output_path, snapshot())
        solver.clear_interrupt()
    out = snapshot()
    atomic_write(output_path, out)
    solver.delete()
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbits", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shell-mode", choices=("affine", "full"), default="full")
    parser.add_argument("--seconds-per-case", type=float, default=30.0)
    parser.add_argument("--seed-order", type=int, default=0)
    args = parser.parse_args()
    out = run(
        args.orbits,
        args.output,
        args.shell_mode,
        args.seconds_per_case,
        args.seed_order,
    )
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
