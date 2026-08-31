#!/usr/bin/env python3
"""Incremental full-shell SAT shard for arbitrary p=5 boundary orbits.

One native-cardinality solver holds the exact 21-edge model, all 260 full
eigenshell score cuts, vertex-degree XOR outputs, and Paley-product XOR.
Boundary orbit cases in ``[start, stop)`` are supplied as assumptions so
learned clauses persist within a shard.  SAT witnesses are independently
audited; timeout is always reported as UNKNOWN.
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

from p5_size_four_incremental_pysat import add_xor_output  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry, verify_witness  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run(
    source_path: Path,
    output_path: Path,
    start: int,
    stop: int | None,
    seconds_per_case: float,
    seed_order: int,
) -> dict:
    from pysat.solvers import Solver

    source = json.loads(source_path.read_text())
    if int(source["p"]) != 5:
        raise ValueError("this solver is scoped to p=5")
    all_orbits = list(source["orbits"])
    stop = len(all_orbits) if stop is None else min(int(stop), len(all_orbits))
    start = max(0, int(start))
    if not start <= stop:
        raise ValueError("invalid orbit range")
    source_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    c_h = int(source["c_H"])
    infinity_value = int(source["infinity_value"])
    data = geometry(5, "full")
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
            edge_literals[index]
            for index, (a, b) in enumerate(edges)
            if vertex in (a, b)
        ]
        output, top_id = add_xor_output(solver, incident, top_id)
        boundary_outputs.append(output)
    negative = [
        edge_literals[index] for index, sign in enumerate(signs) if sign == -1
    ]
    product_output, top_id = add_xor_output(solver, negative, top_id)

    score_constraints = 0
    for eps in (-1, 1):
        for feature in data["features"][eps]:
            bad = [
                edge_literals[index]
                for index in np.flatnonzero(eps * feature < 0).tolist()
            ]
            solver.add_atmost(bad, 9)
            score_constraints += 1

    scope = list(range(start, stop))
    if scope:
        offset = seed_order % len(scope)
        scope = scope[offset:] + scope[:offset]
    prior_rows = []
    if output_path.exists():
        prior = json.loads(output_path.read_text())
        if (
            prior.get("source_sha256") != source_sha256
            or int(prior["start_orbit"]) != start
            or int(prior["stop_orbit"]) != stop
        ):
            raise ValueError("existing output belongs to another shard")
        prior_rows = list(prior.get("rows", []))
    rows = {int(row["orbit_index"]): row for row in prior_rows}
    started = time.time()

    def snapshot() -> dict:
        ordered = [rows[index] for index in sorted(rows)]
        status_counts = {
            status: sum(row["solver_status"] == status for row in ordered)
            for status in ("UNSATISFIABLE", "SATISFIABLE", "UNKNOWN")
        }
        return {
            "experiment": "p5_full_shell_pysat_shard",
            "status": "incremental_exact_full_shell_boundary_sat",
            "source": str(source_path),
            "source_sha256": source_sha256,
            "p": 5,
            "c_H": c_h,
            "boundary_size": int(source.get("boundary_size", 4)),
            "infinity_value": infinity_value,
            "start_orbit": start,
            "stop_orbit": stop,
            "scope_orbits": stop - start,
            "completed": len(ordered),
            "pending": stop - start - len(ordered),
            "status_counts": status_counts,
            "all_unsatisfiable": len(ordered) == stop - start
            and status_counts["UNSATISFIABLE"] == stop - start,
            "edge_variables": len(edges),
            "total_sat_variables": top_id,
            "score_constraints": score_constraints,
            "seconds_per_case": seconds_per_case,
            "seed_order": seed_order,
            "elapsed_seconds_this_run": time.time() - started,
            "rows": ordered,
        }

    atomic_write(output_path, snapshot())
    for index in scope:
        if index in rows and rows[index]["solver_status"] != "UNKNOWN":
            continue
        orbit = all_orbits[index]
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
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
                "SATISFIABLE"
                if result is True
                else "UNSATISFIABLE"
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
                "full",
            )
            if not row["witness_audit"]["valid"]:
                raise AssertionError("SAT witness failed independent audit")
        rows[index] = row
        atomic_write(output_path, snapshot())
        solver.clear_interrupt()
    result = snapshot()
    atomic_write(output_path, result)
    solver.delete()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--seconds-per-case", type=float, default=10.0)
    parser.add_argument("--seed-order", type=int, default=0)
    args = parser.parse_args()
    result = run(
        args.source,
        args.output,
        args.start,
        args.stop,
        args.seconds_per_case,
        args.seed_order,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
