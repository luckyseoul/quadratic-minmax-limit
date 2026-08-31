#!/usr/bin/env python3
"""Exact cutting-plane CP-SAT for hard p=5 fixed-boundary cases."""
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
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import (  # noqa: E402
    add_lex_leader,
    audit_witness,
    boundary_edge_stabilizers,
    parity_vector,
    shell_rows,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def solve(
    source_path: Path,
    orbit_index: int,
    seconds: float,
    workers: int,
    round_seconds: float,
    cut_batch: int,
    seed: int,
    symmetry_breaking: bool,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    source = json.loads(source_path.read_text())
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    boundary_set = set(boundary)
    c_h = int(source["c_H"])
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]

    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 21)
    model.add(selected[edges.index((0, 1))] == 1)
    for vertex in range(26):
        literals = [
            selected[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        if vertex in boundary_set:
            model.add_bool_xor(literals)
        else:
            model.add_bool_xor([~literals[0], *literals[1:]])
    negative = [
        selected[index]
        for index, (a, b) in enumerate(edges)
        if int(C[a, b]) == -1
    ]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])

    stabilizers = boundary_edge_stabilizers(boundary) if symmetry_breaking else ()
    for index, permutation in enumerate(stabilizers):
        add_lex_leader(model, selected, permutation, index)

    normalized_rows = []
    for eps in (-1, 1):
        representatives, normalized = shell_rows(eps)
        parity = parity_vector(representatives, eps, c_h, boundary)
        for row_index, row in enumerate(normalized):
            bad_indices = np.flatnonzero(row < 0).tolist()
            bad_literals = [selected[index] for index in bad_indices]
            # Add every exact parity immediately; cardinality caps are the
            # separated cuts.
            if int(parity[row_index]) == 0:
                model.add_bool_xor(bad_literals)
            else:
                model.add_bool_xor([~bad_literals[0], *bad_literals[1:]])
            normalized_rows.append((eps, row_index, bad_indices))

    added: set[tuple[int, int]] = set()
    rounds = []
    final_status = "UNKNOWN"
    chosen_edges = None
    while time.time() - started < seconds:
        remaining = seconds - (time.time() - started)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = min(float(round_seconds), remaining)
        solver.parameters.num_search_workers = int(workers)
        solver.parameters.random_seed = int(seed + len(rounds))
        solver.parameters.cp_model_presolve = True
        solver.parameters.symmetry_level = 3
        status = solver.solve(model)
        status_name = solver.status_name(status)
        record = {
            "round": len(rounds) + 1,
            "solver_status": status_name,
            "cuts_before": len(added),
            "conflicts": int(solver.num_conflicts),
            "branches": int(solver.num_branches),
            "wall_time_seconds": float(solver.wall_time),
        }
        rounds.append(record)
        if status_name == "INFEASIBLE":
            final_status = "INFEASIBLE"
            break
        if status_name not in {"OPTIMAL", "FEASIBLE"}:
            final_status = status_name
            break
        values = np.asarray([solver.value(variable) for variable in selected], dtype=np.int16)
        violations = []
        for eps, row_index, bad_indices in normalized_rows:
            count = int(values[bad_indices].sum())
            if count > 9:
                violations.append((count - 9, eps, row_index, bad_indices))
        record["violated_rows"] = len(violations)
        record["maximum_excess"] = max((row[0] for row in violations), default=0)
        if not violations:
            final_status = "FEASIBLE"
            chosen_edges = [list(edge) for edge, value in zip(edges, values) if value]
            break
        violations.sort(reverse=True, key=lambda item: item[0])
        cuts_added = 0
        for _excess, eps, row_index, bad_indices in violations:
            key = (eps, row_index)
            if key in added:
                continue
            model.add(sum(selected[index] for index in bad_indices) <= 9)
            added.add(key)
            cuts_added += 1
            if cuts_added >= cut_batch:
                break
        record["cuts_added"] = cuts_added
        if cuts_added == 0:
            final_status = "SEPARATION_STALLED"
            break

    result = {
        "experiment": "p5_full_shell_fixed_boundary_cut",
        "status": "exact_native_xor_full_shell_cutting_plane",
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "orbit_index": orbit_index,
        "orbit_size": int(orbit["size"]),
        "boundary": list(boundary),
        "c_H": c_h,
        "solver_status": final_status,
        "finite_infeasibility_certificate": final_status == "INFEASIBLE",
        "feasible": final_status == "FEASIBLE",
        "cuts_final": len(added),
        "rounds": rounds,
        "workers": workers,
        "round_seconds": round_seconds,
        "cut_batch": cut_batch,
        "symmetry_breaking": symmetry_breaking,
        "boundary_stabilizer_size": len(stabilizers) if symmetry_breaking else None,
        "elapsed_seconds": time.time() - started,
    }
    if chosen_edges is not None:
        result["chosen_edges_H"] = chosen_edges
        result["witness_audit"] = audit_witness(data, c_h, boundary, chosen_edges)
        if not result["witness_audit"]["valid"]:
            raise AssertionError("cutting-plane witness failed independent audit")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--round-seconds", type=float, default=10.0)
    parser.add_argument("--cut-batch", type=int, default=16)
    parser.add_argument("--seed", type=int, default=15657001)
    parser.add_argument("--symmetry-breaking", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = solve(
        args.source,
        args.orbit_index,
        args.seconds,
        args.workers,
        args.round_seconds,
        args.cut_batch,
        args.seed,
        args.symmetry_breaking,
    )
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != "chosen_edges_H"}, indent=2))


if __name__ == "__main__":
    main()
