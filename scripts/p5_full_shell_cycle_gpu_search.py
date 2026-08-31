#!/usr/bin/env python3
"""GPU four-cycle local search for a p=5 full-shell residual witness.

Every move toggles a four-cycle containing exactly two selected edges.  It
therefore preserves edge count and every vertex-degree parity.  Restricting
to cycles with even Paley-negative parity also preserves ``c_H``; cycles
through the distinguished edge are omitted.  The search minimizes complete
eigenshell score violations in large GPU batches.  Finding penalty zero is
a genuine audited witness.  Failure is heuristic only.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_fixed_boundary_cpsat import audit_witness, shell_rows  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def admissible_cycles(edges: list[tuple[int, int]], signs: np.ndarray) -> np.ndarray:
    edge_index = {edge: index for index, edge in enumerate(edges)}
    distinguished = edge_index[(0, 1)]
    cycles = set()
    for vertices in itertools.combinations(range(26), 4):
        a, b, c, d = vertices
        for ordered in ((a, b, c, d), (a, b, d, c), (a, c, b, d)):
            row = tuple(
                sorted(
                    edge_index[tuple(sorted((ordered[index], ordered[(index + 1) % 4])))]
                    for index in range(4)
                )
            )
            if distinguished in row:
                continue
            if sum(int(signs[index] == -1) for index in row) & 1:
                continue
            cycles.add(row)
    return np.asarray(sorted(cycles), dtype=np.int32)


def initial_graph(boundary: tuple[int, ...], c_h: int, seed: int) -> np.ndarray:
    """Find one graph satisfying only cardinality and parity invariants."""
    from ortools.sat.python import cp_model

    data = geometry(5, "full")
    edges = data["edges"]
    signs = data["edge_signs"]
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"edge_{a}_{b}") for a, b in edges]
    model.add(sum(selected) == 21)
    model.add(selected[edges.index((0, 1))] == 1)
    boundary_set = set(boundary)
    for vertex in range(26):
        incident = [
            selected[index]
            for index, edge in enumerate(edges)
            if vertex in edge
        ]
        if vertex in boundary_set:
            model.add_bool_xor(incident)
        else:
            model.add_bool_xor([~incident[0], *incident[1:]])
    negative = [selected[index] for index, sign in enumerate(signs) if sign == -1]
    if c_h == -1:
        model.add_bool_xor(negative)
    else:
        model.add_bool_xor([~negative[0], *negative[1:]])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 10.0
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = int(seed)
    solver.parameters.randomize_search = True
    status = solver.solve(model)
    if solver.status_name(status) not in {"OPTIMAL", "FEASIBLE"}:
        raise RuntimeError("failed to construct parity-feasible initial graph")
    return np.asarray([solver.value(variable) for variable in selected], dtype=np.int8)


def penalty(counts) -> object:
    import cupy as cp

    excess = cp.maximum(counts - 9, 0).astype(cp.int32)
    return cp.sum(100 * excess * excess + excess, axis=-1)


def search(
    source_path: Path,
    orbit_index: int,
    seconds: float,
    batch_size: int,
    seed: int,
    restarts: int,
) -> dict:
    import cupy as cp

    started = time.time()
    source = json.loads(source_path.read_text())
    orbit = source["orbits"][orbit_index]
    boundary = tuple(sorted(int(value) for value in orbit["representative_vertices"]))
    c_h = int(source["c_H"])
    data = geometry(5, "full")
    edges = data["edges"]
    signs = data["edge_signs"]
    bad = np.concatenate([(shell_rows(eps)[1] < 0).astype(np.int8) for eps in (-1, 1)])
    cycles_np = admissible_cycles(edges, signs)
    bad_gpu = cp.asarray(bad, dtype=cp.int16)
    cycles = cp.asarray(cycles_np)
    rng = cp.random.RandomState(seed)
    best_penalty = math.inf
    best_x = None
    best_counts = None
    iterations = 0
    accepted = 0
    completed_restarts = 0

    for restart in range(restarts):
        if time.time() - started >= seconds:
            break
        x = cp.asarray(initial_graph(boundary, c_h, seed + restart), dtype=cp.int16)
        counts = bad_gpu @ x
        current = int(penalty(counts).item())
        if current < best_penalty:
            best_penalty = current
            best_x = cp.asnumpy(x).astype(np.int8)
            best_counts = cp.asnumpy(counts).astype(np.int16)
        completed_restarts += 1
        stagnant = 0
        while time.time() - started < seconds and stagnant < 3000:
            iterations += 1
            sample_indices = rng.randint(0, len(cycles_np), size=batch_size, dtype=cp.int32)
            sample = cycles[sample_indices]
            old = x[sample]
            valid = cp.sum(old, axis=1) == 2
            deltas = cp.sum(
                bad_gpu[:, sample].transpose(1, 0, 2)
                * (1 - 2 * old)[:, None, :],
                axis=2,
            )
            trial_counts = counts[None, :] + deltas
            scores = penalty(trial_counts)
            scores = cp.where(valid, scores, cp.iinfo(cp.int32).max)
            position = int(cp.argmin(scores).item())
            candidate = int(scores[position].item())
            temperature = max(0.25, 20.0 * math.exp(-stagnant / 500.0))
            accept = candidate <= current
            if not accept and candidate < cp.iinfo(cp.int32).max:
                accept = bool(cp.random.random() < math.exp((current - candidate) / temperature))
            if accept:
                move = sample[position]
                x[move] = 1 - x[move]
                counts = trial_counts[position]
                current = candidate
                accepted += 1
                stagnant = 0 if candidate < best_penalty else stagnant + 1
                if candidate < best_penalty:
                    best_penalty = candidate
                    best_x = cp.asnumpy(x).astype(np.int8)
                    best_counts = cp.asnumpy(counts).astype(np.int16)
                    if candidate == 0:
                        chosen = [list(edge) for edge, bit in zip(edges, best_x) if bit]
                        witness = audit_witness(data, c_h, boundary, chosen)
                        if not witness["valid"]:
                            raise AssertionError("GPU zero-penalty graph failed audit")
                        return {
                            "experiment": "p5_full_shell_cycle_gpu_search",
                            "found": True,
                            "source": str(source_path),
                            "orbit_index": orbit_index,
                            "boundary": list(boundary),
                            "c_H": c_h,
                            "best_penalty": 0,
                            "chosen_edges_H": chosen,
                            "witness_audit": witness,
                            "admissible_four_cycles": len(cycles_np),
                            "iterations": iterations,
                            "accepted_moves": accepted,
                            "restarts": completed_restarts,
                            "elapsed_seconds": time.time() - started,
                        }
            else:
                stagnant += 1

    return {
        "experiment": "p5_full_shell_cycle_gpu_search",
        "found": False,
        "source": str(source_path),
        "orbit_index": orbit_index,
        "boundary": list(boundary),
        "c_H": c_h,
        "best_penalty": int(best_penalty),
        "best_max_bad_count": int(best_counts.max()) if best_counts is not None else None,
        "best_violated_shell_rows": (
            int(np.count_nonzero(best_counts > 9)) if best_counts is not None else None
        ),
        "best_edges_H": (
            [list(edge) for edge, bit in zip(edges, best_x) if bit]
            if best_x is not None
            else []
        ),
        "admissible_four_cycles": len(cycles_np),
        "iterations": iterations,
        "accepted_moves": accepted,
        "restarts": completed_restarts,
        "elapsed_seconds": time.time() - started,
        "not_a_certificate": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--seed", type=int, default=15657001)
    parser.add_argument("--restarts", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = search(
        args.source,
        args.orbit_index,
        args.seconds,
        args.batch_size,
        args.seed,
        args.restarts,
    )
    atomic_write(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key not in {"chosen_edges_H", "best_edges_H"}}, indent=2))


if __name__ == "__main__":
    main()
