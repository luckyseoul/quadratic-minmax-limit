#!/usr/bin/env python3
"""Enumerate exact mass-ten quadratic lift vectors on J(7,4)."""
from __future__ import annotations

import argparse
import itertools
import json
import time
from collections import Counter
from pathlib import Path

from p7_quadratic_lift_classify import MASS, N, PAIRS, POINTS


def enumerate_lifts(
    seconds: float,
    limit: int,
    positive_values: tuple[int, ...] | None,
    save_vectors: Path | None,
) -> dict:
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()
    values = [model.new_int_var(0, 2, f"B_{a}") for a in range(len(POINTS))]
    model.add(sum(values) == MASS)
    vertex_masses = []
    for i in range(N):
        u = model.new_int_var(0, MASS, f"U_{i}")
        model.add(u == sum(values[a] for a, X in enumerate(POINTS) if i in X))
        vertex_masses.append(u)
    pair_masses = {}
    for i, j in PAIRS:
        t = model.new_int_var(0, MASS, f"T_{i}_{j}")
        model.add(
            t == sum(values[a] for a, X in enumerate(POINTS) if i in X and j in X)
        )
        pair_masses[i, j] = t
    for a, X in enumerate(POINTS):
        model.add(
            6 * values[a]
            == 2 * sum(pair_masses[pair] for pair in itertools.combinations(X, 2))
            - 3 * sum(vertex_masses[i] for i in X)
            + 36
        )
    if positive_values is not None:
        counts = Counter(positive_values)
        if sum(positive_values) != MASS:
            raise ValueError("histogram values must sum to ten")
        for value in (1, 2):
            bits = []
            for a, B in enumerate(values):
                bit = model.new_bool_var(f"is_{value}_{a}")
                model.add(B == value).only_enforce_if(bit)
                model.add(B != value).only_enforce_if(~bit)
                bits.append(bit)
            model.add(sum(bits) == counts[value])

    class Collector(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.count = 0
            self.histograms = Counter()
            self.vectors = []

        def on_solution_callback(self):
            vector = tuple(self.value(B) for B in values)
            self.count += 1
            self.histograms[tuple(sorted(v for v in vector if v))] += 1
            if save_vectors is not None:
                self.vectors.append(vector)
            if limit and self.count >= limit:
                self.stop_search()

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.symmetry_level = 0
    collector = Collector()
    started = time.time()
    status = solver.solve(model, collector)
    if save_vectors is not None:
        import numpy as np

        np.savez_compressed(save_vectors, values=np.asarray(collector.vectors, dtype=np.int8))
    return {
        "experiment": "p7_quadratic_lift_enumerate",
        "positive_values": list(positive_values) if positive_values is not None else None,
        "solver_status": solver.status_name(status),
        "complete": status in (cp_model.OPTIMAL, cp_model.INFEASIBLE),
        "solution_count": collector.count,
        "limit": limit,
        "histogram_counts": {
            ",".join(map(str, key)): value
            for key, value in sorted(collector.histograms.items())
        },
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "vectors_file": str(save_vectors) if save_vectors is not None else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--positive-values", type=int, nargs="+")
    parser.add_argument("--save-vectors", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = enumerate_lifts(
        args.seconds,
        args.limit,
        tuple(sorted(args.positive_values)) if args.positive_values else None,
        args.save_vectors,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
