#!/usr/bin/env python3
"""Classify mass-ten nonnegative integral quadratics on J(7,4).

For a degree-at-most-two function on the middle slice, write

    B(X) = sum_{ij subset X} c_ij.

The 35-by-21 pair-incidence matrix has full column rank.  If ``sum_X B=10``
then ``sum_ij c_ij=1``.  Put

    U_i  = sum_{X contains i} B(X),
    T_ij = sum_{X contains {i,j}} B(X).

Inverting the pair-incidence Gram matrix gives

    6 c_ij = 2 T_ij - U_i - U_j + 6,

and hence the exact reconstruction test

    6 B(X) = 2 sum_{ij subset X} T_ij
             - 3 sum_{i in X} U_i + 36.

Thus bounded integer variables for the 35 values, together with these
identities, exactly describe the nonnegative integer-valued quadratic lifts
of exceptional mass 10 in the balanced p=7 residual profile.

This is a finite classification utility.  It is not by itself an exclusion
of the residual edge profile.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from fractions import Fraction
from pathlib import Path


N = 7
K = 4
MASS = 10
POINTS = tuple(itertools.combinations(range(N), K))
PAIRS = tuple(itertools.combinations(range(N), 2))


def integer_partitions(total: int, length: int, minimum: int = 1):
    """Yield nondecreasing positive partitions of total with fixed length."""
    if length == 0:
        if total == 0:
            yield ()
        return
    for first in range(minimum, total // length + 1):
        for tail in integer_partitions(total - first, length - 1, first):
            yield (first, *tail)


def solve_histogram(
    positive_values: tuple[int, ...] | None,
    seconds: float,
    workers: int,
    seed: int,
    optimize: str | None = None,
) -> dict:
    from ortools.sat.python import cp_model

    started = time.time()
    model = cp_model.CpModel()
    values = [model.new_int_var(0, MASS, f"B_{''.join(map(str, X))}") for X in POINTS]
    model.add(sum(values) == MASS)

    vertex_masses = []
    for i in range(N):
        u = model.new_int_var(0, MASS, f"U_{i}")
        model.add(u == sum(values[a] for a, X in enumerate(POINTS) if i in X))
        vertex_masses.append(u)

    pair_masses = []
    for i, j in PAIRS:
        t = model.new_int_var(0, MASS, f"T_{i}_{j}")
        model.add(
            t
            == sum(
                values[a]
                for a, X in enumerate(POINTS)
                if i in X and j in X
            )
        )
        pair_masses.append(t)

    pair_index = {pair: a for a, pair in enumerate(PAIRS)}
    for a, X in enumerate(POINTS):
        model.add(
            6 * values[a]
            == 2
            * sum(pair_masses[pair_index[pair]] for pair in itertools.combinations(X, 2))
            - 3 * sum(vertex_masses[i] for i in X)
            + 36
        )

    nonzero = []
    for a, value in enumerate(values):
        used = model.new_bool_var(f"used_{a}")
        model.add(value >= 1).only_enforce_if(used)
        model.add(value == 0).only_enforce_if(~used)
        nonzero.append(used)
    support = model.new_int_var(0, MASS, "support")
    model.add(support == sum(nonzero))

    if positive_values is not None:
        if not positive_values or sum(positive_values) != MASS:
            raise ValueError("positive values must be a positive partition of 10")
        model.add(support == len(positive_values))
        counts = {v: positive_values.count(v) for v in set(positive_values)}
        for v in range(1, MASS + 1):
            equal = []
            for a, value in enumerate(values):
                bit = model.new_bool_var(f"is_{v}_{a}")
                model.add(value == v).only_enforce_if(bit)
                model.add(value != v).only_enforce_if(~bit)
                equal.append(bit)
            model.add(sum(equal) == counts.get(v, 0))

    if optimize == "min-support":
        model.minimize(support)
    elif optimize == "max-support":
        model.maximize(support)
    elif optimize is not None:
        raise ValueError(f"unknown optimization mode: {optimize}")

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    status = solver.solve(model)
    feasible = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    out = {
        "experiment": "p7_quadratic_lift_classify",
        "positive_values": list(positive_values) if positive_values is not None else None,
        "optimize": optimize,
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "optimal": status == cp_model.OPTIMAL,
        "best_objective_bound": solver.best_objective_bound,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if feasible:
        b = [solver.value(value) for value in values]
        u = [solver.value(value) for value in vertex_masses]
        t = [solver.value(value) for value in pair_masses]
        coefficients = {
            f"{i}{j}": str(Fraction(2 * t[a] - u[i] - u[j] + 6, 6))
            for a, (i, j) in enumerate(PAIRS)
        }
        out.update(
            {
                "support": sum(value > 0 for value in b),
                "value_histogram": {
                    str(value): b.count(value) for value in sorted(set(b))
                },
                "nonzero_values": [
                    {"point": list(X), "value": b[a]}
                    for a, X in enumerate(POINTS)
                    if b[a]
                ],
                "vertex_masses": u,
                "pair_masses": {
                    f"{i}{j}": t[a] for a, (i, j) in enumerate(PAIRS)
                },
                "pair_coefficients": coefficients,
            }
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-values", type=int, nargs="+")
    parser.add_argument("--list-partitions", action="store_true")
    parser.add_argument("--scan-support", action="store_true")
    parser.add_argument("--support", type=int)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15649)
    parser.add_argument("--optimize", choices=("min-support", "max-support"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.list_partitions:
        lengths = [args.support] if args.support is not None else range(1, MASS + 1)
        print(json.dumps([list(v) for length in lengths for v in integer_partitions(MASS, length)]))
        return

    if args.scan_support:
        if args.support is None:
            raise ValueError("--scan-support requires --support")
        rows = [
            solve_histogram(
                partition,
                args.seconds,
                args.workers,
                args.seed + index,
                None,
            )
            for index, partition in enumerate(integer_partitions(MASS, args.support))
        ]
        rendered = json.dumps(
            {
                "experiment": "p7_quadratic_lift_histogram_scan",
                "support": args.support,
                "rows": rows,
                "all_decided": all(row["solver_status"] != "UNKNOWN" for row in rows),
            },
            indent=2,
        )
        print(rendered, flush=True)
        if args.output is not None:
            args.output.write_text(rendered + "\n")
        return

    positive_values = (
        tuple(sorted(args.positive_values)) if args.positive_values is not None else None
    )
    out = solve_histogram(
        positive_values,
        args.seconds,
        args.workers,
        args.seed,
        args.optimize,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
