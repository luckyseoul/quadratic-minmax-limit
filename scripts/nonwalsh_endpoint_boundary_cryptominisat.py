#!/usr/bin/env python3
"""Native-XOR endpoint-boundary diagnostic.

Test whether a noncollinear set ``S`` of ``p-2`` points in ``AG(2,p)`` can
have, in every projective direction, either one or ``p-2`` affine fibres
meeting ``S`` oddly.  A noncollinear triple is normalized to
``(0,0),(1,0),(0,1)`` under ``AGL(2,p)``.

This is the same finite-geometry question as
``nonwalsh_endpoint_boundary_cpsat.py``, but it gives the fibre equations to
CryptoMiniSat as native XOR clauses.  Guarded cardinality encodings choose
between the two allowed parity weights.  UNSAT is a finite-prime certificate;
SAT is only a boundary witness and says nothing by itself about a residual
edge set.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)


def solve_endpoint_case(
    p: int,
    seconds: float,
    threads: int,
) -> dict[str, object]:
    from pycryptosat import Solver
    if p < 7 or p % 2 == 0:
        raise ValueError("p must be odd and at least seven")

    started = time.time()
    size = p - 2
    directions = projective_directions(p)

    # Allocate every semantic variable before cardinality encodings introduce
    # auxiliaries.  Variables are one-indexed for DIMACS/PySAT.
    next_id = 0

    def new_var() -> int:
        nonlocal next_id
        next_id += 1
        return next_id

    point = [new_var() for _ in range(p * p)]
    parity = [[new_var() for _ in range(p)] for _ in directions]
    high = [new_var() for _ in directions]
    semantic_variables = next_id

    solver = Solver(verbose=0, threads=max(1, int(threads)))
    clause_count = 0
    xor_count = 0
    cardinality_count = 0

    def add_clauses(clauses: list[list[int]]) -> None:
        nonlocal clause_count
        if clauses:
            solver.add_clauses(clauses)
            clause_count += len(clauses)

    def add_unit(literal: int) -> None:
        add_clauses([[int(literal)]])

    def add_exact(literals: list[int], bound: int, guard: int | None = None) -> None:
        """Guarded exact cardinality via a deterministic count automaton.

        State ``overflow=bound+1`` absorbs every larger count.  One-hot state
        rows plus the two transitions for each input literal make the final
        exact-count unit clause equivalent to ``sum(literals)==bound``.  This
        avoids a runtime dependency on PySAT in the system Python that owns
        the working CryptoMiniSat extension.
        """
        nonlocal cardinality_count
        bound = int(bound)
        if not 0 <= bound <= len(literals):
            add_clauses([[]] if guard is None else [[int(guard)]])
            cardinality_count += 1
            return
        overflow = bound + 1
        states = [
            [new_var() for _ in range(overflow + 1)]
            for _ in range(len(literals) + 1)
        ]

        def guarded(clause: list[int]) -> list[int]:
            return clause if guard is None else [int(guard), *clause]

        clauses: list[list[int]] = []
        for row in states:
            clauses.append(guarded(list(row)))
            for first in range(len(row)):
                for second in range(first + 1, len(row)):
                    clauses.append(guarded([-row[first], -row[second]]))
        clauses.append(guarded([states[0][0]]))
        clauses.append(guarded([states[-1][bound]]))
        for index, literal in enumerate(literals, start=1):
            previous = states[index - 1]
            current = states[index]
            for count in range(overflow + 1):
                incremented = min(overflow, count + 1)
                clauses.append(
                    guarded([-previous[count], int(literal), current[count]])
                )
                clauses.append(
                    guarded(
                        [-previous[count], -int(literal), current[incremented]]
                    )
                )
        add_clauses(clauses)
        cardinality_count += 1

    add_exact(point, size)
    for u in (0, 1, p):
        add_unit(point[u])

    records = []
    for index, direction in enumerate(directions):
        eps, labels = field_direction_data(p, direction)
        fibres = [
            [point[u] for u in range(p * p) if labels[u] == fibre]
            for fibre in range(p)
        ]
        for fibre in range(p):
            # XOR(points in fibre) == parity[index][fibre].
            solver.add_xor_clause(
                [*fibres[fibre], parity[index][fibre]],
                False,
            )
            xor_count += 1

        # high=0 activates weight one; high=1 activates weight p-2.
        add_exact(parity[index], 1, guard=high[index])
        add_exact(parity[index], size, guard=-high[index])

        # In the high branch, parity weight equals |S|, hence every occupied
        # fibre is a singleton.  These redundant guarded clauses materially
        # strengthen propagation before the XOR system is fully assigned.
        for fibre_points in fibres:
            for first in range(len(fibre_points)):
                for second in range(first + 1, len(fibre_points)):
                    add_clauses(
                        [[-high[index], -fibre_points[first], -fibre_points[second]]]
                    )

        # Each normalized pair collides in one of these three directions, so
        # those directions cannot be injective.
        if index in (0, 1, p):
            add_unit(-high[index])
        records.append((direction, eps, labels))

    build_seconds = time.time() - started
    satisfiable, assignment = solver.solve(time_limit=float(seconds))
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    result: dict[str, object] = {
        "experiment": "nonwalsh_endpoint_boundary_cryptominisat",
        "p": p,
        "finite_boundary_points": size,
        "normalization": [0, 1, p],
        "normalization_coordinates": [[0, 0], [1, 0], [0, 1]],
        "allowed_odd_fibre_counts": [1, size],
        "solver": "cryptominisat-native-xor",
        "solver_status": status,
        "feasible_noncollinear_endpoint_set": satisfiable is True,
        "finite_infeasibility_only": satisfiable is False,
        "semantic_variables": semantic_variables,
        "total_variables": next_id,
        "clauses": clause_count,
        "native_xor_constraints": xor_count,
        "cardinality_constraints": cardinality_count,
        "cardinality_encoding": "guarded_deterministic_count_automaton",
        "threads": threads,
        "build_seconds": build_seconds,
        "solve_seconds": time.time() - started - build_seconds,
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        chosen = [u for u, literal in enumerate(point) if assignment[literal]]
        direction_rows = []
        valid = len(chosen) == size and {0, 1, p} <= set(chosen)
        for index, (direction, eps, labels) in enumerate(records):
            counts = [0] * p
            for u in chosen:
                counts[labels[u]] += 1
            odd_fibres = [fibre for fibre, count in enumerate(counts) if count & 1]
            valid = valid and len(odd_fibres) in (1, size)
            direction_rows.append(
                {
                    "direction": list(direction),
                    "eps": eps,
                    "odd_fibre_count": len(odd_fibres),
                    "odd_fibres": odd_fibres,
                }
            )
        if not valid:
            raise AssertionError("CryptoMiniSat endpoint witness failed audit")
        result["points"] = [[u % p, u // p] for u in chosen]
        result["directions"] = direction_rows
        result["witness_audit_valid"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    result = solve_endpoint_case(
        args.p,
        args.seconds,
        args.threads,
    )
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
