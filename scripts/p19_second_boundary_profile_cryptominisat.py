#!/usr/bin/env python3
"""Native-XOR realizability test for the fourteen p=19 survivors.

This is the exact parity-profile problem underlying Proposition 15.689.
It uses both affine Radon equations

    r = A x,             x = A^T r  (mod 2),

where the second identity follows from ``A^T A=I+J`` and ``|x|=16``.
Direction weights and their phase histograms are imposed by guarded exact
cardinality automata.  A pair in a phase-zero b=0 direction is normalized
to field elements 0 and 1.

UNSAT excludes the boundary profile itself.  SAT gives an audited affine
boundary but does not by itself realize a residual edge lift.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import functools
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
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles  # noqa: E402


P = 19
SIZE = 16


def survivor_profiles() -> list[dict[str, object]]:
    rows = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 16
    ]
    if len(rows) != 14:
        raise ArithmeticError("the p=19 high-slack remainder changed")
    return rows


def solve(profile_index: int, seconds: float, threads: int) -> dict[str, object]:
    from pycryptosat import Solver

    profiles = survivor_profiles()
    if not 0 <= profile_index < len(profiles):
        raise ValueError("profile index must be in 0..13")
    profile = profiles[profile_index]
    directions = projective_directions(P)
    started = time.time()

    next_id = 0

    def new_var() -> int:
        nonlocal next_id
        next_id += 1
        return next_id

    point = [new_var() for _ in range(P * P)]
    parity = [[new_var() for _ in range(P)] for _ in directions]
    selectors: list[dict[int, int]] = []
    phases = []
    records = []
    normalized_index = None
    for index, direction in enumerate(directions):
        eps, labels = field_direction_data(P, direction)
        phase = 0 if eps == 1 else 1  # c_H=-1
        allowed = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"][str(phase)].items()
        }
        selectors.append({b: new_var() for b in allowed})
        phases.append(phase)
        records.append((direction, eps, labels))
        if labels[0] == labels[1]:
            if normalized_index is not None or eps != 1:
                raise ArithmeticError("normalized direction changed")
            normalized_index = index
    if normalized_index is None:
        raise ArithmeticError("normalized direction missing")
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
        """Guarded exact cardinality using a deterministic count automaton.

        The clauses are active when ``guard`` is false.  Counts above half
        are encoded through complemented literals to keep the automaton small.
        """
        nonlocal cardinality_count
        literals = [int(literal) for literal in literals]
        bound = int(bound)
        if bound > len(literals) // 2:
            literals = [-literal for literal in literals]
            bound = len(literals) - bound
        if not 0 <= bound <= len(literals):
            add_clauses([[]] if guard is None else [[int(guard)]])
            cardinality_count += 1
            return
        overflow = bound + 1
        states = [
            [new_var() for _ in range(overflow + 1)]
            for _ in range(len(literals) + 1)
        ]

        def gated(clause: list[int]) -> list[int]:
            return clause if guard is None else [int(guard), *clause]

        clauses: list[list[int]] = []
        for row in states:
            clauses.append(gated(list(row)))
            for first in range(len(row)):
                for second in range(first + 1, len(row)):
                    clauses.append(gated([-row[first], -row[second]]))
        clauses.append(gated([states[0][0]]))
        clauses.append(gated([states[-1][bound]]))
        for index, literal in enumerate(literals, start=1):
            previous = states[index - 1]
            current = states[index]
            for count in range(overflow + 1):
                incremented = min(overflow, count + 1)
                clauses.append(gated([-previous[count], literal, current[count]]))
                clauses.append(
                    gated([-previous[count], -literal, current[incremented]])
                )
        add_clauses(clauses)
        cardinality_count += 1

    add_exact(point, SIZE)
    add_unit(point[0])
    add_unit(point[1])

    incident_parities: list[list[int]] = [[] for _ in point]
    for index, (_direction, _eps, labels) in enumerate(records):
        fibres = []
        for fibre in range(P):
            fibre_points = [
                point[v] for v, label in enumerate(labels) if label == fibre
            ]
            fibres.append(fibre_points)
            solver.add_xor_clause(
                [*fibre_points, parity[index][fibre]],
                False,
            )
            xor_count += 1
        for v, label in enumerate(labels):
            incident_parities[v].append(parity[index][label])

        direction_selectors = list(selectors[index].values())
        add_clauses([direction_selectors])
        for first in range(len(direction_selectors)):
            for second in range(first + 1, len(direction_selectors)):
                add_clauses(
                    [[-direction_selectors[first], -direction_selectors[second]]]
                )
        for b, selector in selectors[index].items():
            add_exact(parity[index], b, guard=-selector)
            if b == SIZE:
                # With sixteen selected points, parity weight sixteen means
                # sixteen singleton fibres.  Expose that implication in CNF.
                for fibre_points in fibres:
                    for first in range(len(fibre_points)):
                        for second in range(first + 1, len(fibre_points)):
                            add_clauses(
                                [[-selector, -fibre_points[first], -fibre_points[second]]]
                            )

    for v, rows in enumerate(incident_parities):
        if len(rows) != P + 1:
            raise ArithmeticError("affine point degree changed")
        solver.add_xor_clause([point[v], *rows], False)
        xor_count += 1

    # Exact direction-weight histograms within each Paley phase.
    for phase in (0, 1):
        phase_profile = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"][str(phase)].items()
        }
        phase_indices = [i for i, value in enumerate(phases) if value == phase]
        if len(phase_indices) != 10:
            raise ArithmeticError("quadratic direction split changed")
        for b, target in phase_profile.items():
            add_exact([selectors[i][b] for i in phase_indices], target)

    add_unit(selectors[normalized_index][0])

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
        "experiment": "p19_second_boundary_profile_cryptominisat",
        "profile_index": profile_index,
        "pair_slack": int(profile["pair_slack"]),
        "phase_profiles_b": profile["phase_profiles_b"],
        "normalization": {
            "selected_points": [0, 1],
            "direction": list(directions[normalized_index]),
            "phase": 0,
            "b": 0,
            "c_H": -1,
        },
        "solver": "cryptominisat-native-xor",
        "solver_status": status,
        "feasible_boundary_profile": satisfiable is True,
        "finite_infeasibility_only": satisfiable is False,
        "semantic_variables": semantic_variables,
        "total_variables": next_id,
        "clauses": clause_count,
        "native_xor_constraints": xor_count,
        "cardinality_constraints": cardinality_count,
        "threads": threads,
        "build_seconds": build_seconds,
        "solve_seconds": time.time() - started - build_seconds,
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        chosen = [v for v, literal in enumerate(point) if assignment[literal]]
        observed = {0: {}, 1: {}}
        direction_rows = []
        valid = len(chosen) == SIZE and {0, 1} <= set(chosen)
        for index, (direction, eps, labels) in enumerate(records):
            counts = [0] * P
            for v in chosen:
                counts[labels[v]] += 1
            b = sum(count & 1 for count in counts)
            phase = phases[index]
            observed[phase][b] = observed[phase].get(b, 0) + 1
            direction_rows.append(
                {"direction": list(direction), "eps": eps, "phase": phase, "b": b}
            )
        expected = {
            phase: {
                int(b): int(count)
                for b, count in profile["phase_profiles_b"][str(phase)].items()
            }
            for phase in (0, 1)
        }
        valid = valid and observed == expected
        if not valid:
            raise AssertionError("CryptoMiniSat p=19 witness failed audit")
        result["boundary"] = chosen
        result["boundary_coordinates"] = [[v % P, v // P] for v in chosen]
        result["direction_rows"] = direction_rows
        result["witness_audit_valid"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--profile", type=int)
    group.add_argument("--all", action="store_true")
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.all:
        worker = functools.partial(
            solve, seconds=args.seconds, threads=args.threads
        )
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=max(1, int(args.jobs))
        ) as executor:
            result = list(executor.map(worker, range(14)))
    else:
        result = solve(args.profile, args.seconds, args.threads)
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
