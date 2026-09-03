#!/usr/bin/env python3
"""Exact p=31 decision model for the grouped branch-number lemma.

This is a single-prime attack on a new structural lemma, not a prime sweep.
It asks whether a nonzero binary word ``f`` can satisfy

    wt(f) + #{direction groups where M^T f is nonzero} <= 31.

Odd point weight is already covered by the exact radial-injection argument:
every silent direction has a distinct supported kernel class.  This model
therefore imposes even point weight and attacks only the genuinely open case.
An infeasibility result proves the desired lower bound 32 at p=31.  A feasible
result is printed with a support that is independently rechecked.
"""

from __future__ import annotations

import argparse
import json
import time
from itertools import product

from ortools.sat.python import cp_model

try:
    from e1_gmin_m4_symmetric_halved_mod2 import (
        _antipodal_classes,
        _evaluate,
        _directions,
    )
except ModuleNotFoundError:
    # Keep the decision model directly runnable on an isolated mesh node.
    def _antipodal_classes(p: int) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(
                {
                    min(point, ((-point[0]) % p, (-point[1]) % p))
                    for point in product(range(p), repeat=2)
                    if point != (0, 0)
                }
            )
        )

    def _directions(p: int) -> tuple[tuple[int, int], ...]:
        return tuple((1, slope) for slope in range(p)) + ((0, 1),)

    def _evaluate(
        p: int, functional: tuple[int, int], point: tuple[int, int]
    ) -> int:
        return (
            functional[0] * point[0] + functional[1] * point[1]
        ) % p


P = 31
H = (P - 1) // 2
D = P + 1
N = D * H


def block_layout() -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    """Return point blocks grouped by projective functional."""
    points = _antipodal_classes(P)
    squares = tuple(sorted({value * value % P for value in range(1, P)}))
    groups: list[tuple[tuple[int, ...], ...]] = []
    for functional in _directions(P):
        blocks = tuple(
            tuple(
                index
                for index, point in enumerate(points)
                if _evaluate(P, functional, point) ** 2 % P == square
            )
            for square in squares
        )
        if tuple(sorted(map(len, blocks))) != (P,) * H:
            raise ArithmeticError("affine block layout changed")
        groups.append(blocks)
    if len(points) != N or len(groups) != D:
        raise ArithmeticError("p=31 design dimensions changed")
    return tuple(groups), points


def verify_support(
    support: tuple[int, ...],
    groups: tuple[tuple[tuple[int, ...], ...], ...],
) -> dict[str, object]:
    selected = set(support)
    active_groups: list[int] = []
    syndromes: list[list[int]] = []
    for direction, blocks in enumerate(groups):
        syndrome = [sum(point in selected for point in block) & 1 for block in blocks]
        syndromes.append(syndrome)
        if any(syndrome):
            active_groups.append(direction)
    branch_weight = len(support) + len(active_groups)
    return {
        "support_indices": list(support),
        "support_weight": len(support),
        "active_direction_indices": active_groups,
        "active_direction_count": len(active_groups),
        "silent_direction_count": D - len(active_groups),
        "group_branch_weight": branch_weight,
        "counterexample": bool(support) and branch_weight <= P,
        "verified": bool(support),
    }


def solve(
    time_limit: float,
    workers: int,
    log_search: bool,
    min_point_weight: int,
    max_point_weight: int,
    symmetry_break: str,
) -> dict[str, object]:
    groups, _points = block_layout()
    model = cp_model.CpModel()
    point_bits = [model.new_bool_var(f"f_{index}") for index in range(N)]
    block_bits: list[list[cp_model.IntVar]] = []
    active_bits: list[cp_model.IntVar] = []

    if symmetry_break == "independent-pair":
        # A counterexample cannot be supported on one vector line: its unique
        # annihilating direction is silent and every other direction separates
        # its antipodal magnitudes, giving branch weight at least 32.  Hence it
        # contains an independent pair.  GL(2,31) moves that ordered pair to
        # [(0,1)],[(1,0)], indices 0,15 in this layout.
        model.add(point_bits[0] == 1)
        model.add(point_bits[15] == 1)
    for direction, blocks in enumerate(groups):
        direction_bits: list[cp_model.IntVar] = []
        for block_index, block in enumerate(blocks):
            parity = model.new_bool_var(f"y_{direction}_{block_index}")
            # AddBoolXOr enforces odd parity.  Appending not(parity) says
            # XOR(points in block) == parity.
            model.add_bool_xor([point_bits[index] for index in block] + [parity.Not()])
            direction_bits.append(parity)
        active = model.new_bool_var(f"active_{direction}")
        model.add_max_equality(active, direction_bits)
        block_bits.append(direction_bits)
        active_bits.append(active)

    if symmetry_break == "silent-triple":
        # A counterexample of even weight s>=2 has z>=s+1>=3 silent
        # directions.  PGL(2,31) is triply transitive on the projective line
        # and the grouped design is equivariant, so any ordered silent triple
        # may be moved to functionals (1,0),(1,1),(0,1), indices 0,1,31.
        model.add(active_bits[0] == 0)
        model.add(active_bits[1] == 0)
        model.add(active_bits[31] == 0)

    point_weight = sum(point_bits)
    active_weight = sum(active_bits)
    # The odd-weight case has a short symbolic proof, so spending solver time
    # on it would duplicate settled work.
    model.add_modulo_equality(0, point_weight, 2)
    model.add(point_weight >= min_point_weight)
    model.add(point_weight <= max_point_weight)
    model.add(point_weight + active_weight <= P)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.log_search_progress = log_search
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    started = time.monotonic()
    status = solver.solve(model)
    elapsed = time.monotonic() - started

    out: dict[str, object] = {
        "p": P,
        "decision_problem": "nonzero f with wt(f)+active_groups<=31",
        "scope": "even point weight only; odd weight is symbolically closed",
        "symmetry_break": symmetry_break,
        "symmetry_justification": (
            "a counterexample contains an independent support pair"
            if symmetry_break == "independent-pair"
            else "a counterexample has at least three silent directions; "
            "PGL(2,31) is triply transitive"
        ),
        "status": solver.status_name(status),
        "wall_seconds": elapsed,
        "solver_wall_seconds": solver.wall_time,
        "workers": workers,
        "time_limit_seconds": time_limit,
        "point_weight_range": [min_point_weight, max_point_weight],
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "response_stats": solver.response_stats(),
        "group_branch_at_least_32_proved": status == cp_model.INFEASIBLE,
    }
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        support = tuple(
            index for index, bit in enumerate(point_bits) if solver.value(bit)
        )
        check = verify_support(support, groups)
        if not check["counterexample"]:
            raise ArithmeticError("solver witness failed independent verification")
        out["witness"] = check
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--log-search", action="store_true")
    parser.add_argument("--min-point-weight", type=int, default=2)
    parser.add_argument("--max-point-weight", type=int, default=30)
    parser.add_argument(
        "--symmetry-break",
        choices=("independent-pair", "silent-triple"),
        default="silent-triple",
    )
    args = parser.parse_args()
    if args.time_limit <= 0 or args.workers <= 0:
        parser.error("time limit and workers must be positive")
    if not 2 <= args.min_point_weight <= args.max_point_weight <= 30:
        parser.error("require 2 <= min point weight <= max point weight <= 30")
    print(
        json.dumps(
            solve(
                args.time_limit,
                args.workers,
                args.log_search,
                args.min_point_weight,
                args.max_point_weight,
                args.symmetry_break,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
