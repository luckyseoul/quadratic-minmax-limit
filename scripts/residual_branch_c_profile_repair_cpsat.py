#!/usr/bin/env python3
"""Repair a near-feasible p=31 Mobius parallel design with minimum half edits.

The input is the final JSON object emitted by
``p31_mobius_parallel_profile_search.py``.  Its ``goal_raw`` vector is held
exactly, auxiliary projective directions remain all different, and CP-SAT
minimizes the number of target halves whose (auxiliary direction, relative
scale) option changes.  This deliberately searches for a new exact design;
it does not optimize centers or claim a common-graph completion.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    _parallel_formula,
    paley_direction_sign,
)
from io_atomic import write_json_atomic  # noqa: E402


P = 31


def _scale(vector: tuple[int, int], scalar: int) -> tuple[int, int]:
    return scalar * vector[0] % P, scalar * vector[1] % P


def _last_json(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    starts = [index for index, char in enumerate(text) if char == "{"]
    for start in reversed(starts):
        try:
            value = json.loads(text[start:])
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise ValueError(f"no trailing JSON object in {path}")


def solve(
    near: dict[str, object], *, seconds: float, workers: int, seed: int
) -> dict[str, object]:
    directions = tuple(projective_functionals(P))
    signs = tuple(paley_direction_sign(P, direction) for direction in directions)
    hard = tuple(index for index, sign in enumerate(signs) if sign == 1)
    hinted_aux = tuple(int(value) for value in near["auxiliary_direction_indices"])
    hinted_scales = tuple(int(value) for value in near["relative_scales"])
    goal = tuple(int(value) for value in near["goal_raw"])
    if not (len(hard) == len(hinted_aux) == len(hinted_scales) == 16):
        raise ValueError("the near design must provide sixteen target options")
    if len(goal) != 32:
        raise ValueError("goal_raw must have 32 entries")

    options: list[list[tuple[int, int, tuple[int, ...]]]] = []
    for target_index in hard:
        target = directions[target_index]
        rows = []
        for auxiliary_index, auxiliary_direction in enumerate(directions):
            if auxiliary_index == target_index:
                continue
            for relative_scale in range(1, P):
                auxiliary = _scale(auxiliary_direction, relative_scale)
                profile = tuple(
                    _parallel_formula(P, target, auxiliary, row)
                    for row in directions
                )
                rows.append((auxiliary_index, relative_scale, profile))
        options.append(rows)

    model = cp_model.CpModel()
    choose: list[list[cp_model.IntVar]] = []
    by_auxiliary: list[list[cp_model.IntVar]] = [[] for _ in directions]
    by_direction: list[list[tuple[int, cp_model.IntVar]]] = [
        [] for _ in directions
    ]
    unchanged = []
    for target_pos, rows in enumerate(options):
        variables = [
            model.new_bool_var(f"x_{target_pos}_{option_pos}")
            for option_pos in range(len(rows))
        ]
        choose.append(variables)
        model.add_exactly_one(variables)
        for variable, (auxiliary_index, relative_scale, profile) in zip(
            variables, rows, strict=True
        ):
            by_auxiliary[auxiliary_index].append(variable)
            for direction_index, value in enumerate(profile):
                if value:
                    by_direction[direction_index].append((value, variable))
            if (
                auxiliary_index == hinted_aux[target_pos]
                and relative_scale == hinted_scales[target_pos]
            ):
                unchanged.append(variable)
                model.add_hint(variable, 1)
    if len(unchanged) != 16:
        raise ArithmeticError("one hinted option was not found for every target")
    for variables in by_auxiliary:
        model.add_at_most_one(variables)
    for direction_index, target_value in enumerate(goal):
        model.add(
            sum(
                value * variable
                for value, variable in by_direction[direction_index]
            )
            == target_value
        )
    model.maximize(sum(unchanged))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    status = solver.solve(model)
    result: dict[str, object] = {
        "status": solver.status_name(status),
        "wall_seconds": solver.wall_time,
        "branches": solver.num_branches,
        "conflicts": solver.num_conflicts,
        "seed": seed,
        "goal_raw": list(goal),
        "near_score": near.get("score"),
        "objective_bound_unchanged_halves": solver.best_objective_bound,
        "common_graph_constructed": False,
        "residual_ii_closed": False,
    }
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return result

    selected = []
    raw = [0] * 32
    auxiliary_indices = []
    scales = []
    for target_pos, rows in enumerate(options):
        option_pos = next(
            index
            for index, variable in enumerate(choose[target_pos])
            if solver.value(variable)
        )
        auxiliary_index, relative_scale, profile = rows[option_pos]
        auxiliary_indices.append(auxiliary_index)
        scales.append(relative_scale)
        for direction_index, value in enumerate(profile):
            raw[direction_index] += value
        selected.append(
            {
                "target_index": hard[target_pos],
                "target_functional": list(directions[hard[target_pos]]),
                "auxiliary_direction_index": auxiliary_index,
                "relative_scale": relative_scale,
                "auxiliary_functional": list(
                    _scale(directions[auxiliary_index], relative_scale)
                ),
                "parallel_profile": list(profile),
            }
        )
    if raw != list(goal) or len(set(auxiliary_indices)) != 16:
        raise ArithmeticError("the repaired exact profile failed replay")
    unchanged_count = sum(
        (auxiliary_indices[index], scales[index])
        == (hinted_aux[index], hinted_scales[index])
        for index in range(16)
    )
    result["witness"] = {
        "unchanged_half_count": unchanged_count,
        "changed_half_count": 16 - unchanged_count,
        "auxiliary_direction_indices": auxiliary_indices,
        "relative_scales": scales,
        "raw_parallel_profile": raw,
        "selected_halves": selected,
        "exact_parallel_profile": True,
        "auxiliary_directions_distinct": True,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(
        _last_json(args.input),
        seconds=args.seconds,
        workers=args.workers,
        seed=args.seed,
    )
    if args.output:
        write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
