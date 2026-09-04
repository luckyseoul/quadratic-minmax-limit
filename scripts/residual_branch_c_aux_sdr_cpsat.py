#!/usr/bin/env python3
"""Exact CP-SAT probe of the branch-C paired auxiliary-SDR gate.

This script starts at the open condition in
``e1_gmin_m4_adaptive_mobius_pairing``.  Fix an opposite direction ``F``
and normalize the other projective directions to the affine chart
``L_z=G+zF``.  For prescribed nonzero values

    alpha_i=(L_i/j_i)(x_0)

on the ``m=(p+1)/2`` hard targets, a choice of fixed-edge scale ``c``, a
pair of targets, and two singleton signs uniquely force two auxiliary
directions.  A feasible member of this deliberately paired complementary-SDR
ansatz must choose ``m/2`` such options so that

* every hard target occurs exactly once;
* every auxiliary direction occurs at most once; and
* the auxiliaries contain exactly ``m-2`` hard and two opposite directions.

Those are encoded together, not checked in separate assignments.  A SAT
witness is replayed through ``prescribed_auxiliary_assignment_criterion``.
The calculation is an exploratory target-sensitive probe.  SAT for sampled
centre profiles is not a residual-(ii) proof, and UNSAT excludes only the
stated profile, fixed direction, and paired complementary-SDR ansatz.  It
does not exclude unpaired half assignments or higher-multiplicity physical
overlap patterns.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_adaptive_mobius_pairing import (  # noqa: E402
    forced_affine_auxiliary_pair,
    prescribed_auxiliary_assignment_criterion,
    quadratic_character,
)
from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Point,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
)
from io_atomic import write_json_atomic  # noqa: E402


@dataclass(frozen=True)
class Chart:
    p: int
    fixed_direction_index: int
    fixed_direction: Point
    x0: Point
    base_direction: Point
    direction_by_coordinate: tuple[Point, ...]
    type_by_coordinate: tuple[int, ...]
    hard_coordinates: tuple[int, ...]
    opposite_coordinates: tuple[int, ...]


@dataclass(frozen=True)
class Option:
    first_target: int
    second_target: int
    first_sign: int
    second_sign: int
    first_auxiliary: int
    second_auxiliary: int

    @property
    def auxiliary_pair(self) -> tuple[int, int]:
        return tuple(sorted((self.first_auxiliary, self.second_auxiliary)))


def _inverse(p: int, value: int) -> int:
    value %= p
    if value == 0:
        raise ZeroDivisionError("cannot invert zero")
    return pow(value, -1, p)


def _scale(p: int, scalar: int, vector: Point) -> Point:
    return scalar * vector[0] % p, scalar * vector[1] % p


def _determinant(p: int, first: Point, second: Point) -> int:
    return (first[0] * second[1] - first[1] * second[0]) % p


def _affine_coefficient(p: int, base: Point, fixed: Point, value: Point) -> int:
    """Return z in value=base+z*fixed, checking the normalization."""
    determinant = _determinant(p, base, fixed)
    if determinant == 0:
        raise ArithmeticError("chart basis is dependent")
    base_coefficient = _determinant(p, value, fixed) * _inverse(p, determinant) % p
    fixed_coefficient = _determinant(p, base, value) * _inverse(p, determinant) % p
    if base_coefficient != 1:
        raise ArithmeticError("a normalized direction left the affine chart")
    return fixed_coefficient


def build_chart(p: int, fixed_direction_index: int | None = None) -> Chart:
    """Build the exact F-affine chart with ``F`` of opposite Paley type."""
    directions = projective_functionals(p)
    types = tuple(paley_direction_sign(p, direction) for direction in directions)
    hard_indices = tuple(index for index, value in enumerate(types) if value == 1)
    opposite_indices = tuple(index for index, value in enumerate(types) if value == -1)
    m = (p + 1) // 2
    if len(hard_indices) != m or len(opposite_indices) != m:
        raise ArithmeticError("the projective Paley types lost balance")
    if fixed_direction_index is None:
        fixed_direction_index = opposite_indices[0]
    if fixed_direction_index not in opposite_indices:
        raise ValueError("the fixed direction must have opposite Paley type")

    fixed = directions[fixed_direction_index]
    # This nonzero vector spans ker(F).
    x0 = (fixed[1] % p, -fixed[0] % p)
    if x0 == (0, 0):
        raise ArithmeticError("the fixed direction produced a zero kernel vector")
    base_index = hard_indices[0]
    base_raw = directions[base_index]
    base_value = (base_raw[0] * x0[0] + base_raw[1] * x0[1]) % p
    base = _scale(p, _inverse(p, base_value), base_raw)

    by_coordinate: list[Point | None] = [None] * p
    type_by_coordinate = [0] * p
    for index, direction in enumerate(directions):
        if index == fixed_direction_index:
            continue
        evaluation = (direction[0] * x0[0] + direction[1] * x0[1]) % p
        if evaluation == 0:
            raise ArithmeticError("a second direction annihilated x0")
        normalized = _scale(p, _inverse(p, evaluation), direction)
        coordinate = _affine_coefficient(p, base, fixed, normalized)
        if by_coordinate[coordinate] is not None:
            raise ArithmeticError("the affine chart repeated a coordinate")
        by_coordinate[coordinate] = normalized
        type_by_coordinate[coordinate] = types[index]
    if any(direction is None for direction in by_coordinate):
        raise ArithmeticError("the affine chart missed a coordinate")
    hard_coordinates = tuple(
        coordinate for coordinate, value in enumerate(type_by_coordinate) if value == 1
    )
    opposite_coordinates = tuple(
        coordinate for coordinate, value in enumerate(type_by_coordinate) if value == -1
    )
    if len(hard_coordinates) != m or len(opposite_coordinates) != m - 1:
        raise ArithmeticError("the punctured projective type counts changed")
    return Chart(
        p=p,
        fixed_direction_index=fixed_direction_index,
        fixed_direction=fixed,
        x0=x0,
        base_direction=base,
        direction_by_coordinate=tuple(by_coordinate),  # type: ignore[arg-type]
        type_by_coordinate=tuple(type_by_coordinate),
        hard_coordinates=hard_coordinates,
        opposite_coordinates=opposite_coordinates,
    )


def generate_options(chart: Chart, alphas: Sequence[int], c: int) -> tuple[Option, ...]:
    p = chart.p
    targets = chart.hard_coordinates
    if len(alphas) != len(targets) or any(alpha % p == 0 for alpha in alphas):
        raise ValueError("need one nonzero alpha for every hard target")
    unique: dict[tuple[int, int, int, int], Option] = {}
    for first in range(len(targets)):
        for second in range(first + 1, len(targets)):
            for first_sign in (-1, 1):
                for second_sign in (-1, 1):
                    try:
                        row = forced_affine_auxiliary_pair(
                            p,
                            targets[first],
                            targets[second],
                            int(alphas[first]),
                            int(alphas[second]),
                            c,
                            first_sign,
                            second_sign,
                        )
                    except ValueError:
                        continue
                    first_auxiliary = int(row["first_auxiliary_coordinate_U"])
                    second_auxiliary = int(row["second_auxiliary_coordinate_V"])
                    if first_auxiliary == second_auxiliary:
                        raise ArithmeticError("a valid option repeated its auxiliary")
                    option = Option(
                        first,
                        second,
                        first_sign,
                        second_sign,
                        first_auxiliary,
                        second_auxiliary,
                    )
                    # If two sign choices induce the same target/auxiliary block,
                    # either representative suffices for this existential gate.
                    key = (first, second, *option.auxiliary_pair)
                    unique.setdefault(key, option)
    return tuple(unique[key] for key in sorted(unique))


def solve_scale(
    chart: Chart,
    alphas: Sequence[int],
    c: int,
    seconds: float,
    workers: int,
    random_seed: int,
) -> dict[str, object]:
    from ortools.sat.python import cp_model

    options = generate_options(chart, alphas, c)
    model = cp_model.CpModel()
    selected = tuple(model.new_bool_var(f"x_{index}") for index in range(len(options)))
    target_rows: list[list[int]] = [[] for _ in chart.hard_coordinates]
    auxiliary_rows: list[list[int]] = [[] for _ in range(chart.p)]
    for index, option in enumerate(options):
        target_rows[option.first_target].append(index)
        target_rows[option.second_target].append(index)
        auxiliary_rows[option.first_auxiliary].append(index)
        auxiliary_rows[option.second_auxiliary].append(index)
    if any(not row for row in target_rows):
        return {
            "c": c,
            "status": "INFEASIBLE_BY_EMPTY_TARGET_ROW",
            "option_count": len(options),
            "wall_seconds": 0.0,
        }
    for row in target_rows:
        model.add_exactly_one(selected[index] for index in row)
    for row in auxiliary_rows:
        if row:
            model.add_at_most_one(selected[index] for index in row)
    hard_terms = []
    for index, option in enumerate(options):
        hard_count = (
            int(chart.type_by_coordinate[option.first_auxiliary] == 1)
            + int(chart.type_by_coordinate[option.second_auxiliary] == 1)
        )
        if hard_count:
            hard_terms.append(hard_count * selected[index])
    model.add(sum(hard_terms) == len(chart.hard_coordinates) - 2)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = random_seed
    solver.parameters.log_search_progress = False
    started = time.monotonic()
    status = solver.solve(model)
    wall = time.monotonic() - started
    name = solver.status_name(status)
    result: dict[str, object] = {
        "c": c,
        "status": name,
        "option_count": len(options),
        "wall_seconds": wall,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
    }
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return result

    chosen = tuple(options[index] for index in range(len(options)) if solver.value(selected[index]))
    m = len(chart.hard_coordinates)
    if len(chosen) != m // 2:
        raise ArithmeticError("the selected target matching has wrong size")
    auxiliaries = [0] * m
    partners = [-1] * m
    signs = [0] * m
    for option in chosen:
        i, k = option.first_target, option.second_target
        auxiliaries[i] = option.first_auxiliary
        auxiliaries[k] = option.second_auxiliary
        partners[i] = k
        partners[k] = i
        signs[i] = option.first_sign
        signs[k] = option.second_sign
    replay = prescribed_auxiliary_assignment_criterion(
        chart.p,
        chart.hard_coordinates,
        tuple(int(alpha) % chart.p for alpha in alphas),
        auxiliaries,
        partners,
        tuple(range(m)),
        c,
        signs,
    )
    hard_auxiliaries = sum(chart.type_by_coordinate[value] == 1 for value in auxiliaries)
    if not replay["pair_coherent_distinct_auxiliary_assignment"]:
        raise ArithmeticError("the exact paired-SDR witness replay failed")
    if hard_auxiliaries != m - 2:
        raise ArithmeticError("the selected witness missed the Paley type quota")
    result["witness"] = {
        "target_coordinates": list(chart.hard_coordinates),
        "partners": partners,
        "signs": signs,
        "auxiliary_coordinates": auxiliaries,
        "hard_auxiliary_count": hard_auxiliaries,
        "opposite_auxiliary_count": m - hard_auxiliaries,
        "chosen_options": [
            {
                "targets": [option.first_target, option.second_target],
                "signs": [option.first_sign, option.second_sign],
                "auxiliaries": [option.first_auxiliary, option.second_auxiliary],
            }
            for option in chosen
        ],
        "exact_assignment_replay": True,
    }
    return result


def _profile(mode: str, p: int, m: int, seed: int) -> tuple[int, ...]:
    if mode == "constant":
        return (1,) * m
    if mode == "ramp":
        return tuple(index % (p - 1) + 1 for index in range(m))
    if mode == "alternating":
        return tuple(1 if index % 2 == 0 else p - 1 for index in range(m))
    if mode != "random":
        raise ValueError(f"unknown alpha mode {mode}")
    generator = random.Random(seed)
    return tuple(generator.randrange(1, p) for _ in range(m))


def solve_profile(
    chart: Chart,
    alphas: Sequence[int],
    seconds_per_scale: float,
    workers: int,
    seed: int,
) -> dict[str, object]:
    scale_rows = []
    for offset in range(chart.p - 1):
        c = 1 + (seed + offset) % (chart.p - 1)
        row = solve_scale(chart, alphas, c, seconds_per_scale, workers, seed + offset)
        scale_rows.append(row)
        if row["status"] in ("OPTIMAL", "FEASIBLE"):
            break
    statuses = [str(row["status"]) for row in scale_rows]
    alpha_bytes = b"".join(
        int(value).to_bytes(4, byteorder="little", signed=False)
        for value in alphas
    )
    return {
        "alphas": list(alphas),
        "alpha_sha256": hashlib.sha256(alpha_bytes).hexdigest(),
        "scales_tested": len(scale_rows),
        "status": "SAT" if any(value in ("OPTIMAL", "FEASIBLE") for value in statuses) else (
            "UNSAT" if all(value.startswith("INFEASIBLE") for value in statuses) else "UNKNOWN"
        ),
        "scale_results": scale_rows,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    chart = build_chart(args.p, args.fixed_direction_index)
    modes = tuple(args.alpha_modes.split(","))
    jobs: list[tuple[str, int]] = []
    for mode in modes:
        count = args.trials if mode == "random" else 1
        for trial in range(count):
            global_trial = args.shard_index + args.shard_count * trial
            jobs.append((mode, args.seed + global_trial))
    records = []
    started = time.monotonic()
    for mode, seed in jobs:
        alphas = _profile(mode, args.p, len(chart.hard_coordinates), seed)
        record = solve_profile(chart, alphas, args.seconds_per_scale, args.workers, seed)
        record["alpha_mode"] = mode
        record["seed"] = seed
        records.append(record)
    return {
        "schema": "residual_branch_c_aux_sdr_cpsat_v1",
        "scope": (
            "target-sensitive paired complementary auxiliary SDR; sampled center "
            "profiles and paired ansatz only"
        ),
        "paired_ansatz_exhaustive_for_residual_ii": False,
        "p": args.p,
        "fixed_direction_index": chart.fixed_direction_index,
        "fixed_direction": list(chart.fixed_direction),
        "x0": list(chart.x0),
        "base_direction": list(chart.base_direction),
        "hard_coordinates": list(chart.hard_coordinates),
        "opposite_coordinates": list(chart.opposite_coordinates),
        "hard_target_count": len(chart.hard_coordinates),
        "required_auxiliary_types": {
            "hard": len(chart.hard_coordinates) - 2,
            "opposite": 2,
        },
        "shard": [args.shard_index, args.shard_count],
        "profile_count": len(records),
        "sat_count": sum(row["status"] == "SAT" for row in records),
        "unsat_count": sum(row["status"] == "UNSAT" for row in records),
        "unknown_count": sum(row["status"] == "UNKNOWN" for row in records),
        "profiles": records,
        "elapsed_seconds": time.monotonic() - started,
        "interpretation": (
            "SAT constructs only the paired complementary auxiliary-SDR layer; "
            "UNSAT excludes only that paired ansatz for the recorded alpha profile "
            "and fixed direction. Unpaired assignments and higher-multiplicity "
            "overlaps remain outside scope; neither outcome closes residual (ii)."
        ),
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=31)
    parser.add_argument("--fixed-direction-index", type=int)
    parser.add_argument("--alpha-modes", default="constant,ramp,alternating,random")
    parser.add_argument("--trials", type=int, default=8)
    parser.add_argument("--seed", type=int, default=15766)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seconds-per-scale", type=float, default=2.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    result = run(args)
    if args.output is None:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        write_json_atomic(args.output, result)
        print(
            f"wrote {args.output}: profiles={result['profile_count']} "
            f"SAT={result['sat_count']} UNSAT={result['unsat_count']} "
            f"UNKNOWN={result['unknown_count']}"
        )


if __name__ == "__main__":
    main()
