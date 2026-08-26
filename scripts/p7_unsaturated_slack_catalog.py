#!/usr/bin/env python3
"""Enumerate exact p=7 degree-two slack catalogs at prescribed means.

For a canonical odd-fibre set ``B={0,...,b-1}``, write every admissible
slack on J(7,4) as ``A=parity+2L``.  Fixed scaled mean ``a=14 E[A]`` fixes
the total lift mass.  The 14 primitive left-kernel equations of the
degree-at-most-two evaluation map are then a finite integer feasibility
problem.  Exhaustive CP-SAT enumeration returns the complete catalog.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_size_four_slack_classify import (  # noqa: E402
    _primitive_left_kernel_rows,
    johnson_space,
)


@functools.lru_cache(maxsize=1)
def _no_linear_interpolation_data() -> dict:
    """Exact interpolation data for the constant-plus-pairs gauge.

    On ``J(7,4)`` one has ``sum z_s=1`` and hence

        l_s z_s = l_s + sum_{t != s} l_s z_s z_t.

    Thus every degree-two target has a representation with no linear
    coefficients.  We first interpolate in a fixed rank-21 gauge and then
    apply this integral identity.  The returned inverse has denominator four,
    so conversion of thousands of catalog rows is exact integer arithmetic.
    """
    from sympy import Matrix, ilcm

    points, _monomials, _evaluation, _left_kernel = johnson_space()
    monomials = (
        ((),)
        + tuple((i,) for i in range(7))
        + tuple(itertools.combinations(range(7), 2))
    )
    evaluation = Matrix(
        [
            [
                1
                if not monomial
                else (2 * int(monomial[0] in point) - 1)
                if len(monomial) == 1
                else (2 * int(monomial[0] in point) - 1)
                * (2 * int(monomial[1] in point) - 1)
                for monomial in monomials
            ]
            for point in points
        ]
    )
    _rref, pivot_columns = evaluation.rref()
    pivot_evaluation = evaluation[:, list(pivot_columns)]
    _row_rref, pivot_rows = pivot_evaluation.T.rref()
    inverse = pivot_evaluation[list(pivot_rows), :].inv()
    denominator = 1
    for value in inverse:
        denominator = int(ilcm(denominator, value.q))
    numerator = tuple(
        tuple(int(inverse[i, j] * denominator) for j in range(21))
        for i in range(21)
    )
    return {
        "points": points,
        "pairs": tuple(itertools.combinations(range(7), 2)),
        "evaluation": evaluation,
        "pivot_columns": tuple(int(value) for value in pivot_columns),
        "pivot_rows": tuple(int(value) for value in pivot_rows),
        "inverse_numerator": numerator,
        "inverse_denominator": denominator,
    }


def no_linear_target_row(slack_values: tuple[int, ...]) -> tuple[int, ...]:
    """Return exact ``(constant, 21 pair coefficients)`` for ``3+2A``."""
    data = _no_linear_interpolation_data()
    if len(slack_values) != 35:
        raise ValueError("need one slack value at every point of J(7,4)")
    target = tuple(3 + 2 * int(value) for value in slack_values)
    sampled = tuple(target[index] for index in data["pivot_rows"])
    denominator = int(data["inverse_denominator"])
    coefficient_numerators = tuple(
        sum(row[j] * sampled[j] for j in range(21))
        for row in data["inverse_numerator"]
    )
    if any(value % denominator for value in coefficient_numerators):
        raise AssertionError("target has no integral canonical coefficients")
    pivot_coefficients = tuple(
        value // denominator for value in coefficient_numerators
    )
    full = [0] * 29
    for index, value in zip(data["pivot_columns"], pivot_coefficients):
        full[index] = value
    constant = full[0] + sum(full[1:8])
    pair_coefficients = tuple(
        full[8 + index] + full[1 + s] + full[1 + t]
        for index, (s, t) in enumerate(data["pairs"])
    )
    row = (constant, *pair_coefficients)
    reconstructed = []
    for point in data["points"]:
        reconstructed.append(
            constant
            + sum(
                pair_coefficients[index]
                * (1 if ((s in point) == (t in point)) else -1)
                for index, (s, t) in enumerate(data["pairs"])
            )
        )
    if tuple(reconstructed) != target:
        raise AssertionError("constant-plus-pairs target failed reconstruction")
    return row


def enumerate_catalog(
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
    max_solutions: int | None = None,
    seconds: float | None = None,
    include_values: bool = False,
) -> dict:
    from ortools.sat.python import cp_model

    if not 0 <= odd_fibres <= 7 or phase not in (0, 1):
        raise ValueError("need b in 0..7 and phase in {0,1}")
    if scaled_mean < 0 or scaled_mean % 2:
        raise ValueError("scaled mean must be a nonnegative even integer")
    points, _monomials, evaluation, _left_kernel = johnson_space()
    B = set(range(odd_fibres))
    parity = tuple(
        (sum(index in point for index in B) + phase) & 1 for point in points
    )
    parity_mass = sum(parity)
    numerator = 5 * scaled_mean - 2 * parity_mass
    if numerator < 0 or numerator % 4:
        return {
            "experiment": "p7_unsaturated_slack_catalog",
            "status": "parity_mean_incompatible",
            "complete": True,
            "odd_fibres": odd_fibres,
            "phase": phase,
            "scaled_mean": scaled_mean,
            "parity_mass": parity_mass,
            "solution_count": 0,
        }
    lift_mass = numerator // 4
    model = cp_model.CpModel()
    lifts = [
        model.new_int_var(0, lift_mass, f"lift_{index}")
        for index in range(len(points))
    ]
    model.add(sum(lifts) == lift_mass)
    kernel_rows = _primitive_left_kernel_rows()
    for row in kernel_rows:
        model.add(
            sum(
                row[index] * (parity[index] + 2 * lifts[index])
                for index in range(len(points))
            )
            == 0
        )

    solutions: list[tuple[int, ...]] = []
    started = time.time()

    class Collector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            solutions.append(tuple(self.value(variable) for variable in lifts))
            if max_solutions is not None and len(solutions) >= max_solutions:
                self.stop_search()

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    if seconds is not None:
        solver.parameters.max_time_in_seconds = float(seconds)
    status = solver.solve(model, Collector())
    status_name = solver.status_name(status)
    complete = status_name == "OPTIMAL"
    unique = sorted(set(solutions))
    out = {
        "experiment": "p7_unsaturated_slack_catalog",
        "status": "complete_exact_enumeration" if complete else "partial_enumeration",
        "complete": complete,
        "solver_status": status_name,
        "odd_fibres": odd_fibres,
        "phase": phase,
        "scaled_mean": scaled_mean,
        "johnson_points": len(points),
        "degree_two_rank": evaluation.rank(),
        "left_kernel_dimension": len(kernel_rows),
        "parity_mass": parity_mass,
        "total_slack_mass": 5 * scaled_mean // 2,
        "lift_mass": lift_mass,
        "solution_count": len(solutions),
        "unique_solution_count": len(unique),
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": time.time() - started,
    }
    if include_values:
        out["lift_values"] = unique
    return out


@functools.lru_cache(maxsize=None)
def exact_slack_catalog_values(
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
) -> tuple[tuple[int, ...], ...]:
    """Return the complete exact slack catalog for a canonical odd set.

    The three 1,764-element catalogs are not enumerated independently.  Let
    ``E`` be a nonnegative, even, degree-two function with scaled mean eight.
    For ``b=2`` the minimum slack is exactly its pointwise parity baseline,
    so every next-mean slack is uniquely ``A_min+E``.  This proves a bijection
    with the ``b=0, phase=0, mean=8`` excess catalog and explains the repeated
    count ``1764=36*7^2``.
    """
    if not 0 <= odd_fibres <= 7 or phase not in (0, 1):
        raise ValueError("need b in 0..7 and phase in {0,1}")
    points, _monomials, _evaluation, _left_kernel = johnson_space()
    B = set(range(odd_fibres))
    parity = tuple(
        (sum(index in point for index in B) + phase) & 1 for point in points
    )

    unique_floor: tuple[int, ...] | None = None
    if odd_fibres == 0 and phase == 0 and scaled_mean == 0:
        unique_floor = (0,) * len(points)
    elif odd_fibres == 0 and phase == 1 and scaled_mean == 14:
        unique_floor = (1,) * len(points)
    elif odd_fibres in (1, 2, 5, 6) and scaled_mean == (
        8 if phase == 0 else 6
    ):
        unique_floor = parity
    elif odd_fibres in (3, 4) and phase == 0 and scaled_mean == 8:
        unique_floor = tuple(
            (sum(index in point for index in B) - 2) ** 2 for point in points
        )
    if unique_floor is not None:
        return (unique_floor,)

    # Universal even excess catalog.  Its lift vectors are the slack divided
    # by two because the parity baseline is identically zero.
    if (odd_fibres, phase, scaled_mean) in {
        (0, 0, 8),
        *((b, 0, 16) for b in (1, 2, 5, 6)),
        *((b, 1, 14) for b in (1, 2, 5, 6)),
    }:
        excess = enumerate_catalog(0, 0, 8, include_values=True)
        if not excess["complete"] or excess["unique_solution_count"] != 1764:
            raise AssertionError("universal p=7 excess catalog is incomplete")
        excess_values = tuple(
            tuple(2 * int(value) for value in lift)
            for lift in excess["lift_values"]
        )
        if odd_fibres == 0:
            rows = excess_values
        else:
            rows = tuple(
                tuple(parity[index] + excess_row[index] for index in range(35))
                for excess_row in excess_values
            )
    else:
        direct = enumerate_catalog(
            odd_fibres,
            phase,
            scaled_mean,
            include_values=True,
        )
        if not direct["complete"]:
            raise AssertionError("requested p=7 slack catalog is incomplete")
        rows = tuple(
            tuple(
                parity[index] + 2 * int(lift[index]) for index in range(35)
            )
            for lift in direct["lift_values"]
        )

    unique = tuple(sorted(set(rows)))
    if any(2 * sum(row) != 5 * scaled_mean for row in unique):
        raise AssertionError("catalog row has the wrong scaled mean")
    return unique


@functools.lru_cache(maxsize=None)
def exact_target_catalog_rows(
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
) -> tuple[tuple[int, ...], ...]:
    """Complete no-linear target catalog in canonical fibre coordinates."""
    rows = tuple(
        no_linear_target_row(values)
        for values in exact_slack_catalog_values(odd_fibres, phase, scaled_mean)
    )
    unique = tuple(sorted(set(rows)))
    if len(unique) != len(rows):
        raise AssertionError("target interpolation collapsed distinct slacks")
    return unique


def mapped_target_catalog_rows(
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
    B: set[int],
) -> tuple[tuple[int, ...], ...]:
    """Map a canonical target catalog to an arbitrary odd-fibre set ``B``."""
    if len(B) != odd_fibres:
        raise ValueError("B has the wrong cardinality")
    canonical = exact_target_catalog_rows(odd_fibres, phase, scaled_mean)
    canonical_pairs = tuple(itertools.combinations(range(7), 2))
    pair_index = {pair: index for index, pair in enumerate(canonical_pairs)}
    actual_B = sorted(B)
    actual_complement = sorted(set(range(7)) - B)
    permutation = dict(zip(range(odd_fibres), actual_B)) | dict(
        zip(range(odd_fibres, 7), actual_complement)
    )
    mapped = []
    for row in canonical:
        pairs = [0] * 21
        for index, (s, t) in enumerate(canonical_pairs):
            endpoints = tuple(sorted((permutation[s], permutation[t])))
            pairs[pair_index[endpoints]] = int(row[1 + index])
        mapped.append((int(row[0]), *pairs))
    unique = tuple(sorted(set(mapped)))
    if len(unique) != len(canonical):
        raise AssertionError("fibre relabeling collapsed a target")
    return unique


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--b", type=int, choices=(0, 2, 4), required=True)
    parser.add_argument("--phase", type=int, choices=(0, 1), required=True)
    parser.add_argument("--mean", type=int, required=True)
    parser.add_argument("--max-solutions", type=int)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--include-values", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = enumerate_catalog(
        args.b,
        args.phase,
        args.mean,
        args.max_solutions,
        args.seconds,
        args.include_values,
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
