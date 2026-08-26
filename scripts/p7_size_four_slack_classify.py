#!/usr/bin/env python3
"""Classify saturated p=7 slacks for infinity plus three points.

For ``c_H=+1`` every direction has scaled mean eight.  If the three finite
boundary points occupy three odd fibres, the slack parity has 16 mandatory
ones on ``J(7,4)`` and the total slack mass is 20.  Thus a saturated slack
is the parity vector plus either four at one point or two at two points.
This script checks all 630 possibilities against the exact degree-at-most-two
Johnson evaluation space.  The unique survivor adds four at the complement
of the three odd fibres, namely ``A(X)=(|X cap B|-2)^2``.

For the no-infinity branch it also catalogs the two saturated four-odd-fibre
cases.  Phase zero again has one slack, ``A=(t-2)^2``.  Phase one has exactly
36 degree-two integer slacks; their canonical integer coefficient vectors
are returned for the exact edge model.
"""
from __future__ import annotations

import itertools
import json
import math
from functools import lru_cache, reduce


@lru_cache(maxsize=1)
def johnson_space():
    from sympy import Matrix

    p = 7
    points = tuple(itertools.combinations(range(p), 4))
    monomials = (
        ((),)
        + tuple((i,) for i in range(p))
        + tuple(itertools.combinations(range(p), 2))
    )
    evaluation = Matrix(
        [
            [int(set(monomial) <= set(point)) for monomial in monomials]
            for point in points
        ]
    )
    return points, monomials, evaluation, evaluation.T.nullspace()


def classify_three_odd_fibres() -> dict:
    from sympy import Matrix

    p = 7
    points, monomials, evaluation, left_kernel = johnson_space()
    rank = evaluation.rank()
    B = {0, 1, 2}
    parity = Matrix([sum(i in point for i in B) & 1 for point in points])
    mandatory_mass = sum(int(value) for value in parity)
    candidates = []
    for first in range(len(points)):
        for second in range(first, len(points)):
            values = parity.copy()
            values[first] += 2
            values[second] += 2
            if all((vector.dot(values) == 0) for vector in left_kernel):
                candidates.append(
                    {
                        "first_index": first,
                        "second_index": second,
                        "first_point": points[first],
                        "second_point": points[second],
                        "values": tuple(int(value) for value in values),
                    }
                )
    expected_values = tuple(
        (sum(i in point for i in B) - 2) ** 2 for point in points
    )
    unique_expected = bool(
        len(candidates) == 1
        and candidates[0]["first_point"] == (3, 4, 5, 6)
        and candidates[0]["second_point"] == (3, 4, 5, 6)
        and candidates[0]["values"] == expected_values
    )
    return {
        "experiment": "p7_size_four_slack_classify",
        "status": "complete_exact_sparse_correction_enumeration",
        "johnson_points": len(points),
        "degree_at_most_two_feature_columns": len(monomials),
        "degree_at_most_two_rank": rank,
        "left_kernel_dimension": len(left_kernel),
        "mandatory_parity_mass": mandatory_mass,
        "total_saturated_mass": 20,
        "sparse_correction_candidates": len(points) * (len(points) + 1) // 2,
        "survivor_count": len(candidates),
        "survivors": candidates,
        "unique_formula": "A(X)=(|X cap B|-2)^2",
        "proved": unique_expected,
    }


def classify_four_odd_fibres_phase_zero() -> dict:
    """The saturated phase-zero b=4 slack is uniquely (t-2)^2."""
    from sympy import Matrix

    points, monomials, evaluation, left_kernel = johnson_space()
    B = {0, 1, 2, 3}
    parity = Matrix([sum(i in point for i in B) & 1 for point in points])
    candidates = []
    for first in range(len(points)):
        for second in range(first, len(points)):
            values = parity.copy()
            values[first] += 2
            values[second] += 2
            if all(vector.dot(values) == 0 for vector in left_kernel):
                candidates.append(
                    {
                        "first_point": points[first],
                        "second_point": points[second],
                        "values": tuple(int(value) for value in values),
                    }
                )
    expected = tuple((sum(i in point for i in B) - 2) ** 2 for point in points)
    proved = bool(
        evaluation.rank() == 21
        and len(candidates) == 1
        and candidates[0]["first_point"] == (0, 1, 2, 3)
        and candidates[0]["second_point"] == (0, 1, 2, 3)
        and candidates[0]["values"] == expected
    )
    return {
        "phase": 0,
        "b": 4,
        "mandatory_parity_mass": sum(int(value) for value in parity),
        "total_saturated_mass": 20,
        "sparse_correction_candidates": 630,
        "survivor_count": len(candidates),
        "unique_formula": "A(X)=(|X cap B|-2)^2",
        "proved": proved,
    }


def _primitive_left_kernel_rows() -> tuple[tuple[int, ...], ...]:
    from sympy import ilcm

    _points, _monomials, _evaluation, left_kernel = johnson_space()
    rows = []
    for vector in left_kernel:
        denominator = 1
        for value in vector:
            denominator = int(ilcm(denominator, value.q))
        row = [int(value * denominator) for value in vector]
        divisor = reduce(math.gcd, (abs(value) for value in row))
        rows.append(tuple(value // divisor for value in row))
    return tuple(rows)


@lru_cache(maxsize=1)
def classify_four_odd_fibres_phase_one() -> dict:
    """Enumerate all saturated phase-one b=4 degree-two slacks exactly."""
    from ortools.sat.python import cp_model
    from sympy import Matrix

    points, _monomials, _evaluation, _left_kernel = johnson_space()
    B = {0, 1, 2, 3}
    parity = tuple((sum(i in point for i in B) + 1) & 1 for point in points)
    model = cp_model.CpModel()
    lifts = [model.new_int_var(0, 8, f"lift_{i}") for i in range(len(points))]
    model.add(sum(lifts) == 8)
    for row in _primitive_left_kernel_rows():
        model.add(
            sum(
                row[i] * (parity[i] + 2 * lifts[i])
                for i in range(len(points))
            )
            == 0
        )

    lift_solutions: list[tuple[int, ...]] = []

    class Collector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            lift_solutions.append(tuple(self.value(variable) for variable in lifts))

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    status = solver.solve(model, Collector())

    z_monomials = (
        ((),)
        + tuple((i,) for i in range(7))
        + tuple(itertools.combinations(range(7), 2))
    )
    z_evaluation = Matrix(
        [
            [
                1
                if not monomial
                else (2 * int(monomial[0] in point) - 1)
                if len(monomial) == 1
                else (2 * int(monomial[0] in point) - 1)
                * (2 * int(monomial[1] in point) - 1)
                for monomial in z_monomials
            ]
            for point in points
        ]
    )
    _rref, pivots = z_evaluation.rref()
    pivot_matrix = z_evaluation[:, list(pivots)]
    catalog = []
    for lift in sorted(lift_solutions):
        slack_values = tuple(parity[i] + 2 * lift[i] for i in range(len(points)))
        target = Matrix([3 + 2 * value for value in slack_values])
        pivot_coefficients = pivot_matrix.gauss_jordan_solve(target)[0]
        coefficients = [0] * len(z_monomials)
        for pivot, value in zip(pivots, pivot_coefficients):
            if value.q != 1:
                raise AssertionError("canonical target coefficient is not integral")
            coefficients[pivot] = int(value)
        if z_evaluation * Matrix(coefficients) != target:
            raise AssertionError("canonical target coefficients do not reconstruct")
        catalog.append(
            {
                "lift_values": lift,
                "slack_values": slack_values,
                "target_coefficients": tuple(coefficients),
            }
        )
    proved = bool(
        solver.status_name(status) == "OPTIMAL"
        and len(catalog) == 36
        and len(set(row["slack_values"] for row in catalog)) == 36
        and len(pivots) == 21
    )
    return {
        "phase": 1,
        "b": 4,
        "mandatory_parity_mass": sum(parity),
        "total_saturated_mass": 35,
        "total_lift_mass": 8,
        "johnson_degree_two_rank": len(pivots),
        "left_kernel_dimension": len(_primitive_left_kernel_rows()),
        "survivor_count": len(catalog),
        "catalog": catalog,
        "proved": proved,
    }


def main() -> None:
    print(
        json.dumps(
            {
                "three_odd_phase_zero": classify_three_odd_fibres(),
                "four_odd_phase_zero": classify_four_odd_fibres_phase_zero(),
                "four_odd_phase_one": classify_four_odd_fibres_phase_one(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
