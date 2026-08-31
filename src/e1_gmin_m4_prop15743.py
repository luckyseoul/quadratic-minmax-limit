#!/usr/bin/env python3
r"""Prop. 15.743 -- full translated-cut energy closes ``p=17,t=3``.

Proposition 15.739 leaves one generic fourth-shell branch at ``p=17``.
There are nine hard and nine opposite affine directions.  The hard
quotients satisfy ``sum k_L=12`` and ``k_L>=1``; the opposite directions
all have parallel count three.

The first new point is genuinely global.  For a hard row with parallel
count ``P`` and quotient ``k``, local coefficient comparison gives

    sum_a q_L(a) = 17(P-3)-18k.

On the other hand, let ``h`` be the hard direction sign and ``T`` the signed
global edge total.  Before choosing a representative of an exact star, glue
the two unspecialized sums at ``k=1``:

    17(P-3)-18 = hT-P,

so ``hT=18P-69``.  The same ``hT`` occurs in every direction, hence all
exact stars have one common ``P``.  There are at least six, so ``6P<=75``;
the isolated-chart literal congruence ``P=5 mod 8`` now forces ``P=5``.
Only after that normalization is fixed do we obtain ``hT=21`` and the
exact-star row ``q=(2)^8``.  This agrees with the independent split into 48
hard and 27 opposite displacements.  Every hard-sign-normalized row is the
nonzero part of the same difference-Radon transform and therefore

    sum_a q_L(a) = hT-P = 21-P.

Equating the two forces ``P=4+k``.  Local cells with another ``P`` do not
lift to one common graph.

Six exact hard stars force the global moments ``M_2=M_4=0``.  For each
nonexact row, aggregate the 17 coefficients in each of the eight cyclic
distance classes.  Translation-summing all middle-slice cut inequalities
gives 698 distinct linear inequalities.  Exact one-worker CP-SAT models,
using only row sum, l1, the two modular moments, and those 698 upper cuts,
give

* hard excess one: infeasible;
* hard excess two: energy at most 70;
* hard excess three: energy at most 119;
* opposite: energy at most 72.

The first two hard quotient partitions contain an excess-one row and are
therefore empty.  In the remaining partition ``3``, the nonexact energy is
at most ``119+9*72=767``.  The common difference-Radon Parseval identity is
instead ``1211+34C>=1211``.  This contradiction closes ``p=17,k=74``.

No prior row-energy cap, entry alphabet, lower cut bound, row parity, full
136-entry coefficient matrix, or complete-domain Boolean model is used.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15734 import BRANCH_P1_LAST, baseline_coefficient_rules
from e1_gmin_m4_prop15739 import (
    generic_higher_even_moment_reduction,
    p17_conditioned_cut_reduction,
)


ROOT = Path(__file__).resolve().parents[1]
P = 17
M = 9
Q = 8
H_EDGE_COUNT = 75
HARD_DIRECTION_COUNT = 9
OPPOSITE_DIRECTION_COUNT = 9
HARD_EDGE_COUNT = 48
OPPOSITE_EDGE_COUNT = 27
HARD_SIGN_TIMES_T = 21
DISTANCES = tuple(range(1, Q + 1))
HARD_EXCESS_PARTITIONS = ((1, 1, 1), (2, 1), (3,))

EXPECTED_CUT_VECTOR_COUNT = 698
EXPECTED_CUT_CATALOG_SHA256 = (
    "a8ac7349cb601db5163ef1526949587c766914d774fe26858fe93eac1d940708"
)
EXPECTED_ROW_MAXIMA = {
    "hard_e2": 70,
    "hard_e3": 119,
    "opposite": 72,
}
EXPECTED_MAXIMIZERS = {
    "hard_e2": (1, -2, 5, 3, -1, 2, 5, 1),
    "hard_e3": (6, 4, -1, -3, -2, -1, 4, 6),
    "opposite": (-3, -3, -3, -3, -3, -3, -3, -3),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _digest(rows: Iterable[tuple[int, ...]]) -> str:
    payload = ";".join(",".join(map(str, row)) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def translated_cut_vector(subset: Iterable[int]) -> tuple[int, ...]:
    """Return ``(|X triangle (X+a)|)_(a=1..8)`` in ``F_17``."""
    values = tuple(sorted(int(value) for value in subset))
    if len(values) != M or len(set(values)) != M or not all(
        0 <= value < P for value in values
    ):
        raise ValueError("need nine distinct elements of F_17")
    chosen = set(values)
    return tuple(
        sum(
            (value in chosen) != ((value + distance) % P in chosen)
            for value in range(P)
        )
        for distance in DISTANCES
    )


@lru_cache(maxsize=1)
def translated_cut_vectors() -> tuple[tuple[int, ...], ...]:
    """Generate the exact catalog of distinct translated-cut vectors."""
    return tuple(
        sorted(
            {
                translated_cut_vector(subset)
                for subset in combinations(range(P), M)
            }
        )
    )


@lru_cache(maxsize=1)
def p17_cut_catalog_certificate() -> dict[str, object]:
    """Audit the 698 full translated-cut coefficient vectors."""
    vectors = translated_cut_vectors()
    digest = _digest(vectors)
    proved = bool(
        len(vectors) == EXPECTED_CUT_VECTOR_COUNT
        and digest == EXPECTED_CUT_CATALOG_SHA256
        and all(
            len(vector) == Q
            and sum(vector) == 72
            and all(value % 2 == 0 and 0 <= value <= 16 for value in vector)
            for vector in vectors
        )
    )
    _require(proved, "the p17 translated-cut catalog changed")
    return {
        "p": P,
        "slice_size": M,
        "middle_slice_point_count": 24_310,
        "distinct_translated_cut_vectors": len(vectors),
        "every_vector_sum": 72,
        "every_entry_even_between_zero_and_sixteen": True,
        "catalog_sha256": digest,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def two_source_hard_normalization_certificate() -> dict[str, object]:
    """Glue local coefficient normalization to one global signed graph."""
    moment = generic_higher_even_moment_reduction(P)
    conditioned = p17_conditioned_cut_reduction()
    literal_rule = baseline_coefficient_rules(P)[BRANCH_P1_LAST]
    exact_star_count = int(moment["exact_hard_star_directions_at_least"])

    # Do not insert q=(2)^8 here: before P is fixed, a slice-kernel gauge
    # changes both the zero bin and every nonzero distance bin.  At k=1 the
    # unspecialized local and common sums instead give
    #
    #   17(P_L-3)-18 = hT-P_L,  hence  hT=18P_L-69.
    #
    # The global value hT is common to all directions, and the affine map on
    # the right is injective over the integers.  Thus all exact stars share
    # one P_L before any normalized distance row is used.
    exact_star_quotient = 1
    exact_star_hT_slope = P + 1
    exact_star_hT_intercept = -(3 * P + (P + 1) * exact_star_quotient)
    exact_stars_have_common_parallel_count = exact_star_hT_slope != 0
    exact_star_parallel_upper = H_EDGE_COUNT // exact_star_count
    isolated_chart_I = 0
    exact_star_parallel_candidates = [
        parallel_count
        for parallel_count in range(exact_star_parallel_upper + 1)
        if (
            isolated_chart_I
            + parallel_count
            - int(literal_rule["offset"])
        )
        % Q
        == 0
    ]
    exact_star_parallel_count = exact_star_parallel_candidates[0]
    hard_sign_times_T_from_stars = (
        exact_star_hT_slope * exact_star_parallel_count
        + exact_star_hT_intercept
    )

    # P=5 has now killed the slice-kernel gauge.  Only at this point is the
    # literal-star distance row canonically (2)^8.
    exact_star = {
        tuple(
            sum(
                int(
                    min((center - label) % P, (label - center) % P)
                    == distance
                )
                for label in range(P)
                if label != center
            )
            for distance in DISTANCES
        )
        for center in range(P)
    }
    exact_star_off_bin_sum = sum(next(iter(exact_star)))

    opposite_parallel_count = int(conditioned["opposite_parallel_count_Q"])
    opposite_edge_count = OPPOSITE_DIRECTION_COUNT * opposite_parallel_count
    hard_edge_count = H_EDGE_COUNT - opposite_edge_count
    hard_sign_times_T_from_edge_split = hard_edge_count - opposite_edge_count
    hard_sign_times_T = hard_sign_times_T_from_stars
    rows: list[dict[str, object]] = []
    for k in range(1, 5):
        matches = []
        for parallel_count in range(H_EDGE_COUNT + 1):
            local_sum = P * (parallel_count - 3) - (P + 1) * k
            common_sum = hard_sign_times_T - parallel_count
            if local_sum == common_sum:
                matches.append(parallel_count)
        forced = 4 + k
        rows.append(
            {
                "hard_quotient_k": k,
                "local_coefficient_sum_formula": "17*(P_L-3)-18*k_L",
                "common_difference_sum_formula": "hT-P_L=21-P_L",
                "matching_parallel_counts_in_0_through_75": matches,
                "forced_parallel_count": forced,
                "off_bin_sum": hard_sign_times_T - forced,
                "cellwise_cut_upper_bound": M * (forced - 3 - k),
                "translated_cut_upper_bound": P * M * (forced - 3 - k),
                "l1_bound": H_EDGE_COUNT - forced,
            }
        )

    proved = bool(
        moment["proved"]
        and conditioned["proved"]
        and moment["m"] == HARD_DIRECTION_COUNT
        and moment["hard_quotient_identity"] == "sum k_d=12"
        and exact_star_count == 6
        and moment["global_even_moments_forced_identically_zero"] == [2, 4]
        and literal_rule["offset"] == 5
        and literal_rule["congruence"] == "8 divides I+P-5"
        and exact_star_quotient == 1
        and exact_star_hT_slope == 18
        and exact_star_hT_intercept == -69
        and exact_stars_have_common_parallel_count
        and exact_star_parallel_upper == 12
        and isolated_chart_I == 0
        and exact_star_parallel_candidates == [5]
        and hard_sign_times_T_from_stars == HARD_SIGN_TIMES_T == 21
        and exact_star == {(2,) * Q}
        and exact_star_off_bin_sum == 16
        and exact_star_off_bin_sum
        == hard_sign_times_T_from_stars - exact_star_parallel_count
        and opposite_parallel_count == 3
        and conditioned["sum_W"] == -24
        and opposite_edge_count == OPPOSITE_EDGE_COUNT == 27
        and hard_edge_count == HARD_EDGE_COUNT == 48
        and hard_edge_count + opposite_edge_count == H_EDGE_COUNT
        and hard_sign_times_T_from_edge_split == HARD_SIGN_TIMES_T
        and all(row["matching_parallel_counts_in_0_through_75"] == [4 + k]
                for k, row in zip(range(1, 5), rows))
        and all(row["cellwise_cut_upper_bound"] == M for row in rows)
        and all(row["translated_cut_upper_bound"] == P * M for row in rows)
    )
    _require(proved, "the p17 two-source hard normalization changed")
    return {
        "p": P,
        "exact_star_count_lower_bound": exact_star_count,
        "normalization_order": [
            "glue unspecialized exact-row sums",
            "deduce one common exact-star P_L",
            "apply six-star edge bound and isolated-chart congruence",
            "deduce P_L=5 and hT=21",
            "only then identify q=(2)^8",
        ],
        "exact_star_unspecialized_glue": (
            "17*(P_L-3)-18=hT-P_L"
        ),
        "exact_star_hT_affine_identity": "hT=18*P_L-69",
        "exact_star_common_parallel_count_reason": (
            "common hT and hT=18*P_L-69 force one common exact-star P_L"
        ),
        "exact_star_edge_bound": "6*P_L<=75",
        "exact_literal_parallel_congruence": literal_rule["congruence"],
        "exact_literal_isolated_chart_I": isolated_chart_I,
        "exact_star_parallel_count_upper_bound": exact_star_parallel_upper,
        "exact_star_parallel_candidates": exact_star_parallel_candidates,
        "forced_exact_star_parallel_count": exact_star_parallel_count,
        "hard_sign_times_global_T_from_exact_stars": (
            "hT=18*5-69=21"
        ),
        "exact_star_distance_row_used_to_force_parallel_count": False,
        "opposite_parallel_count_Q": opposite_parallel_count,
        "opposite_edge_count_identity": "9*Q=9*3=27",
        "opposite_edge_count": opposite_edge_count,
        "hard_edge_count_identity": "75-27=48",
        "hard_edge_count": hard_edge_count,
        "hard_sign_times_global_T_identity": "hT=48-27=21",
        "hard_sign_times_global_T": hard_sign_times_T,
        "opposite_normalized_sum_W": conditioned["sum_W"],
        "local_source": "sum q=17*(P_L-3)-18*k_L",
        "common_graph_source": "sum q=hT-P_L=21-P_L",
        "glued_identity": "P_L=4+k_L",
        "hard_rows_k_1_through_4": rows,
        "exact_star_distance_row": [2] * Q,
        "exact_star_energy": 4 * Q,
        "hard_quotient_sum": 12,
        "hard_excess_sum": 3,
        "hard_excess_partitions": [list(row) for row in HARD_EXCESS_PARTITIONS],
        "forced_moment_degrees": [2, 4],
        "local_P_not_equal_4_plus_k_cells_lift_to_common_graph": False,
        "proved": proved,
    }


def _row_parameters(kind: str) -> dict[str, int | None]:
    if kind == "opposite":
        return {
            "excess": None,
            "parallel_count": 3,
            "total": -(HARD_SIGN_TIMES_T + 3),
            "l1_bound": 72,
            "cut_upper": -204,
        }
    if kind.startswith("hard_e"):
        excess = int(kind.removeprefix("hard_e"))
        if excess not in (1, 2, 3):
            raise ValueError(f"unknown hard excess: {excess}")
        k = 1 + excess
        parallel_count = 4 + k
        return {
            "excess": excess,
            "parallel_count": parallel_count,
            "total": HARD_SIGN_TIMES_T - parallel_count,
            "l1_bound": H_EDGE_COUNT - parallel_count,
            "cut_upper": P * M,
        }
    raise ValueError(f"unknown row kind: {kind}")


def _add_common_linear_constraints(
    model: cp_model.CpModel,
    q_values: list[cp_model.IntVar],
    q_abs: list[cp_model.IntVar],
    kind: str,
) -> None:
    parameters = _row_parameters(kind)
    l1_bound = int(parameters["l1_bound"])
    model.Add(sum(q_values) == int(parameters["total"])).WithName("row_sum")
    model.Add(sum(q_abs) <= l1_bound).WithName("l1_bound")

    for degree in (2, 4):
        coefficients = [pow(distance, degree, P) for distance in DISTANCES]
        quotient_bound = max(coefficients) * l1_bound // P + 1
        quotient = model.NewIntVar(
            -quotient_bound,
            quotient_bound,
            f"M{degree}_quotient",
        )
        model.Add(
            sum(
                coefficient * value
                for coefficient, value in zip(coefficients, q_values)
            )
            == P * quotient
        ).WithName(f"M{degree}_mod_{P}")

    cut_upper = int(parameters["cut_upper"])
    for index, cut in enumerate(translated_cut_vectors()):
        model.Add(
            sum(
                coefficient * value
                for coefficient, value in zip(cut, q_values)
            )
            <= cut_upper
        ).WithName(f"translated_cut_{index}")


def _multiplication_energy_model(
    kind: str,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar], cp_model.LinearExpr]:
    """Build the broad-domain optimization encoding without an energy cap."""
    parameters = _row_parameters(kind)
    l1_bound = int(parameters["l1_bound"])
    model = cp_model.CpModel()
    q_values = [
        model.NewIntVar(-l1_bound, l1_bound, f"q_{distance}")
        for distance in DISTANCES
    ]
    q_abs = [
        model.NewIntVar(0, l1_bound, f"q_abs_{distance}")
        for distance in DISTANCES
    ]
    q_square = [
        model.NewIntVar(0, l1_bound * l1_bound, f"q_square_{distance}")
        for distance in DISTANCES
    ]
    for value, absolute, square in zip(q_values, q_abs, q_square):
        model.AddAbsEquality(absolute, value)
        model.AddMultiplicationEquality(square, [value, value])
    _add_common_linear_constraints(model, q_values, q_abs, kind)
    energy = sum(q_square)
    model.Maximize(energy)
    return model, q_values, energy


def _table_threshold_model(
    kind: str,
    forbidden_energy_floor: int | None,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar]]:
    """Independently encode abs/squares by allowed-assignment tables."""
    parameters = _row_parameters(kind)
    l1_bound = int(parameters["l1_bound"])
    model = cp_model.CpModel()
    q_values = [
        model.NewIntVar(-l1_bound, l1_bound, f"q_{distance}")
        for distance in DISTANCES
    ]
    q_abs = [
        model.NewIntVar(0, l1_bound, f"q_abs_{distance}")
        for distance in DISTANCES
    ]
    q_square = [
        model.NewIntVar(0, l1_bound * l1_bound, f"q_square_{distance}")
        for distance in DISTANCES
    ]
    value_table = [
        (value, abs(value), value * value)
        for value in range(-l1_bound, l1_bound + 1)
    ]
    for value, absolute, square in zip(q_values, q_abs, q_square):
        model.AddAllowedAssignments([value, absolute, square], value_table)
    _add_common_linear_constraints(model, q_values, q_abs, kind)
    if forbidden_energy_floor is not None:
        model.Add(sum(q_square) >= forbidden_energy_floor).WithName(
            "forbidden_energy_floor"
        )
    return model, q_values


def _solve_one_worker(model: cp_model.CpModel) -> tuple[cp_model.CpSolver, int]:
    validation = model.Validate()
    _require(not validation, f"invalid CP-SAT model: {validation}")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    return solver, status


def _model_metadata(model: cp_model.CpModel) -> dict[str, object]:
    validation = model.Validate()
    return {
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "workers": 1,
        "seed": 0,
        "model_validation": validation,
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "model_proto_sha256": hashlib.sha256(
            str(model.Proto()).encode("utf-8")
        ).hexdigest(),
    }


def _row_satisfies_constraints(kind: str, row: tuple[int, ...]) -> bool:
    parameters = _row_parameters(kind)
    return bool(
        len(row) == Q
        and sum(row) == parameters["total"]
        and sum(abs(value) for value in row) <= parameters["l1_bound"]
        and all(
            sum(
                pow(distance, degree, P) * value
                for distance, value in zip(DISTANCES, row)
            )
            % P
            == 0
            for degree in (2, 4)
        )
        and all(
            sum(coefficient * value for coefficient, value in zip(cut, row))
            <= parameters["cut_upper"]
            for cut in translated_cut_vectors()
        )
    )


@lru_cache(maxsize=None)
def row_energy_certificate(kind: str) -> dict[str, object]:
    """Optimize one row and replay its sharp threshold independently."""
    if kind not in ("hard_e1", "hard_e2", "hard_e3", "opposite"):
        raise ValueError(f"unknown row kind: {kind}")
    parameters = _row_parameters(kind)

    optimization_model, optimization_q, _ = _multiplication_energy_model(kind)
    optimization_solver, optimization_status = _solve_one_worker(
        optimization_model
    )
    optimization_meta = _model_metadata(optimization_model)
    optimization_meta["status"] = optimization_solver.StatusName(
        optimization_status
    )
    optimization_meta["prior_energy_upper_constraint_used"] = False
    optimization_meta["entry_bounds_used"] = False
    optimization_meta["lower_cut_bounds_used"] = False

    if kind == "hard_e1":
        expected_status = cp_model.INFEASIBLE
        _require(
            optimization_status == expected_status,
            "hard excess-one row unexpectedly became feasible",
        )
        threshold_model, _ = _table_threshold_model(kind, None)
        threshold_solver, threshold_status = _solve_one_worker(threshold_model)
        maximum = None
        maximizer = None
        forbidden_floor = None
    else:
        expected_status = cp_model.OPTIMAL
        expected_maximum = EXPECTED_ROW_MAXIMA[kind]
        _require(
            optimization_status == expected_status
            and int(round(optimization_solver.ObjectiveValue()))
            == expected_maximum,
            f"the {kind} sharp row energy changed",
        )
        maximum = expected_maximum
        maximizer = EXPECTED_MAXIMIZERS[kind]
        _require(
            _row_satisfies_constraints(kind, maximizer)
            and sum(value * value for value in maximizer) == maximum,
            f"the explicit {kind} maximizer changed",
        )
        forbidden_floor = maximum + 1
        threshold_model, _ = _table_threshold_model(kind, forbidden_floor)
        threshold_solver, threshold_status = _solve_one_worker(threshold_model)

    _require(
        threshold_status == cp_model.INFEASIBLE,
        f"the independent {kind} threshold model is no longer infeasible",
    )
    threshold_meta = _model_metadata(threshold_model)
    threshold_meta["status"] = threshold_solver.StatusName(threshold_status)
    threshold_meta["encoding"] = "allowed-assignment tables for (q,abs(q),q^2)"
    threshold_meta["forbidden_energy_floor"] = forbidden_floor
    threshold_meta["prior_energy_upper_constraint_used"] = False
    threshold_meta["entry_bounds_used"] = False
    threshold_meta["lower_cut_bounds_used"] = False

    proved = bool(
        optimization_status == expected_status
        and threshold_status == cp_model.INFEASIBLE
        and (
            kind == "hard_e1"
            or (
                maximum == EXPECTED_ROW_MAXIMA[kind]
                and maximizer == EXPECTED_MAXIMIZERS[kind]
            )
        )
    )
    _require(proved, f"the {kind} row certificate failed")
    cauchy_equality = None
    if kind == "opposite":
        cauchy_lower = int(parameters["total"]) ** 2 // Q
        cauchy_equality = {
            "fixed_sum_energy_lower": cauchy_lower,
            "identity": "8*sum(q_a^2)-(sum q_a)^2=sum_(a<b)(q_a-q_b)^2",
            "equality_requires_all_coordinates_equal": True,
            "forced_equal_coordinate": int(parameters["total"]) // Q,
            "optimization_upper_equals_lower": maximum == cauchy_lower,
            "unique_feasible_row": list(EXPECTED_MAXIMIZERS["opposite"]),
        }
        _require(
            cauchy_lower == 72
            and maximum == cauchy_lower
            and EXPECTED_MAXIMIZERS["opposite"] == (-3,) * Q,
            "the opposite Cauchy equality certificate changed",
        )
    return {
        "kind": kind,
        "parameters": parameters,
        "constraint_scope": [
            "integer eight-bin row",
            "exact row sum",
            "l1 upper bound",
            "M2=0 modulo 17",
            "M4=0 modulo 17",
            "all 698 translated upper-cut inequalities",
        ],
        "constraints_deliberately_not_used": [
            "prior energy upper bound",
            "entry alphabet or aggregate coordinate bounds",
            "lower cut bound",
            "row parity",
            "full 136-entry coefficient matrix",
            "complete-domain Boolean values",
        ],
        "feasible": kind != "hard_e1",
        "sharp_energy_maximum": maximum,
        "explicit_maximizer": list(maximizer) if maximizer is not None else None,
        "fixed_sum_cauchy_equality": cauchy_equality,
        "optimization_model": optimization_meta,
        "independent_threshold_model": threshold_meta,
        "proved": proved,
    }


def _point_scale(scalar: int, point: tuple[int, int]) -> tuple[int, int]:
    return scalar * point[0] % P, scalar * point[1] % P


def _normalize_displacement(point: tuple[int, int]) -> tuple[int, int]:
    negative = ((-point[0]) % P, (-point[1]) % P)
    return min(point, negative)


def _projective_direction(point: tuple[int, int]) -> tuple[int, int]:
    x, y = point
    if x:
        return 1, y * pow(x, -1, P) % P
    if y:
        return 0, 1
    raise ValueError("zero has no projective direction")


def _linear_value(
    functional: tuple[int, int], point: tuple[int, int]
) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % P


def _projected_distance(value: int) -> int:
    return min(value % P, (-value) % P)


@lru_cache(maxsize=1)
def p17_difference_radon_certificate() -> dict[str, object]:
    """Verify the p17 difference-Radon Gram and three Parseval ledgers."""
    points = tuple(product(range(P), repeat=2))
    directions = tuple((1, slope) for slope in range(P)) + ((0, 1),)
    direction_index = {direction: index for index, direction in enumerate(directions)}
    displacements = tuple(
        sorted(
            {
                _normalize_displacement(point)
                for point in points
                if point != (0, 0)
            }
        )
    )
    displacement_directions = tuple(
        direction_index[_projective_direction(displacement)]
        for displacement in displacements
    )
    buckets = tuple(
        tuple(
            _projected_distance(_linear_value(functional, displacement))
            for functional in directions
        )
        for displacement in displacements
    )
    gram_checks = 0
    for left, left_buckets in enumerate(buckets):
        for right, right_buckets in enumerate(buckets):
            actual = sum(
                left_value == right_value
                for left_value, right_value in zip(left_buckets, right_buckets)
            )
            expected = (
                P + 1
                if left == right
                else 1
                if displacement_directions[left] == displacement_directions[right]
                else 2
            )
            _require(actual == expected, "the p17 difference-Radon Gram changed")
            gram_checks += 1

    partition_rows: dict[str, dict[str, object]] = {}
    expected_bases = {
        (1, 1, 1): 1287,
        (2, 1): 1251,
        (3,): 1211,
    }
    for excesses in HARD_EXCESS_PARTITIONS:
        exact_count = HARD_DIRECTION_COUNT - len(excesses)
        hard_parallel_counts = [5] * exact_count + [5 + excess for excess in excesses]
        parallel_counts = hard_parallel_counts + [3] * OPPOSITE_DIRECTION_COUNT
        parallel_square_sum = sum(value * value for value in parallel_counts)
        all_off_bin_base = (
            P * H_EDGE_COUNT
            + 2 * HARD_SIGN_TIMES_T**2
            - 2 * parallel_square_sum
        )
        exact_star_energy = exact_count * 4 * Q
        nonexact_base = all_off_bin_base - exact_star_energy
        key = "+".join(map(str, excesses))
        partition_rows[key] = {
            "hard_excess_partition": list(excesses),
            "hard_quotients": [1] * exact_count
            + [1 + excess for excess in excesses],
            "hard_parallel_counts": hard_parallel_counts,
            "opposite_parallel_counts": [3] * OPPOSITE_DIRECTION_COUNT,
            "parallel_square_sum": parallel_square_sum,
            "exact_hard_star_count": exact_count,
            "exact_hard_star_energy": exact_star_energy,
            "all_off_bin_energy": f"{all_off_bin_base}+34*C",
            "nonexact_off_bin_energy": f"{nonexact_base}+34*C",
            "nonexact_parseval_base": nonexact_base,
        }
        _require(
            sum(hard_parallel_counts) == HARD_EDGE_COUNT
            and sum(parallel_counts) == H_EDGE_COUNT
            and nonexact_base == expected_bases[excesses],
            "a p17 partition Parseval ledger changed",
        )

    proved = bool(
        len(displacements) == 144
        and len(directions) == 18
        and gram_checks == 144 * 144
        and set(partition_rows) == {"1+1+1", "2+1", "3"}
    )
    _require(proved, "the p17 difference-Radon certificate failed")
    return {
        "p": P,
        "difference_class_count": len(displacements),
        "projective_direction_count": len(directions),
        "row_index": "(projective L,a), a in F_17/+/-={0,...,8}",
        "Gram_formula": "B^T*B=17*I+2*J-G_parallel",
        "Gram_entry_values": {
            "same_column": 18,
            "distinct_same_direction": 1,
            "different_directions": 2,
        },
        "Gram_entry_checks": gram_checks,
        "off_bin_parseval": (
            "sum_(L,a>0)q_L(a)^2="
            "17*sum_delta m_delta^2+2*T^2-2*sum_L P_L^2"
        ),
        "collision_parameter": "C=sum_delta binom(m_delta,2)>=0",
        "displacement_square_sum": "sum_delta m_delta^2=75+2*C",
        "partition_ledgers": partition_rows,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15743() -> dict[str, object]:
    """Package the p17 aggregate certificate and exact scope."""
    normalization = two_source_hard_normalization_certificate()
    cuts = p17_cut_catalog_certificate()
    hard_e1 = row_energy_certificate("hard_e1")
    hard_e2 = row_energy_certificate("hard_e2")
    hard_e3 = row_energy_certificate("hard_e3")
    opposite = row_energy_certificate("opposite")
    radon = p17_difference_radon_certificate()

    partition_111_excluded = not hard_e1["feasible"]
    partition_21_excluded = not hard_e1["feasible"]
    partition_3_upper = int(hard_e3["sharp_energy_maximum"]) + (
        OPPOSITE_DIRECTION_COUNT * int(opposite["sharp_energy_maximum"])
    )
    partition_3_parseval_lower = int(
        radon["partition_ledgers"]["3"]["nonexact_parseval_base"]
    )
    partition_3_gap = partition_3_parseval_lower - partition_3_upper
    proved = bool(
        normalization["proved"]
        and cuts["proved"]
        and hard_e1["proved"]
        and hard_e2["proved"]
        and hard_e3["proved"]
        and opposite["proved"]
        and radon["proved"]
        and partition_111_excluded
        and partition_21_excluded
        and partition_3_upper == 767
        and partition_3_parseval_lower == 1211
        and partition_3_gap == 444
    )
    _require(proved, "Proposition 15.743 failed")
    return {
        "prop": "15.743",
        "title": "full translated-cut energy closes the p17 fourth shell",
        "result_status": "exhaustive finite certificate",
        "p": P,
        "layer_index_t": 3,
        "original_k": 4 * P + 6,
        "dependencies": {
            "15.735": "isolated chart and generic branch-B hard/opposite ledger",
            "15.739": "six exact stars force M2=M4=0 and p17 cut reduction",
            "15.741": "difference-Radon common-graph transform, specialized here to p17",
        },
        "two_source_hard_normalization": normalization,
        "translated_cut_catalog": cuts,
        "row_certificates": {
            "hard_excess_one": hard_e1,
            "hard_excess_two": hard_e2,
            "hard_excess_three": hard_e3,
            "opposite": opposite,
        },
        "difference_radon": radon,
        "partition_exclusions": {
            "1+1+1": {
                "reason": "contains a hard excess-one row",
                "excluded": partition_111_excluded,
            },
            "2+1": {
                "reason": "contains a hard excess-one row",
                "excluded": partition_21_excluded,
            },
            "3": {
                "hard_excess_three_energy_upper": hard_e3[
                    "sharp_energy_maximum"
                ],
                "nine_opposite_energy_upper": OPPOSITE_DIRECTION_COUNT
                * int(opposite["sharp_energy_maximum"]),
                "nonexact_energy_upper": partition_3_upper,
                "common_parseval": "1211+34*C",
                "collision_parameter_lower_bound": 0,
                "common_parseval_lower": partition_3_parseval_lower,
                "gap": partition_3_gap,
                "excluded": partition_3_upper < partition_3_parseval_lower,
            },
        },
        "discarded_spectral_cap_or_full_solution_counts_used": False,
        "generic_p17_t3_branch_closed": True,
        "p17_k_eq_74_closed": True,
        "generic_p_ge_17_t3_branch_closed": False,
        "k_eq_4p_plus_6_shell_closed_for_all_primes": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "critical p=5,7; p=11 at k>=50; p=13 at k>=60; p=17 at "
            "k>=76; generic branch-B t=3 from p>=29 and all p>=17 layers "
            "t>=4; multi-level Type I; and the limit"
        ),
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    """Write the deterministic p17 certificate payload."""
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15743.json"
    target.write_text(json.dumps(proposition_15743(), indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    theorem = proposition_15743()
    target = write_evidence()
    print(
        "Prop. 15.743: p17 nonexact energy 767 contradicts "
        "common-graph energy >=1211"
    )
    print(f"result status={theorem['result_status']}")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
