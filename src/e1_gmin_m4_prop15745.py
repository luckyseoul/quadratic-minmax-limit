#!/usr/bin/env python3
r"""Prop. 15.745 -- close the ``u=0`` branch of ``p=13,t=4``.

At ``p=13,k=60`` the odd flip graph has 61 edges.  This proposition treats
only the phase-one common-residue branch ``u=0``.  Its seven hard quotient
variables satisfy ``k_L>=1`` and ``sum k_L=11``.  Gluing an exact literal
star to the common difference transform, before normalizing its distance
row, forces

    hT=17,  P_L=4+k_L.

Thus the hard excess partitions are the five partitions of four, and the
seven opposite parallel counts are ``3^6,4``.  Three exact literal stars
force ``M_2=0``; five force ``M_2=M_4=0``.

All 74 translation-averaged seven-set cuts are imposed on broad six-bin
integer rows.  No prior energy bound, row parity, entry alphabet, lower cut,
or full coefficient matrix is used.  With ``M_2=M_4=0``, an opposite
``Q=3`` row is infeasible, killing the three partitions with at least five
stars.  Under ``M_2=0`` alone the sharp energies are

    hard e=1: 31,  hard e=2: 96,  opposite Q=3: 76,  opposite Q=4: 111.

The ``1^4`` excess partition has upper energy 691, below its exact
difference-Radon baseline ``721+26C``.  For ``(2,1,1)``, the initial upper
bound 725 and baseline ``693+26C`` force ``C<=1``.  But its ``P=7`` hard
direction distributes seven parallel edges among six displacement classes,
so ``C>=1``.  Equality makes that block ``(2,1,1,1,1,1)`` and every other
displacement multiplicity zero or one.  The duplicate is in the zero bin of
the ``e=2`` row.  Each nonzero bin contains one class from every other
projective direction: six have the hard sign and seven the opposite sign.
Consequently every coordinate of that row lies in ``[-7,6]``.  Adding only
this proved equality bound drops its sharp energy to 66, and

    66 + 2*31 + 6*76 + 111 = 695 < 693 + 26 = 719.

Hence the ``u=0`` branch is empty.  This proposition alone does not treat
``u=3,4,6``.  Proposition 15.744 already closes ``u=3``, so the canonical
post-15.745 remainder at ``p=13,k=60`` is exactly ``u=4,6``.
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

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15598 import field_ctx
from e1_gmin_m4_prop15721 import signed_relative_flip_transport
from e1_gmin_m4_prop15734 import (
    BRANCH_P1_LAST,
    baseline_coefficient_rules,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15740 import translated_cut_vector_catalog
from e1_gmin_m4_prop15741 import (
    difference_radon_gram_certificate,
    exact_star_moment_certificate,
    quartic_root_rank_certificate,
)
from e1_gmin_m4_prop15744 import proposition_15744


ROOT = Path(__file__).resolve().parents[1]
P = 13
M = 7
Q = 6
H_EDGE_COUNT = 61
HARD_DIRECTION_COUNT = 7
OPPOSITE_DIRECTION_COUNT = 7
HARD_SIGN_TIMES_T = 17
DISTANCES = tuple(range(1, Q + 1))
HARD_EXCESS_PARTITIONS = (
    (4,),
    (3, 1),
    (2, 2),
    (2, 1, 1),
    (1, 1, 1, 1),
)


ROW_SPECS: dict[str, dict[str, object]] = {
    "opposite_q3_m24": {
        "row_kind": "opposite_q3",
        "total": -20,
        "l1_bound": 58,
        "cut_upper": -130,
        "moment_degrees": (2, 4),
        "coordinate_bounds": None,
        "expected_maximum": None,
        "expected_maximizer": None,
        "optimization_hash": "22d3ff515d20a65de2dea6951a352d31cccec029c3c77ac43c48383e8ea37996",
        "replay_hash": "14c8c33f65c4fff3cc31c5fc1038b51b294da090dcff23bf7dab71566baed814",
    },
    "hard_e1_m2": {
        "row_kind": "hard_e1",
        "total": 11,
        "l1_bound": 55,
        "cut_upper": 91,
        "moment_degrees": (2,),
        "coordinate_bounds": None,
        "expected_maximum": 31,
        "expected_maximizer": (0, 3, 1, 4, 1, 2),
        "optimization_hash": "59525206ec4e18b9ee2662411da86317c6e8cefb26b502d53823ef27bdb08808",
        "replay_hash": "0656dae37434029c854f08d90726b9cb9a9cf546142254e2b8b3f2ccbaff333e",
    },
    "hard_e2_m2": {
        "row_kind": "hard_e2",
        "total": 10,
        "l1_bound": 54,
        "cut_upper": 91,
        "moment_degrees": (2,),
        "coordinate_bounds": None,
        "expected_maximum": 96,
        "expected_maximizer": (5, 0, -3, -2, 3, 7),
        "optimization_hash": "7833446727af6000c5b6d86bff9ba9eb19af995fd44c3a8db1206143c25a9310",
        "replay_hash": "4b8caf7e887b603b61839a474509f3373072841a050beb6f0d1684a3a98568b0",
    },
    "opposite_q3_m2": {
        "row_kind": "opposite_q3",
        "total": -20,
        "l1_bound": 58,
        "cut_upper": -130,
        "moment_degrees": (2,),
        "coordinate_bounds": None,
        "expected_maximum": 76,
        "expected_maximizer": (-5, -3, -2, -5, -2, -3),
        "optimization_hash": "09f908e1d8b1b85132f2f1104955f3e1bf6ab79afd596e0316889fab6e44a96a",
        "replay_hash": "9c099e244b6b341b262fd923186ceb3afcfdfaf92d869b239aa639521679829c",
    },
    "opposite_q4_m2": {
        "row_kind": "opposite_q4",
        "total": -21,
        "l1_bound": 57,
        "cut_upper": -130,
        "moment_degrees": (2,),
        "coordinate_bounds": None,
        "expected_maximum": 111,
        "expected_maximizer": (-4, -6, -1, -3, -7, 0),
        "optimization_hash": "5cf5041f2e1be1cdebd0a2ae7dbb2addb0fa77363bbe36fec04e4cf6c529b969",
        "replay_hash": "4eb57e83ba21f99be79f4e09cd34f7983d9837d6846c42c05014f22022b9badc",
    },
    "hard_e2_c1_m2": {
        "row_kind": "hard_e2",
        "total": 10,
        "l1_bound": 54,
        "cut_upper": 91,
        "moment_degrees": (2,),
        "coordinate_bounds": (-7, 6),
        "expected_maximum": 66,
        "expected_maximizer": (4, 1, -3, 0, 2, 6),
        "optimization_hash": "12d799f142d2be3fbc3c9285e327bb89a07f17e153c3cacfe4bd427c5601dc0e",
        "replay_hash": "b72b57d38790a0eeab87686c3be347771f15adf7263d0ed3eecb8b277703c7d8",
    },
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _proto_hash(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


def _mod_rank(matrix: list[list[int]]) -> int:
    work = [[value % P for value in row] for row in matrix]
    if not work:
        return 0
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, P)
        work[rank] = [(value * inverse) % P for value in work[rank]]
        for row in range(len(work)):
            if row == rank or work[row][column] == 0:
                continue
            scalar = work[row][column]
            work[row] = [
                (left - scalar * right) % P
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == len(work):
            break
    return rank


@lru_cache(maxsize=1)
def moment_root_rank_audit() -> dict[str, object]:
    """Check exactly that three M2 roots and five M4 roots force zero."""
    points = tuple((1, slope) for slope in range(P)) + ((0, 1),)
    checks: dict[int, dict[str, int]] = {}
    for degree, root_count in ((2, 3), (4, 5)):
        ranks = []
        for selected in combinations(points, root_count):
            evaluation = [
                [
                    pow(x, degree - monomial, P) * pow(y, monomial, P) % P
                    for monomial in range(degree + 1)
                ]
                for x, y in selected
            ]
            ranks.append(_mod_rank(evaluation))
        expected_count = 364 if degree == 2 else 2002
        _require(
            len(ranks) == expected_count and set(ranks) == {degree + 1},
            f"the degree-{degree} projective root rank changed",
        )
        checks[degree] = {
            "homogeneous_degree": degree,
            "root_count": root_count,
            "root_count_exceeds_degree": root_count > degree,
            "subsets_checked": len(ranks),
            "common_evaluation_rank": degree + 1,
        }
    return {
        "field": "F_13",
        "degree_2": checks[2],
        "degree_4": checks[4],
        "three_exact_stars_force_M2_zero": True,
        "five_exact_stars_force_M4_zero": True,
        "root_count_rule": (
            "a nonzero homogeneous binary form of degree d has at most d "
            "distinct projective roots"
        ),
        "proved": True,
    }


@lru_cache(maxsize=1)
def translated_cut_vectors() -> tuple[tuple[int, ...], ...]:
    """Import the vetted full 74-vector catalog in its canonical order."""
    catalog = translated_cut_vector_catalog()
    _require(catalog["proved"], "the p13 translated-cut catalog failed")
    vectors = tuple(tuple(int(value) for value in row) for row in catalog["vectors"])
    _require(len(vectors) == 74, "the p13 translated-cut count changed")
    return vectors


@lru_cache(maxsize=1)
def p13_t4_u0_ledger() -> dict[str, object]:
    """Derive the signed normalization, quotient, profile, and moment data."""
    transport = signed_relative_flip_transport()
    floors = residual_even_floor_table(P)
    lift = sharp_integral_quadratic_lift_floor(P)
    literal = baseline_coefficient_rules(P)[BRANCH_P1_LAST]
    stars = exact_star_moment_certificate()
    ranks = quartic_root_rank_certificate()
    root_audit = moment_root_rank_audit()

    ambient_vertices = P * P + 1
    maximum_nonisolated = 2 * H_EDGE_COUNT
    guaranteed_isolated = ambient_vertices - maximum_nonisolated
    type_budget = 2 * M * (M + 4)
    hard_quotient_sum = type_budget // (P + 1)
    hard_excess = hard_quotient_sum - HARD_DIRECTION_COUNT
    phase_one_floors = {
        int(b): int(value) for b, value in floors["phase_one_floors"].items()
    }
    exact_mean = P + 1
    floor_compatible_b_at_exact_mean = tuple(
        b for b, value in phase_one_floors.items() if value <= exact_mean
    )
    b2_excess_above_XNOR = exact_mean - phase_one_floors[2]
    exact_literal_b = P - 1

    # An exact literal star has k=1.  Before inserting its normalized row,
    # local and common sums give 13(P0-3)-14=hT-P0, hence hT=14P0-53.
    # A general hard row then has P_L=P0+k_L-1.  Summing its seven parallel
    # counts and requiring nonnegative opposite edge count bounds P0 by 8.
    exact_parallel_upper = (
        H_EDGE_COUNT - hard_quotient_sum + HARD_DIRECTION_COUNT
    ) // HARD_DIRECTION_COUNT
    exact_parallel_candidates = tuple(
        value
        for value in range(exact_parallel_upper + 1)
        if (value - int(literal["offset"])) % Q == 0
    )
    exact_parallel = exact_parallel_candidates[0]
    hT = 14 * exact_parallel - 53
    hard_parallel_total = (
        HARD_DIRECTION_COUNT * (exact_parallel - 1) + hard_quotient_sum
    )
    opposite_parallel_total = H_EDGE_COUNT - hard_parallel_total

    # For opposite Q, a=14Q-22.  Q=2 has mean six, below both the least
    # nonzero phase-zero parity floor and the b=0 integral-lift floor.
    minimum_opposite_parallel = 3
    opposite_profile = (3,) * 6 + (4,)

    partition_rows = []
    for partition in HARD_EXCESS_PARTITIONS:
        exact_count = HARD_DIRECTION_COUNT - len(partition)
        forced_moments = (2, 4) if exact_count >= 5 else (2,)
        partition_rows.append(
            {
                "hard_excess_partition": list(partition),
                "hard_quotients": [1] * exact_count
                + [1 + value for value in partition],
                "hard_parallel_counts": [5] * exact_count
                + [5 + value for value in partition],
                "exact_literal_star_count": exact_count,
                "forced_global_moment_degrees": list(forced_moments),
            }
        )

    proved = bool(
        transport["proved"]
        and ambient_vertices == 170
        and maximum_nonisolated == 122
        and guaranteed_isolated == 48 > 0
        and type_budget == 154
        and hard_quotient_sum == 11
        and hard_excess == 4
        and min(phase_one_floors.values()) == 12 > 0
        and floor_compatible_b_at_exact_mean == (2, 12)
        and b2_excess_above_XNOR == 2 < lift["sharp_scaled_floor"] == 10
        and exact_literal_b == 12
        and literal["offset"] == 5
        and literal["congruence"] == "6 divides I+P-5"
        and exact_parallel_upper == 8
        and exact_parallel_candidates == (5,)
        and exact_parallel == 5
        and hT == HARD_SIGN_TIMES_T
        and hard_parallel_total == 39
        and opposite_parallel_total == 22
        and floors["least_nonzero_phase_zero_floor"] == 12
        and lift["sharp_scaled_floor"] == 10
        and 14 * 2 - 22 == 6 < 10 < 12
        and minimum_opposite_parallel == 3
        and sum(opposite_profile) == opposite_parallel_total
        and stars["proved"]
        and stars["all_M2_T3_M4_U4_zero"]
        and ranks["proved"]
        and ranks["degree_2_four_roots_force_zero"]
        and ranks["degree_4_four_root_kernel_dimension"] == 1
        and root_audit["proved"]
        and root_audit["three_exact_stars_force_M2_zero"]
        and root_audit["five_exact_stars_force_M4_zero"]
        and [row["exact_literal_star_count"] for row in partition_rows]
        == [6, 5, 5, 4, 3]
        and [row["forced_global_moment_degrees"] for row in partition_rows]
        == [[2, 4], [2, 4], [2, 4], [2], [2]]
    )
    _require(proved, "the p13 t4 u0 ledger changed")
    return {
        "p": P,
        "layer_index_t": 4,
        "original_k": 60,
        "H_edge_count": H_EDGE_COUNT,
        "isolated_chart": {
            "ambient_projective_vertices": ambient_vertices,
            "maximum_nonisolated_vertices": maximum_nonisolated,
            "guaranteed_isolated_vertices": guaranteed_isolated,
            "transported_infinity_degree_I": 0,
        },
        "hard_phase": 1,
        "common_residue_u": 0,
        "hard_type_budget": type_budget,
        "hard_quotient_identity": "k_L>=1 and sum_L k_L=11",
        "hard_k_positive_reason": (
            "phase-one floor is at least 12, so residue-zero mean 14*k "
            "cannot have k=0"
        ),
        "hard_excess_sum": hard_excess,
        "exact_mean_14_floor_compatible_b_values": list(
            floor_compatible_b_at_exact_mean
        ),
        "b2_candidate_excess_above_XNOR": b2_excess_above_XNOR,
        "b2_candidate_excluded_by_integral_lift_floor": True,
        "forced_exact_literal_b": exact_literal_b,
        "exact_literal_unspecialized_identity": "hT=14*P0-53",
        "general_hard_parallel_identity": "P_L=P0+k_L-1",
        "exact_literal_parallel_congruence": literal["congruence"],
        "exact_parallel_upper_from_total_edges": exact_parallel_upper,
        "exact_parallel_candidates": list(exact_parallel_candidates),
        "exact_parallel_count": exact_parallel,
        "hard_sign_times_global_T": hT,
        "general_hard_parallel_count": "P_L=4+k_L=5+e_L",
        "hard_parallel_edge_total": hard_parallel_total,
        "opposite_parallel_edge_total": opposite_parallel_total,
        "minimum_opposite_parallel_count": minimum_opposite_parallel,
        "opposite_parallel_profile": list(opposite_profile),
        "opposite_mean_formula": "a=14*Q-22",
        "hard_cell_formula": "A=7-cut_W",
        "opposite_cell_formula": "A=-10-cut_W; phase zero gives B=A/2",
        "row_formulas": {
            "hard_excess_e": {
                "sum_q": "12-e",
                "l1_bound": "56-e",
                "translated_cut_upper": 91,
            },
            "opposite_Q": {
                "sum_q": "-17-Q",
                "l1_bound": "61-Q",
                "translated_cut_upper": -130,
            },
        },
        "hard_excess_partitions": partition_rows,
        "moment_root_rank_audit": root_audit,
        "exact_literal_distance_row": [2] * Q,
        "exact_literal_distance_energy": 24,
        "proved": proved,
    }


def _build_row_model(
    key: str,
    *,
    table_encoding: bool,
    threshold: int | None,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar], cp_model.LinearExpr]:
    try:
        spec = ROW_SPECS[key]
    except KeyError as exc:
        raise ValueError(f"unknown row key: {key}") from exc
    l1_bound = int(spec["l1_bound"])
    coordinate_bounds = spec["coordinate_bounds"]
    if coordinate_bounds is None:
        lower, upper = -l1_bound, l1_bound
    else:
        lower, upper = (int(value) for value in coordinate_bounds)

    model = cp_model.CpModel()
    q_values = [
        model.NewIntVar(lower, upper, f"q{distance}")
        for distance in DISTANCES
    ]
    q_abs = [
        model.NewIntVar(0, l1_bound, f"a{distance}")
        for distance in DISTANCES
    ]
    q_square = [
        model.NewIntVar(0, l1_bound * l1_bound, f"s{distance}")
        for distance in DISTANCES
    ]
    value_table = [
        (value, abs(value), value * value)
        for value in range(lower, upper + 1)
    ]
    for value, absolute, square in zip(q_values, q_abs, q_square):
        if table_encoding:
            model.AddAllowedAssignments([value, absolute, square], value_table)
        else:
            model.AddAbsEquality(absolute, value)
            model.AddMultiplicationEquality(square, [value, value])

    model.Add(sum(q_values) == int(spec["total"]))
    model.Add(sum(q_abs) <= l1_bound)
    for degree in spec["moment_degrees"]:
        coefficients = [pow(distance, int(degree), P) for distance in DISTANCES]
        quotient_bound = max(coefficients) * l1_bound // P + 2
        quotient = model.NewIntVar(
            -quotient_bound,
            quotient_bound,
            f"m{degree}",
        )
        model.Add(
            sum(
                coefficient * value
                for coefficient, value in zip(coefficients, q_values)
            )
            == P * quotient
        )
    for cut in translated_cut_vectors():
        model.Add(
            sum(
                coefficient * value
                for coefficient, value in zip(cut, q_values)
            )
            <= int(spec["cut_upper"])
        )
    energy = sum(q_square)
    if threshold is None:
        model.Maximize(energy)
    else:
        model.Add(energy >= threshold)
    return model, q_values, energy


def _solve(model: cp_model.CpModel) -> tuple[cp_model.CpSolver, int]:
    validation = model.Validate()
    _require(not validation, f"invalid CP-SAT model: {validation}")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    return solver, status


def _model_metadata(
    model: cp_model.CpModel,
    solver: cp_model.CpSolver,
    status: int,
) -> dict[str, object]:
    return {
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "workers": 1,
        "seed": 0,
        "status": solver.StatusName(status),
        "model_validation": model.Validate(),
        "variables": len(model.Proto().variables),
        "constraints": len(model.Proto().constraints),
        "model_proto_sha256": _proto_hash(model),
    }


def _row_satisfies(key: str, row: Iterable[int]) -> bool:
    spec = ROW_SPECS[key]
    values = tuple(int(value) for value in row)
    coordinate_bounds = spec["coordinate_bounds"]
    return bool(
        len(values) == Q
        and sum(values) == spec["total"]
        and sum(abs(value) for value in values) <= spec["l1_bound"]
        and (
            coordinate_bounds is None
            or all(
                int(coordinate_bounds[0]) <= value <= int(coordinate_bounds[1])
                for value in values
            )
        )
        and all(
            sum(
                pow(distance, int(degree), P) * value
                for distance, value in zip(DISTANCES, values)
            )
            % P
            == 0
            for degree in spec["moment_degrees"]
        )
        and all(
            sum(coefficient * value for coefficient, value in zip(cut, values))
            <= spec["cut_upper"]
            for cut in translated_cut_vectors()
        )
    )


@lru_cache(maxsize=None)
def row_energy_certificate(key: str) -> dict[str, object]:
    """Optimize by tables and replay infeasibility by abs/multiplication."""
    if key not in ROW_SPECS:
        raise ValueError(f"unknown row key: {key}")
    spec = ROW_SPECS[key]
    expected_maximum = spec["expected_maximum"]

    optimization_model, optimization_q, optimization_energy = _build_row_model(
        key,
        table_encoding=True,
        threshold=None,
    )
    optimization_solver, optimization_status = _solve(optimization_model)
    optimization_meta = _model_metadata(
        optimization_model,
        optimization_solver,
        optimization_status,
    )

    if expected_maximum is None:
        _require(
            optimization_status == cp_model.INFEASIBLE,
            f"{key} unexpectedly became feasible",
        )
        maximum = None
        maximizer = None
        replay_threshold = None
    else:
        maximum = int(expected_maximum)
        _require(
            optimization_status == cp_model.OPTIMAL
            and int(round(optimization_solver.ObjectiveValue())) == maximum,
            f"{key} sharp maximum changed",
        )
        maximizer = tuple(int(value) for value in spec["expected_maximizer"])
        _require(
            tuple(optimization_solver.Value(value) for value in optimization_q)
            == maximizer,
            f"{key} deterministic optimizer changed",
        )
        _require(
            _row_satisfies(key, maximizer)
            and sum(value * value for value in maximizer) == maximum,
            f"{key} explicit maximizer failed",
        )
        replay_threshold = maximum + 1

    replay_model, _replay_q, _replay_energy = _build_row_model(
        key,
        table_encoding=False,
        threshold=replay_threshold,
    )
    replay_solver, replay_status = _solve(replay_model)
    replay_meta = _model_metadata(replay_model, replay_solver, replay_status)
    _require(
        replay_status == cp_model.INFEASIBLE,
        f"{key} independent replay is no longer infeasible",
    )
    _require(
        optimization_meta["model_proto_sha256"] == spec["optimization_hash"]
        and replay_meta["model_proto_sha256"] == spec["replay_hash"],
        f"{key} model hash changed",
    )

    return {
        "key": key,
        "row_kind": spec["row_kind"],
        "parameters": {
            "sum": spec["total"],
            "l1_bound": spec["l1_bound"],
            "translated_cut_upper": spec["cut_upper"],
            "moment_degrees_mod_13": list(spec["moment_degrees"]),
            "coordinate_bounds": (
                list(spec["coordinate_bounds"])
                if spec["coordinate_bounds"] is not None
                else None
            ),
        },
        "constraint_scope": [
            "integer six-bin row",
            "exact row sum",
            "l1 upper bound",
            *[
                f"M{degree}=0 modulo 13"
                for degree in spec["moment_degrees"]
            ],
            "all 74 translated upper-cut inequalities",
            *(
                ["C=1-derived coordinate interval [-7,6]"]
                if spec["coordinate_bounds"] is not None
                else []
            ),
        ],
        "constraints_deliberately_not_used": [
            "prior energy upper bound",
            "row parity",
            "local coefficient entry alphabet",
            "lower cut bounds",
            "full 78-entry coefficient matrix",
            "binary midpoint variables",
        ],
        "feasible": expected_maximum is not None,
        "sharp_energy_maximum": maximum,
        "explicit_maximizer": list(maximizer) if maximizer is not None else None,
        "optimization_model": optimization_meta,
        "independent_replay_model": {
            **replay_meta,
            "encoding": "AddAbsEquality and AddMultiplicationEquality",
            "forbidden_energy_floor": replay_threshold,
        },
        "proved": True,
    }


@lru_cache(maxsize=1)
def difference_radon_partition_ledger() -> dict[str, object]:
    """Compute each exact nonstar energy baseline at 61 edges."""
    gram = difference_radon_gram_certificate()
    ledger = p13_t4_u0_ledger()
    rows: dict[str, dict[str, object]] = {}
    expected_bases = {
        (4,): 625,
        (3, 1): 661,
        (2, 2): 665,
        (2, 1, 1): 693,
        (1, 1, 1, 1): 721,
    }
    for partition in HARD_EXCESS_PARTITIONS:
        exact_count = HARD_DIRECTION_COUNT - len(partition)
        hard_parallel = [5] * exact_count + [5 + value for value in partition]
        opposite_parallel = [3] * 6 + [4]
        all_parallel = hard_parallel + opposite_parallel
        parallel_square_sum = sum(value * value for value in all_parallel)
        all_base = (
            P * H_EDGE_COUNT
            + 2 * HARD_SIGN_TIMES_T**2
            - 2 * parallel_square_sum
        )
        exact_energy = exact_count * 24
        nonexact_base = all_base - exact_energy
        key = "+".join(str(value) for value in partition)
        rows[key] = {
            "hard_excess_partition": list(partition),
            "hard_parallel_counts": hard_parallel,
            "opposite_parallel_counts": opposite_parallel,
            "parallel_square_sum": parallel_square_sum,
            "exact_star_count": exact_count,
            "exact_star_energy": exact_energy,
            "all_off_bin_energy": f"{all_base}+26*C",
            "nonexact_off_bin_energy": f"{nonexact_base}+26*C",
            "nonexact_parseval_base": nonexact_base,
        }
        _require(
            sum(all_parallel) == H_EDGE_COUNT
            and nonexact_base == expected_bases[partition],
            "a p13 t4 Radon partition baseline changed",
        )

    proved = bool(
        gram["proved"]
        and gram["Gram_formula"] == "B^T*B=13*I+2*J-G_parallel"
        and ledger["proved"]
        and {key: row["nonexact_parseval_base"] for key, row in rows.items()}
        == {"4": 625, "3+1": 661, "2+2": 665, "2+1+1": 693, "1+1+1+1": 721}
    )
    _require(proved, "the p13 t4 difference-Radon ledger failed")
    return {
        "p": P,
        "Gram_formula": gram["Gram_formula"],
        "collision_parameter": "C=sum_delta binom(m_delta,2)>=0",
        "displacement_square_sum": "sum_delta m_delta^2=61+2*C",
        "off_bin_parseval_formula": (
            "13*61+2*17^2-2*sum_L(P_L^2)+26*C"
        ),
        "partition_ledgers": rows,
        "proved": proved,
    }


def _weak_compositions(total: int, parts: int) -> Iterable[tuple[int, ...]]:
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in _weak_compositions(total - first, parts - 1):
            yield (first,) + tail


def _normalize_displacement(point: tuple[int, int]) -> tuple[int, int]:
    if point == (0, 0):
        raise ValueError("zero is not a displacement class")
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
    functional: tuple[int, int],
    point: tuple[int, int],
) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % P


def _projected_distance(value: int) -> int:
    value %= P
    return min(value, (-value) % P)


@lru_cache(maxsize=1)
def finite_field_bucket_sign_audit() -> dict[str, object]:
    """Check the delicate 6-positive/7-negative nonzero-bin incidence."""
    _q, _mul, _add, chi, _frob, _norm, ia, ib = field_ctx(P)
    points = tuple(product(range(P), repeat=2))
    projective = tuple((1, slope) for slope in range(P)) + ((0, 1),)
    displacements = tuple(
        sorted(
            {
                _normalize_displacement(point)
                for point in points
                if point != (0, 0)
            }
        )
    )

    def field_element(point: tuple[int, int]) -> int:
        return point[0] + P * point[1]

    direction_signs = {
        direction: int(chi(field_element(direction)))
        for direction in projective
    }
    sign_profile = {
        sign: sum(value == sign for value in direction_signs.values())
        for sign in (-1, 1)
    }
    chi_minus_one = int(chi(P - 1))
    base_field_nonzero_characters = tuple(int(chi(value)) for value in range(1, P))
    sign_well_defined_modulo_plus_minus = all(
        int(chi(field_element(point)))
        == int(
            chi(
                field_element(((-point[0]) % P, (-point[1]) % P))
            )
        )
        for point in displacements
    )

    bucket_records = []
    payload_rows = []
    for functional in projective:
        kernel = _projective_direction(
            (functional[1] % P, (-functional[0]) % P)
        )
        epsilon = direction_signs[kernel]
        transverse = set(projective) - {kernel}
        for distance in DISTANCES:
            bucket = tuple(
                point
                for point in displacements
                if _projected_distance(_linear_value(functional, point))
                == distance
            )
            bucket_directions = tuple(
                _projective_direction(point) for point in bucket
            )
            normalized_signs = tuple(
                epsilon * int(chi(field_element(point)))
                for point in bucket
            )
            record = {
                "bucket_size": len(bucket),
                "distinct_transverse_directions": len(set(bucket_directions)),
                "direction_set_is_exactly_transverse": set(bucket_directions)
                == transverse,
                "normalized_plus": normalized_signs.count(1),
                "normalized_minus": normalized_signs.count(-1),
            }
            bucket_records.append(record)
            payload_rows.append(
                f"{functional}:{distance}:"
                + ",".join(
                    f"{point[0]}.{point[1]}.{sign}"
                    for point, sign in zip(bucket, normalized_signs)
                )
            )

    proved = bool(
        len(displacements) == 84
        and sign_profile == {-1: 7, 1: 7}
        and chi_minus_one == 1
        and set(base_field_nonzero_characters) == {1}
        and sign_well_defined_modulo_plus_minus
        and len(bucket_records) == 84
        and all(
            record
            == {
                "bucket_size": 13,
                "distinct_transverse_directions": 13,
                "direction_set_is_exactly_transverse": True,
                "normalized_plus": 6,
                "normalized_minus": 7,
            }
            for record in bucket_records
        )
    )
    _require(proved, "the p13 finite-field bucket/sign incidence changed")
    return {
        "field": "F_13^2",
        "field_model": f"alpha^2={ia}*alpha+{ib}",
        "unoriented_nonzero_displacement_classes": len(displacements),
        "projective_direction_sign_profile": {
            str(sign): count for sign, count in sign_profile.items()
        },
        "chi_minus_one": chi_minus_one,
        "every_nonzero_F13_scalar_is_a_square_in_F13_squared": True,
        "character_sign_well_defined_modulo_plus_minus": (
            sign_well_defined_modulo_plus_minus
        ),
        "hard_normalized_nonzero_buckets_checked": len(bucket_records),
        "classes_per_nonzero_bucket": 13,
        "one_class_from_each_transverse_direction": True,
        "normalized_hard_sign_classes_per_bucket": 6,
        "normalized_opposite_sign_classes_per_bucket": 7,
        "incidence_sign_sha256": hashlib.sha256(
            ";".join(payload_rows).encode("ascii")
        ).hexdigest(),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def collision_one_coordinate_bound() -> dict[str, object]:
    """Turn ``C=1`` into the exact ``[-7,6]`` e2 coordinate interval."""
    bucket_audit = finite_field_bucket_sign_audit()
    collision_rows = [
        (row, sum(value * (value - 1) // 2 for value in row))
        for row in _weak_compositions(7, 6)
    ]
    minimum_collision = min(value for _row, value in collision_rows)
    minimizers = {
        tuple(sorted(row))
        for row, value in collision_rows
        if value == minimum_collision
    }
    hard_sign_directions = int(
        bucket_audit["projective_direction_sign_profile"]["1"]
    )
    opposite_sign_directions = int(
        bucket_audit["projective_direction_sign_profile"]["-1"]
    )
    nonzero_bucket_size = int(bucket_audit["classes_per_nonzero_bucket"])
    plus_classes = int(bucket_audit["normalized_hard_sign_classes_per_bucket"])
    minus_classes = int(
        bucket_audit["normalized_opposite_sign_classes_per_bucket"]
    )
    coordinate_interval = (-minus_classes, plus_classes)
    proved = bool(
        len(collision_rows) == 792
        and minimum_collision == 1
        and minimizers == {(1, 1, 1, 1, 1, 2)}
        and bucket_audit["proved"]
        and bucket_audit["chi_minus_one"] == 1
        and hard_sign_directions == opposite_sign_directions == 7
        and plus_classes == 6
        and minus_classes == 7
        and plus_classes + minus_classes == nonzero_bucket_size
        and coordinate_interval == (-7, 6)
    )
    _require(proved, "the C=1 coordinate bound changed")
    return {
        "parallel_count_in_e2_direction": 7,
        "parallel_displacement_class_count": 6,
        "weak_composition_count_checked": len(collision_rows),
        "minimum_collision_contribution": minimum_collision,
        "unique_sorted_equality_profile": [1, 1, 1, 1, 1, 2],
        "global_C_upper_from_broad_energy": 1,
        "therefore_global_C": 1,
        "all_displacement_multiplicities_outside_e2_parallel_block_are_0_or_1": True,
        "duplicate_class_is_in_e2_zero_bin": True,
        "finite_field_bucket_sign_audit": bucket_audit,
        "each_e2_nonzero_bin_has_one_class_from_each_other_direction": True,
        "hard_sign_classes_per_nonzero_bin": plus_classes,
        "opposite_sign_classes_per_nonzero_bin": minus_classes,
        "hard_normalized_coordinate_interval": list(coordinate_interval),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15745() -> dict[str, object]:
    """Package the complete aggregate close of the p13 t4 u0 branch."""
    ledger = p13_t4_u0_ledger()
    cuts = translated_cut_vector_catalog()
    radon = difference_radon_partition_ledger()
    q3_m24 = row_energy_certificate("opposite_q3_m24")
    h1 = row_energy_certificate("hard_e1_m2")
    h2 = row_energy_certificate("hard_e2_m2")
    q3 = row_energy_certificate("opposite_q3_m2")
    q4 = row_energy_certificate("opposite_q4_m2")
    collision = collision_one_coordinate_bound()
    h2_c1 = row_energy_certificate("hard_e2_c1_m2")
    prior_u3 = proposition_15744()

    first_three_partitions = ((4,), (3, 1), (2, 2))
    first_three_killed = bool(
        not q3_m24["feasible"]
        and all(
            HARD_DIRECTION_COUNT - len(partition) >= 5
            for partition in first_three_partitions
        )
    )
    three_star_upper = 4 * 31 + 6 * 76 + 111
    three_star_base = int(
        radon["partition_ledgers"]["1+1+1+1"]["nonexact_parseval_base"]
    )
    four_star_broad_upper = 96 + 2 * 31 + 6 * 76 + 111
    four_star_base = int(
        radon["partition_ledgers"]["2+1+1"]["nonexact_parseval_base"]
    )
    collision_upper = (four_star_broad_upper - four_star_base) // 26
    four_star_c1_upper = 66 + 2 * 31 + 6 * 76 + 111
    four_star_c1_lower = four_star_base + 26

    proved = bool(
        ledger["proved"]
        and cuts["proved"]
        and cuts["distinct_translated_cut_vector_count"] == 74
        and radon["proved"]
        and first_three_killed
        and h1["sharp_energy_maximum"] == 31
        and h2["sharp_energy_maximum"] == 96
        and q3["sharp_energy_maximum"] == 76
        and q4["sharp_energy_maximum"] == 111
        and three_star_upper == 691 < three_star_base == 721
        and four_star_broad_upper == 725
        and four_star_base == 693
        and collision_upper == 1
        and collision["proved"]
        and collision["minimum_collision_contribution"] == 1
        and h2_c1["sharp_energy_maximum"] == 66
        and four_star_c1_upper == 695 < four_star_c1_lower == 719
        and prior_u3["proved"]
        and prior_u3["p13_t4_u3_branch_closed"]
        and prior_u3["remaining_p13_t4_residues"] == [0, 4, 6]
    )
    _require(proved, "Proposition 15.745 failed")
    return {
        "prop": "15.745",
        "result_status": "exhaustive finite aggregate certificate",
        "scope": "p=13,t=4,k=60, phase-one common residue u=0 only",
        "ledger": ledger,
        "translated_cut_catalog": {
            "distinct_vectors": cuts["distinct_translated_cut_vector_count"],
            "catalog_sha256": cuts["catalog_sha256"],
            "all_74_vectors_used_in_every_row_model": True,
        },
        "row_certificates": {
            key: row
            for key, row in (
                ("opposite_q3_m24", q3_m24),
                ("hard_e1_m2", h1),
                ("hard_e2_m2", h2),
                ("opposite_q3_m2", q3),
                ("opposite_q4_m2", q4),
                ("hard_e2_c1_m2", h2_c1),
            )
        },
        "difference_radon": radon,
        "partition_close": {
            "five_or_more_exact_stars": {
                "partitions": [list(row) for row in first_three_partitions],
                "M2_M4_forced_zero": True,
                "opposite_Q3_row_infeasible": True,
                "closed": first_three_killed,
            },
            "three_exact_stars_1+1+1+1": {
                "row_energy_upper": "4*31+6*76+111=691",
                "exact_Radon_energy": "721+26*C>=721",
                "closed": three_star_upper < three_star_base,
            },
            "four_exact_stars_2+1+1": {
                "initial_row_energy_upper": "96+2*31+6*76+111=725",
                "exact_Radon_energy": "693+26*C",
                "initial_consequence": "C<=1",
                "parallel_P7_consequence": "C>=1",
                "forced_collision_value": 1,
                "collision_one_coordinate_certificate": collision,
                "tightened_row_energy_upper": "66+2*31+6*76+111=695",
                "tightened_exact_Radon_energy": "693+26=719",
                "closed": four_star_c1_upper < four_star_c1_lower,
            },
        },
        "p13_t4_u0_closed": True,
        "p13_k_eq_60_closed": False,
        "not_addressed_by_prop_15745": [3, 4, 6],
        "prior_prop_15744_u3_closed": True,
        "remaining_p13_t4_residues": [4, 6],
        "proved": proved,
    }


def main() -> None:
    result = proposition_15745()
    output = ROOT / "evidence" / "e1_gmin_m4_prop15745.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
