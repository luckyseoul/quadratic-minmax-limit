#!/usr/bin/env python3
r"""Prop. 15.744 -- close the exceptional ``p=13,t=4,u=3`` branch.

At ``p=13,t=4`` the phase-one budget is 154.  In residue ``u=3`` the
seven hard quotients have sum eight.  The exact floor table and the sharp
lift floor therefore force the profile ``1^6 2``: six exact
complement-triple cells of mean 20 and one elevated cell of mean 34.

For an exact hard row with parallel count ``P``, local coefficient
normalization gives

    20 = 14P-hT-39,       so hT=14P-59.

The corrected complement-triple offset is two, hence ``6 | P-2``.  The
elevated hard row has parallel count ``R=P+1``.  Edge accounting leaves
only ``P=2`` and ``P=8`` and gives opposite parallel-count sums 46 and 4.
The opposite mean is ``a(Q)=14Q+hT-39``.  Thus the two ledgers force a
mass-14 opposite cell at ``Q=6`` or ``Q=0``.

Proposition 15.738 cannot be imported verbatim here: its height-four
coefficient models used the old edge count 59 and the bound
``sum |W| <= 59-Q``.  This module rebuilds those models at the changed edge
count 61, with the necessary relaxed bound ``sum |W| <= 61-Q``.  Both
one-worker models are still infeasible.  The edge-count-independent Boolean
support-462 classification from Proposition 15.738 then leaves only
``B=x_i*x_j``.

The six exact hard complement triples are roots of the homogeneous binary
quartic ``G=2hM_4-M_2^2``.  Six is greater than its degree four, so
``G`` vanishes identically.  The selected-pair opposite cell instead gives
``G=-3(i-j)^4 != 0`` in ``F_13``.  This closes exactly the ``u=3`` branch;
the other ``p=13,t=4`` residues and residual (ii) remain open.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    baseline_coefficient_rules,
)
from e1_gmin_m4_prop15738 import (
    MASS14_SUPPORT_SIZE,
    exact_mass14_boolean_classification,
    mass14_boolean_catalog_certificate,
    mass14_height_dichotomy,
    middle_slice_points,
    pair_coordinates,
    selected_third_difference_identities,
    selected_pair_moment_certificate,
    third_difference_rank_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
P = 13
M = 7
Q_MODULUS = 6
LAYER_INDEX = 4
ORIGINAL_K = 60
H_EDGE_COUNT = 61
HARD_DIRECTION_COUNT = 7
OPPOSITE_DIRECTION_COUNT = 7
HARD_RESIDUE_U = 3
EXACT_HARD_MEAN = 20
ELEVATED_HARD_MEAN = 34
EXACT_HARD_COUNT = 6
EXACT_HARD_PARALLEL_COUNTS = (2, 8)
FORCED_OPPOSITE_PARALLEL_COUNTS = (0, 6)
SEARCH_WORKERS = 1

EXPECTED_HEIGHT_FOUR_MODEL_SHA256 = {
    0: "70313e414ca6da2cf6694c11bdd7c7ee8ee985ca05bd30802aa2b6b96353d3d3",
    6: "a94796122b2c1a115b1efec4094031726f118e05fd6825f2e805a406b4f2b9dd",
}
EXPECTED_B10_CONTACT_LAYER_MATRIX_SHA256 = (
    "996269be45189565eaf8717f97f71f2e2f22ad33c8116da5fd2e154ec8eaf695"
)
EXPECTED_B10_PUNCTURED_LIFT_MODEL_SHA256 = (
    "b0d1956f0a173f7c4ce94d7f588af92311d42f5017dd72d814d361e964b6bcd4"
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _model_textproto_sha256(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


def _modular_rank(matrix: list[list[int]], modulus: int) -> tuple[int, list[int]]:
    """Return deterministic row rank and pivot columns over ``F_modulus``."""
    work = [[value % modulus for value in row] for row in matrix]
    if not work:
        return 0, []
    rank = 0
    pivots: list[int] = []
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, modulus)
        work[rank] = [value * inverse % modulus for value in work[rank]]
        for row in range(rank + 1, len(work)):
            scalar = work[row][column]
            if scalar:
                work[row] = [
                    (left - scalar * right) % modulus
                    for left, right in zip(work[row], work[rank])
                ]
        pivots.append(column)
        rank += 1
        if rank == len(work[0]):
            break
    return rank, pivots


@lru_cache(maxsize=1)
def b10_contact_layer_restriction_certificate() -> dict[str, object]:
    """Prove that the three positive contact layers determine a slice quadratic.

    On ``J(13,7)``, pair monomials span every function of degree at most two:
    constants and linears follow from the fixed-weight identities
    ``sum_(i<j)x_i*x_j=21`` and ``sum_(j!=i)x_i*x_j=6*x_i``.  Full column
    rank after restricting the 78 pair monomials to ``r=1,2,3`` therefore
    makes that restriction injective over the rationals.
    """
    complement = frozenset((0, 1, 2))
    points = middle_slice_points()
    pairs = pair_coordinates()
    layer_counts = Counter(len(set(point) & complement) for point in points)
    rows = [
        [int(i in point and j in point) for i, j in pairs]
        for point in points
        if len(set(point) & complement) in (1, 2, 3)
    ]
    modulus = 101
    rank, pivot_columns = _modular_rank(rows, modulus)
    payload = ";".join("".join(map(str, row)) for row in rows).encode("ascii")
    matrix_hash = hashlib.sha256(payload).hexdigest()
    proved = bool(
        len(points) == comb(P, M) == 1716
        and len(pairs) == comb(P, 2) == 78
        and layer_counts == Counter({0: 120, 1: 630, 2: 756, 3: 210})
        and len(rows) == 1596
        and modulus == 101
        and rank == len(pairs) == 78
        and pivot_columns == list(range(78))
        and matrix_hash == EXPECTED_B10_CONTACT_LAYER_MATRIX_SHA256
    )
    _require(proved, "the b=10 contact-layer restriction rank changed")
    return {
        "slice": "J(13,7)",
        "three_point_complement": [0, 1, 2],
        "intersection_layer_counts": {
            str(layer): layer_counts[layer] for layer in range(4)
        },
        "positive_contact_layers": [1, 2, 3],
        "restricted_point_count": len(rows),
        "pair_monomial_count": len(pairs),
        "pair_monomials_span_degree_at_most_two": True,
        "fixed_weight_spanning_identities": [
            "sum_(i<j)x_i*x_j=21",
            "sum_(j!=i)x_i*x_j=6*x_i",
        ],
        "rank_modulus": modulus,
        "restricted_evaluation_rank": rank,
        "pivot_columns": pivot_columns,
        "matrix_sha256": matrix_hash,
        "full_rank_mod_prime_implies_full_rank_over_Q": True,
        "vanishing_on_contact_layers_forces_zero_globally": True,
        "proved": proved,
    }


def _build_b10_punctured_lift_model() -> tuple[cp_model.CpModel, dict[str, int]]:
    """Build the exact two-unit lift model around ``A_0=(2-r)^2``.

    Writing ``A=A_0+2B``, parity and nonnegativity give ``B>=0`` on the
    three contact layers but only ``B>=-2`` on ``r=0``.  Thus Proposition
    15.688's globally nonnegative lift theorem is inapplicable; this model
    handles precisely the missing punctured case.
    """
    complement = frozenset((0, 1, 2))
    points = middle_slice_points()
    identities, _descriptors, _examined = selected_third_difference_identities()
    target_sum = 66  # 4p E[B]=2 on |J(13,7)|=1716.
    omitted_layer_size = sum(
        len(set(point) & complement) == 0 for point in points
    )
    upper_bound = target_sum + 2 * omitted_layer_size

    model = cp_model.CpModel()
    values = []
    for index, point in enumerate(points):
        intersection = len(set(point) & complement)
        lower_bound = -2 if intersection == 0 else 0
        values.append(
            model.NewIntVar(lower_bound, upper_bound, f"B_{index}")
        )
    for row in identities:
        model.Add(sum(sign * values[index] for index, sign in row) == 0)
    model.Add(sum(values) == target_sum)
    return model, {
        "value_variable_count": len(values),
        "third_difference_identity_count": len(identities),
        "omitted_layer_size": omitted_layer_size,
        "target_sum": target_sum,
        "safe_upper_bound": upper_bound,
    }


@lru_cache(maxsize=1)
def b10_floor_plus_two_exclusion() -> dict[str, object]:
    """Exclude the phase-one ``b=10`` cell of scaled mean 22 exactly."""
    rank = b10_contact_layer_restriction_certificate()
    degree_two = third_difference_rank_certificate()
    model, data = _build_b10_punctured_lift_model()
    validation = model.Validate()
    _require(not validation, f"invalid b=10 punctured-lift model: {validation}")
    model_hash = _model_textproto_sha256(model)
    _require(
        model_hash == EXPECTED_B10_PUNCTURED_LIFT_MODEL_SHA256,
        "the b=10 punctured-lift model hash changed",
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = SEARCH_WORKERS
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    proved = bool(
        rank["proved"]
        and degree_two["proved"]
        and degree_two["selected_identity_count"] == 1638
        and degree_two["exact_real_rank"] == 1638
        and degree_two["exact_real_nullity"] == 78
        and data
        == {
            "value_variable_count": 1716,
            "third_difference_identity_count": 1638,
            "omitted_layer_size": 120,
            "target_sum": 66,
            "safe_upper_bound": 306,
        }
        and infeasible
    )
    _require(proved, f"b=10 floor-plus-two model is {status_name}")
    return {
        "p": P,
        "original_phase_and_boundary": {"phase": 1, "b": 10},
        "complement_reduction": {"phase": 0, "b": 3},
        "exact_baseline": "A_0=(2-r)^2",
        "difference_variable": "B=(A-A_0)/2",
        "scaled_mean_assumption": "2p*E[A]=22",
        "scaled_baseline_mean": 20,
        "scaled_difference_identity": "4p*E[B]=2",
        "slice_point_count": data["value_variable_count"],
        "integer_value_sum": data["target_sum"],
        "contact_layer_lower_bound": "B>=0 for r=1,2,3",
        "omitted_layer_lower_bound": "B>=-2 for r=0",
        "omitted_layer_size": data["omitted_layer_size"],
        "safe_coordinate_upper_bound": data["safe_upper_bound"],
        "global_nonnegative_lift_theorem_used": False,
        "third_difference_dependency": {
            "proposition": "15.738",
            "identity_count": data["third_difference_identity_count"],
            "identity_rank": degree_two["exact_real_rank"],
            "degree_two_nullity": degree_two["exact_real_nullity"],
            "proved": degree_two["proved"],
        },
        "contact_layer_restriction_dependency": rank,
        "model_validation": validation,
        "model_variable_count": len(model.Proto().variables),
        "model_constraint_count": len(model.Proto().constraints),
        "model_textproto_sha256": model_hash,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "search_workers": SEARCH_WORKERS,
        "random_seed": 0,
        "symmetry_level": 3,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "punctured_lift_model_infeasible": infeasible,
        "result_status": "exhaustive finite certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def t4_phase_one_baseline_dependencies() -> dict[str, object]:
    """Certify the three exact low baselines used by the all-residue sieve."""
    rules = baseline_coefficient_rules(P)
    # A phase-one b=10 row is indexed by the three-point complement C.
    # Since |X|=7, its parity at r=|X intersect C| is
    # (|X intersect B|+1) mod 2=(7-r+1) mod 2=r mod 2.  It is therefore
    # exactly the phase-zero b=3 problem certified by Proposition 15.652.
    triple_quadrature = parity_floor_certificate(P, 3, 0)
    triple_nodes = tuple(triple_quadrature["quadrature_nodes"])
    triple_weights = tuple(triple_quadrature["quadrature_weights"])
    triple_coefficients = tuple(triple_quadrature["coefficients"])
    triple_restriction = b10_contact_layer_restriction_certificate()
    triple_complement_parity_checks = [
        ((M - r + 1) & 1) == (r & 1) for r in range(4)
    ]
    triple_positive_quadrature_rigidity = bool(
        triple_quadrature["exact_positive_quadrature_certificate"]
        and int(triple_quadrature["scaled_floor"]) == EXACT_HARD_MEAN
        and triple_coefficients
        == (Fraction(1), Fraction(-4), Fraction(4))
        and triple_nodes == (1, 2, 3)
        and all(weight > 0 for weight in triple_weights)
        and all(triple_complement_parity_checks)
        and triple_restriction["proved"]
        and triple_restriction[
            "vanishing_on_contact_layers_forces_zero_globally"
        ]
    )
    triple_scaled_mean_numerator = sum(
        comb(3, r) * comb(10, M - r) * (2 - r) ** 2
        for r in range(4)
        if 0 <= M - r <= 10
    )
    domain_size = comb(P, M)
    triple_scaled_mean = (
        2 * P * triple_scaled_mean_numerator // domain_size
    )
    triple_mean_divides_exactly = (
        (2 * P * triple_scaled_mean_numerator) % domain_size == 0
    )
    triple_target_checks = []
    for bits in product((0, 1), repeat=3):
        r = sum(bits)
        z = tuple(2 * bit - 1 for bit in bits)
        triple_target_checks.append(
            3 + 2 * (2 - r) ** 2
            == 5
            - sum(z)
            + z[0] * z[1]
            + z[0] * z[2]
            + z[1] * z[2]
        )
    proved = bool(
        rules["proved"]
        and rules[BRANCH_B2]["baseline"] == "A=(1-x_i-x_j)^2"
        and rules[BRANCH_P1_LAST]["baseline"] == "A=x_j"
        and triple_positive_quadrature_rigidity
        and triple_mean_divides_exactly
        and triple_scaled_mean == 20
        and all(triple_target_checks)
    )
    _require(proved, "the p13 phase-one baseline dependencies failed")
    return {
        "p": P,
        "b2_exact_baseline": rules[BRANCH_B2]["baseline"],
        "b2_exact_scaled_mean": 12,
        "b2_positive_quadrature_rigidity": rules[
            "b2_phase_one_equality_is_pointwise_XNOR"
        ],
        "b12_exact_baseline": rules[BRANCH_P1_LAST]["baseline"],
        "b12_exact_scaled_mean": 14,
        "b12_positive_quadrature_rigidity": rules[
            "b_p_minus_one_phase_one_equality_is_pointwise_literal"
        ],
        "b10_exact_baseline": "A=(2-r)^2 for r=|X intersect C|, |C|=3",
        "b10_exact_scaled_mean": triple_scaled_mean,
        "b10_target_boolean_checks": len(triple_target_checks),
        "b10_complement_parity_checks": len(triple_complement_parity_checks),
        "b10_positive_quadrature": {
            "proposition": "15.652",
            "reduced_boundary_size": 3,
            "reduced_phase": 0,
            "coefficients": [int(value) for value in triple_coefficients],
            "contact_nodes": list(triple_nodes),
            "weights": [str(weight) for weight in triple_weights],
            "all_weights_strictly_positive": all(
                weight > 0 for weight in triple_weights
            ),
            "positive_weights_force_pointwise_contact_on_all_three_layers": True,
            "contact_layers_determine_slice_quadratic": triple_restriction[
                "vanishing_on_contact_layers_forces_zero_globally"
            ],
            "exact_positive_quadrature_certificate": bool(
                triple_quadrature["exact_positive_quadrature_certificate"]
            ),
        },
        "b10_contact_layer_restriction": triple_restriction,
        "b10_positive_quadrature_rigidity": (
            triple_positive_quadrature_rigidity
        ),
        "positive_quadrature_dependency": "Proposition 15.652",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def t4_all_residue_sieve() -> dict[str, object]:
    """Apply the exact p13 floor/lift sieve to all seven ``t=4`` residues."""
    floors = {b: scaled_direction_floor(P, b, 1) for b in range(0, P, 2)}
    baselines = t4_phase_one_baseline_dependencies()
    b10_floor_plus_two = b10_floor_plus_two_exclusion()
    lift = sharp_integral_quadratic_lift_floor(P)
    lift_floor = int(lift["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for u in range(M):
        quotient_sum = M + LAYER_INDEX - u
        if u <= LAYER_INDEX:
            low_quotient = 1
            low_mean = 2 * u + (P + 1)
            excess_units = quotient_sum - M
            forced_low_count = M - excess_units
            k0_mean = 2 * u
            k0_below_every_floor = k0_mean < min(floors.values())
        else:
            low_quotient = 0
            low_mean = 2 * u
            forced_low_count = M - quotient_sum
            k0_mean = low_mean
            k0_below_every_floor = low_mean < min(floors.values())

        low_cells: list[dict[str, object]] = []
        for b, floor in floors.items():
            if floor > low_mean:
                continue
            excess = low_mean - floor
            ordinary_nonnegative_lift = b in (2, 12)
            forbidden_by_sharp_lift = bool(
                ordinary_nonnegative_lift and 0 < excess < lift_floor
            )
            forbidden_by_punctured_lift = bool(
                b == 10 and excess == 2 and b10_floor_plus_two["proved"]
            )
            forbidden_lift = (
                forbidden_by_sharp_lift or forbidden_by_punctured_lift
            )
            if excess == 0:
                status = "exact baseline"
            elif forbidden_by_punctured_lift:
                status = "punctured floor-plus-two model infeasible"
            elif forbidden_by_sharp_lift:
                status = "forbidden lift below 10"
            elif excess == lift_floor:
                status = "sharp lift equality survives"
            else:
                status = "unclassified positive lift survives"
            low_cells.append(
                {
                    "b": b,
                    "floor": floor,
                    "excess": excess,
                    "status": status,
                    "survives": not forbidden_lift,
                }
            )
        surviving_low_cells = [
            cell for cell in low_cells if bool(cell["survives"])
        ]
        excluded = not surviving_low_cells
        rows.append(
            {
                "u": u,
                "common_residue": 2 * u,
                "quotient_sum": quotient_sum,
                "k0_mean": k0_mean,
                "k0_below_every_phase_one_floor": k0_below_every_floor,
                "forced_low_quotient": low_quotient,
                "forced_low_direction_count_at_least": forced_low_count,
                "forced_low_mean": low_mean,
                "low_cells_at_or_below_mean": low_cells,
                "surviving_low_cells": surviving_low_cells,
                "excluded": excluded,
            }
        )

    survivors = [int(row["u"]) for row in rows if not row["excluded"]]
    expected_surviving_cells = {
        0: [(12, 14, 0)],
        3: [(10, 20, 0)],
        4: [(2, 12, 10)],
        6: [(2, 12, 0)],
    }
    actual_surviving_cells = {
        int(row["u"]): [
            (int(cell["b"]), int(cell["floor"]), int(cell["excess"]))
            for cell in row["surviving_low_cells"]
        ]
        for row in rows
        if not row["excluded"]
    }
    proved = bool(
        floors == {0: 26, 2: 12, 4: 26, 6: 24, 8: 26, 10: 20, 12: 14}
        and lift["proved"]
        and baselines["proved"]
        and b10_floor_plus_two["proved"]
        and lift_floor == 10
        and survivors == [0, 3, 4, 6]
        and actual_surviving_cells == expected_surviving_cells
        and [row["forced_low_direction_count_at_least"] for row in rows[:5]]
        == [3, 4, 5, 6, 7]
        and all(row["k0_below_every_phase_one_floor"] for row in rows[:5])
        and rows[5]["forced_low_quotient"] == 0
        and rows[5]["forced_low_direction_count_at_least"] == 1
        and rows[5]["forced_low_mean"] == 10
        and rows[6]["forced_low_quotient"] == 0
        and rows[6]["forced_low_direction_count_at_least"] == 2
        and rows[6]["forced_low_mean"] == 12
    )
    _require(proved, "the p=13,t=4 all-residue sieve changed")
    return {
        "p": P,
        "layer_index_t": LAYER_INDEX,
        "hard_mean_form": "a=2u+14k",
        "hard_quotient_identity": "sum k=11-u",
        "phase_one_even_b_floors": {str(b): value for b, value in floors.items()},
        "minimum_phase_one_even_b_floor": min(floors.values()),
        "exact_low_baseline_dependencies": baselines,
        "sharp_nonzero_integral_lift_floor": lift_floor,
        "sharp_lift_dependency": "Proposition 15.688 for b=2,12 only",
        "sharp_lift_certificate_called": True,
        "b10_floor_plus_two_exclusion": b10_floor_plus_two,
        "residue_rows": rows,
        "surviving_residues_before_prop_15744": survivors,
        "excluded_residues_before_prop_15744": [1, 2, 5],
        "proved": proved,
    }


@lru_cache(maxsize=1)
def t4_u3_residue_ledger() -> dict[str, object]:
    """Derive the exact ``1^6 2`` hard profile and both edge ledgers."""
    all_residues = t4_all_residue_sieve()
    type_budget = 2 * M * (M + LAYER_INDEX)
    quotient_sum = M + LAYER_INDEX - HARD_RESIDUE_U
    phase_one_floors = {
        b: scaled_direction_floor(P, b, 1) for b in range(0, P, 2)
    }
    lift_floor = int(sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"])
    low_mean = 2 * HARD_RESIDUE_U + (P + 1)

    low_rows = []
    for b, floor in phase_one_floors.items():
        if floor > low_mean:
            continue
        excess = low_mean - floor
        if b in (2, P - 1) and 0 < excess < lift_floor:
            status = "forbidden lift below 10"
        elif b == 10 and excess == 0:
            status = "exact complement-triple cell"
        else:
            status = "not retained"
        low_rows.append({"b": b, "floor": floor, "excess": excess, "status": status})

    quotient_profile = [1] * 6 + [2]
    target_checks = []
    for bits in product((0, 1), repeat=3):
        r = sum(bits)
        z = tuple(2 * bit - 1 for bit in bits)
        target_checks.append(
            3 + 2 * (2 - r) ** 2
            == 5
            - sum(z)
            + z[0] * z[1]
            + z[0] * z[2]
            + z[1] * z[2]
        )

    coefficient_offset = 5 - 3
    candidate_rows: dict[str, dict[str, object]] = {}
    for exact_parallel_count in EXACT_HARD_PARALLEL_COUNTS:
        hard_sign_times_T = 14 * exact_parallel_count - 59
        elevated_parallel_count = (
            ELEVATED_HARD_MEAN + hard_sign_times_T + 39
        ) // 14
        hard_edge_count = (
            EXACT_HARD_COUNT * exact_parallel_count + elevated_parallel_count
        )
        opposite_parallel_sum = H_EDGE_COUNT - hard_edge_count

        def opposite_mean(parallel_count: int) -> int:
            return 14 * parallel_count + hard_sign_times_T - 39

        if exact_parallel_count == 2:
            q5_offset = 3
            q5_compatible = (5 - q5_offset) % Q_MODULUS == 0
            minimum_allowed_q = 6
            forced_q = 6
            minimum_count = 1
            forcing_reason = (
                "Q=5,a=0 has offset 3 and is incompatible modulo 6; "
                "sum Q=46<7*7 forces a Q=6 row"
            )
        else:
            q5_compatible = None
            minimum_allowed_q = 0
            forced_q = 0
            minimum_count = OPPOSITE_DIRECTION_COUNT - opposite_parallel_sum
            forcing_reason = "sum Q=4 across seven nonnegative rows forces a Q=0 row"

        row = {
            "exact_hard_parallel_count_P": exact_parallel_count,
            "hard_sign_times_global_T": hard_sign_times_T,
            "exact_hard_identity": (
                f"20=14*{exact_parallel_count}-({hard_sign_times_T})-39"
            ),
            "elevated_hard_parallel_count_R": elevated_parallel_count,
            "elevated_identity": "34=14*R-hT-39",
            "R_equals_P_plus_one": elevated_parallel_count
            == exact_parallel_count + 1,
            "hard_edge_count": hard_edge_count,
            "hard_edge_count_identity": (
                f"6*{exact_parallel_count}+{elevated_parallel_count}"
            ),
            "opposite_parallel_count_sum": opposite_parallel_sum,
            "opposite_mean_formula": (
                f"a(Q)=14*Q+({hard_sign_times_T})-39"
            ),
            "Q5_zero_cell_signed_target": (
                "epsilon*S_H=3" if exact_parallel_count == 2 else None
            ),
            "Q5_zero_cell_coefficient_offset": (
                3 if exact_parallel_count == 2 else None
            ),
            "Q5_required_coefficient_congruence": (
                "6 divides Q-3" if exact_parallel_count == 2 else None
            ),
            "Q5_zero_cell_coefficient_compatible": q5_compatible,
            "minimum_allowed_opposite_Q": minimum_allowed_q,
            "forced_mass14_parallel_count_Q": forced_q,
            "forced_mass14_mean": opposite_mean(forced_q),
            "directions_at_forced_Q_at_least": minimum_count,
            "forcing_reason": forcing_reason,
        }
        row["proved"] = bool(
            row["R_equals_P_plus_one"]
            and hard_edge_count + opposite_parallel_sum == H_EDGE_COUNT
            and opposite_mean(forced_q) == 14
            and minimum_count >= 1
            and (q5_compatible is not True)
        )
        _require(bool(row["proved"]), "a t=4,u=3 parallel ledger changed")
        candidate_rows[str(exact_parallel_count)] = row

    possible_p = [
        value
        for value in range(H_EDGE_COUNT + 1)
        if (value - coefficient_offset) % Q_MODULUS == 0
        and EXACT_HARD_COUNT * value + (value + 1) <= H_EDGE_COUNT
    ]
    proved = bool(
        type_budget == 154
        and quotient_sum == 8
        and phase_one_floors
        == {0: 26, 2: 12, 4: 26, 6: 24, 8: 26, 10: 20, 12: 14}
        and low_mean == EXACT_HARD_MEAN
        and low_rows
        == [
            {"b": 2, "floor": 12, "excess": 8, "status": "forbidden lift below 10"},
            {"b": 10, "floor": 20, "excess": 0, "status": "exact complement-triple cell"},
            {"b": 12, "floor": 14, "excess": 6, "status": "forbidden lift below 10"},
        ]
        and quotient_profile == [1, 1, 1, 1, 1, 1, 2]
        and all(target_checks)
        and coefficient_offset == 2
        and possible_p == [2, 8]
        and all_residues["surviving_residues_before_prop_15744"]
        == [0, 3, 4, 6]
        and all(bool(row["proved"]) for row in candidate_rows.values())
    )
    _require(proved, "the p=13,t=4,u=3 residue ledger failed")
    return {
        "p": P,
        "layer_index_t": LAYER_INDEX,
        "original_k": ORIGINAL_K,
        "H_edge_count": H_EDGE_COUNT,
        "hard_type": "epsilon_d=c_H",
        "hard_phase": 1,
        "opposite_phase": 0,
        "type_budget": type_budget,
        "type_budget_identity": "2*m*(m+t)=2*7*11=154",
        "hard_residue_u": HARD_RESIDUE_U,
        "all_residue_sieve_dependency": {
            "surviving_residues_before_prop_15744": all_residues[
                "surviving_residues_before_prop_15744"
            ],
            "proved": all_residues["proved"],
        },
        "hard_mean_formula": "a_L=6+14*k_L",
        "hard_quotient_sum": quotient_sum,
        "phase_one_even_b_floors": {
            str(b): value for b, value in phase_one_floors.items()
        },
        "low_mean_rows": low_rows,
        "sharp_nonzero_lift_floor": lift_floor,
        "hard_quotient_profile": quotient_profile,
        "exact_complement_triple_count": EXACT_HARD_COUNT,
        "elevated_hard_count": 1,
        "exact_hard_mean": EXACT_HARD_MEAN,
        "elevated_hard_mean": ELEVATED_HARD_MEAN,
        "exact_hard_target": (
            "epsilon*S_H=5-sum_(i in C)z_i+"
            "sum_({i,j} subset C)z_i*z_j"
        ),
        "exact_hard_target_boolean_checks": len(target_checks),
        "exact_hard_coefficient_offset": coefficient_offset,
        "exact_hard_parallel_congruence": "6 divides P-2",
        "possible_exact_hard_parallel_counts": possible_p,
        "parallel_ledgers": candidate_rows,
        "proved": proved,
    }


def _build_height_four_model(
    parallel_count: int,
) -> tuple[cp_model.CpModel, dict[str, object]]:
    """Build the changed-premise mass-14 height-four model at ``|H|=61``."""
    if parallel_count not in FORCED_OPPOSITE_PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 0 or 6")
    points = middle_slice_points()
    pairs = pair_coordinates()
    total_w = P * parallel_count - 53
    l1_budget = H_EDGE_COUNT - parallel_count

    model = cp_model.CpModel()
    values = [model.NewIntVar(0, 4, f"B_{index}") for index in range(len(points))]
    weights = {
        pair: model.NewIntVar(-l1_budget, l1_budget, f"W_{pair[0]}_{pair[1]}")
        for pair in pairs
    }
    absolute_weights = {
        pair: model.NewIntVar(0, l1_budget, f"absW_{pair[0]}_{pair[1]}")
        for pair in pairs
    }

    model.Add(sum(values) == MASS14_SUPPORT_SIZE)
    model.Add(values[0] == 4)
    model.Add(sum(weights.values()) == total_w)
    for pair in pairs:
        model.AddAbsEquality(absolute_weights[pair], weights[pair])
    model.Add(sum(absolute_weights.values()) <= l1_budget)
    row_halves = [
        model.NewIntVar(-l1_budget, l1_budget, f"row_half_{s}") for s in range(P)
    ]
    for s in range(P):
        model.Add(
            sum(value for pair, value in weights.items() if s in pair)
            == 2 * row_halves[s]
        )
    for index, point in enumerate(points):
        point_set = set(point)
        cut = sum(
            value
            for (s, t), value in weights.items()
            if (s in point_set) != (t in point_set)
        )
        model.Add(
            4 * values[index]
            == parallel_count - 3 + total_w - 2 * cut
        )
    return model, {
        "coefficient_sum": total_w,
        "l1_budget": l1_budget,
        "value_variable_count": len(values),
        "coefficient_variable_count": len(weights),
    }


@lru_cache(maxsize=None)
def h61_height_four_exclusion(parallel_count: int) -> dict[str, object]:
    """Exclude height four after rebuilding the model with ``l1<=61-Q``."""
    model, data = _build_height_four_model(parallel_count)
    validation = model.Validate()
    _require(not validation, f"invalid H=61 model: {validation}")
    model_hash = _model_textproto_sha256(model)
    _require(
        model_hash == EXPECTED_HEIGHT_FOUR_MODEL_SHA256[parallel_count],
        "the H=61 height-four model hash changed",
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = SEARCH_WORKERS
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    _require(infeasible, f"H=61,Q={parallel_count} model is {status_name}")
    return {
        "p": P,
        "H_edge_count": H_EDGE_COUNT,
        "parallel_count_Q": parallel_count,
        "scaled_mass_4pE_B": 14,
        "value_sum": MASS14_SUPPORT_SIZE,
        "height_four_orbit_anchor": "B(first lexicographic 7-set)=4",
        "height_four_orbit_anchor_is_wlog": True,
        "coefficient_sum": data["coefficient_sum"],
        "coefficient_sum_formula": "sum W=13Q-53",
        "cut_identity": "4B(X)=Q-3+sum(W)-2*cut_W(X)",
        "every_coefficient_row_sum_even": True,
        "transverse_edge_count": data["l1_budget"],
        "l1_budget": data["l1_budget"],
        "l1_inequality": f"sum |W_st|<={H_EDGE_COUNT}-{parallel_count}",
        "old_H59_l1_infeasibility_imported": False,
        "changed_premise_model_rebuilt": True,
        "value_variable_count": data["value_variable_count"],
        "coefficient_variable_count": data["coefficient_variable_count"],
        "integer_variable_count": len(model.Proto().variables),
        "model_constraint_count": len(model.Proto().constraints),
        "model_textproto_sha256": model_hash,
        "model_validation": validation,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "search_workers": SEARCH_WORKERS,
        "random_seed": 0,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "height_four_model_infeasible": infeasible,
        "result_status": "exhaustive finite certificate",
        "proved": infeasible,
    }


@lru_cache(maxsize=1)
def h61_mass14_cell_classification() -> dict[str, object]:
    """Classify the forced H=61 mass-14 cell at ``Q=0`` or ``Q=6``."""
    floors = {b: scaled_direction_floor(P, b, 0) for b in range(0, P, 2)}
    lift_floor = int(sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"])
    b2_offset = 4
    b2_compatible = {
        q: (q - b2_offset) % Q_MODULUS == 0
        for q in FORCED_OPPOSITE_PARALLEL_COUNTS
    }
    height_rows = {
        str(q): h61_height_four_exclusion(q)
        for q in FORCED_OPPOSITE_PARALLEL_COUNTS
    }

    height = mass14_height_dichotomy()
    classification = exact_mass14_boolean_classification()
    catalog = mass14_boolean_catalog_certificate()
    survivors = catalog["families_surviving_offset_mod_6"]
    moments = selected_pair_moment_certificate()
    proved = bool(
        floors == {0: 0, 2: 14, 4: 20, 6: 26, 8: 24, 10: 26, 12: 12}
        and lift_floor == 10
        and b2_compatible == {0: False, 6: False}
        and 0 < 14 - floors[12] < lift_floor
        and all(floors[b] > 14 for b in (4, 6, 8, 10))
        and all(bool(row["proved"]) for row in height_rows.values())
        and height["proved"]
        and height["height_dichotomy"] == [1, 4]
        and classification["catalog_exhaustive_at_support_462"]
        and survivors == {"0": ["selected_pair"], "6": ["selected_pair"]}
        and moments["proved"]
    )
    _require(proved, "the H=61 mass-14 cell classification failed")
    return {
        "p": P,
        "H_edge_count": H_EDGE_COUNT,
        "parallel_counts_Q": list(FORCED_OPPOSITE_PARALLEL_COUNTS),
        "scaled_mean": 14,
        "phase_zero_even_b_floors": {str(b): value for b, value in floors.items()},
        "b2_exact_target": "epsilon*S_H=4-z_i*z_j",
        "b2_coefficient_offset": b2_offset,
        "b2_compatible_at_Q": {
            str(q): value for q, value in b2_compatible.items()
        },
        "b12_floor_plus_two_excess": 14 - floors[12],
        "sharp_nonzero_lift_floor": lift_floor,
        "b12_floor_plus_two_excluded": True,
        "b4_b6_b8_b10_floors_above_mean": True,
        "remaining_cell_before_height_audit": "b=0, A=2B, 4p*E[B]=14",
        "height_dichotomy_dependency": {
            "proposition": "15.738",
            "live_certificate_called": True,
            "height_dichotomy": height["height_dichotomy"],
            "proved": height["proved"],
        },
        "H61_height_four_exclusions": height_rows,
        "boolean_catalog_edge_count_dependency": False,
        "boolean_catalog_result_status": classification["result_status"],
        "boolean_catalog_model_sha256": classification["model_textproto_sha256"],
        "boolean_catalog_exhaustive_at_support_462": classification[
            "catalog_exhaustive_at_support_462"
        ],
        "boolean_catalog_support_sha256": catalog["known_support_catalog_sha256"],
        "catalog_survivors_after_offset_mod_6": survivors,
        "forced_form": "B=x_i*x_j",
        "selected_pair_moments": {
            "degree_two": "(i-j)^2",
            "degree_four": "(i-j)^4",
            "Q0_and_Q6_differ_only_by_complete_graph_gauge": moments[
                "Q0_and_Q6_patterns_differ_by_complete_graph"
            ],
        },
        "result_status": "exhaustive finite certificate",
        "proved": proved,
    }


def _normalized_exceptional_hard_pattern(
    triple: tuple[int, int, int], parallel_count: int
) -> dict[tuple[int, int], int]:
    if parallel_count not in EXACT_HARD_PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 2 or 8")
    triple_set = set(triple)
    kernel_scalar = Fraction(parallel_count - 2, 12)
    pattern: dict[tuple[int, int], int] = {}
    for s, t in pair_coordinates():
        value = (
            int(s in triple_set and t in triple_set)
            + 2 * kernel_scalar
            - int(s in triple_set)
            - int(t in triple_set)
        )
        _require(value.denominator == 1, "hard coefficient is nonintegral")
        pattern[(s, t)] = int(value)
    return pattern


def _normalized_selected_pair_pattern(
    i: int, j: int, parallel_count: int
) -> dict[tuple[int, int], int]:
    if parallel_count not in FORCED_OPPOSITE_PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 0 or 6")
    kernel_scalar = Fraction(parallel_count - 6, 12)
    pattern: dict[tuple[int, int], int] = {}
    for s, t in pair_coordinates():
        value = (
            int((s, t) == (i, j))
            + 2 * kernel_scalar
            + int(s in (i, j))
            + int(t in (i, j))
        )
        _require(value.denominator == 1, "selected-pair coefficient is nonintegral")
        pattern[(s, t)] = int(value)
    return pattern


def _pattern_moment(pattern: dict[tuple[int, int], int], degree: int) -> int:
    return sum(
        coefficient * pow(s - t, degree, P)
        for (s, t), coefficient in pattern.items()
    ) % P


@lru_cache(maxsize=1)
def six_root_quartic_contradiction() -> dict[str, object]:
    """Check both signs and both gauges in the six-root quartic argument."""
    sign_rows: dict[str, dict[str, object]] = {}
    gauge_histograms = {
        str(q): {
            str(value): count
            for value, count in sorted(
                Counter(
                    _normalized_exceptional_hard_pattern((0, 1, 2), q).values()
                ).items()
            )
        }
        for q in EXACT_HARD_PARALLEL_COUNTS
    }
    for hard_sign in (-1, 1):
        hard_values = []
        for parallel_count in EXACT_HARD_PARALLEL_COUNTS:
            for triple in combinations(range(P), 3):
                pattern = _normalized_exceptional_hard_pattern(triple, parallel_count)
                s2 = _pattern_moment(pattern, 2)
                s4 = _pattern_moment(pattern, 4)
                m2 = hard_sign * s2
                m4 = hard_sign * s4
                hard_values.append((2 * hard_sign * m4 - m2 * m2) % P)

        opposite_sign = -hard_sign
        opposite_values = []
        for parallel_count in FORCED_OPPOSITE_PARALLEL_COUNTS:
            for i, j in pair_coordinates():
                pattern = _normalized_selected_pair_pattern(i, j, parallel_count)
                s2 = _pattern_moment(pattern, 2)
                s4 = _pattern_moment(pattern, 4)
                m2 = opposite_sign * s2
                m4 = opposite_sign * s4
                opposite_values.append((2 * hard_sign * m4 - m2 * m2) % P)
        expected = {
            (-3 * pow(i - j, 4, P)) % P for i, j in pair_coordinates()
        }
        row = {
            "hard_sign_h": hard_sign,
            "opposite_sign": opposite_sign,
            "hard_G_value_set": sorted(set(hard_values)),
            "opposite_G_value_set": sorted(set(opposite_values)),
            "expected_opposite_nonzero_value_set": sorted(expected),
            "every_opposite_value_nonzero": all(opposite_values),
        }
        row["proved"] = bool(
            set(hard_values) == {0}
            and set(opposite_values) == expected
            and all(opposite_values)
        )
        _require(bool(row["proved"]), "a sign-safe quartic check failed")
        sign_rows[str(hard_sign)] = row

    root_count = EXACT_HARD_COUNT
    quartic_degree = 4
    proved = bool(
        gauge_histograms == {
            "2": {"-1": 33, "0": 45},
            "8": {"0": 33, "1": 45},
        }
        and root_count > quartic_degree
        and all(bool(row["proved"]) for row in sign_rows.values())
    )
    _require(proved, "the six-root quartic contradiction failed")
    return {
        "field": "F_13",
        "global_even_moments": (
            "M_d(L)=sum_({u,v} in H)chi(u-v)*(L(u)-L(v))^d"
        ),
        "homogeneous_quartic": "G=2*h*M_4-M_2^2",
        "quartic_degree": quartic_degree,
        "distinct_exact_hard_projective_roots": root_count,
        "root_count_exceeds_degree": root_count > quartic_degree,
        "hard_gauge_histograms": gauge_histograms,
        "both_hard_gauges_have_triangle_even_moments": True,
        "sign_checks": sign_rows,
        "opposite_evaluation_formula": "G=-3*(i-j)^4",
        "minus_three_nonzero_mod_13": (-3) % P != 0,
        "selected_pair_contradicts_identically_zero_G": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15744() -> dict[str, object]:
    """Package the exact exclusion of ``p=13,t=4,u=3``."""
    all_residues = t4_all_residue_sieve()
    ledger = t4_u3_residue_ledger()
    cells = h61_mass14_cell_classification()
    quartic = six_root_quartic_contradiction()
    remaining_residues = [
        u
        for u in all_residues["surviving_residues_before_prop_15744"]
        if u != HARD_RESIDUE_U
    ]
    proved = bool(
        all_residues["proved"]
        and ledger["proved"]
        and cells["proved"]
        and quartic["proved"]
        and remaining_residues == [0, 4, 6]
    )
    _require(proved, "Proposition 15.744 certificate failed")
    return {
        "prop": "15.744",
        "title": "Six-root quartic closes the p=13,t=4,u=3 branch",
        "result_status": "proved branch theorem",
        "statement": (
            "the residual-(ii) branch p=13,t=4,k=60,u=3 is empty"
        ),
        "changed_premise": (
            "the b=10 equality uses a full-rank contact-layer restriction, "
            "its floor-plus-two case uses a separate punctured-lift model, "
            "and the mass-14 height-four cells are rebuilt at |H|=61"
        ),
        "all_t4_residue_sieve": all_residues,
        "residue_and_parallel_ledger": ledger,
        "H61_mass14_cell": cells,
        "six_root_quartic": quartic,
        "p13_t4_u3_branch_closed": proved,
        "p13_k_eq_60_closed": False,
        "remaining_p13_t4_residues": remaining_residues,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "finite_configuration_search_used": True,
        "finite_search_scope": (
            "one 1716-variable b=10 punctured-lift model and two necessary "
            "78-coefficient height-four cells at Q=0,6; the Boolean catalog "
            "is imported only after the H=61 reruns"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic Proposition 15.744 evidence artifact."""
    target = ROOT / "evidence" / "e1_gmin_m4_prop15744.json"
    target.write_text(json.dumps(proposition_15744(), indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    result = proposition_15744()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.744 failed")
    path = write_evidence()
    print("Prop 15.744 p=13,t=4,u=3: excluded")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
