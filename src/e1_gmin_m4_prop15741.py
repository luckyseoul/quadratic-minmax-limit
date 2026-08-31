#!/usr/bin/env python3
r"""Prop. 15.741 -- four-star mixed moments and the common-graph frontier.

Proposition 15.740 leaves one generic ``p=13,t=3`` hard quotient pattern,
namely ``1^4 2^3``.  A hypothetical realization has four exact positive
stars, three elevated hard cells, and seven opposite cells, all induced by
one 59-edge graph.

Besides the even difference moments, introduce the symmetric endpoint forms

    T_3(L) = sum_e chi(delta_e) (L(u)+L(v)) (L(u)-L(v))^2,
    U_4(L) = sum_e chi(delta_e) (L(u)+L(v))^2 (L(u)-L(v))^2.

An exact star kills ``M_2,T_3,M_4,U_4``.  Four distinct exact directions
therefore force ``M_2=T_3=0`` globally and make ``M_4,U_4`` proportional.
Proposition 15.740 rules out ``M_4=0``, so there is a unique affine-label
invariant ``lambda`` with ``U_4=lambda*M_4`` and ``M_4`` is nonzero in all
ten nonexact directions.

This is a proved open reduction, not a branch close.  Explicit elevated and
opposite one-direction coefficient cells satisfy all their cut, parity,
sum, and l1 conditions with the same ``lambda=7``.  They do not construct a
common graph.  The module also records the exact signed-collision identity
forced by each star and a coefficientwise midpoint/displacement
parameterization of the common graph.
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
from typing import Iterable, Mapping

from e1_gmin_m4_prop15740 import (
    MOMENT_CANDIDATE_COUNT,
    translated_cut_nine_vector_certificate,
    translated_cut_vector,
)


ROOT = Path(__file__).resolve().parents[1]
P = 13
M = 7
H_EDGE_COUNT = 59
HARD_EDGE_COUNT = 38
OPPOSITE_EDGE_COUNT = 21
EXACT_HARD_DIRECTION_COUNT = 4
ELEVATED_HARD_DIRECTION_COUNT = 3
OPPOSITE_DIRECTION_COUNT = 7
PROJECTIVE_DIRECTION_COUNT = P + 1
DISTANCES = tuple(range(1, (P + 1) // 2))

Pair = tuple[int, int]
Pattern = dict[Pair, int]
Point = tuple[int, int]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _pair(s: int, t: int) -> Pair:
    if s == t:
        raise ValueError("a coefficient pair must have distinct labels")
    return (s, t) if s < t else (t, s)


def _add_coefficient(pattern: Pattern, s: int, t: int, value: int) -> None:
    pair = _pair(s, t)
    pattern[pair] = pattern.get(pair, 0) + value
    if pattern[pair] == 0:
        del pattern[pair]


def endpoint_moments(pattern: Mapping[Pair, int]) -> dict[str, int]:
    """Return ``M2,T3,M4,U4`` for one normalized fibre matrix."""
    totals = {"M2": 0, "T3": 0, "M4": 0, "U4": 0}
    for (s, t), coefficient in pattern.items():
        delta = (s - t) % P
        sigma = (s + t) % P
        delta2 = delta * delta % P
        totals["M2"] += coefficient * delta2
        totals["T3"] += coefficient * sigma * delta2
        totals["M4"] += coefficient * delta2 * delta2
        totals["U4"] += coefficient * sigma * sigma * delta2
    return {name: value % P for name, value in totals.items()}


def row_sums(pattern: Mapping[Pair, int]) -> dict[int, int]:
    return {
        vertex: sum(
            coefficient
            for pair, coefficient in pattern.items()
            if vertex in pair
        )
        for vertex in range(P)
    }


def cut_value(pattern: Mapping[Pair, int], subset: Iterable[int]) -> int:
    chosen = set(subset)
    return sum(
        coefficient
        for (s, t), coefficient in pattern.items()
        if (s in chosen) != (t in chosen)
    )


def distance_aggregates(pattern: Mapping[Pair, int]) -> tuple[int, ...]:
    """Aggregate the thirteen pairs of each cyclic distance ``1..6``."""
    return tuple(
        sum(
            pattern.get(_pair(label, (label + distance) % P), 0)
            for label in range(P)
        )
        for distance in DISTANCES
    )


def exact_positive_star(center: int) -> Pattern:
    if not 0 <= center < P:
        raise ValueError("the star center must lie in F_13")
    return {
        _pair(center, label): 1
        for label in range(P)
        if label != center
    }


@lru_cache(maxsize=1)
def exact_star_moment_certificate() -> dict[str, object]:
    """Check all four endpoint contractions on every translated star."""
    checks = {
        center: endpoint_moments(exact_positive_star(center))
        for center in range(P)
    }
    power_sums = {
        degree: sum(pow(value, degree, P) for value in range(1, P)) % P
        for degree in (2, 3, 4)
    }
    distance_rows = {
        center: distance_aggregates(exact_positive_star(center))
        for center in range(P)
    }
    proved = bool(
        power_sums == {2: 0, 3: 0, 4: 0}
        and all(
            moments == {"M2": 0, "T3": 0, "M4": 0, "U4": 0}
            for moments in checks.values()
        )
        and all(row == (2, 2, 2, 2, 2, 2) for row in distance_rows.values())
    )
    _require(proved, "the exact-star endpoint moments changed")
    return {
        "p": P,
        "exact_star_count_checked": P,
        "power_sums_mod_13": power_sums,
        "star_moments_by_center": checks,
        "star_distance_aggregate": list(distance_rows[0]),
        "star_distance_energy": sum(value * value for value in distance_rows[0]),
        "all_M2_T3_M4_U4_zero": True,
        "proved": proved,
    }


def _mod_rank(rows: Iterable[Iterable[int]], modulus: int = P) -> int:
    matrix = [[value % modulus for value in row] for row in rows]
    if not matrix:
        return 0
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, modulus)
        matrix[rank] = [(inverse * value) % modulus for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % modulus
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def _homogeneous_evaluation(
    coefficients: Iterable[int], point: Point
) -> int:
    values = tuple(coefficients)
    degree = len(values) - 1
    x, y = point
    return sum(
        coefficient * pow(x, degree - index, P) * pow(y, index, P)
        for index, coefficient in enumerate(values)
    ) % P


def _polynomial_product(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % P
    return tuple(result)


def _root_factor(point: Point) -> tuple[int, int]:
    x, y = point
    return y % P, (-x) % P


def projective_points() -> tuple[Point, ...]:
    return tuple((1, slope) for slope in range(P)) + ((0, 1),)


@lru_cache(maxsize=1)
def quartic_root_rank_certificate() -> dict[str, object]:
    """Verify the degree 2/3/4 root ranks for every four-point set."""
    points = projective_points()
    rows_for_digest: list[str] = []
    quartet_count = 0
    for roots in combinations(range(len(points)), EXACT_HARD_DIRECTION_COUNT):
        selected = [points[index] for index in roots]
        ranks = {}
        for degree in (2, 3, 4):
            evaluation = [
                [
                    pow(x, degree - monomial, P) * pow(y, monomial, P) % P
                    for monomial in range(degree + 1)
                ]
                for x, y in selected
            ]
            ranks[degree] = _mod_rank(evaluation)
        root_product = (1,)
        for point in selected:
            root_product = _polynomial_product(root_product, _root_factor(point))
        product_vanishes = all(
            _homogeneous_evaluation(root_product, point) == 0
            for point in selected
        )
        _require(
            ranks == {2: 3, 3: 4, 4: 4}
            and len(root_product) == 5
            and any(root_product)
            and product_vanishes,
            "four-point binary-form root rank failed",
        )
        rows_for_digest.append(
            ",".join(map(str, (*roots, *root_product)))
        )
        quartet_count += 1
    digest = hashlib.sha256(";".join(rows_for_digest).encode("ascii")).hexdigest()
    proved = quartet_count == comb(PROJECTIVE_DIRECTION_COUNT, 4)
    _require(proved, "the projective quartet count changed")
    return {
        "field": "F_13",
        "projective_direction_count": PROJECTIVE_DIRECTION_COUNT,
        "four_direction_sets_checked": quartet_count,
        "degree_2_evaluation_rank": 3,
        "degree_3_evaluation_rank": 4,
        "degree_4_evaluation_rank": 4,
        "degree_2_four_roots_force_zero": True,
        "degree_3_four_roots_force_zero": True,
        "degree_4_four_root_kernel_dimension": 1,
        "degree_4_kernel_generator": "product of the four root linear factors",
        "quartet_product_sha256": digest,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def endpoint_contraction_basis_certificate() -> dict[str, object]:
    """Classify symmetric diagonal-zero endpoint forms through degree four."""
    bases = {
        2: ((1, -2, 1),),
        3: ((1, -1, -1, 1),),
        4: ((1, -4, 6, -4, 1), (1, 0, -2, 0, 1)),
    }
    expected_dimensions = {2: 1, 3: 1, 4: 2}
    checks: dict[int, dict[str, object]] = {}
    for degree, rows in bases.items():
        symmetric = all(
            tuple(value % P for value in row)
            == tuple(value % P for value in reversed(row))
            for row in rows
        )
        diagonal_zero = all(sum(row) % P == 0 for row in rows)
        rank = _mod_rank(rows)
        dimension = degree // 2 + 1 - 1
        checks[degree] = {
            "basis_coefficients": [list(row) for row in rows],
            "symmetric": symmetric,
            "diagonal_zero": diagonal_zero,
            "basis_rank": rank,
            "space_dimension": dimension,
        }
        _require(
            symmetric
            and diagonal_zero
            and rank == dimension == expected_dimensions[degree],
            f"degree-{degree} endpoint basis failed",
        )
    return {
        "classification": (
            "a symmetric homogeneous f(s,t) with f(s,s)=0 is divisible "
            "by (s-t)^2"
        ),
        "degrees": checks,
        "named_bases": {
            "degree_2": ["(s-t)^2"],
            "degree_3": ["(s+t)(s-t)^2"],
            "degree_4": ["(s-t)^4", "(s+t)^2(s-t)^2"],
        },
        "orientation_independent_endpoint_contractions_through_degree_4_exhausted": True,
        "proved": True,
    }


@lru_cache(maxsize=1)
def affine_label_invariance_certificate() -> dict[str, object]:
    """Check the four moment transformation laws coefficientwise."""
    checks = 0
    for scale in range(1, P):
        for shift in range(P):
            for s, t in combinations(range(P), 2):
                delta = (s - t) % P
                sigma = (s + t) % P
                new_s = (scale * s + shift) % P
                new_t = (scale * t + shift) % P
                new_delta = (new_s - new_t) % P
                new_sigma = (new_s + new_t) % P
                m2 = delta * delta % P
                t3 = sigma * m2 % P
                m4 = m2 * m2 % P
                u4 = sigma * sigma * m2 % P
                _require(new_delta * new_delta % P == scale * scale * m2 % P, "M2 transform failed")
                _require(
                    new_sigma * new_delta * new_delta % P
                    == (
                        pow(scale, 3, P) * t3
                        + 2 * shift * scale * scale * m2
                    )
                    % P,
                    "T3 transform failed",
                )
                _require(
                    pow(new_delta, 4, P) == pow(scale, 4, P) * m4 % P,
                    "M4 transform failed",
                )
                _require(
                    new_sigma * new_sigma * new_delta * new_delta % P
                    == (
                        pow(scale, 4, P) * u4
                        + 4 * shift * pow(scale, 3, P) * t3
                        + 4 * shift * shift * scale * scale * m2
                    )
                    % P,
                    "U4 transform failed",
                )
                checks += 1
    return {
        "affine_relabeling": "s -> a*s+c with a nonzero",
        "M2_transform": "M2' = a^2*M2",
        "T3_transform": "T3' = a^3*T3+2*c*a^2*M2",
        "M4_transform": "M4' = a^4*M4",
        "U4_transform": "U4' = a^4*U4+4*c*a^3*T3+4*c^2*a^2*M2",
        "coefficientwise_checks": checks,
        "lambda_U4_over_M4_invariant_when_M2_T3_zero": True,
        "proved": checks == (P - 1) * P * comb(P, 2),
    }


def elevated_lambda_seven_pattern() -> Pattern:
    support = (0, 1, 4, 9, 12)
    pattern: Pattern = {}
    for s, t in combinations(support, 2):
        _add_coefficient(pattern, s, t, 1)
    _add_coefficient(pattern, 4, 9, 1)
    return pattern


@lru_cache(maxsize=1)
def elevated_lambda_seven_certificate() -> dict[str, object]:
    """Validate the elevated hard one-direction coefficient cell."""
    pattern = elevated_lambda_seven_pattern()
    rows = row_sums(pattern)
    odd_rows = [vertex for vertex, value in rows.items() if value % 2]
    cuts: Counter[int] = Counter()
    values: Counter[int] = Counter()
    for subset in combinations(range(P), M):
        cut = cut_value(pattern, subset)
        cuts[cut] += 1
        values[7 - cut] += 1
    moments = endpoint_moments(pattern)
    aggregates = distance_aggregates(pattern)
    value_sum = sum(value * count for value, count in values.items())
    scaled_mean = Fraction(2 * P * value_sum, comb(P, M))
    coefficient_sum = sum(pattern.values())
    l1_norm = sum(abs(value) for value in pattern.values())
    cancellation_pairs = (H_EDGE_COUNT - 6 - l1_norm) // 2
    common_lambda = moments["U4"] * pow(moments["M4"], -1, P) % P
    proved = bool(
        pattern[(4, 9)] == 2
        and coefficient_sum == 11
        and l1_norm == 11
        and odd_rows == [4, 9]
        and max(cuts) == 7
        and min(values) == 0
        and scaled_mean == 28
        and moments == {"M2": 0, "T3": 0, "M4": 7, "U4": 10}
        and aggregates == (2, 1, 2, 2, 4, 0)
        and cancellation_pairs == 21
        and common_lambda == 7
    )
    _require(proved, "the elevated lambda-seven cell changed")
    return {
        "cell_type": "elevated hard one-direction coefficient relaxation",
        "hard_parallel_count_P": 6,
        "support_description": (
            "unit K5 on {0,1,4,9,12}, plus one additional unit on {4,9}"
        ),
        "nonzero_coefficients": {
            f"{s},{t}": value for (s, t), value in sorted(pattern.items())
        },
        "coefficient_sum": coefficient_sum,
        "l1_norm": l1_norm,
        "available_nonparallel_edge_count": H_EDGE_COUNT - 6,
        "sign_cancelling_padding_pairs": cancellation_pairs,
        "row_sums": rows,
        "odd_rows": odd_rows,
        "directional_b": len(odd_rows),
        "cut_histogram": dict(sorted(cuts.items())),
        "A_formula": "A(X)=7-cut_W(X)",
        "A_value_histogram": dict(sorted(values.items())),
        "all_1716_middle_cuts_nonnegative": min(values) >= 0,
        "scaled_mean_2pE_A": int(scaled_mean),
        "moments_mod_13": moments,
        "lambda_U4_over_M4": common_lambda,
        "distance_aggregates": list(aggregates),
        "distance_energy": sum(value * value for value in aggregates),
        "constructs_common_59_edge_graph": False,
        "proved": proved,
    }


def opposite_lambda_seven_pattern() -> Pattern:
    pattern: Pattern = {}
    for label in range(1, P):
        _add_coefficient(pattern, 0, label, -1)
    for label in sorted(set(range(P)) - {0, 1, 3, 6}):
        _add_coefficient(pattern, 1, label, -1)
    _add_coefficient(pattern, 3, 6, 1)
    return pattern


@lru_cache(maxsize=1)
def opposite_lambda_seven_certificate() -> dict[str, object]:
    """Validate the opposite one-direction coefficient cell."""
    pattern = opposite_lambda_seven_pattern()
    rows = row_sums(pattern)
    cuts: Counter[int] = Counter()
    values: Counter[int] = Counter()
    analytic_formula_checks = 0
    ordinary = set(range(P)) - {0, 1, 3, 6}
    for subset_tuple in combinations(range(P), M):
        subset = set(subset_tuple)
        cut_w = cut_value(pattern, subset)
        cut_g = -cut_w
        indicator_0 = int(0 in subset)
        indicator_1 = int(1 in subset)
        indicator_3 = int(3 in subset)
        indicator_6 = int(6 in subset)
        ordinary_count = len(subset & ordinary)
        _require(
            indicator_0
            + indicator_1
            + indicator_3
            + indicator_6
            + ordinary_count
            == M,
            "opposite cell partition failed",
        )
        xor_36 = indicator_3 ^ indicator_6
        if indicator_1 == 0:
            formula = 14 - 2 * indicator_0 - indicator_3 - indicator_6 - xor_36
        else:
            formula = 10 + indicator_3 + indicator_6 - xor_36
        _require(cut_g == formula, "opposite analytic cut formula failed")
        cuts[cut_w] += 1
        value = (cut_g - 10) // 2
        _require(2 * value == cut_g - 10, "opposite cell lost integrality")
        values[value] += 1
        analytic_formula_checks += 1
    moments = endpoint_moments(pattern)
    aggregates = distance_aggregates(pattern)
    value_sum = sum(value * count for value, count in values.items())
    scaled_mean = Fraction(4 * P * value_sum, comb(P, M))
    coefficient_sum = sum(pattern.values())
    l1_norm = sum(abs(value) for value in pattern.values())
    cancellation_pairs = (H_EDGE_COUNT - 3 - l1_norm) // 2
    common_lambda = moments["U4"] * pow(moments["M4"], -1, P) % P
    proved = bool(
        coefficient_sum == -20
        and l1_norm == 22
        and all(value % 2 == 0 for value in rows.values())
        and max(cuts) == -10
        and min(values) == 0
        and scaled_mean == 20
        and moments == {"M2": 0, "T3": 0, "M4": 8, "U4": 4}
        and aggregates == (-3, -3, -3, -4, -3, -4)
        and cancellation_pairs == 17
        and common_lambda == 7
        and analytic_formula_checks == comb(P, M)
    )
    _require(proved, "the opposite lambda-seven cell changed")
    return {
        "cell_type": "opposite one-direction coefficient relaxation",
        "opposite_parallel_count_Q": 3,
        "support_description": (
            "w_0x=-1 for x!=0; w_1x=-1 outside {0,1,3,6}; w_36=+1"
        ),
        "nonzero_coefficients": {
            f"{s},{t}": value for (s, t), value in sorted(pattern.items())
        },
        "coefficient_sum": coefficient_sum,
        "l1_norm": l1_norm,
        "available_nonparallel_edge_count": H_EDGE_COUNT - 3,
        "sign_cancelling_padding_pairs": cancellation_pairs,
        "row_sums": rows,
        "every_row_sum_even": True,
        "analytic_cut_formula": {
            "indicator_1_zero": "cut_(-W)=14-2*A-R-S-(R xor S)",
            "indicator_1_one": "cut_(-W)=10+R+S-(R xor S)",
        },
        "analytic_cut_formula_checks": analytic_formula_checks,
        "cut_W_histogram": dict(sorted(cuts.items())),
        "B_formula": "B(X)=-5-cut_W(X)/2=(cut_(-W)(X)-10)/2",
        "B_value_histogram": dict(sorted(values.items())),
        "all_1716_middle_cuts_nonnegative": min(values) >= 0,
        "scaled_mean_4pE_B": int(scaled_mean),
        "moments_mod_13": moments,
        "lambda_U4_over_M4": common_lambda,
        "distance_aggregates": list(aggregates),
        "distance_energy": sum(value * value for value in aggregates),
        "constructs_common_59_edge_graph": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def previous_elevated_witness_cubic_obstruction() -> dict[str, object]:
    """Show that the older M2-only local witness fails the new cubic."""
    support = (0, 1, 2, 3, 5)
    pattern: Pattern = {}
    for s, t in combinations(support, 2):
        _add_coefficient(pattern, s, t, 1)
    _add_coefficient(pattern, 0, 11, 1)
    moments = endpoint_moments(pattern)
    proved = moments == {"M2": 0, "T3": 4, "M4": 5, "U4": 4}
    _require(proved, "the previous elevated witness moments changed")
    return {
        "previous_witness": "K5 on {0,1,2,3,5}, plus {0,11}",
        "moments_mod_13": moments,
        "excluded_by_global_T3_zero": moments["T3"] != 0,
        "result_status": "counterexample retired by stronger invariant",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def exact_star_collision_certificate() -> dict[str, object]:
    """Derive the signed edge-pair collision sum in an exact projection."""
    algebra_checks = 0
    for positive in range(55):
        for negative in range(55):
            left = comb(positive, 2) + comb(negative, 2) - positive * negative
            numerator = (positive - negative) ** 2 - positive - negative
            _require(numerator % 2 == 0 and left == numerator // 2, "collision identity failed")
            algebra_checks += 1
    positive_nonparallel = HARD_EDGE_COUNT - 5
    negative_nonparallel = OPPOSITE_EDGE_COUNT
    nonparallel_total = positive_nonparallel + negative_nonparallel
    star_support_cells = P - 1
    one_direction_collision = (star_support_cells - nonparallel_total) // 2
    four_direction_collision = EXACT_HARD_DIRECTION_COUNT * one_direction_collision
    proved = bool(
        HARD_EDGE_COUNT == 4 * 5 + 3 * 6
        and OPPOSITE_EDGE_COUNT == 7 * 3
        and positive_nonparallel == 33
        and negative_nonparallel == 21
        and nonparallel_total == 54
        and one_direction_collision == -21
        and four_direction_collision == -84
    )
    _require(proved, "the exact-star collision ledger changed")
    return {
        "hard_type_edge_count": HARD_EDGE_COUNT,
        "opposite_type_edge_count": OPPOSITE_EDGE_COUNT,
        "exact_direction_parallel_hard_edges": 5,
        "positive_nonparallel_edges": positive_nonparallel,
        "negative_nonparallel_edges": negative_nonparallel,
        "cell_equation": (
            "a_c-b_c=1 on the 12 star cells and 0 on every other fibre-pair cell"
        ),
        "cell_collision_identity": (
            "C(a_c,2)+C(b_c,2)-a_c*b_c=((a_c-b_c)^2-a_c-b_c)/2"
        ),
        "integer_algebra_checks": algebra_checks,
        "signed_collision_sum_per_exact_direction": one_direction_collision,
        "signed_collision_sum_over_four_exact_directions": four_direction_collision,
        "line_signed_degree_identity": (
            "sum_(v:L_i(v)=s)(deg_+(v)-deg_-(v))="
            "1+11*1_(s=j_i)+2*p_(i,s), with sum_s p_(i,s)=5"
        ),
        "aggregate_identity_alone_couples_the_four_directions": False,
        "common_graph_leverage_requires_tracking_same_edge_pairs": True,
        "proved": proved,
    }


def _point_add(left: Point, right: Point) -> Point:
    return (left[0] + right[0]) % P, (left[1] + right[1]) % P


def _point_sub(left: Point, right: Point) -> Point:
    return (left[0] - right[0]) % P, (left[1] - right[1]) % P


def _point_scale(scalar: int, point: Point) -> Point:
    return scalar * point[0] % P, scalar * point[1] % P


def _linear_value(functional: Point, point: Point) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % P


def _normalize_displacement(displacement: Point) -> Point:
    if displacement == (0, 0):
        raise ValueError("a displacement class cannot be zero")
    negative = ((-displacement[0]) % P, (-displacement[1]) % P)
    return min(displacement, negative)


def _projected_distance(value: int) -> int:
    value %= P
    return min(value, (-value) % P)


@lru_cache(maxsize=1)
def midpoint_displacement_certificate() -> dict[str, object]:
    """Certify the binary midpoint/displacement coordinates coefficientwise."""
    inverse_two = pow(2, -1, P)
    points = tuple(product(range(P), repeat=2))
    displacement_classes = tuple(
        sorted(
            {
                _normalize_displacement(point)
                for point in points
                if point != (0, 0)
            }
        )
    )
    edge_coordinates: dict[tuple[Point, Point], tuple[Point, Point]] = {}
    for u, v in combinations(points, 2):
        midpoint = _point_scale(inverse_two, _point_add(u, v))
        displacement = _normalize_displacement(_point_sub(u, v))
        key = midpoint, displacement
        _require(key not in edge_coordinates, "midpoint/displacement coordinates collided")
        half_displacement = _point_scale(inverse_two, displacement)
        recovered = tuple(
            sorted(
                (
                    _point_add(midpoint, half_displacement),
                    _point_sub(midpoint, half_displacement),
                )
            )
        )
        _require(recovered == tuple(sorted((u, v))), "midpoint inverse failed")
        edge_coordinates[key] = (u, v)

    coefficientwise_checks = 0
    for (midpoint, displacement), (u, v) in edge_coordinates.items():
        for functional in projective_points():
            lu = _linear_value(functional, u)
            lv = _linear_value(functional, v)
            lm = _linear_value(functional, midpoint)
            ld = _linear_value(functional, displacement)
            _require((lu + lv) % P == 2 * lm % P, "midpoint linear sum failed")
            _require((lu - lv) ** 2 % P == ld * ld % P, "displacement square failed")
            _require((lu - lv) ** 4 % P == pow(ld, 4, P), "displacement quartic failed")
            coefficientwise_checks += 1

    bucket_rows = {}
    for index, functional in enumerate(projective_points()):
        counts = Counter(
            _projected_distance(_linear_value(functional, displacement))
            for displacement in displacement_classes
        )
        expected = {0: 6, **{distance: 13 for distance in DISTANCES}}
        _require(dict(sorted(counts.items())) == expected, "difference bucket sizes failed")
        bucket_rows[index] = dict(sorted(counts.items()))

    edge_variable_count = len(edge_coordinates)
    direction_parallel_counts = [5] * 4 + [6] * 3 + [3] * 7
    proved = bool(
        len(points) == P * P
        and len(displacement_classes) == (P * P - 1) // 2
        and edge_variable_count == comb(P * P, 2)
        and coefficientwise_checks == edge_variable_count * (P + 1)
        and sum(direction_parallel_counts) == H_EDGE_COUNT
    )
    _require(proved, "the midpoint/displacement certificate failed")
    return {
        "affine_point_count": len(points),
        "nonzero_displacement_classes_modulo_sign": len(displacement_classes),
        "binary_edge_variables_n_m_delta": edge_variable_count,
        "edge_parameterization": "{u,v}={m+delta/2,m-delta/2}, delta modulo +/-",
        "parameterization_is_bijective": True,
        "coefficientwise_functional_checks": coefficientwise_checks,
        "difference_bucket_sizes_for_every_projective_L": bucket_rows[0],
        "difference_aggregate_definition": (
            "m_delta=sum_m n_(m,delta); q_L(a)=epsilon_L*"
            "sum_(|L(delta)|=a) chi(delta)*m_delta"
        ),
        "certified_contractions": {
            "M2": "sum_(m,delta) chi(delta)*n_(m,delta)*L(delta)^2",
            "T3": "2*sum_(m,delta) chi(delta)*n_(m,delta)*L(m)*L(delta)^2",
            "M4": "sum_(m,delta) chi(delta)*n_(m,delta)*L(delta)^4",
            "U4": "4*sum_(m,delta) chi(delta)*n_(m,delta)*L(m)^2*L(delta)^2",
        },
        "zeroth_midpoint_moment_seen_by_difference_aggregates": "m_delta",
        "new_first_midpoint_moment_seen_by_T3": True,
        "new_second_midpoint_moment_seen_by_U4": True,
        "direction_parallel_count_multiset": direction_parallel_counts,
        "direction_parallel_counts_sum": sum(direction_parallel_counts),
        "midpoint_subcertificate_asserts_difference_Gram_or_inverse": False,
        "constructs_common_59_edge_graph": False,
        "proved": proved,
    }


def _projective_direction(vector: Point) -> Point:
    """Return the normalized projective direction of a nonzero vector."""
    x, y = vector
    if x % P:
        inverse = pow(x % P, -1, P)
        return 1, y * inverse % P
    if y % P:
        return 0, 1
    raise ValueError("zero has no projective direction")


@lru_cache(maxsize=1)
def difference_radon_gram_certificate() -> dict[str, object]:
    """Certify the 84-column difference-Radon Gram and exact inverse.

    Rows are ``(L,a)`` with fourteen projective functionals and
    ``a in F_13/+/-={0,...,6}``; columns are nonzero displacement classes
    modulo sign.  The zero bin is essential in the stated Gram matrix.
    """
    points = tuple(product(range(P), repeat=2))
    directions = projective_points()
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
    column_directions = [
        direction_index[_projective_direction(displacement)]
        for displacement in displacements
    ]
    rows = tuple(
        (functional, distance)
        for functional in directions
        for distance in range((P + 1) // 2)
    )
    incidence = [
        [
            int(
                _projected_distance(_linear_value(functional, displacement))
                == distance
            )
            for displacement in displacements
        ]
        for functional, distance in rows
    ]
    row_sums = [sum(row) for row in incidence]
    gram = [
        [
            sum(incidence[row][left] * incidence[row][right] for row in range(len(rows)))
            for right in range(len(displacements))
        ]
        for left in range(len(displacements))
    ]
    gram_checks = 0
    for left in range(len(displacements)):
        for right in range(len(displacements)):
            same_direction = column_directions[left] == column_directions[right]
            expected = (
                14
                if left == right
                else 1
                if same_direction
                else 2
            )
            _require(gram[left][right] == expected, "difference-Radon Gram failed")
            gram_checks += 1

    # Check the displayed inverse as a coefficient identity.  Choose one
    # hard/opposite sign split; the cancellation is independent of that split.
    direction_signs = [1] * 7 + [-1] * 7
    inverse_checks = 0
    for left, left_direction in enumerate(column_directions):
        left_sign = direction_signs[left_direction]
        for right, right_direction in enumerate(column_directions):
            right_sign = direction_signs[right_direction]
            same_direction = int(left_direction == right_direction)
            coefficient = (
                same_direction
                - 2 * left_sign * right_sign
                + left_sign * right_sign * gram[left][right]
            )
            _require(
                coefficient == (P if left == right else 0),
                "difference-Radon inverse coefficient failed",
            )
            inverse_checks += 1

    parallel_counts = [5] * 4 + [6] * 3 + [3] * 7
    signed_total = HARD_EDGE_COUNT - OPPOSITE_EDGE_COUNT
    parallel_square_sum = sum(value * value for value in parallel_counts)
    fractional_q_values = {
        "exact_hard_P5": Fraction(signed_total - 5, 6),
        "elevated_hard_P6": Fraction(signed_total - 6, 6),
        "opposite_P3": Fraction(-signed_total - 3, 6),
    }
    interval_cut_coefficient_sum = 42
    fractional_translated_cut_values = {
        name: interval_cut_coefficient_sum * value
        for name, value in fractional_q_values.items()
    }
    fractional_point_checked = bool(
        fractional_q_values
        == {
            "exact_hard_P5": Fraction(2),
            "elevated_hard_P6": Fraction(11, 6),
            "opposite_P3": Fraction(-10, 3),
        }
        and fractional_translated_cut_values["elevated_hard_P6"] < 91
        and fractional_translated_cut_values["opposite_P3"] < -130
    )
    constant_off_bin_energy = (
        2 * signed_total * signed_total - 2 * parallel_square_sum
    )
    base_displacement_square_sum = H_EDGE_COUNT
    total_off_bin_energy_at_no_collisions = (
        P * base_displacement_square_sum + constant_off_bin_energy
    )
    exact_star_energy = 6 * 2 * 2
    four_star_energy = EXACT_HARD_DIRECTION_COUNT * exact_star_energy
    nonexact_base_energy = total_off_bin_energy_at_no_collisions - four_star_energy
    matrix_payload = ";".join(
        "".join(map(str, row)) for row in incidence
    ).encode("ascii")
    proved = bool(
        len(rows) == 98
        and len(displacements) == 84
        and row_sums.count(6) == 14
        and row_sums.count(13) == 84
        and signed_total == 17
        and parallel_square_sum == 271
        and constant_off_bin_energy == 36
        and total_off_bin_energy_at_no_collisions == 803
        and four_star_energy == 96
        and nonexact_base_energy == 707
        and fractional_point_checked
    )
    _require(proved, "the difference-Radon energy ledger failed")
    return {
        "difference_class_count": len(displacements),
        "row_index": "(projective L,a), a in F_13/+/-={0,...,6}",
        "row_count_including_zero_bins": len(rows),
        "zero_bin_row_size": 6,
        "nonzero_bin_row_size": 13,
        "incidence_matrix_sha256": hashlib.sha256(matrix_payload).hexdigest(),
        "Gram_formula": "B^T*B=13*I+2*J-G_parallel",
        "G_parallel_definition": (
            "fourteen diagonal J_6 blocks indexed by displacement direction"
        ),
        "Gram_entry_values": {
            "same_column": 14,
            "distinct_same_direction": 1,
            "different_directions": 2,
        },
        "Gram_entry_checks": gram_checks,
        "signed_transform": (
            "q_L(a)=epsilon_L*sum_(|L(delta)|=a)chi(delta)*m_delta; "
            "q_L(0)=P_L"
        ),
        "exact_inverse": (
            "13*m_delta=P_r-2*epsilon_r*T+epsilon_r*"
            "sum_L epsilon_L*q_L(|L(delta)|)"
        ),
        "hard_sign_normalized_inverse": (
            "13*m_delta=P_r-34*sigma_r+sigma_r*"
            "sum_L sigma_L*q_L(|L(delta)|), sigma=+1 hard and -1 opposite"
        ),
        "inverse_coefficient_checks": inverse_checks,
        "branch_signed_total_T_over_h": signed_total,
        "parallel_count_multiset": parallel_counts,
        "parallel_square_sum": parallel_square_sum,
        "off_bin_parseval": "sum_(L,a>0)q_L(a)^2=13*sum_delta m_delta^2+36",
        "collision_parameter": "C=sum_delta binom(m_delta,2)",
        "displacement_square_sum": "sum_delta m_delta^2=59+2*C",
        "all_off_bin_energy": "803+26*C",
        "four_exact_star_off_bin_energy": four_star_energy,
        "three_elevated_plus_seven_opposite_off_bin_energy": "707+26*C",
        "uniform_fractional_assignment": "m_(r,b)=P_r/6",
        "uniform_fractional_q_values": {
            name: str(value) for name, value in fractional_q_values.items()
        },
        "uniform_fractional_translated_interval_cut_values": {
            name: str(value)
            for name, value in fractional_translated_cut_values.items()
        },
        "uniform_fractional_point_checked": fractional_point_checked,
        "uniform_fractional_point_satisfies_integrality": False,
        "uniform_fractional_point_tests_quartic_or_midpoint_constraints": False,
        "quartic_value_code": (
            "sum_(a=1..6) a^4*q_L(a)=epsilon_L*M4(L)"
        ),
        "rational_interior_warning": (
            "the bare direction-count plus four-exact-row transform has a "
            "strictly positive fractional point satisfying the translated "
            "cuts; it has M4=0 and does not test the live quartic or midpoint "
            "constraints"
        ),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p13_opposite_entry_alphabet_certificate() -> dict[str, object]:
    """Derive the integral opposite entry, degree, and pair bounds."""
    total = -20
    cut_upper = -10

    # Hypergeometric conditional-cut coefficients.  In each line the
    # unspecified vertices are sampled uniformly to complete a seven-set.
    vertex_inside_other_cross = Fraction(
        2 * (M - 1) * (P - M), (P - 1) * (P - 2)
    )
    vertex_inside_incident_cross = Fraction(P - M, P - 1)
    vertex_inside_S_coefficient = vertex_inside_other_cross
    vertex_inside_d_coefficient = (
        vertex_inside_incident_cross - vertex_inside_other_cross
    )
    vertex_outside_other_cross = Fraction(
        2 * M * (P - 1 - M), (P - 1) * (P - 2)
    )
    vertex_outside_incident_cross = Fraction(M, P - 1)
    vertex_outside_S_coefficient = vertex_outside_other_cross
    vertex_outside_d_coefficient = (
        vertex_outside_incident_cross - vertex_outside_other_cross
    )
    pair_inside_common_cross = Fraction(
        2 * (M - 2) * (P - M), (P - 2) * (P - 3)
    )
    _require(
        pair_inside_common_cross == Fraction(P - M, P - 2),
        "pair-inside crossing probabilities diverged",
    )
    pair_inside_S_coefficient = pair_inside_common_cross
    pair_inside_w_coefficient = -pair_inside_common_cross
    oriented_pair_other_cross = Fraction(
        2 * (M - 1) * (P - M - 1), (P - 2) * (P - 3)
    )
    oriented_pair_i_cross = Fraction(P - M - 1, P - 2)
    oriented_pair_j_cross = Fraction(M - 1, P - 2)
    oriented_pair_S_coefficient = oriented_pair_other_cross
    oriented_pair_di_coefficient = (
        oriented_pair_i_cross - oriented_pair_other_cross
    )
    oriented_pair_w_coefficient = (
        1
        - oriented_pair_i_cross
        - oriented_pair_j_cross
        + oriented_pair_other_cross
    )
    pair_outside_other_cross = Fraction(
        2 * M * (P - 2 - M), (P - 2) * (P - 3)
    )
    pair_outside_incident_cross = Fraction(M, P - 2)
    pair_outside_S_coefficient = pair_outside_other_cross
    pair_outside_degree_coefficient = (
        pair_outside_incident_cross - pair_outside_other_cross
    )
    pair_outside_w_coefficient = (
        -2 * pair_outside_incident_cross + pair_outside_other_cross
    )

    row_lower = Fraction(
        cut_upper - vertex_inside_S_coefficient * total,
        vertex_inside_d_coefficient,
    )
    row_upper = Fraction(
        cut_upper - vertex_outside_S_coefficient * total,
        vertex_outside_d_coefficient,
    )
    entry_lower = Fraction(
        cut_upper - pair_inside_S_coefficient * total,
        pair_inside_w_coefficient,
    )
    entry_upper = Fraction(
        11 * cut_upper - 6 * total + row_upper,
        6,
    )
    pair_outside_even_bound = Fraction(55 * cut_upper - 28 * total, 7)
    integer_lower = min(value for value in range(-10, 11) if value >= entry_lower)
    integer_upper = max(value for value in range(-10, 11) if value <= entry_upper)
    integer_row_lower = max(
        min(value for value in range(-30, 31, 2) if value >= row_lower),
        (P - 1) * integer_lower,
    )
    integer_row_upper = max(
        value for value in range(-30, 31, 2) if value <= row_upper
    )
    integral_pair_bound = max(
        value
        for value in range(-30, 31, 2)
        if value <= pair_outside_even_bound
    )
    proved = bool(
        vertex_inside_S_coefficient == Fraction(6, 11)
        and vertex_inside_d_coefficient == Fraction(-1, 22)
        and vertex_outside_S_coefficient == Fraction(35, 66)
        and vertex_outside_d_coefficient == Fraction(7, 132)
        and pair_inside_S_coefficient == Fraction(6, 11)
        and pair_inside_w_coefficient == Fraction(-6, 11)
        and oriented_pair_S_coefficient == Fraction(6, 11)
        and oriented_pair_di_coefficient == Fraction(-1, 11)
        and oriented_pair_w_coefficient == Fraction(6, 11)
        and pair_outside_S_coefficient == Fraction(28, 55)
        and pair_outside_degree_coefficient == Fraction(7, 55)
        and pair_outside_w_coefficient == Fraction(-42, 55)
        and row_lower == -20
        and entry_lower == Fraction(-5, 3)
        and row_upper == Fraction(80, 7)
        and entry_upper == Fraction(25, 7)
        and pair_outside_even_bound == Fraction(10, 7)
        and integer_lower == -1
        and integer_upper == 3
        and integer_row_lower == -12
        and integer_row_upper == 10
        and integral_pair_bound == 0
    )
    _require(proved, "the p=13 opposite conditioned bounds changed")
    return {
        "p": P,
        "opposite_total_coefficient_sum": total,
        "balanced_cut_upper_bound": cut_upper,
        "conditional_cut_means": {
            "vertex_inside": "(12*S-d_i)/22",
            "vertex_outside": "(70*S+7*d_i)/132",
            "pair_inside": "6*(S-w_ij)/11",
            "i_inside_j_outside": "(6*S-d_i+6*w_ij)/11",
            "pair_outside": "(28*S+7*(d_i+d_j)-42*w_ij)/55",
        },
        "vertex_inside_rational_row_lower_bound": str(row_lower),
        "pair_inside_rational_lower_bound": str(entry_lower),
        "vertex_outside_rational_row_upper_bound": str(row_upper),
        "oriented_pair_rational_entry_upper_bound": str(entry_upper),
        "pair_outside_rational_bound_on_d_i_plus_d_j_minus_6w": str(
            pair_outside_even_bound
        ),
        "integer_entry_lower_bound": integer_lower,
        "integer_entry_upper_bound": integer_upper,
        "entry_alphabet": list(range(integer_lower, integer_upper + 1)),
        "even_row_sum_bounds": [integer_row_lower, integer_row_upper],
        "pair_outside_even_integral_bound": integral_pair_bound,
        "pair_outside_integral_inequality": "d_i+d_j<=6*w_ij",
        "rational_interior_warning": (
            "-5/3 and 25/7 are strict rational bounds; {-1,...,3} follows "
            "only after using coefficient integrality"
        ),
        "proved": proved,
    }


def _rational_rank(rows: Iterable[Iterable[int]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


@lru_cache(maxsize=1)
def six_dilate_cut_energy_certificate() -> dict[str, object]:
    """Use the six interval dilates to prove the collision bound ``C<=11``."""
    multiplicative_order = (1, 2, 4, 5, 3, 6)
    interval = tuple(range(M))
    interval_vector_natural_order = translated_cut_vector(interval)
    first_row = tuple(
        interval_vector_natural_order[distance - 1]
        for distance in multiplicative_order
    )
    inverse_two = pow(2, -1, P)
    multiplier = 1
    dilate_multipliers = []
    for _ in DISTANCES:
        dilate_multipliers.append(_projected_distance(multiplier))
        multiplier = multiplier * inverse_two % P
    dilated_intervals = [
        tuple(sorted(multiplier * value % P for value in interval))
        for multiplier in dilate_multipliers
    ]
    matrix = [
        [
            translated_cut_vector(subset)[distance - 1]
            for distance in multiplicative_order
        ]
        for subset in dilated_intervals
    ]
    expected_circulant = [
        [first_row[(column + row) % 6] for column in range(6)]
        for row in range(6)
    ]
    _require(
        interval_vector_natural_order == (2, 4, 6, 8, 10, 12)
        and first_row == (2, 4, 8, 10, 6, 12)
        and tuple(dilate_multipliers) == (1, 6, 3, 5, 4, 2)
        and matrix == expected_circulant,
        "the six interval dilates no longer give the claimed circulant",
    )
    gram = [
        [sum(matrix[row][left] * matrix[row][right] for row in range(6)) for right in range(6)]
        for left in range(6)
    ]
    eigenvalue_multiplicities = {1764: 1, 100: 1, 84: 2, 76: 2}
    nullities = {
        eigenvalue: 6
        - _rational_rank(
            [
                [
                    gram[row][column]
                    - (eigenvalue if row == column else 0)
                    for column in range(6)
                ]
                for row in range(6)
            ]
        )
        for eigenvalue in eigenvalue_multiplicities
    }
    _require(
        nullities == eigenvalue_multiplicities,
        "six-dilate squared singular values changed",
    )
    _require(
        all(sum(row) == 42 for row in matrix)
        and all(sum(gram[row]) == 1764 for row in range(6)),
        "six-dilate constant eigenspace changed",
    )

    def row_energy_bound(
        *, coefficient_sum: int, cut_sum_upper: int, parity: int
    ) -> dict[str, object]:
        # r=Cq has sum 42*S.  Put y=cut_sum_upper-r >=0.
        slack_sum = 6 * cut_sum_upper - 42 * coefficient_sum
        _require(slack_sum >= 0, "negative six-dilate slack sum")
        if parity == 1:
            _require(slack_sum % 2 == 0 and slack_sum >= 6, "odd slack failed")
            reduced_sum = (slack_sum - 6) // 2
            maximum_slack_square_sum = (
                (2 * reduced_sum + 1) ** 2 + 5
            )
            extremal_slacks = [2 * reduced_sum + 1] + [1] * 5
        else:
            _require(slack_sum % 2 == 0, "even slack failed")
            reduced_sum = slack_sum // 2
            maximum_slack_square_sum = (2 * reduced_sum) ** 2
            extremal_slacks = [2 * reduced_sum] + [0] * 5
        maximum_image_energy = (
            6 * cut_sum_upper * cut_sum_upper
            - 2 * cut_sum_upper * slack_sum
            + maximum_slack_square_sum
        )
        rational_bound = (
            Fraction(maximum_image_energy, 76)
            - Fraction(1764 - 76, 76 * 6) * coefficient_sum * coefficient_sum
        )
        integer_bound = rational_bound.numerator // rational_bound.denominator
        return {
            "coefficient_sum": coefficient_sum,
            "six_cut_sum_upper": cut_sum_upper,
            "image_sum": 42 * coefficient_sum,
            "slack_sum": slack_sum,
            "slack_parity": "odd positive" if parity else "even nonnegative",
            "extremal_slack_vector_up_to_permutation": extremal_slacks,
            "maximum_slack_square_sum": maximum_slack_square_sum,
            "maximum_image_energy": maximum_image_energy,
            "raw_rational_q_energy_bound": str(rational_bound),
            "integer_q_energy_bound": integer_bound,
        }

    elevated = row_energy_bound(coefficient_sum=11, cut_sum_upper=91, parity=1)
    opposite = row_energy_bound(coefficient_sum=-20, cut_sum_upper=-130, parity=0)
    nonstar_energy_upper = (
        ELEVATED_HARD_DIRECTION_COUNT * elevated["integer_q_energy_bound"]
        + OPPOSITE_DIRECTION_COUNT * opposite["integer_q_energy_bound"]
    )
    collision_upper = (nonstar_energy_upper - 707) // 26
    matched_local_energy = (
        ELEVATED_HARD_DIRECTION_COUNT * 29 + OPPOSITE_DIRECTION_COUNT * 68
    )
    proved = bool(
        elevated["raw_rational_q_energy_bound"] == "4952/57"
        and elevated["integer_q_energy_bound"] == 86
        and opposite["raw_rational_q_energy_bound"] == "6050/57"
        and opposite["integer_q_energy_bound"] == 106
        and nonstar_energy_upper == 1000
        and collision_upper == 11
        and matched_local_energy == 563
        and matched_local_energy < 707
    )
    _require(proved, "the six-dilate collision bound failed")
    return {
        "base_interval_seven_set": list(interval),
        "interval_cut_vector_natural_order": list(interval_vector_natural_order),
        "multiplicative_distance_order": list(multiplicative_order),
        "dilate_multipliers_generated_by_inverse_two_modulo_sign": dilate_multipliers,
        "dilated_interval_seven_sets": [list(row) for row in dilated_intervals],
        "interval_cut_vector_in_that_order": list(first_row),
        "six_by_six_circulant": matrix,
        "squared_singular_value_multiplicities": {
            str(value): multiplicity
            for value, multiplicity in eigenvalue_multiplicities.items()
        },
        "spectral_floor_off_constants": 76,
        "elevated_row": elevated,
        "opposite_row": opposite,
        "rational_interior_warning": (
            "the spectral bounds are 4952/57 and 6050/57, not 86 and 106; "
            "the integer bounds follow only because sum_a q_L(a)^2 is integral"
        ),
        "nonstar_energy_upper_bound": nonstar_energy_upper,
        "combined_with_parseval_identity": "707+26*C<=1000",
        "collision_parameter_upper_bound": collision_upper,
        "matched_lambda_seven_local_rows_total_energy": matched_local_energy,
        "matched_local_rows_fail_common_energy_identity": True,
        "branch_excluded": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def four_star_moment_theorem() -> dict[str, object]:
    """Package the root-count implication and its Proposition 15.740 hinge."""
    stars = exact_star_moment_certificate()
    roots = quartic_root_rank_certificate()
    prop15740_certificate = translated_cut_nine_vector_certificate()
    dependency_ok = bool(
        prop15740_certificate["proved"]
        and prop15740_certificate["moment_degrees"] == [2, 4]
        and prop15740_certificate["candidate_count_after_sum_l1_moments"]
        == MOMENT_CANDIDATE_COUNT
        and prop15740_certificate["remaining_after_nine_vectors"] == 0
        and prop15740_certificate["pure_integer_enumeration_infeasible"]
    )
    proved = bool(stars["proved"] and roots["proved"] and dependency_ok)
    _require(proved, "the four-star moment theorem failed")
    return {
        "hypothesis": (
            "one common p=13 graph has four distinct exact positive-star directions"
        ),
        "global_forms": {
            "M2": "sum_e chi(delta_e)*(L(u)-L(v))^2",
            "T3": "sum_e chi(delta_e)*(L(u)+L(v))*(L(u)-L(v))^2",
            "M4": "sum_e chi(delta_e)*(L(u)-L(v))^4",
            "U4": "sum_e chi(delta_e)*(L(u)+L(v))^2*(L(u)-L(v))^2",
        },
        "homogeneous_degrees": {"M2": 2, "T3": 3, "M4": 4, "U4": 4},
        "forced_global_identities": ["M2=0", "T3=0", "U4=lambda*M4"],
        "quartic_span_rank_at_most_one": True,
        "M4_nonzero_dependency": {
            "proposition": "15.740",
            "reason": (
                "if M4=0, then M2=M4=0 and the nine-vector opposite-cell "
                "certificate is infeasible"
            ),
            "live_certificate_proved": prop15740_certificate["proved"],
            "moment_degrees": prop15740_certificate["moment_degrees"],
            "candidate_count": prop15740_certificate[
                "candidate_count_after_sum_l1_moments"
            ],
            "remaining_after_nine_vectors": prop15740_certificate[
                "remaining_after_nine_vectors"
            ],
        },
        "M4_nonzero": True,
        "M4_projective_zero_set": "exactly the four exact hard directions",
        "three_elevated_and_seven_opposite_M4_values_nonzero": True,
        "unique_lambda_in_F13": True,
        "local_normalized_equations": [
            "sum_(s<t) W_st*(s-t)^2=0",
            "sum_(s<t) W_st*(s+t)*(s-t)^2=0",
            "sum_(s<t) W_st*((s+t)^2-lambda*(s-t)^2)*(s-t)^2=0",
            "sum_(s<t) W_st*(s-t)^4!=0 in every nonexact direction",
        ],
        "global_quartic_value_code": (
            "sum_(a=1..6) a^4*q_L(a)=epsilon_L*M4(L), with "
            "M4 a nonzero scalar times the product of the four root factors"
        ),
        "exact_star_certificate": stars,
        "root_rank_certificate": roots,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15741() -> dict[str, object]:
    """Package the proved common-graph moment reduction and sharp barrier."""
    theorem = four_star_moment_theorem()
    basis = endpoint_contraction_basis_certificate()
    affine = affine_label_invariance_certificate()
    elevated = elevated_lambda_seven_certificate()
    opposite = opposite_lambda_seven_certificate()
    previous = previous_elevated_witness_cubic_obstruction()
    collision = exact_star_collision_certificate()
    midpoint = midpoint_displacement_certificate()
    difference_radon = difference_radon_gram_certificate()
    opposite_alphabet = p13_opposite_entry_alphabet_certificate()
    six_dilate = six_dilate_cut_energy_certificate()
    common_lambda = elevated["lambda_U4_over_M4"]
    proved = bool(
        theorem["proved"]
        and basis["proved"]
        and affine["proved"]
        and elevated["proved"]
        and opposite["proved"]
        and previous["proved"]
        and collision["proved"]
        and midpoint["proved"]
        and difference_radon["proved"]
        and opposite_alphabet["proved"]
        and six_dilate["proved"]
        and common_lambda == opposite["lambda_U4_over_M4"] == 7
    )
    _require(proved, "Proposition 15.741 package failed")
    return {
        "prop": "15.741",
        "title": "four-star mixed moments and the common-graph frontier",
        "result_status": "open reduction",
        "p": P,
        "layer_index_t": 3,
        "original_k": 4 * P + 6,
        "remaining_hard_quotient_partition": [1, 1, 1, 1, 2, 2, 2],
        "moment_theorem": theorem,
        "endpoint_contraction_basis": basis,
        "affine_label_invariance": affine,
        "matched_local_lambda": common_lambda,
        "elevated_local_cell": elevated,
        "opposite_local_cell": opposite,
        "previous_witness_cubic_obstruction": previous,
        "exact_star_collision": collision,
        "midpoint_displacement_formulation": midpoint,
        "difference_radon_transform": difference_radon,
        "opposite_entry_alphabet": opposite_alphabet,
        "six_dilate_cut_energy": six_dilate,
        "strict_advance_over_15_740": (
            "adds the global cubic T3=0, quartic rank one, midpoint moments, "
            "and exact four-star collision identities"
        ),
        "method_barrier": (
            "matched lambda=7 elevated and opposite one-direction cells show "
            "that the independent cellwise scalar consequences M2=T3=0 and "
            "U4=lambda*M4 do not locally exclude either cell type; they do "
            "not witness one common quartic or difference transform"
        ),
        "constructs_common_59_edge_graph": False,
        "p13_generic_four_exact_partition_closed": False,
        "p13_generic_t3_branch_closed": False,
        "p13_k_eq_58_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_exact_gate": (
            "first exclude the 84-class nonnegative integer difference system "
            "with C<=11, exact rows, the full nonzero quartic value code, and "
            "the ten cut/parity lifts; only if it survives, exclude its "
            "14196-variable binary midpoint lift using T3=0, "
            "U4=lambda*M4, simplicity, and the exact fibre-pair equations"
        ),
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15741.json"
    target.write_text(json.dumps(proposition_15741(), indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    theorem = proposition_15741()
    target = write_evidence()
    print(
        "Prop. 15.741: M2=T3=0 and U4=lambda*M4; "
        "matched local cells leave the common graph open"
    )
    print(f"matched local lambda={theorem['matched_local_lambda']}")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
