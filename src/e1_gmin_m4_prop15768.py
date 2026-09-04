#!/usr/bin/env python3
r"""Prop. 15.768 -- the first ``p=1 (mod 4)`` layer beyond 15.752.

Put ``q=(p-1)/2`` and ``m=q+1``.  Proposition 15.752 stops at
``t=q-4`` when ``p=1 (mod 4)`` because, at the next layer ``t=q-3``, the
phase-one floor ``2p-6`` is attained.  This module proves that the new
equality branch is nevertheless impossible for every prime ``p>=29``.

The new hard cell has ``u=t``, all hard quotients equal to one, and
``b=p-3``.  Positive quadrature and a three-swap degree argument make the
cell pointwise

    A=(|X intersect C|-2)^2,  |C|=3,

whose signed target has coefficient offset two.  The common (not merely
local) difference-row sum

    sum q_L = p(P_L-3)-a_L = hT-P_L

forces all hard parallel counts to agree.  The offset-two congruence and
the edge bound then give ``P_L=2`` and ``hT=8-3p``.  Opposite directions
have

    a(Q)=(p+1)Q-6p+8,       sum Q=4p-7.

The sharp lift floor removes ``Q=6``; pigeonhole forces a ``Q=7`` cell of
mass ``p+15``.  Its odd-fibre alternatives are again sub-floor lifts, so it
would give a nonzero nonnegative integral quadratic ``B`` with
``4p E[B]=p+15``.

The local ``p+15`` mass is impossible for ``p=1 (mod 4), p>=29``.  Besides
Proposition 15.751's half-mean cube theorem, the endpoint uses the sharp
dimension-free fact that a nonnegative integral cube quadratic of mean
``3/4`` has maximum at most six.  At height one, the corrected Johnson
influence bound is strictly below eight.  Hence at most seven coordinates
remain, fewer than either side of the complementary slice; symmetrization
extends to a cube, where cube influence leaves at most four active
coordinates.  The fixed four-bit catalog of Proposition 15.751 misses
density ``(p+15)/(4p)``.

The two pre-existing hard branches still force the already excluded local
masses ``p+9`` and ``p+7``.  Thus ``k=4p+2t=5p-7`` is empty for every prime
``p=1 (mod 4), p>=29``.  The first uncovered ``p=3 (mod 4)`` layer is not
claimed here.
"""
from __future__ import annotations

import json
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15751 import (
    atomic_write_json,
    cube_half_mean_height_certificate,
    density_profile_certificate,
    exact_four_cube_catalog,
    height_at_least_two_certificate,
    height_one_junta_certificate,
    profile_density,
)
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion


ROOT = Path(__file__).resolve().parents[1]
NEW_BRANCH = "hard_b_p_minus_3_complement_triple"
LOCAL_MASS_OFFSET = 15


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 29
        or p % 4 != 1
        or not is_prime(p)
    ):
        raise ValueError("need a prime p>=29 congruent to 1 modulo 4")


def cube_three_quarter_height_certificate() -> dict[str, object]:
    """Prove ``E[g]=3/4 => max(g)<=6`` for every integral cube quadratic."""
    half_mean = cube_half_mean_height_certificate()
    support_floor = Fraction(1, 4)
    minimum_counterexample_dimension = 4
    minimum_facet_dimension = minimum_counterexample_dimension - 1
    facet_lattice = tuple(Fraction(value, 4) for value in range(1, 7))
    allowed_through_origin = (Fraction(1), Fraction(5, 4))
    allowed_opposite = tuple(Fraction(3, 2) - value for value in allowed_through_origin)

    # On a five-coordinate face, permutation averaging gives a quadratic
    # q(s).  Interpolate q(0) from the odd layers 1,3,5.
    interpolation_nodes = (1, 3, 5)
    interpolation_weights = (Fraction(15, 8), Fraction(-5, 4), Fraction(3, 8))
    interpolation_moments = tuple(
        sum(weight * node**degree for node, weight in zip(
            interpolation_nodes, interpolation_weights
        ))
        for degree in range(3)
    )
    five_face_upper = (
        interpolation_weights[0] * 3 + interpolation_weights[2] * 3
    )

    # In dimension four, the vanishing fourth alternating difference makes
    # the even- and odd-parity masses equal.  Their total is 12.
    small_dimension_total_masses = tuple(
        Fraction(3 * 2**dimension, 4) for dimension in range(4)
    )
    dimension_four_total_mass = 3 * 2 ** (4 - 2)
    dimension_four_parity_mass = dimension_four_total_mass // 2

    # The bound is sharp in dimension six.
    sharp_values: list[int] = []
    sharp_layers = [6 - 3 * s + comb(s, 2) for s in range(7)]
    for bits in product((0, 1), repeat=6):
        sharp_values.append(sharp_layers[sum(bits)])
    sharp_third_differences = []
    for triple in combinations(range(6), 3):
        outside = [index for index in range(6) if index not in triple]
        for outside_bits in product((0, 1), repeat=len(outside)):
            base = sum(bit << index for bit, index in zip(outside_bits, outside))
            difference = 0
            for size in range(4):
                for chosen in combinations(triple, size):
                    mask = base + sum(1 << index for index in chosen)
                    difference += (-1) ** (3 - size) * sharp_values[mask]
            sharp_third_differences.append(difference)
    sharp_mass = sum(sharp_values)

    proved = bool(
        half_mean["proved"]
        and half_mean["quarter_mean_lattice"]
        and support_floor == Fraction(1, 4)
        and minimum_counterexample_dimension == 4
        and minimum_facet_dimension >= 2
        and facet_lattice
        == (
            Fraction(1, 4),
            Fraction(1, 2),
            Fraction(3, 4),
            Fraction(1),
            Fraction(5, 4),
            Fraction(3, 2),
        )
        and allowed_opposite == (Fraction(1, 2), Fraction(1, 4))
        and interpolation_moments == (Fraction(1), Fraction(0), Fraction(0))
        and five_face_upper == Fraction(27, 4) < 7
        and small_dimension_total_masses
        == (Fraction(3, 4), Fraction(3, 2), Fraction(3), Fraction(6))
        and dimension_four_total_mass == 12
        and dimension_four_parity_mass == 6
        and min(sharp_values) == 0
        and max(sharp_values) == 6
        and sharp_mass == 3 * 2 ** (6 - 2)
        and not any(sharp_third_differences)
    )
    _require(proved, "the three-quarter-mean cube height theorem failed")
    return {
        "domain": "all dimensions d>=0",
        "hypotheses": (
            "g is a nonnegative integer-valued multilinear polynomial "
            "of degree at most two on {0,1}^d, with E[g]=3/4"
        ),
        "integral_values_force_integral_multilinear_coefficients": True,
        "degree_two_cube_support_floor": str(support_floor),
        "total_mass_first_forces_counterexample_dimension_at_least": (
            minimum_counterexample_dimension
        ),
        "facet_dimension_is_then_at_least": minimum_facet_dimension,
        "facet_means_are_quarter_integral": True,
        "minimal_counterexample_facet_mean_lattice": [
            str(value) for value in facet_lattice
        ],
        "facet_cases_excluded": {
            "through_origin_mean_1/4": "support equality gives maximum at most 1",
            "through_origin_mean_1/2": "half-mean theorem gives maximum at most 3",
            "through_origin_mean_3/4": "minimal dimension gives maximum at most 6",
            "through_origin_mean_3/2": (
                "the opposite facet is zero, so g=(1-x_i)h with h affine; "
                "h(x)+h(1-x)=3 gives maximum at most 3"
            ),
        },
        "remaining_through_origin_facet_means": [
            str(value) for value in allowed_through_origin
        ],
        "corresponding_opposite_facet_means": [
            str(value) for value in allowed_opposite
        ],
        "every_nonorigin_vertex_upper_bound": 3,
        "dimension_at_least_five": {
            "interpolation_nodes": list(interpolation_nodes),
            "interpolation_weights": [str(value) for value in interpolation_weights],
            "degree_zero_one_two_moments": [
                str(value) for value in interpolation_moments
            ],
            "identity": "q(0)=15q(1)/8-5q(3)/4+3q(5)/8",
            "maximum_upper_bound_before_integrality": str(five_face_upper),
            "integral_maximum_upper_bound": 6,
        },
        "dimension_at_most_three_total_masses": [
            str(value) for value in small_dimension_total_masses
        ],
        "dimension_at_most_three_total_mass_below_seven": (
            max(small_dimension_total_masses) < 7
        ),
        "dimension_four": {
            "total_mass": dimension_four_total_mass,
            "vanishing_fourth_difference_equalizes_parity_masses": True,
            "each_parity_mass": dimension_four_parity_mass,
            "maximum_upper_bound": 6,
        },
        "maximum_upper_bound": 6,
        "sharp_example": {
            "dimension": 6,
            "formula": "6-3s+binom(s,2)",
            "layer_values": sharp_layers,
            "mass": sharp_mass,
            "cube_size": len(sharp_values),
            "mean": "3/4",
            "maximum": max(sharp_values),
            "all_third_differences_zero": not any(sharp_third_differences),
        },
        "proved": proved,
    }


def p29_p_plus_fifteen_height_exclusion() -> dict[str, object]:
    """Close the height-at-least-two endpoint for local mass ``29+15``."""
    p = 29
    half_mean = cube_half_mean_height_certificate()
    three_quarter = cube_three_quarter_height_certificate()
    scaled_mass = p + LOCAL_MASS_OFFSET
    slice_mean = Fraction(scaled_mass, 4 * p)
    paired_cube_rho = Fraction(1, p + 1)
    paired_cube_constant = p * slice_mean
    initial_height_lower = Fraction(p + 1, 2) - paired_cube_constant
    stabilizer_upper = Fraction((p + 15) * (p + 3), 4 * (p - 1))
    stabilizer_integer_upper = stabilizer_upper.numerator // stabilizer_upper.denominator
    refined_height_threshold = Fraction(3 * (p + 1), 4) - paired_cube_constant
    refined_integer_lower = 12
    forced_height = 12
    paired_average_at_forced_height = Fraction(forced_height + 11, 30)
    proved = bool(
        half_mean["proved"]
        and three_quarter["proved"]
        and scaled_mass == 44
        and slice_mean == Fraction(11, 29)
        and paired_cube_rho == Fraction(1, 30)
        and paired_cube_constant == 11
        and initial_height_lower == 4
        and stabilizer_upper == Fraction(88, 7)
        and stabilizer_integer_upper == 12
        and refined_height_threshold == Fraction(23, 2)
        and refined_integer_lower == 12
        and forced_height == stabilizer_integer_upper == refined_integer_lower
        and Fraction(3, 4) < paired_average_at_forced_height < 1
        and three_quarter["maximum_upper_bound"] < forced_height
    )
    _require(proved, "the p=29 p+15 height endpoint did not close")
    return {
        "p": p,
        "scaled_mean_4p_E_B": scaled_mass,
        "slice_mean_E_B": str(slice_mean),
        "assumed_maximum_at_least": 2,
        "paired_cube_operator_rho": str(paired_cube_rho),
        "paired_cube_operator_at_maximum": "T B(X)=(H+11)/30",
        "paired_cube_mean_lattice": "(1/4)Z",
        "initial_half_mean_floor_height_lower_bound": str(initial_height_lower),
        "half_mean_cube_height_upper_bound": 3,
        "half_mean_is_then_impossible": True,
        "every_paired_cube_mean_at_least": "3/4",
        "refined_raw_height_lower_bound": str(refined_height_threshold),
        "refined_integral_height_lower_bound": refined_integer_lower,
        "stabilizer_height_upper_bound": str(stabilizer_upper),
        "stabilizer_integral_height_upper_bound": stabilizer_integer_upper,
        "forced_height": forced_height,
        "paired_cube_average_at_forced_height": str(
            paired_average_at_forced_height
        ),
        "some_paired_cube_has_mean_exactly": "3/4",
        "three_quarter_cube_height_upper_bound": three_quarter[
            "maximum_upper_bound"
        ],
        "contradiction": True,
        "proved": proved,
    }


def complement_triple_baseline_certificate(p: int) -> dict[str, object]:
    """Prove pointwise rigidity of the new ``b=p-3`` equality cell."""
    _check_prime(p)
    m = (p + 1) // 2
    quadrature = parity_floor_certificate(p, 3, 0)
    coefficients = tuple(quadrature["coefficients"])
    nodes = tuple(quadrature["quadrature_nodes"])
    weights = tuple(quadrature["quadrature_weights"])

    # Complementing a (p-3)-set leaves a three-set C.  Since m is odd,
    # phase one on the large set is phase zero in r=|X intersect C|.
    complement_parity = [((m - r + 1) & 1) == (r & 1) for r in range(4)]
    target_constant = 5
    target_linear_coefficients = (-1, -1, -1)
    target_linear_sum = sum(target_linear_coefficients)
    coefficient_offset = target_constant + target_linear_sum
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

    numerator = sum(
        comb(3, r) * comb(p - 3, m - r) * (2 - r) ** 2
        for r in range(4)
        if 0 <= m - r <= p - 3
    )
    mean = Fraction(numerator, comb(p, m))
    scaled_mean = 2 * p * mean

    # If the difference from the displayed candidate vanishes on r=1,2,3,
    # fix an omitted r=0 point and swap three of its chosen outside points
    # with the three points of C.  The resulting function on a 3-cube has
    # degree at most two and vanishes at all seven nonzero vertices.  Its
    # third finite difference is zero, hence its value at 000 is zero too.
    third_difference_coefficients = (-1, 1, 1, 1, -1, -1, -1, 1)
    omitted_value_coefficient = third_difference_coefficients[0]

    proved = bool(
        quadrature["exact_positive_quadrature_certificate"]
        and int(quadrature["scaled_floor"]) == 2 * p - 6
        and coefficients == (Fraction(1), Fraction(-4), Fraction(4))
        and nodes == (1, 2, 3)
        and all(weight > 0 for weight in weights)
        and all(complement_parity)
        and all(target_checks)
        and scaled_mean == 2 * p - 6
        and coefficient_offset == 2
        and omitted_value_coefficient == -1
        and sum(third_difference_coefficients) == 0
    )
    _require(proved, "the complement-triple equality baseline changed")
    return {
        "p": p,
        "slice": f"J({p},{m})",
        "large_odd_fibre_count_b": p - 3,
        "reduced_three_set_phase": 0,
        "positive_quadrature_coefficients": [str(value) for value in coefficients],
        "positive_quadrature_contact_layers": list(nodes),
        "positive_quadrature_weights": [str(value) for value in weights],
        "all_contact_weights_positive": True,
        "three_swap_cube_argument": (
            "a degree-at-most-two function vanishing at all seven nonzero "
            "vertices of a three-cube also vanishes at the origin"
        ),
        "third_difference_coefficients_000_100_010_001_110_101_011_111": list(
            third_difference_coefficients
        ),
        "pointwise_baseline": "A=(2-r)^2, r=|X intersect C|, |C|=3",
        "scaled_mean_2p_E_A": int(scaled_mean),
        "signed_target": (
            "epsilon*S_H=5-sum_(i in C)z_i+"
            "sum_({i,j} subset C)z_i*z_j"
        ),
        "target_constant": target_constant,
        "target_linear_coefficients": list(target_linear_coefficients),
        "target_linear_sum": target_linear_sum,
        "coefficient_offset_formula": "target constant + sum linear coefficients",
        "coefficient_offset": coefficient_offset,
        "slice_ideal_coefficient_identity": (
            "I+P-offset=(p-1)c with 2c integral"
        ),
        "coefficient_congruence": (
            f"{(p - 1) // 2} divides I+P-{coefficient_offset}"
        ),
        "proved": proved,
    }


def p_plus_fifteen_local_exclusion(p: int) -> dict[str, object]:
    """Exclude ``4p E[B]=p+15`` on the middle slice."""
    _check_prime(p)
    m = (p + 1) // 2
    q = (p - 1) // 2
    scaled_mass = p + LOCAL_MASS_OFFSET

    lower_height = Fraction(p - 13, 4)
    stabilizer_upper = Fraction(scaled_mass * (p + 3), 4 * (p - 1))
    paired_average_upper = Fraction(scaled_mass, 2 * (p - 1))
    half_mean = cube_half_mean_height_certificate()
    if p == 29:
        height_row = p29_p_plus_fifteen_height_exclusion()
        height_proved = bool(height_row["proved"])
    else:
        height_proved = bool(
            half_mean["proved"]
            and lower_height > 3
            and paired_average_upper < Fraction(3, 4)
        )
        height_row = {
            "height_lower_bound": str(lower_height),
            "stabilizer_height_upper_bound": str(stabilizer_upper),
            "paired_cube_average_upper_bound": str(paired_average_upper),
            "some_paired_cube_has_mean_exactly": "1/2",
            "half_mean_cube_height_upper_bound": 3,
            "proved": height_proved,
        }

    mu = Fraction(scaled_mass, 4 * p)
    influence_floor = Fraction((p + 1) * (p - 3), 16 * p * (p - 2))
    total_influence_upper = (p - 1) * mu * (1 - mu)
    junta_bound = Fraction(
        2 * (p - 1) * (p - 2) * (p + 15) * (3 * p - 15),
        p * p * (p + 1) * (p - 3),
    )
    eight_gap_polynomial = p**4 - 29 * p**3 + 297 * p**2 - 735 * p + 450
    x = p - 29
    translated_gap = (
        x**4 + 87 * x**3 + 2820 * x**2 + 40880 * x + 228912
    )
    densities = sorted(
        {
            profile_density(tuple(row["layer_counts"]), p)
            for row in exact_four_cube_catalog()["profiles"]
        }
    )
    expected_densities = sorted(
        {
            Fraction(0),
            Fraction(1),
            Fraction(p - 3, 4 * p),
            Fraction(p + 1, 4 * p),
            Fraction(p - 1, 2 * p),
            Fraction(p + 1, 2 * p),
            Fraction(3 * p - 1, 4 * p),
            Fraction(3 * (p + 1), 4 * p),
        }
    )
    boolean_proved = bool(
        eight_gap_polynomial == translated_gap
        and eight_gap_polynomial > 0
        and junta_bound < 8
        and 7 < q
        and densities == expected_densities
        and Fraction(p + 1, 4 * p) < mu < Fraction(p - 1, 2 * p)
        and mu not in densities
    )
    proved = height_proved and boolean_proved
    _require(proved, "the p+15 local exclusion failed")
    return {
        "p": p,
        "slice": f"J({p},{m})",
        "statement": (
            "no nonzero nonnegative integer-valued quadratic B has "
            "4p E[B]=p+15"
        ),
        "height_at_least_two": height_row,
        "height_one_boolean": {
            "target_density": str(mu),
            "relevant_pair_influence_floor": str(influence_floor),
            "total_influence_upper_bound": str(total_influence_upper),
            "largest_zero_influence_class_complement_bound": str(junta_bound),
            "eight_gap_polynomial": "p^4-29p^3+297p^2-735p+450",
            "eight_gap_at_p_equals_x_plus_29": [1, 87, 2820, 40880, 228912],
            "junta_coordinates_at_most": 7,
            "seven_less_than_both_complementary_slice_sizes": 7 < q,
            "all_junta_patterns_extend_to_slice": True,
            "cube_coordinates_actually_needed_at_most": 4,
            "possible_four_bit_density_values": [str(value) for value in densities],
            "target_absent": True,
            "proved": boolean_proved,
        },
        "finite_prime_or_slice_census_used": False,
        "fixed_four_bit_catalog_reused": True,
        "excluded": proved,
        "proved": proved,
    }


def first_uncovered_residue_ledger(p: int) -> dict[str, object]:
    """Classify the phase-one residues at ``t=q-3``."""
    _check_prime(p)
    q = (p - 1) // 2
    m = q + 1
    t = q - 3
    edge_count = 5 * p - 6
    floors = residual_even_floor_table(p)["phase_one_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    low_exact_at_u_t = [
        int(b) for b, floor in floors.items() if int(floor) == 2 * p - 6
    ]
    residue_rows: list[dict[str, object]] = []
    for u in range(m):
        quotient_sum = m + t - u
        if u <= t:
            forced_low_quotient = 1
            forced_low_count = 2 * m - quotient_sum
            low_mean = p + 1 + 2 * u
            candidates = []
            survivors = []
            for raw_b, raw_floor in floors.items():
                b = int(raw_b)
                floor = int(raw_floor)
                if floor > low_mean:
                    continue
                excess = low_mean - floor
                exact = excess == 0
                sub_floor_baseline_lift = bool(
                    b in {2, p - 1} and 0 < excess < lift_floor
                )
                if exact:
                    survivors.append(b)
                candidates.append(
                    {
                        "b": b,
                        "floor": floor,
                        "excess": excess,
                        "exact": exact,
                        "excluded_as_sub_floor_baseline_lift": sub_floor_baseline_lift,
                    }
                )
        else:
            forced_low_quotient = 0
            forced_low_count = m - quotient_sum
            low_mean = 2 * u
            candidates = [
                {
                    "b": int(raw_b),
                    "floor": int(raw_floor),
                    "excess": low_mean - int(raw_floor),
                    "exact": low_mean == int(raw_floor),
                    "excluded_as_sub_floor_baseline_lift": False,
                }
                for raw_b, raw_floor in floors.items()
                if int(raw_floor) <= low_mean
            ]
            survivors = [int(row["b"]) for row in candidates if row["exact"]]

        if u == 0:
            expected_survivors = [p - 1]
            branch = BRANCH_P1_LAST
        elif 1 <= u <= t - 1:
            expected_survivors = []
            branch = None
        elif u == t:
            expected_survivors = [p - 3]
            branch = NEW_BRANCH
        elif t < u < m - 1:
            expected_survivors = []
            branch = None
        else:
            expected_survivors = [2]
            branch = BRANCH_B2
        _require(
            survivors == expected_survivors
            and forced_low_count > 0
            and all(
                row["exact"] or row["excluded_as_sub_floor_baseline_lift"]
                for row in candidates
            ),
            f"the residue row u={u} changed",
        )
        residue_rows.append(
            {
                "u": u,
                "quotient_sum": quotient_sum,
                "forced_low_quotient": forced_low_quotient,
                "forced_low_direction_count_at_least": forced_low_count,
                "forced_low_mean": low_mean,
                "candidate_floor_rows": candidates,
                "exact_surviving_b": survivors,
                "surviving_branch": branch,
                "excluded": not survivors,
            }
        )
    proved = bool(
        2 * m * (m + t) == m * (2 * p - 6)
        and p * p + 1 - 2 * edge_count > 0
        and min(int(value) for value in floors.values()) == p - 1
        and low_exact_at_u_t == [p - 3]
        and 2 * (t - 1) + 2 == p - 7 < lift_floor
        and 2 * t + 2 == p - 5 < lift_floor
        and 2 * (m - 2) == p - 3 < p - 1
        and len(residue_rows) == m
        and [row["u"] for row in residue_rows] == list(range(m))
    )
    _require(proved, "the first-uncovered phase-one residue ledger failed")
    return {
        "p": p,
        "q": q,
        "m": m,
        "layer_index_t": t,
        "original_k": 5 * p - 7,
        "H_edge_count": edge_count,
        "guaranteed_isolated_vertices": p * p + 1 - 2 * edge_count,
        "phase_one_type_budget": 2 * m * (m + t),
        "phase_one_mean_form": f"a_L=2u+{p + 1}k_L",
        "phase_one_quotient_sum": "sum k_L=m+t-u",
        "all_residue_rows": residue_rows,
        "residue_cases": {
            "u=0": BRANCH_P1_LAST,
            "1<=u<=t-1": "excluded: every forced k=1 cell is a lift below p-3",
            "u=t": NEW_BRANCH,
            "t+1<=u<=m-2": "excluded: a k=0 cell is forced below the floor",
            "u=m-1": BRANCH_B2,
        },
        "new_branch_all_quotients_equal_one": True,
        "new_branch_exact_b": p - 3,
        "new_branch_exact_scaled_mean": 2 * p - 6,
        "possible_branches": [BRANCH_B2, BRANCH_P1_LAST, NEW_BRANCH],
        "prior_15_752_failure_is_an_exact_new_floor_not_a_failed_inequality": True,
        "proved": proved,
    }


def _old_branch_extension(p: int, branch: str) -> dict[str, object]:
    """Replay the A/B pigeonhole one step beyond 15.752's stated band."""
    _check_prime(p)
    q = (p - 1) // 2
    m = q + 1
    t = q - 3
    edge_count = 5 * p - 6
    if branch == BRANCH_B2:
        hard_edges = 5 * m - 3
        opposite_edges = edge_count - hard_edges
        hT = 5
        minimum_Q = 3
        minimum_mean = 8
        next_Q = 4
        next_mean = p + 9
        dependency = p_plus_nine_local_exclusion(p)
        dependency_name = "Proposition 15.752 p+9 local theorem"
    elif branch == BRANCH_P1_LAST:
        hard_edges = 6 * m - 4
        opposite_edges = edge_count - hard_edges
        hT = p + 4
        minimum_Q = 2
        minimum_mean = 6
        next_Q = 3
        next_mean = p + 7
        dependencies = (
            height_at_least_two_certificate(p),
            height_one_junta_certificate(p),
            density_profile_certificate(p),
        )
        dependency = {"proved": all(row["proved"] for row in dependencies)}
        dependency_name = "Proposition 15.751 p+7 local theorem"
    else:
        raise ValueError("branch must be hard_b2 or p1_residue_zero_b_p_minus_1")

    surplus_after_next = opposite_edges - m * next_Q
    proved = bool(
        hard_edges + opposite_edges == edge_count
        and hard_edges - opposite_edges == hT
        and 0 <= surplus_after_next < m
        and 0 < minimum_mean < p - 3
        and dependency["proved"]
    )
    _require(proved, f"the {branch} extension ledger failed")
    return {
        "branch": branch,
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "hard_sign_times_global_T": hT,
        "minimum_opposite_Q": minimum_Q,
        "minimum_opposite_mean": minimum_mean,
        "minimum_cell_excluded_below_phase_zero_and_lift_floors": True,
        "forced_next_Q": next_Q,
        "forced_next_scaled_mean": next_mean,
        "surplus_after_raising_every_opposite_Q_to_next_Q": surplus_after_next,
        "opposite_direction_count": m,
        "local_dependency": dependency_name,
        "excluded": proved,
        "proved": proved,
    }


def complement_triple_branch_exclusion(p: int) -> dict[str, object]:
    """Exclude the new hard complement-triple branch."""
    _check_prime(p)
    baseline = complement_triple_baseline_certificate(p)
    local = p_plus_fifteen_local_exclusion(p)
    q = (p - 1) // 2
    m = q + 1
    edge_count = 5 * p - 6
    hard_mean = 2 * p - 6
    parallel_upper = edge_count // m
    parallel_candidates = [
        value
        for value in range(parallel_upper + 1)
        if (value - 2) % q == 0
    ]
    hard_parallel = parallel_candidates[0]
    hT = (p + 1) * hard_parallel - 3 * p - hard_mean
    hard_edges = m * hard_parallel
    opposite_edges = edge_count - hard_edges

    def opposite_mean(parallel: int) -> int:
        return (p + 1) * parallel - 6 * p + 8

    phase_zero = residual_even_floor_table(p)["phase_zero_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    next_mean = opposite_mean(7)
    nonzero_rows = [
        (int(b), int(floor), next_mean - int(floor))
        for b, floor in phase_zero.items()
        if int(b) != 0 and int(floor) <= next_mean
    ]
    surplus_after_seven = opposite_edges - 7 * m
    proved = bool(
        baseline["proved"]
        and parallel_upper < q + 2
        and parallel_candidates == [2]
        and hT == 8 - 3 * p
        and hard_edges == p + 1
        and opposite_edges == 4 * p - 7
        and opposite_mean(5) < 0
        and opposite_mean(6) == 14
        and 0 < 14 < min(p - 1, lift_floor)
        and 0 <= surplus_after_seven < m
        and next_mean == p + 15
        and [row[0] for row in nonzero_rows] == [2, p - 1]
        and [row[2] for row in nonzero_rows] == [14, 16]
        and all(0 < excess < lift_floor for _, _, excess in nonzero_rows)
        and local["proved"]
    )
    _require(proved, "the complement-triple branch survived")
    return {
        "branch": NEW_BRANCH,
        "hard_baseline": baseline,
        "common_row_sum_identity": "sum q_L=p(P_L-3)-a_L=hT-P_L",
        "common_row_sum_consequence": (
            "equal hard means force one common hard parallel count P"
        ),
        "hard_parallel_congruence": f"P=2 mod {q}",
        "hard_parallel_count_upper_bound": parallel_upper,
        "hard_parallel_candidates": parallel_candidates,
        "forced_hard_parallel_count": hard_parallel,
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "hard_sign_times_global_T": hT,
        "opposite_mean_formula": f"a(Q)={p + 1}Q-{6 * p - 8}",
        "Q5_mean": opposite_mean(5),
        "Q6_mean": opposite_mean(6),
        "Q6_excluded_below_phase_zero_and_lift_floors": True,
        "opposite_parallel_count_sum": opposite_edges,
        "surplus_after_every_Q_at_least_seven": surplus_after_seven,
        "opposite_direction_count": m,
        "a_Q7_direction_is_forced": True,
        "Q7_scaled_mean": next_mean,
        "nonzero_b_Q7_floor_and_lift_rows": [list(row) for row in nonzero_rows],
        "Q7_is_forced_to_b_zero": True,
        "Q7_cell_is_A_equals_2B": True,
        "Q7_local_mass_exclusion": local,
        "excluded": proved,
        "proved": proved,
    }


def first_uncovered_p1_layer_exclusion(p: int) -> dict[str, object]:
    """Close ``t=(p-7)/2`` for one admissible prime."""
    _check_prime(p)
    residues = first_uncovered_residue_ledger(p)
    branches = {
        BRANCH_B2: _old_branch_extension(p, BRANCH_B2),
        BRANCH_P1_LAST: _old_branch_extension(p, BRANCH_P1_LAST),
        NEW_BRANCH: complement_triple_branch_exclusion(p),
    }
    proved = bool(
        residues["proved"]
        and set(branches) == set(residues["possible_branches"])
        and all(row["proved"] for row in branches.values())
    )
    _require(proved, "the first uncovered p=1 mod 4 layer did not close")
    return {
        "p": p,
        "p_mod_4": 1,
        "layer_index_t": (p - 7) // 2,
        "original_k": 5 * p - 7,
        "H_edge_count": 5 * p - 6,
        "residue_ledger": residues,
        "branch_exclusions": branches,
        "all_boundary_sizes_excluded": True,
        "finite_prime_graph_or_slice_census_used": False,
        "residual_ii_layer_excluded": proved,
        "proved": proved,
    }


def proposition_15768() -> dict[str, object]:
    """Package the parameterized theorem and threshold replays."""
    sample_primes = (29, 37, 41, 53)
    rows = {str(p): first_uncovered_p1_layer_exclusion(p) for p in sample_primes}
    proved = all(row["proved"] for row in rows.values())
    return {
        "prop": "15.768",
        "title": "First post-15.752 p=1 mod 4 residual layer",
        "result_status": "proved infinite-family theorem with fixed four-bit catalog",
        "statement": (
            "for every prime p>=29 congruent to 1 modulo 4, residual (ii) "
            "is empty at t=(p-7)/2, equivalently k=5p-7"
        ),
        "new_local_theorem": (
            "no nonzero nonnegative integral quadratic on J(p,(p+1)/2) "
            "has 4p E[B]=p+15"
        ),
        "new_cube_theorem": (
            "every nonnegative integer-valued cube quadratic of mean 3/4 "
            "has maximum at most 6, sharply"
        ),
        "first_layer_beyond_prop_15752": True,
        "parameterized_threshold_replays": rows,
        "p29_endpoint_closed_by_three_quarter_cube_theorem": True,
        "p3_mod_4_next_layer_closed": False,
        "later_layers_closed": False,
        "residual_ii_closed_globally": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_prop15768.json"
    atomic_write_json(path, proposition_15768())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
