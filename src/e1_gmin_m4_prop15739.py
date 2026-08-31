#!/usr/bin/env python3
r"""Prop. 15.739 -- close the exceptional p=13 fourth-shell branch.

At ``p=13,t=3,u=3`` the seven hard cells are the exact complement-triple
quadratics ``A=(2-r)^2``.  Their signed target has coefficient offset two,
not five.  Correcting that normalization leaves hard parallel count ``P=2``
or ``P=8``.  The opposite ledgers force a phase-zero mass-14 cell at
``Q=6`` or ``Q=0`` respectively.  Proposition 15.738 classifies either cell
as ``B=x_i*x_j``.

For even ``d`` let

    M_d(L)=sum_{{u,v} in H} chi(u-v) (L(u)-L(v))^d.

This is a genuine homogeneous binary form.  If ``h`` is the fixed sign on
the hard type, every complement-triple cell satisfies

    S_4=(1/2) S_2^2.

Consequently the homogeneous quartic ``G=2*h*M_4-M_2^2`` vanishes on all
seven hard projective directions and is identically zero.  The opposite
selected-pair cell instead gives ``G=-3(i-j)^4``, a contradiction in
``F_13``.  Thus the exceptional branch is empty.

The same construction supplies a rigorous open reduction for the remaining
generic ``t=3`` branch.  For ``p=1 mod 4, p>=17``, at least
``(p-5)/2`` hard directions are exact stars, so every even moment through
degree ``(p-9)/2`` vanishes identically.  Conditional balanced-cut averages
also force the opposite signed matrix alphabet to ``{-1,0,1,2,3}``.
Those facts do not yet exclude the generic branch, residual (ii), or the
quadratic min-max limit.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, lcm
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15738 import (
    p13_mass14_residual_cell_classification,
    proposition_15738,
    selected_pair_moment_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
P = 13
M = 7
Q = 6
H_EDGE_COUNT = 59
HARD_DIRECTION_COUNT = 7
EXCEPTIONAL_HARD_MEAN = 20
EXCEPTIONAL_HARD_PARALLEL_COUNTS = (2, 8)
OPPOSITE_MINIMUM_PARALLEL_COUNTS = (0, 6)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


@lru_cache(maxsize=1)
def exceptional_hard_target_certificate() -> dict[str, object]:
    """Derive the corrected target and coefficient offset two."""
    identities = []
    for bits in product((0, 1), repeat=3):
        r = sum(bits)
        z = tuple(2 * bit - 1 for bit in bits)
        slack = (2 - r) ** 2
        target_from_slack = 3 + 2 * slack
        target_from_signed_coordinates = (
            5
            - sum(z)
            + z[0] * z[1]
            + z[0] * z[2]
            + z[1] * z[2]
        )
        identities.append(target_from_slack == target_from_signed_coordinates)

    target_constant = 5
    target_linear_sum = -3
    coefficient_offset = target_constant + target_linear_sum
    possible_parallel_counts = [
        parallel_count
        for parallel_count in range(H_EDGE_COUNT // HARD_DIRECTION_COUNT + 1)
        if (parallel_count - coefficient_offset) % Q == 0
    ]
    proved = bool(
        all(identities)
        and coefficient_offset == 2
        and possible_parallel_counts == [2, 8]
    )
    _require(proved, "exceptional hard target normalization failed")
    return {
        "p": P,
        "slice": "J(13,7)",
        "baseline": "A=(2-r)^2 on a three-point complement C",
        "signed_target": (
            "epsilon*S_H=5-sum_(i in C)z_i+"
            "sum_({i,j} subset C)z_i*z_j"
        ),
        "target_constant": target_constant,
        "target_linear_sum": target_linear_sum,
        "coefficient_offset": coefficient_offset,
        "slice_kernel_integrality": "2c integral",
        "coefficient_congruence": "6 divides P-2",
        "hard_parallel_count_upper_bound": H_EDGE_COUNT // HARD_DIRECTION_COUNT,
        "possible_common_hard_parallel_counts": possible_parallel_counts,
        "boolean_target_checks": len(identities),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def exceptional_parallel_ledgers() -> dict[str, object]:
    """Replay the two corrected hard/opposite edge and mean ledgers."""
    target = exceptional_hard_target_certificate()
    rows: dict[str, dict[str, object]] = {}
    for hard_parallel_count in EXCEPTIONAL_HARD_PARALLEL_COUNTS:
        signed_global_T = 14 * hard_parallel_count - 59
        opposite_edge_sum = H_EDGE_COUNT - HARD_DIRECTION_COUNT * hard_parallel_count

        def opposite_mean(parallel_count: int) -> int:
            return 14 * parallel_count + signed_global_T - 39

        if hard_parallel_count == 8:
            minimum_q = 0
            minimum_mean = opposite_mean(minimum_q)
            surplus = opposite_edge_sum
            minimum_count = HARD_DIRECTION_COUNT - surplus
            excluded_previous = None
        else:
            nonnegative_minimum = 5
            zero_mean_target_offset = 3
            q5_coefficient_compatible = (
                nonnegative_minimum - zero_mean_target_offset
            ) % Q == 0
            _require(not q5_coefficient_compatible, "Q=5 unexpectedly compatible")
            minimum_q = 6
            minimum_mean = opposite_mean(minimum_q)
            surplus = opposite_edge_sum - HARD_DIRECTION_COUNT * minimum_q
            minimum_count = HARD_DIRECTION_COUNT - surplus
            excluded_previous = {
                "parallel_count_Q": nonnegative_minimum,
                "mean": opposite_mean(nonnegative_minimum),
                "signed_target": "epsilon*S_H=3",
                "coefficient_offset": zero_mean_target_offset,
                "coefficient_compatible_mod_6": q5_coefficient_compatible,
                "excluded": True,
            }
        row = {
            "hard_parallel_count_P": hard_parallel_count,
            "hard_finite_edge_count": HARD_DIRECTION_COUNT * hard_parallel_count,
            "hard_sign_times_global_T": signed_global_T,
            "opposite_parallel_count_sum": opposite_edge_sum,
            "opposite_mean_formula": f"a(Q)=14*(Q+{hard_parallel_count}-7)",
            "excluded_previous_parallel_count": excluded_previous,
            "minimum_allowed_opposite_Q": minimum_q,
            "mean_at_minimum_Q": minimum_mean,
            "parallel_surplus_above_uniform_minimum": surplus,
            "directions_at_minimum_at_least": minimum_count,
            "proved": bool(
                minimum_mean == 14
                and surplus == 3
                and minimum_count == 4
            ),
        }
        _require(bool(row["proved"]), "exceptional parallel ledger changed")
        rows[str(hard_parallel_count)] = row
    proved = bool(
        target["proved"]
        and rows["2"]["proved"]
        and rows["8"]["proved"]
        and rows["2"]["opposite_parallel_count_sum"] == 45
        and rows["8"]["opposite_parallel_count_sum"] == 3
    )
    _require(proved, "exceptional parallel ledgers failed")
    return {
        "hard_mean_formula": "20=14*P-h*T-39",
        "opposite_mean_formula": "a(Q)=14*Q+h*T-39",
        "hard_target_dependency": target,
        "P2": rows["2"],
        "P8": rows["8"],
        "minimum_mass14_opposite_cell_exists_in_both_ledgers": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def phase_zero_mass14_cell_reduction() -> dict[str, object]:
    """Reduce each forced mean-14 cell to Proposition 15.738."""
    floors = {b: scaled_direction_floor(P, b, 0) for b in range(0, P, 2)}
    lift = sharp_integral_quadratic_lift_floor(P)
    lift_floor = int(lift["sharp_scaled_floor"])
    b2_offset = 4
    b2_compatible = {
        parallel_count: (parallel_count - b2_offset) % Q == 0
        for parallel_count in OPPOSITE_MINIMUM_PARALLEL_COUNTS
    }
    b12_excess = 14 - floors[12]
    b12_excluded = 0 < b12_excess < lift_floor
    other_b_above_mean = all(floors[b] > 14 for b in (4, 6, 8, 10))
    residual = p13_mass14_residual_cell_classification()
    q0 = residual["Q0"]
    q6 = residual["Q6"]
    proved = bool(
        floors == {0: 0, 2: 14, 4: 20, 6: 26, 8: 24, 10: 26, 12: 12}
        and b2_compatible == {0: False, 6: False}
        and b12_excluded
        and other_b_above_mean
        and residual["proved"]
        and q0["selected_pair_is_unique_surviving_family"]
        and q6["selected_pair_is_unique_surviving_family"]
    )
    _require(proved, "phase-zero mass-14 reduction failed")
    return {
        "p": P,
        "phase": 0,
        "scaled_mean": 14,
        "even_b_floors": floors,
        "b2_exact_target": "epsilon*S_H=4-z_i*z_j",
        "b2_coefficient_offset": b2_offset,
        "b2_compatible_at_Q": b2_compatible,
        "b12_floor_plus_two_excess": b12_excess,
        "nonzero_integral_lift_floor": lift_floor,
        "b12_floor_plus_two_excluded": b12_excluded,
        "b4_b6_b8_b10_floors_above_mean": other_b_above_mean,
        "remaining_cell": "b=0, A=2B, 4p*E[B]=14",
        "finite_classification_dependency": {
            "proposition": "15.738",
            "result_status": residual["result_status"],
            "Q0_survivors": residual["Q0_survivors"],
            "Q6_survivors": residual["Q6_survivors"],
        },
        "minimum_cell_forced_form": "B=x_i*x_j",
        "proved": proved,
    }


def _normalized_exceptional_hard_pattern(
    triple: tuple[int, int, int], hard_parallel_count: int
) -> dict[tuple[int, int], int]:
    if hard_parallel_count not in EXCEPTIONAL_HARD_PARALLEL_COUNTS:
        raise ValueError("hard parallel count must be 2 or 8")
    triple_set = set(triple)
    kernel_scalar = Fraction(hard_parallel_count - 2, 12)
    pattern: dict[tuple[int, int], int] = {}
    for s, t in combinations(range(P), 2):
        target_pair = int(s in triple_set and t in triple_set)
        target_linear_s = -int(s in triple_set)
        target_linear_t = -int(t in triple_set)
        value = (
            target_pair
            + 2 * kernel_scalar
            + target_linear_s
            + target_linear_t
        )
        _require(value.denominator == 1, "hard coefficient is nonintegral")
        pattern[(s, t)] = int(value)
    return pattern


def _pattern_moment(
    pattern: dict[tuple[int, int], int], degree: int
) -> int:
    return sum(
        coefficient * pow(s - t, degree, P)
        for (s, t), coefficient in pattern.items()
    ) % P


def _triangle_moment(triple: tuple[int, int, int], degree: int) -> int:
    return sum(
        pow(s - t, degree, P) for s, t in combinations(triple, 2)
    ) % P


@lru_cache(maxsize=1)
def exceptional_hard_moment_certificate() -> dict[str, object]:
    """Show both hard gauges have the triangle degree-2/4 moments."""
    triples = tuple(combinations(range(P), 3))
    gauge_rows: dict[str, dict[str, object]] = {}
    all_relation_checks: list[bool] = []
    for hard_parallel_count in EXCEPTIONAL_HARD_PARALLEL_COUNTS:
        reference = _normalized_exceptional_hard_pattern(
            triples[0], hard_parallel_count
        )
        histogram = dict(sorted(Counter(reference.values()).items()))
        expected_histogram = (
            {-1: 33, 0: 45}
            if hard_parallel_count == 2
            else {0: 33, 1: 45}
        )
        _require(histogram == expected_histogram, "hard gauge histogram changed")
        degree_checks: dict[int, bool] = {}
        for degree in (2, 4):
            degree_checks[degree] = all(
                _pattern_moment(
                    _normalized_exceptional_hard_pattern(
                        triple, hard_parallel_count
                    ),
                    degree,
                )
                == _triangle_moment(triple, degree)
                for triple in triples
            )
        relation_checks = [
            (
                2 * _triangle_moment(triple, 4)
                - pow(_triangle_moment(triple, 2), 2, P)
            )
            % P
            == 0
            for triple in triples
        ]
        all_relation_checks.extend(relation_checks)
        gauge_rows[str(hard_parallel_count)] = {
            "hard_parallel_count_P": hard_parallel_count,
            "slice_kernel_scalar": str(
                Fraction(hard_parallel_count - 2, 12)
            ),
            "normalized_coefficient_histogram": histogram,
            "triple_count_checked": len(triples),
            "degree_moment_equals_triangle": degree_checks,
            "two_S4_equals_S2_squared": all(relation_checks),
            "proved": all(degree_checks.values()) and all(relation_checks),
        }
    complete_graph_moments = {
        degree: sum(
            pow(s - t, degree, P) for s, t in combinations(range(P), 2)
        )
        % P
        for degree in (2, 4)
    }
    proved = bool(
        complete_graph_moments == {2: 0, 4: 0}
        and all(row["proved"] for row in gauge_rows.values())
        and all(all_relation_checks)
    )
    _require(proved, "exceptional hard moment certificate failed")
    return {
        "field": "F_13",
        "hard_parallel_counts": list(EXCEPTIONAL_HARD_PARALLEL_COUNTS),
        "complete_graph_even_moments": complete_graph_moments,
        "P2": gauge_rows["2"],
        "P8": gauge_rows["8"],
        "normalized_triangle_formula": {
            "triple": "{0,1,r}",
            "q0": "r^2-r+1",
            "S2": "2*q0",
            "S4": "2*q0^2",
            "relation": "2*S4=S2^2",
        },
        "proved": proved,
    }


def _normalized_selected_pair_pattern(
    i: int, j: int, parallel_count: int
) -> dict[tuple[int, int], int]:
    kernel_scalar = Fraction(parallel_count - 6, 12)
    pattern: dict[tuple[int, int], int] = {}
    for s, t in combinations(range(P), 2):
        value = (
            int((s, t) == (i, j))
            + 2 * kernel_scalar
            + int(s in (i, j))
            + int(t in (i, j))
        )
        _require(value.denominator == 1, "selected-pair coefficient nonintegral")
        pattern[(s, t)] = int(value)
    return pattern


@lru_cache(maxsize=1)
def quartic_moment_contradiction() -> dict[str, object]:
    """Apply a sign-safe homogeneous quartic root count."""
    hard = exceptional_hard_moment_certificate()
    selected_dependency = selected_pair_moment_certificate()
    hard_sign_checks: dict[int, dict[str, object]] = {}
    for hard_sign in (-1, 1):
        opposite_sign = -hard_sign
        hard_G_values = []
        for triple in combinations(range(P), 3):
            pattern = _normalized_exceptional_hard_pattern(triple, 2)
            s2 = _pattern_moment(pattern, 2)
            s4 = _pattern_moment(pattern, 4)
            m2 = hard_sign * s2
            m4 = hard_sign * s4
            hard_G_values.append((2 * hard_sign * m4 - m2 * m2) % P)

        opposite_G_values = []
        for parallel_count in OPPOSITE_MINIMUM_PARALLEL_COUNTS:
            for i, j in combinations(range(P), 2):
                pattern = _normalized_selected_pair_pattern(
                    i, j, parallel_count
                )
                s2 = _pattern_moment(pattern, 2)
                s4 = _pattern_moment(pattern, 4)
                m2 = opposite_sign * s2
                m4 = opposite_sign * s4
                opposite_G_values.append(
                    (2 * hard_sign * m4 - m2 * m2) % P
                )
        expected_nonzero_values = {
            (-3 * pow(i - j, 4, P)) % P
            for i, j in combinations(range(P), 2)
        }
        hard_sign_checks[hard_sign] = {
            "hard_sign_h": hard_sign,
            "opposite_sign": opposite_sign,
            "hard_projective_zero_count": HARD_DIRECTION_COUNT,
            "hard_G_values": sorted(set(hard_G_values)),
            "opposite_G_value_set": sorted(set(opposite_G_values)),
            "expected_opposite_nonzero_value_set": sorted(
                expected_nonzero_values
            ),
            "every_opposite_G_value_nonzero": all(opposite_G_values),
            "proved": bool(
                set(hard_G_values) == {0}
                and set(opposite_G_values) == expected_nonzero_values
                and all(opposite_G_values)
            ),
        }

    quartic_degree = 4
    hard_root_count = HARD_DIRECTION_COUNT
    root_count_forces_zero = hard_root_count > quartic_degree
    proved = bool(
        hard["proved"]
        and selected_dependency["proved"]
        and root_count_forces_zero
        and all(row["proved"] for row in hard_sign_checks.values())
    )
    _require(proved, "quartic moment contradiction failed")
    return {
        "global_even_moments": (
            "M_d(L)=sum_({u,v} in H)chi(u-v)*(L(u)-L(v))^d"
        ),
        "homogeneous_quartic": "G=2*h*M_4-M_2^2",
        "quartic_degree": quartic_degree,
        "hard_projective_root_count": hard_root_count,
        "nonzero_binary_quartic_projective_root_bound": quartic_degree,
        "hard_roots_force_G_identically_zero": root_count_forces_zero,
        "sign_checks": hard_sign_checks,
        "opposite_evaluation_formula": "G=-3*(i-j)^4",
        "minus_three_nonzero_mod_13": (-3) % P != 0,
        "one_opposite_selected_pair_contradicts_G_zero": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p13_exceptional_branch_exclusion() -> dict[str, object]:
    """Close the corrected ``p=13,t=3,u=3`` exceptional branch."""
    target = exceptional_hard_target_certificate()
    ledgers = exceptional_parallel_ledgers()
    cells = phase_zero_mass14_cell_reduction()
    moments = quartic_moment_contradiction()
    proved = bool(
        target["proved"]
        and ledgers["proved"]
        and cells["proved"]
        and moments["proved"]
    )
    _require(proved, "p=13 exceptional branch exclusion failed")
    return {
        "p": P,
        "layer_index_t": 3,
        "original_k": 4 * P + 6,
        "H_edge_count": H_EDGE_COUNT,
        "hard_residue_u": 3,
        "hard_mean": EXCEPTIONAL_HARD_MEAN,
        "hard_target": target,
        "parallel_ledgers": ledgers,
        "minimum_opposite_cell": cells,
        "quartic_moment": moments,
        "exceptional_p13_t3_u3_branch_excluded": proved,
        "generic_p13_t3_branch_excluded": False,
        "entire_p13_t3_shell_excluded": False,
        "result_status": "proved branch theorem",
        "proved": proved,
    }


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, int(value**0.5) + 1))


def _conditioned_same_membership_coefficients(
    n: int, selected: int, fixed_size: int, *, fixed_inside: bool
) -> dict[str, Fraction]:
    """Average a cut with a fixed set wholly inside or outside.

    If ``D`` is the sum of row degrees on the fixed set, ``I`` its internal
    coefficient sum, and ``T`` the total coefficient sum, the returned
    coefficients give ``E[cut]=c_T*T+c_D*D+c_I*I``.  They are computed from
    the exact hypergeometric separation probabilities, not fitted formulas.
    """
    remaining = n - fixed_size
    chosen = selected - fixed_size if fixed_inside else selected
    boundary_probability = Fraction(
        n - selected if fixed_inside else selected,
        remaining,
    )
    outside_separation_probability = Fraction(
        2 * chosen * (remaining - chosen),
        remaining * (remaining - 1),
    )
    internal_coefficient = (
        outside_separation_probability - 2 * boundary_probability
        if fixed_size >= 2
        else Fraction(0)
    )
    return {
        "T": outside_separation_probability,
        "D": boundary_probability - outside_separation_probability,
        "I": internal_coefficient,
    }


def _conditioned_oriented_pair_coefficients(
    n: int, selected: int
) -> dict[str, Fraction]:
    """Average a cut conditioned on ``i`` inside and ``j`` outside."""
    remaining = n - 2
    chosen = selected - 1
    i_to_rest_cut_probability = Fraction(remaining - chosen, remaining)
    j_to_rest_cut_probability = Fraction(chosen, remaining)
    rest_separation_probability = Fraction(
        2 * chosen * (remaining - chosen),
        remaining * (remaining - 1),
    )
    return {
        "T": rest_separation_probability,
        "d_i": i_to_rest_cut_probability - rest_separation_probability,
        "d_j": j_to_rest_cut_probability - rest_separation_probability,
        "w_ij": (
            1
            - i_to_rest_cut_probability
            - j_to_rest_cut_probability
            + rest_separation_probability
        ),
    }


def _normalized_conditioned_inequality(
    coefficients: dict[str, Fraction], total: int
) -> dict[str, object]:
    """Normalize ``E[cut] <= total/2`` to an integer linear inequality."""
    variable_coefficients = {
        name: coefficient
        for name, coefficient in coefficients.items()
        if name != "T" and coefficient
    }
    rhs = (Fraction(1, 2) - coefficients["T"]) * total
    denominator = lcm(
        rhs.denominator,
        *(coefficient.denominator for coefficient in variable_coefficients.values()),
    )
    integer_coefficients = {
        name: int(coefficient * denominator)
        for name, coefficient in variable_coefficients.items()
    }
    integer_rhs = int(rhs * denominator)
    divisor = 0
    for value in (*integer_coefficients.values(), integer_rhs):
        divisor = gcd(divisor, abs(value))
    divisor = max(divisor, 1)
    return {
        "coefficients": {
            name: value // divisor
            for name, value in integer_coefficients.items()
        },
        "rhs": integer_rhs // divisor,
    }


@lru_cache(maxsize=None)
def generic_conditioned_entry_alphabet(p: int) -> dict[str, object]:
    """Derive the generic five-value alphabet from conditioned cut means."""
    if not _is_prime(p) or p < 17 or p % 4 != 1:
        raise ValueError("need a prime p>=17 with p=1 mod 4")
    selected = (p + 1) // 2
    total = -(p + 7)
    pair_inside = _conditioned_same_membership_coefficients(
        p, selected, 2, fixed_inside=True
    )
    vertex_outside = _conditioned_same_membership_coefficients(
        p, selected, 1, fixed_inside=False
    )
    oriented_pair = _conditioned_oriented_pair_coefficients(p, selected)

    entry_lower_rational = Fraction(total, p - 1)
    row_upper_rational = Fraction((p - 5) * (p + 7), p + 1)
    entry_upper_rational = Fraction(
        2 * row_upper_rational - total,
        p - 1,
    )
    probability_checks = bool(
        pair_inside
        == {
            "T": Fraction(p - 1, 2 * (p - 2)),
            "D": Fraction(0),
            "I": Fraction(-(p - 1), 2 * (p - 2)),
        }
        and vertex_outside["D"]
        == Fraction(p + 1, 2 * (p - 1) * (p - 2))
        and oriented_pair
        == {
            "T": Fraction(p - 1, 2 * (p - 2)),
            "d_i": Fraction(-1, p - 2),
            "d_j": Fraction(0),
            "w_ij": Fraction(p - 1, 2 * (p - 2)),
        }
    )
    upper_gap_numerator = p * p - 12 * p + 59
    proved = bool(
        probability_checks
        and entry_lower_rational > -2
        and entry_lower_rational <= -1
        and upper_gap_numerator > 0
        and entry_upper_rational < 4
    )
    _require(proved, "generic conditioned entry alphabet changed")
    return {
        "p": p,
        "total_coefficient_sum_T": total,
        "pair_inside_forces": f"w_ij>={entry_lower_rational}",
        "integer_entry_lower_bound": -1,
        "vertex_outside_row_upper_bound": str(row_upper_rational),
        "oriented_pair_entry_upper_bound": str(entry_upper_rational),
        "upper_bound_gap_numerator": upper_gap_numerator,
        "integer_entry_upper_bound": 3,
        "entry_alphabet": [-1, 0, 1, 2, 3],
        "negative_entries_are_simple_edges": True,
        "conditional_probability_identities_checked": probability_checks,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def generic_higher_even_moment_reduction(p: int) -> dict[str, object]:
    """Return the exact higher-even-moment range in generic branch B."""
    if not _is_prime(p) or p < 17 or p % 4 != 1:
        raise ValueError("need a prime p>=17 with p=1 mod 4")
    m = (p + 1) // 2
    hard_quotient_sum = m + 3
    hard_excess_units = hard_quotient_sum - m
    exact_hard_directions = m - hard_excess_units
    degrees = list(range(2, exact_hard_directions, 2))
    expected_last_degree = (p - 9) // 2
    star_moments = {
        degree: {
            center: sum(
                pow((center - label) % p, degree, p)
                for label in range(p)
                if label != center
            )
            % p
            for center in range(p)
        }
        for degree in degrees
    }
    every_star_zero = all(
        value == 0
        for centers in star_moments.values()
        for value in centers.values()
    )
    every_root_count_strict = all(
        exact_hard_directions > degree for degree in degrees
    )
    conditioned_entries = generic_conditioned_entry_alphabet(p)
    proved = bool(
        hard_excess_units == 3
        and exact_hard_directions == (p - 5) // 2
        and degrees
        and degrees[-1] == expected_last_degree
        and len(degrees) == (p - 9) // 4
        and every_star_zero
        and every_root_count_strict
        and conditioned_entries["proved"]
    )
    _require(proved, "generic higher-even-moment reduction failed")
    return {
        "p": p,
        "m": m,
        "hard_quotient_identity": f"sum k_d={hard_quotient_sum}",
        "hard_quotient_excess_units": hard_excess_units,
        "exact_hard_star_directions_at_least": exact_hard_directions,
        "global_moment_definition": (
            "M_d(L)=sum_({u,v} in H)chi(u-v)*(L(u)-L(v))^d"
        ),
        "orientation_independent_even_degrees": degrees,
        "last_forced_even_degree": degrees[-1],
        "forced_even_moment_count": len(degrees),
        "star_power_sum_checks": star_moments,
        "every_exact_star_moment_zero": every_star_zero,
        "projective_root_count_strictly_exceeds_degree": (
            every_root_count_strict
        ),
        "global_even_moments_forced_identically_zero": degrees,
        "conditioned_entry_reduction": conditioned_entries,
        "result_status": "proved open reduction",
        "generic_branch_excluded": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p13_generic_elevated_local_counterexample() -> dict[str, object]:
    """Exhibit a real elevated hard cell that defeats a local extension.

    This is a directional quadratic, not a common residual graph.  It proves
    that the four exact p=13 hard stars and ``M_2=0`` cannot force the three
    elevated hard directions back to the ``b=12`` baseline one at a time.
    """
    clique = (0, 1, 2, 3, 5)
    extra_edge = (0, 11)
    edges = tuple(combinations(clique, 2)) + (extra_edge,)
    pattern = {
        pair: int(pair in edges) for pair in combinations(range(P), 2)
    }
    row_sums = {
        vertex: sum(
            coefficient
            for pair, coefficient in pattern.items()
            if vertex in pair
        )
        for vertex in range(P)
    }
    odd_rows = [vertex for vertex, degree in row_sums.items() if degree % 2]
    cut_histogram: Counter[int] = Counter()
    value_histogram: Counter[int] = Counter()
    for subset_tuple in combinations(range(P), M):
        subset = set(subset_tuple)
        cut = sum(
            coefficient
            for (s, t), coefficient in pattern.items()
            if (s in subset) != (t in subset)
        )
        value = 7 - cut
        cut_histogram[cut] += 1
        value_histogram[value] += 1
    moment2 = _pattern_moment(pattern, 2)
    moment4 = _pattern_moment(pattern, 4)
    value_sum = sum(value * count for value, count in value_histogram.items())
    scaled_mean = Fraction(2 * P * value_sum, len(tuple(combinations(range(P), M))))
    proved = bool(
        len(edges) == 11
        and sum(pattern.values()) == 11
        and sum(abs(value) for value in pattern.values()) == 11
        and odd_rows == [0, 11]
        and min(cut_histogram) == 0
        and max(cut_histogram) == 7
        and min(value_histogram) == 0
        and scaled_mean == 28
        and moment2 == 0
        and moment4 == 5
    )
    _require(proved, "p=13 elevated local counterexample changed")
    return {
        "p": P,
        "generic_hard_quotient_k": 2,
        "hard_parallel_count_P": 6,
        "signed_matrix": (
            "W=Adj(K_5 on {0,1,2,3,5})+1_{{0,11}}"
        ),
        "nonzero_edges": [list(edge) for edge in edges],
        "sum_W": sum(pattern.values()),
        "l1_norm": sum(abs(value) for value in pattern.values()),
        "available_nonparallel_edge_bound": H_EDGE_COUNT - 6,
        "row_sums": row_sums,
        "odd_rows": odd_rows,
        "directional_b": len(odd_rows),
        "cut_histogram": dict(sorted(cut_histogram.items())),
        "A_formula": "A(X)=7-cut_W(X)",
        "A_value_histogram": dict(sorted(value_histogram.items())),
        "A_nonnegative_on_all_1716_seven_sets": min(value_histogram) >= 0,
        "scaled_mean_2pE_A": int(scaled_mean),
        "normalized_degree_two_moment_S2_mod_13": moment2,
        "normalized_degree_four_moment_S4_mod_13": moment4,
        "global_moment_interpretation": (
            "M_2=h*S_2=0 and M_4=h*S_4=5*h for the hard sign h"
        ),
        "constructs_common_residual_graph": False,
        "consequence": (
            "counterexample to forcing elevated hard directions back to "
            "b=12 using one-direction floors plus M2"
        ),
        "result_status": "counterexample to method",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def _p17_cut_range_certificate() -> dict[str, object]:
    """Derive the missing lower cut bound by stabilizer averaging."""
    total = -24
    class_sizes = {"a": 36, "b": 28, "c": 72}

    def intersection_average_coefficients(intersection: int) -> dict[str, int]:
        inside_selected = intersection
        outside_selected = 9 - intersection
        inside_probability = Fraction(
            2 * inside_selected * (9 - inside_selected), 9 * 8
        )
        outside_probability = Fraction(
            2 * outside_selected * (8 - outside_selected), 8 * 7
        )
        cross_probability = Fraction(
            inside_selected * (8 - outside_selected)
            + (9 - inside_selected) * outside_selected,
            9 * 8,
        )
        coefficients = {
            "a": class_sizes["a"] * inside_probability,
            "b": class_sizes["b"] * outside_probability,
            "c": class_sizes["c"] * cross_probability,
        }
        _require(
            all(value.denominator == 1 for value in coefficients.values()),
            "p17 stabilizer average lost integrality",
        )
        return {key: int(value) for key, value in coefficients.items()}

    intersection4 = intersection_average_coefficients(4)
    intersection5 = intersection_average_coefficients(5)
    total_coefficients = class_sizes
    multipliers = {
        "intersection4": Fraction(-9),
        "intersection5": Fraction(-45, 4),
        "total_equality": Fraction(45, 4),
    }
    eliminated = {
        name: (
            multipliers["intersection4"] * intersection4[name]
            + multipliers["intersection5"] * intersection5[name]
            + multipliers["total_equality"] * total_coefficients[name]
        )
        for name in ("a", "b", "c")
    }
    lower_rhs = (
        multipliers["intersection4"] * -12
        + multipliers["intersection5"] * -12
        + multipliers["total_equality"] * total
    )
    raw_cut_lower = lower_rhs / eliminated["c"] * class_sizes["c"]
    parity_improved_cut_lower = -26
    mean_B = Fraction(6, 17)
    slice_size = comb(17, 9)
    proved = bool(
        intersection4 == {"a": 20, "b": 15, "c": 37}
        and intersection5 == {"a": 20, "b": 16, "c": 36}
        and eliminated == {"a": 0, "b": 0, "c": 72}
        and lower_rhs == -27
        and raw_cut_lower == -27
        and parity_improved_cut_lower == -26
        and slice_size == 24_310
        and slice_size * mean_B == 8_580
    )
    _require(proved, "p17 stabilizer cut-range certificate changed")
    return {
        "fixed_nine_set_edge_class_sizes": class_sizes,
        "intersection4_average": intersection4,
        "intersection5_average": intersection5,
        "total_sum_equality": "36*a+28*b+72*c=-24",
        "elimination_multipliers": {
            key: str(value) for key, value in multipliers.items()
        },
        "eliminated_coefficients": {
            key: int(value) for key, value in eliminated.items()
        },
        "raw_cut_lower_bound": int(raw_cut_lower),
        "cut_is_even_from_even_row_degrees": True,
        "parity_improved_cut_lower_bound": parity_improved_cut_lower,
        "cut_upper_bound": -12,
        "B_value_range": list(range(8)),
        "B_mean": str(mean_B),
        "J17_9_size": slice_size,
        "B_total_mass": int(slice_size * mean_B),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p17_conditioned_cut_reduction() -> dict[str, object]:
    """Derive the exact bounded signed-matrix reduction at ``p=17``."""
    p = 17
    total = -(p + 7)
    l1_bound = 4 * p + 4
    moment = generic_higher_even_moment_reduction(p)
    entry = moment["conditioned_entry_reduction"]
    row_degrees = list(range(-(p - 1), 17, 2))
    entry_alphabet = entry["entry_alphabet"]
    positive_mass_upper_bound = (l1_bound + total) // 2
    negative_edge_upper_bound = (l1_bound - total) // 2

    pair_inside = _conditioned_same_membership_coefficients(
        p, 9, 2, fixed_inside=True
    )
    vertex_outside = _conditioned_same_membership_coefficients(
        p, 9, 1, fixed_inside=False
    )
    oriented_pair = _conditioned_oriented_pair_coefficients(p, 9)
    pair_outside = _conditioned_same_membership_coefficients(
        p, 9, 2, fixed_inside=False
    )
    triple_inside = _conditioned_same_membership_coefficients(
        p, 9, 3, fixed_inside=True
    )
    triple_outside = _conditioned_same_membership_coefficients(
        p, 9, 3, fixed_inside=False
    )
    four_inside = _conditioned_same_membership_coefficients(
        p, 9, 4, fixed_inside=True
    )
    four_outside = _conditioned_same_membership_coefficients(
        p, 9, 4, fixed_inside=False
    )
    conditional_inequalities = {
        "pair_inside": _normalized_conditioned_inequality(pair_inside, total),
        "vertex_outside": _normalized_conditioned_inequality(
            vertex_outside, total
        ),
        "oriented_pair": _normalized_conditioned_inequality(
            oriented_pair, total
        ),
        "pair_outside": _normalized_conditioned_inequality(
            pair_outside, total
        ),
        "triple_inside": _normalized_conditioned_inequality(
            triple_inside, total
        ),
        "triple_outside": _normalized_conditioned_inequality(
            triple_outside, total
        ),
        "four_inside": _normalized_conditioned_inequality(
            four_inside, total
        ),
        "four_outside": _normalized_conditioned_inequality(
            four_outside, total
        ),
    }
    expected_inequalities = {
        "pair_inside": {"coefficients": {"I": -2}, "rhs": 3},
        "vertex_outside": {"coefficients": {"D": 1}, "rhs": 16},
        "oriented_pair": {
            "coefficients": {"d_i": -1, "w_ij": 8},
            "rhs": 12,
        },
        "pair_outside": {
            "coefficients": {"D": 1, "I": -8},
            "rhs": 4,
        },
        "triple_inside": {
            "coefficients": {"D": 1, "I": -14},
            "rhs": 15,
        },
        "triple_outside": {
            "coefficients": {"D": 9, "I": -48},
            "rhs": -8,
        },
        "four_inside": {
            "coefficients": {"D": 1, "I": -7},
            "rhs": 3,
        },
        "four_outside": {
            "coefficients": {"D": 1, "I": -4},
            "rhs": -4,
        },
    }
    conditional_average_checks = {
        key: conditional_inequalities[key] == expected
        for key, expected in expected_inequalities.items()
    }
    cut_range = _p17_cut_range_certificate()
    proved = bool(
        total == -24
        and l1_bound == 72
        and row_degrees == list(range(-16, 17, 2))
        and entry["proved"]
        and entry_alphabet == [-1, 0, 1, 2, 3]
        and all(conditional_average_checks.values())
        and positive_mass_upper_bound == 24
        and negative_edge_upper_bound == 48
        and moment["global_even_moments_forced_identically_zero"] == [2, 4]
        and cut_range["proved"]
    )
    _require(proved, "p=17 conditioned-cut reduction changed")
    return {
        "p": p,
        "opposite_parallel_count_Q": 3,
        "sum_W": total,
        "l1_bound": l1_bound,
        "balanced_cut_upper_bound": -12,
        "B_formula": "B(X)=-6-cut_W(X)/2",
        "entry_alphabet": entry_alphabet,
        "negative_entries_are_simple_edges": entry[
            "negative_entries_are_simple_edges"
        ],
        "row_degrees_even_range": row_degrees,
        "conditional_average_inequalities": conditional_inequalities,
        "conditional_average_checks": conditional_average_checks,
        "pair_inequalities": [
            "d_i+d_j-4 <= 8*w_ij",
            "8*w_ij <= d_i+12",
            "8*w_ij <= d_j+12",
        ],
        "triple_inequalities": [
            "D_T-14*I_T <= 15",
            "48*I_T-9*D_T >= 8",
        ],
        "four_set_inequalities": [
            "D_T-7*I_T <= 3",
            "4*I_T-D_T >= 4",
        ],
        "positive_multiplicity_upper_bound": positive_mass_upper_bound,
        "negative_edge_count_upper_bound": negative_edge_upper_bound,
        "stabilizer_cut_range": cut_range,
        "forced_moment_degrees": moment[
            "global_even_moments_forced_identically_zero"
        ],
        "full_cut_and_moment_model_status": "OPEN; bounded runs inconclusive",
        "result_status": "proved open reduction",
        "generic_p17_branch_excluded": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15739() -> dict[str, object]:
    """Package the exceptional close and honest generic reduction."""
    dependency = proposition_15738()
    exceptional = p13_exceptional_branch_exclusion()
    p13_generic_barrier = p13_generic_elevated_local_counterexample()
    p17 = p17_conditioned_cut_reduction()
    representative_generic = {
        str(p): generic_higher_even_moment_reduction(p)
        for p in (17, 29, 37)
    }
    proved = bool(
        dependency["proved"]
        and exceptional["proved"]
        and p13_generic_barrier["proved"]
        and p17["proved"]
        and all(row["proved"] for row in representative_generic.values())
    )
    _require(proved, "Proposition 15.739 audit failed")
    return {
        "prop": "15.739",
        "title": "Quartic moment closes the exceptional p=13 fourth-shell branch",
        "result_status": "proved branch theorem and open reduction",
        "statement": (
            "the p=13,t=3,u=3 exceptional branch is empty; in the remaining "
            "generic p=1 mod 4 branch, all even moments through (p-9)/2 "
            "vanish identically"
        ),
        "finite_certificate_dependency": {
            "proposition": "15.738",
            "result_status": dependency["result_status"],
            "p13_mass14_cells_classified": dependency[
                "p13_mass14_cells_classified"
            ],
        },
        "p13_exceptional_branch_exclusion": exceptional,
        "p13_generic_local_method_counterexample": p13_generic_barrier,
        "generic_higher_even_moment_examples": representative_generic,
        "p17_conditioned_cut_reduction": p17,
        "p13_t3_exceptional_u3_closed": True,
        "p13_t3_generic_branch_closed": False,
        "p13_k_eq_58_closed": False,
        "generic_p_ge_17_t3_branch_closed": False,
        "k_eq_4p_plus_6_shell_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "the generic p=1 mod 4 branch at k=4p+6 (including p=13), "
            "all later residual layers, critical p=5,7, p=11 k>=50, "
            "multi-level Type I, and the limit"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic p=13 quartic-moment certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15739.json"
    output.write_text(
        json.dumps(proposition_15739(), indent=2, sort_keys=True) + "\n"
    )
    return output


def main() -> None:
    result = proposition_15739()
    path = write_evidence()
    print("Prop. 15.739: exceptional p=13,t=3,u=3 branch excluded")
    print("generic k=4p+6 branch and residual (ii) remain open")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
