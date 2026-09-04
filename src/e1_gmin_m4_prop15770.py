#!/usr/bin/env python3
r"""Prop. 15.770 -- the next residual layers after 15.768--15.769.

Put ``q=(p-1)/2`` and ``m=q+1``.  Proposition 15.768 closes
``t=q-3`` for ``p=1 (mod 4)`` and Proposition 15.769 closes ``t=q-2``
for ``p=3 (mod 4)``.  The changed premise here is that those newly proved
equality branches can be advanced by one quotient unit without changing the
common signed row sum.

For ``p=1 (mod 4)``, the next layer is ``t=q-2``.  The complement-triple
branch becomes ``m-1`` low rows of parallel count two and one high row of
parallel count three; it still forces local mass ``p+15``.  A new all-low
XNOR branch is a sharp ``p-3`` Boolean lift.  The fixed four-bit catalog
leaves exactly the omitted-pair and all-equal-triple families, of coefficient
offsets three and five, and both force local mass ``p+13``.  At ``p=29`` the
height-at-least-two endpoint of that local theorem is closed by Proposition
15.768's sharp mean-three-quarters cube theorem.  Hence this layer is empty
for every prime ``p=1 (mod 4), p>=29``.

For ``p=3 (mod 4)``, the next layer is ``t=q-1``.  The four sharp families
from Proposition 15.769 become ``m-1`` low rows of parallel count ``P`` and
one high row of count ``P+1``; their signed row sum is unchanged and they
again force local mass ``p+13``.  The only additional residue is a lift of
scaled mass ``p-1``.  Proposition 15.688 makes it Boolean, and the corrected
Johnson/cube influence reduction sends it to the fixed four-bit catalog,
which misses density ``(p-1)/(4p)``.  Thus this layer is empty for every
prime ``p=3 (mod 4), p>=31``.  At the exceptional prime ``p=23``, the
mass-``p-1`` branch is excluded by the same fixed catalog, while the carried
branch has eleven low hard roots.  Eleven still exceeds degrees four and
eight, so the equality/moment certificate of Proposition 15.769 closes
``p=23,t=10,k=112`` as well.

These are two infinite one-layer extensions.  They do not close residual
(ii) globally.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15751 import (
    cube_half_mean_height_certificate,
    density_profile_certificate,
    exact_four_cube_catalog,
    height_at_least_two_certificate,
    height_one_junta_certificate,
    profile_density,
)
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion
from e1_gmin_m4_prop15768 import (
    complement_triple_baseline_certificate,
    cube_three_quarter_height_certificate,
    p_plus_fifteen_local_exclusion,
)
from e1_gmin_m4_prop15769 import (
    hard_family_catalog as p3_hard_family_catalog,
    p_plus_thirteen_local_exclusion as p3_p_plus_thirteen_local_exclusion,
    sharp_p_minus_three_boolean_classification as p3_sharp_classification,
    sharp_p_minus_three_four_bit_catalog,
)
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
P1_CARRIED_BRANCH = "carried_complement_triple"
P1_NEW_SHARP_BRANCH = "all_low_XNOR_sharp_p_minus_3"
P3_CARRIED_BRANCH = "carried_sharp_p_minus_3"
P3_NEW_LOCAL_BRANCH = "all_low_p_minus_one_lift"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int, congruence: int) -> None:
    threshold = 29 if congruence == 1 else 31
    if (
        congruence not in (1, 3)
        or not isinstance(p, int)
        or isinstance(p, bool)
        or p < threshold
        or p % 4 != congruence
        or not is_prime(p)
    ):
        raise ValueError(
            f"need a prime p>={threshold} congruent to {congruence} modulo 4"
        )


def _four_bit_densities(p: int) -> list[Fraction]:
    return sorted(
        {
            profile_density(tuple(row["layer_counts"]), p)
            for row in exact_four_cube_catalog()["profiles"]
        }
    )


def sharp_p_minus_three_classification_all_odd(p: int) -> dict[str, object]:
    """Extend 15.769's sharp ``p-3`` classification to both congruences."""
    congruence = p % 4
    _check_prime(p, congruence)
    q = (p - 1) // 2
    density = Fraction(p - 3, 4 * p)
    influence_floor = Fraction(
        (p + 1) * (p - 3), 16 * p * (p - 2)
    )
    total_influence_upper = (p - 1) * density * (1 - density)
    junta_bound = Fraction(6 * (p - 1) * (p - 2), p * p)
    catalog = sharp_p_minus_three_four_bit_catalog()
    matching_profiles = sorted(
        tuple(int(value) for value in row["layer_counts"])
        for row in exact_four_cube_catalog()["profiles"]
        if profile_density(tuple(row["layer_counts"]), p) == density
    )
    expected_profiles = [(0, 0, 1, 2, 1), (1, 1, 0, 1, 1)]

    prior_agreement = True
    if congruence == 3:
        prior = p3_sharp_classification(p)
        prior_agreement = bool(
            prior["proved"]
            and prior["matching_four_bit_layer_profiles"]
            == [list(row) for row in expected_profiles]
            and prior["four_bit_classification"]["selected_tables_sha256"]
            == catalog["selected_tables_sha256"]
        )

    proved = bool(
        sharp_integral_quadratic_lift_floor(p)["proved"]
        and junta_bound < 6
        and 5 < q
        and matching_profiles == expected_profiles
        and catalog["proved"]
        and catalog["selected_pair_table_count"] == 6
        and catalog["all_equal_triple_table_count"] == 4
        and prior_agreement
    )
    _require(proved, "the all-congruence sharp p-3 classification failed")
    return {
        "p": p,
        "p_mod_4": congruence,
        "slice": f"J({p},{(p + 1) // 2})",
        "density": str(density),
        "H_at_least_two_scaled_floor": int(
            sharp_integral_quadratic_lift_floor(p)[
                "H_at_least_two_scaled_floor"
            ]
        ),
        "therefore_height_one_boolean_at_mass_p_minus_three": True,
        "relevant_pair_influence_floor": str(influence_floor),
        "total_influence_upper_bound": str(total_influence_upper),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "junta_coordinates_at_most": 5,
        "five_less_than_both_slice_sides": 5 < q,
        "cube_active_coordinates_at_most": 4,
        "matching_four_bit_layer_profiles": [
            list(row) for row in matching_profiles
        ],
        "selected_pair_table_count": catalog["selected_pair_table_count"],
        "all_equal_triple_table_count": catalog[
            "all_equal_triple_table_count"
        ],
        "selected_tables_sha256": catalog["selected_tables_sha256"],
        "original_slice_families": ["omitted_pair", "all_equal_triple"],
        "agrees_with_prop_15769_when_p_is_3_mod_4": prior_agreement,
        "proved": proved,
    }


def p1_sharp_family_catalog(p: int) -> dict[str, object]:
    """Classify the new all-low XNOR sharp lifts for ``p=1 (mod 4)``."""
    _check_prime(p, 1)
    baseline = parity_floor_certificate(p, 2, 1)
    classification = sharp_p_minus_three_classification_all_odd(p)
    families = [
        {
            "baseline": "XNOR",
            "lift": "omitted_pair",
            "baseline_offset": 4,
            "lift_offset_increment": -1,
            "coefficient_offset": 3,
        },
        {
            "baseline": "XNOR",
            "lift": "all_equal_triple",
            "baseline_offset": 4,
            "lift_offset_increment": 1,
            "coefficient_offset": 5,
        },
    ]
    q = (p - 1) // 2
    proved = bool(
        baseline["exact_positive_quadrature_certificate"]
        and int(baseline["scaled_floor"]) == p - 1
        and all(weight > 0 for weight in baseline["quadrature_weights"])
        and classification["proved"]
        and [row["coefficient_offset"] for row in families] == [3, 5]
    )
    _require(proved, "the p=1 sharp-lift family catalog failed")
    return {
        "p": p,
        "baseline": "A_0=(1-x_i-x_j)^2",
        "baseline_signed_target": "3+2A_0=4+z_i*z_j",
        "baseline_scaled_mean": p - 1,
        "difference_lift": "B=(A-A_0)/2",
        "difference_scaled_mass": p - 3,
        "sharp_classification": classification,
        "lift_signed_targets": {
            "omitted_pair": "4B=1-z_i-z_j+z_i*z_j",
            "all_equal_triple": "4B=1+z_i*z_j+z_i*z_k+z_j*z_k",
        },
        "coefficient_congruence_modulus": q,
        "families": families,
        "proved": proved,
    }


def p1_p_plus_thirteen_local_exclusion(p: int) -> dict[str, object]:
    """Exclude ``4p E[C]=p+13`` for ``p=1 (mod 4), p>=29``."""
    _check_prime(p, 1)
    mass = p + 13
    half_mean = cube_half_mean_height_certificate()
    three_quarter = cube_three_quarter_height_certificate()
    lower_height = Fraction(p - 11, 4)
    stabilizer_upper = Fraction(mass * (p + 3), 4 * (p - 1))
    paired_average_upper = Fraction(mass, 2 * (p - 1))

    if p == 29:
        half_mean_integer_height = 5
        refined_three_quarter_height = 12
        forced_height = 12
        paired_average = Fraction(4 * forced_height + mass, 4 * (p + 1))
        height_row = {
            "endpoint": True,
            "initial_raw_height_lower_bound": str(lower_height),
            "initial_integral_height_lower_bound": half_mean_integer_height,
            "half_mean_cube_maximum_upper_bound": 3,
            "refined_three_quarter_raw_height_lower_bound": str(
                Fraction(3 * (p + 1) - mass, 4)
            ),
            "refined_integral_height_lower_bound": refined_three_quarter_height,
            "stabilizer_height_upper_bound": str(stabilizer_upper),
            "forced_height": forced_height,
            "paired_cube_average_at_forced_height": str(paired_average),
            "every_paired_cube_has_mean_three_quarters": True,
            "three_quarter_cube_maximum_upper_bound": three_quarter[
                "maximum_upper_bound"
            ],
            "contradiction": True,
        }
        height_proved = bool(
            half_mean["proved"]
            and three_quarter["proved"]
            and lower_height == Fraction(9, 2)
            and half_mean_integer_height == 5 > 3
            and Fraction(3 * (p + 1) - mass, 4) == 12
            and stabilizer_upper == 12
            and paired_average == Fraction(3, 4)
            and int(three_quarter["maximum_upper_bound"]) < forced_height
        )
    else:
        height_row = {
            "endpoint": False,
            "height_lower_bound": str(lower_height),
            "stabilizer_height_upper_bound": str(stabilizer_upper),
            "paired_cube_average_upper_bound": str(paired_average_upper),
            "some_paired_cube_has_mean_exactly_one_half": True,
            "half_mean_cube_maximum_upper_bound": 3,
            "contradiction": True,
        }
        height_proved = bool(
            half_mean["proved"]
            and lower_height > 3
            and paired_average_upper < Fraction(3, 4)
        )

    density = Fraction(mass, 4 * p)
    influence_floor = Fraction(
        (p + 1) * (p - 3), 16 * p * (p - 2)
    )
    total_influence_upper = (p - 1) * density * (1 - density)
    junta_bound = Fraction(
        2 * (p - 1) * (p - 2) * (p + 13) * (3 * p - 13),
        p * p * (p + 1) * (p - 3),
    )
    eight_gap = p**4 - 25 * p**3 + 229 * p**2 - 559 * p + 338
    translated = p - 29
    translated_gap = (
        translated**4
        + 91 * translated**3
        + 3100 * translated**2
        + 47204 * translated
        + 274272
    )
    fixed_catalog = exact_four_cube_catalog()
    densities = _four_bit_densities(p)
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
        eight_gap == translated_gap
        and eight_gap > 0
        and junta_bound < 8
        and 7 < (p - 1) // 2
        and fixed_catalog["proved"]
        and densities == expected_densities
        and Fraction(p + 1, 4 * p) < density < Fraction(p - 1, 2 * p)
        and density not in densities
    )
    proved = height_proved and boolean_proved
    _require(proved, "the p=1 mod 4 p+13 local theorem failed")
    return {
        "p": p,
        "slice": f"J({p},{(p + 1) // 2})",
        "statement": "no nonzero nonnegative integral quadratic has 4p E[C]=p+13",
        "height_at_least_two": {**height_row, "proved": height_proved},
        "height_one_boolean": {
            "density": str(density),
            "relevant_pair_influence_floor": str(influence_floor),
            "total_influence_upper_bound": str(total_influence_upper),
            "largest_zero_influence_class_complement_bound": str(junta_bound),
            "eight_gap_at_p_equals_x_plus_29": [1, 91, 3100, 47204, 274272],
            "junta_coordinates_at_most": 7,
            "cube_active_coordinates_at_most": 4,
            "four_bit_density_values": [str(value) for value in densities],
            "target_density_absent": True,
            "proved": boolean_proved,
        },
        "p29_uses_three_quarter_cube_endpoint": p == 29,
        "finite_prime_or_slice_census_used": False,
        "proved": proved,
    }


def p_minus_one_local_exclusion(p: int) -> dict[str, object]:
    """Exclude a nonzero lift with scaled mass ``p-1`` for p=3 mod 4."""
    _check_prime(p, 3)
    lift = sharp_integral_quadratic_lift_floor(p)
    mass = p - 1
    density = Fraction(mass, 4 * p)
    influence_floor = Fraction(
        (p + 1) * (p - 3), 16 * p * (p - 2)
    )
    total_influence_upper = (p - 1) * density * (1 - density)
    junta_bound = Fraction(
        2 * (p - 1) ** 2 * (p - 2) * (3 * p + 1),
        p * p * (p + 1) * (p - 3),
    )
    six_gap = 5 * p**3 - 20 * p**2 + p + 2
    translated = p - 31
    translated_gap = (
        5 * translated**3
        + 445 * translated**2
        + 13176 * translated
        + 129768
    )
    fixed_catalog = exact_four_cube_catalog()
    densities = _four_bit_densities(p)
    proved = bool(
        lift["proved"]
        and int(lift["H_at_least_two_scaled_floor"]) == p + 1 > mass
        and six_gap == translated_gap
        and six_gap > 0
        and junta_bound < 6
        and 5 < (p - 1) // 2
        and fixed_catalog["proved"]
        and Fraction(p - 3, 4 * p) < density < Fraction(p + 1, 4 * p)
        and density not in densities
    )
    _require(proved, "the p-1 local exclusion failed")
    return {
        "p": p,
        "scaled_mass": mass,
        "H_at_least_two_scaled_floor": p + 1,
        "therefore_height_one_boolean": True,
        "density": str(density),
        "relevant_pair_influence_floor": str(influence_floor),
        "total_influence_upper_bound": str(total_influence_upper),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "six_gap_polynomial": "5p^3-20p^2+p+2",
        "six_gap_at_p_equals_x_plus_31": [5, 445, 13176, 129768],
        "junta_coordinates_at_most": 5,
        "cube_active_coordinates_at_most": 4,
        "four_bit_density_values": [str(value) for value in densities],
        "target_density_absent": True,
        "finite_prime_or_slice_census_used": False,
        "proved": proved,
    }


def _residue_ledger(p: int, congruence: int) -> dict[str, object]:
    """Audit every hard-type residue at the claimed next layer."""
    _check_prime(p, congruence)
    m = (p + 1) // 2
    t = m - (3 if congruence == 1 else 2)
    floors = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(p)[
            "phase_one_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"]
    )
    if congruence == 1:
        expected_live = {
            0: [(p - 1, "exact")],
            m - 4: [(p - 3, "exact")],
            m - 3: [(2, "sharp_p_minus_3")],
            m - 1: [(2, "exact")],
        }
    else:
        expected_live = {
            m - 3: [(2, "sharp_p_minus_3"), (p - 1, "sharp_p_minus_3")],
            m - 2: [(2, "p_minus_one"), (p - 1, "p_minus_one")],
            m - 1: [(2, "exact"), (p - 1, "exact")],
        }

    rows: list[dict[str, object]] = []
    for residue in range(m):
        quotient_sum = m + t - residue
        if quotient_sum >= m:
            low_quotient = 1
            forced_low_count = 2 * m - quotient_sum
            low_mean = p + 1 + 2 * residue
        else:
            low_quotient = 0
            forced_low_count = m - quotient_sum
            low_mean = 2 * residue

        candidates = []
        live = []
        for boundary, floor in floors.items():
            if floor > low_mean:
                continue
            excess = low_mean - floor
            if excess == 0:
                classification = "exact"
            elif excess < lift_floor:
                classification = "excluded_sub_sharp_lift"
            elif excess == lift_floor:
                classification = "sharp_p_minus_3"
            elif congruence == 3 and residue == m - 2 and excess == p - 1:
                classification = "p_minus_one"
            else:
                raise ArithmeticError(
                    f"unclassified phase-one excess at p={p},u={residue},b={boundary}"
                )
            candidates.append(
                {
                    "b": boundary,
                    "floor": floor,
                    "excess": excess,
                    "classification": classification,
                }
            )
            if classification != "excluded_sub_sharp_lift":
                live.append((boundary, classification))

        _require(
            live == expected_live.get(residue, []) and forced_low_count > 0,
            f"the p={p},u={residue} next-layer residue row changed",
        )
        rows.append(
            {
                "u": residue,
                "quotient_sum": quotient_sum,
                "forced_low_quotient": low_quotient,
                "forced_low_direction_count_at_least": forced_low_count,
                "forced_low_mean": low_mean,
                "candidate_rows": candidates,
                "live_rows": [
                    {"b": boundary, "classification": classification}
                    for boundary, classification in live
                ],
                "excluded_arithmetically": not live,
            }
        )

    edge_count = 4 * p + 2 * t + 1
    proved = bool(
        edge_count == (5 * p - 4 if congruence == 1 else 5 * p - 2)
        and p * p + 1 - 2 * edge_count > 0
        and [row["u"] for row in rows] == list(range(m))
    )
    _require(proved, "the next-layer residue ledger failed")
    return {
        "p": p,
        "p_mod_4": congruence,
        "m": m,
        "t": t,
        "k": 4 * p + 2 * t,
        "H_edge_count": edge_count,
        "guaranteed_isolated_vertices": p * p + 1 - 2 * edge_count,
        "phase_one_mean_form": f"a_L=2u+{p + 1}k_L",
        "phase_one_quotient_sum": "sum k_L=m+t-u",
        "sharp_lift_floor": lift_floor,
        "arithmetic_surviving_residues": sorted(expected_live),
        "rows": rows,
        "proved": proved,
    }


def p1_next_residue_ledger(p: int) -> dict[str, object]:
    return _residue_ledger(p, 1)


def p3_next_residue_ledger(p: int) -> dict[str, object]:
    return _residue_ledger(p, 3)


def _opposite_forcing(
    p: int,
    edge_count: int,
    hard_mean: int,
    hard_parallel: int,
    forbidden_Q: int,
    forbidden_mass: int,
    next_Q: int,
    next_mass: int,
    local_dependency: dict[str, object],
) -> dict[str, object]:
    """Apply the common-row identity and opposite-direction pigeonhole."""
    m = (p + 1) // 2
    h_times_T = (p + 1) * hard_parallel - 3 * p - hard_mean
    hard_edges = (edge_count + h_times_T) // 2
    opposite_edges = edge_count - hard_edges
    phase_zero = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(p)[
            "phase_zero_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"]
    )

    def opposite_mean(parallel: int) -> int:
        return (p + 1) * parallel + h_times_T - 3 * p

    nonzero_rows = [
        (boundary, floor, next_mass - floor)
        for boundary, floor in phase_zero.items()
        if boundary != 0 and floor <= next_mass
    ]
    surplus = opposite_edges - m * next_Q
    proved = bool(
        2 * hard_edges == edge_count + h_times_T
        and opposite_mean(forbidden_Q) == forbidden_mass
        and 0 < forbidden_mass < min(
            min(value for boundary, value in phase_zero.items() if boundary != 0),
            lift_floor,
        )
        and opposite_mean(next_Q) == next_mass
        and 0 <= surplus < m
        and [boundary for boundary, _floor, _excess in nonzero_rows]
        == [2, p - 1]
        and all(
            boundary in {2, p - 1} and 0 < excess < lift_floor
            for boundary, _floor, excess in nonzero_rows
        )
        and local_dependency["proved"]
    )
    _require(proved, "the opposite-direction forcing ledger failed")
    return {
        "hard_mean_on_low_rows": hard_mean,
        "low_row_parallel_count": hard_parallel,
        "hard_sign_times_global_T": h_times_T,
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "opposite_mean_formula": f"a(Q)={p + 1}Q+({h_times_T})-{3 * p}",
        "forbidden_Q": forbidden_Q,
        "forbidden_scaled_mean": forbidden_mass,
        "forbidden_below_phase_zero_and_lift_floors": True,
        "forced_next_Q": next_Q,
        "forced_next_scaled_mean": next_mass,
        "surplus_after_every_opposite_Q_at_least_next": surplus,
        "directions_at_next_Q_at_least": m - surplus,
        "nonzero_boundary_floor_and_lift_rows": [
            list(row) for row in nonzero_rows
        ],
        "next_Q_forced_to_boundary_zero": True,
        "local_dependency_proved": True,
        "proved": proved,
    }


def _parallel_candidates(p: int, edge_count: int, offset: int) -> list[int]:
    q = (p - 1) // 2
    m = (p + 1) // 2
    return [
        parallel
        for parallel in range(edge_count // m + 1)
        if (parallel - offset) % q == 0
    ]


def p1_next_layer_exclusion(p: int) -> dict[str, object]:
    """Close ``t=q-2``, equivalently ``k=5p-5``, for p=1 mod 4."""
    _check_prime(p, 1)
    m = (p + 1) // 2
    edge_count = 5 * p - 4
    residues = p1_next_residue_ledger(p)
    p_plus_seven = {
        "proved": all(
            row["proved"]
            for row in (
                height_at_least_two_certificate(p),
                height_one_junta_certificate(p),
                density_profile_certificate(p),
            )
        )
    }
    p_plus_nine = p_plus_nine_local_exclusion(p)
    p_plus_fifteen = p_plus_fifteen_local_exclusion(p)
    p_plus_thirteen = p1_p_plus_thirteen_local_exclusion(p)

    old_literal = _opposite_forcing(
        p, edge_count, p + 1, 5, 2, 6, 3, p + 7, p_plus_seven
    )
    old_XNOR = _opposite_forcing(
        p, edge_count, p - 1, 4, 3, 8, 4, p + 9, p_plus_nine
    )

    triple_baseline = complement_triple_baseline_certificate(p)
    carried = _opposite_forcing(
        p, edge_count, 2 * p - 6, 2, 6, 14, 7, p + 15, p_plus_fifteen
    )
    carried.update(
        {
            "branch": P1_CARRIED_BRANCH,
            "changed_premise": (
                "the all-low complement-triple branch closed in 15.768 "
                "becomes m-1 low P=2 rows and one high P=3 row"
            ),
            "low_direction_count": m - 1,
            "unique_high_direction_count": 1,
            "forced_high_direction_parallel_count": 3,
            "hard_parallel_candidates": _parallel_candidates(
                p, edge_count, 2
            ),
            "hard_baseline_dependency_proved": triple_baseline["proved"],
            "signed_row_sum_unchanged_from_prop_15768": (
                carried["hard_sign_times_global_T"] == 8 - 3 * p
            ),
        }
    )

    family_catalog = p1_sharp_family_catalog(p)
    new_family_rows = []
    for family in family_catalog["families"]:
        offset = int(family["coefficient_offset"])
        candidates = _parallel_candidates(p, edge_count, offset)
        _require(candidates == [offset], "the p=1 sharp offset no longer fixes P")
        row = _opposite_forcing(
            p,
            edge_count,
            2 * p - 4,
            offset,
            8 - offset,
            12,
            9 - offset,
            p + 13,
            p_plus_thirteen,
        )
        row.update({**family, "hard_parallel_candidates": candidates})
        _require(
            row["hard_edge_count"] == m * offset,
            "the all-low p=1 sharp hard-edge count changed",
        )
        new_family_rows.append(row)

    old_literal.update(
        {
            "branch": BRANCH_P1_LAST,
            "hard_parallel_candidates": _parallel_candidates(p, edge_count, 5),
        }
    )
    old_XNOR.update(
        {
            "branch": BRANCH_B2,
            "hard_parallel_candidates": _parallel_candidates(p, edge_count, 4),
        }
    )
    branches = {
        BRANCH_P1_LAST: old_literal,
        BRANCH_B2: old_XNOR,
        P1_CARRIED_BRANCH: carried,
        P1_NEW_SHARP_BRANCH: {
            "family_catalog": family_catalog,
            "family_ledgers": new_family_rows,
            "common_forced_local_mass": p + 13,
            "local_dependency": p_plus_thirteen,
            "proved": all(row["proved"] for row in new_family_rows),
        },
    }
    proved = bool(
        residues["proved"]
        and residues["arithmetic_surviving_residues"]
        == [0, m - 4, m - 3, m - 1]
        and all(row["proved"] for row in branches.values())
        and carried["hard_baseline_dependency_proved"]
        and carried["signed_row_sum_unchanged_from_prop_15768"]
        and carried["hard_parallel_candidates"] == [2]
        and carried["hard_edge_count"] == 2 * (m - 1) + 3
        and old_literal["hard_parallel_candidates"] == [5]
        and old_XNOR["hard_parallel_candidates"] == [4]
    )
    _require(proved, "the p=1 next residual layer survived")
    return {
        "p": p,
        "p_mod_4": 1,
        "layer_index_t": (p - 5) // 2,
        "original_k": 5 * p - 5,
        "H_edge_count": edge_count,
        "residue_ledger": residues,
        "branch_exclusions": branches,
        "all_boundary_sizes_excluded": True,
        "finite_prime_graph_or_slice_census_used": False,
        "residual_ii_layer_excluded": proved,
        "proved": proved,
    }


def p3_next_layer_exclusion(p: int) -> dict[str, object]:
    """Close ``t=q-1``, equivalently ``k=5p-3``, for p=3 mod 4."""
    _check_prime(p, 3)
    m = (p + 1) // 2
    edge_count = 5 * p - 2
    residues = p3_next_residue_ledger(p)
    p_plus_nine = p_plus_nine_local_exclusion(p)
    p_plus_thirteen = p3_p_plus_thirteen_local_exclusion(p)
    p_minus_one = p_minus_one_local_exclusion(p)

    old_XNOR = _opposite_forcing(
        p, edge_count, p - 1, 4, 3, 8, 4, p + 9, p_plus_nine
    )
    old_literal = _opposite_forcing(
        p, edge_count, p - 1, 3, 4, 8, 5, p + 9, p_plus_nine
    )
    old_XNOR.update(
        {
            "branch": BRANCH_B2,
            "hard_parallel_candidates": _parallel_candidates(p, edge_count, 4),
        }
    )
    old_literal.update(
        {
            "branch": BRANCH_P3_LAST,
            "hard_parallel_candidates": _parallel_candidates(p, edge_count, 3),
        }
    )

    prior_catalog = p3_hard_family_catalog(p)
    carried_rows = []
    for family in prior_catalog["families"]:
        offset = int(family["coefficient_offset"])
        candidates = _parallel_candidates(p, edge_count, offset)
        _require(candidates == [offset], "the carried p=3 offset no longer fixes P")
        row = _opposite_forcing(
            p,
            edge_count,
            2 * p - 4,
            offset,
            8 - offset,
            12,
            9 - offset,
            p + 13,
            p_plus_thirteen,
        )
        row.update(
            {
                **family,
                "hard_parallel_candidates": candidates,
                "low_direction_count": m - 1,
                "unique_high_direction_count": 1,
                "forced_high_direction_parallel_count": offset + 1,
                "signed_row_sum_unchanged_from_prop_15769": (
                    row["hard_sign_times_global_T"]
                    == (p + 1) * offset - 5 * p + 4
                ),
            }
        )
        _require(
            row["hard_edge_count"] == m * offset + 1,
            "the carried p=3 hard-edge count changed",
        )
        carried_rows.append(row)

    branches = {
        BRANCH_B2: old_XNOR,
        BRANCH_P3_LAST: old_literal,
        P3_CARRIED_BRANCH: {
            "changed_premise": (
                "each all-low sharp family closed in 15.769 becomes m-1 "
                "low P rows and one high P+1 row"
            ),
            "prior_family_catalog_proved": prior_catalog["proved"],
            "family_ledgers": carried_rows,
            "common_forced_local_mass": p + 13,
            "proved": all(row["proved"] for row in carried_rows),
        },
        P3_NEW_LOCAL_BRANCH: {
            "residue": m - 2,
            "all_hard_quotients_equal_one": True,
            "baseline_boundary_values": [2, p - 1],
            "difference_scaled_mass": p - 1,
            "parity_baseline_rigidity_dependency_proved": prior_catalog[
                "proved"
            ],
            "local_exclusion": p_minus_one,
            "proved": bool(prior_catalog["proved"] and p_minus_one["proved"]),
        },
    }
    proved = bool(
        residues["proved"]
        and residues["arithmetic_surviving_residues"]
        == [m - 3, m - 2, m - 1]
        and all(row["proved"] for row in branches.values())
        and old_XNOR["hard_parallel_candidates"] == [4]
        and old_literal["hard_parallel_candidates"] == [3]
        and all(
            row["signed_row_sum_unchanged_from_prop_15769"]
            for row in carried_rows
        )
    )
    _require(proved, "the p=3 next residual layer survived")
    return {
        "p": p,
        "p_mod_4": 3,
        "layer_index_t": (p - 3) // 2,
        "original_k": 5 * p - 3,
        "H_edge_count": edge_count,
        "residue_ledger": residues,
        "branch_exclusions": branches,
        "all_boundary_sizes_excluded": True,
        "finite_prime_graph_or_slice_census_used": False,
        "residual_ii_layer_excluded": proved,
        "proved": proved,
    }


def proposition_15770() -> dict[str, object]:
    """Package both parameterized next-layer theorems."""
    # Keep the finite exceptional companion independently replayable.
    from e1_gmin_m4_p23_second_post_band_moment_close import (
        p23_second_post_band_moment_close,
    )

    p1_samples = (29, 37, 41, 53)
    p3_samples = (31, 43, 47, 59)
    p1_rows = {str(p): p1_next_layer_exclusion(p) for p in p1_samples}
    p3_rows = {str(p): p3_next_layer_exclusion(p) for p in p3_samples}
    exceptional_p23 = p23_second_post_band_moment_close()
    proved = bool(
        all(row["proved"] for row in p1_rows.values())
        and all(row["proved"] for row in p3_rows.values())
        and exceptional_p23["proved"]
        and exceptional_p23["p23_k112_closed"]
    )
    _require(proved, "Proposition 15.770 failed")
    return {
        "prop": "15.770",
        "title": "Next post-band residual layers by one-row carry",
        "result_status": (
            "proved two infinite one-layer extensions plus the p23 endpoint"
        ),
        "changed_premise": (
            "Props. 15.768--15.769 close the equality branches that become "
            "m-1 low rows plus one high row at the following layer"
        ),
        "p1_mod_4_statement": (
            "for every prime p>=29 congruent to 1 modulo 4, residual (ii) "
            "is empty at t=(p-5)/2, equivalently k=5p-5"
        ),
        "p3_mod_4_statement": (
            "for every prime p>=31 congruent to 3 modulo 4, residual (ii) "
            "is empty at t=(p-3)/2, equivalently k=5p-3"
        ),
        "p1_parameterized_threshold_replays": p1_rows,
        "p3_parameterized_threshold_replays": p3_rows,
        "p23_exceptional_eleven_root_certificate": exceptional_p23,
        "p23_same_layer_closed": True,
        "p23_original_k": 112,
        "fixed_four_bit_catalog_sha256": sharp_p_minus_three_four_bit_catalog()[
            "selected_tables_sha256"
        ],
        "finite_prime_graph_or_slice_census_used": False,
        "later_layers_closed": False,
        "residual_ii_closed_globally": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_prop15770.json"
    write_json_atomic(path, proposition_15770())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
