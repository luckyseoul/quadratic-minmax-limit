#!/usr/bin/env python3
r"""Prop. 15.769 -- the first ``p=3 (mod 4)`` layer beyond 15.752.

Put ``q=(p-1)/2`` and ``m=q+1``.  Proposition 15.752 stops at
``t=q-3`` when ``p=3 (mod 4)``.  At the next layer ``t=q-2``, the
phase-one endpoint baselines can carry a lift of sharp mass ``p-3``.
This module proves that the new branch is nevertheless impossible for
every prime ``p=3 (mod 4)``, ``p>=31``.  The exceptional prime ``p=23`` is
then closed separately by the equality-globalization and moment certificate
in :mod:`e1_gmin_m4_p23_post_band_moment_close`.

At the new residue ``u=t=m-3``, every hard quotient equals one and every
hard cell has scaled mean ``2p-4``.  Its parity baseline is either the XNOR
quadratic ``(1-x_i-x_j)^2`` or the complementary literal ``1-x_j``.  The
difference is twice a nonzero nonnegative integral quadratic ``B`` with

    4p E[B]=p-3.

Sharpness in Proposition 15.688 makes ``B`` Boolean.  A corrected
transposition-influence argument leaves at most five slice coordinates,
and cube influence leaves at most four active coordinates.  The fixed
four-bit catalog has exactly ten tables at density ``(p-3)/(4p)``: six
selected-pair tables and four all-equal-triple tables.  Complementing the
slice identifies the corresponding lifts as omitted pairs or all-equal
triples.

Adding either lift to either parity baseline gives coefficient offsets
``2,3,4,5``.  The common difference-row sum first forces one common hard
parallel count ``P``.  The offset congruence and edge bound then force
``P`` to equal its offset, so the four families cannot mix.  In every
family

    hT=(p+1)P-5p+4,
    a(Q)=(p+1)Q+hT-3p.

The row ``Q=8-P`` has mass twelve and is excluded by the phase-zero and
sharp-lift floors.  Exact edge accounting gives

    sum_L (Q_L-(9-P))=m-9,

so at least nine opposite directions have ``Q=9-P`` and scaled mean
``p+13``.  Their nonzero-boundary alternatives are lifts of mass twelve,
below ``p-3``; hence one obtains a nonzero nonnegative integral quadratic
``C`` with ``4p E[C]=p+13``.

That local mass is impossible for ``p>=31``.  At height at least two,
paired-cube averaging produces a half-mean cube and Proposition 15.751
bounds the height by three, contradicting ``H>=(p-11)/4``.  At height one,
the Johnson influence bound is below eight, all remaining patterns extend
to the complementary slice, cube influence leaves four coordinates, and
the same fixed catalog misses density ``(p+13)/(4p)``.

The two old endpoint branches still force Proposition 15.752's forbidden
mass ``p+9``.  Therefore ``k=5p-5`` is empty for every prime
``p=3 (mod 4)``, ``p>=31``.  At ``p=23`` the local quadratic
``3-2r+binom(r,2)`` has mass ``p+13``, so the parameterized local exclusion
does not apply.  Proposition 15.769's exceptional companion instead
classifies and globalizes all equality forms and excludes the sole compatible
five-coordinate form with a quartic/octic moment certificate.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P3_LAST,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15751 import (
    cube_half_mean_height_certificate,
    exact_four_cube_catalog,
    profile_density,
)
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
NEW_BRANCH = "hard_sharp_p_minus_3_lift"
LOCAL_MASS_OFFSET = 13


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p>=31 congruent to 3 modulo 4")


def _four_bit_anf(table: int) -> tuple[int, ...]:
    """Return integer multilinear coefficients in subset-mask order."""
    coefficients = [(table >> mask) & 1 for mask in range(16)]
    for bit in range(4):
        for mask in range(16):
            if mask & (1 << bit):
                coefficients[mask] -= coefficients[mask ^ (1 << bit)]
    return tuple(coefficients)


def _four_bit_layer_counts(table: int) -> tuple[int, ...]:
    return tuple(
        sum((table >> mask) & 1 for mask in range(16) if mask.bit_count() == weight)
        for weight in range(5)
    )


@lru_cache(maxsize=1)
def sharp_p_minus_three_four_bit_catalog() -> dict[str, object]:
    """Classify the fixed four-bit tables at density ``(p-3)/(4p)``."""
    live_catalog = exact_four_cube_catalog()
    target_profiles = {(0, 0, 1, 2, 1), (1, 1, 0, 1, 1)}
    catalog_profiles = {
        tuple(int(value) for value in row["layer_counts"])
        for row in live_catalog["profiles"]
    }
    selected: list[dict[str, object]] = []
    pair_count = 0
    triple_count = 0
    for table in range(1 << 16):
        coefficients = _four_bit_anf(table)
        if any(
            value
            for mask, value in enumerate(coefficients)
            if mask.bit_count() > 2
        ):
            continue
        profile = _four_bit_layer_counts(table)
        if profile not in target_profiles:
            continue
        nonzero = {
            mask: value for mask, value in enumerate(coefficients) if value
        }
        pair_masks = [mask for mask in nonzero if mask.bit_count() == 2]
        linear_masks = [mask for mask in nonzero if mask.bit_count() == 1]
        if (
            len(nonzero) == 1
            and len(pair_masks) == 1
            and nonzero[pair_masks[0]] == 1
        ):
            family = "selected_pair_on_complementary_slice"
            pair_count += 1
        else:
            support = {
                mask.bit_length() - 1
                for mask in linear_masks
                if nonzero[mask] == -1
            }
            expected_pairs = {
                (1 << left) | (1 << right)
                for left, right in combinations(sorted(support), 2)
            }
            _require(
                nonzero.get(0) == 1
                and len(support) == 3
                and set(pair_masks) == expected_pairs
                and all(nonzero[mask] == 1 for mask in pair_masks)
                and len(nonzero) == 7,
                "an unclassified sharp-density four-bit table appeared",
            )
            family = "all_equal_triple"
            triple_count += 1
        selected.append(
            {
                "table": table,
                "table_hex": f"{table:04x}",
                "layer_counts": list(profile),
                "family": family,
                "anf_coefficients": list(coefficients),
            }
        )

    payload = json.dumps(
        selected, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    digest = hashlib.sha256(payload).hexdigest()
    proved = bool(
        live_catalog["proved"]
        and target_profiles <= catalog_profiles
        and len(selected) == 10
        and pair_count == 6
        and triple_count == 4
    )
    _require(proved, "the sharp p-3 four-bit catalog changed")
    return {
        "source_catalog_valid_tables": live_catalog["valid_tables"],
        "source_catalog_sha256": live_catalog["valid_table_signature_sha256"],
        "target_symbolic_density": "(p-3)/(4p)",
        "target_layer_profiles": [list(row) for row in sorted(target_profiles)],
        "selected_table_count": len(selected),
        "selected_pair_table_count": pair_count,
        "all_equal_triple_table_count": triple_count,
        "complement_back_to_original_slice": {
            "selected_pair": "(1-x_i)(1-x_j), an omitted-pair lift",
            "all_equal_triple": "all-equal is invariant under complementation",
        },
        "selected_tables_sha256": digest,
        "tables": selected,
        "proved": proved,
    }


def sharp_p_minus_three_boolean_classification(p: int) -> dict[str, object]:
    """Reduce every sharp ``p-3`` Boolean lift to two fixed families."""
    _check_prime(p)
    q = (p - 1) // 2
    mu = Fraction(p - 3, 4 * p)
    influence_floor = Fraction(
        (p + 1) * (p - 3), 16 * p * (p - 2)
    )
    total_influence_upper = (p - 1) * mu * (1 - mu)
    junta_bound = Fraction(6 * (p - 1) * (p - 2), p * p)
    catalog = sharp_p_minus_three_four_bit_catalog()
    matching_profiles = [
        tuple(int(value) for value in row["layer_counts"])
        for row in exact_four_cube_catalog()["profiles"]
        if profile_density(tuple(row["layer_counts"]), p) == mu
    ]
    expected_profiles = [(0, 0, 1, 2, 1), (1, 1, 0, 1, 1)]
    proved = bool(
        junta_bound < 6
        and 5 < q
        and sorted(matching_profiles) == expected_profiles
        and catalog["proved"]
    )
    _require(proved, "the sharp p-3 Boolean classification failed")
    return {
        "p": p,
        "complementary_slice": f"J({p},{q})",
        "density": str(mu),
        "relevant_pair_influence_floor": str(influence_floor),
        "total_influence_upper_bound": str(total_influence_upper),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "junta_coordinates_at_most": 5,
        "five_less_than_both_slice_sides": 5 < q,
        "all_junta_patterns_extend_to_slice": True,
        "symmetrized_representative_extends_to_boolean_cube": True,
        "cube_total_influence_upper_bound": 2,
        "cube_relevant_coordinate_influence_floor": "1/2",
        "cube_active_coordinates_at_most": 4,
        "matching_four_bit_layer_profiles": [
            list(row) for row in matching_profiles
        ],
        "four_bit_classification": catalog,
        "original_slice_families": ["omitted_pair", "all_equal_triple"],
        "proved": proved,
    }


def hard_family_catalog(p: int) -> dict[str, object]:
    """Combine the two parity baselines with the two sharp lift families."""
    _check_prime(p)
    b2 = parity_floor_certificate(p, 2, 1)
    literal = parity_floor_certificate(p, 1, 1)
    classification = sharp_p_minus_three_boolean_classification(p)

    baseline_checks = []
    for left, right in product((0, 1), repeat=2):
        z_left, z_right = 2 * left - 1, 2 * right - 1
        baseline_checks.append(
            3 + 2 * (1 - left - right) ** 2
            == 4 + z_left * z_right
        )
    literal_checks = [
        3 + 2 * (1 - value) == 4 - (2 * value - 1)
        for value in (0, 1)
    ]
    omitted_checks = []
    for left, right in product((0, 1), repeat=2):
        z_left, z_right = 2 * left - 1, 2 * right - 1
        omitted_checks.append(
            4 * (1 - left) * (1 - right)
            == 1 - z_left - z_right + z_left * z_right
        )
    triple_checks = []
    for bits in product((0, 1), repeat=3):
        z = tuple(2 * bit - 1 for bit in bits)
        value = 1 - sum(bits) + sum(
            bits[left] * bits[right] for left, right in combinations(range(3), 2)
        )
        triple_checks.append(
            4 * value == 1 + sum(
                z[left] * z[right] for left, right in combinations(range(3), 2)
            )
        )

    rows = []
    for baseline_name, baseline_b, baseline_offset in (
        ("complement_literal", p - 1, 3),
        ("XNOR", 2, 4),
    ):
        for lift_name, lift_offset in (
            ("omitted_pair", -1),
            ("all_equal_triple", 1),
        ):
            offset = baseline_offset + lift_offset
            rows.append(
                {
                    "baseline": baseline_name,
                    "baseline_b": baseline_b,
                    "baseline_scaled_mean": p - 1,
                    "lift": lift_name,
                    "lift_scaled_mass": p - 3,
                    "total_scaled_mean": 2 * p - 4,
                    "baseline_offset": baseline_offset,
                    "lift_offset_increment": lift_offset,
                    "coefficient_offset": offset,
                    "coefficient_congruence": (
                        f"{(p - 1) // 2} divides I+P-{offset}"
                    ),
                }
            )
    proved = bool(
        b2["exact_positive_quadrature_certificate"]
        and literal["exact_positive_quadrature_certificate"]
        and all(weight > 0 for weight in b2["quadrature_weights"])
        and all(weight > 0 for weight in literal["quadrature_weights"])
        and all(baseline_checks + literal_checks + omitted_checks + triple_checks)
        and classification["proved"]
        and sorted(int(row["coefficient_offset"]) for row in rows)
        == [2, 3, 4, 5]
    )
    _require(proved, "the sharp-lift hard family catalog failed")
    return {
        "p": p,
        "parity_baselines": {
            "b=2": "A_0=(1-x_i-x_j)^2, target 4+z_i*z_j, offset 4",
            "b=p-1": "A_0=1-x_j, target 4-z_j, offset 3",
        },
        "difference_lift": "B=(A-A_0)/2",
        "difference_is_nonzero_nonnegative_integral": True,
        "difference_scaled_mass": "4p E[B]=p-3",
        "sharp_lift_dependency": "Proposition 15.688",
        "sharp_lift_forces_boolean": True,
        "boolean_classification": classification,
        "families": rows,
        "proved": proved,
    }


def p_plus_thirteen_local_exclusion(p: int) -> dict[str, object]:
    """Exclude ``4p E[C]=p+13`` on ``J(p,(p+1)/2)``."""
    _check_prime(p)
    q = (p - 1) // 2
    m = q + 1
    scaled_mass = p + LOCAL_MASS_OFFSET
    lower_height = Fraction(p - 11, 4)
    stabilizer_upper = Fraction(scaled_mass, 4)
    paired_average_upper = Fraction(scaled_mass, 2 * (p + 1))
    half_mean = cube_half_mean_height_certificate()
    height_proved = bool(
        half_mean["proved"]
        and lower_height > 3
        and paired_average_upper < Fraction(3, 4)
    )

    mu = Fraction(scaled_mass, 4 * p)
    influence_floor = Fraction(
        (p + 1) * (p - 3), 16 * p * (p - 2)
    )
    total_influence_upper = (p - 1) * mu * (1 - mu)
    junta_bound = Fraction(
        2 * (p - 1) * (p - 2) * (p + 13) * (3 * p - 13),
        p * p * (p + 1) * (p - 3),
    )
    eight_gap = p**4 - 25 * p**3 + 229 * p**2 - 559 * p + 338
    x = p - 31
    translated_gap = (
        x**4 + 99 * x**3 + 3670 * x**2 + 60728 * x + 381824
    )
    catalog = exact_four_cube_catalog()
    densities = sorted(
        {
            profile_density(tuple(row["layer_counts"]), p)
            for row in catalog["profiles"]
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
        eight_gap == translated_gap
        and eight_gap > 0
        and junta_bound < 8
        and 7 < q
        and catalog["proved"]
        and densities == expected_densities
        and Fraction(p + 1, 4 * p) < mu < Fraction(p - 1, 2 * p)
        and mu not in densities
    )
    proved = height_proved and boolean_proved
    _require(proved, "the p+13 local exclusion failed")
    return {
        "p": p,
        "slice": f"J({p},{m})",
        "statement": (
            "no nonzero nonnegative integer-valued quadratic C has "
            "4p E[C]=p+13"
        ),
        "height_at_least_two": {
            "height_lower_bound": str(lower_height),
            "stabilizer_height_upper_bound": str(stabilizer_upper),
            "paired_cube_average_upper_bound": str(paired_average_upper),
            "some_paired_cube_has_mean_exactly": "1/2",
            "half_mean_cube_height_upper_bound": 3,
            "proved": height_proved,
        },
        "height_one_boolean": {
            "target_density": str(mu),
            "relevant_pair_influence_floor": str(influence_floor),
            "total_influence_upper_bound": str(total_influence_upper),
            "largest_zero_influence_class_complement_bound": str(junta_bound),
            "eight_gap_polynomial": "p^4-25p^3+229p^2-559p+338",
            "eight_gap_at_p_equals_x_plus_31": [1, 99, 3670, 60728, 381824],
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


def p23_local_threshold_witness() -> dict[str, object]:
    """Record why the new local theorem cannot include ``p=23``."""
    p = 23
    m = 12
    values = [3 - 2 * r + comb(r, 2) for r in range(5)]
    numerator = sum(
        comb(4, r) * comb(p - 4, m - r) * values[r]
        for r in range(5)
    )
    mean = Fraction(numerator, comb(p, m))
    proved = values == [3, 1, 0, 0, 1] and 4 * p * mean == p + 13
    _require(proved, "the p=23 p+13 local witness changed")
    return {
        "p": p,
        "formula": "C=3-2r+binom(r,2), r=|X intersect R|, |R|=4",
        "layer_values": values,
        "mean": str(mean),
        "scaled_mass_4p_E_C": p + 13,
        "is_only_a_local_quadratic_not_a_residual_graph": True,
        "proved": proved,
    }


def first_uncovered_p3_residue_ledger(p: int) -> dict[str, object]:
    """Classify the phase-one residues at ``t=q-2``."""
    _check_prime(p)
    q = (p - 1) // 2
    m = q + 1
    t = q - 2
    edge_count = 5 * p - 4
    phase_one = residual_even_floor_table(p)["phase_one_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    low_b = sorted(
        int(b) for b, floor in phase_one.items() if int(floor) < 2 * p
    )
    residue_rows: list[dict[str, object]] = []
    for u in range(m):
        quotient_sum = m + t - u
        if u < t:
            status = "excluded: k=0 is sub-floor and k=1 is a nonzero lift below p-3"
        elif u == t:
            status = NEW_BRANCH
        elif u == m - 2:
            status = "excluded: k=0 is sub-floor but sum k<m"
        elif u == m - 1:
            status = "old endpoint A/C dichotomy"
        else:  # pragma: no cover - t=m-3 leaves no other interval
            raise ArithmeticError("unexpected residue interval")
        residue_rows.append(
            {
                "u": u,
                "quotient_sum": quotient_sum,
                "k_zero_mean": 2 * u,
                "k_one_mean": p + 1 + 2 * u,
                "status": status,
            }
        )
    proved = bool(
        t == m - 3
        and edge_count == 4 * p + 2 * t + 1
        and p * p + 1 - 2 * edge_count > 0
        and low_b == [2, p - 1]
        and all(int(phase_one[b]) == p - 1 for b in low_b)
        and p + 1 + 2 * t == 2 * p - 4
        and (2 * p - 4) - (p - 1) == lift_floor == p - 3
        and m + t - t == m
        and m + t - (m - 2) == m - 1
    )
    _require(proved, "the first-uncovered p=3 residue ledger failed")
    return {
        "p": p,
        "q": q,
        "m": m,
        "layer_index_t": t,
        "original_k": 5 * p - 5,
        "H_edge_count": edge_count,
        "guaranteed_isolated_vertices": p * p + 1 - 2 * edge_count,
        "phase_one_mean_form": f"a_L=2u+{p + 1}k_L",
        "phase_one_quotient_sum": "sum k_L=m+t-u",
        "phase_one_low_b_values": low_b,
        "residue_rows": residue_rows,
        "surviving_residues": [t, m - 1],
        "new_branch_all_quotients_equal_one": True,
        "new_branch_scaled_mean": 2 * p - 4,
        "new_branch_baseline_scaled_mean": p - 1,
        "new_branch_sharp_lift_mass": p - 3,
        "possible_branches": [BRANCH_B2, BRANCH_P3_LAST, NEW_BRANCH],
        "proved": proved,
    }


def _old_endpoint_branch_extension(p: int, branch: str) -> dict[str, object]:
    """Replay the old A/C branch one step beyond 15.752's band."""
    _check_prime(p)
    m = (p + 1) // 2
    edge_count = 5 * p - 4
    if branch == BRANCH_B2:
        hard_edges = 5 * m - 2
        hT = 5
        minimum_Q = 3
        next_Q = 4
    elif branch == BRANCH_P3_LAST:
        hard_edges = 4 * m - 2
        hT = 4 - p
        minimum_Q = 4
        next_Q = 5
    else:
        raise ValueError("branch must be hard_b2 or p3_all_low_b_p_minus_1")
    opposite_edges = edge_count - hard_edges
    minimum_mean = (p + 1) * minimum_Q + hT - 3 * p
    next_mean = (p + 1) * next_Q + hT - 3 * p
    surplus_after_next = opposite_edges - m * next_Q
    dependency = p_plus_nine_local_exclusion(p)
    proved = bool(
        hard_edges - opposite_edges == hT
        and minimum_mean == 8
        and next_mean == p + 9
        and surplus_after_next == m - 7
        and 0 <= surplus_after_next < m
        and dependency["proved"]
    )
    _require(proved, f"the {branch} first-outband extension failed")
    return {
        "branch": branch,
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "hard_sign_times_global_T": hT,
        "minimum_opposite_Q": minimum_Q,
        "minimum_opposite_mean": minimum_mean,
        "minimum_cell_excluded": True,
        "forced_next_Q": next_Q,
        "forced_next_scaled_mean": next_mean,
        "surplus_after_raising_every_Q": surplus_after_next,
        "directions_at_next_Q_at_least": 7,
        "local_dependency": "Proposition 15.752 p+9 exclusion",
        "excluded": proved,
        "proved": proved,
    }


def sharp_lift_branch_exclusion(p: int) -> dict[str, object]:
    """Exclude all four new sharp-lift hard families."""
    _check_prime(p)
    q = (p - 1) // 2
    m = q + 1
    edge_count = 5 * p - 4
    hard_mean = 2 * p - 4
    family_catalog = hard_family_catalog(p)
    local = p_plus_thirteen_local_exclusion(p)
    phase_zero = residual_even_floor_table(p)["phase_zero_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    parallel_upper = edge_count // m
    rows: list[dict[str, object]] = []
    for family in family_catalog["families"]:
        offset = int(family["coefficient_offset"])
        parallel_candidates = [
            value
            for value in range(parallel_upper + 1)
            if (value - offset) % q == 0
        ]
        P = parallel_candidates[0]
        hard_edges = m * P
        opposite_edges = edge_count - hard_edges
        hT = (p + 1) * P - 3 * p - hard_mean

        def opposite_mean(Q: int) -> int:
            return (p + 1) * Q + hT - 3 * p

        Q0 = 8 - P
        Q1 = 9 - P
        surplus_after_Q1 = opposite_edges - m * Q1
        next_mean = opposite_mean(Q1)
        nonzero_rows = [
            (int(b), int(floor), next_mean - int(floor))
            for b, floor in phase_zero.items()
            if int(b) != 0 and int(floor) <= next_mean
        ]
        proved = bool(
            parallel_candidates == [offset]
            and P == offset
            and hT == (p + 1) * P - 5 * p + 4
            and hard_edges - opposite_edges == hT
            and opposite_mean(Q0 - 1) == 11 - p < 0
            and opposite_mean(Q0) == 12
            and 0 < 12 < min(p + 1, lift_floor)
            and surplus_after_Q1 == m - 9
            and 0 <= surplus_after_Q1 < m
            and next_mean == p + 13
            and [row[0] for row in nonzero_rows] == [2, p - 1]
            and [row[2] for row in nonzero_rows] == [12, 12]
            and all(0 < excess < lift_floor for _, _, excess in nonzero_rows)
            and local["proved"]
        )
        _require(proved, "a sharp-lift hard family survived")
        rows.append(
            {
                **family,
                "hard_parallel_upper_bound": parallel_upper,
                "hard_parallel_candidates": parallel_candidates,
                "forced_hard_parallel_count": P,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "hard_sign_times_global_T": hT,
                "opposite_mean_formula": f"a(Q)={p + 1}Q+({hT})-{3 * p}",
                "last_negative_Q": Q0 - 1,
                "last_negative_mean": opposite_mean(Q0 - 1),
                "forbidden_mass_twelve_Q": Q0,
                "mass_twelve_Q_mean": 12,
                "forced_low_Q": Q1,
                "surplus_after_every_Q_at_least_forced_low_Q": surplus_after_Q1,
                "directions_at_forced_low_Q_at_least": 9,
                "forced_low_Q_scaled_mean": next_mean,
                "nonzero_b_floor_and_lift_rows": [list(row) for row in nonzero_rows],
                "forced_low_Q_is_b_zero": True,
                "forced_local_cell": "A=2C, 4p E[C]=p+13",
                "excluded": proved,
                "proved": proved,
            }
        )

    proved = bool(
        family_catalog["proved"]
        and sorted(int(row["forced_hard_parallel_count"]) for row in rows)
        == [2, 3, 4, 5]
        and all(row["proved"] for row in rows)
    )
    _require(proved, "the sharp-lift branch exclusion failed")
    return {
        "branch": NEW_BRANCH,
        "hard_family_catalog": family_catalog,
        "common_row_sum_identity": (
            "sum q_L=p*P_L-3p-(2p-4)=hT-P_L"
        ),
        "equal_hard_means_force_one_common_P": True,
        "different_offsets_cannot_mix_mod_q": True,
        "family_ledgers": rows,
        "common_forced_opposite_local_mass": "p+13",
        "local_mass_exclusion": local,
        "excluded": proved,
        "proved": proved,
    }


def first_uncovered_p3_layer_exclusion(p: int) -> dict[str, object]:
    """Close ``t=(p-5)/2`` for one admissible prime."""
    _check_prime(p)
    residues = first_uncovered_p3_residue_ledger(p)
    branches = {
        BRANCH_B2: _old_endpoint_branch_extension(p, BRANCH_B2),
        BRANCH_P3_LAST: _old_endpoint_branch_extension(p, BRANCH_P3_LAST),
        NEW_BRANCH: sharp_lift_branch_exclusion(p),
    }
    proved = bool(
        residues["proved"]
        and set(branches) == set(residues["possible_branches"])
        and all(row["proved"] for row in branches.values())
    )
    _require(proved, "the first uncovered p=3 mod4 layer did not close")
    return {
        "p": p,
        "p_mod_4": 3,
        "layer_index_t": (p - 5) // 2,
        "original_k": 5 * p - 5,
        "H_edge_count": 5 * p - 4,
        "residue_ledger": residues,
        "branch_exclusions": branches,
        "all_boundary_sizes_excluded": True,
        "finite_prime_graph_or_slice_census_used": False,
        "residual_ii_layer_excluded": proved,
        "proved": proved,
    }


def proposition_15769() -> dict[str, object]:
    """Package the parameterized theorem and threshold replays."""
    # Imported lazily so the exceptional certificate remains independently
    # replayable and the main parameterized module keeps a shallow import DAG.
    from e1_gmin_m4_p23_post_band_moment_close import p23_post_band_moment_close

    sample_primes = (31, 43, 47, 59)
    rows = {str(p): first_uncovered_p3_layer_exclusion(p) for p in sample_primes}
    exceptional_p23 = p23_post_band_moment_close()
    proved = bool(
        all(row["proved"] for row in rows.values())
        and exceptional_p23["proved"]
        and exceptional_p23["p23_k110_closed"]
    )
    return {
        "prop": "15.769",
        "title": "First post-15.752 p=3 mod 4 residual layer",
        "result_status": "proved infinite-family theorem with fixed four-bit catalog",
        "statement": (
            "for every prime p>=31 congruent to 3 modulo 4, residual (ii) "
            "is empty at t=(p-5)/2, equivalently k=5p-5"
        ),
        "new_local_theorem": (
            "no nonzero nonnegative integral quadratic on J(p,(p+1)/2) "
            "has 4p E[C]=p+13"
        ),
        "first_layer_beyond_prop_15752": True,
        "sharp_p_minus_three_four_bit_catalog": sharp_p_minus_three_four_bit_catalog(),
        "parameterized_threshold_replays": rows,
        "p23_local_threshold_witness": p23_local_threshold_witness(),
        "p23_exceptional_equality_moment_certificate": exceptional_p23,
        "p23_same_layer_closed": True,
        "p23_original_k": 110,
        "later_layers_closed": False,
        "prop_15274_slope_scope": (
            "only the dual-bad/two-level subcase; it does not apply to this "
            "multi-level isolated-chart equality branch"
        ),
        "residual_ii_closed_globally": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_prop15769.json"
    write_json_atomic(path, proposition_15769())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
