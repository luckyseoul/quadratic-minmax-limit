#!/usr/bin/env python3
r"""Close the exceptional ``p=23, k=110`` first post-band endpoint.

The parameterized first-post-band argument for ``p=3 (mod 4)`` reduces the
new branch to a nonzero nonnegative integral quadratic ``C`` on
``J(23,12)`` with

    92 E[C] = 36.

At this exceptional prime the local mass is attainable, so a local mass
exclusion is false.  The endpoint nevertheless closes because equality is
rigid.  The height is three and every paired cube through a maximizer has
mean one half.  An all-dimensional equality theorem classifies every such
cube restriction as

    F_r(s) = 3 - 2s + binom(s,2),    r in {4,5},

up to dummy coordinates.  Compatibility between different pairings then
globalizes the slice quadratic to the same four- or five-coordinate form.
The coefficient congruence selects only ``F_5`` in the hard ``P=4`` branch.

Every hard coefficient row in that branch is a triangle minus a full star.
Full-star power sums vanish over ``F_23`` in degrees below 22, while triangle
moments satisfy the homogeneous identities

    G4 = 2 h M4 - M2^2 = 0,
    G8 = 24 h M8 - 32 M2 M6 + 5 M2^4 = 0.

Twelve hard projective directions make both binary forms identically zero.
On an opposite ``F_5`` row their values become

    G4 = -2 S4 - S2^2,
    G8 = -24 S8 - 32 S2 S6 + 5 S2^4,

where ``Sd`` is the degree-``d`` edge power sum of a five-set in ``F_23``.
An exact scan of all ``binom(23,5)=33,649`` five-sets (also compressed into
69 affine orbits) finds 1,518 zeros of ``G4``, 2,024 zeros of ``G8``, and no
simultaneous zero.  This is an exhaustive finite coefficient certificate,
not a residual graph census and not a global residual-(ii) closure.
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

from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
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
from e1_gmin_m4_prop15768 import cube_three_quarter_height_certificate
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
P = 23
Q = 11
M = 12
LAYER_INDEX = 9
ORIGINAL_K = 110
H_EDGE_COUNT = 111
LOCAL_SCALED_MASS = 36
NEW_BRANCH = "hard_sharp_p_minus_3_lift"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _f_value(s: int) -> int:
    """The common value polynomial for the sharp four/five-coordinate forms."""
    return 3 - 2 * s + comb(s, 2)


@lru_cache(maxsize=1)
def half_mean_height_three_cube_classification() -> dict[str, object]:
    """Classify all half-mean integral cube quadratics with maximum three.

    If ``g(0)=3`` and ``E[g]=1/2``, a facet through the origin has mean
    ``1/2`` or ``3/4``.  Induction handles a half-mean facet.  An ``F4``
    facet has only the dummy extension or the extension ``ell=-2+s`` to
    ``F5``; an ``F5`` facet has only the dummy extension.  If every origin
    facet has mean ``3/4``, the opposite facets have mean ``1/4`` and are
    Boolean.  The coefficient sums force dimension four and then force all
    linear coefficients to be ``-2`` and all pair coefficients to be one.
    """
    half = cube_half_mean_height_certificate()
    f4_layers = tuple(_f_value(s) for s in range(5))
    f5_layers = tuple(_f_value(s) for s in range(6))
    f4_mass = sum(comb(4, s) * f4_layers[s] for s in range(5))
    f5_mass = sum(comb(5, s) * f5_layers[s] for s in range(6))

    # For an F4 half-mean facet, nonnegativity on its six zero two-sets
    # forces all averaged affine differences there to vanish.  Hence the
    # difference is c+u*s with c=-2u.  The endpoint value bounds and
    # integrality leave exactly these two possibilities.
    f4_extension_candidates = [
        (c, u)
        for c in range(-3, 1)
        for u in range(-3, 4)
        if c == -2 * u and -3 <= c <= 0
    ]
    f4_extension_checks = {
        "dummy": all(
            0 <= _f_value(s) <= 3 for s in range(5)
        ),
        "to_F5": all(
            _f_value(s) + (-2 + s) == _f_value(s + 1)
            for s in range(5)
        ),
    }

    # In the F5 case the empty and full vertices both have value three.
    # Their averaged affine differences are nonpositive and sum to zero;
    # the zero two-layers then force every active coefficient to vanish.
    f5_endpoint_values = (f5_layers[0], f5_layers[-1])
    f5_pair_zero = all(f5_layers[s] == 0 for s in (2, 3))

    # If every origin facet has mean 3/4, write
    # g=3+sum a_i*x_i+sum b_ij*x_i*x_j.  The displayed equations give
    # A=-10+d/2 with a_i in {-3,-2}; hence d is even and d<=4.  Total mass
    # gives d>=3, so d=4, A=-8, and every row of b has three entries >=1
    # summing to three.
    all_three_quarter_dimensions = [
        d
        for d in range(1, 20)
        if d % 2 == 0
        and -10 + Fraction(d, 2) <= -2 * d
        and 2 ** (d - 1) >= 3
    ]
    d = 4
    A = -10 + Fraction(d, 2)
    B = -10 - 2 * A
    row_sum = Fraction(-1) - 2 * Fraction(-2)

    proved = bool(
        half["proved"]
        and f4_layers == (3, 1, 0, 0, 1)
        and f5_layers == (3, 1, 0, 0, 1, 3)
        and f4_mass == 1 << 3
        and f5_mass == 1 << 4
        and f4_extension_candidates == [(-2, 1), (0, 0)]
        and all(f4_extension_checks.values())
        and f5_endpoint_values == (3, 3)
        and f5_pair_zero
        and all_three_quarter_dimensions == [4]
        and A == -8
        and B == 6
        and row_sum == 3
    )
    _require(proved, "the half-mean height-three equality classification changed")
    return {
        "hypotheses": (
            "g is a nonnegative integer-valued polynomial of degree at most "
            "two on a Boolean cube, E[g]=1/2, and g(0)=max(g)=3"
        ),
        "origin_facet_mean_options": ["1/2", "3/4"],
        "excluded_origin_facet_means": {
            "1/4": "support-floor equality would make a facet containing 3 Boolean",
            "1": (
                "the zero opposite facet makes the origin facet affine; mean one "
                "would give coefficient sum -4 but nonnegativity gives at least -3"
            ),
        },
        "induction_half_mean_facet": {
            "F4_extension_affine_pairs_c_u": [list(row) for row in f4_extension_candidates],
            "F4_extensions": ["dummy", "F5 via ell=-2+s"],
            "F5_empty_and_full_values": list(f5_endpoint_values),
            "F5_two_and_three_layers_are_zero": f5_pair_zero,
            "F5_extensions": ["dummy"],
        },
        "all_three_quarter_facets": {
            "global_mean_equation": "2A+B=-10",
            "facet_equations": "2a_i+sum_(j!=i)b_ij=-1",
            "summed_facet_equation": "A+B=-d/2",
            "therefore_A": "-10+d/2",
            "opposite_facets_boolean_implies": "a_i in {-3,-2}",
            "possible_dimensions": all_three_quarter_dimensions,
            "forced_dimension": d,
            "forced_linear_sum": int(A),
            "forced_pair_sum": int(B),
            "forced_pair_row_sum": int(row_sum),
            "forced_coefficients": "a_i=-2 and b_ij=1",
        },
        "classified_forms": {
            "F4": {
                "formula": "3-2s+binom(s,2)",
                "active_coordinates": 4,
                "layer_values": list(f4_layers),
                "mean": str(Fraction(f4_mass, 1 << 4)),
            },
            "F5": {
                "formula": "3-2s+binom(s,2)",
                "active_coordinates": 5,
                "layer_values": list(f5_layers),
                "mean": str(Fraction(f5_mass, 1 << 5)),
            },
        },
        "dummy_coordinates_allowed": True,
        "classification_exhaustive_in_every_dimension": True,
        "proved": proved,
    }


def _four_bit_anf(table: int) -> tuple[int, ...]:
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
def p23_sharp_hard_family_catalog() -> dict[str, object]:
    """Classify the sharp ``p-3`` hard lift and its four target offsets."""
    target_density = Fraction(P - 3, 4 * P)
    target_profiles = {(0, 0, 1, 2, 1), (1, 1, 0, 1, 1)}
    catalog = exact_four_cube_catalog()
    selected: list[dict[str, object]] = []
    family_counts: Counter[str] = Counter()
    for table in range(1 << 16):
        coefficients = _four_bit_anf(table)
        if any(
            value for mask, value in enumerate(coefficients) if mask.bit_count() > 2
        ):
            continue
        profile = _four_bit_layer_counts(table)
        if profile not in target_profiles:
            continue
        nonzero = {mask: value for mask, value in enumerate(coefficients) if value}
        pair_masks = [mask for mask in nonzero if mask.bit_count() == 2]
        linear_masks = [mask for mask in nonzero if mask.bit_count() == 1]
        if (
            len(nonzero) == 1
            and len(pair_masks) == 1
            and nonzero[pair_masks[0]] == 1
        ):
            family = "selected_pair_on_complementary_slice"
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
                "an unclassified sharp p23 four-bit table appeared",
            )
            family = "all_equal_triple"
        family_counts[family] += 1
        selected.append(
            {
                "table": table,
                "table_hex": f"{table:04x}",
                "layer_counts": list(profile),
                "family": family,
            }
        )

    selected_digest = hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    junta_bound = Fraction(6 * (P - 1) * (P - 2), P * P)
    matching_profiles = sorted(
        tuple(int(value) for value in row["layer_counts"])
        for row in catalog["profiles"]
        if profile_density(tuple(row["layer_counts"]), P) == target_density
    )

    # Combine the two phase-one equality baselines with the two sharp lift
    # families.  The signed offset is target constant plus all linear terms.
    b2 = parity_floor_certificate(P, 2, 1)
    literal = parity_floor_certificate(P, 1, 1)
    families = []
    for baseline, baseline_b, baseline_offset in (
        ("complement_literal", P - 1, 3),
        ("XNOR", 2, 4),
    ):
        for lift, lift_increment in (("omitted_pair", -1), ("all_equal_triple", 1)):
            families.append(
                {
                    "baseline": baseline,
                    "baseline_b": baseline_b,
                    "lift": lift,
                    "baseline_offset": baseline_offset,
                    "lift_offset_increment": lift_increment,
                    "coefficient_offset": baseline_offset + lift_increment,
                }
            )

    proved = bool(
        catalog["proved"]
        and target_density == Fraction(5, 23)
        and junta_bound == Fraction(2772, 529) < 6
        and matching_profiles == sorted(target_profiles)
        and len(selected) == 10
        and family_counts
        == Counter(
            {
                "selected_pair_on_complementary_slice": 6,
                "all_equal_triple": 4,
            }
        )
        and b2["exact_positive_quadrature_certificate"]
        and literal["exact_positive_quadrature_certificate"]
        and sorted(int(row["coefficient_offset"]) for row in families)
        == [2, 3, 4, 5]
    )
    _require(proved, "the p23 sharp hard-family catalog changed")
    return {
        "p": P,
        "sharp_scaled_mass": P - 3,
        "sharp_floor_equality_forces_boolean": True,
        "complementary_slice_density": str(target_density),
        "corrected_transposition_junta_bound": str(junta_bound),
        "junta_coordinates_at_most": 5,
        "cube_active_coordinates_at_most": 4,
        "matching_four_bit_profiles": [list(row) for row in matching_profiles],
        "selected_table_count": len(selected),
        "selected_pair_table_count": family_counts[
            "selected_pair_on_complementary_slice"
        ],
        "all_equal_triple_table_count": family_counts["all_equal_triple"],
        "selected_tables_sha256": selected_digest,
        "original_slice_lift_families": ["omitted_pair", "all_equal_triple"],
        "hard_families": families,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_local_height_equality() -> dict[str, object]:
    """Pin the forced local mass-36 quadratic to height three."""
    half = cube_half_mean_height_certificate()
    three_quarter = cube_three_quarter_height_certificate()
    density = Fraction(LOCAL_SCALED_MASS, 4 * P)
    influence_floor = Fraction((P + 1) * (P - 3), 16 * P * (P - 2))
    total_influence_upper = (P - 1) * density * (1 - density)
    junta_bound = Fraction(
        2
        * (P - 1)
        * (P - 2)
        * (P + 13)
        * (3 * P - 13),
        P * P * (P + 1) * (P - 3),
    )
    densities = sorted(
        {
            profile_density(tuple(row["layer_counts"]), P)
            for row in exact_four_cube_catalog()["profiles"]
        }
    )
    boolean_excluded = bool(
        junta_bound == Fraction(19404, 2645)
        and junta_bound < 8
        and 7 < Q
        and density == Fraction(9, 23)
        and density not in densities
    )

    stabilizer_upper = Fraction(P + 13, 4)
    paired_average_at_three = Fraction(3 + 9, P + 1)
    paired_average_at_nine = Fraction(9 + 9, P + 1)
    height_at_least_two_proved = bool(
        half["proved"]
        and three_quarter["proved"]
        and stabilizer_upper == 9
        and paired_average_at_three == Fraction(1, 2)
        and paired_average_at_nine == Fraction(3, 4)
        and int(half["maximum_upper_bound"]) == 3
        and int(three_quarter["maximum_upper_bound"]) == 6
    )
    proved = bool(boolean_excluded and height_at_least_two_proved)
    _require(proved, "the p23 local height equality changed")
    return {
        "p": P,
        "slice": "J(23,12)",
        "scaled_mass_4p_E_C": LOCAL_SCALED_MASS,
        "mean": str(density),
        "height_one_exclusion": {
            "relevant_pair_influence_floor": str(influence_floor),
            "total_influence_upper_bound": str(total_influence_upper),
            "largest_zero_class_complement_bound": str(junta_bound),
            "junta_coordinates_at_most": 7,
            "four_bit_density_values": [str(value) for value in densities],
            "target_density_absent": True,
            "excluded": boolean_excluded,
        },
        "height_at_least_two": {
            "paired_cube_operator": "T C(X)=(C(X)+23*E[C])/24=(H+9)/24",
            "paired_cube_mean_lattice": "(1/4)Z",
            "every_paired_cube_mean_at_least": "1/2",
            "initial_height_lower_bound": 3,
            "stabilizer_height_upper_bound": str(stabilizer_upper),
            "if_H_at_least_four": (
                "half-mean cubes are impossible, so every cube has mean at "
                "least 3/4 and H>=9"
            ),
            "H_equals_nine_would_force_every_cube_mean": "3/4",
            "three_quarter_cube_maximum_upper_bound": int(
                three_quarter["maximum_upper_bound"]
            ),
            "therefore_H_at_least_four_is_impossible": True,
        },
        "forced_height": 3,
        "paired_cube_average_at_forced_height": str(paired_average_at_three),
        "every_paired_cube_through_a_maximizer_has_mean": "1/2",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_slice_half_mean_classification() -> dict[str, object]:
    """Globalize all paired-cube equality forms on ``J(23,12)``."""
    cube = half_mean_height_three_cube_classification()
    local = p23_local_height_equality()

    f012 = tuple(_f_value(s) for s in range(3))
    rectangle_rows = []
    for diagonal in product((0, 1), repeat=2):
        for off_diagonal in product((0, 1), repeat=2):
            same_value = _f_value(sum(diagonal)) == _f_value(sum(off_diagonal))
            same_sum = sum(diagonal) == sum(off_diagonal)
            rectangle_rows.append(same_value == same_sum)

    # A binary 2-by-2 matrix with additive determinant zero is necessarily
    # constant along its rows or along its columns.  This is the local step
    # in the dimension-free classification: if one row of a larger matrix is
    # nonconstant, every other row equals it; otherwise all rows are constant.
    additive_binary_rectangles = []
    for a, b, c, d in product((0, 1), repeat=4):
        if a + d != b + c:
            continue
        row_constant = a == b and c == d
        column_constant = a == c and b == d
        additive_binary_rectangles.append(
            {
                "entries_row_major": [a, b, c, d],
                "row_constant": row_constant,
                "column_constant": column_constant,
            }
        )
    binary_rectangle_classification = bool(
        len(additive_binary_rectangles) == 6
        and all(
            row["row_constant"] or row["column_constant"]
            for row in additive_binary_rectangles
        )
    )

    # A twelve-set Z differs from X in d=0,...,11 points because the
    # complement of X has size eleven.  Retain one of the 12-d common points,
    # match the d removed points to the d added points, and extend that partial
    # bijection on the remaining equal-size sets.  This certifies that every Z
    # lies in at least one eleven-dimensional paired cube through X.
    coverage_rows = []
    for swaps in range(Q + 1):
        retained = M - swaps
        unused_rows = M - 1 - swaps
        unused_columns = Q - swaps
        coverage_rows.append(
            {
                "swap_count": swaps,
                "common_points_available": retained,
                "unused_rows_after_retaining_one": unused_rows,
                "unused_columns": unused_columns,
                "partial_matching_extends": bool(
                    retained >= 1 and unused_rows == unused_columns >= 0
                ),
            }
        )
    every_slice_point_covered = all(
        row["partial_matching_extends"] for row in coverage_rows
    )

    column_counts = [count for count in range(12) if count in (4, 5)]
    row_counts = []
    for count in range(13):
        near_perfect_counts = set()
        if count:
            near_perfect_counts.add(count - 1)
        if count < 12:
            near_perfect_counts.add(count)
        if near_perfect_counts and near_perfect_counts <= {4, 5}:
            row_counts.append(count)

    forms = []
    for active in (4, 5):
        values = tuple(_f_value(s) for s in range(active + 1))
        numerator = sum(
            comb(active, s) * comb(P - active, M - s) * values[s]
            for s in range(active + 1)
        )
        mean = Fraction(numerator, comb(P, M))
        constant = 15 - 4 * active + comb(active, 2)
        linear = active - 5
        offset = constant + active * linear
        linear_term = (
            "-sum_(i in R)z_i"
            if linear == -1
            else "0"
        )
        forms.append(
            {
                "name": f"F{active}",
                "active_coordinates": active,
                "layer_values": list(values),
                "slice_mean": str(mean),
                "scaled_mass_4p_E": int(4 * P * mean),
                "signed_target_constant": constant,
                "signed_target_linear_coefficient": linear,
                "signed_target_pair_coefficients": 1,
                "coefficient_offset": offset,
                "signed_target": (
                    f"3+4C={constant}+({linear_term})+"
                    "sum_({i,j} subset R)z_i*z_j"
                ),
            }
        )

    compatibility = []
    for hard_parallel in range(2, 6):
        opposite_Q = 9 - hard_parallel
        for form in forms:
            offset = int(form["coefficient_offset"])
            if (opposite_Q - offset) % Q == 0:
                compatibility.append(
                    {
                        "hard_P": hard_parallel,
                        "opposite_Q": opposite_Q,
                        "form": form["name"],
                        "offset": offset,
                    }
                )

    proved = bool(
        cube["proved"]
        and local["proved"]
        and f012 == (3, 1, 0)
        and all(rectangle_rows)
        and binary_rectangle_classification
        and column_counts == [4, 5]
        and row_counts == [5]
        and every_slice_point_covered
        and [form["slice_mean"] for form in forms] == ["9/23", "9/23"]
        and [form["coefficient_offset"] for form in forms] == [1, 5]
        and compatibility
        == [{"hard_P": 4, "opposite_Q": 5, "form": "F5", "offset": 5}]
    )
    _require(proved, "the p23 slice half-mean classification changed")
    return {
        "p": P,
        "slice": "J(23,12)",
        "cube_equality_dependency": cube,
        "cross_activity_matrix": (
            "D_xy=1 iff C(X-x+y)=1 for x in X and y outside X"
        ),
        "same_double_swap_two_pairings": (
            "F(D_xu+D_yv)=F(D_xv+D_yu); F(0),F(1),F(2)=(3,1,0) "
            "are distinct"
        ),
        "all_additive_two_by_two_minors_vanish": True,
        "additive_binary_rectangle_catalog": additive_binary_rectangles,
        "additive_binary_rectangle_classification_proved": (
            binary_rectangle_classification
        ),
        "binary_additive_matrix_classification": ["column-only", "row-only"],
        "active_column_counts": column_counts,
        "active_row_counts": row_counts,
        "paired_cube_coverage_by_swap_count": coverage_rows,
        "every_slice_point_lies_in_a_paired_cube_through_X": (
            every_slice_point_covered
        ),
        "row_only_symmetry": "F(5-r)=F(r) for r=0,...,5",
        "global_slice_forms": forms,
        "opposite_coefficient_congruence": "11 divides Q-offset",
        "compatible_hard_opposite_rows": compatibility,
        "all_forced_local_equality_forms_enumerated": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_first_post_band_ledger() -> dict[str, object]:
    """Replay all branches at ``p=23,t=9,k=110``."""
    floors = residual_even_floor_table(P)
    phase_one = floors["phase_one_floors"]
    phase_zero = floors["phase_zero_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"])
    low_b = sorted(int(b) for b, floor in phase_one.items() if int(floor) < 2 * P)

    residue_rows = []
    for u in range(M):
        quotient_sum = M + LAYER_INDEX - u
        if u < LAYER_INDEX:
            status = "excluded below phase-one or sharp-lift floor"
        elif u == LAYER_INDEX:
            status = NEW_BRANCH
        elif u == M - 2:
            status = "excluded because the quotient sum is below twelve"
        elif u == M - 1:
            status = "old endpoint A/C dichotomy"
        else:  # pragma: no cover
            raise ArithmeticError("unexpected p23 residue interval")
        residue_rows.append(
            {
                "u": u,
                "quotient_sum": quotient_sum,
                "k_zero_mean": 2 * u,
                "k_one_mean": P + 1 + 2 * u,
                "status": status,
            }
        )

    old_branches = []
    for branch, hard_edges, hT, minimum_Q in (
        (BRANCH_B2, 5 * M - 2, 5, 3),
        (BRANCH_P3_LAST, 4 * M - 2, 4 - P, 4),
    ):
        opposite_edges = H_EDGE_COUNT - hard_edges
        next_Q = minimum_Q + 1
        minimum_mean = (P + 1) * minimum_Q + hT - 3 * P
        next_mean = (P + 1) * next_Q + hT - 3 * P
        surplus = opposite_edges - M * next_Q
        old_branches.append(
            {
                "branch": branch,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "hard_sign_times_global_T": hT,
                "minimum_Q": minimum_Q,
                "minimum_mean": minimum_mean,
                "forced_next_Q": next_Q,
                "forced_next_scaled_mean": next_mean,
                "surplus_after_next_Q": surplus,
                "directions_at_next_Q_at_least": M - surplus,
                "dependency": "Proposition 15.752 p+9 local exclusion",
                "excluded": bool(
                    minimum_mean == 8
                    and next_mean == P + 9
                    and 0 <= surplus < M
                    and p_plus_nine_local_exclusion(P)["proved"]
                ),
            }
        )

    hard_catalog = p23_sharp_hard_family_catalog()
    slice_forms = p23_slice_half_mean_classification()
    new_rows = []
    for family in hard_catalog["hard_families"]:
        hard_parallel = int(family["coefficient_offset"])
        hard_edges = M * hard_parallel
        opposite_edges = H_EDGE_COUNT - hard_edges
        hT = (P + 1) * hard_parallel - 5 * P + 4
        mass_twelve_Q = 8 - hard_parallel
        forced_Q = 9 - hard_parallel
        forced_mean = (P + 1) * forced_Q + hT - 3 * P
        surplus = opposite_edges - M * forced_Q
        nonzero_floor_rows = [
            [int(b), int(floor), forced_mean - int(floor)]
            for b, floor in phase_zero.items()
            if int(b) != 0 and int(floor) <= forced_mean
        ]
        compatible_forms = [
            row["form"]
            for row in slice_forms["compatible_hard_opposite_rows"]
            if int(row["hard_P"]) == hard_parallel
            and int(row["opposite_Q"]) == forced_Q
        ]
        new_rows.append(
            {
                **family,
                "forced_hard_parallel_count": hard_parallel,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "hard_sign_times_global_T": hT,
                "mass_twelve_Q": mass_twelve_Q,
                "mass_twelve_excluded": bool(0 < 12 < min(P + 1, lift_floor)),
                "forced_low_Q": forced_Q,
                "forced_low_scaled_mean": forced_mean,
                "surplus_after_every_Q_at_least_forced_low_Q": surplus,
                "directions_at_forced_low_Q_at_least": M - surplus,
                "nonzero_b_floor_and_lift_rows": nonzero_floor_rows,
                "forced_low_Q_is_b_zero": True,
                "forced_local_quadratic": "A=2C and 92*E[C]=36",
                "compatible_slice_forms": compatible_forms,
                "excluded_by_local_offset": not compatible_forms,
            }
        )

    surviving_new_rows = [row for row in new_rows if row["compatible_slice_forms"]]
    proved = bool(
        floors["proved"]
        and low_b == [2, P - 1]
        and lift_floor == P - 3 == 20
        and P * P + 1 - 2 * H_EDGE_COUNT == 308
        and [row["u"] for row in residue_rows if row["status"] in (NEW_BRANCH, "old endpoint A/C dichotomy")]
        == [9, 11]
        and all(bool(row["excluded"]) for row in old_branches)
        and [int(row["forced_hard_parallel_count"]) for row in new_rows]
        == [2, 4, 3, 5]
        and all(int(row["forced_low_scaled_mean"]) == 36 for row in new_rows)
        and all(
            int(row["surplus_after_every_Q_at_least_forced_low_Q"]) == 3
            for row in new_rows
        )
        and all(
            int(row["directions_at_forced_low_Q_at_least"]) == 9
            for row in new_rows
        )
        and all(
            row["nonzero_b_floor_and_lift_rows"] == [[2, 24, 12], [22, 24, 12]]
            for row in new_rows
        )
        and len(surviving_new_rows) == 1
        and int(surviving_new_rows[0]["forced_hard_parallel_count"]) == 4
        and surviving_new_rows[0]["compatible_slice_forms"] == ["F5"]
    )
    _require(proved, "the p23 first post-band ledger changed")
    return {
        "p": P,
        "q": Q,
        "m": M,
        "layer_index_t": LAYER_INDEX,
        "original_k": ORIGINAL_K,
        "H_edge_count": H_EDGE_COUNT,
        "guaranteed_isolated_vertices": P * P + 1 - 2 * H_EDGE_COUNT,
        "phase_one_mean_form": "a_L=2u+24*k_L",
        "phase_one_quotient_sum": "sum k_L=21-u",
        "phase_one_low_b_values": low_b,
        "sharp_integral_lift_floor": lift_floor,
        "residue_rows": residue_rows,
        "surviving_residues": [9, 11],
        "old_branch_exclusions": old_branches,
        "new_sharp_branch_family_ledgers": new_rows,
        "unique_new_survivor_before_moments": {
            "hard_P": 4,
            "opposite_Q": 5,
            "hard_family": "complement-literal plus all-equal-triple",
            "opposite_form": "F5",
            "opposite_Q5_directions_at_least": 9,
        },
        "proved": proved,
    }


def _edge_power_sum(vertices: tuple[int, ...], degree: int) -> int:
    return sum(
        pow((left - right) % P, degree, P)
        for left, right in combinations(vertices, 2)
    ) % P


def _k5_g_values(vertices: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]]:
    moments = tuple(_edge_power_sum(vertices, degree) for degree in (2, 4, 6, 8))
    s2, s4, s6, s8 = moments
    g4 = (-2 * s4 - s2 * s2) % P
    g8 = (-24 * s8 - 32 * s2 * s6 + 5 * pow(s2, 4, P)) % P
    return g4, g8, moments


@lru_cache(maxsize=1)
def p23_hard_moment_root_certificate() -> dict[str, object]:
    """Prove the hard triangle-minus-star rows annihilate ``G4`` and ``G8``."""
    star_sums = {
        degree: sum(pow(value, degree, P) for value in range(1, P)) % P
        for degree in (2, 4, 6, 8)
    }
    triangle_checks = []
    for x in range(P):
        for y in range(P):
            moments = tuple(
                (pow(x, degree, P) + pow(y, degree, P) + pow(x - y, degree, P))
                % P
                for degree in (2, 4, 6, 8)
            )
            s2, s4, s6, s8 = moments
            triangle_checks.append(
                (
                    (2 * s4 - s2 * s2) % P,
                    (24 * s8 - 32 * s2 * s6 + 5 * pow(s2, 4, P)) % P,
                )
            )

    hard_direction_count = M
    g4_degree = 4
    g8_degree = 8
    proved = bool(
        star_sums == {2: 0, 4: 0, 6: 0, 8: 0}
        and len(triangle_checks) == P * P
        and set(triangle_checks) == {(0, 0)}
        and hard_direction_count > g4_degree
        and hard_direction_count > g8_degree
    )
    _require(proved, "the p23 hard moment root certificate changed")
    return {
        "field": "F_23",
        "global_even_moments": (
            "M_d(L)=sum_({u,v} in H)chi(u-v)*(L(u)-L(v))^d"
        ),
        "hard_target": (
            "5-z_j+sum_({a,b} subset T)z_a*z_b = "
            "4 + triangle(T) - full_star(j) on sum z_i=1"
        ),
        "full_star_power_sums_degrees_2_4_6_8": {
            str(degree): value for degree, value in star_sums.items()
        },
        "full_star_vanishes_because": (
            "sum_(a in F_23^*) a^d=0 for 0<d<22"
        ),
        "triangle_pairs_checked": len(triangle_checks),
        "triangle_G4_identity": "2*S4-S2^2=0",
        "triangle_G8_identity": "24*S8-32*S2*S6+5*S2^4=0",
        "hard_sign": "h in {+1,-1}, common on the twelve hard directions",
        "homogeneous_forms": {
            "G4": "2*h*M4-M2^2",
            "G8": "24*h*M8-32*M2*M6+5*M2^4",
        },
        "form_degrees": {"G4": g4_degree, "G8": g8_degree},
        "distinct_hard_projective_roots": hard_direction_count,
        "both_forms_identically_zero": True,
        "opposite_F5_evaluations": {
            "G4": "-2*S4-S2^2",
            "G8": "-24*S8-32*S2*S6+5*S2^4",
        },
        "proved": proved,
    }


def _affine_image(vertices: tuple[int, ...], scale: int, shift: int) -> tuple[int, ...]:
    return tuple(sorted((scale * value + shift) % P for value in vertices))


def _rows_sha256(rows: list[tuple[int, ...]]) -> str:
    payload = "\n".join(",".join(str(value) for value in row) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def p23_k5_moment_sieve() -> dict[str, object]:
    """Exact all-five-set and affine-orbit replay of the two moment forms."""
    all_sets = list(combinations(range(P), 5))
    g4_zero: list[tuple[int, ...]] = []
    g8_zero: list[tuple[int, ...]] = []
    for vertices in all_sets:
        g4, g8, _ = _k5_g_values(vertices)
        if g4 == 0:
            g4_zero.append(vertices)
        if g8 == 0:
            g8_zero.append(vertices)
    intersection = sorted(set(g4_zero) & set(g8_zero))

    remaining = set(all_sets)
    orbit_rows = []
    while remaining:
        representative = min(remaining)
        orbit = {
            _affine_image(representative, scale, shift)
            for scale in range(1, P)
            for shift in range(P)
        }
        canonical = min(orbit)
        _require(canonical == representative, "the affine orbit representative changed")
        g4, g8, moments = _k5_g_values(representative)
        orbit_rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(orbit),
                "moments_S2_S4_S6_S8": list(moments),
                "G4": g4,
                "G8": g8,
            }
        )
        remaining.difference_update(orbit)

    orbit_size_histogram = Counter(int(row["orbit_size"]) for row in orbit_rows)
    g4_zero_orbits = [
        row for row in orbit_rows if int(row["G4"]) == 0
    ]
    g8_zero_orbits = [
        row for row in orbit_rows if int(row["G8"]) == 0
    ]
    orbit_digest_payload = "\n".join(
        ",".join(str(value) for value in row["representative"])
        + f":{row['orbit_size']}"
        for row in orbit_rows
    )
    orbit_digest = hashlib.sha256(orbit_digest_payload.encode("ascii")).hexdigest()

    expected_g4_reps = [
        [0, 1, 2, 3, 12],
        [0, 1, 2, 4, 15],
        [0, 1, 2, 7, 17],
    ]
    expected_g8_reps = [
        [0, 1, 2, 3, 10],
        [0, 1, 2, 4, 17],
        [0, 1, 2, 4, 18],
        [0, 1, 2, 7, 10],
    ]
    proved = bool(
        len(all_sets) == comb(P, 5) == 33649
        and len(g4_zero) == 1518
        and len(g8_zero) == 2024
        and not intersection
        and _rows_sha256(g4_zero)
        == "82460f67f3414a1f461b24605c108861d215f970063c0d0af82772de21240c1a"
        and _rows_sha256(g8_zero)
        == "733bc62c7ad8d0d7083388480d307ad7298d56b4f9e1fcd12562848350c8d6c7"
        and len(orbit_rows) == 69
        and orbit_size_histogram == Counter({506: 64, 253: 5})
        and sum(int(row["orbit_size"]) for row in orbit_rows) == len(all_sets)
        and orbit_digest
        == "34eeb59b625d24907758658f78c0f966291728a72cebb0426a3d4a883fb2022a"
        and [row["representative"] for row in g4_zero_orbits] == expected_g4_reps
        and [row["representative"] for row in g8_zero_orbits] == expected_g8_reps
        and sum(int(row["orbit_size"]) for row in g4_zero_orbits) == len(g4_zero)
        and sum(int(row["orbit_size"]) for row in g8_zero_orbits) == len(g8_zero)
    )
    _require(proved, "the p23 K5 moment sieve changed")
    return {
        "field": "F_23",
        "five_sets_checked": len(all_sets),
        "moment_definition": "S_d(R)=sum_({i,j} subset R)(i-j)^d mod 23",
        "G4": "-2*S4-S2^2",
        "G8": "-24*S8-32*S2*S6+5*S2^4",
        "G4_zero_count": len(g4_zero),
        "G8_zero_count": len(g8_zero),
        "simultaneous_zero_count": len(intersection),
        "G4_zero_sets_sha256": _rows_sha256(g4_zero),
        "G8_zero_sets_sha256": _rows_sha256(g8_zero),
        "affine_group": "AGL(1,23), x maps to a*x+b with a nonzero",
        "affine_orbit_count": len(orbit_rows),
        "affine_orbit_size_histogram": {
            str(size): count for size, count in sorted(orbit_size_histogram.items())
        },
        "affine_orbit_representatives_sha256": orbit_digest,
        "G4_zero_orbits": g4_zero_orbits,
        "G8_zero_orbits": g8_zero_orbits,
        "zero_orbit_representatives_disjoint": not (
            {tuple(row["representative"]) for row in g4_zero_orbits}
            & {tuple(row["representative"]) for row in g8_zero_orbits}
        ),
        "independent_accelerator_replay": {
            "script": "scripts/p23_k5_moment_gpu.py",
            "script_sha256": (
                "4afa5d397ccf38dc6e61f8b006c7a697c533ff067b9c7af944441f35c986ee1d"
            ),
            "reported_backends": [
                "NUKA gfx1201 OpenCL",
                "Jellyfin Arc A380 OpenCL",
                "Soulkiller V100 CUDA",
            ],
            "each_reported_counts": [33649, 1518, 2024, 0],
            "authoritative_certificate": False,
            "purpose": "independent implementation replay",
        },
        "classification": "exhaustive finite coefficient certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_post_band_moment_close() -> dict[str, object]:
    """Package the complete exceptional endpoint proof."""
    cube = half_mean_height_three_cube_classification()
    height = p23_local_height_equality()
    slice_classification = p23_slice_half_mean_classification()
    ledger = p23_first_post_band_ledger()
    roots = p23_hard_moment_root_certificate()
    sieve = p23_k5_moment_sieve()
    proved = bool(
        cube["proved"]
        and height["proved"]
        and slice_classification["proved"]
        and ledger["proved"]
        and roots["proved"]
        and sieve["proved"]
        and int(sieve["simultaneous_zero_count"]) == 0
    )
    _require(proved, "the p23 post-band endpoint did not close")
    return {
        "title": "p23 first post-band equality and K5 moment close",
        "result_status": (
            "proved endpoint theorem with exhaustive finite coefficient certificate"
        ),
        "statement": (
            "the residual-(ii) isolated-chart branch at p=23,t=9,k=110 "
            "is empty for every boundary size"
        ),
        "changed_premise": (
            "all half-mean height-three equality forms are classified and "
            "globalized before the simultaneous degree-four/eight K5 sieve"
        ),
        "half_mean_cube_equality_classification": cube,
        "local_mass_36_height_equality": height,
        "slice_equality_globalization": slice_classification,
        "endpoint_residue_and_parallel_ledger": ledger,
        "hard_triangle_minus_star_root_certificate": roots,
        "opposite_K5_moment_sieve": sieve,
        "p23_k110_closed": proved,
        "all_boundary_sizes_excluded": proved,
        "finite_graph_or_residual_configuration_census_used": False,
        "fixed_five_set_coefficient_certificate_used": True,
        "later_layers_closed": False,
        "residual_ii_closed_globally": False,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_p23_post_band_moment_close.json"
    write_json_atomic(path, p23_post_band_moment_close())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
