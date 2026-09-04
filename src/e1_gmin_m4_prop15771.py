#!/usr/bin/env python3
r"""Prop. 15.771 -- third exceptional ``p=23`` post-band endpoint close.

The general-slice equality, covering swap-cube, and phase-zero mass-32
bridges are explicit and independently reviewed. See the standalone proof
NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md and its two local proof notes.
The four-node arithmetic checks supplement, and do not replace, that proof.

At ``p=23,t=11,k=114`` the isolated-chart phase-one ledger has only three
non-arithmetic residues.  The carried sharp ``u=9`` residue still has ten
low roots, so the fixed quartic/octic five-set certificate used at the two
preceding endpoints excludes it.  The ``u=10`` residue contains a forbidden
scaled-mass-22 lift.

The genuinely new residue is ``u=11``.  All twelve hard rows have scaled
mean 46.  Exact positive contact quadratures and the degree-two slice ideal
classify every equality cell:

If a quotient-zero row occurs, it is one of the two old exact parity
baselines and already fixes the common row ledger.  Otherwise every hard
quotient equals one, and:

* ``b=0`` is constant one;
* ``b=2,22`` is an exact parity baseline plus a Boolean mass-24 lift;
* ``b=4`` has the two ``4000``/``2200`` equality types;
* ``b=6,8,...,18`` is impossible by the injective even-half evaluation of
  degree-two functions on a cube of dimension at least five;
* ``b=20`` has four three-bit equality types.

Their signed coefficient offsets are exactly ``4,5,6,7,8``.  Common-row
normalization makes the hard parallel count ``P`` common.  Since ``12P``
is at most 115, the slice-kernel congruence forces ``P`` to equal its offset,
so different offset classes cannot mix.  For every remaining class,

    hT=24P-115,
    sum Q=115-12P,
    a(Q)=24(P+Q)-184.

The mass-eight row ``Q=8-P`` is impossible.  Assigning the next value
``Q=9-P`` to all twelve opposite rows leaves surplus seven, hence at least
five such rows occur.  They have scaled mass 32, exactly the ``p+9`` local
cell excluded by Proposition 15.752.  This closes the endpoint for every
boundary size.  It is not a graph census, a later-layer theorem, or a global
closure of residual (ii).
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

from e1_gmin_m4_p23_post_band_moment_close import (
    p23_hard_moment_root_certificate,
    p23_k5_moment_sieve,
    p23_sharp_hard_family_catalog,
    p23_slice_half_mean_classification,
)
from e1_gmin_m4_p23_second_post_band_moment_close import (
    p23_p_minus_one_local_exclusion,
)
from e1_gmin_m4_prop15632 import hypergeometric_weights
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P3_LAST,
    baseline_coefficient_rules,
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
P = 23
Q = 11
M = 12
LAYER_INDEX = 11
ORIGINAL_K = 114
H_EDGE_COUNT = 115
MODULUS = 1_000_003


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _rank_mod(rows: list[list[int]], modulus: int = MODULUS) -> int:
    """Exact row rank over a fixed prime field, for the tiny cube checks."""
    if not rows:
        return 0
    matrix = [[value % modulus for value in row] for row in rows]
    columns = len(matrix[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], modulus - 2, modulus)
        matrix[rank] = [value * inverse % modulus for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % modulus
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == columns:
            break
    return rank


def _degree_two_even_half_rank(dimension: int) -> tuple[int, int]:
    """Rank degree-at-most-two cube monomials on even-parity vertices."""
    monomials = [0]
    monomials.extend(1 << index for index in range(dimension))
    monomials.extend(
        (1 << left) | (1 << right)
        for left, right in combinations(range(dimension), 2)
    )
    rows = [
        [int(mask & monomial == monomial) for monomial in monomials]
        for mask in range(1 << dimension)
        if mask.bit_count() % 2 == 0
    ]
    return len(monomials), _rank_mod(rows)


@lru_cache(maxsize=1)
def p23_third_post_band_residue_ledger() -> dict[str, object]:
    """Classify every phase-one residue at ``p=23,t=11``."""
    floors = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(P)[
            "phase_one_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"]
    )
    expected_quotient_one_live = {
        M - 3: [(2, "sharp_p_minus_3"), (P - 1, "sharp_p_minus_3")],
        M - 2: [(2, "p_minus_one"), (P - 1, "p_minus_one")],
        M - 1: [
            (0, "exact_mean_46"),
            (2, "p_plus_one"),
            *[(boundary, "exact_mean_46") for boundary in range(4, 22, 2)],
            (P - 1, "p_plus_one"),
        ],
    }
    expected_quotient_zero_live = {
        M - 1: [(2, "exact"), (P - 1, "exact")],
    }

    rows = []
    for residue in range(M):
        quotient_sum = M + LAYER_INDEX - residue
        low_quotient = 1
        forced_low_count_if_no_zero = 2 * M - quotient_sum
        low_mean = P + 1 + 2 * residue
        zero_mean = 2 * residue
        zero_candidates = []
        zero_live = []
        for boundary, floor in floors.items():
            if floor > zero_mean:
                continue
            excess = zero_mean - floor
            if excess == 0:
                classification = "exact"
            elif 0 < excess < lift_floor:
                classification = "excluded_sub_sharp_lift"
            else:  # pragma: no cover - guarded by the p23 floor table
                raise ArithmeticError(
                    f"unclassified p23 quotient-zero excess at "
                    f"u={residue},b={boundary}"
                )
            zero_candidates.append(
                {
                    "b": boundary,
                    "floor": floor,
                    "excess": excess,
                    "classification": classification,
                }
            )
            if classification != "excluded_sub_sharp_lift":
                zero_live.append((boundary, classification))
        candidates = []
        live = []
        for boundary, floor in floors.items():
            if floor > low_mean:
                continue
            excess = low_mean - floor
            if excess == 0:
                classification = "exact_mean_46"
            elif 0 < excess < lift_floor:
                classification = "excluded_sub_sharp_lift"
            elif excess == lift_floor:
                classification = "sharp_p_minus_3"
            elif excess == P - 1:
                classification = "p_minus_one"
            elif excess == P + 1:
                classification = "p_plus_one"
            else:  # pragma: no cover - guarded by the exact floor table
                raise ArithmeticError(
                    f"unclassified p23 excess at u={residue},b={boundary}"
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
            live == expected_quotient_one_live.get(residue, [])
            and zero_live == expected_quotient_zero_live.get(residue, [])
            and forced_low_count_if_no_zero > 0,
            f"the p23 third-post-band residue u={residue} changed",
        )
        rows.append(
            {
                "u": residue,
                "quotient_sum": quotient_sum,
                "quotient_zero_mean": zero_mean,
                "quotient_zero_candidate_rows": zero_candidates,
                "quotient_zero_live_rows": [
                    {"b": boundary, "classification": classification}
                    for boundary, classification in zero_live
                ],
                "forced_low_quotient": low_quotient,
                "forced_quotient_one_count_if_no_quotient_zero": (
                    forced_low_count_if_no_zero
                ),
                "forced_low_mean": low_mean,
                "candidate_rows": candidates,
                "live_rows": [
                    {"b": boundary, "classification": classification}
                    for boundary, classification in live
                ],
                "excluded_arithmetically": not live,
            }
        )

    isolated = P * P + 1 - 2 * H_EDGE_COUNT
    proved = bool(
        H_EDGE_COUNT == 4 * P + 2 * LAYER_INDEX + 1
        and isolated == 300 > 0
        and lift_floor == P - 3 == 20
        and [row["u"] for row in rows if row["live_rows"]] == [9, 10, 11]
    )
    _require(proved, "the p23 third-post-band residue ledger failed")
    return {
        "p": P,
        "q": Q,
        "m": M,
        "layer_index_t": LAYER_INDEX,
        "original_k": ORIGINAL_K,
        "H_edge_count": H_EDGE_COUNT,
        "guaranteed_isolated_vertices": isolated,
        "phase_one_mean_form": "a_L=2u+24*k_L",
        "phase_one_quotient_sum": "sum k_L=23-u",
        "u11_dichotomy": (
            "either an exact quotient-zero b=2/22 row occurs, or all twelve "
            "hard quotients equal one"
        ),
        "sharp_lift_floor": lift_floor,
        "arithmetic_surviving_residues": [9, 10, 11],
        "rows": rows,
        "proved": proved,
    }


# Each row is an exact degree-two quadrature for the hypergeometric law.
# All listed nodes have phase-one parity one.  For b=4,...,18 every even
# layer has positive weight; at b=20 the two positive layers become r=2,0
# after passing to the three-point complement.
_CONTACT_QUADRATURES: dict[int, tuple[tuple[int, ...], tuple[Fraction, ...]]] = {
    4: ((0, 2, 4), (Fraction(2, 23), Fraction(18, 23), Fraction(3, 23))),
    6: (
        (0, 2, 4, 6),
        (Fraction(1, 92), Fraction(39, 92), Fraction(51, 92), Fraction(1, 92)),
    ),
    8: (
        (0, 2, 4, 6, 8),
        (Fraction(3, 115), Fraction(3, 115), Fraction(93, 115), Fraction(13, 115), Fraction(3, 115)),
    ),
    10: (
        (0, 2, 4, 6, 8, 10),
        (Fraction(3, 368), Fraction(3, 368), Fraction(3, 8), Fraction(109, 184), Fraction(3, 368), Fraction(3, 368)),
    ),
    12: (
        (2, 4, 6, 8, 10, 12),
        (Fraction(3, 184), Fraction(3, 184), Fraction(77, 92), Fraction(9, 92), Fraction(3, 184), Fraction(3, 184)),
    ),
    14: (
        (4, 6, 8, 10, 12),
        (Fraction(3, 230), Fraction(83, 230), Fraction(3, 5), Fraction(3, 230), Fraction(3, 230)),
    ),
    16: (
        (6, 8, 10, 12),
        (Fraction(1, 23), Fraction(18, 23), Fraction(3, 23), Fraction(1, 23)),
    ),
    18: (
        (8, 10, 12),
        (Fraction(15, 46), Fraction(15, 23), Fraction(1, 46)),
    ),
    20: ((10, 12), (Fraction(18, 23), Fraction(5, 23))),
}


@lru_cache(maxsize=1)
def mean_46_contact_quadratures() -> dict[str, object]:
    """Verify positive contact quadratures for every broad floor-46 row."""
    rows = []
    for boundary, (nodes, weights) in _CONTACT_QUADRATURES.items():
        distribution = hypergeometric_weights(P, boundary)
        source_moments = tuple(
            sum(weight * intersection**degree for intersection, weight in distribution.items())
            for degree in range(3)
        )
        quadrature_moments = tuple(
            sum(weight * node**degree for node, weight in zip(nodes, weights))
            for degree in range(3)
        )
        exact = bool(
            len(nodes) == len(weights)
            and all(weight > 0 for weight in weights)
            and all(node in distribution and node % 2 == 0 for node in nodes)
            and source_moments == quadrature_moments
            and source_moments[0] == 1
        )
        _require(exact, f"the b={boundary} contact quadrature changed")
        rows.append(
            {
                "b": boundary,
                "nodes": list(nodes),
                "weights": [str(weight) for weight in weights],
                "moments_0_1_2": [str(value) for value in source_moments],
                "equality_forces_A_equals_one_pointwise_on_nodes": True,
                "proved": exact,
            }
        )
    return {
        "p": P,
        "phase": 1,
        "candidate": "q(t)=1",
        "rows": rows,
        "all_weights_strictly_positive": True,
        "proved": True,
    }


def _truth_table_anf(values: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    coefficients = list(values)
    for bit in range(dimension):
        for mask in range(1 << dimension):
            if mask & (1 << bit):
                coefficients[mask] -= coefficients[mask ^ (1 << bit)]
    return tuple(coefficients)


def _signed_target(values: tuple[int, ...], dimension: int) -> dict[str, object]:
    """Convert ``3+2A`` from 0/1 variables to signed ``z=2x-1`` form."""
    anf = _truth_table_anf(values, dimension)
    _require(
        all(
            coefficient == 0
            for mask, coefficient in enumerate(anf)
            if mask.bit_count() > 2
        ),
        "an equality table ceased to be quadratic",
    )
    constant = 3 + 2 * anf[0]
    linear_binary = [2 * anf[1 << index] for index in range(dimension)]
    pair_binary = {
        (left, right): 2 * anf[(1 << left) | (1 << right)]
        for left, right in combinations(range(dimension), 2)
    }
    signed_constant = constant + Fraction(sum(linear_binary), 2) + Fraction(
        sum(pair_binary.values()), 4
    )
    signed_linear = [
        Fraction(linear_binary[index], 2)
        + Fraction(
            sum(
                value
                for pair, value in pair_binary.items()
                if index in pair
            ),
            4,
        )
        for index in range(dimension)
    ]
    signed_pairs = {
        f"{left},{right}": str(Fraction(value, 4))
        for (left, right), value in pair_binary.items()
        if value
    }
    _require(
        signed_constant.denominator == 1
        and all(value.denominator == 1 for value in signed_linear),
        "the signed equality coefficients ceased to be integral",
    )
    offset = signed_constant + sum(signed_linear)
    _require(offset.denominator == 1, "the equality offset ceased to be integral")
    return {
        "signed_constant": int(signed_constant),
        "signed_linear_coefficients": [int(value) for value in signed_linear],
        "signed_pair_coefficients": signed_pairs,
        "coefficient_offset": int(offset),
    }


def _expected_value_on_slice(values: tuple[int, ...], active: int) -> Fraction:
    distribution = hypergeometric_weights(P, active)
    total = Fraction(0)
    for weight, probability in distribution.items():
        layer = [
            values[mask]
            for mask in range(1 << active)
            if mask.bit_count() == weight
        ]
        total += probability * Fraction(sum(layer), len(layer))
    return total


@lru_cache(maxsize=1)
def mean_46_small_support_equality_catalog() -> dict[str, object]:
    r"""Classify general-slice ``b=0,4,20`` floor-equality cells.

    For a small side S of size3 or4, put F=A-1. Contact r=0 makes the
    pure outside polynomial q(y) vanish on J(23-|S|,12). The quadratic
    fixed-weight kernel identity q=(sum y-12)(lambda+sum mu_j y_j)
    eliminates that polynomial using sum y=12-sum x. At r=2, the mixed
    coefficient pair sums D_ij+D_lj are constant in outside index j.
    Since |S|>=3, every D_ij is itself constant in j, and outside weight
    substitution leaves a quadratic in S alone. All bit patterns extend.
    Thus the catalog is exhaustive for arbitrary slice quadratics, not an
    assumed small junta. The elementary kernel proof is in the associated
    SMALL_BOUNDARY_EQUALITY_PROOF note.
    """
    quadratures = mean_46_contact_quadratures()

    b4_rows = []
    for e_values in product(range(0, 5, 2), repeat=4):
        if sum(e_values) != 4:
            continue
        values = []
        for mask in range(16):
            weight = mask.bit_count()
            if weight in (0, 2, 4):
                value = 1
            elif weight == 1:
                value = e_values[(mask & -mask).bit_length() - 1]
            else:
                omitted = next(index for index in range(4) if not mask & (1 << index))
                value = e_values[omitted]
            values.append(value)
        values_tuple = tuple(values)
        signed = _signed_target(values_tuple, 4)
        _require(
            all(value >= 0 for value in values_tuple)
            and all(
                value % 2 == ((mask.bit_count() + 1) & 1)
                for mask, value in enumerate(values_tuple)
            )
            and _expected_value_on_slice(values_tuple, 4) == 1
            and signed["signed_constant"] == 5
            and signed["signed_linear_coefficients"] == [0, 0, 0, 0]
            and signed["coefficient_offset"] == 5,
            "a b=4 equality form changed",
        )
        b4_rows.append(
            {
                "singleton_values": list(e_values),
                "orbit_type": "4000" if 4 in e_values else "2200",
                **signed,
            }
        )

    b20_rows = []
    for e_values in product(range(0, 5, 2), repeat=3):
        endpoint = 4 - sum(e_values)
        if endpoint < 0 or endpoint % 2:
            continue
        values = []
        for mask in range(8):
            weight = mask.bit_count()
            if weight in (0, 2):
                value = 1
            elif weight == 1:
                value = e_values[(mask & -mask).bit_length() - 1]
            else:
                value = endpoint
            values.append(value)
        values_tuple = tuple(values)
        signed = _signed_target(values_tuple, 3)
        half_sum = sum(e_values) // 2
        expected_offset = 8 - 2 * half_sum
        _require(
            all(value >= 0 for value in values_tuple)
            and all(
                value % 2 == ((mask.bit_count() + 1) & 1)
                for mask, value in enumerate(values_tuple)
            )
            and _expected_value_on_slice(values_tuple, 3) == 1
            and signed["signed_constant"] == 5
            and signed["coefficient_offset"] == expected_offset,
            "a b=20 equality form changed",
        )
        partition = tuple(sorted(e_values, reverse=True))
        b20_rows.append(
            {
                "singleton_values": list(e_values),
                "triple_value": endpoint,
                "orbit_type": f"{''.join(map(str, partition))};{endpoint}",
                **signed,
            }
        )

    b4_orbits = Counter(row["orbit_type"] for row in b4_rows)
    b20_orbits = Counter(row["orbit_type"] for row in b20_rows)
    proved = bool(
        quadratures["proved"]
        and len(b4_rows) == 10
        and b4_orbits == Counter({"2200": 6, "4000": 4})
        and len(b20_rows) == 10
        and b20_orbits
        == Counter({"220;0": 3, "400;0": 3, "200;2": 3, "000;4": 1})
        and {row["coefficient_offset"] for row in b20_rows} == {4, 6, 8}
    )
    _require(proved, "the small-support mean-46 equality catalog failed")
    return {
        "b=0": {
            "pointwise_form": "A=1",
            "signed_target": "3+2A=5",
            "coefficient_offset": 5,
        },
        "slice_ideal_reduction": (
            "the r=0 contact kills the pure outside polynomial; constancy on "
            "every r=2 contact makes each cross-coefficient column constant, "
            "so the equality cell depends only on the small parity side"
        ),
        "general_slice_globalization": {
            "outside_sizes": [19, 20],
            "zero_contact_outside_weight": 12,
            "pair_contact_outside_weight": 10,
            "fixed_weight_quadratic_kernel": "q(y)=(sum y-12)*(lambda+sum mu_j*y_j)",
            "cross_pair_sums_force_each_row_constant": True,
            "all_small_side_bit_patterns_extend": True,
            "proof_note": "evidence/NOTE_2026-09-04_P23_SMALL_BOUNDARY_EQUALITY_PROOF.md",
        },
        "cross_column_pair_sum_rank": {
            "b=4_active_columns": 4,
            "b=4_pair_sum_matrix_rank": _rank_mod(
                [
                    [int(index in pair) for index in range(4)]
                    for pair in combinations(range(4), 2)
                ]
            ),
            "b=20_complement_columns": 3,
            "b=20_pair_sum_matrix_rank": _rank_mod(
                [
                    [int(index in pair) for index in range(3)]
                    for pair in combinations(range(3), 2)
                ]
            ),
        },
        "b=4": {
            "labeled_form_count": len(b4_rows),
            "orbit_histogram": dict(sorted(b4_orbits.items())),
            "forms": b4_rows,
            "coefficient_offsets": [5],
        },
        "b=20": {
            "labeled_form_count": len(b20_rows),
            "orbit_histogram": dict(sorted(b20_orbits.items())),
            "forms": b20_rows,
            "coefficient_offsets": [4, 6, 8],
        },
        "proved": proved,
    }


@lru_cache(maxsize=1)
def middle_boundary_equality_exclusion() -> dict[str, object]:
    """Exclude exact mean 46 for ``b=6,8,...,18``."""
    quadratures = mean_46_contact_quadratures()
    rows = []
    for boundary in range(6, 20, 2):
        active = min(boundary, P - boundary)
        columns, rank = _degree_two_even_half_rank(active)
        row = next(item for item in quadratures["rows"] if item["b"] == boundary)
        forced_original_layers = list(row["nodes"])
        covering_cases = []
        for selected_on_small_side in range(active + 1):
            # Every selected small-side point needs an unselected outside
            # partner; every unselected one needs a selected outside partner.
            outside_selected = M - selected_on_small_side
            outside_unselected = P - active - outside_selected
            covering_cases.append({
                "selected_on_small_side": selected_on_small_side,
                "outside_selected": outside_selected,
                "outside_unselected": outside_unselected,
                "selected_partner_slack": outside_selected - (active - selected_on_small_side),
                "unselected_partner_slack": outside_unselected - selected_on_small_side,
            })
        original_intersections = [
            weight if boundary <= Q else M - weight
            for weight in range(0, active + 1, 2)
        ]
        covers = all(
            case["selected_partner_slack"] == M - active >= 0
            and case["unselected_partner_slack"] == Q - active >= 0
            for case in covering_cases
        )
        proved = bool(
            5 <= active <= Q
            and rank == columns == 1 + active + comb(active, 2)
            and all(pattern % 2 == 0 for pattern in forced_original_layers)
            and sorted(original_intersections) == forced_original_layers
            and covers
        )
        _require(proved, f"the b={boundary} even-half rank failed")
        rows.append(
            {
                "b": boundary,
                "smaller_parity_side_size": active,
                "all_patterns_on_smaller_side_extend_to_J(23,12)": True,
                "every_slice_point_has_cross_boundary_swap_cube": covers,
                "covering_capacity_cases": covering_cases,
                "fixed_selected_outside_count": M - active,
                "cube_original_intersection": "r" if boundary <= Q else "12-r",
                "degree_at_most_two_dimension": columns,
                "even_half_evaluation_rank_mod_1000003": rank,
                "contact_layers_in_original_b_coordinate": forced_original_layers,
                "conclusion": (
                    "A-1 vanishes on the even half, hence identically; A=1 "
                    "then contradicts phase-one parity on the odd half"
                ),
                "excluded": proved,
            }
        )
    return {
        "rows": rows,
        "fourier_reason": (
            "on an even-parity d-cube, characters S and S^c agree; for d>=5 "
            "no two subsets of size at most two are complements, so their "
            "restrictions are linearly independent"
        ),
        "covering_proof": (
            "pair each smaller-side coordinate with a distinct outside coordinate "
            "of opposite membership in the prescribed 12-set; capacities have "
            "slacks 12-d and11-d. Fix the remaining12-d selected coordinates. "
            "Choosing one point from every pair gives a d-cube through the set. "
            "Original boundary intersection has the parity of the cube weight."
        ),
        "proved": all(row["excluded"] for row in rows),
    }


def middle_boundary_swap_cube(boundary: int, selected: tuple[int, ...]) -> dict[str, object]:
    """Construct a full cross-boundary cube through a prescribed 12-set."""
    if boundary not in range(6, 20, 2):
        raise ValueError("need b in {6,8,10,12,14,16,18}")
    chosen = set(selected)
    if len(selected) != M or len(chosen) != M or not chosen <= set(range(P)):
        raise ValueError("selected must contain twelve distinct coordinates in0..22")
    boundary_set = set(range(boundary))
    small = boundary_set if boundary <= Q else set(range(P)) - boundary_set
    outside = set(range(P)) - small
    selected_partners = iter(sorted(outside & chosen))
    unselected_partners = iter(sorted(outside - chosen))
    pairs = tuple(
        (point, next(unselected_partners if point in chosen else selected_partners))
        for point in sorted(small)
    )
    paired = {point for pair in pairs for point in pair}
    fixed = tuple(sorted(chosen - paired))
    original_mask = sum(1 << index for index, pair in enumerate(pairs) if pair[0] in chosen)
    recovered = set(fixed) | {
        pair[0] if original_mask & (1 << index) else pair[1]
        for index, pair in enumerate(pairs)
    }
    _require(
        len(paired) == 2 * len(small)
        and len(fixed) == M - len(small)
        and recovered == chosen,
        "the covering cube lost its base point",
    )
    return {
        "b": boundary,
        "dimension": len(small),
        "pairs": pairs,
        "fixed_selected": fixed,
        "base_point_mask": original_mask,
    }


def _four_bit_anf(table: int) -> tuple[int, ...]:
    values = tuple((table >> mask) & 1 for mask in range(16))
    return _truth_table_anf(values, 4)


def _four_bit_layer_counts(table: int) -> tuple[int, ...]:
    return tuple(
        sum((table >> mask) & 1 for mask in range(16) if mask.bit_count() == weight)
        for weight in range(5)
    )


def _table_from_function(function) -> int:
    return sum(int(function(mask)) << mask for mask in range(16))


@lru_cache(maxsize=1)
def p23_mass_24_lift_catalog() -> dict[str, object]:
    """Classify every nonnegative integral lift with ``92 E[L]=24``."""
    mass = P + 1
    density = Fraction(mass, 4 * P)
    lift = sharp_integral_quadratic_lift_floor(P)
    half = cube_half_mean_height_certificate()

    paired_lower_height = Fraction(2 * (P + 1) - mass, 4)
    stabilizer_upper_height = Fraction(mass, 4)
    forced_height = int(paired_lower_height)
    paired_mean = Fraction(forced_height + mass // 4, P + 1)
    height_at_least_two_excluded = bool(
        lift["proved"]
        and half["proved"]
        and paired_lower_height == stabilizer_upper_height == 6
        and paired_mean == Fraction(1, 2)
        and int(half["maximum_upper_bound"]) == 3 < forced_height
    )

    junta_bound = Fraction(
        2 * (P - 1) * (P - 2) * (3 * P - 1),
        P * P * (P - 3),
    )
    target_tables = []
    for table in range(1 << 16):
        coefficients = _four_bit_anf(table)
        if any(
            coefficient
            for mask, coefficient in enumerate(coefficients)
            if mask.bit_count() > 2
        ):
            continue
        profile = _four_bit_layer_counts(table)
        if profile_density(profile, P) == density:
            target_tables.append(table)

    omitted_pair_tables = {
        _table_from_function(
            lambda mask, left=left, right=right: not (mask & (1 << left))
            and not (mask & (1 << right))
        )
        for left, right in combinations(range(4), 2)
    }
    oriented_pair_tables = {
        _table_from_function(
            lambda mask, left=left, right=right: bool(mask & (1 << left))
            and not (mask & (1 << right))
        )
        for left in range(4)
        for right in range(4)
        if left != right
    }
    compact_triangle_tables = {
        _table_from_function(
            lambda mask, left=left, right=right, apex=apex: (
                int(bool(mask & (1 << left)))
                * int(bool(mask & (1 << right)))
                + int(bool(mask & (1 << apex)))
                * (
                    1
                    - int(bool(mask & (1 << left)))
                    - int(bool(mask & (1 << right)))
                )
            )
        )
        for active in combinations(range(4), 3)
        for apex in active
        for left, right in [tuple(index for index in active if index != apex)]
    }
    classified = {
        "selected_pair_on_original_slice": omitted_pair_tables,
        "oriented_pair_on_original_slice": oriented_pair_tables,
        "compact_triangle": compact_triangle_tables,
    }
    union = set().union(*classified.values())

    # The catalog variables are complementary-slice bits y=1-x.  For 4L,
    # changing from signed w=2y-1 to z=2x-1=-w reverses every linear term.
    family_rows = []
    for family, tables in classified.items():
        increments = set()
        for table in tables:
            anf = _four_bit_anf(table)
            binary_linear = [4 * anf[1 << index] for index in range(4)]
            binary_pairs = {
                (left, right): 4 * anf[(1 << left) | (1 << right)]
                for left, right in combinations(range(4), 2)
            }
            signed_constant = 4 * anf[0] + Fraction(sum(binary_linear), 2) + Fraction(
                sum(binary_pairs.values()), 4
            )
            signed_linear_in_y = [
                Fraction(binary_linear[index], 2)
                + Fraction(
                    sum(value for pair, value in binary_pairs.items() if index in pair),
                    4,
                )
                for index in range(4)
            ]
            increments.add(int(signed_constant - sum(signed_linear_in_y)))
        _require(len(increments) == 1, "a mass-24 family lost its common offset")
        family_rows.append(
            {
                "family": family,
                "table_count": len(tables),
                "four_L_offset_increment_after_slice_complement": next(iter(increments)),
            }
        )

    digest = hashlib.sha256(
        json.dumps(sorted(target_tables), separators=(",", ":")).encode("ascii")
    ).hexdigest()
    proved = bool(
        height_at_least_two_excluded
        and density == Fraction(6, 23)
        and junta_bound == Fraction(15708, 2645) < 6
        and 5 < Q
        and exact_four_cube_catalog()["proved"]
        and len(target_tables) == 30
        and target_tables == sorted(union)
        and [(row["family"], row["table_count"], row["four_L_offset_increment_after_slice_complement"]) for row in family_rows]
        == [
            ("selected_pair_on_original_slice", 6, 3),
            ("oriented_pair_on_original_slice", 12, 1),
            ("compact_triangle", 12, 1),
        ]
    )
    _require(proved, "the p23 mass-24 lift catalog failed")
    return {
        "p": P,
        "scaled_mass_4p_E_L": mass,
        "height_at_least_two": {
            "paired_bound_forces_H_at_least": str(paired_lower_height),
            "stabilizer_bound_forces_H_at_most": str(stabilizer_upper_height),
            "therefore_H": forced_height,
            "every_paired_cube_through_a_maximizer_has_mean": str(paired_mean),
            "half_mean_cube_maximum_upper_bound": int(half["maximum_upper_bound"]),
            "excluded": height_at_least_two_excluded,
        },
        "therefore_boolean": True,
        "density": str(density),
        "corrected_johnson_junta_bound": str(junta_bound),
        "slice_coordinates_at_most": 5,
        "cube_active_coordinates_at_most": 4,
        "target_table_count": len(target_tables),
        "target_tables_sha256": digest,
        "families": family_rows,
        "fixed_four_bit_catalog_reused": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def mean_46_hard_family_catalog() -> dict[str, object]:
    """Collect every exact mean-46 family and its signed offset."""
    small = mean_46_small_support_equality_catalog()
    middle = middle_boundary_equality_exclusion()
    mass24 = p23_mass_24_lift_catalog()
    rows = [
        {"b": 0, "family": "constant_one", "coefficient_offset": 5},
    ]
    rows.extend(
        {
            "b": 2,
            "family": row["family"],
            "baseline_offset": 4,
            "lift_offset_increment": row[
                "four_L_offset_increment_after_slice_complement"
            ],
            "coefficient_offset": 4
            + int(row["four_L_offset_increment_after_slice_complement"]),
        }
        for row in mass24["families"]
    )
    rows.extend(
        {
            "b": 4,
            "family": orbit,
            "coefficient_offset": 5,
        }
        for orbit in sorted(small["b=4"]["orbit_histogram"])
    )
    rows.extend(
        {
            "b": 20,
            "family": orbit,
            "coefficient_offset": next(
                int(form["coefficient_offset"])
                for form in small["b=20"]["forms"]
                if form["orbit_type"] == orbit
            ),
        }
        for orbit in sorted(small["b=20"]["orbit_histogram"])
    )
    rows.extend(
        {
            "b": 22,
            "family": row["family"],
            "baseline_offset": 3,
            "lift_offset_increment": row[
                "four_L_offset_increment_after_slice_complement"
            ],
            "coefficient_offset": 3
            + int(row["four_L_offset_increment_after_slice_complement"]),
        }
        for row in mass24["families"]
    )
    offsets = sorted({int(row["coefficient_offset"]) for row in rows})
    covered = sorted({int(row["b"]) for row in rows} | {6, 8, 10, 12, 14, 16, 18})
    proved = bool(
        small["proved"]
        and middle["proved"]
        and mass24["proved"]
        and offsets == [4, 5, 6, 7, 8]
        and covered == list(range(0, 23, 2))
    )
    _require(proved, "the mean-46 hard family catalog is incomplete")
    return {
        "exact_family_rows": rows,
        "excluded_boundary_values": [6, 8, 10, 12, 14, 16, 18],
        "possible_coefficient_offsets": offsets,
        "all_even_boundary_values_accounted_for": covered,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_u9_two_unit_carry_exclusion() -> dict[str, object]:
    """Exclude the ``u=9`` residue using ten low hard roots."""
    catalog = p23_sharp_hard_family_catalog()
    slice_forms = p23_slice_half_mean_classification()
    root_identities = p23_hard_moment_root_certificate()
    sieve = p23_k5_moment_sieve()
    phase_zero = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(P)[
            "phase_zero_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"]
    )

    quotient_excess = 2
    low_hard_count = M - quotient_excess
    family_rows = []
    for family in catalog["hard_families"]:
        offset = int(family["coefficient_offset"])
        parallel_candidates = [
            value
            for value in range(H_EDGE_COUNT // M + 1)
            if (value - offset) % Q == 0
        ]
        _require(parallel_candidates == [offset], "a u=9 offset changed")
        hard_edges = M * offset + quotient_excess
        opposite_edges = H_EDGE_COUNT - hard_edges
        h_times_T = (P + 1) * offset - 5 * P + 4

        def opposite_mean(parallel: int) -> int:
            return (P + 1) * parallel + h_times_T - 3 * P

        forbidden_Q = 8 - offset
        forced_Q = 9 - offset
        surplus = opposite_edges - M * forced_Q
        nonzero_rows = [
            [boundary, floor, opposite_mean(forced_Q) - floor]
            for boundary, floor in phase_zero.items()
            if boundary and floor <= opposite_mean(forced_Q)
        ]
        compatible_forms = [
            form["name"]
            for form in slice_forms["global_slice_forms"]
            if (forced_Q - int(form["coefficient_offset"])) % Q == 0
        ]
        family_rows.append(
            {
                **family,
                "hard_parallel_candidates": parallel_candidates,
                "hard_quotient_excess_above_all_low": quotient_excess,
                "low_hard_direction_count_at_least": low_hard_count,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "hard_sign_times_global_T": h_times_T,
                "forbidden_Q": forbidden_Q,
                "forbidden_scaled_mass": opposite_mean(forbidden_Q),
                "forced_low_Q": forced_Q,
                "forced_low_scaled_mass": opposite_mean(forced_Q),
                "surplus_after_forced_low_Q": surplus,
                "directions_at_forced_low_Q_at_least": M - surplus,
                "nonzero_boundary_floor_and_lift_rows": nonzero_rows,
                "compatible_slice_forms": compatible_forms,
            }
        )

    survivors = [row for row in family_rows if row["compatible_slice_forms"]]
    maximum_form_degree = max(
        int(value) for value in root_identities["form_degrees"].values()
    )
    roots_force_zero = low_hard_count > maximum_form_degree
    excluded = bool(
        len(survivors) == 1
        and int(survivors[0]["coefficient_offset"]) == 4
        and int(survivors[0]["forced_low_Q"]) == 5
        and survivors[0]["compatible_slice_forms"] == ["F5"]
        and roots_force_zero
        and int(sieve["simultaneous_zero_count"]) == 0
    )
    proved = bool(
        catalog["proved"]
        and slice_forms["proved"]
        and root_identities["proved"]
        and sieve["proved"]
        and [int(row["coefficient_offset"]) for row in family_rows]
        == [2, 4, 3, 5]
        and all(int(row["forbidden_scaled_mass"]) == 12 for row in family_rows)
        and all(int(row["forced_low_scaled_mass"]) == 36 for row in family_rows)
        and all(int(row["surplus_after_forced_low_Q"]) == 5 for row in family_rows)
        and all(int(row["directions_at_forced_low_Q_at_least"]) == 7 for row in family_rows)
        and all(
            row["nonzero_boundary_floor_and_lift_rows"]
            == [[2, 24, 12], [22, 24, 12]]
            for row in family_rows
        )
        and lift_floor == 20
        and excluded
    )
    _require(proved, "the p23 u=9 two-unit carry survived")
    return {
        "residue_u": 9,
        "changed_premise": (
            "the t=9 all-low sharp family now has total quotient excess two, "
            "hence at least ten low hard rows"
        ),
        "family_ledgers": family_rows,
        "unique_survivor_before_moments": {
            "hard_P": 4,
            "opposite_Q": 5,
            "opposite_form": "F5",
        },
        "low_triangle_minus_star_projective_roots_at_least": low_hard_count,
        "maximum_common_form_degree": maximum_form_degree,
        "roots_force_G4_and_G8_identically_zero": roots_force_zero,
        "opposite_K5_simultaneous_zero_count": sieve["simultaneous_zero_count"],
        "fixed_five_set_certificate_reused": True,
        "excluded": excluded,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_phase_zero_mass_exclusions() -> dict[str, object]:
    r"""Exclude opposite scaled means 8 and 32 for every even boundary.

    Write a=2p E[A], with A nonnegative, integral, quadratic and of
    phase-zero parity |B intersect X| on J(23,12). The exact floors leave
    b=0 only at a=8, and b=0,2,22 at a=32. At b=0, A=2L. At b=2 the
    parity minimum is (x_i-x_j)^2; at b=22 it is x_j for the omitted point.
    Both minima are Boolean quadratics of scaled mean24. Hence subtracting
    them and dividing by two gives a NONNEGATIVE integral quadratic L,
    not merely a formal difference. Its mass4p E[L] is the excess.

    Positive mass8 contradicts Proposition15.688's sharp mass20 floor;
    mass32 at b=0 contradicts the local, band-independent Proposition15.752.
    """
    floor_certificate = residual_even_floor_table(P)
    floors = {
        int(boundary): int(value)
        for boundary, value in floor_certificate["phase_zero_floors"].items()
    }
    sharp = sharp_integral_quadratic_lift_floor(P)
    sharp_mass = int(sharp["sharp_scaled_floor"])
    local = p_plus_nine_local_exclusion(P)
    baseline_tables = {
        2: [(left - right) ** 2 for left, right in product((0, 1), repeat=2)],
        22: [0, 1],
    }
    baseline_parities = {
        2: [(left + right) % 2 for left, right in product((0, 1), repeat=2)],
        22: [(M - omitted) % 2 for omitted in (0, 1)],
    }
    baseline_means = {
        2: Fraction(2 * M * (P - M), P * (P - 1)),
        22: Fraction(M, P),
    }
    _require(
        all(baseline_tables[b] == baseline_parities[b] for b in (2, 22))
        and all(2 * P * baseline_means[b] == 24 for b in (2, 22)),
        "the phase-zero pointwise parity baselines changed",
    )
    rows = []
    for mass in (8, 32):
        boundary_rows = []
        for boundary, floor in floors.items():
            if floor > mass:
                reason = "scaled mean below exact phase-zero floor"
                lift_mass = None
                excluded = True
            else:
                _require(boundary in (0, 2, 22), "unclassified phase-zero boundary")
                lift_mass = mass - floor
                if 0 < lift_mass < sharp_mass:
                    reason = "nonzero nonnegative integral quadratic lift below mass20"
                    excluded = bool(sharp["proved"])
                elif boundary == 0 and lift_mass == P + 9:
                    reason = "b=0 gives A=2L, excluded by local p+9 theorem15.752"
                    excluded = bool(local["proved"])
                else:
                    raise ArithmeticError("unclassified phase-zero lift mass")
            boundary_rows.append({
                "b": boundary,
                "phase_zero_floor": floor,
                "lift_mass": lift_mass,
                "reason": reason,
                "excluded": excluded,
            })
        rows.append({
            "scaled_mean": mass,
            "floor_surviving_boundaries": [b for b, floor in floors.items() if floor <= mass],
            "boundary_rows": boundary_rows,
            "excluded": all(row["excluded"] for row in boundary_rows),
        })
    proved = bool(
        floor_certificate["proved"] and sharp["proved"] and local["proved"]
        and sharp_mass == 20
        and rows[0]["floor_surviving_boundaries"] == [0]
        and rows[1]["floor_surviving_boundaries"] == [0, 2, 22]
        and all(row["excluded"] for row in rows)
    )
    _require(proved, "an opposite mass8 or mass32 boundary survived")
    return {
        "p": P,
        "phase": 0,
        "scaled_mean_definition": "a=2p E[A]",
        "lift_scaled_mass_definition": "a-floor=4p E[L]",
        "pointwise_baselines": {"b=0": "0", "b=2": "(x_i-x_j)^2", "b=22": "x_j"},
        "baseline_subtraction_is_nonnegative_by_parity": True,
        "sharp_positive_lift_mass_floor": sharp_mass,
        "local_p_plus_nine_dependency": local,
        "rows": rows,
        "excluded_scaled_means": [8, 32],
        "all_even_boundaries_excluded": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_u11_all_one_common_row_exclusion() -> dict[str, object]:
    """Exclude the no-quotient-zero, all-mean-46 hard branch."""
    catalog = mean_46_hard_family_catalog()
    local = p23_phase_zero_mass_exclusions()
    rows = []
    for offset in catalog["possible_coefficient_offsets"]:
        offset = int(offset)
        candidates = [
            parallel
            for parallel in range(H_EDGE_COUNT // M + 1)
            if (parallel - offset) % Q == 0
        ]
        _require(candidates == [offset], "a mean-46 offset lost rigidity")
        hard_parallel = offset
        h_times_T = (P + 1) * hard_parallel - 3 * P - 2 * P
        hard_edges = M * hard_parallel
        opposite_edges = H_EDGE_COUNT - hard_edges

        def opposite_mean(parallel: int) -> int:
            return (P + 1) * parallel + h_times_T - 3 * P

        forbidden_Q = 8 - hard_parallel
        forced_Q = 9 - hard_parallel
        surplus = opposite_edges - M * forced_Q
        rows.append(
            {
                "coefficient_offset": offset,
                "hard_parallel_candidates": candidates,
                "hard_parallel_count": hard_parallel,
                "hard_sign_times_global_T": h_times_T,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "forbidden_Q": forbidden_Q,
                "forbidden_scaled_mass": opposite_mean(forbidden_Q),
                "forced_Q": forced_Q,
                "forced_scaled_mass": opposite_mean(forced_Q),
                "opposite_surplus": surplus,
                "directions_at_forced_Q_at_least": M - surplus,
                "excluded_by_Proposition_15_752": bool(local["proved"]),
            }
        )
    proved = bool(
        catalog["proved"]
        and local["proved"]
        and [int(row["coefficient_offset"]) for row in rows] == [4, 5, 6, 7, 8]
        and all(int(row["forbidden_scaled_mass"]) == 8 for row in rows)
        and all(int(row["forced_scaled_mass"]) == P + 9 == 32 for row in rows)
        and all(int(row["opposite_surplus"]) == 7 for row in rows)
        and all(int(row["directions_at_forced_Q_at_least"]) == 5 for row in rows)
    )
    _require(proved, "a p23 mean-46 common-row branch survived")
    return {
        "residue_u": 11,
        "branch": "no quotient-zero row, hence all twelve quotients equal one",
        "hard_scaled_mean": 46,
        "hard_family_catalog": catalog,
        "common_row_identity": "hT=24P-115",
        "opposite_mean_identity": "a(Q)=24(P+Q)-184",
        "offset_ledgers": rows,
        "all_boundary_opposite_mass_bridge": local,
        "dependency": "phase-zero floors, pointwise parity baselines, and15.688/15.752",
        "excluded": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_u11_zero_quotient_exclusion() -> dict[str, object]:
    """Exclude ``u=11`` when one old exact quotient-zero row occurs."""
    rules = baseline_coefficient_rules(P)
    local = p23_phase_zero_mass_exclusions()
    baseline_rows = (
        (2, BRANCH_B2, int(rules[BRANCH_B2]["offset"])),
        (P - 1, BRANCH_P3_LAST, int(rules[BRANCH_P3_LAST]["offset"])),
    )
    rows = []
    for boundary, branch, offset in baseline_rows:
        # For any hard row, a=22+24k.  Common hT therefore makes P-k one
        # common integer c.  Here k=0, and sum P=12c+sum k=12c+12.
        base_candidates = [
            base
            for base in range((H_EDGE_COUNT - M) // M + 1)
            if (base - offset) % Q == 0
        ]
        _require(base_candidates == [offset], "a quotient-zero offset changed")
        base = offset
        h_times_T = (P + 1) * base - 3 * P - (P - 1)
        hard_edges = M * base + M
        opposite_edges = H_EDGE_COUNT - hard_edges

        def opposite_mean(parallel: int) -> int:
            return (P + 1) * parallel + h_times_T - 3 * P

        forbidden_Q = 7 - base
        forced_Q = 8 - base
        surplus = opposite_edges - M * forced_Q
        rows.append(
            {
                "quotient_zero_boundary_b": boundary,
                "baseline_branch": branch,
                "baseline_coefficient_offset": offset,
                "common_base_candidates_for_P_minus_k": base_candidates,
                "common_base_c": base,
                "hard_parallel_formula": "P_L=c+k_L",
                "hard_sign_times_global_T": h_times_T,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "opposite_mean_identity": f"a(Q)=24*({base}+Q)-160",
                "forbidden_Q": forbidden_Q,
                "forbidden_scaled_mass": opposite_mean(forbidden_Q),
                "forced_Q": forced_Q,
                "forced_scaled_mass": opposite_mean(forced_Q),
                "opposite_surplus": surplus,
                "directions_at_forced_Q_at_least": M - surplus,
                "excluded_by_Proposition_15_752": bool(local["proved"]),
            }
        )
    proved = bool(
        rules["proved"]
        and local["proved"]
        and [(row["quotient_zero_boundary_b"], row["common_base_c"]) for row in rows]
        == [(2, 4), (22, 3)]
        and all(int(row["forbidden_scaled_mass"]) == 8 for row in rows)
        and all(int(row["forced_scaled_mass"]) == P + 9 == 32 for row in rows)
        and all(int(row["opposite_surplus"]) == 7 for row in rows)
        and all(int(row["directions_at_forced_Q_at_least"]) == 5 for row in rows)
    )
    _require(proved, "a p23 u=11 quotient-zero branch survived")
    return {
        "residue_u": 11,
        "branch": "at least one quotient-zero exact b=2/22 row",
        "common_row_identity": "hT=24c-91 and P_L=c+k_L",
        "opposite_mean_identity": "a(Q)=24(c+Q)-160",
        "baseline_ledgers": rows,
        "all_boundary_opposite_mass_bridge": local,
        "dependency": "phase-zero floors, pointwise parity baselines, and15.688/15.752",
        "excluded": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_u11_common_row_exclusion() -> dict[str, object]:
    """Package the exhaustive quotient-zero/all-one ``u=11`` split."""
    zero = p23_u11_zero_quotient_exclusion()
    all_one = p23_u11_all_one_common_row_exclusion()
    proved = bool(zero["proved"] and all_one["proved"])
    _require(proved, "the exhaustive p23 u=11 split failed")
    return {
        "residue_u": 11,
        "exhaustive_quotient_profile_dichotomy": (
            "sum k_L=12: either some k_L=0, or every one of the twelve "
            "nonnegative integer quotients equals one"
        ),
        "quotient_zero_present": zero,
        "no_quotient_zero": all_one,
        "excluded": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15771() -> dict[str, object]:
    """Package the reviewed all-boundary ``p=23,t=11,k=114`` theorem."""
    residues = p23_third_post_band_residue_ledger()
    u9 = p23_u9_two_unit_carry_exclusion()
    u10 = p23_p_minus_one_local_exclusion()
    u11 = p23_u11_common_row_exclusion()
    proved = bool(
        residues["proved"]
        and u9["proved"]
        and u10["proved"]
        and u11["proved"]
        and residues["arithmetic_surviving_residues"] == [9, 10, 11]
    )
    _require(proved, "Proposition 15.771 failed")
    return {
        "prop": "15.771",
        "status": "PROVED all-boundary exceptional endpoint theorem",
        "certificate_checks_passed": proved,
        "proof_review_complete": True,
        "pending_proof_bridges": [],
        "completed_proof_bridges": [
            "general-slice equality reduction for b=4 and b=20",
            "covering swap cubes for the middle-boundary rank certificates",
            "phase-zero mass-32 floor and p+9 exclusion bridge",
        ],
        "proof_notes": [
            "evidence/NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md",
            "evidence/NOTE_2026-09-04_P23_SMALL_BOUNDARY_EQUALITY_PROOF.md",
            "evidence/NOTE_2026-09-04_P23_MIDDLE_BOUNDARY_CUBE_PROOF.md",
        ],
        "independent_mesh_replay": "evidence/p23_third_post_band_mesh_replay.json",
        "statement": (
            "the residual-(ii) isolated-chart branch at p=23,t=11,k=114 "
            "is empty for every boundary size"
        ),
        "residue_ledger": residues,
        "branch_exclusions": {
            "u=9": u9,
            "u=10": {
                "forced_low_direction_count_at_least": 11,
                "forced_low_scaled_mean": 44,
                "difference_scaled_mass": 22,
                "local_exclusion": u10,
                "proved": u10["proved"],
            },
            "u=11": u11,
        },
        "p23_k114_closed": proved,
        "all_boundary_sizes_excluded": proved,
        "new_graph_or_residual_configuration_census_used": False,
        "fixed_four_bit_boolean_catalog_reused": True,
        "fixed_p23_five_set_coefficient_certificate_reused_for_u9": True,
        "later_p23_layers_closed": False,
        "residual_ii_closed_globally": False,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_prop15771.json"
    write_json_atomic(path, proposition_15771())
    return path


def main() -> None:
    path = write_evidence()
    row = proposition_15771()
    print(json.dumps({"proved": row["proved"], "status": row["status"],
                      "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
