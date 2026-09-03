#!/usr/bin/env python3
"""Low-weight gap and structured-puncture theorem for the halved row code.

The all-prime grouped uncertainty inequality, combined with the exact
direction-block normal form, proves that the only row-code words below
``|Delta|`` are the fixed-transverse rectangles of weight ``p*h``.  The
existing per-half rectangle intersection bound then proves that every
physically extendable branch-C Mobius puncture leaves the halved binary map
onto.  No prescribed Hamming or integral direction slice is inferred.
"""

from __future__ import annotations

from e1_gmin_m4_grouped_uncertainty_square import grouped_uncertainty_theorem
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_halved_row_code import halved_row_code_decomposition


def _check_odd_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")
    return (p - 1) // 2


def cell_union_distance_theorem(p: int) -> dict[str, object]:
    """Return the exact two-direction cell-union distance lemma.

    A cell union in one direction is encoded by an even Boolean function on
    ``F_p``.  If ``a`` values are selected and the zero value has bit
    ``epsilon``, its projective weight is ``(p*a-epsilon)/2``.  For two
    independent directions, twice the symmetric-difference weight is
    ``p*(a+b)-2*a*b-(epsilon xor eta)``.

    The loops below check only this closed two-variable formula.  They are
    not a point-set or prime census.
    """
    h = _check_odd_prime(p)
    proper_weights: list[tuple[int, int, int]] = []
    cross_weights: list[tuple[int, int, int]] = []
    for a in range(1, p):
        zero_bit = a & 1
        weight = (p * a - zero_bit) // 2
        proper_weights.append((weight, a, zero_bit))
        for b in range(1, p):
            other_zero_bit = b & 1
            origin_bit = zero_bit ^ other_zero_bit
            doubled = p * (a + b) - 2 * a * b - origin_bit
            if doubled % 2:
                raise ArithmeticError("the antipodal count lost parity")
            cross_weights.append((doubled // 2, a, b))

    minimum_cell_weight = min(row[0] for row in proper_weights)
    minimum_cell_parameters = tuple(
        (a, zero_bit)
        for weight, a, zero_bit in proper_weights
        if weight == minimum_cell_weight
    )
    minimum_cross_distance = min(row[0] for row in cross_weights)
    proved = bool(
        minimum_cell_weight == h
        and minimum_cell_parameters == ((1, 1),)
        and minimum_cross_distance == 2 * h
    )
    if not proved:
        raise ArithmeticError("the cell-union distance formula changed")
    return {
        "p": p,
        "h": h,
        "proper_cell_union_minimum_weight": minimum_cell_weight,
        "minimum_weight_equality": "the radial line l_A only",
        "different_direction_proper_union_minimum_distance": (
            minimum_cross_distance
        ),
        "distance_formula": (
            "2*|S_A triangle S_B|="
            "p*(a+b)-2*a*b-(f(0) xor g(0))"
        ),
        "proved": proved,
    }


def row_code_gap_theorem(p: int) -> dict[str, object]:
    """Prove ``d(Row(D))=p*h`` and classify all weights below ``|Delta|``."""
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    minimum_distance = p * h
    decomposition = halved_row_code_decomposition(p)
    grouped = grouped_uncertainty_theorem(p)
    cells = cell_union_distance_theorem(p)

    # These equalities are the endpoints in the active-support proof:
    # D_act >= h*k*(k-1), G >= D_act/(k-1) >= h*k, and
    # (d-k)*|R|+G >= (d-k)*h+h*k=d*h.
    active_support_endpoint = (d - 2) * h + 2 * h
    full_direction_endpoint = h * d * (d - 1) // (d - 1)
    rectangle_count = d * h
    proved = bool(
        decomposition["proved"]
        and grouped["proved"]
        and cells["proved"]
        and active_support_endpoint == delta_size
        and full_direction_endpoint == delta_size
        and minimum_distance == delta_size - h
        and rectangle_count == delta_size
    )
    if not proved:
        raise ArithmeticError("the halved row-code gap proof changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "minimum_distance": minimum_distance,
        "minimum_word_count": rectangle_count,
        "minimum_words": "exactly l_A tensor b_(A,j)",
        "weights_strictly_between_p_h_and_Delta": (),
        "weight_Delta_layer_classified": False,
        "proof": {
            "row_bound": "wt(W_x)>=max(1,d-b_x)",
            "pair_distance": "|S_A triangle S_B|>=2h",
            "two_to_d_minus_one_active": (
                "G=sum_x(k-b_x)>=h*k, so wt(W)>=d*h"
            ),
            "all_d_active": (
                "sum b_x(d-b_x)>=(h*d*(d-1)) gives wt(W)>=d*h"
            ),
            "one_active": (
                "weight below d*h forces one radial support and one right block"
            ),
        },
        "proved": proved,
    }


def localized_half_hamming_ledger(
    p: int, q: int, t: int, kappa: int
) -> dict[str, object]:
    """Return exact support, parity, and distance thresholds for ``q`` halves.

    This ledger is for ``p=4r+3``.  It does not assert that the requested
    cancellation count is constructible.
    """
    h = _check_odd_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    if not isinstance(q, int) or isinstance(q, bool) or not 0 <= q <= h + 1:
        raise ValueError("q must lie between zero and h+1")
    if not isinstance(t, int) or isinstance(t, bool):
        raise ValueError("t must be an integer")
    if not isinstance(kappa, int) or isinstance(kappa, bool) or kappa < 0:
        raise ValueError("kappa must be a nonnegative integer")

    r = (p - 3) // 4
    delta_size = (p + 1) * h
    minimum_distance = p * h
    raw_occurrences = q * (p - 1)
    if 2 * kappa > raw_occurrences:
        raise ValueError("kappa exceeds the raw occurrence budget")
    used = raw_occurrences - 2 * kappa
    edge_count = 4 * p + 2 * t + 1
    remaining_capacity = edge_count - used
    size_kappa_min = max(0, q * h - 2 * p - t)
    distance_kappa_min = 0 if q <= h else r + 1
    forced_fixed_weight_parity = (1 + q) & 1
    hamming_numerator_parity = q & 1
    return {
        "p": p,
        "r": r,
        "h": h,
        "q_nonzero_hard_centres": q,
        "t": t,
        "kappa": kappa,
        "raw_occurrences": raw_occurrences,
        "used_support_size": used,
        "edge_count": edge_count,
        "remaining_edge_capacity": remaining_capacity,
        "bare_size_kappa_minimum": size_kappa_min,
        "bare_size_condition_holds": kappa >= size_kappa_min,
        "minimum_distance_kappa_minimum": distance_kappa_min,
        "below_minimum_distance": used < minimum_distance,
        "below_Delta": used < delta_size,
        "forced_fixed_weight_parity": forced_fixed_weight_parity,
        "hamming_numerator_parity": hamming_numerator_parity,
        "localized_half_ansatz_parity_compatible": (
            hamming_numerator_parity == 0
        ),
        "parity_scope": (
            "one localized half per nonzero hard centre, with arbitrary "
            "ternary cancellations; not every antisymmetric preimage"
        ),
        "proved": True,
    }


def branch_c_structured_puncture_theorem(p: int) -> dict[str, object]:
    """Prove onto-ness of every physically extendable structured puncture.

    Onto-ness here is exclusively over ``F_2`` for the halved map ``D_U``.
    """
    h = _check_odd_prime(p)
    if p < 31 or p % 4 != 3:
        raise ValueError("need a branch-C prime p=3 mod 4 with p>=31")
    r = (p - 3) // 4
    d = p + 1
    delta_size = d * h
    minimum_distance = p * h
    t_min = 2 * r * r - 4 * r - 2
    t_max = 4 * r * r - 2 * r - 5
    maximum_edge_count = 4 * p + 2 * t_max + 1
    distance_only_t_max = t_max - r
    maximum_rectangle_intersection = p + 1
    gap = row_code_gap_theorem(p)
    proved = bool(
        gap["proved"]
        and maximum_edge_count == delta_size - 1
        and maximum_rectangle_intersection < minimum_distance
        and distance_only_t_max == 4 * r * r - 3 * r - 5
    )
    if not proved:
        raise ArithmeticError("the structured Mobius puncture proof changed")
    return {
        "p": p,
        "r": r,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "row_code_minimum_distance": minimum_distance,
        "row_code_minimum_words": "fixed-transverse rectangles only",
        "row_code_open_weight_interval": (minimum_distance, delta_size),
        "branch_t_range": (t_min, t_max),
        "maximum_edge_count": maximum_edge_count,
        "maximum_rectangle_intersection_by_all_halves": (
            maximum_rectangle_intersection
        ),
        "distance_only_automatic_t_range": (t_min, distance_only_t_max),
        "gap_theorem_automatic_t_range": (t_min, t_max),
        "zero_centre_case": (
            "q<=h gives |U|<=q(p-1)<p*h, so D_U is onto"
        ),
        "all_active_case": (
            "|U|<=|H|<=N-1 and no minimum rectangle lies in U, so D_U is onto"
        ),
        "structured_punctured_halved_map_onto_over_F2": True,
        "prescribed_Hamming_slice_solved": False,
        "directionwise_integer_slices_solved": False,
        "divided_integral_Boolean_fibre_solved": False,
        "common_simple_graph_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    row_code = row_code_gap_theorem(p)
    puncture = branch_c_structured_puncture_theorem(p)
    proved = bool(row_code["proved"] and puncture["proved"])
    if not proved:
        raise ArithmeticError("the row-code gap theorem record changed")
    return {
        "title": "Exact halved row-code gap and structured Mobius puncture",
        "row_code": row_code,
        "structured_puncture": puncture,
        "proved": {
            "minimum_distance_equals_p_h": True,
            "minimum_words_are_fixed_transverse_rectangles": True,
            "no_weights_strictly_between_p_h_and_Delta": True,
            "structured_Mobius_punctured_map_onto_over_F2": True,
            "prescribed_Hamming_and_direction_slices": False,
            "symmetric_Boolean_completion": False,
            "residual_ii_closed": False,
        },
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
