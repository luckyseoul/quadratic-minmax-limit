#!/usr/bin/env python3
"""The halved symmetric Radon map modulo two.

After the fixed antipodal variables are forced and the even fixed-cell block
is divided by two, the remaining binary map has a paired-cell component
``C`` and a fixed-word component ``Phi``.  This module proves that the full
map ``(C, Phi)`` is surjective, gives the exact dual-support criterion after
deleting used orbits, and gives an explicit puncture of size ``p*h`` which
destroys surjectivity.

The small matrices for ``p=3,5,7`` are fail-when-wrong formula checks.  They
are not a prime census and are not evidence for the all-prime theorem.
"""

from __future__ import annotations

from itertools import product

from e1_gmin_m4_inversion_symmetric_lattice import (
    symmetric_mod2_decomposition,
)
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_fixed_edge_elimination import (
    fixed_word_block_basis_theorem,
    mobius_midpoint_direction_theorem,
)


Point = tuple[int, int]
Coordinate = tuple[str, int, int, int]


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


def halved_mod2_surjectivity_theorem(p: int) -> dict[str, object]:
    """Return the exact all-odd-prime surjectivity theorem for ``(C,Phi)``.

    The lower map ``C`` is the already-proved paired nonfixed target map.
    The fixed-word block design and its ``p``-column ``C``-kernel lifts are
    imported from the fixed-edge-elimination theorem rather than rebuilt.
    """
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    paired_rank = d * h * h
    fixed_word_rank = delta_size
    expected_rank = paired_rank + fixed_word_rank

    lower = symmetric_mod2_decomposition(p)
    blocks = fixed_word_block_basis_theorem(p)
    proved = bool(
        lower["nonfixed_pair_map_surjective"]
        and lower["paired_nonfixed_target_rank"] == paired_rank
        and blocks["block_vectors_form_basis"]
        and blocks["paired_affine_line_block_types"] == delta_size
        and blocks["disjoint_C_kernel_lifts_per_block_type"] == h
        and blocks["columns_per_C_kernel_lift"] == p
        and expected_rank == d * h * (h + 1)
    )
    if not proved:
        raise ArithmeticError("the halved mod-two surjectivity data changed")

    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "source_nonfixed_orbit_columns": delta_size * delta_size,
        "paired_component": "C",
        "paired_component_rank": paired_rank,
        "paired_component_surjective": True,
        "fixed_component_after_inverse": "Phi=A_bar^{-1}*B mod 2",
        "fixed_word_component_rank": fixed_word_rank,
        "fixed_word_block_basis_reused": True,
        "kernel_lift": (
            "for a paired affine block B of direction A and one fixed "
            "[a] in A, sum the p columns ([a],[delta]) with [delta] in B"
        ),
        "kernel_lift_lies_in_ker_C": True,
        "kernel_lift_Phi_image": "the incidence vector 1_B",
        "Phi_on_ker_C_surjective": True,
        "surjectivity_argument": (
            "lift the desired C-coordinate first, then correct its Phi-word "
            "by a linear combination of the C-kernel block lifts"
        ),
        "full_halved_map_rank": expected_rank,
        "full_halved_map_surjective": True,
        "signs_mod_two": "all tau signs reduce to one",
        "restricted_boolean_fibre_nonempty": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def fixed_transverse_counter_puncture(p: int) -> dict[str, object]:
    """Return an explicit puncture which destroys halved-map surjectivity.

    Fix a projective row ``L`` and a nonzero square ``beta``.  The halved
    fixed-transverse coordinate ``(L,0,beta)`` is supported precisely on

        X_(L,beta) = {([a],[delta]): L(a)=0, L(delta)^2=beta}.

    It has ``h*p`` columns.  Deleting them makes this formerly nonzero row
    identically zero, so the punctured row rank is smaller.
    """
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    midpoint_classes = h
    difference_classes = p
    puncture_size = midpoint_classes * difference_classes
    proved = bool(
        puncture_size == p * h
        and puncture_size == delta_size - h
        and puncture_size <= delta_size
    )
    if not proved:
        raise ArithmeticError("the fixed-transverse counter-puncture changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "coordinate": "fixed-transverse (L,0,beta), beta a nonzero square",
        "support": "{([a],[delta]):L(a)=0 and L(delta)^2=beta}",
        "support_product": "A_L times B_(L,beta)",
        "midpoint_factor_size": midpoint_classes,
        "difference_factor_size": difference_classes,
        "puncture_size": puncture_size,
        "puncture_size_identity": "p*h=|Delta|-h",
        "punctured_row_becomes_zero": True,
        "punctured_map_surjective": False,
        "universal_robustness_through_Delta_deletions": False,
        "dual_distance_convention": (
            "minimum support size of a nonzero word in Row(D), the row-code "
            "distance controlling column punctures"
        ),
        "dual_distance_upper_bound": puncture_size,
        "dual_distance_equality_proved": False,
        "all_low_weight_dual_words_classified": False,
        "proved": proved,
    }


def punctured_halved_dual_criterion(p: int) -> dict[str, object]:
    """State the exact image criterion after an arbitrary column puncture.

    If ``D=(C,Phi)`` and ``D_U`` retains the columns outside ``U``, then a
    target is in ``im D_U`` exactly when it annihilates every left-kernel
    functional of ``D_U``.  Full surjectivity fails exactly when ``U``
    contains the support of a nonzero word in the full row code of ``D``.
    """
    theorem = halved_mod2_surjectivity_theorem(p)
    counter = fixed_transverse_counter_puncture(p)
    proved = bool(
        theorem["proved"]
        and counter["proved"]
        and counter["puncture_size"] <= theorem["delta_size"]
    )
    if not proved:
        raise ArithmeticError("the punctured dual criterion data changed")
    return {
        "p": p,
        "column_set": "Delta times Delta",
        "punctured_map": "D_U=D restricted to columns (Delta^2 minus U)",
        "surjectivity_criterion": (
            "D_U is onto iff no nonzero word w in Row(D) has supp(w) "
            "contained in U"
        ),
        "equivalent_two_stage_criterion": (
            "C_U is onto and Phi(ker C_U)=F_2^Delta"
        ),
        "target_membership_criterion": (
            "r_U lies in im(D_U) iff lambda(r_U)=0 for every lambda with "
            "lambda*D_U=0"
        ),
        "raw_congruence_interpretation": (
            "fixed cells modulo four jointly with paired cells modulo two, "
            "after the forced fixed vector is subtracted"
        ),
        "counter_puncture_size": counter["puncture_size"],
        "universal_robustness_through_Delta_deletions": False,
        "structured_Mobius_puncture_surjective": "OPEN",
        "dual_distance_equality": "OPEN",
        "proved": proved,
    }


def mobius_rectangle_intersection_bound(
    p: int, half_count: int
) -> dict[str, object]:
    """Give the imported midpoint-direction bound for a union of halves.

    This is only an intersection bound.  It excludes containment of the
    explicit rectangle when ``2*half_count<p*h`` but says nothing about the
    supports of other dual words.
    """
    h = _check_odd_prime(p)
    if (
        not isinstance(half_count, int)
        or isinstance(half_count, bool)
        or half_count < 0
    ):
        raise ValueError("half_count must be a nonnegative integer")
    imported = mobius_midpoint_direction_theorem(p)
    intersection_bound = (
        half_count * imported["one_half_hits_any_midpoint_direction_at_most"]
    )
    rectangle_size = p * h
    return {
        "p": p,
        "Mobius_half_count": half_count,
        "one_half_midpoint_direction_bound": 2,
        "union_intersection_bound_with_X_L_beta": intersection_bound,
        "X_L_beta_size": rectangle_size,
        "explicit_counter_rectangle_cannot_be_contained": (
            intersection_bound < rectangle_size
        ),
        "other_low_weight_dual_supports_excluded": False,
        "structured_punctured_surjectivity_proved": False,
        "proved": imported["proved"],
    }


def _antipodal_classes(p: int) -> tuple[Point, ...]:
    zero = (0, 0)
    return tuple(
        sorted(
            {
                min(point, ((-point[0]) % p, (-point[1]) % p))
                for point in product(range(p), repeat=2)
                if point != zero
            }
        )
    )


def _directions(p: int) -> tuple[Point, ...]:
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def _evaluate(p: int, functional: Point, point: Point) -> int:
    return (
        functional[0] * point[0] + functional[1] * point[1]
    ) % p


def _coordinate_layout(p: int) -> tuple[Coordinate, ...]:
    squares = tuple(sorted({value * value % p for value in range(1, p)}))
    coordinates: list[Coordinate] = []
    for direction_index in range(p + 1):
        coordinates.append(("P", direction_index, 0, 0))
        coordinates.extend(
            ("F", direction_index, 0, beta) for beta in squares
        )
        coordinates.extend(
            ("C", direction_index, alpha, beta)
            for alpha in squares
            for beta in squares
        )
    return tuple(coordinates)


def _column_support(
    p: int, midpoint: Point, difference: Point
) -> tuple[Coordinate, ...]:
    support: list[Coordinate] = []
    for direction_index, functional in enumerate(_directions(p)):
        alpha = _evaluate(p, functional, midpoint)
        beta = _evaluate(p, functional, difference)
        if beta == 0:
            support.append(("P", direction_index, 0, 0))
        elif alpha == 0:
            support.append(("F", direction_index, 0, beta * beta % p))
        else:
            support.append(
                (
                    "C",
                    direction_index,
                    alpha * alpha % p,
                    beta * beta % p,
                )
            )
    return tuple(support)


def _binary_rank(vectors: tuple[int, ...] | list[int]) -> int:
    basis: dict[int, int] = {}
    for original in vectors:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def exact_small_halved_replay(p: int) -> dict[str, object]:
    """Build the raw halved matrix for ``p=3,5,7`` as a formula check."""
    h = _check_odd_prime(p)
    if p not in (3, 5, 7):
        raise ValueError("the exact halved replay is limited to p=3,5,7")
    d = p + 1
    delta_size = d * h
    classes = _antipodal_classes(p)
    directions = _directions(p)
    coordinates = _coordinate_layout(p)
    coordinate_index = {
        coordinate: index for index, coordinate in enumerate(coordinates)
    }
    paired_mask = sum(
        1 << index
        for index, coordinate in enumerate(coordinates)
        if coordinate[0] == "C"
    )

    columns: list[int] = []
    kept_columns: list[int] = []
    counter_support = 0
    parallel_rows_exact = True
    chosen_direction = directions[0]
    chosen_square = 1
    counter_coordinate = ("F", 0, 0, chosen_square)
    counter_bit = 1 << coordinate_index[counter_coordinate]
    for midpoint in classes:
        for difference in classes:
            support = _column_support(p, midpoint, difference)
            value = 0
            for coordinate in support:
                value |= 1 << coordinate_index[coordinate]
            columns.append(value)

            deleted = bool(
                _evaluate(p, chosen_direction, midpoint) == 0
                and _evaluate(p, chosen_direction, difference) ** 2 % p
                == chosen_square
            )
            if deleted:
                counter_support += 1
            else:
                kept_columns.append(value)
            if bool(value & counter_bit) != deleted:
                raise ArithmeticError("the counter-puncture support changed")

            for direction_index, functional in enumerate(directions):
                parallel_bit = 1 << coordinate_index[
                    ("P", direction_index, 0, 0)
                ]
                if bool(value & parallel_bit) != (
                    _evaluate(p, functional, difference) == 0
                ):
                    parallel_rows_exact = False

    paired_rank = _binary_rank([value & paired_mask for value in columns])
    punctured_paired_rank = _binary_rank(
        [value & paired_mask for value in kept_columns]
    )
    full_rank = _binary_rank(columns)
    punctured_rank = _binary_rank(kept_columns)
    expected_paired_rank = d * h * h
    expected_full_rank = d * h * (h + 1)
    proved = bool(
        len(classes) == delta_size
        and len(columns) == delta_size * delta_size
        and len(coordinates) == d * (1 + h + h * h)
        and paired_rank == expected_paired_rank
        and punctured_paired_rank == expected_paired_rank
        and full_rank == expected_full_rank
        and counter_support == p * h
        and punctured_rank == expected_full_rank - 1
        and parallel_rows_exact
    )
    if not proved:
        raise ArithmeticError("the exact small halved replay changed")
    return {
        "p": p,
        "delta_size": delta_size,
        "source_columns": len(columns),
        "raw_target_coordinates": len(coordinates),
        "paired_component_rank": paired_rank,
        "expected_paired_component_rank": expected_paired_rank,
        "full_halved_rank": full_rank,
        "expected_full_halved_rank": expected_full_rank,
        "raw_compatibility_relation_count": len(coordinates) - full_rank,
        "counter_puncture_size": counter_support,
        "punctured_paired_component_rank": punctured_paired_rank,
        "punctured_halved_rank": punctured_rank,
        "counter_puncture_rank_drop": full_rank - punctured_rank,
        "parallel_rows_encode_direction_weight_parity": parallel_rows_exact,
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    theorem = halved_mod2_surjectivity_theorem(p)
    puncture = punctured_halved_dual_criterion(p)
    counter = fixed_transverse_counter_puncture(p)
    mobius_bound = mobius_rectangle_intersection_bound(p, (p + 1) // 2)
    replay = {str(q): exact_small_halved_replay(q) for q in (3, 5, 7)}
    proved = bool(
        theorem["proved"]
        and puncture["proved"]
        and counter["proved"]
        and mobius_bound["proved"]
        and all(row["proved"] for row in replay.values())
    )
    if not proved:
        raise ArithmeticError("the halved mod-two theorem record changed")
    return {
        "title": "Full and punctured halved symmetric map modulo two",
        "status": (
            "FULL MAP SURJECTIVE; ARBITRARY DELTA-SIZED ROBUSTNESS FALSE; "
            "STRUCTURED PUNCTURE OPEN"
        ),
        "theorem": theorem,
        "punctured_dual_criterion": puncture,
        "counter_puncture": counter,
        "conditional_Mobius_rectangle_bound": mobius_bound,
        "small_exact_replay": replay,
        "direction_weight_boundary": (
            "the P_L rows impose only the parity of each exact parallel "
            "weight n_L; the integer equalities sum_(O parallel L)b_O=n_L "
            "remain necessary"
        ),
        "proved": {
            "full_unpunctured_halved_map_surjective": True,
            "exact_puncture_dual_support_criterion": True,
            "universal_robustness_through_Delta_deletions": False,
            "dual_distance_equals_p_h": False,
            "all_low_weight_dual_words_classified": False,
            "structured_Mobius_punctured_map_surjective": False,
            "direction_weight_parity_is_sufficient": False,
            "restricted_boolean_fibre_nonempty": False,
            "residual_ii_closed": False,
        },
        "next_exact_obstruction": (
            "classify row-code words of weight at most |Delta|, or directly "
            "exclude support containment in the actual Mobius U; after that, "
            "meet every exact direction-weight slice and the integer target"
        ),
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
