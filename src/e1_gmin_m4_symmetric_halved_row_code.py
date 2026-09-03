#!/usr/bin/env python3
"""Structure and low-weight boundary of the halved symmetric row code.

The paired-affine-line incidence matrix splits the point space into one
orthogonal block space per projective direction.  In those terms the full
halved row code is a boundary summand plus the algebra of direction-block
diagonal matrices.  This gives an exact normal form and exhibits weight
``|Delta|`` words which are not fixed-transverse rectangles.

The result does not determine the minimum distance.  In particular, the
lower bound ``p*h`` and the classification of its equality words remain
open here.
"""

from __future__ import annotations

from itertools import combinations

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_fixed_edge_elimination import (
    fixed_word_block_basis_theorem,
)
from e1_gmin_m4_symmetric_halved_mod2 import (
    _antipodal_classes,
    _binary_rank,
    _column_support,
    _coordinate_layout,
    _directions,
    _evaluate,
)


Point = tuple[int, int]


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


def halved_row_code_decomposition(p: int) -> dict[str, object]:
    """Return the exact boundary-plus-block-diagonal row-code theorem."""
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    design = fixed_word_block_basis_theorem(p)
    boundary_rank = delta_size
    diagonal_rank = d * h * h
    total_rank = boundary_rank + diagonal_rank
    proved = bool(
        design["proved"]
        and design["block_vectors_form_basis"]
        and design["paired_affine_line_block_types"] == delta_size
        and total_rank == d * h * (h + 1)
    )
    if not proved:
        raise ArithmeticError("the halved row-code decomposition changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "point_space": "H=F_2^Delta",
        "direction_block_space": (
            "B_A=span of the h disjoint p-point nonorigin affine blocks "
            "of direction A"
        ),
        "orthogonal_point_decomposition": "H=direct_sum_A B_A",
        "row_code_decomposition": (
            "Row(D)=(<1_Delta> tensor H) direct_sum "
            "direct_sum_A(B_A tensor B_A)"
        ),
        "boundary_rank": boundary_rank,
        "direction_block_diagonal_rank": diagonal_rank,
        "total_rank": total_rank,
        "generator_conversion": {
            "paired_row": "1_I tensor 1_J for I,J in the same B_A",
            "fixed_transverse_row": (
                "1_(line A) tensor 1_J = 1 tensor 1_J + "
                "sum_(I in B_A) 1_I tensor 1_J"
            ),
            "parallel_row": "1 tensor 1_(line A)",
            "boundary_block": (
                "1 tensor 1_J is fixed-transverse plus all paired rows "
                "with right block J"
            ),
        },
        "matrix_normal_form": (
            "after conjugating by the orthogonal block incidence matrix M, "
            "a word is 1*q^T plus a matrix block-diagonal by direction"
        ),
        "minimum_distance_p_h_proved": False,
        "minimum_weight_words_classified": False,
        "proved": proved,
    }


def block_incidence_branch_number_theorem(p: int) -> dict[str, object]:
    """Prove ``wt(x)+wt(M^T x)>=p+1`` for every nonzero binary word.

    If ``S=supp(x)``, put ``s=|S|`` and let ``n_B=|S intersect B|``.
    A point lies on ``p`` blocks; two noncollinear point classes have two
    common blocks, while two distinct collinear classes have none.  If ``c``
    counts collinear pairs, then

        sum_B n_B = p*s,
        sum_B n_B^2 = p*s + 4*(binom(s,2)-c).

    Since ``1_(n odd)>=2*n-n^2``, for ``t=wt(M^T x)`` one gets
    ``t>=s*(p-2*s+2)+4*c``.  This implies ``s+t>=p+1`` when
    ``s<=h``.  If ``s>h``, either ``t>h`` or the same argument applies to
    ``M^T x`` and uses ``M*M^T=I``.
    """
    h = _check_odd_prime(p)
    design = fixed_word_block_basis_theorem(p)
    lower_branch = p + 1
    singleton_weight = 1
    singleton_transform_weight = p
    block_weight = p
    block_transform_weight = 1
    proved = bool(
        design["proved"]
        and design["binary_gram_identity"] == "M*M^T=I, hence M^T*M=I"
        and singleton_weight + singleton_transform_weight == lower_branch
        and block_weight + block_transform_weight == lower_branch
    )
    if not proved:
        raise ArithmeticError("the affine-block branch number changed")
    return {
        "p": p,
        "h": h,
        "transform": "x maps to M^T*x",
        "branch_inequality": "wt(x)+wt(M^T*x)>=p+1 for x nonzero",
        "branch_number": lower_branch,
        "second_moment_identity": (
            "sum_B n_B^2=p*s+4*(binom(s,2)-c_collinear)"
        ),
        "odd_block_lower_bound": (
            "wt(M^T*x)>=s*(p-2*s+2)+4*c_collinear"
        ),
        "small_side_range": "1<=s<=h",
        "large_side_reason": (
            "if s>=h+1 and t>=h+1 the sum bound is immediate; otherwise "
            "apply the small-side bound to M^T*x and invert with M*M^T=I"
        ),
        "equality_witnesses": (
            "one point has transform weight p; one affine block has "
            "transform weight one"
        ),
        "minimum_halved_row_code_distance_deduced": False,
        "proved": proved,
    }


def low_weight_counterfamilies(p: int) -> dict[str, object]:
    """Record exact row-code families at weights ``p*h`` and ``|Delta|``."""
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    rectangle_weight = p * h
    vertical_count = delta_size
    scalar_graph_count = h
    proved = bool(
        rectangle_weight == delta_size - h
        and vertical_count == delta_size
        and scalar_graph_count == (p - 1) // 2
    )
    if not proved:
        raise ArithmeticError("the low-weight counterfamily counts changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "fixed_transverse_rectangles": {
            "support": "line_A times affine_block_(A,beta)",
            "weight": rectangle_weight,
            "count": delta_size,
        },
        "vertical_fibres": {
            "support": "Delta times {[delta_0]}",
            "weight": delta_size,
            "count": vertical_count,
            "row_code_reason": "the boundary summand <1> tensor H",
        },
        "scalar_graphs": {
            "support": "{([a],[delta]):[a]=[c*delta]}",
            "weight": delta_size,
            "count": scalar_graph_count,
            "row_code_reason": (
                "scalar multiplication permutes the h affine blocks inside "
                "each B_A, so its permutation matrix is block-diagonal"
            ),
        },
        "all_words_of_weight_at_most_Delta_are_rectangles": False,
        "minimum_distance_equals_p_h": "OPEN",
        "equality_only_fixed_transverse_rectangles": "OPEN",
        "complete_weight_at_most_Delta_classification": "OPEN",
        "proved": proved,
    }


def _affine_blocks(p: int) -> tuple[tuple[Point, tuple[Point, ...]], ...]:
    classes = _antipodal_classes(p)
    blocks: list[tuple[Point, tuple[Point, ...]]] = []
    for functional in _directions(p):
        squares = sorted({value * value % p for value in range(1, p)})
        for square in squares:
            block = tuple(
                point
                for point in classes
                if _evaluate(p, functional, point) ** 2 % p == square
            )
            blocks.append((functional, block))
    return tuple(blocks)


def exact_small_row_code_replay(p: int) -> dict[str, object]:
    """Check the normal form and explicit families for ``p=3,5,7``."""
    h = _check_odd_prime(p)
    if p not in (3, 5, 7):
        raise ValueError("the exact row-code replay is limited to p=3,5,7")
    d = p + 1
    classes = _antipodal_classes(p)
    delta_size = len(classes)
    class_index = {point: index for index, point in enumerate(classes)}
    pairs = tuple((left, right) for left in classes for right in classes)
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    blocks = _affine_blocks(p)

    block_vectors = tuple(
        sum(1 << class_index[point] for point in block)
        for _functional, block in blocks
    )
    gram_identity = all(
        ((left & right).bit_count() & 1) == (i == j)
        for i, left in enumerate(block_vectors)
        for j, right in enumerate(block_vectors)
    )

    coordinates = _coordinate_layout(p)
    raw_rows: list[int] = []
    for coordinate in coordinates:
        value = 0
        for index, (midpoint, difference) in enumerate(pairs):
            if coordinate in _column_support(p, midpoint, difference):
                value |= 1 << index
        raw_rows.append(value)

    structured_rows: list[int] = []
    for difference in classes:
        structured_rows.append(
            sum(
                1 << pair_index[(midpoint, difference)]
                for midpoint in classes
            )
        )
    by_direction: dict[Point, list[tuple[Point, ...]]] = {}
    for functional, block in blocks:
        by_direction.setdefault(functional, []).append(block)
    for direction_blocks in by_direction.values():
        for left_block in direction_blocks:
            for right_block in direction_blocks:
                structured_rows.append(
                    sum(
                        1 << pair_index[(midpoint, difference)]
                        for midpoint in left_block
                        for difference in right_block
                    )
                )

    raw_rank = _binary_rank(raw_rows)
    structured_rank = _binary_rank(structured_rows)
    combined_rank = _binary_rank(raw_rows + structured_rows)

    vertical_words = structured_rows[:delta_size]
    vertical_weights = {word.bit_count() for word in vertical_words}
    scalar_words: list[int] = []
    for scalar in range(1, h + 1):
        scalar_words.append(
            sum(
                1
                << pair_index[
                    (
                        min(
                            (
                                scalar * difference[0] % p,
                                scalar * difference[1] % p,
                            ),
                            (
                                -scalar * difference[0] % p,
                                -scalar * difference[1] % p,
                            ),
                        ),
                        difference,
                    )
                ]
                for difference in classes
            )
        )
    scalar_rank_unchanged = (
        _binary_rank(raw_rows + scalar_words) == raw_rank
    )
    scalar_weights = {word.bit_count() for word in scalar_words}

    branch_witnesses_exact = all(
        1
        + sum(
            1
            for block_vector in block_vectors
            if block_vector & (1 << point_index)
        )
        == p + 1
        for point_index in range(delta_size)
    ) and all(
        block_vector.bit_count() + 1 == p + 1
        for block_vector in block_vectors
    )

    chosen_points = classes[: min(3, delta_size)]
    chosen_set = set(chosen_points)
    collinear_pairs = sum(
        1
        for left, right in combinations(chosen_points, 2)
        if (left[0] * right[1] - left[1] * right[0]) % p == 0
    )
    intersection_counts = tuple(
        len(chosen_set.intersection(block)) for _functional, block in blocks
    )
    second_moment_exact = bool(
        sum(value * value for value in intersection_counts)
        == p * len(chosen_points)
        + 4
        * (
            len(chosen_points) * (len(chosen_points) - 1) // 2
            - collinear_pairs
        )
    )

    expected_rank = d * h * (h + 1)
    proved = bool(
        delta_size == d * h
        and len(blocks) == delta_size
        and gram_identity
        and raw_rank == expected_rank
        and structured_rank == expected_rank
        and combined_rank == expected_rank
        and vertical_weights == {delta_size}
        and len(scalar_words) == h
        and scalar_weights == {delta_size}
        and scalar_rank_unchanged
        and branch_witnesses_exact
        and second_moment_exact
    )
    if not proved:
        raise ArithmeticError("the exact small row-code replay changed")
    return {
        "p": p,
        "delta_size": delta_size,
        "block_Gram_identity": gram_identity,
        "raw_row_code_rank": raw_rank,
        "structured_row_code_rank": structured_rank,
        "combined_rank": combined_rank,
        "vertical_fibre_count": len(vertical_words),
        "vertical_fibre_weights": sorted(vertical_weights),
        "scalar_graph_count": len(scalar_words),
        "scalar_graph_weights": sorted(scalar_weights),
        "scalar_graphs_lie_in_raw_row_code": scalar_rank_unchanged,
        "branch_number_point_and_block_witnesses": branch_witnesses_exact,
        "second_moment_formula_check": second_moment_exact,
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    decomposition = halved_row_code_decomposition(p)
    branch = block_incidence_branch_number_theorem(p)
    families = low_weight_counterfamilies(p)
    replay = {str(q): exact_small_row_code_replay(q) for q in (3, 5, 7)}
    proved = bool(
        decomposition["proved"]
        and branch["proved"]
        and families["proved"]
        and all(row["proved"] for row in replay.values())
    )
    if not proved:
        raise ArithmeticError("the halved row-code theorem record changed")
    return {
        "title": "Direction-block normal form of the halved row code",
        "status": (
            "EXACT NORMAL FORM AND WEIGHT-|Delta| COUNTERFAMILIES; "
            "MINIMUM DISTANCE OPEN"
        ),
        "decomposition": decomposition,
        "block_incidence_branch_number": branch,
        "low_weight_families": families,
        "small_exact_replay": replay,
        "proved": {
            "row_code_normal_form": True,
            "block_incidence_branch_number_p_plus_1": True,
            "nonrectangle_words_of_weight_Delta_exist": True,
            "minimum_distance_equals_p_h": False,
            "minimum_words_are_only_rectangles": False,
            "weight_at_most_Delta_classified": False,
            "structured_Mobius_punctured_surjectivity": False,
            "residual_ii_closed": False,
        },
        "next_exact_obstruction": (
            "bound the Hamming weight of 1*q^T+T when T is direction-block "
            "diagonal, or classify its words through weight |Delta|; then "
            "compare their supports with the actual Mobius puncture"
        ),
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
