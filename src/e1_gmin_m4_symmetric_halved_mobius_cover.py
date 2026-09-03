#!/usr/bin/env python3
"""Block projections forced by a Mobius-contained halved dual word.

The halved row code is considered on columns ``([a], [delta])`` in
``Delta x Delta``.  Projecting a row word through one affine difference
block gives a union of the parallel midpoint cells.  If the word is
contained in the actual support of ``(p+1)/2`` direction-localized Mobius
halves, this projection forces a nearly saturated all-halves cover.

The theorem is symbolic.  ``exact_small_cover_replay`` checks the incidence
and Mobius formulas only at ``p=3,7``; it is not theorem evidence and does
not search over choices of halves or dual words.
"""

from __future__ import annotations

from collections import Counter
from itertools import product

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    Point,
    _functional_value,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_prop15721 import is_prime


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


def _check_paley_prime(p: int, minimum: int = 3) -> int:
    h = _check_odd_prime(p)
    if p % 4 != 3 or p < minimum:
        raise ValueError(f"need a prime p=3 mod 4 with p>={minimum}")
    return h


def dual_block_projection_theorem(p: int) -> dict[str, object]:
    """Return the exact affine-block projection theorem for ``Row(D)``.

    For a nonorigin difference block ``C=B_(K,beta)``, put

        r_C([a]) = sum_([delta] in C) w([a],[delta])  (mod 2).

    If ``w`` belongs to the full halved row code, then ``r_C`` is a union
    of cells in the partition ``A_K, B_(K,alpha)`` of ``Delta``.  If every
    projection uses only ``A_K``, the word is a disjoint sum of the known
    fixed-transverse rectangles.
    """
    h = _check_odd_prime(p)
    d = p + 1
    delta_size = d * h
    block_count = delta_size
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "difference_block_count": block_count,
        "difference_block_size": p,
        "origin_midpoint_cell_size": h,
        "affine_midpoint_cell_size": p,
        "projection": (
            "r_C([a])=sum_[delta in C] w([a],[delta]) mod 2"
        ),
        "projection_partition": (
            "for C=B_(K,beta), r_C is a disjoint union of A_K and some "
            "B_(K,alpha)"
        ),
        "projection_weight": "|r_C|=epsilon_C*h+q_C*p",
        "generator_images": {
            "paired_C_row": (
                "B_(L,alpha) when C=B_(L,beta), and zero otherwise"
            ),
            "fixed_F_row": (
                "A_L when C=B_(L,beta), and zero otherwise"
            ),
            "parallel_P_row": (
                "zero when L=K and all of Delta when L!=K"
            ),
        },
        "generator_reason": (
            "distinct affine blocks have even intersection and each block "
            "has odd self-intersection because M^T*M=I over F_2"
        ),
        "projections_determine_word": True,
        "reconstruction_reason": (
            "for each midpoint [a], the vector of block parities is its "
            "difference fibre multiplied by the invertible matrix M"
        ),
        "no_affine_projection_characterization": (
            "if q_C=0 for every C, w is exactly the disjoint XOR of those "
            "rectangles A_K x B_(K,beta) for which r_C=A_K"
        ),
        "fixed_transverse_rectangles_pairwise_disjoint": True,
        "proved": True,
    }


def mobius_block_intersection_theorem(
    p: int, half_count: int
) -> dict[str, object]:
    """Return the exact cover bounds for a union of localized halves.

    The application only needs Paley primes.  Orientations and ternarity are
    irrelevant here: if ``U`` is the actual support after cancellations,
    then merely ``U`` being a subset of the union of the half supports is
    used.  Cancellations can only make the cover smaller.
    """
    h = _check_paley_prime(p)
    if (
        not isinstance(half_count, int)
        or isinstance(half_count, bool)
        or half_count < 0
    ):
        raise ValueError("half_count must be a nonnegative integer")
    d = p + 1
    rectangle_size = p * h
    return {
        "p": p,
        "h": h,
        "d": d,
        "half_count": half_count,
        "assumption_on_actual_support": (
            "U is the actual post-cancellation support and "
            "U is contained in the union of the half supports"
        ),
        "ternarity_used": False,
        "orientation_used": False,
        "midpoint_parameter": "z=t+1 in F_p^*",
        "midpoint_conic": "L(a)*(L-M)(a)=j^2/4",
        "difference_fibre_formula": (
            "if K=x*L+y*M, then "
            "K(delta_z)=j*((x+y)*(z-2)+y/z)/2"
        ),
        "one_half_midpoints_on_origin_line_at_most": 1,
        "one_half_midpoints_in_affine_block_at_most": 2,
        "one_half_columns_over_difference_block_at_most": 4,
        "union_columns_over_difference_block_at_most": 4 * half_count,
        "one_half_intersection_with_fixed_transverse_rectangle_at_most": 2,
        "union_intersection_with_fixed_transverse_rectangle_at_most": (
            2 * half_count
        ),
        "fixed_transverse_rectangle_size": rectangle_size,
        "fixed_transverse_rectangle_cannot_be_contained": (
            2 * half_count < rectangle_size
        ),
        "proof_of_four": (
            "K(delta_z)^2=beta splits into two nonzero quadratic "
            "equations in z"
        ),
        "proof_of_midpoint_bounds": (
            "the nonzero conic meets an origin line in at most one "
            "antipodal class and a paired nonorigin affine line in at "
            "most two"
        ),
        "proved": True,
    }


def common_block_resultant_theorem(p: int) -> dict[str, object]:
    """Return the exact per-half locus for two common-block midpoint classes.

    Write ``K=x*L+y*M``, ``A=x+y``, and choose square roots
    ``r^2=alpha``, ``s^2=beta``.  Eliminating the Mobius parameter from
    ``K(a)^2=alpha`` and ``K(delta)^2=beta`` gives one quadratic resultant.
    The intrinsic form describes the two midpoint classes without reference
    to the coefficients ``x,y``.
    """
    _check_odd_prime(p)
    return {
        "p": p,
        "coefficient_setup": "K=x*L+y*M and A=x+y",
        "square_root_setup": "r^2=alpha and s^2=beta, with r,s nonzero",
        "signed_resultant": (
            "alpha-beta-A*x*j^2=2*eta*A*j*s for eta in {+1,-1}"
        ),
        "sign_free_resultant": (
            "(alpha-beta-A*x*j^2)^2-4*A^2*j^2*beta=0"
        ),
        "parameter_candidates": (
            "z_epsilon=(A*j+eta*s+epsilon*r)/(A*j), "
            "epsilon in {+1,-1}"
        ),
        "two_distinct_class_locus": "A*y*(A*j+eta*s)!=0",
        "degenerate_cases": (
            "A=0 or y=0 gives at most one midpoint class; "
            "A*j+eta*s=0 identifies the two candidates antipodally"
        ),
        "intrinsic_hypotheses": (
            "K and L independent; represent two distinct classes by "
            "a1,a2 with K(a1)=K(a2)=r and put l_k=L(a_k)"
        ),
        "intrinsic_pair_criterion": (
            "l1*l2*(l1+l2)!=0 and, for one eta, "
            "(r-eta*s)*l1-(r+eta*s)*l2=r*j"
        ),
        "unique_auxiliary": (
            "M(a_k)=l_k-j^2/(4*l_k), k=1,2"
        ),
        "corresponding_parameters": "z1=2*l1/j and z2=-2*l2/j",
        "resulting_difference_values": "K(delta_z1)=K(delta_z2)=eta*s",
        "proved": True,
    }


def branch_c_all_halves_cover_obstruction(p: int) -> dict[str, object]:
    """State the forced cover profile for ``m=(p+1)/2`` Mobius halves.

    This is a necessary condition for a nonzero halved dual support to be
    contained in the actual used support.  It deliberately does not assert
    that the forced all-halves cover is impossible.
    """
    h = _check_paley_prime(p, minimum=31)
    d = p + 1
    half_count = d // 2
    projection = dual_block_projection_theorem(p)
    cover = mobius_block_intersection_theorem(p, half_count)
    delta_size = d * h
    raw_half_occurrences = half_count * (p - 1)
    projection_capacity = 4 * half_count
    proved = bool(
        half_count == h + 1
        and raw_half_occurrences == delta_size
        and projection_capacity == 2 * p + 2
        and 2 * half_count < p * h
        and 2 * p + h > projection_capacity
        and projection["proved"]
        and cover["proved"]
    )
    if not proved:
        raise ArithmeticError("the branch-C cover inequalities changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "Mobius_half_count": half_count,
        "raw_half_occurrences": raw_half_occurrences,
        "delta_size": delta_size,
        "actual_support_size_at_most_delta": True,
        "hypothesis": (
            "0!=w in Row(D), with supp(w) contained in the actual support "
            "U of the chosen halves"
        ),
        "first_conclusion": (
            "some difference block C has an affine midpoint block in r_C"
        ),
        "reason_all_origin_cells_are_impossible": (
            "otherwise w is a disjoint sum of fixed-transverse rectangles, "
            "one of which would have to be contained in U"
        ),
        "projection_capacity": projection_capacity,
        "possible_nonzero_projection_profiles": [
            "A_K",
            "B_(K,alpha)",
            "A_K union B_(K,alpha)",
            "B_(K,alpha_1) union B_(K,alpha_2)",
        ],
        "at_most_two_affine_cells": True,
        "two_affine_cells_exclude_origin_cell": True,
        "all_halves_cover_conclusion": (
            "for every affine block appearing in this r_C, all h+1 halves "
            "must supply a midpoint class in that block over C"
        ),
        "halves_supplying_two_classes_at_least": h,
        "two_affine_block_total_incidence_slack_at_most": 2,
        "ternarity_needed_for_cover_theorem": False,
        "actual_support_note": (
            "orientations and cancellations only shrink U relative to the "
            "union, so the necessary cover conclusions remain valid"
        ),
        "all_halves_cover_impossible_proved": False,
        "structured_punctured_map_surjective_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def _antipodal_class(p: int, point: Point) -> Point:
    point = (point[0] % p, point[1] % p)
    if point == (0, 0):
        raise ValueError("zero has no antipodal nonzero class")
    negative = ((-point[0]) % p, (-point[1]) % p)
    return min(point, negative)


def _antipodal_classes(p: int) -> tuple[Point, ...]:
    return tuple(
        sorted(
            {
                _antipodal_class(p, point)
                for point in product(range(p), repeat=2)
                if point != (0, 0)
            }
        )
    )


def _square_values(p: int) -> tuple[int, ...]:
    return tuple(sorted({value * value % p for value in range(1, p)}))


def _origin_cell(
    p: int, functional: Functional, classes: tuple[Point, ...]
) -> frozenset[Point]:
    return frozenset(
        point
        for point in classes
        if _functional_value(p, functional, point) == 0
    )


def _affine_block(
    p: int,
    functional: Functional,
    beta: int,
    classes: tuple[Point, ...],
) -> frozenset[Point]:
    return frozenset(
        point
        for point in classes
        if _functional_value(p, functional, point) ** 2 % p == beta % p
    )


def _midpoint_difference(p: int, edge: Edge) -> tuple[Point, Point]:
    inverse_two = pow(2, -1, p)
    midpoint = (
        (edge[0][0] + edge[1][0]) * inverse_two % p,
        (edge[0][1] + edge[1][1]) * inverse_two % p,
    )
    difference = (
        (edge[1][0] - edge[0][0]) * inverse_two % p,
        (edge[1][1] - edge[0][1]) * inverse_two % p,
    )
    return _antipodal_class(p, midpoint), _antipodal_class(p, difference)


def _functional_from_two_values(
    p: int,
    first_point: Point,
    first_value: int,
    second_point: Point,
    second_value: int,
) -> Functional:
    determinant = (
        first_point[0] * second_point[1]
        - first_point[1] * second_point[0]
    ) % p
    if determinant == 0:
        raise ValueError("the two interpolation points must be independent")
    inverse = pow(determinant, -1, p)
    first_coefficient = (
        first_value * second_point[1]
        - first_point[1] * second_value
    ) * inverse % p
    second_coefficient = (
        first_point[0] * second_value
        - first_value * second_point[0]
    ) * inverse % p
    return first_coefficient, second_coefficient


def _basis_coefficients(
    p: int, first: Functional, second: Functional, target: Functional
) -> tuple[int, int]:
    determinant = (
        first[0] * second[1] - first[1] * second[0]
    ) % p
    if determinant == 0:
        raise ValueError("the functionals must be independent")
    inverse = pow(determinant, -1, p)
    first_coefficient = (
        target[0] * second[1] - target[1] * second[0]
    ) * inverse % p
    second_coefficient = (
        first[0] * target[1] - first[1] * target[0]
    ) * inverse % p
    return first_coefficient, second_coefficient


def distinct_direction_saturated_cover_counterexample(
    p: int,
) -> dict[str, object]:
    """Construct an all-halves common-block cover with distinct targets.

    This refutes only an obstruction based on distinctness of the target
    directions.  It chooses the star centers and auxiliaries.  It does not
    impose the prescribed branch-C centers, mutual ternarity, or containment
    of an entire dual support.

    The explicit witness uses ``K=(1,0)``, ``J=(0,1)``, and
    ``alpha=beta=1``.  The block ``B_(K,1)`` is represented by the points
    ``a(q)=(1,q)``.  We cover the field by ``h+1`` pairs, choose distinct
    slopes ``u_i`` outside at most four forbidden values per pair, and put
    ``L_i=J+u_i*K``.  The intrinsic criterion uniquely supplies ``M_i``.
    """
    h = _check_paley_prime(p, minimum=11)
    inverse_two = pow(2, -1, p)
    inverse_four = pow(4, -1, p)
    K: Functional = (1, 0)
    J: Functional = (0, 1)
    pairs = [(2 * index, 2 * index + 1) for index in range(h)]
    pairs.append((p - 1, 0))
    used_slopes: set[int] = set()
    covered: set[Point] = set()
    records: list[dict[str, object]] = []

    for first_q, second_q in pairs:
        forbidden = {
            (-first_q) % p,
            (-second_q) % p,
            (-(first_q + second_q) * inverse_two) % p,
        }
        slope = next(
            (
                candidate
                for candidate in range(p)
                if candidate not in forbidden and candidate not in used_slopes
            ),
            None,
        )
        if slope is None:
            raise ArithmeticError("the greedy distinct-direction choice failed")
        used_slopes.add(slope)

        first_point = (1, first_q)
        second_point = (1, second_q)
        direction: Functional = (slope, 1)
        first_l = (first_q + slope) % p
        second_l = (second_q + slope) % p
        center = (-2 * second_l) % p
        if (
            first_l == 0
            or second_l == 0
            or (first_l + second_l) % p == 0
            or center == 0
        ):
            raise ArithmeticError("a forbidden value entered the construction")

        first_m = (
            first_l - center * center * inverse_four * pow(first_l, -1, p)
        ) % p
        second_m = (
            second_l - center * center * inverse_four * pow(second_l, -1, p)
        ) % p
        auxiliary = _functional_from_two_values(
            p,
            first_point,
            first_m,
            second_point,
            second_m,
        )
        determinant = (
            direction[0] * auxiliary[1]
            - direction[1] * auxiliary[0]
        ) % p
        if determinant == 0:
            raise ArithmeticError("the constructed auxiliary became dependent")

        x, y = _basis_coefficients(p, direction, auxiliary, K)
        A = (x + y) % p
        signed_resultant = bool(
            (-A * x * center * center) % p
            == (2 * A * center) % p
        )
        nondegenerate = bool(
            A != 0 and y != 0 and (A * center + 1) % p != 0
        )

        edges = mobius_parameter_edges(
            p, direction, auxiliary, center
        )
        hits: set[Point] = set()
        for edge in edges.values():
            midpoint, difference = _midpoint_difference(p, edge)
            if (
                _functional_value(p, K, midpoint) ** 2 % p == 1
                and _functional_value(p, K, difference) ** 2 % p == 1
            ):
                hits.add(midpoint)
        assigned = {
            _antipodal_class(p, first_point),
            _antipodal_class(p, second_point),
        }
        if hits != assigned or not signed_resultant or not nondegenerate:
            raise ArithmeticError("the constructed common-block pair changed")
        covered.update(hits)
        records.append(
            {
                "assigned_q_values": [first_q, second_q],
                "target_direction": list(direction),
                "center": center,
                "auxiliary": list(auxiliary),
                "resultant_holds": signed_resultant,
                "two_class_locus_nondegenerate": nondegenerate,
                "common_block_midpoint_classes": [
                    list(point) for point in sorted(hits)
                ],
            }
        )

    classes = _antipodal_classes(p)
    target_block = _affine_block(p, K, 1, classes)
    incidence_count = sum(
        len(record["common_block_midpoint_classes"]) for record in records
    )
    proved = bool(
        len(records) == h + 1
        and len(used_slopes) == h + 1
        and covered == set(target_block)
        and incidence_count == p + 1
        and incidence_count - len(covered) == 1
    )
    if not proved:
        raise ArithmeticError("the saturated common-block cover changed")
    return {
        "p": p,
        "h": h,
        "target_block": "B_((1,0),1)",
        "difference_block": "B_((1,0),1)",
        "distinct_target_direction_count": len(used_slopes),
        "Mobius_half_count": len(records),
        "covered_midpoint_classes": len(covered),
        "target_block_size": p,
        "common_block_incidences": incidence_count,
        "duplicate_incidences": incidence_count - len(covered),
        "all_target_directions_distinct": True,
        "all_auxiliaries_independent": True,
        "all_halves_supply_two_classes": True,
        "records": records,
        "counterexample_scope": (
            "distinct target directions alone do not obstruct a saturated "
            "common B over C cover"
        ),
        "prescribed_centers_respected": False,
        "mutual_ternarity_proved": False,
        "full_dual_support_containment_proved": False,
        "structured_punctured_map_surjective_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def exact_small_cover_replay(p: int) -> dict[str, object]:
    """Check the incidence and one-half formulas at ``p=3,7`` only."""
    h = _check_paley_prime(p)
    if p not in (3, 7):
        raise ValueError("the exact cover replay is limited to p=3 or p=7")
    classes = _antipodal_classes(p)
    directions = projective_functionals(p)
    squares = _square_values(p)
    origin_cells = {
        direction: _origin_cell(p, direction, classes)
        for direction in directions
    }
    blocks = {
        (direction, beta): _affine_block(
            p, direction, beta, classes
        )
        for direction in directions
        for beta in squares
    }

    partitions_exact = all(
        origin_cells[direction].isdisjoint(blocks[direction, beta])
        and len(origin_cells[direction]) == h
        and len(blocks[direction, beta]) == p
        for direction in directions
        for beta in squares
    ) and all(
        origin_cells[direction].union(
            *(blocks[direction, beta] for beta in squares)
        )
        == frozenset(classes)
        for direction in directions
    )
    gram_exact = all(
        len(first & second) % 2 == int(first_key == second_key)
        for first_key, first in blocks.items()
        for second_key, second in blocks.items()
    )
    parallel_projection_exact = all(
        len(origin_cells[row] & block) % 2
        == int(row != block_key[0])
        for row in directions
        for block_key, block in blocks.items()
    )
    origin_cells_disjoint = all(
        not (origin_cells[first] & origin_cells[second])
        for first_index, first in enumerate(directions)
        for second in directions[first_index + 1 :]
    )

    chosen_direction = (1, 0)
    chosen_auxiliary = (0, 1)
    edges = mobius_parameter_edges(
        p, chosen_direction, chosen_auxiliary, 1
    )
    columns = tuple(
        _midpoint_difference(p, edge) for edge in edges.values()
    )
    midpoint_counts = Counter(midpoint for midpoint, _difference in columns)
    conic_exact = all(
        _functional_value(p, chosen_direction, midpoint)
        * (
            _functional_value(p, chosen_direction, midpoint)
            - _functional_value(p, chosen_auxiliary, midpoint)
        )
        % p
        == pow(4, -1, p)
        for midpoint in midpoint_counts
    )
    max_origin_midpoints = max(
        len(set(midpoint_counts) & cell) for cell in origin_cells.values()
    )
    max_affine_midpoints = max(
        len(set(midpoint_counts) & block) for block in blocks.values()
    )
    max_columns_over_difference_block = max(
        sum(difference in block for _midpoint, difference in columns)
        for block in blocks.values()
    )
    proved = bool(
        len(classes) == (p + 1) * h
        and partitions_exact
        and gram_exact
        and parallel_projection_exact
        and origin_cells_disjoint
        and len(columns) == p - 1
        and set(midpoint_counts.values()) == {2}
        and conic_exact
        and max_origin_midpoints <= 1
        and max_affine_midpoints <= 2
        and max_columns_over_difference_block <= 4
    )
    if not proved:
        raise ArithmeticError("the exact small cover replay changed")
    return {
        "p": p,
        "delta_size": len(classes),
        "parallel_partitions_exact": partitions_exact,
        "affine_block_gram_identity_exact": gram_exact,
        "parallel_generator_projection_exact": parallel_projection_exact,
        "distinct_origin_cells_disjoint": origin_cells_disjoint,
        "Mobius_half_columns": len(columns),
        "midpoint_class_multiplicity_set": sorted(set(midpoint_counts.values())),
        "midpoint_conic_exact": conic_exact,
        "observed_max_origin_line_midpoints": max_origin_midpoints,
        "observed_max_affine_block_midpoints": max_affine_midpoints,
        "observed_max_columns_over_difference_block": (
            max_columns_over_difference_block
        ),
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    projection = dual_block_projection_theorem(p)
    cover = mobius_block_intersection_theorem(p, (p + 1) // 2)
    resultant = common_block_resultant_theorem(p)
    obstruction = branch_c_all_halves_cover_obstruction(p)
    distinct_direction_cover = distinct_direction_saturated_cover_counterexample(p)
    replays = {str(q): exact_small_cover_replay(q) for q in (3, 7)}
    proved = bool(
        projection["proved"]
        and cover["proved"]
        and resultant["proved"]
        and obstruction["proved"]
        and distinct_direction_cover["proved"]
        and all(replay["proved"] for replay in replays.values())
    )
    if not proved:
        raise ArithmeticError("the Mobius cover theorem record changed")
    return {
        "title": "Affine-block projections of Mobius-contained dual words",
        "status": (
            "ALL-F CASE EXCLUDED; ANY CONTAINED DUAL WORD FORCES AN "
            "ALL-HALVES AFFINE-BLOCK COVER; IMPOSSIBILITY OPEN"
        ),
        "projection_theorem": projection,
        "Mobius_intersection_theorem": cover,
        "per_half_common_block_resultant": resultant,
        "branch_C_cover_obstruction": obstruction,
        "distinct_direction_saturated_cover": distinct_direction_cover,
        "small_exact_replays": replays,
        "proved_all_claimed_statements": proved,
        "contained_nonzero_dual_word_excluded": False,
        "prescribed_center_cover_excluded": False,
        "mutual_ternarity_of_saturated_cover_proved": False,
        "structured_punctured_map_surjective": False,
        "residual_ii_closed": False,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
