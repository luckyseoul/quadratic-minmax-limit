#!/usr/bin/env python3
"""Prescribed-center obstruction for an equal-square common-block cover.

This module records a symbolic consequence of the intrinsic pair criterion
from ``e1_gmin_m4_symmetric_halved_mobius_cover``.  It concerns only a fixed
choice of common midpoint/difference direction ``K`` and the slice
``alpha=beta``.  It does not exclude a cover for another ``K`` or for
``alpha!=beta``.

The p=31 record is a hard-coded exact witness, not a search or census.
"""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Hashable, Iterable
from itertools import combinations

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Functional,
    Point,
    _functional_value,
    _point_from_coordinates,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import paley_direction_sign
from e1_gmin_m4_prop15721 import is_prime


def _check_odd_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")


def _canonical_functional(p: int, functional: Functional) -> Functional:
    first, second = (coordinate % p for coordinate in functional)
    if first:
        return 1, second * pow(first, -1, p) % p
    if second:
        return 0, 1
    raise ValueError("a projective functional must be nonzero")


def alpha_equal_beta_anchor_hall_theorem(p: int) -> dict[str, object]:
    """Return the exact anchor and Hall criterion for fully doubled covers.

    Fix ``K,L`` independent, ``r!=0``, and identify the p classes in
    ``B_(K,r^2)`` by their representatives satisfying ``K(a)=r``.  For a
    prescribed nonzero center ``j``, let ``z,a^+,a^-`` be the unique points
    on this affine line at which ``L`` has values ``0,j/2,-j/2``.

    The intrinsic criterion says that a half supplies two distinct classes
    over the same block ``C=B_(K,r^2)`` exactly when its pair consists of one
    anchor in ``{a^+,a^-}`` and one point outside
    ``{z,a^+,a^-}``.

    For m=(p+1)/2 fully doubled halves, create an anchor slot and a free slot
    for every half.  A saturated cover has p+1 incidences on p points, hence
    is equivalent to a perfect matching after duplicating one right vertex.
    The returned Hall inequality is therefore necessary and sufficient.
    """
    _check_odd_prime(p)
    return {
        "p": p,
        "half_count": (p + 1) // 2,
        "slice": "alpha=beta=r^2 with r!=0",
        "anchor_values": "L(a)=+j/2 or L(a)=-j/2",
        "zero_value": "L(z)=0",
        "valid_double_pair": (
            "one point in A_i={a_i^+,a_i^-} and one point in "
            "X\\T_i, where T_i={z_i,a_i^+,a_i^-}"
        ),
        "eta_plus_reduction": "L(a_2)=-j/2",
        "eta_minus_reduction": "L(a_1)=+j/2",
        "slot_neighborhoods": (
            "N(i,anchor)=A_i and N(i,free)=X\\T_i"
        ),
        "fully_doubled_hall_criterion": (
            "there is d in X such that for every P,Q subset I, "
            "|union_(i in P) A_i union (X\\intersection_(i in Q) T_i)|"
            "+1_[d lies in that union] >= |P|+|Q|; for Q=empty omit "
            "the X\\intersection term"
        ),
        "criterion_is_necessary_and_sufficient": True,
        "anchor_necessary_condition": (
            "|union_(i in P) A_i| >= |P|-1 for every P subset I"
        ),
        "deficiency_two_obstruction": (
            "if |union_(i in P) A_i| <= |P|-2 for some P, then neither "
            "the all-doubled profile nor the one-single profile can cover X"
        ),
        "one_single_reason": (
            "at most one member of P can be the single half; the remaining "
            "|P|-1 doubled halves already repeat an anchor, while the total "
            "incidence count is exactly |X|"
        ),
        "proved": True,
    }


def anchor_graph_pseudoforest_profile(
    edges: Iterable[tuple[Hashable, Hashable]],
) -> dict[str, object]:
    """Classify a finite simple anchor graph component by component.

    An edge-to-incident-vertex system has an SDR exactly when every graph
    component has at most one cycle, equivalently ``edge_count<=vertex_count``
    in every component.  This helper checks the finite graph identity used by
    ``dependent_anchor_line_cover_theorem``; it is not proof by enumeration.
    """
    edge_list = tuple(edges)
    if any(first == second for first, second in edge_list):
        raise ValueError("anchor edges must have two distinct endpoints")
    frozen_edges = tuple(frozenset(edge) for edge in edge_list)
    if len(set(frozen_edges)) != len(frozen_edges):
        raise ValueError("the anchor graph must be simple")

    adjacency: defaultdict[Hashable, set[Hashable]] = defaultdict(set)
    for first, second in edge_list:
        adjacency[first].add(second)
        adjacency[second].add(first)

    seen: set[Hashable] = set()
    components: list[dict[str, int | bool]] = []
    for start in adjacency:
        if start in seen:
            continue
        queue: deque[Hashable] = deque([start])
        seen.add(start)
        vertices: set[Hashable] = set()
        degree_sum = 0
        while queue:
            vertex = queue.popleft()
            vertices.add(vertex)
            degree_sum += len(adjacency[vertex])
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        edge_count = degree_sum // 2
        vertex_count = len(vertices)
        components.append(
            {
                "vertex_count": vertex_count,
                "edge_count": edge_count,
                "cycle_rank": edge_count - vertex_count + 1,
                "at_most_unicyclic": edge_count <= vertex_count,
            }
        )

    pseudoforest = all(
        bool(component["at_most_unicyclic"]) for component in components
    )
    return {
        "edge_count": len(edge_list),
        "vertex_count": len(adjacency),
        "component_count": len(components),
        "components": components,
        "edge_to_incident_vertex_SDR_exists": pseudoforest,
        "pseudoforest": pseudoforest,
        "proved_graph_equivalence": True,
    }


def dependent_anchor_line_cover_theorem(p: int) -> dict[str, object]:
    """Give a one-way all-prime construction from a good anchor line.

    Assume every prescribed hard center ``j_i`` is nonzero.  For a chosen
    hard target ``(L_0,j_0)``, take ``K=L_0`` and
    ``alpha=beta=(j_0/2)^2``.  On the representative affine line
    ``ell={L_0=j_0/2}``, the dependent half supplies an arbitrary singleton
    at parameter ``z=1``.  The other ``h=(p-1)/2`` halves have the anchor and
    free-slot description from ``alpha_equal_beta_anchor_hall_theorem``.

    If their simple anchor graph on ``ell`` is a pseudoforest, its edges have
    an SDR.  Hall then always matches the h free slots into h of the remaining
    h+1 points; the dependent singleton fills the last point.

    This is deliberately one-way.  It does not prove that one of the
    prescribed anchor lines always has a pseudoforest link graph.
    """
    _check_odd_prime(p)
    if p < 31 or p % 4 != 3:
        raise ValueError("need a branch prime p=3 mod 4 with p>=31")
    h = (p - 1) // 2
    hard_count = h + 1
    if h <= 3:
        raise ArithmeticError("the all-h triple-intersection bound changed")
    return {
        "p": p,
        "h": h,
        "hard_half_count": hard_count,
        "hypothesis": "every prescribed hard center j_i is nonzero",
        "chosen_base_target": "(L_0,j_0) with j_0!=0",
        "common_direction": "K=L_0",
        "common_slice": "alpha=beta=(j_0/2)^2",
        "representative_anchor_line": "ell={a:L_0(a)=j_0/2}",
        "dependent_half_parameter": "z=1 (equivalently t=0)",
        "dependent_half_unique_common_parameter": True,
        "dependent_half_uniqueness_reason": (
            "K(a_z)^2=alpha gives z^2=1, while "
            "K(delta_z)^2=beta gives (z-2)^2=1; their only common "
            "solution in odd characteristic is z=1"
        ),
        "dependent_half_edge": "{0,2a}, with midpoint a and difference -a",
        "dependent_half_freedom": (
            "every a in ell is available by choosing an auxiliary M_0 "
            "with M_0(a)=0"
        ),
        "dependent_half_is_single": True,
        "other_half_count": h,
        "other_half_anchor_sets": (
            "A_i=ell intersection {L_i=+j_i/2 or L_i=-j_i/2}"
        ),
        "other_half_forbidden_triples": (
            "T_i=ell intersection {L_i=0,+j_i/2,-j_i/2}"
        ),
        "zero_points_distinct": True,
        "zero_points_reason": (
            "a common zero would be a nonzero vector annihilated by two "
            "distinct projective functionals"
        ),
        "forbidden_triples_are_distinct_APs": True,
        "distinct_AP_reason": (
            "T_i is a nondegenerate three-term affine progression with "
            "unique center z_i in characteristic other than 3"
        ),
        "anchor_graph_simple": True,
        "anchor_SDR_iff_pseudoforest": True,
        "pseudoforest_definition": (
            "every connected component is a tree or is unicyclic"
        ),
        "remaining_point_count_after_anchor_SDR": h + 1,
        "free_slot_neighborhood": "R\\T_i",
        "free_slot_Hall_bounds": {
            "q_at_most_h_minus_2": (
                "|union N_i|>=h+1-3=h-2>=q"
            ),
            "q_equal_h_minus_1": (
                "distinct triples have common intersection size at most 2, "
                "so |union N_i|>=h-1"
            ),
            "q_equal_h": (
                "at most three distinct nondegenerate three-term APs can "
                "contain a fixed pair; h>3, so the common intersection has "
                "size at most 1 and |union N_i|>=h"
            ),
        },
        "free_slot_matching_exists": True,
        "cover_profile": "h doubled halves and one dependent singleton",
        "covered_midpoint_class_count": p,
        "conditional_conclusion": (
            "if one prescribed anchor line has a pseudoforest link graph, "
            "then a saturated equal-square common-block incidence cover exists"
        ),
        "some_good_base_line_proved": False,
        "converse_claimed": False,
        "mutual_ternarity_proved": False,
        "full_dual_support_containment_proved": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def p31_prescribed_center_anchor_obstruction() -> dict[str, object]:
    """Return an exact p=31 deficiency-two hard-center witness.

    Put ``K=(1,4)``, ``J=(1,0)``, and represent the common affine block by
    ``a(q)`` with ``K(a(q))=1`` and ``J(a(q))=q``.  The six edges of the
    four-point set ``{0,2,11,14}`` have distinct midpoints
    ``{1,7,8,21,22,28}``.  Their annihilator directions are all Paley-hard.
    Prescribing ``j=2L(a(u))`` makes ``a(u),a(v)`` the two anchors for the
    edge ``{u,v}``.  Thus six anchor sets have union four, a deficiency-two
    Hall obstruction.

    The other ten Paley-hard centers are set to one.  This gives a complete
    prescribed nonzero hard-center list, but the conclusion concerns only
    the displayed fixed ``K`` and equal-square slice.
    """
    p = 31
    inverse_two = pow(2, -1, p)
    K: Functional = (1, 4)
    J: Functional = (1, 0)
    vertices = (0, 2, 11, 14)
    edges = tuple(combinations(vertices, 2))
    expected_midpoints = {1, 7, 8, 21, 22, 28}
    hard_directions = tuple(
        direction
        for direction in projective_functionals(p)
        if paley_direction_sign(p, direction) == 1
    )
    centers = {direction: 1 for direction in hard_directions}
    records: list[dict[str, object]] = []

    def block_point(q: int) -> Point:
        return _point_from_coordinates(p, K, J, 1, q % p)

    for first_q, second_q in edges:
        midpoint_q = (first_q + second_q) * inverse_two % p
        midpoint = block_point(midpoint_q)
        direction = _canonical_functional(
            p, ((-midpoint[1]) % p, midpoint[0])
        )
        first_value = _functional_value(
            p, direction, block_point(first_q)
        )
        second_value = _functional_value(
            p, direction, block_point(second_q)
        )
        center = 2 * first_value % p
        if (
            midpoint_q not in expected_midpoints
            or direction not in hard_directions
            or paley_direction_sign(p, direction) != 1
            or _functional_value(p, direction, midpoint) != 0
            or first_value == 0
            or second_value != -first_value % p
            or center == 0
        ):
            raise ArithmeticError("the prescribed-center witness changed")
        centers[direction] = center
        records.append(
            {
                "anchor_q_values": [first_q, second_q],
                "zero_q_value": midpoint_q,
                "hard_direction": list(direction),
                "center": center,
                "anchor_L_values": [first_value, second_value],
            }
        )

    anchor_union = {
        q
        for record in records
        for q in record["anchor_q_values"]
    }
    midpoint_set = {int(record["zero_q_value"]) for record in records}
    directions_used = {
        tuple(record["hard_direction"]) for record in records
    }
    all_centers_nonzero = all(value % p for value in centers.values())
    proved = bool(
        paley_direction_sign(p, K) == -1
        and len(hard_directions) == 16
        and len(centers) == 16
        and all_centers_nonzero
        and len(records) == 6
        and len(directions_used) == 6
        and midpoint_set == expected_midpoints
        and anchor_union == set(vertices)
        and len(records) - len(anchor_union) == 2
    )
    if not proved:
        raise ArithmeticError("the p=31 Hall obstruction changed")
    return {
        "p": p,
        "common_direction_K": list(K),
        "coordinate_functional_J": list(J),
        "alpha": 1,
        "beta": 1,
        "K_is_Paley_opposite": True,
        "hard_direction_count": len(hard_directions),
        "all_prescribed_hard_centers_nonzero": all_centers_nonzero,
        "prescribed_centers": [
            {"hard_direction": list(direction), "center": centers[direction]}
            for direction in hard_directions
        ],
        "four_anchor_q_values": list(vertices),
        "six_hard_half_records": records,
        "six_distinct_zero_q_values": sorted(midpoint_set),
        "anchor_family_size": len(records),
        "anchor_union_size": len(anchor_union),
        "anchor_deficiency": len(records) - len(anchor_union),
        "fixed_K_equal_square_saturated_cover_impossible": True,
        "covers_both_incidence_profiles": (
            "all 16 halves doubled, or exactly 15 doubled and one single"
        ),
        "another_K_or_unequal_square_cover_excluded": False,
        "full_dual_support_excluded": False,
        "residual_ii_closed": False,
        "witness_role": "direct exact arithmetic, not a prime scan",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Combine the general criterion with the fixed p=31 witness."""
    theorem = alpha_equal_beta_anchor_hall_theorem(p)
    dependent = dependent_anchor_line_cover_theorem(p)
    out: dict[str, object] = {
        "p": p,
        "anchor_hall_theorem": theorem,
        "dependent_anchor_line_theorem": dependent,
        "proved_all_claimed_statements": bool(
            theorem["proved"] and dependent["proved"]
        ),
        "arbitrary_prescribed_centers_always_cover_fixed_K": False,
        "existence_after_varying_K_excluded": False,
        "residual_ii_closed": False,
    }
    if p == 31:
        witness = p31_prescribed_center_anchor_obstruction()
        out["p31_witness"] = witness
        out["proved_all_claimed_statements"] = bool(
            theorem["proved"] and dependent["proved"] and witness["proved"]
        )
    return out


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
