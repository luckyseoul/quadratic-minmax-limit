#!/usr/bin/env python3
"""The exact common-origin ledger at the p=31 branch-C top endpoint.

This module combines the already proved Mobius parallel-parity identity with
the elementary geometry of the ``t=0`` edge of a localized half.  At
``p=31, t=177`` the ``j=0`` Hamming endpoint has one fixed antipodal edge,
no unused doubled orbit, and one cancellation unit among the sixteen
localized halves.

The conclusion is narrow: all sixteen auxiliary projective directions are
distinct, so no cancellation can occur on a common-origin orbit.  The graph
contains sixteen origin edges, and their Paley-signed degree is determined
by the sign of the fixed-edge direction.  This is a necessary condition for
the localized-half top construction, not a common-graph construction and not
a closure of residual (ii).
"""

from __future__ import annotations

import json
from typing import Literal

from e1_gmin_m4_inversion_antisymmetric_radon import projective_functionals
from e1_gmin_m4_mobius_half_symmetric import paley_direction_sign


Point = tuple[int, int]
FixedType = Literal["hard", "opposite"]

P = 31
M = 16
T_TOP = 177
GRAPH_EDGES = 479
RAW_HALF_EDGES = M * (P - 1)
USED_NONFIXED_ORBITS = GRAPH_EDGES - 1


def _negative(point: Point, p: int = P) -> Point:
    return (-point[0] % p, -point[1] % p)


def _antipodal_class(point: Point, p: int = P) -> Point:
    if point == (0, 0):
        raise ValueError("zero has no nonzero antipodal class")
    return min(point, _negative(point, p))


def kernel_representative(functional: Point, p: int = P) -> Point:
    """Return a nonzero vector spanning the kernel of ``functional``."""
    a, b = functional
    if a % p == 0 and b % p == 0:
        raise ValueError("a projective functional must be nonzero")
    return (b % p, -a % p)


def origin_line_injectivity_certificate(p: int = P) -> dict[str, object]:
    """Check the implication needed to separate the actual origin edges.

    In the normalized localized-half parameterization, the ``t=0`` edge is
    ``{0,u}`` with ``u`` a nonzero vector in ``ker(M_i)``.  Hence two such
    edges agreeing modulo central inversion forces their auxiliary
    projective directions to agree.  The converse is not asserted: the
    actual magnitude of ``u`` along a fixed kernel line also depends on the
    target, center, and auxiliary scale.  Also ``Q(u)=Q(M_i)`` up to a
    square, so the edge's Paley sign is the hard/opposite sign of ``M_i``.
    """
    directions = projective_functionals(p)
    orbit_by_direction = tuple(
        _antipodal_class(kernel_representative(direction, p), p)
        for direction in directions
    )
    edge_signs = tuple(
        paley_direction_sign(p, kernel_representative(direction, p))
        for direction in directions
    )
    auxiliary_signs = tuple(
        paley_direction_sign(p, direction) for direction in directions
    )
    distinct_kernel_directions = len(set(orbit_by_direction)) == len(directions)
    sign_match = edge_signs == auxiliary_signs
    proved = bool(
        len(directions) == p + 1
        and distinct_kernel_directions
        and sign_match
    )
    if not proved:
        raise ArithmeticError("the auxiliary/origin-orbit geometry changed")
    return {
        "p": p,
        "projective_auxiliary_directions": len(directions),
        "distinct_projective_kernel_directions": len(set(orbit_by_direction)),
        "auxiliary_to_kernel_direction_bijection": distinct_kernel_directions,
        "equal_origin_orbits_force_equal_auxiliary_direction": (
            distinct_kernel_directions
        ),
        "origin_edge_sign_equals_auxiliary_sign": sign_match,
        "proved": proved,
    }


def p31_top_origin_endpoint(fixed_edge_type: FixedType) -> dict[str, object]:
    """Return the two exact signed-origin cases at ``p=31,t=177``.

    Let ``v=P mod 2 + 1_hard`` be the parity target from the established
    parallel-slice theorem and let ``c`` be the parity vector of the sixteen
    auxiliary directions.  The top endpoint has a single fixed direction
    ``F`` and no doubled column, hence ``c=v+e_F`` over ``F_2``.
    """
    if fixed_edge_type not in ("hard", "opposite"):
        raise ValueError("fixed_edge_type must be 'hard' or 'opposite'")

    # Hard quotas are 14^14,15^2, so v is supported on the fourteen even
    # hard quotas.  Opposite quotas are 15^3,16^13, so v is supported on the
    # three odd opposite quotas.
    hard_v_support = 14
    opposite_v_support = 3
    v_weight = hard_v_support + opposite_v_support

    # If F were outside supp(v), wt(v+e_F)=18, exceeding the number of
    # auxiliary occurrences.  Thus F is in supp(v), and deleting it leaves
    # a parity vector of weight sixteen.  Sixteen occurrences realizing a
    # parity vector of weight sixteen must occupy its directions once each.
    fixed_in_v_support = True
    auxiliary_parity_weight = v_weight - 1
    auxiliary_occurrences = M
    auxiliaries_distinct = auxiliary_parity_weight == auxiliary_occurrences

    if fixed_edge_type == "hard":
        positive_origin_edges = hard_v_support - 1
        negative_origin_edges = opposite_v_support
        fixed_quota = 14
    else:
        positive_origin_edges = hard_v_support
        negative_origin_edges = opposite_v_support - 1
        fixed_quota = 15

    origin_degree = positive_origin_edges + negative_origin_edges
    raw_signed_origin_degree = positive_origin_edges - negative_origin_edges
    opposite_minus_hard_origin_degree = -raw_signed_origin_degree
    cancellation_units = (RAW_HALF_EDGES - USED_NONFIXED_ORBITS) // 2

    geometry = origin_line_injectivity_certificate(P)
    origin_cancellation_excluded = bool(
        auxiliaries_distinct
        and geometry["equal_origin_orbits_force_equal_auxiliary_direction"]
    )
    fixed_edge_avoids_origin = P % 2 == 1
    proved = bool(
        T_TOP == 177
        and GRAPH_EDGES == 4 * P + 2 * T_TOP + 1
        and RAW_HALF_EDGES == 480
        and USED_NONFIXED_ORBITS == 478
        and cancellation_units == 1
        and v_weight == 17
        and fixed_in_v_support
        and auxiliary_parity_weight == 16
        and auxiliaries_distinct
        and origin_cancellation_excluded
        and fixed_edge_avoids_origin
        and origin_degree == 16
        and raw_signed_origin_degree in (10, 12)
    )
    if not proved:
        raise ArithmeticError("the p31 top common-origin ledger changed")

    return {
        "p": P,
        "t": T_TOP,
        "localized_halves": M,
        "raw_half_edge_occurrences": RAW_HALF_EDGES,
        "graph_edges": GRAPH_EDGES,
        "used_nonfixed_orbits": USED_NONFIXED_ORBITS,
        "fixed_edges": 1,
        "unused_doubled_orbits": 0,
        "cancellation_units": cancellation_units,
        "hard_quota_multiset": {"14": 14, "15": 2},
        "opposite_quota_multiset": {"15": 3, "16": 13},
        "v_support_by_type": {
            "hard": hard_v_support,
            "opposite": opposite_v_support,
        },
        "v_weight": v_weight,
        "fixed_edge_type": fixed_edge_type,
        "fixed_edge_quota": fixed_quota,
        "fixed_direction_lies_in_v_support": fixed_in_v_support,
        "auxiliary_parity_weight": auxiliary_parity_weight,
        "auxiliary_direction_occurrences": auxiliary_occurrences,
        "auxiliary_directions_are_distinct": auxiliaries_distinct,
        "auxiliary_sign_counts": {
            "hard": positive_origin_edges,
            "opposite": negative_origin_edges,
        },
        "origin_cancellation_excluded": origin_cancellation_excluded,
        "sole_cancellation_is_nonorigin": origin_cancellation_excluded,
        "fixed_antipodal_edge_avoids_origin": fixed_edge_avoids_origin,
        "origin_unsigned_degree": origin_degree,
        "origin_raw_signed_degree_hard_minus_opposite": raw_signed_origin_degree,
        "origin_signed_degree_opposite_minus_hard": opposite_minus_hard_origin_degree,
        "scope": (
            "necessary condition for the all-active localized-Mobius p31 "
            "top endpoint; not a common-graph construction"
        ),
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    hard = p31_top_origin_endpoint("hard")
    opposite = p31_top_origin_endpoint("opposite")
    return {
        "title": "p31 top localized-Mobius common-origin endpoint",
        "geometry": origin_line_injectivity_certificate(P),
        "fixed_edge_cases": {
            "hard": hard,
            "opposite": opposite,
        },
        "origin_unsigned_degree": 16,
        "origin_raw_signed_degree_cases": [10, 12],
        "origin_cancellation_possible": False,
        "residual_ii_closed": False,
        "proved": hard["proved"] and opposite["proved"],
    }


def main() -> dict[str, object]:
    result = theorem_record()
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
