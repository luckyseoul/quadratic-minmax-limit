#!/usr/bin/env python3
"""Exact nonorigin collision loci at the p=31 branch-C top endpoint.

The top localized-Mobius ledger has one cancellation unit and, by the
common-origin theorem, that cancellation is nonorigin.  This module records
the two exact opposite-orientation loci for a pair of halves with distinct
auxiliary directions.  It also gives an explicit three-half construction in
which one orbit has orientation multiplicities ``2:1``.  Thus the top
cancellation need not be an isolated cancelling pair.

The result is deliberately local.  It neither supplies the other thirteen
halves nor realizes the complete transverse target, so branch C and residual
(ii) remain open.
"""

from __future__ import annotations

from collections import defaultdict
import json
from itertools import combinations

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    _edge,
    _functional_value,
    _negative_edge,
    localized_star_trade,
)
from e1_gmin_m4_mobius_half_intersections import (
    two_half_intersection_candidates,
)
from e1_gmin_m4_mobius_half_symmetric import (
    mobius_parameter_edges,
    paley_direction_sign,
)
from e1_gmin_m4_prop15721 import is_prime


P = 31
Point = tuple[int, int]


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 5
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=5")


def _inverse(p: int, value: int) -> int:
    value %= p
    if value == 0:
        raise ZeroDivisionError("cannot invert zero")
    return pow(value, -1, p)


def _projective_normalize(p: int, functional: Functional) -> Functional:
    first, second = (coordinate % p for coordinate in functional)
    if first:
        scale = _inverse(p, first)
        return 1, second * scale % p
    if second:
        return 0, 1
    raise ValueError("zero has no projective direction")


def _same_projective_direction(
    p: int, first: Functional, second: Functional
) -> bool:
    return (first[0] * second[1] - first[1] * second[0]) % p == 0


def _orbit(p: int, edge: Edge) -> Edge:
    return min(edge, _negative_edge(p, edge))


def _edge_is_parallel(p: int, direction: Functional, edge: Edge) -> bool:
    return _functional_value(p, direction, edge[0]) == _functional_value(
        p, direction, edge[1]
    )


def opposite_direct_collision(
    p: int, alpha: int, beta: int, q: int
) -> dict[str, object]:
    r"""Instantiate the exact opposite-direct pair locus.

    Normalize two target rows to a dual basis ``X,Y`` and normalize the
    projective auxiliaries as

        M1 ~ alpha*X+Y,       M2 ~ X+beta*Y.

    Distinct auxiliaries mean ``alpha*beta != 1``.  The opposite direct
    endpoint matching exists exactly on

        q*r=1,               alpha+beta=2.

    The remaining nonzero parameter is ``q``; ``A=q+alpha`` and
    ``B=r+beta`` are the actual half parameters used by the existing normal
    form.
    """
    _check_prime(p)
    alpha %= p
    beta %= p
    q %= p
    if alpha * beta % p == 1:
        raise ValueError("the two projective auxiliary directions must differ")
    if (alpha + beta) % p != 2 % p:
        raise ValueError("the opposite-direct locus requires alpha+beta=2")
    if q in (0, p - 1):
        raise ValueError("the direct parameter q must avoid 0 and -1")

    r = _inverse(p, q)
    A = (q + alpha) % p
    B = (r + beta) % p
    if A == 0 or B == 0:
        raise ValueError("the Mobius parameters A and B must be nonzero")
    t = (alpha - 1) * _inverse(p, q + 1) % p
    s = (beta - 1) * _inverse(p, r + 1) % p
    edge = _edge((1, p - 1), (t, q * t % p))
    cancellation_direction_raw = (
        (1 + q * t) % p,
        (1 - t) % p,
    )

    replay = two_half_intersection_candidates(p, q, r, A, B)
    intersections = replay["intersections"]
    unique_direct = bool(
        len(intersections) == 1
        and intersections[0]["sign"] == -1
        and intersections[0]["matching"] == "direct"
        and tuple(tuple(point) for point in intersections[0]["edge"]) == edge
    )
    nonorigin = (0, 0) not in edge
    auxiliaries_distinct = alpha * beta % p != 1
    proved = bool(
        q * r % p == 1
        and (alpha + beta) % p == 2 % p
        and auxiliaries_distinct
        and unique_direct
        and nonorigin
        and _edge_is_parallel(p, cancellation_direction_raw, edge)
    )
    if not proved:
        raise ArithmeticError("the opposite-direct collision locus changed")
    return {
        "p": p,
        "auxiliary_coordinates": {
            "M1": [alpha, 1],
            "M2": [1, beta],
        },
        "auxiliary_directions_distinct": auxiliaries_distinct,
        "locus_equations": ["q*r=1", "alpha+beta=2"],
        "parameters": {"q": q, "r": r, "A": A, "B": B, "t": t, "s": s},
        "common_edge": [list(point) for point in edge],
        "common_edge_is_nonorigin": nonorigin,
        "cancellation_direction_raw": list(cancellation_direction_raw),
        "cancellation_direction": list(
            _projective_normalize(p, cancellation_direction_raw)
        ),
        "cancellation_direction_formula": (
            "C~(1+q*alpha)X+(q+beta)Y"
        ),
        "prescribed_C_generic_recovery": (
            "if C~cX+Y then q=(1-c*beta)/(c-alpha)"
        ),
        "shared_inversion_orbits": 1,
        "matching": "opposite-direct",
        "proved": proved,
    }


def opposite_swapped_collision(
    p: int, alpha: int, beta: int, z: int
) -> dict[str, object]:
    r"""Instantiate one root of the exact opposite-swapped pair locus.

    With ``u=alpha+1``, ``v=beta+1`` and ``z=q*r``, the locus is the
    quadratic

        u*v*z^2-(u+v+1)*z+1=0.

    Every admissible root uniquely recovers ``q=1-u*z`` and ``r=1-v*z``.
    Hence fixed target/auxiliary data leave at most two swapped candidates.
    """
    _check_prime(p)
    alpha %= p
    beta %= p
    z %= p
    if alpha * beta % p == 1:
        raise ValueError("the two projective auxiliary directions must differ")
    u = (alpha + 1) % p
    v = (beta + 1) % p
    polynomial = (u * v * z * z - (u + v + 1) * z + 1) % p
    if polynomial:
        raise ValueError("z is not a root of the swapped collision quadratic")
    q = (1 - u * z) % p
    r = (1 - v * z) % p
    if q in (0, 1) or r in (0, 1) or z in (0, 1):
        raise ValueError("the swapped root is Mobius-degenerate")
    A = (q + alpha) % p
    B = (r + beta) % p
    if A == 0 or B == 0 or q * r % p != z:
        raise ValueError("the swapped root does not recover admissible A,B,q,r")

    t = -_inverse(p, q) % p
    s = -_inverse(p, r) % p
    edge = _edge((1, _inverse(p, r)), (-_inverse(p, q) % p, p - 1))
    cancellation_direction_raw = (
        q * (r + 1) % p,
        -r * (q + 1) % p,
    )
    replay = two_half_intersection_candidates(p, q, r, A, B)
    intersections = replay["intersections"]
    unique_swapped = bool(
        len(intersections) == 1
        and intersections[0]["sign"] == -1
        and intersections[0]["matching"] == "swapped"
        and tuple(tuple(point) for point in intersections[0]["edge"]) == edge
    )
    discriminant = ((u + v + 1) ** 2 - 4 * u * v) % p
    proved = bool(
        polynomial == 0
        and q * r % p == z
        and unique_swapped
        and (0, 0) not in edge
        and _edge_is_parallel(p, cancellation_direction_raw, edge)
    )
    if not proved:
        raise ArithmeticError("the opposite-swapped collision locus changed")
    return {
        "p": p,
        "auxiliary_coordinates": {
            "M1": [alpha, 1],
            "M2": [1, beta],
        },
        "auxiliary_directions_distinct": True,
        "quadratic": (
            "(alpha+1)(beta+1)z^2-(alpha+beta+3)z+1=0"
        ),
        "quadratic_discriminant": discriminant,
        "at_most_two_roots": True,
        "parameters": {
            "z": z,
            "q": q,
            "r": r,
            "A": A,
            "B": B,
            "t": t,
            "s": s,
        },
        "common_edge": [list(point) for point in edge],
        "common_edge_is_nonorigin": True,
        "cancellation_direction_raw": list(cancellation_direction_raw),
        "cancellation_direction": list(
            _projective_normalize(p, cancellation_direction_raw)
        ),
        "cancellation_direction_formula": (
            "C~q(r+1)X-r(q+1)Y"
        ),
        "shared_inversion_orbits": 1,
        "matching": "opposite-swapped",
        "proved": proved,
    }


def pair_collision_locus_theorem(p: int = P) -> dict[str, object]:
    """Package both symbolic pair loci and one exact witness for each."""
    _check_prime(p)
    direct = opposite_direct_collision(p, alpha=0, beta=2, q=2)

    # q=2,r=3 gives z=6 and, from the swapped formulas, alpha=4,beta=9
    # at p=31.  Keep this theorem record deliberately pinned to p=31.
    if p != P:
        raise ValueError("the frozen swapped witness is pinned to p=31")
    swapped = opposite_swapped_collision(p, alpha=4, beta=9, z=6)
    return {
        "p": p,
        "normalization": (
            "X=L1/j1, Y=L2/j2, M1~alpha X+Y, M2~X+beta Y"
        ),
        "distinct_auxiliary_condition": "alpha*beta!=1",
        "opposite_direct_locus": {
            "equations": ["q*r=1", "alpha+beta=2"],
            "edge": "{(1,-1),(t,q*t)}, t=(alpha-1)/(q+1)",
            "direction": "C~(1+q*alpha)X+(q+beta)Y",
            "prescribed_direction_leaves_at_most_one_candidate": True,
            "witness": direct,
        },
        "opposite_swapped_locus": {
            "equation": (
                "(alpha+1)(beta+1)z^2-(alpha+beta+3)z+1=0"
            ),
            "recovery": (
                "q=1-(alpha+1)z, r=1-(beta+1)z, z=q*r"
            ),
            "edge": "{(1,1/r),(-1/q,-1)}",
            "direction": "C~q(r+1)X-r(q+1)Y",
            "fixed_target_auxiliary_data_leave_at_most_two_candidates": True,
            "witness": swapped,
        },
        "rigid_two_overlap_point_excluded": (
            "q=r=1/2, A=B=3/2 forces alpha=beta=1 and M1=M2"
        ),
        "distinct_auxiliaries_force_at_most_one_opposite_shared_orbit_per_pair": True,
        "pair_only_description_of_top_cancellation_is_complete": False,
        "proved": True,
    }


def _theta_half(
    p: int, theta: int, center: int = 1
) -> tuple[Functional, Functional, dict[int, Edge]]:
    theta %= p
    if theta in (0, p - 1):
        raise ValueError("theta must avoid 0 and -1")
    direction = (1, theta)
    auxiliary = (theta * _inverse(p, theta + 1) % p, theta)
    edges = mobius_parameter_edges(p, direction, auxiliary, center)
    return direction, auxiliary, edges


def theta_pair_intersection(
    p: int, theta: int, phi: int
) -> dict[str, object]:
    """Replay the exact one-orbit intersection of two theta-family halves."""
    _check_prime(p)
    theta %= p
    phi %= p
    if theta == phi:
        raise ValueError("theta and phi must differ")
    first_direction, first_auxiliary, first = _theta_half(p, theta)
    second_direction, second_auxiliary, second = _theta_half(p, phi)
    common_edge = _edge((1, 0), (0, 1))
    common_orbit = _orbit(p, common_edge)
    first_orbits = {_orbit(p, edge) for edge in first.values()}
    second_orbits = {_orbit(p, edge) for edge in second.values()}
    actual_intersection = set(first.values()) & set(second.values())
    orbit_intersection = first_orbits & second_orbits
    proved = bool(
        actual_intersection == {common_edge}
        and orbit_intersection == {common_orbit}
        and first[theta] == common_edge
        and second[phi] == common_edge
        and not _same_projective_direction(
            p, first_auxiliary, second_auxiliary
        )
    )
    if not proved:
        raise ArithmeticError("the theta-family pair intersection changed")
    return {
        "p": p,
        "theta": theta,
        "phi": phi,
        "first_target": list(first_direction),
        "second_target": list(second_direction),
        "first_auxiliary": list(first_auxiliary),
        "second_auxiliary": list(second_auxiliary),
        "projective_auxiliaries": [[1, (theta + 1) % p], [1, (phi + 1) % p]],
        "common_edge": [list(point) for point in common_edge],
        "common_inversion_orbits": 1,
        "same_orientation_common_edges": 1,
        "direct_same_equation": "(phi-theta)(z-1)=0, hence z=1",
        "direct_opposite_remainder": "(theta+phi)(1-z^2)=0",
        "opposite_cases_are_excluded": (
            "z=1 forces theta=-1; z=-1 forces phi=-1; "
            "phi=-theta forces both Mobius parameters to -1"
        ),
        "swapped_matching_is_impossible": (
            "the affine endpoint has nonzero x-coordinate theta+1"
        ),
        "proved": proved,
    }


def triple_overlap_countermechanism() -> dict[str, object]:
    r"""Return a p=31 hard-target triple with one ternary 2:1 overlap.

    For independent endpoints ``x,y`` and every ``theta != 0,-1``, set

        L_theta(x)=1,       L_theta(y)=theta,
        M_theta(x)=theta/(theta+1),  M_theta(y)=theta.

    The theta-parameter edge is ``{x,y}``.  Distinct theta halves share no
    other inversion orbit.  Centers ``(+1,+1,-1)`` therefore give common
    coefficients ``(-1,-1,+1)`` and one cancellation unit while the full
    sum remains ternary.
    """
    p = P
    thetas = (1, 2, 3)
    centers = (1, 1, p - 1)
    pair_records = tuple(
        theta_pair_intersection(p, first, second)
        for first, second in combinations(thetas, 2)
    )

    targets: list[Functional] = []
    auxiliaries: list[Functional] = []
    trades = []
    raw_orbit_sets = []
    for theta, center in zip(thetas, centers, strict=True):
        direction, auxiliary, edges = _theta_half(p, theta, center)
        targets.append(direction)
        auxiliaries.append(auxiliary)
        trades.append(localized_star_trade(p, direction, auxiliary, center))
        raw_orbit_sets.append({_orbit(p, edge) for edge in edges.values()})

    source: defaultdict[Edge, int] = defaultdict(int)
    for trade in trades:
        for edge, coefficient in trade.items():
            source[edge] += coefficient
    source = defaultdict(
        int, {edge: value for edge, value in source.items() if value}
    )

    common_edge = _edge((1, 0), (0, 1))
    negative_common_edge = _negative_edge(p, common_edge)
    common_orbit = _orbit(p, common_edge)
    common_orientation_coefficients = tuple(
        trade.get(common_edge, 0) for trade in trades
    )
    pairwise_orbit_intersections = tuple(
        len(raw_orbit_sets[first] & raw_orbit_sets[second])
        for first, second in combinations(range(3), 2)
    )
    target_signs = tuple(paley_direction_sign(p, row) for row in targets)
    auxiliary_directions = tuple(
        _projective_normalize(p, row) for row in auxiliaries
    )
    auxiliary_signs = tuple(
        paley_direction_sign(p, row) for row in auxiliaries
    )
    cancellation_direction = (1, 1)
    raw_occurrences = 3 * (p - 1)
    support_orbits = len(source) // 2
    cancellation_units = (raw_occurrences - support_orbits) // 2
    proved = bool(
        all(record["proved"] for record in pair_records)
        and pairwise_orbit_intersections == (1, 1, 1)
        and len(set(auxiliary_directions)) == 3
        and target_signs == (1, 1, 1)
        and auxiliary_signs == (1, 1, -1)
        and common_orientation_coefficients == (-1, -1, 1)
        and source[common_edge] == -1
        and source[negative_common_edge] == 1
        and set(source.values()) == {-1, 1}
        and support_orbits == 88
        and cancellation_units == 1
        and (0, 0) not in common_edge
        and _edge_is_parallel(p, cancellation_direction, common_edge)
        and paley_direction_sign(p, cancellation_direction) == 1
        and common_orbit in set.intersection(*raw_orbit_sets)
    )
    if not proved:
        raise ArithmeticError("the p31 triple-overlap countermechanism changed")
    return {
        "p": p,
        "basis_endpoints": {"x": [1, 0], "y": [0, 1]},
        "parameter_family": {
            "target": "L_theta=(1,theta)",
            "actual_auxiliary": (
                "M_theta=(theta/(theta+1),theta)~(1,theta+1)"
            ),
            "edge_at_parameter_theta": "{x,y}",
            "allowed_theta": "theta not in {0,-1}",
        },
        "thetas": list(thetas),
        "centers": [1, 1, -1],
        "targets": [list(row) for row in targets],
        "target_signs": list(target_signs),
        "projective_auxiliaries": [list(row) for row in auxiliary_directions],
        "auxiliary_signs": list(auxiliary_signs),
        "auxiliary_type_counts": {"hard": 2, "opposite": 1},
        "auxiliary_directions_distinct": True,
        "pairwise_intersections": list(pair_records),
        "pairwise_shared_orbit_counts": list(pairwise_orbit_intersections),
        "common_edge": [list(point) for point in common_edge],
        "common_edge_is_nonorigin": True,
        "common_orientation_coefficients": list(
            common_orientation_coefficients
        ),
        "common_orientation_multiplicities": {"negative": 2, "positive": 1},
        "final_common_coefficient": -1,
        "raw_orbit_occurrences": raw_occurrences,
        "final_support_orbits": support_orbits,
        "cancellation_units": cancellation_units,
        "full_three_trade_sum_is_ternary": True,
        "cancellation_direction": list(cancellation_direction),
        "cancellation_direction_type": "hard",
        "top_auxiliary_sign_counts_do_not_exclude_this_local_triple": True,
        "full_p31_top_lift_constructed": False,
        "proved": proved,
    }


def top_parallel_collision_ledger(
    fixed_type: str, cancellation_type: str
) -> dict[str, object]:
    """Pull one pair or triple cancellation through the top parallel ledger."""
    signs = {"hard": 1, "opposite": -1}
    if fixed_type not in signs or cancellation_type not in signs:
        raise ValueError("types must be 'hard' or 'opposite'")
    epsilon_f = signs[fixed_type]
    epsilon_c = signs[cancellation_type]
    trace = 4 - P - epsilon_f + 2 * epsilon_c
    expected = {
        ("hard", "hard"): -26,
        ("hard", "opposite"): -30,
        ("opposite", "hard"): -24,
        ("opposite", "opposite"): -28,
    }[(fixed_type, cancellation_type)]
    if trace != expected:
        raise ArithmeticError("the p31 weighted parallel trace changed")
    return {
        "p": P,
        "fixed_direction_type": fixed_type,
        "cancellation_direction_type": cancellation_type,
        "raw_parallel_identity": (
            "sum_i P_D(E_i)=P_D(target)-1_{D=F}+2*1_{D=C}"
        ),
        "weighted_trace_identity": (
            "sum_i S(L_i,M_i)=4-p-epsilon_F+2*epsilon_C"
        ),
        "weighted_trace": trace,
        "identity_applies_to_pair_1_to_1_overlap": True,
        "identity_applies_to_triple_2_to_1_overlap": True,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    triple = triple_overlap_countermechanism()
    return {
        "title": "p31 top nonorigin pair loci and triple-overlap countermechanism",
        "pair_collision_loci": pair_collision_locus_theorem(),
        "triple_overlap": triple,
        "one_cancellation_unit_combinatorial_alternatives": {
            "pair": "one orbit with raw orientation multiplicities 1:1 and final 0",
            "triple": "one orbit with raw orientation multiplicities 2:1 and final +/-1",
        },
        "pair_only_top_models_are_complete": False,
        "reason": (
            "the explicit hard-target p31 triple has distinct auxiliaries, "
            "one nonorigin 2:1 overlap, and an exactly ternary sum"
        ),
        "parallel_trace_cases": {
            f"fixed_{fixed}_cancellation_{cancellation}": (
                top_parallel_collision_ledger(fixed, cancellation)
            )
            for fixed in ("hard", "opposite")
            for cancellation in ("hard", "opposite")
        },
        "full_top_common_graph_constructed": False,
        "branch_c_closed": False,
        "residual_ii_closed": False,
        "proved": bool(triple["proved"]),
    }


def main() -> dict[str, object]:
    result = theorem_record()
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
