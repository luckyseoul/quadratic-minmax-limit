#!/usr/bin/env python3
"""Exact symmetric ledger for one direction-localized Mobius trade.

This module is deliberately not a finite configuration census.  It reuses
the all-prime localized trade from ``e1_gmin_m4_inversion_antisymmetric_radon``
and records closed formulas for one half, its Paley signs, its forced
inversion-symmetric source chain, and an exact two-trade origin cancellation.
The coupled symmetric Boolean lift remains open.
"""
from __future__ import annotations

from collections import defaultdict

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    Point,
    _edge,
    _functional_value,
    _negative_edge,
    _point_from_coordinates,
    edge_radon_image,
    hard_star_chain,
    localized_star_trade,
    projective_functionals,
)
from e1_gmin_m4_prop15721 import is_prime


def _check_paley_prime(p: int, minimum: int = 3) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < minimum
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError(f"need a prime p=3 mod 4 with p>={minimum}")


def _legendre(p: int, value: int) -> int:
    value %= p
    if value == 0:
        return 0
    result = pow(value, (p - 1) // 2, p)
    if result == 1:
        return 1
    if result == p - 1:
        return -1
    raise ArithmeticError("Euler's criterion returned a nonsign")


def _add_functionals(
    p: int,
    first: Functional,
    first_scale: int,
    second: Functional,
    second_scale: int,
) -> Functional:
    return (
        (first_scale * first[0] + second_scale * second[0]) % p,
        (first_scale * first[1] + second_scale * second[1]) % p,
    )


def _relative_coefficients(
    p: int,
    direction: Functional,
    auxiliary: Functional,
    row: Functional,
) -> tuple[int, int]:
    """Return the unique ``(a,b)`` with ``row=a*direction+b*auxiliary``."""
    determinant = (
        direction[0] * auxiliary[1]
        - direction[1] * auxiliary[0]
    ) % p
    if determinant == 0:
        raise ValueError("direction and auxiliary must be independent")
    inverse = pow(determinant, -1, p)
    a = (row[0] * auxiliary[1] - row[1] * auxiliary[0]) * inverse % p
    b = (direction[0] * row[1] - direction[1] * row[0]) * inverse % p
    return a, b


def paley_direction_sign(p: int, functional: Functional) -> int:
    """Return ``epsilon_L=eta(a^2+b^2)`` for ``L=(a,b)``."""
    _check_paley_prime(p)
    a, b = functional
    if a % p == 0 and b % p == 0:
        raise ValueError("a projective functional must be nonzero")
    sign = _legendre(p, a * a + b * b)
    if sign == 0:
        raise ArithmeticError("x^2+y^2 must be anisotropic for p=3 mod 4")
    return sign


def paley_edge_sign(p: int, edge: Edge) -> int:
    """Return the Paley column sign of a nonloop edge."""
    _check_paley_prime(p)
    dx = (edge[0][0] - edge[1][0]) % p
    dy = (edge[0][1] - edge[1][1]) % p
    if dx == 0 and dy == 0:
        raise ValueError("a source edge cannot be a loop")
    sign = _legendre(p, dx * dx + dy * dy)
    if sign == 0:
        raise ArithmeticError("the anisotropic norm vanished")
    return sign


def mobius_parameter_edges(
    p: int,
    direction: Functional,
    auxiliary: Functional,
    center: int,
) -> dict[int, Edge]:
    """Return ``t -> {u_t,v_t}`` for ``t in F_p minus {-1}``."""
    _check_paley_prime(p)
    center %= p
    if center == 0:
        raise ValueError("the nontrivial Mobius half needs a nonzero center")
    out: dict[int, Edge] = {}
    for parameter in range(p):
        if parameter == p - 1:
            continue
        fraction = parameter * pow(parameter + 1, -1, p) % p
        first = _point_from_coordinates(
            p,
            direction,
            auxiliary,
            center,
            center * fraction % p,
        )
        second = _point_from_coordinates(
            p,
            direction,
            auxiliary,
            center * parameter % p,
            center * parameter % p,
        )
        out[parameter] = _edge(first, second)
    trade = localized_star_trade(p, direction, auxiliary, center)
    if set(out.values()) != {edge for edge, value in trade.items() if value == -1}:
        raise ArithmeticError("the parameterized half disagrees with the trade")
    return out


def _parallel_formula(
    p: int,
    direction: Functional,
    auxiliary: Functional,
    row: Functional,
) -> int:
    a, b = _relative_coefficients(p, direction, auxiliary, row)
    if b == 0:  # the chosen row L
        return 1
    if a == 0:  # the auxiliary row M
        return 1
    multiplier = b * pow(a, -1, p) % p
    if multiplier == p - 1:  # L-M
        return 0
    return 1 + _legendre(p, 1 + multiplier)


def _row_total(
    image: dict[tuple[object, ...], int], direction_index: int
) -> int:
    return sum(
        value for key, value in image.items()
        if int(key[1]) == direction_index
    )


def _normalized_image(
    p: int, signed_source: dict[Edge, int]
) -> dict[tuple[object, ...], int]:
    ordinary = edge_radon_image(p, signed_source)
    directions = projective_functionals(p)
    return {
        key: paley_direction_sign(p, directions[int(key[1])]) * value
        for key, value in ordinary.items()
        if value
    }


def _centrally_symmetric_rows(
    p: int,
    image: dict[tuple[object, ...], int],
    excluded_direction_index: int | None = None,
) -> bool:
    for key, value in image.items():
        if key[0] != "K" or int(key[1]) == excluded_direction_index:
            continue
        direction_index = int(key[1])
        left, right = int(key[2]), int(key[3])
        negative = tuple(sorted(((-left) % p, (-right) % p)))
        if image.get(("K", direction_index, *negative), 0) != value:
            return False
    return True


def ordinary_half_certificate(
    p: int,
    direction_index: int = 0,
    auxiliary: Functional = (0, 1),
    center: int = 1,
) -> dict[str, object]:
    """Replay the exact ``P_N`` formula and centrality of ``R_N(E)``."""
    _check_paley_prime(p)
    directions = projective_functionals(p)
    if not 0 <= direction_index < len(directions):
        raise ValueError("direction index is out of range")
    direction = directions[direction_index]
    parameter_edges = mobius_parameter_edges(
        p, direction, auxiliary, center
    )
    image = edge_radon_image(
        p, {edge: 1 for edge in parameter_edges.values()}
    )
    actual_parallel = tuple(
        image.get(("P", index), 0) for index in range(p + 1)
    )
    formula_parallel = tuple(
        _parallel_formula(p, direction, auxiliary, row)
        for row in directions
    )
    chosen_transverse = {
        key[2:]: value
        for key, value in image.items()
        if key[0] == "K" and int(key[1]) == direction_index
    }
    expected_chosen = {
        tuple(sorted((center % p, center * parameter % p))): 1
        for parameter in range(p)
        if parameter not in (1, p - 1)
    }
    proved = bool(
        actual_parallel == formula_parallel
        and sum(actual_parallel) == p - 1
        and actual_parallel[direction_index] == 1
        and chosen_transverse == expected_chosen
        and _centrally_symmetric_rows(
            p, image, excluded_direction_index=direction_index
        )
        and all(_row_total(image, index) == p - 1 for index in range(p + 1))
    )
    if not proved:
        raise ArithmeticError("the ordinary Mobius-half formulas changed")
    return {
        "p": p,
        "direction": list(direction),
        "auxiliary": list(auxiliary),
        "center": center % p,
        "edge_count": p - 1,
        "parallel_counts": list(actual_parallel),
        "parallel_count_formula": {
            "L": 1,
            "M": 1,
            "L-M": 0,
            "L+mM_other": "1+eta(1+m)",
        },
        "sum_parallel_counts": sum(actual_parallel),
        "chosen_K_row": "S_j minus {j,-j}",
        "all_other_K_rows_centrally_symmetric": True,
        "proved": proved,
    }


def paley_tau_certificate(
    p: int,
    direction_index: int = 0,
    auxiliary: Functional = (0, 1),
    center: int = 1,
) -> dict[str, object]:
    r"""Certify ``tau_t=eta(Q(e1-t^2(e1+e2)))`` and its trace."""
    _check_paley_prime(p)
    directions = projective_functionals(p)
    direction = directions[direction_index]
    edges = mobius_parameter_edges(p, direction, auxiliary, center)
    e1 = _point_from_coordinates(p, direction, auxiliary, 1, 0)
    e2 = _point_from_coordinates(p, direction, auxiliary, 0, 1)

    direct: dict[int, int] = {}
    formula: dict[int, int] = {}
    for parameter, edge in edges.items():
        direct[parameter] = paley_edge_sign(p, edge)
        square = parameter * parameter % p
        vector = (
            (e1[0] - square * (e1[0] + e2[0])) % p,
            (e1[1] - square * (e1[1] + e2[1])) % p,
        )
        formula[parameter] = _legendre(
            p, vector[0] * vector[0] + vector[1] * vector[1]
        )

    complete_trace = 0
    for parameter in range(p):
        square = parameter * parameter % p
        vector = (
            (e1[0] - square * (e1[0] + e2[0])) % p,
            (e1[1] - square * (e1[1] + e2[1])) % p,
        )
        complete_trace += _legendre(
            p, vector[0] * vector[0] + vector[1] * vector[1]
        )
    domain_sum = sum(direct.values())
    epsilon_l = paley_direction_sign(p, direction)
    even_checks = [
        direct[parameter] == direct[-parameter % p]
        for parameter in range(p)
        if parameter not in (1, p - 1)
    ]
    proved = bool(
        direct == formula
        and direct[1] == epsilon_l
        and domain_sum == complete_trace - epsilon_l
        and all(even_checks)
    )
    if not proved:
        raise ArithmeticError("the Paley quartic trace identity changed")
    return {
        "p": p,
        "direction": list(direction),
        "auxiliary": list(auxiliary),
        "dual_basis_e1_e2": [list(e1), list(e2)],
        "tau_by_parameter": {str(key): value for key, value in direct.items()},
        "formula": "eta(Q(e1-t^2(e1+e2)))",
        "tau_t_equals_tau_minus_t_when_both_parameters_are_in_domain": True,
        "tau_1": direct[1],
        "epsilon_L": epsilon_l,
        "complete_quartic_character_sum": complete_trace,
        "domain_sum_S": domain_sum,
        "trace_identity": "S=sum_Fp eta(Q(e1-t^2(e1+e2)))-epsilon_L",
        "proved": proved,
    }


def forced_symmetric_certificate(
    p: int,
    direction_index: int = 0,
    auxiliary: Functional = (0, 1),
    center: int = 1,
) -> dict[str, object]:
    """Return the normalized forced pair-total and selected-half ledgers."""
    _check_paley_prime(p)
    directions = projective_functionals(p)
    direction = directions[direction_index]
    if paley_direction_sign(p, direction) != 1:
        raise ValueError("the selected target direction must be hard")
    edges = mobius_parameter_edges(p, direction, auxiliary, center)
    tau = {parameter: paley_edge_sign(p, edge) for parameter, edge in edges.items()}

    symmetric_source: dict[Edge, int] = {}
    selected_source: dict[Edge, int] = {}
    for parameter, edge in edges.items():
        negative = _negative_edge(p, edge)
        sign = tau[parameter]
        symmetric_source[edge] = sign
        symmetric_source[negative] = sign
        chosen = edge if sign == -1 else negative
        selected_source[chosen] = sign

    symmetric_image = _normalized_image(p, symmetric_source)
    selected_image = _normalized_image(p, selected_source)
    parallel_formula = tuple(
        _parallel_formula(p, direction, auxiliary, row)
        for row in directions
    )
    domain_sum = sum(tau.values())
    symmetric_parallel = tuple(
        symmetric_image.get(("P", index), 0) for index in range(p + 1)
    )
    selected_parallel = tuple(
        selected_image.get(("P", index), 0) for index in range(p + 1)
    )
    expected_antisymmetric = {
        ("K", direction_index, *cell): value
        for cell, value in hard_star_chain(p, center).items()
    }
    actual_antisymmetric: defaultdict[tuple[object, ...], int] = defaultdict(int)
    for key, value in selected_image.items():
        if key[0] != "K":
            continue
        actual_antisymmetric[key] += value
        left, right = int(key[2]), int(key[3])
        negative = tuple(sorted(((-left) % p, (-right) % p)))
        actual_antisymmetric[("K", int(key[1]), *negative)] -= value
    actual_antisymmetric = defaultdict(
        int,
        {key: value for key, value in actual_antisymmetric.items() if value},
    )
    proved = bool(
        len(symmetric_source) == 2 * (p - 1)
        and len(selected_source) == p - 1
        and symmetric_parallel == tuple(2 * value for value in parallel_formula)
        and selected_parallel == parallel_formula
        and _centrally_symmetric_rows(p, symmetric_image)
        and _centrally_symmetric_rows(
            p, selected_image, excluded_direction_index=direction_index
        )
        and dict(actual_antisymmetric) == expected_antisymmetric
        and all(
            _row_total(symmetric_image, index)
            == 2 * paley_direction_sign(p, row) * domain_sum
            for index, row in enumerate(directions)
        )
        and all(
            _row_total(selected_image, index)
            == paley_direction_sign(p, row) * domain_sum
            for index, row in enumerate(directions)
        )
    )
    if not proved:
        raise ArithmeticError("the forced symmetric Mobius-half ledger changed")
    return {
        "p": p,
        "used_inversion_orbits": p - 1,
        "selected_graph_edges": p - 1,
        "domain_tau_sum_S": domain_sum,
        "selected_parallel_counts": list(selected_parallel),
        "forced_pair_total_parallel_counts": list(symmetric_parallel),
        "selected_row_total_formula": "epsilon_N*S",
        "forced_symmetric_row_total_formula": "2*epsilon_N*S",
        "forced_symmetric_K_rows_central": True,
        "selected_nonchosen_K_rows_central": True,
        "selected_antisymmetric_image": "A_j in L and zero elsewhere",
        "coupled_symmetric_target_realized": False,
        "proved": proved,
    }


def two_trade_origin_cancellation_certificate(
    p: int,
    first_direction_index: int,
    second_direction_index: int,
    first_center: int,
    second_center: int,
) -> dict[str, object]:
    """Construct two localized trades sharing exactly one cancelled orbit."""
    _check_paley_prime(p)
    directions = projective_functionals(p)
    if first_direction_index == second_direction_index:
        raise ValueError("the target directions must be distinct")
    first_direction = directions[first_direction_index]
    second_direction = directions[second_direction_index]
    first_center %= p
    second_center %= p
    if first_center == 0 or second_center == 0:
        raise ValueError("both star centers must be nonzero")

    first_auxiliary = _add_functionals(
        p,
        first_direction,
        1,
        second_direction,
        first_center * pow(second_center, -1, p),
    )
    second_auxiliary = _add_functionals(
        p,
        first_direction,
        second_center * pow(first_center, -1, p),
        second_direction,
        1,
    )
    first = localized_star_trade(
        p, first_direction, first_auxiliary, first_center
    )
    second = localized_star_trade(
        p, second_direction, second_auxiliary, second_center
    )
    common_edges = set(first) & set(second)
    same_sign = {edge for edge in common_edges if first[edge] == second[edge]}
    opposite_sign = common_edges - same_sign
    source = {
        edge: first.get(edge, 0) + second.get(edge, 0)
        for edge in set(first) | set(second)
    }
    source = {edge: value for edge, value in source.items() if value}
    expected: dict[tuple[object, ...], int] = {}
    for index, center in (
        (first_direction_index, first_center),
        (second_direction_index, second_center),
    ):
        expected.update({
            ("K", index, *cell): value
            for cell, value in hard_star_chain(p, center).items()
        })
    proved = bool(
        len(common_edges) == 2
        and not same_sign
        and len(opposite_sign) == 2
        and len(source) == 4 * (p - 1) - 4
        and set(source.values()) == {-1, 1}
        and edge_radon_image(p, source) == expected
    )
    if not proved:
        raise ArithmeticError("the exact origin-orbit cancellation changed")
    return {
        "p": p,
        "first_direction": list(first_direction),
        "second_direction": list(second_direction),
        "first_center": first_center,
        "second_center": second_center,
        "first_auxiliary": list(first_auxiliary),
        "second_auxiliary": list(second_auxiliary),
        "shared_actual_edges": len(common_edges),
        "shared_inversion_orbits": 1,
        "same_sign_shared_orbits": 0,
        "opposite_sign_shared_orbits": 1,
        "nonzero_orbits_after_cancellation": 2 * (p - 1) - 2,
        "source_stays_ternary": True,
        "both_direction_targets_stay_exact": True,
        "proved": proved,
    }


def branch_c_capacity_ledger(p: int) -> dict[str, object]:
    """Record the exact all-active support gap across the branch-C ray.

    The deliberately disjoint family has ``N=m(p-1)`` used inversion
    orbits.  If opposite-sign overlaps reduce its actual ternary support to
    ``|U|=N-2*kappa``, extension to an ``H``-edge graph first requires
    ``|U|<=|H|``.  This pins the minimum cancellation count at every ``t``;
    it does not assert that a family with that many compatible overlaps
    exists.
    """
    _check_paley_prime(p, minimum=31)
    r = (p - 3) // 4
    hard_count = (p + 1) // 2
    lower_t = 2 * r * r - 4 * r - 2
    upper_t = 4 * r * r - 2 * r - 5
    lower_edges = 4 * p + 2 * lower_t + 1
    upper_edges = 4 * p + 2 * upper_t + 1
    disjoint_used = hard_count * (p - 1)
    one_pair_used = disjoint_used - 2
    lower_gap = disjoint_used - lower_edges
    upper_gap = disjoint_used - upper_edges
    lower_required_cancellations = (lower_gap + 1) // 2
    upper_required_cancellations = (upper_gap + 1) // 2
    proved = bool(
        lower_edges == 4 * r * r + 8 * r + 9
        and upper_edges == 8 * r * r + 12 * r + 3
        and disjoint_used == upper_edges + 1
        and one_pair_used == upper_edges - 1
        and lower_gap == 2 * (upper_t - lower_t) + 1
        and upper_gap == 1
        and lower_required_cancellations == upper_t - lower_t + 1
        and upper_required_cancellations == 1
    )
    if not proved:
        raise ArithmeticError("the branch-C capacity identity changed")
    return {
        "p": p,
        "r": r,
        "hard_direction_count": hard_count,
        "t_interval": [lower_t, upper_t],
        "H_edge_count_interval": [lower_edges, upper_edges],
        "all_centers_nonzero_disjoint_trade_edges": disjoint_used,
        "disjoint_support_excess_over_H_interval": [lower_gap, upper_gap],
        "minimum_opposite_sign_cancellations_interval": [
            lower_required_cancellations,
            upper_required_cancellations,
        ],
        "minimum_cancellations_at_t": "t_max-t+1",
        "support_after_kappa_cancellations": "m(p-1)-2*kappa",
        "remaining_edge_capacity_after_kappa": (
            "2*(kappa-(t_max-t))-1"
        ),
        "disjoint_lift_extendable_anywhere_on_ray": False,
        "reason_disjoint_lift_is_not_extendable": (
            "m(p-1)-|H|=2(t_max-t)+1>0"
        ),
        "minimum_cancellation_forces": (
            "one fixed antipodal edge and zero unused double orbits"
        ),
        "forced_fixed_edge_weight_parity_after_any_cancellations": "odd",
        "disjoint_minus_H_upper": 1,
        "one_cancelled_pair_trade_edges": one_pair_used,
        "one_pair_minus_H_upper": -1,
        "universal_support_lower_bound_proved": False,
        "required_multi_overlap_family_constructed": False,
        "capacity_contradiction_proved": False,
        "reason_not_an_obstruction": (
            "two arbitrary nonzero target stars can cancel one origin orbit"
        ),
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Return the exact new reduction and explicit OPEN flags."""
    half = ordinary_half_certificate(p)
    tau = paley_tau_certificate(p)
    forced = forced_symmetric_certificate(p)
    overlap = two_trade_origin_cancellation_certificate(p, 0, 1, 1, 2)
    capacity = branch_c_capacity_ledger(p)
    proved = all(
        record["proved"] for record in (half, tau, forced, overlap, capacity)
    )
    return {
        "title": "Mobius-half symmetric image and overlap gate",
        "method": "closed formulas and exact identity replay; no census",
        "ordinary_half": half,
        "paley_tau": tau,
        "forced_symmetric_half": forced,
        "two_trade_overlap": overlap,
        "branch_C_capacity": capacity,
        "remaining_exact_gate": (
            "realize the prescribed hard compact and opposite AE+compact "
            "central rows in the restricted unused-orbit Boolean box"
        ),
        "signed_boolean_lift_proved": False,
        "residual_ii_closed": False,
        "E1_closed": False,
        "L_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
