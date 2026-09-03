#!/usr/bin/env python3
"""Exact symbolic ledger for the inversion-antisymmetric edge-Radon block.

This module performs no configuration census.  It records the closed-form
rank, kernel, odd-moment cokernel, hard-star row, direction-localized Mobius
trade, and ternary lift proved in
NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations

from e1_gmin_m4_prop15721 import is_prime


Cell = tuple[int, int]
Point = tuple[int, int]
Edge = tuple[Point, Point]
Functional = tuple[int, int]


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")


def antisymmetric_dimensions(p: int) -> dict[str, int | bool]:
    """Return all characteristic-zero dimensions and half-direction ranks."""
    _check_prime(p)
    h = (p - 1) // 2
    d = p + 1
    difference_classes = d * h
    source_fixed_edges = difference_classes
    source_edges = p * p * difference_classes
    source_minus = (source_edges - source_fixed_edges) // 2
    target_minus_per_direction = h * h
    target_minus = d * target_minus_per_direction
    kernel_minus = source_minus - target_minus
    half_direction_count = d // 2
    half_rank = half_direction_count * target_minus_per_direction
    odd_cokernel_rank_direct = sum(
        (degree // 2) * (p - degree)
        for degree in range(3, p - 1, 2)
    )
    odd_cokernel_rank_closed = h * (h - 1) * (h + 1) // 3
    proved = bool(
        difference_classes == (p * p - 1) // 2
        and source_minus == difference_classes * difference_classes
        and target_minus == d * h * h
        and kernel_minus == d * p * h * h
        and 2 * half_rank == target_minus
        and odd_cokernel_rank_direct == odd_cokernel_rank_closed
    )
    if not proved:
        raise ArithmeticError("antisymmetric dimension identity changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "difference_classes": difference_classes,
        "source_fixed_edges": source_fixed_edges,
        "source_minus_rank": source_minus,
        "target_minus_rank_per_direction": target_minus_per_direction,
        "target_minus_rank": target_minus,
        "kernel_minus_rank": kernel_minus,
        "hard_direction_count": half_direction_count,
        "opposite_direction_count": half_direction_count,
        "hard_restriction_rank": half_rank,
        "opposite_restriction_rank": half_rank,
        "odd_moment_cokernel_rank_direct": odd_cokernel_rank_direct,
        "odd_moment_cokernel_rank_closed": odd_cokernel_rank_closed,
        "proved": proved,
    }


def _cell(p: int, u: int, v: int) -> Cell:
    inv2 = pow(2, -1, p)
    alpha = ((u + v) * inv2) % p
    delta = ((u - v) * inv2) % p
    return alpha, (delta * delta) % p


def hard_star_difference_direct(p: int, j: int) -> dict[Cell, int]:
    """Return S_{-j}-S_j by direct signed edge collection."""
    _check_prime(p)
    j %= p
    out: defaultdict[Cell, int] = defaultdict(int)
    for other in range(p):
        if other != (-j) % p:
            out[_cell(p, -j, other)] += 1
        if other != j:
            out[_cell(p, j, other)] -= 1
    return {cell: value for cell, value in out.items() if value}


def hard_star_difference_formula(p: int, j: int) -> dict[Cell, int]:
    """Return the indicator formula in equation (10) of the note."""
    _check_prime(p)
    j %= p
    squares = {pow(value, 2, p) for value in range(1, p)}
    out: dict[Cell, int] = {}
    for alpha in range(p):
        for beta in squares:
            value = int(pow(alpha + j, 2, p) == beta) - int(
                pow(alpha - j, 2, p) == beta
            )
            if value:
                out[(alpha, beta)] = value
    return out


def hard_star_boundary(p: int, j: int) -> dict[int, int]:
    """Return the signed vertex boundary of S_{-j}-S_j."""
    _check_prime(p)
    j %= p
    if j == 0:
        return {}
    return {(-j) % p: p - 2, j: -(p - 2)}


def star_moment_contraction(
    p: int, j: int, degree: int, channel: int
) -> int:
    """Directly replay one field-sum identity for a unit star."""
    _check_prime(p)
    if not 2 <= degree <= p - 2:
        raise ValueError("need 2 <= degree <= p-2")
    if not 0 <= channel < degree // 2:
        raise ValueError("invalid moment channel")
    j %= p
    total = 0
    for other in range(p):
        if other == j:
            continue
        total += (
            pow(j - other, 2, p)
            * pow(j * other, channel, p)
            * pow(
                j + other,
                degree - 2 - 2 * channel,
                p,
            )
        )
    return total % p


def ternary_defect(values: Iterable[int]) -> int:
    """Nonnegative integer defect vanishing exactly on {-1,0,1}."""
    total = 0
    for value in values:
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError("ternary defect requires integers")
        total += value * value * (value * value - 1) // 2
    return total


def _functional_value(p: int, functional: Functional, point: Point) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % p


def _edge(first: Point, second: Point) -> Edge:
    return tuple(sorted((first, second)))  # type: ignore[return-value]


def _negative_point(p: int, point: Point) -> Point:
    return (-point[0] % p, -point[1] % p)


def _negative_edge(p: int, edge: Edge) -> Edge:
    return _edge(_negative_point(p, edge[0]), _negative_point(p, edge[1]))


def _independent(p: int, first: Functional, second: Functional) -> bool:
    return (first[0] * second[1] - first[1] * second[0]) % p != 0


def _point_from_coordinates(
    p: int,
    first: Functional,
    second: Functional,
    first_value: int,
    second_value: int,
) -> Point:
    determinant = (first[0] * second[1] - first[1] * second[0]) % p
    if determinant == 0:
        raise ValueError("the two functionals must be independent")
    inverse = pow(determinant, -1, p)
    x = (
        second[1] * first_value - first[1] * second_value
    ) * inverse % p
    y = (
        -second[0] * first_value + first[0] * second_value
    ) * inverse % p
    return x, y


def projective_functionals(p: int) -> tuple[Functional, ...]:
    """Return one deterministic representative of every direction row."""
    _check_prime(p)
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def hard_star_chain(p: int, j: int) -> dict[tuple[int, int], int]:
    """Return ``A_j=S_{-j}-S_j`` in unordered fibre-edge coordinates."""
    _check_prime(p)
    j %= p
    out: defaultdict[tuple[int, int], int] = defaultdict(int)
    for other in range(p):
        if other != -j % p:
            out[tuple(sorted((-j % p, other)))] += 1
        if other != j:
            out[tuple(sorted((j, other)))] -= 1
    return {cell: value for cell, value in out.items() if value}


def localized_star_trade(
    p: int,
    direction: Functional,
    auxiliary: Functional,
    center: int,
) -> dict[Edge, int]:
    r"""Return a ternary anti-chain supported over exactly ``p-1`` orbits.

    In ``(direction,auxiliary)`` coordinates, put

    ``u_t=j*(1,t/(t+1))`` and ``v_t=j*(t,t)`` for ``t!=-1``.

    The returned chain is ``1_{-E}-1_E``.  Its ordinary edge--Radon image is
    ``A_j`` in ``direction`` and zero in every other direction.
    """
    _check_prime(p)
    center %= p
    if center == 0:
        return {}
    if not _independent(p, direction, auxiliary):
        raise ValueError("direction and auxiliary must be independent")

    source: defaultdict[Edge, int] = defaultdict(int)
    positive_half: set[Edge] = set()
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
        edge = _edge(first, second)
        negative = _negative_edge(p, edge)
        if first == second:
            raise ArithmeticError("the Mobius trade created a loop")
        if edge == negative:
            raise ArithmeticError("the Mobius trade created an antipodal edge")
        if edge in positive_half or negative in positive_half:
            raise ArithmeticError("the Mobius trade repeated an edge orbit")
        positive_half.add(edge)
        source[edge] -= 1
        source[negative] += 1

    out = {edge: value for edge, value in source.items() if value}
    if (
        len(positive_half) != p - 1
        or len(out) != 2 * (p - 1)
        or set(out.values()) != {-1, 1}
        or any(out.get(_negative_edge(p, edge)) != -value for edge, value in out.items())
    ):
        raise ArithmeticError("the Mobius trade lost ternary antisymmetry")
    return out


def edge_radon_image(
    p: int, source: dict[Edge, int]
) -> dict[tuple[object, ...], int]:
    """Return the ordinary unsigned edge--Radon image of a sparse chain."""
    _check_prime(p)
    target: defaultdict[tuple[object, ...], int] = defaultdict(int)
    for edge, coefficient in source.items():
        if edge[0] == edge[1]:
            raise ValueError("source chains cannot contain loops")
        for direction_index, functional in enumerate(projective_functionals(p)):
            first = _functional_value(p, functional, edge[0])
            second = _functional_value(p, functional, edge[1])
            if first == second:
                key: tuple[object, ...] = ("P", direction_index)
            else:
                key = ("K", direction_index, *sorted((first, second)))
            target[key] += coefficient
    return {key: value for key, value in target.items() if value}


def localized_star_trade_certificate(
    p: int,
    direction_index: int = 0,
    auxiliary: Functional = (0, 1),
    center: int = 1,
) -> dict[str, object]:
    """Replay one direct Mobius trade and its all-direction cancellation."""
    directions = projective_functionals(p)
    if not 0 <= direction_index < len(directions):
        raise ValueError("direction index is out of range")
    direction = directions[direction_index]
    source = localized_star_trade(p, direction, auxiliary, center)
    image = edge_radon_image(p, source)
    expected = {
        ("K", direction_index, *cell): value
        for cell, value in hard_star_chain(p, center).items()
    }
    proved = image == expected
    if not proved:
        raise ArithmeticError("the direction-localized star trade changed")
    return {
        "p": p,
        "direction": list(direction),
        "auxiliary": list(auxiliary),
        "center": center % p,
        "parameter_domain": "F_p minus {-1}",
        "source_inversion_orbits": 0 if center % p == 0 else p - 1,
        "source_actual_edges": len(source),
        "target_nonzero_direction_count": 0 if center % p == 0 else 1,
        "target": "A_j=S_-j-S_j",
        "proved": proved,
    }


def mobius_pairing_certificate(p: int) -> dict[str, object]:
    """Replay the three projective cancellation involutions symbolically mod p."""
    _check_prime(p)
    domain = tuple(value for value in range(p) if value != p - 1)
    generic_checks = []
    for multiplier in range(1, p):
        if multiplier == p - 1:
            continue
        c = multiplier + 1
        for parameter in domain:
            fraction = parameter * pow(parameter + 1, -1, p) % p
            first = (1 + multiplier * fraction) % p
            second = c * parameter % p
            paired = -(c * parameter + 1) * pow(
                c * (parameter + 1), -1, p
            ) % p
            paired_fraction = paired * pow(paired + 1, -1, p) % p
            generic_checks.append(
                paired != p - 1
                and c * paired % p == -first % p
                and (1 + multiplier * paired_fraction) % p == -second % p
            )

    minus_one_checks = []
    auxiliary_checks = []
    for parameter in domain:
        paired = (-parameter - 2) % p
        minus_one_checks.append(
            paired != p - 1
            and pow(paired + 1, -1, p) == -pow(parameter + 1, -1, p) % p
            and (-paired - 2) % p == parameter
        )

        fraction = parameter * pow(parameter + 1, -1, p) % p
        auxiliary_pair = -fraction % p
        auxiliary_pair_fraction = auxiliary_pair * pow(
            auxiliary_pair + 1, -1, p
        ) % p
        auxiliary_checks.append(
            auxiliary_pair != p - 1
            and auxiliary_pair_fraction == -parameter % p
            and (-auxiliary_pair_fraction) % p == parameter
        )

    proved = all(generic_checks + minus_one_checks + auxiliary_checks)
    if not proved:
        raise ArithmeticError("a Mobius cancellation involution changed")
    return {
        "p": p,
        "generic_projective_rows_checked": len(generic_checks),
        "exceptional_multiplier_minus_one_checked": len(minus_one_checks),
        "auxiliary_direction_checked": len(auxiliary_checks),
        "all_involutions_preserve_Fp_minus_minus_one": True,
        "fixed_points_give_centrally_fixed_cells": True,
        "proved": proved,
    }


def greedy_auxiliary_margin(p: int) -> dict[str, int | bool]:
    """Return the exact support-avoidance count for at most half the rows."""
    _check_prime(p)
    maximum_trades = (p + 1) // 2
    auxiliary_choices = p * (p - 1)
    forbidden_per_prior_trade = 2 * p - 1
    maximum_prior_trades = maximum_trades - 1
    maximum_forbidden = maximum_prior_trades * forbidden_per_prior_trade
    margin = auxiliary_choices - maximum_forbidden
    proved = margin == (p - 1) // 2 and margin > 0
    if not proved:
        raise ArithmeticError("the greedy auxiliary margin changed")
    return {
        "p": p,
        "maximum_target_directions": maximum_trades,
        "auxiliary_functional_choices": auxiliary_choices,
        "forbidden_choices_per_prior_trade_upper_bound": (
            forbidden_per_prior_trade
        ),
        "maximum_prior_trades": maximum_prior_trades,
        "maximum_forbidden_choices": maximum_forbidden,
        "guaranteed_remaining_auxiliaries": margin,
        "proved": proved,
    }


def simultaneous_localized_star_lift(
    p: int, targets: dict[int, int]
) -> dict[Edge, int]:
    """Construct disjoint ternary trades for at most ``(p+1)/2`` rows.

    ``targets`` maps projective-direction indices to the desired center.  A
    zero center is skipped because ``A_0=0``.  The deterministic greedy
    search is an implementation of the counting proof, not a SAT search.
    """
    directions = projective_functionals(p)
    if len(targets) > (p + 1) // 2:
        raise ValueError("at most (p+1)/2 target directions are allowed")
    if any(not 0 <= index < len(directions) for index in targets):
        raise ValueError("target direction index is out of range")

    used_orbits: set[Edge] = set()
    source: dict[Edge, int] = {}
    for direction_index in sorted(targets):
        center = targets[direction_index] % p
        if center == 0:
            continue
        direction = directions[direction_index]
        selected: dict[Edge, int] | None = None
        for auxiliary in (
            (x, y) for x in range(p) for y in range(p)
            if _independent(p, direction, (x, y))
        ):
            candidate = localized_star_trade(
                p, direction, auxiliary, center
            )
            candidate_orbits = {
                min(edge, _negative_edge(p, edge)) for edge in candidate
            }
            if candidate_orbits.isdisjoint(used_orbits):
                selected = candidate
                used_orbits.update(candidate_orbits)
                break
        if selected is None:
            raise ArithmeticError("the proved greedy auxiliary choice failed")
        source.update(selected)

    expected_orbits = sum(center % p != 0 for center in targets.values()) * (p - 1)
    if (
        len(used_orbits) != expected_orbits
        or len(source) != 2 * expected_orbits
        or any(value not in (-1, 1) for value in source.values())
    ):
        raise ArithmeticError("the simultaneous ternary supports overlapped")
    return source


def theorem_record(p: int = 31) -> dict[str, object]:
    """Return the exact reduction and explicit OPEN flags."""
    dimensions = antisymmetric_dimensions(p)
    star_checks = {
        str(j): (
            hard_star_difference_direct(p, j)
            == hard_star_difference_formula(p, j)
        )
        for j in (0, 1, p // 2)
    }
    odd_star_checks = {
        f"d{degree}k{channel}": star_moment_contraction(
            p, 1, degree, channel
        )
        for degree in range(3, min(p - 1, 10), 2)
        for channel in range(degree // 2)
    }
    pairing = mobius_pairing_certificate(p)
    localized = localized_star_trade_certificate(p)
    greedy = greedy_auxiliary_margin(p)
    proved = bool(
        dimensions["proved"]
        and all(star_checks.values())
        and all(value == 0 for value in odd_star_checks.values())
        and ternary_defect((-1, 0, 1)) == 0
        and ternary_defect((2,)) == 6
        and pairing["proved"]
        and localized["proved"]
        and greedy["proved"]
    )
    if not proved:
        raise ArithmeticError("antisymmetric theorem replay changed")
    return {
        "title": "Central-inversion antisymmetric edge-Radon reduction",
        "method": "symbolic rank and exact field-sum identities; no census",
        "dimensions": dimensions,
        "integral_cokernel": (
            "(Z/pZ)^[h(h-1)(h+1)/3], exactly the odd moment rows"
        ),
        "central_opposite_consequence": (
            "hard antisymmetric target has an unrestricted integral lift "
            "when all odd global forms vanish"
        ),
        "hard_star_target": (
            "A_j(alpha,beta)=1[(alpha+j)^2=beta]"
            "-1[(alpha-j)^2=beta]"
        ),
        "star_formula_checks": star_checks,
        "odd_star_moment_checks": odd_star_checks,
        "direction_localized_trade": localized,
        "mobius_pairing": pairing,
        "disjoint_support_greedy_bound": greedy,
        "antisymmetric_hard_star_ternary_lift_proved": True,
        "remaining_exact_gate": (
            "the coupled symmetric half, with s_e=1 on every nonzero "
            "Mobius-trade orbit and s_e in {0,2} elsewhere"
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
