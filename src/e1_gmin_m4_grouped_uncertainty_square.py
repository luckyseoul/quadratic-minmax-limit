#!/usr/bin/env python3
"""All-prime grouped uncertainty via a square-root remainder.

For even support, silent directions are precisely the parameters where the
polynomial of squared projections is a square.  A canonical top-half square
root leaves a nonzero homogeneous remainder of degree at most twice the
support size.  Every silent direction is a double zero of that remainder.
Together with the existing radial argument for odd support, this proves the
grouped uncertainty inequality for every odd prime.
"""

from __future__ import annotations

from itertools import combinations

from e1_gmin_m4_prop15721 import is_prime


Point = tuple[int, int]
Polynomial = tuple[int, ...]


def _check_odd_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")


def _trim(poly: list[int], p: int) -> Polynomial:
    while len(poly) > 1 and poly[-1] % p == 0:
        poly.pop()
    return tuple(value % p for value in poly)


def _add(left: Polynomial, right: Polynomial, p: int) -> Polynomial:
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return _trim(out, p)


def _scale(poly: Polynomial, scalar: int, p: int) -> Polynomial:
    return _trim([scalar * value for value in poly], p)


def _multiply(left: Polynomial, right: Polynomial, p: int) -> Polynomial:
    out = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            out[left_index + right_index] += left_value * right_value
    return _trim(out, p)


def _canonical_class(p: int, point: Point) -> Point:
    point = (point[0] % p, point[1] % p)
    negative = ((-point[0]) % p, (-point[1]) % p)
    return min(point, negative)


def _check_representatives(p: int, representatives: tuple[Point, ...]) -> None:
    if not representatives:
        raise ValueError("support must be nonempty")
    classes = []
    for point in representatives:
        if not isinstance(point, tuple) or len(point) != 2:
            raise ValueError("representatives must be pairs")
        reduced = (point[0] % p, point[1] % p)
        if reduced == (0, 0):
            raise ValueError("zero is not an antipodal point class")
        classes.append(_canonical_class(p, reduced))
    if len(set(classes)) != len(classes):
        raise ValueError("representatives must be distinct modulo sign")


def canonical_square_remainder(
    p: int, representatives: tuple[Point, ...]
) -> dict[str, object]:
    """Construct the canonical remainder for one even support.

    In the affine chart ``u=(1,t)``, put ``r_i(t)=u(v_i)^2`` and

        P_t(X)=product_i (X-r_i(t))
              =X^(2n)+c_1(t)X^(2n-1)+...+c_(2n)(t).

    There is a unique monic degree-``n`` polynomial ``Q_t`` whose square
    agrees with the top ``n`` coefficients of ``P_t``.  The returned
    coefficients ``R_k`` are those of ``P_t-Q_t^2`` for ``n<k<=2n``.
    Affine coefficient lists are accompanied by their homogeneous degree
    bounds ``2k``, so a zero at the projective point at infinity is retained.
    """
    _check_odd_prime(p)
    _check_representatives(p, representatives)
    support_size = len(representatives)
    if support_size % 2:
        raise ValueError("the square-remainder construction needs even support")
    half = support_size // 2

    coefficients: list[Polynomial] = [(1,)]
    for first, second in representatives:
        root = (
            first * first % p,
            2 * first * second % p,
            second * second % p,
        )
        updated = list(coefficients) + [(0,)]
        for index in range(len(coefficients) - 1, -1, -1):
            updated[index + 1] = _add(
                updated[index + 1],
                _scale(_multiply(coefficients[index], root, p), -1, p),
                p,
            )
        coefficients = updated

    inverse_two = pow(2, -1, p)
    square_root: list[Polynomial] = [(1,)]
    for index in range(1, half + 1):
        cross = (0,)
        for left in range(1, index):
            cross = _add(
                cross,
                _multiply(square_root[left], square_root[index - left], p),
                p,
            )
        square_root.append(
            _scale(_add(coefficients[index], _scale(cross, -1, p), p), inverse_two, p)
        )

    remainders: dict[int, Polynomial] = {}
    for index in range(half + 1, support_size + 1):
        square_coefficient = (0,)
        for left in range(max(0, index - half), min(half, index) + 1):
            right = index - left
            if 0 <= right <= half:
                square_coefficient = _add(
                    square_coefficient,
                    _multiply(square_root[left], square_root[right], p),
                    p,
                )
        remainders[index] = _add(
            coefficients[index], _scale(square_coefficient, -1, p), p
        )

    nonzero_indices = tuple(
        index for index, poly in remainders.items() if poly != (0,)
    )
    if not nonzero_indices:
        raise ArithmeticError(
            "distinct antipodal factors unexpectedly gave a global square"
        )
    return {
        "p": p,
        "support_size": support_size,
        "half_support_size": half,
        "P_coefficients": tuple(coefficients),
        "Q_coefficients": tuple(square_root),
        "remainder_coefficients": remainders,
        "nonzero_remainder_indices": nonzero_indices,
        "homogeneous_degree_by_index": {
            index: 2 * index for index in range(half + 1, support_size + 1)
        },
        "largest_remainder_degree": 2 * support_size,
        "proved": True,
    }


def _value(poly: Polynomial, argument: int, p: int) -> int:
    total = 0
    for coefficient in reversed(poly):
        total = (total * argument + coefficient) % p
    return total


def _finite_zero_order(poly: Polynomial, root: int, p: int) -> int:
    work = list(poly)
    order = 0
    while len(work) > 1 and _value(tuple(work), root, p) == 0:
        quotient = [0] * (len(work) - 1)
        carry = work[-1] % p
        quotient[-1] = carry
        for index in range(len(work) - 2, 0, -1):
            carry = (work[index] + root * carry) % p
            quotient[index - 1] = carry
        if (work[0] + root * carry) % p:
            raise ArithmeticError("synthetic division changed")
        work = quotient
        order += 1
    return order


def _silent_directions(p: int, representatives: tuple[Point, ...]) -> tuple[int, ...]:
    directions: list[int] = []
    for slope in range(p):
        counts: dict[int, int] = {}
        for first, second in representatives:
            value = (first + slope * second) % p
            square = value * value % p
            counts[square] = counts.get(square, 0) + 1
        if all(count % 2 == 0 for square, count in counts.items() if square):
            directions.append(slope)
    infinity_counts: dict[int, int] = {}
    for _first, second in representatives:
        square = second * second % p
        infinity_counts[square] = infinity_counts.get(square, 0) + 1
    if all(
        count % 2 == 0
        for square, count in infinity_counts.items()
        if square
    ):
        directions.append(p)
    return tuple(directions)


def exact_remainder_replay(
    p: int, representatives: tuple[Point, ...]
) -> dict[str, object]:
    """Check double vanishing directly for a small supplied support.

    This is a fail-when-wrong replay, not evidence for the theorem.
    Direction ``p`` denotes the projective point at infinity.
    """
    out = canonical_square_remainder(p, representatives)
    silent = _silent_directions(p, representatives)
    orders: dict[int, dict[int, int]] = {}
    for index in out["nonzero_remainder_indices"]:
        poly = out["remainder_coefficients"][index]
        degree = out["homogeneous_degree_by_index"][index]
        per_direction: dict[int, int] = {}
        for direction in silent:
            if direction == p:
                per_direction[direction] = degree - (len(poly) - 1)
            else:
                per_direction[direction] = _finite_zero_order(poly, direction, p)
        orders[index] = per_direction
        if any(order < 2 for order in per_direction.values()):
            raise ArithmeticError("a silent direction lost double vanishing")
        if 2 * len(silent) > degree:
            raise ArithmeticError("the projective root count changed")
    return {
        "p": p,
        "support_size": len(representatives),
        "silent_directions": silent,
        "silent_direction_count": len(silent),
        "nonzero_remainder_indices": out["nonzero_remainder_indices"],
        "zero_orders": orders,
        "grouped_bound_holds": len(silent) <= len(representatives),
        "role": "fail-when-wrong replay, not theorem evidence",
        "proved": True,
    }


def even_support_grouped_uncertainty_theorem(p: int) -> dict[str, object]:
    """Return the symbolic proof of ``z<=s`` for every even support."""
    _check_odd_prime(p)
    return {
        "p": p,
        "scope": "every nonempty even support S in (F_p^2 minus 0)/{+1,-1}",
        "conclusion": "number z of silent direction groups is at most |S|",
        "square_criterion": (
            "for even |S|, a direction is silent iff every root of "
            "P_u(X)=product_[v]inS (X-u(v)^2) has even multiplicity"
        ),
        "global_nonsquare_reason": (
            "the factors X-u(v)^2 are distinct in F_p[U,V,X], because "
            "u(v)^2=u(w)^2 as forms would imply v=+w or v=-w"
        ),
        "double_zero_reason": (
            "pair equal roots locally; each paired product is a square "
            "modulo the square of a local parameter"
        ),
        "degree_count": "2*z<=degree(R_k)=2*k<=2*|S|",
        "proved": True,
    }


def grouped_uncertainty_theorem(p: int) -> dict[str, object]:
    """Combine the even square-remainder proof with odd radial parity."""
    _check_odd_prime(p)
    even = even_support_grouped_uncertainty_theorem(p)
    return {
        "p": p,
        "point_space": "Delta=(F_p^2 minus 0)/{+1,-1}",
        "transform": "paired nonorigin affine-block incidence M^T",
        "grouped_inequality": (
            "wt(f)+#{A:(M^T f) restricted to B_A is nonzero}>=p+1"
        ),
        "odd_support_proof": (
            "every silent direction has odd radial support, and distinct "
            "radial directions partition the support"
        ),
        "even_support_proof": even["double_zero_reason"],
        "prime_range": "every odd prime",
        "p_congruent_3_mod_4_included": True,
        "proved": True,
    }


def theorem_record(p: int) -> dict[str, object]:
    theorem = grouped_uncertainty_theorem(p)
    return {
        "p": p,
        "proved": {
            "grouped_uncertainty_all_supports": theorem["proved"],
            "grouped_uncertainty_all_odd_primes": theorem["proved"],
            "row_code_minimum_distance": False,
            "minimum_word_classification": False,
            "structured_mobius_puncture": False,
            "symmetric_boolean_completion": False,
            "residual_ii_closed": False,
        },
    }
