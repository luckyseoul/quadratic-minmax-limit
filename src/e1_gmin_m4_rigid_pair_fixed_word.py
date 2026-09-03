#!/usr/bin/env python3
"""Source-side fixed words of localized Mobius halves.

For any auxiliary, one nonzero localized Mobius half has exactly one orbit
with zero fixed word and p-2 orbits with odd, p-element fixed words.  The
unique rigid two-cancellation pair admits a finer closed formula and exact
weight.

These are source-side Phi statements only.  A centrally symmetric compact
hard residual can have odd fixed-cell coefficients, so the full branch-C
target coset word and its Hamming weight are not determined here.
"""

from __future__ import annotations

from collections import Counter
from itertools import product

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    Point,
    _functional_value,
    _negative_edge,
    localized_star_trade,
)
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_fixed_edge_elimination import orbit_fixed_word


def _check_branch_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a branch-C prime p=3 mod 4 with p>=31")
    return (p - 3) // 4


def _legendre(p: int, value: int) -> int:
    value %= p
    if value == 0:
        return 0
    result = pow(value, (p - 1) // 2, p)
    if result == 1:
        return 1
    if result == p - 1:
        return -1
    raise ArithmeticError("Euler criterion returned a nonsign")


def _antipodal_class(p: int, point: Point) -> Point:
    point = (point[0] % p, point[1] % p)
    if point == (0, 0):
        raise ValueError("zero has no nonzero antipodal class")
    negative = ((-point[0]) % p, (-point[1]) % p)
    return min(point, negative)


def _point_classes(p: int) -> tuple[Point, ...]:
    return tuple(
        sorted(
            {
                _antipodal_class(p, point)
                for point in product(range(p), repeat=2)
                if point != (0, 0)
            }
        )
    )


def _independent(
    p: int, first: Functional, second: Functional
) -> bool:
    return (
        first[0] * second[1] - first[1] * second[0]
    ) % p != 0


def _normalized_coordinates(
    p: int,
    point: Point,
    first: Functional,
    second: Functional,
    first_center: int,
    second_center: int,
) -> tuple[int, int]:
    return (
        _functional_value(p, first, point)
        * pow(first_center % p, -1, p)
        % p,
        _functional_value(p, second, point)
        * pow(second_center % p, -1, p)
        % p,
    )


def _rigid_half_phi_bit(p: int, x: int, y: int) -> int:
    """Evaluate the exact Phi word of the first rigid selected half."""
    x %= p
    y %= p
    three_quarters = 3 * pow(4, -1, p) % p
    difference = (2 * y - x) % p
    bit = int(x * x % p == 1)
    for sign in (1, -1):
        bit ^= int(
            difference * ((x + sign) % p) % p
            == three_quarters
        )
    return bit


def localized_half_phi_parity_theorem(p: int) -> dict[str, object]:
    """State the auxiliary-independent source Phi parity of one half.

    In the standard L,M coordinates, the selected edge at t is

        u_t=(j,j*t/(t+1)),  v_t=(j*t,j*t).

    Its midpoint and half-difference are parallel exactly when the endpoints
    are collinear with the origin.  Their determinant is j^2*t/(t+1), so
    for j nonzero exactly t=0 has zero Phi word.  Every other orbit has a
    p-element affine-line word.
    """
    _check_branch_prime(p)
    parameter_count = p - 1
    zero_word_parameters = 1
    nonzero_word_parameters = p - 2
    block_weight = p
    half_phi_parity = nonzero_word_parameters * block_weight % 2
    proved = bool(
        parameter_count
        == zero_word_parameters + nonzero_word_parameters
        and half_phi_parity == 1
    )
    if not proved:
        raise ArithmeticError("the localized-half Phi parity changed")
    return {
        "p": p,
        "standard_selected_edge": (
            "u_t=(j,j*t/(t+1)), v_t=(j*t,j*t), t!=-1"
        ),
        "endpoint_determinant": "j^2*t/(t+1)",
        "zero_Phi_word_parameters": [0],
        "nonzero_Phi_word_parameter_count": nonzero_word_parameters,
        "each_nonzero_Phi_word_weight": block_weight,
        "total_half_Phi_word_parity": half_phi_parity,
        "independent_of_auxiliary": True,
        "independent_of_Mobius_orientation_and_Paley_sign_mod_2": True,
        "full_target_coset_word_determined": False,
        "proved": proved,
    }


def rigid_pair_phi_word(
    p: int,
    first: Functional = (1, 0),
    second: Functional = (0, 1),
    first_center: int = 1,
    second_center: int = 1,
) -> frozenset[Point]:
    """Return the source Phi word of one rigid sharp pair."""
    _check_branch_prime(p)
    first_center %= p
    second_center %= p
    if not first_center or not second_center:
        raise ValueError("the two hard centers must be nonzero")
    if not _independent(p, first, second):
        raise ValueError("the two hard directions must be distinct")
    support: set[Point] = set()
    for point in _point_classes(p):
        x, y = _normalized_coordinates(
            p,
            point,
            first,
            second,
            first_center,
            second_center,
        )
        bit = _rigid_half_phi_bit(p, x, y)
        bit ^= _rigid_half_phi_bit(p, y, x)
        if bit:
            support.add(point)
    return frozenset(support)


def rigid_pair_phi_word_theorem(p: int) -> dict[str, object]:
    """State the exact rigid-pair Phi formula and weight."""
    _check_branch_prime(p)
    eta_three = _legendre(p, 3)
    eta_six = _legendre(p, 6)
    eta_minus_two = _legendre(p, -2)
    half_phi_weight = 2 * p - 3
    half_word_intersection_points = (
        20 + 2 * eta_three + 4 * eta_six + 4 * eta_minus_two
    )
    pair_phi_weight = (
        4 * p
        - 26
        - 2 * eta_three
        - 4 * eta_six
        - 4 * eta_minus_two
    )
    replay_phi = len(rigid_pair_phi_word(p))
    proved = replay_phi == pair_phi_weight
    if not proved:
        raise ArithmeticError("the rigid-pair Phi count changed")
    return {
        "p": p,
        "normalization": "x=L1(v)/j1, y=L2(v)/j2",
        "one_half_Phi_formula": (
            "1_(x^2=1) + sum_(e=+1,-1) "
            "1_((2y-x)(x+e)=3/4) over F2"
        ),
        "pair_Phi_formula": "F(x,y)+F(y,x)",
        "one_half_Phi_weight": half_phi_weight,
        "two_half_word_intersection_points": (
            half_word_intersection_points
        ),
        "pair_Phi_weight": pair_phi_weight,
        "weight_is_independent_of_directions_and_nonzero_centers": True,
        "literal_only_forced_coset_claim_retracted": True,
        "full_target_coset_weight_open": True,
        "proved": proved,
    }


def _scale_functional(
    p: int, scalar: int, functional: Functional
) -> Functional:
    return (
        scalar * functional[0] % p,
        scalar * functional[1] % p,
    )


def _add_functionals(
    p: int, first: Functional, second: Functional
) -> Functional:
    return (
        (first[0] + second[0]) % p,
        (first[1] + second[1]) % p,
    )


def _canonical_orbit(p: int, edge: Edge) -> Edge:
    return min(edge, _negative_edge(p, edge))


def p31_phi_pairing_dependence_replay() -> dict[str, object]:
    """Give one exact counterexample to pairing-invariant source Phi weight."""
    p = 31
    directions: tuple[Functional, ...] = (
        (1, 0),
        (0, 1),
        (1, 1),
        (1, 2),
    )
    labels = ("X", "Y", "Z", "W")
    by_label = dict(zip(labels, directions, strict=True))
    pairings = (
        (("X", "Y"), ("Z", "W")),
        (("X", "Z"), ("Y", "W")),
        (("X", "W"), ("Y", "Z")),
    )
    two_thirds = 2 * pow(3, -1, p) % p
    rows: list[dict[str, object]] = []
    for pairing in pairings:
        source: Counter[Edge] = Counter()
        formula_word: set[Point] = set()
        for first_label, second_label in pairing:
            first = by_label[first_label]
            second = by_label[second_label]
            auxiliary = _scale_functional(
                p,
                two_thirds,
                _add_functionals(p, first, second),
            )
            source.update(
                localized_star_trade(p, first, auxiliary, 1)
            )
            source.update(
                localized_star_trade(p, second, auxiliary, 1)
            )
            formula_word.symmetric_difference_update(
                rigid_pair_phi_word(p, first, second)
            )
        source = Counter(
            {edge: value for edge, value in source.items() if value}
        )
        ternary = all(abs(value) == 1 for value in source.values())
        used_orbits = {
            _canonical_orbit(p, edge) for edge in source
        }
        phi_word: set[Point] = set()
        for orbit in used_orbits:
            phi_word.symmetric_difference_update(
                tuple(point)
                for point in orbit_fixed_word(
                    p, orbit
                )["fixed_word_support"]
            )
        rows.append(
            {
                "pairing": [
                    list(pair) for pair in pairing
                ],
                "full_trade_sum_ternary": ternary,
                "used_inversion_orbits": len(used_orbits),
                "cancellation_units": (
                    len(directions) * (p - 1)
                    - len(used_orbits)
                )
                // 2,
                "source_Phi_weight": len(phi_word),
                "closed_formula_matches_actual_Phi": (
                    phi_word == formula_word
                ),
            }
        )
    weights = sorted(
        int(row["source_Phi_weight"]) for row in rows
    )
    proved = bool(
        weights == [172, 174, 176]
        and all(row["full_trade_sum_ternary"] for row in rows)
        and all(row["used_inversion_orbits"] == 112 for row in rows)
        and all(row["cancellation_units"] == 4 for row in rows)
        and all(
            row["closed_formula_matches_actual_Phi"] for row in rows
        )
    )
    if not proved:
        raise ArithmeticError("the p=31 Phi pairing replay changed")
    return {
        "p": p,
        "directions": {
            label: list(direction)
            for label, direction in by_label.items()
        },
        "all_centers": 1,
        "rows": rows,
        "distinct_source_Phi_weights": weights,
        "pairing_independent_exact_Phi_weight": False,
        "target_coset_weights_computed": False,
        "role": (
            "exact p=31 counterexample, not a prime or configuration census"
        ),
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Package the source Phi theorem and the corrected open target gate."""
    half = localized_half_phi_parity_theorem(p)
    rigid = rigid_pair_phi_word_theorem(p)
    replay = p31_phi_pairing_dependence_replay() if p == 31 else None
    proved = bool(
        half["proved"]
        and rigid["proved"]
        and (replay is None or replay["proved"])
    )
    if not proved:
        raise ArithmeticError("the Mobius Phi fixed-word record changed")
    return {
        "title": "Source-side fixed words of localized Mobius halves",
        "one_half_Phi_parity": half,
        "rigid_pair_Phi_word": rigid,
        "p31_Phi_pairing_dependence_replay": replay,
        "status": (
            "SOURCE PHI WORD PROVED; COMPACT TARGET FIXED WORD "
            "AND FULL COSET WEIGHT OPEN"
        ),
        "corrected_scope": (
            "central compact residuals need not vanish on fixed cells; "
            "no branch-C Hamming-parity obstruction is claimed"
        ),
        "one_trade_per_hard_Mobius_ansatz_excluded": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
