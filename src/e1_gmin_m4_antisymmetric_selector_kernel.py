#!/usr/bin/env python3
"""Boundary-selector behavior on the antisymmetric edge--Radon kernel.

There are two different boundary maps here, and conflating them gives a
false extension of the adaptive Mobius obstruction.

* The mod-two boundary of the *full signed edge chain* is fixed by the
  ordinary edge--Radon target.  Indeed every Type-P and Type-K ridge has
  zero mod-two boundary, and odd-prime ridge p-saturation then proves this
  on the complete integral kernel.
* The boundary of the *physical positive half* of a ternary antisymmetric
  chain is not fixed.  A Type-K circuit below has zero antisymmetric Radon
  target, while its positive half has projective-kernel selector signature
  ``(0,1,...,1)`` and hence odd aggregate parity.

The countercircuit works for every odd prime for the unsigned Radon map.  At
``p=31`` both difference classes have Paley sign ``+1``, so it also survives
unchanged in the normalized residual convention.  This blocks promotion of
the localized-Mobius selector obstruction to arbitrary antisymmetric lifts;
it does not construct a residual-(ii) graph.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from itertools import combinations
from typing import Iterable

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Point,
    _edge,
    _functional_value,
    _negative_edge,
    _point_from_coordinates,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import paley_edge_sign
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


def _translated_edge(p: int, midpoint: Point, difference: Point) -> Edge:
    return _edge(
        (
            (midpoint[0] - difference[0]) % p,
            (midpoint[1] - difference[1]) % p,
        ),
        (
            (midpoint[0] + difference[0]) % p,
            (midpoint[1] + difference[1]) % p,
        ),
    )


def mod_two_boundary(edges: Iterable[Edge]) -> frozenset[Point]:
    """Return the odd-degree vertex set of an edge iterable."""
    boundary: set[Point] = set()
    for first, second in edges:
        boundary.symmetric_difference_update((first, second))
    return frozenset(boundary)


def kernel_selector(p: int, direction_index: int) -> frozenset[Point]:
    """Choose zero and the lexicographic member of each pair in ``ker L``."""
    _check_odd_prime(p)
    directions = projective_functionals(p)
    if not 0 <= direction_index < len(directions):
        raise ValueError("direction index is out of range")
    functional = directions[direction_index]
    selector = {(0, 0)}
    for x in range(p):
        for y in range(p):
            point = (x, y)
            if point == (0, 0) or _functional_value(p, functional, point):
                continue
            negative = (-x % p, -y % p)
            selector.add(min(point, negative))
    expected_size = (p + 1) // 2
    if len(selector) != expected_size:
        raise ArithmeticError("a projective-kernel selector changed size")
    return frozenset(selector)


def boundary_selector_signature(
    p: int, boundary: Iterable[Point]
) -> tuple[int, ...]:
    """Pair a centrally symmetric boundary with all kernel selectors."""
    _check_odd_prime(p)
    support = frozenset(boundary)
    for point in support:
        if (-point[0] % p, -point[1] % p) not in support:
            raise ValueError("selector pairing requires a central boundary")
    return tuple(
        len(support & kernel_selector(p, index)) % 2
        for index in range(p + 1)
    )


def antisymmetric_transverse_countercircuit(
    p: int,
) -> tuple[dict[Edge, int], frozenset[Edge]]:
    r"""Return ``z=1_H-1_{-H}`` in ``ker_Z R cap E^-``.

    Take ``L(x,y)=x``, half-differences ``delta_1=(1,0)`` and
    ``delta_2=(1,1)``, and the odd ridge profile
    ``g=1_{1}-1_{-1}``.  The positive half is

    ``H={(a,delta_1):a=(1,y)} union {(a,delta_2):a=(-1,y)}``.

    The two differences have the same nonzero squared L-projection, so this
    is a Type-K edge--Radon circuit.  Its displayed construction also makes
    inversion antisymmetry literal, without choosing orbit coordinates.
    """
    _check_odd_prime(p)
    first_difference = (1, 0)
    second_difference = (1, 1)
    positive: set[Edge] = set()
    for y in range(p):
        positive.add(_translated_edge(p, (1, y), first_difference))
        positive.add(_translated_edge(p, (-1 % p, y), second_difference))
    negative = {_negative_edge(p, edge) for edge in positive}
    if len(positive) != 2 * p or positive & negative:
        raise ArithmeticError("the positive half lost inversion orbits")
    source: defaultdict[Edge, int] = defaultdict(int)
    for edge in positive:
        source[edge] += 1
        source[_negative_edge(p, edge)] -= 1
    out = {edge: value for edge, value in source.items() if value}
    if (
        len(out) != 4 * p
        or set(out.values()) != {-1, 1}
        or any(out.get(_negative_edge(p, edge)) != -value for edge, value in out.items())
    ):
        raise ArithmeticError("the countercircuit lost ternary antisymmetry")
    return out, frozenset(positive)


def _antisymmetric_chain_from_half(
    p: int, positive: Iterable[Edge]
) -> dict[Edge, int]:
    """Return ``1_H-1_{-H}``, rejecting repeated inversion orbits."""
    half = frozenset(positive)
    negative = frozenset(_negative_edge(p, edge) for edge in half)
    if half & negative:
        raise ValueError("a physical half may use only one edge per inversion orbit")
    source: defaultdict[Edge, int] = defaultdict(int)
    for edge in half:
        source[edge] += 1
        source[_negative_edge(p, edge)] -= 1
    out = {edge: value for edge, value in source.items() if value}
    if len(out) != 2 * len(half) or set(out.values()) != {-1, 1}:
        raise ArithmeticError("the physical half did not give a ternary anti-chain")
    return out


def _coordinate_point(
    p: int,
    first_direction_index: int,
    second_direction_index: int,
    first_value: int,
    second_value: int,
) -> Point:
    directions = projective_functionals(p)
    if not 0 <= first_direction_index < p + 1 or not 0 <= second_direction_index < p + 1:
        raise ValueError("direction index is out of range")
    if first_direction_index == second_direction_index:
        raise ValueError("the two projective directions must be distinct")
    return _point_from_coordinates(
        p,
        directions[first_direction_index],
        directions[second_direction_index],
        first_value % p,
        second_value % p,
    )


def transverse_countercircuit_half(
    p: int, first_direction_index: int, second_direction_index: int
) -> frozenset[Edge]:
    """Return the Type-K positive half in arbitrary ``(L,M)`` coordinates."""
    _check_odd_prime(p)
    delta_one = _coordinate_point(
        p, first_direction_index, second_direction_index, 1, 0
    )
    delta_two = _coordinate_point(
        p, first_direction_index, second_direction_index, 1, 1
    )
    positive = set()
    for value in range(p):
        first_midpoint = _coordinate_point(
            p, first_direction_index, second_direction_index, 1, value
        )
        second_midpoint = _coordinate_point(
            p, first_direction_index, second_direction_index, -1, value
        )
        positive.add(_translated_edge(p, first_midpoint, delta_one))
        positive.add(_translated_edge(p, second_midpoint, delta_two))
    if len(positive) != 2 * p:
        raise ArithmeticError("the transverse positive half changed size")
    source = _antisymmetric_chain_from_half(p, positive)
    if edge_radon_image(p, source):
        raise ArithmeticError("the transverse positive half left the anti-kernel")
    return frozenset(positive)


def clique_star_kernel_half(
    p: int, first_direction_index: int, second_direction_index: int
) -> frozenset[Edge]:
    r"""Return the unsigned clique--star anti-kernel half.

    In ``(L,M)`` coordinates put ``c=1/2``.  Take the clique on
    ``{(c,y):y != 0}`` and the star from ``(c,0)`` to
    ``{(-c,y):y != 0}``.

    In direction ``L`` its projection is a parallel clique plus repeated
    copies of the inversion-fixed pair ``{c,-c}``; in direction ``M`` it is
    the complete graph on ``F_p``.  In a direction ``aL+bM`` with
    ``a*b != 0``, it contains every off-diagonal label pair except
    ``{ac,-ac}``, together with one parallel edge.  Every projected row is
    therefore invariant under label negation, proving ``R(H)=R(-H)``.
    """
    _check_odd_prime(p)
    inverse_two = pow(2, -1, p)
    positive_line = tuple(
        _coordinate_point(
            p, first_direction_index, second_direction_index, inverse_two, value
        )
        for value in range(1, p)
    )
    center = _coordinate_point(
        p, first_direction_index, second_direction_index, inverse_two, 0
    )
    negative_line = tuple(
        _coordinate_point(
            p, first_direction_index, second_direction_index, -inverse_two, value
        )
        for value in range(1, p)
    )
    positive = {
        _edge(first, second)
        for first, second in combinations(positive_line, 2)
    }
    positive.update(_edge(center, endpoint) for endpoint in negative_line)
    expected_size = p * (p - 1) // 2
    if len(positive) != expected_size:
        raise ArithmeticError("the clique--star half changed size")
    source = _antisymmetric_chain_from_half(p, positive)
    if edge_radon_image(p, source):
        raise ArithmeticError("the clique--star half left the anti-kernel")
    return frozenset(positive)


def unsigned_unit_selector_kernel_half(
    p: int, unit_direction_index: int, auxiliary_direction_index: int | None = None
) -> frozenset[Edge]:
    """Return a ternary unsigned anti-kernel half with a one-hot signature.

    The wanted unit direction is ``M``.  Any different direction can serve
    as ``L``; by default the first available canonical direction is used.
    The disjoint union of the transverse half and the clique--star half has
    selector word exactly ``e_M``.
    """
    _check_odd_prime(p)
    if not 0 <= unit_direction_index < p + 1:
        raise ValueError("unit direction index is out of range")
    if auxiliary_direction_index is None:
        auxiliary_direction_index = 0 if unit_direction_index != 0 else 1
    if (
        not 0 <= auxiliary_direction_index < p + 1
        or auxiliary_direction_index == unit_direction_index
    ):
        raise ValueError("the auxiliary direction must differ from the unit direction")
    transverse = transverse_countercircuit_half(
        p, auxiliary_direction_index, unit_direction_index
    )
    clique_star = clique_star_kernel_half(
        p, auxiliary_direction_index, unit_direction_index
    )
    if transverse & clique_star:
        raise ArithmeticError("the two kernel halves share a physical edge")
    positive = transverse | clique_star
    source = _antisymmetric_chain_from_half(p, positive)
    if edge_radon_image(p, source):
        raise ArithmeticError("the unit-word half left the anti-kernel")
    signature = boundary_selector_signature(p, mod_two_boundary(positive))
    expected = tuple(int(index == unit_direction_index) for index in range(p + 1))
    if signature != expected:
        raise ArithmeticError("the unit-word selector signature changed")
    return positive


def unsigned_selector_code_theorem(p: int) -> dict[str, object]:
    """Record the exact full unsigned selector-code theorem.

    This theorem deliberately does not assert the analogous statement after
    Paley source signing.  That signing preserves the anti-Radon kernel but
    can reverse the physical member selected from individual edge orbits.
    """
    _check_odd_prime(p)
    base_unit_direction = p
    positive = unsigned_unit_selector_kernel_half(p, base_unit_direction, 0)
    source = _antisymmetric_chain_from_half(p, positive)
    signature = boundary_selector_signature(p, mod_two_boundary(positive))
    expected_edge_count = p * (p + 3) // 2
    proved = bool(
        not edge_radon_image(p, source)
        and len(positive) == expected_edge_count
        and len(source) == 2 * expected_edge_count
        and signature == (0,) * p + (1,)
    )
    if not proved:
        raise ArithmeticError("the unsigned unit-selector theorem changed")
    return {
        "p": p,
        "map": "physical-half boundary followed by p+1 kernel selectors",
        "coordinate_convention": "ordinary unsigned anti-Radon map",
        "transverse_half_edges": 2 * p,
        "clique_star_half_edges": p * (p - 1) // 2,
        "unit_half_edges": len(positive),
        "unit_anti_chain_support": len(source),
        "base_unit_signature": list(signature),
        "every_projective_unit_realized": True,
        "transport": "choose the desired M and any independent L",
        "ternary_signature_set_contains_standard_basis": True,
        "induced_linear_code": "F_2^(p+1)",
        "induced_linear_code_dimension": p + 1,
        "minimum_nonzero_codeword_weight": 1,
        "normalized_Paley_code_fullness_proved": False,
        "normalized_Paley_minimum_weight": "OPEN",
        "p31_normalized_Type_K_barrier_retained": p == 31,
        "residual_ii_closed": False,
        "proved": proved,
    }


def full_chain_boundary_invariance_theorem(p: int) -> dict[str, object]:
    """Record the symbolic full-kernel boundary proof for an odd prime."""
    _check_odd_prime(p)
    return {
        "p": p,
        "type_P_boundary_mod_2": (
            "zero: a->a-delta and a->a+delta are two permutations of "
            "each midpoint fibre when L(delta)=0"
        ),
        "type_K_boundary_mod_2": (
            "zero: L(delta_1)^2=L(delta_2)^2 gives the same unordered "
            "pair of endpoint fibres for both difference classes"
        ),
        "ridge_p_saturation": "p*ker_Z(R) is contained in the ridge lattice",
        "odd_prime_reduction": "boundary(p*z)=boundary(z) over F_2",
        "conclusion": "boundary(z)=0 mod 2 for every z in ker_Z(R)",
        "applies_to_antisymmetric_kernel": True,
        "full_chain_boundary_target_determined": True,
        "full_chain_selector_signature_target_determined": True,
        "physical_half_signature_target_determined": False,
        "proved": True,
    }


def antisymmetric_selector_countercircuit_certificate(p: int) -> dict[str, object]:
    """Replay the all-odd-prime countercircuit and its selector word."""
    source, positive = antisymmetric_transverse_countercircuit(p)
    image = edge_radon_image(p, source)
    positive_boundary = mod_two_boundary(positive)
    expected_boundary = frozenset(
        {(2 % p, y) for y in range(p)}
        | {(-2 % p, y) for y in range(p)}
    )
    full_boundary = mod_two_boundary(source)
    signature = boundary_selector_signature(p, positive_boundary)
    expected_signature = (0,) + (1,) * p
    serialized_positive = json.dumps(
        [[list(first), list(second)] for first, second in sorted(positive)],
        separators=(",", ":"),
    ).encode()
    p31_signs: tuple[int, ...] | None = None
    if p == 31:
        p31_signs = tuple(sorted({paley_edge_sign(p, edge) for edge in source}))
    proved = bool(
        not image
        and not full_boundary
        and positive_boundary == expected_boundary
        and signature == expected_signature
        and sum(signature) % 2 == 1
        and (p != 31 or p31_signs == (1,))
    )
    if not proved:
        raise ArithmeticError("the antisymmetric selector countercircuit changed")
    return {
        "p": p,
        "ordinary_edge_Radon_image": {},
        "antisymmetric": True,
        "ternary": True,
        "difference_classes": [[1, 0], [1, 1]],
        "ridge_direction": [1, 0],
        "ridge_profile": "1_{x=1}-1_{x=-1}",
        "positive_half_edge_count": len(positive),
        "anti_chain_support": len(source),
        "positive_half_boundary": "1_{x=2}+1_{x=-2}",
        "positive_half_boundary_weight": len(positive_boundary),
        "full_chain_boundary_weight": len(full_boundary),
        "selector_signature": list(signature),
        "selector_signature_weight": sum(signature),
        "selector_signature_aggregate_mod_2": sum(signature) % 2,
        "zero_target_empty_lift_signature": [0] * (p + 1),
        "physical_half_signature_target_determined": False,
        "aggregate_selector_parity_target_determined": False,
        "positive_half_sha256": hashlib.sha256(serialized_positive).hexdigest(),
        "p31_paley_edge_sign_set": None if p31_signs is None else list(p31_signs),
        "p31_normalized_residual_countercircuit_unchanged": p == 31,
        "blocks_arbitrary_lift_extension_of_adaptive_Mobius_selector": True,
        "produces_residual_witness": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Return the positive full-chain result and negative half-chain result."""
    return {
        "full_chain": full_chain_boundary_invariance_theorem(p),
        "physical_half": antisymmetric_selector_countercircuit_certificate(p),
        "unsigned_physical_half_code": unsigned_selector_code_theorem(p),
    }


if __name__ == "__main__":
    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
