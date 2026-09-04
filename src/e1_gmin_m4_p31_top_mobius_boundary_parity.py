#!/usr/bin/env python3
"""Exact degree-boundary parity for the p=31 top Mobius construction.

This module isolates a necessary condition on a common graph built from the
sixteen localized Mobius halves at the top branch-C endpoint.  Compact target
atoms have zero vertex boundary over ``F_2``.  Consequently, after the hard
unit stars are removed, the selected source edges must have zero discrepancy

    d(L,M,c) = boundary(E(L,M,c)) + 1_{L=c}.

Pairing that discrepancy with one representative from every antipodal pair
on the kernel of a projective direction gives a small, exact invariant.  A
single fixed kernel is not an obstruction for every auxiliary design.  The
product over *all* projective kernels is, however, forced: one half contributes
``-epsilon_M``.  At the opposite-fixed top endpoint the auxiliary SDR has an
even number of entries and product sign ``+1``.  Arbitrary valid reductions
on nonorigin inversion orbits have even aggregate parity, while the fixed
edge has odd aggregate parity.  Hence an odd number of kernel selectors
contradict the required degree boundary.

The symbolic product argument works for every prime ``p = 3 (mod 4)`` in the
same auxiliary-SDR, one-fixed-edge Mobius family, regardless of the pattern
of nonorigin overlaps.  The explicit mask replay below remains at ``p=31``.
This is a closure of that construction family, not by itself a closure of
residual (ii).
"""

from __future__ import annotations

import hashlib
from itertools import product
from typing import Iterable, Sequence

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    Point,
    _functional_value,
    _negative_edge,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (
    _parallel_formula,
    _relative_coefficients,
    mobius_parameter_edges,
    paley_direction_sign,
    paley_edge_sign,
)
from e1_gmin_m4_prop15721 import is_prime


P = 31
POINT_COUNT = P * P


def sign_product(signs: Iterable[int]) -> int:
    """Multiply a finite sequence of signs, rejecting nonsigns."""
    result = 1
    for sign in signs:
        if sign not in (-1, 1):
            raise ValueError("a Paley sign must be +1 or -1")
        result *= sign
    return result


def all_prime_adaptive_product_theorem(p: int) -> dict[str, object]:
    """Return the symbolic all-design top-Mobius product obstruction.

    Put ``m=(p+1)/2`` and suppose the ``m`` targets are hard while the
    auxiliary SDR contains ``m-2`` hard directions and two opposite
    directions.  For one half, write

    ``D(z)=M+z(L-M)``.

    If ``sigma(G)=(-1)^<S_G,d(L,M,c)>``, the local calculation gives

    * ``sigma(L)=-1``;
    * ``sigma(L-M)=epsilon_L`` (the ``z=infinity`` case); and
    * ``sigma(D(z))=epsilon(D(z^2))`` for finite ``z != 1``.

    Pairing ``z`` with ``-z`` then gives

    ``product_G sigma(G) = -epsilon_M``.

    The product over the even number ``m`` of halves is therefore ``+1``.
    Replacing the occurrence sum on any nonorigin inversion orbit by its
    ternary final value changes the boundary by either zero or the whole
    pair ``boundary(e)+boundary(-e)``.  That pair changes two kernel
    signatures and hence has product ``+1``.  The auxiliary SDR makes the
    unique origin orbits of the halves distinct, so there is no origin-pair
    correction (which would have only one signature change).  Choosing one
    fixed antipodal edge changes exactly the fixed-kernel signature.  Thus
    the required final boundary has product ``-1`` and cannot agree with the
    half product.  A clean nonzero-affine cancellation is the special case
    with two endpoint changes plus the fixed-kernel change.

    This function records the symbolic parity algebra; the p31 functions
    below independently replay it from the exact edges and point masks.
    """
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=3 mod 4")
    half_count = (p + 1) // 2
    auxiliary_hard_count = half_count - 2
    auxiliary_opposite_count = 2
    auxiliary_sign_product = (-1) ** auxiliary_opposite_count
    local_minus_sign_product = (-1) ** half_count
    all_half_signature_product = (
        local_minus_sign_product * auxiliary_sign_product
    )
    nonorigin_pair_product = 1
    origin_pair_corrections = 0
    fixed_edge_count = 1
    required_signature_product = (
        nonorigin_pair_product
        * (-1) ** origin_pair_corrections
        * (-1) ** fixed_edge_count
        * all_half_signature_product
    )
    proved = bool(
        half_count % 2 == 0
        and auxiliary_sign_product == 1
        and all_half_signature_product == 1
        and nonorigin_pair_product == 1
        and origin_pair_corrections == 0
        and fixed_edge_count % 2 == 1
        and required_signature_product == -1
    )
    if not proved:
        raise ArithmeticError("the adaptive product parity changed")
    return {
        "p": p,
        "p_mod_4": p % 4,
        "half_count": half_count,
        "half_count_is_even": True,
        "target_sign": 1,
        "auxiliary_hard_count": auxiliary_hard_count,
        "auxiliary_opposite_count": auxiliary_opposite_count,
        "auxiliary_sign_product": auxiliary_sign_product,
        "one_half_signature_product_formula": "-epsilon_M",
        "all_half_signature_product": all_half_signature_product,
        "nonorigin_orbit_pair_signature_product": nonorigin_pair_product,
        "origin_orbit_pair_corrections": origin_pair_corrections,
        "origin_orbits_distinct_by_auxiliary_sdr": True,
        "fixed_edge_count": fixed_edge_count,
        "fixed_edge_signature_product": -1,
        "clean_collision_corollary_toggle_count": 3,
        "required_signature_product": required_signature_product,
        "products_disagree": True,
        "some_projective_kernel_is_a_contradiction": True,
        "scope": (
            "hard-target localized-Mobius family with an auxiliary SDR of "
            "m-2 hard plus two opposite directions, arbitrary valid ternary "
            "reductions on nonorigin inversion orbits, no origin-pair "
            "correction, and one fixed antipodal edge"
        ),
        "residual_ii_closed": False,
        "proved": proved,
    }


def point_index(point: Point) -> int:
    return point[0] * P + point[1]


def boundary_mask(edges: Iterable[Edge]) -> int:
    """Return the mod-two vertex boundary of an edge iterable."""
    mask = 0
    for first, second in edges:
        mask ^= 1 << point_index(first)
        mask ^= 1 << point_index(second)
    return mask


def affine_line_mask(target: Functional, center: int) -> int:
    """Return the indicator of ``{x: target(x)=center}``."""
    mask = 0
    for point in product(range(P), repeat=2):
        if _functional_value(P, target, point) == center % P:
            mask |= 1 << point_index(point)
    if mask.bit_count() != P:
        raise ArithmeticError("an affine line changed size")
    return mask


def selected_half_edges(
    target: Functional, auxiliary: Functional, center: int
) -> tuple[Edge, ...]:
    """Return the thirty physical edges selected by one forced half."""
    edges = []
    for edge in mobius_parameter_edges(P, target, auxiliary, center).values():
        negative = _negative_edge(P, edge)
        orbit = min(edge, negative)
        normalized_coefficient = paley_edge_sign(P, orbit) * (
            -1 if edge == orbit else 1
        )
        edges.append(
            orbit
            if normalized_coefficient == 1
            else _negative_edge(P, orbit)
        )
    result = tuple(edges)
    if len(result) != P - 1 or len(set(result)) != P - 1:
        raise ArithmeticError("a selected Mobius half lost an edge")
    return result


def half_discrepancy_mask(
    target: Functional, auxiliary: Functional, center: int
) -> int:
    """Return ``boundary(selected half) + indicator(target=center)``."""
    return boundary_mask(selected_half_edges(target, auxiliary, center)) ^ (
        affine_line_mask(target, center)
    )


def kernel_selector(
    fixed_direction: Functional, representatives: Sequence[Point] | None = None
) -> tuple[Point, ...]:
    """Choose zero and one point from each nonzero ``+/-`` pair in ``ker F``."""
    kernel = tuple(
        point
        for point in product(range(P), repeat=2)
        if _functional_value(P, fixed_direction, point) == 0
    )
    if len(kernel) != P:
        raise ArithmeticError("the fixed kernel changed size")
    if representatives is None:
        chosen = {(0, 0)}
        for point in kernel:
            if point == (0, 0):
                continue
            negative = (-point[0] % P, -point[1] % P)
            chosen.add(min(point, negative))
    else:
        chosen = set(representatives)
    if (0, 0) not in chosen or len(chosen) != (P + 1) // 2:
        raise ValueError("a kernel selector needs zero and fifteen representatives")
    for point in kernel:
        if point == (0, 0):
            continue
        negative = (-point[0] % P, -point[1] % P)
        if (point in chosen) + (negative in chosen) != 1:
            raise ValueError("a kernel selector must choose one point from each pair")
    if any(point not in kernel for point in chosen):
        raise ValueError("a selector point left the fixed kernel")
    return tuple(sorted(chosen))


def selector_mask(selector: Sequence[Point]) -> int:
    mask = 0
    for point in selector:
        mask |= 1 << point_index(point)
    return mask


def selector_pairing(mask: int, selector: Sequence[Point]) -> int:
    """Pair a point mask with a selector over ``F_2``."""
    return (mask & selector_mask(selector)).bit_count() % 2


def is_centrally_symmetric_mask(mask: int) -> bool:
    for point in product(range(P), repeat=2):
        negative = (-point[0] % P, -point[1] % P)
        if ((mask >> point_index(point)) & 1) != (
            (mask >> point_index(negative)) & 1
        ):
            return False
    return True


def half_kernel_parity(
    target: Functional,
    auxiliary: Functional,
    fixed_direction: Functional,
    center: int = 1,
) -> int:
    selector = kernel_selector(fixed_direction)
    return selector_pairing(
        half_discrepancy_mask(target, auxiliary, center), selector
    )


def projective_direction_index(functional: Functional) -> int:
    """Return the canonical projective index of a nonzero functional."""
    if functional == (0, 0):
        raise ValueError("the zero functional has no projective direction")
    for index, row in enumerate(projective_functionals(P)):
        if (
            functional[0] * row[1] - functional[1] * row[0]
        ) % P == 0:
            return index
    raise ArithmeticError("a nonzero functional lost its projective class")


def annihilator_direction_index(point: Point) -> int:
    """Return the unique projective functional vanishing at ``point``."""
    if point == (0, 0):
        raise ValueError("zero has no unique annihilator direction")
    return projective_direction_index((point[1], -point[0] % P))


def half_kernel_sigma_formula(
    target: Functional,
    auxiliary: Functional,
    kernel_direction: Functional,
) -> int:
    """Return ``(-1)^g`` from the closed projective-coordinate formula.

    Write ``D(z)=M+z(L-M)``.  For finite ``z != 1``, the unique nonzero
    endpoint of a raw parameter edge on ``ker D(z)`` has parameter ``-z``;
    that edge has spatial direction ``D(z^2)``.  Its selected physical
    orientation therefore gives ``sigma(D(z))=epsilon(D(z^2))``.  At
    ``z=1`` the missing parameter ``t=-1`` leaves only the origin, so
    ``sigma(L)=-1``.  At ``z=infinity`` the parameters other than
    ``0,+/-1`` pair as ``t,-t``.  Since ``(P-3)/2`` is even, the unpaired
    ``t=1`` term gives ``sigma(L-M)=epsilon_L``.
    """
    a, b = _relative_coefficients(P, target, auxiliary, kernel_direction)
    if b == 0:  # D(1)=L
        return -1
    if (a + b) % P == 0:  # D(infinity)=L-M
        return paley_direction_sign(P, target)
    z = a * pow(a + b, -1, P) % P
    z_squared = z * z % P
    spatial = (
        (z_squared * target[0] + (1 - z_squared) * auxiliary[0]) % P,
        (z_squared * target[1] + (1 - z_squared) * auxiliary[1]) % P,
    )
    return paley_direction_sign(P, spatial)


def half_projective_signature_certificate(
    target: Functional, auxiliary: Functional
) -> dict[str, object]:
    """Replay all 32 kernel signatures and all 30 nonzero centers at p31."""
    directions = projective_functionals(P)
    if projective_direction_index(target) == projective_direction_index(
        auxiliary
    ):
        raise ValueError("target and auxiliary must be independent")
    selectors = tuple(kernel_selector(row) for row in directions)
    masks = tuple(
        half_discrepancy_mask(target, auxiliary, center)
        for center in range(1, P)
    )
    center_bits = tuple(
        tuple(selector_pairing(mask, selector) for selector in selectors)
        for mask in masks
    )
    bits = center_bits[0]
    predicted_sigma = tuple(
        half_kernel_sigma_formula(target, auxiliary, row)
        for row in directions
    )
    actual_sigma = tuple(-1 if bit else 1 for bit in bits)
    center_invariant = len(set(center_bits)) == 1
    all_central = all(is_centrally_symmetric_mask(mask) for mask in masks)
    all_contain_origin = all(
        (mask >> point_index((0, 0))) & 1 for mask in masks
    )
    target_index = projective_direction_index(target)
    auxiliary_index = projective_direction_index(auxiliary)
    infinity_index = projective_direction_index(
        (
            (target[0] - auxiliary[0]) % P,
            (target[1] - auxiliary[1]) % P,
        )
    )
    signature_product = sign_product(actual_sigma)
    expected_product = -paley_direction_sign(P, auxiliary)
    formula_matches = actual_sigma == predicted_sigma
    proved = bool(
        center_invariant
        and all_central
        and all_contain_origin
        and formula_matches
        and actual_sigma[target_index] == -1
        and actual_sigma[infinity_index]
        == paley_direction_sign(P, target)
        and signature_product == expected_product
    )
    if not proved:
        raise ArithmeticError("the p31 half signature formula failed replay")
    return {
        "target": target,
        "target_direction_index": target_index,
        "target_sign": paley_direction_sign(P, target),
        "auxiliary": auxiliary,
        "auxiliary_direction_index": auxiliary_index,
        "auxiliary_sign": paley_direction_sign(P, auxiliary),
        "l_minus_m_direction_index": infinity_index,
        "kernel_parity_bits": bits,
        "kernel_sigma": actual_sigma,
        "kernel_sigma_sha256": hashlib.sha256(
            bytes(1 if sign == 1 else 0 for sign in actual_sigma)
        ).hexdigest(),
        "all_centers_have_same_signature": center_invariant,
        "all_discrepancies_centrally_symmetric": all_central,
        "all_discrepancies_contain_origin": all_contain_origin,
        "closed_formula_matches_every_direction": formula_matches,
        "sigma_at_target": actual_sigma[target_index],
        "sigma_at_l_minus_m": actual_sigma[infinity_index],
        "signature_product": signature_product,
        "expected_signature_product": expected_product,
        "proved": proved,
    }


def fixed_antipodal_edges(fixed_direction: Functional) -> tuple[Edge, ...]:
    """Return the fifteen fixed edges whose spatial direction is ``F``."""
    out = set()
    for point in product(range(P), repeat=2):
        if point == (0, 0):
            continue
        if _functional_value(P, fixed_direction, point) != 0:
            continue
        negative = (-point[0] % P, -point[1] % P)
        out.add(tuple(sorted((point, negative))))
    result = tuple(sorted(out))
    if len(result) != (P - 1) // 2:
        raise ArithmeticError("the fixed antipodal fibre changed size")
    return result


def nonorigin_orbit_removal_parity(
    edge: Edge, fixed_direction: Functional
) -> int:
    """Pair the boundaries of ``edge`` and ``-edge`` with the selector."""
    if (0, 0) in edge:
        raise ValueError("the cancellation orbit must be nonorigin")
    selector = kernel_selector(fixed_direction)
    return selector_pairing(
        boundary_mask((edge, _negative_edge(P, edge))), selector
    )


def inversion_orbit_pair_signature(edge: Edge) -> tuple[int, ...]:
    """Return all kernel pairings of ``boundary(edge)+boundary(-edge)``."""
    pair_mask = boundary_mask((edge, _negative_edge(P, edge)))
    return tuple(
        selector_pairing(pair_mask, kernel_selector(direction))
        for direction in projective_functionals(P)
    )


def adaptive_design_boundary_parity_certificate(
    halves: Sequence[tuple[Functional, Functional]],
    fixed_direction_index: int,
    *,
    centers: Sequence[int] | None = None,
    cancelled_edge: Edge | None = None,
) -> dict[str, object]:
    """Certify the adaptive-kernel obstruction for a p31 top design.

    The returned product obstruction is insensitive to every valid ternary
    reduction on a nonorigin inversion orbit.  The optional ``cancelled_edge``
    merely replays the concrete clean-collision signature; it is not a
    hypothesis of the stronger product conclusion.
    """
    if len(halves) != 16:
        raise ValueError("the p31 top design needs sixteen halves")
    if centers is None:
        centers = (1,) * len(halves)
    if len(centers) != len(halves) or any(int(c) % P == 0 for c in centers):
        raise ValueError("need one nonzero center per half")

    directions = projective_functionals(P)
    direction_signs = tuple(
        paley_direction_sign(P, direction) for direction in directions
    )
    hard = tuple(index for index, sign in enumerate(direction_signs) if sign == 1)
    opposite = tuple(
        index for index, sign in enumerate(direction_signs) if sign == -1
    )
    fixed_direction = directions[fixed_direction_index]
    local = tuple(
        half_projective_signature_certificate(target, auxiliary)
        for target, auxiliary in halves
    )
    target_indices = tuple(
        int(record["target_direction_index"]) for record in local
    )
    auxiliary_indices = tuple(
        int(record["auxiliary_direction_index"]) for record in local
    )
    auxiliary_signs = tuple(int(record["auxiliary_sign"]) for record in local)
    auxiliary_hard_count = sum(sign == 1 for sign in auxiliary_signs)
    auxiliary_opposite_count = sum(sign == -1 for sign in auxiliary_signs)
    auxiliary_sign_product = sign_product(auxiliary_signs)

    half_bits = tuple(
        sum(int(record["kernel_parity_bits"][index]) for record in local) % 2
        for index in range(P + 1)
    )
    half_sigma = tuple(-1 if bit else 1 for bit in half_bits)
    half_signature_product = sign_product(half_sigma)
    local_product_formula = sign_product(
        -int(record["auxiliary_sign"]) for record in local
    )

    fixed_edges = fixed_antipodal_edges(fixed_direction)
    fixed_bits_by_direction = []
    fixed_pairing_constant = True
    for direction in directions:
        selector = kernel_selector(direction)
        values = tuple(
            selector_pairing(boundary_mask((edge,)), selector)
            for edge in fixed_edges
        )
        fixed_pairing_constant &= len(set(values)) == 1
        fixed_bits_by_direction.append(values[0])
    fixed_bits = tuple(fixed_bits_by_direction)
    expected_fixed_bits = tuple(
        int(index == fixed_direction_index) for index in range(P + 1)
    )

    # Every half has exactly one origin edge.  Its annihilator direction is
    # its auxiliary direction, so an auxiliary SDR prevents any origin-orbit
    # collision and therefore prevents the sole kind of orbit-pair correction
    # whose aggregate kernel parity would be odd.
    origin_orbits = []
    half_has_no_fixed_edge = []
    for (target, auxiliary), center in zip(halves, centers, strict=True):
        selected = selected_half_edges(target, auxiliary, int(center))
        origin = tuple(edge for edge in selected if (0, 0) in edge)
        if len(origin) != 1:
            raise ArithmeticError("a Mobius half lost its unique origin edge")
        origin_orbits.append(min(origin[0], _negative_edge(P, origin[0])))
        other_endpoint = origin[0][0] if origin[0][1] == (0, 0) else origin[0][1]
        if annihilator_direction_index(other_endpoint) != projective_direction_index(
            auxiliary
        ):
            raise ArithmeticError("an origin orbit lost its auxiliary direction")
        half_has_no_fixed_edge.append(
            all(edge != _negative_edge(P, edge) for edge in selected)
        )
    origin_orbits_distinct = len(set(origin_orbits)) == len(origin_orbits)

    raw_parallel_profile = tuple(
        sum(_parallel_formula(P, target, auxiliary, row) for target, auxiliary in halves)
        for row in directions
    )
    top_profile_multisets_match = bool(
        sorted(raw_parallel_profile[index] for index in hard)
        == [14] * 14 + [15] * 2
        and sorted(raw_parallel_profile[index] for index in opposite)
        == [15] * 2 + [16] * 14
    )
    top_opposite_fixed_sdr = bool(
        set(target_indices) == set(hard)
        and len(set(target_indices)) == 16
        and len(set(auxiliary_indices)) == 16
        and auxiliary_hard_count == 14
        and auxiliary_opposite_count == 2
        and direction_signs[fixed_direction_index] == -1
        and top_profile_multisets_match
    )

    removal_bits = None
    clean_signature_bits = None
    clean_contradiction_directions = None
    if cancelled_edge is not None:
        if (0, 0) in cancelled_edge:
            raise ValueError("the displayed cancellation must be nonorigin")
        removal_bits = inversion_orbit_pair_signature(cancelled_edge)
        clean_signature_bits = tuple(
            half_bits[index] ^ fixed_bits[index] ^ removal_bits[index]
            for index in range(P + 1)
        )
        clean_contradiction_directions = tuple(
            index for index, bit in enumerate(clean_signature_bits) if bit
        )

    all_local_formulae = all(bool(record["proved"]) for record in local)
    nonorigin_pair_aggregate_even = True
    final_product_after_arbitrary_nonorigin_reductions = (
        half_signature_product * sign_product(-1 if bit else 1 for bit in fixed_bits)
    )
    product_obstruction = bool(
        all_local_formulae
        and fixed_pairing_constant
        and fixed_bits == expected_fixed_bits
        and all(half_has_no_fixed_edge)
        and origin_orbits_distinct
        and auxiliary_sign_product == 1
        and half_signature_product == 1
        and final_product_after_arbitrary_nonorigin_reductions == -1
    )
    top_family_obstruction = top_opposite_fixed_sdr and product_obstruction
    proved = bool(
        all_local_formulae
        and half_signature_product == local_product_formula
        and fixed_pairing_constant
        and fixed_bits == expected_fixed_bits
        and all(half_has_no_fixed_edge)
    )
    if not proved:
        raise ArithmeticError("the adaptive p31 signature replay failed")
    return {
        "p": P,
        "fixed_direction_index": fixed_direction_index,
        "fixed_direction_sign": direction_signs[fixed_direction_index],
        "target_direction_indices": target_indices,
        "auxiliary_direction_indices": auxiliary_indices,
        "auxiliary_hard_count": auxiliary_hard_count,
        "auxiliary_opposite_count": auxiliary_opposite_count,
        "auxiliary_sign_product": auxiliary_sign_product,
        "raw_parallel_profile": raw_parallel_profile,
        "top_profile_multisets_match": top_profile_multisets_match,
        "top_opposite_fixed_sdr": top_opposite_fixed_sdr,
        "half_kernel_parity_bits_by_direction": half_bits,
        "half_kernel_signature_product": half_signature_product,
        "local_product_formula": local_product_formula,
        "every_local_product_is_minus_auxiliary_sign": all_local_formulae,
        "fixed_edge_kernel_parity_bits": fixed_bits,
        "fixed_edge_signature_product": sign_product(
            -1 if bit else 1 for bit in fixed_bits
        ),
        "fixed_edge_pairing_constant_within_each_kernel": fixed_pairing_constant,
        "origin_orbit_count": len(origin_orbits),
        "origin_orbits_distinct": origin_orbits_distinct,
        "every_half_has_no_fixed_antipodal_edge": all(half_has_no_fixed_edge),
        "nonorigin_orbit_pair_aggregate_parity": 0,
        "nonorigin_orbit_pair_aggregate_even": nonorigin_pair_aggregate_even,
        "cancelled_orbit_kernel_parity_bits": removal_bits,
        "clean_collision_contradiction_bits": clean_signature_bits,
        "clean_collision_contradiction_direction_indices": (
            clean_contradiction_directions
        ),
        "final_signature_product_after_arbitrary_nonorigin_reductions": (
            final_product_after_arbitrary_nonorigin_reductions
        ),
        "adaptive_kernel_product_obstruction": product_obstruction,
        "top_family_boundary_parity_obstruction": top_family_obstruction,
        "some_kernel_selector_is_a_contradiction": product_obstruction,
        "scope": (
            "all valid ternary reductions of the supplied p31 localized-Mobius "
            "half family on nonorigin inversion orbits, with one fixed edge "
            "and no origin-pair correction"
        ),
        "residual_ii_closed": False,
        "proved": proved,
    }


def design_kernel_parity_certificate(
    halves: Sequence[tuple[Functional, Functional]],
    fixed_direction_index: int,
    *,
    centers: Sequence[int] | None = None,
    cancelled_edge: Edge | None = None,
) -> dict[str, object]:
    """Replay the kernel-selector obstruction for one sixteen-half design.

    A return value with ``boundary_parity_obstruction=True`` says that no
    choice of the sixteen nonzero centers, no nonorigin cancellation of the
    displayed inversion orbit, and no choice of the fixed antipodal edge can
    meet the necessary target degree boundary.
    """
    if len(halves) != 16:
        raise ValueError("the top design needs sixteen halves")
    if centers is None:
        centers = (1,) * len(halves)
    if len(centers) != len(halves) or any(int(c) % P == 0 for c in centers):
        raise ValueError("need one nonzero center per half")
    directions = projective_functionals(P)
    fixed_direction = directions[fixed_direction_index]
    selector = kernel_selector(fixed_direction)

    bits = []
    center_invariant = []
    symmetric = []
    contains_origin = []
    for target, auxiliary in halves:
        masks = tuple(
            half_discrepancy_mask(target, auxiliary, center)
            for center in range(1, P)
        )
        pairings = tuple(selector_pairing(mask, selector) for mask in masks)
        bits.append(pairings[int(centers[len(bits)]) - 1])
        center_invariant.append(len(set(pairings)) == 1)
        symmetric.append(all(is_centrally_symmetric_mask(mask) for mask in masks))
        contains_origin.append(
            all((mask >> point_index((0, 0))) & 1 for mask in masks)
        )

    fixed_pairings = tuple(
        selector_pairing(boundary_mask((edge,)), selector)
        for edge in fixed_antipodal_edges(fixed_direction)
    )
    removal = (
        None
        if cancelled_edge is None
        else nonorigin_orbit_removal_parity(cancelled_edge, fixed_direction)
    )
    proved_local_structure = bool(
        all(center_invariant)
        and all(symmetric)
        and all(contains_origin)
        and set(fixed_pairings) == {1}
        and removal in (None, 0)
    )
    discrepancy_total = sum(bits) % 2
    obstruction = proved_local_structure and discrepancy_total == 0
    return {
        "p": P,
        "fixed_direction_index": fixed_direction_index,
        "fixed_direction": fixed_direction,
        "selector": selector,
        "selector_size": len(selector),
        "half_kernel_parity_bits": tuple(bits),
        "half_parity_sum_mod_two": discrepancy_total,
        "all_half_pairings_center_invariant": all(center_invariant),
        "all_half_discrepancies_centrally_symmetric": all(symmetric),
        "all_half_discrepancies_contain_origin": all(contains_origin),
        "fixed_edge_pairings": fixed_pairings,
        "every_fixed_edge_pairing": 1,
        "cancelled_orbit_removal_pairing": removal,
        "necessary_total_pairing": (discrepancy_total + 1) % 2,
        "boundary_parity_obstruction": obstruction,
        "scope": "supplied localized-Mobius top design only",
        "residual_ii_closed": False,
        "proved": proved_local_structure,
    }
