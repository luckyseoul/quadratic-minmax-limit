#!/usr/bin/env python3
"""Projective collision keys for the p=31 localized-Mobius halves.

For independent functionals ``L,M`` and a nonorigin Mobius parameter
``t != 0,-1``, the unit-centre edge has endpoint-annihilator directions

``U=[tL-(t+1)M]`` and ``V=[L-M]``

and spatial direction

``D=[t^2 L+(1-t^2)M]``.

The triple ``(D,{U,V})`` is independent of the centre.  Two such skeleton
edges can be centred onto the same inversion orbit if and only if their
triples agree.  This module records that exact algebra and replays the
hard-fixed sixteen-half witness for which a brute-force 108,000 centre-pair
audit had found no cancellation in the required direction.

This is a physical-collision prefilter for one localized-Mobius family.  It
does not assert that the hard-fixed branch, branch C, or residual (ii) is
closed.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from itertools import combinations
from typing import Iterable

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    Point,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (
    _parallel_formula,
    mobius_parameter_edges,
    paley_direction_sign,
)
from e1_gmin_m4_p31_top_mobius_boundary_parity import (
    half_kernel_sigma_formula,
)
from e1_gmin_m4_prop15721 import is_prime


P = 31
DIRECTIONS = tuple(projective_functionals(P))
FIXED_DIRECTION_INDEX = 0
REQUIRED_CANCELLATION_DIRECTION_INDEX = 0
FROZEN_CORRECTION_SUPPORT = (2, 23)
SOURCE_WITNESS_SHA256 = (
    "9e241aeefee9d91132833b83d0d2e6cf2e4f875b531ddd039551fafb0b5cea91"
)
OPTION_CATALOG_SHA256 = (
    "cec2298fc950ed6a97a60b203756c6caddcb15f239cc78e6ab2de79b6cdea5e0"
)
FROZEN_LABELLED_KEY_CATALOG_SHA256 = (
    "f3512567d9655427212ec0383b6d1e0f2f15b364031fa348980a44cb5724d80b"
)
FULL_HARD_FIXED_KEY_MEMBERSHIP_SHA256 = (
    "2775cd4ad86a1844834dba1981a8396a92fb2f6da9b840a739e9099abef32897"
)
FROZEN_HALF_CHOICES = (
    (0, 22, 12),
    (1, 27, 27),
    (2, 15, 11),
    (3, 28, 7),
    (7, 16, 9),
    (9, 4, 5),
    (10, 7, 20),
    (15, 26, 16),
    (16, 24, 28),
    (21, 29, 15),
    (22, 9, 19),
    (24, 10, 15),
    (28, 21, 21),
    (29, 31, 27),
    (30, 2, 1),
    (31, 3, 25),
)
EXPECTED_RAW_PROFILE = (
    15, 15, 14, 14, 15, 16, 16, 14,
    16, 14, 14, 16, 16, 16, 16, 14,
    14, 16, 16, 16, 16, 14, 14, 16,
    14, 16, 15, 15, 14, 14, 15, 14,
)
EXPECTED_AGGREGATE_SIGNATURE = 0x00800005
EXPECTED_CORRECTION_SIGNATURE = 0x00800004
EXPECTED_SHARED_KEYS = (
    (4, 3, 9),
    (9, 16, 31),
    (20, 28, 29),
    (21, 27, 31),
    (29, 7, 16),
)


@dataclass(frozen=True, order=True)
class CollisionKey:
    """Centre-free class of a nonorigin Mobius skeleton edge."""

    spatial_direction: int
    endpoint_directions: tuple[int, int]

    def __post_init__(self) -> None:
        first, second = self.endpoint_directions
        if first >= second:
            raise ValueError("endpoint directions must be sorted and distinct")
        if self.spatial_direction in self.endpoint_directions:
            raise ValueError("spatial and endpoint directions must be distinct")

    def as_tuple(self) -> tuple[int, int, int]:
        return (self.spatial_direction, *self.endpoint_directions)


def determinant(p: int, first: Functional, second: Functional) -> int:
    """Return the determinant of two row functionals modulo ``p``."""
    return (first[0] * second[1] - first[1] * second[0]) % p


def _check_paley_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=3 mod 4")


def _check_independent(p: int, target: Functional, auxiliary: Functional) -> None:
    if determinant(p, target, auxiliary) == 0:
        raise ValueError("target and auxiliary functionals must be independent")


def projective_direction_index(p: int, functional: Functional) -> int:
    """Return the deterministic projective index of a nonzero functional."""
    functional = functional[0] % p, functional[1] % p
    if functional == (0, 0):
        raise ValueError("zero has no projective direction")
    for index, row in enumerate(projective_functionals(p)):
        if determinant(p, functional, row) == 0:
            return index
    raise ArithmeticError("a nonzero functional lost its projective class")


def annihilator_direction_index(p: int, point: Point) -> int:
    """Return the projective functional direction annihilating ``point``."""
    point = point[0] % p, point[1] % p
    if point == (0, 0):
        raise ValueError("the origin has no unique annihilator direction")
    return projective_direction_index(p, (point[1], -point[0] % p))


def _scaled_auxiliary(
    p: int, auxiliary_index: int, relative_scale: int
) -> Functional:
    row = projective_functionals(p)[auxiliary_index]
    return (
        relative_scale * row[0] % p,
        relative_scale * row[1] % p,
    )


def nonorigin_collision_key(
    p: int,
    target: Functional,
    auxiliary: Functional,
    parameter: int,
) -> CollisionKey:
    """Return ``(D,{U,V})`` from the closed Mobius formulas.

    Here ``U=[tL-(t+1)M]``, ``V=[L-M]``, and
    ``D=[t^2 L+(1-t^2)M]``.  The determinant identities

    ``det(U,V)=det(L,M)``, ``det(V,D)=det(L,M)``, and
    ``det(U,D)=t(t+1)det(L,M)``

    prove that all three directions are distinct for ``t != 0,-1``.
    """
    _check_paley_prime(p)
    _check_independent(p, target, auxiliary)
    parameter %= p
    if parameter in (0, p - 1):
        raise ValueError("a nonorigin collision parameter must avoid 0 and -1")
    t2 = parameter * parameter % p
    endpoint_u = (
        parameter * target[0] - (parameter + 1) * auxiliary[0],
        parameter * target[1] - (parameter + 1) * auxiliary[1],
    )
    endpoint_v = (
        target[0] - auxiliary[0],
        target[1] - auxiliary[1],
    )
    spatial = (
        t2 * target[0] + (1 - t2) * auxiliary[0],
        t2 * target[1] + (1 - t2) * auxiliary[1],
    )
    endpoint_u = endpoint_u[0] % p, endpoint_u[1] % p
    endpoint_v = endpoint_v[0] % p, endpoint_v[1] % p
    spatial = spatial[0] % p, spatial[1] % p

    delta = determinant(p, target, auxiliary)
    if (
        determinant(p, endpoint_u, endpoint_v) != delta
        or determinant(p, endpoint_v, spatial) != delta
        or determinant(p, endpoint_u, spatial)
        != parameter * (parameter + 1) * delta % p
    ):
        raise ArithmeticError("the collision-key determinant identities changed")

    endpoint_indices = tuple(
        sorted(
            (
                projective_direction_index(p, endpoint_u),
                projective_direction_index(p, endpoint_v),
            )
        )
    )
    key = CollisionKey(
        projective_direction_index(p, spatial), endpoint_indices
    )
    if key.spatial_direction in key.endpoint_directions:
        raise ArithmeticError("a nonorigin collision key lost transversality")
    return key


def physical_edge_collision_key(p: int, edge: Edge) -> CollisionKey:
    """Read the same collision key directly from a physical nonorigin edge."""
    first, second = edge
    if first == (0, 0) or second == (0, 0):
        raise ValueError("the origin edge has only one endpoint direction")
    endpoint_directions = tuple(
        sorted(
            (
                annihilator_direction_index(p, first),
                annihilator_direction_index(p, second),
            )
        )
    )
    difference = (
        (first[0] - second[0]) % p,
        (first[1] - second[1]) % p,
    )
    return CollisionKey(
        annihilator_direction_index(p, difference), endpoint_directions
    )


def mobius_half_collision_keys(
    p: int, target: Functional, auxiliary: Functional
) -> dict[int, CollisionKey]:
    """Return all ``p-2`` nonorigin keys of one unit-centre half."""
    edges = mobius_parameter_edges(p, target, auxiliary, center=1)
    out: dict[int, CollisionKey] = {}
    for parameter, edge in edges.items():
        if parameter == 0:
            continue
        formula = nonorigin_collision_key(p, target, auxiliary, parameter)
        physical = physical_edge_collision_key(p, edge)
        if formula != physical:
            raise ArithmeticError("the formula and physical collision keys disagree")
        out[parameter] = formula
    if len(out) != p - 2:
        raise ArithmeticError("one half lost a nonorigin collision parameter")
    return out


def _point_on_endpoint_direction(
    p: int, edge: Edge
) -> dict[int, Point]:
    out = {
        annihilator_direction_index(p, point): point
        for point in edge
    }
    if len(out) != 2:
        raise ArithmeticError("a nonorigin edge lost independent endpoint rays")
    return out


def _ray_scale(p: int, first: Point, second: Point) -> int:
    """Return the nonzero ``scale`` for ``second = scale*first``."""
    if first[0] % p:
        scale = second[0] * pow(first[0], -1, p) % p
    elif first[1] % p:
        scale = second[1] * pow(first[1], -1, p) % p
    else:
        raise ValueError("the origin has no ray scale")
    if scale == 0 or (
        scale * first[0] % p,
        scale * first[1] % p,
    ) != (second[0] % p, second[1] % p):
        raise ValueError("the points are not on one nonzero ray")
    return scale


def collision_homothety_scalar(
    p: int,
    first_target: Functional,
    first_auxiliary: Functional,
    first_parameter: int,
    second_target: Functional,
    second_auxiliary: Functional,
    second_parameter: int,
) -> int | None:
    """Return ``lambda`` with ``E2(1)=lambda E1(1)``, or ``None``.

    This is the necessary-and-sufficient centre theorem.  Equality of the
    two collision keys matches the two independent endpoint rays.  If their
    ray scales are ``a,b``, equality of the spatial direction implies
    ``a*x-b*y`` is parallel to ``x-y``; independence gives ``a=b``.  Thus
    equal keys are equivalent to a global homothety of the unit-centre
    edges.  If the returned scalar is ``lambda``, centred edges lie on the
    same inversion orbit exactly when

    ``c1 = lambda*c2`` or ``c1 = -lambda*c2``.

    The second equation makes the two physical edges negatives and hence
    gives opposite normalized antisymmetric coefficients (cancellation).
    """
    first_key = nonorigin_collision_key(
        p, first_target, first_auxiliary, first_parameter
    )
    second_key = nonorigin_collision_key(
        p, second_target, second_auxiliary, second_parameter
    )
    if first_key != second_key:
        return None

    first_edge = mobius_parameter_edges(
        p, first_target, first_auxiliary, center=1
    )[first_parameter % p]
    second_edge = mobius_parameter_edges(
        p, second_target, second_auxiliary, center=1
    )[second_parameter % p]
    first_points = _point_on_endpoint_direction(p, first_edge)
    second_points = _point_on_endpoint_direction(p, second_edge)
    scales = {
        _ray_scale(p, first_points[index], second_points[index])
        for index in first_key.endpoint_directions
    }
    if len(scales) != 1:
        raise ArithmeticError("equal collision keys were not homothetic")
    return scales.pop()


def centred_edges_share_orbit(
    p: int,
    first_target: Functional,
    first_auxiliary: Functional,
    first_parameter: int,
    first_center: int,
    second_target: Functional,
    second_auxiliary: Functional,
    second_parameter: int,
    second_center: int,
) -> dict[str, object]:
    """Apply the exact centre relation to two nonorigin skeleton edges."""
    first_center %= p
    second_center %= p
    if first_center == 0 or second_center == 0:
        raise ValueError("Mobius half centres must be nonzero")
    scale = collision_homothety_scalar(
        p,
        first_target,
        first_auxiliary,
        first_parameter,
        second_target,
        second_auxiliary,
        second_parameter,
    )
    if scale is None:
        return {
            "keys_equal": False,
            "homothety_scalar": None,
            "same_inversion_orbit": False,
            "opposite_physical_edges": False,
            "cancellation_orientation": False,
        }
    same = first_center == scale * second_center % p
    opposite = first_center == -scale * second_center % p
    return {
        "keys_equal": True,
        "homothety_scalar": scale,
        "same_inversion_orbit": same or opposite,
        "opposite_physical_edges": opposite,
        "cancellation_orientation": opposite,
    }


def prescribed_endpoint_candidate(
    p: int,
    target: Functional,
    auxiliary: Functional,
    endpoint_directions: Iterable[int],
) -> tuple[int, CollisionKey] | None:
    """Return the unique parameter/key compatible with a prescribed pair.

    One endpoint direction of every edge in the half is fixed as
    ``V=[L-M]``.  If the requested pair is ``{V,H}``, then the other endpoint
    equation has the unique solution

    ``t = det(M,H)/det(L-M,H)``.

    Therefore a fixed two-bit boundary correction leaves at most one
    collision parameter, and hence at most one spatial-direction bucket, in
    each labelled half.
    """
    endpoints = tuple(sorted(set(endpoint_directions)))
    if len(endpoints) != 2 or any(not 0 <= index <= p for index in endpoints):
        raise ValueError("need two distinct projective endpoint indices")
    _check_independent(p, target, auxiliary)
    difference = (
        (target[0] - auxiliary[0]) % p,
        (target[1] - auxiliary[1]) % p,
    )
    fixed_endpoint = projective_direction_index(p, difference)
    if fixed_endpoint not in endpoints:
        return None
    other_endpoint = endpoints[0] if endpoints[1] == fixed_endpoint else endpoints[1]
    other = projective_functionals(p)[other_endpoint]
    denominator = determinant(p, difference, other)
    if denominator == 0:
        raise ArithmeticError("distinct endpoint directions became parallel")
    parameter = determinant(p, auxiliary, other) * pow(denominator, -1, p) % p
    if parameter in (0, p - 1):
        return None
    key = nonorigin_collision_key(p, target, auxiliary, parameter)
    if key.endpoint_directions != endpoints:
        raise ArithmeticError("the prescribed-endpoint recovery changed")
    return parameter, key


def target_options_for_collision_key(
    p: int,
    target_index: int,
    key: CollisionKey,
) -> tuple[tuple[int, int, int, int], ...]:
    """Invert one key at a fixed canonical target in at most two ways.

    Each returned row is
    ``(auxiliary_index, relative_scale, parameter, fixed_endpoint)``.
    The two possible rows correspond to deciding which member of the
    unordered endpoint pair is ``V=[L-M]``.

    For an ordered assignment ``(K,H)=(V,U)``, express ``L`` in the basis
    ``K,H`` and compare ``D`` with

    ``D = L+(t^2-1)(L-M) = t*l_K*K+l_H*H``.

    This gives the projectively invariant formula for ``t``

    ``t = det(D,H)det(K,L)/(det(K,D)det(L,H))``

    and, for the chosen representative ``K``, ``L-M=sK`` with

    ``s = det(L,H)/(det(K,H)(t+1))``.

    Hence every fixed target has at most two labelled options containing a
    prescribed collision key.  This is the sharp variable reduction used by
    a key-sharded simultaneous solver.
    """
    _check_paley_prime(p)
    directions = projective_functionals(p)
    if not 0 <= target_index < len(directions):
        raise ValueError("target index left the projective line")
    if (
        not 0 <= key.spatial_direction < len(directions)
        or any(not 0 <= index < len(directions) for index in key.endpoint_directions)
    ):
        raise ValueError("collision key left the projective line")
    target = directions[target_index]
    spatial = directions[key.spatial_direction]
    out = []
    for fixed_endpoint in key.endpoint_directions:
        moving_endpoint = (
            key.endpoint_directions[0]
            if fixed_endpoint == key.endpoint_directions[1]
            else key.endpoint_directions[1]
        )
        fixed = directions[fixed_endpoint]
        moving = directions[moving_endpoint]
        # A valid nonorigin edge has L transverse to both endpoint kernels.
        if determinant(p, target, fixed) == 0 or determinant(p, target, moving) == 0:
            continue
        denominator = (
            determinant(p, fixed, spatial)
            * determinant(p, target, moving)
        ) % p
        if denominator == 0:
            raise ArithmeticError("a valid collision key lost an inverse denominator")
        parameter = (
            determinant(p, spatial, moving)
            * determinant(p, fixed, target)
            * pow(denominator, -1, p)
        ) % p
        if parameter in (0, p - 1):
            continue
        s = (
            determinant(p, target, moving)
            * pow(determinant(p, fixed, moving), -1, p)
            * pow(parameter + 1, -1, p)
        ) % p
        if s == 0:
            raise ArithmeticError("the reconstructed L-M scale vanished")
        auxiliary = (
            (target[0] - s * fixed[0]) % p,
            (target[1] - s * fixed[1]) % p,
        )
        _check_independent(p, target, auxiliary)
        auxiliary_index = projective_direction_index(p, auxiliary)
        canonical_auxiliary = directions[auxiliary_index]
        relative_scale = (
            auxiliary[0] if canonical_auxiliary[0] else auxiliary[1]
        )
        if relative_scale == 0 or _scaled_auxiliary(
            p, auxiliary_index, relative_scale
        ) != auxiliary:
            raise ArithmeticError("the reconstructed auxiliary scale changed")
        if nonorigin_collision_key(p, target, auxiliary, parameter) != key:
            raise ArithmeticError("the inverse collision-key formula changed")
        out.append(
            (
                auxiliary_index,
                relative_scale,
                parameter,
                fixed_endpoint,
            )
        )
    result = tuple(sorted(set(out)))
    if len(result) > 2:
        raise ArithmeticError("one target acquired more than two key options")
    return result


def full_p31_inverse_catalog_replay() -> dict[str, object]:
    """Exhaustively compare the forward and inverse hard-fixed key catalogs.

    This optional, roughly 24-second audit covers all 14,430 labelled choices
    allowed when fixed auxiliary direction zero is omitted, and all 29
    nonorigin parameters of each choice.  It is kept out of the default unit
    test so the focused replay remains fast.
    """
    signs = tuple(paley_direction_sign(P, row) for row in DIRECTIONS)
    hard = tuple(index for index, sign in enumerate(signs) if sign == 1)
    forward = set()
    for target_index in hard:
        target = DIRECTIONS[target_index]
        for auxiliary_index, _auxiliary_row in enumerate(DIRECTIONS):
            if auxiliary_index in (target_index, FIXED_DIRECTION_INDEX):
                continue
            for relative_scale in range(1, P):
                auxiliary = _scaled_auxiliary(
                    P, auxiliary_index, relative_scale
                )
                for parameter in range(1, P - 1):
                    key = nonorigin_collision_key(
                        P, target, auxiliary, parameter
                    )
                    forward.add(
                        (
                            target_index,
                            auxiliary_index,
                            relative_scale,
                            parameter,
                            *key.as_tuple(),
                        )
                    )

    inverse = set()
    for spatial_direction in range(P + 1):
        for first_endpoint in range(P + 1):
            for second_endpoint in range(first_endpoint + 1, P + 1):
                if spatial_direction in (first_endpoint, second_endpoint):
                    continue
                key = CollisionKey(
                    spatial_direction, (first_endpoint, second_endpoint)
                )
                for target_index in hard:
                    for (
                        auxiliary_index,
                        relative_scale,
                        parameter,
                        _fixed_endpoint,
                    ) in target_options_for_collision_key(P, target_index, key):
                        if auxiliary_index in (
                            target_index,
                            FIXED_DIRECTION_INDEX,
                        ):
                            continue
                        inverse.add(
                            (
                                target_index,
                                auxiliary_index,
                                relative_scale,
                                parameter,
                                *key.as_tuple(),
                            )
                        )

    ordered_forward = sorted(forward)
    digest = hashlib.sha256(
        json.dumps(ordered_forward, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    proved = bool(
        len(forward) == 14_430 * 29
        and forward == inverse
        and digest == FULL_HARD_FIXED_KEY_MEMBERSHIP_SHA256
    )
    if not proved:
        raise ArithmeticError("the forward/inverse collision catalog changed")
    return {
        "schema": "p31_hard_fixed_collision_key_inverse_catalog_v1",
        "labelled_choice_count": 14_430,
        "nonorigin_parameters_per_choice": 29,
        "forward_membership_count": len(forward),
        "inverse_membership_count": len(inverse),
        "forward_minus_inverse_count": len(forward - inverse),
        "inverse_minus_forward_count": len(inverse - forward),
        "catalogs_equal": True,
        "membership_catalog_sha256": digest,
        "centre_enumeration_used": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def _frozen_design() -> tuple[tuple[int, int, int, Functional, Functional], ...]:
    out = []
    for target_index, auxiliary_index, relative_scale in FROZEN_HALF_CHOICES:
        out.append(
            (
                target_index,
                auxiliary_index,
                relative_scale,
                DIRECTIONS[target_index],
                _scaled_auxiliary(P, auxiliary_index, relative_scale),
            )
        )
    return tuple(out)


def frozen_witness_collision_key_replay() -> dict[str, object]:
    """Replay the frozen hard-fixed witness by a centre-free class join."""
    design = _frozen_design()
    signs = tuple(paley_direction_sign(P, row) for row in DIRECTIONS)
    hard = tuple(index for index, sign in enumerate(signs) if sign == 1)
    raw_profile = tuple(
        sum(_parallel_formula(P, target, auxiliary, row) for *_, target, auxiliary in design)
        for row in DIRECTIONS
    )
    aggregate_signature = 0
    for *_, target, auxiliary in design:
        for kernel_index, kernel in enumerate(DIRECTIONS):
            if half_kernel_sigma_formula(target, auxiliary, kernel) == -1:
                aggregate_signature ^= 1 << kernel_index
    correction_signature = aggregate_signature ^ (1 << FIXED_DIRECTION_INDEX)

    owners: defaultdict[CollisionKey, list[tuple[int, int]]] = defaultdict(list)
    catalog_rows = []
    for half_index, (*_, target, auxiliary) in enumerate(design):
        for parameter, key in mobius_half_collision_keys(P, target, auxiliary).items():
            owners[key].append((half_index, parameter))
            catalog_rows.append((half_index, parameter, *key.as_tuple()))
    catalog_sha256 = hashlib.sha256(
        json.dumps(catalog_rows, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    shared_records = []
    pair_class_count = 0
    for key in sorted(owners):
        entries = sorted(owners[key])
        distinct_halves = {half_index for half_index, _ in entries}
        if len(distinct_halves) < 2:
            continue
        if len(entries) != len(distinct_halves):
            raise ArithmeticError("one half repeated a collision key")
        pair_class_count += len(tuple(combinations(entries, 2)))
        shared_records.append(
            {
                "key": key.as_tuple(),
                "owners_half_index_parameter": tuple(entries),
            }
        )

    prescribed_records = []
    for half_index, (
        target_index,
        auxiliary_index,
        relative_scale,
        target,
        auxiliary,
    ) in enumerate(design):
        candidate = prescribed_endpoint_candidate(
            P, target, auxiliary, FROZEN_CORRECTION_SUPPORT
        )
        if candidate is None:
            continue
        parameter, key = candidate
        prescribed_records.append(
            {
                "half_index": half_index,
                "target_index": target_index,
                "auxiliary_index": auxiliary_index,
                "relative_scale": relative_scale,
                "parameter": parameter,
                "key": key.as_tuple(),
            }
        )

    expected_prescribed = (
        {
            "half_index": 13,
            "target_index": 29,
            "auxiliary_index": 31,
            "relative_scale": 27,
            "parameter": 12,
            "key": (8, 2, 23),
        },
    )
    shared_keys = tuple(record["key"] for record in shared_records)
    auxiliary_indices = tuple(row[1] for row in design)
    input_replayed = bool(
        tuple(row[0] for row in design) == hard
        and len(set(auxiliary_indices)) == len(auxiliary_indices)
        and raw_profile == EXPECTED_RAW_PROFILE
        and aggregate_signature == EXPECTED_AGGREGATE_SIGNATURE
        and correction_signature == EXPECTED_CORRECTION_SIGNATURE
        and tuple(
            index for index in range(P + 1) if correction_signature >> index & 1
        )
        == FROZEN_CORRECTION_SUPPORT
    )
    frozen_family_excluded = bool(
        input_replayed
        and catalog_sha256 == FROZEN_LABELLED_KEY_CATALOG_SHA256
        and shared_keys == EXPECTED_SHARED_KEYS
        and tuple(prescribed_records) == expected_prescribed
        and len(prescribed_records) < 2
        and all(
            record["key"][0] != REQUIRED_CANCELLATION_DIRECTION_INDEX
            for record in prescribed_records
        )
    )
    if not frozen_family_excluded:
        raise ArithmeticError("the frozen collision-key obstruction changed")

    # Equal keys give two centre ratios per nonzero second centre: one same
    # physical edge and one negative edge.  Exactly one of them cancels.
    derived_shared_orbit_center_incidences = pair_class_count * 2 * (P - 1)
    derived_cancelling_center_incidences = pair_class_count * (P - 1)
    return {
        "schema": "p31_mobius_collision_key_v1",
        "classification": (
            "proved collision-key theorem and exhaustive finite obstruction "
            "for one frozen sixteen-half family"
        ),
        "p": P,
        "source_witness_sha256": SOURCE_WITNESS_SHA256,
        "option_catalog_sha256": OPTION_CATALOG_SHA256,
        "half_choices_target_auxiliary_scale": FROZEN_HALF_CHOICES,
        "raw_parallel_profile": raw_profile,
        "aggregate_signature_hex": f"{aggregate_signature:08x}",
        "correction_signature_hex": f"{correction_signature:08x}",
        "correction_signature_support": FROZEN_CORRECTION_SUPPORT,
        "required_cancellation_direction": REQUIRED_CANCELLATION_DIRECTION_INDEX,
        "nonorigin_keys_per_half": P - 2,
        "total_labelled_nonorigin_keys": len(catalog_rows),
        "labelled_key_catalog_sha256": catalog_sha256,
        "shared_collision_keys": tuple(shared_records),
        "shared_collision_key_count": len(shared_records),
        "distinct_half_pair_key_count": pair_class_count,
        "derived_shared_orbit_center_incidences": (
            derived_shared_orbit_center_incidences
        ),
        "derived_cancelling_center_incidences": (
            derived_cancelling_center_incidences
        ),
        "prescribed_endpoint_candidates": tuple(prescribed_records),
        "prescribed_endpoint_candidate_count": len(prescribed_records),
        "at_least_two_prescribed_endpoint_owners": False,
        "required_direction_prescribed_collision_exists": False,
        "frozen_sixteen_half_family_physically_excluded": frozen_family_excluded,
        "centre_enumeration_needed": False,
        "scope": (
            "the displayed hard-fixed sixteen-half family only; equal keys are "
            "necessary and sufficient for one chosen pair of nonorigin skeleton "
            "edges to share an inversion orbit, but do not by themselves solve "
            "the simultaneous sixteen-centre or compact-atom equations"
        ),
        "hard_fixed_branch_closed": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    """Return the symbolic statement together with the frozen exact replay."""
    replay = frozen_witness_collision_key_replay()
    return {
        "p": P,
        "collision_key": "(D,{U,V})",
        "U": "[tL-(t+1)M]",
        "V": "[L-M]",
        "D": "[t^2 L+(1-t^2)M]",
        "parameter_domain": "t != 0,-1",
        "directions_pairwise_distinct": True,
        "equal_keys_iff_unit_edges_are_homothetic": True,
        "same_orbit_center_condition": "c1=+/-lambda*c2",
        "cancellation_center_condition": "c1=-lambda*c2",
        "fixed_endpoint_prefilter": (
            "V=[L-M] must be in {A,B}; then "
            "t=det(M,H)/det(L-M,H) is unique"
        ),
        "frozen_replay": replay,
        "scope": replay["scope"],
        "residual_ii_closed": False,
        "proved": True,
    }
