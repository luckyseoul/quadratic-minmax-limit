#!/usr/bin/env python3
"""Exact direct one-half parallel design at the p=31 branch-C top row.

The certificate below uses sixteen independently scaled localized Mobius
halves.  It replays every parallel count from the thirty physical edges of
the half; no complementary-profile identity is assumed in the replay.

The raw ledger has hard profile ``15^2 14^14`` and opposite profile
``15^2 16^14``.  Taking the fixed antipodal edge and the unique cancelled
nonfixed orbit in the same opposite spatial direction 5 changes the latter
to ``15^3 16^13``, exactly the top target.

The center-one seed has disjoint inversion-orbit supports and certifies only
the parallel ledger.  A second exact certificate below keeps the same scaled
functionals, changes their nonzero centers, and realizes the required clean
cancellation.  That centered 479-edge graph still fails a necessary
transverse compact-atom budget, so neither certificate is a full common graph.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Functional,
    _functional_value,
    _negative_edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (
    _parallel_formula,
    mobius_parameter_edges,
    paley_direction_sign,
    paley_edge_sign,
)


P = 31
FIXED_DIRECTION_INDEX = 5
CANCELLATION_DIRECTION_INDEX = 5

PHYSICAL_FIXED_DIRECTION_INDEX = 5
PHYSICAL_CANCELLATION_DIRECTION_INDEX = 5
PHYSICAL_CENTERS = (27, 9, 28, 6, 19, 21, 28, 9, 19, 16, 9, 3, 4, 1, 30, 13)
PHYSICAL_FIXED_POINT = (5, 30)

# Actual scaled functionals (L,M), not merely projective directions.
HALVES: tuple[tuple[Functional, Functional], ...] = (
    ((8, 0), (20, 25)),
    ((9, 4), (9, 12)),
    ((10, 8), (13, 14)),
    ((12, 16), (29, 30)),
    ((14, 24), (1, 28)),
    ((16, 1), (4, 23)),
    ((17, 5), (16, 26)),
    ((18, 9), (12, 24)),
    ((20, 17), (10, 16)),
    ((21, 21), (0, 24)),
    ((25, 6), (26, 4)),
    ((26, 10), (8, 8)),
    ((0, 30), (5, 13)),
    ((1, 3), (27, 22)),
    ((5, 19), (20, 11)),
    ((6, 23), (8, 25)),
)


def _determinant(first: Functional, second: Functional) -> int:
    return (first[0] * second[1] - first[1] * second[0]) % P


def _projective_index(functional: Functional) -> int:
    directions = projective_functionals(P)
    matches = [
        index
        for index, direction in enumerate(directions)
        if _determinant(functional, direction) == 0
    ]
    if len(matches) != 1:
        raise ArithmeticError("functional lost its projective direction")
    return matches[0]


def _physical_parallel_profile(
    direction: Functional, auxiliary: Functional
) -> tuple[int, ...]:
    directions = projective_functionals(P)
    edges = mobius_parameter_edges(P, direction, auxiliary, center=1)
    counts = []
    for row in directions:
        counts.append(
            sum(
                _functional_value(P, row, edge[0])
                == _functional_value(P, row, edge[1])
                for edge in edges.values()
            )
        )
    formula = tuple(
        _parallel_formula(P, direction, auxiliary, row)
        for row in directions
    )
    if tuple(counts) != formula or sum(counts) != P - 1:
        raise ArithmeticError("physical edges disagree with the half formula")
    return formula


def _oriented_orbit_coefficients(
    direction: Functional, auxiliary: Functional, center: int = 1
) -> dict[Edge, int]:
    """Return the normalized antisymmetric coefficient on each orbit."""
    out: dict[Edge, int] = {}
    for edge in mobius_parameter_edges(P, direction, auxiliary, center).values():
        negative = _negative_edge(P, edge)
        orbit = min(edge, negative)
        z_value = -1 if edge == orbit else 1
        out[orbit] = paley_edge_sign(P, orbit) * z_value
    if len(out) != P - 1 or set(out.values()) != {-1, 1}:
        raise ArithmeticError("one half lost its thirty oriented orbits")
    return out


def direct_parallel_design_certificate() -> dict[str, object]:
    directions = projective_functionals(P)
    signs = tuple(paley_direction_sign(P, row) for row in directions)
    hard = tuple(index for index, sign in enumerate(signs) if sign == 1)
    opposite = tuple(index for index, sign in enumerate(signs) if sign == -1)

    records = []
    raw = [0] * len(directions)
    orbit_total: Counter[Edge] = Counter()
    for direction, auxiliary in HALVES:
        if _determinant(direction, auxiliary) == 0:
            raise ArithmeticError("an auxiliary became dependent")
        target_index = _projective_index(direction)
        auxiliary_index = _projective_index(auxiliary)
        profile = _physical_parallel_profile(direction, auxiliary)
        for index, value in enumerate(profile):
            raw[index] += value
        orbit_total.update(_oriented_orbit_coefficients(direction, auxiliary))
        records.append(
            {
                "target": direction,
                "target_direction_index": target_index,
                "auxiliary": auxiliary,
                "auxiliary_direction_index": auxiliary_index,
                "parallel_profile": profile,
            }
        )

    target_indices = tuple(row["target_direction_index"] for row in records)
    auxiliary_indices = tuple(row["auxiliary_direction_index"] for row in records)
    if set(target_indices) != set(hard) or len(set(auxiliary_indices)) != 16:
        raise ArithmeticError("the target/auxiliary SDR replay failed")

    final = list(raw)
    final[CANCELLATION_DIRECTION_INDEX] -= 2
    final[FIXED_DIRECTION_INDEX] += 1
    hard_values = sorted(final[index] for index in hard)
    opposite_values = sorted(final[index] for index in opposite)

    repeated_orbits = {
        orbit: value for orbit, value in orbit_total.items()
        if abs(value) != 1
    }
    proved = bool(
        sorted(raw[index] for index in hard) == [14] * 14 + [15] * 2
        and sorted(raw[index] for index in opposite) == [15] * 2 + [16] * 14
        and hard_values == [14] * 14 + [15] * 2
        and opposite_values == [15] * 3 + [16] * 13
        and signs[FIXED_DIRECTION_INDEX] == -1
        and FIXED_DIRECTION_INDEX == CANCELLATION_DIRECTION_INDEX
        and len(set(auxiliary_indices)) == 16
        and sum(signs[index] == 1 for index in auxiliary_indices) == 14
        and sum(signs[index] == -1 for index in auxiliary_indices) == 2
        and not repeated_orbits
    )
    if not proved:
        raise ArithmeticError("the direct p31 profile certificate changed")
    return {
        "p": P,
        "half_count": len(HALVES),
        "edges_per_half": P - 1,
        "target_direction_indices": target_indices,
        "auxiliary_direction_indices": auxiliary_indices,
        "auxiliary_hard_count": 14,
        "auxiliary_opposite_count": 2,
        "fixed_direction_index": FIXED_DIRECTION_INDEX,
        "fixed_direction": directions[FIXED_DIRECTION_INDEX],
        "fixed_direction_type": "opposite",
        "cancellation_direction_index": CANCELLATION_DIRECTION_INDEX,
        "cancellation_direction": directions[CANCELLATION_DIRECTION_INDEX],
        "raw_parallel_profile": tuple(raw),
        "final_parallel_profile": tuple(final),
        "raw_hard_multiset": tuple(sorted(raw[index] for index in hard)),
        "raw_opposite_multiset": tuple(sorted(raw[index] for index in opposite)),
        "final_hard_multiset": tuple(hard_values),
        "final_opposite_multiset": tuple(opposite_values),
        "raw_common_orbit_count": len(repeated_orbits),
        "physical_cancellation_realized_by_frozen_halves": False,
        "profiles_replayed_from_physical_edges": True,
        "halves": tuple(records),
        "scope": "parallel profile only; centre-dependent physical cancellation remains open",
        "proved": True,
    }


def _spatial_direction_index(edge: Edge) -> int:
    difference = (
        (edge[0][0] - edge[1][0]) % P,
        (edge[0][1] - edge[1][1]) % P,
    )
    annihilators = [
        index
        for index, direction in enumerate(projective_functionals(P))
        if _functional_value(P, direction, difference) == 0
    ]
    if len(annihilators) != 1:
        raise ArithmeticError("edge lost its unique spatial direction")
    return annihilators[0]


def centered_physical_parallel_design_certificate() -> dict[str, object]:
    """Replay the clean one-cancellation physical realization.

    Keeping the frozen scaled functionals in ``HALVES``, use the displayed
    nonzero centre list.  Halves 2 and 13 then meet in the single inversion
    orbit ``{(2,25),(29,1)}`` with opposite normalized coefficients.  All
    other pairs are disjoint.  The cancelled orbit and the fixed antipodal
    edge both have spatial direction 5.
    """
    directions = projective_functionals(P)
    direction_signs = tuple(
        paley_direction_sign(P, direction) for direction in directions
    )
    hard = tuple(index for index, sign in enumerate(direction_signs) if sign == 1)
    opposite = tuple(
        index for index, sign in enumerate(direction_signs) if sign == -1
    )

    half_orbits = tuple(
        _oriented_orbit_coefficients(direction, auxiliary, center)
        for (direction, auxiliary), center in zip(
            HALVES, PHYSICAL_CENTERS, strict=True
        )
    )
    total: Counter[Edge] = Counter()
    for orbit_map in half_orbits:
        total.update(orbit_map)

    pair_intersections = []
    for first in range(len(half_orbits)):
        for second in range(first + 1, len(half_orbits)):
            for orbit in set(half_orbits[first]) & set(half_orbits[second]):
                pair_intersections.append(
                    {
                        "halves": (first, second),
                        "orbit": orbit,
                        "first_coefficient": half_orbits[first][orbit],
                        "second_coefficient": half_orbits[second][orbit],
                        "spatial_direction_index": _spatial_direction_index(orbit),
                        "contains_origin": (0, 0) in orbit,
                    }
                )

    cancelled = tuple(orbit for orbit, value in total.items() if value == 0)
    surviving = {orbit: value for orbit, value in total.items() if value}
    if any(abs(value) != 1 for value in surviving.values()):
        raise ArithmeticError("the centered realization lost ternarity")

    # The coefficient stored above is the normalized inversion difference;
    # choose the canonical edge exactly when that coefficient is +1.
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in surviving.items()
    ]
    negative_fixed_point = (
        -PHYSICAL_FIXED_POINT[0] % P,
        -PHYSICAL_FIXED_POINT[1] % P,
    )
    fixed_edge = tuple(sorted((PHYSICAL_FIXED_POINT, negative_fixed_point)))
    if fixed_edge in graph:
        raise ArithmeticError("the fixed edge collided with a nonfixed orbit")
    graph.append(fixed_edge)  # type: ignore[arg-type]
    graph = sorted(graph)

    normalized_source = {
        edge: paley_edge_sign(P, edge) for edge in graph
    }
    image = edge_radon_image(P, normalized_source)
    final_profile = tuple(
        direction_signs[index] * image.get(("P", index), 0)
        for index in range(P + 1)
    )
    graph_bytes = json.dumps(graph, separators=(",", ":")).encode()
    origin_edges = tuple(edge for edge in graph if (0, 0) in edge)

    expected_collision = ((2, 25), (29, 1))
    proved = bool(
        len(pair_intersections) == 1
        and pair_intersections[0]["halves"] == (2, 13)
        and pair_intersections[0]["orbit"] == expected_collision
        and pair_intersections[0]["first_coefficient"]
        == -pair_intersections[0]["second_coefficient"]
        and pair_intersections[0]["spatial_direction_index"]
        == PHYSICAL_CANCELLATION_DIRECTION_INDEX
        and not pair_intersections[0]["contains_origin"]
        and cancelled == (expected_collision,)
        and len(surviving) == 478
        and len(graph) == 479
        and _spatial_direction_index(fixed_edge)
        == PHYSICAL_FIXED_DIRECTION_INDEX
        and PHYSICAL_FIXED_DIRECTION_INDEX
        == PHYSICAL_CANCELLATION_DIRECTION_INDEX
        and len(origin_edges) == 16
        and sorted(final_profile[index] for index in hard)
        == [14] * 14 + [15] * 2
        and sorted(final_profile[index] for index in opposite)
        == [15] * 3 + [16] * 13
    )
    if not proved:
        raise ArithmeticError("the centered physical p31 certificate changed")
    return {
        "p": P,
        "centers": PHYSICAL_CENTERS,
        "unique_pair_intersection_count": len(pair_intersections),
        "pair_intersections": tuple(pair_intersections),
        "cancelled_orbit_count": len(cancelled),
        "cancelled_orbit": cancelled[0],
        "cancellation_is_nonorigin": True,
        "cancellation_direction_index": PHYSICAL_CANCELLATION_DIRECTION_INDEX,
        "cancellation_direction": directions[
            PHYSICAL_CANCELLATION_DIRECTION_INDEX
        ],
        "surviving_nonfixed_orbit_count": len(surviving),
        "fixed_edge": fixed_edge,
        "fixed_direction_index": PHYSICAL_FIXED_DIRECTION_INDEX,
        "fixed_direction": directions[PHYSICAL_FIXED_DIRECTION_INDEX],
        "fixed_direction_type": "opposite",
        "graph_edge_count": len(graph),
        "origin_edge_count": len(origin_edges),
        "final_parallel_profile": final_profile,
        "final_hard_multiset": tuple(
            sorted(final_profile[index] for index in hard)
        ),
        "final_opposite_multiset": tuple(
            sorted(final_profile[index] for index in opposite)
        ),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "exact_physical_edge_radon_replay": True,
        "same_sign_overlap_count": 0,
        "triple_overlap_count": 0,
        "scope": (
            "one explicit nonzero center profile; parallel quotas and the "
            "antisymmetric support only, not transverse compact atoms"
        ),
        "proved": True,
    }


def _canonical_center(
    target: Functional, center: int
) -> tuple[int, int]:
    """Return ``(direction_index, center)`` in the canonical row scaling."""
    directions = projective_functionals(P)
    direction_index = _projective_index(target)
    canonical = directions[direction_index]
    if target[0] % P:
        multiplier = canonical[0] * pow(target[0], -1, P) % P
    else:
        multiplier = canonical[1] * pow(target[1], -1, P) % P
    if (
        multiplier * target[0] % P,
        multiplier * target[1] % P,
    ) != canonical:
        raise ArithmeticError("canonical target scaling failed")
    return direction_index, multiplier * center % P


def centered_physical_graph() -> dict[str, object]:
    """Expose the deterministic 479-edge graph and canonical star centers.

    The returned edges are canonical unordered point pairs.  To form the
    normalized source for a full row replay, give each edge coefficient
    ``paley_edge_sign(31, edge)`` and multiply each output row by its Paley
    direction sign.
    """
    total: Counter[Edge] = Counter()
    for (direction, auxiliary), center in zip(
        HALVES, PHYSICAL_CENTERS, strict=True
    ):
        total.update(
            _oriented_orbit_coefficients(direction, auxiliary, center)
        )
    surviving = {orbit: value for orbit, value in total.items() if value}
    if len(surviving) != 478 or any(
        abs(value) != 1 for value in surviving.values()
    ):
        raise ArithmeticError("the exposed graph lost its ternary support")
    graph = [
        orbit if value == 1 else _negative_edge(P, orbit)
        for orbit, value in surviving.items()
    ]
    negative_fixed_point = (
        -PHYSICAL_FIXED_POINT[0] % P,
        -PHYSICAL_FIXED_POINT[1] % P,
    )
    fixed_edge = tuple(sorted((PHYSICAL_FIXED_POINT, negative_fixed_point)))
    graph.append(fixed_edge)  # type: ignore[arg-type]
    edges = tuple(sorted(graph))
    if len(edges) != 479 or len(set(edges)) != 479:
        raise ArithmeticError("the exposed physical graph repeated an edge")
    graph_bytes = json.dumps(edges, separators=(",", ":")).encode()
    graph_sha256 = hashlib.sha256(graph_bytes).hexdigest()
    expected_sha256 = "c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d"
    if graph_sha256 != expected_sha256:
        raise ArithmeticError("the exposed graph hash changed")

    directions = projective_functionals(P)
    center_records = []
    for half_index, ((target, auxiliary), center) in enumerate(
        zip(HALVES, PHYSICAL_CENTERS, strict=True)
    ):
        direction_index, canonical_center = _canonical_center(target, center)
        center_records.append(
            {
                "half_index": half_index,
                "target_functional": target,
                "target_center_in_displayed_scaling": center,
                "target_direction_index": direction_index,
                "canonical_direction": directions[direction_index],
                "canonical_target_center": canonical_center,
                "auxiliary_functional": auxiliary,
                "auxiliary_direction_index": _projective_index(auxiliary),
            }
        )
    if len({row["target_direction_index"] for row in center_records}) != 16:
        raise ArithmeticError("the canonical target-center map repeated a row")
    return {
        "p": P,
        "edges": edges,
        "edge_count": len(edges),
        "graph_sha256": graph_sha256,
        "fixed_edge": fixed_edge,
        "fixed_direction_index": PHYSICAL_FIXED_DIRECTION_INDEX,
        "hard_target_centers": tuple(center_records),
        "normalized_source_rule": (
            "coefficient(edge)=paley_edge_sign(p,edge), then multiply row "
            "image by paley_direction_sign(p,row)"
        ),
        "proved": True,
    }


def transverse_compact_l1_diagnostic() -> dict[str, object]:
    """Test this frozen graph against the necessary compact-atom l1 budget.

    A hard row is ``-S_j`` plus ``P-3`` compact three-edge atoms.  An
    opposite row has six all-positive three-edge atoms plus ``Q-9`` compact
    atoms, again ``Q-3`` atoms total.  Therefore, after subtracting the hard
    literal star, every row must have transverse l1 norm at most
    ``3*(parallel-3)``.  This diagnostic concerns only the frozen graph; it
    is not an obstruction to other Mobius profile designs.
    """
    exposed = centered_physical_graph()
    directions = projective_functionals(P)
    signs = tuple(paley_direction_sign(P, row) for row in directions)
    centers = {
        int(record["target_direction_index"]): int(
            record["canonical_target_center"]
        )
        for record in exposed["hard_target_centers"]
    }
    image = edge_radon_image(
        P,
        {
            edge: paley_edge_sign(P, edge)
            for edge in exposed["edges"]
        },
    )
    rows = []
    for direction_index, sign in enumerate(signs):
        parallel = sign * image.get(("P", direction_index), 0)
        center = centers.get(direction_index)
        coefficients = []
        for left in range(P):
            for right in range(left + 1, P):
                value = sign * image.get(
                    ("K", direction_index, left, right), 0
                )
                if center is not None and center in (left, right):
                    # Subtract the target baseline -S_center.
                    value += 1
                coefficients.append(value)
        transverse_l1 = sum(abs(value) for value in coefficients)
        budget = 3 * (parallel - 3)
        rows.append(
            {
                "direction_index": direction_index,
                "direction_sign": sign,
                "parallel_count": parallel,
                "canonical_literal_center": center,
                "transverse_residual_l1": transverse_l1,
                "necessary_compact_atom_l1_budget": budget,
                "excess": transverse_l1 - budget,
                "within_budget": transverse_l1 <= budget,
            }
        )
    excesses = tuple(int(row["excess"]) for row in rows)
    proved = bool(
        len(rows) == P + 1
        and all(excess > 0 for excess in excesses)
        and min(excesses) == 122
        and max(excesses) == 194
    )
    if not proved:
        raise ArithmeticError("the frozen graph l1 diagnostic changed")
    return {
        "row_count": len(rows),
        "violating_row_count": sum(not row["within_budget"] for row in rows),
        "minimum_l1_excess": min(excesses),
        "maximum_l1_excess": max(excesses),
        "rows": tuple(rows),
        "frozen_graph_fails_necessary_compact_l1_budget": True,
        "scope": "this graph only; not a profile-family obstruction",
        "proved": True,
    }


if __name__ == "__main__":
    print(
        json.dumps(
            {
                "profile_ledger": direct_parallel_design_certificate(),
                "centered_physical_ledger": (
                    centered_physical_parallel_design_certificate()
                ),
            },
            indent=2,
        )
    )
