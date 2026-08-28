#!/usr/bin/env python3
"""Prop. 15.683 -- close the p=41 next all-finite endpoint.

Proposition 15.681 removes every positive phase-zero residue at
``p=41,s=34``.  Exact enumeration of the residue-zero row leaves nine
phase-labelled profiles: seven 34-arcs and two sets with exactly one
3-secant.  Every profile has fourteen phase-zero directions with seventeen
secants and twenty phase-one directions with sixteen secants.  The other
eight directions carry only three floor-secants in total.

The new ingredient is Segre's tangent envelope in the polynomial form of
Ball--Lavrauw, *Planar arcs*, Theorem 11.  For an arc of size ``q+2-t`` in
odd order, when ``|A|>2t+2``, there is a degree-``2t`` polynomial in the
dual plane whose restriction to each point-pencil is the square of the
tangent polynomial.

For a slack-zero profile, the eight exceptional directions contain at
least 28 tangents each.  Their eight dual pencil lines therefore divide the
degree-18 envelope twice, leaving a conic.  The three exceptional secants
touch at least three arc points.  At each such point at least two of the
remaining tangent points give distinct double zeros on its dual line, so
that line must divide the residual conic.  Three distinct line components
cannot fit in degree two.

For a slack-four profile, choose a point of the unique 3-secant whose
deletion preserves the unique exceptional floor-secant.  (If that
floor-secant is the triple, any deletion leaves a pair; otherwise at least
one of the three triple points is outside its ordinary pair.)  The result is
a 33-arc with seven directions containing 33 tangents and the eighth
exceptional direction containing 31.  Their squares leave a quartic from
the degree-20 envelope.  Each endpoint of the surviving exceptional pair
has three remaining tangents, forcing both point-pencil lines out of the
quartic.  The residual conic would then have to contain the point-pencil
line of every other arc point, again impossible.

This closes only ``p=41,s=34``.  The same second finite boundary remains
open at ``p=17,19,23``; later all-finite sizes, the infinity-present
remainder, residual (ii), R1, global QVAR, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15675 import first_even_survivor
from e1_gmin_m4_prop15681 import endpoint_residue_ledger, pair_slack_divisibility


ROOT = Path(__file__).resolve().parents[1]
P = 41
M = 21
S = first_even_survivor(P) + 2
PERIOD = P + 1
PAIR_DEFICIT_BUDGET = S * (S - 1)


def _histogram(profile: tuple[int, ...]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(profile).items())}


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate every exact p=41 direction profile under a tight cap."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<21")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor_value = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor_value
            if excess >= 0 and excess != 2:
                options.append((quotient, S - b, b))

    infinity = deficit_cap + S * M + 1
    completion_minimum = [[infinity] * (target + 1) for _ in range(M + 1)]
    completion_minimum[0][0] = 0
    for count in range(1, M + 1):
        for quotient_sum in range(target + 1):
            completion_minimum[count][quotient_sum] = min(
                (
                    added
                    + completion_minimum[count - 1][quotient_sum - quotient]
                    for quotient, added, _b in options
                    if quotient <= quotient_sum
                ),
                default=infinity,
            )

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for count in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                remaining_count = M - count - 1
                remaining_sum = target - new_used
                if (
                    new_used <= target
                    and new_deficit <= deficit_cap
                    and new_deficit
                    + completion_minimum[remaining_count][remaining_sum]
                    <= deficit_cap
                ):
                    next_states.add(
                        (new_used, new_deficit, tuple(sorted(profile + (b,))))
                    )
        states = next_states
    return tuple(
        sorted(
            (deficit, profile)
            for used, deficit, profile in states
            if used == target
        )
    )


def p41_residue_zero_profiles() -> dict[str, object]:
    """Classify all residue-zero profiles left by the exact pair ledger."""
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 640)
    phase_one = _profile_rows(1, 20, PAIR_DEFICIT_BUDGET - 476)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            pair_slack = PAIR_DEFICIT_BUDGET - deficit_zero - deficit_one
            if pair_slack < 0 or pair_slack % 4:
                continue
            secants = Counter((S - b) // 2 for b in profile_zero + profile_one)
            candidates.append(
                {
                    "phase_profiles_b": {
                        "0": _histogram(profile_zero),
                        "1": _histogram(profile_one),
                    },
                    "phase_deficits": {
                        "0": deficit_zero,
                        "1": deficit_one,
                    },
                    "total_deficit": deficit_zero + deficit_one,
                    "pair_slack": pair_slack,
                    "arc": pair_slack == 0,
                    "global_secant_distribution": {
                        str(key): value for key, value in sorted(secants.items())
                    },
                    "undetermined_directions": secants[0],
                }
            )

    shapes = {
        (int(row["pair_slack"]), tuple(row["global_secant_distribution"].items()))
        for row in candidates
    }
    expected_shapes = {
        (4, tuple({"0": 7, "1": 1, "16": 20, "17": 14}.items())),
        (0, tuple({"0": 7, "3": 1, "16": 20, "17": 14}.items())),
        (
            0,
            tuple({"0": 6, "1": 1, "2": 1, "16": 20, "17": 14}.items()),
        ),
        (0, tuple({"0": 5, "1": 3, "16": 20, "17": 14}.items())),
    }
    arcs = [row for row in candidates if row["pair_slack"] == 0]
    near_arcs = [row for row in candidates if row["pair_slack"] == 4]
    if (
        S != 34
        or len(phase_zero) != 7
        or len(phase_one) != 4
        or len(candidates) != 9
        or len(arcs) != 7
        or len(near_arcs) != 2
        or shapes != expected_shapes
    ):
        raise ArithmeticError("p=41 residue-zero profile census changed")
    if not all(row["undetermined_directions"] >= 5 for row in candidates):
        raise ArithmeticError("p=41 undetermined-direction floor changed")
    return {
        "p": P,
        "s": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_zero_rows_within_cap": [
            {"deficit": deficit, "profile": _histogram(profile)}
            for deficit, profile in phase_zero
        ],
        "phase_one_rows_within_cap": [
            {"deficit": deficit, "profile": _histogram(profile)}
            for deficit, profile in phase_one
        ],
        "pair_slack_divisibility": pair_slack_divisibility(),
        "profiles": candidates,
        "profile_count": len(candidates),
        "distinct_global_shape_count": len(shapes),
        "arc_profile_count": len(arcs),
        "near_arc_profile_count": len(near_arcs),
        "minimum_undetermined_directions": min(
            int(row["undetermined_directions"]) for row in candidates
        ),
        "slack_zero_geometry": "34-arc",
        "slack_four_geometry": (
            "exactly one 3-secant and every other line has occupancy at most two"
        ),
        "proved": True,
    }


def tangent_envelope_input() -> dict[str, object]:
    """Record the external theorem used by both geometric contradictions."""
    return {
        "external_dependency": True,
        "source": (
            "S. Ball and M. Lavrauw, Planar arcs, Journal of Combinatorial "
            "Theory Series A 160 (2018), 261-287"
        ),
        "doi": "10.1016/j.jcta.2018.06.015",
        "arxiv": "1705.10940v4",
        "theorem": 11,
        "odd_order_statement": (
            "if an arc A has size q+2-t and |A|>2t+2, a nonzero degree-2t "
            "dual polynomial restricts on every point-pencil to the square "
            "of that point's tangent polynomial"
        ),
        "component_argument": (
            "if a direction pencil contains b>2t tangents it divides the "
            "envelope; the double-root restriction then leaves b zeros on "
            "the degree-(2t-1) quotient, so b>2t-1 forces its square"
        ),
        "proved": True,
    }


def _minimum_vertices_for_distinct_edges(edge_count: int) -> int:
    vertices = 0
    while math.comb(vertices, 2) < edge_count:
        vertices += 1
    return vertices


def p41_arc_envelope_exclusion() -> dict[str, object]:
    """Exclude all seven slack-zero profiles with one envelope argument."""
    census = p41_residue_zero_profiles()
    arcs = [row for row in census["profiles"] if row["pair_slack"] == 0]
    rows = []
    for profile in arcs:
        secants = {
            int(key): int(value)
            for key, value in profile["global_secant_distribution"].items()
        }
        exceptional = {
            key: value for key, value in secants.items() if key not in (16, 17)
        }
        rows.append(
            {
                "exceptional_direction_count": sum(exceptional.values()),
                "exceptional_secant_edges": sum(
                    key * value for key, value in exceptional.items()
                ),
                "minimum_exceptional_tangents": min(
                    S - 2 * key for key in exceptional
                ),
            }
        )
    if not all(
        row["exceptional_direction_count"] == 8
        and row["exceptional_secant_edges"] == 3
        and row["minimum_exceptional_tangents"] >= 28
        for row in rows
    ):
        raise ArithmeticError("p=41 arc envelope inputs changed")

    arc_size = 34
    tangent_number = P + 2 - arc_size
    envelope_degree = 2 * tangent_number
    forced_double_lines = 8
    residual_degree = envelope_degree - 2 * forced_double_lines
    incident_point_floor = _minimum_vertices_for_distinct_edges(3)
    minimum_low_tangents_at_incident_point = 2
    forced_root_multiplicity = 2 * minimum_low_tangents_at_incident_point
    proved = (
        arc_size > 2 * tangent_number + 2
        and 28 > envelope_degree
        and 28 > envelope_degree - 1
        and residual_degree == 2
        and incident_point_floor == 3
        and forced_root_multiplicity > residual_degree
        and incident_point_floor > residual_degree
    )
    if not proved:
        raise ArithmeticError("p=41 arc envelope contradiction changed")
    return {
        "profile_count": len(arcs),
        "arc_size": arc_size,
        "tangents_per_point": tangent_number,
        "envelope_degree": envelope_degree,
        "exceptional_directions": forced_double_lines,
        "minimum_tangents_per_exceptional_direction": 28,
        "forced_double_direction_component_degree": 2 * forced_double_lines,
        "residual_curve_degree": residual_degree,
        "exceptional_secant_edges": 3,
        "minimum_incident_arc_points": incident_point_floor,
        "minimum_low_tangents_at_each_incident_point": (
            minimum_low_tangents_at_incident_point
        ),
        "double_zero_multiplicity_on_point_pencil": forced_root_multiplicity,
        "forced_point_pencil_components": incident_point_floor,
        "contradiction": "at least three distinct line components in a conic",
        "proved": True,
    }


def p41_near_arc_envelope_exclusion() -> dict[str, object]:
    """Exclude both one-triple profiles after one lossless deletion."""
    census = p41_residue_zero_profiles()
    near_arcs = [row for row in census["profiles"] if row["pair_slack"] == 4]
    expected = {"0": 7, "1": 1, "16": 20, "17": 14}
    if not all(row["global_secant_distribution"] == expected for row in near_arcs):
        raise ArithmeticError("p=41 near-arc global profile changed")

    arc_size_after_deletion = 33
    tangent_number = P + 2 - arc_size_after_deletion
    envelope_degree = 2 * tangent_number
    no_secant_directions = 7
    exceptional_pair_directions = 1
    high_direction_count = no_secant_directions + exceptional_pair_directions
    high_tangent_counts = [33] * no_secant_directions + [31]
    residual_quartic_degree = envelope_degree - 2 * high_direction_count
    surviving_high_pair_endpoints = 2
    low_tangents_at_high_pair_endpoint = 3
    first_forced_multiplicity = 2 * low_tangents_at_high_pair_endpoint
    residual_conic_degree = (
        residual_quartic_degree - surviving_high_pair_endpoints
    )
    low_tangents_at_other_point = 2
    second_forced_multiplicity = 2 * low_tangents_at_other_point
    other_point_count = arc_size_after_deletion - surviving_high_pair_endpoints
    proved = (
        arc_size_after_deletion > 2 * tangent_number + 2
        and min(high_tangent_counts) > envelope_degree
        and min(high_tangent_counts) > envelope_degree - 1
        and residual_quartic_degree == 4
        and first_forced_multiplicity > residual_quartic_degree
        and residual_conic_degree == 2
        and second_forced_multiplicity > residual_conic_degree
        and other_point_count >= 3
    )
    if not proved:
        raise ArithmeticError("p=41 near-arc envelope contradiction changed")
    return {
        "profile_count": len(near_arcs),
        "unique_triple_deletion": {
            "deleted_points": 1,
            "choice": (
                "preserve the unique exceptional floor-secant: if it is "
                "the triple it becomes a pair; otherwise choose a triple "
                "point outside its ordinary pair"
            ),
            "resulting_arc_size": arc_size_after_deletion,
            "seven_no_secant_direction_tangents": 33,
            "surviving_exceptional_pair_direction_tangents": 31,
        },
        "tangents_per_point": tangent_number,
        "envelope_degree": envelope_degree,
        "forced_double_direction_components": high_direction_count,
        "forced_double_direction_component_degree": 2 * high_direction_count,
        "residual_curve_degree": residual_quartic_degree,
        "surviving_high_pair_endpoints": surviving_high_pair_endpoints,
        "low_tangents_at_each_high_pair_endpoint": (
            low_tangents_at_high_pair_endpoint
        ),
        "first_forced_point_pencil_components": surviving_high_pair_endpoints,
        "residual_after_two_point_pencils": residual_conic_degree,
        "other_arc_points": other_point_count,
        "low_tangents_at_each_other_point": low_tangents_at_other_point,
        "contradiction": (
            "three ordinary point-pencil lines would divide the residual conic"
        ),
        "proved": True,
    }


def p41_endpoint_theorem() -> dict[str, object]:
    """Combine Proposition 15.681's lift with both envelope exclusions."""
    ledger = endpoint_residue_ledger(P)
    census = p41_residue_zero_profiles()
    arc = p41_arc_envelope_exclusion()
    near = p41_near_arc_envelope_exclusion()
    if not ledger["positive_residues_all_excluded"]:
        raise ArithmeticError("Proposition 15.681 no longer removes positive residues")
    if [row["u0"] for row in ledger["residue_zero_rows"]] != [0]:
        raise ArithmeticError("p=41 residue-zero handoff changed")
    proved = (
        bool(census["proved"])
        and int(census["profile_count"])
        == int(arc["profile_count"]) + int(near["profile_count"])
        and bool(arc["proved"])
        and bool(near["proved"])
    )
    if not proved:
        raise ArithmeticError("p=41 endpoint proof did not close every profile")
    return {
        "p": P,
        "s": S,
        "positive_residues_excluded_by_prop15681": True,
        "residue_zero_profile_count": census["profile_count"],
        "arc_profiles_excluded": arc["profile_count"],
        "near_arc_profiles_excluded": near["profile_count"],
        "second_all_finite_endpoint_closed": True,
        "remaining_same_boundary_primes": [17, 19, 23],
        "later_boundaries_open": True,
        "top_level_gates_changed": False,
        "proved": True,
    }


def build_evidence() -> dict[str, object]:
    return {
        "proposition": "15.683",
        "title": "the p=41 next all-finite endpoint is impossible",
        "tangent_envelope_input": tangent_envelope_input(),
        "residue_zero_profiles": p41_residue_zero_profiles(),
        "arc_envelope_exclusion": p41_arc_envelope_exclusion(),
        "near_arc_envelope_exclusion": p41_near_arc_envelope_exclusion(),
        "theorem": p41_endpoint_theorem(),
    }


def main() -> None:
    output = ROOT / "evidence" / "e1_gmin_m4_prop15683.json"
    output.write_text(json.dumps(build_evidence(), indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
