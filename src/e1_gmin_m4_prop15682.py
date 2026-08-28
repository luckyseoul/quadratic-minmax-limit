#!/usr/bin/env python3
"""Prop. 15.682 -- close the p=31 next all-finite endpoint.

Proposition 15.681's integral paired-cube lift removes every positive
phase-zero residue at ``p=31,s=26``.  The sole residue-zero row has phase
minima

    10*b0 + b2 + 5*b26,   deficit 284,
    15*b2 + b26,          deficit 360,

against pair budget 650.  Pair-slack divisibility leaves fourteen
phase-labelled profiles.  Globally they are either 26-arcs with at least
three undetermined directions, or sets with one 3-secant and five
undetermined directions.  Deleting one triple point in the latter case
gives a 25-arc.

Coolsaet's exhaustive classification of complete arcs in ``PG(2,31)`` has
no complete arc of size 23 through 31.  Every finite arc extends to a
complete arc, so every 27- or 28-arc extends to size 32 and is contained in
the resulting conic.  Adjoining pairs from three undetermined infinity
points to the arcs above then forces one conic to contain three collinear
points, a contradiction.

This closes only ``p=31,s=26``.  At the same boundary the remaining smaller
endpoints are ``p=17,19,23,41``; Proposition 15.681 has reduced ``p=41`` to
residue zero.  Later all-finite sizes, the infinity-present remainder,
residual (ii), R1, global QVAR, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15681 import (
    endpoint_residue_ledger,
    pair_slack_divisibility,
)


ROOT = Path(__file__).resolve().parents[1]
P = 31
M = 16
S = 26
PERIOD = 32
PAIR_DEFICIT_BUDGET = S * (S - 1)


def _histogram(profile: tuple[int, ...]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(profile).items())}


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate every exact p=31 profile within a tight deficit cap."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<16")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor_value = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor_value
            if excess >= 0 and excess != 2:
                options.append((quotient, S - b, b))

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for _ in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                if new_used <= target and new_deficit <= deficit_cap:
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


def p31_residue_zero_profiles() -> dict[str, object]:
    """Classify every residue-zero profile allowed by the pair budget."""
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 360)
    phase_one = _profile_rows(1, 15, PAIR_DEFICIT_BUDGET - 284)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            total = deficit_zero + deficit_one
            slack = PAIR_DEFICIT_BUDGET - total
            if slack < 0 or slack % 4:
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
                    "total_deficit": total,
                    "pair_slack": slack,
                    "arc": slack == 0,
                    "global_secant_distribution": {
                        str(key): value for key, value in sorted(secants.items())
                    },
                    "undetermined_directions": secants[0],
                }
            )

    global_shapes = {
        tuple(sorted(row["global_secant_distribution"].items()))
        for row in candidates
    }
    expected_shapes = {
        tuple(sorted({"0": 5, "1": 1, "12": 16, "13": 10}.items())),
        tuple(sorted({"0": 5, "3": 1, "12": 16, "13": 10}.items())),
        tuple(
            sorted({"0": 4, "1": 1, "2": 1, "12": 16, "13": 10}.items())
        ),
        tuple(sorted({"0": 3, "1": 3, "12": 16, "13": 10}.items())),
        tuple(sorted({"0": 6, "12": 15, "13": 11}.items())),
        tuple(sorted({"0": 5, "2": 1, "12": 15, "13": 11}.items())),
        tuple(sorted({"0": 4, "1": 2, "12": 15, "13": 11}.items())),
    }
    if len(candidates) != 14 or global_shapes != expected_shapes:
        raise ArithmeticError("p=31 endpoint profile classification changed")
    near_arcs = [row for row in candidates if int(row["pair_slack"]) == 4]
    arcs = [row for row in candidates if int(row["pair_slack"]) == 0]
    if len(near_arcs) != 3 or len(arcs) != 11:
        raise ArithmeticError("p=31 arc/near-arc split changed")
    if min(int(row["undetermined_directions"]) for row in near_arcs) != 5:
        raise ArithmeticError("p=31 near-arc direction count changed")
    if min(int(row["undetermined_directions"]) for row in arcs) != 3:
        raise ArithmeticError("p=31 arc direction floor changed")
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
        "distinct_global_shapes": [dict(shape) for shape in sorted(global_shapes)],
        "near_arc_profile_count": len(near_arcs),
        "near_arc_geometry": (
            "exactly one affine line has occupancy three, every other line "
            "has occupancy at most two, and five directions are undetermined"
        ),
        "near_arc_delete_one_triple_point": (
            "a 25-arc preserving all five undetermined directions"
        ),
        "arc_profile_count": len(arcs),
        "arc_minimum_undetermined_directions": min(
            int(row["undetermined_directions"]) for row in arcs
        ),
        "proved": True,
    }


def p31_complete_arc_classification() -> dict[str, object]:
    """Record the exhaustive external complete-arc input."""
    return {
        "external_dependency": True,
        "source": (
            "K. Coolsaet, The Complete Arcs of PG(2,31), "
            "J. Combin. Des. 23 (2015), 522-533"
        ),
        "doi": "10.1002/jcd.21410",
        "classification_scope": "all complete arcs in PG(2,31)",
        "no_complete_arc_sizes": list(range(23, 32)),
        "largest_nonconic_complete_arc_size": 22,
        "complete_32_arcs": "nondegenerate conics by Segre's odd-order theorem",
        "finite_extension_argument": (
            "greedily extending any 27- or 28-arc must terminate at a "
            "complete arc; the classification forces size 32"
        ),
        "consequence": "every 27- and 28-arc in PG(2,31) is conic-contained",
        "proved_conditional_on_external_classification": True,
    }


def p31_geometric_exclusion() -> dict[str, object]:
    """Exclude every p=31 residue-zero profile by conic extensions."""
    profiles = p31_residue_zero_profiles()
    classification = p31_complete_arc_classification()
    proved = bool(
        profiles["proved"]
        and int(profiles["arc_minimum_undetermined_directions"]) >= 3
        and int(profiles["near_arc_profile_count"]) == 3
        and classification["proved_conditional_on_external_classification"]
    )
    if not proved:
        raise ArithmeticError("p=31 geometric exclusion audit failed")
    return {
        "arc_case": {
            "starting_size": 26,
            "minimum_undetermined_infinity_points": 3,
            "adjoin_two_size": 28,
            "classification_consequence": "contained in a conic",
        },
        "near_arc_case": {
            "starting_size": 26,
            "unique_3_secant": True,
            "delete_one_point_from_triple_size": 25,
            "result_is_arc": True,
            "undetermined_infinity_points_preserved": 5,
            "adjoin_two_size": 27,
            "classification_consequence": "contained in a conic",
        },
        "common_contradiction": (
            "choose three undetermined infinity points; conics through two "
            "extensions share at least 25 affine arc points and coincide, "
            "forcing one nondegenerate conic through three collinear points"
        ),
        "classification": classification,
        "excluded": True,
    }


def theorem_record() -> dict[str, object]:
    residue = endpoint_residue_ledger(P)
    survivors = [int(row["u0"]) for row in residue["pair_survivors"]]
    if survivors != [0, 2, 3, 4, 5, 6]:
        raise ArithmeticError("p=31 residue ledger changed")
    if residue["positive_residues_all_excluded"] is not True:
        raise ArithmeticError("p=31 positive residues returned")
    profiles = p31_residue_zero_profiles()
    geometry = p31_geometric_exclusion()
    proved = bool(
        residue["residue_zero_remains"] is True
        and profiles["proved"]
        and geometry["excluded"]
    )
    return {
        "prop": "15.682",
        "title": "The p=31 next all-finite endpoint is impossible",
        "proved": proved,
        "theorem": {
            "p31_s26_next_all_finite_endpoint": "EXCLUDED",
            "remaining_smaller_endpoints": [17, 19, 23, 41],
            "p41_status": "ONLY_RESIDUE_ZERO_REMAINS",
            "later_all_finite_sizes": "OPEN",
            "infinity_present_remainder": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "p31_residue_ledger": residue,
        "p31_profiles": profiles,
        "p31_geometric_exclusion": geometry,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.682 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15682.json"
    destination.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    print("Prop 15.682 p=31,s=26 next all-finite endpoint: excluded")
    print("  p=41: only residue zero remains")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
