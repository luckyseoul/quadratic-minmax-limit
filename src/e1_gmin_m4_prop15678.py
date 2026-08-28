#!/usr/bin/env python3
"""Prop. 15.678 -- close the exceptional p=17 all-finite endpoint.

At ``p=17`` the first even all-finite size above ``3(p-1)/4`` is
``s=14``.  Proposition 15.677 deliberately left this endpoint open because
phase zero has an additional same-type residue ``u_0=0``.

The exact floor/mean ledger leaves only ``u_1=8`` in phase one.  In phase
zero, ``u_0=2`` is excluded by the six-unit quadratic-lift floor and
``u_0=3`` by the coefficient ``l1`` bound; all residues at least four exceed
the pair budget.  For ``u_0=0``, the fact that pair slack is divisible by
four leaves exactly two profiles.  In secant notation both have the same
global distribution

    six directions with 7 secants,
    eight directions with 6 secants,
    one direction with 1 secant,
    three directions with 0 secants.

Their total deficit is the pair budget, so the fourteen affine points form
an arc.  Adjoin any two of the three undetermined points on the line at
infinity to obtain a 16-arc in ``PG(2,17)``.  Sticker's exhaustive
classification records exactly one PGL class of 16-arcs; since deleting two
points from a conic gives such an arc, every 16-arc is conic-contained.

If ``S`` is the original 14-arc, the conic has four points outside ``S``.
The third undetermined infinity point is off the conic, because the line at
infinity already contains the two adjoined conic points.  An external point
of a conic in odd order lies on eight conic secants and an internal point on
nine.  Deleting four conic points destroys at most four of these secants, so
the third point lies on at least four secants of ``S`` -- a contradiction.

The finite classification is an explicit external dependency, not a local
re-enumeration.  This proposition closes only the exceptional first-survivor
endpoint.  Later all-finite sizes, residual (ii), R1, global QVAR, Type I,
and the limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15642 import nonbaseline_scaled_cost_floor
from e1_gmin_m4_prop15669 import full_symbolic_floor


ROOT = Path(__file__).resolve().parents[1]
P = 17
S = 14
M = 9
PERIOD = 18
PAIR_DEFICIT_BUDGET = S * (S - 1)


def _histogram(profile: tuple[int, ...]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(profile).items())}


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int = PAIR_DEFICIT_BUDGET
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate the exact same-type floor relaxation at one residue.

    Means are ``2u+18k`` and their quotient sum is ``9-u``.  The exclusion
    of floor plus two is Proposition 15.642's nonzero-lift consequence.
    Profiles are sorted tuples so duplicate allocations collapse exactly.
    """
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<9")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor
            if excess >= 0 and excess != 2:
                options.append((quotient, S - b, b))

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for _ in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added_deficit, b in options:
                new_used = used + quotient
                new_deficit = deficit + added_deficit
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


def type_residue_ledger() -> dict[str, object]:
    """Return every relaxed residue minimum at the p=17 endpoint."""
    rows: dict[str, dict[str, object]] = {}
    for phase in (0, 1):
        phase_rows: dict[str, object] = {}
        for u in range(M):
            profiles = _profile_rows(phase, u)
            if not profiles:
                continue
            minimum = profiles[0][0]
            phase_rows[str(u)] = {
                "minimum_deficit": minimum,
                "minimizing_profiles": [
                    _histogram(profile)
                    for deficit, profile in profiles
                    if deficit == minimum
                ],
                "profile_count_below_pair_budget": len(profiles),
            }
        rows[str(phase)] = phase_rows

    expected_zero = {
        "0": 84,
        "2": 82,
        "3": 84,
        "4": 96,
        "5": 98,
        "6": 110,
        "7": 112,
        "8": 112,
    }
    observed_zero = {
        u: int(row["minimum_deficit"])
        for u, row in rows["0"].items()  # type: ignore[union-attr]
    }
    observed_one = {
        u: int(row["minimum_deficit"])
        for u, row in rows["1"].items()  # type: ignore[union-attr]
    }
    if observed_zero != expected_zero or observed_one != {"8": 96}:
        raise ArithmeticError("p=17 residue minima changed")
    return {
        "p": P,
        "s": S,
        "m": M,
        "period": PERIOD,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_residue_rows": rows,
        "phase_one_only_residue": 8,
        "proved": True,
    }


def u2_lift_exclusion() -> dict[str, object]:
    """Exclude phase-zero u=2 using the exact p=17 lift floor."""
    u = 2
    quotient_sum = M - u
    scaled_mean = 2 * u
    minimum_nonzero_b_floor = min(
        full_symbolic_floor(P, b, 0) for b in range(2, S + 1, 2)
    )
    lift_floor = nonbaseline_scaled_cost_floor(P)
    excluded = (
        quotient_sum < M
        and scaled_mean < minimum_nonzero_b_floor
        and lift_floor > scaled_mean
    )
    if not excluded or lift_floor != 6:
        raise ArithmeticError("p=17 u=2 lift contradiction changed")
    return {
        "u": u,
        "quotient_sum": quotient_sum,
        "direction_count": M,
        "zero_quotient_direction_forced": True,
        "zero_quotient_scaled_mean": scaled_mean,
        "minimum_phase_zero_floor_at_nonzero_even_b": minimum_nonzero_b_floor,
        "forced_b": 0,
        "factorization": "A=2B with B nonzero and nonnegative",
        "nonzero_B_scaled_cost_floor": lift_floor,
        "excluded": True,
    }


def _balanced_base_l1(scalar: int, total: int) -> dict[str, int]:
    """Minimum of sum_{s<t}|scalar-n_s-n_t| at fixed integer sum."""
    low, high_count = divmod(total, P)
    high = low + 1
    low_count = P - high_count
    value = (
        comb(low_count, 2) * abs(scalar - 2 * low)
        + low_count * high_count * abs(scalar - low - high)
        + comb(high_count, 2) * abs(scalar - 2 * high)
    )
    return {
        "low": low,
        "low_count": low_count,
        "high": high,
        "high_count": high_count,
        "base_l1_minimum": value,
    }


def u3_coefficient_exclusion() -> dict[str, object]:
    """Exclude phase-zero u=3 by the pre-lift coefficient l1 ledger."""
    q = (P - 1) // 2
    u = 3
    j_candidates = [j for j in range(1, 8) if (u - j) % q == 0]
    if j_candidates != [3]:
        raise ArithmeticError("p=17 u=3 baseline offset changed")
    j = j_candidates[0]
    rows = []
    for ell in range(0, 9):
        phase_zero_finite = M * j - u
        phase_one_finite = M * ell + 1
        finite = phase_zero_finite + phase_one_finite
        infinity = 4 * P + 1 - finite
        if infinity < 0 or infinity % 2:
            continue
        if infinity > S + 2 * finite:
            continue
        if (infinity + ell - 4) % q:
            continue
        scalar = (infinity + ell - 4) // q
        balanced = _balanced_base_l1(scalar, infinity)
        lower = max(0, balanced["base_l1_minimum"] - 1)
        capacity = finite - ell
        rows.append(
            {
                "ell": ell,
                "I": infinity,
                "E": finite,
                "xnor_scalar": scalar,
                "balanced_base": balanced,
                "xnor_l1_lower_bound": lower,
                "transverse_edge_capacity": capacity,
                "excluded_by_l1": lower > capacity,
            }
        )
    expected = [
        (0, 44, 25, 65, 25),
        (2, 26, 43, 63, 41),
        (4, 8, 61, 63, 57),
    ]
    observed = [
        (
            int(row["ell"]),
            int(row["I"]),
            int(row["E"]),
            int(row["xnor_l1_lower_bound"]),
            int(row["transverse_edge_capacity"]),
        )
        for row in rows
    ]
    if observed != expected or not all(row["excluded_by_l1"] for row in rows):
        raise ArithmeticError("p=17 u=3 coefficient exclusion changed")
    return {
        "u": u,
        "q": q,
        "j_box": "1<=j<=7",
        "congruence": "q divides u-j",
        "j_candidates": j_candidates,
        "rows": rows,
        "excluded": True,
    }


def pair_slack_divisibility() -> dict[str, object]:
    """Record why pair-budget slack is a nonnegative multiple of four."""
    return {
        "line_contribution": "2*(C(n,2)-floor(n/2))",
        "n=2r": "4*r*(r-1)",
        "n=2r+1": "4*r*r",
        "global_slack_modulus": 4,
        "zero_slack_iff": "every affine line occupancy is at most two",
        "zero_slack_consequence": "the finite point set is an arc",
        "proved": True,
    }


def endpoint_profiles() -> dict[str, object]:
    """Derive the exact two u0=0 profiles after pair slack divisibility."""
    # Phase one costs at least 96 and phase zero u=0 costs at least 84.
    phase_zero = _profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 96)
    phase_one = _profile_rows(1, 8, PAIR_DEFICIT_BUDGET - 84)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            total = deficit_zero + deficit_one
            slack = PAIR_DEFICIT_BUDGET - total
            if slack < 0 or slack % 4:
                continue
            by_phase = {
                "0": _histogram(profile_zero),
                "1": _histogram(profile_one),
            }
            secants = Counter()
            for profile in (profile_zero, profile_one):
                secants.update((S - b) // 2 for b in profile)
            candidates.append(
                {
                    "phase_profiles_b": by_phase,
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
    expected = [
        {
            "0": {"0": 6, "14": 3},
            "1": {"2": 8, "12": 1},
        },
        {
            "0": {"0": 6, "12": 1, "14": 2},
            "1": {"2": 8, "14": 1},
        },
    ]
    observed = [row["phase_profiles_b"] for row in candidates]
    if observed != expected:
        raise ArithmeticError("p=17 endpoint profile list changed")
    expected_secants = {"0": 3, "1": 1, "6": 8, "7": 6}
    if not all(
        row["total_deficit"] == PAIR_DEFICIT_BUDGET
        and row["arc"] is True
        and row["global_secant_distribution"] == expected_secants
        for row in candidates
    ):
        raise ArithmeticError("p=17 endpoint arc consequence changed")
    return {
        "p": P,
        "s": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "pair_slack_divisibility": pair_slack_divisibility(),
        "profiles": {
            "A": candidates[0],
            "B": candidates[1],
        },
        "common_global_secant_distribution": expected_secants,
        "common_undetermined_directions": 3,
        "all_profiles_are_arcs": True,
        "proved": True,
    }


def p17_arc_classification_ledger() -> dict[str, object]:
    """Record the external finite classification used in the endpoint proof."""
    classes = {"14": 4, "15": 1, "16": 1, "17": 1, "18": 1}
    return {
        "external_dependency": True,
        "source": (
            "H. Sticker, Classification of Arcs in Small Desarguesian "
            "Projective Planes, PhD thesis, Ghent University, 2012"
        ),
        "source_url": (
            "https://cage.ugent.be/geometry/Theses/57/PhDHeideSticker.pdf"
        ),
        "location": "Section 5.3, printed page 119 (PDF page 129)",
        "classification_scope": (
            "PGL-inequivalent (k,2)-arcs in PG(2,q), not necessarily complete"
        ),
        "pgl_classes_in_pg2_17": classes,
        "independent_generation_check": (
            "the thesis reports the orbit-stabilizer/double-count consistency "
            "check for every (k,2)-arc classification with q<=27"
        ),
        "known_representative": (
            "a nondegenerate conic has 18 points; deleting two gives a 16-arc"
        ),
        "unique_class_consequence": (
            "every 16-arc in PG(2,17) is projectively equivalent to conic-minus-two"
        ),
        "every_16_arc_is_conic_contained": classes["16"] == 1,
        "imported_result_accepted": True,
    }


def conic_secant_survival_ledger() -> dict[str, object]:
    """Count conic secants through an off-conic point after four deletions."""
    q = P
    removed = 4
    external_secants = (q - 1) // 2
    internal_secants = (q + 1) // 2
    external_remaining = external_secants - removed
    internal_remaining = internal_secants - removed
    if (external_remaining, internal_remaining) != (4, 5):
        raise ArithmeticError("p=17 conic secant lower bound changed")
    return {
        "q": q,
        "conic_size": q + 1,
        "points_removed_from_conic": removed,
        "external_point": {
            "tangents": 2,
            "conic_secants": external_secants,
            "remaining_S_secants_at_least": external_remaining,
        },
        "internal_point": {
            "tangents": 0,
            "conic_secants": internal_secants,
            "remaining_S_secants_at_least": internal_remaining,
        },
        "reason_deleting_one_point_destroys_at_most_one_secant": (
            "each deleted conic point lies on one line through the fixed point"
        ),
        "off_conic_point_is_never_undetermined": True,
        "proved": True,
    }


def three_undetermined_direction_contradiction() -> dict[str, object]:
    """Close either endpoint profile using the unique 16-arc class."""
    profiles = endpoint_profiles()
    classification = p17_arc_classification_ledger()
    secants = conic_secant_survival_ledger()
    excluded = bool(
        profiles["common_undetermined_directions"] == 3
        and classification["every_16_arc_is_conic_contained"] is True
        and secants["off_conic_point_is_never_undetermined"] is True
    )
    if not excluded:
        raise ArithmeticError("three-direction conic contradiction changed")
    return {
        "initial_arc_size": 14,
        "undetermined_points_on_line_at_infinity": 3,
        "extension": (
            "adjoin any two undetermined infinity points; the result is a 16-arc"
        ),
        "classification": classification,
        "extension_lies_on_conic": True,
        "conic_minus_original_arc_size": 4,
        "line_at_infinity_meets_conic_in_the_two_adjoined_points": True,
        "third_undetermined_point_is_off_conic": True,
        "secant_survival": secants,
        "contradiction": (
            "the third point has at least four S-secants but was undetermined"
        ),
        "excluded": True,
    }


def theorem_record() -> dict[str, object]:
    residues = type_residue_ledger()
    u2 = u2_lift_exclusion()
    u3 = u3_coefficient_exclusion()
    profiles = endpoint_profiles()
    geometry = three_undetermined_direction_contradiction()
    phase_zero_rows = residues["phase_residue_rows"]["0"]  # type: ignore[index]
    later_residues_over_budget = all(
        int(phase_zero_rows[str(u)]["minimum_deficit"]) + 96
        > PAIR_DEFICIT_BUDGET
        for u in range(4, M)
    )
    proved = bool(
        u2["excluded"] is True
        and u3["excluded"] is True
        and later_residues_over_budget
        and profiles["all_profiles_are_arcs"] is True
        and geometry["excluded"] is True
    )
    return {
        "prop": "15.678",
        "title": "Exceptional p=17 first all-finite survivor exclusion",
        "proved": proved,
        "theorem": {
            "p17_first_all_finite_survivor": "EXCLUDED_HERE",
            "all_odd_primes_p_at_least_17": "FIRST_SURVIVOR_EXCLUDED",
            "later_all_finite_boundary_sizes": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "residue_ledger": residues,
        "u2_lift_exclusion": u2,
        "u3_coefficient_exclusion": u3,
        "later_phase_zero_residues_over_pair_budget": later_residues_over_budget,
        "endpoint_profiles": profiles,
        "geometry_exclusion": geometry,
        "external_dependency_is_explicit": True,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.678 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15678.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.678 p=17 first all-finite survivor: EXCLUDED")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
