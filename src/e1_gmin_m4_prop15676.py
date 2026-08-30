#!/usr/bin/env python3
"""Prop. 15.676 -- close pair-deficit equality in infinity plus p.

Let the odd-degree boundary contain infinity and ``p`` finite points.  The
pair-deficit inequality is

    sum_d (p-b_d) <= p(p-1).

At equality the finite set is a ``p``-arc.  Segre's odd-order theorem puts
it on a nondegenerate conic.  A secant line at infinity leaves only ``p-1``
affine conic points.  A tangent line gives the profile

    p directions with b=1, one direction with b=p,

and an external line, after deleting one of the ``p+1`` affine conic
points, gives

    m+1 directions with b=1, m-1 directions with b=3,

where ``m=(p+1)/2``.

The external case fails the exact type budgets immediately.  In phase zero
one ``b=3`` direction already exceeds a type budget; in phase one each type
can contain at most one, whereas there are ``m-1>=8`` globally.  In the
tangent case, baseline ``b=1`` coefficient congruences force both baseline
parallel-edge counts to vanish.  The resulting zero or two finite edges
violate the boundary-support inequality.  Thus pair-deficit equality is
impossible for both signs and every prime ``p>=17``.  This does not close
the strict pair-deficit branch of the infinity-plus-p shell.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15723 import floor_excess_admissible


ROOT = Path(__file__).resolve().parents[1]


def _check_parameters(p: int, phase: int | None = None) -> None:
    if p < 17 or p % 2 == 0:
        raise ValueError("need an odd prime p>=17")
    if phase is not None and phase not in (0, 1):
        raise ValueError("phase must be zero or one")


def legendre(value: int, p: int) -> int:
    value %= p
    if value == 0:
        return 0
    return 1 if pow(value, (p - 1) // 2, p) == 1 else -1


def odd_fibre_profile(p: int, points: list[tuple[int, int]]) -> dict[str, object]:
    """Return the exact directional odd-fibre profile of an affine set."""
    rows = []
    for direction in projective_directions(p):
        eps, labels = field_direction_data(p, direction)
        occupancy = [0] * p
        for x, y in points:
            occupancy[labels[y * p + x]] += 1
        b = sum(value & 1 for value in occupancy)
        rows.append(
            {
                "direction": list(direction),
                "eps": eps,
                "b": b,
                "deficit": len(points) - b,
                "maximum_fibre_size": max(occupancy),
            }
        )
    return {
        "size": len(points),
        "profile": {
            str(b): count
            for b, count in sorted(Counter(row["b"] for row in rows).items())
        },
        "pair_deficit": sum(int(row["deficit"]) for row in rows),
        "maximum_fibre_size": max(int(row["maximum_fibre_size"]) for row in rows),
        "rows": rows,
    }


def canonical_conic_profile_audit(p: int) -> dict[str, object]:
    """Construct one tangent and one external-line conic representative."""
    _check_parameters(p)
    m = (p + 1) // 2

    # y=x^2 has one point at infinity, so all p affine conic points remain.
    tangent_points = [(x, x * x % p) for x in range(p)]
    tangent = odd_fibre_profile(p, tangent_points)

    # x^2-nu*y^2=z^2 has no point on z=0 when nu is a nonsquare.  Its
    # affine norm-one conic has p+1 points; delete (1,0) to obtain a p-arc.
    nonsquare = next(value for value in range(2, p) if legendre(value, p) == -1)
    full_external = [
        (x, y)
        for y in range(p)
        for x in range(p)
        if (x * x - nonsquare * y * y - 1) % p == 0
    ]
    if len(full_external) != p + 1 or (1, 0) not in full_external:
        raise ArithmeticError("canonical external conic count changed")
    external_points = [point for point in full_external if point != (1, 0)]
    external = odd_fibre_profile(p, external_points)

    expected_tangent = {"1": p, str(p): 1}
    expected_external = {"1": m + 1, "3": m - 1}
    pair_budget = p * (p - 1)
    if tangent["profile"] != expected_tangent:
        raise ArithmeticError("tangent-conic odd-fibre profile changed")
    if external["profile"] != expected_external:
        raise ArithmeticError("external-conic odd-fibre profile changed")
    if tangent["pair_deficit"] != pair_budget:
        raise ArithmeticError("tangent conic did not attain pair equality")
    if external["pair_deficit"] != pair_budget:
        raise ArithmeticError("punctured external conic did not attain equality")

    return {
        "p": p,
        "m": m,
        "pair_deficit_budget": pair_budget,
        "tangent_line_at_infinity": {
            "profile": tangent["profile"],
            "pair_deficit": tangent["pair_deficit"],
            "maximum_fibre_size": tangent["maximum_fibre_size"],
        },
        "external_line_at_infinity": {
            "nonsquare": nonsquare,
            "full_affine_conic_size": len(full_external),
            "deleted_point": [1, 0],
            "profile": external["profile"],
            "pair_deficit": external["pair_deficit"],
            "maximum_fibre_size": external["maximum_fibre_size"],
        },
        "proved_for_representatives": True,
    }


def conic_classification_ledger() -> dict[str, object]:
    """Record the projective classification used after Segre's theorem."""
    return {
        "pair_deficit_equality": (
            "every affine fibre has size at most two, so the p finite "
            "points form a p-arc"
        ),
        "Segre": "every odd-order p-arc is contained in a conic",
        "secant_infinity": "only p-1 affine conic points; impossible",
        "tangent_infinity": {
            "affine_points": "p, hence the boundary is the full affine conic",
            "profile": "p copies of b=1 and one copy of b=p",
        },
        "external_infinity": {
            "affine_points": "p+1, hence the boundary deletes one point",
            "direction_types": (
                "the external infinity line has m internal and m external "
                "points; deletion leaves m+1 copies of b=1 and m-1 of b=3"
            ),
            "character_sum_check": (
                "in x^2-nu*y^2=z^2, sum_r chi(r^2-nu)=-1 and the infinity "
                "direction contributes +1, so the two direction classes split m,m"
            ),
        },
        "exhaustive": True,
    }


def external_conic_floor_exclusion(p: int, phase: int) -> dict[str, object]:
    """Exclude the conic-minus-point profile using only type floor budgets."""
    _check_parameters(p, phase)
    period = p + 1
    m = period // 2
    b1_floor = full_symbolic_floor(p, 1, phase)
    b3_floor = full_symbolic_floor(p, 3, phase)
    if phase == 0:
        if (b1_floor, b3_floor) != (period, 2 * period - 8):
            raise ArithmeticError("phase-zero b=1,3 floors changed")
        maximum_b3_per_type = 0
        type_floor_with_t_b3 = "mP+t(P-8), so t=0"
    else:
        if (b1_floor, b3_floor) != (period - 2, 2 * period - 2):
            raise ArithmeticError("phase-one b=1,3 floors changed")
        maximum_b3_per_type = 1
        type_floor_with_t_b3 = "mP+(t-1)P, so t<=1"
    global_b3 = m - 1
    capacity = 2 * maximum_b3_per_type
    excluded = global_b3 > capacity
    if not excluded:
        raise ArithmeticError("external conic entered the type budgets")
    return {
        "p": p,
        "phase": phase,
        "P": period,
        "m": m,
        "floors": {"b=1": b1_floor, "b=3": b3_floor},
        "type_floor_formula": type_floor_with_t_b3,
        "global_b3_directions": global_b3,
        "maximum_b3_per_type": maximum_b3_per_type,
        "two_type_capacity": capacity,
        "excluded": True,
    }


def tangent_type_residue_audit(
    p: int, phase: int, contains_bp_direction: bool
) -> dict[str, object]:
    """Enumerate exact common residues for one tangent-conic type.

    The audit uses only the exact type equation and the theorem that a
    nonzero lift cannot have scaled cost two.  It records the minimum number
    of genuine baseline ``b=1`` directions in every feasible allocation.
    """
    _check_parameters(p, phase)
    period = p + 1
    m = period // 2
    profile = [1] * (m - int(contains_bp_direction))
    if contains_bp_direction:
        profile.append(p)
    baseline_mean = period if phase == 0 else period - 2
    rows = []
    for u in range(m):
        residue = 2 * u
        quotient_sum = m - u
        states = {(0, 0)}
        for b in profile:
            floor = full_symbolic_floor(p, b, phase)
            options = []
            for k in range(quotient_sum + 1):
                mean = residue + period * k
                excess = mean - floor
                if floor_excess_admissible(p, b, phase, excess):
                    options.append((k, int(b == 1 and mean == baseline_mean)))
            states = {
                (used + k, baselines + is_baseline)
                for used, baselines in states
                for k, is_baseline in options
                if used + k <= quotient_sum
            }
        baseline_counts = sorted(
            baselines for used, baselines in states if used == quotient_sum
        )
        if baseline_counts:
            rows.append(
                {
                    "u": u,
                    "residue": residue,
                    "minimum_b1_baselines": min(baseline_counts),
                    "maximum_b1_baselines": max(baseline_counts),
                }
            )
    expected_u = 0 if phase == 0 else m - 1
    if p % 4 == 1 and phase == 0 and contains_bp_direction:
        if rows:
            raise ArithmeticError("over-budget tangent type became feasible")
        minimum_baselines = None
    else:
        if [row["u"] for row in rows] != [expected_u]:
            raise ArithmeticError("tangent type acquired another residue")
        minimum_baselines = int(rows[0]["minimum_b1_baselines"])
        if minimum_baselines < m - 2:
            raise ArithmeticError("tangent type lost its b=1 baseline")
    return {
        "p": p,
        "phase": phase,
        "contains_b=p": contains_bp_direction,
        "profile": {"b=1": profile.count(1), "b=p": profile.count(p)},
        "feasible_residue_rows": rows,
        "minimum_b1_baselines": minimum_baselines,
        "aggregate_parallel_offset_from_baseline": (
            None if not rows else phase
        ),
    }


def tangent_conic_exclusion(p: int, phase: int) -> dict[str, object]:
    """Exclude the full affine tangent conic by coefficient arithmetic."""
    _check_parameters(p, phase)
    period = p + 1
    m = period // 2
    q = (p - 1) // 2
    sigma = 1 if phase == 0 else -1
    b1_floor = full_symbolic_floor(p, 1, phase)
    bp_floor = full_symbolic_floor(p, p, phase)
    expected_b1 = period if phase == 0 else period - 2
    if b1_floor != expected_b1:
        raise ArithmeticError("tangent b=1 floor changed")
    ordinary_type = tangent_type_residue_audit(p, phase, False)
    exceptional_type = tangent_type_residue_audit(p, phase, True)

    # For p=1 mod 4 in phase zero, the type containing b=p is already over
    # budget: (m-1)P+(2P-2)=mP+(P-2).
    if p % 4 == 1 and phase == 0:
        floor_sum = (m - 1) * b1_floor + bp_floor
        excluded = floor_sum > m * period
        if bp_floor != 2 * period - 2 or not excluded:
            raise ArithmeticError("tangent floor contradiction changed")
        return {
            "p": p,
            "phase": phase,
            "branch": "p=1 mod 4, phase zero",
            "floors": {"b=1": b1_floor, "b=p": bp_floor},
            "type_containing_b=p_floor_sum": floor_sum,
            "type_budget": m * period,
            "ordinary_type_residue_audit": ordinary_type,
            "exceptional_type_residue_audit": exceptional_type,
            "excluded": True,
            "method": "type floor exceeds budget",
        }

    expected_bp = 0 if (p % 4, phase) in ((3, 0), (1, 1)) else 2 * period - 2
    if bp_floor != expected_bp:
        raise ArithmeticError("tangent b=p floor changed")

    # Exact common residues leave a b=1 baseline in each type.  Relative to
    # that baseline, the sum of P_d offsets in a type is zero in phase zero
    # and one in phase one, directly from the exact type mean sum mP.
    offset_per_type = phase
    total_offset = 2 * offset_per_type
    candidates = []
    for x in range(8):
        for y in range(8 - x):
            finite_edges = m * (x + y) + total_offset
            infinity_edges = 4 * p + 1 - finite_edges
            if infinity_edges < 1:
                continue
            x_congruence = (infinity_edges + x - (4 + sigma)) % q == 0
            y_congruence = (infinity_edges + y - (4 + sigma)) % q == 0
            if x_congruence and y_congruence:
                candidates.append(
                    {
                        "x": x,
                        "y": y,
                        "E": finite_edges,
                        "I": infinity_edges,
                        "support_upper": p + 2 * finite_edges,
                        "support_contradiction": infinity_edges > p + 2 * finite_edges,
                    }
                )
    if len(candidates) != 1 or candidates[0]["x"] or candidates[0]["y"]:
        raise ArithmeticError("tangent coefficient candidates changed")
    if not candidates[0]["support_contradiction"]:
        raise ArithmeticError("tangent support contradiction disappeared")
    return {
        "p": p,
        "phase": phase,
        "branch": f"p={p % 4} mod 4, phase {phase}",
        "floors": {"b=1": b1_floor, "b=p": bp_floor},
        "same_type_residue_ledger": {
            "baseline_mean": expected_b1,
            "baseline_b1_directions_per_type_at_least": m - 2,
            "parallel_count_offset_per_type": offset_per_type,
            "finite_edge_formula": f"E=m(x+y)+{total_offset}",
            "proof": (
                "the common-residue equation and forbidden two-unit lift leave "
                "u=0 in phase zero or u=m-1 in phase one; the exact mean sum "
                "then fixes the aggregate P_d offset"
            ),
        },
        "ordinary_type_residue_audit": ordinary_type,
        "exceptional_type_residue_audit": exceptional_type,
        "baseline_congruence": "q divides I+P_d-(4+sigma)",
        "substituted_congruences": "q divides x and q divides y",
        "x_plus_y_upper_bound": 7,
        "candidates": candidates,
        "excluded": True,
        "method": "coefficient divisibility plus I<=p+2E",
    }


def theorem_record() -> dict[str, object]:
    sample_primes = (17, 19, 23, 29, 31, 37, 41, 101)
    profile_samples = {
        str(p): canonical_conic_profile_audit(p) for p in sample_primes[:-1]
    }
    branches = {
        str(p): {
            str(phase): {
                "tangent": tangent_conic_exclusion(p, phase),
                "external": external_conic_floor_exclusion(p, phase),
            }
            for phase in (0, 1)
        }
        for p in sample_primes
    }
    proved = all(
        branch["excluded"]
        for by_phase in branches.values()
        for row in by_phase.values()
        for branch in row.values()
    )
    return {
        "prop": "15.676",
        "title": "Pair-deficit equality exclusion in infinity plus p",
        "proved": proved,
        "theorem": {
            "boundary": "infinity plus p finite points",
            "pair_deficit_equality": "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS",
            "all_odd_primes_p_at_least_17": True,
            "strict_pair_deficit_branch": "OPEN",
            "whole_infinity_plus_p_shell": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "classification": conic_classification_ledger(),
        "canonical_profile_samples": profile_samples,
        "branch_samples": branches,
        "literature_scope": (
            "Segre supplies conic containment; the affine profiles are derived "
            "directly, and no searched source combines them with the Paley budgets"
        ),
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.676 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15676.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.676 infinity-plus-p pair-deficit equality: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
