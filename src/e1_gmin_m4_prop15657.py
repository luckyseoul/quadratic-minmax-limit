#!/usr/bin/env python3
"""Prop. 15.657 -- exclude every six-vertex boundary for p>=11.

At residual size ``|H|=4p+1``, Proposition 15.632 gives total affine
slack budget ``(p+1)^2`` and half that budget in each quadratic direction
type.  If a direction has finite boundary-fibre multiplicities ``n_i``,
put ``b_d=#{i:n_i odd}``.  Then

    s-b_d = 2 sum_i floor(n_i/2) <= 2 sum_i binom(n_i,2),

where ``s`` is the number of finite boundary points.  Every finite pair
collides in exactly one projective direction, so

    sum_d (s-b_d) <= s(s-1).

For a six-vertex boundary, ``s=6`` without infinity and ``s=5`` with
infinity.  Thus the respective total pair-deficit budgets are only 30 and
20.  Exact positive quadrature extends Proposition 15.652's parity-floor
table through six odd fibres.  Combining those floors with the pair-deficit
budgets exceeds the total affine slack budget for all odd ``p>=13`` and for
every infinity-present case already at ``p=11``.  In the remaining finite
``p=11`` case, the two quadratic types have opposite phases; their separate
budgets force deficits at least 20 and 18, contradicting the total bound 30.

Hence every six-point residual boundary is impossible for every odd prime
``p>=11``.  The exceptional primes ``p=5,7``, larger boundaries,
residual (ii), R1, and the limit remain open.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15632 import (
    eval_quadratic,
    hypergeometric_weights,
    scaled_direction_floor,
)
from e1_gmin_m4_prop15652 import small_boundary_floor_table

ROOT = Path(__file__).resolve().parents[1]


def _large_floor_quadrature_data(
    p: int, b: int, phase: int
) -> tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[int, int, int],
    tuple[Fraction, Fraction, Fraction],
]:
    """Candidate and exact dual quadrature for ``b in {5,6}``.

    The quadrature has the same moments through degree two as the relevant
    hypergeometric law and is supported on contacts of the candidate.
    Its nonnegative weights therefore certify optimality of the candidate.
    """
    if p < 11 or p % 2 == 0:
        raise ValueError("p must be odd and at least 11")
    if b not in (5, 6) or phase not in (0, 1):
        raise ValueError("need b in {5,6} and phase in {0,1}")

    if phase == 0:
        coefficients = (Fraction(0), Fraction(0), Fraction(1))
        nodes = (1, 3, 5)
        if b == 5:
            weights = (
                Fraction(5 * (p - 5), 16 * p),
                Fraction(5 * (p + 3), 8 * p),
                Fraction(p - 5, 16 * p),
            )
        else:
            weights = (
                Fraction(3 * (p - 9), 16 * p),
                Fraction(5 * (p + 3), 8 * p),
                Fraction(3 * (p - 1), 16 * p),
            )
        return coefficients, nodes, weights

    if p <= 15:
        coefficients = (Fraction(1), Fraction(-6), Fraction(9))
        nodes = (2, 3, 4)
        if b == 5:
            weights = (
                Fraction(p - 5, p),
                Fraction(15 - p, 2 * p),
                Fraction(p - 5, 2 * p),
            )
        else:
            weights = (
                Fraction(3 * (p - 7), 4 * p),
                Fraction(15 - p, 2 * p),
                Fraction(3 * (p - 3), 4 * p),
            )
        return coefficients, nodes, weights

    coefficients = (Fraction(0), Fraction(0), Fraction(1))
    nodes = (0, 2, 4)
    if b == 5:
        weights = (
            Fraction(p - 15, 16 * p),
            Fraction(5 * (p + 1), 8 * p),
            Fraction(5 * (p + 1), 16 * p),
        )
    else:
        weights = (
            Fraction(p - 15, 16 * p),
            Fraction(3 * (p + 1), 8 * p),
            Fraction(9 * (p + 1), 16 * p),
        )
    return coefficients, nodes, weights


def large_boundary_floor_certificate(p: int, b: int, phase: int) -> dict:
    """Verify one ``b=5,6`` symbolic floor by exact rational arithmetic."""
    coefficients, nodes, quadrature = _large_floor_quadrature_data(p, b, phase)
    distribution = hypergeometric_weights(p, b)
    support = tuple(distribution)
    candidate_expectation = sum(
        weight * eval_quadratic(coefficients, t)
        for t, weight in distribution.items()
    )
    quadrature_expectation = sum(
        weight * eval_quadratic(coefficients, node)
        for node, weight in zip(nodes, quadrature)
    )
    moment_match = all(
        sum(weight * Fraction(t**degree) for t, weight in distribution.items())
        == sum(
            weight * Fraction(node**degree)
            for node, weight in zip(nodes, quadrature)
        )
        for degree in range(3)
    )
    positive_nodes_are_contacts = all(
        weight == 0
        or (
            node in support
            and eval_quadratic(coefficients, node) == ((node + phase) & 1)
        )
        for node, weight in zip(nodes, quadrature)
    )
    majorizes = all(
        eval_quadratic(coefficients, t) >= ((t + phase) & 1)
        for t in support
    )
    expected = (
        Fraction(1)
        if phase == 0 or p >= 15
        else Fraction(3 * (p - 5), 2 * p)
    )
    exact = bool(
        all(weight >= 0 for weight in quadrature)
        and sum(quadrature) == 1
        and moment_match
        and positive_nodes_are_contacts
        and majorizes
        and candidate_expectation == quadrature_expectation == expected
    )
    return {
        "p": p,
        "b": b,
        "phase": phase,
        "coefficients": coefficients,
        "quadrature_nodes": nodes,
        "quadrature_weights": quadrature,
        "expectation": candidate_expectation,
        "scaled_floor": scaled_direction_floor(p, b, phase),
        "exact_positive_quadrature_certificate": exact,
    }


def size_six_floor_table(p: int) -> dict[int, tuple[int, int]]:
    """Return exact phase-zero/phase-one floors for ``b=0,...,6``."""
    if p < 11 or p % 2 == 0:
        raise ValueError("p must be odd and at least 11")
    small = small_boundary_floor_table(p)
    large = {
        (b, phase): large_boundary_floor_certificate(p, b, phase)
        for b in (5, 6)
        for phase in (0, 1)
    }
    if not all(
        row["exact_positive_quadrature_certificate"] for row in large.values()
    ):
        raise AssertionError("invalid large-boundary floor certificate")
    return {
        **small,
        5: (int(large[5, 0]["scaled_floor"]), int(large[5, 1]["scaled_floor"])),
        6: (int(large[6, 0]["scaled_floor"]), int(large[6, 1]["scaled_floor"])),
    }


def _integer_partitions(total: int, maximum: int | None = None):
    """Yield integer partitions in nonincreasing order."""
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for tail in _integer_partitions(total - first, first):
            yield (first, *tail)


def finite_fibre_partition_rows(s: int) -> tuple[dict, ...]:
    """Audit the local pair-deficit inequality for every partition of ``s``."""
    if s < 1:
        raise ValueError("s must be positive")
    rows = []
    for partition in _integer_partitions(s):
        b = sum(part & 1 for part in partition)
        deficit = s - b
        pair_collisions = sum(math.comb(part, 2) for part in partition)
        rows.append(
            {
                "partition": partition,
                "b": b,
                "deficit": deficit,
                "pair_collisions": pair_collisions,
                "local_bound": 2 * pair_collisions,
                "inequality_verified": deficit <= 2 * pair_collisions,
            }
        )
    return tuple(rows)


def pair_deficit_budget(s: int) -> dict:
    """Global deficit budget from unique pair directions in the affine plane."""
    rows = finite_fibre_partition_rows(s)
    return {
        "finite_boundary_points": s,
        "pair_count": math.comb(s, 2),
        "total_deficit_upper_bound": s * (s - 1),
        "all_local_partition_inequalities_verified": all(
            row["inequality_verified"] for row in rows
        ),
        "every_pair_has_one_projective_direction": True,
    }


def _count_vectors(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in _count_vectors(total - first, length - 1):
            yield (first, *tail)


def minimum_type_deficit(p: int, phase: int, s: int) -> dict:
    """Exact tiny count-profile LP for one quadratic direction type."""
    if p < 11 or p % 2 == 0 or phase not in (0, 1):
        raise ValueError("need odd p>=11 and phase in {0,1}")
    allowed_b = tuple(range(s & 1, s + 1, 2))
    floors = size_six_floor_table(p)
    direction_count = (p + 1) // 2
    budget = (p + 1) ** 2 // 2
    feasible = []
    for counts in _count_vectors(direction_count, len(allowed_b)):
        cost = sum(
            count * floors[b][phase] for b, count in zip(allowed_b, counts)
        )
        if cost <= budget:
            deficit = sum(
                count * (s - b) for b, count in zip(allowed_b, counts)
            )
            feasible.append((deficit, cost, counts))
    minimum = min(deficit for deficit, _cost, _counts in feasible)
    minimizers = [
        {"deficit": deficit, "cost": cost, "counts": dict(zip(allowed_b, counts))}
        for deficit, cost, counts in feasible
        if deficit == minimum
    ]
    return {
        "p": p,
        "phase": phase,
        "finite_boundary_points": s,
        "direction_count": direction_count,
        "type_budget": budget,
        "allowed_odd_fibre_counts": allowed_b,
        "feasible_count_profiles": len(feasible),
        "minimum_deficit": minimum,
        "minimizers": minimizers,
    }


def no_infinity_size_six_exclusion(p: int, c_h: int) -> dict:
    """Pair-deficit contradiction for six finite boundary points."""
    if p < 11 or p % 2 == 0 or c_h not in (-1, 1):
        raise ValueError("need odd p>=11 and c_h in {+-1}")
    deficit_budget = pair_deficit_budget(6)["total_deficit_upper_bound"]
    if p == 11:
        phase_zero = minimum_type_deficit(11, 0, 6)
        phase_one = minimum_type_deficit(11, 1, 6)
        required = phase_zero["minimum_deficit"] + phase_one["minimum_deficit"]
        return {
            "p": p,
            "c_H": c_h,
            "infinity_in_boundary": False,
            "phase_rule": "phase=1 iff eps_d*c_H=+1",
            "phase_zero_minimum_type_deficit": phase_zero["minimum_deficit"],
            "phase_one_minimum_type_deficit": phase_one["minimum_deficit"],
            "required_total_deficit": required,
            "pair_deficit_budget": deficit_budget,
            "contradiction_gap": required - deficit_budget,
            "excluded": required > deficit_budget,
        }

    total_budget = (p + 1) ** 2
    baseline_cost = (p + 1) * (2 * p - 2)
    # Relative to b=6, a deficit unit can save at most (p-1)/3:
    # b=4 saves <=4 for deficit 2, b=2 saves <=p-1 for deficit 4,
    # and b=0 saves <=2p-2 for deficit 6.  For p>=13 the last ratio
    # dominates, so the pair budget 30 saves at most 10(p-1).
    maximum_saving = 10 * (p - 1)
    lower_bound = baseline_cost - maximum_saving
    gap = lower_bound - total_budget
    return {
        "p": p,
        "c_H": c_h,
        "infinity_in_boundary": False,
        "minimum_floor_bounds": {0: 0, 2: p - 1, 4: 2 * p - 6, 6: 2 * p - 2},
        "pair_deficit_budget": deficit_budget,
        "baseline_cost": baseline_cost,
        "maximum_pair_deficit_saving": maximum_saving,
        "total_cost_lower_bound": lower_bound,
        "total_budget": total_budget,
        "contradiction_gap": gap,
        "gap_polynomial": "p^2-12p+7",
        "gap_increases_by": "4p-20 when p is replaced by p+2",
        "excluded": gap > 0,
    }


def infinity_size_six_exclusion(p: int, c_h: int) -> dict:
    """Pair-deficit contradiction for infinity plus five finite points."""
    if p < 11 or p % 2 == 0 or c_h not in (-1, 1):
        raise ValueError("need odd p>=11 and c_h in {+-1}")
    deficit_budget = pair_deficit_budget(5)["total_deficit_upper_bound"]
    total_budget = (p + 1) ** 2
    baseline_cost = (p + 1) * (2 * p - 4)
    # Relative to b=5, b=3 saves <=2 for deficit 2 and b=1 saves
    # <=p-3 for deficit 4.  The latter ratio dominates for p>=11.
    maximum_saving = 5 * (p - 3)
    lower_bound = baseline_cost - maximum_saving
    gap = lower_bound - total_budget
    return {
        "p": p,
        "c_H": c_h,
        "infinity_in_boundary": True,
        "phase": 1 if c_h == -1 else 0,
        "phase_independent_of_direction_type": True,
        "minimum_floor_bounds": {1: p - 1, 3: 2 * p - 6, 5: 2 * p - 4},
        "pair_deficit_budget": deficit_budget,
        "baseline_cost": baseline_cost,
        "maximum_pair_deficit_saving": maximum_saving,
        "total_cost_lower_bound": lower_bound,
        "total_budget": total_budget,
        "contradiction_gap": gap,
        "gap_polynomial": "p^2-9p+10",
        "gap_increases_by": "4p-14 when p is replaced by p+2",
        "excluded": gap > 0,
    }


def theorem_size_six_boundary_pge11() -> dict:
    """Machine-readable theorem with the exceptional-prime scope explicit."""
    sample_primes = (11, 13, 17, 19, 23, 29, 31, 37)
    quadrature_certificates = all(
        large_boundary_floor_certificate(p, b, phase)[
            "exact_positive_quadrature_certificate"
        ]
        for p in sample_primes
        for b in (5, 6)
        for phase in (0, 1)
    )
    partition_certificates = all(
        pair_deficit_budget(s)["all_local_partition_inequalities_verified"]
        for s in (5, 6)
    )
    p11_type_split = all(
        no_infinity_size_six_exclusion(11, c_h)["excluded"]
        for c_h in (-1, 1)
    )
    analytic_base_cases = bool(
        no_infinity_size_six_exclusion(13, 1)["contradiction_gap"] == 20
        and infinity_size_six_exclusion(11, 1)["contradiction_gap"] == 32
    )
    odd_step_monotonicity = bool(4 * 13 - 20 > 0 and 4 * 11 - 14 > 0)
    sample_audit = all(
        no_infinity_size_six_exclusion(p, c_h)["excluded"]
        and infinity_size_six_exclusion(p, c_h)["excluded"]
        for p in sample_primes
        for c_h in (-1, 1)
    )
    proved = bool(
        quadrature_certificates
        and partition_certificates
        and p11_type_split
        and analytic_base_cases
        and odd_step_monotonicity
        and sample_audit
    )
    return {
        "proved": proved,
        "exact_floor_formulas_through_six_odd_fibres": quadrature_certificates,
        "pair_deficit_inequality": "sum_d(s-b_d)<=s(s-1)",
        "p11_no_infinity_type_split_closed": p11_type_split,
        "pge13_no_infinity_gap": "p^2-12p+7>0",
        "pge11_with_infinity_gap": "p^2-9p+10>0",
        "six_point_boundary_all_odd_primes_p_at_least_11": "CLOSED",
        "depends_on_affine_budget_prop": "15.632",
        "depends_on_small_floor_prop": "15.652",
        "p5_size_six": "OPEN",
        "p7_size_six": "OPEN",
        "larger_boundary_sizes": "OPEN",
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    return value


def main() -> dict:
    theorem = theorem_size_six_boundary_pge11()
    out = {
        "prop": "15.657",
        "title": "Six-point residual boundary exclusion for p>=11",
        "proved": theorem["proved"],
        "theorem": theorem,
        "large_floor_certificates": {
            f"p={p},b={b},phase={phase}": _jsonable(
                large_boundary_floor_certificate(p, b, phase)
            )
            for p in (11, 13, 17)
            for b in (5, 6)
            for phase in (0, 1)
        },
        "partition_rows": {
            str(s): _jsonable(finite_fibre_partition_rows(s)) for s in (5, 6)
        },
        "p11_type_profiles": {
            f"s={s},phase={phase}": _jsonable(minimum_type_deficit(11, phase, s))
            for s in (5, 6)
            for phase in (0, 1)
        },
        "no_infinity_cases": {
            f"p={p},c_H={c_h}": no_infinity_size_six_exclusion(p, c_h)
            for p in (11, 13, 17)
            for c_h in (-1, 1)
        },
        "with_infinity_cases": {
            f"p={p},c_H={c_h}": infinity_size_six_exclusion(p, c_h)
            for p in (11, 13, 17)
            for c_h in (-1, 1)
        },
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15657.json"
    destination.write_text(json.dumps(_jsonable(out), indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
