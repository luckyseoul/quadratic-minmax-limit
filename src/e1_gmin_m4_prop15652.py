#!/usr/bin/env python3
"""Prop. 15.652 -- exclude every four-vertex boundary for p>=11.

At residual size ``|H|=4p+1``, Proposition 15.632 gives each quadratic
direction type the exact slack budget ``(p+1)^2/2``.  For an odd prime
``p>=7``, the parity-majorant floors needed by a boundary of size four are

=================  ========  ========
odd finite fibres  phase 0   phase 1
=================  ========  ========
0                 0         2p
1,2               p+1       p-1
3,4               2p-6      2p
=================  ========  ========

The formulas have exact positive quadrature certificates, recorded below;
they are not a numerical extrapolation of the LP tables.

If infinity is absent from the four-vertex boundary, the four finite points
create six unordered pairs.  Every pair collides in exactly one projective
direction.  The expensive quadratic type needs at least ``(p-1)/2``
directions with two odd fibres.  Six pairs exclude ``p>=17`` immediately;
the exact leftover accounting excludes ``p=11,13``.

If infinity is present, the three finite points create only three pair
collisions, so at most three directions have one odd fibre.  The floors
above exclude the negative edge-product sign for ``p>=7`` and the positive
sign for ``p>=11``.  Hence every size-four boundary is impossible for every
odd prime ``p>=11``.  Together with Propositions 15.632 and 15.650--15.651,
the first still-open boundary size there is at least six.  This does not
close the exceptional p=5,7 size-four cases, larger boundaries, residual
(ii), R1, or the limit.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15632 import (
    eval_quadratic,
    hypergeometric_weights,
    scaled_direction_floor,
)

ROOT = Path(__file__).resolve().parents[1]


def _quadrature_data(
    p: int, b: int, phase: int
) -> tuple[
    tuple[Fraction, Fraction, Fraction],
    tuple[int, ...],
    tuple[Fraction, ...],
]:
    """Candidate majorant and exact positive degree-two quadrature.

    The quadrature nodes are parity contacts of the candidate whenever the
    corresponding weight is positive.  Therefore every feasible quadratic
    has expectation at least the candidate's expectation.
    """
    if p < 7 or p % 2 == 0:
        raise ValueError("p must be odd and at least 7")
    if b not in range(5) or phase not in (0, 1):
        raise ValueError("need b in 0..4 and phase in {0,1}")

    if b == 0:
        return (
            (Fraction(0), Fraction(0), Fraction(phase)),
            (0,),
            (Fraction(1),),
        )
    if b == 1:
        coefficients = (
            (Fraction(0), Fraction(1), Fraction(0))
            if phase == 0
            else (Fraction(0), Fraction(-1), Fraction(1))
        )
        return coefficients, (0, 1), (
            Fraction(p - 1, 2 * p),
            Fraction(p + 1, 2 * p),
        )
    if b == 2:
        coefficients = (
            (Fraction(-1), Fraction(2), Fraction(0))
            if phase == 0
            else (Fraction(1), Fraction(-2), Fraction(1))
        )
        weights = hypergeometric_weights(p, 2)
        return coefficients, (0, 1, 2), tuple(weights[t] for t in range(3))
    if b == 3 and phase == 0:
        return (
            (Fraction(1), Fraction(-4), Fraction(4)),
            (1, 2, 3),
            (
                Fraction(3 * (p - 3), 4 * p),
                Fraction(3, p),
                Fraction(p - 3, 4 * p),
            ),
        )
    if b == 3 and phase == 1:
        return (
            (Fraction(1), Fraction(-2), Fraction(1)),
            (0, 1, 2),
            (
                Fraction(p - 3, 4 * p),
                Fraction(0),
                Fraction(3 * (p + 1), 4 * p),
            ),
        )
    if b == 4 and phase == 0:
        return (
            (Fraction(1), Fraction(-4), Fraction(4)),
            (1, 2, 3),
            (
                Fraction(p - 5, 2 * p),
                Fraction(3, p),
                Fraction(p - 1, 2 * p),
            ),
        )
    return (
        (Fraction(0), Fraction(0), Fraction(1)),
        (0, 2, 4),
        (
            Fraction(p - 7, 8 * p),
            Fraction(3 * (p + 1), 4 * p),
            Fraction(p + 1, 8 * p),
        ),
    )


def parity_floor_certificate(p: int, b: int, phase: int) -> dict:
    """Verify one symbolic floor formula by exact rational arithmetic."""
    coefficients, nodes, quadrature = _quadrature_data(p, b, phase)
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
        == sum(weight * Fraction(node**degree) for node, weight in zip(nodes, quadrature))
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
    expected = {
        0: (Fraction(0), Fraction(1)),
        1: (Fraction(p + 1, 2 * p), Fraction(p - 1, 2 * p)),
        2: (Fraction(p + 1, 2 * p), Fraction(p - 1, 2 * p)),
        3: (Fraction(p - 3, p), Fraction(1)),
        4: (Fraction(p - 3, p), Fraction(1)),
    }[b][phase]
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


def small_boundary_floor_table(p: int) -> dict[int, tuple[int, int]]:
    """Return the exact phase-zero/phase-one floors for b=0,...,4."""
    certificates = {
        (b, phase): parity_floor_certificate(p, b, phase)
        for b in range(5)
        for phase in (0, 1)
    }
    if not all(
        row["exact_positive_quadrature_certificate"]
        for row in certificates.values()
    ):
        raise AssertionError("invalid parity-floor quadrature certificate")
    return {
        b: (
            int(certificates[b, 0]["scaled_floor"]),
            int(certificates[b, 1]["scaled_floor"]),
        )
        for b in range(5)
    }


def four_finite_partition_rows() -> tuple[dict, ...]:
    """Odd-fibre and collision counts for partitions of four points."""
    return (
        {"partition": (1, 1, 1, 1), "b": 4, "pair_collisions": 0},
        {"partition": (2, 1, 1), "b": 2, "pair_collisions": 1},
        {"partition": (2, 2), "b": 0, "pair_collisions": 2},
        {"partition": (3, 1), "b": 2, "pair_collisions": 3},
        {"partition": (4,), "b": 0, "pair_collisions": 6},
    )


def three_finite_partition_rows() -> tuple[dict, ...]:
    """Odd-fibre and collision counts for partitions of three points."""
    return (
        {"partition": (1, 1, 1), "b": 3, "pair_collisions": 0},
        {"partition": (2, 1), "b": 1, "pair_collisions": 1},
        {"partition": (3,), "b": 1, "pair_collisions": 3},
    )


def no_infinity_size_four_exclusion(p: int) -> dict:
    """All-prime pair-direction contradiction when D has four finite points."""
    if p < 11 or p % 2 == 0:
        raise ValueError("this exclusion applies to odd p>=11")
    m = (p + 1) // 2
    q = (p - 1) // 2
    budget = m * (p + 1)
    expensive = 2 * p
    discounted = p - 1
    required_b2 = q
    assert (
        required_b2 * discounted + (m - required_b2) * expensive
        == budget
    )

    if p >= 15:
        reason = "required_b2_directions_exceed_six_pair_directions"
        contradiction_gap = required_b2 - 6
    elif p == 13:
        # Six required b=2 directions consume all six pairs.  Every good-type
        # direction then has b=4 and phase-zero cost 2p-6.
        good_type_lower_bound = m * (2 * p - 6)
        reason = "all_pairs_consumed_then_good_type_exceeds_budget"
        contradiction_gap = good_type_lower_bound - budget
    elif p == 11:
        # Five required b=2 directions leave at most one collision outside
        # them.  A good-type collision can only be a single 2+1+1 partition,
        # so at least five of its six directions have b=4.
        good_type_lower_bound = (m - 1) * (2 * p - 6) + (p + 1)
        reason = "at_most_one_good_collision_then_good_type_exceeds_budget"
        contradiction_gap = good_type_lower_bound - budget
    else:  # pragma: no cover - guarded above for odd p
        raise AssertionError("unreachable")
    return {
        "p": p,
        "infinity_in_boundary": False,
        "total_pair_collisions": 6,
        "required_bad_type_b2_directions": required_b2,
        "reason": reason,
        "contradiction_gap": contradiction_gap,
        "excluded": contradiction_gap > 0,
    }


def infinity_size_four_exclusion(p: int, c_h: int) -> dict:
    """Pair-direction contradiction for infinity plus three finite points."""
    if p < 7 or p % 2 == 0 or c_h not in (-1, 1):
        raise ValueError("need odd p>=7 and c_h in {+-1}")
    m = (p + 1) // 2
    budget = m * (p + 1)
    if c_h == -1:
        # Phase one: b=1 costs p-1 and b=3 costs 2p.  Each type can
        # therefore contain at most one b=3 direction.
        required_b1 = p - 1
        actual_b1_upper_bound = 3
        return {
            "p": p,
            "c_H": c_h,
            "phase": 1,
            "total_pair_collisions": 3,
            "required_b1_directions": required_b1,
            "actual_b1_upper_bound": actual_b1_upper_bound,
            "contradiction_gap": required_b1 - actual_b1_upper_bound,
            "excluded": required_b1 > actual_b1_upper_bound,
        }

    # Phase zero: b=1 costs p+1 and b=3 costs 2p-6.  For p>=11 every
    # b=3 direction is strictly more expensive, while the budget is exactly
    # the all-b=1 cost.  Hence all p+1 directions would have to collide.
    required_b1 = p + 1
    actual_b1_upper_bound = 3
    return {
        "p": p,
        "c_H": c_h,
        "phase": 0,
        "total_pair_collisions": 3,
        "required_b1_directions": required_b1 if p >= 11 else None,
        "actual_b1_upper_bound": actual_b1_upper_bound,
        "strict_b3_surcharge": p - 7,
        "contradiction_gap": (
            required_b1 - actual_b1_upper_bound if p >= 11 else None
        ),
        "excluded": p >= 11 and required_b1 > actual_b1_upper_bound,
        "p7_exception_remains": p == 7,
        "type_budget": budget,
    }


def theorem_size_four_boundary() -> dict:
    """Machine-readable statement, with exceptional-prime scope explicit."""
    sample_primes = (7, 11, 13, 17, 19, 23, 29, 31)
    floor_certificates = all(
        small_boundary_floor_table(p)
        == {
            0: (0, 2 * p),
            1: (p + 1, p - 1),
            2: (p + 1, p - 1),
            3: (2 * p - 6, 2 * p),
            4: (2 * p - 6, 2 * p),
        }
        for p in sample_primes
    )
    all_cases = all(
        no_infinity_size_four_exclusion(p)["excluded"]
        and infinity_size_four_exclusion(p, -1)["excluded"]
        and infinity_size_four_exclusion(p, 1)["excluded"]
        for p in sample_primes
        if p >= 11
    )
    return {
        "proved": bool(floor_certificates and all_cases),
        "exact_floor_formulas_all_odd_p_at_least_7": True,
        "four_point_boundary_all_odd_primes_p_at_least_11": "CLOSED",
        "infinity_plus_three_finite_c_H_minus_all_odd_p_at_least_7": "CLOSED",
        "first_open_boundary_size_for_p_at_least_11": 6,
        "depends_on_empty_boundary_prop": "15.632",
        "depends_on_two_point_boundary_props": "15.650--15.651",
        "p5_size_four": "OPEN",
        "p7_size_four": (
            "OPEN (the infinity-present c_H=-1 subbranch alone is closed)"
        ),
        "closes_larger_boundary_shapes": False,
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
    theorem = theorem_size_four_boundary()
    out = {
        "prop": "15.652",
        "title": "Four-point residual boundary exclusion for p>=11",
        "proved": theorem["proved"],
        "theorem": theorem,
        "floor_certificates": {
            str(p): {
                f"b={b},phase={phase}": _jsonable(
                    parity_floor_certificate(p, b, phase)
                )
                for b in range(5)
                for phase in (0, 1)
            }
            for p in (7, 11, 13)
        },
        "no_infinity_cases": {
            str(p): no_infinity_size_four_exclusion(p)
            for p in (11, 13, 17)
        },
        "with_infinity_cases": {
            f"p={p},c_H={c_h}": infinity_size_four_exclusion(p, c_h)
            for p in (7, 11, 13)
            for c_h in (-1, 1)
        },
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15652.json"
    destination.write_text(json.dumps(_jsonable(out), indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
