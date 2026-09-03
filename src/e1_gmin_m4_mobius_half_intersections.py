#!/usr/bin/env python3
"""Exact intersections of two direction-localized Mobius halves.

The classification is symbolic.  After normalizing two distinct target
directions to coordinates X,Y, every allowed pair of auxiliaries has four
parameters q,r,A,B with A and B nonzero.  Each possible common inversion
orbit is forced by one of four sign-and-endpoint matchings, so no search over
the field is needed.

The resulting pairwise cancellation bound is sharp, but by itself does not
solve the fixed-edge objective or the restricted symmetric Boolean fibre.
"""

from __future__ import annotations

from typing import Literal

from e1_gmin_m4_prop15721 import is_prime


Point = tuple[int, int]
Edge = tuple[Point, Point]
Matching = Literal["direct", "swapped"]


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 5
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=5")


def _inverse(p: int, value: int) -> int:
    value %= p
    if value == 0:
        raise ZeroDivisionError("cannot invert zero")
    return pow(value, -1, p)


def _point_scale(p: int, sign: int, point: Point) -> Point:
    return (sign * point[0] % p, sign * point[1] % p)


def _edge(first: Point, second: Point) -> Edge:
    return tuple(sorted((first, second)))  # type: ignore[return-value]


def _first_half_edge(
    p: int, q: int, A: int, parameter: int
) -> Edge:
    """The first selected half in normalized X,Y coordinates."""
    t = parameter % p
    if t == p - 1:
        raise ValueError("the Mobius parameter cannot be -1")
    return _edge(
        (
            1,
            (q - A * _inverse(p, t + 1)) % p,
        ),
        (t, q * t % p),
    )


def _second_half_edge(
    p: int, r: int, B: int, parameter: int
) -> Edge:
    """The second selected half in normalized X,Y coordinates."""
    s = parameter % p
    if s == p - 1:
        raise ValueError("the Mobius parameter cannot be -1")
    return _edge(
        (
            (r - B * _inverse(p, s + 1)) % p,
            1,
        ),
        (r * s % p, s),
    )


def _candidate_parameters(
    p: int,
    q: int,
    r: int,
    A: int,
    B: int,
    sign: int,
    matching: Matching,
) -> tuple[int, int] | None:
    """Return the uniquely forced parameters for one matching, if defined."""
    if matching == "direct":
        q_denominator = (q - sign) % p
        r_denominator = (r - sign) % p
        if q_denominator == 0 or r_denominator == 0:
            return None
        return (
            (A * _inverse(p, q_denominator) - 1) % p,
            (B * _inverse(p, r_denominator) - 1) % p,
        )
    if matching == "swapped":
        if q % p == 0 or r % p == 0:
            return None
        return (
            sign * _inverse(p, q) % p,
            sign * _inverse(p, r) % p,
        )
    raise ValueError("matching must be direct or swapped")


def two_half_intersection_candidates(
    p: int,
    q: int,
    r: int,
    A: int,
    B: int,
) -> dict[str, object]:
    """Classify every common inversion orbit of two Mobius halves.

    The normalized auxiliary equations are

        M1/j1 = (1-q/A) X + (1/A) Y,
        M2/j2 = (1/B) X + (1-r/B) Y.

    A shared orbit has E1(t)=sign*E2(s), where sign is +1 or -1.
    For each sign there are only direct and swapped endpoint matchings.
    Two coordinates force t,s uniquely in each matching; the remaining two
    coordinates merely accept or reject that candidate.
    """
    _check_prime(p)
    q %= p
    r %= p
    A %= p
    B %= p
    if A == 0 or B == 0:
        raise ValueError("A and B must be nonzero")

    accepted: dict[tuple[int, Edge], dict[str, object]] = {}
    attempted: list[dict[str, object]] = []
    for sign in (1, -1):
        for matching in ("direct", "swapped"):
            parameters = _candidate_parameters(
                p, q, r, A, B, sign, matching
            )
            record: dict[str, object] = {
                "sign": sign,
                "matching": matching,
                "defined": parameters is not None,
                "accepted": False,
            }
            if parameters is None:
                attempted.append(record)
                continue
            t, s = parameters
            record["t"] = t
            record["s"] = s
            if t == p - 1 or s == p - 1:
                record["reason"] = "forced parameter equals excluded value -1"
                attempted.append(record)
                continue
            first = _first_half_edge(p, q, A, t)
            second = _second_half_edge(p, r, B, s)
            signed_second = _edge(
                _point_scale(p, sign, second[0]),
                _point_scale(p, sign, second[1]),
            )
            if first != signed_second:
                record["reason"] = "remaining endpoint equations reject"
                attempted.append(record)
                continue
            record["accepted"] = True
            record["edge"] = [list(point) for point in first]
            attempted.append(record)
            key = (sign, first)
            if key not in accepted:
                accepted[key] = record

    intersections = tuple(
        sorted(
            accepted.values(),
            key=lambda item: (
                int(item["sign"]),
                str(item["matching"]),
                int(item["t"]),
                int(item["s"]),
            ),
        )
    )
    same = tuple(item for item in intersections if item["sign"] == 1)
    opposite = tuple(item for item in intersections if item["sign"] == -1)
    if len(intersections) > 4 or len(same) > 2 or len(opposite) > 2:
        raise ArithmeticError("the four-candidate intersection bound changed")
    return {
        "p": p,
        "normalized_parameters": {"q": q, "r": r, "A": A, "B": B},
        "normalized_auxiliaries": (
            "M1/j1=(1-q/A)X+(1/A)Y; "
            "M2/j2=(1/B)X+(1-r/B)Y"
        ),
        "first_half": (
            "E1(t)={{(1,q-A/(t+1)),(t,q*t)}} for t!=-1"
        ),
        "second_half": (
            "E2(s)={{(r-B/(s+1),1),(r*s,s)}} for s!=-1"
        ),
        "candidate_formulas": {
            "direct": (
                "t=A/(q-sign)-1, s=B/(r-sign)-1"
            ),
            "swapped": "t=sign/q, s=sign/r",
        },
        "attempted_candidates": attempted,
        "intersections": list(intersections),
        "shared_inversion_orbits": len(intersections),
        "same_orientation_shared_orbits": len(same),
        "opposite_orientation_shared_orbits": len(opposite),
        "uniform_total_bound": 4,
        "uniform_opposite_orientation_bound": 2,
        "two_trade_sum_is_ternary": len(same) == 0,
        "proved": True,
    }


def sharp_two_cancellation_witness(p: int) -> dict[str, object]:
    """Return an all-prime pair attaining two opposite-orientation overlaps."""
    _check_prime(p)
    half = _inverse(p, 2)
    three_halves = 3 * half % p
    result = two_half_intersection_candidates(
        p,
        q=half,
        r=half,
        A=three_halves,
        B=three_halves,
    )
    opposite = result["opposite_orientation_shared_orbits"]
    same = result["same_orientation_shared_orbits"]
    parameters = {
        (int(item["t"]), int(item["s"]))
        for item in result["intersections"]
        if item["sign"] == -1
    }
    expected_parameters = {(0, 0), ((-2) % p, (-2) % p)}
    proved = bool(
        opposite == 2
        and same == 0
        and parameters == expected_parameters
        and result["two_trade_sum_is_ternary"]
    )
    if not proved:
        raise ArithmeticError("the sharp two-cancellation witness changed")
    return {
        "p": p,
        "q": half,
        "r": half,
        "A": three_halves,
        "B": three_halves,
        "normalized_common_auxiliary": "(2/3)*(X+Y)",
        "opposite_orientation_parameters": [
            [0, 0],
            [(-2) % p, (-2) % p],
        ],
        "cancelled_inversion_orbits": 2,
        "same_orientation_shared_orbits": 0,
        "two_trade_sum_is_ternary": True,
        "pairwise_bound_is_sharp": True,
        "proved": proved,
    }


def two_cancellation_locus_theorem(p: int) -> dict[str, object]:
    """Solve the locus where both opposite-orientation candidates occur.

    The opposite direct matching leaves

        t_d + r*s_d = 0,  q*t_d + s_d = 0.

    The opposite swapped matching requires q,r nonzero and not one, and

        A=(q-1)(q*r-1)/(q*r),
        B=(r-1)(q*r-1)/(q*r).

    Hence q*r is not one.  The direct equations force t_d=s_d=0, so
    A=q+1 and B=r+1.  Equating the formulas gives q=r=1-2*q*r, hence
    (2*q-1)(q+1)=0.  The root q=-1 makes A zero and is inadmissible.
    The locus is therefore the single normalized point
    q=r=1/2, A=B=3/2.
    """
    _check_prime(p)
    half = _inverse(p, 2)
    three_halves = 3 * half % p
    result = two_half_intersection_candidates(
        p,
        q=half,
        r=half,
        A=three_halves,
        B=three_halves,
    )
    accepted_opposite = tuple(
        item for item in result["intersections"] if item["sign"] == -1
    )
    proved = bool(
        len(accepted_opposite) == 2
        and result["same_orientation_shared_orbits"] == 0
        and {
            (int(item["t"]), int(item["s"]))
            for item in accepted_opposite
        }
        == {(0, 0), ((-2) % p, (-2) % p)}
    )
    if not proved:
        raise ArithmeticError("the rigid two-cancellation locus changed")
    return {
        "p": p,
        "opposite_direct_conditions": (
            "t_d+r*s_d=0 and q*t_d+s_d=0"
        ),
        "opposite_swapped_conditions": (
            "A=(q-1)(q*r-1)/(q*r), "
            "B=(r-1)(q*r-1)/(q*r)"
        ),
        "combined_factorization": "(2*q-1)*(q+1)=0 with q=r",
        "inadmissible_root": "q=r=-1 forces A=B=0",
        "unique_normalized_parameters": {
            "q": half,
            "r": half,
            "A": three_halves,
            "B": three_halves,
        },
        "unique_normalized_auxiliaries": (
            "M1/j1=M2/j2=(2/3)*(X+Y)"
        ),
        "same_orientation_matches_at_unique_point": 0,
        "free_parameter_after_two_cancellations": False,
        "greedy_pairing_from_free_locus_available": False,
        "proved": proved,
    }


def known_one_origin_witness(p: int) -> dict[str, object]:
    """Return the previously used common-origin auxiliary choice."""
    _check_prime(p)
    result = two_half_intersection_candidates(
        p, q=0, r=0, A=1, B=1
    )
    intersections = result["intersections"]
    proved = bool(
        result["opposite_orientation_shared_orbits"] == 1
        and result["same_orientation_shared_orbits"] == 0
        and len(intersections) == 1
        and intersections[0]["t"] == 0
        and intersections[0]["s"] == 0
    )
    if not proved:
        raise ArithmeticError("the known one-origin overlap changed")
    return {
        "p": p,
        "normalized_common_auxiliary": "X+Y",
        "cancelled_inversion_orbits": 1,
        "parameters": [0, 0],
        "two_trade_sum_is_ternary": True,
        "proved": proved,
    }


def branch_c_pairwise_cancellation_bound(p: int, t: int) -> dict[str, object]:
    """Compare the pairwise bound with the branch-C cancellation demand."""
    _check_prime(p)
    if p % 4 != 3 or p < 31:
        raise ValueError("need a branch-C prime p=3 mod 4 with p>=31")
    r0 = (p - 3) // 4
    hard_stars = 2 * r0 + 2
    t_min = 2 * r0 * r0 - 4 * r0 - 2
    t_max = 4 * r0 * r0 - 2 * r0 - 5
    if not t_min <= t <= t_max:
        raise ValueError("t is outside the balanced branch-C interval")

    raw_occurrences = hard_stars * (p - 1)
    target_edges = 4 * p + 2 * t + 1
    required_cancellation = t_max - t + 1
    pair_count = hard_stars * (hard_stars - 1) // 2
    pairwise_upper_bound = 2 * pair_count
    weakest_support_lower_bound = max(
        0, raw_occurrences - 2 * pairwise_upper_bound
    )
    objective_parity = (target_edges - raw_occurrences) % 2
    proved = bool(
        raw_occurrences == (p * p - 1) // 2
        and target_edges
        == raw_occurrences - (2 * (t_max - t) + 1)
        and required_cancellation == (raw_occurrences - target_edges + 1) // 2
        and objective_parity == 1
    )
    if not proved:
        raise ArithmeticError("the branch-C cancellation identities changed")
    return {
        "p": p,
        "r": r0,
        "t": t,
        "t_interval": [t_min, t_max],
        "hard_star_count": hard_stars,
        "raw_selected_orbit_occurrences": raw_occurrences,
        "target_edge_count": target_edges,
        "required_cancellation_units": required_cancellation,
        "opposite_overlap_bound_per_trade_pair": 2,
        "trade_pair_count": pair_count,
        "pairwise_cancellation_upper_bound": pairwise_upper_bound,
        "support_lower_bound_from_pairwise_theorem": weakest_support_lower_bound,
        "pairwise_bound_rules_out_required_cancellation": (
            pairwise_upper_bound < required_cancellation
        ),
        "exact_feasibility_objective": (
            "|U|+|a_Y+sum_(O in U)Phi(O)| <= |H|"
        ),
        "forced_fixed_weight_must_be_odd_if_feasible": True,
        "fixed_edge_objective_evaluated_by_intersection_theorem": False,
        "reason_fixed_edge_objective_remains_open": (
            "intersection data determines |U| but not the binary affine-line "
            "sum a_Y+sum Phi(O)"
        ),
        "restricted_symmetric_fibre_nonempty_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Package the exact theorem, sharpness, and deliberately open gate."""
    _check_prime(p)
    if p % 4 != 3 or p < 31:
        raise ValueError("the theorem record uses a branch-C prime p>=31")
    r0 = (p - 3) // 4
    t_min = 2 * r0 * r0 - 4 * r0 - 2
    lower = branch_c_pairwise_cancellation_bound(p, t_min)
    upper = branch_c_pairwise_cancellation_bound(
        p, int(lower["t_interval"][1])
    )
    return {
        "title": "Exact two-Mobius-half intersection classification",
        "known_origin_overlap": known_one_origin_witness(p),
        "sharp_two_cancellation_witness": sharp_two_cancellation_witness(p),
        "two_cancellation_locus": two_cancellation_locus_theorem(p),
        "branch_c_lower_endpoint": lower,
        "branch_c_upper_endpoint": upper,
        "status": (
            "PAIRWISE INTERSECTION THEOREM PROVED; "
            "FIXED-EDGE OBJECTIVE AND SYMMETRIC FIBRE OPEN"
        ),
        "proved": True,
        "fixed_edge_objective_closed": False,
        "residual_ii_closed": False,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
