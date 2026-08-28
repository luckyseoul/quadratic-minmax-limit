#!/usr/bin/env python3
"""Prop. 15.669 -- uniform boundary-range exclusion from parity floors.

For odd p>=17 and 5<=b<=p-5, the exact degree-two hypergeometric
parity-majorant floor is 2p in both phases. The proof is an explicit positive
quadrature: the hypergeometric first two moments lie between the lower and
upper convex envelopes of the parity-one integer nodes.

At residual size |H|=4p+1, each quadratic direction type has
m=(p+1)/2 directions and budget m(p+1). With no infinity in the boundary,
the two types have opposite phases. The middle floor, the exact low-b floors
of Proposition 15.652, and the pair-deficit inequality

    sum_d (s-b_d) <= s(s-1)

exclude every even finite boundary size 6<=s<=3(p-1)/4. With infinity
present, every odd number s of finite boundary points with 5<=s<=p-4 is
excluded. Exact small count-profile dynamic programs additionally exclude
infinity plus seven points at p=11, and at p=13 exclude eight finite points
as well as infinity plus seven or nine points.

This is a range theorem, not a closure of residual (ii). Larger boundary
sizes, the surviving small-prime profiles, R1, QVAR, Type I, and the limit
remain open.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15632 import hypergeometric_weights, scaled_direction_floor


ROOT = Path(__file__).resolve().parents[1]


def hypergeometric_moments(p: int, b: int) -> tuple[Fraction, Fraction, Fraction]:
    """Return E[t], E[t^2], and Var(t) exactly."""
    mean = Fraction(b * (p + 1), 2 * p)
    second = Fraction(b * (b + 1) * (p + 1), 4 * p)
    variance = Fraction((p + 1) * b * (p - b), 4 * p * p)
    if second - mean * mean != variance:
        raise ArithmeticError("hypergeometric moment formulas disagree")
    return mean, second, variance


def _add_weight(target: dict[int, Fraction], node: int, weight: Fraction) -> None:
    if weight:
        target[node] = target.get(node, Fraction()) + weight


def _two_point_mean_measure(
    left: int, right: int, mean: Fraction
) -> dict[int, Fraction]:
    """Probability measure on two endpoints with prescribed mean."""
    if left == right:
        if mean != left:
            raise ArithmeticError("degenerate endpoint measure has wrong mean")
        return {left: Fraction(1)}
    if not Fraction(left) <= mean <= Fraction(right):
        raise ArithmeticError("mean lies outside endpoint interval")
    width = right - left
    return {
        left: Fraction(right - mean, width),
        right: Fraction(mean - left, width),
    }


def middle_floor_quadrature(p: int, b: int, phase: int) -> dict[str, object]:
    """Explicit positive quadrature proving M(p,b,phase)=1."""
    if p < 17 or p % 2 == 0:
        raise ValueError("need odd p>=17")
    if not 5 <= b <= p - 5 or phase not in (0, 1):
        raise ValueError("need 5<=b<=p-5 and phase in {0,1}")

    m = (p + 1) // 2
    complemented = b > (p - 1) // 2
    b0 = p - b if complemented else b
    phase0 = phase ^ ((m & 1) if complemented else 0)
    mean, second, variance = hypergeometric_moments(p, b0)

    contact_parity = 1 - phase0
    left = contact_parity
    right = b0 if (b0 & 1) == contact_parity else b0 - 1
    if not left <= mean <= right:
        raise ArithmeticError("target mean is outside the contact hull")

    offset_steps = math.floor(Fraction(mean - left, 2))
    lower_left = left + 2 * offset_steps
    if mean == lower_left:
        lower_measure = {lower_left: Fraction(1)}
        lower_variance = Fraction()
    else:
        lower_right = lower_left + 2
        lower_measure = _two_point_mean_measure(lower_left, lower_right, mean)
        lower_variance = (mean - lower_left) * (lower_right - mean)

    upper_measure = _two_point_mean_measure(left, right, mean)
    upper_variance = (mean - left) * (right - mean)
    if not lower_variance <= variance <= upper_variance:
        raise ArithmeticError("target second moment is outside contact hull")
    mix_upper = (
        Fraction()
        if upper_variance == lower_variance
        else Fraction(variance - lower_variance, upper_variance - lower_variance)
    )

    reduced_weights: dict[int, Fraction] = {}
    for node, weight in lower_measure.items():
        _add_weight(reduced_weights, node, (1 - mix_upper) * weight)
    for node, weight in upper_measure.items():
        _add_weight(reduced_weights, node, mix_upper * weight)

    original_weights: dict[int, Fraction] = {}
    for node, weight in reduced_weights.items():
        original_node = m - node if complemented else node
        _add_weight(original_weights, original_node, weight)
    original_weights = dict(sorted(original_weights.items()))

    distribution = hypergeometric_weights(p, b)
    moment_match = all(
        sum(weight * Fraction(t**degree) for t, weight in distribution.items())
        == sum(
            weight * Fraction(node**degree)
            for node, weight in original_weights.items()
        )
        for degree in range(3)
    )
    contacts = all(
        node in distribution and ((node + phase) & 1) == 1
        for node in original_weights
    )
    exact = bool(
        all(weight >= 0 for weight in original_weights.values())
        and sum(original_weights.values(), Fraction()) == 1
        and moment_match
        and contacts
    )
    if not exact:
        raise ArithmeticError("constructed middle-floor quadrature failed")

    return {
        "p": p,
        "b": b,
        "phase": phase,
        "complemented": complemented,
        "reduced_b": b0,
        "reduced_phase": phase0,
        "contact_parity": contact_parity,
        "mean": mean,
        "second_moment": second,
        "variance": variance,
        "lower_variance": lower_variance,
        "upper_variance": upper_variance,
        "upper_envelope_mixing_weight": mix_upper,
        "quadrature_weights": original_weights,
        "candidate_coefficients": (Fraction(), Fraction(), Fraction(1)),
        "candidate_expectation": Fraction(1),
        "scaled_floor": 2 * p,
        "exact_positive_quadrature_certificate": True,
    }


def symbolic_floor_margin_ledger(p: int) -> dict[str, object]:
    """Finite algebraic ledger proving every contact-hull margin is positive."""
    if p < 17 or p % 2 == 0:
        raise ValueError("need odd p>=17")
    return {
        "reduced_interval": (5, (p - 1) // 2),
        "upper_elementary_lower_bounds": {
            "U1_b_ge_1": p + 1,
            "U2_b_ge_3": 2 * p + 6,
            "U3_b_ge_3": 2 * (p + 1),
            "U4_b_ge_4": p + 5,
        },
        "lower_elementary_lower_bounds": {
            "A_p_minus_b_minus_3": Fraction(p - 5, 2),
            "D_p_minus_b_times_b_minus_3": p + 1,
        },
        "B_concave_endpoint_numerators": (
            2 * p - 30,
            Fraction(p * p - 14 * p - 3, 4),
        ),
        "C_concave_endpoint_numerators": (
            2 * p - 30,
            Fraction(p * p - 14 * p + 1, 4),
        ),
        "all_strictly_positive": bool(
            2 * p - 30 > 0
            and p * p - 14 * p - 3 > 0
            and p * p - 14 * p + 1 > 0
        ),
    }


def low_floor(p: int, b: int, phase: int) -> int:
    """Exact Proposition 15.652 floors for b<=4."""
    table = {
        0: (0, 2 * p),
        1: (p + 1, p - 1),
        2: (p + 1, p - 1),
        3: (2 * p - 6, 2 * p),
        4: (2 * p - 6, 2 * p),
    }
    return table[b][phase]


def full_symbolic_floor(p: int, b: int, phase: int) -> int:
    """All exact floors for p>=17 from low, middle, and complement."""
    if p < 17 or p % 2 == 0 or not 0 <= b <= p or phase not in (0, 1):
        raise ValueError("need odd p>=17, 0<=b<=p, phase in {0,1}")
    if b <= 4:
        return low_floor(p, b, phase)
    if b <= p - 5:
        return 2 * p
    m = (p + 1) // 2
    return low_floor(p, p - b, phase ^ (m & 1))


def phase_zero_deficit_lower_bound(p: int, s: int) -> Fraction:
    """Required deficit in the phase-zero type for an even boundary."""
    if p < 17 or p % 2 == 0 or s < 6 or s % 2:
        raise ValueError("need odd p>=17 and even s>=6")
    return Fraction(phase_zero_deficit_ledger(p, s)["uniform_lower_bound"])


def phase_zero_deficit_ledger(p: int, s: int) -> dict[str, object]:
    """Exact two-branch knapsack dual behind the phase-zero bound."""
    if p < 17 or p % 2 == 0 or s < 6 or s % 2:
        raise ValueError("need odd p>=17 and even s>=6")
    required_saving = (p * p - 1) // 2
    zero_count = required_saving // (2 * p)
    remainder = required_saving - 2 * p * zero_count

    # b=0 has the best saving/deficit ratio. After its count is fixed,
    # b=2 dominates b=4. If at least zero_count+1 zero-fibre directions
    # occur, the first branch applies. Otherwise the second expression is
    # decreasing in their count and is minimized at zero_count.
    zero_beats_two_margin = s * (p + 1) - 4 * p
    two_beats_four_margin = s * (p - 7) - 4 * (p - 4)
    branch_extra_zero = Fraction((zero_count + 1) * s)
    branch_remainder = Fraction(zero_count * s) + Fraction(
        remainder * (s - 2), p - 1
    )
    uniform = Fraction((p + 1) * s, 4) - 1
    exact_residue_formula = (
        "p=1 mod 4: remainder=(p-1)/2; "
        "p=3 mod 4: remainder=(3p-1)/2"
    )
    valid = bool(
        zero_beats_two_margin > 0
        and two_beats_four_margin > 0
        and branch_extra_zero >= uniform
        and branch_remainder >= uniform
    )
    if not valid:
        raise ArithmeticError("phase-zero deficit ledger failed")
    return {
        "required_type_saving": required_saving,
        "full_zero_fibre_directions": zero_count,
        "remaining_saving": remainder,
        "zero_beats_two_ratio_margin": zero_beats_two_margin,
        "two_beats_four_ratio_margin": two_beats_four_margin,
        "extra_zero_branch_lower_bound": branch_extra_zero,
        "remainder_branch_lower_bound": branch_remainder,
        "uniform_lower_bound": uniform,
        "residue_formula": exact_residue_formula,
        "proved": True,
    }


def phase_one_deficit_lower_bound(p: int, s: int) -> int:
    """Required deficit in the phase-one type for an even boundary."""
    if p < 17 or p % 2 == 0 or s < 6 or s % 2:
        raise ValueError("need odd p>=17 and even s>=6")
    return (p - 1) * (s - 2) // 2


def no_infinity_range_exclusion(p: int, s: int) -> dict[str, object]:
    """Exclude an even all-finite boundary in the uniform range."""
    if p < 17 or p % 2 == 0 or s < 6 or s % 2:
        raise ValueError("need odd p>=17 and even s>=6")
    if 4 * s > 3 * (p - 1):
        raise ValueError("s lies beyond the proved no-infinity range")
    if s > p - 5:
        raise ArithmeticError("uniform range should lie inside middle floors")
    d0 = phase_zero_deficit_lower_bound(p, s)
    d1 = phase_one_deficit_lower_bound(p, s)
    required = d0 + d1
    pair_budget = s * (s - 1)
    gap = required - pair_budget
    numerator = s * (3 * p + 3 - 4 * s) - 4 * p
    if gap != Fraction(numerator, 4):
        raise ArithmeticError("deficit-gap simplification failed")
    return {
        "p": p,
        "finite_boundary_points": s,
        "infinity_in_boundary": False,
        "phase_zero_required_deficit_lower_bound": d0,
        "phase_one_required_deficit_lower_bound": d1,
        "required_total_deficit_lower_bound": required,
        "pair_deficit_budget": pair_budget,
        "contradiction_gap": gap,
        "gap_numerator": numerator,
        "excluded": gap > 0,
    }


def infinity_range_exclusion(p: int, s: int, phase: int) -> dict[str, object]:
    """Exclude infinity plus s finite points in the uniform range."""
    if p < 17 or p % 2 == 0 or s < 5 or s % 2 != 1 or phase not in (0, 1):
        raise ValueError("need odd p>=17, odd s>=5, phase in {0,1}")
    if s > p - 4:
        raise ValueError("s lies beyond the proved infinity-present range")
    type_ledger = infinity_type_deficit_ledger(p, s, phase)
    required = 2 * int(type_ledger["required_type_deficit_lower_bound"])
    pair_budget = s * (s - 1)
    return {
        "p": p,
        "finite_boundary_points": s,
        "total_boundary_points": s + 1,
        "infinity_in_boundary": True,
        "common_phase": phase,
        "type_deficit_ledger": type_ledger,
        "required_total_deficit_lower_bound": required,
        "pair_deficit_budget": pair_budget,
        "contradiction_gap": required - pair_budget,
        "excluded": required > pair_budget,
    }


def infinity_type_deficit_ledger(p: int, s: int, phase: int) -> dict[str, object]:
    """One-type saving/deficit lemma for an infinity-present boundary.

    Phase zero needs every direction to have ``b=1``.  Phase one needs at
    least ``m-1=(p-1)/2`` such directions.  The latter statement accounts
    explicitly for the complemented endpoint ``b=p-4``: when ``p=1 mod 4``
    it also saves six units, but even two non-``b=1`` directions leave a
    saving shortfall of at least ``p-11``.
    """
    if p < 17 or p % 2 == 0 or s < 5 or s % 2 != 1 or phase not in (0, 1):
        raise ValueError("need odd p>=17, odd s>=5, phase in {0,1}")
    if s > p - 4:
        raise ValueError("s lies beyond the proved infinity-present range")

    m = (p + 1) // 2
    required_saving = (p * p - 1) // 2
    allowed_b = tuple(range(1, s + 1, 2))
    savings = {b: 2 * p - full_symbolic_floor(p, b, phase) for b in allowed_b}
    b1_saving = p - 1 if phase == 0 else p + 1
    if savings[1] != b1_saving:
        raise ArithmeticError("b=1 saving disagrees with the symbolic floor")
    maximum_other_saving = max(
        (saving for b, saving in savings.items() if b != 1), default=0
    )
    if maximum_other_saving > 6:
        raise ArithmeticError("a non-b=1 direction saves more than six")

    if phase == 0:
        minimum_b1_count = m
        maximum_saving_with_one_replacement = (
            (m - 1) * (p - 1) + maximum_other_saving
        )
        saving_shortfall = required_saving - maximum_saving_with_one_replacement
        symbolic_shortfall_lower_bound = p - 7
    else:
        minimum_b1_count = m - 1
        maximum_saving_with_one_replacement = (
            (m - 2) * (p + 1) + 2 * maximum_other_saving
        )
        saving_shortfall = required_saving - maximum_saving_with_one_replacement
        symbolic_shortfall_lower_bound = p - 11

    if saving_shortfall < symbolic_shortfall_lower_bound or saving_shortfall <= 0:
        raise ArithmeticError("the forced b=1 count was not certified")
    required_deficit = minimum_b1_count * (s - 1)
    return {
        "direction_count": m,
        "required_type_saving": required_saving,
        "allowed_odd_fibre_counts": allowed_b,
        "savings_from_middle_baseline": savings,
        "b1_saving": b1_saving,
        "maximum_other_saving": maximum_other_saving,
        "maximum_saving_below_forced_b1_count": (
            maximum_saving_with_one_replacement
        ),
        "saving_shortfall": saving_shortfall,
        "symbolic_shortfall_lower_bound": symbolic_shortfall_lower_bound,
        "minimum_forced_b1_directions": minimum_b1_count,
        "required_type_deficit_lower_bound": required_deficit,
        "proved": True,
    }


def minimum_type_deficit_dp(p: int, s: int, phase: int) -> dict[str, object]:
    """Exact floor-only count-profile DP for one quadratic direction type."""
    if p < 11 or p % 2 == 0 or s < 1 or phase not in (0, 1):
        raise ValueError("need odd p>=11, s>=1, phase in {0,1}")
    direction_count = (p + 1) // 2
    type_budget = (p + 1) ** 2 // 2
    allowed_b = tuple(range(s & 1, min(s, p) + 1, 2))
    floors = {
        b: (
            full_symbolic_floor(p, b, phase)
            if p >= 17
            else scaled_direction_floor(p, b, phase)
        )
        for b in allowed_b
    }
    states: dict[int, int] = {0: 0}
    for _ in range(direction_count):
        next_states: dict[int, int] = {}
        for cost, deficit in states.items():
            for b in allowed_b:
                new_cost = cost + floors[b]
                if new_cost > type_budget:
                    continue
                new_deficit = deficit + s - b
                old = next_states.get(new_cost)
                if old is None or new_deficit < old:
                    next_states[new_cost] = new_deficit
        states = next_states
    if not states:
        raise ArithmeticError("type floor profile is itself infeasible")
    return {
        "p": p,
        "s": s,
        "phase": phase,
        "direction_count": direction_count,
        "type_budget": type_budget,
        "allowed_b": allowed_b,
        "floors": floors,
        "reachable_costs": len(states),
        "minimum_deficit": min(states.values()),
    }


def abstract_no_infinity_gap(p: int, s: int) -> int:
    """Exact floor-plus-pair gap for an all-finite count profile."""
    return (
        int(minimum_type_deficit_dp(p, s, 0)["minimum_deficit"])
        + int(minimum_type_deficit_dp(p, s, 1)["minimum_deficit"])
        - s * (s - 1)
    )


def abstract_infinity_gap(p: int, s: int, phase: int) -> int:
    """Exact floor-plus-pair gap for an infinity-present count profile."""
    return (
        2 * int(minimum_type_deficit_dp(p, s, phase)["minimum_deficit"])
        - s * (s - 1)
    )


def largest_even_in_general_range(p: int) -> int:
    value = 3 * (p - 1) // 4
    return value if value % 2 == 0 else value - 1


def range_gap_concavity_ledger(p: int) -> dict[str, object]:
    """Prove the no-infinity gap is positive on the whole real interval."""
    if p < 17 or p % 2 == 0:
        raise ValueError("need odd p>=17")
    # Four times the gap is h(s)=s(3p+3-4s)-4p, a concave quadratic.
    # Its minimum on [6,3(p-1)/4] occurs at an endpoint.
    lower_endpoint_numerator = 14 * p - 126
    upper_endpoint_numerator = Fraction(p - 9, 2)
    return {
        "gap_numerator": "h(s)=s(3p+3-4s)-4p",
        "second_derivative": -8,
        "lower_endpoint_s": 6,
        "lower_endpoint_numerator": lower_endpoint_numerator,
        "upper_endpoint_s": Fraction(3 * (p - 1), 4),
        "upper_endpoint_numerator": upper_endpoint_numerator,
        "positive_on_full_interval": bool(
            lower_endpoint_numerator > 0 and upper_endpoint_numerator > 0
        ),
    }


def small_prime_extensions() -> dict[str, object]:
    """Exact new p=11,13 cases beyond Proposition 15.657."""
    p11_inf7 = {phase: abstract_infinity_gap(11, 7, phase) for phase in (0, 1)}
    p13_no8 = abstract_no_infinity_gap(13, 8)
    p13_inf = {
        s: {phase: abstract_infinity_gap(13, s, phase) for phase in (0, 1)}
        for s in (7, 9)
    }
    first_survivors = {
        "p11_8_finite_gap": abstract_no_infinity_gap(11, 8),
        "p11_infinity_plus_9_gaps": {
            phase: abstract_infinity_gap(11, 9, phase) for phase in (0, 1)
        },
        "p13_10_finite_gap": abstract_no_infinity_gap(13, 10),
        "p13_infinity_plus_11_gaps": {
            phase: abstract_infinity_gap(13, 11, phase) for phase in (0, 1)
        },
    }
    return {
        "p11_infinity_plus_7_gaps": p11_inf7,
        "p11_infinity_plus_7_excluded_both_phases": all(
            gap > 0 for gap in p11_inf7.values()
        ),
        "p13_8_finite_gap": p13_no8,
        "p13_8_finite_excluded": p13_no8 > 0,
        "p13_infinity_gaps": p13_inf,
        "p13_infinity_plus_7_or_9_excluded_both_phases": all(
            gap > 0 for rows in p13_inf.values() for gap in rows.values()
        ),
        "first_floor_pair_survivors": first_survivors,
        "first_floor_pair_survivors_verified": all(
            gap <= 0
            for key, value in first_survivors.items()
            for gap in (value.values() if isinstance(value, dict) else (value,))
        ),
    }


def theorem_record() -> dict[str, object]:
    sample_primes = (17, 19, 23, 29, 31, 37, 47)
    floor_ledgers = {p: symbolic_floor_margin_ledger(p) for p in sample_primes}
    gap_ledgers = {p: range_gap_concavity_ledger(p) for p in sample_primes}
    quadrature_samples = {
        f"p={p},b={b},phase={phase}": middle_floor_quadrature(p, b, phase)
        for p in (17, 19, 23)
        for b in sorted({5, 6, (p - 1) // 2, p - 6, p - 5})
        for phase in (0, 1)
    }
    range_endpoints = {}
    for p in sample_primes:
        endpoint = largest_even_in_general_range(p)
        next_even = endpoint + 2
        range_endpoints[str(p)] = {
            "last_excluded_even_s": endpoint,
            "symbolic_gap_at_endpoint": no_infinity_range_exclusion(
                p, endpoint
            )["contradiction_gap"],
            "exact_floor_pair_gap_at_endpoint": abstract_no_infinity_gap(
                p, endpoint
            ),
            "first_floor_pair_survivor_s": next_even,
            "exact_floor_pair_gap_at_first_survivor": abstract_no_infinity_gap(
                p, next_even
            ),
            "last_excluded_infinity_finite_s": p - 4,
            "symbolic_infinity_endpoint_exclusions": {
                phase: infinity_range_exclusion(p, p - 4, phase)
                for phase in (0, 1)
            },
            "infinity_endpoint_gaps": {
                phase: abstract_infinity_gap(p, p - 4, phase)
                for phase in (0, 1)
            },
            "infinity_first_floor_pair_survivor_s": p - 2,
            "infinity_first_survivor_gaps": {
                phase: abstract_infinity_gap(p, p - 2, phase)
                for phase in (0, 1)
            },
        }
    small = small_prime_extensions()
    proved = bool(
        all(row["all_strictly_positive"] for row in floor_ledgers.values())
        and all(row["positive_on_full_interval"] for row in gap_ledgers.values())
        and all(
            row["exact_positive_quadrature_certificate"]
            for row in quadrature_samples.values()
        )
        and all(
            row["symbolic_gap_at_endpoint"] > 0
            and row["exact_floor_pair_gap_at_endpoint"] > 0
            and row["exact_floor_pair_gap_at_first_survivor"] <= 0
            and all(
                endpoint["excluded"]
                and endpoint["type_deficit_ledger"]["proved"]
                for endpoint in row[
                    "symbolic_infinity_endpoint_exclusions"
                ].values()
            )
            and all(gap > 0 for gap in row["infinity_endpoint_gaps"].values())
            and all(
                gap <= 0 for gap in row["infinity_first_survivor_gaps"].values()
            )
            for row in range_endpoints.values()
        )
        and small["p11_infinity_plus_7_excluded_both_phases"]
        and small["p13_8_finite_excluded"]
        and small["p13_infinity_plus_7_or_9_excluded_both_phases"]
        and small["first_floor_pair_survivors_verified"]
    )
    return {
        "prop": "15.669",
        "title": "Uniform non-Walsh boundary-range exclusion",
        "proved": proved,
        "theorem": {
            "middle_parity_floors_both_phases": (
                "M(p,b,phase)=1 and scaled floor=2p for odd p>=17, "
                "5<=b<=p-5"
            ),
            "no_infinity_general": (
                "every even s with 6<=s<=3(p-1)/4 is excluded for odd p>=17"
            ),
            "infinity_present_general": (
                "every odd finite s with 5<=s<=p-4 is excluded for odd p>=17"
            ),
            "small_prime_extensions": small,
            "pair_deficit_inequality": "sum_d(s-b_d)<=s(s-1)",
            "general_residual_ii": False,
            "all_non_Walsh_multilevel": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "symbolic_floor_margin_ledgers": floor_ledgers,
        "symbolic_range_gap_ledgers": gap_ledgers,
        "quadrature_samples": quadrature_samples,
        "range_endpoint_audit": range_endpoints,
        "remaining": {
            "p11": (
                "eight finite points or infinity plus nine, and larger, remain"
            ),
            "p13": "ten finite points or infinity plus eleven, and larger, remain",
            "p_at_least_17_no_infinity": (
                "smallest floor-plus-pair survivor is the first even integer "
                "strictly above 3(p-1)/4"
            ),
            "p_at_least_17_with_infinity": (
                "smallest floor-plus-pair survivor has p-2 finite points"
            ),
            "floor_pair_survivor_is_not_an_actual_graph": True,
        },
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


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.669 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15669.json"
    destination.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    print("Prop 15.669 uniform boundary-range exclusion: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
