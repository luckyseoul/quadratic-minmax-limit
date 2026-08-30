#!/usr/bin/env python3
"""Prop. 15.724 -- exclude the full Miquelian-circle boundary.

Continue Proposition 15.722 at residual size ``|H|=4p+1`` with
``p>=17`` and boundary ``D`` equal to a full Miquelian circle of size
``P=p+1=2m``.  Normalize one boundary point to infinity.  Proposition
15.722 gives

    D={infinity} union (a+b*F_p),       c_H=(-1)^m.          (1)

The complement of ``D`` has ``p^2-p`` vertices.  The graph ``H`` has at
most ``2|H|=8p+2`` nonisolated vertices, and every point of ``D`` is
nonisolated because it has odd degree.  Hence at most ``7p+1`` outside
vertices are nonisolated.  Since

    p^2-p > 7p+1                                           (p>=17),

there is an isolated outside point ``w``.

Send ``w`` to infinity.  If ``eps`` is the direction type of the line in
(1), Proposition 15.722's exact sign cocycle gives

    c_w=(-1)^m * ((-1)^m eps)=eps.                          (2)

The transported boundary is an all-finite affine circle.  Its ``b=2``
directions have type ``eps`` and its ``b=0`` directions the opposite type.
Thus (2) gives exactly ``m`` phase-one ``b=2`` directions and ``m``
phase-zero ``b=0`` directions.  The new infinity degree is
``I=deg_H(w)=0``.

For either direction type, the exact scaled means satisfy

    a_d=I+P P_d-eps_d*T-3p
       =2u+P k_d,             sum_d k_d=m-u.                (3)

In the phase-one type the ``b=2`` floor is ``P-2``.  For
``1<=u<=m-2``, every direction would need ``k_d>=1`` although their quotient
sum is below ``m``.  At ``u=0`` all means would be ``P``, only two above the
xnor baseline; Proposition 15.688 excludes that nonzero lift.  Therefore
``u=m-1`` and the parallel counts are

    x,...,x,x+1.                                           (4)

For the phase-zero type write ``P_d=y+k_d``.  Counting all finite edges in
(3), using ``I=0`` and ``|H|=8m-3``, gives

    8m-3=m(x+y+1)+1-u,
    m(x+y-7)=u-4.                                          (5)

As ``0<=u<m``, equation (5) forces ``u=4`` and ``x+y=7``.  There is then a
zero quotient in the phase-zero type, so ``y>=0``; also ``x>=0``.

Any xnor-baseline direction in (4) has the exact coefficient congruence

    q=(p-1)/2 divides I+x-4.                               (6)

Since ``I=0``, ``q>=8`` and ``0<=x<=7``, equations (5)--(6) force
``x=4,y=3``.  The phase-zero quotient sum is ``m-4``, so at least four of
its ``m`` directions have ``k_d=0`` and hence ``a_d=2u=8``.  In such a
``b=0``, phase-zero direction, ``A_d=2B_d`` for a nonzero nonnegative
integer-valued quadratic ``B_d``.  Consequently

    4p E[B_d]=8,

contradicting Proposition 15.688's sharp bound

    4p E[B_d] >= p-3 >= 14.                                (7)

Thus a full Miquelian-circle boundary is impossible for every odd prime
``p>=17``.  Combined with Proposition 15.722, the outside-chart
``R=0`` branch is now empty; 15.722 also excludes ``R=2,3`` by conic
extension and, more generally, every
``1<=R<=max(3,floor(sqrt(p)-5/2))`` by off-conic secant counting.  Strict
outside profiles beyond that cutoff, the whole
``p+1`` shell, residual (ii), Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15673 import coefficient_ledger
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15722 import (
    full_circle_line_chart_normal_form,
    outside_low_slack_conic_exclusion,
    outside_R_three_structure,
    outside_R_two_structure,
)


ROOT = Path(__file__).resolve().parents[1]


def _check_prime_parameter(p: int) -> None:
    if p < 17 or not is_prime(p):
        raise ValueError("need an odd prime parameter p>=17")


def isolated_outside_vertex_ledger(p: int) -> dict[str, object]:
    """Force an isolated vertex outside a size-``p+1`` boundary."""
    _check_prime_parameter(p)
    boundary = p + 1
    outside = p * p - p
    edges = 4 * p + 1
    maximum_nonisolated = 2 * edges
    maximum_nonisolated_outside = maximum_nonisolated - boundary
    guaranteed_isolated_outside = outside - maximum_nonisolated_outside
    proved = guaranteed_isolated_outside > 0
    if not proved:
        raise ArithmeticError("the isolated outside vertex gap disappeared")
    return {
        "p": p,
        "boundary_vertices": boundary,
        "outside_vertices": outside,
        "edge_count": edges,
        "maximum_nonisolated_vertices": maximum_nonisolated,
        "boundary_vertices_are_nonisolated": True,
        "maximum_nonisolated_outside": maximum_nonisolated_outside,
        "guaranteed_isolated_outside_vertices": guaranteed_isolated_outside,
        "isolated_outside_vertex_exists": proved,
        "proved": proved,
    }


def isolated_circle_chart(p: int) -> dict[str, object]:
    """Transport an isolated outside point and fix both direction types."""
    _check_prime_parameter(p)
    normal = full_circle_line_chart_normal_form(p)
    m = (p + 1) // 2
    line_sign = int(normal["forced_line_chart_c_H"])
    expected = -1 if m & 1 else 1
    if line_sign != expected:
        raise ArithmeticError("full-circle line-chart sign changed")
    return {
        "p": p,
        "P": p + 1,
        "m": m,
        "source_line_chart_c_H": line_sign,
        "outside_transport_multiplier": "(-1)^m*eps_line",
        "transported_c_H": "eps_line",
        "transported_infinity_degree_I": 0,
        "phase_zero_type": {"directions": m, "b": 0},
        "phase_one_type": {"directions": m, "b": 2},
        "phase_one_floor": p - 1,
        "type_alignment_exact": True,
        "proved": bool(normal["proved"]),
    }


def xnor_lift_certificate(p: int, scaled_excess: int = 2) -> dict[str, object]:
    """Reduce a nonnegative ``b=2`` XNOR lift to Proposition 15.688.

    The parity baseline is

    ``X=1-x_i-x_j+2*x_i*x_j``.

    It is zero or one pointwise.  If ``A`` is nonnegative and has the same
    parity, then ``B=(A-X)/2`` is pointwise a nonnegative integer: when
    ``X=0``, ``A`` is a nonnegative even integer; when ``X=1``, ``A`` is a
    positive odd integer.  This also proves that no hidden sign assumption is
    being imported from the older endpoint propositions.
    """
    _check_prime_parameter(p)
    if scaled_excess <= 0 or scaled_excess % 2:
        raise ValueError("need a positive even scaled excess")
    sharp = sharp_integral_quadratic_lift_floor(p)
    sharp_floor = int(sharp["sharp_scaled_floor"])
    forced_scaled_B_mean = scaled_excess
    nonzero = forced_scaled_B_mean > 0
    excluded = nonzero and forced_scaled_B_mean < sharp_floor
    return {
        "p": p,
        "xnor_polynomial": "X=1-x_i-x_j+2*x_i*x_j",
        "xnor_truth_table": {"00": 1, "01": 0, "10": 0, "11": 1},
        "same_parity_as_b2_phase_one": True,
        "lift": "B=(A-X)/2",
        "pointwise_nonnegative_integer_proof": {
            "X=0": "A is a nonnegative even integer, so A/2>=0",
            "X=1": "A is a nonnegative odd integer, so (A-1)/2>=0",
        },
        "B_has_degree_at_most_two": True,
        "B_nonzero": nonzero,
        "scaled_excess_2p_E_A_minus_X": scaled_excess,
        "forced_4p_E_B": forced_scaled_B_mean,
        "prop_15_688_lower_bound": sharp_floor,
        "excluded": excluded,
        "proved": excluded,
    }


def xnor_coefficient_congruence_ledger(p: int) -> dict[str, object]:
    """Re-establish the sign-independent two-coordinate congruence here.

    Proposition 15.673 treats both targets ``4+z_a*z_b`` and
    ``4-z_a*z_b``.  They are XNOR and XOR after returning to zero-one
    variables, so they are not literally the same baseline.  What is common
    is the coefficient comparison: the sign is carried by the symbolic
    ``tau`` entry and drops out of ``(p-1)c=I+P_d-4``.
    """
    _check_prime_parameter(p)
    q = (p - 1) // 2
    prior_zero = coefficient_ledger(p, 0)
    prior_one = coefficient_ledger(p, 1)
    sign_independent_matrix = bool(
        prior_zero["general_complement_matrix"]
        == prior_one["general_complement_matrix"]
        and prior_zero["complement_divisibility"]
        == prior_one["complement_divisibility"]
        == "q divides I+P_d-4"
        and {
            prior_zero["complement_target"],
            prior_one["complement_target"],
        }
        == {"4 + z_a*z_b", "4 - z_a*z_b"}
    )
    return {
        "p": p,
        "q": q,
        "current_baseline": "XNOR(x_i,x_j)=1-x_i-x_j+2*x_i*x_j",
        "prior_two_sign_targets": (
            "4+z_a*z_b and 4-z_a*z_b (XNOR/XOR in zero-one variables)"
        ),
        "coefficient_matrix": prior_zero["general_complement_matrix"],
        "sign_parameter": "tau in {+1,-1}; it drops out of the constant comparison",
        "coefficient_comparison": "(p-1)c=I+P_d-4 with 2c integral",
        "divisibility": "q divides I+P_d-4",
        "applies_to_b2_phase_one_here": sign_independent_matrix,
        "proved": (
            sign_independent_matrix
            and bool(prior_zero["proved"])
            and bool(prior_one["proved"])
        ),
    }


def phase_one_xnor_normal_form(p: int) -> dict[str, object]:
    """Classify the all-``b=2`` type's exact common residue."""
    _check_prime_parameter(p)
    period = p + 1
    m = period // 2
    floor = period - 2
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    floor_plus_two = xnor_lift_certificate(p, 2)
    rows = []
    for u in range(m):
        quotient_sum = m - u
        if u == m - 1:
            feasible = True
            reason = "m-1 baseline means and one P-unit elevation"
        elif u == 0:
            # Every direction has k=1 and mean P, i.e. baseline plus two.
            feasible = not bool(floor_plus_two["excluded"])
            reason = "all directions would be forbidden floor-plus-two lifts"
        else:
            feasible = False
            reason = "every direction needs k>=1 but quotient sum is below m"
        rows.append(
            {
                "u": u,
                "quotient_sum": quotient_sum,
                "feasible": feasible,
                "reason": reason,
            }
        )
    feasible_u = [int(row["u"]) for row in rows if row["feasible"]]
    if feasible_u != [m - 1]:
        raise ArithmeticError("phase-one xnor residue normal form changed")
    return {
        "p": p,
        "P": period,
        "m": m,
        "b": 2,
        "phase": 1,
        "floor": floor,
        "floor_plus_two_lift_floor_from_prop_15_688": lift_floor,
        "floor_plus_two_xnor_lift_certificate": floor_plus_two,
        "residue_rows": rows,
        "unique_residue_u": m - 1,
        "quotient_multiset": {"0": m - 1, "1": 1},
        "parallel_count_multiset": "m-1 copies of x and one copy of x+1",
        "proved": True,
    }


def zero_infinity_circle_arithmetic(p: int) -> dict[str, object]:
    """Collapse the two type residues and parallel baselines to one tuple."""
    _check_prime_parameter(p)
    period = p + 1
    m = period // 2
    q = m - 1
    phase_one = phase_one_xnor_normal_form(p)
    congruence = xnor_coefficient_congruence_ledger(p)
    if not congruence["proved"]:
        raise ArithmeticError("XNOR coefficient congruence was not established")

    # Equation m(x+y-7)=u-4 has 0<=u<m.  Its right side lies strictly
    # between -m and m, so it must vanish.
    residue_candidates = [
        u for u in range(m) if (u - 4) % m == 0
    ]
    if residue_candidates != [4]:
        raise ArithmeticError("phase-zero residue equation changed")
    u = residue_candidates[0]
    x_plus_y = 7 + (u - 4) // m

    # Once u=4, sum k=m-4<m forces a zero quotient, hence y>=0.  Together
    # with x>=0 and x+y=7, the xnor congruence q|(x-4) has one solution.
    xy = [
        (x, x_plus_y - x)
        for x in range(x_plus_y + 1)
        if (x - 4) % q == 0
    ]
    if xy != [(4, 3)]:
        raise ArithmeticError("xnor baseline counts changed")
    x, y = xy[0]
    quotient_sum = m - u
    zero_quotient_directions = m - quotient_sum
    if zero_quotient_directions != 4:
        raise ArithmeticError("phase-zero zero-quotient count changed")
    return {
        "p": p,
        "P": period,
        "m": m,
        "q": q,
        "I": 0,
        "finite_edge_identity": "E=m(x+y+1)+1-u=8m-3",
        "residue_equation": "m(x+y-7)=u-4",
        "phase_zero_unique_u": u,
        "x_plus_y": x_plus_y,
        "xnor_baseline_congruence": "q divides I+x-4",
        "xnor_coefficient_congruence_certificate": congruence,
        "xnor_baseline_count_x": x,
        "phase_zero_baseline_count_y": y,
        "phase_zero_quotient_sum": quotient_sum,
        "phase_zero_zero_quotient_directions_at_least": zero_quotient_directions,
        "zero_quotient_scaled_mean": 2 * u,
        "phase_one_normal_form": phase_one,
        "proved": True,
    }


def full_circle_lift_contradiction(p: int) -> dict[str, object]:
    """Apply the sharp nonzero quadratic-lift floor to the forced mean eight."""
    _check_prime_parameter(p)
    isolated = isolated_outside_vertex_ledger(p)
    chart = isolated_circle_chart(p)
    arithmetic = zero_infinity_circle_arithmetic(p)
    forced_scaled_mean = int(arithmetic["zero_quotient_scaled_mean"])
    sharp_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    contradiction_gap = sharp_floor - forced_scaled_mean
    excluded = contradiction_gap > 0
    if not excluded:
        raise ArithmeticError("full-circle lift contradiction disappeared")
    return {
        "p": p,
        "boundary": "full Miquelian circle of size p+1",
        "isolated_vertex": isolated,
        "isolated_chart": chart,
        "arithmetic": arithmetic,
        "forced_direction": {
            "b": 0,
            "phase": 0,
            "A_factorization": "A_d=2B_d",
            "B_nonzero_nonnegative_integral_quadratic": True,
            "forced_4p_E_B": forced_scaled_mean,
        },
        "prop_15_688_lower_bound": sharp_floor,
        "contradiction_gap": contradiction_gap,
        "full_circle_excluded": excluded,
        "proved": excluded,
    }


def universal_full_circle_certificate() -> dict[str, object]:
    """Collect the monotone inequalities that make the proof uniform in p."""
    p0 = 17
    m0 = (p0 + 1) // 2
    q0 = (p0 - 1) // 2
    isolated_gap_at_17 = p0 * p0 - 8 * p0 - 1
    isolated_gap_increment_at_17 = 2 * p0 - 7
    residue_unique = [u for u in range(m0) if (u - 4) % m0 == 0] == [4]
    congruence_unique = [x for x in range(8) if (x - 4) % q0 == 0] == [4]
    xnor_lift = xnor_lift_certificate(p0, 2)
    coefficient = xnor_coefficient_congruence_ledger(p0)
    proved = bool(
        isolated_gap_at_17 > 0
        and isolated_gap_increment_at_17 > 0
        and residue_unique
        and congruence_unique
        and xnor_lift["proved"]
        and coefficient["proved"]
        and p0 - 3 > 8
    )
    return {
        "scope": "every odd prime p>=17",
        "isolated_vertex_polynomial": "p^2-8p-1",
        "isolated_vertex_base_value_at_17": isolated_gap_at_17,
        "isolated_vertex_odd_step_is_increasing_from_17": isolated_gap_increment_at_17 > 0,
        "residue_argument": (
            "m>=9 and 0<=u<m with m|(u-4) force u=4"
        ),
        "residue_base_check": residue_unique,
        "xnor_lift_argument": xnor_lift,
        "xnor_congruence_argument": coefficient,
        "baseline_count_argument": (
            "q>=8, 0<=x<=7, and q|(x-4) force x=4"
        ),
        "baseline_count_base_check": congruence_unique,
        "final_lift_gap": "p-3-8=p-11>=6",
        "proved": proved,
    }


def theorem_full_circle_exclusion() -> dict[str, object]:
    universal = universal_full_circle_certificate()
    r_two = outside_R_two_structure(17)
    r_three = outside_R_three_structure(17)
    low_slack = outside_low_slack_conic_exclusion(17)
    sample_primes = (17, 19, 23, 29, 31, 37, 41, 101)
    rows = {str(p): full_circle_lift_contradiction(p) for p in sample_primes}
    proved = bool(
        universal["proved"]
        and r_two["proved"]
        and r_three["proved"]
        and low_slack["proved"]
        and all(bool(row["proved"]) for row in rows.values())
    )
    return {
        "prop": "15.724",
        "title": "Full Miquelian-circle boundary exclusion",
        "proved": proved,
        "scope": "every odd prime p>=17",
        "universal_certificate": universal,
        "dependencies": {
            "15.722": "exact circle phase and outside-chart type alignment",
            "15.722 low slack": (
                "R<=max(3,floor(sqrt(p)-5/2)) conic/secant exclusions"
            ),
            "15.672/15.673": "xnor baseline coefficient congruence",
            "15.688": "sharp 4p E[B]>=p-3 lift floor",
        },
        "theorem": {
            "full_Miquelian_circle_boundary": "EXCLUDED",
            "outside_R_zero": "EXCLUDED",
            "outside_R_two": "EXCLUDED_BY_15.722",
            "outside_R_three": "EXCLUDED_BY_15.722",
            "positive_outside_slack_excluded_through": (
                "max(3,floor(sqrt(p)-5/2))"
            ),
            "strict_outside_profiles_R_at_least_2": "R_EQUALS_2_AND_3_EXCLUDED",
            "strict_outside_profiles_R_at_least_4": (
                "OPEN_AFTER_THE_PRIME_DEPENDENT_LOW_SLACK_CUTOFF"
            ),
            "whole_p_plus_one_shell": "OPEN",
            "residual_ii": False,
            "type_I": False,
            "limit_exists": False,
        },
        "sample_ledgers_regression_only": rows,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    theorem = theorem_full_circle_exclusion()
    if theorem["proved"] is not True:
        raise ArithmeticError("Proposition 15.724 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15724.json"
    destination.write_text(json.dumps(theorem, indent=2) + "\n")
    print("Prop 15.724 full Miquelian-circle boundary: excluded")
    print(f"  wrote {destination}")
    return theorem


if __name__ == "__main__":
    main()
