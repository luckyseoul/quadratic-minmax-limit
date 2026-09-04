#!/usr/bin/env python3
r"""Proposition 15.766 -- a two-type collision closes an even bridge band.

Let ``C`` be the Paley conference matrix of order ``p^2+1``, let ``H``
have even size

    h=4p+2t,

and suppose that both signed Boolean eigenshells have the level-four floor

    T_H^eps(y)=eps*sum_(uv in H) C_uv y_u y_v >= 4.       (1)

This is exactly the branch left by Proposition 15.764 when an even
minimal-four-gap set has no Type-I shell bridge.  In the ranges below an
isolated vertex can be transported to infinity, so every edge of ``H`` is
finite.  Put ``M=p+1``, ``m=M/2`` and

    tau=sum_(uv in H) C_uv.

For a projective ``F_p``-direction ``d`` of quadratic type ``eps``, let
``P_d`` count the edges of ``H`` parallel to ``d`` and put, on the affine
middle-slice chart,

    A_d=(T_H^eps-4)/2,       a_d=2p E_d[A_d].              (2)

Exact middle-slice moments and the split of the ``p+1`` directions into
``m`` directions of each type give

    a_d=M P_d-eps*tau-4p,
    sum_(d:type eps) a_d=M t,
    a_d == 4-eps*tau (mod M).                              (3)

Every ``a_d`` is a nonnegative even integer.  Boundary parity, the central
Krawtchouk bound of Proposition 15.750, and the integral quadratic mass
bound of Proposition 15.681 sharpen this to

    a_d=0  or  a_d>=beta_p,                                (4)

where ``beta_p=(p+1)/2`` for ``p=3 mod 4`` and
``beta_p=(p-1)/2`` for ``p=1 mod 4``.  Indeed, a nonconstant parity costs
at least ``p-1``; a constantly odd chart costs at least ``2p``; and on a
constantly even chart ``A_d=2B_d``, so 15.681 gives
``a_d=4p E[B_d]>=beta_p``.

If ``a_d=0``, then ``A_d`` vanishes pointwise.  Johnson-slice constancy
makes all signed off-fibre block sums equal to one integer ``kappa`` and

    P_d-(p-1)kappa/2=4.

The edge-count inequality

    h-P_d >= binom(p,2)*abs(kappa)

forces ``kappa=0`` in the theorem range.  Hence ``P_d=4`` and (3) gives

    eps*tau=4.                                             (5)

The two direction types now collide.  If neither type has a zero, let
``r_+`` and ``r_-`` be the least residues in ``[0,M)`` from (3).  Since
the typewise average is ``2t<M``, (4) gives

    beta_p <= r_+,r_- <= 2t,
    r_+ + r_- == 8 (mod M).

As ``2 beta_p>8``, this requires ``r_++r_-=M+8``.  Below the displayed
threshold this is impossible, so one type has a zero.  Equation (5) makes
the other type's residue exactly eight.  For ``p>=19``, ``beta_p>8``, so
every positive entry in that type is at least ``M+8``, again above its
average.  For ``p=11,13,17`` the same last step holds when ``t<=3`` because
``2t<8``.

Consequently (1) is impossible in the following all-prime bands:

* ``p in {11,13,17}`` and ``1<=t<=3``;
* ``p>=19``, ``p=3 mod 4``, and ``1<=t<=(p+5)/4``;
* ``p>=19``, ``p=1 mod 4``, and ``1<=t<=(p+7)/4``.

In particular, ``h=4p+2`` is excluded for every prime ``p>=11``.  This is
a genuine enlargement of the Proposition 15.764 bridge, not a global
closure: larger even sizes and ``p=5,7`` remain open.

The same collision does not transfer to an odd residual-(ii) separator
``h=4p+2t+1``.  With baseline three the exact typewise average becomes
``M+2t``, rather than ``2t``.  ``odd_residual_nozero_profile`` gives, for
every ``t>=2``, an exact no-zero scalar/parallel-count profile satisfying
those budgets.  It is a method barrier only and is not asserted to be
realized by one Paley graph.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15681 import paired_cube_integral_quadratic_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15750 import parity_bias_theorem
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15766.json"


def _validate_prime(p: int) -> None:
    if type(p) is not int or p < 11 or not is_prime(p):
        raise ValueError("p must be a prime at least eleven")


def beta_p(p: int) -> int:
    """Least universal positive directional excess in (4)."""
    _validate_prime(p)
    beta = (p + 1) // 2 if p % 4 == 3 else (p - 1) // 2
    if beta % 2:
        raise ArithmeticError("the positive directional floor must be even")
    lift = paired_cube_integral_quadratic_floor(p)
    if Fraction(lift["universal_scaled_mass_floor"]) != beta:
        raise ArithmeticError("Proposition 15.681's mass floor changed")
    if not parity_bias_theorem()["proved"]:
        raise ArithmeticError("Proposition 15.750's parity bound changed")
    return beta


def closed_t_max(p: int) -> int:
    """Largest ``t`` proved impossible by the two-type collision."""
    _validate_prime(p)
    if p in (11, 13, 17):
        return 3
    if p % 4 == 3:
        return (p + 5) // 4
    return (p + 7) // 4


def directional_identity_ledger(p: int, t: int, eps: int) -> dict[str, object]:
    """Record the exact moment, type-sum, and congruence identities (3)."""
    _validate_prime(p)
    if type(t) is not int or t < 1 or eps not in (-1, 1):
        raise ValueError("need t>=1 and eps in {-1,+1}")
    M = p + 1
    m = M // 2
    h = 4 * p + 2 * t
    return {
        "p": p,
        "t": t,
        "H_size": h,
        "direction_type": eps,
        "directions_of_type": m,
        "A_d": "(eps*S_H-4)/2",
        "a_d": "2p*E_d[A_d]",
        "pointwise_formula": "a_d=(p+1)P_d-eps*tau-4p",
        "type_sum": M * t,
        "type_average": 2 * t,
        "common_residue_mod_p_plus_1": "4-eps*tau",
        "a_d_even": True,
        "proved": True,
    }


def zero_direction_rigidity(p: int, t: int) -> dict[str, object]:
    """Audit that a zero direction forces ``P_d=4`` and ``eps*tau=4``."""
    _validate_prime(p)
    if type(t) is not int or t < 1:
        raise ValueError("t must be positive")
    q = (p - 1) // 2
    h = 4 * p + 2 * t
    positive_kappa_edge_floor = 4 + q * (p + 1)
    negative_kappa_forces_P_at_most = 4 - q
    return {
        "p": p,
        "t": t,
        "H_size": h,
        "Johnson_constant_block_sum": "kappa",
        "constant_equation": "P_d-q*kappa=4",
        "nonparallel_edge_bound": "h-P_d>=binom(p,2)*abs(kappa)",
        "kappa_ge_one_forces_at_least_edges": positive_kappa_edge_floor,
        "strict_edge_gap": positive_kappa_edge_floor - h,
        "kappa_le_minus_one_forces_P_at_most": negative_kappa_forces_P_at_most,
        "kappa_zero": bool(
            positive_kappa_edge_floor > h
            and negative_kappa_forces_P_at_most < 0
        ),
        "P_d": 4,
        "zero_direction_consequence": "eps*tau=4",
        "proved": bool(
            positive_kappa_edge_floor > h
            and negative_kappa_forces_P_at_most < 0
        ),
    }


def even_band_row(p: int, t: int) -> dict[str, object]:
    """Audit every strict inequality used in the even-band contradiction."""
    _validate_prime(p)
    if type(t) is not int or not 1 <= t <= closed_t_max(p):
        raise ValueError("t lies outside the proved band")
    M = p + 1
    m = M // 2
    h = 4 * p + 2 * t
    beta = beta_p(p)
    isolation_margin = p * p + 1 - 2 * h
    rigidity = zero_direction_rigidity(p, t)

    nozero_both_required_sum = M + 8
    nozero_both_maximum_sum = 4 * t
    if beta > 8:
        opposite_least_positive = M + 8
    else:
        opposite_least_positive = 8

    proved = bool(
        isolation_margin > 0
        and bool(rigidity["proved"])
        and 2 * t < M
        and 2 * beta > 8
        and nozero_both_maximum_sum < nozero_both_required_sum
        and 2 * t < opposite_least_positive
    )
    if not proved:
        raise ArithmeticError("the even-band collision inequality failed")
    return {
        "p": p,
        "t": t,
        "H_size": h,
        "vertices": p * p + 1,
        "isolated_vertex_margin": isolation_margin,
        "M": M,
        "directions_per_type": m,
        "beta_p": beta,
        "type_average": 2 * t,
        "type_average_below_modulus": 2 * t < M,
        "if_both_types_no_zero_required_residue_sum": nozero_both_required_sum,
        "maximum_possible_residue_sum_from_average": nozero_both_maximum_sum,
        "therefore_some_type_has_zero": True,
        "zero_type_forces_opposite_residue": 8,
        "opposite_least_positive": opposite_least_positive,
        "opposite_type_average": 2 * t,
        "zero_direction_rigidity": rigidity,
        "contradiction": True,
        "proved": proved,
    }


def odd_residual_nozero_profile(p: int, t: int) -> dict[str, object]:
    """Exact scalar obstruction to reusing the collision for residual (ii).

    This realizes the moment, congruence, and parallel-count ledgers only.
    It does not assert that the rows are restrictions of one Paley graph.
    """
    _validate_prime(p)
    if type(t) is not int or t < 2:
        raise ValueError("the displayed no-zero profile needs t>=2")
    M = p + 1
    m = M // 2
    h = 4 * p + 2 * t + 1
    tau = 1

    plus_q = [t - 1] + [0] * (m - 1)
    minus_q = [t - 2] + [0] * (m - 1)
    plus_a = [M + 2 + M * q for q in plus_q]
    minus_a = [M + 4 + M * q for q in minus_q]
    plus_P = [4 + q for q in plus_q]
    minus_P = [4 + q for q in minus_q]

    target_type_sum = m * (M + 2 * t)
    target_plus_parallel = (h + tau) // 2
    target_minus_parallel = (h - tau) // 2
    beta = beta_p(p)
    proved = bool(
        sum(plus_a) == sum(minus_a) == target_type_sum
        and sum(plus_P) == target_plus_parallel
        and sum(minus_P) == target_minus_parallel
        and all(value > 0 and value >= beta for value in plus_a + minus_a)
        and all((value - 2) % M == 0 for value in plus_a)
        and all((value - 4) % M == 0 for value in minus_a)
    )
    if not proved:
        raise ArithmeticError("the odd residual no-zero profile changed")
    return {
        "p": p,
        "t": t,
        "H_size": h,
        "tau": tau,
        "baseline": 3,
        "type_average": M + 2 * t,
        "target_type_sum": target_type_sum,
        "plus_type": {
            "residue": 2,
            "a_values": plus_a,
            "parallel_counts": plus_P,
            "quotient_sum_above_M_plus_2": t - 1,
        },
        "minus_type": {
            "residue": 4,
            "a_values": minus_a,
            "parallel_counts": minus_P,
            "quotient_sum_above_M_plus_4": t - 2,
        },
        "all_directions_nonzero": True,
        "Paley_graph_realizability_claimed": False,
        "purpose": "the zero-direction collision alone cannot close residual (ii)",
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    sample_primes = (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    rows = {
        str(p): [even_band_row(p, t) for t in range(1, closed_t_max(p) + 1)]
        for p in sample_primes
    }
    odd_barrier = {
        str(p): odd_residual_nozero_profile(p, 2)
        for p in (11, 13, 17, 19, 29)
    }
    proved = bool(
        all(row["proved"] for prime_rows in rows.values() for row in prime_rows)
        and all(row["proved"] for row in odd_barrier.values())
    )
    return {
        "prop": "15.766",
        "title": "Two-direction-type collision closes an even minimal-gap-four band",
        "proved": {
            "directional_moment_and_type_sum_identities": proved,
            "positive_directional_excess_floor": proved,
            "zero_direction_rigidity": proved,
            "even_bridge_band": proved,
            "all_even_sizes": False,
            "p5_p7": False,
            "residual_ii": False,
            "minimal_gap4_shell_bridge_closed_general": False,
            "e1_closed_general": False,
        },
        "closed_ranges": {
            "p_11_13_17": "1<=t<=3",
            "p_ge_19_p_3_mod_4": "1<=t<=(p+5)/4",
            "p_ge_19_p_1_mod_4": "1<=t<=(p+7)/4",
            "H_size": "4p+2t",
        },
        "sample_rows": rows,
        "odd_residual_method_barrier": odd_barrier,
        "remaining": (
            "Even sizes above the displayed band and p=5,7 remain open; "
            "the odd residual-(ii) type budget has average p+1+2t and admits "
            "the displayed no-zero scalar profiles."
        ),
    }


def main() -> dict[str, object]:
    theorem = theorem_record()
    write_json_atomic(EV, theorem)
    print("Prop. 15.766 even bridge band: proved")
    print("  residual (ii), full bridge, and E1: OPEN")
    print(f"  wrote {EV}")
    return theorem


if __name__ == "__main__":
    main()
