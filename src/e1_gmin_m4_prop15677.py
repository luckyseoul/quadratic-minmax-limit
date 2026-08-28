#!/usr/bin/env python3
"""Prop. 15.677 -- close the first finite survivor for every prime p>=19.

Proposition 15.675 excludes the first even all-finite boundary size above
``3(p-1)/4`` when ``p=3,5 (mod 8)``.  This proposition treats the two
remaining residue classes for ``p>=23``.  Together, they close the first
survivor for every prime ``p>=19``.

Write ``m=(p+1)/2`` and let ``u_0,u_1`` be the common half-residues of the
two exact directional-mean types.  The phase-one quotient sum forces
``u_1=m-1``; all but one phase-one direction are therefore exact ``b=2``
xnor baselines.  A complete quotient/deficit calculation leaves

* ``u_0=2`` in both remaining classes; and
* the additional ``u_0=3`` row only when ``p=1 (mod 8)``.

Let ``j-1`` be the phase-zero baseline parallel count and ``l`` the xnor
baseline parallel count.  The xnor coefficient congruence forces ``j=2``
in the first row and ``j=3`` in the second.  Exact inter-fibre l1 bounds
then eliminate every ``u_0=3`` row and, for ``p=7 (mod 8)``, the odd
``l=1,3,5`` rows.  The only pre-lift arithmetic regimes are

    l=2: I=2p,   E=2p+1,
    l=4: I=p-1,  E=3p+2.

They are not actual survivors.  Since ``sum k_d=m-u_0<m``, every retained
phase-zero row has a zero-quotient direction.  Its mean is ``2u_0`` and its
floor forces ``b=0``.  Hence its nonnegative slack is even pointwise,
``A_d=2B_d``, where ``B_d`` is a nonzero nonnegative integer-valued
quadratic on the middle slice.  Proposition 15.642 gives

    2u_0 = 4p E[B_d] >= nonbaseline_scaled_cost_floor(p).

The right side is greater than six for every ``p>=23``.  Thus ``u_0=2`` is
always impossible and the extra ``u_0=3`` row is impossible whenever it can
occur (the first in-scope prime is ``p=41``).

For completeness, in every one of the ``m-1`` xnor directions, if ``n_s`` counts
infinity neighbours in fibre ``s``, the coefficient capacity forces

* ``l=2``: ``(n_s)`` is all twos or one 1, one 3, and all other entries 2;
* ``l=4``: ``(n_s)`` is ``p-1`` ones and one zero, or one 2, ``p-3``
  ones, and two zeros.

Together with Proposition 15.675, the first all-finite survivor is excluded
for every odd prime ``p>=19``.  The exceptional smaller endpoint ``p=17`` is
not claimed here: its exact residue ledger has an additional ``u_0=0`` row
that the uniform argument does not remove.  The near-perfect profiles remain
useful as an independent pre-lift normal form, but no in-scope graph reaches
them.  The ``p=17`` endpoint, later all-finite sizes, residual (ii), R1,
QVAR, Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15642 import nonbaseline_scaled_cost_floor


ROOT = Path(__file__).resolve().parents[1]


def parameters(p: int) -> tuple[int, int, int]:
    if p < 23 or p % 8 not in (1, 7):
        raise ValueError("need p>=23 with p=1 or 7 mod 8")
    m = (p + 1) // 2
    q = (p - 1) // 2
    s = {
        1: (3 * p + 5) // 4,
        7: (3 * p + 3) // 4,
    }[p % 8]
    if s % 2 or not 4 * s > 3 * (p - 1):
        raise ArithmeticError("first-survivor size formula changed")
    return q, m, s


def zero_sum_additive_lower_bound(length: int, positive_mass: int) -> int:
    """Lower bound for sum_{i<j}|a_i+a_j| at fixed positive mass.

    This is Proposition 15.645's elementary sign-counting inequality, with
    the vector length made explicit.
    """
    if length < 3 or positive_mass < 0:
        raise ValueError("invalid additive-vector parameters")
    if positive_mass == 0:
        return 0
    if 2 * positive_mass <= length:
        return 2 * positive_mass * (length - positive_mass - 1)
    return (length - 2) * positive_mass


def phase_residue_reduction(p: int) -> dict[str, object]:
    """Complete residue reduction before edge-count arithmetic."""
    _q, m, s = parameters(p)
    if m & 1:
        d2 = ((m + 3) // 2) * s - 2
    else:
        d2 = ((m + 2) // 2) * s
    d1 = (m - 1) * (s - 2)
    pair_budget = s * (s - 1)
    slack = pair_budget - d2 - d1
    expected_slack = (p - 1) // 4 if p % 8 == 1 else (p - 7) // 4
    if slack != expected_slack:
        raise ArithmeticError("first-survivor pair slack changed")

    # For u>=2 and u<=m-5, quotient weights 0,1,2 occur at b=0,2,s.
    # The first two increments are therefore exact.  Once u>=m-4, at
    # most four directions can have positive quotient, giving the displayed
    # endpoint lower bound.  These two ranges overlap at the needed barrier.
    u3_increment = 2 if m & 1 else s - 2
    u4_increment = s
    endpoint_lower = (m - 4) * s
    endpoint_covers = endpoint_lower >= d2 + s

    # At u=0, at most floor(m/3) directions can use the largest b=s option.
    # This coarse lower bound is already strictly beyond the available slack.
    u0_lower = ((2 * m + 2) // 3) * s - 4
    u0_excluded = u0_lower > d2 + slack
    u1_infeasible = True  # every direction needs k>=1 but sum k=m-1
    u3_survives = u3_increment <= slack
    candidates = [2] + ([3] if u3_survives else [])
    proved = bool(
        u0_excluded
        and u1_infeasible
        and endpoint_covers
        and u4_increment > slack
        and (u3_survives == (p % 8 == 1))
    )
    return {
        "p": p,
        "p_mod_8": p % 8,
        "m": m,
        "s": s,
        "phase_one_only_residue": m - 1,
        "phase_one_baseline_directions": m - 1,
        "phase_one_special_directions": 1,
        "phase_zero_u2_deficit": d2,
        "phase_one_deficit": d1,
        "pair_budget": pair_budget,
        "available_slack": slack,
        "u0_coarse_lower_bound": u0_lower,
        "u0_excluded": u0_excluded,
        "u1_infeasible": u1_infeasible,
        "u3_increment": u3_increment,
        "u4_and_later_increment_lower_bound": u4_increment,
        "endpoint_range_lower_bound": endpoint_lower,
        "endpoint_range_covers_u4_barrier": endpoint_covers,
        "phase_zero_candidate_residues": candidates,
        "proved": proved,
    }


def _balanced_base_l1(p: int, scalar: int, total: int) -> dict[str, int]:
    """Minimum of sum_{s<t}|scalar-n_s-n_t| at fixed integer sum."""
    low, high_count = divmod(total, p)
    high = low + 1
    low_count = p - high_count
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


def arithmetic_rows(p: int, u0: int) -> dict[str, object]:
    """Enumerate all baseline counts after the xnor congruence."""
    q, m, s = parameters(p)
    residue = phase_residue_reduction(p)
    if u0 not in residue["phase_zero_candidate_residues"]:
        raise ValueError("u0 is not a surviving phase-zero residue")
    # I>=0 gives j+l<=7, while q>=8.  Substitution in the baseline
    # coefficient congruence gives q|(u0-j), so this short box has one row.
    j_candidates = [j0 for j0 in range(1, 8) if (u0 - j0) % q == 0]
    if j_candidates != [u0]:
        raise ArithmeticError("phase-zero baseline offset is no longer unique")
    j = j_candidates[0]
    rows = []
    for l in range(0, 9):
        n0 = m * j - u0
        n1 = m * l + 1
        finite = n0 + n1
        infinity = 4 * p + 1 - finite
        if infinity < 0 or infinity % 2:
            continue
        if infinity > s + 2 * finite:
            continue
        if (infinity + l - 4) % q:
            continue
        scalar = (infinity + l - 4) // q
        capacity = finite - l
        balanced = _balanced_base_l1(p, scalar, infinity)
        # The distinguished +1 entry can reduce the base norm by at most one.
        lower = max(0, balanced["base_l1_minimum"] - 1)
        rows.append(
            {
                "l": l,
                "phase_zero_baseline_offset_j": j,
                "phase_zero_finite_edges": n0,
                "phase_one_finite_edges": n1,
                "E": finite,
                "I": infinity,
                "xnor_scalar": scalar,
                "transverse_edge_capacity": capacity,
                "balanced_base": balanced,
                "xnor_l1_lower_bound": lower,
                "excluded_by_l1": lower > capacity,
            }
        )
    return {
        "p": p,
        "u0": u0,
        "q": q,
        "m": m,
        "j_box": "1<=j<=7",
        "j_candidates": j_candidates,
        "forced_j": j,
        "congruence": "q divides I+l-4, hence q divides u0-j",
        "rows": rows,
    }


def remaining_regimes(p: int) -> dict[str, object]:
    """Return the exact two pre-lift regimes left after l1 exclusion."""
    reduction = phase_residue_reduction(p)
    ledgers = {
        str(u0): arithmetic_rows(p, u0)
        for u0 in reduction["phase_zero_candidate_residues"]
    }
    survivors = [
        row
        for ledger in ledgers.values()
        for row in ledger["rows"]
        if not row["excluded_by_l1"]
    ]
    expected = [
        {"u0": 2, "l": 2, "I": 2 * p, "E": 2 * p + 1},
        {"u0": 2, "l": 4, "I": p - 1, "E": 3 * p + 2},
    ]
    observed = [
        {"u0": 2, "l": row["l"], "I": row["I"], "E": row["E"]}
        for row in survivors
    ]
    if observed != expected:
        raise ArithmeticError("remaining first-survivor normal form changed")
    return {
        "p": p,
        "residue_reduction": reduction,
        "arithmetic_ledgers": ledgers,
        "survivors": survivors,
        "normal_form": expected,
        "status": "PRE_LIFT_ONLY",
        "proved": True,
    }


def zero_quotient_lift_exclusion(p: int, u0: int) -> dict[str, object]:
    """Apply Proposition 15.642 to a forced phase-zero k=0 direction."""
    _q, m, _s = parameters(p)
    if not 1 <= u0 < m:
        raise ValueError("u0 must force quotient sum below the direction count")
    scaled_mean = 2 * u0
    cost = nonbaseline_scaled_cost_floor(p)
    next_even_b_floor = p + 1
    if not scaled_mean < next_even_b_floor:
        raise ArithmeticError("a zero quotient no longer forces b=0")
    return {
        "p": p,
        "u0": u0,
        "quotient_sum": m - u0,
        "direction_count": m,
        "zero_quotient_direction_forced": m - u0 < m,
        "zero_quotient_mean": scaled_mean,
        "zero_quotient_b": 0,
        "next_even_b_floor": next_even_b_floor,
        "mean_below_every_nonzero_even_b_floor": True,
        "pointwise_factorization": "A_d=2B_d",
        "B_d": "nonzero nonnegative integer-valued quadratic on J(p,m)",
        "prop15642_nonbaseline_scaled_cost": cost,
        "excluded": cost > scaled_mean,
    }


def first_survivor_exclusion(p: int) -> dict[str, object]:
    """Close the p=1,7 mod 8 first-survivor branch."""
    reduction = phase_residue_reduction(p)
    pre_lift = remaining_regimes(p)
    branches = []
    for u0 in reduction["phase_zero_candidate_residues"]:
        lift = zero_quotient_lift_exclusion(p, u0)
        arithmetic = arithmetic_rows(p, u0)
        all_coefficient_rows_excluded = all(
            row["excluded_by_l1"] for row in arithmetic["rows"]
        )
        excluded = bool(lift["excluded"] or all_coefficient_rows_excluded)
        branches.append(
            {
                "u0": u0,
                "lift": lift,
                "all_coefficient_rows_excluded": all_coefficient_rows_excluded,
                "method": (
                    "Prop. 15.642 nonzero-lift mass"
                    if lift["excluded"]
                    else "xnor inter-fibre l1"
                ),
                "excluded": excluded,
            }
        )
    return {
        "p": p,
        "s": reduction["s"],
        "p_mod_8": p % 8,
        "residue_reduction": reduction,
        "pre_lift_normal_form": pre_lift,
        "branches": branches,
        "excluded": all(row["excluded"] for row in branches),
    }


def fibre_profile_classification(p: int, l: int) -> dict[str, object]:
    """Classify infinity-neighbour counts in every xnor baseline direction."""
    parameters(p)
    if l == 2:
        infinity = 2 * p
        capacity = 2 * p - 1
        # Put a_s=n_s-2.  The target +1 changes l1 by at most one, so the
        # additive zero-sum matrix has norm at most 2p.  Positive mass >=2
        # already exceeds this for p>=7.
        first_forbidden = zero_sum_additive_lower_bound(p, 2)
        excluded = first_forbidden > capacity + 1
        profiles = [
            {"2": p},
            {"1": 1, "2": p - 2, "3": 1},
        ]
        proof = "zero-sum deviations have positive mass at most one"
    elif l == 4:
        infinity = p - 1
        capacity = 3 * p - 2
        # Put a_s=n_s-1 and append one ghost coordinate +1.  This is a
        # zero-sum vector of length p+1.  Its additive norm is the base
        # coefficient norm plus p-1, hence at most 4p-2.  Positive mass
        # >=3 is already impossible.
        first_forbidden = zero_sum_additive_lower_bound(p + 1, 3)
        excluded = first_forbidden > 4 * p - 2
        profiles = [
            {"0": 1, "1": p - 1},
            {"0": 2, "1": p - 3, "2": 1},
        ]
        proof = "ghost-augmented zero-sum deviations have positive mass at most two"
    else:
        raise ValueError("l must be 2 or 4")
    return {
        "p": p,
        "l": l,
        "I": infinity,
        "transverse_edge_capacity": capacity,
        "first_forbidden_additive_l1": first_forbidden,
        "forbidden_mass_exceeds_capacity": excluded,
        "allowed_fibre_count_histograms": profiles,
        "proof": proof,
        "proved": excluded,
    }


def theorem_record() -> dict[str, object]:
    sample_primes = (23, 31, 41, 47, 71, 73, 79, 89)
    samples = {str(p): first_survivor_exclusion(p) for p in sample_primes}
    profiles = {
        str(p): {
            str(l): fibre_profile_classification(p, l) for l in (2, 4)
        }
        for p in (23, 41, 73)
    }
    proved = bool(
        all(row["excluded"] for row in samples.values())
        and all(
            row["proved"]
            for by_l in profiles.values()
            for row in by_l.values()
        )
    )
    return {
        "prop": "15.677",
        "title": "Complete close of the first all-finite survivor",
        "proved": proved,
        "theorem": {
            "scope": "first even all-finite s>3(p-1)/4, prime p>=19",
            "p_mod_8_3_5": "CLOSED_BY_15.675",
            "p_mod_8_1_7_p_at_least_23": "CLOSED_HERE",
            "all_odd_primes_p_at_least_19": "FIRST_SURVIVOR_EXCLUDED",
            "p17_endpoint": "OPEN_ADDITIONAL_U0_ZERO_ROW",
            "pre_lift_l2": {"I": "2p", "E": "2p+1"},
            "pre_lift_l4": {"I": "p-1", "E": "3p+2"},
            "simultaneous_near_perfect_direction_geometry": (
                "VACUOUS_AFTER_LIFT_EXCLUSION"
            ),
            "later_all_finite_sizes": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "fibre_profile_classification": profiles,
        "samples": samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.677 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15677.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.677 first all-finite survivor: closed for every prime p>=19")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
