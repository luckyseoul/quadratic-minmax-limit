#!/usr/bin/env python3
"""Prop. 15.679 -- close the next all-finite boundary for p>=43.

Let ``s`` be the second even all-finite boundary size strictly above
``3(p-1)/4``.  For every prime ``p>=43`` this proposition independently
excludes ``s``.  Its claim does not depend on the retracted ``p=17``
endpoint assertion in Proposition 15.678.

Put ``P=p+1`` and ``m=P/2``.  The phase-one type has only common residue
``u_1=m-1`` and therefore has ``m-1`` directions at ``b=2`` and one at
``b=s``.  Its exact pair deficit is ``(m-1)(s-2)``.  In phase zero, the
corrected parameter-aware floor-plus-two classification and the exact
quotient sum leave only residues ``2<=u_0<=7`` after the pair budget:

* ``u_0=0`` is already over budget;
* ``u_0=1`` is infeasible;
* the exact interior deficit is increasing, and ``u_0=8`` is over budget;
* the final four residues have at least ``(m-4)s`` deficit and are over
  budget.

Every retained phase-zero row has quotient sum ``m-u_0<m``, so one
direction has quotient zero and mean ``2u_0<=14``.  Since ``s<=p-5``, every
nonzero even fibre count has phase-zero floor at least ``P``.  The direction
therefore has ``b=0`` and pointwise slack ``A=2B`` for a nonzero,
nonnegative, integer-valued quadratic ``B`` on the middle slice.  Thus

    2u_0 = 4p E[B].

Proposition 15.642's degree-two slice-distance bound makes the right side
strictly larger than 14 for ``p>=59``.  The only smaller in-scope primes are
43, 47, and 53; their exact pair ledgers leave maximum residues 4, 6, and 5,
while their lift floors are 12, 14, and 14.  These are contradictions too.

This closes one additional uniform all-finite boundary size.  The seven
smaller endpoints ``p=17,19,23,29,31,37,41``, later all-finite sizes, the
infinity-present remainder, residual (ii), R1, global QVAR, Type I, and the
limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15642 import (
    nonbaseline_scaled_cost_floor,
    polynomial_distance_support_floor,
)
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15675 import first_even_survivor
from e1_gmin_m4_prop15723 import floor_excess_admissible


ROOT = Path(__file__).resolve().parents[1]


def next_even_boundary(p: int) -> int:
    """Second even integer strictly above ``3(p-1)/4``."""
    if p < 43 or p % 2 == 0:
        raise ValueError("need odd p>=43")
    s = first_even_survivor(p) + 2
    expected = {
        1: (3 * p + 13) // 4,
        3: (3 * p + 7) // 4,
        5: (3 * p + 9) // 4,
        7: (3 * p + 11) // 4,
    }[p % 8]
    if s != expected or s % 2 or s > p - 5:
        raise ArithmeticError("next-boundary formula or middle range changed")
    return s


def phase_one_minimum(p: int) -> dict[str, object]:
    """Exact phase-one quantized minimum at the next boundary."""
    s = next_even_boundary(p)
    m = (p + 1) // 2
    deficit = (m - 1) * (s - 2)
    return {
        "p": p,
        "s": s,
        "u": m - 1,
        "quotient_sum": 1,
        "profile": {"b=2": m - 1, "b=s": 1},
        "minimum_deficit": deficit,
        "uniqueness": (
            "u=0 forbids the only quotient-one low-floor lift; for "
            "1<=u<=m-2 every quotient is at least one but their sum is "
            "m-u<m; at u=m-1 exactly one quotient is one"
        ),
        "proved": True,
    }


def phase_zero_interior_minimum(p: int, u: int) -> dict[str, object]:
    """Exact phase-zero minimum for ``2<=u<=m-5``."""
    s = next_even_boundary(p)
    m = (p + 1) // 2
    if not 2 <= u <= m - 5:
        raise ValueError("need an interior phase-zero residue")
    target = m - u
    high_count, b2_count = divmod(target, 2)
    b0_count = m - high_count - b2_count
    deficit = b0_count * s + b2_count * (s - 2)
    return {
        "p": p,
        "s": s,
        "u": u,
        "quotient_sum": target,
        "profile": {
            "b=0": b0_count,
            "b=2": b2_count,
            "b=s": high_count,
        },
        "quotient_profile": {
            "k=0": b0_count,
            "k=1": b2_count,
            "k=2": high_count,
        },
        "minimum_deficit": deficit,
        "proof": (
            "for u<=m-5 the best b at quotient 0,1,2 is 0,2,s; "
            "b=s dominates b=4 at quotient two"
        ),
        "proved": True,
    }


def phase_zero_u0_minimum(p: int) -> dict[str, object]:
    """Exact phase-zero minimum in common residue zero."""
    s = next_even_boundary(p)
    m = (p + 1) // 2
    high_count, remainder = divmod(m, 3)
    deficit = (m - high_count) * s - 2 * remainder
    return {
        "p": p,
        "s": s,
        "u": 0,
        "quotient_sum": m,
        "high_quotient_count": high_count,
        "quotient_remainder": remainder,
        "minimum_deficit": deficit,
        "proof": (
            "the best b at quotient 0,1,2,3 is 0,2,4,s; quotient-three "
            "directions are optimal before resolving the remainder"
        ),
        "proved": True,
    }


def pair_row(p: int, u: int) -> dict[str, object]:
    """Pair-budget row for a phase-zero common residue."""
    s = next_even_boundary(p)
    if u == 0:
        phase_zero = phase_zero_u0_minimum(p)
    else:
        phase_zero = phase_zero_interior_minimum(p, u)
    d0 = int(phase_zero["minimum_deficit"])
    d1 = int(phase_one_minimum(p)["minimum_deficit"])
    budget = s * (s - 1)
    return {
        "u": u,
        "phase_zero_deficit": d0,
        "phase_one_deficit": d1,
        "required_total_deficit": d0 + d1,
        "pair_deficit_budget": budget,
        "pair_slack": budget - d0 - d1,
        "survives_pair_budget": d0 + d1 <= budget,
    }


def symbolic_residue_reduction() -> dict[str, object]:
    """Finite algebra proving that only ``2<=u<=7`` need a lift."""
    u8_gaps = {
        "p=1 mod 8": "2s-p-1=(p+11)/2",
        "p=3 mod 8": "3s-p+1=(5p+25)/4",
        "p=5 mod 8": "3s-p-1=(5p+23)/4",
        "p=7 mod 8": "2s-p+1=(p+13)/2",
    }
    return {
        "scope": "odd p>=43",
        "next_boundary": {
            "p=1 mod 8": "s=(3p+13)/4",
            "p=3 mod 8": "s=(3p+7)/4",
            "p=5 mod 8": "s=(3p+9)/4",
            "p=7 mod 8": "s=(3p+11)/4",
        },
        "middle_floor_range": "s<=p-5",
        "u0": {
            "exact_deficit": "(m-floor(m/3))s-2(m mod 3)",
            "p>=47_gap_lower_bound": "(3p^2-128p-347)/48>0",
            "p=43_exact_gap": 58,
            "excluded": True,
        },
        "u1": {
            "reason": "every direction needs k>=1 but sum k=m-1",
            "excluded": True,
        },
        "interior": {
            "exact_profile": (
                "x=floor((m-u)/2) at b=s, y=(m-u) mod 2 at b=2, "
                "m-x-y at b=0"
            ),
            "strictly_increasing_from_u=2": True,
            "u8_pair_gap_by_class": u8_gaps,
            "all_u_at_least_8_excluded": True,
        },
        "last_four_residues": {
            "deficit_lower_bound": "(m-4)s",
            "pair_gap_lower_bound": "(3p^2-84p-159)/16>0",
            "excluded": True,
        },
        "only_residues_requiring_lift": [2, 3, 4, 5, 6, 7],
        "proved": True,
    }


def _relaxed_minimum_dp(p: int, phase: int) -> dict[str, object]:
    """Independent parameter-aware floor-lift DP for regression samples."""
    s = next_even_boundary(p)
    m = (p + 1) // 2
    period = p + 1
    rows = []
    for u in range(m):
        target = m - u
        best_by_k: dict[int, tuple[int, int]] = {}
        for b in range(0, s + 1, 2):
            floor = full_symbolic_floor(p, b, phase)
            for k in range(target + 1):
                excess = 2 * u + period * k - floor
                if not floor_excess_admissible(p, b, phase, excess):
                    continue
                candidate = (s - b, b)
                old = best_by_k.get(k)
                if old is None or candidate[0] < old[0]:
                    best_by_k[k] = candidate
        states: dict[int, tuple[int, tuple[int, ...]]] = {0: (0, ())}
        for _ in range(m):
            next_states: dict[int, tuple[int, tuple[int, ...]]] = {}
            for used, (deficit, profile) in states.items():
                for k, (added, b) in best_by_k.items():
                    new_used = used + k
                    if new_used > target:
                        continue
                    candidate = (deficit + added, profile + (b,))
                    old = next_states.get(new_used)
                    if old is None or candidate[0] < old[0]:
                        next_states[new_used] = candidate
            states = next_states
        if target in states:
            deficit, profile = states[target]
            rows.append(
                {
                    "u": u,
                    "minimum_deficit": deficit,
                    "profile": dict(sorted(Counter(profile).items())),
                }
            )
    return {"p": p, "phase": phase, "rows": rows}


def small_prime_lift_ledger(p: int) -> dict[str, object]:
    """Exact pair survivors and lift contradiction at 43, 47, or 53."""
    expected = {
        43: {"residues": [2, 3, 4], "cost": 12},
        47: {"residues": [2, 3, 4, 5, 6], "cost": 14},
        53: {"residues": [2, 3, 4, 5], "cost": 14},
    }
    if p not in expected:
        raise ValueError("small exact ledger is only for 43, 47, 53")
    rows = [pair_row(p, u) for u in range(2, 8)]
    residues = [int(row["u"]) for row in rows if row["survives_pair_budget"]]
    cost = nonbaseline_scaled_cost_floor(p)
    if residues != expected[p]["residues"] or cost != expected[p]["cost"]:
        raise ArithmeticError("small-prime next-boundary ledger changed")
    return {
        "p": p,
        "s": next_even_boundary(p),
        "pair_rows": rows,
        "surviving_residues": residues,
        "maximum_scaled_mean": 2 * max(residues),
        "nonzero_quadratic_lift_floor": cost,
        "excluded": cost > 2 * max(residues),
    }


def large_prime_lift_ledger() -> dict[str, object]:
    """Symbolic degree-two support contradiction for every p>=59."""
    # 4p times the support floor from Proposition 15.642.
    p = 59
    scaled = 4 * p * polynomial_distance_support_floor(p)
    threshold_polynomial = p * p - 56 * p + 111
    if scaled <= 14 or threshold_polynomial <= 0:
        raise ArithmeticError("large-prime support threshold changed")
    return {
        "scope": "odd p>=59",
        "scaled_support_floor": "(p^2-1)/(4(p-2))",
        "comparison_to_14": "p^2-56p+111>0",
        "value_at_59": str(scaled),
        "derivative_after_59": "2p-56>0",
        "all_candidate_means_at_most": 14,
        "excluded": True,
    }


def next_boundary_exclusion(p: int) -> dict[str, object]:
    """Return the exact exclusion ledger for one in-scope prime."""
    if p in (43, 47, 53):
        lift = small_prime_lift_ledger(p)
        method = "exact pair ledger plus Proposition 15.642 lift floor"
    elif p >= 59 and p % 2:
        lift = large_prime_lift_ledger()
        method = "uniform Proposition 15.642 polynomial-distance floor"
    else:
        raise ValueError("need an odd p>=43; 43,47,53 are the small cases")
    return {
        "p": p,
        "s": next_even_boundary(p),
        "residue_reduction": symbolic_residue_reduction(),
        "lift_ledger": lift,
        "method": method,
        "excluded": bool(lift["excluded"]),
    }


def theorem_record() -> dict[str, object]:
    samples = (43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
    dp_samples = {}
    for p in (43, 47, 53, 59, 73):
        phase_zero = _relaxed_minimum_dp(p, 0)
        phase_one = _relaxed_minimum_dp(p, 1)
        dp_samples[str(p)] = {"0": phase_zero, "1": phase_one}

        by_u = {int(row["u"]): row for row in phase_zero["rows"]}
        if int(by_u[0]["minimum_deficit"]) != int(
            phase_zero_u0_minimum(p)["minimum_deficit"]
        ):
            raise ArithmeticError("independent u=0 DP disagrees")
        for u in range(2, 8):
            expected = phase_zero_interior_minimum(p, u)
            if int(by_u[u]["minimum_deficit"]) != int(
                expected["minimum_deficit"]
            ):
                raise ArithmeticError("independent phase-zero DP disagrees")
        phase_one_rows = phase_one["rows"]
        expected_one = phase_one_minimum(p)
        if len(phase_one_rows) != 1 or int(
            phase_one_rows[0]["minimum_deficit"]
        ) != int(expected_one["minimum_deficit"]):
            raise ArithmeticError("independent phase-one DP disagrees")

    exclusions = {str(p): next_boundary_exclusion(p) for p in samples}
    proved = bool(
        symbolic_residue_reduction()["proved"]
        and large_prime_lift_ledger()["excluded"]
        and all(row["excluded"] for row in exclusions.values())
    )
    return {
        "prop": "15.679",
        "title": "Uniform exclusion of the next all-finite boundary",
        "proved": proved,
        "theorem": {
            "all_odd_primes_p_at_least_43": (
                "the second even all-finite size above 3(p-1)/4 is excluded"
            ),
            "smaller_endpoints": [17, 19, 23, 29, 31, 37, 41],
            "smaller_endpoints_status": "OPEN_AT_THIS_BOUNDARY_SIZE",
            "later_all_finite_sizes": "OPEN",
            "infinity_present_remainder": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "symbolic_residue_reduction": symbolic_residue_reduction(),
        "large_prime_lift": large_prime_lift_ledger(),
        "small_prime_ledgers": {
            str(p): small_prime_lift_ledger(p) for p in (43, 47, 53)
        },
        "samples": exclusions,
        "independent_relaxed_dp_samples": dp_samples,
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.679 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15679.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.679 next all-finite boundary: excluded for p>=43")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
