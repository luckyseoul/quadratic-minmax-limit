#!/usr/bin/env python3
"""Parallel-parity obstruction at the all-active Mobius support endpoint.

This module records a symbolic all-prime consequence of three exact inputs:

* the balanced branch-C parallel quotas;
* one localized Mobius half has odd parallel count only in its target and
  auxiliary directions; and
* cancellation changes every direction count by an even integer.

It does not enumerate primes, auxiliaries, supports, or target cells.  The
result is a necessary lower bound on the cancellation offset beyond the bare
support-size floor.  Passing the bound does not construct the remaining
symmetric Boolean fibre.
"""

from __future__ import annotations

import json

from e1_gmin_m4_prop15721 import is_prime


def _check_branch_c_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=3 mod 4 with p>=31")
    return (p - 3) // 4


def _balanced(total: int, count: int) -> tuple[int, ...]:
    quotient, remainder = divmod(total, count)
    return tuple(
        quotient + int(index < remainder) for index in range(count)
    )


def _closed_parity_weight(m: int, residue: int) -> int:
    """Return the closed value of ``wt(P mod 2 + 1_hard)``.

    Here ``residue=(t+1) mod (2m)`` and the theorem is used only for
    ``m>=16``.  The six cases are written symmetrically around the two
    transition intervals.
    """
    if not 0 <= residue < 2 * m:
        raise ValueError("residue must lie in [0,2m)")
    if 5 <= residue <= m:
        return m + 5
    if residue in (4, m + 1):
        return m + 3
    if residue in (3, m + 2):
        return m + 1
    if residue in (2, m + 3):
        return m - 1
    if residue in (1, m + 4):
        return m - 3
    return m - 5


def balanced_parallel_parity_profile(p: int, t: int) -> dict[str, object]:
    """Compute the exact all-active balanced branch-C parity profile.

    There are ``m=(p+1)/2`` hard rows and ``m`` opposite rows.  The hard
    compact counts are the balanced allocation of ``t+1`` and their physical
    parallel quotas are ``3+e_L``.  The opposite quotas are the balanced
    allocation of ``10r+6+t``.
    """
    r = _check_branch_c_prime(p)
    m = 2 * r + 2
    t_min = 2 * r * r - 4 * r - 2
    t_max = 4 * r * r - 2 * r - 5
    if not isinstance(t, int) or isinstance(t, bool) or not t_min <= t <= t_max:
        raise ValueError(f"need {t_min}<=t<={t_max}")

    hard_compact_total = t + 1
    opposite_quota_total = 10 * r + 6 + t
    hard_excesses = _balanced(hard_compact_total, m)
    hard_quotas = tuple(3 + value for value in hard_excesses)
    opposite_quotas = _balanced(opposite_quota_total, m)

    # Adding the hard-direction indicator cancels the odd constant 3.
    base_direction_parity = tuple(value & 1 for value in hard_excesses) + tuple(
        value & 1 for value in opposite_quotas
    )
    direct_weight = sum(base_direction_parity)
    residue = hard_compact_total % (2 * m)
    closed_weight = _closed_parity_weight(m, residue)

    quotient, alpha = divmod(hard_compact_total, m)
    if alpha >= 5:
        derived_weight = m + 5 if quotient % 2 == 0 else m - 5
    elif quotient % 2 == 0:
        derived_weight = m + 2 * alpha - 5
    else:
        derived_weight = m + 5 - 2 * alpha

    # The parity of m auxiliary directions has even weight at most m.  Since
    # direct_weight is odd, the sharp aggregate lower bound is never zero.
    minimum_fixed_edges = max(1, direct_weight - m)
    minimum_extra_cancellations = (minimum_fixed_edges - 1) // 2

    graph_edge_count = 4 * p + 2 * t + 1
    raw_half_occurrences = m * (p - 1)
    size_floor_cancellations = t_max - t + 1
    if raw_half_occurrences - graph_edge_count != 2 * size_floor_cancellations - 1:
        raise ArithmeticError("the all-active support-size ledger changed")
    if opposite_quota_total != hard_compact_total + 5 * m - 5:
        raise ArithmeticError("the two balanced quota totals lost their shift")
    if not direct_weight == derived_weight == closed_weight:
        raise ArithmeticError("the closed parallel-parity weight changed")

    return {
        "p": p,
        "r": r,
        "m_hard_rows": m,
        "direction_count": 2 * m,
        "t": t,
        "t_interval": [t_min, t_max],
        "hard_compact_total_E": hard_compact_total,
        "opposite_quota_total": opposite_quota_total,
        "opposite_minus_hard_total": 5 * m - 5,
        "hard_compact_counts": list(hard_excesses),
        "hard_parallel_quotas": list(hard_quotas),
        "opposite_parallel_quotas": list(opposite_quotas),
        "base_direction_parity": list(base_direction_parity),
        "residue_s_equals_E_mod_p_plus_1": residue,
        "base_direction_parity_weight_w0": direct_weight,
        "minimum_forced_fixed_edge_weight": minimum_fixed_edges,
        "minimum_extra_cancellations_beyond_size_floor": (
            minimum_extra_cancellations
        ),
        "size_floor_cancellations": size_floor_cancellations,
        "strengthened_cancellation_lower_bound": (
            size_floor_cancellations + minimum_extra_cancellations
        ),
        "graph_edge_count": graph_edge_count,
        "raw_Mobius_half_occurrences": raw_half_occurrences,
        "proved": True,
    }


def cancellation_offset_consequence(p: int, t: int, j: int) -> dict[str, object]:
    """Apply the parity theorem at cancellation offset ``j``.

    If ``kappa=t_max-t+1+j``, then the physical capacity left after the
    actual Mobius support is ``1+2j``.  The exact Hamming equation makes this
    an upper bound for the forced fixed-edge weight.  The parity theorem
    supplies the lower bound returned by :func:`balanced_parallel_parity_profile`.
    """
    profile = balanced_parallel_parity_profile(p, t)
    if not isinstance(j, int) or isinstance(j, bool) or j < 0:
        raise ValueError("j must be a nonnegative integer")

    minimum_j = int(profile["minimum_extra_cancellations_beyond_size_floor"])
    minimum_fixed = int(profile["minimum_forced_fixed_edge_weight"])
    remaining_capacity = 1 + 2 * j
    excluded = j < minimum_j
    fixed_only_forced = bool(j == minimum_j and minimum_fixed == remaining_capacity)

    if j == 0:
        endpoint = (
            "excluded by parallel parity"
            if excluded
            else "not excluded by parallel parity; T_U=A*a(T_U) remains target-dependent"
        )
    elif j == 1:
        if excluded:
            endpoint = "both three-fixed and one-fixed/one-column branches excluded"
        elif fixed_only_forced:
            endpoint = "three fixed edges and zero divided columns forced"
        else:
            endpoint = (
                "parallel parity leaves the three-fixed branch and the exact "
                "one-fixed/one-column target test"
            )
    else:
        endpoint = "necessary cancellation-offset test only"

    return {
        **profile,
        "cancellation_offset_j": j,
        "remaining_physical_edge_capacity": remaining_capacity,
        "maximum_forced_fixed_edge_weight_from_Hamming": remaining_capacity,
        "excluded_by_parallel_parity": excluded,
        "all_remaining_capacity_forced_to_fixed_edges_at_bound": fixed_only_forced,
        "conditional_unused_double_orbits_at_bound": 0 if fixed_only_forced else None,
        "endpoint_consequence": endpoint,
        "symmetric_Boolean_completion_constructed": False,
        "residual_ii_closed": False,
    }


def theorem_record(p: int = 31, t: int | None = None) -> dict[str, object]:
    r = _check_branch_c_prime(p)
    if t is None:
        t = 2 * r * r - 4 * r - 2
    sample = balanced_parallel_parity_profile(p, t)
    m = int(sample["m_hard_rows"])
    return {
        "title": "All-active Mobius parallel-parity endpoint obstruction",
        "scope": "balanced branch C, p=4r+3>=31, every hard center nonzero",
        "exact_input": {
            "hard_quotas": "P_L=3+e_L, balanced sum e_L=t+1",
            "opposite_quotas": "Q_L balanced with sum 10r+6+t",
            "one_half_mod_2": "parallel word e_L+e_M",
            "cancellation": "changes each parallel support count by an even integer",
        },
        "closed_residue_table": {
            "5<=s<=m": {"w0": "m+5", "minimum_j": 2},
            "s in {4,m+1}": {"w0": "m+3", "minimum_j": 1},
            "s in {3,m+2}": {"w0": "m+1", "minimum_j": 0},
            "s in {2,m+3}": {"w0": "m-1", "minimum_j": 0},
            "s in {1,m+4}": {"w0": "m-3", "minimum_j": 0},
            "remaining residues": {"w0": "m-5", "minimum_j": 0},
        },
        "residue_definition": f"s=(t+1) mod (p+1), m={m}",
        "general_bound": (
            "|a(T_U)|>=max(1,w0-m), hence "
            "j>=max(0,(w0-m-1)/2)"
        ),
        "j0_excluded_residues": "4<=s<=m+1",
        "j1_fully_excluded_residues": "5<=s<=m",
        "j1_fixed_only_residues": "s in {4,m+1}",
        "sample": sample,
        "proved": True,
        "integral_transverse_fibre_solved": False,
        "residual_ii_closed": False,
    }


def main() -> None:
    print(json.dumps(theorem_record(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
