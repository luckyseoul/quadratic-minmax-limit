#!/usr/bin/env python3
"""Cardinality barrier for the direction-sliced halved Boolean fibre.

This is a negative theorem about what punctured mod-two surjectivity can
prove.  It does not decide the particular branch-C target syndrome.

Put ``p=2h+1``, ``d=p+1`` and ``N=d*h=|Delta|``.  The halved symmetric map
has codomain dimension

    R = d*h*(h+1).

Its columns split into ``d`` groups according to their unique parallel
direction.  Under the fixed-edge inverse, the parallel functionals are the
sums on the ``d`` disjoint direction blocks which partition ``Delta``; they
are therefore independent.  If the punctured map is onto, then fixing these
``d`` bits leaves exactly ``2**(R-d)`` target syndromes.  On the other hand, a
product of exact direction-weight slices of total weight
``s <= (N-1)/2`` contains strictly fewer than ``2**(R-d)`` Boolean points
for every ``h>=7``.  Consequently onto plus the scalar quota bounds cannot
imply that every compatible target has a point in its prescribed slice.
"""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15721 import is_prime


def _check_branch_c_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a branch-C prime p=3 mod 4 with p>=31")
    return (p - 1) // 2


def quota_cardinality_barrier(
    p: int,
    quotas: tuple[int, ...] | list[int],
    unused_capacities: tuple[int, ...] | list[int],
) -> dict[str, object]:
    """Prove that one feasible low-total quota slice cannot cover a parity fibre.

    ``unused_capacities[L]`` is the number of retained columns in parallel
    direction ``L`` and ``quotas[L]`` is their prescribed Hamming weight.
    The only map hypothesis used in the conclusion is that the punctured
    halved map ``D_U`` is onto its full ``R``-dimensional codomain.

    The function computes the exact slice cardinality for auditability, but
    the proof uses only Vandermonde and the symbolic upper bounds recorded in
    the return value.  It builds no Radon matrix and performs no prime census.
    """
    h = _check_branch_c_prime(p)
    d = p + 1
    n_delta = d * h
    codomain_rank = d * h * (h + 1)

    quotas = tuple(quotas)
    unused_capacities = tuple(unused_capacities)
    if len(quotas) != d or len(unused_capacities) != d:
        raise ValueError("quota and capacity vectors must have p+1 entries")
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for value in quotas + unused_capacities
    ):
        raise ValueError("quotas and capacities must be nonnegative integers")
    if any(
        quota > capacity
        for quota, capacity in zip(quotas, unused_capacities, strict=True)
    ):
        raise ValueError("every quota must be at most its unused capacity")

    total_weight = sum(quotas)
    if 2 * total_weight > n_delta - 1:
        raise ValueError("the theorem requires total quota at most (|Delta|-1)/2")
    if sum(unused_capacities) > n_delta * n_delta:
        raise ValueError("unused capacities cannot exceed the full column count")

    exact_slice_size = 1
    for capacity, quota in zip(unused_capacities, quotas, strict=True):
        exact_slice_size *= comb(capacity, quota)

    fixed_parallel_parity_fibre_size = 1 << (codomain_rank - d)
    # For h>=7, N=2h(h+1)<2^h (base h=7, then induction).  If s>=1,
    # |slice| <= C(N^2,s) <= (N^2)^s < 2^(2hs)
    #          <= 2^(h(N-1)) < 2^(R-d).
    # If s=0, the slice is the singleton zero word and the final strict
    # comparison is immediate.
    delta_below_two_to_h = n_delta < (1 << h)
    symbolic_slice_exponent = h * (n_delta - 1)
    parity_fibre_exponent = codomain_rank - d
    exponent_gap = parity_fibre_exponent - symbolic_slice_exponent
    exact_gap = fixed_parallel_parity_fibre_size - exact_slice_size

    proved = bool(
        h >= 7
        and delta_below_two_to_h
        and exponent_gap > 0
        and exact_slice_size < fixed_parallel_parity_fibre_size
    )
    if not proved:
        raise ArithmeticError("the quota cardinality barrier changed")

    return {
        "p": p,
        "h": h,
        "direction_count_d": d,
        "Delta_size_N": n_delta,
        "full_halved_codomain_rank_R": codomain_rank,
        "parallel_rows_fixed": d,
        "fixed_parallel_parity_fibre_exponent": parity_fibre_exponent,
        "fixed_parallel_parity_fibre_size": fixed_parallel_parity_fibre_size,
        "quota_total_s": total_weight,
        "quota_total_bound": "2s<=N-1",
        "unused_column_total": sum(unused_capacities),
        "exact_quota_slice_size": exact_slice_size,
        "symbolic_chain": (
            "if s=0 the slice has size one; if s>=1, "
            "prod_L binom(c_L,n_L) <= binom(sum c_L,s) <= "
            "binom(N^2,s) <= (N^2)^s < 2^(2hs) <= "
            "2^(h(N-1)) < 2^(R-d)"
        ),
        "N_below_2_to_h": delta_below_two_to_h,
        "symbolic_slice_exponent_upper": symbolic_slice_exponent,
        "parity_fibre_exponent_gap": exponent_gap,
        "exact_missing_target_count_lower_bound": exact_gap,
        "conditional_map_hypothesis": "D_U is onto",
        "conclusion": (
            "for these exact feasible quotas, some target syndrome with the "
            "same parallel-row parity has no Boolean preimage in the quota slice"
        ),
        "onto_plus_quota_bounds_implies_slice_nonempty_for_every_target": False,
        "actual_branch_C_target_excluded": False,
        "proved": proved,
    }


def branch_c_uniform_information_barrier(p: int) -> dict[str, object]:
    """State the uniform branch-C consequence without choosing actual quotas.

    Along the balanced branch-C ray the target graph has at most ``N-1``
    edges.  Fixed-edge elimination therefore gives total unused-pair quota
    ``s<=(N-1)/2``.  For every feasible directionwise quota vector in this
    range, the theorem above applies.  This wrapper records the symbolic
    comparison only; it does not assert that the actual syndrome is one of
    the missing syndromes.
    """
    h = _check_branch_c_prime(p)
    d = p + 1
    n_delta = d * h
    codomain_rank = d * h * (h + 1)
    maximum_quota_total = (n_delta - 1) // 2
    symbolic_slice_exponent = h * (n_delta - 1)
    parity_fibre_exponent = codomain_rank - d
    proved = bool(
        n_delta == (p * p - 1) // 2
        and n_delta < (1 << h)
        and symbolic_slice_exponent < parity_fibre_exponent
    )
    if not proved:
        raise ArithmeticError("the uniform branch-C information barrier changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "N": n_delta,
        "branch_C_H_max": n_delta - 1,
        "maximum_unused_pair_quota_total": maximum_quota_total,
        "R": codomain_rank,
        "R_minus_d": parity_fibre_exponent,
        "slice_log2_strict_upper_bound": symbolic_slice_exponent,
        "exponent_gap": parity_fibre_exponent - symbolic_slice_exponent,
        "scope": (
            "all feasible quota vectors arising after fixed-edge elimination; "
            "target syndrome otherwise arbitrary"
        ),
        "methodological_conclusion": (
            "punctured surjectivity and scalar direction-capacity inequalities "
            "alone cannot establish the prescribed-weight Boolean lift"
        ),
        "additional_needed_input": (
            "specific transverse-target structure, a target-sensitive kernel "
            "exchange theorem, or a direct target obstruction"
        ),
        "actual_target_status": "OPEN",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    return branch_c_uniform_information_barrier(p)


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
