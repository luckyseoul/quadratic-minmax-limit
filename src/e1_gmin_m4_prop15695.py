#!/usr/bin/env python3
"""Prop. 15.695 -- exclude the two p=19 slack-twenty b=14 profiles.

In either profile the phase-one type has nine ``b=2`` directions and one
``b=14`` direction.  Their floor sum is the complete type budget

    9*18 + 38 = 200,

so the b=14 direction attains its exact floor.  Its symmetrized quadratic
has expectation one.  The positive quadrature for ``(p,b,phase)=(19,14,1)``
is supported with positive weights on intersection layers ``t=6,8,10``.
All values on those layers are positive odd integers, hence equality forces
the original directional slack to equal one pointwise there.

On ``J(19,10)``, pair monomials span every degree-at-most-two function.  The
fixed 171-row inclusion submatrix below, drawn from the three forced layers,
has rank 171 modulo 101.  It therefore has full rank over Q.  A quadratic
equal to one on those layers is consequently identically one on the entire
slice, contradicting the required even parity on the nonempty layer t=5.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15669 import middle_floor_quadrature
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
P = 19
M = 10
B_SIZE = 14
PHASE = 1
MODULUS = 101

# A square full-rank minor of the pair-vs.-10-subset inclusion matrix.
# Masks 0..13 are the b=14 odd-fibre set.  The rows use 91 masks from
# t=10, 75 from t=8, and 5 from t=6.
RANK_WITNESS_MASKS = (
    1023, 1535, 2559, 4607, 8703, 1791, 2815, 4863, 8959, 3327, 5375,
    9471, 6399, 10495, 12543, 1919, 2943, 4991, 9087, 3455, 3711, 1983,
    3007, 5055, 9151, 3519, 3775, 3903, 2015, 3039, 5087, 9183, 3551,
    3807, 3935, 3999, 2031, 3055, 5103, 9199, 3567, 3823, 3951, 4015,
    4047, 2039, 3063, 5111, 9207, 3575, 3831, 3959, 4023, 4055, 4071,
    2043, 3067, 5115, 9211, 3579, 3835, 3963, 4027, 4059, 4075, 4083,
    2045, 3069, 5117, 9213, 3581, 3837, 3965, 4029, 4061, 4077, 4085,
    4089, 2046, 3070, 5118, 9214, 3582, 3838, 3966, 4030, 4062, 4078,
    4086, 4090, 4092, 49407, 82175, 147711, 278783, 98559, 164095,
    295167, 196863, 327935, 393471, 49535, 82303, 147839, 278911, 98687,
    49791, 82559, 148095, 279167, 98943, 50303, 83071, 148607, 279679,
    99455, 51327, 84095, 149631, 280703, 100479, 53375, 86143, 151679,
    282751, 102527, 57471, 90239, 155775, 286847, 106623, 49599, 82367,
    147903, 278975, 98751, 49631, 82399, 147935, 279007, 98783, 49647,
    82415, 147951, 279023, 98799, 49655, 82423, 147959, 279031, 98807,
    49659, 82427, 147963, 279035, 98811, 49661, 82429, 147965, 279037,
    98813, 49662, 82430, 147966, 279038, 98814, 245823, 376895, 442431,
    475199, 491583,
)


def _rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    """Exact row rank by finite-field elimination."""
    matrix = [[value % prime for value in row] for row in rows]
    rank = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def p19_b14_layer_rank_certificate() -> dict[str, object]:
    """Full-rank certificate for the t=6,8,10 determining layers."""
    pairs = tuple(combinations(range(P), 2))
    rows = [
        [int((mask >> a) & 1 and (mask >> b) & 1) for a, b in pairs]
        for mask in RANK_WITNESS_MASKS
    ]
    weights = [mask.bit_count() for mask in RANK_WITNESS_MASKS]
    layers = [
        (mask & ((1 << B_SIZE) - 1)).bit_count()
        for mask in RANK_WITNESS_MASKS
    ]
    layer_histogram = dict(sorted(Counter(layers).items()))
    rank = _rank_mod_prime(rows, MODULUS)
    dimension = comb(P, 2)
    if len(rows) != dimension or set(weights) != {M}:
        raise ArithmeticError("rank witness left J(19,10)")
    if layer_histogram != {6: 5, 8: 75, 10: 91}:
        raise ArithmeticError("rank witness layer counts changed")
    if rank != dimension:
        raise ArithmeticError("three-layer quadratic evaluation lost rank")
    return {
        "slice": "J(19,10)",
        "degree_at_most_two_dimension": dimension,
        "pair_monomial_column_count": len(pairs),
        "witness_row_count": len(rows),
        "witness_layer_histogram": layer_histogram,
        "finite_field_modulus": MODULUS,
        "finite_field_rank": rank,
        "therefore_rational_rank": dimension,
        "pair_monomials_span_degree_at_most_two_on_slice": True,
        "proved": True,
    }


def p19_b14_floor_equality_ledger() -> dict[str, object]:
    """Positive-quadrature equality and its pointwise consequences."""
    quadrature = middle_floor_quadrature(P, B_SIZE, PHASE)
    weights = quadrature["quadrature_weights"]
    if set(weights) != {6, 8, 10} or not all(weight > 0 for weight in weights.values()):
        raise ArithmeticError("b=14 equality quadrature changed")
    if int(quadrature["scaled_floor"]) != 38:
        raise ArithmeticError("b=14 phase-one floor changed")
    if 9 * 18 + 38 != (P + 1) ** 2 // 2:
        raise ArithmeticError("phase-one floor saturation changed")
    return {
        "p": P,
        "b": B_SIZE,
        "phase": PHASE,
        "phase_profile": {2: 9, 14: 1},
        "type_budget": (P + 1) ** 2 // 2,
        "type_floor_sum": 9 * 18 + 38,
        "b14_scaled_floor": 38,
        "symmetrized_minimizer": "q(t)=1",
        "positive_quadrature_weights": weights,
        "forced_pointwise_one_layers": sorted(weights),
        "contradicting_even_parity_layer": 5,
        "contradicting_layer_nonempty": comb(B_SIZE, 5) * comb(P - B_SIZE, 5) > 0,
        "proved": True,
    }


def p19_slack_twenty_b14_exclusion() -> dict[str, object]:
    """Proposition 15.695."""
    slack_twenty = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) == 20
    ]
    excluded = [
        row
        for row in slack_twenty
        if row["phase_profiles_b"]["1"] == {2: 9, 14: 1}
    ]
    remaining = [row for row in slack_twenty if row not in excluded]
    if len(slack_twenty) != 4 or len(excluded) != 2 or len(remaining) != 2:
        raise ArithmeticError("slack-twenty b=14 profile split changed")
    if Counter(int(row["undetermined_directions"]) for row in excluded) != {4: 2}:
        raise ArithmeticError("excluded profiles are no longer the two t=4 rows")
    high_slack = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 20
    ]
    remaining_histogram = dict(
        sorted(
            Counter(
                int(row["pair_slack"])
                for row in high_slack
                if row not in excluded
            ).items()
        )
    )
    if len(high_slack) != 7 or remaining_histogram != {20: 2, 24: 1, 28: 1, 32: 1}:
        raise ArithmeticError("post-15.695 p=19 remainder changed")
    return {
        "proposition": "15.695",
        "p": P,
        "boundary_size": 16,
        "floor_equality": p19_b14_floor_equality_ledger(),
        "three_layer_rank_certificate": p19_b14_layer_rank_certificate(),
        "quadratic_identity_contradiction": (
            "A=1 on t=6,8,10 forces A=1 on J(19,10), "
            "but phase-one parity makes A even on nonempty t=5"
        ),
        "slack_twenty_profiles_before": len(slack_twenty),
        "excluded_b14_profiles": len(excluded),
        "excluded_profiles": excluded,
        "slack_twenty_profiles_after": len(remaining),
        "remaining_slack_twenty_profiles": remaining,
        "p19_profiles_before": 7,
        "p19_profiles_after": 5,
        "remaining_slack_histogram": remaining_histogram,
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def main() -> None:
    theorem = p19_slack_twenty_b14_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15695.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True, default=str) + "\n")
    print(
        "Prop. 15.695: excluded both p=19 slack-twenty b=14 profiles; "
        "five p=19 profiles remain"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
