#!/usr/bin/env python3
"""Prop. 15.696 -- exclude the p=19 slack-twenty b=16 profile.

The phase-one floors ``9*18+38`` saturate the exact 200-unit type budget.
For the unique b=16 direction, equality fixes the slack on intersection
layers t=7,8,10.  Their quadratic evaluation map on J(19,10) has rank 169:
a fixed 169-row minor has rank 169 modulo two, while the two-dimensional
kernel consists of

    (alpha dot y_C) * (2-|y_C|),  sum(alpha)=0,

where C is the set of three empty fibres.  Integrality and parity on t=9
then leave, up to relabelling C={z,u,v}, two forms

    A_022 = 1-z+u+v-2uv,
    A_400 = 1+3z-u-v-2zu-2zv+2uv.

Coefficient comparison converts this and the nine rigid b=2 directions
into exact linear identities for an affine edge lift.  Aggregate signed
capacity and parallel-edge accounting leave ten possible infinity degrees.
The checked-in CP-SAT model exhausts every edge subset in each shard; all
twenty corrected shape/degree shards are INFEASIBLE. Finite-field differences
are formed componentwise in the fixed ``F_p`` basis and regression-tested
against the canonical Paley conference matrix. A nonsquare anti-isometry
transfers c_H=-1 to c_H=1, so the computation excludes both signs.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path

from e1_gmin_m4_prop15632 import eval_quadratic, parity_majorant_floor
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
P = 19
M = 10
B_SIZE = 16
PHASE = 1
H_SIZE = 77
PAIR_COLUMNS = tuple(combinations(range(P), 2))
EMPTY = (16, 17, 18)

# (infinity degree, phase-one gauge, phase-one parallel-edge count).
AGGREGATE_ROWS = (
    (2, 0, 2),
    (8, 1, 5),
    (10, 1, 3),
    (12, 1, 1),
    (18, 2, 4),
    (20, 2, 2),
    (28, 3, 3),
    (30, 3, 1),
    (38, 4, 2),
    (48, 5, 1),
)

# Greedy independent rows over F_2 from t=10, then t=8, then t=7.
# Together with the explicit two-dimensional Q-kernel this proves rank 169.
RANK_WITNESS_MASKS = (
    1023, 1535, 2559, 4607, 8703, 16895, 33279, 1791, 2815, 4863,
    8959, 17151, 33535, 3327, 5375, 9471, 17663, 34047, 6399, 10495,
    18687, 35071, 12543, 20735, 37119, 24831, 41215, 49407, 196863,
    327935, 393471, 1919, 2943, 4991, 9087, 17279, 33663, 3455,
    196991, 328063, 3711, 197247, 328319, 197759, 328831, 198783,
    329855, 200831, 331903, 204927, 335999, 213119, 344191, 229503,
    360575, 458879, 1983, 3007, 5055, 9151, 17343, 33727, 3519,
    197055, 328127, 3775, 458943, 3903, 459071, 459327, 459839,
    460863, 462911, 467007, 475199, 491583, 2015, 3039, 5087, 9183,
    17375, 33759, 3551, 197087, 328159, 3807, 458975, 3935, 3999,
    2031, 3055, 5103, 9199, 17391, 33775, 3567, 197103, 328175,
    3823, 458991, 3951, 4015, 4047, 2039, 3063, 5111, 9207, 17399,
    33783, 3575, 197111, 328183, 3831, 458999, 3959, 4023, 4055,
    4071, 2043, 3067, 5115, 9211, 17403, 33787, 3579, 197115,
    328187, 3835, 459003, 3963, 4027, 4059, 4075, 4083, 2045,
    3069, 5117, 9213, 17405, 33789, 3581, 197117, 328189, 3837,
    459005, 3965, 4029, 4061, 4077, 4085, 4089, 2046, 3070, 5118,
    9214, 17406, 33790, 3582, 197118, 328190, 3838, 459006, 3966,
    4030, 4062, 4078, 4086, 4090, 4092,
)


def _pair_row_bits(mask: int) -> int:
    row = 0
    for column, (left, right) in enumerate(PAIR_COLUMNS):
        if (mask >> left) & 1 and (mask >> right) & 1:
            row |= 1 << column
    return row


def _rank_mod_two(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def p19_b16_floor_and_kernel_certificate() -> dict[str, object]:
    """Exact floor equality, rank 169, and integral t=9 classification."""
    floor = parity_majorant_floor(P, B_SIZE, PHASE)
    coefficients = floor["coefficients"]
    assert isinstance(coefficients, tuple)
    q_values = {
        t: eval_quadratic(coefficients, t) for t in floor["support"]
    }
    expected_q = {7: Fraction(0), 8: Fraction(1), 9: Fraction(4, 3), 10: Fraction(1)}
    if q_values != expected_q or floor["contacts"] != (7, 8, 10):
        raise ArithmeticError("b=16 phase-one floor polynomial changed")

    masks = list(RANK_WITNESS_MASKS)
    layers = [sum((mask >> i) & 1 for i in range(B_SIZE)) for mask in masks]
    if len(masks) != 169 or any(mask.bit_count() != M for mask in masks):
        raise ArithmeticError("rank witness left J(19,10)")
    rank = _rank_mod_two([_pair_row_bits(mask) for mask in masks])
    if rank != 169:
        raise ArithmeticError("forced-layer rank witness lost rank")

    # Two independent rational kernel elements.  On the forced layers the
    # number k of selected empty fibres is 3,2,0; hence either 2-k=0 or the
    # alpha sum/linear form is zero.
    alphas = ((1, -1, 0), (1, 0, -1))
    for chosen in product((0, 1), repeat=3):
        k = sum(chosen)
        if k not in (0, 2, 3):
            continue
        for alpha in alphas:
            value = sum(a * y for a, y in zip(alpha, chosen)) * (2 - k)
            if value:
                raise ArithmeticError("displayed kernel failed a forced layer")

    t9_values = tuple(
        values
        for values in product(range(0, 6, 2), repeat=3)
        if sum(values) == 4
    )
    expected_t9 = set()
    for values in ((0, 2, 2), (4, 0, 0)):
        expected_t9.update(permutations(values))
    if set(t9_values) != expected_t9:
        raise ArithmeticError("integral t=9 classification changed")

    # In the first orbit z has value zero; in the second it has value four.
    normal_values_022 = {}
    normal_values_400 = {}
    for z, u, v in product((0, 1), repeat=3):
        k = z + u + v
        normal_values_022.setdefault(k, set()).add(1 - z + u + v - 2 * u * v)
        normal_values_400.setdefault(k, set()).add(
            1 + 3 * z - u - v - 2 * z * u - 2 * z * v + 2 * u * v
        )
    expected_layers = {0: {1}, 1: {0, 2}, 2: {1}, 3: {0}}
    expected_layers_400 = {0: {1}, 1: {0, 4}, 2: {1}, 3: {0}}
    if normal_values_022 != expected_layers or normal_values_400 != expected_layers_400:
        raise ArithmeticError("integral b=16 normal forms changed")

    def canonical_targets(
        linear: dict[int, int], quadratic: dict[tuple[int, int], int]
    ) -> tuple[int, ...]:
        row_sums = {
            s: sum(
                quadratic.get((min(s, t), max(s, t)), 0)
                for t in range(P) if t != s
            )
            for s in range(P)
        }
        shifts = {s: row_sums[s] + 2 * linear.get(s, 0) for s in range(P)}
        targets = []
        for s, t in PAIR_COLUMNS:
            numerator = quadratic.get((s, t), 0) + shifts[s] + shifts[t]
            if numerator & 1:
                raise ArithmeticError("normal form has nonintegral pair target")
            targets.append(numerator // 2)
        return tuple(targets)

    target_b2 = canonical_targets({0: -1, 1: -1}, {(0, 1): 2})
    target_022 = canonical_targets(
        {16: -1, 17: 1, 18: 1}, {(17, 18): -2}
    )
    target_400 = canonical_targets(
        {16: 3, 17: -1, 18: -1},
        {(16, 17): -2, (16, 18): -2, (17, 18): 2},
    )
    target_histograms = {
        "b2": dict(sorted(Counter(target_b2).items())),
        "022": dict(sorted(Counter(target_022).items())),
        "400": dict(sorted(Counter(target_400).items())),
    }
    if target_histograms != {
        "b2": {0: 170, 1: 1},
        "022": {-1: 19, 0: 152},
        "400": {-1: 35, 0: 120, 1: 16},
    }:
        raise ArithmeticError("normal-form pair targets changed")

    return {
        "p": P,
        "b": B_SIZE,
        "phase": PHASE,
        "phase_profile": {2: 9, 16: 1},
        "type_budget": 200,
        "type_floor_sum": 9 * 18 + 38,
        "floor_polynomial_coefficients": coefficients,
        "floor_polynomial_values": q_values,
        "forced_pointwise_layers": {7: 0, 8: 1, 10: 1},
        "pair_monomial_dimension": len(PAIR_COLUMNS),
        "rank_witness_rows": len(masks),
        "rank_witness_layer_histogram": dict(sorted(Counter(layers).items())),
        "rank_mod_two": rank,
        "explicit_rational_kernel_dimension": len(alphas),
        "kernel": "(alpha dot y_C)(2-|y_C|), sum(alpha)=0",
        "therefore_rational_rank": 169,
        "t9_nonnegative_even_value_triples": t9_values,
        "normal_forms": {
            "022": "A=1-z+u+v-2uv",
            "400": "A=1+3z-u-v-2zu-2zv+2uv",
        },
        "normal_form_values_by_selected_empty_count": {
            "022": {key: sorted(value) for key, value in normal_values_022.items()},
            "400": {key: sorted(value) for key, value in normal_values_400.items()},
        },
        "canonical_pair_target_histograms": target_histograms,
        "canonical_pair_target_sums": {
            "b2": sum(target_b2),
            "022": sum(target_022),
            "400": sum(target_400),
        },
        "proved": True,
    }


def _phase_zero_parallel_options(infinity_degree: int) -> dict[tuple[int, int], tuple[int, ...]]:
    """Parallel counts surviving the exact signed-cross capacity bound."""
    finite_edges = H_SIZE - infinity_degree
    options: dict[tuple[int, int], tuple[int, ...]] = {}
    for role, base_mean in ((0, 0), (2, 20), (16, 40)):
        for elevated in (0, 1):
            mean = base_mean + 20 * elevated
            allowed = []
            for parallel in range(finite_edges + 1):
                signed_cross = infinity_degree + P * parallel - 3 * P - mean
                if abs(signed_cross) <= finite_edges - parallel:
                    allowed.append(parallel)
            options[(role, elevated)] = tuple(allowed)
    return options


def _phase_zero_parallel_sum_possible(
    infinity_degree: int, target_sum: int
) -> bool:
    options = _phase_zero_parallel_options(infinity_degree)
    # State is (number b0, number b2, number b16, elevations, parallel sum).
    states = {(0, 0, 0, 0, 0)}
    limits = {0: 5, 2: 1, 16: 4}
    index = {0: 0, 2: 1, 16: 2}
    for _ in range(10):
        next_states = set()
        for counts in states:
            for role in (0, 2, 16):
                updated = list(counts[:3])
                updated[index[role]] += 1
                if updated[index[role]] > limits[role]:
                    continue
                for elevated in (0, 1):
                    if counts[3] + elevated > 1:
                        continue
                    for parallel in options[(role, elevated)]:
                        total = counts[4] + parallel
                        if total <= target_sum:
                            next_states.add((*updated, counts[3] + elevated, total))
        states = next_states
    return (5, 1, 4, 1, target_sum) in states


def p19_b16_aggregate_degree_certificate() -> dict[str, object]:
    """Recompute the complete ten-row aggregate shard table."""
    rows = []
    for infinity_degree in range(H_SIZE + 1):
        # Infinity is not in the boundary, while the normalized infinity-zero
        # edge is selected.
        if infinity_degree == 0 or infinity_degree & 1:
            continue
        for gauge in range(-H_SIZE, H_SIZE + 1):
            parallel = 4 + 9 * gauge - infinity_degree
            cross_edges = 73 - 9 * gauge
            if parallel < 0 or cross_edges < 0:
                continue
            # Sum the 171 phase-one coefficient identities.  The b=2 target
            # has total +1; the b=16 target has total -19.
            b2_signed = 1 + 171 * gauge - 18 * infinity_degree
            b16_signed = -19 + 171 * gauge - 18 * infinity_degree
            if abs(b2_signed) > cross_edges or abs(b16_signed) > cross_edges:
                continue
            finite_edges = H_SIZE - infinity_degree
            phase_zero_target = finite_edges - 10 * parallel
            if phase_zero_target < 0:
                continue
            if _phase_zero_parallel_sum_possible(infinity_degree, phase_zero_target):
                rows.append((infinity_degree, gauge, parallel))
    if tuple(rows) != AGGREGATE_ROWS:
        raise ArithmeticError("aggregate infinity-degree table changed")
    return {
        "normalization": "infinity-zero edge selected",
        "infinity_outside_boundary": True,
        "phase_one_parallel_identity": "P_d=4+9*g_d-I",
        "phase_one_signed_sums": {
            "b2": "1+171*g_d-18*I",
            "b16": "-19+171*g_d-18*I",
        },
        "phase_zero_profile": {0: 5, 2: 1, 16: 4},
        "phase_zero_elevation_count": 1,
        "admissible_rows": rows,
        "admissible_infinity_degrees": [row[0] for row in rows],
        "proved": True,
    }


def p19_b16_solver_shard_certificate() -> dict[str, object]:
    """Audit both normal-form orbits across all ten exact shards."""
    directory = ROOT / "evidence" / "p19_slack20_b16_lift_shards"
    expected = [row[0] for row in AGGREGATE_ROWS]
    shards = []
    for shape in ("022", "400"):
        for infinity_degree in expected:
            raw_shards = (
                tuple(
                    (
                        f"p19_slack20_b16_022_i28_role{role}_correct.json",
                        role,
                    )
                    for role in (0, 2, 16)
                )
                if shape == "022" and infinity_degree == 28
                else (
                    (
                        f"p19_slack20_b16_{shape}_i{infinity_degree}_correct.json",
                        None,
                    ),
                )
            )
            for suffix, elevated_role in raw_shards:
                path = directory / suffix
                raw = path.read_bytes()
                payload = json.loads(raw)
                required = {
                    "experiment": "p19_slack20_b16_lift_cpsat",
                    "status": "exact_affine_edge_lift_model",
                    "p": P,
                    "c_H": 1,
                    "fixed_infinity_degree": infinity_degree,
                    "phase_zero_profile": {"0": 5, "2": 1, "16": 4},
                    "phase_one_profile": {"2": 9, "16": 1},
                    "b16_shape": shape,
                    "edge_variables": 65341,
                    "rigid_coefficient_identities": 1720,
                    "solver_status": "INFEASIBLE",
                    "feasible": False,
                    "finite_infeasibility_certificate": True,
                }
                if elevated_role is not None:
                    required["fixed_phase_zero_elevated_role"] = elevated_role
                for key, value in required.items():
                    if payload.get(key) != value:
                        raise ArithmeticError(
                            f"solver shard {shape}/{infinity_degree} failed {key}"
                        )
                shards.append(
                    {
                        "shape": shape,
                        "infinity_degree": infinity_degree,
                        "phase_zero_elevated_role": elevated_role,
                        "file": suffix,
                        "sha256": hashlib.sha256(raw).hexdigest(),
                        "solver_status": payload["solver_status"],
                        "wall_time_seconds": payload["wall_time_seconds"],
                        "branches": payload["branches"],
                        "conflicts": payload["conflicts"],
                    }
                )
    return {
        "solver": "OR-Tools CP-SAT 9.15.6755",
        "python_runtime_recorded_in_raw_shards": False,
        "machines": [
            "soulkiller x86_64",
            "jellyfin x86_64",
        ],
        "model": "scripts/p19_slack20_b16_lift_cpsat.py",
        "finite_field_sign_convention": (
            "componentwise subtraction in the fixed F_p basis; full edge "
            "tables regression-tested against paley_conference_prime_power"
        ),
        "supersedes_original_raw_shards": True,
        "normal_form_orbits": ["022", "400"],
        "shards_per_orbit": len(expected),
        "shard_count": 2 * len(expected),
        "raw_shard_count": len(shards),
        "split_logical_shard": {
            "shape": "022",
            "infinity_degree": 28,
            "exhaustive_phase_zero_elevated_roles": [0, 2, 16],
        },
        "all_statuses": "INFEASIBLE",
        "shards": shards,
        "c_h_plus_one_excluded": True,
        "c_h_minus_one_transfer": (
            "multiplication by a nonsquare fixes infinity and zero, flips "
            "both direction type and c_H, and preserves every phase/profile"
        ),
        "both_c_h_signs_excluded": True,
        "proved_computationally": True,
    }


def p19_slack_twenty_b16_exclusion() -> dict[str, object]:
    """Proposition 15.696."""
    previous = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 20
        and not (
            int(row["pair_slack"]) == 20
            and row["phase_profiles_b"]["1"] == {2: 9, 14: 1}
        )
    ]
    excluded = [
        row for row in previous
        if row["phase_profiles_b"]["1"] == {2: 9, 16: 1}
    ]
    remaining = [row for row in previous if row not in excluded]
    histogram = dict(sorted(Counter(int(row["pair_slack"]) for row in remaining).items()))
    if len(previous) != 5 or len(excluded) != 1 or histogram != {20: 1, 24: 1, 28: 1, 32: 1}:
        raise ArithmeticError("post-15.696 p=19 remainder changed")
    return {
        "proposition": "15.696",
        "p": P,
        "boundary_size": 16,
        "floor_kernel_certificate": p19_b16_floor_and_kernel_certificate(),
        "aggregate_degree_certificate": p19_b16_aggregate_degree_certificate(),
        "solver_shard_certificate": p19_b16_solver_shard_certificate(),
        "excluded_profile": excluded[0],
        "p19_profiles_before": len(previous),
        "p19_profiles_after": len(remaining),
        "remaining_slack_histogram": histogram,
        "remaining_profiles": remaining,
        "p19_second_all_finite_endpoint_closed": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved_computationally": True,
    }


def main() -> None:
    theorem = p19_slack_twenty_b16_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15696.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True, default=str) + "\n")
    print(
        "Prop. 15.696: excluded the p=19 slack-twenty b=16 profile; "
        "four p=19 profiles remain"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
