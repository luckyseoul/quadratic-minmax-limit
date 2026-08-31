#!/usr/bin/env python3
r"""Prop. 15.737 -- close the first three ``p=11`` residual layers.

Proposition 15.736 gives a self-contained exhaustive catalog of the sharp
Boolean quadratic lifts on ``J(11,6)``.  For ``t in {0,1,2}``, consider

    k=44+2t,                 |H|=45+2t.

There are respectively ``32,28,24`` isolated projective vertices beyond the
support bound.  Signed PSL transport therefore gives an all-finite chart with
``I=0`` and every directional ``b`` even.

The phase-one type has exact budget ``12(6+t)``.  Writing
``a_d=2u+12k_d`` gives ``sum k_d=6+t-u``.  The residues ``u=0,...,t`` force
a mean-``12+2u`` direction, a nonzero lift of excess at most six below the
sharp floor eight.  Residues ``t<u<5`` have too few quotient units.  Thus
``u=5``: at least ``5-t`` directions have exact mean ten.  Positive
quadrature fixes their baselines to ``b=2`` or ``b=10``, and coefficient
offsets four and three prevent mixing.

The 15.736 catalog excludes the hard-``b=2`` branch.  In the hard-``b=10``
branch, at least ``4-t`` opposite directions have the all-equal-triple target

    eps_L S_H = 4 + z_i*z_j + z_i*z_k + z_j*z_k.        (1)

There is a short cross-direction obstruction to (1).  Work in the isolated
chart of Proposition 15.734, so every selected edge is finite.  For an
``F_11``-linear fibre functional ``L`` define

    M_H(L) = sum_({u,v} in H) chi(u-v) (L(u)-L(v))^2.   (2)

This is a homogeneous binary quadratic in the two coefficients of ``L``.
If ``K^L_st`` is the signed selected-edge sum between the fibres ``s,t``,
then

    M_H(L) = sum_(s<t) K^L_st (s-t)^2,                 (3)

because edges parallel to ``ker L`` contribute zero.

The hard branch has at least ``5-t>=3`` low directions with ``P=3`` and exact target
``eps_L S_H=4-z_j``.  Coefficient comparison on ``sum z_s=1`` gives

    eps_L K^L_st = -1  if exactly one of s,t equals j,
                       0  otherwise.                  (4)

Consequently

    M_H(L) = -eps_L sum_(t != j) (j-t)^2 = 0           (5)

in ``F_11``.  These low directions are distinct projective linear forms.
A nonzero homogeneous binary quadratic has at most two projective
zeros, so (5) forces ``M_H`` to be the zero quadratic.

For an opposite all-equal target (1), ``P=4`` and coefficient comparison
instead gives ``eps_L K^L_st=1`` on the triangle ``{i,j,k}`` and zero off
it.  Thus ``M_H=0`` would require

    (i-j)^2 + (i-k)^2 + (j-k)^2 = 0.                  (6)

For three distinct fibres translate and scale to ``(i,j,k)=(0,1,r)``.
The left side is ``2(r^2-r+1)``.  Its discriminant is ``-3=8``, a
nonsquare modulo 11, so (6) has no solution.  Even one all-equal target is
impossible, contradicting the at-least-``4-t>=2`` opposite targets.

Therefore residual (ii) is empty at ``p=11`` for ``k=44,46,48``.  This is a
proved theorem whose only finite-certificate dependency is Proposition
15.736's exact Boolean catalog.  It stops at ``t=3``: there are then only two
forced hard stars, and excess eight reaches the equality floor.  It does not
treat ``p=5,7``, ``k>=50`` at ``p=11``, or residual (ii) as a whole.
"""
from __future__ import annotations

import json
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from typing import Sequence

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import signed_relative_flip_transport
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P3_LAST
from e1_gmin_m4_prop15736 import proposition_15736, residual_p11_consequence


ROOT = Path(__file__).resolve().parents[1]
P = 11
Q = 5
M = 6
CRITICAL_K = 44
H_EDGE_COUNT = 45
CLOSED_LAYER_INDICES = (0, 1, 2)


def _check_layer_index(t: int) -> None:
    if not isinstance(t, int) or isinstance(t, bool) or t not in CLOSED_LAYER_INDICES:
        raise ValueError("need layer index t in {0,1,2}")


@lru_cache(maxsize=None)
def p11_isolated_layer_chart(t: int) -> dict[str, object]:
    """Replay the isolated-chart argument at ``k=44+2t``."""
    _check_layer_index(t)
    edge_count = H_EDGE_COUNT + 2 * t
    ambient_vertices = P * P + 1
    maximum_nonisolated = 2 * edge_count
    isolated_gap = ambient_vertices - maximum_nonisolated
    transport = signed_relative_flip_transport()
    proved = bool(
        isolated_gap == 32 - 4 * t
        and isolated_gap > 0
        and transport["proved"]
        and transport["flip_set_size_preserved"]
        and transport["odd_degree_boundary_is_permuted"]
        and transport["both_separation_inequalities_preserved"]
    )
    if not proved:
        raise ArithmeticError("the p=11 isolated-layer chart failed")
    return {
        "p": P,
        "layer_index_t": t,
        "original_k": CRITICAL_K + 2 * t,
        "H_edge_count": edge_count,
        "ambient_vertex_count": ambient_vertices,
        "maximum_nonisolated_vertices": maximum_nonisolated,
        "guaranteed_isolated_vertices": isolated_gap,
        "isolated_vertex_is_outside_odd_degree_boundary": True,
        "signed_PSL_transport_dependency": transport,
        "transported_infinity_degree_I": 0,
        "transported_boundary_is_all_finite": True,
        "boundary_size_even_by_handshake": True,
        "every_directional_odd_fibre_count_b_is_even": True,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def p11_phase_one_residue_replay(t: int) -> dict[str, object]:
    """Derive the hard low-baseline dichotomy directly at p=11."""
    _check_layer_index(t)
    floors = {b: scaled_direction_floor(P, b, 1) for b in range(0, P, 2)}
    expected_floors = {0: 22, 2: 10, 4: 22, 6: 18, 8: 22, 10: 10}
    if floors != expected_floors:
        raise ArithmeticError("the exact p=11 phase-one floor ledger changed")
    lift = sharp_integral_quadratic_lift_floor(P)
    lift_floor = int(lift["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for u in range(M):
        quotient_sum = M + t - u
        if u <= t:
            low_mean = P + 1 + 2 * u
            extra_units_above_one = t - u
            low_direction_count = M - extra_units_above_one
            compatible = [b for b, floor in floors.items() if floor <= low_mean]
            excesses = {b: low_mean - floors[b] for b in compatible}
            excluded = bool(
                compatible == [2, 10]
                and low_direction_count >= 1
                and all(0 < excess < lift_floor for excess in excesses.values())
            )
            reason = (
                "a k=1 direction is forced and every compatible cell is a "
                "nonzero lift of excess 2+2u below the sharp floor 8"
            )
        elif u < M - 1:
            low_mean = 2 * u
            low_direction_count = 0
            compatible = []
            excesses = {}
            excluded = quotient_sum < M and low_mean < min(floors.values())
            reason = "k=0 is below every floor, but sum k is less than six"
        else:
            low_mean = P - 1
            low_direction_count = M - (t + 1)
            compatible = [b for b, floor in floors.items() if floor <= low_mean]
            excesses = {b: low_mean - floors[b] for b in compatible}
            excluded = False
            reason = (
                "sum k=t+1, so at least 5-t directions have k=0 and exact "
                "mean ten"
            )
        rows.append(
            {
                "u": u,
                "common_residue": 2 * u,
                "quotient_sum": quotient_sum,
                "forced_low_mean": low_mean,
                "forced_low_direction_count": low_direction_count,
                "floor_compatible_even_b": compatible,
                "excess_above_floor": excesses,
                "excluded": excluded,
                "reason": reason,
            }
        )

    b2_quadrature = parity_floor_certificate(P, 2, 1)
    complementary_b1_quadrature = parity_floor_certificate(P, 1, 1)
    positive_quadrature_rigidity = bool(
        b2_quadrature["exact_positive_quadrature_certificate"]
        and complementary_b1_quadrature["exact_positive_quadrature_certificate"]
        and all(weight > 0 for weight in b2_quadrature["quadrature_weights"])
        and all(
            weight > 0
            for weight in complementary_b1_quadrature["quadrature_weights"]
        )
    )
    feasible_u = [int(row["u"]) for row in rows if not row["excluded"]]
    endpoint = rows[M - 1]
    proved = bool(
        p11_isolated_layer_chart(t)["proved"]
        and 12 * (M + t) == 2 * M * (M + t)
        and feasible_u == [M - 1]
        and endpoint["floor_compatible_even_b"] == [2, 10]
        and int(endpoint["forced_low_direction_count"]) == 5 - t
        and positive_quadrature_rigidity
        and lift_floor == 8
    )
    if not proved:
        raise ArithmeticError("the p=11 hard residue replay failed")
    b2_quadrature_summary = {
        "p": P,
        "b": 2,
        "phase": 1,
        "scaled_floor": int(b2_quadrature["scaled_floor"]),
        "quadrature_nodes": list(b2_quadrature["quadrature_nodes"]),
        "quadrature_weights": [
            str(weight) for weight in b2_quadrature["quadrature_weights"]
        ],
        "all_quadrature_weights_strictly_positive": all(
            weight > 0 for weight in b2_quadrature["quadrature_weights"]
        ),
        "exact_positive_quadrature_certificate": bool(
            b2_quadrature["exact_positive_quadrature_certificate"]
        ),
    }
    complementary_quadrature_summary = {
        "p": P,
        "b": 1,
        "phase": 1,
        "scaled_floor": int(complementary_b1_quadrature["scaled_floor"]),
        "quadrature_nodes": list(complementary_b1_quadrature["quadrature_nodes"]),
        "quadrature_weights": [
            str(weight)
            for weight in complementary_b1_quadrature["quadrature_weights"]
        ],
        "all_quadrature_weights_strictly_positive": all(
            weight > 0
            for weight in complementary_b1_quadrature["quadrature_weights"]
        ),
        "exact_positive_quadrature_certificate": bool(
            complementary_b1_quadrature["exact_positive_quadrature_certificate"]
        ),
    }
    return {
        "p": P,
        "layer_index_t": t,
        "H_edge_count": H_EDGE_COUNT + 2 * t,
        "phase_exponent": 21 + t,
        "phase_one_type": "eps_d=(-1)^t*c_H",
        "hard_type_budget": 12 * (M + t),
        "hard_type_budget_formula": "2m(m+t)",
        "same_type_mean_form": "a_d=2u+12*k_d",
        "quotient_identity": "sum_d k_d=6+t-u",
        "exact_phase_one_even_floors": floors,
        "residue_rows": rows,
        "feasible_u": feasible_u,
        "endpoint_u": M - 1,
        "endpoint_low_mean": P - 1,
        "endpoint_low_direction_count_at_least": 5 - t,
        "endpoint_low_b_candidates": [2, 10],
        "b2_exact_baseline": "A=(1-x_i-x_j)^2; eps*S_H=4+z_i*z_j",
        "b10_exact_baseline": "A=1-x_j; eps*S_H=4-z_j",
        "b2_coefficient_offset": 4,
        "b10_coefficient_offset": 3,
        "equal_mean_forces_equal_parallel_count": True,
        "offsets_differ_modulo_q_so_baselines_cannot_mix": True,
        "positive_quadrature_dependency": {
            "proposition": "15.652",
            "b2": b2_quadrature_summary,
            "complementary_b1": complementary_quadrature_summary,
        },
        "nonzero_integral_lift_floor": lift_floor,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def p11_branch_parallel_replay(t: int, branch: str) -> dict[str, object]:
    """Replay the exact hard/opposite counts before using the catalog."""
    _check_layer_index(t)
    if branch not in (BRANCH_B2, BRANCH_P3_LAST):
        raise ValueError("branch must be hard_b2 or p3_all_low_b_p_minus_1")
    if branch == BRANCH_B2:
        offset, expected_p, minimum_q, mean_constant, edge_delta = 4, 4, 3, 9, t
        surviving_family: str | None = None
    else:
        offset, expected_p, minimum_q, mean_constant, edge_delta = 3, 3, 4, 7, t + 1
        surviving_family = "all_equal_triple"
    parameter_rows: list[dict[str, object]] = []
    for parallel_count in range(9):
        numerator = parallel_count - offset
        congruent = numerator % Q == 0
        rho = numerator // Q if congruent else None
        nonnegative_rho = rho is not None and rho >= 0
        s = parallel_count + int(rho) if nonnegative_rho else None
        opposite_edges = Q * (8 - int(s)) + edge_delta if s is not None else None
        feasible = bool(
            congruent
            and nonnegative_rho
            and opposite_edges is not None
            and opposite_edges >= 0
        )
        parameter_rows.append(
            {
                "P": parallel_count,
                "rho": rho,
                "s": s,
                "opposite_edges": opposite_edges,
                "feasible": feasible,
            }
        )
    feasible_rows = [row for row in parameter_rows if row["feasible"]]
    hard_edges = M * expected_p + t + 1
    total_edges = H_EDGE_COUNT + 2 * t
    opposite_edges = total_edges - hard_edges
    previous_mean = (
        (P - 1) * expected_p
        + (P + 1) * (minimum_q - 1)
        + mean_constant
        - 7 * P
    )
    minimum_mean = previous_mean + P + 1
    surplus = opposite_edges - M * minimum_q
    minimum_count = M - surplus
    catalog_offsets = {"omitted_pair": 2, "all_equal_triple": 4}
    catalog_survivors = [
        name
        for name, catalog_offset in catalog_offsets.items()
        if (minimum_q - catalog_offset) % Q == 0
    ]
    expected_survivors = [] if surviving_family is None else [surviving_family]
    proved = bool(
        p11_phase_one_residue_replay(t)["proved"]
        and [(row["P"], row["rho"], row["s"]) for row in feasible_rows]
        == [(expected_p, 0, expected_p)]
        and hard_edges + opposite_edges == total_edges
        and previous_mean == -4
        and minimum_mean == 8
        and surplus == t + 2
        and minimum_count == 4 - t
        and catalog_survivors == expected_survivors
    )
    if not proved:
        raise ArithmeticError("the p=11 branch parallel replay failed")
    return {
        "p": P,
        "layer_index_t": t,
        "branch": branch,
        "coefficient_offset": offset,
        "parameter_rows_P_0_through_8": parameter_rows,
        "forced_P": expected_p,
        "forced_rho": 0,
        "forced_s": expected_p,
        "hard_finite_edge_count": hard_edges,
        "opposite_parallel_count_sum": opposite_edges,
        "minimum_opposite_Q": minimum_q,
        "mean_at_Q_minus_1": previous_mean,
        "mean_at_minimum_Q": minimum_mean,
        "parallel_surplus_above_minimum": surplus,
        "directions_at_minimum_at_least": minimum_count,
        "sharp_catalog_offsets_mod_5": catalog_offsets,
        "catalog_forms_with_offset_congruent_to_Q": catalog_survivors,
        "branch_excluded_by_catalog_offsets": not catalog_survivors,
        "all_equal_triple_is_only_catalog_survivor": (
            catalog_survivors == ["all_equal_triple"]
        ),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p11_sharp_equality_dependency() -> dict[str, object]:
    """Import the 15.688 Boolean rigidity used by the 15.736 catalog."""
    lift = sharp_integral_quadratic_lift_floor(P)
    phase_zero_floors = {
        b: scaled_direction_floor(P, b, 0) for b in range(0, P, 2)
    }
    least_nonzero_b_floor = min(
        floor for b, floor in phase_zero_floors.items() if b != 0
    )
    catalog = proposition_15736()
    proved = bool(
        lift["proved"]
        and lift["sharp_scaled_floor"] == 8
        and "B is Boolean" in str(lift["equality_rigidity"])
        and least_nonzero_b_floor == 12
        and catalog["proved"]
        and catalog["sharp_boolean_catalog_certified"]
    )
    if not proved:
        raise ArithmeticError("the p=11 sharp-equality dependency failed")
    return {
        "p": P,
        "phase_zero_even_floors": phase_zero_floors,
        "least_nonzero_b_floor": least_nonzero_b_floor,
        "sharp_b0_integral_lift_floor": lift["sharp_scaled_floor"],
        "equality_rigidity": lift["equality_rigidity"],
        "equality_forces_B_boolean": True,
        "boolean_catalog_dependency": {
            "proposition": "15.736",
            "result_status": catalog["result_status"],
            "certified": catalog["sharp_boolean_catalog_certified"],
        },
        "proved": proved,
    }


def _normalized_pair_coefficients(
    *,
    p: int,
    parallel_count: int,
    target_constant: int,
    target_linear: Sequence[int],
    target_pairs: dict[tuple[int, int], int],
) -> dict[str, object]:
    """Compare a target with ``eps*S=P+eps*sum K_st z_s z_t``.

    If their difference vanishes on the middle slice, its multilinear form
    is ``(sum z_s-1)(c+sum a_s z_s)``.  Constant and linear comparison gives

        a_s=c+target_linear_s,
        (p-1)c=P-target_constant-sum target_linear.

    Pair comparison then gives

        eps*K_st=target_pair_st+2c+target_linear_s+target_linear_t.

    The two uses below have ``c=0``.  Returning ``eps*K`` makes the sign
    convention independent of which Paley type is called hard.
    """
    if len(target_linear) != p:
        raise ValueError("target_linear must have one entry per fibre")
    if any(not (0 <= s < t < p) for s, t in target_pairs):
        raise ValueError("target pair lies outside the fibre set")
    numerator = parallel_count - target_constant - sum(target_linear)
    if numerator % (p - 1):
        raise ArithmeticError("slice-kernel scalar is not integral")
    kernel_scalar = numerator // (p - 1)
    a = [kernel_scalar + int(value) for value in target_linear]
    normalized = {
        (s, t): int(target_pairs.get((s, t), 0)) + a[s] + a[t]
        for s, t in combinations(range(p), 2)
    }
    constant_check = sum(a) - kernel_scalar
    return {
        "p": p,
        "parallel_count": parallel_count,
        "target_constant": target_constant,
        "target_linear_sum": sum(target_linear),
        "kernel_scalar_numerator": numerator,
        "kernel_scalar": kernel_scalar,
        "kernel_linear_coefficients": a,
        "constant_coefficient_reconstructed": constant_check,
        "constant_coefficient_expected": parallel_count - target_constant,
        "normalized_pair_coefficients_eps_times_K": normalized,
        "proved": constant_check == parallel_count - target_constant,
    }


@lru_cache(maxsize=1)
def exact_coefficient_patterns() -> dict[str, object]:
    """Derive the hard-star and opposite-triangle signed cell patterns."""
    hard_center = 0
    hard_linear = [-1 if s == hard_center else 0 for s in range(P)]
    hard = _normalized_pair_coefficients(
        p=P,
        parallel_count=3,
        target_constant=4,
        target_linear=hard_linear,
        target_pairs={},
    )
    hard_pattern = hard["normalized_pair_coefficients_eps_times_K"]
    assert isinstance(hard_pattern, dict)
    hard_support = [pair for pair, value in hard_pattern.items() if value]

    triple = (0, 1, 2)
    triple_pairs = {pair: 1 for pair in combinations(triple, 2)}
    opposite = _normalized_pair_coefficients(
        p=P,
        parallel_count=4,
        target_constant=4,
        target_linear=[0] * P,
        target_pairs=triple_pairs,
    )
    opposite_pattern = opposite["normalized_pair_coefficients_eps_times_K"]
    assert isinstance(opposite_pattern, dict)
    opposite_support = [pair for pair, value in opposite_pattern.items() if value]

    hard_ok = bool(
        hard["proved"]
        and hard["kernel_scalar"] == 0
        and len(hard_support) == P - 1
        and all(hard_center in pair for pair in hard_support)
        and {hard_pattern[pair] for pair in hard_support} == {-1}
    )
    opposite_ok = bool(
        opposite["proved"]
        and opposite["kernel_scalar"] == 0
        and opposite_support == list(combinations(triple, 2))
        and {opposite_pattern[pair] for pair in opposite_support} == {1}
    )
    if not hard_ok or not opposite_ok:
        raise ArithmeticError("the exact p=11 signed cell patterns changed")
    return {
        "sign_convention": (
            "K^L_st=sum chi(u-v) over selected edges between fibres; "
            "eps_L*S_H=P+eps_L*sum K^L_st*z_s*z_t"
        ),
        "hard_low": {
            "target": "eps_L*S_H=4-z_j",
            "parallel_count": 3,
            "normalized_pattern": "eps_L*K_st=-1 on the star at j, else 0",
            "support_size": len(hard_support),
            "kernel_scalar": hard["kernel_scalar"],
            "proved": hard_ok,
        },
        "opposite_minimum": {
            "target": "eps_L*S_H=4+z_i*z_j+z_i*z_k+z_j*z_k",
            "parallel_count": 4,
            "normalized_pattern": (
                "eps_L*K_st=1 on the triangle {i,j,k}, else 0"
            ),
            "support_size": len(opposite_support),
            "kernel_scalar": opposite["kernel_scalar"],
            "proved": opposite_ok,
        },
        "proved": hard_ok and opposite_ok,
    }


def star_square_moment(p: int, center: int) -> int:
    """Return ``sum_(t != center) (center-t)^2`` modulo ``p``."""
    if p < 3 or p % 2 == 0 or not 0 <= center < p:
        raise ValueError("need an odd p and a fibre center in 0..p-1")
    return sum((center - t) ** 2 for t in range(p) if t != center) % p


@lru_cache(maxsize=1)
def hard_star_moment_certificate() -> dict[str, object]:
    """Check that every hard-star normalized second moment vanishes."""
    values = [star_square_moment(P, center) for center in range(P)]
    nonzero_square_sum = sum(t * t for t in range(1, P))
    proved = all(value == 0 for value in values) and nonzero_square_sum % P == 0
    if not proved:
        raise ArithmeticError("the hard-star second moment stopped vanishing")
    return {
        "field": "F_11",
        "raw_nonzero_square_sum": nonzero_square_sum,
        "nonzero_square_sum_mod_11": nonzero_square_sum % P,
        "translation_identity": (
            "sum_(t!=j)(j-t)^2=sum_(u in F_11^*)u^2"
        ),
        "center_values": values,
        "every_hard_star_moment_is_zero": True,
        "proved": proved,
    }


def projective_quadratic_zero_count(
    p: int, coefficient_r2: int, coefficient_rs: int, coefficient_s2: int
) -> int:
    """Count zeros of a binary quadratic on ``P^1(F_p)`` exactly."""
    coefficients = (
        coefficient_r2 % p,
        coefficient_rs % p,
        coefficient_s2 % p,
    )
    if coefficients == (0, 0, 0):
        return p + 1
    representatives = [(1, t) for t in range(p)] + [(0, 1)]
    a, b, c = coefficients
    return sum((a * r * r + b * r * s + c * s * s) % p == 0 for r, s in representatives)


@lru_cache(maxsize=1)
def binary_quadratic_projective_root_certificate() -> dict[str, object]:
    """Audit the elementary at-most-two projective-root bound over F_11."""
    histogram: dict[int, int] = {}
    maximum = 0
    for a, b, c in product(range(P), repeat=3):
        if (a, b, c) == (0, 0, 0):
            continue
        count = projective_quadratic_zero_count(P, a, b, c)
        histogram[count] = histogram.get(count, 0) + 1
        maximum = max(maximum, count)
    proved = maximum == 2 and sum(histogram.values()) == P**3 - 1
    if not proved:
        raise ArithmeticError("binary quadratic projective root bound changed")
    return {
        "field": "F_11",
        "nonzero_binary_quadratic_count": P**3 - 1,
        "projective_point_count": P + 1,
        "zero_count_histogram": dict(sorted(histogram.items())),
        "maximum_projective_zero_count": maximum,
        "symbolic_reason": (
            "dehomogenize away from infinity; if infinity is a root, the "
            "remaining polynomial has degree at most one"
        ),
        "proved": proved,
    }


def triangle_square_moment(p: int, triple: Sequence[int]) -> int:
    """Return the triangle second moment modulo ``p``."""
    if len(triple) != 3 or len(set(triple)) != 3:
        raise ValueError("need three distinct fibre labels")
    if any(not isinstance(value, int) or not 0 <= value < p for value in triple):
        raise ValueError("triple labels must lie in 0..p-1")
    return sum((a - b) ** 2 for a, b in combinations(triple, 2)) % p


@lru_cache(maxsize=1)
def all_equal_triangle_nondegeneracy() -> dict[str, object]:
    """Prove that no three distinct F_11 labels have zero square moment."""
    histogram: dict[int, int] = {}
    zero_triples: list[tuple[int, int, int]] = []
    for triple in combinations(range(P), 3):
        value = triangle_square_moment(P, triple)
        histogram[value] = histogram.get(value, 0) + 1
        if value == 0:
            zero_triples.append(triple)
    discriminant = (-3) % P
    nonzero_squares = {value * value % P for value in range(1, P)}
    proved = bool(
        not zero_triples
        and discriminant == 8
        and discriminant not in nonzero_squares
        and sum(histogram.values()) == 165
        and set(histogram) == set(range(1, P))
    )
    if not proved:
        raise ArithmeticError("an all-equal triangle acquired zero moment")
    return {
        "field": "F_11",
        "distinct_triple_count": 165,
        "normal_form": "(i,j,k)=(0,1,r)",
        "normalized_moment": "2*(r^2-r+1)",
        "quadratic_discriminant_mod_11": discriminant,
        "nonzero_quadratic_residues_mod_11": sorted(nonzero_squares),
        "discriminant_is_nonsquare": True,
        "moment_value_histogram": dict(sorted(histogram.items())),
        "zero_moment_triples": zero_triples,
        "every_all_equal_triangle_has_nonzero_moment": True,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def p11_binary_moment_exclusion(t: int) -> dict[str, object]:
    """Exclude both sharp branches in one of the first three p=11 layers."""
    _check_layer_index(t)
    chart = p11_isolated_layer_chart(t)
    residues = p11_phase_one_residue_replay(t)
    equality = p11_sharp_equality_dependency()
    branch_a = p11_branch_parallel_replay(t, BRANCH_B2)
    branch_c = p11_branch_parallel_replay(t, BRANCH_P3_LAST)
    consequence = residual_p11_consequence() if t == 0 else None
    patterns = exact_coefficient_patterns()
    stars = hard_star_moment_certificate()
    roots = binary_quadratic_projective_root_certificate()
    triangles = all_equal_triangle_nondegeneracy()

    hard_low_direction_count = int(residues["endpoint_low_direction_count_at_least"])
    minimum_triple_direction_count = int(branch_c["directions_at_minimum_at_least"])
    root_bound = int(roots["maximum_projective_zero_count"])
    global_moment_forced_zero = hard_low_direction_count > root_bound
    triple_contradiction = bool(
        global_moment_forced_zero
        and minimum_triple_direction_count >= 1
        and triangles["every_all_equal_triangle_has_nonzero_moment"]
    )
    proved = bool(
        chart["proved"]
        and residues["proved"]
        and equality["proved"]
        and branch_a["proved"]
        and branch_a["branch_excluded_by_catalog_offsets"]
        and branch_c["proved"]
        and branch_c["all_equal_triple_is_only_catalog_survivor"]
        and (
            t != 0
            or (
                consequence is not None
                and consequence["proved_reduction"]
                and consequence["hard_b2_branch"]["excluded"]
                and consequence["hard_b_p_minus_1_branch"][
                    "all_equal_triple_survives"
                ]
            )
        )
        and patterns["proved"]
        and stars["proved"]
        and roots["proved"]
        and triangles["proved"]
        and triple_contradiction
    )
    if not proved:
        raise ArithmeticError("the p=11 binary-moment exclusion failed")
    return {
        "p": P,
        "layer_index_t": t,
        "original_k": CRITICAL_K + 2 * t,
        "H_edge_count": H_EDGE_COUNT + 2 * t,
        "isolated_chart": chart,
        "phase_one_residue_replay": residues,
        "sharp_equality_dependency": equality,
        "hard_b2_parallel_replay": branch_a,
        "hard_b10_parallel_replay": branch_c,
        "global_moment": (
            "M_H(L)=sum_{e={u,v} in H}chi(u-v)*(L(u)-L(v))^2"
        ),
        "global_moment_degree": 2,
        "hard_low_direction_count": hard_low_direction_count,
        "hard_low_projective_zeros": hard_low_direction_count,
        "nonzero_binary_quadratic_projective_root_bound": root_bound,
        "global_moment_forced_identically_zero": global_moment_forced_zero,
        "minimum_all_equal_triple_direction_count": minimum_triple_direction_count,
        "one_triple_already_contradicts_zero_moment": True,
        "hard_b2_branch_excluded_by_15_736": True,
        "hard_b_p_minus_1_branch_excluded_by_binary_moment": triple_contradiction,
        "p11_layer_excluded": proved,
        "result_status": "proved theorem",
        "proved": proved,
    }


def proposition_15737() -> dict[str, object]:
    """Package the exact first-three-layer p=11 close."""
    dependency = proposition_15736()
    exclusions = {str(t): p11_binary_moment_exclusion(t) for t in CLOSED_LAYER_INDICES}
    proved = bool(
        dependency["proved"]
        and dependency["result_status"] == "exhaustive finite certificate"
        and all(row["proved"] for row in exclusions.values())
    )
    return {
        "prop": "15.737",
        "title": "Binary quadratic moment closes the first three p=11 layers",
        "result_status": "proved theorem",
        "statement": (
            "residual (ii) at p=11 is empty for k=44,46,48; at least three "
            "hard star baselines force the global binary quadratic moment "
            "to vanish, contradicting every all-equal triple target"
        ),
        "finite_certificate_dependency": {
            "proposition": "15.736",
            "result_status": dependency["result_status"],
            "sharp_boolean_catalog_certified": dependency[
                "sharp_boolean_catalog_certified"
            ],
        },
        "coefficient_patterns": exact_coefficient_patterns(),
        "hard_star_second_moment": hard_star_moment_certificate(),
        "binary_quadratic_projective_roots": (
            binary_quadratic_projective_root_certificate()
        ),
        "all_equal_triangle_nondegeneracy": all_equal_triangle_nondegeneracy(),
        "p11_layer_exclusions": exclusions,
        "closed_layer_indices_t": list(CLOSED_LAYER_INDICES),
        "closed_even_k": [44, 46, 48],
        "critical_p11_k_eq_44_closed": proved,
        "p11_k_eq_46_closed": proved,
        "p11_k_eq_48_closed": proved,
        "critical_p5_closed": False,
        "critical_p7_closed": False,
        "p11_k_at_least_50_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "critical p=5,7; p=11 even k>=50; the later p>=13 shells not "
            "covered by 15.734-15.735; multi-level Type I; and the limit"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic p=11 binary-moment certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15737.json"
    payload = json.dumps(proposition_15737(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15737()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.737 audit failed")
    path = write_evidence()
    print("Prop. 15.737: p=11 residual (ii) excluded at k=44,46,48")
    print("p=5,7, p=11 k>=50, and full residual (ii) remain open")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
