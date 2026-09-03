#!/usr/bin/env python3
"""Exact certificate for the live irreducible-conic odd-Radon branch.

The proof is in evidence/NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md.
This module records its exact constants and independently replays the atom
witness at p=31.  It does not assert simultaneous vanishing at degrees six
and eight, a common finite-field lift, or residual-(ii) closure.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations

from e1_gmin_m4_compact_ray_moment_gate import (
    all_equal_moment,
    compact_moment,
)


P31_AE_ATOMS = (
    (1, 3, 13),
    (2, 7, 8),
    (4, 18, 26),
    (5, 20, 23),
    (6, 14, 28),
    (24, 25, 30),
)

# Entries are (sorted triple, distinguished label).  The compact boundary is
# +ab-ac-bc when c is distinguished.
P31_COMPACT_ATOMS = (
    ((0, 2, 12), 2),
    ((0, 12, 19), 12),
    ((4, 19, 22), 4),
    ((9, 12, 16), 16),
    ((10, 14, 17), 14),
    ((11, 19, 22), 19),
    ((12, 16, 20), 20),
)


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def _check_p(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not _is_prime(p)
    ):
        raise ValueError("need a prime p=4r+3 with p>=31")
    return (p - 3) // 4


def conic_reduction_constants(p: int, b: int) -> dict[str, object]:
    """Return the exact inequalities in the peeling/classification proof."""
    r = _check_p(p)
    if not isinstance(b, int) or isinstance(b, bool) or not 0 <= b <= r:
        raise ValueError("need 0<=b<=r")
    h = 2 * r + 1
    m = h - 2
    occurrence_budget = 3 * (r + b - 1)

    # Exact squared forms of the two strict square-root inequalities.
    no_constant_rhs = 3 * p - 31
    one_constant_rhs = 2 * p - 26
    no_constant_strict = (
        no_constant_rhs > 0
        and no_constant_rhs * no_constant_rhs > 36 * p
    )
    one_constant_strict = (
        one_constant_rhs > 0
        and one_constant_rhs * one_constant_rhs > 16 * p
    )

    outside_upper_bound = occurrence_budget - (2 * m + 2)
    return {
        "p": p,
        "r": r,
        "h": h,
        "m": m,
        "b": b,
        "occurrence_budget": occurrence_budget,
        "support_upper_bound": occurrence_budget,
        "support_at_most_3m": occurrence_budget <= 3 * m,
        "conic_points_required": 2 * m + 2,
        "outside_support_upper_bound": outside_upper_bound,
        "outside_nonzero_dual_minimum_support": m,
        "peeling_contradiction": outside_upper_bound <= m - 2,
        "no_constant_character_bound_strict": no_constant_strict,
        "one_constant_character_bound_strict": one_constant_strict,
        "forced_normal_form": "U=u*z^2, D=d*(z-1)^2 with u,d squares",
        "normal_form_Omega_points": p - 2,
        "dual_weight": "W(z)=c*z*(z-1)*(A*z+B)",
        "orbit_difference": "n(z)=alpha+beta/(z-1)",
        "beta_nonzero_least_l1": h * (h - 1),
        "beta_nonzero_excluded_by_l1": h * (h - 1) > occurrence_budget,
        "constant_branch_requires": f"3*b >= {r + 4}",
        "constant_branch_numerically_possible": 3 * b >= r + 4,
        "proved": (
            occurrence_budget <= 3 * m
            and outside_upper_bound <= m - 2
            and no_constant_strict
            and one_constant_strict
            and h * (h - 1) > occurrence_budget
        ),
    }


def _edge(p: int, a: int, b: int) -> tuple[int, int]:
    a %= p
    b %= p
    if a == b:
        raise ValueError("loop is not an atom edge")
    return tuple(sorted((a, b)))


def _negate_edge(p: int, edge: tuple[int, int]) -> tuple[int, int]:
    return _edge(p, -edge[0], -edge[1])


def _orbit_coordinate(
    p: int, edge: tuple[int, int]
) -> tuple[tuple[int, int] | None, int]:
    negative = _negate_edge(p, edge)
    if edge == negative:
        return None, 0
    representative = min(edge, negative)
    return representative, 1 if edge == representative else -1


def tangent_conic_target(p: int, k: int) -> dict[tuple[int, int], int]:
    """Return the constant-n tangent-conic edge-orbit word."""
    _check_p(p)
    k %= p
    if not k:
        raise ValueError("k must be nonzero")
    inverse_two = pow(2, -1, p)
    target: dict[tuple[int, int], int] = {}
    for z in range(p):
        if z in (0, 1):
            continue
        a = (((1 + k) * z - k) * inverse_two) % p
        b = (((1 - k) * z + k) * inverse_two) % p
        edge = _edge(p, a, b)
        representative, sign = _orbit_coordinate(p, edge)
        if representative is None or representative in target:
            raise AssertionError("tangent-conic parameterization is not injective")
        target[representative] = sign
    if len(target) != p - 2:
        raise AssertionError("wrong tangent-conic support")
    return target


def _add_orbit_edge(
    chain: dict[tuple[int, int], int],
    p: int,
    a: int,
    b: int,
    coefficient: int,
) -> None:
    representative, sign = _orbit_coordinate(p, _edge(p, a, b))
    if representative is not None:
        chain[representative] += coefficient * sign


def _p31_witness_imbalance() -> dict[tuple[int, int], int]:
    p = 31
    chain: dict[tuple[int, int], int] = defaultdict(int)
    for triangle in P31_AE_ATOMS:
        for a, b in combinations(triangle, 2):
            _add_orbit_edge(chain, p, a, b, 1)
    for triple, distinguished in P31_COMPACT_ATOMS:
        a, b = (value for value in triple if value != distinguished)
        _add_orbit_edge(chain, p, a, b, 1)
        _add_orbit_edge(chain, p, a, distinguished, -1)
        _add_orbit_edge(chain, p, b, distinguished, -1)
    return {edge: value for edge, value in chain.items() if value}


def _p31_moment_vector(degree: int) -> tuple[int, ...]:
    p = 31
    values = []
    for channel in range(degree // 2):
        value = sum(
            all_equal_moment(p, *triangle, degree, channel)
            for triangle in P31_AE_ATOMS
        )
        for triple, distinguished in P31_COMPACT_ATOMS:
            a, b = (entry for entry in triple if entry != distinguished)
            value += compact_moment(
                p, a, b, distinguished, degree, channel
            )
        values.append(value % p)
    return tuple(values)


def _multiplicative_order(value: int, p: int) -> int:
    product = 1
    for order in range(1, p):
        product = product * value % p
        if product == 1:
            return order
    raise AssertionError("nonzero field element has no multiplicative order")


def _q_cycle_triangles(p: int, q: int) -> tuple[tuple[int, int, int], ...]:
    centre = pow(2, -1, p)
    unused = set(range(p)) - {centre}
    cycles = []
    while unused:
        first = min(unused)
        cycle = []
        value = first
        while value not in cycle:
            cycle.append(value)
            value = (centre + q * (value - centre)) % p
        if value != first or len(cycle) != 3:
            raise AssertionError("expected order-three affine cycles")
        cycles.append(tuple(sorted(cycle)))
        unused.difference_update(cycle)
    return tuple(sorted(cycles))


def _occurrence_alignment(
    p: int,
    target: dict[tuple[int, int], int],
    a: int,
    b: int,
    coefficient: int,
) -> int:
    representative, sign = _orbit_coordinate(p, _edge(p, a, b))
    if representative is None or representative not in target:
        return 0
    return coefficient * sign * target[representative]


def _compact_alignment_score(
    p: int,
    target: dict[tuple[int, int], int],
    triple: tuple[int, int, int],
    distinguished: int,
) -> int:
    a, b = (value for value in triple if value != distinguished)
    return (
        _occurrence_alignment(p, target, a, b, 1)
        + _occurrence_alignment(p, target, a, distinguished, -1)
        + _occurrence_alignment(p, target, b, distinguished, -1)
    )


def nonequianharmonic_score_three_compact_candidates(
    p: int, k: int
) -> tuple[tuple[tuple[int, int, int], int], ...]:
    """Return the at-most-two score-three compact atoms from formula (11b)."""
    _check_p(p)
    k %= p
    if k in (0, 1, p - 1):
        raise ValueError("need a nonstar tangent conic with k nonzero")
    q = (1 - k) * pow(1 + k, -1, p) % p
    if pow(q, 3, p) == 1:
        raise ValueError("need the nonequianharmonic branch q^3 != 1")
    denominator = (q * q + q + 1) % p
    inverse_q = pow(q, -1, p)
    inverse_denominator = pow(denominator, -1, p)
    inverse_two = pow(2, -1, p)

    # Formula (11b) is written in X=2x coordinates.
    raw_x_atoms = (
        (
            (-(q + 2) * inverse_q) % p,
            (-(2 * q + 1)) % p,
            (-3) % p,
        ),
        (
            ((q * q - q - 1) * inverse_denominator) % p,
            (-(q * q + q - 1) * inverse_denominator) % p,
            (-(q * q - q + 1) * inverse_denominator) % p,
        ),
    )
    target = tangent_conic_target(p, k)
    candidates = []
    for x_a, x_b, x_c in raw_x_atoms:
        a, b, distinguished = (
            x_a * inverse_two % p,
            x_b * inverse_two % p,
            x_c * inverse_two % p,
        )
        triple = tuple(sorted((a, b, distinguished)))
        if len(set(triple)) < 3:
            continue
        candidate = (triple, distinguished)
        if (
            _compact_alignment_score(
                p, target, triple, distinguished
            )
            == 3
            and candidate not in candidates
        ):
            candidates.append(candidate)
    return tuple(candidates)


def nonequianharmonic_constant_fiber_no_go(
    p: int, b: int, k: int
) -> dict[str, object]:
    """Certify that a nonstar constant atom fiber forces q^3=1."""
    r = _check_p(p)
    if not isinstance(b, int) or isinstance(b, bool) or not 0 <= b <= r:
        raise ValueError("need 0<=b<=r")
    k %= p
    if k in (0, 1, p - 1):
        raise ValueError("need a nonstar tangent conic with k nonzero")
    q = (1 - k) * pow(1 + k, -1, p) % p
    if pow(q, 3, p) == 1:
        raise ValueError("equianharmonic fibers are the surviving branch")
    candidates = nonequianharmonic_score_three_compact_candidates(p, k)
    total_atoms = r - 1 + b
    score_upper_bound = 2 * total_atoms + 2
    target_score = p - 2
    return {
        "p": p,
        "r": r,
        "b": b,
        "k": k,
        "q": q,
        "q_cubed": pow(q, 3, p),
        "all_equal_atom_score_upper_bound": 2,
        "ordinary_compact_atom_score_upper_bound": 2,
        "score_three_compact_candidate_count": len(candidates),
        "score_three_compact_candidates": [
            {"triple": list(triple), "distinguished": distinguished}
            for triple, distinguished in candidates
        ],
        "score_three_compact_candidate_bound": 2,
        "total_atom_count": total_atoms,
        "total_atom_count_upper_bound": 2 * r - 1,
        "defect_lower_bound": total_atoms - 2,
        "score_upper_bound": score_upper_bound,
        "target_score": target_score,
        "strict_score_contradiction": score_upper_bound < target_score,
        "conclusion": "no exact constant tangent-conic atom fiber",
        "proved": len(candidates) <= 2 and score_upper_bound < target_score,
    }


def p31_equianharmonic_witness_certificate() -> dict[str, object]:
    """Replay the exact noncentral odd-zero row at p=31,b=7."""
    p = 31
    k = 11
    q = (1 - k) * pow(1 + k, -1, p) % p
    target = tangent_conic_target(p, k)
    imbalance = _p31_witness_imbalance()
    odd_degrees = tuple(range(3, p - 1, 2))
    odd_vectors = {degree: _p31_moment_vector(degree) for degree in odd_degrees}
    degree_six = _p31_moment_vector(6)
    degree_eight = _p31_moment_vector(8)
    cycles = _q_cycle_triangles(p, q)
    all_odd_zero = all(not any(vector) for vector in odd_vectors.values())
    edge_replay_exact = imbalance == target
    ae_are_cycles = all(tuple(atom) in cycles for atom in P31_AE_ATOMS)
    return {
        "p": p,
        "r": 7,
        "b": 7,
        "k": k,
        "k_squared_mod_p": k * k % p,
        "minus_three_mod_p": (-3) % p,
        "q": q,
        "q_order": _multiplicative_order(q, p),
        "q_satisfies_equianharmonic_polynomial": (q * q + q + 1) % p == 0,
        "ae_atoms": [list(atom) for atom in P31_AE_ATOMS],
        "compact_atoms": [
            {"triple": list(triple), "distinguished": distinguished}
            for triple, distinguished in P31_COMPACT_ATOMS
        ],
        "ae_atoms_are_q_cycle_triangles": ae_are_cycles,
        "target_support": len(target),
        "target_l1": sum(abs(value) for value in target.values()),
        "edge_orbit_replay_exact": edge_replay_exact,
        "odd_degrees_checked": list(odd_degrees),
        "odd_channel_count": sum(degree // 2 for degree in odd_degrees),
        "all_odd_channels_zero": all_odd_zero,
        "degree_six": list(degree_six),
        "degree_eight": list(degree_eight),
        "degree_six_and_eight_both_zero": (
            not any(degree_six) and not any(degree_eight)
        ),
        "central_signed_chain": not imbalance,
        "proved": (
            k * k % p == (-3) % p
            and _multiplicative_order(q, p) == 3
            and ae_are_cycles
            and edge_replay_exact
            and all_odd_zero
            and degree_six == (11, 19, 10)
            and degree_eight == (12, 11, 23, 6)
        ),
    }


def scaled_family_exceptional_row_obstruction(p: int) -> dict[str, object]:
    """Record the unsigned projective obstruction for one scaled row family."""
    _check_p(p)
    identity_degree = 24
    return {
        "p": p,
        "projective_direction_count": p + 1,
        "identity_degree": identity_degree,
        "point_count_forces_A4_equals_B3": p + 1 > identity_degree,
        "ufd_factorization": "A=c*Q^3, B=d*Q^4 for a binary quadratic Q",
        "constant_nonzero_quadratic_character_on_P1_possible": False,
        "conclusion": (
            "under unsigned common-form interpolation, one fixed nonzero "
            "degree-6/8 row cannot be scaled into every projective direction"
        ),
        "actual_signed_Paley_exceptional_row_proved": False,
        "paley_half_exception": (
            "on one Paley type, Q can be a scalar multiple of the anisotropic "
            "norm form and the obstruction becomes a compatible half-system"
        ),
        "proved": p + 1 > identity_degree,
    }


def theorem_record() -> dict[str, object]:
    witness = p31_equianharmonic_witness_certificate()
    return {
        "title": "Irreducible-conic odd-Radon dichotomy",
        "status": "PROVED REDUCTION WITH EXACT LIVE COUNTEREXAMPLE",
        "proved": {
            "conic_containing_low_weight_word_is_fully_conic_supported": True,
            "high_intersection_irreducible_conic_is_triangle_tangent": True,
            "nonconstant_affine_dual_weight_is_excluded_by_integer_l1": True,
            "only_constant_plus_or_minus_one_conic_word_survives": True,
            "star_constant_branch_is_excluded_by_quotient_parity": True,
            "nonequianharmonic_constant_branch_is_excluded": True,
            "nonstar_constant_branch_is_excluded": False,
            "constant_branch_forces_q_cubed_equals_one": True,
            "constant_branch_forces_p_congruent_7_mod_12": True,
            "p31_b7_equianharmonic_odd_zero_atom_witness_exists": witness["proved"],
            "p31_witness_degree_six_and_eight_both_zero": False,
            "unsigned_one_scaled_nonzero_row_family_works_on_every_direction": False,
            "actual_signed_Paley_exceptional_row_proved": False,
            "paley_half_norm_form_coordination_is_algebraically_compatible": True,
            "common_Fp_atom_lift_constructed": False,
            "Boolean_lift_constructed": False,
            "residual_ii_closed": False,
        },
        "p31_witness": witness,
        "remaining_gate": (
            "classify equianharmonic constant tangent-conic atom fibers and match "
            "the exceptional opposite rows and all hard rows to common "
            "degree-six/eight forms, then solve the integral/Boolean lift"
        ),
        "duplicate_work_guard": (
            "Do not claim all low-weight odd-zero chains are central: the "
            "p31,b=7,k=11 witness is an exact counterexample."
        ),
        "L_status": "OPEN",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), sort_keys=True, indent=2))
