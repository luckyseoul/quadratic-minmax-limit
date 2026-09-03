#!/usr/bin/env python3
"""Adaptive center coherence for complementary Mobius halves.

This module proves one bounded intermediate theorem at the balanced branch-C
endpoint.  After an opposite fixed-edge direction is chosen, arbitrary
nonzero hard-star centers admit a magnitude of that fixed edge, singleton
signs, and a perfect matching of the hard target directions for which every
matched pair has a nonsquare complementary parameter.

The theorem is deliberately *not* an endpoint construction.  Once a target
pair and its two singleton signs are fixed, its two auxiliary projective
directions are forced.  The argument below does not prove that the auxiliary
directions from different pairs are distinct, nor that their Paley types are
the quota-required ``m-2`` hard and two opposite directions.  Those coupled
block conditions are the remaining SDR/induced-matching gate.
"""

from __future__ import annotations

from collections.abc import Sequence

from e1_gmin_m4_prop15721 import is_prime


def _check_branch_prime(p: int) -> tuple[int, int]:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=4r+3 with p>=31")
    h = (p - 1) // 2
    return h, h + 1


def _inverse(p: int, value: int) -> int:
    value %= p
    if value == 0:
        raise ValueError("cannot invert zero")
    return pow(value, -1, p)


def quadratic_character(p: int, value: int) -> int:
    """Return the quadratic character in ``F_p`` as ``-1, 0, 1``."""
    value %= p
    if value == 0:
        return 0
    residue = pow(value, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def endpoint_color_profile(p: int, a: int) -> dict[str, object]:
    """Classify the singleton signs available at one normalized target.

    Here ``a=X(x)`` for ``X=L/j``.  A sign ``eps`` is valid precisely when
    ``a != 2*eps``.  Its color is

        chi(eps*(a-2*eps)).

    Two endpoints form a complementary pair exactly when their selected
    colors are opposite.
    """
    _check_branch_prime(p)
    a %= p
    if a == 0:
        raise ValueError("the fixed direction must differ from every target")
    colors: dict[int, int] = {}
    for eps in (1, -1):
        denominator = (a - 2 * eps) % p
        if denominator:
            colors[eps] = quadratic_character(p, eps * denominator)
    if not colors:
        raise ArithmeticError("both singleton signs cannot be singular")
    distinct = set(colors.values())
    classification = "flexible" if len(distinct) == 2 else "monochrome"
    return {
        "a": a,
        "sign_colors": colors,
        "valid_signs": tuple(colors),
        "classification": classification,
        "monochrome_color": next(iter(distinct)) if len(distinct) == 1 else None,
        "a_squared_minus_four_character": quadratic_character(p, a * a - 4),
        "proved": True,
    }


def complementary_pair_parameters(
    p: int,
    a: int,
    b: int,
    eps: int,
    eps_prime: int,
) -> dict[str, object]:
    """Solve the exact scaled complementary-pair singleton equations.

    Normalize the prescribed targets by ``X=L/j`` and ``Y=L'/j'`` and put
    ``rho=j'/j``.  If ``N=M/j`` and ``N'=M'/j'`` are the normalized
    auxiliaries, complementary Mobius halves obey

        rho*Y = nu*X + (1-nu)*N,
        N' = (nu/rho)*(X-N).

    The inputs are ``a=X(x)``, ``b=Y(x)`` and the requested singleton
    evaluations ``N(x)=2*eps``, ``N'(x)=2*eps_prime``.
    """
    _check_branch_prime(p)
    if eps not in (-1, 1) or eps_prime not in (-1, 1):
        raise ValueError("singleton signs must be +1 or -1")
    a %= p
    b %= p
    if a == 0 or b == 0:
        raise ValueError("the fixed direction must differ from both targets")
    left_denominator = (a - 2 * eps) % p
    right_denominator = (b - 2 * eps_prime) % p
    if left_denominator == 0:
        raise ValueError("a=2*eps makes the first complementary scale singular")
    if right_denominator == 0:
        raise ValueError("b=2*eps_prime makes rho singular")

    rho = 2 * eps * _inverse(p, right_denominator) % p
    nu = (
        4
        * eps
        * eps_prime
        * _inverse(p, left_denominator * right_denominator)
    ) % p
    if nu == 1:
        raise ValueError("nu=1 collapses the complementary pair")
    one_minus_nu_inverse = _inverse(p, 1 - nu)
    first_auxiliary_coefficients = (
        -nu * one_minus_nu_inverse % p,
        rho * one_minus_nu_inverse % p,
    )
    second_auxiliary_coefficients = (
        nu * _inverse(p, rho) * one_minus_nu_inverse % p,
        -nu * one_minus_nu_inverse % p,
    )

    first_evaluation = (
        first_auxiliary_coefficients[0] * a
        + first_auxiliary_coefficients[1] * b
    ) % p
    second_evaluation = (
        second_auxiliary_coefficients[0] * a
        + second_auxiliary_coefficients[1] * b
    ) % p
    first_color = quadratic_character(p, eps * left_denominator)
    second_color = quadratic_character(p, eps_prime * right_denominator)
    proved = bool(
        rho
        and nu
        and first_evaluation == 2 * eps % p
        and second_evaluation == 2 * eps_prime % p
        and quadratic_character(p, nu) == first_color * second_color
    )
    if not proved:
        raise ArithmeticError("the complementary singleton formulas changed")
    return {
        "a": a,
        "b": b,
        "eps": eps,
        "eps_prime": eps_prime,
        "rho_definition": "j_prime/j",
        "rho": rho,
        "nu": nu,
        "nu_character": quadratic_character(p, nu),
        "first_endpoint_color": first_color,
        "second_endpoint_color": second_color,
        "first_auxiliary_coefficients_in_X_Y": first_auxiliary_coefficients,
        "second_auxiliary_coefficients_in_X_Y": second_auxiliary_coefficients,
        "first_auxiliary_projective_direction": "[rho*Y-nu*X]",
        "second_auxiliary_projective_direction": "[X-rho*Y]",
        "first_auxiliary_evaluation": first_evaluation,
        "second_auxiliary_evaluation": second_evaluation,
        "complementary_nonsquare": quadratic_character(p, nu) == -1,
        "proved": True,
    }


def color_count_theorem(p: int) -> dict[str, object]:
    """Count the nonzero normalized evaluations of each endpoint kind."""
    h, m = _check_branch_prime(p)
    profiles = tuple(endpoint_color_profile(p, a) for a in range(1, p))
    flexible = sum(row["classification"] == "flexible" for row in profiles)
    mono_plus = sum(row["monochrome_color"] == 1 for row in profiles)
    mono_minus = sum(row["monochrome_color"] == -1 for row in profiles)
    eta_two = quadratic_character(p, 2)
    expected_plus = (h + eta_two) // 2
    expected_minus = (h + 2 - eta_two) // 2
    proved = bool(
        flexible == h - 1
        and mono_plus == expected_plus
        and mono_minus == expected_minus
        and mono_plus + mono_minus == m
        and flexible + m == p - 1
    )
    if not proved:
        raise ArithmeticError("the endpoint-color count changed")
    return {
        "p": p,
        "h": h,
        "m": m,
        "eta_2": eta_two,
        "flexible_nonzero_evaluations": flexible,
        "monochrome_plus_evaluations": mono_plus,
        "monochrome_minus_evaluations": mono_minus,
        "total_monochrome_nonzero_evaluations": m,
        "generic_flexible_criterion": "chi(a^2-4)=+1",
        "generic_monochrome_criterion": "chi(a^2-4)=-1",
        "exceptional_a=+2": "only eps=-1, monochrome color -1",
        "exceptional_a=-2": "only eps=+1, monochrome color -1",
        "proved": True,
    }


def _sign_for_color(profile: dict[str, object], color: int) -> int:
    sign_colors = dict(profile["sign_colors"])
    for sign, available_color in sign_colors.items():
        if available_color == color:
            return int(sign)
    raise ArithmeticError("a requested endpoint color is unavailable")


def adaptive_center_pairing_certificate(
    p: int, alpha_values: Sequence[int]
) -> dict[str, object]:
    """Construct the all-center target perfect matching.

    Fix a nonzero ``x0`` in the kernel of an opposite direction ``F`` and
    put ``alpha_i=(L_i/j_i)(x0)``.  All ``alpha_i`` are nonzero.  Replacing
    the fixed edge representative by ``x=c*x0`` changes the normalized
    evaluation to ``a_i=c*alpha_i``.

    Averaging over ``c`` finds a magnitude with at most ``m/2`` monochrome
    vertices.  The feasible-pair graph is complete except inside each of
    its two monochrome color classes, so the explicit construction below
    gives a perfect matching.
    """
    h, m = _check_branch_prime(p)
    if len(alpha_values) != m:
        raise ValueError("need one alpha value for every hard target")
    alphas = tuple(value % p for value in alpha_values)
    if any(value == 0 for value in alphas):
        raise ValueError("all hard centers and evaluations must be nonzero")

    scale_rows: list[dict[str, object]] = []
    for c in range(1, p):
        profiles = tuple(endpoint_color_profile(p, c * alpha) for alpha in alphas)
        mono_plus = tuple(
            index for index, row in enumerate(profiles)
            if row["monochrome_color"] == 1
        )
        mono_minus = tuple(
            index for index, row in enumerate(profiles)
            if row["monochrome_color"] == -1
        )
        scale_rows.append({
            "c": c,
            "profiles": profiles,
            "mono_plus": mono_plus,
            "mono_minus": mono_minus,
            "monochrome_total": len(mono_plus) + len(mono_minus),
        })

    total_monochrome_incidence = sum(
        int(row["monochrome_total"]) for row in scale_rows
    )
    expected_incidence = m * m
    candidates = [
        row for row in scale_rows if int(row["monochrome_total"]) <= m // 2
    ]
    if total_monochrome_incidence != expected_incidence or not candidates:
        raise ArithmeticError("the magnitude-averaging argument changed")
    witness = min(candidates, key=lambda row: (int(row["monochrome_total"]), int(row["c"])))
    c = int(witness["c"])
    profiles = tuple(witness["profiles"])
    plus = list(witness["mono_plus"])
    minus = list(witness["mono_minus"])
    flexible = [
        index for index, row in enumerate(profiles)
        if row["classification"] == "flexible"
    ]

    colored_pairs: list[tuple[int, int, int, int]] = []
    while plus and minus:
        colored_pairs.append((plus.pop(), 1, minus.pop(), -1))
    if plus:
        while plus:
            colored_pairs.append((plus.pop(), 1, flexible.pop(), -1))
    elif minus:
        while minus:
            colored_pairs.append((minus.pop(), -1, flexible.pop(), 1))
    if len(flexible) % 2:
        raise ArithmeticError("the unmatched flexible count must be even")
    while flexible:
        colored_pairs.append((flexible.pop(), 1, flexible.pop(), -1))

    labeled_pairs: list[dict[str, object]] = []
    used_targets: list[int] = []
    for first, first_color, second, second_color in colored_pairs:
        eps = _sign_for_color(profiles[first], first_color)
        eps_prime = _sign_for_color(profiles[second], second_color)
        a = c * alphas[first] % p
        b = c * alphas[second] % p
        parameters = complementary_pair_parameters(p, a, b, eps, eps_prime)
        if not parameters["complementary_nonsquare"]:
            raise ArithmeticError("an opposite-color edge was not complementary")
        used_targets.extend((first, second))
        labeled_pairs.append({
            "first_target": first,
            "second_target": second,
            "first_color": first_color,
            "second_color": second_color,
            "first_sign": eps,
            "second_sign": eps_prime,
            "parameters": parameters,
        })

    proved = bool(
        len(labeled_pairs) == m // 2
        and sorted(used_targets) == list(range(m))
        and int(witness["monochrome_total"]) <= m // 2
        and total_monochrome_incidence == m * m
        and m * m < (p - 1) * (m // 2 + 1)
    )
    if not proved:
        raise ArithmeticError("the adaptive target pairing certificate changed")
    return {
        "p": p,
        "h": h,
        "m": m,
        "alpha_values": alphas,
        "chosen_fixed_edge_scale_c": c,
        "chosen_scale_monochrome_total": witness["monochrome_total"],
        "chosen_scale_monochrome_plus": witness["mono_plus"],
        "chosen_scale_monochrome_minus": witness["mono_minus"],
        "sum_over_scales_of_monochrome_vertices": total_monochrome_incidence,
        "sum_over_scales_formula": "m^2",
        "strict_average_bound": "m^2/(p-1)<m/2+1",
        "matched_target_pairs": tuple(labeled_pairs),
        "target_perfect_matching_proved": True,
        "all_pair_parameters_nonsquare": True,
        "auxiliary_directions_evaluated": False,
        "auxiliary_direction_SDR_proved": False,
        "auxiliary_type_quota_proved": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def forced_affine_auxiliary_pair(
    p: int,
    z_i: int,
    z_k: int,
    alpha_i: int,
    alpha_k: int,
    c: int,
    eps_i: int,
    eps_k: int,
) -> dict[str, object]:
    """Evaluate the two forced auxiliary directions in the ``F``-affine chart.

    Every direction other than ``F`` is normalized to take value one at
    ``x0`` and has affine coordinate ``z``.  Put

        w_i=2*eps_i/alpha_i,  mu_i=w_i/(w_i-c).

    The two auxiliary coordinates are the unique two-cycle of the affine
    homotheties centered at ``z_i,z_k`` with multipliers ``mu_i,mu_k``.
    """
    _check_branch_prime(p)
    if eps_i not in (-1, 1) or eps_k not in (-1, 1):
        raise ValueError("singleton signs must be +1 or -1")
    z_i %= p
    z_k %= p
    alpha_i %= p
    alpha_k %= p
    c %= p
    if z_i == z_k:
        raise ValueError("target directions must be distinct")
    if alpha_i == 0 or alpha_k == 0 or c == 0:
        raise ValueError("alpha values and fixed-edge scale must be nonzero")
    w_i = 2 * eps_i * _inverse(p, alpha_i) % p
    w_k = 2 * eps_k * _inverse(p, alpha_k) % p
    if c == w_i:
        raise ValueError("a_i=2*eps_i")
    if c == w_k:
        raise ValueError("a_k=2*eps_k")
    denominator = (w_i + w_k - c) % p
    if denominator == 0:
        raise ValueError("nu=1 gives no finite two-cycle")
    inverse_denominator = _inverse(p, denominator)
    u = (
        w_k * z_i + (w_i - c) * z_k
    ) * inverse_denominator % p
    v = (
        (w_k - c) * z_i + w_i * z_k
    ) * inverse_denominator % p
    mu_i = w_i * _inverse(p, w_i - c) % p
    mu_k = w_k * _inverse(p, w_k - c) % p
    nu = mu_i * mu_k % p
    if quadratic_character(p, nu) != -1:
        raise ValueError("the selected signs do not give a complementary pair")

    h_i_u = (z_i + mu_i * (u - z_i)) % p
    h_k_v = (z_k + mu_k * (v - z_k)) % p
    separated = len({u, v, z_i, z_k}) == 4
    proved = bool(
        h_i_u == v
        and h_k_v == u
        and separated
        and (u - v) % p
        == c * (z_i - z_k) * inverse_denominator % p
    )
    if not proved:
        raise ArithmeticError("the forced affine auxiliary map changed")
    return {
        "z_i": z_i,
        "z_k": z_k,
        "alpha_i": alpha_i,
        "alpha_k": alpha_k,
        "c": c,
        "eps_i": eps_i,
        "eps_k": eps_k,
        "w_i": w_i,
        "w_k": w_k,
        "mu_i": mu_i,
        "mu_k": mu_k,
        "nu": nu,
        "nu_character": -1,
        "first_auxiliary_coordinate_U": u,
        "second_auxiliary_coordinate_V": v,
        "cycle_equations": "V=H_i(U), U=H_k(V)",
        "within_pair_auxiliaries_distinct": True,
        "auxiliaries_avoid_both_paired_targets": True,
        "cross_pair_auxiliary_distinctness_proved": False,
        "auxiliary_type_quota_proved": False,
        "proved": True,
    }


def prescribed_auxiliary_assignment_criterion(
    p: int,
    target_coordinates: Sequence[int],
    alpha_values: Sequence[int],
    auxiliary_coordinates: Sequence[int],
    target_partner: Sequence[int],
    sigma: Sequence[int],
    c: int,
    signs: Sequence[int],
) -> dict[str, object]:
    """Check the exact paired-SDR equations and their square invariant.

    ``target_partner`` is a fixed-point-free involution ``tau`` on target
    indices.  ``sigma[i]`` labels the auxiliary occurrence assigned to
    target ``i``.  For an actual SDR, ``sigma`` and the auxiliary coordinates
    must both be injective; accepting occurrences here makes collisions
    visible without hiding the exact algebra.
    """
    _check_branch_prime(p)
    size = len(target_coordinates)
    if (
        size == 0
        or size % 2
        or len(alpha_values) != size
        or len(auxiliary_coordinates) != size
        or len(target_partner) != size
        or len(sigma) != size
        or len(signs) != size
    ):
        raise ValueError("all assignment lists must have one common positive even size")
    targets = tuple(value % p for value in target_coordinates)
    alphas = tuple(value % p for value in alpha_values)
    auxiliaries = tuple(value % p for value in auxiliary_coordinates)
    tau = tuple(int(value) for value in target_partner)
    sigma = tuple(int(value) for value in sigma)
    signs = tuple(int(value) for value in signs)
    c %= p
    if len(set(targets)) != size:
        raise ValueError("target directions must be distinct")
    if any(value == 0 for value in alphas) or c == 0:
        raise ValueError("alpha values and c must be nonzero")
    if sorted(sigma) != list(range(size)):
        raise ValueError("sigma must permute the labeled auxiliary occurrences")
    if any(index < 0 or index >= size for index in tau):
        raise ValueError("a target partner index is out of range")
    if any(tau[index] == index or tau[tau[index]] != index for index in range(size)):
        raise ValueError("target_partner must be a fixed-point-free involution")
    if any(sign not in (-1, 1) for sign in signs):
        raise ValueError("all singleton signs must be +1 or -1")

    endpoint_equations: list[bool] = []
    for index in range(size):
        own_auxiliary = auxiliaries[sigma[index]]
        partner_auxiliary = auxiliaries[sigma[tau[index]]]
        left = c * alphas[index] * (
            partner_auxiliary - targets[index]
        ) % p
        right = 2 * signs[index] * (
            partner_auxiliary - own_auxiliary
        ) % p
        endpoint_equations.append(left == right)

    nonsquare_pairs: list[bool] = []
    pair_indices: list[tuple[int, int]] = []
    for index in range(size):
        partner = tau[index]
        if index > partner:
            continue
        pair_indices.append((index, partner))
        try:
            parameters = complementary_pair_parameters(
                p,
                c * alphas[index],
                c * alphas[partner],
                signs[index],
                signs[partner],
            )
        except ValueError:
            nonsquare_pairs.append(False)
        else:
            nonsquare_pairs.append(bool(parameters["complementary_nonsquare"]))

    # Cross-assign each auxiliary endpoint to the target at the other end of
    # its pair: phi(sigma(tau(i)))=i.
    phi = [-1] * size
    for index in range(size):
        phi[sigma[tau[index]]] = index
    g = tuple(
        alphas[phi[auxiliary_index]]
        * (auxiliaries[auxiliary_index] - targets[phi[auxiliary_index]])
        % p
        for auxiliary_index in range(size)
    )
    pair_invariants: list[dict[str, object]] = []
    inverse_four = _inverse(p, 4)
    for index, partner in pair_indices:
        first_auxiliary_index = sigma[index]
        second_auxiliary_index = sigma[partner]
        chord = (
            auxiliaries[first_auxiliary_index]
            - auxiliaries[second_auxiliary_index]
        ) % p
        first_g = g[first_auxiliary_index]
        second_g = g[second_auxiliary_index]
        equal_squares = first_g * first_g % p == second_g * second_g % p
        common_ratio = None
        if first_g:
            common_ratio = (
                chord
                * chord
                * _inverse(p, first_g * first_g)
            ) % p
        pair_invariants.append({
            "target_pair": (index, partner),
            "auxiliary_occurrence_pair": (
                first_auxiliary_index,
                second_auxiliary_index,
            ),
            "g_squares_equal": equal_squares,
            "chord_square_over_g_square": common_ratio,
            "expected_common_ratio": c * c * inverse_four % p,
        })

    equations_hold = all(endpoint_equations)
    invariants_hold = all(
        row["g_squares_equal"]
        and row["chord_square_over_g_square"] == row["expected_common_ratio"]
        for row in pair_invariants
    )
    occurrence_bijection = len(set(auxiliaries)) == size
    exact_assignment = bool(
        equations_hold
        and all(nonsquare_pairs)
        and invariants_hold
        and occurrence_bijection
    )
    return {
        "endpoint_equations": tuple(endpoint_equations),
        "all_endpoint_equations_hold": equations_hold,
        "all_pair_parameters_nonsquare": all(nonsquare_pairs),
        "cross_assignment_phi": tuple(phi),
        "g_values": g,
        "pair_square_invariants": tuple(pair_invariants),
        "all_pair_square_invariants_hold": invariants_hold,
        "auxiliary_occurrences_are_distinct_directions": occurrence_bijection,
        "pair_coherent_distinct_auxiliary_assignment": exact_assignment,
        "paley_type_quota_checked": False,
        "full_endpoint_assignment_proved": False,
        "necessary_multiset_condition": (
            "the cross-assigned g^2 values pair equally, with "
            "(U-V)^2/g^2=c^2/4 for every pair"
        ),
        "proved": True,
    }


def exact_magnitude_scope_guard(p: int) -> dict[str, object]:
    """Record why a bad exact magnitude is not a fixed-direction obstruction."""
    _check_branch_prime(p)
    profile = endpoint_color_profile(p, 2)
    eps = -1
    mu = 2 * _inverse(p, 4) % p
    nu = mu * mu % p
    proved = bool(
        profile["valid_signs"] == (-1,)
        and mu == _inverse(p, 2)
        and quadratic_character(p, nu) == 1
    )
    if not proved:
        raise ArithmeticError("the exact-magnitude scope guard changed")
    return {
        "aligned_evaluation_a": 2,
        "only_valid_sign": eps,
        "forced_mu": mu,
        "pair_nu": nu,
        "pair_nu_is_square": True,
        "pairing_at_this_exact_magnitude": False,
        "rescaling_x_to_c_times_x_changes_a_to": "2*c",
        "whole_fixed_direction_obstructed": False,
        "role": "scope guard, not an obstruction",
        "proved": True,
    }


def paired_square_invariant_counterexample(p: int) -> dict[str, object]:
    """Show that even square fibres do not imply one common chord ratio.

    This is an abstract invariant counterexample, not a branch-C target.
    """
    _check_branch_prime(p)
    points = (0, 1, 2, 4)
    matchings = (
        ((0, 1), (2, 4)),
        ((0, 2), (1, 4)),
        ((0, 4), (1, 2)),
    )
    squared_lengths = tuple(
        tuple((first - second) ** 2 % p for first, second in matching)
        for matching in matchings
    )
    if not all(left != right for left, right in squared_lengths):
        raise ArithmeticError("the paired-square limitation witness changed")
    return {
        "p": p,
        "auxiliary_points": points,
        "g_squared_values": (1, 1, 1, 1),
        "square_fibre_polynomial": "(T-1)^4",
        "perfect_matching_squared_chords": squared_lengths,
        "one_common_chord_ratio_exists": False,
        "paired_square_polynomial_condition_sufficient": False,
        "scope": "abstract invariant counterexample, not a branch-C target",
        "proved": True,
    }


def theorem_record(
    p: int = 31, alpha_values: Sequence[int] | None = None
) -> dict[str, object]:
    """Return the proved center-pairing theorem and explicit open frontier."""
    _, m = _check_branch_prime(p)
    if alpha_values is None:
        alpha_values = tuple(range(1, m + 1))
    counts = color_count_theorem(p)
    pairing = adaptive_center_pairing_certificate(p, alpha_values)
    guard = exact_magnitude_scope_guard(p)
    invariant_barrier = paired_square_invariant_counterexample(p)
    proved = bool(
        counts["proved"]
        and pairing["proved"]
        and guard["proved"]
        and invariant_barrier["proved"]
    )
    return {
        "title": "Adaptive complementary Mobius target pairing",
        "status": "PROVED CENTER-COHERENT TARGET MATCHING; AUXILIARY SDR OPEN",
        "color_counts": counts,
        "pairing_certificate": pairing,
        "exact_magnitude_scope_guard": guard,
        "paired_square_invariant_barrier": invariant_barrier,
        "required_global_auxiliary_set": "m distinct directions: m-2 hard and 2 opposite",
        "target_perfect_matching_proved": True,
        "cross_pair_auxiliary_distinctness_proved": False,
        "auxiliary_type_quota_proved": False,
        "full_parallel_quota_identity_constructed": False,
        "fixed_word_singleton_constructed": False,
        "common_graph_constructed": False,
        "residual_ii_closed": False,
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    from pprint import pprint

    pprint(theorem_record(), sort_dicts=True)
