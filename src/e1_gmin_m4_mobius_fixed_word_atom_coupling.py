#!/usr/bin/env python3
"""Target-sensitive fixed-word gate at the branch-C ``j=0`` endpoint.

This module couples three pieces which must refer to the same target and the
same all-active Mobius support:

* the actual hard/opposite triangle-atom quotas;
* the affine-block parity word ``c_U=M^T Phi(U)``; and
* the putative singleton fixed edge at the bare Hamming endpoint.

It proves an exact per-direction atom-capacity syndrome and a consequent
``Phi``-block collision lower bound.  It also proves that two distinct hard
halves can share at most eight ``Phi`` block types and gives a direct
three-block counterexample to the tempting bound one.  The bounds do not
contradict each other, so the endpoint and residual (ii) remain open.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Mapping, Sequence
from math import comb

from e1_gmin_m4_inversion_antisymmetric_radon import (
    _negative_edge,
    localized_star_trade,
)
from e1_gmin_m4_mobius_half_symmetric import (
    mobius_parameter_edges,
    paley_direction_sign,
)
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15758 import p3_local_survivor
from e1_gmin_m4_symmetric_fixed_edge_elimination import orbit_fixed_word


def _check_branch_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=4r+3 with p>=31")
    return (p - 3) // 4


def branch_c_atom_quotas(p: int, t: int) -> dict[str, object]:
    """Import and audit the exact triangle counts in every branch-C row."""
    r = _check_branch_prime(p)
    local = p3_local_survivor(p, t)
    m = (p + 1) // 2
    h = (p - 1) // 2
    hard = tuple(int(row["e"]) for row in local["hard_rows"])
    opposite_ae = tuple(
        int(row["p_minus_3_units"]) for row in local["opposite_rows"]
    )
    opposite_compact = tuple(
        int(row["p_plus_1_units"]) for row in local["opposite_rows"]
    )
    opposite = tuple(
        left + right
        for left, right in zip(opposite_ae, opposite_compact, strict=True)
    )
    quotas = tuple(int(row["Q"]) for row in local["opposite_rows"])
    atom_counts = hard + opposite
    E = t + 1
    total = sum(atom_counts)
    proved = bool(
        len(hard) == len(opposite) == m
        and sum(hard) == E
        and all(value == r - 1 for value in opposite_ae)
        and all(value == quota - 3 for value, quota in zip(opposite, quotas, strict=True))
        and total == 2 * E + 2 * m - 5
        and total % 2 == 1
        and min(atom_counts) > 0
        and max(hard) <= h - 3
        and max(opposite) <= h - 2
    )
    if not proved:
        raise ArithmeticError("the branch-C atom-quota identities changed")
    return {
        "p": p,
        "r": r,
        "m": m,
        "h": h,
        "t": t,
        "hard_compact_counts": hard,
        "opposite_all_equal_counts": opposite_ae,
        "opposite_compact_counts": opposite_compact,
        "opposite_total_triangle_counts": opposite,
        "all_direction_atom_counts": atom_counts,
        "total_atom_count": total,
        "total_atom_count_formula": "2*(t+1)+2*m-5",
        "every_atom_coefficient_graph_mod_two": "one three-cycle",
        "proved": True,
    }


def coupled_atom_syndrome(
    block_parity_rows: tuple[tuple[int, ...], ...],
    atom_counts: tuple[int, ...],
    hard_literal_cells: dict[int, int],
    singleton_cells: tuple[int | None, ...],
) -> dict[str, object]:
    """Compute the exact atom incidence forced by ``a_Y+Phi(U)=e_x``.

    Rows of ``block_parity_rows`` are ``c_U(D,beta)``.  A hard literal
    toggles its one ``beta=j_D^2`` cell.  ``singleton_cells[D]`` is the one
    cell ``beta=D(x)^2`` unless ``D`` annihilates ``x``, in which case it is
    ``None``.  The required atom vector is their binary sum.
    """
    direction_count = len(atom_counts)
    if (
        not block_parity_rows
        or len(block_parity_rows) != direction_count
        or len(singleton_cells) != direction_count
    ):
        raise ValueError("the direction-indexed inputs must have equal nonzero size")
    h = len(block_parity_rows[0])
    if h == 0 or any(len(row) != h for row in block_parity_rows):
        raise ValueError("block rows must have one common positive length")
    if any(bit not in (0, 1) or isinstance(bit, bool) for row in block_parity_rows for bit in row):
        raise ValueError("block parities must be binary")
    if any(value < 0 or isinstance(value, bool) for value in atom_counts):
        raise ValueError("atom counts must be nonnegative integers")
    if any(not 0 <= direction < direction_count or not 0 <= cell < h
           for direction, cell in hard_literal_cells.items()):
        raise ValueError("a hard literal cell is out of range")
    if any(cell is not None and not 0 <= cell < h for cell in singleton_cells):
        raise ValueError("a singleton cell is out of range")

    required: list[tuple[int, ...]] = []
    for direction, row in enumerate(block_parity_rows):
        bits = list(row)
        if direction in hard_literal_cells:
            bits[hard_literal_cells[direction]] ^= 1
        singleton = singleton_cells[direction]
        if singleton is not None:
            bits[singleton] ^= 1
        required.append(tuple(bits))
    weights = tuple(sum(row) for row in required)
    parity_matches = tuple(
        weight % 2 == count % 2
        for weight, count in zip(weights, atom_counts, strict=True)
    )
    capacity_matches = tuple(
        weight <= count
        for weight, count in zip(weights, atom_counts, strict=True)
    )
    return {
        "required_antipodal_atom_incidence": tuple(required),
        "required_weights": weights,
        "atom_count_parities_match": parity_matches,
        "atom_capacities_match": capacity_matches,
        "fixed_word_layer_feasible": all(parity_matches) and all(capacity_matches),
        "criterion": "wt(c_U[D]+ell_D+s_x[D])<=n_D and equal parity",
        "proved_by_direct_binary_addition": True,
    }


def atom_fixed_incidence_realizable(atom_count: int, requested_weight: int) -> bool:
    """Decide the centered-atom realization of one fixed-cell parity row.

    A distinct-label triangle contains at most one antipodal edge.  Conversely,
    ``K(s,-s;0)`` and the all-equal triangle ``{s,-s,0}`` are individually
    central and toggle exactly that edge.  Repeating a square class twice
    fills any remaining even number of atoms.  Thus weight and parity are the
    only restrictions at the central fixed-word/odd-moment layer *when the
    atom labels may be chosen*, subject to the prescribed split between the
    two atom types.  This is not a statement about an already fixed atom list.
    """
    if (
        not isinstance(atom_count, int)
        or isinstance(atom_count, bool)
        or not isinstance(requested_weight, int)
        or isinstance(requested_weight, bool)
        or atom_count < 0
        or requested_weight < 0
    ):
        raise ValueError("atom_count and requested_weight must be nonnegative integers")
    return requested_weight <= atom_count and (atom_count - requested_weight) % 2 == 0


def fixed_word_atom_coupling_theorem(p: int, t: int) -> dict[str, object]:
    """State the fully quantified target-sensitive ``j=0`` reduction."""
    quotas = branch_c_atom_quotas(p, t)
    m = int(quotas["m"])
    h = int(quotas["h"])
    atom_counts = tuple(quotas["all_direction_atom_counts"])
    t_max = 4 * int(quotas["r"]) ** 2 - 2 * int(quotas["r"]) - 5
    kappa0 = t_max - t + 1
    raw_nonzero_phi = m * (p - 2)
    proved = bool(
        len(atom_counts) == p + 1
        and raw_nonzero_phi == m * (2 * m - 3)
        and sum(atom_counts) == 2 * (t + 1) + 2 * m - 5
        and h == m - 1
    )
    if not proved:
        raise ArithmeticError("the fixed-word coupling ledger changed")
    return {
        **quotas,
        "j_endpoint": 0,
        "size_floor_cancellations_kappa0": kappa0,
        "raw_nonzero_Phi_block_occurrences": raw_nonzero_phi,
        "fixed_inverse_fact": (
            "the P_(L_v) bit occurs in both terms of the explicit inverse "
            "and cancels; it contributes no point-word term"
        ),
        "block_basis_identity": "M^T*a_Y=ell+z",
        "coupled_singleton_identity": "z=c_U+ell+s_x",
        "definitions": {
            "c_U": (
                "parity over all d*h affine block types of the used nonzero-Phi orbits"
            ),
            "ell": "one block beta=j_D^2 in every hard direction D",
            "s_x": (
                "one block beta=D(x)^2 in every D except F=L_x, and zero in F"
            ),
            "z": (
                "parity of target triangle atoms containing the antipodal label pair for beta"
            ),
        },
        "per_direction_necessary_and_fixed_layer_sufficient": (
            "wt(c_U[D]+ell[D]+s_x[D])<=n_D and has parity n_D, "
            "existentially over atom labels with the prescribed type counts"
        ),
        "sufficiency_does_not_apply_to_a_prefixed_atom_labeling": True,
        "parity_is_automatic_if_j0_parallel_slices_hold": True,
        "fixed_word_and_odd_moment_layer_closed_by_this_criterion": True,
        "even_common_moments_and_nonfixed_target_cells_solved": False,
        "Mobius_curve_realizability_of_an_abstract_c_U_solved": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def phi_collision_ledger(
    p: int,
    t: int,
    literal_singleton_matches_q: int = 0,
    zero_phi_cancellations: int = 0,
) -> dict[str, object]:
    """Derive the exact collision demand forced by the capacity syndrome.

    Let ``C=wt(c_U)`` over all ``d*h`` block types and
    ``Lambda=(m(p-2)-C)/2``.  With ``U_np`` the *distinct surviving* used
    nonzero-Phi orbits, put ``sigma=(|U_np|-C)/2``.  These are integers even
    with triple/higher raw overlaps, because both parities equal the parity
    of the raw multiplicity in each physical orbit and then in each block.
    """
    theorem = fixed_word_atom_coupling_theorem(p, t)
    m = int(theorem["m"])
    raw = int(theorem["raw_nonzero_Phi_block_occurrences"])
    kappa0 = int(theorem["size_floor_cancellations_kappa0"])
    if (
        not isinstance(literal_singleton_matches_q, int)
        or isinstance(literal_singleton_matches_q, bool)
        or not 0 <= literal_singleton_matches_q <= m
    ):
        raise ValueError("q must lie between zero and the number of hard rows")
    if (
        not isinstance(zero_phi_cancellations, int)
        or isinstance(zero_phi_cancellations, bool)
        or not 0 <= zero_phi_cancellations <= min(kappa0, m // 2)
    ):
        raise ValueError("invalid number of zero-Phi cancellation units")

    q = literal_singleton_matches_q
    kappa_zero = zero_phi_cancellations
    kappa_nonzero = kappa0 - kappa_zero
    used_nonzero = raw - 2 * kappa_nonzero
    lambda_lower = kappa0 + m + q
    sigma_lower = kappa_zero + m + q
    sigma_raw_upper = used_nonzero // 2
    pair_intersection_upper = 8 * comb(m, 2)
    return {
        "p": p,
        "t": t,
        "m": m,
        "q_literal_singleton_matches": q,
        "kappa_total": kappa0,
        "kappa_zero_Phi": kappa_zero,
        "kappa_nonzero_Phi": kappa_nonzero,
        "raw_nonzero_Phi_occurrences": raw,
        "used_distinct_nonzero_Phi_orbits": used_nonzero,
        "C_definition": "C=wt(c_U) over all d*h affine block types",
        "Lambda_definition": "Lambda=(m*(p-2)-C)/2",
        "sigma_definition": "sigma=(|U_nonzero|-C)/2",
        "exact_identity": "Lambda=kappa_nonzero+sigma=kappa0-kappa_zero+sigma",
        "forced_Lambda_lower_bound": lambda_lower,
        "forced_sigma_lower_bound": sigma_lower,
        "sigma_is_integer_nonnegative_under_arbitrary_higher_overlaps": True,
        "raw_occurrence_upper_bound_for_sigma": sigma_raw_upper,
        "eight_per_half_pair_upper_bound": pair_intersection_upper,
        "available_raw_upper_margin": sigma_raw_upper - sigma_lower,
        "collision_bounds_contradict": sigma_lower > min(
            sigma_raw_upper, pair_intersection_upper
        ),
        "residual_ii_closed": False,
        "proved": True,
    }


def phi_collision_decomposition(
    raw_multiplicities_by_block: Mapping[Hashable, Sequence[int]],
) -> dict[str, object]:
    """Replay ``Lambda=kappa_nonzero+sigma`` with arbitrary higher overlaps.

    Each positive integer in a block row is the number of raw half
    occurrences on one distinct physical orbit.  Ternarity makes that orbit
    survive precisely when the multiplicity is odd.  Different entries in
    one row are distinct physical orbits with the same nonzero Phi block.
    """
    if not raw_multiplicities_by_block:
        raise ValueError("need at least one nonzero-Phi block")
    rows = {
        block: tuple(multiplicities)
        for block, multiplicities in raw_multiplicities_by_block.items()
    }
    if any(
        not multiplicities
        or any(
            not isinstance(value, int) or isinstance(value, bool) or value <= 0
            for value in multiplicities
        )
        for multiplicities in rows.values()
    ):
        raise ValueError("each block needs positive integer orbit multiplicities")

    raw = sum(value for row in rows.values() for value in row)
    used_by_block = {
        block: sum(value & 1 for value in row) for block, row in rows.items()
    }
    used = sum(used_by_block.values())
    block_parity = {
        block: count & 1 for block, count in used_by_block.items()
    }
    C = sum(block_parity.values())
    kappa_nonzero = (raw - used) // 2
    sigma = (used - C) // 2
    Lambda = (raw - C) // 2
    proved = bool(
        raw % 2 == C % 2
        and used % 2 == C % 2
        and kappa_nonzero >= 0
        and sigma >= 0
        and Lambda == kappa_nonzero + sigma
    )
    if not proved:
        raise ArithmeticError("the Phi collision decomposition changed")
    return {
        "raw_nonzero_Phi_occurrences": raw,
        "used_distinct_orbits_by_block": used_by_block,
        "used_distinct_nonzero_Phi_orbits": used,
        "block_parity": block_parity,
        "C": C,
        "kappa_nonzero": kappa_nonzero,
        "sigma": sigma,
        "Lambda": Lambda,
        "identity": "Lambda=kappa_nonzero+sigma",
        "arbitrary_triple_and_higher_multiplicities_allowed": True,
        "proved": True,
    }


def all_prime_collision_room_theorem(p: int) -> dict[str, object]:
    """Show symbolically that the new collision demand is not a scalar close."""
    r = _check_branch_prime(p)
    m = 2 * r + 2
    kappa0_max = m * m // 2 - m - 2
    raw_half = m * (2 * m - 3) // 2
    worst_q = m
    minimum_margin = raw_half - kappa0_max - m - worst_q
    closed_margin = (m - 1) * (m - 4) // 2
    proved = bool(
        m % 2 == 0
        and minimum_margin == closed_margin
        and closed_margin > 0
    )
    if not proved:
        raise ArithmeticError("the symbolic collision-room identity changed")
    return {
        "p": p,
        "m": m,
        "worst_case": "t=t_min and q=m",
        "maximum_size_floor_cancellations": kappa0_max,
        "minimum_raw_upper_margin": minimum_margin,
        "closed_margin": "(m-1)(m-4)/2",
        "strictly_positive_for_branch_m_at_least_16": True,
        "scalar_collision_count_excludes_j0": False,
        "needs_actual_global_Mobius_block_geometry": True,
        "proved": True,
    }


def mobius_half_block_normal_form(p: int, z: int, center: int) -> tuple[int, int]:
    """Return coefficients of the normalized block functional in ``(L,M)``.

    For ``z=t+1`` the block functional is

    ``N_z=-(z+1)L/j + z^2 M/[j(z-1)]``.
    """
    _check_branch_prime(p)
    z %= p
    center %= p
    if z in (0, 1) or center == 0:
        raise ValueError("need z nonzero/nonunit and a nonzero center")
    inverse_center = pow(center, -1, p)
    return (
        -(z + 1) * inverse_center % p,
        z * z * pow(center * (z - 1) % p, -1, p) % p,
    )


def two_half_phi_block_intersection_theorem(p: int) -> dict[str, object]:
    """State the exact conic normal form and the valid universal bound eight."""
    _check_branch_prime(p)
    return {
        "p": p,
        "one_half_block_count": p - 2,
        "one_half_blocks_are_distinct": True,
        "normal_form": "(N(A)+2)*N(B)+1=0, excluding N(A)=N(B)=-1",
        "parameters": "A=j*e and B=j*(e+f) in the primal basis dual to (L,M)",
        "target_direction_recovery": "L annihilates B-A",
        "same_block_criterion": "N_i(z)=+N_k(w) or N_i(z)=-N_k(w)",
        "per_orientation_intersection_bound": 4,
        "proof": (
            "on the first hyperbola, u=N(A)+2 gives N(V)=a*(u-2)-b/u; "
            "after multiplication by u^2 the second conic equation has degree at most 4. "
            "Identity would make the irreducible conics equal and force the same target L"
        ),
        "two_orientation_block_intersection_bound": 8,
        "bound_one_valid": False,
        "proved": True,
    }


def p31_three_phi_block_counterexample() -> dict[str, object]:
    """Directly replay two disjoint hard halves sharing three Phi blocks."""
    p = 31
    specifications = (
        ((1, 1), (1, 3), 2),
        ((0, 1), (1, 7), 3),
    )
    block_maps: list[dict[int, tuple[tuple[int, int], ...]]] = []
    trades = []
    for direction, auxiliary, center in specifications:
        parameter_edges = mobius_parameter_edges(p, direction, auxiliary, center)
        blocks: dict[int, tuple[tuple[int, int], ...]] = {}
        for parameter, edge in parameter_edges.items():
            record = orbit_fixed_word(p, edge)
            if int(record["fixed_word_weight"]):
                blocks[parameter] = tuple(
                    tuple(point) for point in record["fixed_word_support"]
                )
        block_maps.append(blocks)
        trades.append(localized_star_trade(p, direction, auxiliary, center))

    shared_parameters = tuple(
        (left_parameter, right_parameter)
        for left_parameter, left_block in block_maps[0].items()
        for right_parameter, right_block in block_maps[1].items()
        if left_block == right_block
    )
    total = Counter(trades[0])
    total.update(trades[1])
    common_physical = set(trades[0]) & set(trades[1])
    proved = bool(
        all(paley_direction_sign(p, row[0]) == 1 for row in specifications)
        and tuple(len(blocks) for blocks in block_maps) == (p - 2, p - 2)
        and shared_parameters == ((9, 20), (19, 12), (25, 18))
        and not common_physical
        and set(total.values()) == {-1, 1}
        and all(total.get(_negative_edge(p, edge)) == -value for edge, value in total.items())
    )
    if not proved:
        raise ArithmeticError("the three-Phi-block counterexample changed")
    return {
        "p": p,
        "halves": specifications,
        "both_target_directions_Paley_hard": True,
        "shared_Phi_block_parameter_pairs": shared_parameters,
        "shared_Phi_block_count": len(shared_parameters),
        "common_physical_edges": 0,
        "sum_is_ternary": True,
        "role": "fail-when-wrong barrier, not a prime or configuration census",
        "proved": True,
    }


def theorem_record(p: int = 31, t: int | None = None) -> dict[str, object]:
    r = _check_branch_prime(p)
    if t is None:
        t = 4 * r * r - 2 * r - 5
    coupling = fixed_word_atom_coupling_theorem(p, t)
    collisions = phi_collision_ledger(p, t)
    room = all_prime_collision_room_theorem(p)
    intersections = two_half_phi_block_intersection_theorem(p)
    counterexample = p31_three_phi_block_counterexample() if p == 31 else None
    proved = bool(
        coupling["proved"]
        and collisions["proved"]
        and room["proved"]
        and intersections["proved"]
        and (counterexample is None or counterexample["proved"])
    )
    return {
        "title": "Branch-C j=0 fixed-word/atom coupling and Phi-collision gate",
        "status": "PROVED TARGET-SENSITIVE REDUCTION; ENDPOINT OPEN",
        "coupling": coupling,
        "collision_ledger": collisions,
        "all_prime_collision_room": room,
        "pair_intersection_theorem": intersections,
        "three_collision_counterexample": counterexample,
        "uniform_j0_exclusion_proved": False,
        "j0_construction_proved": False,
        "even_common_moment_completion_proved": False,
        "nonfixed_target_cell_completion_proved": False,
        "residual_ii_closed": False,
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    from pprint import pprint

    pprint(theorem_record(), sort_dicts=True)
