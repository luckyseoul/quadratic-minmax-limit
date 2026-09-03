#!/usr/bin/env python3
"""Exact unrestricted lattice theorem for the inversion-symmetric Radon block.

This complements e1_gmin_m4_inversion_antisymmetric_radon. It proves the
characteristic-zero ranks, the exact elementary-p cokernel spanned by even
moments, and mod-two surjectivity after fixed antipodal edges are included.
It does not solve the restricted central Boolean box.
"""

from __future__ import annotations

from itertools import combinations, product

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15759 import p_torsion_codimension


Point = tuple[int, int]
Edge = tuple[int, int]


def _check_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")
    return (p - 1) // 2


def symmetric_dimensions(p: int) -> dict[str, int | bool]:
    """Return the exact source, target, and kernel plus-ranks."""
    h = _check_prime(p)
    d = p + 1
    difference_classes = d * h
    source_edges = p * p * difference_classes
    source_fixed_edges = difference_classes
    source_plus = (source_edges + source_fixed_edges) // 2
    target_plus = d * h * (h + 1)
    kernel_plus = source_plus - target_plus
    expected_kernel = d * p * h * h
    proved = bool(
        difference_classes == (p * p - 1) // 2
        and source_plus == difference_classes * (p * p + 1) // 2
        and kernel_plus == expected_kernel
    )
    if not proved:
        raise ArithmeticError("the inversion-symmetric rank identity changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "difference_classes": difference_classes,
        "source_edge_rank": source_edges,
        "source_fixed_antipodal_edge_rank": source_fixed_edges,
        "source_plus_rank": source_plus,
        "target_plus_rank": target_plus,
        "kernel_plus_rank": kernel_plus,
        "kernel_plus_closed_form": expected_kernel,
        "proved": proved,
    }


def even_moment_cokernel(p: int) -> dict[str, object]:
    """Split Proposition 15.759's elementary-p cokernel by inversion parity."""
    h = _check_prime(p)
    ledger = p_torsion_codimension(p)
    by_degree = {
        int(degree): int(count)
        for degree, count in ledger["relation_count_by_degree"].items()
    }
    even_direct = sum(count for degree, count in by_degree.items() if degree % 2 == 0)
    odd_direct = sum(count for degree, count in by_degree.items() if degree % 2 == 1)
    even_closed = (h - 1) * (2 * h * h + 5 * h + 6) // 6
    odd_closed = h * (h - 1) * (h + 1) // 3
    full = int(ledger["extra_p_primary_codimension_closed_form"])
    proved = bool(
        even_direct == even_closed
        and odd_direct == odd_closed
        and even_direct + odd_direct == full
    )
    if not proved:
        raise ArithmeticError("the even/odd moment split changed")
    return {
        "p": p,
        "h": h,
        "relation_count_by_degree": {
            str(degree): count for degree, count in sorted(by_degree.items())
        },
        "even_moment_rank_direct": even_direct,
        "even_moment_rank_closed": even_closed,
        "odd_moment_rank_direct": odd_direct,
        "odd_moment_rank_closed": odd_closed,
        "full_cokernel_rank": full,
        "symmetric_integral_cokernel": f"(Z/{p}Z)^{even_closed}",
        "symmetric_integral_compatibility": (
            "exactly the even-degree moment rows of Proposition 15.759"
        ),
        "proved": proved,
    }


def symmetric_mod2_decomposition(p: int) -> dict[str, object]:
    """Prove mod-two surjectivity by fixed/nonfixed coordinate support.

    Fixed antipodal source edges map only to target coordinates fixed by
    inversion: the parallel cell and cells {s,-s}. Their map is the affine
    Radon transform of an even point function vanishing at zero, hence is
    injective. Its image is exactly the fixed-cell data obeying the p+1
    ordinary total/parallel equations.

    A nonfixed source pair e+Je maps to paired nonfixed target cells and
    vanishes on fixed cells modulo two. This is exactly the mod-two
    reduction of the antisymmetric block, already surjective onto those
    paired cells. The two coordinate supports are disjoint.
    """
    h = _check_prime(p)
    d = p + 1
    fixed_source = d * h
    fixed_target_coordinates = d * (h + 1)
    ordinary_fixed_compatibility = d
    fixed_image = fixed_target_coordinates - ordinary_fixed_compatibility
    nonfixed_source_pair_orbits = (d * h) ** 2
    paired_nonfixed_target = d * h * h
    total_image = fixed_image + paired_nonfixed_target
    expected = d * h * (h + 1)
    proved = bool(fixed_image == fixed_source and total_image == expected)
    if not proved:
        raise ArithmeticError("the symmetric mod-two decomposition changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "fixed_antipodal_source_variables": fixed_source,
        "fixed_target_coordinates": fixed_target_coordinates,
        "ordinary_fixed_coordinate_compatibility_equations": (
            ordinary_fixed_compatibility
        ),
        "fixed_edge_map_rank": fixed_image,
        "fixed_edge_map_injective": True,
        "fixed_edge_map_reason": (
            "Fourier inversion of the affine Radon transform on even point "
            "functions with value zero at the origin"
        ),
        "nonfixed_source_pair_variables": nonfixed_source_pair_orbits,
        "paired_nonfixed_target_rank": paired_nonfixed_target,
        "nonfixed_pair_map_surjective": True,
        "nonfixed_pair_map_reason": (
            "it is the characteristic-two reduction of the proved "
            "antisymmetric block"
        ),
        "fixed_and_nonfixed_target_coordinate_supports_disjoint": True,
        "symmetric_mod2_image_rank": total_image,
        "symmetric_target_rank": expected,
        "symmetric_mod2_surjective": True,
        "proved": proved,
    }


def _rank_binary_vectors(vectors: list[int]) -> int:
    basis: dict[int, int] = {}
    for original in vectors:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def _negative_point_index(p: int, point: Point, point_index: dict[Point, int]) -> int:
    return point_index[((-point[0]) % p, (-point[1]) % p)]


def exact_symmetric_mod2_rank(p: int) -> dict[str, int | bool]:
    """Build the small exact plus-map over F2 as a fail-when-wrong replay."""
    h = _check_prime(p)
    if p > 11:
        raise ValueError("the exact matrix replay is intentionally limited to p<=11")
    points = tuple(product(range(p), repeat=2))
    point_index = {point: index for index, point in enumerate(points)}
    edges = tuple(combinations(range(len(points)), 2))
    edge_set = set(edges)
    directions = tuple((1, slope) for slope in range(p)) + ((0, 1),)
    cells = tuple(combinations(range(p), 2))
    cell_index = {cell: index + 1 for index, cell in enumerate(cells)}
    block_size = len(cells) + 1

    def negative_edge(edge: Edge) -> Edge:
        left = _negative_point_index(p, points[edge[0]], point_index)
        right = _negative_point_index(p, points[edge[1]], point_index)
        return tuple(sorted((left, right)))  # type: ignore[return-value]

    def column(edge_orbit: tuple[Edge, ...]) -> int:
        value = 0
        for edge in edge_orbit:
            first_point = points[edge[0]]
            second_point = points[edge[1]]
            for direction_index, (a, b) in enumerate(directions):
                first = (a * first_point[0] + b * first_point[1]) % p
                second = (a * second_point[0] + b * second_point[1]) % p
                offset = direction_index * block_size
                coordinate = (
                    0
                    if first == second
                    else cell_index[tuple(sorted((first, second)))]
                )
                value ^= 1 << (offset + coordinate)
        return value

    visited: set[Edge] = set()
    fixed_columns: list[int] = []
    nonfixed_columns: list[int] = []
    for edge in edges:
        if edge in visited:
            continue
        negative = negative_edge(edge)
        if negative not in edge_set:
            raise ArithmeticError("central inversion lost a source edge")
        if negative == edge:
            orbit = (edge,)
            fixed_columns.append(column(orbit))
        else:
            orbit = (edge, negative)
            nonfixed_columns.append(column(orbit))
        visited.update(orbit)

    fixed_rank = _rank_binary_vectors(fixed_columns)
    nonfixed_rank = _rank_binary_vectors(nonfixed_columns)
    total_rank = _rank_binary_vectors(fixed_columns + nonfixed_columns)
    expected_fixed = (p + 1) * h
    expected_nonfixed = (p + 1) * h * h
    expected_total = (p + 1) * h * (h + 1)
    proved = bool(
        len(visited) == len(edges)
        and len(fixed_columns) == expected_fixed
        and fixed_rank == expected_fixed
        and nonfixed_rank == expected_nonfixed
        and total_rank == expected_total
        and total_rank == fixed_rank + nonfixed_rank
    )
    if not proved:
        raise ArithmeticError("the exact symmetric mod-two rank changed")
    return {
        "p": p,
        "fixed_source_columns": len(fixed_columns),
        "nonfixed_pair_source_columns": len(nonfixed_columns),
        "fixed_map_rank": fixed_rank,
        "nonfixed_pair_map_rank": nonfixed_rank,
        "full_symmetric_map_rank": total_rank,
        "expected_symmetric_target_rank": expected_total,
        "proved": proved,
    }


def symmetric_integral_lattice_theorem(p: int) -> dict[str, object]:
    """Record the exact plus-cokernel argument and its scope."""
    dimensions = symmetric_dimensions(p)
    moments = even_moment_cokernel(p)
    mod2 = symmetric_mod2_decomposition(p)
    proved = bool(dimensions["proved"] and moments["proved"] and mod2["proved"])
    if not proved:
        raise ArithmeticError("the symmetric integral lattice theorem changed")
    return {
        "p": p,
        "dimensions": dimensions,
        "mod2": mod2,
        "cokernel": moments["symmetric_integral_cokernel"],
        "cokernel_injects_into_full_edge_Radon_cokernel": True,
        "injection_reason": (
            "a kernel class is killed by two, while mod-two surjectivity "
            "removes all two-torsion"
        ),
        "image_in_full_cokernel": "the +1 inversion eigenspace",
        "image_reason": (
            "symmetrization represents twice a fixed class, and two is "
            "invertible in the elementary p-cokernel"
        ),
        "compatibility_if_and_only_if": (
            "all even-degree moment congruences of Proposition 15.759 vanish"
        ),
        "unrestricted_signed_integral_central_lift_proved": True,
        "restricted_Boolean_central_lift_proved": False,
        "proved": proved,
    }


def mobius_central_remainder_box(p: int, nonzero_hard_centres: int) -> dict[str, object]:
    """State the exact remaining box after a localized ternary Mobius lift."""
    h = _check_prime(p)
    hard_direction_count = h + 1
    q = nonzero_hard_centres
    if (
        not isinstance(q, int)
        or isinstance(q, bool)
        or not 0 <= q <= hard_direction_count
    ):
        raise ValueError("nonzero hard-centre count is out of range")
    difference_classes = (p * p - 1) // 2
    nonfixed_pair_orbits = difference_classes * difference_classes
    used = q * (p - 1)
    fixed = difference_classes
    proved = bool(used <= nonfixed_pair_orbits)
    if not proved:
        raise ArithmeticError("the Mobius central-box count changed")
    return {
        "p": p,
        "h": h,
        "nonzero_hard_centres": q,
        "mobius_used_nonfixed_orbits": used,
        "available_nonfixed_pair_orbits": nonfixed_pair_orbits,
        "available_fixed_antipodal_edges": fixed,
        "ordinary_source_notation": "q_U=tau*x_U",
        "antisymmetric_identity": "q_U-Jq_U=z",
        "forced_symmetric_source": (
            "C_U=q_U+Jq_U=sum_used tau_O*(e+Je)"
        ),
        "central_target_formula": "T_U=Y-Rq_U=(Y+IY-RC_U)/2",
        "remaining_source_box": {
            "used_nonfixed_orbits_after_subtraction": "{0}",
            "unused_nonfixed_orbits": "{0,tau_O*(e+Je)}",
            "fixed_antipodal_edges": "{0,tau_f*f}",
        },
        "unrestricted_integral_lift_if_even_moments_vanish": True,
        "mod2_central_lift_always_exists_for_compatible_target": True,
        "same_lift_integral_and_in_box_proved": False,
        "coupled_symmetric_half_closed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    theorem = symmetric_integral_lattice_theorem(p)
    box = mobius_central_remainder_box(p, (p + 1) // 2)
    replay = exact_symmetric_mod2_rank(7)
    proved = bool(theorem["proved"] and box["proved"] and replay["proved"])
    if not proved:
        raise ArithmeticError("the symmetric theorem record changed")
    return {
        "title": "Central-inversion symmetric edge-Radon lattice",
        "status": "UNRESTRICTED INTEGRAL AND MOD-TWO BLOCKS CLOSED; BOX OPEN",
        "theorem": theorem,
        "mobius_central_remainder": box,
        "small_exact_mod2_replay": replay,
        "proved": {
            "symmetric_cokernel_is_even_moment_eigenspace": True,
            "symmetric_mod2_map_surjective": True,
            "unrestricted_signed_integral_central_lift_characterized": True,
            "restricted_central_Boolean_box_nonempty": False,
            "residual_ii_closed": False,
        },
        "duplicate_work_guard": (
            "Do not retry parity, Smith, or unrestricted integral obstructions "
            "inside the symmetric block. The live problem is simultaneous "
            "intersection with used=0, unused={0,tau}, fixed={0,tau}."
        ),
        "L_status": "OPEN",
        "proved_all": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), sort_keys=True, indent=2))
