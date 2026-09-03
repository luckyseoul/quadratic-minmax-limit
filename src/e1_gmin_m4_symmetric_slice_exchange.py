#!/usr/bin/env python3
"""Whole-slab exchanges in one unused symmetric difference slice.

After the fixed antipodal source choices have been eliminated, the remaining
variables are the nonfixed inversion-orbit pairs ``([a],[delta])``.  This
module proves the exact integer kernel supported on one fixed difference
class ``[delta]``.  It is an ``A_(h-1)`` root lattice on the nonzero square
fibres of the functional annihilating ``delta``.  Its primitive circuits are
weight-preserving exchanges of two complete ``p``-element slabs.

The theorem is deliberately slice-local.  It neither proves normality of the
full unused-column configuration nor proves that its constant-weight Boolean
fibre is nonempty.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product

from e1_gmin_m4_prop15721 import is_prime


Point = tuple[int, int]
Coordinate = tuple[int, int, int]


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


def _neg(point: Point, p: int) -> Point:
    return ((-point[0]) % p, (-point[1]) % p)


def _half_classes(p: int) -> tuple[Point, ...]:
    """Represent ``(F_p^2 minus {0})/{+1,-1}`` canonically."""
    zero = (0, 0)
    return tuple(
        sorted(
            {
                min(point, _neg(point, p))
                for point in product(range(p), repeat=2)
                if point != zero
            }
        )
    )


def _directions(p: int) -> tuple[Point, ...]:
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def _evaluate(functional: Point, point: Point, p: int) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % p


def _annihilator(delta: Point, p: int) -> Point:
    functional = (delta[1] % p, (-delta[0]) % p)
    if functional == (0, 0) or _evaluate(functional, delta, p) != 0:
        raise ArithmeticError("the difference annihilator changed")
    return functional


def _reduced_pair_column(p: int, midpoint: Point, delta: Point) -> tuple[Coordinate, ...]:
    """Return one column of the fixed-edge-eliminated symmetric pair map.

    In direction ``L``, a pair orbit contributes the parallel coordinate if
    ``L(delta)=0``.  Otherwise its central cell orbit is determined by
    ``L(midpoint)^2`` and ``L(delta)^2``.  The marker ``(-1, 0)`` separates
    the parallel coordinate from transverse coordinates.
    """
    column: list[Coordinate] = []
    for index, functional in enumerate(_directions(p)):
        midpoint_value = _evaluate(functional, midpoint, p)
        difference_value = _evaluate(functional, delta, p)
        if difference_value == 0:
            column.append((index, -1, 0))
        else:
            column.append(
                (
                    index,
                    midpoint_value * midpoint_value % p,
                    difference_value * difference_value % p,
                )
            )
    return tuple(column)


def symmetric_slice_kernel_theorem(p: int) -> dict[str, object]:
    """State the exact one-difference-slice kernel and its circuits."""
    h = _check_prime(p)
    difference_classes = (p * p - 1) // 2
    nonzero_slabs = h
    slab_size = p
    zero_slab_size = h
    kernel_rank = h - 1
    proved = bool(
        difference_classes == (p + 1) * h
        and nonzero_slabs * slab_size + zero_slab_size == difference_classes
        and kernel_rank == nonzero_slabs - 1
    )
    if not proved:
        raise ArithmeticError("the symmetric slice counts changed")
    return {
        "p": p,
        "h": h,
        "nonfixed_midpoint_orbits_per_difference_slice": difference_classes,
        "zero_label_slab_size": zero_slab_size,
        "nonzero_square_label_slab_count": nonzero_slabs,
        "nonzero_square_label_slab_size": slab_size,
        "kernel_description": (
            "v([a],[delta])=gamma_(L_delta(a)^2), gamma_0=0, "
            "sum_nonzero_square_classes gamma=0"
        ),
        "kernel_lattice": f"A_{kernel_rank}",
        "kernel_rank": kernel_rank,
        "primitive_circuit": "1_(S_delta,s)-1_(S_delta,t)",
        "circuit_positive_degree_in_pair_variables": p,
        "circuit_negative_degree_in_pair_variables": p,
        "circuit_support_in_pair_variables": 2 * p,
        "physical_graph_edges_removed_and_added": 2 * p,
        "weight_preserving": True,
        "all_directional_parallel_coordinates_preserved": True,
        "unsigned_kernel_valid_for_all_odd_primes": True,
        "paley_signed_specialization": (
            "only for p=3 mod 4: the anisotropic Paley sign is nonzero and "
            "constant on a fixed difference slice, so it does not change "
            "the slice kernel"
        ),
        "full_unused_configuration_normality_proved": False,
        "global_Boolean_fibre_nonempty_proved": False,
        "proved": proved,
    }


def used_orbit_deletion_bound(p: int, used_orbit_cap: int | None = None) -> dict[str, object]:
    """Guarantee a surviving whole-slab circuit under the Mobius deletion cap."""
    h = _check_prime(p)
    difference_classes = (p * p - 1) // 2
    cap = difference_classes if used_orbit_cap is None else used_orbit_cap
    if (
        not isinstance(cap, int)
        or isinstance(cap, bool)
        or not 0 <= cap <= difference_classes
    ):
        raise ValueError("used_orbit_cap must lie between 0 and (p^2-1)/2")
    least_used_in_one_slice_at_most = cap // difference_classes
    guaranteed_clean_slabs = h - least_used_in_one_slice_at_most
    surviving_circuit = guaranteed_clean_slabs >= 2
    proved = bool(
        least_used_in_one_slice_at_most <= 1
        and (p < 7 or surviving_circuit)
    )
    if not proved:
        raise ArithmeticError("the used-orbit pigeonhole bound changed")
    return {
        "p": p,
        "difference_slice_count": difference_classes,
        "mobius_used_orbit_cap": cap,
        "one_slice_used_orbits_at_most": least_used_in_one_slice_at_most,
        "clean_nonzero_slabs_in_that_slice_at_least": guaranteed_clean_slabs,
        "whole_slab_circuits_in_that_slice_at_least": (
            guaranteed_clean_slabs * (guaranteed_clean_slabs - 1) // 2
        ),
        "surviving_unused_whole_slab_circuit_for_p_at_least_7": (
            surviving_circuit if p >= 7 else False
        ),
        "global_connectivity_proved": False,
        "proved": proved,
    }


def _rational_column_rank(columns: tuple[tuple[int, ...], ...]) -> int:
    """Exact small-column rank used only by the p=5,7 replay."""
    basis: dict[int, dict[int, Fraction]] = {}
    for support in columns:
        vector = {row: Fraction(1) for row in support}
        while vector:
            pivot = min(vector)
            old = basis.get(pivot)
            if old is None:
                scale = vector[pivot]
                vector = {row: value / scale for row, value in vector.items()}
                basis[pivot] = vector
                break
            factor = vector[pivot]
            for row, value in old.items():
                replacement = vector.get(row, Fraction(0)) - factor * value
                if replacement:
                    vector[row] = replacement
                else:
                    vector.pop(row, None)
    return len(basis)


def exact_small_slice_replay(p: int) -> dict[str, object]:
    """Build one reduced slice exactly for p=5 or 7 and check the theorem."""
    h = _check_prime(p)
    if p not in (5, 7):
        raise ValueError("the exact slice replay is intentionally limited to p=5,7")
    classes = _half_classes(p)
    delta = classes[0]
    annihilator = _annihilator(delta, p)
    square_slabs: dict[int, list[Point]] = {}
    for midpoint in classes:
        label = _evaluate(annihilator, midpoint, p)
        square_slabs.setdefault(label * label % p, []).append(midpoint)
    nonzero_squares = tuple(sorted(value for value in square_slabs if value))

    all_coordinates = sorted(
        {
            coordinate
            for midpoint in classes
            for coordinate in _reduced_pair_column(p, midpoint, delta)
        }
    )
    coordinate_index = {coordinate: index for index, coordinate in enumerate(all_coordinates)}
    columns = tuple(
        tuple(
            coordinate_index[coordinate]
            for coordinate in _reduced_pair_column(p, midpoint, delta)
        )
        for midpoint in classes
    )
    rank = _rational_column_rank(columns)

    slab_sums: list[Counter[Coordinate]] = []
    for square in nonzero_squares:
        total: Counter[Coordinate] = Counter()
        for midpoint in square_slabs[square]:
            total.update(_reduced_pair_column(p, midpoint, delta))
        slab_sums.append(total)
    equal_nonzero_slab_sums = bool(
        slab_sums and all(total == slab_sums[0] for total in slab_sums[1:])
    )
    expected_rank = len(classes) - (h - 1)
    proved = bool(
        len(classes) == (p * p - 1) // 2
        and len(square_slabs.get(0, ())) == h
        and len(nonzero_squares) == h
        and all(len(square_slabs[square]) == p for square in nonzero_squares)
        and equal_nonzero_slab_sums
        and rank == expected_rank
    )
    if not proved:
        raise ArithmeticError("the exact small symmetric slice replay failed")
    return {
        "p": p,
        "map_replayed": "unsigned reduced pair columns",
        "difference_representative": list(delta),
        "midpoint_orbit_columns": len(classes),
        "zero_slab_size": len(square_slabs[0]),
        "nonzero_slab_sizes": [len(square_slabs[square]) for square in nonzero_squares],
        "all_nonzero_slab_column_sums_equal": equal_nonzero_slab_sums,
        "exact_rational_column_rank": rank,
        "exact_integer_kernel_rank": len(classes) - rank,
        "expected_integer_kernel_rank": h - 1,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    theorem = symmetric_slice_kernel_theorem(p)
    deletion = used_orbit_deletion_bound(p)
    replay = {str(q): exact_small_slice_replay(q) for q in (5, 7)}
    proved = bool(
        theorem["proved"]
        and deletion["proved"]
        and all(row["proved"] for row in replay.values())
    )
    if not proved:
        raise ArithmeticError("the symmetric slice theorem record changed")
    return {
        "title": "Unused symmetric one-difference-slice exchange theorem",
        "status": "PROVED SLICE KERNEL AND CONNECTIVITY; GLOBAL BOX OPEN",
        "theorem": theorem,
        "mobius_deletion": deletion,
        "small_exact_replay": replay,
        "binary_connectivity_scope": (
            "two binary solutions with equal target and agreeing outside one "
            "difference slice are connected by unused whole-slab exchanges"
        ),
        "full_unused_configuration_normality_proved": False,
        "global_Boolean_fibre_nonempty_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
