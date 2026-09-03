#!/usr/bin/env python3
r"""Explicit midpoint-ridge moves in the integral edge--Radon kernel.

Write an edge of ``V=F_p^2`` as ``(a,[delta])`` where ``a`` is its
midpoint and ``[delta]`` is its nonzero half-difference modulo sign.  For a
projective functional ``L`` and a zero-sum integer function ``g:F_p->Z``
there are two elementary families.

* If ``L(delta)=0``, the parallel ridge has coefficient ``g(L(a))`` on
  ``(a,[delta])`` and zero on every other difference class.
* If ``[delta_1] != [delta_2]`` and
  ``L(delta_1)^2=L(delta_2)^2 != 0``, the transverse ridge has coefficient
  ``g(L(a))`` on ``[delta_1]`` and ``-g(L(a))`` on ``[delta_2]``.

Both families lie in ``ker_Z R``.  More strongly, if ``K_ridge`` is the
lattice they generate, then

    p ker_Z R  subset  K_ridge  subset  ker_Z R.              (1)

The proof of the first inclusion is the affine-plane Radon inversion

    p f(a) = sum_L sum_(b:L(b)=L(a)) f(b)                     (2)

for every zero-total integer function ``f`` on ``V``.  An edge--Radon
kernel vector has zero total over the midpoints separately for every
``[delta]``: after forgetting midpoints its totals lie in the kernel of the
injective pure-difference map ``S`` from Proposition 15.760.  Apply (2) to
each difference class.  In a fixed direction, the kernel equations say
that the fibre-sum functions add to zero among the ``p`` difference
classes with the same nonzero square.  The directional summand is therefore
an integral sum of the two ridge families above.

Thus the ridge moves span the rational kernel.  If ``m=(p-1)/2`` and
``d=p+1``, the exact remaining quotient is

    ker_Z R/K_ridge = (Z/pZ)^[d*p*m^2+m(m-1)(4m+1)/6].

The elementary two-fibre Type-P moves are Graver circuits for every odd
prime, and the four-fibre Type-K moves are Graver circuits for ``p>=5``.
The ridge lattice is nevertheless proper, so these are not the complete
Graver basis and do not prove a signed Boolean lift.  The executable
routines below audit the identities for supplied primes and supplied
kernel vectors; the all-prime proof is in the accompanying evidence note.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Sequence


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def _prime(value: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or not _is_prime(value):
        raise ValueError("p must be a prime")
    if value == 2:
        raise ValueError("p must be odd")
    return value


def projective_functionals(p: int) -> tuple[tuple[int, int], ...]:
    """Return ``x+s*y`` for ``s in F_p``, followed by ``y``."""
    p = _prime(p)
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def half_difference_classes(p: int) -> tuple[tuple[int, int], ...]:
    """Return one deterministic representative of each nonzero ``+/-`` class."""
    p = _prime(p)
    representatives = []
    for y in range(p):
        for x in range(p):
            if x == 0 and y == 0:
                continue
            point = (x, y)
            negative = ((-x) % p, (-y) % p)
            if point < negative:
                representatives.append(point)
    return tuple(representatives)


def _linear_value(functional: tuple[int, int], point: tuple[int, int], p: int) -> int:
    return (functional[0] * point[0] + functional[1] * point[1]) % p


def _integer_vector(values: Iterable[int], expected: int, name: str) -> tuple[int, ...]:
    out = tuple(values)
    if len(out) != expected:
        raise ValueError(f"{name} must have length {expected}")
    if any(not isinstance(value, int) or isinstance(value, bool) for value in out):
        raise ValueError(f"{name} must contain integers")
    return out


def _source_index(delta_index: int, point_index: int, p: int) -> int:
    return delta_index * p * p + point_index


def _point(point_index: int, p: int) -> tuple[int, int]:
    return point_index % p, point_index // p


def edge_radon_image(p: int, source: Sequence[int]) -> dict[tuple[object, ...], int]:
    """Return the nonzero rows of the unsigned edge--Radon image.

    Parallel keys are ``("P",L)``.  Off-diagonal keys are
    ``("K",L,alpha,beta)`` in midpoint/squared-half-difference coordinates.
    """
    p = _prime(p)
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    vector = _integer_vector(source, len(differences) * p * p, "source")
    target: defaultdict[tuple[object, ...], int] = defaultdict(int)
    for delta_index, delta in enumerate(differences):
        projected_differences = tuple(
            _linear_value(functional, delta, p) for functional in directions
        )
        for point_index in range(p * p):
            coefficient = vector[_source_index(delta_index, point_index, p)]
            if coefficient == 0:
                continue
            point = _point(point_index, p)
            for direction_index, functional in enumerate(directions):
                projected_delta = projected_differences[direction_index]
                if projected_delta == 0:
                    key: tuple[object, ...] = ("P", direction_index)
                else:
                    key = (
                        "K",
                        direction_index,
                        _linear_value(functional, point, p),
                        projected_delta * projected_delta % p,
                    )
                target[key] += coefficient
    return {key: value for key, value in target.items() if value}


def _ridge_source(
    p: int,
    direction_index: int,
    first_delta_index: int,
    g: Sequence[int],
    second_delta_index: int | None,
) -> tuple[int, ...]:
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    if direction_index < 0 or direction_index >= len(directions):
        raise ValueError("direction index is out of range")
    if first_delta_index < 0 or first_delta_index >= len(differences):
        raise ValueError("first difference index is out of range")
    if second_delta_index is not None and (
        second_delta_index < 0 or second_delta_index >= len(differences)
    ):
        raise ValueError("second difference index is out of range")
    profile = _integer_vector(g, p, "g")
    if sum(profile):
        raise ValueError("g must have sum zero")

    functional = directions[direction_index]
    source = [0] * (len(differences) * p * p)
    for point_index in range(p * p):
        point = _point(point_index, p)
        coefficient = profile[_linear_value(functional, point, p)]
        source[_source_index(first_delta_index, point_index, p)] = coefficient
        if second_delta_index is not None:
            source[_source_index(second_delta_index, point_index, p)] = -coefficient
    return tuple(source)


def parallel_ridge_move(
    p: int, direction_index: int, delta_index: int, g: Sequence[int]
) -> tuple[int, ...]:
    """Construct a Type-P ridge and prove its returned Radon image is zero."""
    p = _prime(p)
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    if direction_index < 0 or direction_index >= len(directions):
        raise ValueError("direction index is out of range")
    if delta_index < 0 or delta_index >= len(differences):
        raise ValueError("difference index is out of range")
    if _linear_value(directions[direction_index], differences[delta_index], p):
        raise ValueError("the difference class must be parallel to the direction")
    source = _ridge_source(p, direction_index, delta_index, g, None)
    if edge_radon_image(p, source):
        raise ArithmeticError("the parallel ridge left the edge-Radon kernel")
    return source


def transverse_ridge_move(
    p: int,
    direction_index: int,
    first_delta_index: int,
    second_delta_index: int,
    g: Sequence[int],
) -> tuple[int, ...]:
    """Construct a Type-K same-square ridge and prove its Radon image is zero."""
    p = _prime(p)
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    if direction_index < 0 or direction_index >= len(directions):
        raise ValueError("direction index is out of range")
    if first_delta_index == second_delta_index:
        raise ValueError("the two difference classes must be distinct")
    if min(first_delta_index, second_delta_index) < 0 or max(
        first_delta_index, second_delta_index
    ) >= len(differences):
        raise ValueError("difference index is out of range")
    functional = directions[direction_index]
    first_value = _linear_value(functional, differences[first_delta_index], p)
    second_value = _linear_value(functional, differences[second_delta_index], p)
    if first_value == 0 or first_value * first_value % p != second_value * second_value % p:
        raise ValueError("the two classes must have the same nonzero projected square")
    source = _ridge_source(
        p, direction_index, first_delta_index, g, second_delta_index
    )
    if edge_radon_image(p, source):
        raise ArithmeticError("the transverse ridge left the edge-Radon kernel")
    return source


def canonical_ridge_basis(p: int) -> tuple[tuple[int, ...], ...]:
    """Materialize the canonical Type-P/Type-K basis from equation (12).

    This is an exact audit helper, not a claim that materializing the matrix
    is practical for large primes.
    """
    p = _prime(p)
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    profiles = []
    for alpha in range(1, p):
        profile = [0] * p
        profile[0] = -1
        profile[alpha] = 1
        profiles.append(tuple(profile))

    columns = []
    for direction_index, functional in enumerate(directions):
        parallel = [
            delta_index
            for delta_index, delta in enumerate(differences)
            if _linear_value(functional, delta, p) == 0
        ]
        for delta_index in parallel:
            for profile in profiles:
                columns.append(
                    parallel_ridge_move(
                        p, direction_index, delta_index, profile
                    )
                )

        by_square: defaultdict[int, list[int]] = defaultdict(list)
        for delta_index, delta in enumerate(differences):
            value = _linear_value(functional, delta, p)
            if value:
                by_square[value * value % p].append(delta_index)
        for square in sorted(by_square):
            delta_indices = by_square[square]
            if len(delta_indices) != p:
                raise ArithmeticError("a nonzero square fibre does not have p classes")
            base = delta_indices[0]
            for delta_index in delta_indices[1:]:
                for profile in profiles:
                    columns.append(
                        transverse_ridge_move(
                            p,
                            direction_index,
                            delta_index,
                            base,
                            profile,
                        )
                    )

    expected = (p + 1) * ((p - 1) // 2) * p * (p - 1)
    if len(columns) != expected:
        raise ArithmeticError("the canonical ridge basis has the wrong size")
    return tuple(columns)


def _column_rank_mod_prime(columns: Sequence[Sequence[int]], prime: int) -> int:
    """Return the exact rank of supplied columns over the prime field."""
    prime = _prime(prime)
    if not columns:
        return 0
    width = len(columns[0])
    if any(len(column) != width for column in columns):
        raise ValueError("all columns must have the same length")
    pivots: dict[int, dict[int, int]] = {}
    for column in columns:
        vector = {
            index: int(value) % prime
            for index, value in enumerate(column)
            if int(value) % prime
        }
        while vector:
            pivot = min(vector)
            previous = pivots.get(pivot)
            if previous is None:
                inverse = pow(vector[pivot], -1, prime)
                pivots[pivot] = {
                    index: value * inverse % prime
                    for index, value in vector.items()
                }
                break
            factor = vector[pivot]
            for index, value in previous.items():
                reduced = (vector.get(index, 0) - factor * value) % prime
                if reduced:
                    vector[index] = reduced
                else:
                    vector.pop(index, None)
    return len(pivots)


def ridge_mod_p_dependency_certificate(p: int) -> dict[str, object]:
    """Materialize B_p and verify the closed dependency-nullity formula."""
    p = _prime(p)
    columns = canonical_ridge_basis(p)
    rank = _column_rank_mod_prime(columns, p)
    nullity = len(columns) - rank
    m = (p - 1) // 2
    smith_defect = m * (m - 1) * (4 * m + 1) // 6
    expected = (p + 1) * p * m * m + smith_defect
    if nullity != expected:
        raise ArithmeticError("the ridge mod-p dependency formula failed")
    return {
        "p": p,
        "edge_rows": len(columns[0]),
        "ridge_columns": len(columns),
        "rank_mod_p": rank,
        "dependency_nullity_nu_p": nullity,
        "closed_formula": expected,
        "matches_uniform_quotient_dimension": True,
        "finite_identity_check_not_target_census": True,
        "proved_arithmetic": True,
    }


def affine_radon_inversion(p: int, values: Sequence[int]) -> dict[str, object]:
    """Audit (2) for one zero-total integer function on ``F_p^2``."""
    p = _prime(p)
    function = _integer_vector(values, p * p, "values")
    if sum(function):
        raise ValueError("values must have total zero")
    directions = projective_functionals(p)
    fibre_sums = []
    for functional in directions:
        row = [0] * p
        for point_index, value in enumerate(function):
            row[_linear_value(functional, _point(point_index, p), p)] += value
        fibre_sums.append(tuple(row))
    reconstruction = tuple(
        sum(
            fibre_sums[direction_index][
                _linear_value(functional, _point(point_index, p), p)
            ]
            for direction_index, functional in enumerate(directions)
        )
        for point_index in range(p * p)
    )
    proved = reconstruction == tuple(p * value for value in function)
    if not proved:
        raise ArithmeticError("the affine Radon inversion identity failed")
    return {
        "p": p,
        "fibre_sums": [list(row) for row in fibre_sums],
        "reconstruction": list(reconstruction),
        "equals_p_times_input": True,
        "proved": True,
    }


def ridge_psaturation_certificate(p: int, source: Sequence[int]) -> dict[str, object]:
    """Decompose ``p*source`` into ridge moves for a supplied kernel vector.

    The routine verifies, rather than assumes, that the source is in the
    edge--Radon kernel and that its midpoint total is zero in every
    half-difference class.  The latter follows uniformly from injectivity of
    Proposition 15.760's pure-difference map ``S``.
    """
    p = _prime(p)
    directions = projective_functionals(p)
    differences = half_difference_classes(p)
    vector = _integer_vector(source, len(differences) * p * p, "source")
    image = edge_radon_image(p, vector)
    if image:
        raise ValueError("source must lie in the edge-Radon kernel")

    midpoint_totals = tuple(
        sum(
            vector[_source_index(delta_index, point_index, p)]
            for point_index in range(p * p)
        )
        for delta_index in range(len(differences))
    )
    if any(midpoint_totals):
        raise ArithmeticError(
            "a kernel vector survived the injective pure-difference quotient"
        )

    decomposition = [0] * len(vector)
    piece_count = 0
    for direction_index, functional in enumerate(directions):
        profiles = []
        for delta_index in range(len(differences)):
            profile = [0] * p
            for point_index in range(p * p):
                profile[
                    _linear_value(functional, _point(point_index, p), p)
                ] += vector[_source_index(delta_index, point_index, p)]
            if sum(profile):
                raise ArithmeticError("a ridge profile lost its zero-sum property")
            profiles.append(tuple(profile))

        parallel = [
            delta_index
            for delta_index, delta in enumerate(differences)
            if _linear_value(functional, delta, p) == 0
        ]
        for delta_index in parallel:
            profile = profiles[delta_index]
            if not any(profile):
                continue
            piece = parallel_ridge_move(p, direction_index, delta_index, profile)
            decomposition = [left + right for left, right in zip(decomposition, piece)]
            piece_count += 1

        by_square: defaultdict[int, list[int]] = defaultdict(list)
        for delta_index, delta in enumerate(differences):
            value = _linear_value(functional, delta, p)
            if value:
                by_square[value * value % p].append(delta_index)
        for square, delta_indices in by_square.items():
            if len(delta_indices) != p:
                raise ArithmeticError("a nonzero square fibre does not have p classes")
            for alpha in range(p):
                if sum(profiles[index][alpha] for index in delta_indices):
                    raise ArithmeticError(
                        f"the K row at square {square} did not cancel"
                    )
            base = delta_indices[0]
            for delta_index in delta_indices[1:]:
                profile = profiles[delta_index]
                if not any(profile):
                    continue
                piece = transverse_ridge_move(
                    p, direction_index, delta_index, base, profile
                )
                decomposition = [
                    left + right for left, right in zip(decomposition, piece)
                ]
                piece_count += 1

    expected = tuple(p * value for value in vector)
    proved = tuple(decomposition) == expected
    if not proved:
        raise ArithmeticError("the ridge decomposition did not reconstruct p*source")
    return {
        "p": p,
        "source_support": sum(value != 0 for value in vector),
        "ridge_piece_count": piece_count,
        "midpoint_total_zero_for_every_difference_class": True,
        "ridge_decomposition_equals_p_times_source": True,
        "p_kernel_contained_in_ridge_lattice": True,
        "proved": True,
    }


def ridge_kernel_theorem_record(p: int) -> dict[str, object]:
    """Return the exact ranks, shortest displayed supports, and live scope."""
    p = _prime(p)
    m = (p - 1) // 2
    directions = p + 1
    difference_classes = directions * m
    edges = p * p * difference_classes
    radon_rank = directions * p * m
    kernel_rank = edges - radon_rank
    parallel_basis_count = directions * m * (p - 1)
    transverse_basis_count = directions * m * (p - 1) * (p - 1)
    ridge_basis_count = parallel_basis_count + transverse_basis_count
    midpoint_smith_defect = m * (m - 1) * (4 * m + 1) // 6
    ridge_quotient_dimension = directions * p * m * m + midpoint_smith_defect
    if ridge_basis_count != kernel_rank:
        raise ArithmeticError("the ridge count failed to match the Radon nullity")
    return {
        "p": p,
        "directions": directions,
        "difference_classes": difference_classes,
        "edge_source_rank": edges,
        "edge_radon_rank": radon_rank,
        "integer_kernel_rank": kernel_rank,
        "canonical_parallel_ridges": parallel_basis_count,
        "canonical_transverse_ridges": transverse_basis_count,
        "canonical_ridge_total": ridge_basis_count,
        "midpoint_smith_defect_S0": midpoint_smith_defect,
        "ridge_quotient_dimension_nu_p": ridge_quotient_dimension,
        "ridge_quotient_order": f"{p}^{ridge_quotient_dimension}",
        "shortest_displayed_parallel_support": 2 * p,
        "shortest_displayed_parallel_squared_norm": 2 * p,
        "shortest_displayed_transverse_support": 4 * p,
        "shortest_displayed_transverse_squared_norm": 4 * p,
        "exact_one_step_saturation": (
            "ker_Z R = {v in Z^E : p*v lies in the ridge lattice}"
        ),
        "explicit_quotient_invariant": (
            "ker_Z R / K_ridge is isomorphic to ker_Fp(B_p mod p), "
            "where B_p is the canonical ridge-basis matrix"
        ),
        "closed_quotient_formula": (
            "nu_p=(p+1)*p*((p-1)/2)^2+"
            "m*(m-1)*(4*m+1)/6, with m=(p-1)/2"
        ),
        "minimum_nonridge_graver_elements": 2 * ridge_quotient_dimension,
        "exact_fibre_parameterization": (
            "z=z0+B_p*q+sum_j a_j*(B_p*c_j/p), "
            "0<=a_j<p, for a mod-p dependency basis c_j"
        ),
        "proved": {
            "ridge_moves_lie_in_integer_kernel": True,
            "elementary_parallel_ridges_are_graver": True,
            "elementary_transverse_ridges_are_graver_for_p_at_least_5": (
                p >= 5
            ),
            "p_times_integer_kernel_is_in_ridge_lattice": True,
            "ridge_moves_span_rational_kernel": True,
            "kernel_mod_ridge_is_finite_elementary_p_torsion": True,
            "ridge_quotient_dimension_formula": True,
            "ridge_lattice_is_proper": True,
            "one_p_saturation_recovers_full_integer_kernel": True,
            "mod_p_dependencies_parameterize_saturating_moves": True,
            "displayed_ridges_are_not_the_complete_graver_basis": True,
            "complete_graver_basis_has_at_least_2nu_nonridge_elements": True,
            "compact_target_has_boolean_lift": False,
            "residual_ii_closed": False,
        },
        "scope_warning": (
            "No Type-K Graver assertion is made at p=3.  The compact-target "
            "and residual-closure false entries are unproved claims, not "
            "negative mathematical statements."
        ),
    }
