#!/usr/bin/env python3
r"""Prop. 15.738 -- exact p=13 mass-14 Boolean residual catalog.

This module repairs the exceptional ``p=13,t=3`` arithmetic row left by
Proposition 15.735.  It supplies three exact ingredients.

* The 78 pair monomials have full column rank on ``J(13,7)``.  A
  deterministic family of 1,638 independent third-difference identities
  therefore cuts out exactly the degree-at-most-two evaluation space.
* Proposition 15.688 narrows a nonnegative integral quadratic of scaled mass
  14 to height one or four.  For each residual parallel count ``Q=0,6``, an
  exact coefficient/cut/row-parity/l1 CP-SAT relaxation with a safe height-4
  orbit anchor is infeasible.  Thus every residual mass-14 cell is Boolean.
* An anchored no-good model exhausts the support-462 Boolean quadratics:
  78 selected pairs, 156 oriented mixed pairs, and 858 sign-flipped
  all-equal triples.  Their coefficient offsets are respectively 6, 4, 4,
  so only selected pairs survive when ``Q=0`` or ``Q=6`` modulo 6.

For a surviving selected pair ``{i,j}``, exact coefficient normalization at
both parallel counts gives the finite-field moments

    M_2=(i-j)^2,   M_4=(i-j)^4  in F_13.

This is an exhaustive finite certificate and a residual cell
classification.  The cross-direction moment contradiction is deliberately
left to the next proposition.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor


ROOT = Path(__file__).resolve().parents[1]
P = 13
N = 13
K = 7
M = 7
Q_MODULUS = 6
RANK_MODULUS = 101
DOMAIN_SIZE = 1716
PAIR_COLUMN_COUNT = 78
IDENTITY_RANK = 1638
MASS14_SUPPORT_SIZE = 462
SELECTED_PAIR_COUNT = 78
ORIENTED_MIXED_PAIR_COUNT = 156
MIXED_SIGNED_TRIPLE_COUNT = 858
KNOWN_SUPPORT_COUNT = 1092
ANCHORED_KNOWN_SUPPORT_COUNT = 294
RESIDUAL_PARALLEL_COUNTS = (0, 6)
RESIDUAL_EDGE_COUNT = 59
SEARCH_WORKERS = 32

Point = tuple[int, ...]
SparseRow = tuple[tuple[int, int], ...]
CubeDescriptor = tuple[tuple[int, ...], tuple[tuple[int, int], ...]]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


@lru_cache(maxsize=1)
def middle_slice_points() -> tuple[Point, ...]:
    """Return ``J(13,7)`` in deterministic lexicographic order."""
    points = tuple(combinations(range(N), K))
    _require(len(points) == DOMAIN_SIZE, "J(13,7) domain size changed")
    return points


@lru_cache(maxsize=1)
def pair_coordinates() -> tuple[tuple[int, int], ...]:
    """Return the 78 pair monomials spanning degree at most two."""
    pairs = tuple(combinations(range(N), 2))
    _require(len(pairs) == PAIR_COLUMN_COUNT, "pair count changed")
    return pairs


def _modular_rank(
    rows: Iterable[Sequence[int]], modulus: int = RANK_MODULUS
) -> int:
    """Compute deterministic row rank over a prime field."""
    basis: dict[int, list[int]] = {}
    width: int | None = None
    for source in rows:
        row = [int(value) % modulus for value in source]
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise ValueError("rank rows must have one width")
        while True:
            pivot = next((index for index, value in enumerate(row) if value), None)
            if pivot is None:
                break
            if pivot not in basis:
                inverse = pow(row[pivot], -1, modulus)
                basis[pivot] = [(value * inverse) % modulus for value in row]
                break
            factor = row[pivot]
            old = basis[pivot]
            row = [
                (value - factor * old_value) % modulus
                for value, old_value in zip(row, old)
            ]
    return len(basis)


@lru_cache(maxsize=1)
def degree_two_space_certificate() -> dict[str, object]:
    """Certify exact dimension 78 of the quadratic evaluation space."""
    points = middle_slice_points()
    pairs = pair_coordinates()
    rows = [
        [int(i in point and j in point) for i, j in pairs]
        for point in points
    ]
    rank = _modular_rank(rows)
    proved = bool(
        rank == PAIR_COLUMN_COUNT
        and DOMAIN_SIZE - rank == IDENTITY_RANK
    )
    _require(proved, "degree-two evaluation rank changed")
    return {
        "slice": "J(13,7)",
        "point_count": len(points),
        "pair_monomial_count": len(pairs),
        "evaluation_matrix_shape": [len(points), len(pairs)],
        "rank_modulus": RANK_MODULUS,
        "rank_mod_101": rank,
        "real_rank_lower_bound_from_modular_minor": rank,
        "real_rank_upper_bound_from_column_count": PAIR_COLUMN_COUNT,
        "linear_recovery_identity": "sum_(j!=i) x_i*x_j=6*x_i",
        "constant_recovery_identity": "sum_(i<j) x_i*x_j=21",
        "pair_monomials_span_degree_at_most_two": True,
        "exact_real_dimension": PAIR_COLUMN_COUNT,
        "exact_annihilator_dimension": IDENTITY_RANK,
        "proved": proved,
    }


def _perfect_matchings(
    items: tuple[int, ...],
) -> Iterator[tuple[tuple[int, int], ...]]:
    if not items:
        yield ()
        return
    first = items[0]
    for position in range(1, len(items)):
        second = items[position]
        remaining = items[1:position] + items[position + 1 :]
        for tail in _perfect_matchings(remaining):
            yield ((first, second),) + tail


def _cube_descriptors() -> Iterator[CubeDescriptor]:
    universe = tuple(range(N))
    for paired_elements in combinations(universe, 6):
        complement = tuple(
            value for value in universe if value not in paired_elements
        )
        for matching in _perfect_matchings(tuple(paired_elements)):
            for base in combinations(complement, K - 3):
                yield base, matching


def _cube_row(
    descriptor: CubeDescriptor, point_index: dict[Point, int]
) -> SparseRow:
    base, matching = descriptor
    entries: list[tuple[int, int]] = []
    for bits in range(8):
        selected = tuple(
            pair[(bits >> axis) & 1]
            for axis, pair in enumerate(matching)
        )
        point = tuple(sorted(base + selected))
        sign = -1 if bits.bit_count() % 2 else 1
        entries.append((point_index[point], sign))
    row = tuple(sorted(entries))
    _require(
        len(row) == 8 and len({index for index, _sign in row}) == 8,
        "third-difference cube lost a vertex",
    )
    return row


def _reduce_sparse_modular_row(
    source: SparseRow,
    basis: dict[int, dict[int, int]],
    modulus: int = RANK_MODULUS,
) -> bool:
    row = {index: value % modulus for index, value in source if value % modulus}
    while row:
        pivot = min(row)
        if pivot not in basis:
            inverse = pow(row[pivot], -1, modulus)
            basis[pivot] = {
                index: (value * inverse) % modulus
                for index, value in row.items()
                if (value * inverse) % modulus
            }
            return True
        factor = row[pivot]
        for index, value in basis[pivot].items():
            updated = (row.get(index, 0) - factor * value) % modulus
            if updated:
                row[index] = updated
            else:
                row.pop(index, None)
    return False


def _row_annihilates_pairs(row: SparseRow) -> bool:
    points = middle_slice_points()
    return all(
        sum(
            sign * int(i in points[index] and j in points[index])
            for index, sign in row
        )
        == 0
        for i, j in pair_coordinates()
    )


def _sparse_row_digest(rows: Sequence[SparseRow]) -> str:
    payload = ";".join(
        ",".join(f"{index}:{sign}" for index, sign in row)
        for row in rows
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def selected_third_difference_identities() -> tuple[
    tuple[SparseRow, ...], tuple[CubeDescriptor, ...], int
]:
    """Retain 1,638 independent third differences modulo 101."""
    point_index = {
        point: index for index, point in enumerate(middle_slice_points())
    }
    basis: dict[int, dict[int, int]] = {}
    selected_rows: list[SparseRow] = []
    selected_descriptors: list[CubeDescriptor] = []
    examined = 0
    for descriptor in _cube_descriptors():
        examined += 1
        row = _cube_row(descriptor, point_index)
        if _reduce_sparse_modular_row(row, basis):
            selected_rows.append(row)
            selected_descriptors.append(descriptor)
            if len(selected_rows) == IDENTITY_RANK:
                break
    _require(
        len(selected_rows) == IDENTITY_RANK,
        "third differences did not span the annihilator",
    )
    return tuple(selected_rows), tuple(selected_descriptors), examined


@lru_cache(maxsize=1)
def third_difference_rank_certificate() -> dict[str, object]:
    """Certify the exact real third-difference annihilator."""
    rows, descriptors, examined = selected_third_difference_identities()
    all_annihilate = all(_row_annihilates_pairs(row) for row in rows)
    space = degree_two_space_certificate()
    candidate_count = comb(N, 6) * 15 * comb(N - 6, K - 3)
    proved = bool(
        len(rows) == IDENTITY_RANK
        and all_annihilate
        and int(space["exact_real_dimension"]) + len(rows) == DOMAIN_SIZE
    )
    _require(proved, "third-difference rank certificate failed")
    return {
        "identity_family": "base 4-set plus three disjoint swap pairs",
        "identity_term_count": 8,
        "candidate_descriptor_count": candidate_count,
        "candidate_rows_examined": examined,
        "selected_identity_count": len(rows),
        "selected_identity_sha256": _sparse_row_digest(rows),
        "rank_modulus": RANK_MODULUS,
        "selected_rank_mod_101": len(rows),
        "every_identity_annihilates_all_78_pair_monomials": all_annihilate,
        "real_rank_lower_bound_from_modular_independence": len(rows),
        "real_rank_upper_bound_from_78_dimensional_nullspace": IDENTITY_RANK,
        "exact_real_rank": IDENTITY_RANK,
        "exact_real_nullity": PAIR_COLUMN_COUNT,
        "nullspace_equals_degree_at_most_two_evaluation_space": True,
        "first_selected_cube": {
            "base": list(descriptors[0][0]),
            "pairs": [list(pair) for pair in descriptors[0][1]],
        },
        "last_selected_cube": {
            "base": list(descriptors[-1][0]),
            "pairs": [list(pair) for pair in descriptors[-1][1]],
        },
        "proved": proved,
    }


@lru_cache(maxsize=1)
def mass14_height_dichotomy() -> dict[str, object]:
    """Use Proposition 15.688 to reduce height to one or four."""
    lift = sharp_integral_quadratic_lift_floor(P)
    scaled_mass = 14
    stabilizer_coefficient = Fraction(
        lift["H_at_least_two_stabilizer_coefficient"]
    )
    candidates = [
        height
        for height in range(2, 1 + scaled_mass)
        if scaled_mass >= 2 * (P + 1) - 4 * height
        and scaled_mass >= stabilizer_coefficient * height
    ]
    proved = bool(
        stabilizer_coefficient == 3
        and candidates == [4]
        and int(lift["sharp_scaled_floor"]) == 10
    )
    _require(proved, "mass-14 height dichotomy changed")
    return {
        "p": P,
        "dependency": "Proposition 15.688",
        "scaled_mass_4pE_B": scaled_mass,
        "height_one_is_boolean": True,
        "H_at_least_two_paired_bound": "14>=28-4H",
        "H_at_least_two_stabilizer_bound": "14>=3H",
        "H_at_least_two_integer_candidates": candidates,
        "height_dichotomy": [1, 4],
        "only_nonboolean_case_to_exclude": 4,
        "proved": proved,
    }


def _model_textproto_sha256(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


@lru_cache(maxsize=None)
def residual_height_four_exclusion(parallel_count: int) -> dict[str, object]:
    """Exclude a residual height-four mass-14 cell at ``Q=0`` or ``6``.

    Put ``W_st=epsilon*K_st``.  For ``b=0``, ``A=2B`` and coefficient
    evaluation on a slice point ``X`` gives

        4 B(X) = Q-3 + sum(W) - 2 cut_W(X).

    Averaging fixes ``sum(W)=13Q-53``.  A selected transverse edge can
    contribute only one unit to one signed cell, so ``sum |W|<=59-Q``.
    Empty odd-fibre boundary gives even row sums.  These are necessary
    conditions, hence infeasibility of this relaxation is an exclusion.
    """
    if parallel_count not in RESIDUAL_PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 0 or 6")
    points = middle_slice_points()
    pairs = pair_coordinates()
    total_w = P * parallel_count - 53
    l1_budget = RESIDUAL_EDGE_COUNT - parallel_count

    model = cp_model.CpModel()
    values = [model.NewIntVar(0, 4, f"B_{index}") for index in range(DOMAIN_SIZE)]
    weights = {
        pair: model.NewIntVar(-l1_budget, l1_budget, f"W_{pair[0]}_{pair[1]}")
        for pair in pairs
    }
    absolute_weights = {
        pair: model.NewIntVar(0, l1_budget, f"absW_{pair[0]}_{pair[1]}")
        for pair in pairs
    }

    model.Add(sum(values) == MASS14_SUPPORT_SIZE)
    # If the only non-Boolean height is four, transitivity on J(13,7)
    # permits a maximum point to be moved to the first lexicographic set.
    model.Add(values[0] == 4)
    model.Add(sum(weights.values()) == total_w)
    for pair in pairs:
        model.AddAbsEquality(absolute_weights[pair], weights[pair])
    model.Add(sum(absolute_weights.values()) <= l1_budget)
    row_halves = [
        model.NewIntVar(-l1_budget, l1_budget, f"row_half_{s}")
        for s in range(P)
    ]
    for s in range(P):
        model.Add(
            sum(value for pair, value in weights.items() if s in pair)
            == 2 * row_halves[s]
        )
    for index, point in enumerate(points):
        point_set = set(point)
        cut = sum(
            value
            for (s, t), value in weights.items()
            if (s in point_set) != (t in point_set)
        )
        model.Add(
            4 * values[index]
            == parallel_count - 3 + total_w - 2 * cut
        )

    validation = model.Validate()
    _require(not validation, f"invalid residual model: {validation}")
    model_hash = _model_textproto_sha256(model)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = SEARCH_WORKERS
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    _require(infeasible, f"Q={parallel_count} height-four model is {status_name}")
    return {
        "p": P,
        "parallel_count_Q": parallel_count,
        "residual_edge_count": RESIDUAL_EDGE_COUNT,
        "transverse_edge_count": l1_budget,
        "scaled_mass_4pE_B": 14,
        "value_sum": MASS14_SUPPORT_SIZE,
        "height_four_orbit_anchor": "B(first lexicographic 7-set)=4",
        "height_four_orbit_anchor_is_wlog": True,
        "coefficient_variable_count": len(weights),
        "coefficient_sum": total_w,
        "coefficient_sum_formula": "sum W=13Q-53",
        "cut_identity": "4B(X)=Q-3+sum(W)-2*cut_W(X)",
        "every_coefficient_row_sum_even": True,
        "l1_budget": l1_budget,
        "l1_inequality": "sum |W_st|<=59-Q",
        "relaxation_is_necessary_for_a_residual_cell": True,
        "integer_variable_count": len(model.Proto().variables),
        "model_constraint_count": len(model.Proto().constraints),
        "model_textproto_sha256": model_hash,
        "model_validation": validation,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "search_workers": SEARCH_WORKERS,
        "random_seed": 0,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "solver_wall_time_seconds": solver.WallTime(),
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "height_four_model_infeasible": infeasible,
        "result_status": "exhaustive finite certificate",
        "proved": infeasible,
    }


@lru_cache(maxsize=1)
def known_mass14_boolean_supports() -> tuple[
    tuple[tuple[int, ...], ...], tuple[dict[str, object], ...]
]:
    """Return the 1,092 candidate support-462 Boolean quadratics."""
    points = middle_slice_points()
    supports: list[tuple[int, ...]] = []
    forms: list[dict[str, object]] = []

    for i, j in combinations(range(N), 2):
        support = tuple(
            index
            for index, point in enumerate(points)
            if i in point and j in point
        )
        supports.append(support)
        forms.append(
            {
                "family": "selected_pair",
                "coordinates": [i, j],
                "polynomial": "x_i*x_j",
                "signed_target": "3+4*f=4+z_i+z_j+z_i*z_j",
                "coefficient_offset": 6,
                "support_size": len(support),
            }
        )

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            support = tuple(
                index
                for index, point in enumerate(points)
                if i in point and j not in point
            )
            supports.append(support)
            forms.append(
                {
                    "family": "oriented_mixed_pair",
                    "coordinates": [i, j],
                    "polynomial": "x_i*(1-x_j)",
                    "signed_target": "3+4*f=4+z_i-z_j-z_i*z_j",
                    "coefficient_offset": 4,
                    "support_size": len(support),
                }
            )

    for triple in combinations(range(N), 3):
        for opposite in triple:
            i, j = tuple(value for value in triple if value != opposite)
            k = opposite
            support = tuple(
                index
                for index, point in enumerate(points)
                if (i in point) == (j in point) != (k in point)
            )
            supports.append(support)
            forms.append(
                {
                    "family": "mixed_all_equal_signed_triple",
                    "coordinates": [i, j, k],
                    "opposite_coordinate": k,
                    "polynomial": "x_i*x_j+x_k-x_i*x_k-x_j*x_k",
                    "signed_target": "3+4*f=4+z_i*z_j-z_i*z_k-z_j*z_k",
                    "coefficient_offset": 4,
                    "support_size": len(support),
                }
            )

    _require(len(supports) == KNOWN_SUPPORT_COUNT, "catalog count changed")
    _require(
        len(set(supports)) == KNOWN_SUPPORT_COUNT,
        "catalog supports are not distinct",
    )
    _require(
        all(len(support) == MASS14_SUPPORT_SIZE for support in supports),
        "catalog support size changed",
    )
    return tuple(supports), tuple(forms)


def _support_digest(supports: Sequence[tuple[int, ...]]) -> str:
    payload = ";".join(",".join(map(str, support)) for support in supports)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def mass14_boolean_catalog_certificate() -> dict[str, object]:
    """Check candidate counts, targets, offsets, and quadratic identities."""
    supports, forms = known_mass14_boolean_supports()
    rows, _descriptors, _examined = selected_third_difference_identities()
    support_sets = tuple(map(set, supports))
    all_satisfy_identities = all(
        sum(sign * int(index in support) for index, sign in row) == 0
        for support in support_sets
        for row in rows
    )
    family_counts = dict(Counter(str(form["family"]) for form in forms))
    offsets = {
        family: sorted(
            {
                int(form["coefficient_offset"])
                for form in forms
                if form["family"] == family
            }
        )
        for family in family_counts
    }
    target_identities = {
        "selected_pair": all(
            3 + 4 * xi * xj
            == 4
            + (2 * xi - 1)
            + (2 * xj - 1)
            + (2 * xi - 1) * (2 * xj - 1)
            for xi in (0, 1)
            for xj in (0, 1)
        ),
        "oriented_mixed_pair": all(
            3 + 4 * xi * (1 - xj)
            == 4
            + (2 * xi - 1)
            - (2 * xj - 1)
            - (2 * xi - 1) * (2 * xj - 1)
            for xi in (0, 1)
            for xj in (0, 1)
        ),
        "mixed_all_equal_signed_triple": all(
            3 + 4 * (xi * xj + xk - xi * xk - xj * xk)
            == 4
            + (2 * xi - 1) * (2 * xj - 1)
            - (2 * xi - 1) * (2 * xk - 1)
            - (2 * xj - 1) * (2 * xk - 1)
            for xi in (0, 1)
            for xj in (0, 1)
            for xk in (0, 1)
        ),
    }
    survivors = {
        parallel_count: sorted(
            family
            for family, family_offsets in offsets.items()
            if (parallel_count - family_offsets[0]) % Q_MODULUS == 0
        )
        for parallel_count in RESIDUAL_PARALLEL_COUNTS
    }
    anchored_indices = [
        index for index, support in enumerate(supports) if 0 in support
    ]
    anchored_family_counts = dict(
        Counter(str(forms[index]["family"]) for index in anchored_indices)
    )
    expected_counts = {
        "selected_pair": SELECTED_PAIR_COUNT,
        "oriented_mixed_pair": ORIENTED_MIXED_PAIR_COUNT,
        "mixed_all_equal_signed_triple": MIXED_SIGNED_TRIPLE_COUNT,
    }
    expected_anchored_counts = {
        "selected_pair": 21,
        "oriented_mixed_pair": 42,
        "mixed_all_equal_signed_triple": 231,
    }
    proved = bool(
        family_counts == expected_counts
        and anchored_family_counts == expected_anchored_counts
        and len(anchored_indices) == ANCHORED_KNOWN_SUPPORT_COUNT
        and offsets
        == {
            "selected_pair": [6],
            "oriented_mixed_pair": [4],
            "mixed_all_equal_signed_triple": [4],
        }
        and survivors == {0: ["selected_pair"], 6: ["selected_pair"]}
        and all(target_identities.values())
        and all_satisfy_identities
    )
    _require(proved, "mass-14 Boolean catalog verification failed")
    return {
        "support_count": len(supports),
        "distinct_support_count": len(set(supports)),
        "support_size": MASS14_SUPPORT_SIZE,
        "density": "462/1716=7/26",
        "family_counts": family_counts,
        "coefficient_offsets": offsets,
        "coefficient_offset_definition": "target_constant+sum(target_linear)",
        "families_surviving_offset_mod_6": {
            str(key): value for key, value in survivors.items()
        },
        "signed_target_normalization": "target=3+4*f",
        "signed_target_identities_verified": target_identities,
        "every_catalog_support_satisfies_all_1638_identities": (
            all_satisfy_identities
        ),
        "known_support_catalog_sha256": _support_digest(supports),
        "support_point_anchor_index": 0,
        "anchored_catalog_support_count": len(anchored_indices),
        "anchored_family_counts": anchored_family_counts,
        "anchored_support_catalog_sha256": _support_digest(
            [supports[index] for index in anchored_indices]
        ),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def exact_mass14_boolean_classification() -> dict[str, object]:
    """Exclude every support-462 Boolean quadratic outside the catalog.

    A support point can be relabelled to the first point because ``S_13`` is
    transitive on ``J(13,7)``.  The catalog is itself invariant under this
    action.  After imposing ``f_0=1``, only catalog supports containing point
    zero can possibly equal the unknown support, so only those no-goods are
    needed.
    """
    rank = third_difference_rank_certificate()
    rows, _descriptors, _examined = selected_third_difference_identities()
    supports, forms = known_mass14_boolean_supports()
    catalog = mass14_boolean_catalog_certificate()
    anchored = [
        (support, form)
        for support, form in zip(supports, forms)
        if 0 in support
    ]

    model = cp_model.CpModel()
    values = [model.NewBoolVar(f"f_{index}") for index in range(DOMAIN_SIZE)]
    for row in rows:
        model.Add(sum(sign * values[index] for index, sign in row) == 0)
    model.Add(sum(values) == MASS14_SUPPORT_SIZE)
    model.Add(values[0] == 1)
    for support, _form in anchored:
        model.Add(
            sum(values[index] for index in support)
            <= MASS14_SUPPORT_SIZE - 1
        )

    validation = model.Validate()
    _require(not validation, f"invalid Boolean catalog model: {validation}")
    model_hash = _model_textproto_sha256(model)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = SEARCH_WORKERS
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    _require(infeasible, f"Boolean classification status is {status_name}")
    return {
        "slice": "J(13,7)",
        "boolean_variable_count": DOMAIN_SIZE,
        "third_difference_equality_count": len(rows),
        "support_size": MASS14_SUPPORT_SIZE,
        "support_point_orbit_anchor": "f(first lexicographic 7-set)=1",
        "support_point_orbit_anchor_is_wlog": True,
        "catalog_is_invariant_under_S13": True,
        "full_known_support_count": len(supports),
        "anchored_known_support_nogood_count": len(anchored),
        "anchored_nogood_filter_is_exact": True,
        "known_support_catalog_sha256": _support_digest(supports),
        "anchored_support_catalog_sha256": catalog[
            "anchored_support_catalog_sha256"
        ],
        "model_constraint_count": len(model.Proto().constraints),
        "model_textproto_sha256": model_hash,
        "model_validation": validation,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "search_workers": SEARCH_WORKERS,
        "random_seed": 0,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "solver_wall_time_seconds": solver.WallTime(),
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "integer_model_infeasible": infeasible,
        "catalog_exhaustive_at_support_462": infeasible,
        "catalog_certificate": catalog,
        "linear_rank_certificate": rank,
        "result_status": "exhaustive finite certificate",
        "proved": infeasible,
    }


def _normalized_selected_pair_pattern(
    i: int, j: int, parallel_count: int
) -> dict[tuple[int, int], int]:
    if not 0 <= i < j < P:
        raise ValueError("need 0<=i<j<13")
    if parallel_count not in RESIDUAL_PARALLEL_COUNTS:
        raise ValueError("parallel_count must be 0 or 6")
    # The target 4+z_i+z_j+z_i*z_j has offset 6.  The slice-kernel
    # scalar c=(Q-6)/12 may be half-integral, while all normalized pair
    # coefficients target_pair+2c+linear_s+linear_t are integral.
    kernel_scalar = Fraction(parallel_count - 6, 12)
    pattern: dict[tuple[int, int], int] = {}
    for s, t in pair_coordinates():
        value = (
            int((s, t) == (i, j))
            + 2 * kernel_scalar
            + int(s in (i, j))
            + int(t in (i, j))
        )
        _require(value.denominator == 1, "normalized coefficient not integral")
        pattern[(s, t)] = int(value)
    return pattern


@lru_cache(maxsize=1)
def selected_pair_moment_certificate() -> dict[str, object]:
    """Certify the degree-two and degree-four moments at ``Q=0,6``."""
    checks: list[dict[str, object]] = []
    expected_histograms = {
        0: {-1: 55, 0: 22, 2: 1},
        6: {0: 55, 1: 22, 3: 1},
    }
    for parallel_count in RESIDUAL_PARALLEL_COUNTS:
        reference_pattern = _normalized_selected_pair_pattern(
            0, 1, parallel_count
        )
        histogram = dict(sorted(Counter(reference_pattern.values()).items()))
        _require(
            histogram == expected_histograms[parallel_count],
            "selected-pair coefficient histogram changed",
        )
        for degree in (2, 4):
            values: list[int] = []
            for i, j in pair_coordinates():
                pattern = _normalized_selected_pair_pattern(
                    i, j, parallel_count
                )
                actual = sum(
                    coefficient * pow(s - t, degree, P)
                    for (s, t), coefficient in pattern.items()
                ) % P
                expected = pow(i - j, degree, P)
                _require(actual == expected, "selected-pair moment changed")
                values.append(actual)
            checks.append(
                {
                    "parallel_count_Q": parallel_count,
                    "slice_kernel_scalar": str(
                        Fraction(parallel_count - 6, 12)
                    ),
                    "normalized_coefficient_histogram": histogram,
                    "degree": degree,
                    "pair_count_checked": len(values),
                    "zero_moment_count": values.count(0),
                    "normalized_moment_formula": (
                        f"M_{degree}(i,j)=(i-j)^{degree} mod 13"
                    ),
                }
            )

    complete_graph_moments = {
        degree: sum(
            pow(s - t, degree, P) for s, t in pair_coordinates()
        )
        % P
        for degree in (2, 4)
    }
    q0 = _normalized_selected_pair_pattern(0, 1, 0)
    q6 = _normalized_selected_pair_pattern(0, 1, 6)
    patterns_differ_by_complete_graph = all(
        q6[pair] - q0[pair] == 1 for pair in pair_coordinates()
    )
    proved = bool(
        complete_graph_moments == {2: 0, 4: 0}
        and patterns_differ_by_complete_graph
        and all(int(check["zero_moment_count"]) == 0 for check in checks)
    )
    _require(proved, "selected-pair moment certificate failed")
    return {
        "field": "F_13",
        "parallel_counts": list(RESIDUAL_PARALLEL_COUNTS),
        "degrees": [2, 4],
        "complete_graph_even_moments": complete_graph_moments,
        "Q0_and_Q6_patterns_differ_by_complete_graph": (
            patterns_differ_by_complete_graph
        ),
        "checks": checks,
        "every_selected_pair_second_moment_nonzero": True,
        "every_selected_pair_fourth_moment_nonzero": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p13_mass14_residual_cell_classification() -> dict[str, object]:
    """Classify mass-14 residual cells at the two relevant Q values."""
    height = mass14_height_dichotomy()
    catalog = exact_mass14_boolean_classification()
    moments = selected_pair_moment_certificate()
    offsets = mass14_boolean_catalog_certificate()[
        "families_surviving_offset_mod_6"
    ]
    q_rows: dict[str, dict[str, object]] = {}
    for parallel_count in RESIDUAL_PARALLEL_COUNTS:
        exclusion = residual_height_four_exclusion(parallel_count)
        survivors = list(offsets[str(parallel_count)])
        q_rows[str(parallel_count)] = {
            "parallel_count_Q": parallel_count,
            "height_four_exclusion": exclusion,
            "height_four_excluded": bool(exclusion["proved"]),
            "mass14_cell_forced_boolean": bool(exclusion["proved"]),
            "catalog_survivors_after_offset_mod_6": survivors,
            "selected_pair_is_unique_surviving_family": (
                survivors == ["selected_pair"]
            ),
            "selected_pair_signed_target": (
                "epsilon*S_H=4+z_i+z_j+z_i*z_j"
            ),
            "proved": bool(
                exclusion["proved"] and survivors == ["selected_pair"]
            ),
        }
    proved = bool(
        height["proved"]
        and catalog["proved"]
        and moments["proved"]
        and all(row["proved"] for row in q_rows.values())
    )
    _require(proved, "p=13 mass-14 residual classification failed")
    return {
        "p": P,
        "slice": "J(13,7)",
        "scaled_mass_4pE_B": 14,
        "height_dichotomy": height,
        "boolean_catalog_dependency": {
            "result_status": catalog["result_status"],
            "catalog_exhaustive_at_support_462": catalog[
                "catalog_exhaustive_at_support_462"
            ],
        },
        "Q0": q_rows["0"],
        "Q6": q_rows["6"],
        "Q0_survivors": q_rows["0"][
            "catalog_survivors_after_offset_mod_6"
        ],
        "Q6_survivors": q_rows["6"][
            "catalog_survivors_after_offset_mod_6"
        ],
        "selected_pair_moment_dependency": moments,
        "result_status": "exhaustive finite certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15738() -> dict[str, object]:
    """Package the exact p=13 mass-14 residual cell classification."""
    space = degree_two_space_certificate()
    identities = third_difference_rank_certificate()
    catalog = mass14_boolean_catalog_certificate()
    classification = exact_mass14_boolean_classification()
    residual = p13_mass14_residual_cell_classification()
    moments = selected_pair_moment_certificate()
    proved = bool(
        space["proved"]
        and identities["proved"]
        and catalog["proved"]
        and classification["proved"]
        and residual["proved"]
        and moments["proved"]
    )
    _require(proved, "Proposition 15.738 certificate failed")
    return {
        "prop": "15.738",
        "title": "Exact p=13 mass-14 Boolean residual catalog",
        "result_status": "exhaustive finite certificate",
        "statement": (
            "at Q=0 or Q=6, every p=13 residual mass-14 cell is Boolean "
            "and its unique catalog family is x_i*x_j"
        ),
        "changed_premise": (
            "the corrected p=13,t=3 exceptional ledger has hard parallel "
            "counts P=2 or 8 and minimum opposite counts Q=6 or 0"
        ),
        "degree_two_space": space,
        "third_difference_annihilator": identities,
        "mass14_height_dichotomy": mass14_height_dichotomy(),
        "mass14_boolean_catalog": catalog,
        "exact_boolean_classification": classification,
        "residual_cell_classification": residual,
        "selected_pair_moments": moments,
        "p13_mass14_cells_classified": True,
        "p13_t3_exceptional_branch_closed_here": False,
        "cross_direction_moment_step_deferred_to": "Proposition 15.739",
        "residual_ii_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "proved": proved,
    }


def write_evidence() -> Path:
    output = ROOT / "evidence" / "e1_gmin_m4_prop15738.json"
    output.write_text(
        json.dumps(proposition_15738(), indent=2, sort_keys=True) + "\n"
    )
    return output


def main() -> None:
    result = proposition_15738()
    path = write_evidence()
    classification = result["exact_boolean_classification"]
    residual = result["residual_cell_classification"]
    print(
        "Prop. 15.738: 1,092 support-462 Boolean quadratics exhaust "
        "J(13,7); only x_i*x_j survives at Q=0,6"
    )
    print(
        f"catalog CP-SAT {classification['solver_status']} in "
        f"{classification['solver_wall_time_seconds']:.3f}s; "
        f"Q=0/Q=6 height-four models "
        f"{residual['Q0']['height_four_exclusion']['solver_status']}/"
        f"{residual['Q6']['height_four_exclusion']['solver_status']}"
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
