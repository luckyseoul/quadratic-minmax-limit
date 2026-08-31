#!/usr/bin/env python3
r"""Prop. 15.736 -- exact ``p=11`` sharp Boolean-quadratic catalog.

This proposition replaces the conditional restriction--extension input in the
``p=11`` critical residual reduction by a self-contained exhaustive finite
certificate.  It does *not* close the remaining residual branch.

Let ``Omega=J(11,6)``.  Pair-monomial evaluations ``x_i x_j`` span the
degree-at-most-two functions on ``Omega``: on the slice,

    sum_(j != i) x_i x_j = 5 x_i,
    sum_(i<j) x_i x_j = 15.

Their ``462 x 55`` evaluation matrix has rank 55 modulo 101, hence also over
the reals.  For a base 3-set and three disjoint swap pairs, the alternating
sum over the resulting eight 6-sets annihilates every quadratic.  A
deterministic greedy elimination retains 407 independent such identities.
They all annihilate the 55-dimensional evaluation space, so their real rank
is exactly ``462-55=407`` and their nullspace is exactly the quadratic
evaluation space.

The remaining classification is an exact Boolean feasibility problem.  It
has one variable ``f_X`` for every ``X in Omega``, the 407 identities,
``sum_X f_X=84``, and no-goods for the 220 known supports:

* 55 omitted-pair forms ``(1-x_i)(1-x_j)``;
* 165 all-equal triple forms
  ``1-x_i-x_j-x_k+x_i*x_j+x_i*x_k+x_j*x_k``.

CP-SAT proves the resulting integer model infeasible.  Consequently these
220 forms exhaust the Boolean quadratics on ``J(11,6)`` with sharp support
84 (density ``2/11``).  This is an **exhaustive finite certificate**, not a
symbolic classification theorem.

For the residual application, Proposition 15.688 supplies the missing bridge
to this Boolean model: at ``p=11`` a sharp nonnegative integral lift has
scaled mass 8, while every lift with maximum at least two has scaled mass at
least 12.  Hence equality forces maximum one, so the lift is Boolean; its
mass ``2/11`` on 462 slice points gives support 84.

For the residual ledger the associated signed targets are

    4-z_i-z_j+z_i*z_j                         (offset 2),
    4+z_i*z_j+z_i*z_k+z_j*z_k                 (offset 4).

Thus the hard-``b=2`` branch, whose minimum opposite parallel count is
``Q=3``, is impossible modulo ``q=5``.  In the hard-``b=p-1`` branch
``Q=4``: omitted-pair forms are impossible, but all-equal triples remain.
At least four minimum opposite directions must have that triple form.  Their
simultaneous compatibility is open here, so neither ``p=11`` nor residual
(ii) is claimed closed.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P3_LAST,
    p11_equality_obstruction,
)
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor


ROOT = Path(__file__).resolve().parents[1]
N = 11
K = 6
MODULUS = 101
DOMAIN_SIZE = 462
PAIR_COLUMN_COUNT = 55
IDENTITY_RANK = 407
SHARP_SUPPORT_SIZE = 84
OMITTED_PAIR_COUNT = 55
ALL_EQUAL_TRIPLE_COUNT = 165
KNOWN_SUPPORT_COUNT = 220

Point = tuple[int, ...]
SparseRow = tuple[tuple[int, int], ...]
CubeDescriptor = tuple[tuple[int, ...], tuple[tuple[int, int], ...]]


@lru_cache(maxsize=1)
def middle_slice_points() -> tuple[Point, ...]:
    """Return ``J(11,6)`` in deterministic lexicographic order."""
    points = tuple(combinations(range(N), K))
    if len(points) != DOMAIN_SIZE:
        raise ArithmeticError("J(11,6) domain size changed")
    return points


@lru_cache(maxsize=1)
def pair_coordinates() -> tuple[tuple[int, int], ...]:
    """The 55 pair monomials, which span degree at most two on the slice."""
    pairs = tuple(combinations(range(N), 2))
    if len(pairs) != PAIR_COLUMN_COUNT:
        raise ArithmeticError("pair-coordinate count changed")
    return pairs


def _modular_rank(rows: Iterable[Sequence[int]], modulus: int = MODULUS) -> int:
    """Compute row rank over ``F_modulus`` by deterministic elimination."""
    basis: dict[int, list[int]] = {}
    width: int | None = None
    for source in rows:
        row = [int(value) % modulus for value in source]
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise ValueError("rank rows must have one common width")
        while True:
            pivot = next((i for i, value in enumerate(row) if value), None)
            if pivot is None:
                break
            if pivot not in basis:
                inverse = pow(row[pivot], -1, modulus)
                row = [(value * inverse) % modulus for value in row]
                basis[pivot] = row
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
    """Certify that the quadratic evaluation space has exact dimension 55."""
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
    if not proved:
        raise ArithmeticError("degree-two evaluation rank changed")
    return {
        "slice": "J(11,6)",
        "point_count": len(points),
        "pair_monomial_count": len(pairs),
        "evaluation_matrix_shape": [len(points), len(pairs)],
        "rank_modulus": MODULUS,
        "rank_mod_101": rank,
        "real_rank_lower_bound_from_modular_minor": rank,
        "pair_monomials_span_degree_at_most_two": True,
        "linear_recovery_identity": "sum_(j!=i) x_i*x_j=5*x_i",
        "constant_recovery_identity": "sum_(i<j) x_i*x_j=15",
        "exact_real_dimension": PAIR_COLUMN_COUNT,
        "exact_annihilator_dimension": IDENTITY_RANK,
        "proved": proved,
    }


def _perfect_matchings(items: tuple[int, ...]) -> Iterator[tuple[tuple[int, int], ...]]:
    """Yield each perfect matching once, in lexicographic order."""
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
    """Generate base-three-set/three-swap-pair cubes deterministically."""
    universe = tuple(range(N))
    # Choosing the six swapped coordinates first distributes the fixed base
    # sets across the slice early.  This deterministic order reaches full
    # annihilator rank after 8,321 of the 69,300 possible descriptors.
    for paired_elements in combinations(universe, 6):
        complement = tuple(
            value for value in universe if value not in paired_elements
        )
        for matching in _perfect_matchings(tuple(paired_elements)):
            for base in combinations(complement, 3):
                yield (base, matching)


def _cube_row(
    descriptor: CubeDescriptor,
    point_index: dict[Point, int],
) -> SparseRow:
    """Return the signed eight-term third-difference row of one cube."""
    base, pairs = descriptor
    entries: list[tuple[int, int]] = []
    for bits in range(8):
        swapped = tuple(
            pair[(bits >> axis) & 1] for axis, pair in enumerate(pairs)
        )
        point = tuple(
            sorted(base + swapped)
        )
        sign = -1 if bits.bit_count() % 2 else 1
        entries.append((point_index[point], sign))
    row = tuple(sorted(entries))
    if len(row) != 8 or len({index for index, _sign in row}) != 8:
        raise ArithmeticError("a third-difference cube lost a vertex")
    return row


def _reduce_sparse_modular_row(
    source: SparseRow,
    basis: dict[int, dict[int, int]],
    modulus: int = MODULUS,
) -> bool:
    """Insert a sparse row into a normalized echelon basis if independent."""
    row = {index: value % modulus for index, value in source if value % modulus}
    while row:
        pivot = min(row)
        if pivot not in basis:
            inverse = pow(row[pivot], -1, modulus)
            normalized = {
                index: (value * inverse) % modulus
                for index, value in row.items()
                if (value * inverse) % modulus
            }
            basis[pivot] = normalized
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


def _identity_digest(rows: Sequence[SparseRow]) -> str:
    payload = ";".join(
        ",".join(f"{index}:{sign}" for index, sign in row)
        for row in rows
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def selected_third_difference_identities() -> tuple[
    tuple[SparseRow, ...], tuple[CubeDescriptor, ...], int
]:
    """Greedily retain exactly 407 independent cube identities modulo 101."""
    point_index = {point: index for index, point in enumerate(middle_slice_points())}
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
    if len(selected_rows) != IDENTITY_RANK:
        raise ArithmeticError("third-difference cubes did not span the annihilator")
    return tuple(selected_rows), tuple(selected_descriptors), examined


@lru_cache(maxsize=1)
def third_difference_rank_certificate() -> dict[str, object]:
    """Verify annihilation and exact real rank of the selected identities."""
    rows, descriptors, examined = selected_third_difference_identities()
    dense_rows = [
        [dict(row).get(index, 0) for index in range(DOMAIN_SIZE)]
        for row in rows
    ]
    rank = _modular_rank(dense_rows)
    all_annihilate = all(_row_annihilates_pairs(row) for row in rows)
    space = degree_two_space_certificate()
    proved = bool(
        len(rows) == IDENTITY_RANK
        and len(descriptors) == IDENTITY_RANK
        and rank == IDENTITY_RANK
        and all_annihilate
        and int(space["exact_real_dimension"]) + rank == DOMAIN_SIZE
    )
    if not proved:
        raise ArithmeticError("third-difference rank certificate failed")
    return {
        "identity_family": "base 3-set plus three disjoint swap pairs",
        "identity_term_count": 8,
        "candidate_rows_examined": examined,
        "selected_identity_count": len(rows),
        "selected_identity_sha256": _identity_digest(rows),
        "rank_modulus": MODULUS,
        "selected_rank_mod_101": rank,
        "every_identity_annihilates_all_55_pair_monomials": all_annihilate,
        "real_rank_lower_bound_from_modular_minor": rank,
        "real_rank_upper_bound_from_55_dimensional_nullspace": (
            DOMAIN_SIZE - int(space["exact_real_dimension"])
        ),
        "exact_real_rank": rank,
        "exact_real_nullity": DOMAIN_SIZE - rank,
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
def known_sharp_supports() -> tuple[
    tuple[tuple[int, ...], ...], tuple[dict[str, object], ...]
]:
    """Return the 220 known sharp supports and their form descriptors."""
    points = middle_slice_points()
    supports: list[tuple[int, ...]] = []
    forms: list[dict[str, object]] = []
    for i, j in combinations(range(N), 2):
        support = tuple(
            index for index, point in enumerate(points)
            if i not in point and j not in point
        )
        supports.append(support)
        forms.append(
            {
                "family": "omitted_pair",
                "coordinates": [i, j],
                "polynomial": "(1-x_i)(1-x_j)",
                "signed_target": "4-z_i-z_j+z_i*z_j",
                "coefficient_offset": 2,
                "support_size": len(support),
            }
        )
    for i, j, k in combinations(range(N), 3):
        support = tuple(
            index for index, point in enumerate(points)
            if sum(value in point for value in (i, j, k)) in (0, 3)
        )
        supports.append(support)
        forms.append(
            {
                "family": "all_equal_triple",
                "coordinates": [i, j, k],
                "polynomial": (
                    "1-x_i-x_j-x_k+x_i*x_j+x_i*x_k+x_j*x_k"
                ),
                "signed_target": "4+z_i*z_j+z_i*z_k+z_j*z_k",
                "coefficient_offset": 4,
                "support_size": len(support),
            }
        )
    unique_supports = set(supports)
    if (
        len(supports) != KNOWN_SUPPORT_COUNT
        or len(unique_supports) != KNOWN_SUPPORT_COUNT
        or any(len(support) != SHARP_SUPPORT_SIZE for support in supports)
    ):
        raise ArithmeticError("known sharp-support catalog changed")
    return tuple(supports), tuple(forms)


def _support_catalog_digest(supports: Sequence[tuple[int, ...]]) -> str:
    payload = ";".join(
        ",".join(str(index) for index in support) for support in supports
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def p11_sharp_lift_equality_is_boolean() -> dict[str, object]:
    """Bridge a residual equality lift to the Boolean support-84 model."""
    lift = sharp_integral_quadratic_lift_floor(11)
    equality_scaled_mass = int(lift["sharp_scaled_floor"])
    h_at_least_two_floor = int(lift["H_at_least_two_scaled_floor"])
    sharp_mass = lift["sharp_mass_floor"]
    if not isinstance(sharp_mass, Fraction):
        raise ArithmeticError("the sharp lift mass lost exact rational form")
    support_size = sharp_mass * DOMAIN_SIZE
    h_at_least_two_excluded = equality_scaled_mass < h_at_least_two_floor
    forced_maximum = 1 if h_at_least_two_excluded and sharp_mass > 0 else None
    integral_range_is_boolean = forced_maximum == 1
    proved = bool(
        equality_scaled_mass == 8
        and h_at_least_two_floor == 12
        and h_at_least_two_excluded
        and forced_maximum == 1
        and integral_range_is_boolean
        and support_size == SHARP_SUPPORT_SIZE
    )
    if not proved:
        raise ArithmeticError("p=11 equality no longer forces Boolean support 84")
    return {
        "p": 11,
        "dependency": "Proposition 15.688",
        "input_class": "nonzero nonnegative integer-valued quadratic B",
        "equality_scaled_mass_4pE_B": equality_scaled_mass,
        "H_at_least_two_scaled_floor": h_at_least_two_floor,
        "H_at_least_two_excluded": h_at_least_two_excluded,
        "forced_maximum_H": forced_maximum,
        "integer_values_between_zero_and_H_are_boolean": integral_range_is_boolean,
        "sharp_mass": str(sharp_mass),
        "slice_point_count": DOMAIN_SIZE,
        "forced_support_size": int(support_size),
        "enters_boolean_support_84_model": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def sharp_support_catalog_certificate() -> dict[str, object]:
    """Verify the size, distinctness, offsets, and identities of the catalog."""
    supports, forms = known_sharp_supports()
    rows, _descriptors, _examined = selected_third_difference_identities()
    support_sets = [set(support) for support in supports]
    all_satisfy_identities = all(
        sum(sign * int(index in support) for index, sign in row) == 0
        for support in support_sets
        for row in rows
    )
    omitted_count = sum(form["family"] == "omitted_pair" for form in forms)
    triple_count = sum(form["family"] == "all_equal_triple" for form in forms)
    offsets = {
        str(form["family"]): int(form["coefficient_offset"])
        for form in forms
    }
    pair_target_identity = all(
        3 + 4 * (1 - x_i) * (1 - x_j)
        == 4 - (2 * x_i - 1) - (2 * x_j - 1)
        + (2 * x_i - 1) * (2 * x_j - 1)
        for x_i in (0, 1)
        for x_j in (0, 1)
    )
    triple_target_identity = all(
        3
        + 4
        * (
            1
            - x_i
            - x_j
            - x_k
            + x_i * x_j
            + x_i * x_k
            + x_j * x_k
        )
        == 4
        + (2 * x_i - 1) * (2 * x_j - 1)
        + (2 * x_i - 1) * (2 * x_k - 1)
        + (2 * x_j - 1) * (2 * x_k - 1)
        for x_i in (0, 1)
        for x_j in (0, 1)
        for x_k in (0, 1)
    )
    proved = bool(
        len(supports) == KNOWN_SUPPORT_COUNT
        and len(set(supports)) == KNOWN_SUPPORT_COUNT
        and all(len(support) == SHARP_SUPPORT_SIZE for support in supports)
        and omitted_count == OMITTED_PAIR_COUNT
        and triple_count == ALL_EQUAL_TRIPLE_COUNT
        and offsets == {"omitted_pair": 2, "all_equal_triple": 4}
        and pair_target_identity
        and triple_target_identity
        and all_satisfy_identities
    )
    if not proved:
        raise ArithmeticError("sharp-support catalog verification failed")
    return {
        "support_count": len(supports),
        "distinct_support_count": len(set(supports)),
        "support_size": SHARP_SUPPORT_SIZE,
        "omitted_pair_support_count": omitted_count,
        "all_equal_triple_support_count": triple_count,
        "coefficient_offsets": offsets,
        "signed_target_normalization": "target=3+4*f",
        "omitted_pair_signed_target_verified": pair_target_identity,
        "all_equal_triple_signed_target_verified": triple_target_identity,
        "every_catalog_support_satisfies_all_407_identities": (
            all_satisfy_identities
        ),
        "known_support_catalog_sha256": _support_catalog_digest(supports),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def exact_boolean_classification() -> dict[str, object]:
    """Run the exact CP-SAT exclusion of every support outside the catalog."""
    rank = third_difference_rank_certificate()
    rows, _descriptors, _examined = selected_third_difference_identities()
    supports, _forms = known_sharp_supports()
    catalog = sharp_support_catalog_certificate()

    model = cp_model.CpModel()
    values = [model.NewBoolVar(f"f_{index}") for index in range(DOMAIN_SIZE)]
    for row in rows:
        model.Add(sum(sign * values[index] for index, sign in row) == 0)
    model.Add(sum(values) == SHARP_SUPPORT_SIZE)
    for support in supports:
        model.Add(
            sum(values[index] for index in support) <= SHARP_SUPPORT_SIZE - 1
        )

    model_proto = model.Proto()
    model_sha256 = hashlib.sha256(str(model_proto).encode("utf-8")).hexdigest()

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 32
    solver.parameters.random_seed = 0
    solver.parameters.symmetry_level = 3
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    if not infeasible:
        raise ArithmeticError(
            f"sharp Boolean classification was not certified: {status_name}"
        )

    return {
        "slice": "J(11,6)",
        "boolean_variable_count": DOMAIN_SIZE,
        "third_difference_equality_count": len(rows),
        "sharp_support_size": SHARP_SUPPORT_SIZE,
        "sharp_density": "84/462=2/11",
        "known_support_nogood_count": len(supports),
        "omitted_pair_support_count": catalog["omitted_pair_support_count"],
        "all_equal_triple_support_count": catalog["all_equal_triple_support_count"],
        "known_support_catalog_sha256": _support_catalog_digest(supports),
        "model_constraint_count": len(model_proto.constraints),
        "model_textproto_sha256": model_sha256,
        "orbit_anchor_used": False,
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "search_workers": 32,
        "interleave_search": False,
        "random_seed": 0,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "wall_time_seconds": solver.WallTime(),
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "integer_model_infeasible": infeasible,
        "catalog_exhaustive_at_support_84": infeasible,
        "catalog_certificate": catalog,
        "linear_rank_certificate": rank,
        "result_status": "exhaustive finite certificate",
        "proved": infeasible,
    }


def residual_p11_consequence() -> dict[str, object]:
    """Apply only the certified catalog information to the two p=11 branches."""
    prior = p11_equality_obstruction()
    equality_bridge = p11_sharp_lift_equality_is_boolean()
    q = 5
    direction_count = 6
    branch_a_s = int(prior["branch_reductions"][BRANCH_B2]["forced_s"])
    branch_c_s = int(prior["branch_reductions"][BRANCH_P3_LAST]["forced_s"])
    branch_a_edges = q * (8 - branch_a_s)
    branch_c_edges = q * (8 - branch_c_s) + 1
    branch_a_minimum = 3
    branch_c_minimum = 4
    branch_a_previous_mean = (
        10 * branch_a_s + 12 * (branch_a_minimum - 1) + 9 - 77
    )
    branch_c_previous_mean = (
        10 * branch_c_s + 12 * (branch_c_minimum - 1) + 7 - 77
    )
    branch_a_minimum_mean = branch_a_previous_mean + 12
    branch_c_minimum_mean = branch_c_previous_mean + 12
    branch_a_surplus = branch_a_edges - direction_count * branch_a_minimum
    branch_c_surplus = branch_c_edges - direction_count * branch_c_minimum
    branch_a_minimum_count = direction_count - branch_a_surplus
    branch_c_minimum_count = direction_count - branch_c_surplus
    offsets = {"omitted_pair": 2, "all_equal_triple": 4}
    branch_a_survivors = [
        name
        for name, offset in offsets.items()
        if (branch_a_minimum - offset) % q == 0
    ]
    branch_c_survivors = [
        name
        for name, offset in offsets.items()
        if (branch_c_minimum - offset) % q == 0
    ]
    proved_reduction = bool(
        branch_a_survivors == []
        and branch_c_survivors == ["all_equal_triple"]
        and branch_c_minimum_count >= 4
        and branch_a_previous_mean < 0
        and branch_c_previous_mean < 0
        and branch_a_minimum_mean == 8
        and branch_c_minimum_mean == 8
        and prior["small_mean_equals_lift_floor"]
        and equality_bridge["proved"]
    )
    if not proved_reduction:
        raise ArithmeticError("p=11 residual consequence changed")
    return {
        "p": 11,
        "q": q,
        "catalog_input": "exact_boolean_classification",
        "integral_lift_to_boolean_bridge": equality_bridge,
        "prior_reduction_dependency": "Proposition 15.734",
        "coefficient_comparison_dependency": (
            "Proposition 15.734 slice coefficient congruence"
        ),
        "coefficient_offsets_mod_q": offsets,
        "hard_b2_branch": {
            "forced_s": branch_a_s,
            "opposite_parallel_count_sum": branch_a_edges,
            "opposite_direction_count": direction_count,
            "minimum_Q": branch_a_minimum,
            "mean_at_Q_minus_1": branch_a_previous_mean,
            "mean_at_minimum_Q": branch_a_minimum_mean,
            "surplus": branch_a_surplus,
            "directions_at_minimum_at_least": branch_a_minimum_count,
            "catalog_forms_with_offset_congruent_to_Q": branch_a_survivors,
            "excluded": True,
        },
        "hard_b_p_minus_1_branch": {
            "forced_s": branch_c_s,
            "opposite_parallel_count_sum": branch_c_edges,
            "opposite_direction_count": direction_count,
            "minimum_Q": branch_c_minimum,
            "mean_at_Q_minus_1": branch_c_previous_mean,
            "mean_at_minimum_Q": branch_c_minimum_mean,
            "surplus": branch_c_surplus,
            "directions_at_minimum_at_least": branch_c_minimum_count,
            "catalog_forms_with_offset_congruent_to_Q": branch_c_survivors,
            "omitted_pair_excluded": True,
            "all_equal_triple_survives": True,
            "excluded": False,
        },
        "remaining_branch": (
            "simultaneous compatibility of at least four all-equal triple "
            "targets in minimum opposite directions"
        ),
        "p11_closed": False,
        "residual_ii_closed": False,
        "result_status": "open reduction",
        "proved_reduction": proved_reduction,
    }


def proposition_15736() -> dict[str, object]:
    """Package the exact finite classification and its narrow consequence."""
    space = degree_two_space_certificate()
    identities = third_difference_rank_certificate()
    classification = exact_boolean_classification()
    equality_bridge = p11_sharp_lift_equality_is_boolean()
    consequence = residual_p11_consequence()
    certified = bool(
        space["proved"]
        and identities["proved"]
        and classification["proved"]
        and equality_bridge["proved"]
        and consequence["proved_reduction"]
    )
    return {
        "prop": "15.736",
        "title": "Exact p=11 sharp Boolean-quadratic catalog",
        "result_status": "exhaustive finite certificate",
        "statement": (
            "the 55 omitted-pair and 165 all-equal-triple forms exhaust "
            "Boolean quadratics of support 84 on J(11,6)"
        ),
        "changed_premise": (
            "self-contained 407-row third-difference rank certificate "
            "replaces the unavailable external restriction-extension input"
        ),
        "degree_two_space": space,
        "third_difference_annihilator": identities,
        "exact_boolean_classification": classification,
        "p11_integral_lift_equality_bridge": equality_bridge,
        "residual_p11_consequence": consequence,
        "sharp_boolean_catalog_certified": certified,
        "hard_b2_branch_excluded_p11": certified,
        "hard_b_p_minus_1_branch_excluded_p11": False,
        "simultaneous_all_equal_triple_branch_closed": False,
        "p11_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "simultaneous all-equal triple compatibility in the p=11 hard-"
            "b=p-1 branch, plus the other recorded residual and Type-I fronts"
        ),
        "proved": certified,
    }


def write_evidence() -> Path:
    """Write the certificate summary after exact validation."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15736.json"
    payload = json.dumps(proposition_15736(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15736()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.736 certificate failed")
    path = write_evidence()
    solver = result["exact_boolean_classification"]
    print(
        "Prop. 15.736: 220 sharp Boolean quadratics exhaust J(11,6); "
        f"CP-SAT {solver['solver_status']} in {solver['wall_time_seconds']:.3f}s"
    )
    print("p=11 hard-b=2 excluded; simultaneous all-equal triples remain open")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
