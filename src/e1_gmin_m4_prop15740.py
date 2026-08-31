#!/usr/bin/env python3
r"""Prop. 15.740 -- translation averages split the generic p=13 branch.

In the generic ``p=13,t=3`` branch the seven hard quotient units satisfy

    k_L >= 1,              sum_L k_L = 10.

Up to permutation their only patterns are

    1^6 4,       1^5 2 3,       1^4 2^3.

An exact ``k_L=1`` hard cell is a signed star.  Five such directions force
both global binary forms ``M_2`` and ``M_4`` to vanish identically.  This
module gives an exact finite certificate that no opposite ``Q=3,b=0`` cell
can obey those two moment congruences.

Write ``W`` for its normalized signed coefficient matrix.  The exact local
conditions imply

    sum W = -20,       sum |W| <= 56,
    cut_W(X) <= -10                    (|X|=7),

and conditioning on a pair contained in ``X`` gives ``W_ij >= -1``.  For
``a=1,...,6``, aggregate the thirteen cyclic distance-``+-a`` entries as

    n_a = sum_i W_{i,i+a}.

Then ``-13<=n_a<=18``, ``sum n_a=-20``, ``sum |n_a|<=56``, and

    sum a^2 n_a = sum a^4 n_a = 0 mod 13.

For a seven-set ``X`` put ``c_a(X)=|X triangle (X+a)|``.  Summing the cut
inequality over all thirteen translates of ``X`` gives

    sum_a c_a(X)n_a <= -130.

The 1,716 seven-sets give 74 distinct vectors ``c(X)``.  Exact enumeration
leaves 32,313 aggregate vectors after the sum, l1, and moment conditions;
nine deterministic translated-cut vectors eliminate all of them.  A second
six-variable CP-SAT model independently returns INFEASIBLE.

Therefore the first two quotient patterns above are impossible.  The only
remaining generic p=13 pattern is ``1^4 2^3``: four exact stars and three
elevated hard cells.  This is a proved branch split, not a close of that
last pattern or of residual (ii).

The binary affine-Radon identity used to reconstruct boundary words is
Proposition 15.692 and is imported here; it is not renumbered or reproved as
a new proposition.
"""
from __future__ import annotations

import hashlib
import json
from bisect import bisect_right
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb
from pathlib import Path
from typing import Iterable

from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model

from e1_gmin_m4_prop15692 import affine_binary_radon_isomorphism


ROOT = Path(__file__).resolve().parents[1]
P = 13
M = 7
HARD_DIRECTION_COUNT = 7
HARD_QUOTIENT_SUM = 10
DISTANCES = tuple(range(1, 7))
DISTANCE_EDGE_COUNT = 13
OPPOSITE_PARALLEL_COUNT = 3
OPPOSITE_SCALED_MEAN = 20
OPPOSITE_TOTAL_W = -20
OPPOSITE_L1_BOUND = 56
OPPOSITE_CUT_UPPER_BOUND = -10
TRANSLATED_CUT_UPPER_BOUND = P * OPPOSITE_CUT_UPPER_BOUND
AGGREGATE_LOWER_BOUND = -13
AGGREGATE_UPPER_BOUND = 18
MOMENT_CANDIDATE_COUNT = 32_313
TRANSLATED_CUT_VECTOR_COUNT = 74


EXPECTED_NINE_CUT_VECTORS = (
    (2, 4, 6, 8, 10, 12),
    (8, 8, 6, 8, 8, 4),
    (8, 8, 8, 6, 4, 8),
    (12, 2, 10, 4, 8, 6),
    (4, 8, 10, 8, 6, 6),
    (8, 8, 4, 8, 6, 8),
    (6, 12, 8, 2, 4, 10),
    (10, 6, 4, 12, 2, 8),
    (8, 8, 6, 6, 10, 4),
)
EXPECTED_NINE_REPRESENTATIVES = (
    (0, 1, 2, 3, 4, 5, 6),
    (0, 1, 2, 4, 7, 8, 10),
    (0, 1, 2, 5, 6, 8, 10),
    (0, 1, 3, 5, 7, 9, 11),
    (0, 1, 2, 3, 4, 8, 9),
    (0, 1, 2, 4, 5, 7, 10),
    (0, 1, 2, 5, 6, 9, 10),
    (0, 1, 3, 5, 6, 8, 11),
    (0, 1, 2, 4, 7, 8, 11),
)
EXPECTED_GREEDY_ELIMINATIONS = (14_222, 9_967, 6_087, 1_395, 417, 168, 45, 8, 4)
EXPECTED_GREEDY_REMAINDERS = (18_091, 8_124, 2_037, 642, 225, 57, 12, 4, 0)


Aggregate = tuple[int, int, int, int, int, int]
CutVector = tuple[int, int, int, int, int, int]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _digest(rows: Iterable[tuple[int, ...]]) -> str:
    payload = ";".join(",".join(map(str, row)) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def p13_binary_radon_dependency() -> dict[str, object]:
    """Import Proposition 15.692 and specialize its dimensions to p=13."""
    radon = affine_binary_radon_isomorphism(P)
    hard_direction_count = HARD_DIRECTION_COUNT
    hard_profile_dimension = hard_direction_count * (P - 1)
    proved = bool(
        radon["proved"]
        and radon["inverse"] == "x = A^T r"
        and radon["source_dimension"] == P * P - 1
        and radon["target_dimension"] == P * P - 1
        and hard_profile_dimension == 84
    )
    _require(proved, "the Proposition 15.692 dependency changed")
    return {
        "dependency": "Proposition 15.692",
        "new_radon_proposition_asserted": False,
        "field": "F_2",
        "p": P,
        "radon_isomorphism": radon,
        "opposite_zero_profile_direction_count": HARD_DIRECTION_COUNT,
        "remaining_hard_even_profile_dimension": hard_profile_dimension,
        "boundary_reconstruction": (
            "1_D(x)=sum_(hard L) r_L(L(x)) because all seven opposite "
            "profiles vanish"
        ),
        "exact_hard_profile": "r_L=1+delta_(j_L)",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def generic_p13_hard_partition_split() -> dict[str, object]:
    """Enumerate the three excess-three hard quotient partitions."""
    partitions = tuple(
        row
        for row in combinations_with_replacement(range(1, 5), HARD_DIRECTION_COUNT)
        if sum(row) == HARD_QUOTIENT_SUM
    )
    expected = (
        (1, 1, 1, 1, 1, 1, 4),
        (1, 1, 1, 1, 1, 2, 3),
        (1, 1, 1, 1, 2, 2, 2),
    )
    star_power_sums = {
        degree: {
            center: sum(
                pow((center - label) % P, degree, P)
                for label in range(P)
                if label != center
            )
            % P
            for center in range(P)
        }
        for degree in (2, 4)
    }
    rows = []
    for partition in partitions:
        exact_count = partition.count(1)
        forced_zero_degrees = [
            degree for degree in (2, 4) if exact_count > degree
        ]
        rows.append(
            {
                "hard_quotient_partition": list(partition),
                "exact_hard_star_count": exact_count,
                "global_even_moments_forced_zero": forced_zero_degrees,
                "excluded_by_translation_average": exact_count >= 5,
            }
        )
    proved = bool(
        partitions == expected
        and all(
            value == 0
            for centers in star_power_sums.values()
            for value in centers.values()
        )
        and [row["exact_hard_star_count"] for row in rows] == [6, 5, 4]
        and [row["global_even_moments_forced_zero"] for row in rows]
        == [[2, 4], [2, 4], [2]]
    )
    _require(proved, "the generic p=13 hard partition split changed")
    return {
        "p": P,
        "hard_direction_count": HARD_DIRECTION_COUNT,
        "hard_quotient_constraints": "k_L>=1 and sum_L k_L=10",
        "hard_excess_units": HARD_QUOTIENT_SUM - HARD_DIRECTION_COUNT,
        "partitions": rows,
        "exact_star_power_sum_checks": star_power_sums,
        "homogeneous_binary_root_rule": (
            "more than d distinct projective roots force a degree-d form "
            "to vanish identically"
        ),
        "at_least_five_exact_stars_force_M2_and_M4_zero": True,
        "four_exact_stars_force_M2_but_not_M4_zero": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def opposite_cell_aggregate_reduction() -> dict[str, object]:
    """Derive the six cyclic aggregate variables and their exact bounds."""
    minimum_integer_entry = min(
        value for value in range(-100, 101) if 20 + 12 * value >= 0
    )
    aggregate_lower = DISTANCE_EDGE_COUNT * minimum_integer_entry
    aggregate_positive_mass_bound = (
        OPPOSITE_L1_BOUND + OPPOSITE_TOTAL_W
    ) // 2
    aggregate_upper = aggregate_positive_mass_bound
    proved = bool(
        minimum_integer_entry == -1
        and aggregate_lower == AGGREGATE_LOWER_BOUND
        and aggregate_positive_mass_bound == AGGREGATE_UPPER_BOUND
        and TRANSLATED_CUT_UPPER_BOUND == -130
    )
    _require(proved, "the opposite aggregate bounds changed")
    return {
        "p": P,
        "opposite_parallel_count_Q": OPPOSITE_PARALLEL_COUNT,
        "opposite_scaled_mean": OPPOSITE_SCALED_MEAN,
        "normalized_matrix": "W=epsilon_L*K^L",
        "coefficient_sum": OPPOSITE_TOTAL_W,
        "l1_bound": OPPOSITE_L1_BOUND,
        "row_sums_even": True,
        "B_formula": "B(X)=-5-cut_W(X)/2",
        "balanced_cut_upper_bound": OPPOSITE_CUT_UPPER_BOUND,
        "pair_in_conditional_mean": "E[B|i,j in X]=(20+12*w_ij)/44",
        "entry_lower_bound": minimum_integer_entry,
        "distance_classes": list(DISTANCES),
        "edges_per_distance_class": DISTANCE_EDGE_COUNT,
        "aggregate_definition": (
            "n_a=sum W_ij over the 13 unordered pairs with i-j=+-a"
        ),
        "aggregate_lower_bound": aggregate_lower,
        "aggregate_upper_bound": aggregate_upper,
        "aggregate_l1_bound": OPPOSITE_L1_BOUND,
        "aggregate_sum": OPPOSITE_TOTAL_W,
        "aggregate_moment_congruences": [
            "sum_(a=1)^6 a^2*n_a=0 mod 13",
            "sum_(a=1)^6 a^4*n_a=0 mod 13",
        ],
        "translated_cut_identity": (
            "sum_(t in F_13) cut_W(X+t)=sum_(a=1)^6 "
            "|X triangle (X+a)|*n_a"
        ),
        "translated_cut_upper_bound": TRANSLATED_CUT_UPPER_BOUND,
        "relaxation_status": (
            "necessary aggregate relaxation; feasibility would not imply a "
            "directional cell or common graph"
        ),
        "proved": proved,
    }


def translated_cut_vector(subset: Iterable[int]) -> CutVector:
    """Return ``(|X triangle (X+a)|)_(a=1..6)`` in F_13."""
    values = tuple(sorted(int(value) for value in subset))
    if len(values) != M or len(set(values)) != M or not all(
        0 <= value < P for value in values
    ):
        raise ValueError("need seven distinct elements of F_13")
    chosen = set(values)
    return tuple(
        sum(
            (value in chosen) != ((value + distance) % P in chosen)
            for value in range(P)
        )
        for distance in DISTANCES
    )  # type: ignore[return-value]


@lru_cache(maxsize=1)
def translated_cut_vector_catalog() -> dict[str, object]:
    """Generate all 74 translation-average coefficient vectors."""
    vector_to_representative: dict[CutVector, tuple[int, ...]] = {}
    for subset in combinations(range(P), M):
        vector = translated_cut_vector(subset)
        vector_to_representative.setdefault(vector, subset)
    vectors = tuple(sorted(vector_to_representative))
    proved = bool(
        len(vectors) == TRANSLATED_CUT_VECTOR_COUNT
        and all(
            value % 2 == 0 and 0 <= value <= 12
            for vector in vectors
            for value in vector
        )
        and all(sum(vector) == 42 for vector in vectors)
    )
    _require(proved, "the translated-cut vector catalog changed")
    return {
        "p": P,
        "slice_size": M,
        "seven_set_count": comb(P, M),
        "distinct_translated_cut_vector_count": len(vectors),
        "vectors": [list(vector) for vector in vectors],
        "representatives": {
            ",".join(map(str, vector)): list(vector_to_representative[vector])
            for vector in vectors
        },
        "catalog_sha256": _digest(vectors),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def _moment_candidates() -> tuple[Aggregate, ...]:
    """Enumerate aggregates by an exact modular row-reduction parameterization."""
    values = tuple(range(AGGREGATE_LOWER_BOUND, AGGREGATE_UPPER_BOUND + 1))
    values_by_residue = {
        residue: tuple(value for value in values if value % P == residue)
        for residue in range(P)
    }
    candidates: list[Aggregate] = []
    for n4, n5, n6 in product(values, repeat=3):
        # RREF over F_13 of the sum, second-moment, and fourth-moment rows.
        residues = (
            (9 - 10 * n4 - n5 - 10 * n6) % P,
            (12 - 6 * n4 - 3 * n5 - 2 * n6) % P,
            (11 - 11 * n4 - 10 * n5 - 2 * n6) % P,
        )
        for n1, n2, n3 in product(
            *(values_by_residue[residue] for residue in residues)
        ):
            row: Aggregate = (n1, n2, n3, n4, n5, n6)
            if sum(row) != OPPOSITE_TOTAL_W:
                continue
            if sum(abs(value) for value in row) > OPPOSITE_L1_BOUND:
                continue
            _require(
                all(
                    sum(
                        pow(distance, degree, P) * value
                        for distance, value in zip(DISTANCES, row)
                    )
                    % P
                    == 0
                    for degree in (2, 4)
                ),
                "the modular aggregate parameterization admitted a bad row",
            )
            candidates.append(row)
    result = tuple(candidates)
    _require(
        len(result) == MOMENT_CANDIDATE_COUNT,
        "the aggregate moment candidate count changed",
    )
    return result


def _independent_meet_in_middle_count() -> int:
    """Count the same pre-cut rows without using the RREF parameterization."""
    values = tuple(range(AGGREGATE_LOWER_BOUND, AGGREGATE_UPPER_BOUND + 1))
    right: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for row in product(values, repeat=3):
        key = (
            sum(row),
            sum(
                pow(distance, 2, P) * value
                for distance, value in zip((4, 5, 6), row)
            )
            % P,
            sum(
                pow(distance, 4, P) * value
                for distance, value in zip((4, 5, 6), row)
            )
            % P,
        )
        right[key].append(sum(abs(value) for value in row))
    for l1_values in right.values():
        l1_values.sort()

    count = 0
    for row in product(values, repeat=3):
        left_m2 = sum(
            pow(distance, 2, P) * value
            for distance, value in zip((1, 2, 3), row)
        )
        left_m4 = sum(
            pow(distance, 4, P) * value
            for distance, value in zip((1, 2, 3), row)
        )
        key = (
            OPPOSITE_TOTAL_W - sum(row),
            (-left_m2) % P,
            (-left_m4) % P,
        )
        remaining_l1 = OPPOSITE_L1_BOUND - sum(abs(value) for value in row)
        count += bisect_right(right.get(key, ()), remaining_l1)
    return count


@lru_cache(maxsize=1)
def translated_cut_nine_vector_certificate() -> dict[str, object]:
    """Greedily and deterministically extract nine vectors killing all rows."""
    catalog = translated_cut_vector_catalog()
    vectors = tuple(tuple(row) for row in catalog["vectors"])
    candidates = _moment_candidates()
    remaining = set(range(len(candidates)))
    selected: list[CutVector] = []
    elimination_counts: list[int] = []
    remainder_counts: list[int] = []
    available = list(vectors)
    while remaining:
        best_vector: CutVector | None = None
        best_killed: set[int] = set()
        for vector in available:
            killed = {
                index
                for index in remaining
                if sum(
                    coefficient * value
                    for coefficient, value in zip(vector, candidates[index])
                )
                > TRANSLATED_CUT_UPPER_BOUND
            }
            if len(killed) > len(best_killed):
                best_vector = vector  # sorted order makes ties deterministic
                best_killed = killed
        _require(best_vector is not None and best_killed, "greedy cover stalled")
        selected.append(best_vector)
        elimination_counts.append(len(best_killed))
        remaining.difference_update(best_killed)
        remainder_counts.append(len(remaining))
        available.remove(best_vector)

    representatives = tuple(
        next(
            subset
            for subset in combinations(range(P), M)
            if translated_cut_vector(subset) == vector
        )
        for vector in selected
    )
    independent_count = _independent_meet_in_middle_count()
    proved = bool(
        tuple(selected) == EXPECTED_NINE_CUT_VECTORS
        and representatives == EXPECTED_NINE_REPRESENTATIVES
        and tuple(elimination_counts) == EXPECTED_GREEDY_ELIMINATIONS
        and tuple(remainder_counts) == EXPECTED_GREEDY_REMAINDERS
        and independent_count == MOMENT_CANDIDATE_COUNT
        and not remaining
    )
    _require(proved, "the deterministic nine-vector certificate changed")
    return {
        "aggregate_bounds": [AGGREGATE_LOWER_BOUND, AGGREGATE_UPPER_BOUND],
        "aggregate_sum": OPPOSITE_TOTAL_W,
        "aggregate_l1_bound": OPPOSITE_L1_BOUND,
        "moment_degrees": [2, 4],
        "candidate_count_after_sum_l1_moments": len(candidates),
        "independent_meet_in_middle_candidate_count": independent_count,
        "candidate_catalog_sha256": _digest(candidates),
        "deterministic_selection_rule": (
            "from lexicographically sorted 74-vector catalog, repeatedly "
            "choose the first vector eliminating the most remaining rows"
        ),
        "selected_vector_count": len(selected),
        "selected_vectors": [list(vector) for vector in selected],
        "representative_seven_sets": [list(row) for row in representatives],
        "eliminated_at_each_step": elimination_counts,
        "remaining_after_each_step": remainder_counts,
        "remaining_after_nine_vectors": len(remaining),
        "selected_vectors_sha256": _digest(selected),
        "full_74_vector_catalog_sha256": catalog["catalog_sha256"],
        "pure_integer_enumeration_infeasible": not remaining,
        "proved": proved,
    }


def _model_textproto_sha256(model: cp_model.CpModel) -> str:
    return hashlib.sha256(str(model.Proto()).encode("utf-8")).hexdigest()


@lru_cache(maxsize=1)
def independent_six_variable_cpsat_check() -> dict[str, object]:
    """Independently certify the nine-vector aggregate model infeasible."""
    certificate = translated_cut_nine_vector_certificate()
    selected = tuple(tuple(row) for row in certificate["selected_vectors"])
    model = cp_model.CpModel()
    aggregates = [
        model.NewIntVar(
            AGGREGATE_LOWER_BOUND,
            AGGREGATE_UPPER_BOUND,
            f"n_{distance}",
        )
        for distance in DISTANCES
    ]
    magnitudes = [
        model.NewIntVar(0, max(-AGGREGATE_LOWER_BOUND, AGGREGATE_UPPER_BOUND), f"abs_n_{distance}")
        for distance in DISTANCES
    ]
    for aggregate, magnitude in zip(aggregates, magnitudes):
        model.AddAbsEquality(magnitude, aggregate)
    model.Add(sum(aggregates) == OPPOSITE_TOTAL_W)
    model.Add(sum(magnitudes) <= OPPOSITE_L1_BOUND)
    for degree in (2, 4):
        quotient = model.NewIntVar(-100, 100, f"moment_{degree}_quotient")
        model.Add(
            sum(
                pow(distance, degree, P) * aggregate
                for distance, aggregate in zip(DISTANCES, aggregates)
            )
            == P * quotient
        )
    for index, vector in enumerate(selected):
        model.Add(
            sum(
                coefficient * aggregate
                for coefficient, aggregate in zip(vector, aggregates)
            )
            <= TRANSLATED_CUT_UPPER_BOUND
        ).WithName(f"translated_cut_{index}")

    validation = model.Validate()
    _require(not validation, f"invalid six-variable model: {validation}")
    model_hash = _model_textproto_sha256(model)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.cp_model_presolve = True
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    _require(infeasible, f"the six-variable model is {status_name}")
    return {
        "independent_check": "six-variable exact CP-SAT",
        "solver": "OR-Tools CP-SAT",
        "solver_version": ORTOOLS_VERSION,
        "integer_variable_count": len(model.Proto().variables),
        "model_constraint_count": len(model.Proto().constraints),
        "model_validation": validation,
        "model_textproto_sha256": model_hash,
        "search_workers": 1,
        "random_seed": 0,
        "cp_model_presolve": True,
        "solver_status": status_name,
        "solver_wall_time_seconds": solver.WallTime(),
        "branches": solver.NumBranches(),
        "conflicts": solver.NumConflicts(),
        "infeasible": infeasible,
        "proved": infeasible,
    }


@lru_cache(maxsize=1)
def proposition_15740() -> dict[str, object]:
    """Package the proved p=13 generic partition split."""
    radon = p13_binary_radon_dependency()
    partitions = generic_p13_hard_partition_split()
    reduction = opposite_cell_aggregate_reduction()
    vectors = translated_cut_vector_catalog()
    certificate = translated_cut_nine_vector_certificate()
    independent = independent_six_variable_cpsat_check()
    proved = bool(
        radon["proved"]
        and partitions["proved"]
        and reduction["proved"]
        and vectors["proved"]
        and certificate["proved"]
        and independent["proved"]
    )
    _require(proved, "Proposition 15.740 package failed")
    return {
        "prop": "15.740",
        "title": "translation-average split of the generic p=13 fourth shell",
        "result_status": "proved branch split with exhaustive finite certificate",
        "p": P,
        "layer_index_t": 3,
        "original_k": 4 * P + 6,
        "radon_dependency": radon,
        "hard_partition_split": partitions,
        "opposite_aggregate_reduction": reduction,
        "translated_cut_vector_catalog": vectors,
        "nine_vector_certificate": certificate,
        "independent_solver_check": independent,
        "excluded_hard_quotient_partitions": [
            [1, 1, 1, 1, 1, 1, 4],
            [1, 1, 1, 1, 1, 2, 3],
        ],
        "remaining_hard_quotient_partitions": [[1, 1, 1, 1, 2, 2, 2]],
        "p13_generic_partitions_with_at_least_five_exact_stars_excluded": True,
        "p13_generic_four_exact_partition_closed": False,
        "p13_generic_t3_branch_closed": False,
        "p13_k_eq_58_closed": False,
        "generic_p_ge_17_t3_branch_closed": False,
        "k_eq_4p_plus_6_shell_closed": False,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_exact_gate": (
            "exclude the common 59-edge completion of four exact P=5 stars, "
            "three elevated P=6 hard cells, and seven Q=3,b=0 opposite cells"
        ),
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15740.json"
    target.write_text(json.dumps(proposition_15740(), indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    theorem = proposition_15740()
    target = write_evidence()
    print(
        "Prop. 15.740: generic p=13 five/six-exact partitions excluded; "
        "only 1^4 2^3 remains"
    )
    print(
        f"aggregate candidates={theorem['nine_vector_certificate']['candidate_count_after_sum_l1_moments']}; "
        f"after nine cuts={theorem['nine_vector_certificate']['remaining_after_nine_vectors']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
