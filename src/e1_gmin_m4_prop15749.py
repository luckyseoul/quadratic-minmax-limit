#!/usr/bin/env python3
r"""Prop. 15.749 -- translated-cut moments close ``p=13,t=4,u=4``.

Proposition 15.748 leaves, in the ``P=5`` branch, two exact literal
opposite directions and five excess-one ``Q=4`` directions.  It also gives
336 exact moment-level candidates for each choice of the hard sign.

For a ``Q=4`` direction let ``q=(q_1,...,q_6)`` be its signed cyclic
distance row.  The common signed total ``hT=9`` and nonnegativity of the
cell give

    sum q_a = -13,
    sum |q_a| <= 57,
    c dot q <= -78

for every one of the 74 already-certified translated-cut vectors ``c``.
Exact dual combinations of those cuts imply ``-52/9 <= q_a <= 26/15``;
therefore every integral coordinate lies in ``[-5,1]``.  Exhaustive exact
list recovery in that six-variable box gives 522 rows and 492 moment
triples ``(N2,N4,N6)``.

Reconstructing the full binary forms from each Proposition 15.748 record
shows that their opposite-direction evaluation alphabet has 48 triples.
Its intersection with the 492 admissible ``Q=4`` triples consists of twelve
triples, all with fourth coordinate zero.  Thus every one of the five
``Q=4`` directions would be a root of the common quartic.  Together with
the two literal roots this gives seven projective roots, while every hard
fourth moment is nonzero.  The quartic cannot exist, closing ``P=5`` and,
with Proposition 15.747, all of ``u=4``.

This is an exact aggregate certificate, not a close of ``u=6`` or global
residual (ii).
"""
from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.p13_p5_literal_interpolation import (  # noqa: E402
    P,
    POINTS,
    audit,
    form_value,
    quartic_code_membership_setup,
    root_product_values,
)

from e1_gmin_m4_prop15740 import translated_cut_vector_catalog  # noqa: E402
from e1_gmin_m4_prop15748 import proposition_15748  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402


DISTANCES = tuple(range(1, 7))
Q4_SUM = -13
Q4_L1_BOUND = 57
Q4_TRANSLATED_CUT_UPPER = -78
EXPECTED_Q4_ROW_COUNT = 522
EXPECTED_Q4_MOMENT_COUNT = 492
EXPECTED_EVALUATION_ALPHABET_SIZE = 48
EXPECTED_Q4_ROW_SHA256 = (
    "2e4cf7f733ffd6d85a68a6b37ebd380d93d962e6119e9306673ddd3a1df8cb35"
)
EXPECTED_Q4_MOMENT_SHA256 = (
    "16e529a3ea7263b66f7af61cf6eaa59441747622c62d0dbc9a8267837c76f378"
)
EXPECTED_EVALUATION_ALPHABET_SHA256 = (
    "580f8080d29895e7c320bde39a71b5c804c110fe86657460ccc347294974f561"
)
EXPECTED_INTERSECTION_SHA256 = (
    "5ea8fc294df3cc5ff587b153a50b945834836ba69f037045fbdfc3766ba5fe94"
)
EXPECTED_INTERSECTION = (
    (1, 0, 3),
    (2, 0, 1),
    (3, 0, 3),
    (4, 0, 10),
    (5, 0, 1),
    (6, 0, 1),
    (7, 0, 12),
    (8, 0, 12),
    (9, 0, 3),
    (10, 0, 10),
    (11, 0, 12),
    (12, 0, 10),
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _digest(rows: object) -> str:
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _dot(left: Iterable[int | Fraction], right: Iterable[int | Fraction]) -> Fraction:
    return sum(
        (Fraction(a) * Fraction(b) for a, b in zip(left, right)),
        Fraction(),
    )


def _distance_permutation(multiplier: int) -> tuple[int, ...]:
    """Permutation of ``F_13^*/{+-1}`` induced by multiplication."""
    return tuple(min(multiplier * value % P, -multiplier * value % P) - 1 for value in DISTANCES)


@lru_cache(maxsize=1)
def translated_cut_coordinate_bound_certificate() -> dict[str, object]:
    """Validate exact dual cut certificates giving ``q_a in [-5,1]``."""
    catalog = translated_cut_vector_catalog()
    vectors = tuple(tuple(int(value) for value in row) for row in catalog["vectors"])
    vector_set = set(vectors)
    ones = (1,) * 6

    # e_1 = 19/9*1 - c_0/18 - c_6/6 - c_34/18.
    lower_terms = ((0, Fraction(-1, 18)), (6, Fraction(-1, 6)), (34, Fraction(-1, 18)))
    lower_equality = tuple(
        Fraction(19, 9) + sum(weight * vectors[index][coordinate] for index, weight in lower_terms)
        for coordinate in range(6)
    )
    lower_rhs = Fraction(19, 9) * Q4_SUM + sum(
        weight * Q4_TRANSLATED_CUT_UPPER for _index, weight in lower_terms
    )

    # -e_1 = 29/15*1 - c_63/15 - c_69/30 - c_71/6 - c_73/30.
    upper_terms = (
        (63, Fraction(-1, 15)),
        (69, Fraction(-1, 30)),
        (71, Fraction(-1, 6)),
        (73, Fraction(-1, 30)),
    )
    upper_equality = tuple(
        Fraction(29, 15) + sum(weight * vectors[index][coordinate] for index, weight in upper_terms)
        for coordinate in range(6)
    )
    upper_rhs_for_minus_q = Fraction(29, 15) * Q4_SUM + sum(
        weight * Q4_TRANSLATED_CUT_UPPER for _index, weight in upper_terms
    )

    permutations = tuple(_distance_permutation(multiplier) for multiplier in DISTANCES)
    invariant = all(
        tuple(vector[index] for index in permutation) in vector_set
        for permutation in permutations
        for vector in vectors
    )
    transitive = {permutation[0] for permutation in permutations} == set(range(6))

    proved = bool(
        catalog["proved"]
        and len(vectors) == 74
        and lower_equality == (1, 0, 0, 0, 0, 0)
        and lower_rhs == Fraction(-52, 9)
        and upper_equality == (-1, 0, 0, 0, 0, 0)
        and upper_rhs_for_minus_q == Fraction(-26, 15)
        and invariant
        and transitive
        and -5 >= Fraction(-52, 9)
        and 1 <= Fraction(26, 15)
    )
    _require(proved, "the translated-cut coordinate bound changed")
    return {
        "translated_cut_vector_count": len(vectors),
        "sum_constraint": Q4_SUM,
        "cut_upper": Q4_TRANSLATED_CUT_UPPER,
        "lower_dual_terms": [[index, str(weight)] for index, weight in lower_terms],
        "lower_vector_identity": [str(value) for value in lower_equality],
        "rational_lower_bound": str(lower_rhs),
        "upper_dual_terms": [[index, str(weight)] for index, weight in upper_terms],
        "upper_vector_identity": [str(value) for value in upper_equality],
        "rational_upper_bound": str(-upper_rhs_for_minus_q),
        "multiplicative_distance_permutations": [list(row) for row in permutations],
        "cut_catalog_invariant_under_all_six_permutations": invariant,
        "distance_action_transitive": transitive,
        "integral_coordinate_bounds": [-5, 1],
        "proved": proved,
    }


def _moments(row: tuple[int, ...]) -> tuple[int, int, int]:
    return tuple(
        sum(pow(distance, degree, P) * row[distance - 1] for distance in DISTANCES) % P
        for degree in (2, 4, 6)
    )


@lru_cache(maxsize=1)
def q4_translated_cut_moment_certificate() -> dict[str, object]:
    """Recover every admissible integral six-bin ``Q=4`` row exactly."""
    bounds = translated_cut_coordinate_bound_certificate()
    vectors = tuple(
        tuple(int(value) for value in row)
        for row in translated_cut_vector_catalog()["vectors"]
    )
    rows = tuple(
        row
        for row in product(range(-5, 2), repeat=6)
        if sum(row) == Q4_SUM
        and sum(abs(value) for value in row) <= Q4_L1_BOUND
        and all(_dot(cut, row) <= Q4_TRANSLATED_CUT_UPPER for cut in vectors)
    )
    moments = tuple(sorted({_moments(row) for row in rows}))
    row_digest = _digest(rows)
    moment_digest = _digest(moments)
    proved = bool(
        bounds["proved"]
        and len(rows) == EXPECTED_Q4_ROW_COUNT
        and len(moments) == EXPECTED_Q4_MOMENT_COUNT
        and row_digest == EXPECTED_Q4_ROW_SHA256
        and moment_digest == EXPECTED_Q4_MOMENT_SHA256
        and all(all(-5 <= value <= 1 for value in row) for row in rows)
    )
    _require(proved, "the Q=4 translated-cut moment list changed")
    return {
        "row_definition": {
            "sum": Q4_SUM,
            "l1_upper": Q4_L1_BOUND,
            "translated_cut_upper": Q4_TRANSLATED_CUT_UPPER,
            "coordinate_bounds": [-5, 1],
        },
        "candidate_box_size": 7**6,
        "admissible_row_count": len(rows),
        "admissible_row_sha256": row_digest,
        "moment_degrees": [2, 4, 6],
        "admissible_moment_count": len(moments),
        "admissible_moment_sha256": moment_digest,
        "moments": [list(row) for row in moments],
        "proved": proved,
    }


def _polynomial_product(left: Iterable[int], right: Iterable[int]) -> tuple[int, ...]:
    left = tuple(left)
    right = tuple(right)
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % P
    return tuple(result)


def _root_polynomial(roots: tuple[int, int]) -> tuple[int, ...]:
    result = (1,)
    for index in roots:
        x, y = POINTS[index]
        result = _polynomial_product(result, (y % P, (-x) % P))
    return result


def _m6_quotient_coefficients(
    hard: tuple[int, ...],
    roots: tuple[int, int],
    hard_n6: tuple[int, ...],
) -> tuple[int, ...]:
    root_values = root_product_values(roots)
    quotient_values = tuple(
        value * pow(root_values[index], -1, P) % P
        for value, index in zip(hard_n6, hard)
    )
    evaluation, inverse = quartic_code_membership_setup(hard)
    coefficients = tuple(
        sum(inverse[row][column] * quotient_values[column] for column in range(5)) % P
        for row in range(5)
    )
    _require(
        all(
            sum(coefficient * basis for coefficient, basis in zip(coefficients, row)) % P
            == value
            for row, value in zip(evaluation, quotient_values)
        ),
        "a stored M6 quotient left the quartic evaluation code",
    )
    return coefficients


@lru_cache(maxsize=1)
def survivor_q4_moment_intersection_certificate() -> dict[str, object]:
    """Intersect all stored ``z=2`` evaluations with admissible ``Q=4`` moments."""
    prior = proposition_15748()
    q4 = q4_translated_cut_moment_certificate()
    admissible = {tuple(int(value) for value in row) for row in q4["moments"]}
    raw = audit()
    sign_rows: dict[str, dict[str, object]] = {}
    all_intersections = []
    for hard_sign in (-1, 1):
        source = raw["sign_rows"][str(hard_sign)]
        hard = tuple(int(value) for value in source["hard_direction_indices"])
        opposite = tuple(int(value) for value in source["opposite_direction_indices"])
        evaluation_alphabet: set[tuple[int, int, int]] = set()
        compatible_count_histogram: dict[int, int] = {}
        every_hard_replay = True
        maximum_compatible = 0
        for survivor in source["z2_survivors"]:
            roots = tuple(int(value) for value in survivor["roots"])
            root_polynomial = _root_polynomial(roots)
            m2 = tuple(
                int(survivor["M2_scalar"]) * value % P
                for value in root_polynomial
            )
            m4 = _polynomial_product(
                root_polynomial,
                tuple(int(value) for value in survivor["M4_quotient_coefficients"]),
            )
            m6 = _polynomial_product(
                root_polynomial,
                _m6_quotient_coefficients(
                    hard,
                    roots,
                    tuple(int(value) for value in survivor["hard_N6"]),
                ),
            )
            forms = (m2, m4, m6)
            replayed_hard = tuple(
                tuple(form_value(form, POINTS[index]) for index in hard)
                for form in forms
            )
            every_hard_replay &= replayed_hard == (
                tuple(int(value) for value in survivor["hard_N2"]),
                tuple(int(value) for value in survivor["hard_N4"]),
                tuple(int(value) for value in survivor["hard_N6"]),
            )
            nonroots = tuple(index for index in opposite if index not in roots)
            _require(len(nonroots) == 5, "a z=2 survivor lost a Q=4 direction")
            evaluations = tuple(
                tuple((-form_value(form, POINTS[index])) % P for form in forms)
                for index in nonroots
            )
            evaluation_alphabet.update(evaluations)
            compatible = sum(value in admissible for value in evaluations)
            maximum_compatible = max(maximum_compatible, compatible)
            compatible_count_histogram[compatible] = (
                compatible_count_histogram.get(compatible, 0) + 1
            )

        intersection = tuple(sorted(evaluation_alphabet & admissible))
        evaluation_digest = _digest(tuple(sorted(evaluation_alphabet)))
        intersection_digest = _digest(intersection)
        row_proved = bool(
            source["z2_survivor_count"] == 336
            and every_hard_replay
            and len(evaluation_alphabet) == EXPECTED_EVALUATION_ALPHABET_SIZE
            and evaluation_digest == EXPECTED_EVALUATION_ALPHABET_SHA256
            and intersection == EXPECTED_INTERSECTION
            and intersection_digest == EXPECTED_INTERSECTION_SHA256
            and all(n4 == 0 for _n2, n4, _n6 in intersection)
            and compatible_count_histogram == {0: 252, 1: 42, 2: 42}
            and maximum_compatible == 2 < 5
        )
        _require(row_proved, f"the hard-sign {hard_sign} Q=4 intersection changed")
        sign_rows[str(hard_sign)] = {
            "hard_sign": hard_sign,
            "z2_survivor_count": int(source["z2_survivor_count"]),
            "opposite_nonroot_directions_per_survivor": 5,
            "evaluation_alphabet_size": len(evaluation_alphabet),
            "evaluation_alphabet_sha256": evaluation_digest,
            "admissible_intersection_size": len(intersection),
            "admissible_intersection": [list(row) for row in intersection],
            "admissible_intersection_sha256": intersection_digest,
            "every_intersection_fourth_moment_zero": True,
            "compatible_direction_count_histogram": {
                str(key): value for key, value in sorted(compatible_count_histogram.items())
            },
            "maximum_compatible_Q4_directions_in_one_survivor": maximum_compatible,
            "every_hard_value_vector_replayed": every_hard_replay,
            "proved": row_proved,
        }
        all_intersections.append(intersection)

    proved = bool(
        prior["proved"]
        and prior["P5_excess_partition_reduction"]["forced_z"] == 2
        and prior["P5_excess_partition_reduction"]["positive_opposite_excess_count"] == 5
        and raw["payload_sha256"]
        == prior["exact_literal_interpolation"]["raw_payload_sha256"]
        and all(row["proved"] for row in sign_rows.values())
        and all_intersections[0] == all_intersections[1] == EXPECTED_INTERSECTION
    )
    _require(proved, "the survivor/Q=4 moment intersection failed")
    return {
        "field": "F_13",
        "normalization": (
            "hard records store h*M_(2r); an opposite row stores "
            "q-moment=(-h)*M_(2r), hence the minus sign"
        ),
        "prior_raw_payload_sha256": raw["payload_sha256"],
        "sign_rows": sign_rows,
        "common_intersection": [list(row) for row in EXPECTED_INTERSECTION],
        "common_intersection_has_N4_zero": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p13_t4_u4_close() -> dict[str, object]:
    """Close the last ``P=5`` branch and hence all of ``u=4``."""
    prior = proposition_15748()
    intersection = survivor_q4_moment_intersection_certificate()
    literal_roots = int(prior["P5_excess_partition_reduction"]["forced_z"])
    q4_directions = int(
        prior["P5_excess_partition_reduction"]["positive_opposite_excess_count"]
    )
    total_roots = literal_roots + q4_directions
    quartic_degree = 4
    hard_fourth_never_zero = bool(
        prior["literal_root_and_hard_alphabet"]["hard_fourth_moment_never_zero"]
    )
    proved = bool(
        prior["proved"]
        and intersection["proved"]
        and literal_roots == 2
        and q4_directions == 5
        and total_roots == 7 > quartic_degree
        and intersection["common_intersection_has_N4_zero"]
        and hard_fourth_never_zero
    )
    _require(proved, "the p13 t4 u4 quartic contradiction failed")
    return {
        "p": 13,
        "t": 4,
        "k": 60,
        "u": 4,
        "P3_branch_closed_by_prop_15747": True,
        "P5_literal_root_count": literal_roots,
        "P5_Q4_direction_count": q4_directions,
        "every_admissible_Q4_direction_forces_M4_zero": True,
        "total_forced_projective_M4_roots": total_roots,
        "M4_homogeneous_degree": quartic_degree,
        "root_count_exceeds_degree": total_roots > quartic_degree,
        "M4_forced_identically_zero": True,
        "hard_fourth_moment_never_zero": hard_fourth_never_zero,
        "contradiction": "seven projective roots force M4=0, but every hard N4 is nonzero",
        "P5_branch_closed": True,
        "p13_t4_u4_closed": True,
        "remaining_p13_t4_residues": [6],
        "result_status": "proved branch theorem with exact aggregate certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15749() -> dict[str, object]:
    """Package Proposition 15.749 without changing the global predicate."""
    bounds = translated_cut_coordinate_bound_certificate()
    q4 = q4_translated_cut_moment_certificate()
    intersection = survivor_q4_moment_intersection_certificate()
    close = p13_t4_u4_close()
    proved = bool(
        bounds["proved"]
        and q4["proved"]
        and intersection["proved"]
        and close["proved"]
        and close["remaining_p13_t4_residues"] == [6]
    )
    _require(proved, "Proposition 15.749 failed")
    return {
        "prop": "15.749",
        "title": "Translated-cut moments close p13 t4 u4",
        "result_status": "proved branch theorem with exact aggregate certificate",
        "statement": (
            "the only Q4 moment triples compatible with a z=2 survivor have "
            "N4=0, so two literal plus five Q4 roots contradict nonzero M4"
        ),
        "translated_cut_coordinate_bounds": bounds,
        "Q4_translated_cut_moments": q4,
        "survivor_Q4_moment_intersection": intersection,
        "p13_t4_u4_close": close,
        "p13_t4_u4_closed": True,
        "p13_k_eq_60_closed": False,
        "remaining_p13_t4_residues": [6],
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    """Write the deterministic Proposition 15.749 evidence atomically."""
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15749.json"
    write_json_atomic(target, proposition_15749())
    return target


def main() -> None:
    result = proposition_15749()
    path = write_evidence()
    print(
        json.dumps(
            {
                "prop": result["prop"],
                "result_status": result["result_status"],
                "p13_t4_u4_closed": result["p13_t4_u4_closed"],
                "remaining_p13_t4_residues": result["remaining_p13_t4_residues"],
                "output": str(path),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
