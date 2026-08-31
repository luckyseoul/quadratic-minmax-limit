#!/usr/bin/env python3
"""Exact moment interpolation for the p=13,t=4,u=4 P=5 branch.

Once every minimum Q=3 opposite cell is the b=12 literal, it is a common
projective root of M2, M4, and M6.  This script checks the remaining root
patterns against the complete 69-element local moment alphabet of a hard
baseline pair plus an all-equal triple.
"""
from __future__ import annotations

import hashlib
import json
import sys
from itertools import combinations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402


P = 13
POINTS = tuple((1, slope) for slope in range(P)) + ((0, 1),)
FOURTH_ALPHABET = frozenset({7, 8, 11})


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def projective_direction(point: tuple[int, int]) -> tuple[int, int]:
    x, y = point
    if x:
        return 1, y * pow(x, -1, P) % P
    if y:
        return 0, 1
    raise ValueError("zero has no projective direction")


def direction_signs() -> tuple[int, ...]:
    _q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(P)

    def field_element(point: tuple[int, int]) -> int:
        return point[0] + P * point[1]

    signs = []
    for functional in POINTS:
        x, y = functional
        kernel = projective_direction((y % P, (-x) % P))
        signs.append(int(chi(field_element(kernel))))
    _require(signs.count(-1) == signs.count(1) == 7, "sign split changed")
    return tuple(signs)


def hard_moment_alphabet() -> frozenset[tuple[int, int, int]]:
    rows = set()
    for baseline in combinations(range(P), 2):
        for triple in combinations(range(P), 3):
            edges = (baseline,) + tuple(combinations(triple, 2))
            rows.add(
                tuple(
                    sum((left - right) ** degree for left, right in edges) % P
                    for degree in (2, 4, 6)
                )
            )
    _require(len(rows) == 69, "hard moment alphabet changed")
    return frozenset(rows)


def root_product_values(roots: tuple[int, ...]) -> tuple[int, ...]:
    values = []
    for r, s in POINTS:
        value = 1
        for root in roots:
            x, y = POINTS[root]
            value = value * (y * r - x * s) % P
        values.append(value)
    return tuple(values)


def form_value(coefficients: tuple[int, ...], point: tuple[int, int]) -> int:
    degree = len(coefficients) - 1
    r, s = point
    return sum(
        coefficient * pow(r, degree - index, P) * pow(s, index, P)
        for index, coefficient in enumerate(coefficients)
    ) % P


def invert_matrix(matrix: list[list[int]]) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    work = [
        [value % P for value in row]
        + [int(row_index == column) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if work[row][column] % P
        )
        work[column], work[pivot] = work[pivot], work[column]
        inverse = pow(work[column][column], -1, P)
        work[column] = [value * inverse % P for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                (value - factor * pivot_value) % P
                for value, pivot_value in zip(work[row], work[column])
            ]
    inverse = tuple(tuple(row[size:]) for row in work)
    _require(
        all(
            sum(matrix[i][k] * inverse[k][j] for k in range(size)) % P
            == int(i == j)
            for i in range(size)
            for j in range(size)
        ),
        "matrix inverse failed",
    )
    return inverse


def quartic_code_membership_setup(
    hard: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    evaluation = tuple(
        tuple(
            pow(POINTS[index][0], 4 - monomial, P)
            * pow(POINTS[index][1], monomial, P)
            % P
            for monomial in range(5)
        )
        for index in hard
    )
    inverse = invert_matrix([list(row) for row in evaluation[:5]])
    return evaluation, inverse


def quartic_code_contains(
    values: tuple[int, ...],
    evaluation: tuple[tuple[int, ...], ...],
    inverse: tuple[tuple[int, ...], ...],
) -> bool:
    coefficients = tuple(
        sum(inverse[row][column] * values[column] for column in range(5)) % P
        for row in range(5)
    )
    return all(
        sum(coefficient * basis for coefficient, basis in zip(coefficients, row))
        % P
        == value
        for row, value in zip(evaluation, values)
    )


def z4_and_z3_survivors(
    hard: tuple[int, ...], opposite: tuple[int, ...]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    z4 = []
    for roots in combinations(opposite, 4):
        root_values = root_product_values(roots)
        for scalar in range(1, P):
            values = tuple(scalar * root_values[index] % P for index in hard)
            if set(values) <= FOURTH_ALPHABET:
                z4.append({"roots": roots, "scalar": scalar, "values": values})

    z3 = []
    for roots in combinations(opposite, 3):
        root_values = root_product_values(roots)
        for coefficients in product(range(P), repeat=2):
            if coefficients == (0, 0):
                continue
            values = tuple(
                root_values[index] * form_value(coefficients, POINTS[index]) % P
                for index in hard
            )
            if set(values) <= FOURTH_ALPHABET:
                z3.append(
                    {
                        "roots": roots,
                        "linear_coefficients": coefficients,
                        "values": values,
                    }
                )
    return z4, z3


def z2_survivors(
    hard: tuple[int, ...], opposite: tuple[int, ...]
) -> tuple[list[dict[str, object]], int, int]:
    alphabet = hard_moment_alphabet()
    allowed_24 = {(n2, n4) for n2, n4, _n6 in alphabet}
    n6_by_24 = {
        pair: tuple(sorted(n6 for n2, n4, n6 in alphabet if (n2, n4) == pair))
        for pair in allowed_24
    }
    evaluation, inverse = quartic_code_membership_setup(hard)
    survivors = []
    q2_candidates = 0
    n6_vectors_checked = 0
    for roots in combinations(opposite, 2):
        root_values = root_product_values(roots)
        hard_root_values = tuple(root_values[index] for index in hard)
        _require(all(hard_root_values), "hard and opposite directions overlap")
        for scalar in range(P):
            n2_values = tuple(scalar * value % P for value in hard_root_values)
            for quadratic in product(range(P), repeat=3):
                n4_values = tuple(
                    root_values[index]
                    * form_value(quadratic, POINTS[index])
                    % P
                    for index in hard
                )
                pairs_24 = tuple(zip(n2_values, n4_values))
                if not all(pair in allowed_24 for pair in pairs_24):
                    continue
                q2_candidates += 1
                options = tuple(n6_by_24[pair] for pair in pairs_24)
                for n6_values in product(*options):
                    n6_vectors_checked += 1
                    quotient_values = tuple(
                        n6 * pow(root_value, -1, P) % P
                        for n6, root_value in zip(n6_values, hard_root_values)
                    )
                    if quartic_code_contains(quotient_values, evaluation, inverse):
                        survivors.append(
                            {
                                "roots": roots,
                                "M2_scalar": scalar,
                                "M4_quotient_coefficients": quadratic,
                                "hard_N2": n2_values,
                                "hard_N4": n4_values,
                                "hard_N6": n6_values,
                            }
                        )
    return survivors, q2_candidates, n6_vectors_checked


def audit() -> dict[str, object]:
    signs = direction_signs()
    alphabet = hard_moment_alphabet()
    sign_rows = {}
    payload = []
    for hard_sign in (-1, 1):
        hard = tuple(index for index, sign in enumerate(signs) if sign == hard_sign)
        opposite = tuple(index for index, sign in enumerate(signs) if sign == -hard_sign)
        z4, z3 = z4_and_z3_survivors(hard, opposite)
        z2, q2_candidates, n6_vectors_checked = z2_survivors(hard, opposite)
        row = {
            "hard_sign": hard_sign,
            "hard_direction_indices": hard,
            "opposite_direction_indices": opposite,
            "z4_survivor_count": len(z4),
            "z3_survivor_count": len(z3),
            "z2_M2_M4_candidate_count": q2_candidates,
            "z2_N6_vectors_checked": n6_vectors_checked,
            "z2_survivor_count": len(z2),
            "z4_survivors": z4,
            "z3_survivors": z3,
            "z2_survivors": z2,
        }
        sign_rows[str(hard_sign)] = row
        payload.append(json.dumps(row, sort_keys=True))
    all_empty = all(
        not row["z4_survivors"]
        and not row["z3_survivors"]
        and not row["z2_survivors"]
        for row in sign_rows.values()
    )
    return {
        "field": "F_13",
        "projective_points": [list(point) for point in POINTS],
        "direction_signs": signs,
        "hard_moment_alphabet_size": len(alphabet),
        "hard_moment_alphabet": [list(row) for row in sorted(alphabet)],
        "fourth_moment_alphabet_under_M2_zero": sorted(FOURTH_ALPHABET),
        "sign_rows": sign_rows,
        "all_z2_z3_z4_interpolation_cases_empty": all_empty,
        "payload_sha256": hashlib.sha256("\n".join(payload).encode()).hexdigest(),
        "proved": all_empty,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
