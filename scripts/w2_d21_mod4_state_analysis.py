#!/usr/bin/env python3
"""Analyze local mod-four states after the first d=21 boundary gate stalls."""
from __future__ import annotations

import argparse
import collections
import itertools
import json
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_boundary_gpu import FACTORS, remainder_bits  # noqa: E402
from w2_d21_projective_invariants import (  # noqa: E402
    multiply,
    power,
    reciprocal_to_common,
)


def polynomial_bits(coefficients: list[int]) -> int:
    return sum((value & 1) << index for index, value in enumerate(coefficients))


def common_pair(coefficients: list[int]) -> tuple[int, int]:
    value = polynomial_bits(coefficients)
    return (
        remainder_bits(value, FACTORS[0]),
        reciprocal_to_common(remainder_bits(value, FACTORS[1])),
    )


def quotient_parity(item: dict) -> tuple[int, int]:
    return common_pair([value & 1 for value in item["quotient_mod4"]])


def half_pair(item: dict) -> tuple[int, int] | None:
    value = item["half_lift_mod2"]
    return None if value is None else common_pair(value)


def binary_rank(values: list[int]) -> int:
    basis: list[int] = []
    for value in values:
        reduced = value
        for pivot in basis:
            reduced = min(reduced, reduced ^ pivot)
        if reduced:
            basis.append(reduced)
            basis.sort(reverse=True)
    return len(basis)


def vector_bits(values: tuple[int, ...]) -> int:
    result = 0
    offset = 0
    for value in values:
        result |= value << offset
        offset += 3
    return result


def monomials(max_degree: int) -> list[tuple[int, ...]]:
    return [()] + [
        indices
        for degree in range(1, max_degree + 1)
        for indices in itertools.combinations(range(6), degree)
    ]


def evaluate_boolean_monomials(
    values: tuple[int, ...], terms: list[tuple[int, ...]]
) -> int:
    return sum(
        1 << index
        for index, term in enumerate(terms)
        if all(values[variable] for variable in term)
    )


def support_polynomial_row(
    patterns: set[tuple[int, ...]], degree: int
) -> dict:
    terms = monomials(degree)
    evaluations = [evaluate_boolean_monomials(row, terms) for row in patterns]
    rank = binary_rank(evaluations)
    bad_separable = []
    for left, right in itertools.product(range(2), repeat=2):
        bad = evaluate_boolean_monomials((left, right, 0, 0, 0, 0), terms)
        bad_separable.append(binary_rank(evaluations + [bad]) > rank)
    return {
        "degree": degree,
        "features": len(terms),
        "rank": rank,
        "relation_dimension": len(terms) - rank,
        "all_counterexample_supports_separable": all(bad_separable),
        "counterexample_supports_separable": bad_separable,
    }


def solve_f8(matrix: list[list[int]], rhs: list[int]) -> tuple[list[int] | None, int]:
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    pivot_columns = []
    pivot_row = 0
    for column in range(len(matrix[0])):
        found = next(
            (row for row in range(pivot_row, len(augmented)) if augmented[row][column]),
            None,
        )
        if found is None:
            continue
        augmented[pivot_row], augmented[found] = (
            augmented[found],
            augmented[pivot_row],
        )
        inverse = power(augmented[pivot_row][column], 6)
        augmented[pivot_row] = [
            multiply(inverse, value) for value in augmented[pivot_row]
        ]
        for row in range(len(augmented)):
            if row == pivot_row or not augmented[row][column]:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                left ^ multiply(scale, right)
                for left, right in zip(augmented[row], augmented[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
    if any(not any(row[:-1]) and row[-1] for row in augmented):
        return None, len(pivot_columns)
    solution = [0] * len(matrix[0])
    for row, column in enumerate(pivot_columns):
        solution[column] = augmented[row][-1]
    return solution, len(pivot_columns)


def transition_summary(rows: list[dict], side: int) -> dict:
    outputs: dict[tuple[int, int], set[int]] = collections.defaultdict(set)
    matrix = []
    rhs = []
    for row in rows:
        h1 = half_pair(row["differences"][0])
        assert h1 is not None
        second = quotient_parity(row["differences"][1])
        third = quotient_parity(row["differences"][2])
        outputs[(h1[side], second[side])].add(third[side])
        h_value, second_value = h1[side], second[side]
        matrix.append(
            [
                h_value,
                power(h_value, 2),
                power(h_value, 4),
                second_value,
                power(second_value, 2),
                power(second_value, 4),
            ]
        )
        rhs.append(third[side])
    solution, rank = solve_f8(matrix, rhs)
    return {
        "input_states": len(outputs),
        "ambiguous_input_states": sum(len(values) > 1 for values in outputs.values()),
        "maximum_outputs_per_state": max(map(len, outputs.values())),
        "frobenius_linear_rank": rank,
        "frobenius_linear_solution": solution,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--discovery-stop", type=int, default=5_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    rows = source["rows"]

    prefix_counts = collections.Counter(
        tuple(row["global_zero_differences"]) for row in rows
    )
    first_states: dict[tuple[int, int], list[bool]] = collections.defaultdict(list)
    paired_outputs: dict[tuple[int, ...], set[tuple[int, int]]] = (
        collections.defaultdict(set)
    )
    subfield_violations = []
    support_patterns: dict[str, set[tuple[int, ...]]] = {
        "discovery": set(),
        "holdout": set(),
        "all": set(),
    }
    deep_vectors = []
    for row in rows:
        differences = row["differences"]
        h1 = half_pair(differences[0])
        if h1 is None:
            raise AssertionError(f"p={row['p']}: first gate is not zero")
        second = quotient_parity(differences[1])
        third = quotient_parity(differences[2])
        second_zero = differences[1]["parity_zero"]
        first_states[h1].append(second_zero)
        paired_outputs[h1 + second].add(third)
        support = tuple(int(value != 0) for value in h1 + second + third)
        support_patterns["all"].add(support)
        support_patterns[
            "discovery" if row["p"] <= args.discovery_stop else "holdout"
        ].add(support)
        for difference_index, item in enumerate(differences, 1):
            lifted = half_pair(item)
            if lifted is not None and any(power(value, 8) != value for value in lifted):
                subfield_violations.append([row["p"], difference_index, *lifted])
        if second_zero:
            deep_vectors.append(vector_bits(third))

    affine_origin = deep_vectors[0]
    support_rows = {}
    for split, patterns in support_patterns.items():
        support_rows[split] = {
            "patterns": len(patterns),
            "polynomial_ranks": [
                support_polynomial_row(patterns, degree) for degree in range(1, 4)
            ],
        }
    result = {
        "input": str(args.input),
        "n_rows": len(rows),
        "prefix_pattern_counts": {
            ",".join("1" if value else "0" for value in pattern): count
            for pattern, count in sorted(prefix_counts.items())
        },
        "subfield_violations": subfield_violations,
        "first_half_states": len(first_states),
        "mixed_gate_two_states": sum(
            any(values) and not all(values) for values in first_states.values()
        ),
        "single_side_transitions": [
            transition_summary(rows, side) for side in range(2)
        ],
        "reciprocal_transition": {
            "input_states": len(paired_outputs),
            "ambiguous_input_states": sum(
                len(values) > 1 for values in paired_outputs.values()
            ),
            "maximum_outputs_per_state": max(map(len, paired_outputs.values())),
        },
        "deep_prefix": {
            "rows": len(deep_vectors),
            "third_linear_rank_over_f2": binary_rank(deep_vectors),
            "third_affine_rank_over_f2": binary_rank(
                [value ^ affine_origin for value in deep_vectors[1:]]
            ),
        },
        "support_bits": support_rows,
        "conclusion": (
            "mod-four half-lifts lie in F8, but neither first-state, paired-state, "
            "Frobenius-linear, affine-support, nor quadratic-support laws close gate three"
        ),
    }
    write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
