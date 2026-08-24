#!/usr/bin/env python3
"""Finite F_8 classification forced by the order-21 generator identities.

For p == 29 (mod 84), Frobenius reduces the parity support of the Bose
generator modulo Phi_21 to two subsets S,T of Z/7.  The square component has
classes ``{0}``, S, S over the three residues modulo 3; the nonsquare
component has classes T, empty, T.  The already-proved Phi_3 generator value
forces both |S| and |T| odd.  This script exhausts the remaining 64^2 support
pairs against the diagonal and off-diagonal parts of

    Gamma Gamma^(-1) = e + G + N.

The computation is a finite exact lemma-discovery aid.  Its output records
whether those identities alone force the observed four-state ratio law or
whether another structural input is still required.
"""
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

from io_atomic import write_json_atomic


MODULUS = 0x57
ALPHA = 0x2


def multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> 6:
            left ^= MODULUS
    return result


def power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


LAMBDA = power(ALPHA, 15)  # order 7, CRT class (1 mod 7, 0 mod 3)
U = power(ALPHA, 3)  # order 7, and ALPHA = LAMBDA * ALPHA^7
H = power(U, 4)  # unique square root of U in F_8; the forbidden finite ratio

# Local factor X^3+X+1 of Phi_7 modulo two.  This second F_8 model is used
# for the phase-sensitive order-14 Jacobi quotient.
JACOBI_F8_MODULUS = 0xB
JACOBI_F8_ROOT = 0x2


def trace(value: int) -> int:
    return value ^ power(value, 2) ^ power(value, 4)


def jacobi_f8_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> 3:
            left ^= JACOBI_F8_MODULUS
    return result


def jacobi_f8_power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = jacobi_f8_multiply(result, value)
        value = jacobi_f8_multiply(value, value)
        exponent >>= 1
    return result


def jacobi_f8_trace(value: int) -> int:
    square = jacobi_f8_multiply(value, value)
    fourth = jacobi_f8_multiply(square, square)
    return value ^ square ^ fourth


def jacobi_phase_trace(square_mask: int, nonsquare_mask: int) -> int:
    """First 2-adic trace of the order-14/order-7 Jacobi quotient."""
    word = [0] * 14
    for exponent in range(7):
        word[2 * exponent] = (
            ((square_mask >> exponent) & 1) ^ int(exponent == 0)
        )
        word[2 * exponent + 1] = (nonsquare_mask >> exponent) & 1
    collapsed = [word[index] ^ word[index + 7] for index in range(7)]
    if sum(collapsed) != 1:
        return -1
    support = collapsed.index(1)
    nonsquare_value = 0
    for exponent in range(7):
        if word[2 * exponent + 1]:
            nonsquare_value ^= jacobi_f8_power(JACOBI_F8_ROOT, exponent)
    shifted = jacobi_f8_multiply(
        jacobi_f8_power(JACOBI_F8_ROOT, (4 - 4 * support) % 7),
        nonsquare_value,
    )
    return jacobi_f8_trace(shifted)


def subset_value(mask: int, inverse: bool = False) -> int:
    value = 0
    for exponent in range(7):
        if mask & (1 << exponent):
            value ^= power(LAMBDA, (-exponent if inverse else exponent) % 7)
    return value


def reciprocal_to_common(value: int) -> int:
    alpha_inverse = power(ALPHA, 62)
    result = 0
    term = 1
    for exponent in range(6):
        if value & (1 << exponent):
            result ^= term
        term = multiply(term, alpha_inverse)
    return result


def support_data(square_mask: int, nonsquare_mask: int) -> dict:
    p_value = subset_value(square_mask)
    p_star = subset_value(square_mask, True)
    q_value = subset_value(nonsquare_mask)
    q_star = subset_value(nonsquare_mask, True)
    square = 1 ^ p_value
    square_star = 1 ^ p_star
    nonsquare_over_x = multiply(power(LAMBDA, 6), q_value)
    nonsquare_over_x_star = multiply(LAMBDA, q_star)
    diagonal = (
        multiply(square, square_star)
        ^ multiply(nonsquare_over_x, nonsquare_over_x_star)
    )
    off_diagonal = (
        multiply(square, nonsquare_over_x_star)
        ^ multiply(U, multiply(nonsquare_over_x, square_star))
    )
    ratio = (
        None
        if nonsquare_over_x == 0
        else multiply(square, power(nonsquare_over_x, 6))
    )
    ratio_star = (
        None
        if nonsquare_over_x_star == 0
        else multiply(square_star, power(nonsquare_over_x_star, 6))
    )
    # The diagonal norm excludes ratio H whenever the denominator is nonzero,
    # so this Mobius coordinate is defined on all solutions.  It sends the
    # denominator-zero state to zero and the seven finite states to F_8^*.
    state_parameter = (
        0
        if nonsquare_over_x == 0
        else multiply(
            nonsquare_over_x,
            power(square ^ multiply(H, nonsquare_over_x), 6),
        )
    )
    return {
        "square_mask": square_mask,
        "nonsquare_mask": nonsquare_mask,
        "square_weight": square_mask.bit_count(),
        "nonsquare_weight": nonsquare_mask.bit_count(),
        "square": square,
        "square_star": square_star,
        "nonsquare_over_x": nonsquare_over_x,
        "nonsquare_over_x_star": nonsquare_over_x_star,
        "diagonal": diagonal,
        "off_diagonal": off_diagonal,
        "ratio": ratio,
        "ratio_star": ratio_star,
        "state_parameter": state_parameter,
        "state_trace": trace(multiply(H, state_parameter)),
        "jacobi_phase_trace": jacobi_phase_trace(square_mask, nonsquare_mask),
    }


def encoded_ratio(value: int | None) -> str:
    return "none" if value is None else hex(value)


def solve_binary(matrix: list[list[int]], rhs: list[int]) -> list[int] | None:
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    n_columns = len(matrix[0])
    pivot_row = 0
    pivots = []
    for column in range(n_columns):
        found = next(
            (
                row
                for row in range(pivot_row, len(augmented))
                if augmented[row][column]
            ),
            None,
        )
        if found is None:
            continue
        augmented[pivot_row], augmented[found] = (
            augmented[found],
            augmented[pivot_row],
        )
        for row in range(len(augmented)):
            if row != pivot_row and augmented[row][column]:
                augmented[row] = [
                    left ^ right
                    for left, right in zip(augmented[row], augmented[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
    if any(not any(row[:n_columns]) and row[-1] for row in augmented):
        return None
    solution = [0] * n_columns
    for row, column in enumerate(pivots):
        solution[column] = augmented[row][-1]
    return solution


def full_parity_word(square_mask: int, nonsquare_mask: int) -> list[int]:
    """Expand the Frobenius support form to the 42 multiplicative cosets."""
    word = [0] * 42
    word[0] = 1
    for exponent in range(7):
        if square_mask & (1 << exponent):
            for mod3 in (1, 2):
                order21 = (15 * exponent + 7 * mod3) % 21
                word[2 * order21] = 1
        if nonsquare_mask & (1 << exponent):
            for mod3 in (0, 2):
                order21 = (15 * exponent + 7 * mod3) % 21
                word[2 * order21 + 1] = 1
    return word


def mod4_relative_difference_set_lift(row: dict, k_parity: int) -> bool:
    """Test the exact RDS autocorrelation and six class totals modulo four.

    Write p=84k+29 and each integer cyclotomic count as v_i+2*y_i
    modulo four.  The parity word v is fixed by S,Q; all remaining conditions
    are linear in the 42 unknown bits y_i.
    """
    parity = full_parity_word(row["square_mask"], row["nonsquare_mask"])
    matrix = []
    rhs = []
    # For this congruence class, |C_i|=(p^2-1)/42 and |C_i cap F_p^*|
    # are both 0 mod 4.  R R^- therefore has autocorrelation 1 at shift
    # zero and 0 at every nonzero shift modulo four.
    if sum(parity) % 4 != 1:
        return False
    for shift in range(1, 42):
        base = sum(
            parity[index] * parity[(index + shift) % 42]
            for index in range(42)
        )
        difference = (-base) % 4
        if difference % 2:
            return False
        matrix.append(
            [
                parity[(index + shift) % 42]
                ^ parity[(index - shift) % 42]
                for index in range(42)
            ]
        )
        rhs.append(difference // 2)

    # The six index classes modulo six have totals N=(p+1)/6, except class
    # three (the nonsquare Frobenius-fixed class), whose total is N-1.
    n_mod4 = 1 if k_parity == 0 else 3
    for residue in range(6):
        indices = [index for index in range(42) if index % 6 == residue]
        target = (n_mod4 - (residue == 3)) % 4
        base = sum(parity[index] for index in indices)
        difference = (target - base) % 4
        if difference % 2:
            return False
        matrix.append([int(index in indices) for index in range(42)])
        rhs.append(difference // 2)
    return solve_binary(matrix, rhs) is not None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--generator-evidence", type=Path)
    args = parser.parse_args()

    odd_masks = [mask for mask in range(128) if mask.bit_count() % 2]
    solutions = []
    for square_mask, nonsquare_mask in itertools.product(odd_masks, repeat=2):
        row = support_data(square_mask, nonsquare_mask)
        if row["diagonal"] == 1 and row["off_diagonal"] == 0:
            solutions.append(row)
    jacobi_phase_equals_state_trace = all(
        row["jacobi_phase_trace"] == row["state_trace"] for row in solutions
    )
    if not jacobi_phase_equals_state_trace:
        raise AssertionError("Jacobi local trace and Mobius trace disagree")

    ratio_counts = collections.Counter(
        (row["ratio"], row["ratio_star"]) for row in solutions
    )
    weight_counts = collections.Counter(
        (row["square_weight"], row["nonsquare_weight"]) for row in solutions
    )
    weight_counts_by_trace = {
        trace_value: collections.Counter(
            (row["square_weight"], row["nonsquare_weight"])
            for row in solutions
            if row["state_trace"] == trace_value
        )
        for trace_value in (0, 1)
    }
    nonsquare_zero = [
        row
        for row in solutions
        if row["nonsquare_over_x"] == 0
        or row["nonsquare_over_x_star"] == 0
    ]
    observed_left = {0x0, 0x8, 0x9, 0x1F}
    observed_pairs = {(value, multiply(value, power(U, 6))) for value in observed_left}
    outside_observed = [
        row
        for row in solutions
        if (row["ratio"], row["ratio_star"]) not in observed_pairs
    ]
    trace_counts = collections.Counter(row["state_trace"] for row in solutions)
    jacobi_phase_trace_counts = collections.Counter(
        row["jacobi_phase_trace"] for row in solutions
    )
    ratio_counts_by_trace = {
        trace_value: collections.Counter(
            (row["ratio"], row["ratio_star"])
            for row in solutions
            if row["state_trace"] == trace_value
        )
        for trace_value in (0, 1)
    }
    support_labels = ["constant"] + [
        f"S_{index}" for index in range(7)
    ] + [f"T_{index}" for index in range(7)]
    support_matrix = [
        [1]
        + [(row["square_mask"] >> index) & 1 for index in range(7)]
        + [(row["nonsquare_mask"] >> index) & 1 for index in range(7)]
        for row in solutions
    ]
    support_trace_linear = solve_binary(
        support_matrix, [row["state_trace"] for row in solutions]
    )
    mod4_lift_counts = {
        k_parity: collections.Counter(
            row["state_trace"]
            for row in solutions
            if mod4_relative_difference_set_lift(row, k_parity)
        )
        for k_parity in (0, 1)
    }
    mod4_lift_ratio_counts = {
        k_parity: collections.Counter(
            (row["ratio"], row["ratio_star"])
            for row in solutions
            if mod4_relative_difference_set_lift(row, k_parity)
        )
        for k_parity in (0, 1)
    }
    generator_validation = None
    if args.generator_evidence:
        generator_rows = json.loads(args.generator_evidence.read_text())["rows"]
        validation_rows = []
        failures = []
        for generator_row in generator_rows:
            support = generator_row["component_support"]
            square_groups = [
                {exponent % 7 for exponent in support[0] if exponent % 3 == cls}
                for cls in range(3)
            ]
            nonsquare_groups = [
                {exponent % 7 for exponent in support[1] if exponent % 3 == cls}
                for cls in range(3)
            ]
            form_ok = (
                square_groups[0] == {0}
                and square_groups[1] == square_groups[2]
                and not nonsquare_groups[1]
                and nonsquare_groups[0] == nonsquare_groups[2]
            )
            square_mask = sum(1 << value for value in square_groups[1])
            nonsquare_mask = sum(1 << value for value in nonsquare_groups[0])
            row = support_data(square_mask, nonsquare_mask)
            direct_left = int(
                generator_row["factors"][0]["ratio_hex"], 16
            )
            direct_right = reciprocal_to_common(
                int(generator_row["factors"][1]["ratio_hex"], 16)
            )
            checks = {
                "frobenius_support_form": form_ok,
                "odd_support_weights": (
                    square_mask.bit_count() % 2 == 1
                    and nonsquare_mask.bit_count() % 2 == 1
                ),
                "diagonal_norm": row["diagonal"] == 1,
                "off_diagonal_norm": row["off_diagonal"] == 0,
                "left_ratio": row["ratio"] == direct_left,
                "right_ratio": row["ratio_star"] == direct_right,
                "trace_one": row["state_trace"] == 1,
                "jacobi_phase_trace_one": row["jacobi_phase_trace"] == 1,
            }
            validation_rows.append(
                {
                    "p": generator_row["p"],
                    "square_mask_hex": hex(square_mask),
                    "nonsquare_mask_hex": hex(nonsquare_mask),
                    "ratio_pair_common_hex": [hex(direct_left), hex(direct_right)],
                    "state_parameter_hex": hex(row["state_parameter"]),
                    "state_trace": row["state_trace"],
                    "jacobi_phase_trace": row["jacobi_phase_trace"],
                    "checks": checks,
                }
            )
            if not all(checks.values()):
                failures.append(validation_rows[-1])
        generator_validation = {
            "path": str(args.generator_evidence),
            "rows": len(validation_rows),
            "failures": failures,
            "records": validation_rows,
        }
    result = {
        "field_modulus_hex": hex(MODULUS),
        "alpha_hex": hex(ALPHA),
        "lambda_hex": hex(LAMBDA),
        "u_equals_alpha_cubed_hex": hex(U),
        "h_equals_sqrt_u_hex": hex(H),
        "assumptions": [
            "square support by mod-3 class is ({0}, S, S)",
            "nonsquare support by mod-3 class is (T, empty, T)",
            "|S| and |T| are odd from the Phi_3 generator value",
            "diagonal generator norm is 1",
            "off-diagonal generator norm is 0",
        ],
        "odd_masks": len(odd_masks),
        "candidate_pairs": len(odd_masks) ** 2,
        "solutions": len(solutions),
        "ratio_pair_counts": {
            f"{encoded_ratio(left)}|{encoded_ratio(right)}": count
            for (left, right), count in sorted(ratio_counts.items(), key=str)
        },
        "state_trace_counts": dict(sorted(trace_counts.items())),
        "jacobi_phase_trace_counts": dict(
            sorted(jacobi_phase_trace_counts.items())
        ),
        "jacobi_phase_equals_state_trace": jacobi_phase_equals_state_trace,
        "state_trace_linear_in_support_bits": (
            None
            if support_trace_linear is None
            else [
                label
                for label, coefficient in zip(
                    support_labels, support_trace_linear
                )
                if coefficient
            ]
        ),
        "mod4_rds_lift_counts_by_k_parity_and_trace": {
            str(k_parity): {
                str(trace_value): count
                for trace_value, count in sorted(mod4_lift_counts[k_parity].items())
            }
            for k_parity in (0, 1)
        },
        "mod4_rds_lift_ratio_counts_by_k_parity": {
            str(k_parity): {
                f"{encoded_ratio(left)}|{encoded_ratio(right)}": count
                for (left, right), count in sorted(
                    mod4_lift_ratio_counts[k_parity].items(), key=str
                )
            }
            for k_parity in (0, 1)
        },
        "ratio_pair_counts_by_state_trace": {
            str(trace_value): {
                f"{encoded_ratio(left)}|{encoded_ratio(right)}": count
                for (left, right), count in sorted(
                    ratio_counts_by_trace[trace_value].items(), key=str
                )
            }
            for trace_value in (0, 1)
        },
        "weight_pair_counts": {
            f"{left}|{right}": count
            for (left, right), count in sorted(weight_counts.items())
        },
        "weight_pair_counts_by_state_trace": {
            str(trace_value): {
                f"{left}|{right}": count
                for (left, right), count in sorted(
                    weight_counts_by_trace[trace_value].items()
                )
            }
            for trace_value in (0, 1)
        },
        "nonsquare_zero_solutions": len(nonsquare_zero),
        "outside_observed_four_states": len(outside_observed),
        "nonsquare_zero_examples": nonsquare_zero[:12],
        "outside_observed_examples": outside_observed[:12],
        "generator_validation": generator_validation,
    }
    print(json.dumps(result, indent=2))
    if args.output:
        write_json_atomic(args.output, result)


if __name__ == "__main__":
    main()
