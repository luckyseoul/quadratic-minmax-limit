#!/usr/bin/env python3
"""Test 2-power lifts of the 56 order-21 F_8 support states.

The affine line R={1+t*sigma} is a relative difference set.  If n_i is its
intersection size with the i-th multiplicative coset modulo 42, then all
cyclic autocorrelations of n are known exactly.  Frobenius fixes n_i under
i -> 29*i, and the six index classes modulo six have known totals.  This
script asks whether each parity support admitted by the F_8 norm identities
can lift to counts modulo 2^b satisfying all those necessary conditions.

The constraints are finite b-bit bit-vector equations.  An UNSAT result is an
exact exclusion; SAT is only a surviving necessary-condition test.
"""
from __future__ import annotations

import argparse
import collections
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_f8_support_classify import (  # noqa: E402
    encoded_ratio,
    full_parity_word,
    support_data,
)


def bv_sum(terms: list[str], bits: int) -> str:
    if not terms:
        return f"#b{0:0{bits}b}"
    result = terms[0]
    for term in terms[1:]:
        result = f"(bvadd {result} {term})"
    return result


def target_bits(value: int, bits: int) -> str:
    return f"#b{value % (1 << bits):0{bits}b}"


def smt_problem(row: dict, k_residue: int, bits: int, timeout_ms: int) -> str:
    parity = full_parity_word(row["square_mask"], row["nonsquare_mask"])
    lines = [
        "(set-logic QF_BV)",
        f"(set-option :timeout {timeout_ms})",
    ]
    for index in range(42):
        lines.append(f"(declare-fun n{index} () (_ BitVec {bits}))")
        lines.append(
            f"(assert (= ((_ extract 0 0) n{index}) #b{parity[index]}))"
        )

    # Exact Frobenius invariance of the affine line in this congruence class.
    for index in range(42):
        lines.append(f"(assert (= n{index} n{(29 * index) % 42}))")

    # Six projective character classes.  N=(p+1)/6=14k+5; the omitted kernel
    # direction lies in residue class three.
    modulus = 1 << bits
    p = 84 * k_residue + 29
    n_modulus = ((p + 1) // 6) % modulus
    for residue in range(6):
        terms = [f"n{index}" for index in range(42) if index % 6 == residue]
        target = n_modulus - (residue == 3)
        lines.append(
            f"(assert (= {bv_sum(terms, bits)} {target_bits(target, bits)}))"
        )

    # Exact R R^- autocorrelation reduced modulo the requested 2-power.
    class_size = (p * p - 1) // 42
    scalar_intersection = (p - 1) // 7
    for shift in range(22):
        products = [
            f"(bvmul n{index} n{(index + shift) % 42})"
            for index in range(42)
        ]
        if shift == 0:
            target = p + class_size - scalar_intersection
        elif shift % 6 == 0:
            target = class_size - scalar_intersection
        else:
            target = class_size
        lines.append(
            f"(assert (= {bv_sum(products, bits)} {target_bits(target, bits)}))"
        )
    lines.extend(["(check-sat)", "(exit)"])
    return "\n".join(lines) + "\n"


def solve(z3: Path, problem: str) -> str:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".smt2") as handle:
        handle.write(problem)
        handle.flush()
        completed = subprocess.run(
            [str(z3), handle.name],
            check=False,
            capture_output=True,
            text=True,
        )
    answer = completed.stdout.strip().splitlines()
    if completed.returncode != 0 or not answer:
        raise RuntimeError(
            f"z3 failed rc={completed.returncode}: {completed.stderr.strip()}"
        )
    if answer[0] not in {"sat", "unsat", "unknown"}:
        raise RuntimeError(f"unexpected z3 output: {answer[0]}")
    return answer[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--z3", type=Path, required=True)
    parser.add_argument("--bits", type=int, default=3)
    parser.add_argument("--timeout-ms", type=int, default=10000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.z3.is_file() or args.timeout_ms < 1 or args.bits < 2:
        raise ValueError("--z3 must be a file; timeout/bits must be valid")

    odd_masks = [mask for mask in range(128) if mask.bit_count() % 2]
    solutions = []
    for square_mask in odd_masks:
        for nonsquare_mask in odd_masks:
            row = support_data(square_mask, nonsquare_mask)
            if row["diagonal"] == 1 and row["off_diagonal"] == 0:
                solutions.append(row)

    records = []
    k_period = 1 << (args.bits - 1)
    counts = {
        k_residue: collections.Counter() for k_residue in range(k_period)
    }
    for k_residue in range(k_period):
        for index, row in enumerate(solutions, 1):
            answer = solve(
                args.z3,
                smt_problem(row, k_residue, args.bits, args.timeout_ms),
            )
            counts[k_residue][(row["state_trace"], answer)] += 1
            records.append(
                {
                    "k_residue": k_residue,
                    "square_mask_hex": hex(row["square_mask"]),
                    "nonsquare_mask_hex": hex(row["nonsquare_mask"]),
                    "state_trace": row["state_trace"],
                    "ratio_pair": [
                        encoded_ratio(row["ratio"]),
                        encoded_ratio(row["ratio_star"]),
                    ],
                    "answer": answer,
                }
            )
        print(
            f"k_residue={k_residue} counts={dict(counts[k_residue])}",
            flush=True,
        )
    result = {
        "model": f"42 multiplicative-coset counts modulo {1 << args.bits}",
        "solver": str(args.z3),
        "bits": args.bits,
        "k_period": k_period,
        "timeout_ms": args.timeout_ms,
        "f8_states": len(solutions),
        "systems": len(records),
        "counts_by_k_residue_trace_answer": {
            str(k_residue): {
                f"trace={trace},{answer}": count
                for (trace, answer), count in sorted(counts[k_residue].items())
            }
            for k_residue in range(k_period)
        },
        "records": records,
    }
    write_json_atomic(args.output, result)


if __name__ == "__main__":
    main()
