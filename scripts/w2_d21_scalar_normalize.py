#!/usr/bin/env python3
"""Normalize d=21 boundary residues by the proved nonzero Bose generator.

The signature census stores raw square/nonsquare residues at the reciprocal
sextics of Phi_21.  The local Jacobi-norm lemma proves that the nonsquare
generator residue is nonzero at both factors, so division recovers the actual
three boundary scalars instead of only their projective class.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_generator_ratio import FACTORS, generator_record  # noqa: E402
from w2_d21_signature_scan import multiply_bits, power_bits  # noqa: E402
from w2_d21_projective_invariants import reciprocal_to_common  # noqa: E402


def normalize_row(row: dict) -> dict:
    generator = generator_record(row["p"])
    triples = []
    checks = []
    for factor_index, factor in enumerate(FACTORS):
        generator_factor = generator["factors"][factor_index]
        generator_square = int(generator_factor["square_remainder_hex"], 16)
        generator_nonsquare = int(
            generator_factor["nonsquare_remainder_hex"], 16
        )
        if not generator_nonsquare:
            raise AssertionError(
                f"p={row['p']} factor={factor:#x}: zero generator denominator"
            )
        inverse = power_bits(generator_nonsquare, 62, factor)
        values = []
        for difference_index, (square, nonsquare) in enumerate(
            row["remainders"][factor_index], 1
        ):
            scalar = multiply_bits(nonsquare, inverse, factor)
            reconstructed_square = multiply_bits(
                scalar, generator_square, factor
            )
            checks.append(
                {
                    "factor_index": factor_index,
                    "difference": difference_index,
                    "square_reconstruction": reconstructed_square == square,
                    "scalar_in_f8": power_bits(scalar, 8, factor) == scalar,
                }
            )
            values.append(scalar)
        if factor_index == 1:
            values = [reciprocal_to_common(value) for value in values]
        triples.append(values)
    return {
        "p": row["p"],
        "masks": row["masks"],
        "generator_nonsquare_remainders": [
            int(item["nonsquare_remainder_hex"], 16)
            for item in generator["factors"]
        ],
        "scalar_triples_common_f8": triples,
        "checks": checks,
        "all_checks": all(
            item["square_reconstruction"] and item["scalar_in_f8"]
            for item in checks
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    source = json.loads(args.input.read_text())
    input_rows = source["retained"]
    started = time.perf_counter()
    rows = []
    if args.workers == 1:
        for index, row in enumerate(input_rows, 1):
            rows.append(normalize_row(row))
            if index % 25 == 0 or index == len(input_rows):
                print(f"[{index}/{len(input_rows)}]", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(normalize_row, row): row["p"] for row in input_rows
            }
            for index, future in enumerate(as_completed(futures), 1):
                rows.append(future.result())
                if index % 25 == 0 or index == len(input_rows):
                    print(f"[{index}/{len(input_rows)}]", flush=True)
    rows.sort(key=lambda row: row["p"])
    pair_counts = collections.Counter(
        tuple(tuple(values) for values in row["scalar_triples_common_f8"])
        for row in rows
    )
    result = {
        "input": str(args.input),
        "description": "affine d=21 boundary scalars after generator division",
        "n_rows": len(rows),
        "all_checks": all(row["all_checks"] for row in rows),
        "distinct_scalar_pairs": len(pair_counts),
        "scalar_pair_counts": {
            "|".join(",".join(hex(value) for value in triple) for triple in pair): count
            for pair, count in sorted(pair_counts.items())
        },
        "simultaneous_zero_triples": [
            row
            for row in rows
            if row["scalar_triples_common_f8"] == [[0, 0, 0], [0, 0, 0]]
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done rows={len(rows)} pairs={len(pair_counts)} "
        f"checks={result['all_checks']} "
        f"failures={len(result['simultaneous_zero_triples'])}",
        flush=True,
    )


if __name__ == "__main__":
    main()
