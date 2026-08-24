#!/usr/bin/env python3
"""Cross-check sparse low-order W2 records against dense exact certificates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from w2_low_order_atomic_gpu import cyclotomic_bits
from w2_translated_antipodal_norm_scan import f2_gcd_bits


def load_rows(paths: list[Path]) -> dict[int, dict]:
    rows: dict[int, dict] = {}
    for path in paths:
        with path.open() as source:
            document = json.load(source)
        for row in document["rows"]:
            p = int(row["p"])
            if p in rows:
                raise ValueError(f"duplicate p={p} in inputs")
            rows[p] = row
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dense", type=Path, nargs="+", required=True)
    parser.add_argument("--sparse", type=Path, nargs="+", required=True)
    args = parser.parse_args()

    dense_rows = load_rows(args.dense)
    sparse_rows = load_rows(args.sparse)
    if dense_rows.keys() != sparse_rows.keys():
        missing_dense = sorted(sparse_rows.keys() - dense_rows.keys())
        missing_sparse = sorted(dense_rows.keys() - sparse_rows.keys())
        raise AssertionError(
            f"prime sets differ: missing_dense={missing_dense}, "
            f"missing_sparse={missing_sparse}"
        )

    comparisons = 0
    exceptional = 0
    for p in sorted(dense_rows):
        dense = dense_rows[p]
        sparse = sparse_rows[p]
        expected: dict[tuple[int, int], int] = {}
        for endpoint in dense["endpoints"]:
            endpoint_offset = int(endpoint["offset"])
            if endpoint_offset not in sparse["endpoint_offsets"]:
                continue
            exact_bad = int(endpoint["aut_bad_hex"], 16)
            for order in sparse["orders_tested"]:
                factor = f2_gcd_bits(exact_bad, cyclotomic_bits(order))
                comparisons += 1
                if factor != 1:
                    expected[(endpoint_offset, order)] = factor

        observed = {
            (int(item["endpoint_offset"]), int(item["order"])):
                int(item["aut_bad_hex"], 16)
            for item in sparse["exceptions"]
        }
        if observed != expected:
            raise AssertionError(
                f"p={p} mismatch: expected={expected}, observed={observed}"
            )
        exceptional += len(observed)

    print(
        f"exact match: primes={len(dense_rows)} "
        f"endpoint_order_comparisons={comparisons} exceptions={exceptional}"
    )


if __name__ == "__main__":
    main()
