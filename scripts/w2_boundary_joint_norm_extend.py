#!/usr/bin/env python3
"""Recompute truncated joint-norm exceptions through a full line prefix."""
from __future__ import annotations

import argparse
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
from w2_boundary_line_direct import record  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--prefix", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prefix < 1 or args.workers < 1:
        raise ValueError("prefix and workers must be positive")

    source = json.loads(args.input.read_bytes())
    orders_by_prime: dict[int, list[int]] = {}
    for exception in source["exceptions"]:
        orders_by_prime.setdefault(exception["p"], []).append(exception["order"])
    orders_by_prime = {
        p: sorted(set(orders)) for p, orders in orders_by_prime.items()
    }

    started = time.perf_counter()
    rows = []
    if args.workers == 1:
        for index, (p, orders) in enumerate(orders_by_prime.items(), 1):
            row = record(p, orders, args.prefix, False, True)
            rows.append(row)
            print(f"[{index}/{len(orders_by_prime)}] p={p}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(record, p, orders, args.prefix, False, True): p
                for p, orders in orders_by_prime.items()
            }
            for index, future in enumerate(as_completed(futures), 1):
                row = future.result()
                rows.append(row)
                print(
                    f"[{index}/{len(orders_by_prime)}] p={row['p']}",
                    flush=True,
                )

    rows.sort(key=lambda row: row["p"])
    result = {
        "source": str(args.input),
        "prefix": args.prefix,
        "trace_full_prefix": True,
        "pairs": sum(len(orders) for orders in orders_by_prime.values()),
        "n_primes": len(rows),
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} pairs={result['pairs']} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
