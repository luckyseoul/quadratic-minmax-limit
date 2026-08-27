#!/usr/bin/env python3
"""Split unknown deep p=7 orbits into all exact directional mean allocations."""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_size6_finite_deep_modular_cpsat import solve, source_hash  # noqa: E402


def exact_mean_allocations(direction_rows: list[dict]) -> tuple[dict[int, int], ...]:
    tables = [tuple(int(value) for value in row["allowed_scaled_means"]) for row in direction_rows]
    out = []
    for values in itertools.product(*tables):
        if all(
            sum(values[index] for index, row in enumerate(direction_rows) if int(row["eps"]) == eps) == 32
            for eps in (-1, 1)
        ):
            out.append({index: value for index, value in enumerate(values)})
    unique = tuple(dict(row) for row in {tuple(sorted(row.items())) for row in out})
    return tuple(sorted(unique, key=lambda row: tuple(row[index] for index in range(8))))


def run(
    source: Path,
    initial_shards: tuple[Path, ...],
    moduli: tuple[int, ...],
    shard_index: int,
    shard_count: int,
    seconds: float,
    workers: int,
    seed: int,
) -> dict:
    started = time.time()
    initial_rows = []
    source_digests = set()
    initial_indices = set()
    for path in initial_shards:
        payload = json.loads(path.read_text())
        source_digests.add(payload["source_sha256"])
        for row in payload["rows"]:
            orbit_index = int(row["orbit_index"])
            if orbit_index in initial_indices:
                raise ValueError("duplicate orbit in initial shards")
            initial_indices.add(orbit_index)
            if row["solver_status"] == "UNKNOWN":
                initial_rows.append(row)
    digest = source_hash(source)
    if source_digests != {digest}:
        raise ValueError("initial shards do not match the orbit source")
    selected = [
        row
        for row in initial_rows
        if int(row["orbit_index"]) % shard_count == shard_index
    ]
    orbit_rows = []
    leaf_status_counts = Counter()
    for initial in selected:
        orbit_index = int(initial["orbit_index"])
        allocations = exact_mean_allocations(initial["direction_rows"])
        leaves = []
        for allocation_index, fixed_means in enumerate(allocations):
            result = solve(
                source,
                orbit_index,
                moduli,
                seconds,
                workers,
                seed + 100 * orbit_index + allocation_index,
                fixed_means,
            )
            leaf_status_counts[result["solver_status"]] += 1
            leaves.append(result)
        orbit_rows.append(
            {
                "orbit_index": orbit_index,
                "allocation_count": len(allocations),
                "infeasible_allocations": sum(row["solver_status"] == "INFEASIBLE" for row in leaves),
                "feasible_allocations": sum(row["feasible"] for row in leaves),
                "unknown_allocations": sum(row["solver_status"] == "UNKNOWN" for row in leaves),
                "all_allocations_infeasible": bool(leaves) and all(
                    row["solver_status"] == "INFEASIBLE" for row in leaves
                ),
                "leaves": leaves,
            }
        )
    return {
        "experiment": "p7_size6_finite_deep_allocation_batch",
        "status": "complete_shard_exact_directional_mean_partition",
        "source": str(source),
        "source_sha256": digest,
        "moduli": list(moduli),
        "initial_shards": [str(path) for path in initial_shards],
        "initial_orbits": len(initial_indices),
        "initial_unknown_orbits": len(initial_rows),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "seconds_per_allocation": seconds,
        "workers_per_allocation": workers,
        "processed_unknown_orbits": len(orbit_rows),
        "processed_allocations": sum(row["allocation_count"] for row in orbit_rows),
        "closed_orbits": sum(row["all_allocations_infeasible"] for row in orbit_rows),
        "orbits_with_unknown_allocations": sum(row["unknown_allocations"] > 0 for row in orbit_rows),
        "orbits_with_feasible_allocations": sum(row["feasible_allocations"] > 0 for row in orbit_rows),
        "leaf_status_counts": dict(sorted(leaf_status_counts.items())),
        "elapsed_seconds": time.time() - started,
        "rows": orbit_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--initial-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 7))
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15664001)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run(
        args.source,
        tuple(args.initial_shards),
        tuple(args.moduli),
        args.shard_index,
        args.shard_count,
        args.seconds,
        args.workers,
        args.seed,
    )
    atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
