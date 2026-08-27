#!/usr/bin/env python3
"""Sharded compact modular-catalog sweep of all deep six-finite orbits."""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_size6_finite_deep_modular_cpsat import (  # noqa: E402
    load_source,
    solve,
    source_hash,
)


def run(
    source: Path,
    moduli: tuple[int, ...],
    shard_index: int,
    shard_count: int,
    seconds: float,
    workers: int,
    seed: int,
) -> dict:
    started = time.time()
    payload = load_source(source)
    deep = [
        index
        for index, orbit in enumerate(payload["orbits"])
        if any(int(value) not in (24, 32) for value in orbit["type_costs"].values())
    ]
    selected = [index for index in deep if index % shard_count == shard_index]
    rows = []
    for orbit_index in selected:
        rows.append(
            solve(
                source,
                orbit_index,
                moduli,
                seconds,
                workers,
                seed + orbit_index,
            )
        )
    counts = Counter(row["solver_status"] for row in rows)
    return {
        "experiment": "p7_size6_finite_deep_modular_batch",
        "status": "complete_shard_compact_high_mean_catalog_modular_models",
        "source": str(source),
        "source_sha256": source_hash(source),
        "moduli": list(moduli),
        "deep_orbits_in_source": len(deep),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "seconds_per_orbit": seconds,
        "workers_per_orbit": workers,
        "processed_orbits": len(rows),
        "solver_status_counts": dict(sorted(counts.items())),
        "infeasible_orbits": sum(row["finite_infeasibility_certificate"] for row in rows),
        "feasible_orbits": sum(row["feasible"] for row in rows),
        "unknown_orbits": sum(row["solver_status"] == "UNKNOWN" for row in rows),
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 7))
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=15663001)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    out = run(
        args.source,
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
