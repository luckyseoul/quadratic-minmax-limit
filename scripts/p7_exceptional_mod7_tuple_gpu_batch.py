#!/usr/bin/env python3
"""Resume-safe V100 batch for tractable exceptional mean leaves."""
from __future__ import annotations

import argparse
import gc
import json
import math
import time
from pathlib import Path

import cupy as cp
import numpy as np

from p7_exceptional_projected_join_gpu import atomic_json, choose_partition, run_mod7_pair


def main() -> None:
    parser = argparse.ArgumentParser()
    for index in (1, 2, 3):
        suffix = "" if index == 1 else f"-{index}"
        parser.add_argument(f"--projection-summary{suffix}", type=Path, required=True)
        parser.add_argument(f"--projection-catalog{suffix}", type=Path, required=True)
    parser.add_argument("--mean-batch", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-workload", type=int, default=200_000_000)
    parser.add_argument("--max-build", type=int, default=150_000_000)
    parser.add_argument("--chunk-size", type=int, default=20_000_000)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    started = time.time()
    summaries = [args.projection_summary, args.projection_summary_2, args.projection_summary_3]
    catalogs = [args.projection_catalog, args.projection_catalog_2, args.projection_catalog_3]
    batch = json.loads(args.mean_batch.read_text())
    args.output_dir.mkdir(parents=True, exist_ok=True)

    eligible = []
    deferred = []
    with np.load(catalogs[0], allow_pickle=False) as source:
        for leaf in batch["leaves"]:
            if leaf.get("solver_status") == "INFEASIBLE":
                continue
            sizes = []
            supported = True
            for direction_index, mean in enumerate(leaf["scaled_means_direction_order"]):
                key = f"d{direction_index}_m{mean}_p7_lo"
                if key not in source:
                    supported = False
                    break
                count = len(source[key])
                if count > 1:
                    sizes.append(count)
            if not supported:
                continue
            first, second = choose_partition(tuple(sizes), args.max_build)
            build = math.prod(sizes[index] for index in first)
            probe = math.prod(sizes[index] for index in second)
            row = {"leaf_index": int(leaf["leaf_index"]), "build": build, "probe": probe}
            (eligible if max(build, probe) <= args.max_workload else deferred).append(row)

    results = []
    for number, workload in enumerate(eligible, 1):
        leaf_index = workload["leaf_index"]
        output = args.output_dir / f"leaf{leaf_index:03d}.json"
        if output.exists():
            result = json.loads(output.read_text())
            if result.get("leaf_index") != leaf_index or result.get("status") != "complete_exact_selected_dependency_gpu_join":
                raise ValueError(f"invalid resume output {output}")
        else:
            result = run_mod7_pair(
                summaries[0], catalogs[0], summaries[1], catalogs[1],
                args.mean_batch, leaf_index, args.max_build, args.chunk_size,
                summaries[2], catalogs[2],
            )
            atomic_json(output, result)
        compact = {
            "leaf_index": leaf_index,
            "matches": int(result["exact_projected_matches"]),
            "excluded": bool(result["finite_mean_allocation_exclusion"]),
            "elapsed_seconds": float(result["elapsed_seconds"]),
        }
        results.append(compact)
        print(json.dumps({"progress": [number, len(eligible)], **compact}), flush=True)
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()

    out = {
        "experiment": "p7_exceptional_mod7_tuple_gpu_batch",
        "status": "complete_exact_selected_dependency_gpu_batch",
        "fixed_boundary": batch["fixed_boundary"],
        "eligible_leaf_count": len(eligible),
        "excluded_leaf_count": sum(row["excluded"] for row in results),
        "unresolved_leaf_count": sum(not row["excluded"] for row in results),
        "deferred_supported_leaf_count": len(deferred),
        "max_workload": args.max_workload,
        "results": results,
        "deferred": deferred,
        "elapsed_seconds": time.time() - started,
    }
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(json.dumps({key: value for key, value in out.items() if key not in ("results", "deferred")}, indent=2))


if __name__ == "__main__":
    main()
