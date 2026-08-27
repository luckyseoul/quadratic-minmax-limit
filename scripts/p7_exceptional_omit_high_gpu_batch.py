#!/usr/bin/env python3
"""Resume-safe V100 batch for all exceptional high-mean leaves in one orbit."""
from __future__ import annotations

import argparse
import gc
import json
import math
import time
from pathlib import Path

import cupy as cp
import numpy as np

from p7_exceptional_omit_high_catalogs import OMITTED_DIRECTIONS, high_leaves
from p7_exceptional_projected_join_gpu import (
    atomic_json,
    choose_partition,
    run_mod7_pair,
)


def paths_for(
    projection_dir: Path, orbit_index: int, omitted_direction: int
) -> tuple[list[Path], list[Path]]:
    summaries = [
        projection_dir
        / f"p7_exceptional_omitd{omitted_direction}_p{profile}_all.json"
        for profile in range(3)
    ]
    catalogs = [
        projection_dir
        / (
            f"cminus_exceptional_omitd{omitted_direction}_p{profile}_"
            f"orbit{orbit_index:02d}.npz"
        )
        for profile in range(3)
    ]
    if any(not path.exists() for path in (*summaries, *catalogs)):
        raise FileNotFoundError(
            f"missing omission projections for direction {omitted_direction}"
        )
    return summaries, catalogs


def leaf_omitted_direction(batch: dict, leaf: dict) -> int | None:
    hits = [
        direction_index
        for direction_index, (row, mean) in enumerate(
            zip(batch["direction_rows"], leaf["scaled_means_direction_order"])
        )
        if int(row["b"]) == 0
        and int(row["phase"]) == 0
        and int(mean) > 16
    ]
    if not hits:
        return None
    if len(hits) != 1 or hits[0] not in OMITTED_DIRECTIONS:
        raise AssertionError("unexpected high-mean direction pattern")
    return hits[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--projection-dir", type=Path, required=True)
    parser.add_argument("--mean-batch", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-workload", type=int, default=20_000_000)
    parser.add_argument("--max-build", type=int, default=20_000_000)
    parser.add_argument("--chunk-size", type=int, default=5_000_000)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    started = time.time()
    batch = json.loads(args.mean_batch.read_text())
    if (
        batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
        or int(batch.get("allocation_count", -1)) != 180
    ):
        raise ValueError("unexpected exceptional mean batch")
    selected_by_direction = {
        direction: {
            int(leaf["leaf_index"]): leaf
            for leaf in high_leaves(batch, direction)
        }
        for direction in OMITTED_DIRECTIONS
    }
    all_selected = {
        leaf_index
        for leaves in selected_by_direction.values()
        for leaf_index in leaves
    }
    independently_selected = {
        int(leaf["leaf_index"])
        for leaf in batch["leaves"]
        if leaf_omitted_direction(batch, leaf) is not None
        and leaf.get("solver_status") != "INFEASIBLE"
    }
    if all_selected != independently_selected:
        raise AssertionError("high-leaf partition mismatch")

    projection_paths = {
        direction: paths_for(args.projection_dir, args.orbit_index, direction)
        for direction in OMITTED_DIRECTIONS
    }
    eligible = []
    deferred = []
    for leaf_index in sorted(all_selected):
        leaf = batch["leaves"][leaf_index]
        if int(leaf["leaf_index"]) != leaf_index:
            raise AssertionError("mean leaf index mismatch")
        omitted_direction = leaf_omitted_direction(batch, leaf)
        if omitted_direction is None:
            raise AssertionError("selected leaf has no high direction")
        _summaries, catalogs = projection_paths[omitted_direction]
        with np.load(catalogs[0], allow_pickle=False) as source:
            sizes = []
            for direction_index, mean in enumerate(
                leaf["scaled_means_direction_order"]
            ):
                key = f"d{direction_index}_m{mean}_p7_lo"
                if key not in source:
                    raise KeyError(f"projection catalog lacks {key}")
                count = int(len(source[key]))
                if count > 1:
                    sizes.append(count)
        first, second = choose_partition(tuple(sizes), args.max_build)
        build = math.prod(sizes[index] for index in first)
        probe = math.prod(sizes[index] for index in second)
        row = {
            "leaf_index": leaf_index,
            "omitted_direction": omitted_direction,
            "variable_catalog_sizes": sizes,
            "build": build,
            "probe": probe,
        }
        (eligible if max(build, probe) <= args.max_workload else deferred).append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for number, workload in enumerate(eligible, 1):
        leaf_index = workload["leaf_index"]
        omitted_direction = workload["omitted_direction"]
        summaries, catalogs = projection_paths[omitted_direction]
        output = args.output_dir / f"leaf{leaf_index:03d}.json"
        if output.exists():
            result = json.loads(output.read_text())
            if (
                int(result.get("leaf_index", -1)) != leaf_index
                or result.get("status")
                != "complete_exact_selected_dependency_gpu_join"
                or result.get("projection_mode")
                != "injective_disjoint_mod7_22x3_omitted_high_direction"
                or int(result.get("omitted_high_mean_direction", -1))
                != omitted_direction
            ):
                raise ValueError(f"invalid resume output {output}")
        else:
            result = run_mod7_pair(
                summaries[0],
                catalogs[0],
                summaries[1],
                catalogs[1],
                args.mean_batch,
                leaf_index,
                args.max_build,
                args.chunk_size,
                summaries[2],
                catalogs[2],
            )
            atomic_json(output, result)
        compact = {
            **workload,
            "matches": int(result["exact_projected_matches"]),
            "excluded": bool(result["finite_mean_allocation_exclusion"]),
            "elapsed_seconds": float(result["elapsed_seconds"]),
        }
        results.append(compact)
        print(
            json.dumps({"progress": [number, len(eligible)], **compact}),
            flush=True,
        )
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()

    out = {
        "experiment": "p7_exceptional_omit_high_gpu_batch",
        "status": "complete_exact_high_direction_omission_gpu_batch",
        "p": 7,
        "c_H": -1,
        "orbit_index": args.orbit_index,
        "fixed_boundary": batch["fixed_boundary"],
        "selected_high_leaf_count": len(all_selected),
        "eligible_leaf_count": len(eligible),
        "excluded_leaf_count": sum(row["excluded"] for row in results),
        "unresolved_leaf_count": sum(not row["excluded"] for row in results),
        "deferred_leaf_count": len(deferred),
        "max_workload": args.max_workload,
        "max_build": args.max_build,
        "results": results,
        "deferred": deferred,
        "elapsed_seconds": time.time() - started,
    }
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(
        json.dumps(
            {key: value for key, value in out.items() if key not in ("results", "deferred")},
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
