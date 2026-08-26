#!/usr/bin/env python3
"""Merge deterministic retry shards while auditing their exact coverage."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def key(row: list[int]) -> tuple[int, ...]:
    return tuple(int(value) for value in row)


def run(parent_source: Path, shard_sources: list[Path]) -> dict:
    parent = json.loads(parent_source.read_text())
    parent_unknown = [key(row) for row in parent["unknown_branches"]]
    if len(parent_unknown) != len(set(parent_unknown)):
        raise AssertionError("parent has duplicate unknown branch keys")
    parent_hash = sha256(parent_source)
    shards = [json.loads(path.read_text()) for path in shard_sources]
    shard_count = len(shards)
    if shard_count < 1:
        raise ValueError("at least one shard is required")
    indices = {int(shard["shard_index"]) for shard in shards}
    if indices != set(range(shard_count)):
        raise AssertionError("shard indices do not cover 0..count-1")

    rows = []
    shard_records = []
    covered: set[tuple[int, ...]] = set()
    for path, shard in zip(shard_sources, shards):
        if shard.get("experiment") != "p5_full_shell_circle_branch_retry_scip":
            raise AssertionError(f"wrong experiment in {path}")
        if int(shard["shard_count"]) != shard_count:
            raise AssertionError(f"wrong shard count in {path}")
        if int(shard["all_input_unknown_branch_count"]) != len(parent_unknown):
            raise AssertionError(f"wrong full input count in {path}")
        if shard.get("split_source_sha256") != parent_hash:
            raise AssertionError(f"parent hash mismatch in {path}")
        shard_index = int(shard["shard_index"])
        expected = {
            branch
            for position, branch in enumerate(parent_unknown)
            if position % shard_count == shard_index
        }
        actual = {key(row) for row in shard["input_unknown_branches"]}
        if actual != expected:
            raise AssertionError(f"shard partition mismatch in {path}")
        if covered & actual:
            raise AssertionError("shards overlap")
        covered |= actual
        rows.extend(shard["rows"])
        shard_records.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "shard_index": shard_index,
                "branch_count": int(shard["branch_count"]),
                "infeasible_count": len(shard["infeasible_branches"]),
                "unknown_count": len(shard["unknown_branches"]),
                "feasible_count": len(shard["feasible_branches"]),
            }
        )
    if covered != set(parent_unknown):
        raise AssertionError("merged shards do not cover the parent tail")
    if len(rows) != len(parent_unknown):
        raise AssertionError("merged row count does not equal parent tail")
    rows.sort(key=lambda row: key(row["source_branch_key"]))
    row_keys = [key(row["source_branch_key"]) for row in rows]
    if len(row_keys) != len(set(row_keys)) or set(row_keys) != covered:
        raise AssertionError("merged rows are not a bijection with parent branches")

    infeasible = [
        row["source_branch_key"] for row in rows if row["solver_status"] == "infeasible"
    ]
    feasible = [row["source_branch_key"] for row in rows if row.get("feasible") is True]
    unknown = [
        row["source_branch_key"]
        for row in rows
        if row["solver_status"] != "infeasible" and row.get("feasible") is not True
    ]
    first = shards[0]
    shard_records.sort(key=lambda row: int(row["shard_index"]))
    return {
        "experiment": "p5_full_shell_circle_branch_retry_scip",
        "status": "coverage_audited_merge_of_exact_retry_shards",
        "source": first["source"],
        "source_sha256": first["source_sha256"],
        "split_source": str(parent_source),
        "split_source_sha256": parent_hash,
        "split_experiment": parent["experiment"],
        "orbit_index": int(first["orbit_index"]),
        "boundary": first["boundary"],
        "c_H": int(first["c_H"]),
        "boundary_internal_edges": int(first["boundary_internal_edges"]),
        "all_input_unknown_branch_count": len(parent_unknown),
        "shard_index": 0,
        "shard_count": 1,
        "input_unknown_branches": [list(row) for row in parent_unknown],
        "branch_count": len(rows),
        "infeasible_branches": infeasible,
        "unknown_branches": unknown,
        "feasible_branches": feasible,
        "finite_infeasibility_certificate": len(infeasible) == len(rows),
        "merged_shards": shard_records,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent-source", type=Path, required=True)
    parser.add_argument("--shard-sources", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.parent_source, args.shard_sources)
    atomic_write(args.output, result)
    print(
        json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
