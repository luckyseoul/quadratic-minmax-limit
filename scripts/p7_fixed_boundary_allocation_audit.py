#!/usr/bin/env python3
"""Audit complete exact mean-allocation coverage for one p=7 boundary.

The older CP-SAT batch may contain a superset of the admissible leaves.  This
audit independently reconstructs the exact type-sum and common-residue leaves,
then accepts a leaf only when either the CP-SAT batch certifies it INFEASIBLE or
an exact catalog join has zero simultaneous modular tuples.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def audit(batch_path: Path, join_paths: tuple[Path, ...]) -> dict:
    batch = json.loads(batch_path.read_text())
    if batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch":
        raise ValueError("unexpected CP-SAT batch experiment")
    if batch.get("status") != "complete_exact_mean_allocation_exhaustion":
        raise ValueError("the CP-SAT batch did not finish its enumerated leaves")
    if int(batch.get("p", 0)) != 7 or int(batch.get("c_H", 0)) not in (-1, 1):
        raise ValueError("unexpected CP-SAT batch scope")

    c_h = int(batch["c_H"])
    boundary = tuple(int(value) for value in batch["fixed_boundary"])
    rows = direction_data(c_h, boundary)
    exact_leaves = allocations(rows)
    if not exact_leaves or len(set(exact_leaves)) != len(exact_leaves):
        raise AssertionError("independent exact allocation enumeration failed")

    batch_by_means: dict[tuple[int, ...], dict] = {}
    for row in batch.get("leaves", []):
        means = tuple(int(value) for value in row["scaled_means_direction_order"])
        if means in batch_by_means:
            raise AssertionError("CP-SAT batch contains duplicate mean tuples")
        batch_by_means[means] = row

    joins_by_means: dict[tuple[int, ...], tuple[dict, Path]] = {}
    join_sources = []
    for path in join_paths:
        payload = json.loads(path.read_text())
        means = tuple(int(value) for value in payload.get("fixed_scaled_means", []))
        if means in joins_by_means:
            raise AssertionError("duplicate exact catalog join for one mean tuple")
        if (
            payload.get("experiment") != "p7_fixed_boundary_catalog_join"
            or payload.get("status")
            != "complete_exact_multimodular_catalog_join"
            or int(payload.get("p", 0)) != 7
            or int(payload.get("c_H", 0)) != c_h
            or tuple(int(value) for value in payload.get("fixed_boundary", []))
            != boundary
        ):
            raise ValueError(f"catalog join has the wrong scope: {path}")
        joins_by_means[means] = (payload, path)
        join_sources.append({"path": str(path), "sha256": sha256(path)})

    coverage = []
    missing = []
    malformed = []
    method_counts = {"cp_sat": 0, "catalog_join": 0}
    for leaf_index, means in enumerate(exact_leaves):
        batch_row = batch_by_means.get(means)
        batch_proves = bool(
            batch_row is not None
            and batch_row.get("solver_status") == "INFEASIBLE"
            and batch_row.get("finite_infeasibility_certificate") is True
            and batch_row.get("feasible") is False
        )
        join_item = joins_by_means.get(means)
        join_proves = bool(
            join_item is not None
            and int(join_item[0].get("consistent_catalog_tuples", -1)) == 0
            and join_item[0].get("modularly_infeasible") is True
            and join_item[0].get("finite_mean_allocation_exclusion") is True
        )
        if batch_proves:
            method = "cp_sat"
        elif join_proves:
            method = "catalog_join"
        else:
            method = "missing"
            missing.append(list(means))
        if method != "missing":
            method_counts[method] += 1
        coverage.append(
            {
                "exact_leaf_index": leaf_index,
                "scaled_means_direction_order": list(means),
                "method": method,
                "batch_leaf_index": (
                    int(batch_row["leaf_index"]) if batch_row is not None else None
                ),
                "join_path": str(join_item[1]) if join_item is not None else None,
            }
        )

    exact_leaf_set = set(exact_leaves)
    for means, (payload, path) in joins_by_means.items():
        if means not in exact_leaf_set:
            malformed.append(f"join is not an admissible exact leaf: {path}")
        if not (
            int(payload.get("consistent_catalog_tuples", -1)) == 0
            and payload.get("modularly_infeasible") is True
            and payload.get("finite_mean_allocation_exclusion") is True
        ):
            malformed.append(f"join is not an exclusion: {path}")

    checks = {
        "exact_allocations_pairwise_distinct": len(set(exact_leaves))
        == len(exact_leaves),
        "exact_allocations_have_type_sum_32": all(
            sum(value for value, row in zip(means, rows) if row["eps"] == eps)
            == 32
            for means in exact_leaves
            for eps in (-1, 1)
        ),
        "exact_allocations_have_common_type_residue_mod_8": all(
            len(
                {
                    value % 8
                    for value, row in zip(means, rows)
                    if row["eps"] == eps
                }
            )
            == 1
            for means in exact_leaves
            for eps in (-1, 1)
        ),
        "all_exact_allocations_covered": not missing,
        "all_join_sources_valid": not malformed,
    }
    proved = all(checks.values())
    return {
        "experiment": "p7_fixed_boundary_allocation_audit",
        "status": "complete_independent_exact_allocation_coverage_audit",
        "proved": proved,
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "direction_rows": rows,
        "exact_allocation_count": len(exact_leaves),
        "old_batch_allocation_count": len(batch_by_means),
        "old_batch_extra_nonadmissible_leaves": len(
            set(batch_by_means) - exact_leaf_set
        ),
        "method_counts": method_counts,
        "missing_allocations": missing,
        "malformed_sources": malformed,
        "checks": checks,
        "batch_source": {"path": str(batch_path), "sha256": sha256(batch_path)},
        "join_sources": join_sources,
        "coverage": coverage,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--joins", type=Path, nargs="*", default=())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = audit(args.batch, tuple(args.joins))
    if args.output is not None:
        atomic_write(args.output, out)
    compact = {key: value for key, value in out.items() if key != "coverage"}
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
