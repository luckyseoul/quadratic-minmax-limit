#!/usr/bin/env python3
"""Exact catalog joins for unresolved saturated p=7 conic mean leaves."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_fixed_boundary_catalog_join import run as run_join
from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data
from p7_size8_saturated_mean_batch import saturated_orbits


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def valid_batch(path: Path, boundary: tuple[int, ...]) -> dict:
    payload = json.loads(path.read_text())
    if (
        payload.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
        or payload.get("status") != "complete_exact_mean_allocation_exhaustion"
        or int(payload.get("p", 0)) != 7
        or int(payload.get("c_H", 0)) != -1
        or tuple(payload.get("fixed_boundary", [])) != boundary
        or int(payload.get("allocation_count", -1)) != 24
        or len(payload.get("leaves", [])) != 24
    ):
        raise ValueError(f"mean batch has incompatible metadata: {path}")
    return payload


def valid_join(path: Path, boundary: tuple[int, ...], means: tuple[int, ...]) -> dict | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if (
        payload.get("experiment") == "p7_fixed_boundary_catalog_join"
        and payload.get("status") == "complete_exact_multimodular_catalog_join"
        and int(payload.get("p", 0)) == 7
        and int(payload.get("c_H", 0)) == -1
        and tuple(payload.get("fixed_boundary", [])) == boundary
        and tuple(payload.get("fixed_scaled_means", [])) == means
    ):
        return payload
    raise ValueError(f"existing join has incompatible metadata: {path}")


def join_one(
    orbit_index: int,
    leaf_index: int,
    boundary: tuple[int, ...],
    means: tuple[int, ...],
    cache_path_text: str,
    cache_summary_text: str,
    output_dir_text: str,
) -> dict:
    output_dir = Path(output_dir_text)
    path = output_dir / f"cminus_saturated_orbit{orbit_index:02d}_leaf{leaf_index:02d}.json"
    existing = valid_join(path, boundary, means)
    if existing is not None:
        payload = existing
        reused = True
    else:
        payload = run_join(
            -1,
            boundary,
            means,
            (3, 7),
            (Path(cache_path_text), Path(cache_summary_text)),
        )
        atomic_write(path, payload)
        reused = False
    return {
        "orbit_index": orbit_index,
        "leaf_index": leaf_index,
        "fixed_scaled_means": list(means),
        "path": str(path),
        "reused": reused,
        "consistent_catalog_tuples": int(payload["consistent_catalog_tuples"]),
        "modularly_infeasible": bool(payload["modularly_infeasible"]),
        "elapsed_seconds": float(payload["elapsed_seconds"]),
    }


def run(
    source_path: Path,
    mean_dir: Path,
    output_dir: Path,
    cache_path: Path,
    cache_summary: Path,
    workers: int,
    extra_cp_paths: tuple[Path, ...] = (),
) -> dict:
    started = time.time()
    source = json.loads(source_path.read_text())
    orbits = saturated_orbits(source)
    output_dir.mkdir(parents=True, exist_ok=True)
    extra_cp = {}
    for path in extra_cp_paths:
        payload = json.loads(path.read_text())
        boundary = tuple(int(value) for value in payload.get("fixed_boundary", []))
        fixed = payload.get("fixed_scaled_means", {})
        means = tuple(int(fixed[str(index)]) for index in range(8))
        key = (boundary, means)
        if key in extra_cp:
            raise AssertionError("duplicate extra CP-SAT certificate")
        if (
            payload.get("experiment") != "p7_fixed_boundary_modular_cpsat"
            or payload.get("solver_status") != "INFEASIBLE"
            or payload.get("finite_infeasibility_certificate") is not True
            or payload.get("feasible") is not False
            or int(payload.get("p", 0)) != 7
            or int(payload.get("c_H", 0)) != -1
        ):
            raise ValueError(f"extra CP-SAT file is not an exact exclusion: {path}")
        extra_cp[key] = (payload, path)
    tasks = []
    cp_infeasible = 0
    extra_cp_infeasible = 0
    used_extra_cp = set()
    deferred = []
    for orbit_index, boundary in orbits:
        batch_path = mean_dir / f"cminus_saturated_orbit{orbit_index:02d}_means.json"
        batch = valid_batch(batch_path, boundary)
        exact_leaves = allocations(direction_data(-1, boundary))
        if len(exact_leaves) != 24:
            raise AssertionError("independent allocation count changed")
        batch_by_means = {
            tuple(row["scaled_means_direction_order"]): row for row in batch["leaves"]
        }
        if len(batch_by_means) != 24 or set(batch_by_means) != set(exact_leaves):
            raise AssertionError("mean batch does not equal independent allocation set")
        for leaf_index, means in enumerate(exact_leaves):
            row = batch_by_means[means]
            if (
                row.get("solver_status") == "INFEASIBLE"
                and row.get("finite_infeasibility_certificate") is True
                and row.get("feasible") is False
            ):
                cp_infeasible += 1
            else:
                extra_key = (boundary, means)
                if extra_key in extra_cp:
                    extra_cp_infeasible += 1
                    used_extra_cp.add(extra_key)
                    continue
                high_catalogs = []
                for direction_index, (direction, mean) in enumerate(
                    zip(batch["direction_rows"], means)
                ):
                    b = int(direction["b"])
                    phase = int(direction["phase"])
                    minimum = 8 if b in (1, 2, 5, 6) and phase == 0 else (
                        6 if b in (1, 2, 5, 6) else 0
                    )
                    excess_mean = int(mean) - minimum
                    if phase == 0 and excess_mean > 16:
                        high_catalogs.append(
                            {
                                "direction_index": direction_index,
                                "b": b,
                                "phase": phase,
                                "scaled_mean": int(mean),
                                "zero_parity_excess_mean": excess_mean,
                            }
                        )
                if high_catalogs:
                    deferred.append(
                        {
                            "orbit_index": orbit_index,
                            "leaf_index": leaf_index,
                            "fixed_boundary": list(boundary),
                            "fixed_scaled_means": list(means),
                            "reason": "requires an uncached zero-parity catalog above mean 16",
                            "high_catalogs": high_catalogs,
                        }
                    )
                else:
                    tasks.append((orbit_index, leaf_index, boundary, means))

    if used_extra_cp != set(extra_cp):
        raise AssertionError("an extra CP-SAT certificate did not match an exact open leaf")

    rows = []
    if workers == 1:
        for orbit_index, leaf_index, boundary, means in tasks:
            row = join_one(
                orbit_index,
                leaf_index,
                boundary,
                means,
                str(cache_path),
                str(cache_summary),
                str(output_dir),
            )
            rows.append(row)
            print(json.dumps(row), flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    join_one,
                    orbit_index,
                    leaf_index,
                    boundary,
                    means,
                    str(cache_path),
                    str(cache_summary),
                    str(output_dir),
                ): (orbit_index, leaf_index)
                for orbit_index, leaf_index, boundary, means in tasks
            }
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                print(json.dumps(row), flush=True)
    rows.sort(key=lambda row: (row["orbit_index"], row["leaf_index"]))
    unresolved = [row for row in rows if not row["modularly_infeasible"]]
    return {
        "experiment": "p7_size8_saturated_join_batch",
        "status": "complete_exact_unknown_leaf_catalog_joins",
        "p": 7,
        "c_H": -1,
        "source": str(source_path),
        "saturated_orbit_count": len(orbits),
        "total_exact_allocations": 24 * len(orbits),
        "cp_sat_infeasible_allocations": cp_infeasible,
        "extra_cp_sat_infeasible_allocations": extra_cp_infeasible,
        "catalog_join_allocations": len(rows),
        "catalog_join_infeasible_allocations": sum(
            row["modularly_infeasible"] for row in rows
        ),
        "remaining_modularly_consistent_allocations": len(unresolved),
        "deferred_high_catalog_allocations": len(deferred),
        "all_saturated_orbits_excluded": not unresolved
        and not deferred
        and cp_infeasible + extra_cp_infeasible + len(rows)
        == 24 * len(orbits),
        "workers": workers,
        "elapsed_seconds": time.time() - started,
        "unresolved": unresolved,
        "deferred": deferred,
        "extra_cp_sources": [str(path) for path in extra_cp_paths],
        "joins": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mean-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--catalog-cache", type=Path, required=True)
    parser.add_argument("--catalog-cache-summary", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--extra-cp-results", type=Path, nargs="*", default=())
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(
        args.source,
        args.mean_dir,
        args.output_dir,
        args.catalog_cache,
        args.catalog_cache_summary,
        args.workers,
        tuple(args.extra_cp_results),
    )
    if args.summary is not None:
        atomic_write(args.summary, out)
    compact = {
        key: value
        for key, value in out.items()
        if key not in ("joins", "unresolved", "deferred")
    }
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
