#!/usr/bin/env python3
"""Build exact low-dimensional modular signatures for exceptional p=7 joins."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog  # noqa: E402
from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def atomic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    with temporary.open("wb") as handle:
        np.savez(handle, **arrays)
    os.replace(temporary, path)


def pack_digits(digits: np.ndarray, bits: int, dtype) -> np.ndarray:
    if digits.ndim != 2:
        raise ValueError("projected digits must be a matrix")
    packed = np.zeros(digits.shape[1], dtype=dtype)
    for index, row in enumerate(digits):
        packed |= row.astype(dtype) << (bits * index)
    return packed


def supported_unknown_means(batch: dict) -> tuple[tuple[int, ...], ...]:
    rows = batch["direction_rows"]
    out = []
    for leaf in batch["leaves"]:
        if leaf.get("solver_status") == "INFEASIBLE":
            continue
        means = tuple(int(value) for value in leaf["scaled_means_direction_order"])
        if all(
            not (
                int(row["b"]) == 0
                and int(row["phase"]) == 0
                and mean > 16
            )
            for row, mean in zip(rows, means)
        ):
            out.append(means)
    return tuple(out)


def build_one(
    orbit_index: int,
    boundary: tuple[int, ...],
    mean_path_text: str,
    cache_path_text: str,
    cache_summary_text: str,
    output_dir_text: str,
    rows3: tuple[int, ...],
    rows7: tuple[int, ...],
) -> dict:
    started = time.time()
    mean_path = Path(mean_path_text)
    batch = json.loads(mean_path.read_text())
    if (
        batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
        or tuple(batch.get("fixed_boundary", [])) != boundary
        or int(batch.get("allocation_count", -1)) != 180
    ):
        raise ValueError("exceptional mean batch has incompatible metadata")
    exact = allocations(direction_data(-1, boundary))
    if len(exact) != 180:
        raise AssertionError("independent exceptional allocation count changed")
    supported = supported_unknown_means(batch)
    if not supported:
        raise AssertionError("no supported unknown leaves")

    _matrix, dependencies, _linear_rows = linear_data((3, 7))
    selected = {3: np.asarray(rows3), 7: np.asarray(rows7)}
    arrays: dict[str, np.ndarray] = {}
    global_rhs = np.asarray([29, 1], dtype=np.int64)
    for modulus, bits, dtype in ((3, 2, np.uint32), (7, 4, np.uint64)):
        dependency = dependencies[modulus]
        digits = dependency[selected[modulus], :2] @ global_rhs % modulus
        arrays[f"base_p{modulus}"] = pack_digits(
            digits.reshape(-1, 1), bits, dtype
        )

    cache = (Path(cache_path_text), Path(cache_summary_text))
    rows = direction_rows(-1, boundary)
    catalog_metadata = []
    for direction_index, row in enumerate(rows):
        means = sorted({leaf[direction_index] for leaf in supported})
        for mean in means:
            values = mapped_catalog(
                int(row["b"]),
                int(row["phase"]),
                int(mean),
                set(row["B"]),
                cache,
            )
            bad = 13 - values.astype(np.int64)
            for modulus, bits, dtype in ((3, 2, np.uint32), (7, 4, np.uint64)):
                dependency = dependencies[modulus]
                block = dependency[
                    np.ix_(
                        selected[modulus],
                        np.arange(
                            2 + 35 * direction_index,
                            2 + 35 * (direction_index + 1),
                        ),
                    )
                ]
                digits = block @ (bad.T % modulus) % modulus
                arrays[f"d{direction_index}_m{mean}_p{modulus}"] = pack_digits(
                    digits, bits, dtype
                )
            catalog_metadata.append(
                {
                    "direction_index": direction_index,
                    "eps": int(row["eps"]),
                    "b": int(row["b"]),
                    "phase": int(row["phase"]),
                    "scaled_mean": int(mean),
                    "catalog_rows": int(len(values)),
                }
            )

    output = Path(output_dir_text) / f"cminus_exceptional_orbit{orbit_index:02d}.npz"
    atomic_npz(output, arrays)
    return {
        "orbit_index": orbit_index,
        "fixed_boundary": list(boundary),
        "mean_source": str(mean_path),
        "supported_unknown_leaves": len(supported),
        "output": str(output),
        "sha256": sha256(output),
        "bytes": output.stat().st_size,
        "catalogs": catalog_metadata,
        "elapsed_seconds": time.time() - started,
    }


def run(
    source_path: Path,
    mean_dir: Path,
    cache_path: Path,
    cache_summary: Path,
    output_dir: Path,
    workers: int,
    seed: int,
) -> dict:
    started = time.time()
    source = json.loads(source_path.read_text())
    orbits = exceptional_orbits(source)
    _matrix, dependencies, _linear_rows = linear_data((3, 7))
    rng = np.random.default_rng(seed)
    rows3 = tuple(
        sorted(rng.choice(dependencies[3].shape[0], size=12, replace=False).tolist())
    )
    rows7 = tuple(
        sorted(rng.choice(dependencies[7].shape[0], size=13, replace=False).tolist())
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                build_one,
                orbit_index,
                boundary,
                str(mean_dir / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"),
                str(cache_path),
                str(cache_summary),
                str(output_dir),
                rows3,
                rows7,
            ): orbit_index
            for orbit_index, boundary in orbits
        }
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps({key: value for key, value in row.items() if key != "catalogs"}), flush=True)
    rows.sort(key=lambda row: row["orbit_index"])
    return {
        "experiment": "p7_exceptional_projected_catalogs",
        "status": "complete_exact_selected_dependency_catalog_projection",
        "p": 7,
        "c_H": -1,
        "source": str(source_path),
        "projection_seed": seed,
        "selected_dependency_rows": {"3": list(rows3), "7": list(rows7)},
        "projection_group_order": 3**12 * 7**13,
        "projection_is_necessary_not_sufficient": True,
        "orbit_count": len(rows),
        "total_supported_unknown_leaves": sum(
            row["supported_unknown_leaves"] for row in rows
        ),
        "workers": workers,
        "elapsed_seconds": time.time() - started,
        "orbits": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mean-dir", type=Path, required=True)
    parser.add_argument("--catalog-cache", type=Path, required=True)
    parser.add_argument("--catalog-cache-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--seed", type=int, default=15708301)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(
        args.source,
        args.mean_dir,
        args.catalog_cache,
        args.catalog_cache_summary,
        args.output_dir,
        args.workers,
        args.seed,
    )
    if args.summary is not None:
        atomic_json(args.summary, out)
    compact = {key: value for key, value in out.items() if key != "orbits"}
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
