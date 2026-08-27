#!/usr/bin/env python3
"""Build a direction-aware exact 40-coordinate mod-3 exceptional projection."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_exceptional_projected_catalogs import (  # noqa: E402
    atomic_json, atomic_npz, sha256, supported_unknown_means,
)
from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog  # noqa: E402
from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402


# All 16 rows touching direction 6, all 16 touching direction 7, and the eight
# highest-support remaining rows for direction 5.  The resulting coverage by
# direction is [5,5,20,12,5,13,16,16].  The whole group is injectively encoded
# because 3**40 < 2**64.
ROWS3 = (
    30, 33, 34, 40, 44, 47, 48, 49,
    59, 64, 65, 67, 68, 69, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84,
    94, 99, 100, 102, 103, 104, 109, 110, 112, 113, 114, 115, 116, 117, 118, 119,
)


def pack20(digits: np.ndarray) -> np.ndarray:
    if digits.ndim != 2 or digits.shape[0] > 20:
        raise ValueError("pack20 expects at most twenty rows")
    packed = np.zeros(digits.shape[1], dtype=np.uint64)
    for index, row in enumerate(digits):
        packed |= row.astype(np.uint64) << (2 * index)
    return packed


def store40(arrays: dict[str, np.ndarray], stem: str, digits: np.ndarray) -> None:
    if digits.shape[0] != 40 or np.any(digits < 0) or np.any(digits >= 3):
        raise AssertionError("invalid 40-coordinate mod-3 signature")
    arrays[f"{stem}_lo"] = pack20(digits[:20])
    arrays[f"{stem}_hi"] = pack20(digits[20:])


def build_one(orbit_index, boundary, mean_path_text, cache_path_text, cache_summary_text, output_dir_text):
    started = time.time()
    mean_path = Path(mean_path_text)
    batch = json.loads(mean_path.read_text())
    if (
        batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
        or tuple(batch.get("fixed_boundary", [])) != boundary
        or int(batch.get("allocation_count", -1)) != 180
    ):
        raise ValueError("exceptional mean batch has incompatible metadata")
    if len(allocations(direction_data(-1, boundary))) != 180:
        raise AssertionError("independent exceptional allocation count changed")
    supported = supported_unknown_means(batch)
    if not supported:
        raise AssertionError("no supported unknown leaves")

    _matrix, dependencies, _linear_rows = linear_data((3,))
    dependency = dependencies[3]
    selected = np.asarray(ROWS3)
    arrays: dict[str, np.ndarray] = {}
    base_digits = dependency[selected, :2] @ np.asarray([29, 1], dtype=np.int64) % 3
    store40(arrays, "base_p3", base_digits.reshape(40, 1))
    cache = (Path(cache_path_text), Path(cache_summary_text))
    rows = direction_rows(-1, boundary)
    metadata = []
    for direction_index, row in enumerate(rows):
        means = sorted({leaf[direction_index] for leaf in supported})
        block = dependency[
            np.ix_(selected, np.arange(2 + 35 * direction_index, 2 + 35 * (direction_index + 1)))
        ]
        for mean in means:
            values = mapped_catalog(int(row["b"]), int(row["phase"]), int(mean), set(row["B"]), cache)
            digits = block @ ((13 - values.astype(np.int64)).T % 3) % 3
            stem = f"d{direction_index}_m{mean}_p3"
            store40(arrays, stem, digits)
            metadata.append({
                "direction_index": direction_index,
                "scaled_mean": int(mean),
                "catalog_rows": int(len(values)),
                "distinct_projected_signatures": int(len(set(zip(arrays[f"{stem}_lo"].tolist(), arrays[f"{stem}_hi"].tolist())))),
            })
    output = Path(output_dir_text) / f"cminus_exceptional_p3_40_orbit{orbit_index:02d}.npz"
    atomic_npz(output, arrays)
    return {
        "orbit_index": orbit_index,
        "fixed_boundary": list(boundary),
        "mean_source": str(mean_path),
        "supported_unknown_leaves": len(supported),
        "output": str(output),
        "sha256": sha256(output),
        "bytes": output.stat().st_size,
        "catalogs": metadata,
        "elapsed_seconds": time.time() - started,
    }


def run(args):
    started = time.time()
    orbits = exceptional_orbits(json.loads(args.source.read_text()))
    if args.orbit_index is not None:
        orbits = [row for row in orbits if row[0] == args.orbit_index]
        if not orbits:
            raise ValueError("requested orbit is not exceptional")
    _matrix, dependencies, _rows = linear_data((3,))
    dependency = dependencies[3]
    selected = np.asarray(ROWS3)
    coverage = [
        int(np.sum(np.any(dependency[np.ix_(selected, np.arange(2 + 35*i, 2 + 35*(i+1)))], axis=1)))
        for i in range(8)
    ]
    if coverage != [5, 5, 20, 12, 5, 13, 16, 16]:
        raise AssertionError("directed mod-3 coverage audit failed")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(
            build_one, orbit_index, boundary,
            str(args.mean_dir / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"),
            str(args.catalog_cache), str(args.catalog_cache_summary), str(args.output_dir),
        ) for orbit_index, boundary in orbits]
        for future in as_completed(futures):
            row = future.result(); rows.append(row)
            print(json.dumps({key: value for key, value in row.items() if key != "catalogs"}), flush=True)
    rows.sort(key=lambda row: row["orbit_index"])
    return {
        "experiment": "p7_exceptional_p3_40_catalogs",
        "status": "complete_exact_selected_dependency_catalog_projection",
        "p": 7,
        "c_H": -1,
        "selected_dependency_rows_3": list(ROWS3),
        "selected_row_coverage_by_direction": coverage,
        "projection_group_order": 3**40,
        "projection_encoding": "two_packed_two_bit_words_to_injective_uint64_base3",
        "projection_is_necessary_not_sufficient": True,
        "orbit_count": len(rows),
        "total_supported_unknown_leaves": sum(row["supported_unknown_leaves"] for row in rows),
        "elapsed_seconds": time.time() - started,
        "orbits": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mean-dir", type=Path, required=True)
    parser.add_argument("--catalog-cache", type=Path, required=True)
    parser.add_argument("--catalog-cache-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(args)
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(json.dumps({key: value for key, value in out.items() if key != "orbits"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
