#!/usr/bin/env python3
"""Build a balanced exact 27-coordinate mod-5 exceptional projection."""
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

from p7_exceptional_projected_catalogs import atomic_json, atomic_npz, sha256, supported_unknown_means  # noqa: E402
from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog  # noqa: E402
from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402


# Five rows shared by directions 5 and 7, seven additional high-support rows
# for each, and eight high-support rows for direction 6.  Coverage is balanced
# at 12, 8, 12 over the three large variable catalogs in the hard leaf.
ROWS5 = (
    23, 28, 38, 40, 41, 42, 43,
    58, 63, 73, 74, 75, 76, 77, 78,
    93, 94, 96, 97, 98, 103, 108,
    109, 110, 111, 112, 113,
)


def pack(digits: np.ndarray) -> np.ndarray:
    if digits.ndim != 2 or digits.shape[0] > 14:
        raise ValueError("packed mod-5 half has at most fourteen rows")
    out = np.zeros(digits.shape[1], dtype=np.uint64)
    for index, row in enumerate(digits):
        out |= row.astype(np.uint64) << (3 * index)
    return out


def store27(arrays: dict[str, np.ndarray], stem: str, digits: np.ndarray) -> None:
    if digits.shape[0] != 27 or np.any(digits < 0) or np.any(digits >= 5):
        raise AssertionError("invalid 27-coordinate mod-5 signature")
    arrays[f"{stem}_lo"] = pack(digits[:13])
    arrays[f"{stem}_hi"] = pack(digits[13:])


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

    _matrix, dependencies, _linear_rows = linear_data((5,))
    dependency = dependencies[5]
    selected = np.asarray(ROWS5)
    arrays: dict[str, np.ndarray] = {}
    base = dependency[selected, :2] @ np.asarray([29, 1], dtype=np.int64) % 5
    store27(arrays, "base_p5", base.reshape(27, 1))
    cache = (Path(cache_path_text), Path(cache_summary_text))
    rows = direction_rows(-1, boundary)
    metadata = []
    for direction_index, row in enumerate(rows):
        means = sorted({leaf[direction_index] for leaf in supported})
        block = dependency[
            np.ix_(selected, np.arange(2 + 35*direction_index, 2 + 35*(direction_index + 1)))
        ]
        for mean in means:
            values = mapped_catalog(int(row["b"]), int(row["phase"]), int(mean), set(row["B"]), cache)
            digits = block @ ((13 - values.astype(np.int64)).T % 5) % 5
            stem = f"d{direction_index}_m{mean}_p5"
            store27(arrays, stem, digits)
            metadata.append({
                "direction_index": direction_index,
                "scaled_mean": int(mean),
                "catalog_rows": int(len(values)),
                "distinct_projected_signatures": int(len(set(zip(arrays[f"{stem}_lo"].tolist(), arrays[f"{stem}_hi"].tolist())))),
            })
    output = Path(output_dir_text) / f"cminus_exceptional_p5_27_orbit{orbit_index:02d}.npz"
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
    _matrix, dependencies, _rows = linear_data((5,))
    dependency = dependencies[5]
    selected = np.asarray(ROWS5)
    coverage = [
        int(np.sum(np.any(dependency[np.ix_(selected, np.arange(2 + 35*i, 2 + 35*(i+1)))], axis=1)))
        for i in range(8)
    ]
    if coverage != [5, 5, 5, 5, 5, 12, 8, 12]:
        raise AssertionError("balanced mod-5 coverage audit failed")
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
        "experiment": "p7_exceptional_p5_27_catalogs",
        "status": "complete_exact_selected_dependency_catalog_projection",
        "p": 7,
        "c_H": -1,
        "selected_dependency_rows_5": list(ROWS5),
        "selected_row_coverage_by_direction": coverage,
        "projection_group_order": 5**27,
        "projection_encoding": "two_packed_three_bit_words_to_injective_uint64_base5",
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
