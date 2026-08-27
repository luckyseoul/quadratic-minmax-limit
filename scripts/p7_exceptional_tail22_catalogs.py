#!/usr/bin/env python3
"""Build a strong exact 22-coordinate mod-7 projection for exceptional joins."""
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
    atomic_json,
    atomic_npz,
    sha256,
    supported_unknown_means,
)
from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog  # noqa: E402
from p7_fixed_boundary_mean_allocation_batch import allocations, direction_data  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402


# These are 22 of the 23 mod-7 dependency rows that touch direction 7.  They
# also touch every other direction; row 133 is omitted because it has the
# smallest support in the sparse direction-7 block.  Since 7**22 < 2**62,
# the complete projection has an injective uint64 mixed-radix encoding.
ROWS7 = (
    109, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122,
    123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 134,
)

ROWS7_MID = (
    63, 78, 79, 80, 81, 82, 83, 84, 86, 87, 88,
    89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 133,
)

ROWS7_HEAD = (
    25, 26, 27, 28, 29,
    44, 45, 46, 47, 48, 49, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 64,
)


def pack11(digits: np.ndarray) -> np.ndarray:
    if digits.ndim != 2 or digits.shape[0] > 11:
        raise ValueError("pack11 expects at most eleven rows")
    packed = np.zeros(digits.shape[1], dtype=np.uint64)
    for index, row in enumerate(digits):
        packed |= row.astype(np.uint64) << (4 * index)
    return packed


def store22(arrays: dict[str, np.ndarray], stem: str, digits: np.ndarray) -> None:
    if digits.shape[0] != 22 or np.any(digits < 0) or np.any(digits >= 7):
        raise AssertionError("invalid 22-coordinate mod-7 signature")
    arrays[f"{stem}_lo"] = pack11(digits[:11])
    arrays[f"{stem}_hi"] = pack11(digits[11:])


def build_one(
    orbit_index: int,
    boundary: tuple[int, ...],
    mean_path_text: str,
    cache_path_text: str,
    cache_summary_text: str,
    output_dir_text: str,
    selected_rows: tuple[int, ...],
    profile: str,
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
    if len(allocations(direction_data(-1, boundary))) != 180:
        raise AssertionError("independent exceptional allocation count changed")
    supported = supported_unknown_means(batch)
    if not supported:
        raise AssertionError("no supported unknown leaves")

    _matrix, dependencies, _linear_rows = linear_data((7,))
    dependency = dependencies[7]
    selected = np.asarray(selected_rows)
    arrays: dict[str, np.ndarray] = {}
    base_digits = dependency[selected, :2] @ np.asarray([29, 1], dtype=np.int64) % 7
    store22(arrays, "base_p7", base_digits.reshape(22, 1))

    cache = (Path(cache_path_text), Path(cache_summary_text))
    rows = direction_rows(-1, boundary)
    metadata = []
    for direction_index, row in enumerate(rows):
        means = sorted({leaf[direction_index] for leaf in supported})
        block = dependency[
            np.ix_(selected, np.arange(2 + 35 * direction_index, 2 + 35 * (direction_index + 1)))
        ]
        for mean in means:
            values = mapped_catalog(
                int(row["b"]), int(row["phase"]), int(mean), set(row["B"]), cache
            )
            bad = 13 - values.astype(np.int64)
            digits = block @ (bad.T % 7) % 7
            store22(arrays, f"d{direction_index}_m{mean}_p7", digits)
            metadata.append(
                {
                    "direction_index": direction_index,
                    "scaled_mean": int(mean),
                    "catalog_rows": int(len(values)),
                    "distinct_projected_signatures": int(
                        len(set(zip(arrays[f"d{direction_index}_m{mean}_p7_lo"].tolist(),
                                    arrays[f"d{direction_index}_m{mean}_p7_hi"].tolist())))
                    ),
                }
            )

    output = Path(output_dir_text) / f"cminus_exceptional_{profile}22_orbit{orbit_index:02d}.npz"
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


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    source = json.loads(args.source.read_text())
    orbits = exceptional_orbits(source)
    if args.orbit_index is not None:
        orbits = [row for row in orbits if row[0] == args.orbit_index]
        if not orbits:
            raise ValueError("requested orbit is not exceptional")
    selected_rows = {
        "tail": ROWS7,
        "mid": ROWS7_MID,
        "head": ROWS7_HEAD,
    }[args.profile]
    _matrix, dependencies, _rows = linear_data((7,))
    dependency = dependencies[7]
    coverage = []
    for direction_index in range(8):
        block = dependency[np.ix_(np.asarray(selected_rows), np.arange(2 + 35 * direction_index, 2 + 35 * (direction_index + 1)))]
        coverage.append(int(np.sum(np.any(block, axis=1))))
    expected = {
        "tail": [22, 22, 22, 22, 22, 22, 21, 22],
        "mid": [22, 22, 22, 22, 22, 22, 21, 1],
        "head": [20, 20, 20, 20, 20, 17, 0, 0],
    }[args.profile]
    if coverage != expected:
        raise AssertionError("22-row dependency coverage audit failed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                build_one, orbit_index, boundary,
                str(args.mean_dir / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"),
                str(args.catalog_cache), str(args.catalog_cache_summary), str(args.output_dir),
                selected_rows, args.profile,
            )
            for orbit_index, boundary in orbits
        ]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps({key: value for key, value in row.items() if key != "catalogs"}), flush=True)
    rows.sort(key=lambda row: row["orbit_index"])
    return {
        "experiment": "p7_exceptional_tail22_catalogs",
        "status": "complete_exact_selected_dependency_catalog_projection",
        "p": 7,
        "c_H": -1,
        "projection_profile": args.profile,
        "selected_dependency_rows_7": list(selected_rows),
        "selected_row_coverage_by_direction": coverage,
        "projection_group_order": 7**22,
        "projection_encoding": "two_packed_nibble_words_to_injective_uint64_base7",
        "projection_is_necessary_not_sufficient": True,
        "orbit_count": len(rows),
        "total_supported_unknown_leaves": sum(row["supported_unknown_leaves"] for row in rows),
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
    parser.add_argument("--orbit-index", type=int)
    parser.add_argument("--profile", choices=("tail", "mid", "head"), default="tail")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    out = run(args)
    if args.summary is not None:
        atomic_json(args.summary, out)
    print(json.dumps({key: value for key, value in out.items() if key != "orbits"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
