#!/usr/bin/env python3
"""Exact catalog-pair check after the p=7 negative-infinity affine sieve."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from p7_unsaturated_mod7_batch import contribution_matrix  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


MODULUS = 7


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def run(orbit_source: Path, affine_source: Path) -> dict:
    started = time.time()
    orbit_payload = json.loads(orbit_source.read_text())
    affine = json.loads(affine_source.read_text())
    orbit_sha = sha256(orbit_source)
    if (
        int(orbit_payload["p"]) != 7
        or int(orbit_payload["c_H"]) != -1
        or int(orbit_payload["boundary_size"]) != 6
        or int(orbit_payload["infinity_value"]) != 1
        or int(orbit_payload["orbit_size_sum"])
        != int(orbit_payload["candidate_boundaries"])
    ):
        raise ValueError("orbit source has the wrong or incomplete scope")
    if (
        affine.get("experiment")
        != "p7_size6_negative_infinity_mod7_affine"
        or not affine.get("complete_orbit_enumeration")
        or affine.get("source_sha256") != orbit_sha
        or int(affine["checked_boundaries"]) != int(orbit_payload["orbit_count"])
        or int(affine["checked_boundary_weight"])
        != int(orbit_payload["candidate_boundaries"])
    ):
        raise ValueError("affine source does not certify the supplied orbit catalog")
    cases = affine["surviving_cases"]
    if len(cases) != int(affine["affine_span_surviving_cases"]):
        raise AssertionError("affine survivor case list is incomplete")

    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, MODULUS)
    if rank != 147 or dependencies.shape != (135, 282):
        raise AssertionError("unexpected mod-seven dependency dimensions")
    edge_base = (
        dependencies[:, :2] @ np.asarray([29, 1], dtype=np.int64) % MODULUS
    ).astype(np.int16)
    labels = []
    for direction in projective_directions(7):
        _eps, row = field_direction_data(7, direction)
        labels.append(tuple(int(value) for value in row))

    minimum_cache: dict[tuple[int, int], np.ndarray] = {}
    catalog_cache: dict[tuple[int, int], np.ndarray] = {}

    def fibre_set(mask: int) -> set[int]:
        return {value for value in range(7) if (mask >> value) & 1}

    def minimum(direction: int, mask: int) -> np.ndarray:
        key = (direction, mask)
        if key not in minimum_cache:
            b = mask.bit_count()
            if b not in (1, 5):
                raise AssertionError("mean-six direction does not have b=1 or b=5")
            values = contribution_matrix(
                dependencies, direction, b, 1, 6, fibre_set(mask)
            ).astype(np.int16)
            if values.shape[1] != 1:
                raise AssertionError("mean-six catalog is not unique")
            minimum_cache[key] = values[:, 0]
        return minimum_cache[key]

    def catalog(direction: int, mask: int) -> np.ndarray:
        key = (direction, mask)
        if key not in catalog_cache:
            b = mask.bit_count()
            values = contribution_matrix(
                dependencies, direction, b, 1, 14, fibre_set(mask)
            ).astype(np.int16)
            catalog_cache[key] = values
        return catalog_cache[key]

    total_pairs = 0
    surviving_pairs = []
    checked_cases = []
    for case_index, record in enumerate(cases):
        finite = tuple(int(value) for value in record["finite_boundary"])
        masks = []
        for direction in range(8):
            mask = 0
            for point in finite:
                mask ^= 1 << labels[direction][point]
            masks.append(mask)
        if masks != [int(value) for value in record["direction_masks"]]:
            raise AssertionError("recorded affine case has incorrect direction masks")
        directions = tuple(int(value) for value in record["elevated_directions"])
        if len(directions) != 2:
            raise AssertionError("affine survivor does not have two elevated directions")

        base = edge_base.copy()
        for direction, mask in enumerate(masks):
            if mask.bit_count() in (1, 5):
                base = (base + minimum(direction, mask)) % MODULUS
        catalogs = []
        for direction in directions:
            mask = masks[direction]
            if mask.bit_count() in (1, 5):
                base = (base - minimum(direction, mask)) % MODULUS
            catalogs.append(catalog(direction, mask))
        observed_sizes = [int(values.shape[1]) for values in catalogs]
        if observed_sizes != [int(value) for value in record["catalog_sizes"]]:
            raise AssertionError("exact catalog size differs from affine recording")

        case_survivors = []
        left, right = catalogs
        total_pairs += int(left.shape[1] * right.shape[1])
        right_lookup: dict[bytes, list[int]] = {}
        for right_index in range(right.shape[1]):
            key = np.asarray(right[:, right_index] % MODULUS, dtype=np.uint8).tobytes()
            right_lookup.setdefault(key, []).append(right_index)
        for left_index in range(left.shape[1]):
            target = np.asarray(
                (-base - left[:, left_index]) % MODULUS, dtype=np.uint8
            ).tobytes()
            for right_index in right_lookup.get(target, []):
                case_survivors.append([left_index, right_index])
                surviving_pairs.append(
                    {
                        "case_index": case_index,
                        "finite_boundary": list(finite),
                        "elevated_directions": list(directions),
                        "catalog_indices": [left_index, right_index],
                    }
                )
        checked_cases.append(
            {
                "case_index": case_index,
                "finite_boundary": list(finite),
                "orbit_size": int(record["orbit_size"]),
                "elevated_directions": list(directions),
                "catalog_sizes": observed_sizes,
                "checked_pairs": int(left.shape[1] * right.shape[1]),
                "surviving_pairs": len(case_survivors),
            }
        )

    return {
        "experiment": "p7_size6_negative_infinity_mod7_exact",
        "status": "complete_exact_catalog_pair_syndrome_check",
        "p": 7,
        "c_H": -1,
        "boundary_size": 6,
        "infinity_value": 1,
        "orbit_source": str(orbit_source),
        "orbit_source_sha256": orbit_sha,
        "affine_source": str(affine_source),
        "affine_source_sha256": sha256(affine_source),
        "source_orbits": int(orbit_payload["orbit_count"]),
        "source_boundary_weight": int(orbit_payload["candidate_boundaries"]),
        "affine_rejected_cases": int(affine["elevation_cases"])
        - int(affine["affine_span_surviving_cases"]),
        "exact_cases": len(cases),
        "checked_catalog_pairs": total_pairs,
        "checked_cases": checked_cases,
        "surviving_catalog_pairs": len(surviving_pairs),
        "surviving_pairs": surviving_pairs,
        "all_exact_cases_mod7_infeasible": not surviving_pairs,
        "complete_branch_mod7_infeasible": bool(
            not surviving_pairs
            and affine.get("complete_orbit_enumeration")
            and int(affine["floor_rejected_boundaries"]) == 0
        ),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbit-source", type=Path, required=True)
    parser.add_argument("--affine-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.orbit_source, args.affine_source)
    atomic_write(args.output, out)
    print(
        json.dumps(
            {
                "source_orbits": out["source_orbits"],
                "affine_rejected_cases": out["affine_rejected_cases"],
                "exact_cases": out["exact_cases"],
                "checked_catalog_pairs": out["checked_catalog_pairs"],
                "surviving_catalog_pairs": out["surviving_catalog_pairs"],
                "complete_branch_mod7_infeasible": out[
                    "complete_branch_mod7_infeasible"
                ],
                "elapsed_seconds": out["elapsed_seconds"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
