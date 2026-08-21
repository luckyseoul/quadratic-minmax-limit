#!/usr/bin/env python3
"""Extract depressed p=13,k=7 representatives from a packed PSL orbit.

This recovers the polynomial profile coefficients independently from the
CP-SAT variables.  The small scalar-7 output can be used to ask CP-SAT only
whether a solution exists outside the known free orbit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "evidence"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15588 import _fit_poly_modp  # noqa: E402
from k5_p29_coefficient_sieve import (  # noqa: E402
    homogeneous_matrix,
    kernel_modp,
    square_directions,
)


def unpack_finite_bits(packed: np.ndarray, q: int) -> np.ndarray:
    indices = np.arange(1, q + 1, dtype=np.int64)
    return (
        (packed[:, indices // 64] >> (indices % 64).astype(np.uint64))
        & np.uint64(1)
    ).astype(np.uint8)


def interpolation_inverse(p: int) -> np.ndarray:
    """Return A^-1 for A[s,e]=s^e over F_p."""
    inverse = np.empty((p, p), dtype=np.int64)
    for value in range(p):
        basis = np.zeros(p, dtype=np.int64)
        basis[value] = 1
        inverse[:, value] = _fit_poly_modp(basis, p)
    return inverse


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("orbit", type=Path)
    parser.add_argument("representatives", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--batch-size", type=int, default=32_768)
    args = parser.parse_args()

    p = 13
    q = p * p
    directions = square_directions(p)
    coordinates = [coordinate for coordinate, _form in directions]
    forms = [form for _coordinate, form in directions]
    line_incidence = np.zeros((q, 7 * p), dtype=np.uint8)
    for direction, coordinate in enumerate(coordinates):
        for value in range(p):
            line_incidence[coordinate == value, direction * p + value] = 1
    inverse = interpolation_inverse(p)
    top_kernel = kernel_modp(homogeneous_matrix(forms, 5, p), p)[0]

    packed = np.load(args.orbit, mmap_mode="r")
    scalar_histogram: Counter[int] = Counter()
    depressed_count = 0
    selected = []
    for lo in range(0, len(packed), args.batch_size):
        hi = min(lo + args.batch_size, len(packed))
        batch = np.asarray(packed[lo:hi])
        batch = batch[(batch[:, 0] & 1) == 0]
        if len(batch) == 0:
            continue
        bits = unpack_finite_bits(batch, q)
        negative_counts = (bits @ line_incidence).reshape(-1, 7, p)
        active = np.any(negative_counts != 6, axis=2)
        keep = np.sum(active, axis=1) == 7
        if not np.any(keep):
            continue
        batch = batch[keep]
        negative_counts = negative_counts[keep].astype(np.int64)
        line_sums = p - 2 * negative_counts
        rho = ((line_sums + p - 2) // 2) % p
        coefficients = (rho @ inverse.T) % p
        depressed = np.all(coefficients[:, :, 4] == 0, axis=1)
        batch = batch[depressed]
        leading = coefficients[depressed, :, 5]
        depressed_count += len(batch)
        for scalar in range(1, p):
            scalar_mask = np.all(
                leading == (scalar * top_kernel % p)[None, :], axis=1
            )
            count = int(np.count_nonzero(scalar_mask))
            scalar_histogram[scalar] += count
            if scalar == 7 and count:
                selected.append(batch[scalar_mask])

    representatives = np.concatenate(selected, axis=0)
    representatives = np.unique(representatives, axis=0)
    args.representatives.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.representatives, representatives)
    digest = hashlib.sha256(args.representatives.read_bytes()).hexdigest()
    report = {
        "p": p,
        "source_orbit": str(args.orbit),
        "depressed_k7_representatives": depressed_count,
        "top_scalar_histogram": {
            str(scalar): scalar_histogram[scalar] for scalar in range(1, p)
        },
        "scalar7_representatives": len(representatives),
        "representatives_path": str(args.representatives),
        "representatives_sha256": digest,
        "independent_profile_reconstruction": True,
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    main()
