#!/usr/bin/env python3
"""Measure the common anisotropic rank in sparse harmonic-shell scouts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--relative-spread-floor", type=float, default=1e-8)
    parser.add_argument("--scaled-parity", choices=("all", "even", "odd"), default="all")
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    labels = []
    centered_rows = []
    relative_spreads = []
    for scaled_norm, row in payload["aggregate_by_scaled_norm"].items():
        if args.scaled_parity != "all" and int(scaled_norm) % 2 != (args.scaled_parity == "odd"):
            continue
        values = np.asarray(row["phased_harmonic_sums"], dtype=np.float64)
        centered = values - values.mean()
        relative_spread = float(np.linalg.norm(centered) / max(1.0, abs(values.mean())))
        if relative_spread <= args.relative_spread_floor:
            continue
        labels.append(int(scaled_norm))
        centered_rows.append(centered)
        relative_spreads.append(relative_spread)
    matrix = np.asarray(centered_rows)
    normalized = matrix / np.linalg.norm(matrix, axis=1, keepdims=True)
    singular = np.linalg.svd(normalized, compute_uv=False)
    raw_singular = np.linalg.svd(matrix, compute_uv=False)
    rank_one_residual = float(
        np.linalg.norm(normalized - np.linalg.svd(normalized, full_matrices=False)[0][:, :1]
                       @ np.diag(singular[:1])
                       @ np.linalg.svd(normalized, full_matrices=False)[2][:1])
        / np.linalg.norm(normalized)
    )
    result = {
        "input": str(args.input),
        "scaled_parity": args.scaled_parity,
        "scaled_norms": labels,
        "relative_spreads": relative_spreads,
        "row_normalized_singular_values": singular.tolist(),
        "raw_singular_values": raw_singular.tolist(),
        "rank_one_relative_residual": rank_one_residual,
        "rank_one_to_second_singular_ratio": (
            float(singular[0] / singular[1]) if len(singular) > 1 and singular[1] else None
        ),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
