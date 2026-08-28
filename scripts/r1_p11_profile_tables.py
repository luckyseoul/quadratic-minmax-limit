#!/usr/bin/env python3
"""Exact one-profile tables for the p=11 glue-dual theta counter.

For a value-distribution ``n_y`` and an integral profile ``a_s``, tabulate

    sum_s a_s = b,
    sum_s a_s^2 = b + 2 k,
    sum_s y_s a_s = h (mod 11).

Only thirteen affine output types occur among the 604 quartic histograms.
The dynamic program is run on those thirteen canonical types; all remaining
tables are reconstructed under ``y -> alpha*y + beta``.  Shifting every
profile entry by an integer shows that residue sums ``b=0,...,10`` suffice
for every common sum used by the lattice theta series.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from r1_p11_profile_dual_orbits import P


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def affine_transform(histogram: tuple[int, ...], alpha: int, beta: int) -> tuple[int, ...]:
    transformed = [0] * P
    for value, count in enumerate(histogram):
        transformed[(alpha * value + beta) % P] = count
    return tuple(transformed)


def canonicalize_histograms(
    histograms: np.ndarray,
) -> tuple[list[tuple[int, ...]], list[tuple[int, int, int]]]:
    canonical_types: list[tuple[int, ...]] = []
    canonical_ids: dict[tuple[int, ...], int] = {}
    descriptors: list[tuple[int, int, int]] = []
    for row in histograms:
        histogram = tuple(int(value) for value in row)
        candidates = [
            (affine_transform(histogram, alpha, beta), alpha, beta)
            for alpha in range(1, P)
            for beta in range(P)
        ]
        canonical, alpha, beta = min(candidates)
        identifier = canonical_ids.get(canonical)
        if identifier is None:
            identifier = len(canonical_types)
            canonical_ids[canonical] = identifier
            canonical_types.append(canonical)
        descriptors.append((identifier, alpha, beta))
    return canonical_types, descriptors


def distribution_table(histogram: tuple[int, ...], max_k: int) -> np.ndarray:
    """Return counts with shape ``(11, max_k+1, 11)`` using uint64 exactly."""
    max_energy = (P - 1) + 2 * max_k
    sum_radius = math.isqrt(P * max_energy) + 1
    width = 2 * sum_radius + 1
    values = [value for value, count in enumerate(histogram) for _ in range(count)]
    if len(values) != P:
        raise ArithmeticError("histogram does not contain eleven positions")

    current = np.zeros((width, max_energy + 1, P), dtype=np.uint64)
    current[sum_radius, 0, 0] = 1
    coordinate_bound = math.isqrt(max_energy)
    for phase_value in values:
        following = np.zeros_like(current)
        for coordinate in range(-coordinate_bound, coordinate_bound + 1):
            square = coordinate * coordinate
            if square > max_energy:
                continue
            if coordinate >= 0:
                source_sum = slice(0, width - coordinate)
                target_sum = slice(coordinate, width)
            else:
                source_sum = slice(-coordinate, width)
                target_sum = slice(0, width + coordinate)
            source = current[source_sum, : max_energy + 1 - square, :]
            shifted = np.roll(source, phase_value * coordinate, axis=2)
            following[target_sum, square:, :] += shifted
        current = following

    output = np.zeros((P, max_k + 1, P), dtype=np.uint64)
    for residue_sum in range(P):
        for excess in range(max_k + 1):
            energy = residue_sum + 2 * excess
            output[residue_sum, excess] = current[residue_sum + sum_radius, energy]
    return output


def reconstruct_tables(
    canonical_tables: np.ndarray,
    descriptors: list[tuple[int, int, int]],
) -> np.ndarray:
    table_count, _, max_k_plus_one, _ = canonical_tables.shape
    if table_count == 0:
        raise ArithmeticError("empty canonical table family")
    output = np.zeros((len(descriptors), P, max_k_plus_one, P), dtype=np.uint64)
    for histogram_id, (canonical_id, alpha, beta) in enumerate(descriptors):
        alpha_inverse = pow(alpha, P - 2, P)
        for residue_sum in range(P):
            for original_phase in range(P):
                canonical_phase = (alpha * original_phase + beta * residue_sum) % P
                output[histogram_id, residue_sum, :, original_phase] = canonical_tables[
                    canonical_id, residue_sum, :, canonical_phase
                ]
        # The inverse is recorded implicitly by the exact phase permutation;
        # retaining this calculation catches a malformed nonunit alpha early.
        if alpha * alpha_inverse % P != 1:
            raise ArithmeticError("affine multiplier is not invertible")
    return output


def direct_audit(
    histograms: np.ndarray,
    tables: np.ndarray,
    max_k: int,
    sample_count: int,
) -> list[int]:
    rng = np.random.default_rng(0x11A771)
    fixed = {0, 1, len(histograms) - 1}
    fixed.update(int(value) for value in rng.choice(len(histograms), sample_count, replace=False))
    audited = sorted(fixed)
    for identifier in audited:
        direct = distribution_table(tuple(int(value) for value in histograms[identifier]), max_k)
        if not np.array_equal(direct, tables[identifier]):
            mismatch = np.argwhere(direct != tables[identifier])[0].tolist()
            raise ArithmeticError(
                f"affine table reconstruction failed for histogram {identifier} at {mismatch}"
            )
    return audited


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tuples", type=Path, required=True)
    parser.add_argument("--max-k", type=int, default=21)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--direct-audit-samples", type=int, default=13)
    args = parser.parse_args()
    if args.max_k < 0:
        raise ValueError("max-k must be nonnegative")

    started = time.monotonic()
    with np.load(args.tuples) as archive:
        histograms = np.asarray(archive["histograms"], dtype=np.uint8)
    if histograms.shape != (604, P):
        raise ArithmeticError(f"unexpected histogram table shape {histograms.shape}")
    canonical_types, descriptors = canonicalize_histograms(histograms)
    if len(canonical_types) != 13:
        raise ArithmeticError(f"affine histogram orbit count changed to {len(canonical_types)}")
    canonical_tables = np.stack(
        [distribution_table(histogram, args.max_k) for histogram in canonical_types]
    )
    tables = reconstruct_tables(canonical_tables, descriptors)
    audited = direct_audit(
        histograms, tables, args.max_k, args.direct_audit_samples
    )

    unphased = tables.sum(axis=3, dtype=np.uint64)
    if not np.all(unphased == unphased[0:1]):
        mismatch = np.argwhere(unphased != unphased[0:1])[0].tolist()
        raise ArithmeticError(f"unphased profile count depends on histogram at {mismatch}")
    if np.any(tables[0, :, :, 1:]):
        raise ArithmeticError("zero phase histogram has a nonzero phase count")
    if not np.array_equal(tables[0, :, :, 0], unphased[0]):
        raise ArithmeticError("zero phase histogram does not recover unphased counts")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        counts=tables,
        canonical_histograms=np.asarray(canonical_types, dtype=np.uint8),
        affine_descriptors=np.asarray(descriptors, dtype=np.uint8),
    )
    elapsed = time.monotonic() - started
    report = {
        "experiment": "r1_p11_profile_tables",
        "status": "complete_exact_single_profile_tables",
        "p": P,
        "quartic_histogram_types": int(len(histograms)),
        "affine_output_types": len(canonical_types),
        "max_excess_parameter_k": args.max_k,
        "maximum_profile_energy": P - 1 + 2 * args.max_k,
        "table_shape": list(tables.shape),
        "maximum_table_entry": int(tables.max()),
        "maximum_unphased_entry": int(unphased.max()),
        "direct_histogram_tables_audited": audited,
        "all_unphased_counts_histogram_independent": True,
        "zero_phase_table_audited": True,
        "elapsed_seconds": elapsed,
        "output": str(args.output),
        "output_sha256": sha256(args.output),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
