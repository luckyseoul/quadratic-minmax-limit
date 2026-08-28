#!/usr/bin/env python3
"""Exact one-profile count and square-circle U4 tables at p=11.

For a quartic residue sequence ``y_s`` and an integral profile ``a_s`` the
ordinary table records

    sum a_s=b,  sum a_s^2=b+2k,  sum y_s a_s=h (mod 11).

This refinement also records the nonnegative integer weight

    U4(a)=sum_c (sum_s eta(s-c) a_s)^4.

Input-affine permutations preserve U4.  The 1,007 rich profile types are
therefore reconstructed from only 20 domain/output-affine canonical dynamic
programs, with the output-affine phase transported exactly.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from r1_p11_channel_profile_types import rich_profile_data
from r1_p11_profile_dual_orbits import P


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def legendre(value: int) -> int:
    value %= P
    if value == 0:
        return 0
    return 1 if pow(value, (P - 1) // 2, P) == 1 else -1


LEGENDRE_SHIFTS = np.asarray(
    [[legendre(position - shift) for position in range(P)] for shift in range(P)],
    dtype=np.int64,
)


def distribution_tables(
    sequence: tuple[int, ...], max_k: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact count and U4 tables, each shaped ``(11,max_k+1,11)``."""
    max_energy = (P - 1) + 2 * max_k
    sum_radius = math.isqrt(P * max_energy) + 1
    width = 2 * sum_radius + 1
    coordinate_bound = math.isqrt(max_energy)

    count = np.zeros((width, max_energy + 1, P), dtype=np.uint64)
    count[sum_radius, 0, 0] = 1
    # Axis order: Legendre shift, moment order 1..4, sum, energy, phase.
    moments = np.zeros((P, 4, width, max_energy + 1, P), dtype=np.int64)
    delta_shape = (P, 1, 1, 1)

    for position, phase_value in enumerate(sequence):
        following_count = np.zeros_like(count)
        following_moments = np.zeros_like(moments)
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
            source_count_u = count[source_sum, : max_energy + 1 - square, :]
            source_count = source_count_u.astype(np.int64, copy=False)
            phase = int(phase_value) * coordinate
            following_count[target_sum, square:, :] += np.roll(
                source_count_u, phase, axis=-1
            )

            source = moments[:, :, source_sum, : max_energy + 1 - square, :]
            first = source[:, 0]
            second = source[:, 1]
            third = source[:, 2]
            fourth = source[:, 3]
            delta = (LEGENDRE_SHIFTS[:, position] * coordinate).reshape(delta_shape)
            delta2 = delta * delta
            delta3 = delta2 * delta
            delta4 = delta2 * delta2
            broadcast_count = source_count[None, ...]
            transformed = np.stack(
                (
                    first + delta * broadcast_count,
                    second + 2 * delta * first + delta2 * broadcast_count,
                    third
                    + 3 * delta * second
                    + 3 * delta2 * first
                    + delta3 * broadcast_count,
                    fourth
                    + 4 * delta * third
                    + 6 * delta2 * second
                    + 4 * delta3 * first
                    + delta4 * broadcast_count,
                ),
                axis=1,
            )
            following_moments[:, :, target_sum, square:, :] += np.roll(
                transformed, phase, axis=-1
            )
        count = following_count
        moments = following_moments

    output_count = np.zeros((P, max_k + 1, P), dtype=np.uint64)
    output_u4 = np.zeros_like(output_count)
    summed_fourth = moments[:, 3].sum(axis=0)
    if np.any(summed_fourth < 0):
        raise ArithmeticError("U4 dynamic program produced a negative bucket")
    for residue_sum in range(P):
        for excess in range(max_k + 1):
            energy = residue_sum + 2 * excess
            output_count[residue_sum, excess] = count[
                residue_sum + sum_radius, energy
            ]
            output_u4[residue_sum, excess] = summed_fourth[
                residue_sum + sum_radius, energy
            ].astype(np.uint64)
    return output_count, output_u4


def reconstruct_tables(
    canonical_count: np.ndarray,
    canonical_u4: np.ndarray,
    descriptors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    max_k_plus_one = canonical_count.shape[2]
    shape = (len(descriptors), P, max_k_plus_one, P)
    counts = np.zeros(shape, dtype=np.uint64)
    u4 = np.zeros(shape, dtype=np.uint64)
    for rich_id, descriptor in enumerate(descriptors):
        table_id, alpha, beta = (int(value) for value in descriptor)
        for residue_sum in range(P):
            for original_phase in range(P):
                canonical_phase = (
                    alpha * original_phase + beta * residue_sum
                ) % P
                counts[rich_id, residue_sum, :, original_phase] = canonical_count[
                    table_id, residue_sum, :, canonical_phase
                ]
                u4[rich_id, residue_sum, :, original_phase] = canonical_u4[
                    table_id, residue_sum, :, canonical_phase
                ]
    return counts, u4


def zero_sequence_audit(counts: np.ndarray, u4: np.ndarray, rich_sequences: np.ndarray) -> None:
    zero_id = next(
        index for index, row in enumerate(rich_sequences) if not np.any(row)
    )
    if int(counts[zero_id, 1, 0, 0]) != P:
        raise ArithmeticError("unit-profile count audit failed")
    if np.any(counts[zero_id, 1, 0, 1:]):
        raise ArithmeticError("zero sequence has a nonzero character phase")
    if int(u4[zero_id, 1, 0, 0]) != P * (P - 1):
        raise ArithmeticError("unit-profile U4 audit failed")
    if np.any(u4[zero_id, 1, 0, 1:]):
        raise ArithmeticError("zero-sequence U4 has a nonzero character phase")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--types", type=Path, help="optional saved type archive")
    parser.add_argument("--max-k", type=int, default=30)
    parser.add_argument("--ordinary-profiles", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.max_k < 0:
        raise ValueError("max-k must be nonnegative")

    started = time.monotonic()
    if args.types is None:
        _lookup, rich_sequences, table_sequences, descriptors = rich_profile_data()
    else:
        with np.load(args.types) as archive:
            rich_sequences = np.asarray(archive["rich_sequences"], dtype=np.uint8)
            table_sequences = np.asarray(
                archive["canonical_table_sequences"], dtype=np.uint8
            )
            descriptors = np.asarray(archive["affine_descriptors"], dtype=np.uint8)
    canonical_pairs = [
        distribution_tables(tuple(int(value) for value in row), args.max_k)
        for row in table_sequences
    ]
    canonical_count = np.stack([pair[0] for pair in canonical_pairs])
    canonical_u4 = np.stack([pair[1] for pair in canonical_pairs])
    counts, u4 = reconstruct_tables(canonical_count, canonical_u4, descriptors)
    zero_sequence_audit(counts, u4, rich_sequences)

    ordinary_audit = None
    if args.ordinary_profiles is not None:
        with np.load(args.ordinary_profiles) as archive:
            ordinary = np.asarray(archive["counts"], dtype=np.uint64)
        from r1_p11_profile_dual_tuple_gpu import histogram_lookup

        _ordinary_lookup, ordinary_histograms = histogram_lookup()
        histogram_ids = {
            histogram: index for index, histogram in enumerate(ordinary_histograms)
        }
        mismatches = []
        for rich_id, sequence in enumerate(rich_sequences):
            histogram = tuple(
                int(np.count_nonzero(sequence == value)) for value in range(P)
            )
            histogram_id = histogram_ids[histogram]
            if not np.array_equal(counts[rich_id], ordinary[histogram_id]):
                mismatches.append(rich_id)
                break
        if mismatches:
            raise ArithmeticError(
                f"rich count table disagrees with ordinary table at {mismatches[0]}"
            )
        ordinary_audit = len(rich_sequences)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        counts=counts,
        legendre_fourth=u4,
        rich_sequences=rich_sequences,
        canonical_table_sequences=table_sequences,
        affine_descriptors=descriptors,
    )
    report = {
        "experiment": "r1_p11_channel_profile_tables",
        "status": "complete_exact_count_and_legendre_fourth_tables",
        "p": P,
        "domain_affine_profile_types": int(len(rich_sequences)),
        "canonical_dynamic_programs": int(len(table_sequences)),
        "max_excess_parameter_k": args.max_k,
        "maximum_profile_energy": (P - 1) + 2 * args.max_k,
        "table_shape": list(counts.shape),
        "maximum_count_entry": int(counts.max()),
        "maximum_legendre_fourth_entry": int(u4.max()),
        "zero_sequence_unit_profile_audit": True,
        "ordinary_profile_tables_audited": ordinary_audit,
        "elapsed_seconds": time.monotonic() - started,
        "output": str(args.output),
        "output_sha256": sha256(args.output),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
