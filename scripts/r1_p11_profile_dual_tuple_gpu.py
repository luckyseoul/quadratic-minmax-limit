#!/usr/bin/env python3
"""Aggregate p=11 glue-dual words by six quartic value distributions.

Translation orbits and nonzero scalar conjugacy reduce 11^10 dual words to
21,437,340 representatives.  A V100 kernel maps each representative to the
sorted six-tuple of value-distribution IDs of its direction polynomials.
Radix sorting then aggregates the exact translation-orbit weights.  The
result is independent of the profile sum and can be reused for every theta
coefficient.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from r1_p11_profile_dual_orbits import (
    P,
    column_complement,
    kernel_bases,
    polynomial_histogram,
    projective_representatives,
    translation_image,
)


H4_COUNT = P**7
H3_COUNT = (P + 1) * P**5
H2_COUNT = ((P**3 - 1) // (P - 1)) * P**2
H1_COUNT = (P**4 - 1) // (P - 1)
TOTAL_REPRESENTATIVES = H4_COUNT + H3_COUNT + H2_COUNT + H1_COUNT


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def histogram_lookup() -> tuple[np.ndarray, list[tuple[int, ...]]]:
    identifiers: dict[tuple[int, ...], int] = {}
    lookup = np.empty(P**4, dtype=np.uint16)
    for c4 in range(P):
        for c3 in range(P):
            for c2 in range(P):
                for c1 in range(P):
                    histogram = polynomial_histogram((c1, c2, c3, c4))
                    identifier = identifiers.setdefault(histogram, len(identifiers))
                    packed = c1 + P * c2 + P**2 * c3 + P**3 * c4
                    lookup[packed] = identifier
    if len(identifiers) != 604:
        raise ArithmeticError("quartic value-distribution census changed")
    ordered = [()] * len(identifiers)
    for histogram, identifier in identifiers.items():
        ordered[identifier] = histogram
    return lookup, ordered


def orbit_descriptors() -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    bases = kernel_bases()
    basis_rows = []
    for degree in range(1, 5):
        padded = list(bases[degree]) + [[0] * 6] * (4 - len(bases[degree]))
        basis_rows.append(np.asarray(padded, dtype=np.uint8))
    basis_array = np.stack(basis_rows)

    h3_top = np.asarray(projective_representatives(2), dtype=np.uint8)
    h3_complement = []
    for top in map(tuple, h3_top.tolist()):
        complement = column_complement(translation_image(bases, 3, top), 3)
        if len(complement) != 1:
            raise ArithmeticError("degree-three translation complement changed")
        h3_complement.append(complement[0])

    h2_top = np.asarray(projective_representatives(3), dtype=np.uint8)
    h2_complement = []
    for top in map(tuple, h2_top.tolist()):
        complement = column_complement(translation_image(bases, 2, top), 4)
        if len(complement) != 2:
            raise ArithmeticError("degree-two translation complement changed")
        h2_complement.append(complement)

    h1_top = np.asarray(projective_representatives(4), dtype=np.uint8)
    return (
        basis_array,
        h3_top,
        np.asarray(h3_complement, dtype=np.uint8),
        h2_top,
        np.asarray(h2_complement, dtype=np.uint8),
        h1_top,
    )


def cpu_tuple_key(
    basis: np.ndarray,
    lookup: np.ndarray,
    c1: np.ndarray,
    c2: np.ndarray,
    c3: np.ndarray,
    c4: int,
) -> int:
    """Independent host implementation used to audit sampled GPU entries."""
    identifiers = []
    coordinates = (c1, c2, c3, np.asarray([c4], dtype=np.uint8))
    for direction in range(6):
        coefficients = [
            sum(
                int(coordinate) * int(basis[degree, offset, direction])
                for offset, coordinate in enumerate(block)
            )
            % P
            for degree, block in enumerate(coordinates)
        ]
        packed = coefficients[0] + P * (
            coefficients[1] + P * (coefficients[2] + P * coefficients[3])
        )
        identifiers.append(int(lookup[packed]))
    identifiers.sort()
    return sum(identifier << (10 * index) for index, identifier in enumerate(identifiers))


def cpu_orbit_entry(
    index: int,
    basis: np.ndarray,
    lookup: np.ndarray,
    h3_top: np.ndarray,
    h3_complement: np.ndarray,
    h2_top: np.ndarray,
    h2_complement: np.ndarray,
    h1_top: np.ndarray,
) -> tuple[int, int]:
    """Decode one global representative index without sharing CUDA logic."""
    if index < H4_COUNT:
        value = index
        c1 = np.asarray([(value // P**k) % P for k in range(4)], dtype=np.uint8)
        value //= P**4
        c2 = np.asarray([(value // P**k) % P for k in range(3)], dtype=np.uint8)
        return cpu_tuple_key(basis, lookup, c1, c2, np.zeros(2, dtype=np.uint8), 1), P**2

    index -= H4_COUNT
    if index < H3_COUNT:
        value = index
        c1 = np.asarray([(value // P**k) % P for k in range(4)], dtype=np.uint8)
        value //= P**4
        alpha = value % P
        cls = value // P
        c2 = alpha * h3_complement[cls] % P
        return cpu_tuple_key(basis, lookup, c1, c2, h3_top[cls], 0), P**2

    index -= H3_COUNT
    if index < H2_COUNT:
        alpha = index % P
        beta = (index // P) % P
        cls = index // P**2
        c1 = (alpha * h2_complement[cls, 0] + beta * h2_complement[cls, 1]) % P
        return cpu_tuple_key(
            basis, lookup, c1, h2_top[cls], np.zeros(2, dtype=np.uint8), 0
        ), P**2

    index -= H2_COUNT
    if index >= H1_COUNT:
        raise IndexError("representative index outside all four strata")
    return cpu_tuple_key(
        basis,
        lookup,
        h1_top[index],
        np.zeros(3, dtype=np.uint8),
        np.zeros(2, dtype=np.uint8),
        0,
    ), 1


def audit_generated_entries(
    keys: object,
    weights: object,
    descriptors: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    lookup: np.ndarray,
    samples_per_stratum: int,
) -> int:
    """Compare deterministic random GPU entries with an independent CPU decoder."""
    import cupy as cp

    basis, h3_top, h3_complement, h2_top, h2_complement, h1_top = descriptors
    rng = np.random.default_rng(0x11C0DE)
    starts_and_counts = (
        (0, H4_COUNT),
        (H4_COUNT, H3_COUNT),
        (H4_COUNT + H3_COUNT, H2_COUNT),
        (H4_COUNT + H3_COUNT + H2_COUNT, H1_COUNT),
    )
    sampled: set[int] = set()
    for start, count in starts_and_counts:
        sampled.update((start, start + count // 2, start + count - 1))
        draw_count = min(samples_per_stratum, count)
        sampled.update(start + int(value) for value in rng.choice(count, draw_count, replace=False))
    indices = np.asarray(sorted(sampled), dtype=np.int64)
    gpu_keys = cp.asnumpy(keys[indices])
    gpu_weights = cp.asnumpy(weights[indices])
    for position, index in enumerate(indices):
        expected_key, expected_weight = cpu_orbit_entry(
            int(index), basis, lookup, h3_top, h3_complement, h2_top, h2_complement, h1_top
        )
        if int(gpu_keys[position]) != expected_key or int(gpu_weights[position]) != expected_weight:
            raise ArithmeticError(
                f"GPU orbit entry mismatch at {index}: "
                f"got ({int(gpu_keys[position])}, {int(gpu_weights[position])}), "
                f"expected ({expected_key}, {expected_weight})"
            )
    return len(indices)


CUDA_SOURCE = r"""
extern "C" {

__device__ __forceinline__ unsigned int mod11(unsigned int value) {
    return value % 11u;
}

__device__ __forceinline__ unsigned long long tuple_key(
    const unsigned char *basis,
    const unsigned short *lookup,
    const unsigned int c1[4],
    const unsigned int c2[3],
    const unsigned int c3[2],
    unsigned int c4
) {
    unsigned short ids[6];
    for (int j = 0; j < 6; ++j) {
        unsigned int values[4] = {0u, 0u, 0u, 0u};
        for (int k = 0; k < 4; ++k)
            values[0] += c1[k] * (unsigned int)basis[((0 * 4 + k) * 6) + j];
        for (int k = 0; k < 3; ++k)
            values[1] += c2[k] * (unsigned int)basis[((1 * 4 + k) * 6) + j];
        for (int k = 0; k < 2; ++k)
            values[2] += c3[k] * (unsigned int)basis[((2 * 4 + k) * 6) + j];
        values[3] = c4 * (unsigned int)basis[((3 * 4 + 0) * 6) + j];
        unsigned int packed = mod11(values[0])
            + 11u * mod11(values[1])
            + 121u * mod11(values[2])
            + 1331u * mod11(values[3]);
        ids[j] = lookup[packed];
    }
    for (int i = 1; i < 6; ++i) {
        unsigned short value = ids[i];
        int j = i - 1;
        while (j >= 0 && ids[j] > value) {
            ids[j + 1] = ids[j];
            --j;
        }
        ids[j + 1] = value;
    }
    unsigned long long key = 0ull;
    for (int j = 0; j < 6; ++j)
        key |= ((unsigned long long)ids[j]) << (10 * j);
    return key;
}

__global__ void generate_h4(
    unsigned long long count,
    const unsigned char *basis,
    const unsigned short *lookup,
    unsigned long long *keys,
    unsigned short *weights
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    unsigned long long value = index;
    unsigned int c1[4], c2[3], c3[2] = {0u, 0u};
    for (int k = 0; k < 4; ++k) { c1[k] = value % 11ull; value /= 11ull; }
    for (int k = 0; k < 3; ++k) { c2[k] = value % 11ull; value /= 11ull; }
    keys[index] = tuple_key(basis, lookup, c1, c2, c3, 1u);
    weights[index] = 121u;
}

__global__ void generate_h3(
    unsigned long long count,
    unsigned long long output_offset,
    const unsigned char *basis,
    const unsigned short *lookup,
    const unsigned char *tops,
    const unsigned char *complements,
    unsigned long long *keys,
    unsigned short *weights
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    unsigned long long value = index;
    unsigned int c1[4], c2[3], c3[2];
    for (int k = 0; k < 4; ++k) { c1[k] = value % 11ull; value /= 11ull; }
    unsigned int alpha = value % 11ull; value /= 11ull;
    unsigned int cls = (unsigned int)value;
    for (int k = 0; k < 3; ++k)
        c2[k] = alpha * (unsigned int)complements[cls * 3 + k] % 11u;
    for (int k = 0; k < 2; ++k) c3[k] = tops[cls * 2 + k];
    unsigned long long out = output_offset + index;
    keys[out] = tuple_key(basis, lookup, c1, c2, c3, 0u);
    weights[out] = 121u;
}

__global__ void generate_h2(
    unsigned long long count,
    unsigned long long output_offset,
    const unsigned char *basis,
    const unsigned short *lookup,
    const unsigned char *tops,
    const unsigned char *complements,
    unsigned long long *keys,
    unsigned short *weights
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    unsigned long long value = index;
    unsigned int alpha = value % 11ull; value /= 11ull;
    unsigned int beta = value % 11ull; value /= 11ull;
    unsigned int cls = (unsigned int)value;
    unsigned int c1[4], c2[3], c3[2] = {0u, 0u};
    for (int k = 0; k < 4; ++k)
        c1[k] = (alpha * (unsigned int)complements[(cls * 2 + 0) * 4 + k]
               + beta  * (unsigned int)complements[(cls * 2 + 1) * 4 + k]) % 11u;
    for (int k = 0; k < 3; ++k) c2[k] = tops[cls * 3 + k];
    unsigned long long out = output_offset + index;
    keys[out] = tuple_key(basis, lookup, c1, c2, c3, 0u);
    weights[out] = 121u;
}

__global__ void generate_h1(
    unsigned long long count,
    unsigned long long output_offset,
    const unsigned char *basis,
    const unsigned short *lookup,
    const unsigned char *tops,
    unsigned long long *keys,
    unsigned short *weights
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    unsigned int c1[4], c2[3] = {0u, 0u, 0u}, c3[2] = {0u, 0u};
    for (int k = 0; k < 4; ++k) c1[k] = tops[index * 4 + k];
    unsigned long long out = output_offset + index;
    keys[out] = tuple_key(basis, lookup, c1, c2, c3, 0u);
    weights[out] = 1u;
}

}
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--audit-samples-per-stratum", type=int, default=64)
    args = parser.parse_args()

    import cupy as cp

    cp.cuda.Device(args.device).use()
    lookup, histograms = histogram_lookup()
    descriptors = orbit_descriptors()
    basis, h3_top, h3_complement, h2_top, h2_complement, h1_top = descriptors
    arrays = [
        cp.asarray(value)
        for value in (basis, lookup, h3_top, h3_complement, h2_top, h2_complement, h1_top)
    ]
    d_basis, d_lookup, d_h3_top, d_h3_comp, d_h2_top, d_h2_comp, d_h1_top = arrays
    keys = cp.empty(TOTAL_REPRESENTATIVES, dtype=cp.uint64)
    weights = cp.empty(TOTAL_REPRESENTATIVES, dtype=cp.uint16)
    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++11",))
    kernels = {
        name: module.get_function(name)
        for name in ("generate_h4", "generate_h3", "generate_h2", "generate_h1")
    }
    block = 256
    started = time.monotonic()
    kernels["generate_h4"](
        ((H4_COUNT + block - 1) // block,), (block,),
        (np.uint64(H4_COUNT), d_basis, d_lookup, keys, weights),
    )
    kernels["generate_h3"](
        ((H3_COUNT + block - 1) // block,), (block,),
        (
            np.uint64(H3_COUNT), np.uint64(H4_COUNT), d_basis, d_lookup,
            d_h3_top, d_h3_comp, keys, weights,
        ),
    )
    offset_h2 = H4_COUNT + H3_COUNT
    kernels["generate_h2"](
        ((H2_COUNT + block - 1) // block,), (block,),
        (
            np.uint64(H2_COUNT), np.uint64(offset_h2), d_basis, d_lookup,
            d_h2_top, d_h2_comp, keys, weights,
        ),
    )
    offset_h1 = offset_h2 + H2_COUNT
    kernels["generate_h1"](
        ((H1_COUNT + block - 1) // block,), (block,),
        (
            np.uint64(H1_COUNT), np.uint64(offset_h1), d_basis, d_lookup,
            d_h1_top, keys, weights,
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    generation_seconds = time.monotonic() - started
    audited_entries = audit_generated_entries(
        keys, weights, descriptors, lookup, args.audit_samples_per_stratum
    )

    started = time.monotonic()
    order = cp.argsort(keys)
    sorted_keys = keys[order]
    sorted_weights = weights[order].astype(cp.uint64)
    del keys, weights, order
    change = cp.empty(TOTAL_REPRESENTATIVES, dtype=cp.bool_)
    change[0] = True
    change[1:] = sorted_keys[1:] != sorted_keys[:-1]
    starts = cp.nonzero(change)[0]
    unique_keys = sorted_keys[starts]
    prefix = cp.cumsum(sorted_weights, dtype=cp.uint64)
    ends = cp.concatenate((starts[1:], cp.asarray([TOTAL_REPRESENTATIVES], dtype=starts.dtype)))
    unique_weights = prefix[ends - 1]
    if len(starts) > 1:
        unique_weights[1:] -= prefix[starts[1:] - 1]
    cp.cuda.runtime.deviceSynchronize()
    aggregation_seconds = time.monotonic() - started

    host_keys = cp.asnumpy(unique_keys)
    host_weights = cp.asnumpy(unique_weights)
    representative_weight = int(host_weights.sum(dtype=np.uint64))
    if 1 + (P - 1) * representative_weight != P**10:
        raise ArithmeticError("aggregated orbit weights do not reconstruct 11^10")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        keys=host_keys,
        weights=host_weights,
        histograms=np.asarray(histograms, dtype=np.uint8),
    )
    gpu_name = cp.cuda.runtime.getDeviceProperties(args.device)["name"]
    if isinstance(gpu_name, bytes):
        gpu_name = gpu_name.decode()
    report = {
        "experiment": "r1_p11_profile_dual_tuple_gpu",
        "status": "complete_exact_orbit_tuple_aggregation",
        "p": P,
        "device": args.device,
        "gpu_name": str(gpu_name),
        "dual_codewords": P**10,
        "translation_scalar_representatives": TOTAL_REPRESENTATIVES,
        "unique_sorted_value_distribution_tuples": int(len(host_keys)),
        "translation_orbit_weight_sum_per_scalar_class": representative_weight,
        "codewords_reconstructed": 1 + (P - 1) * representative_weight,
        "quartic_value_distribution_types": len(histograms),
        "cpu_gpu_entries_audited": audited_entries,
        "generation_seconds": generation_seconds,
        "aggregation_seconds": aggregation_seconds,
        "output": str(args.output),
        "output_sha256": sha256(args.output),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
