#!/usr/bin/env python3
"""Aggregate p=11 glue words by channel-sensitive profile types.

This is the square-circle analogue of ``r1_p11_profile_dual_tuple_gpu``.
The ordinary counter stores only six value-distribution IDs.  Here each
quartic profile is instead reduced by affine permutations of its input,
which preserve the Legendre-convolution fourth statistic U4.  There are
1,007 such IDs, so the sorted six-tuple still fits exactly in one 60-bit
integer.  Translation/scalar orbit weights and all sampled GPU entries are
checked independently on the CPU.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from r1_p11_channel_profile_types import rich_profile_data
from r1_p11_profile_dual_orbits import P
from r1_p11_profile_dual_tuple_gpu import (
    CUDA_SOURCE,
    H1_COUNT,
    H2_COUNT,
    H3_COUNT,
    H4_COUNT,
    TOTAL_REPRESENTATIVES,
    audit_generated_entries,
    orbit_descriptors,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--audit-samples-per-stratum", type=int, default=64)
    args = parser.parse_args()

    import cupy as cp

    cp.cuda.Device(args.device).use()
    lookup, rich_sequences, canonical_tables, affine_descriptors = rich_profile_data()
    if int(lookup.max()) >= 1 << 10:
        raise ArithmeticError("rich profile ID no longer fits in ten bits")
    descriptors = orbit_descriptors()
    basis, h3_top, h3_complement, h2_top, h2_complement, h1_top = descriptors
    arrays = [
        cp.asarray(value)
        for value in (
            basis,
            lookup,
            h3_top,
            h3_complement,
            h2_top,
            h2_complement,
            h1_top,
        )
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
        ((H4_COUNT + block - 1) // block,),
        (block,),
        (np.uint64(H4_COUNT), d_basis, d_lookup, keys, weights),
    )
    kernels["generate_h3"](
        ((H3_COUNT + block - 1) // block,),
        (block,),
        (
            np.uint64(H3_COUNT),
            np.uint64(H4_COUNT),
            d_basis,
            d_lookup,
            d_h3_top,
            d_h3_comp,
            keys,
            weights,
        ),
    )
    offset_h2 = H4_COUNT + H3_COUNT
    kernels["generate_h2"](
        ((H2_COUNT + block - 1) // block,),
        (block,),
        (
            np.uint64(H2_COUNT),
            np.uint64(offset_h2),
            d_basis,
            d_lookup,
            d_h2_top,
            d_h2_comp,
            keys,
            weights,
        ),
    )
    offset_h1 = offset_h2 + H2_COUNT
    kernels["generate_h1"](
        ((H1_COUNT + block - 1) // block,),
        (block,),
        (
            np.uint64(H1_COUNT),
            np.uint64(offset_h1),
            d_basis,
            d_lookup,
            d_h1_top,
            keys,
            weights,
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
    ends = cp.concatenate(
        (starts[1:], cp.asarray([TOTAL_REPRESENTATIVES], dtype=starts.dtype))
    )
    unique_weights = prefix[ends - 1]
    if len(starts) > 1:
        unique_weights[1:] -= prefix[starts[1:] - 1]
    cp.cuda.runtime.deviceSynchronize()
    aggregation_seconds = time.monotonic() - started

    host_keys = cp.asnumpy(unique_keys)
    host_weights = cp.asnumpy(unique_weights)
    representative_weight = int(host_weights.sum(dtype=np.uint64))
    if 1 + (P - 1) * representative_weight != P**10:
        raise ArithmeticError("aggregated rich orbit weights do not reconstruct 11^10")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        keys=host_keys,
        weights=host_weights,
        rich_sequences=rich_sequences,
        canonical_table_sequences=canonical_tables,
        affine_descriptors=affine_descriptors,
    )
    gpu_name = cp.cuda.runtime.getDeviceProperties(args.device)["name"]
    if isinstance(gpu_name, bytes):
        gpu_name = gpu_name.decode()
    report = {
        "experiment": "r1_p11_channel_dual_tuple_gpu",
        "status": "complete_exact_channel_tuple_aggregation",
        "p": P,
        "device": args.device,
        "gpu_name": str(gpu_name),
        "dual_codewords": P**10,
        "translation_scalar_representatives": TOTAL_REPRESENTATIVES,
        "domain_affine_profile_types": int(len(rich_sequences)),
        "domain_output_affine_table_types": int(len(canonical_tables)),
        "unique_sorted_channel_profile_tuples": int(len(host_keys)),
        "translation_orbit_weight_sum_per_scalar_class": representative_weight,
        "codewords_reconstructed": 1 + (P - 1) * representative_weight,
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
