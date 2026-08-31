#!/usr/bin/env python3
"""Exact OpenCL companion for the positive p=7, z=7 semigroup case join.

The mathematical construction is delegated unchanged to
``p7_infinity7_positive_z7_semigroup_case_join_gpu.py`` and its CPU source.
This file replaces only the final CUDA Minkowski-support engine with OpenCL.
It is intended for Nuka's AMD RX 9070 XT (gfx1201), but selects any requested
OpenCL GPU that advertises ``cl_khr_int64_base_atomics``.

All six characteristic-three coordinates and an explicit width-five or
width-six subset of the 21 characteristic-seven coordinates are encoded in
one mixed-radix uint64 key.  An OpenCL hash table uses exact uint64 compare-
and-exchange plus linear probing.  Hash collisions are resolved by full-key
comparison; no Bloom filter, fingerprint, or probabilistic membership test is
used.  State, memory, allocation, or probe caps are explicit nondecision skips.

The self-audit compares complete CPU and OpenCL sets, forces a genuine initial
hash-slot collision, exercises present/absent, chunk-invariance, same-row
cross-prime, state-cap, and memory-cap traps, and checks one small real case
elementwise.  Production case sharding is the audited deterministic sharding
from the CUDA companion.  JSON output is replaced atomically.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join as cpu_join  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_join_gpu as case_gpu  # noqa: E402
import p7_infinity7_positive_z7_torsion_support_gpu as cuda_join  # noqa: E402


EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_join_opencl"
SUPPORTED_K = case_gpu.SUPPORTED_K
EMPTY_KEY = np.uint64(np.iinfo(np.uint64).max)
THREADS_PER_GROUP = 256
DEFAULT_PAIR_CHUNK_CAP = case_gpu.DEFAULT_PAIR_CHUNK_CAP
DEFAULT_OPENCL_MEMORY_CAP_MIB = 12_000
MAX_DIRECTION_SUPPORT_WORKERS = 1


OPENCL_SOURCE = r"""
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

#define EMPTY_KEY ((ulong)0xffffffffffffffffUL)

inline ulong mix64(ulong x) {
    x ^= x >> 30;
    x *= (ulong)0xbf58476d1ce4e5b9UL;
    x ^= x >> 27;
    x *= (ulong)0x94d049bb133111ebUL;
    x ^= x >> 31;
    return x;
}

inline ulong mixed_sum_key(
    ulong left,
    ulong right,
    __global const uint *moduli,
    int width
) {
    ulong key = 0UL;
    ulong place = 1UL;
    for (int coordinate = 0; coordinate < width; ++coordinate) {
        const uint modulus = moduli[coordinate];
        const uint left_digit = (uint)(left % (ulong)modulus);
        const uint right_digit = (uint)(right % (ulong)modulus);
        const uint digit = (left_digit + right_digit) % modulus;
        key += place * (ulong)digit;
        place *= (ulong)modulus;
        left /= (ulong)modulus;
        right /= (ulong)modulus;
    }
    return key;
}

__kernel void insert_mixed_sums(
    __global const ulong *states,
    ulong state_count,
    __global const ulong *support,
    ulong support_count,
    __global const uint *moduli,
    int width,
    ulong pair_offset,
    ulong pair_count,
    __global volatile ulong *table,
    ulong table_mask,
    uint state_cap,
    __global volatile uint *counters
) {
    const ulong work_item = (ulong)get_global_id(0);
    if (work_item >= pair_count) return;

    const ulong pair_index = pair_offset + work_item;
    const ulong left_index = pair_index / support_count;
    const ulong right_index = pair_index % support_count;
    if (left_index >= state_count) return;

    const ulong key = mixed_sum_key(
        states[left_index], support[right_index], moduli, width
    );
    ulong slot = mix64(key) & table_mask;
    ulong probes = 0UL;
    while (probes <= table_mask) {
        const ulong previous = atom_cmpxchg(table + slot, EMPTY_KEY, key);
        if (previous == EMPTY_KEY) {
            const uint old_count = atomic_inc(counters + 0);
            if (old_count >= state_cap) atomic_xchg(counters + 1, 1U);
            return;
        }
        if (previous == key) return;
        slot = (slot + 1UL) & table_mask;
        ++probes;
    }
    atomic_xchg(counters + 2, 1U);
}

__kernel void compact_table(
    __global const ulong *table,
    ulong table_size,
    __global ulong *output,
    __global volatile uint *output_count
) {
    const ulong index = (ulong)get_global_id(0);
    if (index >= table_size) return;
    const ulong key = table[index];
    if (key != EMPTY_KEY) {
        const uint destination = atomic_inc(output_count);
        output[destination] = key;
    }
}

__kernel void contains_key(
    __global const ulong *table,
    ulong table_mask,
    ulong key,
    __global volatile uint *result
) {
    if (get_global_id(0) != 0) return;
    ulong slot = mix64(key) & table_mask;
    ulong probes = 0UL;
    while (probes <= table_mask) {
        const ulong observed = table[slot];
        if (observed == key) {
            result[0] = 1U;
            return;
        }
        if (observed == EMPTY_KEY) return;
        slot = (slot + 1UL) & table_mask;
        ++probes;
    }
    result[1] = 1U;
}
"""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def mix64_cpu(value: int) -> int:
    mask = (1 << 64) - 1
    value &= mask
    value ^= value >> 30
    value = value * 0xBF58476D1CE4E5B9 & mask
    value ^= value >> 27
    value = value * 0x94D049BB133111EB & mask
    value ^= value >> 31
    return value & mask


def load_opencl(
    platform_index: int,
    device_index: int,
    device_name_contains: str | None,
) -> dict[str, Any]:
    """Import PyOpenCL lazily and build the exact integer kernels."""
    try:
        import pyopencl as cl  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - host packaging dependent
        raise RuntimeError("pyopencl is required for the OpenCL companion") from exc

    platforms = cl.get_platforms()
    require(bool(platforms), "no OpenCL platform found")
    selected_platform = None
    selected_device = None
    if device_name_contains:
        needle = device_name_contains.casefold()
        matches = [
            (platform, device)
            for platform in platforms
            for device in platform.get_devices(device_type=cl.device_type.GPU)
            if needle in str(device.name).casefold()
        ]
        require(len(matches) == 1, "device-name selector must match exactly one GPU")
        selected_platform, selected_device = matches[0]
    else:
        require(0 <= platform_index < len(platforms), "OpenCL platform index out of range")
        selected_platform = platforms[platform_index]
        devices = selected_platform.get_devices(device_type=cl.device_type.GPU)
        require(bool(devices), "selected OpenCL platform has no GPU")
        require(0 <= device_index < len(devices), "OpenCL GPU index out of range")
        selected_device = devices[device_index]

    extensions = set(str(selected_device.extensions).split())
    require(
        "cl_khr_int64_base_atomics" in extensions,
        "selected device lacks cl_khr_int64_base_atomics",
    )
    context = cl.Context([selected_device])
    queue = cl.CommandQueue(
        context,
        properties=cl.command_queue_properties.PROFILING_ENABLE,
    )
    program = cl.Program(context, OPENCL_SOURCE).build(options=["-cl-std=CL1.2"])
    kernels = {
        name: cl.Kernel(program, name)
        for name in ("insert_mixed_sums", "compact_table", "contains_key")
    }
    maximum_group = int(selected_device.max_work_group_size)
    local_size = 1 << int(math.log2(min(THREADS_PER_GROUP, maximum_group)))
    return {
        "cl": cl,
        "platform": selected_platform,
        "device": selected_device,
        "context": context,
        "queue": queue,
        "program": program,
        "kernels": kernels,
        "local_size": local_size,
        "audit": {
            "pyopencl_version": getattr(cl, "VERSION_TEXT", "unknown"),
            "platform_name": str(selected_platform.name),
            "platform_vendor": str(selected_platform.vendor),
            "device_name": str(selected_device.name),
            "device_vendor": str(selected_device.vendor),
            "device_version": str(selected_device.version),
            "driver_version": str(selected_device.driver_version),
            "opencl_c_version": str(selected_device.opencl_c_version),
            "global_memory_bytes": int(selected_device.global_mem_size),
            "maximum_single_allocation_bytes": int(selected_device.max_mem_alloc_size),
            "maximum_work_group_size": maximum_group,
            "selected_local_size": local_size,
            "cl_khr_int64_base_atomics_advertised": True,
            "opencl_kernel_sha256": hashlib.sha256(
                OPENCL_SOURCE.encode("utf-8")
            ).hexdigest(),
            "hash_collisions_resolved_by_exact_uint64_comparison_and_linear_probing": True,
            "probabilistic_filter_used": False,
        },
    }


def launch_1d(runtime: dict[str, Any], kernel_name: str, work_items: int, *args: Any) -> Any:
    require(work_items > 0, "OpenCL launch must contain work")
    local_size = int(runtime["local_size"])
    global_size = ((work_items + local_size - 1) // local_size) * local_size
    kernel = runtime["kernels"][kernel_name]
    kernel.set_args(*args)
    return runtime["cl"].enqueue_nd_range_kernel(
        runtime["queue"], kernel, (global_size,), (local_size,)
    )


def opencl_skip(
    status: str,
    normalized_sizes: list[int],
    group_size: int,
    state_cap: int,
    completed_state_sizes: list[int],
    completed_factors: int,
    pair_candidates_launched: int,
    **extra: object,
) -> dict:
    result = cuda_join.convolution_skip(
        status,
        normalized_sizes,
        group_size,
        state_cap,
        completed_state_sizes,
        completed_factors,
        pair_candidates_launched,
        **extra,
    )
    result["opencl_hash_set_is_exact_not_probabilistic"] = True
    return result


def opencl_exact_support_convolution(
    runtime: dict[str, Any],
    supports: list[np.ndarray],
    target: np.ndarray,
    moduli: tuple[int, ...],
    state_cap: int,
    pair_chunk_cap: int,
    opencl_memory_cap_bytes: int,
    return_final_codes: bool = False,
) -> tuple[dict, np.ndarray | None]:
    """Compute an exact capped support convolution with an OpenCL hash set."""
    started = time.perf_counter()
    require(state_cap > 0, "state cap must be positive")
    require(state_cap < 2**32, "OpenCL state cap does not fit the exact uint32 counter")
    require(pair_chunk_cap > 0, "pair chunk cap must be positive")
    require(opencl_memory_cap_bytes > 0, "OpenCL memory cap must be positive")
    require(len(target) == len(moduli), "projected target width changed")
    group_size = cuda_join.group_size_for(moduli)
    normalized_codes = [cuda_join.unique_support_codes(rows, moduli) for rows in supports]
    normalized_sizes = [len(codes) for codes in normalized_codes]
    require(all(size > 0 for size in normalized_sizes), "empty exact catalog support")
    target_code = int(
        cuda_join.encode_rows_cpu(
            np.asarray(target, dtype=np.uint8)[None, :], moduli
        )[0]
    )

    cl = runtime["cl"]
    context = runtime["context"]
    queue = runtime["queue"]
    flags = cl.mem_flags
    device = runtime["device"]
    effective_memory_limit = min(
        opencl_memory_cap_bytes, int(0.85 * int(device.global_mem_size))
    )
    maximum_allocation = int(device.max_mem_alloc_size)
    moduli_host = np.ascontiguousarray(moduli, dtype=np.uint32)
    moduli_buffer = cl.Buffer(
        context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=moduli_host
    )
    states_host = np.asarray([0], dtype=np.uint64)
    states_buffer = cl.Buffer(
        context, flags.READ_WRITE | flags.COPY_HOST_PTR, hostbuf=states_host
    )
    state_count = 1
    final_table = None
    final_table_size = 0
    completed_state_sizes: list[int] = []
    factor_audits = []
    pair_candidates_launched = 0
    peak_estimated_bytes = 0

    for factor_index, support_host in enumerate(normalized_codes):
        factor_started = time.perf_counter()
        if final_table is not None:
            final_table.release()
            final_table = None
        support_host = np.ascontiguousarray(support_host, dtype=np.uint64)
        support_count = len(support_host)
        total_pairs = state_count * support_count
        launch_chunk = min(pair_chunk_cap, total_pairs)
        maximum_insertions = min(
            group_size,
            total_pairs,
            state_cap + launch_chunk + 1,
        )
        table_size = cuda_join.next_power_of_two(max(8, 2 * maximum_insertions))
        table_bytes = table_size * 8
        compact_bytes = maximum_insertions * 8
        estimated_peak_bytes = (
            table_bytes
            + compact_bytes
            + (state_count + support_count) * 8
            + len(moduli) * 4
            + 32
        )
        peak_estimated_bytes = max(peak_estimated_bytes, estimated_peak_bytes)
        largest_allocation = max(table_bytes, compact_bytes, state_count * 8, support_count * 8)
        if estimated_peak_bytes > effective_memory_limit or largest_allocation > maximum_allocation:
            result = opencl_skip(
                "skipped_opencl_memory_cap",
                normalized_sizes,
                group_size,
                state_cap,
                completed_state_sizes,
                factor_index,
                pair_candidates_launched,
                required_estimated_peak_bytes=estimated_peak_bytes,
                configured_opencl_memory_cap_bytes=opencl_memory_cap_bytes,
                effective_opencl_memory_limit_bytes=effective_memory_limit,
                maximum_single_allocation_bytes=maximum_allocation,
                required_largest_single_allocation_bytes=largest_allocation,
                skipped_factor_pair_count=total_pairs,
            )
            result["elapsed_seconds"] = time.perf_counter() - started
            return result, None

        try:
            support_buffer = cl.Buffer(
                context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=support_host
            )
            table_buffer = cl.Buffer(context, flags.READ_WRITE, size=table_bytes)
            counters_host = np.zeros(3, dtype=np.uint32)
            counters_buffer = cl.Buffer(
                context, flags.READ_WRITE | flags.COPY_HOST_PTR, hostbuf=counters_host
            )
            cl.enqueue_fill_buffer(
                queue, table_buffer, EMPTY_KEY, 0, table_bytes
            ).wait()
        except cl.MemoryError:
            result = opencl_skip(
                "skipped_opencl_allocation_failure",
                normalized_sizes,
                group_size,
                state_cap,
                completed_state_sizes,
                factor_index,
                pair_candidates_launched,
                required_estimated_peak_bytes=estimated_peak_bytes,
                configured_opencl_memory_cap_bytes=opencl_memory_cap_bytes,
                skipped_factor_pair_count=total_pairs,
            )
            result["elapsed_seconds"] = time.perf_counter() - started
            return result, None

        factor_launched = 0
        for offset in range(0, total_pairs, pair_chunk_cap):
            count = min(pair_chunk_cap, total_pairs - offset)
            event = launch_1d(
                runtime,
                "insert_mixed_sums",
                count,
                states_buffer,
                np.uint64(state_count),
                support_buffer,
                np.uint64(support_count),
                moduli_buffer,
                np.int32(len(moduli)),
                np.uint64(offset),
                np.uint64(count),
                table_buffer,
                np.uint64(table_size - 1),
                np.uint32(state_cap),
                counters_buffer,
            )
            event.wait()
            factor_launched += count
            pair_candidates_launched += count
            cl.enqueue_copy(queue, counters_host, counters_buffer).wait()
            if int(counters_host[2]):
                result = opencl_skip(
                    "skipped_opencl_hash_table_capacity",
                    normalized_sizes,
                    group_size,
                    state_cap,
                    completed_state_sizes,
                    factor_index,
                    pair_candidates_launched,
                    hash_table_entries=table_size,
                    successful_distinct_insertions=int(counters_host[0]),
                    skipped_factor_pair_count=total_pairs,
                    skipped_factor_pairs_launched=factor_launched,
                )
                result["elapsed_seconds"] = time.perf_counter() - started
                return result, None
            if int(counters_host[1]):
                inserted = int(counters_host[0])
                require(inserted >= state_cap + 1, "overflow lacks cap-plus-one witness")
                result = opencl_skip(
                    "skipped_state_cap",
                    normalized_sizes,
                    group_size,
                    state_cap,
                    completed_state_sizes,
                    factor_index,
                    pair_candidates_launched,
                    first_capped_factor_index=factor_index,
                    distinct_states_lower_bound_at_skip=state_cap + 1,
                    successful_distinct_insertions_before_stop=inserted,
                    hash_table_entries=table_size,
                    skipped_factor_pair_count=total_pairs,
                    skipped_factor_pairs_launched=factor_launched,
                )
                result["elapsed_seconds"] = time.perf_counter() - started
                return result, None

        inserted = int(counters_host[0])
        require(inserted <= state_cap, "completed factor exceeded state cap")
        require(inserted <= maximum_insertions, "compact allocation bound failed")
        compact_buffer = cl.Buffer(context, flags.READ_WRITE, size=max(8, inserted * 8))
        compact_count_host = np.zeros(1, dtype=np.uint32)
        compact_count_buffer = cl.Buffer(
            context,
            flags.READ_WRITE | flags.COPY_HOST_PTR,
            hostbuf=compact_count_host,
        )
        compact_event = launch_1d(
            runtime,
            "compact_table",
            table_size,
            table_buffer,
            np.uint64(table_size),
            compact_buffer,
            compact_count_buffer,
        )
        compact_event.wait()
        cl.enqueue_copy(queue, compact_count_host, compact_count_buffer).wait()
        require(int(compact_count_host[0]) == inserted, "OpenCL hash occupancy/count mismatch")

        completed_state_sizes.append(inserted)
        factor_audits.append(
            {
                "factor_index": factor_index,
                "input_state_count": state_count,
                "exact_support_size": support_count,
                "pair_count": total_pairs,
                "pair_candidates_launched": factor_launched,
                "output_distinct_state_count": inserted,
                "hash_table_entries": table_size,
                "estimated_peak_bytes": estimated_peak_bytes,
                "elapsed_seconds": time.perf_counter() - factor_started,
            }
        )
        states_buffer.release()
        states_buffer = compact_buffer
        state_count = inserted
        support_buffer.release()
        counters_buffer.release()
        compact_count_buffer.release()
        final_table = table_buffer
        final_table_size = table_size

        if inserted == group_size:
            completed_state_sizes.extend(
                [group_size] * (len(normalized_codes) - factor_index - 1)
            )
            factor_audits.append(
                {
                    "after_factor_index": factor_index,
                    "exact_full_group_saturation_shortcut": True,
                    "remaining_nonempty_factors": len(normalized_codes) - factor_index - 1,
                }
            )
            break

    require(final_table is not None and final_table_size > 0, "final hash table missing")
    membership_host = np.zeros(2, dtype=np.uint32)
    membership_buffer = cl.Buffer(
        context, flags.READ_WRITE | flags.COPY_HOST_PTR, hostbuf=membership_host
    )
    membership_event = launch_1d(
        runtime,
        "contains_key",
        1,
        final_table,
        np.uint64(final_table_size - 1),
        np.uint64(target_code),
        membership_buffer,
    )
    membership_event.wait()
    cl.enqueue_copy(queue, membership_host, membership_buffer).wait()
    require(not int(membership_host[1]), "final OpenCL membership probe exhausted table")
    present = bool(membership_host[0])

    final_codes_host = np.empty(state_count, dtype=np.uint64)
    cl.enqueue_copy(queue, final_codes_host, states_buffer).wait()
    final_codes = np.sort(final_codes_host)
    require(
        len(final_codes) <= 1 or bool(np.all(final_codes[1:] > final_codes[:-1])),
        "final OpenCL support is not an exact set",
    )
    require(
        not len(final_codes) or int(final_codes[-1]) < group_size,
        "OpenCL hash produced a code outside the projected group",
    )
    index = int(np.searchsorted(final_codes, np.uint64(target_code)))
    host_present = bool(index < len(final_codes) and int(final_codes[index]) == target_code)
    require(host_present == present, "OpenCL hash query and compacted set disagree")

    result = {
        "completed_exact_convolution": True,
        "decision_status": (
            "necessary_only_projected_target_present"
            if present
            else "rigorous_projected_support_rejection"
        ),
        "rigorous_rejection": not present,
        "projected_target_present": present,
        "necessary_only_projected_survivor": present,
        "completed_catalog_factors": len(normalized_codes),
        "state_cap": state_cap,
        "finite_group_size": group_size,
        "exact_support_sizes": normalized_sizes,
        "completed_state_sizes": completed_state_sizes,
        "final_state_count": len(final_codes),
        "final_support_sha256_uint64": cuda_join.array_sha256(final_codes),
        "target_code": target_code,
        "pair_candidates_launched": pair_candidates_launched,
        "factor_audits": factor_audits,
        "peak_estimated_bytes": peak_estimated_bytes,
        "full_group_saturation_shortcut_is_exact": bool(len(final_codes) == group_size),
        "opencl_hash_set_is_exact_not_probabilistic": True,
        "elapsed_seconds": time.perf_counter() - started,
    }
    membership_buffer.release()
    states_buffer.release()
    final_table.release()
    moduli_buffer.release()
    return result, final_codes if return_final_codes else None


def collision_pair(group_size: int, table_mask: int) -> tuple[int, int, int]:
    by_slot: dict[int, int] = {}
    for code in range(group_size):
        slot = mix64_cpu(code) & table_mask
        if slot in by_slot:
            return by_slot[slot], code, slot
        by_slot[slot] = code
    raise AssertionError("could not manufacture a hash-slot collision")


def manufactured_cpu_opencl_self_audit(
    runtime: dict[str, Any], opencl_memory_cap_bytes: int
) -> dict:
    """Cross-check OpenCL sets with CPU brute force and force safety traps."""
    moduli = (3,) * 6 + (7,) * 3
    rng = np.random.default_rng(20_260_830)
    supports = [
        np.column_stack(
            [rng.integers(0, modulus, size=row_count) for modulus in moduli]
        ).astype(np.uint8)
        for row_count in (7, 6, 5, 4)
    ]
    expected_codes = cuda_join.brute_support_codes(supports, moduli)
    expected_set = {int(code) for code in expected_codes}
    present_code = int(expected_codes[len(expected_codes) // 3])
    absent_code = next(
        code for code in range(cuda_join.group_size_for(moduli)) if code not in expected_set
    )
    present_target = cuda_join.decode_code_cpu(present_code, moduli)
    absent_target = cuda_join.decode_code_cpu(absent_code, moduli)

    large_chunk, large_codes = opencl_exact_support_convolution(
        runtime,
        supports,
        present_target,
        moduli,
        state_cap=20_000,
        pair_chunk_cap=257,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    small_chunk, small_codes = opencl_exact_support_convolution(
        runtime,
        supports,
        absent_target,
        moduli,
        state_cap=20_000,
        pair_chunk_cap=17,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    require(large_codes is not None and small_codes is not None, "audit unexpectedly skipped")
    require(
        np.array_equal(large_codes, expected_codes)
        and np.array_equal(small_codes, expected_codes),
        "OpenCL support differs from exhaustive CPU support",
    )
    require(
        large_chunk["projected_target_present"] is True
        and small_chunk["rigorous_rejection"] is True,
        "OpenCL present/absent decision differs from CPU",
    )

    collision_table_size = 8
    first, second, slot = collision_pair(
        cuda_join.group_size_for(moduli), collision_table_size - 1
    )
    collision_rows = np.stack(
        [
            cuda_join.decode_code_cpu(first, moduli),
            cuda_join.decode_code_cpu(second, moduli),
        ]
    ).astype(np.uint8)
    collision_result, collision_codes = opencl_exact_support_convolution(
        runtime,
        [collision_rows],
        cuda_join.decode_code_cpu(second, moduli),
        moduli,
        state_cap=4,
        pair_chunk_cap=2,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    require(collision_codes is not None, "collision trap skipped")
    require(
        collision_result["factor_audits"][0]["hash_table_entries"] == collision_table_size
        and np.array_equal(collision_codes, np.asarray([first, second], dtype=np.uint64))
        and collision_result["projected_target_present"] is True,
        "uint64 collision-resolution trap failed",
    )

    trap = np.zeros((2, len(moduli)), dtype=np.uint8)
    trap[0, 6] = 1
    trap[1, 0] = 1
    trap_result, trap_codes = opencl_exact_support_convolution(
        runtime,
        [trap],
        np.zeros(len(moduli), dtype=np.uint8),
        moduli,
        state_cap=20,
        pair_chunk_cap=2,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    require(
        trap_codes is not None
        and trap_result["rigorous_rejection"] is True
        and 0 not in {int(code) for code in trap_codes},
        "same-row cross-prime trap failed",
    )

    cap_rows = np.stack(
        [cuda_join.decode_code_cpu(code, moduli) for code in range(11)]
    ).astype(np.uint8)
    cap_result, cap_codes = opencl_exact_support_convolution(
        runtime,
        [cap_rows],
        np.zeros(len(moduli), dtype=np.uint8),
        moduli,
        state_cap=10,
        pair_chunk_cap=13,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    require(
        cap_codes is None
        and cap_result["decision_status"] == "skipped_state_cap"
        and cap_result["distinct_states_lower_bound_at_skip"] == 11,
        "OpenCL state-cap trap failed",
    )

    memory_result, memory_codes = opencl_exact_support_convolution(
        runtime,
        [collision_rows],
        cuda_join.decode_code_cpu(first, moduli),
        moduli,
        state_cap=4,
        pair_chunk_cap=2,
        opencl_memory_cap_bytes=1,
        return_final_codes=True,
    )
    require(
        memory_codes is None
        and memory_result["decision_status"] == "skipped_opencl_memory_cap"
        and memory_result["skip_is_explicit_not_approximation"] is True,
        "OpenCL memory-cap trap failed",
    )

    supported_width_audits = []
    width_rng = np.random.default_rng(20_260_831)
    for width in SUPPORTED_K:
        width_moduli = (3,) * 6 + (7,) * width
        width_supports = [
            np.column_stack(
                [
                    width_rng.integers(0, modulus, size=row_count)
                    for modulus in width_moduli
                ]
            ).astype(np.uint8)
            for row_count in (4, 3, 2)
        ]
        width_expected = cuda_join.brute_support_codes(width_supports, width_moduli)
        width_target = cuda_join.decode_code_cpu(
            int(width_expected[len(width_expected) // 2]), width_moduli
        )
        width_result, width_codes = opencl_exact_support_convolution(
            runtime,
            width_supports,
            width_target,
            width_moduli,
            state_cap=1_000,
            pair_chunk_cap=7,
            opencl_memory_cap_bytes=opencl_memory_cap_bytes,
            return_final_codes=True,
        )
        require(
            width_codes is not None
            and np.array_equal(width_codes, width_expected)
            and width_result["projected_target_present"] is True,
            f"OpenCL k={width} width audit differs from CPU",
        )
        supported_width_audits.append(
            {
                "k": width,
                "projected_group": f"F3^6 x F7^{width}",
                "projected_group_size": cuda_join.group_size_for(width_moduli),
                "exact_state_count": len(width_expected),
                "support_sha256_uint64": cuda_join.array_sha256(width_expected),
                "OpenCL_elapsed_seconds": width_result["elapsed_seconds"],
                "OpenCL_peak_estimated_bytes": width_result["peak_estimated_bytes"],
                "CPU_and_OpenCL_complete_support_sets_equal_elementwise": True,
            }
        )

    semantic = {
        "expected_support_sha256_uint64": cuda_join.array_sha256(expected_codes),
        "large_chunk_support_sha256_uint64": large_chunk[
            "final_support_sha256_uint64"
        ],
        "small_chunk_support_sha256_uint64": small_chunk[
            "final_support_sha256_uint64"
        ],
        "forced_collision_codes": [first, second],
        "forced_collision_initial_slot": slot,
        "forced_collision_table_size": collision_table_size,
        "cap_decision_status": cap_result["decision_status"],
        "memory_decision_status": memory_result["decision_status"],
        "supported_width_audits": supported_width_audits,
    }
    return {
        "status": "passed",
        "projected_group": "F3^6 x F7^3",
        "cpu_brute_force_final_state_count": len(expected_codes),
        "cpu_brute_force_support_sha256_uint64": cuda_join.array_sha256(expected_codes),
        "opencl_large_chunk_support_sha256_uint64": large_chunk[
            "final_support_sha256_uint64"
        ],
        "opencl_small_chunk_support_sha256_uint64": small_chunk[
            "final_support_sha256_uint64"
        ],
        "cpu_vs_opencl_exact_set_equality": True,
        "present_target_detected": True,
        "absent_target_rigorously_rejected": True,
        "pair_chunk_invariance_audited": True,
        "forced_distinct_uint64_initial_slot_collision": True,
        "collision_keys": [first, second],
        "collision_initial_slot": slot,
        "collision_resolved_without_loss_or_false_membership": True,
        "same_catalog_row_cross_prime_trap_rejected": True,
        "state_cap_produces_explicit_skip": True,
        "memory_cap_produces_explicit_skip": True,
        "supported_production_widths_k5_k6_crosschecked_against_CPU": True,
        "supported_width_audits": supported_width_audits,
        "probabilistic_membership_used": False,
        "semantic_audit_sha256": json_sha256(semantic),
    }


def manufactured_cross_engine_audit(
    runtime: dict[str, Any], opencl_memory_cap_bytes: int
) -> dict:
    cpu_audit = cpu_join.manufactured_case_join_audit()
    opencl_audit = manufactured_cpu_opencl_self_audit(
        runtime, opencl_memory_cap_bytes
    )
    require(cpu_audit["passed"], "CPU manufactured case-join audit failed")
    require(
        cpu_audit["same_row_cross_prime_false_positive_trap_rejected"]
        and opencl_audit["same_catalog_row_cross_prime_trap_rejected"],
        "CPU/OpenCL same-row trap failed",
    )
    require(
        cpu_audit["state_cap_is_explicit_skip_and_partial_support_is_discarded"]
        and opencl_audit["state_cap_produces_explicit_skip"],
        "CPU/OpenCL cap semantics failed",
    )
    synthetic_targets = [{"case_key": f"synthetic_{index}"} for index in range(11)]
    shard_rows = []
    covered_indices = []
    for shard_index in range(3):
        selected, audit = case_gpu.shard_cases(synthetic_targets, shard_index, 3)
        indices = [index for index, _row in selected]
        require(indices == list(range(shard_index, 11, 3)), "shard rule changed")
        require(audit["all_shards_form_a_disjoint_cover_of_the_51_cases"], "shard audit failed")
        covered_indices.extend(indices)
        shard_rows.append(indices)
    require(sorted(covered_indices) == list(range(11)), "synthetic shards do not cover")
    require(len(covered_indices) == len(set(covered_indices)), "synthetic shards overlap")
    return {
        "status": "passed",
        "CPU_manufactured_case_join": cpu_audit,
        "OpenCL_manufactured_exact_set": opencl_audit,
        "same_row_cross_prime_traps_passed_on_CPU_and_OpenCL": True,
        "state_and_memory_caps_are_explicit_nondecision_skips": True,
        "OpenCL_pair_chunk_invariance_and_collision_resolution_audited": True,
        "deterministic_case_sharding_audit": {
            "rule": "target_index_mod_shard_count_equals_shard_index",
            "shard_count": 3,
            "selected_indices_by_shard": shard_rows,
            "disjoint_complete_cover": True,
        },
    }


def validate_opencl_decision(decision: dict) -> None:
    status = str(decision["decision_status"])
    if status.startswith("skipped_"):
        require(not decision["completed_exact_convolution"], "skip labeled complete")
        require(not decision["rigorous_rejection"], "skip became rejection")
        require(decision["projected_target_present"] is None, "skip guessed target")
        require(decision["skip_is_explicit_not_approximation"], "skip semantics weakened")
    else:
        require(decision["completed_exact_convolution"], "non-skip did not complete")
        require(
            decision["rigorous_rejection"] is (not decision["projected_target_present"]),
            "completed OpenCL decision is inconsistent",
        )
        require(
            decision["opencl_hash_set_is_exact_not_probabilistic"],
            "OpenCL set semantics weakened",
        )


def build_direction_supports(
    context: dict[str, Any],
    projection: tuple[int, ...],
    pair_chunk_cap: int,
    direction_workers: int,
) -> tuple[dict[int, dict[int, np.ndarray]], dict]:
    """Use the current companion's audited serial support builder unchanged."""
    require(
        direction_workers == 1,
        "the current CUDA companion exposes only the serial support builder",
    )
    return case_gpu.build_direction_supports(context, projection, pair_chunk_cap)


def real_case_cpu_opencl_crosscheck(
    runtime: dict[str, Any],
    opencl_memory_cap_bytes: int,
    context: dict[str, Any],
    pair_chunk_cap: int,
    direction_workers: int,
) -> dict:
    """Compare complete CPU/OpenCL sets for one real F3^6 x F7^2 case."""
    projection = case_gpu.DEFAULT_CROSSCHECK_PROJECTION
    table, direction_audit = build_direction_supports(
        context, projection, pair_chunk_cap, direction_workers
    )
    target_row = context["targets"][0]
    support_rows, target, moduli, ordered, case_audit = case_gpu.exact_case_inputs(
        context, target_row, projection, table
    )
    codec = cpu_join.semigroup.MixedRadixCodec(moduli)
    cpu_codes, cpu_convolution = cpu_join.convolve_support_sequence(
        ordered,
        codec,
        state_cap=codec.group_size,
        pair_chunk_cap=min(pair_chunk_cap, codec.group_size),
    )
    require(cpu_codes is not None and cpu_convolution["completed"], "CPU crosscheck capped")
    decision, opencl_codes = opencl_exact_support_convolution(
        runtime,
        support_rows,
        target,
        moduli,
        state_cap=codec.group_size,
        pair_chunk_cap=pair_chunk_cap,
        opencl_memory_cap_bytes=opencl_memory_cap_bytes,
        return_final_codes=True,
    )
    validate_opencl_decision(decision)
    require(opencl_codes is not None, "real OpenCL crosscheck unexpectedly skipped")
    require(np.array_equal(cpu_codes, opencl_codes), "real CPU/OpenCL sets differ")
    target_code = int(codec.encode(target[None, :])[0])
    target_index = int(np.searchsorted(cpu_codes, np.uint64(target_code)))
    cpu_present = bool(
        target_index < len(cpu_codes) and int(cpu_codes[target_index]) == target_code
    )
    require(cpu_present == decision["projected_target_present"], "target decisions differ")
    return {
        "status": "passed",
        "projection_mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "case_input_audit": case_audit,
        "direction_support_audit": direction_audit,
        "CPU_convolution_audit": cpu_convolution,
        "OpenCL_convolution_audit": decision,
        "complete_support_state_count": len(cpu_codes),
        "complete_support_sha256_uint64": cpu_join.array_sha256(
            cpu_codes.astype("<u8", copy=False)
        ),
        "target_present": cpu_present,
        "CPU_and_OpenCL_complete_support_sets_equal_elementwise": True,
        "CPU_and_OpenCL_target_decisions_equal": True,
    }


def expected_engine_memory(k: int, state_cap: int, pair_chunk_cap: int) -> dict:
    moduli = (3,) * 6 + (7,) * k
    group_size = cuda_join.group_size_for(moduli)
    effective_state_cap = min(state_cap, group_size)
    maximum_insertions = min(group_size, effective_state_cap + pair_chunk_cap + 1)
    table_size = cuda_join.next_power_of_two(max(8, 2 * maximum_insertions))
    estimated_bytes = (
        table_size * 8
        + maximum_insertions * 8
        + (effective_state_cap + case_gpu.EXACT_DIRECTION_SUPPORT_CAP) * 8
        + len(moduli) * 4
        + 32
    )
    return {
        "projected_group_size": group_size,
        "effective_state_cap": effective_state_cap,
        "maximum_direction_support_rows": case_gpu.EXACT_DIRECTION_SUPPORT_CAP,
        "hash_table_entries_upper_bound": table_size,
        "conservative_opencl_device_peak_bytes_upper_bound": estimated_bytes,
        "conservative_opencl_device_peak_mib_upper_bound": estimated_bytes / (1 << 20),
        "host_final_code_array_bytes_if_full_group": group_size * 8,
        "uint32_distinct_counter_is_exact_for_group": group_size < 2**32,
    }


def run_production(
    runtime: dict[str, Any],
    context: dict[str, Any],
    context_audit: dict,
    manufactured_audit: dict,
    real_crosscheck: dict | None,
    projections: tuple[tuple[int, ...], ...],
    selected_cases: list[tuple[int, dict]],
    shard_audit: dict,
    state_cap: int,
    pair_chunk_cap: int,
    direction_workers: int,
    opencl_memory_cap_bytes: int,
    output_path: Path,
) -> dict:
    started = time.time()
    projection_tables = {}
    support_audits = []
    for projection in projections:
        table, audit = build_direction_supports(
            context, projection, pair_chunk_cap, direction_workers
        )
        projection_tables[projection] = table
        support_audits.append(audit)

    projection_totals: Counter[str] = Counter()
    case_results = []
    for target_index, target_row in selected_cases:
        projection_results = []
        for projection in projections:
            support_rows, target, moduli, _ordered, case_input_audit = (
                case_gpu.exact_case_inputs(
                    context, target_row, projection, projection_tables[projection]
                )
            )
            decision, _codes = opencl_exact_support_convolution(
                runtime,
                support_rows,
                target,
                moduli,
                state_cap=state_cap,
                pair_chunk_cap=pair_chunk_cap,
                opencl_memory_cap_bytes=opencl_memory_cap_bytes,
            )
            validate_opencl_decision(decision)
            projection_totals[str(decision["decision_status"])] += 1
            projection_results.append(
                {
                    "mod7_coordinates": list(projection),
                    "retained_all_six_mod3_coordinates": True,
                    "case_input_audit": case_input_audit,
                    **decision,
                }
            )

        rejected = any(row["rigorous_rejection"] for row in projection_results)
        skipped = not rejected and any(
            str(row["decision_status"]).startswith("skipped_")
            for row in projection_results
        )
        necessary = not rejected and not skipped
        require(sum((rejected, skipped, necessary)) == 1, "case decision is ambiguous")
        case_key = str(target_row["case_key"])
        row = {
            "target_index_in_audited_51_case_order": target_index,
            "case_key": case_key,
            "catalog_pattern": context["current_by_key"][case_key]["catalog_pattern"],
            "prior_global_join_decision": context["current_by_key"][case_key][
                "decision_status"
            ],
            "projection_results": projection_results,
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
            "decision_status": (
                "rigorous_exact_OpenCL_semigroup_projection_rejection"
                if rejected
                else "explicit_cap_skip_without_negative_decision"
                if skipped
                else "necessary_only_survivor_of_all_completed_OpenCL_projections"
            ),
        }
        row["decision_certificate_sha256"] = json_sha256(row)
        case_results.append(row)

    counts = {
        "selected": len(case_results),
        "rejected": sum(row["rigorously_rejected"] for row in case_results),
        "surviving": sum(row["necessary_only_survivor"] for row in case_results),
        "skipped": sum(row["skipped"] for row in case_results),
    }
    require(sum(counts[key] for key in ("rejected", "surviving", "skipped")) == counts["selected"], "case census failed")
    k = len(projections[0])
    script_path = Path(__file__).resolve()
    return {
        "experiment": EXPERIMENT,
        "status": "complete_sharded_exact_OpenCL_semigroup_case_join",
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "51 orbit0/branch-A grade-three-only representatives",
        "opencl_invoked": True,
        "solver_invoked": False,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "CPU_semigroup_case_join_path": str(Path(cpu_join.__file__).resolve()),
            "CPU_semigroup_case_join_sha256": file_sha256(Path(cpu_join.__file__)),
            "CUDA_semigroup_companion_reused_for_construction_path": str(
                Path(case_gpu.__file__).resolve()
            ),
            "CUDA_semigroup_companion_sha256": file_sha256(Path(case_gpu.__file__)),
            "CUDA_torsion_engine_reused_for_CPU_codec_path": str(
                Path(cuda_join.__file__).resolve()
            ),
            "CUDA_torsion_engine_sha256": file_sha256(Path(cuda_join.__file__)),
        },
        "configuration": {
            "k": k,
            "explicit_mod7_projection_subsets": [list(row) for row in projections],
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "opencl_memory_cap_bytes": opencl_memory_cap_bytes,
            "exact_direction_support_cap": case_gpu.EXACT_DIRECTION_SUPPORT_CAP,
            "direction_support_workers": direction_workers,
        },
        "expected_memory": expected_engine_memory(k, state_cap, pair_chunk_cap),
        "opencl_engine_audit": runtime["audit"],
        "CPU_construction_and_input_audits": context_audit,
        "manufactured_CPU_OpenCL_audit": manufactured_audit,
        "small_real_case_CPU_OpenCL_crosscheck": real_crosscheck,
        "small_real_case_crosscheck_skipped_by_explicit_flag": real_crosscheck is None,
        "case_sharding_audit": shard_audit,
        "projected_direction_support_audits": support_audits,
        "result_counts": counts,
        "projection_decision_census": dict(sorted(projection_totals.items())),
        "case_results_sha256": cpu_join.old_join.canonical_case_digest(case_results),
        "case_results": case_results,
        "logical_semantics": {
            "all_51_cases_were_reconstructed_before_case_sharding": True,
            "all_six_mod3_coordinates_are_retained": True,
            "mod7_coordinates_are_explicitly_selected_from_the_derived_21D_quotient": True,
            "same_Hilbert_or_catalog_row_supplies_mod3_and_mod7_before_deduplication": True,
            "grade_zero_through_three_direction_supports_are_exact_and_direct_catalog_calibrated": True,
            "all_eight_direction_supports_are_convolved_in_every_completed_projection": True,
            "OpenCL_uint64_hash_collisions_are_resolved_by_exact_key_comparison": True,
            "cl_khr_int64_base_atomics_is_required_and_audited": True,
            "exact_full_group_saturation_shortcut_only_proves_target_presence": True,
            "missing_target_after_completed_convolution_is_a_rigorous_rejection": True,
            "target_presence_is_necessary_only": True,
            "all_state_memory_allocation_or_hash_capacity_caps_are_explicit_skips": True,
            "partial_support_after_any_cap_is_never_used": True,
            "probabilistic_membership_used": False,
            "binary_edge_feasibility_claimed": False,
            "positive_z7_closure_claimed": False,
        },
        "positive_z7_excluded": False,
        "full_theorem_claimed": False,
        "output_path": str(output_path.resolve()),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform-index", type=int, default=0)
    parser.add_argument("--device-index", type=int, default=0)
    parser.add_argument("--device-name-contains")
    parser.add_argument("--parent-input", type=Path, default=cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument("--k", type=int, choices=SUPPORTED_K, default=5)
    parser.add_argument(
        "--mod7-projections",
        help="semicolon-separated width-k subsets; default is the width-k prefix",
    )
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument(
        "--state-cap",
        type=int,
        default=0,
        help="global support cap; 0 means the complete projected-group size",
    )
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument("--direction-support-workers", type=int, default=1)
    parser.add_argument(
        "--opencl-memory-cap-mib", type=int, default=DEFAULT_OPENCL_MEMORY_CAP_MIB
    )
    parser.add_argument("--skip-small-real-crosscheck", action="store_true")
    parser.add_argument("--manufactured-only", action="store_true")
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    require(
        1 <= args.direction_support_workers <= MAX_DIRECTION_SUPPORT_WORKERS,
        "--direction-support-workers lies outside the audited range",
    )
    require(args.opencl_memory_cap_mib > 0, "--opencl-memory-cap-mib must be positive")
    require(
        not (args.manufactured_only and args.self_audit_only),
        "manufactured-only and self-audit-only are mutually exclusive",
    )

    projections = case_gpu.parse_projection_subsets(args.mod7_projections, args.k)
    group_size = cuda_join.group_size_for((3,) * 6 + (7,) * args.k)
    require(args.state_cap >= 0, "--state-cap cannot be negative")
    state_cap = group_size if args.state_cap == 0 else min(args.state_cap, group_size)
    opencl_memory_cap_bytes = args.opencl_memory_cap_mib * (1 << 20)
    runtime = load_opencl(
        args.platform_index, args.device_index, args.device_name_contains
    )
    manufactured = manufactured_cross_engine_audit(
        runtime, opencl_memory_cap_bytes
    )

    context = None
    context_audit = None
    real_crosscheck = None
    if not args.manufactured_only:
        context, context_audit = case_gpu.build_cpu_context(
            args.parent_input, args.current_join, args.hilbert_basis
        )
        real_crosscheck = (
            None
            if args.skip_small_real_crosscheck and not args.self_audit_only
            else real_case_cpu_opencl_crosscheck(
                runtime,
                opencl_memory_cap_bytes,
                context,
                args.pair_chunk_cap,
                args.direction_support_workers,
            )
        )

    if args.manufactured_only:
        result = {
            "experiment": f"{EXPERIMENT}_manufactured_self_audit",
            "status": "manufactured_CPU_OpenCL_self_audit_passed",
            "opencl_engine_audit": runtime["audit"],
            "manufactured_CPU_OpenCL_audit": manufactured,
            "small_real_case_crosscheck_run": False,
            "production_cases_processed": 0,
            "output_path": str(args.output.resolve()),
        }
    elif args.self_audit_only:
        require(real_crosscheck is not None, "self-audit omitted the real crosscheck")
        result = {
            "experiment": f"{EXPERIMENT}_self_audit",
            "status": "manufactured_and_real_case_CPU_OpenCL_self_audits_passed",
            "opencl_engine_audit": runtime["audit"],
            "CPU_construction_and_input_audits": context_audit,
            "manufactured_CPU_OpenCL_audit": manufactured,
            "small_real_case_CPU_OpenCL_crosscheck": real_crosscheck,
            "production_cases_processed": 0,
            "output_path": str(args.output.resolve()),
        }
    else:
        require(context is not None and context_audit is not None, "CPU context missing")
        selected_cases, shard_audit = case_gpu.shard_cases(
            context["targets"], args.shard_index, args.shard_count
        )
        result = run_production(
            runtime,
            context,
            context_audit,
            manufactured,
            real_crosscheck,
            projections,
            selected_cases,
            shard_audit,
            state_cap,
            args.pair_chunk_cap,
            args.direction_support_workers,
            opencl_memory_cap_bytes,
            args.output,
        )

    cuda_join.atomic_write(args.output, result)
    summary = {
        "status": result["status"],
        "output": str(args.output),
        "opencl_device": runtime["audit"]["device_name"],
        "manufactured_CPU_OpenCL_audit": manufactured["status"],
        "small_real_case_CPU_OpenCL_crosscheck": (
            None if real_crosscheck is None else real_crosscheck["status"]
        ),
        "processed_cases": result.get("result_counts", {}).get("selected", 0),
        "rigorously_rejected_cases": result.get("result_counts", {}).get("rejected", 0),
        "necessary_only_survivors": result.get("result_counts", {}).get("surviving", 0),
        "skipped_cases": result.get("result_counts", {}).get("skipped", 0),
        "expected_memory": expected_engine_memory(args.k, state_cap, args.pair_chunk_cap),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
