#!/usr/bin/env python3
"""Exact CUDA torsion-support projections for the four z=7 H0_S0_M7 cases.

This is the GPU companion to
``p7_infinity7_positive_z7_torsion_support_projection.py``.  It reuses that
script's audited construction of the pointed orbit-0/A system, its exact
114-dimensional rational dependency space, and its derived characteristic-3
and characteristic-7 quotient bases.  In particular, the effective
dimensions are derived at run time; they are not hard-coded here.

All six mod-3 quotient coordinates and a bounded selection of at least three
mod-7 quotient coordinates are encoded together in one mixed-radix uint64
key.  Thus a catalog row always contributes its mod-3 and mod-7 signatures
with the same row index.  A CUDA open-addressing hash set computes every
Minkowski support union exactly.  Hash collisions are resolved by key
comparison and probing; there are no probabilistic rejections.

``--state-cap`` and ``--gpu-memory-cap-mib`` are proof-safety boundaries.
Crossing either boundary produces an explicit skip and never a rejection.
The manufactured self-audit compares CUDA results with complete CPU brute
force, checks chunk invariance and the shared-row cross-prime trap, and forces
the state-cap path.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MIN_PROJECTED_MOD7_COORDINATES = 3
MAX_PROJECTIONS = 16
DEFAULT_STATE_CAP = 2_000_000
DEFAULT_PAIR_CHUNK_CAP = 1_000_000
DEFAULT_GPU_MEMORY_CAP_MIB = 2_048
THREADS_PER_BLOCK = 256
EMPTY_KEY = np.uint64(np.iinfo(np.uint64).max)


CUDA_SOURCE = r"""
extern "C" {

__device__ __forceinline__ unsigned long long mix64(unsigned long long x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

__device__ __forceinline__ unsigned long long mixed_sum_key(
    unsigned long long left,
    unsigned long long right,
    const unsigned int* moduli,
    int width
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    for (int coordinate = 0; coordinate < width; ++coordinate) {
        const unsigned int modulus = moduli[coordinate];
        const unsigned int left_digit = (unsigned int)(left % modulus);
        const unsigned int right_digit = (unsigned int)(right % modulus);
        const unsigned int digit = (left_digit + right_digit) % modulus;
        key += place * (unsigned long long)digit;
        place *= (unsigned long long)modulus;
        left /= (unsigned long long)modulus;
        right /= (unsigned long long)modulus;
    }
    return key;
}

__global__ void insert_mixed_sums(
    const unsigned long long* states,
    unsigned long long state_count,
    const unsigned long long* support,
    unsigned long long support_count,
    const unsigned int* moduli,
    int width,
    unsigned long long pair_offset,
    unsigned long long pair_count,
    unsigned long long* table,
    unsigned long long table_mask,
    unsigned long long state_cap,
    unsigned long long* distinct_count,
    unsigned int* overflow,
    unsigned int* probe_failure
) {
    const unsigned long long local =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= pair_count || *overflow != 0U || *probe_failure != 0U) return;

    const unsigned long long pair_index = pair_offset + local;
    const unsigned long long left_index = pair_index / support_count;
    const unsigned long long right_index = pair_index % support_count;
    if (left_index >= state_count) return;

    const unsigned long long key = mixed_sum_key(
        states[left_index], support[right_index], moduli, width
    );
    unsigned long long slot = mix64(key) & table_mask;
    unsigned long long probes = 0ULL;
    while (probes <= table_mask) {
        const unsigned long long previous = atomicCAS(
            table + slot, 0xffffffffffffffffULL, key
        );
        if (previous == 0xffffffffffffffffULL) {
            const unsigned long long old_count = atomicAdd(distinct_count, 1ULL);
            if (old_count >= state_cap) atomicExch(overflow, 1U);
            return;
        }
        if (previous == key) return;
        slot = (slot + 1ULL) & table_mask;
        ++probes;
    }
    atomicExch(probe_failure, 1U);
}

}
"""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_cupy(device: int) -> tuple[Any, Any, dict]:
    """Import CuPy lazily so py_compile and --help need no CUDA installation."""
    try:
        import cupy as cp  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - depends on host packaging
        raise RuntimeError("CuPy is required for the CUDA support engine") from exc

    cp.cuda.Device(device).use()
    properties = cp.cuda.runtime.getDeviceProperties(device)
    raw_name = properties.get("name", b"unknown")
    name = raw_name.decode("utf-8", "replace") if isinstance(raw_name, bytes) else str(raw_name)
    free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()
    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++14",))
    kernel = module.get_function("insert_mixed_sums")
    audit = {
        "cupy_version": cp.__version__,
        "cuda_runtime_version": int(cp.cuda.runtime.runtimeGetVersion()),
        "device_index": device,
        "device_name": name,
        "compute_capability": [
            int(properties.get("major", -1)),
            int(properties.get("minor", -1)),
        ],
        "total_global_memory_bytes": int(total_bytes),
        "free_global_memory_bytes_at_start": int(free_bytes),
        "cuda_kernel_sha256": hashlib.sha256(CUDA_SOURCE.encode("utf-8")).hexdigest(),
        "hash_collisions_resolved_by_exact_key_comparison_and_linear_probing": True,
        "probabilistic_filter_used": False,
    }
    return cp, kernel, audit


def load_projection_module() -> Any:
    """Load the audited CPU derivation only for a real four-case run."""
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.path.insert(0, str(ROOT / "src"))
    import p7_infinity7_positive_z7_torsion_support_projection as source

    return source


def group_size_for(moduli: tuple[int, ...]) -> int:
    require(bool(moduli), "empty projected group")
    require(all(int(modulus) >= 2 for modulus in moduli), "invalid projected modulus")
    group_size = math.prod(int(modulus) for modulus in moduli)
    require(group_size < int(EMPTY_KEY), "mixed-radix projected group does not fit uint64")
    return group_size


def encode_rows_cpu(rows: np.ndarray, moduli: tuple[int, ...]) -> np.ndarray:
    source = np.ascontiguousarray(rows, dtype=np.uint8)
    require(source.ndim == 2 and source.shape[1] == len(moduli), "encoding width changed")
    codes = np.zeros(len(source), dtype=np.uint64)
    place = 1
    for coordinate, modulus in enumerate(moduli):
        require(
            not np.any(source[:, coordinate] >= modulus),
            f"coordinate {coordinate} escaped modulus {modulus}",
        )
        codes += source[:, coordinate].astype(np.uint64) * np.uint64(place)
        place *= int(modulus)
    require(place == group_size_for(moduli), "mixed-radix place product changed")
    return np.ascontiguousarray(codes)


def decode_code_cpu(code: int, moduli: tuple[int, ...]) -> np.ndarray:
    require(0 <= code < group_size_for(moduli), "code lies outside projected group")
    value = int(code)
    row = np.empty(len(moduli), dtype=np.uint8)
    for coordinate, modulus in enumerate(moduli):
        row[coordinate] = value % modulus
        value //= modulus
    require(value == 0, "mixed-radix decode left a quotient")
    return row


def unique_support_codes(rows: np.ndarray, moduli: tuple[int, ...]) -> np.ndarray:
    codes = np.unique(encode_rows_cpu(rows, moduli))
    return np.ascontiguousarray(codes, dtype=np.uint64)


def brute_support_codes(
    supports: list[np.ndarray], moduli: tuple[int, ...]
) -> np.ndarray:
    """Complete CPU brute force for manufactured, deliberately tiny supports."""
    width = len(moduli)
    modulus_row = np.asarray(moduli, dtype=np.uint16)
    states = {tuple(0 for _ in moduli)}
    for support in supports:
        rows = np.ascontiguousarray(support, dtype=np.uint8)
        require(rows.ndim == 2 and rows.shape[1] == width, "brute support width changed")
        states = {
            tuple(
                int((left[coordinate] + int(right[coordinate])) % moduli[coordinate])
                for coordinate in range(width)
            )
            for left in states
            for right in rows
        }
        require(all(np.all(np.asarray(row) < modulus_row) for row in states), "bad brute state")
    rows = np.asarray(sorted(states), dtype=np.uint8)
    return np.sort(unique_support_codes(rows, moduli))


def next_power_of_two(value: int) -> int:
    require(value > 0, "power-of-two request must be positive")
    return 1 << (value - 1).bit_length()


def convolution_skip(
    status: str,
    normalized_sizes: list[int],
    group_size: int,
    state_cap: int,
    completed_state_sizes: list[int],
    completed_factors: int,
    pair_candidates_launched: int,
    **extra: object,
) -> dict:
    require(status.startswith("skipped_"), "non-skip status passed to skip constructor")
    return {
        "completed_exact_convolution": False,
        "decision_status": status,
        "rigorous_rejection": False,
        "projected_target_present": None,
        "necessary_only_projected_survivor": None,
        "completed_catalog_factors": completed_factors,
        "first_skipped_factor_index": completed_factors,
        "state_cap": state_cap,
        "finite_group_size": group_size,
        "exact_support_sizes": normalized_sizes,
        "completed_state_sizes": completed_state_sizes,
        "pair_candidates_launched": pair_candidates_launched,
        "skip_is_explicit_not_approximation": True,
        **extra,
    }


def gpu_exact_support_convolution(
    cp: Any,
    kernel: Any,
    supports: list[np.ndarray],
    target: np.ndarray,
    moduli: tuple[int, ...],
    state_cap: int,
    pair_chunk_cap: int,
    gpu_memory_cap_bytes: int,
    return_final_codes: bool = False,
) -> tuple[dict, np.ndarray | None]:
    """Compute an exact capped support convolution using a CUDA hash set."""
    started = time.perf_counter()
    require(state_cap > 0, "state cap must be positive")
    require(pair_chunk_cap > 0, "pair chunk cap must be positive")
    require(gpu_memory_cap_bytes > 0, "GPU memory cap must be positive")
    require(len(target) == len(moduli), "projected target width changed")
    group_size = group_size_for(moduli)
    normalized_codes = [unique_support_codes(rows, moduli) for rows in supports]
    normalized_sizes = [len(codes) for codes in normalized_codes]
    require(all(size > 0 for size in normalized_sizes), "empty exact catalog support")
    target_code = int(encode_rows_cpu(np.asarray(target, dtype=np.uint8)[None, :], moduli)[0])

    moduli_device = cp.asarray(moduli, dtype=cp.uint32)
    states = cp.asarray([0], dtype=cp.uint64)
    completed_state_sizes: list[int] = []
    factor_audits = []
    pair_candidates_launched = 0
    peak_estimated_bytes = 0

    for factor_index, support_host in enumerate(normalized_codes):
        factor_started = time.perf_counter()
        support_device = cp.asarray(support_host, dtype=cp.uint64)
        state_count = int(states.size)
        support_count = int(support_device.size)
        total_pairs = state_count * support_count
        launch_chunk = min(pair_chunk_cap, total_pairs)

        # Before the overflow flag becomes visible, at most one launched chunk
        # can add keys beyond the cap.  Sizing for that entire chunk keeps the
        # open-addressing table below 50% load even under the worst race.
        maximum_insertions = min(group_size, state_cap + launch_chunk + 1)
        table_size = next_power_of_two(max(8, 2 * maximum_insertions))
        require(table_size - 1 <= np.iinfo(np.uint64).max, "hash table mask overflow")
        estimated_peak_bytes = (
            table_size * 24
            + (state_count + support_count) * 8
            + maximum_insertions * 16
            + len(moduli) * 4
        )
        peak_estimated_bytes = max(peak_estimated_bytes, estimated_peak_bytes)
        free_bytes, _total_bytes = cp.cuda.runtime.memGetInfo()
        effective_memory_limit = min(gpu_memory_cap_bytes, int(0.85 * free_bytes))
        if estimated_peak_bytes > effective_memory_limit:
            result = convolution_skip(
                "skipped_gpu_memory_cap",
                normalized_sizes,
                group_size,
                state_cap,
                completed_state_sizes,
                factor_index,
                pair_candidates_launched,
                required_estimated_peak_bytes=estimated_peak_bytes,
                configured_gpu_memory_cap_bytes=gpu_memory_cap_bytes,
                free_gpu_memory_bytes_before_factor=int(free_bytes),
                effective_gpu_memory_limit_bytes=effective_memory_limit,
                skipped_factor_pair_count=total_pairs,
            )
            result["elapsed_seconds"] = time.perf_counter() - started
            return result, None

        try:
            table = cp.full(table_size, EMPTY_KEY, dtype=cp.uint64)
            distinct_count = cp.zeros(1, dtype=cp.uint64)
            overflow = cp.zeros(1, dtype=cp.uint32)
            probe_failure = cp.zeros(1, dtype=cp.uint32)
        except cp.cuda.memory.OutOfMemoryError:
            result = convolution_skip(
                "skipped_gpu_out_of_memory",
                normalized_sizes,
                group_size,
                state_cap,
                completed_state_sizes,
                factor_index,
                pair_candidates_launched,
                required_estimated_peak_bytes=estimated_peak_bytes,
                configured_gpu_memory_cap_bytes=gpu_memory_cap_bytes,
                skipped_factor_pair_count=total_pairs,
            )
            result["elapsed_seconds"] = time.perf_counter() - started
            return result, None

        factor_launched = 0
        for offset in range(0, total_pairs, pair_chunk_cap):
            count = min(pair_chunk_cap, total_pairs - offset)
            blocks = (count + THREADS_PER_BLOCK - 1) // THREADS_PER_BLOCK
            kernel(
                (blocks,),
                (THREADS_PER_BLOCK,),
                (
                    states,
                    np.uint64(state_count),
                    support_device,
                    np.uint64(support_count),
                    moduli_device,
                    np.int32(len(moduli)),
                    np.uint64(offset),
                    np.uint64(count),
                    table,
                    np.uint64(table_size - 1),
                    np.uint64(state_cap),
                    distinct_count,
                    overflow,
                    probe_failure,
                ),
            )
            factor_launched += count
            pair_candidates_launched += count
            if int(probe_failure.get()[0]):
                result = convolution_skip(
                    "skipped_gpu_hash_table_capacity",
                    normalized_sizes,
                    group_size,
                    state_cap,
                    completed_state_sizes,
                    factor_index,
                    pair_candidates_launched,
                    hash_table_entries=table_size,
                    successful_distinct_insertions=int(distinct_count.get()[0]),
                    skipped_factor_pair_count=total_pairs,
                    skipped_factor_pairs_launched=factor_launched,
                )
                result["elapsed_seconds"] = time.perf_counter() - started
                return result, None
            if int(overflow.get()[0]):
                inserted = int(distinct_count.get()[0])
                require(inserted >= state_cap + 1, "overflow flag lacks cap-plus-one witness")
                result = convolution_skip(
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

        cp.cuda.runtime.deviceSynchronize()
        inserted = int(distinct_count.get()[0])
        require(inserted <= state_cap, "completed factor exceeded state cap")
        occupied = table != EMPTY_KEY
        occupied_count = int(cp.count_nonzero(occupied).get())
        require(occupied_count == inserted, "CUDA hash occupancy/count mismatch")
        next_states = cp.sort(table[occupied])
        require(
            bool(cp.all(next_states[1:] > next_states[:-1]).get()) if inserted > 1 else True,
            "CUDA hash compaction retained a duplicate or lost ordering",
        )
        require(
            bool(cp.all(next_states < np.uint64(group_size)).get()),
            "CUDA hash produced a code outside the projected group",
        )
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
        states = next_states
        del table, occupied, next_states, support_device, distinct_count, overflow, probe_failure

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

    final_codes = np.ascontiguousarray(cp.asnumpy(states), dtype=np.uint64)
    require(
        len(final_codes) == len(np.unique(final_codes)),
        "final CUDA support is not a set",
    )
    index = int(np.searchsorted(final_codes, np.uint64(target_code)))
    present = bool(index < len(final_codes) and int(final_codes[index]) == target_code)
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
        "final_support_sha256_uint64": array_sha256(final_codes),
        "target_code": target_code,
        "pair_candidates_launched": pair_candidates_launched,
        "factor_audits": factor_audits,
        "peak_estimated_bytes": peak_estimated_bytes,
        "full_group_saturation_shortcut_is_exact": bool(len(final_codes) == group_size),
        "cuda_hash_set_is_exact_not_probabilistic": True,
        "elapsed_seconds": time.perf_counter() - started,
    }
    return result, final_codes if return_final_codes else None


def manufactured_cpu_gpu_self_audit(
    cp: Any,
    kernel: Any,
    gpu_memory_cap_bytes: int,
) -> dict:
    """Cross-check the CUDA set engine against exhaustive CPU enumeration."""
    moduli = (3,) * 6 + (7,) * 3
    rng = np.random.default_rng(20_260_830)
    supports = []
    for row_count in (7, 6, 5, 4):
        rows = np.column_stack(
            [rng.integers(0, modulus, size=row_count) for modulus in moduli]
        ).astype(np.uint8)
        supports.append(rows)
    expected_codes = brute_support_codes(supports, moduli)
    expected_set = {int(code) for code in expected_codes}
    present_code = int(expected_codes[len(expected_codes) // 3])
    absent_code = next(code for code in range(group_size_for(moduli)) if code not in expected_set)
    present_target = decode_code_cpu(present_code, moduli)
    absent_target = decode_code_cpu(absent_code, moduli)

    large_chunk, large_codes = gpu_exact_support_convolution(
        cp,
        kernel,
        supports,
        present_target,
        moduli,
        state_cap=20_000,
        pair_chunk_cap=257,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    small_chunk, small_codes = gpu_exact_support_convolution(
        cp,
        kernel,
        supports,
        absent_target,
        moduli,
        state_cap=20_000,
        pair_chunk_cap=17,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    require(large_codes is not None and small_codes is not None, "self-audit unexpectedly skipped")
    require(
        np.array_equal(large_codes, expected_codes)
        and np.array_equal(small_codes, expected_codes),
        "CUDA support differs from complete CPU brute force",
    )
    require(
        large_chunk["projected_target_present"] is True
        and small_chunk["rigorous_rejection"] is True,
        "CUDA present/absent decision disagrees with CPU brute force",
    )

    # No row is zero jointly, although one row is zero on the F3 side and the
    # other is zero on the F7 side.  Separate-prime supports would accept zero.
    trap = np.zeros((2, len(moduli)), dtype=np.uint8)
    trap[0, 6] = 1
    trap[1, 0] = 1
    trap_result, trap_codes = gpu_exact_support_convolution(
        cp,
        kernel,
        [trap],
        np.zeros(len(moduli), dtype=np.uint8),
        moduli,
        state_cap=10,
        pair_chunk_cap=2,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    require(
        trap_codes is not None
        and trap_result["rigorous_rejection"] is True
        and np.all(trap[0, :6] == 0)
        and np.all(trap[1, 6:] == 0),
        "same-row cross-prime CUDA trap failed",
    )

    cap_rows = np.stack([decode_code_cpu(code, moduli) for code in range(32)])
    cap_result, cap_codes = gpu_exact_support_convolution(
        cp,
        kernel,
        [cap_rows],
        np.zeros(len(moduli), dtype=np.uint8),
        moduli,
        state_cap=10,
        pair_chunk_cap=13,
        gpu_memory_cap_bytes=gpu_memory_cap_bytes,
        return_final_codes=True,
    )
    require(
        cap_codes is None
        and cap_result["decision_status"] == "skipped_state_cap"
        and cap_result["distinct_states_lower_bound_at_skip"] == 11
        and cap_result["skip_is_explicit_not_approximation"] is True,
        "CUDA state-cap self-audit failed",
    )

    semantic_audit = {
        "moduli": list(moduli),
        "expected_state_count": len(expected_codes),
        "expected_support_sha256_uint64": array_sha256(expected_codes),
        "present_target_code": present_code,
        "absent_target_code": absent_code,
        "large_chunk_support_sha256_uint64": large_chunk[
            "final_support_sha256_uint64"
        ],
        "small_chunk_support_sha256_uint64": small_chunk[
            "final_support_sha256_uint64"
        ],
        "trap_decision_status": trap_result["decision_status"],
        "cap_decision_status": cap_result["decision_status"],
        "cap_distinct_states_lower_bound": cap_result[
            "distinct_states_lower_bound_at_skip"
        ],
    }
    return {
        "status": "passed",
        "projected_group": "F3^6 x F7^3",
        "projected_mod7_coordinate_count_k": 3,
        "cpu_brute_force_final_state_count": len(expected_codes),
        "cpu_brute_force_support_sha256_uint64": array_sha256(expected_codes),
        "cuda_large_chunk_support_sha256_uint64": large_chunk[
            "final_support_sha256_uint64"
        ],
        "cuda_small_chunk_support_sha256_uint64": small_chunk[
            "final_support_sha256_uint64"
        ],
        "cpu_vs_gpu_exact_set_equality": True,
        "present_target_detected": True,
        "absent_target_rigorously_rejected": True,
        "pair_chunk_invariance_audited": True,
        "same_catalog_row_cross_prime_trap_rejected": True,
        "state_cap_produces_explicit_skip": True,
        "hash_collision_false_positive_path_exists": False,
        "semantic_audit_sha256": json_sha256(semantic_audit),
    }


def projected_supports(
    signatures: list[np.ndarray],
    target: np.ndarray,
    q3: int,
    mod7_subset: tuple[int, ...],
) -> tuple[list[np.ndarray], np.ndarray, tuple[int, ...]]:
    columns = tuple(range(q3)) + tuple(q3 + coordinate for coordinate in mod7_subset)
    supports = [np.ascontiguousarray(rows[:, columns], dtype=np.uint8) for rows in signatures]
    projected_target = np.ascontiguousarray(target[list(columns)], dtype=np.uint8)
    moduli = (3,) * q3 + (7,) * len(mod7_subset)
    group_size_for(moduli)
    return supports, projected_target, moduli


def projected_size_score(
    signatures: list[np.ndarray], q3: int, subset: tuple[int, ...]
) -> tuple[int, tuple[int, ...]]:
    zero = np.zeros(signatures[0].shape[1], dtype=np.uint8)
    supports, _target, moduli = projected_supports(signatures, zero, q3, subset)
    sizes = tuple(len(unique_support_codes(rows, moduli)) for rows in supports)
    return math.prod(sizes), sizes


def gpu_projection_subsets(
    signatures: list[np.ndarray],
    q3: int,
    q7: int,
    k: int,
    count: int,
) -> list[dict]:
    """Choose deterministic projections without the CPU prototype's k<=4 bound."""
    require(MIN_PROJECTED_MOD7_COORDINATES <= k <= q7, "invalid mod-7 projection width")
    require(1 <= count <= MAX_PROJECTIONS, "invalid projection count")
    require(group_size_for((3,) * q3 + (7,) * k) < int(EMPTY_KEY), "projection overflow")
    candidates: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    score_cache: dict[tuple[int, ...], tuple[int, tuple[int, ...]]] = {}

    def score(subset: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
        if subset not in score_cache:
            score_cache[subset] = projected_size_score(signatures, q3, subset)
        return score_cache[subset]

    def add(strategy: str, values: tuple[int, ...]) -> None:
        subset = tuple(sorted(int(value) for value in values))
        if len(subset) == k and len(set(subset)) == k and subset not in seen:
            require(all(0 <= value < q7 for value in subset), "projection coordinate escaped q7")
            seen.add(subset)
            candidates.append({"strategy": strategy, "mod7_coordinates": subset})

    add("prefix", tuple(range(k)))
    add("suffix", tuple(range(q7 - k, q7)))
    spread = tuple(sorted({int(index * q7 / k) for index in range(k)}))
    add("evenly_spaced", spread)

    for maximize, name in (
        (False, "greedy_min_support_product"),
        (True, "greedy_max_support_product"),
    ):
        selected: tuple[int, ...] = ()
        while len(selected) < k:
            choices = []
            for coordinate in range(q7):
                if coordinate in selected:
                    continue
                subset = tuple(sorted((*selected, coordinate)))
                support_product, sizes = score(subset)
                choices.append((support_product, coordinate, subset, sizes))
            selected = (
                max(choices, key=lambda row: (row[0], -row[1]))
                if maximize
                else min(choices, key=lambda row: (row[0], row[1]))
            )[2]
        add(name, selected)

    for subset in itertools.combinations(range(q7), k):
        if len(candidates) >= count:
            break
        add("lexicographic_fill", subset)
    rows = candidates[:count]
    require(len(rows) == count, "could not construct requested projection census")
    for row in rows:
        support_product, sizes = score(row["mod7_coordinates"])
        row["individual_support_product_score"] = support_product
        row["individual_projected_support_sizes"] = list(sizes)
    return rows


def maximum_uint64_k(q3: int, q7: int) -> int:
    admissible = [
        k
        for k in range(q7 + 1)
        if 3**q3 * 7**k < int(EMPTY_KEY)
    ]
    require(bool(admissible), "no uint64 projection width exists")
    return max(admissible)


def run_four_cases(
    cp: Any,
    kernel: Any,
    gpu_audit: dict,
    state_cap: int,
    pair_chunk_cap: int,
    gpu_memory_cap_bytes: int,
    k: int,
    projection_count: int,
) -> dict:
    started = time.time()
    gpu_self_audit = manufactured_cpu_gpu_self_audit(
        cp, kernel, gpu_memory_cap_bytes
    )
    source = load_projection_module()
    context = source.symmetry.construct_pointed_systems()
    system = context["systems"][0]["A"]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    require(matrix.shape == (282, 1_225), "orbit0/A pointed matrix shape changed")
    kernel_rows = np.asarray(
        source.johnson._primitive_left_kernel_rows(),  # noqa: SLF001
        dtype=np.int64,
    )
    require(kernel_rows.shape == (14, 35), "primitive Johnson kernel shape changed")
    common, common_audit = source.exact_common_dependency_basis(matrix, kernel_rows)
    quotient_data, quotient_audit = source.derive_torsion_quotients(matrix, common)
    q3 = len(quotient_data[3]["complement"])
    q7 = len(quotient_data[7]["complement"])
    max_k = maximum_uint64_k(q3, q7)
    require(
        MIN_PROJECTED_MOD7_COORDINATES <= k <= max_k,
        f"k must lie in {MIN_PROJECTED_MOD7_COORDINATES}..{max_k} for uint64 encoding",
    )

    leaves_by_orbit, leaf_audit = source.parent.exact_mean_leaves(context["orbits"])
    selected = [
        (leaf_index, leaf)
        for leaf_index, leaf in enumerate(leaves_by_orbit[0])
        if tuple(int(value) for value in leaf["pattern"]) == source.TARGET_PATTERN
    ]
    require(len(selected) == 4, "orbit0 H0_S0_M7 representative census changed")
    require(
        all(
            Counter(leaf["catalog_classes"]) == Counter({"M": 7, "U": 1})
            and not leaf["high_directions"]
            for _leaf_index, leaf in selected
        ),
        "selected leaf is not seven exact M catalogs plus one U",
    )

    hull_kernel, _hull_bases, hull_audit = source.affine.build_hull_audit()
    require(np.array_equal(hull_kernel, kernel_rows), "hull and compact kernel rows disagree")
    source_catalog = source.affine.canonical_catalog(7, 4).astype(np.int64)
    anchors = source.affine.AnchorFactory(kernel_rows, source_catalog[1:] - source_catalog[0])
    signature_cache: dict = {}
    case_results = []
    projection_totals: Counter[str] = Counter()

    for case_index, (leaf_index, leaf) in enumerate(selected):
        case_audit, signatures, target = source.quotient_case_signatures(
            case_index,
            leaf_index,
            leaf,
            context["orbits"][0],
            system,
            common,
            quotient_data,
            anchors,
            signature_cache,
        )
        subset_rows = gpu_projection_subsets(
            signatures, q3, q7, k, projection_count
        )
        projections = []
        for subset_row in subset_rows:
            subset = tuple(subset_row["mod7_coordinates"])
            raw_supports, projected_target, moduli = projected_supports(
                signatures, target, q3, subset
            )
            decision, _final_codes = gpu_exact_support_convolution(
                cp,
                kernel,
                raw_supports,
                projected_target,
                moduli,
                state_cap,
                pair_chunk_cap,
                gpu_memory_cap_bytes,
            )
            projection_totals[decision["decision_status"]] += 1
            projections.append(
                {
                    **subset_row,
                    "retained_all_mod3_quotient_coordinates": list(range(q3)),
                    "retained_mod7_quotient_coordinates": list(subset),
                    "projected_group": f"F3^{q3} x F7^{len(subset)}",
                    "projected_target_sha256_uint8": array_sha256(projected_target),
                    "same_catalog_row_identity_used_across_projected_primes": True,
                    **decision,
                }
            )
            cp.get_default_memory_pool().free_all_blocks()

        rejected = any(row["rigorous_rejection"] for row in projections)
        skipped = any(row["decision_status"].startswith("skipped_") for row in projections)
        case_results.append(
            {
                **case_audit,
                "projection_count": len(projections),
                "projections": projections,
                "rigorously_rejected_by_at_least_one_exact_projection": rejected,
                "any_projection_explicitly_skipped": skipped,
                "case_status": (
                    "rigorously_rejected_by_exact_gpu_projected_support"
                    if rejected
                    else "bounded_gpu_projected_support_survivor_with_skips"
                    if skipped
                    else "necessary_only_survivor_of_all_tested_exact_gpu_projections"
                ),
            }
        )

    rejected_cases = sum(
        row["rigorously_rejected_by_at_least_one_exact_projection"]
        for row in case_results
    )
    return {
        "experiment": "p7_infinity7_positive_z7_torsion_support_gpu",
        "status": "complete_bounded_exact_cuda_torsion_support_projection",
        "p": source.P,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "scope": "four orbit0/A H0_S0_M7 affine-hull survivors",
        "solver_invoked": False,
        "cuda_invoked": True,
        "configuration": {
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
            "gpu_memory_cap_bytes": gpu_memory_cap_bytes,
            "retained_all_mod3_quotient_coordinates": True,
            "selected_mod7_coordinate_count_k": k,
            "projection_subsets_per_case": projection_count,
            "minimum_supported_k": MIN_PROJECTED_MOD7_COORDINATES,
            "maximum_uint64_supported_k_derived": max_k,
        },
        "gpu_engine_audit": gpu_audit,
        "manufactured_cpu_gpu_self_audit": gpu_self_audit,
        "source_provenance": {
            "cpu_projection_module": "p7_infinity7_positive_z7_torsion_support_projection.py",
            "orbit_source": context["orbit_source"],
            "mean_leaf_coverage": leaf_audit,
            "pointed_case": "orbit0/A",
            "pointed_fixed_edge_rows": system["fixed_edge_rows"],
            "pointed_base_rhs_sha256_int64": array_sha256(
                np.asarray(system["base_rhs"], dtype=np.int64)
            ),
            "selected_leaf_indices": [index for index, _leaf in selected],
            "selected_case_keys": [
                f"orbit0_leaf{index}_branchA" for index, _leaf in selected
            ],
        },
        "exact_common_rational_dependency_audit": common_audit,
        "varying_torsion_quotient_audit": quotient_audit,
        "degree_two_zero_mean_hull_audit": hull_audit,
        "processed_cases": len(case_results),
        "rigorously_rejected_cases": rejected_cases,
        "necessary_only_or_skipped_cases": len(case_results) - rejected_cases,
        "projection_decision_census": dict(sorted(projection_totals.items())),
        "case_results": case_results,
        "logical_semantics": {
            "quotient_dimensions_are_derived_not_assumed": True,
            "all_seven_M_catalogs_are_complete_exact_1764_row_sets": True,
            "same_catalog_row_index_supplies_mod3_and_mod7_signatures": True,
            "coordinate_projection_is_an_exact_group_homomorphism": True,
            "cuda_hash_collisions_are_resolved_exactly": True,
            "projected_target_absence_rigorously_rejects_the_full_global_join_case": True,
            "projected_target_presence_is_necessary_only": True,
            "state_or_memory_cap_excess_is_an_explicit_skip_not_a_rejection": True,
            "raw_binary_edge_feasibility_claimed": False,
            "full_z7_closure_claimed": False,
        },
        "all_selected_case_and_projection_audits_passed": True,
        "elapsed_seconds": time.time() - started,
    }


def self_audit_only(cp: Any, kernel: Any, gpu_audit: dict, memory_bytes: int) -> dict:
    started = time.time()
    return {
        "experiment": "p7_infinity7_positive_z7_torsion_support_gpu_self_audit",
        "status": "manufactured_cpu_gpu_self_audit_passed",
        "gpu_engine_audit": gpu_audit,
        "manufactured_cpu_gpu_self_audit": manufactured_cpu_gpu_self_audit(
            cp, kernel, memory_bytes
        ),
        "real_H0_S0_M7_cases_processed": 0,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--state-cap", type=int, default=DEFAULT_STATE_CAP)
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    parser.add_argument("--gpu-memory-cap-mib", type=int, default=DEFAULT_GPU_MEMORY_CAP_MIB)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--projections", type=int, default=4)
    parser.add_argument("--self-audit-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    require(args.state_cap > 0, "--state-cap must be positive")
    require(args.pair_chunk_cap > 0, "--pair-chunk-cap must be positive")
    require(args.gpu_memory_cap_mib > 0, "--gpu-memory-cap-mib must be positive")
    require(1 <= args.projections <= MAX_PROJECTIONS, "--projections lies outside bounds")
    if not args.self_audit_only:
        require(args.output is not None, "--output is required for a real four-case run")

    cp, kernel, gpu_audit = load_cupy(args.device)
    memory_bytes = args.gpu_memory_cap_mib * (1 << 20)
    with cp.cuda.Device(args.device):
        result = (
            self_audit_only(cp, kernel, gpu_audit, memory_bytes)
            if args.self_audit_only
            else run_four_cases(
                cp,
                kernel,
                gpu_audit,
                args.state_cap,
                args.pair_chunk_cap,
                memory_bytes,
                args.k,
                args.projections,
            )
        )
    if args.output is not None:
        atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output) if args.output is not None else None,
                "gpu": result["gpu_engine_audit"]["device_name"],
                "self_audit": result["manufactured_cpu_gpu_self_audit"]["status"],
                "processed_cases": result.get("processed_cases", 0),
                "rigorously_rejected_cases": result.get("rigorously_rejected_cases", 0),
                "projection_decision_census": result.get("projection_decision_census", {}),
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
