#!/usr/bin/env python3
"""Exact V100 sieve for the doubly saturated p=7 size-eight stratum.

For ``c_H=-1``, this script selects every finite eight-point boundary whose
two quadratic direction types both have exact floor sum 32.  Since each type
also has exact mean sum 32, every directional mean is forced to its floor.
All floor catalogs are singletons except ``b=4, phase=1, mean=14``, whose
complete Johnson-slice catalog has only 36 rows.  A floor-(32,32) boundary
has at most one such direction.  Consequently a complete mod-seven
dependency check is a zero-table test or a one-table membership test.

The CUDA pass uses eight deterministic linear combinations of the 135 exact
left dependencies only as a lossless rejection prefilter: every projected
survivor is rechecked against all 135 dependencies on the host.  Therefore
the final exclusion is exact regardless of projection collisions.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import scaled_direction_floor  # noqa: E402
from p7_exceptional_omit_high_catalogs import modular_rank  # noqa: E402
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_conic_orbits import unrank_lex  # noqa: E402
from p7_size8_floor_profile_gpu import (  # noqa: E402
    choose_table,
    direction_tables,
)


PROJECTION_DIMENSION = 8
PROJECTION_SEED = 15_662_032
MAX_CATALOG_ROWS = 36


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def floor_for(eps: int, b: int) -> int:
    phase = int(eps == -1)
    return int(scaled_direction_floor(7, b, phase))


def type_costs(profile: tuple[int, ...], epsilons: tuple[int, ...]) -> tuple[int, int]:
    return tuple(
        sum(floor_for(epsilons[index], b) for index, b in enumerate(profile) if epsilons[index] == eps)
        for eps in (-1, 1)
    )


def expected_scope(source: dict, epsilons: tuple[int, ...]) -> dict:
    if (
        source.get("experiment") != "p7_size8_floor_profile_gpu"
        or source.get("status") != "complete_exact_floor_profile_census"
        or int(source.get("p", 0)) != 7
        or int(source.get("c_H", 0)) != -1
        or int(source.get("checked_boundaries", 0)) != math.comb(49, 8)
    ):
        raise ValueError("source is not the complete c_H=-1 p=7 size-eight floor census")
    profile_count = 0
    boundary_count = 0
    odd_histogram: Counter[int] = Counter()
    for row in source["survivor_ordered_profiles"]:
        profile = tuple(int(value) for value in row["b_by_direction"])
        count = int(row["count"])
        if type_costs(profile, epsilons) == (32, 32):
            profile_count += 1
            boundary_count += count
            odd_histogram[sum(profile)] += count
    if odd_histogram.get(8, 0):
        raise AssertionError("the forced-floor stratum unexpectedly contains a conic profile")
    return {
        "ordered_profile_count": profile_count,
        "boundary_count": boundary_count,
        "odd_secant_histogram": dict(sorted(odd_histogram.items())),
    }


def projection_family(dependency: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(PROJECTION_SEED)
    for _attempt in range(10_000):
        coefficients = rng.integers(
            0, 7, size=(PROJECTION_DIMENSION, len(dependency)), dtype=np.int64
        )
        projected = coefficients @ dependency.astype(np.int64) % 7
        if modular_rank(coefficients) != PROJECTION_DIMENSION:
            continue
        if modular_rank(projected) != PROJECTION_DIMENSION:
            continue
        block_ranks = [
            modular_rank(projected[:, 2 + 35 * d : 2 + 35 * (d + 1)])
            for d in range(8)
        ]
        if min(block_ranks) == PROJECTION_DIMENSION:
            return coefficients.astype(np.uint8), projected.astype(np.uint8)
    raise RuntimeError("failed to construct a full-rank deterministic projection")


def contribution(
    dependency: np.ndarray,
    direction_index: int,
    mask: int,
    eps: int,
) -> np.ndarray:
    b = mask.bit_count()
    phase = int(eps == -1)
    mean = floor_for(eps, b)
    B = {index for index in range(7) if mask & (1 << index)}
    values = mapped_catalog(b, phase, mean, B, None).astype(np.int64)
    bad = 13 - values
    block = dependency[:, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)].astype(np.int64)
    return (block @ (bad.T % 7) % 7).astype(np.uint8)


def build_tables(
    dependency: np.ndarray,
    epsilons: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    base = (
        dependency[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % 7
    ).astype(np.uint8)
    singleton = np.full((8, 128, len(dependency)), 255, dtype=np.uint8)
    variable = np.zeros(
        (8, 128, MAX_CATALOG_ROWS, len(dependency)), dtype=np.uint8
    )
    variable_count = np.zeros((8, 128), dtype=np.uint8)
    catalog_histogram: Counter[int] = Counter()
    for direction_index, eps in enumerate(epsilons):
        for mask in range(128):
            if mask.bit_count() not in (0, 2, 4, 6):
                continue
            rows = contribution(dependency, direction_index, mask, eps)
            catalog_rows = int(rows.shape[1])
            catalog_histogram[catalog_rows] += 1
            if catalog_rows == 1:
                singleton[direction_index, mask] = rows[:, 0]
            elif (
                catalog_rows == MAX_CATALOG_ROWS
                and eps == -1
                and mask.bit_count() == 4
            ):
                variable[direction_index, mask] = rows.T
                variable_count[direction_index, mask] = catalog_rows
            else:
                raise AssertionError(
                    "unexpected forced-floor catalog: "
                    f"d={direction_index} mask={mask} rows={catalog_rows}"
                )
    return base, singleton, variable, variable_count, {
        "catalog_count_histogram": dict(sorted(catalog_histogram.items())),
        "singleton_table_sha256": array_sha256(singleton),
        "variable_table_sha256": array_sha256(variable),
    }


CUDA_SOURCE = r'''
extern "C" __global__
void forced_floor_scan(
    const signed char* labels,
    const signed char* epsilons,
    const unsigned long long* choose,
    const int choose_stride,
    const unsigned char* base,
    const unsigned char* singleton,
    const unsigned char* variable,
    const unsigned char* variable_count,
    const unsigned long long start_rank,
    const unsigned long long stop_rank,
    unsigned long long* checked,
    unsigned long long* forced,
    unsigned long long* projected_pass,
    unsigned long long* odd_histogram,
    unsigned long long* survivor_ranks,
    const unsigned long long survivor_capacity)
{
    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;
    unsigned long long local_checked = 0;
    unsigned long long local_forced = 0;

    for (unsigned long long rank0 = start_rank + logical_thread;
         rank0 < stop_rank;
         rank0 += logical_threads)
    {
        int index[8];
        unsigned long long rank = rank0;
        int next = 0;
        bool valid = true;
        #pragma unroll
        for (int position = 0; position < 8; ++position)
        {
            const int remaining = 7 - position;
            bool selected = false;
            const int last = 49 - (remaining + 1);
            for (int candidate = next; candidate <= last; ++candidate)
            {
                const unsigned long long ways =
                    choose[(49 - candidate - 1) * choose_stride + remaining];
                if (rank < ways)
                {
                    index[position] = candidate;
                    next = candidate + 1;
                    selected = true;
                    break;
                }
                rank -= ways;
            }
            if (!selected) valid = false;
        }
        if (!valid) continue;
        ++local_checked;

        int masks[8];
        int type_minus = 0;
        int type_plus = 0;
        int odd_secants = 0;
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction)
        {
            int mask = 0;
            #pragma unroll
            for (int column = 0; column < 8; ++column)
                mask ^= 1 << (int)labels[49 * direction + index[column]];
            masks[direction] = mask;
            const int b = __popc((unsigned int)mask);
            odd_secants += b;
            const int eps = (int)epsilons[direction];
            int cost;
            if (eps > 0)
                cost = b == 0 ? 0 : 8;
            else
                cost = (b == 0 || b == 4) ? 14 : 6;
            if (eps < 0) type_minus += cost;
            else type_plus += cost;
        }
        if (type_minus != 32 || type_plus != 32) continue;
        ++local_forced;
        atomicAdd(odd_histogram + odd_secants, 1ULL);

        int syndrome[8];
        #pragma unroll
        for (int row = 0; row < 8; ++row) syndrome[row] = (int)base[row];
        int variable_direction = -1;
        int variable_mask = -1;
        bool table_valid = true;
        #pragma unroll
        for (int direction = 0; direction < 8; ++direction)
        {
            const int mask = masks[direction];
            const int count = (int)variable_count[128 * direction + mask];
            if (count)
            {
                if (variable_direction >= 0) table_valid = false;
                variable_direction = direction;
                variable_mask = mask;
            }
            else
            {
                #pragma unroll
                for (int row = 0; row < 8; ++row)
                {
                    const int value = (int)singleton[(128 * direction + mask) * 8 + row];
                    if (value == 255) table_valid = false;
                    syndrome[row] += value;
                }
            }
        }
        if (!table_valid) continue;

        bool passing = false;
        if (variable_direction < 0)
        {
            passing = true;
            #pragma unroll
            for (int row = 0; row < 8; ++row)
                if (syndrome[row] % 7 != 0) passing = false;
        }
        else
        {
            const int offset =
                ((128 * variable_direction + variable_mask) * 36) * 8;
            for (int catalog_row = 0; catalog_row < 36 && !passing; ++catalog_row)
            {
                bool match = true;
                #pragma unroll
                for (int row = 0; row < 8; ++row)
                {
                    const int value = (int)variable[offset + 8 * catalog_row + row];
                    if ((syndrome[row] + value) % 7 != 0) match = false;
                }
                passing = match;
            }
        }
        if (passing)
        {
            const unsigned long long slot = atomicAdd(projected_pass, 1ULL);
            if (slot < survivor_capacity) survivor_ranks[slot] = rank0;
        }
    }
    if (local_checked) atomicAdd(checked, local_checked);
    if (local_forced) atomicAdd(forced, local_forced);
}
'''


def gpu_scan(
    kernel,
    labels_gpu,
    epsilons_gpu,
    choose_gpu,
    base_gpu,
    singleton_gpu,
    variable_gpu,
    variable_count_gpu,
    start_rank: int,
    stop_rank: int,
    blocks: int,
    threads: int,
    survivor_capacity: int,
) -> dict:
    import cupy as cp

    checked = cp.zeros(1, dtype=cp.uint64)
    forced = cp.zeros(1, dtype=cp.uint64)
    passed = cp.zeros(1, dtype=cp.uint64)
    histogram = cp.zeros(57, dtype=cp.uint64)
    ranks = cp.zeros(survivor_capacity, dtype=cp.uint64)
    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            labels_gpu,
            epsilons_gpu,
            choose_gpu,
            np.int32(choose_gpu.shape[1]),
            base_gpu,
            singleton_gpu,
            variable_gpu,
            variable_count_gpu,
            np.uint64(start_rank),
            np.uint64(stop_rank),
            checked,
            forced,
            passed,
            histogram,
            ranks,
            np.uint64(survivor_capacity),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    passed_count = int(cp.asnumpy(passed)[0])
    if passed_count > survivor_capacity:
        raise RuntimeError(
            f"projection survivor capacity exceeded: {passed_count}>{survivor_capacity}"
        )
    histogram_host = cp.asnumpy(histogram)
    return {
        "checked": int(cp.asnumpy(checked)[0]),
        "forced": int(cp.asnumpy(forced)[0]),
        "projected_pass": passed_count,
        "odd_histogram": {
            int(index): int(value)
            for index, value in enumerate(histogram_host)
            if value
        },
        "survivor_ranks": sorted(int(value) for value in cp.asnumpy(ranks[:passed_count])),
        "elapsed_seconds": elapsed,
    }


def masks_for_boundary(
    boundary: tuple[int, ...], labels: np.ndarray
) -> tuple[int, ...]:
    return tuple(
        int(np.bitwise_xor.reduce(1 << labels[direction, list(boundary)]))
        for direction in range(8)
    )


def projected_passes(
    masks: tuple[int, ...],
    base: np.ndarray,
    singleton: np.ndarray,
    variable: np.ndarray,
    variable_count: np.ndarray,
) -> bool:
    syndrome = base.astype(np.int16).copy()
    variable_key = None
    for direction, mask in enumerate(masks):
        if variable_count[direction, mask]:
            if variable_key is not None:
                raise AssertionError("forced boundary has two variable catalogs")
            variable_key = (direction, mask)
        else:
            row = singleton[direction, mask]
            if np.any(row == 255):
                raise AssertionError("missing singleton contribution")
            syndrome += row
    if variable_key is None:
        return bool(np.all(syndrome % 7 == 0))
    direction, mask = variable_key
    values = variable[direction, mask, : int(variable_count[direction, mask])]
    return bool(np.any(np.all((values.astype(np.int16) + syndrome) % 7 == 0, axis=1)))


def cpu_prefix(
    stop_rank: int,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    base: np.ndarray,
    singleton: np.ndarray,
    variable: np.ndarray,
    variable_count: np.ndarray,
) -> dict:
    forced = 0
    passed = []
    histogram: Counter[int] = Counter()
    for rank, boundary in enumerate(itertools.islice(itertools.combinations(range(49), 8), stop_rank)):
        masks = masks_for_boundary(boundary, labels)
        profile = tuple(mask.bit_count() for mask in masks)
        if type_costs(profile, epsilons) != (32, 32):
            continue
        forced += 1
        histogram[sum(profile)] += 1
        if projected_passes(masks, base, singleton, variable, variable_count):
            passed.append(rank)
    return {
        "checked": stop_rank,
        "forced": forced,
        "projected_pass": len(passed),
        "odd_histogram": dict(sorted(histogram.items())),
        "survivor_ranks": passed,
    }


def exact_recheck(
    ranks: list[int],
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    dependency: np.ndarray,
) -> tuple[list[int], list[dict]]:
    base = (
        dependency[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % 7
    ).astype(np.int64)
    cache: dict[tuple[int, int], np.ndarray] = {}
    survivors = []
    rows = []
    for rank in ranks:
        boundary = unrank_lex(rank)
        masks = masks_for_boundary(boundary, labels)
        syndrome = base.copy()
        variable = None
        variable_metadata = None
        for direction, mask in enumerate(masks):
            key = (direction, mask)
            if key not in cache:
                cache[key] = contribution(dependency, direction, mask, epsilons[direction]).astype(np.int64)
            values = cache[key]
            if values.shape[1] == 1:
                syndrome = (syndrome + values[:, 0]) % 7
            else:
                if variable is not None:
                    raise AssertionError("exact forced-floor recheck found two variable catalogs")
                variable = values
                variable_metadata = {
                    "direction_index": direction,
                    "odd_fibre_mask": mask,
                    "catalog_rows": int(values.shape[1]),
                }
        if variable is None:
            matches = int(not np.any(syndrome % 7))
        else:
            matches = int(
                np.count_nonzero(
                    np.all((variable.astype(np.int64) + syndrome[:, None]) % 7 == 0, axis=0)
                )
            )
        row = {
            "rank": rank,
            "boundary_finite_field": list(boundary),
            "boundary_vertices": [value + 1 for value in boundary],
            "odd_fibre_masks": list(masks),
            "variable_catalog": variable_metadata,
            "full_dependency_matches": matches,
        }
        rows.append(row)
        if matches:
            survivors.append(rank)
    return survivors, rows


def run(args: argparse.Namespace) -> dict:
    import cupy as cp

    started = time.time()
    source = json.loads(args.source.read_text())
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    scope = expected_scope(source, epsilons)

    matrix, dependencies, linear_rows = linear_data((7,))
    dependency = dependencies[7].astype(np.uint8)
    if (
        matrix.shape != (282, 1225)
        or dependency.shape != (135, 282)
        or linear_rows != [
            {
                "modulus": 7,
                "rank": 147,
                "left_dependency_dimension": 135,
                "left_null_audit": True,
            }
        ]
    ):
        raise AssertionError("unexpected common mod-seven linear system")
    if np.any(dependency.astype(np.int64) @ (matrix.astype(np.int64) % 7) % 7):
        raise AssertionError("full left-null dependency audit failed")
    coefficients, projected = projection_family(dependency)
    if np.any(projected.astype(np.int64) @ (matrix.astype(np.int64) % 7) % 7):
        raise AssertionError("projected left-null audit failed")
    block_ranks = [
        modular_rank(projected[:, 2 + 35 * d : 2 + 35 * (d + 1)])
        for d in range(8)
    ]
    base, singleton, variable, variable_count, table_audit = build_tables(
        projected, epsilons
    )

    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("CUDA is unavailable")
    kernel = cp.RawKernel(CUDA_SOURCE, "forced_floor_scan", options=("--std=c++11",))
    labels_gpu = cp.asarray(labels, dtype=cp.int8)
    epsilons_gpu = cp.asarray(epsilon_array, dtype=cp.int8)
    choose_gpu = cp.asarray(choose_table(49, 8))
    base_gpu = cp.asarray(base)
    singleton_gpu = cp.asarray(singleton)
    variable_gpu = cp.asarray(variable)
    variable_count_gpu = cp.asarray(variable_count)

    prefix_gpu = gpu_scan(
        kernel,
        labels_gpu,
        epsilons_gpu,
        choose_gpu,
        base_gpu,
        singleton_gpu,
        variable_gpu,
        variable_count_gpu,
        0,
        args.verify_prefix,
        args.blocks,
        args.threads,
        args.survivor_capacity,
    )
    prefix_cpu = cpu_prefix(
        args.verify_prefix,
        labels,
        epsilons,
        base,
        singleton,
        variable,
        variable_count,
    )
    prefix_gpu_comparable = {
        key: value for key, value in prefix_gpu.items() if key != "elapsed_seconds"
    }
    if prefix_gpu_comparable != prefix_cpu:
        raise AssertionError(
            f"independent CPU/GPU prefix mismatch: gpu={prefix_gpu_comparable}, cpu={prefix_cpu}"
        )

    total = math.comb(49, 8)
    full = gpu_scan(
        kernel,
        labels_gpu,
        epsilons_gpu,
        choose_gpu,
        base_gpu,
        singleton_gpu,
        variable_gpu,
        variable_count_gpu,
        0,
        total,
        args.blocks,
        args.threads,
        args.survivor_capacity,
    )
    if full["checked"] != total:
        raise AssertionError("full CUDA rank interval was not exhausted")
    if full["forced"] != scope["boundary_count"]:
        raise AssertionError("CUDA forced-floor count disagrees with the floor-profile source")
    if full["odd_histogram"] != scope["odd_secant_histogram"]:
        raise AssertionError("CUDA forced-floor odd-secant histogram changed")

    exact_survivors, exact_rows = exact_recheck(
        full["survivor_ranks"], labels, epsilons, dependency
    )
    properties = cp.cuda.runtime.getDeviceProperties(0)
    device = properties["name"]
    if isinstance(device, bytes):
        device = device.decode()
    return {
        "experiment": "p7_size8_forced_floor_gpu",
        "status": "complete_exact_doubly_saturated_boundary_exhaustion",
        "p": 7,
        "c_H": -1,
        "finite_boundary_size": 8,
        "source": str(args.source),
        "source_sha256": sha256(args.source),
        "source_scope": scope,
        "mathematical_reduction": {
            "type_floor_sums": [32, 32],
            "exact_type_mean_sums": [32, 32],
            "all_directional_means_forced_to_floor": True,
            "maximum_variable_catalogs_per_boundary": 1,
            "sole_variable_floor_catalog": {
                "odd_fibres": 4,
                "phase": 1,
                "scaled_mean": 14,
                "rows": 36,
            },
        },
        "linear_system": linear_rows,
        "full_dependency_shape": list(dependency.shape),
        "full_dependency_sha256": array_sha256(dependency),
        "projection": {
            "seed": PROJECTION_SEED,
            "dimension": PROJECTION_DIMENSION,
            "coefficient_rank_mod7": modular_rank(coefficients),
            "projected_rank_mod7": modular_rank(projected),
            "direction_block_ranks_mod7": block_ranks,
            "coefficients_mod7": coefficients.astype(int).tolist(),
            "projected_dependency_sha256": array_sha256(projected),
            "used_only_as_rejection_prefilter": True,
        },
        "catalog_tables": table_audit,
        "verification": {
            "method": "independent itertools CPU prefix versus direct-rank CUDA",
            "prefix_checked": args.verify_prefix,
            "all_counts_histograms_and_projected_ranks_match": True,
        },
        "device": str(device),
        "launch": {"blocks": args.blocks, "threads": args.threads},
        "all_boundaries": total,
        "checked_boundaries": full["checked"],
        "forced_floor_boundaries": full["forced"],
        "forced_floor_odd_secant_histogram": {
            str(key): value for key, value in full["odd_histogram"].items()
        },
        "projected_dependency_survivors": full["projected_pass"],
        "projected_survivor_ranks": full["survivor_ranks"],
        "full_135_dependency_recheck_rows": exact_rows,
        "full_dependency_survivors": len(exact_survivors),
        "full_dependency_survivor_ranks": exact_survivors,
        "all_forced_floor_boundaries_mod7_infeasible": not exact_survivors,
        "proved_forced_floor_stratum_exclusion_cminus1": not exact_survivors,
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "elapsed_seconds": time.time() - started,
        "cuda_scan_seconds": full["elapsed_seconds"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--verify-prefix", type=int, default=100_000)
    parser.add_argument("--survivor-capacity", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.verify_prefix <= math.comb(49, 8):
        raise ValueError("verify-prefix is outside the rank interval")
    out = run(args)
    atomic_json(args.output, out)
    compact = {
        key: value
        for key, value in out.items()
        if key
        not in {
            "projected_survivor_ranks",
            "full_135_dependency_recheck_rows",
            "full_dependency_survivor_ranks",
        }
    }
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
