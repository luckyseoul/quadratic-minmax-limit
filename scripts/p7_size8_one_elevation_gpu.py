#!/usr/bin/env python3
"""Exact V100 sieve for the p=7 size-eight four-allocation stratum.

After the conic and doubly forced-floor strata have been removed, almost all
remaining boundaries have one quadratic type with floor sum 24 and the other
with floor sum 32.  There are exactly four admissible mean allocations: raise
one of the four deficient-type directions by eight and leave every other
direction at its floor.

For each possible raised direction this script uses exact mod-seven left
dependencies which vanish on all 35 score columns of that direction.  The
large raised catalog then contributes identically zero.  A complete direct-
rank CUDA scan tests all four allocations of every boundary in this stratum.
The 22-coordinate test is only a necessary-condition prefilter; every
projected survivor is rechecked with all 135 common dependencies on the host.
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
from p7_fixed_boundary_mean_allocation_batch import POINTS, allocations  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_size8_conic_orbits import unrank_lex  # noqa: E402
from p7_size8_floor_profile_gpu import choose_table, direction_tables  # noqa: E402


MODULUS = 7
PROJECTION_ROWS = 22
MAX_FLOOR_CATALOG_ROWS = 36
EXPECTED_ALLOCATION_BOUNDARIES = {
    4: 23_563_806,
    11: 154_056,
    16: 1_194_816,
    24: 1_176,
    44: 69_384,
}
EXPECTED_ALLOCATION_PROFILES = {4: 2_245, 11: 248, 16: 516, 24: 8, 44: 110}
EXPECTED_FOUR_ODD_HISTOGRAM = {
    16: 691_488,
    20: 5_603_640,
    24: 9_190_146,
    28: 5_990_544,
    32: 1_846_908,
    36: 232_848,
    40: 5_880,
    44: 2_352,
}
EXPECTED_FOUR_FLOOR_PAIRS = {(24, 32): 17_298_078, (32, 24): 6_265_728}


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
    return int(scaled_direction_floor(7, b, int(eps == -1)))


def type_costs(profile: tuple[int, ...], epsilons: tuple[int, ...]) -> tuple[int, int]:
    return tuple(
        sum(
            floor_for(epsilons[index], b)
            for index, b in enumerate(profile)
            if epsilons[index] == eps
        )
        for eps in (-1, 1)
    )


def abstract_rows(profile: tuple[int, ...], epsilons: tuple[int, ...]) -> list[dict]:
    rows = []
    for b, eps in zip(profile, epsilons):
        phase = int(eps == -1)
        floor = floor_for(eps, b)
        odd_fibres = set(range(b))
        parity_mass = sum(
            (sum(value in odd_fibres for value in point) + phase) & 1
            for point in POINTS
        )
        allowed = tuple(
            mean
            for mean in range(floor, 33, 2)
            if 5 * mean >= 2 * parity_mass
            and (5 * mean - 2 * parity_mass) % 4 == 0
        )
        rows.append(
            {
                "eps": int(eps),
                "b": int(b),
                "phase": phase,
                "floor": floor,
                "allowed_scaled_means": allowed,
            }
        )
    return rows


def expected_scope(source: dict, epsilons: tuple[int, ...]) -> dict:
    if (
        source.get("experiment") != "p7_size8_floor_profile_gpu"
        or source.get("status") != "complete_exact_floor_profile_census"
        or int(source.get("p", 0)) != 7
        or int(source.get("c_H", 0)) != -1
        or int(source.get("checked_boundaries", 0)) != math.comb(49, 8)
    ):
        raise ValueError("source is not the complete c_H=-1 p=7 size-eight census")

    allocation_boundaries: Counter[int] = Counter()
    allocation_profiles: Counter[int] = Counter()
    four_odd: Counter[int] = Counter()
    four_floor_pairs: Counter[tuple[int, int]] = Counter()
    conic_boundaries = 0
    forced_floor_boundaries = 0
    selected_profiles = 0
    for source_row in source["survivor_ordered_profiles"]:
        profile = tuple(int(value) for value in source_row["b_by_direction"])
        count = int(source_row["count"])
        floors = type_costs(profile, epsilons)
        if sum(profile) == 8:
            conic_boundaries += count
            continue
        if floors == (32, 32):
            forced_floor_boundaries += count
            continue

        rows = abstract_rows(profile, epsilons)
        leaves = allocations(rows)
        allocation_boundaries[len(leaves)] += count
        allocation_profiles[len(leaves)] += 1
        if len(leaves) != 4:
            if floors in EXPECTED_FOUR_FLOOR_PAIRS:
                raise AssertionError("a (24,32)/(32,24) profile lost its four-leaf shape")
            continue
        if floors not in EXPECTED_FOUR_FLOOR_PAIRS:
            raise AssertionError("a four-allocation profile has an unexpected floor pair")

        deficient_eps = -1 if floors == (24, 32) else 1
        floor_vector = tuple(row["floor"] for row in rows)
        elevated = []
        for leaf in leaves:
            differences = tuple(value - floor for value, floor in zip(leaf, floor_vector))
            support = [index for index, value in enumerate(differences) if value]
            if (
                len(support) != 1
                or differences[support[0]] != 8
                or epsilons[support[0]] != deficient_eps
            ):
                raise AssertionError("four-allocation leaf is not one deficient direction +8")
            elevated.append(support[0])
        expected_elevated = sorted(
            index for index, eps in enumerate(epsilons) if eps == deficient_eps
        )
        if sorted(elevated) != expected_elevated:
            raise AssertionError("four-allocation leaves do not cover the deficient type")
        selected_profiles += 1
        four_odd[sum(profile)] += count
        four_floor_pairs[floors] += count

    if dict(allocation_boundaries) != EXPECTED_ALLOCATION_BOUNDARIES:
        raise AssertionError(f"remaining allocation boundary census changed: {allocation_boundaries}")
    if dict(allocation_profiles) != EXPECTED_ALLOCATION_PROFILES:
        raise AssertionError(f"remaining allocation profile census changed: {allocation_profiles}")
    if dict(four_odd) != EXPECTED_FOUR_ODD_HISTOGRAM:
        raise AssertionError(f"four-allocation odd histogram changed: {four_odd}")
    if dict(four_floor_pairs) != EXPECTED_FOUR_FLOOR_PAIRS:
        raise AssertionError(f"four-allocation floor-pair census changed: {four_floor_pairs}")
    if sum(allocation_boundaries.values()) != 24_983_238:
        raise AssertionError("post-forced-floor remainder changed")
    return {
        "conic_boundaries_previously_excluded": conic_boundaries,
        "forced_floor_boundaries_previously_excluded": forced_floor_boundaries,
        "remaining_boundaries_before_this_stratum": sum(allocation_boundaries.values()),
        "allocation_count_boundary_histogram": dict(sorted(allocation_boundaries.items())),
        "allocation_count_ordered_profile_histogram": dict(sorted(allocation_profiles.items())),
        "mean_allocation_count": 4,
        "ordered_profile_count": selected_profiles,
        "boundary_count": allocation_boundaries[4],
        "allocation_leaf_count": 4 * allocation_boundaries[4],
        "type_floor_pair_boundary_histogram": [
            {"type_floor_sums": list(key), "boundaries": value}
            for key, value in sorted(four_floor_pairs.items())
        ],
        "odd_secant_histogram": dict(sorted(four_odd.items())),
        "remaining_boundaries_after_stratum_if_excluded": (
            sum(allocation_boundaries.values()) - allocation_boundaries[4]
        ),
    }


def load_and_validate_tables(
    table_path: Path,
    summary_path: Path,
    matrix: np.ndarray,
    dependency: np.ndarray,
    expected_labels: np.ndarray,
    expected_epsilons: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict]:
    summary = json.loads(summary_path.read_text())
    if (
        summary.get("experiment") != "p7_size8_one_elevation_tables"
        or summary.get("status")
        != "complete_exact_elevated_direction_omission_tables"
        or int(summary.get("projection_rows", 0)) != PROJECTION_ROWS
        or summary.get("output_sha256") != sha256(table_path)
    ):
        raise ValueError("table cache or summary has the wrong identity")
    with np.load(table_path, allow_pickle=False) as handle:
        arrays = {key: handle[key] for key in handle.files}
    expected_shapes = {
        "labels": (8, 49),
        "epsilons": (8,),
        "dependency": (135, 282),
        "selected_coefficients": (8, 22, 135),
        "projected_dependencies": (8, 22, 282),
        "base": (8, 22),
        "singleton": (8, 8, 128, 22),
        "variable": (8, 8, 128, 36, 22),
        "variable_count": (8, 8, 128),
    }
    if set(arrays) != set(expected_shapes):
        raise AssertionError("one-elevation cache array names changed")
    for key, shape in expected_shapes.items():
        if arrays[key].shape != shape:
            raise AssertionError(f"one-elevation cache shape changed for {key}")
        if summary["array_sha256"].get(key) != array_sha256(arrays[key]):
            raise AssertionError(f"one-elevation cache hash mismatch for {key}")
    if not np.array_equal(arrays["labels"], expected_labels):
        raise AssertionError("cached direction labels changed")
    if not np.array_equal(arrays["epsilons"], expected_epsilons):
        raise AssertionError("cached direction signs changed")
    if not np.array_equal(arrays["dependency"], dependency):
        raise AssertionError("cached full dependency basis changed")

    coefficients = arrays["selected_coefficients"].astype(np.int64)
    projected = arrays["projected_dependencies"].astype(np.int64)
    rebuilt = coefficients @ dependency.astype(np.int64) % MODULUS
    if not np.array_equal(rebuilt, projected):
        raise AssertionError("conditioned dependencies do not rebuild from the full basis")
    expected_base = (
        projected[:, :, :2] @ np.asarray([29, 1], dtype=np.int64) % MODULUS
    ).astype(np.uint8)
    if not np.array_equal(expected_base, arrays["base"]):
        raise AssertionError("conditioned base syndromes changed")
    for omitted in range(8):
        columns = slice(2 + 35 * omitted, 2 + 35 * (omitted + 1))
        if np.any(projected[omitted, :, columns]):
            raise AssertionError("conditioned projection touches its omitted direction")
        if np.any(projected[omitted] @ (matrix.astype(np.int64) % 7) % 7):
            raise AssertionError("conditioned projection is not left-null")
        if modular_rank(projected[omitted]) != PROJECTION_ROWS:
            raise AssertionError("conditioned projection lost rank")

    singleton = arrays["singleton"]
    variable_count = arrays["variable_count"]
    for omitted in range(8):
        for direction, eps in enumerate(expected_epsilons):
            for mask in range(128):
                if mask.bit_count() not in (0, 2, 4, 6):
                    continue
                expected_count = (
                    MAX_FLOOR_CATALOG_ROWS
                    if int(eps) == -1 and mask.bit_count() == 4
                    else 0
                )
                if int(variable_count[omitted, direction, mask]) != expected_count:
                    raise AssertionError("floor variable-catalog marker changed")
                has_singleton = not np.any(singleton[omitted, direction, mask] == 255)
                if has_singleton != (expected_count == 0):
                    raise AssertionError("floor singleton-catalog marker changed")
    return arrays, summary


CUDA_SOURCE = r'''
extern "C" __global__
void one_elevation_scan(
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
    unsigned long long* stratum,
    unsigned long long* projected_leaves,
    unsigned long long* projected_boundaries,
    unsigned long long* odd_histogram,
    unsigned long long* floor_pair_histogram,
    unsigned long long* survivor_ranks,
    unsigned char* survivor_directions,
    const unsigned long long survivor_capacity)
{
    const unsigned long long logical_thread =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned long long logical_threads =
        (unsigned long long)gridDim.x * blockDim.x;
    unsigned long long local_checked = 0;
    unsigned long long local_stratum = 0;

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
            const int cost = eps > 0
                ? (b == 0 ? 0 : 8)
                : ((b == 0 || b == 4) ? 14 : 6);
            if (eps < 0) type_minus += cost;
            else type_plus += cost;
        }
        int deficient_eps = 0;
        int floor_pair_index = -1;
        if (type_minus == 24 && type_plus == 32)
        {
            deficient_eps = -1;
            floor_pair_index = 0;
        }
        else if (type_minus == 32 && type_plus == 24)
        {
            deficient_eps = 1;
            floor_pair_index = 1;
        }
        else continue;

        ++local_stratum;
        atomicAdd(odd_histogram + odd_secants, 1ULL);
        atomicAdd(floor_pair_histogram + floor_pair_index, 1ULL);
        bool boundary_pass = false;

        #pragma unroll
        for (int omitted = 0; omitted < 8; ++omitted)
        {
            if ((int)epsilons[omitted] != deficient_eps) continue;
            int syndrome[22];
            #pragma unroll
            for (int row = 0; row < 22; ++row)
                syndrome[row] = (int)base[22 * omitted + row];
            int variable_direction = -1;
            int variable_mask = -1;
            bool table_valid = true;
            #pragma unroll
            for (int direction = 0; direction < 8; ++direction)
            {
                if (direction == omitted) continue;
                const int mask = masks[direction];
                const int key = (omitted * 8 + direction) * 128 + mask;
                const int count = (int)variable_count[key];
                if (count)
                {
                    if (variable_direction >= 0) table_valid = false;
                    variable_direction = direction;
                    variable_mask = mask;
                }
                else
                {
                    const int offset = key * 22;
                    #pragma unroll
                    for (int row = 0; row < 22; ++row)
                    {
                        const int value = (int)singleton[offset + row];
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
                for (int row = 0; row < 22; ++row)
                    if (syndrome[row] % 7 != 0) passing = false;
            }
            else
            {
                const int key =
                    (omitted * 8 + variable_direction) * 128 + variable_mask;
                const int offset = key * 36 * 22;
                for (int catalog_row = 0; catalog_row < 36 && !passing; ++catalog_row)
                {
                    bool match = true;
                    #pragma unroll
                    for (int row = 0; row < 22; ++row)
                    {
                        const int value =
                            (int)variable[offset + catalog_row * 22 + row];
                        if ((syndrome[row] + value) % 7 != 0) match = false;
                    }
                    passing = match;
                }
            }
            if (passing)
            {
                boundary_pass = true;
                const unsigned long long slot = atomicAdd(projected_leaves, 1ULL);
                if (slot < survivor_capacity)
                {
                    survivor_ranks[slot] = rank0;
                    survivor_directions[slot] = (unsigned char)omitted;
                }
            }
        }
        if (boundary_pass) atomicAdd(projected_boundaries, 1ULL);
    }
    if (local_checked) atomicAdd(checked, local_checked);
    if (local_stratum) atomicAdd(stratum, local_stratum);
}
'''


def gpu_scan(
    kernel,
    arrays_gpu: dict,
    choose_gpu,
    start_rank: int,
    stop_rank: int,
    blocks: int,
    threads: int,
    survivor_capacity: int,
) -> dict:
    import cupy as cp

    checked = cp.zeros(1, dtype=cp.uint64)
    stratum = cp.zeros(1, dtype=cp.uint64)
    passed_leaves = cp.zeros(1, dtype=cp.uint64)
    passed_boundaries = cp.zeros(1, dtype=cp.uint64)
    odd_histogram = cp.zeros(57, dtype=cp.uint64)
    floor_pairs = cp.zeros(2, dtype=cp.uint64)
    ranks = cp.zeros(survivor_capacity, dtype=cp.uint64)
    directions = cp.zeros(survivor_capacity, dtype=cp.uint8)
    started = time.perf_counter()
    kernel(
        (blocks,),
        (threads,),
        (
            arrays_gpu["labels"],
            arrays_gpu["epsilons"],
            choose_gpu,
            np.int32(choose_gpu.shape[1]),
            arrays_gpu["base"],
            arrays_gpu["singleton"],
            arrays_gpu["variable"],
            arrays_gpu["variable_count"],
            np.uint64(start_rank),
            np.uint64(stop_rank),
            checked,
            stratum,
            passed_leaves,
            passed_boundaries,
            odd_histogram,
            floor_pairs,
            ranks,
            directions,
            np.uint64(survivor_capacity),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.perf_counter() - started
    passed_count = int(cp.asnumpy(passed_leaves)[0])
    if passed_count > survivor_capacity:
        raise RuntimeError(
            f"projection survivor capacity exceeded: {passed_count}>{survivor_capacity}"
        )
    rank_values = cp.asnumpy(ranks[:passed_count])
    direction_values = cp.asnumpy(directions[:passed_count])
    pairs = sorted(
        (int(rank), int(direction))
        for rank, direction in zip(rank_values, direction_values)
    )
    histogram_host = cp.asnumpy(odd_histogram)
    floor_pair_host = cp.asnumpy(floor_pairs)
    return {
        "checked": int(cp.asnumpy(checked)[0]),
        "stratum": int(cp.asnumpy(stratum)[0]),
        "projected_leaves": passed_count,
        "projected_boundaries": int(cp.asnumpy(passed_boundaries)[0]),
        "odd_histogram": {
            int(index): int(value)
            for index, value in enumerate(histogram_host)
            if value
        },
        "floor_pair_histogram": {
            (24, 32): int(floor_pair_host[0]),
            (32, 24): int(floor_pair_host[1]),
        },
        "survivor_pairs": pairs,
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
    omitted: int,
    base: np.ndarray,
    singleton: np.ndarray,
    variable: np.ndarray,
    variable_count: np.ndarray,
) -> bool:
    syndrome = base[omitted].astype(np.int16).copy()
    variable_key = None
    for direction, mask in enumerate(masks):
        if direction == omitted:
            continue
        count = int(variable_count[omitted, direction, mask])
        if count:
            if variable_key is not None:
                raise AssertionError("four-allocation boundary has two floor variables")
            variable_key = (direction, mask, count)
        else:
            row = singleton[omitted, direction, mask]
            if np.any(row == 255):
                raise AssertionError("missing floor singleton contribution")
            syndrome += row
    if variable_key is None:
        return bool(np.all(syndrome % 7 == 0))
    direction, mask, count = variable_key
    values = variable[omitted, direction, mask, :count]
    return bool(
        np.any(np.all((values.astype(np.int16) + syndrome) % 7 == 0, axis=1))
    )


def cpu_prefix(
    stop_rank: int,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    arrays: dict[str, np.ndarray],
) -> dict:
    stratum = 0
    passed_boundaries = 0
    passed_pairs = []
    odd_histogram: Counter[int] = Counter()
    floor_pairs: Counter[tuple[int, int]] = Counter()
    for rank, boundary in enumerate(
        itertools.islice(itertools.combinations(range(49), 8), stop_rank)
    ):
        masks = masks_for_boundary(boundary, labels)
        profile = tuple(mask.bit_count() for mask in masks)
        floors = type_costs(profile, epsilons)
        if floors not in EXPECTED_FOUR_FLOOR_PAIRS:
            continue
        stratum += 1
        odd_histogram[sum(profile)] += 1
        floor_pairs[floors] += 1
        deficient_eps = -1 if floors == (24, 32) else 1
        boundary_pass = False
        for omitted, eps in enumerate(epsilons):
            if eps != deficient_eps:
                continue
            if projected_passes(
                masks,
                omitted,
                arrays["base"],
                arrays["singleton"],
                arrays["variable"],
                arrays["variable_count"],
            ):
                passed_pairs.append((rank, omitted))
                boundary_pass = True
        passed_boundaries += int(boundary_pass)
    return {
        "checked": stop_rank,
        "stratum": stratum,
        "projected_leaves": len(passed_pairs),
        "projected_boundaries": passed_boundaries,
        "odd_histogram": dict(sorted(odd_histogram.items())),
        "floor_pair_histogram": {
            key: floor_pairs[key] for key in sorted(EXPECTED_FOUR_FLOOR_PAIRS)
        },
        "survivor_pairs": passed_pairs,
    }


def contribution(
    dependency: np.ndarray,
    direction: int,
    mask: int,
    eps: int,
    mean: int,
) -> np.ndarray:
    odd_fibres = {value for value in range(7) if mask & (1 << value)}
    values = mapped_catalog(
        mask.bit_count(), int(eps == -1), mean, odd_fibres, None
    ).astype(np.int64)
    bad = 13 - values
    block = dependency[:, 2 + 35 * direction : 2 + 35 * (direction + 1)].astype(
        np.int64
    )
    return (block @ (bad.T % 7) % 7).astype(np.uint8)


def exact_catalog_matches(
    syndrome: np.ndarray,
    variables: list[np.ndarray],
) -> list[list[int]]:
    if len(variables) == 1:
        hits = np.flatnonzero(
            np.all((variables[0].astype(np.int16) + syndrome[:, None]) % 7 == 0, axis=0)
        )
        return [[int(index)] for index in hits]
    if len(variables) != 2:
        raise AssertionError(f"unexpected exact variable-catalog count {len(variables)}")
    lookup: dict[bytes, list[int]] = {}
    for second_index in range(variables[1].shape[1]):
        key = np.ascontiguousarray(variables[1][:, second_index]).tobytes()
        lookup.setdefault(key, []).append(second_index)
    matches = []
    for first_index in range(variables[0].shape[1]):
        target = (
            -syndrome - variables[0][:, first_index].astype(np.int16)
        ) % MODULUS
        for second_index in lookup.get(np.ascontiguousarray(target.astype(np.uint8)).tobytes(), []):
            matches.append([first_index, second_index])
    return matches


def exact_recheck(
    pairs: list[tuple[int, int]],
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    dependency: np.ndarray,
) -> tuple[list[tuple[int, int]], list[dict]]:
    base = (
        dependency[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % MODULUS
    ).astype(np.int16)
    cache: dict[tuple[int, int, int], np.ndarray] = {}
    survivors = []
    rows = []
    for rank, omitted in pairs:
        boundary = unrank_lex(rank)
        masks = masks_for_boundary(boundary, labels)
        means = [floor_for(eps, mask.bit_count()) for eps, mask in zip(epsilons, masks)]
        means[omitted] += 8
        syndrome = base.copy()
        variables = []
        variable_metadata = []
        for direction, (mask, eps, mean) in enumerate(zip(masks, epsilons, means)):
            key = (direction, mask, mean)
            if key not in cache:
                cache[key] = contribution(dependency, direction, mask, eps, mean)
            values = cache[key]
            if values.shape[1] == 1:
                syndrome = (syndrome + values[:, 0]) % MODULUS
            else:
                variables.append(values)
                variable_metadata.append(
                    {
                        "direction_index": direction,
                        "odd_fibre_mask": mask,
                        "scaled_mean": mean,
                        "catalog_rows": int(values.shape[1]),
                    }
                )
        if not 1 <= len(variables) <= 2:
            raise AssertionError("raised leaf has an unexpected variable-catalog shape")
        matches = exact_catalog_matches(syndrome, variables)
        row = {
            "rank": rank,
            "elevated_direction_index": omitted,
            "boundary_finite_field": list(boundary),
            "boundary_vertices": [value + 1 for value in boundary],
            "odd_fibre_masks": list(masks),
            "scaled_means_direction_order": means,
            "variable_catalogs": variable_metadata,
            "full_dependency_match_count": len(matches),
            "matching_catalog_rows": matches,
        }
        rows.append(row)
        if matches:
            survivors.append((rank, omitted))
    return survivors, rows


def run(args: argparse.Namespace) -> dict:
    import cupy as cp

    started = time.time()
    source = json.loads(args.source.read_text())
    labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)
    scope = expected_scope(source, epsilons)

    matrix, dependencies, linear_rows = linear_data((MODULUS,))
    dependency = dependencies[MODULUS].astype(np.uint8)
    if (
        matrix.shape != (282, 1225)
        or dependency.shape != (135, 282)
        or linear_rows
        != [
            {
                "modulus": 7,
                "rank": 147,
                "left_dependency_dimension": 135,
                "left_null_audit": True,
            }
        ]
    ):
        raise AssertionError("common mod-seven score system changed")
    if np.any(dependency.astype(np.int64) @ (matrix.astype(np.int64) % 7) % 7):
        raise AssertionError("full dependency basis is not left-null")
    arrays, table_summary = load_and_validate_tables(
        args.tables,
        args.table_summary,
        matrix,
        dependency,
        labels,
        epsilon_array,
    )

    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("CUDA is unavailable")
    kernel = cp.RawKernel(CUDA_SOURCE, "one_elevation_scan", options=("--std=c++11",))
    arrays_gpu = {
        key: cp.asarray(arrays[key])
        for key in ("labels", "epsilons", "base", "singleton", "variable", "variable_count")
    }
    choose_gpu = cp.asarray(choose_table(49, 8))

    prefix_gpu = gpu_scan(
        kernel,
        arrays_gpu,
        choose_gpu,
        0,
        args.verify_prefix,
        args.blocks,
        args.threads,
        args.survivor_capacity,
    )
    prefix_cpu = cpu_prefix(args.verify_prefix, labels, epsilons, arrays)
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
        arrays_gpu,
        choose_gpu,
        0,
        total,
        args.blocks,
        args.threads,
        args.survivor_capacity,
    )
    if full["checked"] != total:
        raise AssertionError("full CUDA rank interval was not exhausted")
    if full["stratum"] != scope["boundary_count"]:
        raise AssertionError("CUDA four-allocation count disagrees with source scope")
    if full["odd_histogram"] != scope["odd_secant_histogram"]:
        raise AssertionError("CUDA four-allocation odd histogram changed")
    if full["floor_pair_histogram"] != EXPECTED_FOUR_FLOOR_PAIRS:
        raise AssertionError("CUDA four-allocation floor-pair census changed")

    exact_survivors, exact_rows = exact_recheck(
        full["survivor_pairs"], labels, epsilons, dependency
    )
    properties = cp.cuda.runtime.getDeviceProperties(0)
    device = properties["name"]
    if isinstance(device, bytes):
        device = device.decode()
    return {
        "experiment": "p7_size8_one_elevation_gpu",
        "status": "complete_exact_four_allocation_boundary_exhaustion",
        "p": 7,
        "c_H": -1,
        "finite_boundary_size": 8,
        "source": str(args.source),
        "source_sha256": sha256(args.source),
        "source_scope": scope,
        "mathematical_reduction": {
            "deficient_type_floor_sum": 24,
            "saturated_type_floor_sum": 32,
            "exact_type_mean_sums": [32, 32],
            "allocations_per_boundary": 4,
            "each_allocation": "raise exactly one deficient-type direction by 8",
            "conditioned_dependency_dimension": 112,
            "selected_conditioned_dependencies": PROJECTION_ROWS,
            "elevated_direction_block_is_identically_zero": True,
            "maximum_non_elevated_variable_catalogs_per_leaf": 1,
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
        "conditioned_table_cache": {
            "path": str(args.tables),
            "sha256": sha256(args.tables),
            "summary": str(args.table_summary),
            "summary_sha256": sha256(args.table_summary),
            "projected_dependency_sha256": table_summary["array_sha256"][
                "projected_dependencies"
            ],
            "all_arrays_rehashed_and_algebraically_validated": True,
        },
        "verification": {
            "method": "independent itertools CPU prefix versus direct-rank CUDA",
            "prefix_checked": args.verify_prefix,
            "all_counts_histograms_floor_pairs_and_projected_pairs_match": True,
        },
        "device": str(device),
        "launch": {"blocks": args.blocks, "threads": args.threads},
        "all_boundaries": total,
        "checked_boundaries": full["checked"],
        "four_allocation_boundaries": full["stratum"],
        "four_allocation_leaves": 4 * full["stratum"],
        "four_allocation_odd_secant_histogram": {
            str(key): value for key, value in full["odd_histogram"].items()
        },
        "four_allocation_floor_pair_histogram": [
            {"type_floor_sums": list(key), "boundaries": value}
            for key, value in sorted(full["floor_pair_histogram"].items())
        ],
        "projected_dependency_survivor_leaves": full["projected_leaves"],
        "projected_dependency_survivor_boundaries": full["projected_boundaries"],
        "projected_survivor_rank_direction_pairs": [
            [rank, direction] for rank, direction in full["survivor_pairs"]
        ],
        "full_135_dependency_recheck_rows": exact_rows,
        "full_dependency_survivor_leaves": len(exact_survivors),
        "full_dependency_survivor_rank_direction_pairs": [
            [rank, direction] for rank, direction in exact_survivors
        ],
        "all_four_allocation_leaves_mod7_infeasible": not exact_survivors,
        "proved_four_allocation_stratum_exclusion_cminus1": not exact_survivors,
        "remaining_nonconic_floor_survivors_if_excluded": (
            scope["remaining_boundaries_after_stratum_if_excluded"]
            if not exact_survivors
            else scope["remaining_boundaries_before_this_stratum"]
        ),
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "elapsed_seconds": time.time() - started,
        "cuda_scan_seconds": full["elapsed_seconds"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--table-summary", type=Path, required=True)
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
            "projected_survivor_rank_direction_pairs",
            "full_135_dependency_recheck_rows",
            "full_dependency_survivor_rank_direction_pairs",
        }
    }
    print(json.dumps(compact, indent=2), flush=True)


if __name__ == "__main__":
    main()
