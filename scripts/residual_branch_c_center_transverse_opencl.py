#!/usr/bin/env python3
"""OpenCL coordinate search for the exact p=31 centre-transverse fibre.

This is an independent backend for
``residual_branch_c_center_transverse_gpu.py``.  It deliberately imports the
canonical table builder and exact graph replay instead of reimplementing the
mathematics.  Each OpenCL workgroup evaluates one (chain, centre-option)
candidate with one work-item per Radon row, using integer arithmetic
throughout.

The search is bounded and explores only the frozen one-cancellation
half/auxiliary fibre.  A successful compact-l1 test would still be only a
necessary condition for the common residual graph, and a failed bounded run
is not an exclusion certificate.
"""

from __future__ import annotations

import argparse
import json
import random
import socket
import time
from pathlib import Path

import numpy as np

from scripts.residual_branch_c_center_transverse_gpu import (
    HALVES,
    PAIR_FIRST,
    PAIR_FIRST_CENTER,
    PAIR_SECOND,
    PAIR_SECOND_CENTER,
    build_exact_tables,
    exact_graph_replay,
    validate_published_witness,
)
from io_atomic import write_json_atomic


ROWS = 32
CELLS = 465
CENTRE_OPTIONS = 30

V100_REFERENCE_CENTERS = (8, 9, 28, 12, 7, 9, 11, 6, 3, 15, 12, 3, 15, 1, 10, 11)
V100_REFERENCE_FIXED_INDEX = 2
V100_REFERENCE_TOTAL = 4604
V100_REFERENCE_MAXIMUM = 178
V100_REFERENCE_GRAPH_SHA256 = (
    "46b4a125eee33a6b66a358b7a04bf9cb640a9f1398205e345c882f2d10a02889"
)


KERNEL_SOURCE = r"""
#define ROWS 32
#define CELLS 465

__kernel void score_options(
    __global const short *state,
    __global const short *option_tables,
    const int table_base,
    __global const short *old_choices,
    __global const short *budgets,
    const int option_count,
    __global int *out_total,
    __global int *out_maximum)
{
    const int group = get_group_id(0);
    const int row = get_local_id(0);
    const int chain = group / option_count;
    const int option = group - chain * option_count;
    const int old_option = (int)old_choices[chain];
    const int state_offset = (chain * ROWS + row) * CELLS;
    const int old_offset = table_base + (old_option * ROWS + row) * CELLS;
    const int new_offset = table_base + (option * ROWS + row) * CELLS;
    int l1 = 0;
    for (int cell = 0; cell < CELLS; ++cell) {
        int value = (int)state[state_offset + cell]
                  - (int)option_tables[old_offset + cell]
                  + (int)option_tables[new_offset + cell];
        l1 += abs(value);
    }
    int excess = l1 - (int)budgets[row];
    excess = excess > 0 ? excess : 0;
    __local int row_excess[ROWS];
    row_excess[row] = excess;
    barrier(CLK_LOCAL_MEM_FENCE);
    if (row == 0) {
        int total = 0;
        int maximum = 0;
        for (int index = 0; index < ROWS; ++index) {
            const int value = row_excess[index];
            total += value;
            maximum = value > maximum ? value : maximum;
        }
        out_total[group] = total;
        out_maximum[group] = maximum;
    }
}

__kernel void apply_options(
    __global short *state,
    __global const short *option_tables,
    const int table_base,
    __global const short *old_choices,
    __global const short *new_choices,
    const int element_count)
{
    const int index = get_global_id(0);
    if (index >= element_count)
        return;
    const int per_chain = ROWS * CELLS;
    const int chain = index / per_chain;
    const int inner = index - chain * per_chain;
    const int old_offset = table_base + (int)old_choices[chain] * per_chain + inner;
    const int new_offset = table_base + (int)new_choices[chain] * per_chain + inner;
    state[index] = (short)((int)state[index]
                         - (int)option_tables[old_offset]
                         + (int)option_tables[new_offset]);
}
"""


def _reference_cpu_replays(tables: dict[str, object]) -> dict[str, object]:
    baseline = validate_published_witness(tables)
    v100 = exact_graph_replay(
        list(V100_REFERENCE_CENTERS), V100_REFERENCE_FIXED_INDEX, tables
    )
    if (
        baseline["total_positive_l1_excess"] != 5068
        or baseline["maximum_l1_excess"] != 194
    ):
        raise ArithmeticError("published 5068/194 CPU reference changed")
    if (
        v100["total_positive_l1_excess"] != V100_REFERENCE_TOTAL
        or v100["maximum_l1_excess"] != V100_REFERENCE_MAXIMUM
        or v100["graph_sha256"] != V100_REFERENCE_GRAPH_SHA256
    ):
        raise ArithmeticError("recorded V100 winner failed exact CPU replay")
    return {"published_5068_194": baseline, "v100_4604_178": v100}


def _select_device():
    import pyopencl as cl

    devices = [
        device
        for platform in cl.get_platforms()
        for device in platform.get_devices(device_type=cl.device_type.GPU)
    ]
    if not devices:
        raise RuntimeError("no OpenCL GPU device found")
    device = max(devices, key=lambda item: (item.max_compute_units, item.global_mem_size))
    return cl, device


def _valid_options(
    choices: np.ndarray, variable: int, pair_bad: np.ndarray
) -> np.ndarray:
    valid = np.ones((len(choices), CENTRE_OPTIONS), dtype=np.bool_)
    for other in range(choices.shape[1]):
        if other == variable:
            continue
        for option in range(CENTRE_OPTIONS):
            valid[:, option] &= ~pair_bad[
                variable, other, option, choices[:, other]
            ]
    return valid


def _valid_assignment(
    choices: np.ndarray,
    pair_bad: np.ndarray,
    ignored_pair: tuple[int, int],
) -> bool:
    ignored = frozenset(ignored_pair)
    for first in range(len(choices)):
        for second in range(first + 1, len(choices)):
            if frozenset((first, second)) == ignored:
                continue
            if pair_bad[first, second, choices[first], choices[second]]:
                return False
    return True


def _random_design_assignments(
    count: int,
    pair_bad: np.ndarray,
    seed: int,
    frozen_pair: tuple[int, int],
    frozen_centers: tuple[int, int],
    initial_centers: tuple[int, ...] | None,
) -> np.ndarray:
    """Generate collision-clean centre tuples for an arbitrary exact design."""
    rng = random.Random(seed)
    half_count = pair_bad.shape[0]
    frozen = {
        int(frozen_pair[0]): int(frozen_centers[0]) - 1,
        int(frozen_pair[1]): int(frozen_centers[1]) - 1,
    }
    out: list[np.ndarray] = []
    if initial_centers is not None:
        baseline = np.asarray(
            [int(center) - 1 for center in initial_centers], dtype=np.int16
        )
        if len(baseline) != half_count or not _valid_assignment(
            baseline, pair_bad, frozen_pair
        ):
            raise ArithmeticError("the design's initial center tuple is invalid")
        if any(baseline[index] != value for index, value in frozen.items()):
            raise ArithmeticError("the initial tuple changed a frozen collision center")
        out.append(baseline)
    attempts = 0
    while len(out) < count:
        attempts += 1
        if attempts > max(10_000, 1_000 * count):
            raise ArithmeticError("could not construct enough clean center assignments")
        choices = np.full(half_count, -1, dtype=np.int16)
        for index, value in frozen.items():
            choices[index] = value
        order = [index for index in range(half_count) if index not in frozen]
        rng.shuffle(order)
        for index in order:
            options = list(range(CENTRE_OPTIONS))
            rng.shuffle(options)
            selected = next(
                (
                    option
                    for option in options
                    if all(
                        other == index
                        or choices[other] < 0
                        or not pair_bad[index, other, option, choices[other]]
                        for other in range(half_count)
                    )
                ),
                None,
            )
            if selected is None:
                break
            choices[index] = selected
        if np.all(choices >= 0) and _valid_assignment(
            choices, pair_bad, frozen_pair
        ):
            out.append(choices)
    return np.stack(out)


def _lexicographic_choice(
    totals: np.ndarray,
    maxima: np.ndarray,
    valid: np.ndarray,
    primary: str,
) -> np.ndarray:
    chains, options = totals.shape
    limit = np.iinfo(np.int32).max
    best_first = np.full(chains, limit, dtype=np.int64)
    best_second = np.full(chains, limit, dtype=np.int64)
    best_option = np.zeros(chains, dtype=np.int16)
    first_values, second_values = (
        (totals, maxima) if primary == "total" else (maxima, totals)
    )
    for option in range(options):
        first = first_values[:, option].astype(np.int64, copy=False)
        second = second_values[:, option].astype(np.int64, copy=False)
        improve = valid[:, option] & (
            (first < best_first)
            | ((first == best_first) & (second < best_second))
        )
        best_first[improve] = first[improve]
        best_second[improve] = second[improve]
        best_option[improve] = option
    if np.any(best_first == limit):
        raise ArithmeticError("a chain lost every coordinate option")
    return best_option


def _initial_state(
    choices: np.ndarray,
    fixed_choices: np.ndarray,
    half_tables: np.ndarray,
    fixed_tables: np.ndarray,
    collision_table: np.ndarray,
) -> np.ndarray:
    state = fixed_tables[fixed_choices].copy()
    state -= collision_table[None, :, :]
    for half_index in range(choices.shape[1]):
        state += half_tables[
            half_index, choices[:, half_index], :, :
        ]
    return np.ascontiguousarray(state, dtype=np.int16)


def opencl_coordinate_search(
    tables: dict[str, object],
    chains: int,
    restarts: int,
    sweeps: int,
    seed: int,
    primary: str,
    time_limit: float | None,
    frozen_pair: tuple[int, int] = (PAIR_FIRST, PAIR_SECOND),
    frozen_centers: tuple[int, int] = (
        PAIR_FIRST_CENTER,
        PAIR_SECOND_CENTER,
    ),
    initial_centers: tuple[int, ...] | None = None,
    initial_fixed_index: int | None = None,
) -> dict[str, object]:
    cl, device = _select_device()
    context = cl.Context([device])
    queue = cl.CommandQueue(context)
    program = cl.Program(context, KERNEL_SOURCE).build()
    score_kernel = cl.Kernel(program, "score_options")
    apply_kernel = cl.Kernel(program, "apply_options")
    memory_flags = cl.mem_flags

    half_tables = np.ascontiguousarray(tables["half_tables"], dtype=np.int16)
    fixed_tables = np.ascontiguousarray(tables["fixed_tables"], dtype=np.int16)
    collision_table = np.asarray(tables["collision_table"], dtype=np.int16)
    pair_bad = np.asarray(tables["pair_bad"], dtype=np.bool_)
    budgets = np.ascontiguousarray(tables["budgets"], dtype=np.int16)
    fixed_option_count = len(fixed_tables)

    half_buffer = cl.Buffer(
        context, memory_flags.READ_ONLY | memory_flags.COPY_HOST_PTR, hostbuf=half_tables
    )
    fixed_buffer = cl.Buffer(
        context, memory_flags.READ_ONLY | memory_flags.COPY_HOST_PTR, hostbuf=fixed_tables
    )
    budget_buffer = cl.Buffer(
        context, memory_flags.READ_ONLY | memory_flags.COPY_HOST_PTR, hostbuf=budgets
    )
    state_buffer = cl.Buffer(
        context, memory_flags.READ_WRITE, chains * ROWS * CELLS * np.dtype(np.int16).itemsize
    )
    old_buffer = cl.Buffer(
        context, memory_flags.READ_ONLY, chains * np.dtype(np.int16).itemsize
    )
    new_buffer = cl.Buffer(
        context, memory_flags.READ_ONLY, chains * np.dtype(np.int16).itemsize
    )
    output_count = chains * CENTRE_OPTIONS
    total_buffer = cl.Buffer(
        context, memory_flags.WRITE_ONLY, output_count * np.dtype(np.int32).itemsize
    )
    maximum_buffer = cl.Buffer(
        context, memory_flags.WRITE_ONLY, output_count * np.dtype(np.int32).itemsize
    )
    totals_flat = np.empty(output_count, dtype=np.int32)
    maxima_flat = np.empty(output_count, dtype=np.int32)

    rng = random.Random(seed)
    half_count = half_tables.shape[0]
    if pair_bad.shape[:2] != (half_count, half_count):
        raise ArithmeticError("the pair-intersection table has the wrong shape")
    mutable = [
        index for index in range(half_count)
        if index not in frozen_pair
    ]
    started = time.perf_counter()
    evaluated = 0
    completed = 0
    best: dict[str, object] | None = None
    stop = False

    def evaluate_and_apply(
        table_buffer,
        table_base: int,
        option_count: int,
        old_choices: np.ndarray,
        valid: np.ndarray,
    ) -> np.ndarray:
        nonlocal evaluated
        old_choices = np.ascontiguousarray(old_choices, dtype=np.int16)
        cl.enqueue_copy(queue, old_buffer, old_choices)
        score_kernel(
            queue,
            (chains * option_count * ROWS,),
            (ROWS,),
            state_buffer,
            table_buffer,
            np.int32(table_base),
            old_buffer,
            budget_buffer,
            np.int32(option_count),
            total_buffer,
            maximum_buffer,
        )
        cl.enqueue_copy(queue, totals_flat, total_buffer)
        cl.enqueue_copy(queue, maxima_flat, maximum_buffer).wait()
        totals = totals_flat[: chains * option_count].reshape(chains, option_count)
        maxima = maxima_flat[: chains * option_count].reshape(chains, option_count)
        selected = _lexicographic_choice(totals, maxima, valid, primary)
        selected = np.ascontiguousarray(selected, dtype=np.int16)
        cl.enqueue_copy(queue, new_buffer, selected)
        element_count = chains * ROWS * CELLS
        global_size = ((element_count + 255) // 256) * 256
        apply_kernel(
            queue,
            (global_size,),
            (256,),
            state_buffer,
            table_buffer,
            np.int32(table_base),
            old_buffer,
            new_buffer,
            np.int32(element_count),
        )
        evaluated += int(valid.sum())
        return selected

    for restart in range(restarts):
        if time_limit is not None and restart and time.perf_counter() - started >= time_limit:
            break
        choices = _random_design_assignments(
            chains,
            pair_bad,
            rng.randrange(1 << 62),
            frozen_pair,
            frozen_centers,
            initial_centers,
        )
        fixed_choices = np.asarray(
            [rng.randrange(fixed_option_count) for _ in range(chains)], dtype=np.int16
        )
        if restart == 0 and initial_fixed_index is not None:
            fixed_choices[0] = initial_fixed_index
        state = _initial_state(
            choices, fixed_choices, half_tables, fixed_tables, collision_table
        )
        cl.enqueue_copy(queue, state_buffer, state)

        for _sweep in range(sweeps):
            changed = False
            rng.shuffle(mutable)
            for variable in mutable:
                valid = _valid_options(choices, variable, pair_bad)
                selected = evaluate_and_apply(
                    half_buffer,
                    variable * CENTRE_OPTIONS * ROWS * CELLS,
                    CENTRE_OPTIONS,
                    choices[:, variable],
                    valid,
                )
                changed |= bool(np.any(selected != choices[:, variable]))
                choices[:, variable] = selected
                if time_limit is not None and time.perf_counter() - started >= time_limit:
                    stop = True
                    break
            if stop:
                break
            valid_fixed = np.ones((chains, fixed_option_count), dtype=np.bool_)
            selected_fixed = evaluate_and_apply(
                fixed_buffer, 0, fixed_option_count, fixed_choices, valid_fixed
            )
            changed |= bool(np.any(selected_fixed != fixed_choices))
            fixed_choices = selected_fixed
            if not changed:
                break

        cl.enqueue_copy(queue, state, state_buffer).wait()
        l1 = np.abs(state.astype(np.int32)).sum(axis=2)
        excess = np.maximum(l1 - budgets[None, :], 0)
        totals = excess.sum(axis=1)
        maxima = excess.max(axis=1)
        winner = min(
            range(chains),
            key=lambda index: (
                (int(totals[index]), int(maxima[index]))
                if primary == "total"
                else (int(maxima[index]), int(totals[index]))
            ),
        )
        record = {
            "centers": [int(value) + 1 for value in choices[winner]],
            "fixed_edge_index": int(fixed_choices[winner]),
            "total_positive_l1_excess": int(totals[winner]),
            "maximum_l1_excess": int(maxima[winner]),
            "row_l1": [int(value) for value in l1[winner]],
            "restart": restart,
        }
        record_key = (
            (record["total_positive_l1_excess"], record["maximum_l1_excess"])
            if primary == "total"
            else (record["maximum_l1_excess"], record["total_positive_l1_excess"])
        )
        best_key = None if best is None else (
            (best["total_positive_l1_excess"], best["maximum_l1_excess"])
            if primary == "total"
            else (best["maximum_l1_excess"], best["total_positive_l1_excess"])
        )
        if best_key is None or record_key < best_key:
            best = record
            print(json.dumps({"new_best": best}), flush=True)
        completed += 1
        if record["total_positive_l1_excess"] == 0 or stop:
            break

    queue.finish()
    elapsed = time.perf_counter() - started
    if best is None:
        raise ArithmeticError("the OpenCL search completed no restart")
    return {
        "opencl_platform": device.platform.name,
        "opencl_device": device.name,
        "opencl_device_version": device.version,
        "global_memory_bytes": int(device.global_mem_size),
        "compute_units": int(device.max_compute_units),
        "elapsed_seconds": elapsed,
        "objective_states_evaluated": evaluated,
        "objective_states_per_second": evaluated / elapsed,
        "chains": chains,
        "restarts_requested": restarts,
        "restarts_completed": completed,
        "maximum_sweeps": sweeps,
        "primary_objective": primary,
        "time_limit_seconds": time_limit,
        "frozen_collision_halves": list(frozen_pair),
        "frozen_collision_centers": list(frozen_centers),
        "best": best,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    tables = build_exact_tables()
    references = _reference_cpu_replays(tables)
    published = references["published_5068_194"]
    search = opencl_coordinate_search(
        tables,
        chains=args.chains,
        restarts=args.restarts,
        sweeps=args.sweeps,
        seed=args.seed,
        primary=args.primary,
        time_limit=args.time_limit,
        initial_centers=tuple(int(value) for value in published["centers"]),
        initial_fixed_index=int(published["fixed_edge_index"]),
    )
    best = search["best"]
    replay = exact_graph_replay(
        list(best["centers"]), int(best["fixed_edge_index"]), tables
    )
    for key in (
        "row_l1",
        "total_positive_l1_excess",
        "maximum_l1_excess",
    ):
        if replay[key] != best[key]:
            raise ArithmeticError(f"OpenCL candidate failed exact CPU replay at {key}")
    result = {
        "schema": "residual_branch_c_center_transverse_opencl_v1",
        "classification": "bounded OpenCL search in one fixed p31 half/auxiliary fibre",
        "host": socket.gethostname(),
        "reference_cpu_replays": references,
        "search": search,
        "best_exact_replay": replay,
        "l1_pass_is_sufficient_for_atom_decomposition": False,
        "fibre_exhausted": False,
        "residual_ii_closed": False,
    }
    if args.output:
        write_json_atomic(args.output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chains", type=int, default=256)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--sweeps", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--primary", choices=("total", "max"), default="total")
    parser.add_argument("--time-limit", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.chains < 1 or args.restarts < 1 or args.sweeps < 1:
        parser.error("chains, restarts, and sweeps must be positive")
    if args.time_limit is not None and args.time_limit <= 0:
        parser.error("time limit must be positive")
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
