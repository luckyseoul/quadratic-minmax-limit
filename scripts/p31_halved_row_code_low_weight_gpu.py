#!/usr/bin/env python3
"""V100 heuristic for a p=31 halved-row-code word below weight 465.

This is a single-prime discovery search, not theorem evidence and not a prime
sweep.  It uses the exact normal form

    W = 1 q^T + B X B^T,

where the 480 columns of ``B`` are the paired affine-block basis and ``X`` is
block diagonal with 32 blocks of size 15.  For fixed ``X`` the boundary word
``q`` is optimized exactly, independently in every physical column.  The GPU
therefore minimizes

    sum_delta min(w_delta, 480-w_delta),

where ``w_delta`` is the column weight of ``B X B^T``.

Each CUDA chain starts at a known weight-465 fixed-transverse rectangle.  The
default run diversifies the chains with one through eight forced moves and
then performs deterministic greedy descent; multi-round temperature schedules
remain available.  The two exact code-preserving moves are one coefficient
flip and a 2-by-2 rank-one block flip.  Physical columns are stored as eight
uint64 words.  Every reported best coefficient vector is reconstructed from
scratch on the CPU; a sub-465 result is emitted only after the CPU weight
agrees.

If CuPy is unavailable, the script still emits a machine-readable backend
record and names the needed CUDA-matched package.  It never installs packages.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_symmetric_halved_mod2 import _antipodal_classes  # noqa: E402
from e1_gmin_m4_symmetric_halved_row_code import (  # noqa: E402
    _affine_blocks,
    halved_row_code_decomposition,
)


P = 31
H = 15
D = 32
N = 480
WORDS = 8
VARIABLES = D * H * H
THRESHOLD = P * H
DEFAULT_SEED = 0x31D15761


CUDA_SOURCE = r'''
extern "C" {

__device__ __forceinline__ unsigned int xs32(unsigned int &x) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

__device__ __forceinline__ int folded_weight(int weight, int n) {
    const int complement = n - weight;
    return weight < complement ? weight : complement;
}

__global__ void init_rectangle_chains(
    unsigned long long *state,
    unsigned char *coeff,
    unsigned short *column_weight,
    int *objective,
    int *coeff_weight,
    unsigned int *rng_state,
    int *best_objective,
    unsigned char *best_coeff,
    const unsigned long long *star_masks,
    const int *block_points,
    const int chains,
    const int n,
    const int words,
    const int p,
    const int h,
    const int d,
    const int variables,
    const unsigned int seed,
    const int restart,
    const int reset_best)
{
    const int chain = blockDim.x * blockIdx.x + threadIdx.x;
    if (chain >= chains) return;

    unsigned int rng = seed
        ^ (0x9e3779b9U * (unsigned int)(chain + 1))
        ^ (0x85ebca6bU * (unsigned int)(restart + 1));
    xs32(rng);
    const int direction = (int)(xs32(rng) % (unsigned int)d);
    const int right_block = (int)(xs32(rng) % (unsigned int)h);

    const long long state_base = (long long)chain * n * words;
    const long long coeff_base = (long long)chain * variables;
    const long long column_base = (long long)chain * n;
    for (int index = 0; index < n * words; ++index)
        state[state_base + index] = 0ULL;
    for (int index = 0; index < variables; ++index)
        coeff[coeff_base + index] = 0U;
    for (int column = 0; column < n; ++column)
        column_weight[column_base + column] = 0U;

    const int global_right = direction * h + right_block;
    for (int k = 0; k < p; ++k) {
        const int column = block_points[global_right * p + k];
        for (int word = 0; word < words; ++word)
            state[state_base + column * words + word]
                = star_masks[direction * words + word];
        column_weight[column_base + column] = (unsigned short)(n - h);
    }
    for (int left_block = 0; left_block < h; ++left_block) {
        const int variable = direction * h * h
            + left_block * h + right_block;
        coeff[coeff_base + variable] = 1U;
    }

    objective[chain] = p * h;
    coeff_weight[chain] = h;
    rng_state[chain] = rng == 0U ? 0xa341316cU : rng;
    if (reset_best) {
        best_objective[chain] = p * h;
        const long long best_base = (long long)chain * variables;
        for (int index = 0; index < variables; ++index)
            best_coeff[best_base + index] = coeff[coeff_base + index];
    }
}

__global__ void anneal_block_code(
    unsigned long long *state,
    unsigned char *coeff,
    unsigned short *column_weight,
    int *objective,
    int *coeff_weight,
    unsigned int *rng_state,
    int *best_objective,
    unsigned char *best_coeff,
    const unsigned long long *block_masks,
    const int *block_points,
    const int chains,
    const int n,
    const int words,
    const int p,
    const int h,
    const int d,
    const int variables,
    const int steps,
    const int kick_steps,
    const int move2_percent,
    const float temperature)
{
    const int chain = blockDim.x * blockIdx.x + threadIdx.x;
    if (chain >= chains) return;
    const long long state_base = (long long)chain * n * words;
    const long long coeff_base = (long long)chain * variables;
    const long long column_base = (long long)chain * n;
    unsigned int rng = rng_state[chain];
    int obj = objective[chain];
    int coefficient_count = coeff_weight[chain];
    const int chain_kick_steps = kick_steps > 0
        ? 1 + chain % kick_steps
        : 0;

    for (int step = 0; step < steps; ++step) {
        const bool move2 =
            (int)(xs32(rng) % 100U) < move2_percent;
        const int direction = (int)(xs32(rng) % (unsigned int)d);
        const int first_left = (int)(xs32(rng) % (unsigned int)h);
        const int first_right = (int)(xs32(rng) % (unsigned int)h);
        int second_left = first_left;
        int second_right = first_right;
        int row_count = p;
        int column_blocks = 1;
        if (move2) {
            second_left = (int)(xs32(rng) % (unsigned int)(h - 1));
            if (second_left >= first_left) ++second_left;
            second_right = (int)(xs32(rng) % (unsigned int)(h - 1));
            if (second_right >= first_right) ++second_right;
            row_count = 2 * p;
            column_blocks = 2;
        }

        const int first_global_left = direction * h + first_left;
        const int second_global_left = direction * h + second_left;
        int delta_objective = 0;
        for (int side = 0; side < column_blocks; ++side) {
            const int right = side == 0 ? first_right : second_right;
            const int global_right = direction * h + right;
            for (int k = 0; k < p; ++k) {
                const int column = block_points[global_right * p + k];
                int intersection = 0;
                for (int word = 0; word < words; ++word) {
                    unsigned long long row_mask
                        = block_masks[first_global_left * words + word];
                    if (move2)
                        row_mask ^= block_masks[
                            second_global_left * words + word];
                    intersection += __popcll(
                        state[state_base + column * words + word] & row_mask);
                }
                const int old_weight =
                    (int)column_weight[column_base + column];
                const int new_weight = old_weight + row_count - 2 * intersection;
                delta_objective += folded_weight(new_weight, n)
                    - folded_weight(old_weight, n);
            }
        }

        const int proposed_objective = obj + delta_objective;
        bool accept = proposed_objective > 0 && step < chain_kick_steps;
        if (proposed_objective > 0 && !accept) {
            if (delta_objective <= 0) {
                accept = true;
            } else if (temperature > 0.0f) {
                const float uniform =
                    ((float)xs32(rng) + 1.0f) * 2.3283064365386963e-10f;
                accept = uniform < __expf(
                    -(float)delta_objective / temperature);
            }
        }
        if (!accept) continue;

        for (int side = 0; side < column_blocks; ++side) {
            const int right = side == 0 ? first_right : second_right;
            const int global_right = direction * h + right;
            for (int k = 0; k < p; ++k) {
                const int column = block_points[global_right * p + k];
                int intersection = 0;
                for (int word = 0; word < words; ++word) {
                    unsigned long long row_mask
                        = block_masks[first_global_left * words + word];
                    if (move2)
                        row_mask ^= block_masks[
                            second_global_left * words + word];
                    intersection += __popcll(
                        state[state_base + column * words + word] & row_mask);
                }
                const int old_weight =
                    (int)column_weight[column_base + column];
                column_weight[column_base + column] = (unsigned short)(
                    old_weight + row_count - 2 * intersection);
                for (int word = 0; word < words; ++word) {
                    unsigned long long row_mask
                        = block_masks[first_global_left * words + word];
                    if (move2)
                        row_mask ^= block_masks[
                            second_global_left * words + word];
                    state[state_base + column * words + word] ^= row_mask;
                }
            }
        }

        const int first_variable = direction * h * h
            + first_left * h + first_right;
        coefficient_count += coeff[coeff_base + first_variable] ? -1 : 1;
        coeff[coeff_base + first_variable] ^= 1U;
        if (move2) {
            const int variables4[3] = {
                direction * h * h + first_left * h + second_right,
                direction * h * h + second_left * h + first_right,
                direction * h * h + second_left * h + second_right
            };
            for (int k = 0; k < 3; ++k) {
                const int variable = variables4[k];
                coefficient_count += coeff[coeff_base + variable] ? -1 : 1;
                coeff[coeff_base + variable] ^= 1U;
            }
        }
        obj = proposed_objective;
        if (obj < best_objective[chain]) {
            best_objective[chain] = obj;
            const long long best_base = (long long)chain * variables;
            for (int index = 0; index < variables; ++index)
                best_coeff[best_base + index] = coeff[coeff_base + index];
        }
    }
    objective[chain] = obj;
    coeff_weight[chain] = coefficient_count;
    rng_state[chain] = rng;
}

}
'''


def build_geometry() -> dict[str, Any]:
    """Build the p=31 affine-block basis used by the exact normal form."""
    theorem = halved_row_code_decomposition(P)
    if (
        theorem["delta_size"] != N
        or theorem["h"] != H
        or theorem["d"] != D
        or theorem["total_rank"] != D * H * (H + 1)
    ):
        raise ArithmeticError("the imported p=31 row-code normal form changed")

    classes = _antipodal_classes(P)
    class_index = {point: index for index, point in enumerate(classes)}
    labelled_blocks = _affine_blocks(P)
    if len(classes) != N or len(labelled_blocks) != N:
        raise ArithmeticError("the p=31 point/block count changed")
    block_points = np.array(
        [
            [class_index[point] for point in block]
            for _direction, block in labelled_blocks
        ],
        dtype=np.int32,
    )
    if block_points.shape != (N, P):
        raise ArithmeticError("the affine-block size changed")
    for direction in range(D):
        labels = {
            labelled_blocks[direction * H + index][0]
            for index in range(H)
        }
        if len(labels) != 1:
            raise ArithmeticError("affine blocks are not grouped by direction")
        if len(set(block_points[direction * H].tolist())) != P:
            raise ArithmeticError("an affine block has repeated points")

    block_masks = np.zeros((N, WORDS), dtype=np.uint64)
    for block, points in enumerate(block_points):
        for point in points:
            word, bit = divmod(int(point), 64)
            block_masks[block, word] |= np.uint64(1 << bit)
    if not np.all(
        np.array(
            [sum(int(word).bit_count() for word in mask) for mask in block_masks]
        )
        == P
    ):
        raise ArithmeticError("a packed affine block has the wrong weight")

    star_masks = np.zeros((D, WORDS), dtype=np.uint64)
    for direction in range(D):
        for block in range(H):
            star_masks[direction] ^= block_masks[direction * H + block]
    star_weights = np.array(
        [sum(int(word).bit_count() for word in mask) for mask in star_masks]
    )
    if not np.all(star_weights == N - H):
        raise ArithmeticError("the direction stars have the wrong weight")

    return {
        "theorem": theorem,
        "classes": classes,
        "block_points": np.ascontiguousarray(block_points),
        "block_masks": np.ascontiguousarray(block_masks),
        "star_masks": np.ascontiguousarray(star_masks),
    }


def cpu_reconstruct(
    coefficients: np.ndarray,
    geometry: dict[str, Any],
    reported_weight: int | None = None,
) -> dict[str, Any]:
    """Reconstruct and independently verify one normal-form word on CPU."""
    coefficient_bits = np.asarray(coefficients, dtype=np.uint8).reshape(-1)
    if coefficient_bits.shape != (VARIABLES,) or not np.all(
        np.isin(coefficient_bits, (0, 1))
    ):
        raise ValueError("coefficient vector is not binary of length 7200")

    block_points = geometry["block_points"]
    physical = np.zeros((N, N), dtype=np.bool_)
    support = np.flatnonzero(coefficient_bits)
    triples: list[list[int]] = []
    for variable_raw in support:
        variable = int(variable_raw)
        direction, remainder = divmod(variable, H * H)
        left, right = divmod(remainder, H)
        rows = block_points[direction * H + left]
        columns = block_points[direction * H + right]
        physical[np.ix_(rows, columns)] ^= True
        triples.append([direction, left, right])

    column_weights = np.count_nonzero(physical, axis=0)
    q = column_weights > N // 2
    word = np.logical_xor(physical, q[np.newaxis, :])
    weight = int(np.count_nonzero(word))
    folded_weight = int(np.minimum(column_weights, N - column_weights).sum())
    if weight != folded_weight:
        raise ArithmeticError("CPU boundary optimization changed the weight")
    if reported_weight is not None and weight != reported_weight:
        raise ArithmeticError(
            f"CPU/GPU best-weight mismatch: {weight} != {reported_weight}"
        )
    if weight == 0:
        raise ArithmeticError("the heuristic returned the zero word")

    packed_coefficients = np.packbits(coefficient_bits, bitorder="little")
    packed_word = np.packbits(word.reshape(-1), bitorder="little")
    classes = geometry["classes"]
    q_support = np.flatnonzero(q).tolist()
    return {
        "cpu_verified": True,
        "weight": weight,
        "below_ph_465": weight < THRESHOLD,
        "coefficient_weight": int(len(support)),
        "coefficient_support_direction_left_right": triples,
        "optimized_boundary_q_support_indices": q_support,
        "optimized_boundary_q_support_classes": [
            list(classes[index]) for index in q_support
        ],
        "coefficient_bits_sha256": hashlib.sha256(
            packed_coefficients.tobytes()
        ).hexdigest(),
        "physical_word_bits_sha256": hashlib.sha256(
            packed_word.tobytes()
        ).hexdigest(),
        "normal_form_membership_reason": (
            "coefficients use only B_A tensor B_A and q is a boundary word"
        ),
    }


def baseline_certificate(geometry: dict[str, Any]) -> dict[str, Any]:
    """CPU-check the fixed-transverse weight-465 starting point."""
    coefficients = np.zeros(VARIABLES, dtype=np.uint8)
    for left in range(H):
        coefficients[left * H] = 1
    certificate = cpu_reconstruct(coefficients, geometry, THRESHOLD)
    if certificate["weight"] != THRESHOLD:
        raise ArithmeticError("the fixed-transverse baseline changed")
    return certificate


def emit(result: dict[str, Any], output: Path | None) -> None:
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if output is not None:
        output.write_text(text)
    print(text, end="", flush=True)


def backend_record(
    status: str,
    args: argparse.Namespace,
    baseline: dict[str, Any],
    detail: str,
) -> dict[str, Any]:
    return {
        "experiment": "p31_halved_row_code_low_weight_gpu",
        "status": status,
        "scope": "single-prime p=31 heuristic, not theorem evidence",
        "p": P,
        "ph_threshold": THRESHOLD,
        "seed": args.seed,
        "host": platform.node(),
        "detail": detail,
        "required_package_if_cuda_requested": (
            "CuPy wheel matching the CUDA runtime, for example cupy-cuda12x"
        ),
        "cpu_baseline": baseline,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=lambda value: int(value, 0), default=DEFAULT_SEED)
    parser.add_argument("--chains", type=int, default=2048)
    parser.add_argument("--restarts", type=int, default=16)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--steps-per-round", type=int, default=5000)
    parser.add_argument("--quench-steps", type=int, default=0)
    parser.add_argument(
        "--kick-moves",
        type=int,
        default=8,
        help="maximum forced kick; chains use one through this many moves",
    )
    parser.add_argument("--move2-percent", type=int, default=25)
    parser.add_argument("--temp-start", type=float, default=128.0)
    parser.add_argument("--temp-end", type=float, default=0.0)
    parser.add_argument("--threads", type=int, default=256)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.chains < 1 or args.restarts < 1 or args.rounds < 1:
        parser.error("chains, restarts, and rounds must be positive")
    if args.steps_per_round < 1 or args.quench_steps < 0 or args.kick_moves < 0:
        parser.error("step counts must be nonnegative and rounds nonempty")
    if not 0 <= args.move2_percent <= 100:
        parser.error("move2-percent must lie in 0..100")
    if args.temp_start <= 0 or args.temp_end < 0:
        parser.error("temperatures must be nonnegative with temp-start positive")
    if args.threads < 1 or args.threads > 1024:
        parser.error("threads must lie in 1..1024")
    return args


def main() -> int:
    args = parse_args()
    geometry = build_geometry()
    baseline = baseline_certificate(geometry)
    if args.cpu_only:
        emit(
            backend_record(
                "cpu_model_check_only",
                args,
                baseline,
                "exact normal-form geometry and weight-465 baseline verified",
            ),
            args.output,
        )
        return 0

    try:
        import cupy as cp
    except ImportError as error:
        emit(
            backend_record(
                "cuda_backend_unavailable",
                args,
                baseline,
                f"CuPy import failed: {error}",
            ),
            args.output,
        )
        return 2

    try:
        device = cp.cuda.Device(0)
        device.use()
        properties = cp.cuda.runtime.getDeviceProperties(0)
    except Exception as error:  # pragma: no cover - hardware-specific path
        emit(
            backend_record(
                "cuda_device_unavailable",
                args,
                baseline,
                f"CUDA device initialization failed: {type(error).__name__}: {error}",
            ),
            args.output,
        )
        return 2

    compile_start = time.perf_counter()
    init_kernel = cp.RawKernel(
        CUDA_SOURCE, "init_rectangle_chains", options=("-std=c++11",)
    )
    anneal_kernel = cp.RawKernel(
        CUDA_SOURCE, "anneal_block_code", options=("-std=c++11",)
    )
    init_kernel.compile()
    anneal_kernel.compile()
    compile_seconds = time.perf_counter() - compile_start

    chains = args.chains
    launch_blocks = math.ceil(chains / args.threads)
    state = cp.empty((chains, N, WORDS), dtype=cp.uint64)
    coefficients = cp.empty((chains, VARIABLES), dtype=cp.uint8)
    column_weights = cp.empty((chains, N), dtype=cp.uint16)
    objective = cp.empty(chains, dtype=cp.int32)
    coefficient_weight = cp.empty(chains, dtype=cp.int32)
    rng_state = cp.empty(chains, dtype=cp.uint32)
    best_objective = cp.empty(chains, dtype=cp.int32)
    best_coefficients = cp.empty((chains, VARIABLES), dtype=cp.uint8)
    block_masks = cp.asarray(geometry["block_masks"])
    star_masks = cp.asarray(geometry["star_masks"])
    block_points = cp.asarray(geometry["block_points"])

    kernel_args_common = (
        state,
        coefficients,
        column_weights,
        objective,
        coefficient_weight,
        rng_state,
        best_objective,
        best_coefficients,
    )
    init_geometry_args = (
        star_masks,
        block_points,
        np.int32(chains),
        np.int32(N),
        np.int32(WORDS),
        np.int32(P),
        np.int32(H),
        np.int32(D),
        np.int32(VARIABLES),
    )
    anneal_geometry_args = (
        block_masks,
        block_points,
        np.int32(chains),
        np.int32(N),
        np.int32(WORDS),
        np.int32(P),
        np.int32(H),
        np.int32(D),
        np.int32(VARIABLES),
    )

    search_start = time.perf_counter()
    rounds_completed = 0
    restarts_completed = 0
    quench_calls = 0
    found = False
    for restart in range(args.restarts):
        init_kernel(
            (launch_blocks,),
            (args.threads,),
            kernel_args_common
            + init_geometry_args
            + (
                np.uint32(args.seed),
                np.int32(restart),
                np.int32(restart == 0),
            ),
        )
        cp.cuda.runtime.deviceSynchronize()
        for round_index in range(args.rounds):
            fraction = (
                round_index / (args.rounds - 1) if args.rounds > 1 else 1.0
            )
            temperature = args.temp_start * (
                args.temp_end / args.temp_start
            ) ** fraction
            anneal_kernel(
                (launch_blocks,),
                (args.threads,),
                kernel_args_common
                + anneal_geometry_args
                + (
                    np.int32(args.steps_per_round),
                    np.int32(args.kick_moves if round_index == 0 else 0),
                    np.int32(args.move2_percent),
                    np.float32(temperature),
                ),
            )
            cp.cuda.runtime.deviceSynchronize()
            rounds_completed += 1
            current_best = int(cp.min(best_objective).get())
            print(
                f"restart={restart} round={round_index} "
                f"temperature={temperature:.6g} best={current_best}",
                file=sys.stderr,
                flush=True,
            )
            if 0 < current_best < THRESHOLD:
                found = True
                break
        if not found and args.quench_steps:
            anneal_kernel(
                (launch_blocks,),
                (args.threads,),
                kernel_args_common
                + anneal_geometry_args
                + (
                    np.int32(args.quench_steps),
                    np.int32(0),
                    np.int32(args.move2_percent),
                    np.float32(0.0),
                ),
            )
            cp.cuda.runtime.deviceSynchronize()
            quench_calls += 1
            current_best = int(cp.min(best_objective).get())
            print(
                f"restart={restart} quench best={current_best}",
                file=sys.stderr,
                flush=True,
            )
            found = 0 < current_best < THRESHOLD
        restarts_completed += 1
        if found:
            break

    search_seconds = time.perf_counter() - search_start
    best_values = cp.asnumpy(best_objective)
    current_values = cp.asnumpy(objective)
    current_coefficient_weights = cp.asnumpy(coefficient_weight)
    best_chain = int(np.argmin(best_values))
    gpu_best_weight = int(best_values[best_chain])
    best_coeff = cp.asnumpy(best_coefficients[best_chain])
    certificate = cpu_reconstruct(best_coeff, geometry, gpu_best_weight)
    audit_chain = int(np.argmax(current_values))
    audit_coeff = cp.asnumpy(coefficients[audit_chain])
    terminal_audit = cpu_reconstruct(
        audit_coeff, geometry, int(current_values[audit_chain])
    )
    if certificate["below_ph_465"] != found:
        raise ArithmeticError("GPU stopping status disagrees with CPU certificate")

    device_name_raw = properties["name"]
    device_name = (
        device_name_raw.decode()
        if isinstance(device_name_raw, bytes)
        else str(device_name_raw)
    )
    best_unique, best_counts = np.unique(best_values, return_counts=True)
    result = {
        "experiment": "p31_halved_row_code_low_weight_gpu",
        "status": (
            "counterexample_found_cpu_verified"
            if found
            else "no_subthreshold_word_found_heuristic_not_proof"
        ),
        "scope": "single-prime p=31 heuristic, not theorem evidence",
        "p": P,
        "h": H,
        "d": D,
        "delta_size": N,
        "ph_threshold": THRESHOLD,
        "seed": args.seed,
        "host": platform.node(),
        "device": device_name,
        "compute_capability": str(device.compute_capability),
        "cupy": cp.__version__,
        "configuration": {
            "chains": chains,
            "restarts_requested": args.restarts,
            "restarts_completed": restarts_completed,
            "rounds_requested_per_restart": args.rounds,
            "rounds_completed": rounds_completed,
            "quench_calls": quench_calls,
            "steps_per_round": args.steps_per_round,
            "quench_steps": args.quench_steps,
            "kick_moves": args.kick_moves,
            "kick_schedule": (
                "chain c uses 1+c mod kick_moves forced moves "
                "when kick_moves is positive"
            ),
            "move2_percent": args.move2_percent,
            "temp_start": args.temp_start,
            "temp_end": args.temp_end,
            "threads": args.threads,
            "launch_blocks": launch_blocks,
            "proposals_executed": int(
                chains
                * (
                    rounds_completed * args.steps_per_round
                    + quench_calls * args.quench_steps
                )
            ),
        },
        "normal_form": "W=1*q^T+B*X*B^T with X direction-block-diagonal",
        "q_optimization": "exact per physical column by majority complement",
        "moves": [
            "one direction-block coefficient",
            "one 2-by-2 rank-one cycle inside a direction block",
        ],
        "zero_word_rejected": True,
        "kernel_compile_seconds": compile_seconds,
        "search_seconds": search_seconds,
        "gpu_best_weight": gpu_best_weight,
        "best_chain": best_chain,
        "best_weight_histogram": {
            str(int(value)): int(count)
            for value, count in zip(best_unique, best_counts, strict=True)
        },
        "final_current_objective_summary": {
            "minimum": int(current_values.min()),
            "median": float(np.median(current_values)),
            "maximum": int(current_values.max()),
        },
        "final_current_coefficient_weight_summary": {
            "minimum": int(current_coefficient_weights.min()),
            "median": float(np.median(current_coefficient_weights)),
            "maximum": int(current_coefficient_weights.max()),
        },
        "cpu_terminal_state_audit": {
            "chain": audit_chain,
            **terminal_audit,
        },
        "cpu_certificate": certificate,
        "cpu_baseline": baseline,
        "conclusion": (
            "a nonzero Row(D) word below 465 was found"
            if found
            else "this deterministic heuristic run found no word below 465; "
            "this is not a minimum-distance certificate"
        ),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    emit(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
