#!/usr/bin/env python3
"""Mesh scan for simultaneous degree-six/degree-eight label sums.

The compact-ray ledger proves the degree-six and degree-eight conditions
separately, but deliberately leaves open whether the *same* field labels can
satisfy both.  For a prime ``p`` this program computes the exact-count sumsets

    R_p = {(a^6, a^8) : a in F_p^*},     N R_p subset F_p^2.

It records coverage for every requested ``N`` and reconstructs witnesses for
the two normalized targets used by both hard compact and all-equal rows,
``(0,0)`` and ``(-1/2,-1/2)``.  Repetitions are allowed, matching repeated
atom scales in the local directional model.

The fast backend uses a cyclic two-dimensional FFT.  Each step convolves two
0/1 arrays and therefore has integer coefficients in ``[0, |R_p|]``; the
program rounds only after checking the observed distance to the nearest
integer.  ``--backend direct`` is an independent integer implementation for
selected audit primes.  This is exhaustive finite-field discovery, not by
itself an all-prime theorem or a residual-(ii) close.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from io_atomic import write_json_atomic  # noqa: E402


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def power_pairs(p: int) -> list[tuple[int, int]]:
    if p < 5 or not is_prime(p):
        raise ValueError("p must be an odd prime at least five")
    return sorted({(pow(a, 6, p), pow(a, 8, p)) for a in range(1, p)})


def _array_digest(array: np.ndarray) -> str:
    packed = np.packbits(np.asarray(array, dtype=np.uint8), bitorder="little")
    return hashlib.sha256(packed.tobytes()).hexdigest()


def _fft_step(module: Any, current: Any, base_fft: Any) -> tuple[Any, float, int]:
    raw = module.fft.ifft2(module.fft.fft2(current) * base_fft).real
    rounded = module.rint(raw)
    error = float(module.max(module.abs(raw - rounded)).get()) if module is not np else float(np.max(np.abs(raw - rounded)))
    minimum = int(module.min(rounded).get()) if module is not np else int(np.min(rounded))
    if error >= 1.0e-5 or minimum < 0:
        raise ArithmeticError(
            f"unsafe FFT rounding: max error {error:.3e}, minimum {minimum}"
        )
    return rounded > 0, error, int(module.max(rounded).get()) if module is not np else int(np.max(rounded))


def _direct_step(current: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    out = np.zeros_like(current, dtype=bool)
    for first, second in pairs:
        out |= np.roll(np.roll(current, first, axis=0), second, axis=1)
    return out


def _witness(
    history: list[np.ndarray], pairs: list[tuple[int, int]], target: tuple[int, int], p: int
) -> list[int] | None:
    count = len(history) - 1
    if not bool(history[count][target]):
        return None
    residue_to_scale: dict[tuple[int, int], int] = {}
    for scale in range(1, p):
        residue_to_scale.setdefault((pow(scale, 6, p), pow(scale, 8, p)), scale)
    point = target
    scales: list[int] = []
    for step in range(count, 0, -1):
        for pair in pairs:
            previous = ((point[0] - pair[0]) % p, (point[1] - pair[1]) % p)
            if history[step - 1][previous]:
                scales.append(residue_to_scale[pair])
                point = previous
                break
        else:  # pragma: no cover - guards a corrupt reachability history
            raise ArithmeticError("failed to reconstruct a reported sumset point")
    scales.reverse()
    if (
        sum(pow(a, 6, p) for a in scales) % p,
        sum(pow(a, 8, p) for a in scales) % p,
    ) != target:
        raise ArithmeticError("reconstructed witness does not hit its target")
    return scales


def scan_prime(p: int, max_count: int, backend: str) -> dict[str, Any]:
    pairs = power_pairs(p)
    base = np.zeros((p, p), dtype=bool)
    for first, second in pairs:
        base[first, second] = True
    current_np = np.zeros_like(base)
    current_np[0, 0] = True
    history = [current_np.copy()]
    records: list[dict[str, Any]] = []
    start = time.time()

    if backend == "cupy":
        import cupy as cp

        module: Any = cp
        current: Any = cp.asarray(current_np.astype(np.float64))
        base_fft: Any = cp.fft.fft2(cp.asarray(base.astype(np.float64)))
        device_name = cp.cuda.runtime.getDeviceProperties(0)["name"]
        if isinstance(device_name, bytes):
            device_name = device_name.decode(errors="replace")
    elif backend == "numpy":
        module = np
        current = current_np.astype(np.float64)
        base_fft = np.fft.fft2(base.astype(np.float64))
        device_name = "CPU"
    elif backend == "opencl":
        import pyopencl as cl

        platforms = cl.get_platforms()
        devices = [
            device
            for platform_ in platforms
            for device in platform_.get_devices()
            if device.type & cl.device_type.GPU
        ]
        if not devices:
            raise RuntimeError("no OpenCL GPU device found")
        device = devices[0]
        context = cl.Context([device])
        queue = cl.CommandQueue(context)
        program = cl.Program(
            context,
            r"""
            __kernel void sumset_step(
                __global const uchar *current,
                __global const int2 *pairs,
                __global uchar *output,
                const int pair_count,
                const int p)
            {
                int index = get_global_id(0);
                int total = p*p;
                if (index >= total) return;
                int first = index / p;
                int second = index - first*p;
                uchar hit = 0;
                for (int k=0; k<pair_count; ++k) {
                    int a = first-pairs[k].x;
                    int b = second-pairs[k].y;
                    if (a < 0) a += p;
                    if (b < 0) b += p;
                    if (current[a*p+b]) { hit=1; break; }
                }
                output[index] = hit;
            }
            """,
        ).build()
        pair_array = np.asarray(pairs, dtype=np.int32)
        flags = cl.mem_flags
        pair_buffer = cl.Buffer(
            context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=pair_array
        )
        current = np.asarray(current_np, dtype=np.uint8)
        base_fft = (cl, context, queue, program, pair_buffer)
        module = np
        device_name = device.name.strip()
    elif backend == "direct":
        module = np
        current = current_np
        base_fft = None
        device_name = "CPU exact integer shifts"
    else:
        raise ValueError("backend must be cupy, numpy, opencl, or direct")

    minus_half = (-pow(2, -1, p)) % p
    targets = ((0, 0), (minus_half, minus_half))
    first_full: int | None = None
    for count in range(1, max_count + 1):
        if backend == "direct":
            current_np = _direct_step(np.asarray(current, dtype=bool), pairs)
            current = current_np
            rounding_error = 0.0
            convolution_max = len(pairs)
        elif backend == "opencl":
            cl, context, queue, program, pair_buffer = base_fft
            flags = cl.mem_flags
            current_buffer = cl.Buffer(
                context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=current
            )
            next_array = np.empty_like(current, dtype=np.uint8)
            next_buffer = cl.Buffer(context, flags.WRITE_ONLY, next_array.nbytes)
            program.sumset_step(
                queue,
                (p * p,),
                None,
                current_buffer,
                pair_buffer,
                next_buffer,
                np.int32(len(pairs)),
                np.int32(p),
            )
            cl.enqueue_copy(queue, next_array, next_buffer).wait()
            current = next_array
            current_np = next_array.astype(bool)
            rounding_error = 0.0
            convolution_max = len(pairs)
        else:
            current, rounding_error, convolution_max = _fft_step(module, current, base_fft)
            current_np = module.asnumpy(current) if module is not np else np.asarray(current)
        current_np = np.asarray(current_np, dtype=bool)
        history.append(current_np.copy())
        reached = int(np.count_nonzero(current_np))
        row = {
            "count": count,
            "reached": reached,
            "total": p * p,
            "coverage": reached / (p * p),
            "zero_target": bool(current_np[0, 0]),
            "minus_half_diagonal_target": bool(current_np[minus_half, minus_half]),
            "max_fft_rounding_error": rounding_error,
            "convolution_max": convolution_max,
            "packed_reachability_sha256": _array_digest(current_np),
        }
        records.append(row)
        if reached == p * p and first_full is None:
            first_full = count
            break

    witness_count = first_full if first_full is not None else len(history) - 1
    witnesses: dict[str, Any] = {}
    for name, target in zip(("zero", "minus_half_diagonal"), targets, strict=True):
        witness = _witness(history[: witness_count + 1], pairs, target, p)
        witnesses[name] = {
            "target": list(target),
            "count": witness_count,
            "scales": witness,
        }

    return {
        "p": p,
        "p_mod_4": p % 4,
        "distinct_power_pairs": len(pairs),
        "base_pair_sha256": hashlib.sha256(
            np.asarray(pairs, dtype="<i8").tobytes()
        ).hexdigest(),
        "backend": backend,
        "device": str(device_name),
        "first_full_exact_count": first_full,
        "records": records,
        "witnesses_at_terminal_count": witnesses,
        "elapsed_seconds": time.time() - start,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", type=int, nargs="+", required=True)
    parser.add_argument("--max-count", type=int, default=16)
    parser.add_argument(
        "--backend", choices=("cupy", "numpy", "opencl", "direct"), default="numpy"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_count < 1:
        raise ValueError("max-count must be positive")
    invalid = [p for p in args.primes if p < 5 or not is_prime(p)]
    if invalid:
        raise ValueError(f"invalid primes: {invalid}")
    started = time.time()
    rows = [scan_prime(p, args.max_count, args.backend) for p in args.primes]
    payload = {
        "classification": (
            "exhaustive finite-field paired-power sumset discovery; not an all-prime proof"
        ),
        "host": platform.node(),
        "backend": args.backend,
        "primes": args.primes,
        "max_count": args.max_count,
        "curve": "{(a^6,a^8):a in F_p^*}",
        "repetitions_allowed": True,
        "normalized_targets": ["(0,0)", "(-1/2,-1/2)"],
        "results": rows,
        "elapsed_seconds": time.time() - started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    write_json_atomic(args.output, payload)
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
