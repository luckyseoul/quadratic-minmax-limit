#!/usr/bin/env python3
"""Full p=11 exceptional-projection shell census on an Intel XPU.

The complete epsilon=+1 Max+ array is already ordered by profile activity
``k=1,3,4,5,6``.  This program computes the squared norm of every vector's
projection onto the unique 61-dimensional exceptional constituent, quantizes
the result in the empirically exact unit ``16/15``, and records full and
per-stratum shell histograms.

Environment:
  E1WORK_P11              directory containing the two untracked p=11 arrays
  EXCEPTIONAL_XPU_CHUNK   row batch size (default 8192)
  EXCEPTIONAL_SHELL_OUTPUT output JSON path
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15588 import z_basis  # noqa: E402


P = 11
N = 122
EXCEPTIONAL_DIM = 61
QUANTUM = Fraction(16, 15)
COUNTS = {
    "k1": 2_772,
    "k3": 24_200,
    "k4": 58_080,
    "k5": 1_306_800,
    "k6": 36_065_260,
}


def _clusters(values: np.ndarray, tol: float = 1e-6) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    lo = 0
    for hi in range(1, len(values) + 1):
        if hi == len(values) or values[hi] - values[hi - 1] > tol:
            out.append((lo, hi))
            lo = hi
    return out


def _exceptional_forms(phi: np.ndarray) -> tuple[np.ndarray, float]:
    values, vectors = np.linalg.eigh(phi)
    lo, hi = next(
        block
        for block in _clusters(values)
        if block[1] - block[0] == EXCEPTIONAL_DIM
    )
    forms = np.einsum(
        "tr,tij->rij", vectors[:, lo:hi], z_basis(P), optimize=True
    ).astype(np.float32)
    return forms, float(values[lo:hi].mean())


def _render_histogram(histogram: Counter[int]) -> dict[str, int]:
    return {
        str(QUANTUM * quantum_index): int(count)
        for quantum_index, count in sorted(histogram.items())
    }


def main() -> dict:
    if not hasattr(torch, "xpu") or not torch.xpu.is_available():
        raise RuntimeError("PyTorch XPU is unavailable")

    work = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    rows = np.load(work / "maxplus_p11_eps1.npy", mmap_mode="r")
    phi = np.load(work / "phiZ_p11.npy")
    if rows.shape != (sum(COUNTS.values()), N):
        raise RuntimeError(f"unexpected Max+ shape {rows.shape}")

    started = time.time()
    forms_cpu, eigenvalue = _exceptional_forms(phi)
    forms = torch.as_tensor(forms_cpu, device="xpu")
    chunk = int(os.environ.get("EXCEPTIONAL_XPU_CHUNK", "8192"))
    total_hist: Counter[int] = Counter()
    stratum_hist = {label: Counter() for label in COUNTS}
    bounds = []
    start = 0
    for label, count in COUNTS.items():
        bounds.append((start, start + count, label))
        start += count

    worst_quantization_error = 0.0
    for batch_index, lo in enumerate(range(0, len(rows), chunk), start=1):
        hi = min(len(rows), lo + chunk)
        y = torch.as_tensor(np.asarray(rows[lo:hi], dtype=np.float32), device="xpu")
        # Broadcasted batched GEMM: (1,b,n) @ (r,n,n) -> (r,b,n).
        transformed = torch.matmul(y.unsqueeze(0), forms)
        coordinates = torch.sum(transformed * y.unsqueeze(0), dim=2).T
        norms = torch.sum(coordinates * coordinates, dim=1)
        scaled = norms * (QUANTUM.denominator / QUANTUM.numerator)
        quantum_indices_gpu = torch.round(scaled)
        error = float(torch.max(torch.abs(scaled - quantum_indices_gpu)).cpu())
        worst_quantization_error = max(worst_quantization_error, error)
        quantum_indices = quantum_indices_gpu.to(torch.int64).cpu().numpy()
        del y, transformed, coordinates, norms, scaled, quantum_indices_gpu

        values, frequencies = np.unique(quantum_indices, return_counts=True)
        total_hist.update({int(v): int(c) for v, c in zip(values, frequencies)})
        for bound_lo, bound_hi, label in bounds:
            left, right = max(lo, bound_lo), min(hi, bound_hi)
            if left >= right:
                continue
            part = quantum_indices[left - lo : right - lo]
            values, frequencies = np.unique(part, return_counts=True)
            stratum_hist[label].update(
                {int(v): int(c) for v, c in zip(values, frequencies)}
            )

        if batch_index % 100 == 0 or hi == len(rows):
            torch.xpu.synchronize()
            print(
                f"rows {hi}/{len(rows)} shells={len(total_hist)} "
                f"elapsed={time.time()-started:.1f}s",
                flush=True,
            )

    if sum(total_hist.values()) != len(rows):
        raise RuntimeError("histogram count mismatch")
    if worst_quantization_error >= 0.01:
        raise RuntimeError(f"quantization error too large: {worst_quantization_error}")

    weighted_index_sum = sum(index * count for index, count in total_hist.items())
    mean_norm = QUANTUM * Fraction(weighted_index_sum, len(rows))
    expected_mean = EXCEPTIONAL_DIM * eigenvalue
    result = {
        "p": P,
        "n": N,
        "device": torch.xpu.get_device_name(0),
        "torch": torch.__version__,
        "row_count_eps_plus": len(rows),
        "quantum": str(QUANTUM),
        "worst_quantization_error": worst_quantization_error,
        "lambda_exceptional": eigenvalue,
        "mean_projection_norm_sq": str(mean_norm),
        "mean_projection_norm_sq_float": float(mean_norm),
        "expected_dim_times_lambda": expected_mean,
        "mean_match_error": abs(float(mean_norm) - expected_mean),
        "target_3n": 3 * N,
        "zero_count": int(total_hist.get(0, 0)),
        "min_nonzero_norm_sq": str(
            QUANTUM * min(index for index in total_hist if index > 0)
        ),
        "histogram": _render_histogram(total_hist),
        "strata": {
            label: {
                "count": sum(hist.values()),
                "mean_projection_norm_sq": str(
                    QUANTUM
                    * Fraction(
                        sum(index * count for index, count in hist.items()),
                        sum(hist.values()),
                    )
                ),
                "zero_count": int(hist.get(0, 0)),
                "min_nonzero_norm_sq": str(
                    QUANTUM * min(index for index in hist if index > 0)
                ),
                "histogram": _render_histogram(hist),
            }
            for label, hist in stratum_hist.items()
        },
        "elapsed_seconds": time.time() - started,
    }
    output = Path(
        os.environ.get(
            "EXCEPTIONAL_SHELL_OUTPUT",
            str(Path(__file__).with_suffix(".json")),
        )
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    print(f"wrote {output}", flush=True)
    return result


if __name__ == "__main__":
    main()
