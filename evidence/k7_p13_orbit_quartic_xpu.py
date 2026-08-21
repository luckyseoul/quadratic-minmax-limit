#!/usr/bin/env python3
"""Exact quartic histogram on the signed-PSL orbit of the p=13 witness.

The packed orbit is produced by ``k7_p13_signed_psl_orbit.py`` on the fast
single-thread node.  This second stage uses Torch XPU on Jellyfin's Intel Arc
A380.  Float32 is exact here: every matrix entry and input bit is integral,
and every intermediate integer has magnitude below 2**24.  A direct int64
CPU sample audit is nevertheless included in the certificate.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "evidence"))

from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k5_p29_coefficient_sieve import square_directions  # noqa: E402


def unpack_finite_bits(packed: np.ndarray, q: int) -> np.ndarray:
    """Unpack finite-coordinate negative indicators, skipping bit 0=inf."""
    indices = np.arange(1, q + 1, dtype=np.int64)
    words = indices // 64
    shifts = (indices % 64).astype(np.uint64)
    return ((packed[:, words] >> shifts) & np.uint64(1)).astype(np.float32)


def direct_values(
    bits: np.ndarray, kernel_real: np.ndarray, kernel_imag: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    bits64 = bits.astype(np.int64)
    real = np.einsum(
        "bi,ij,bj->b", bits64, kernel_real.astype(np.int64), bits64,
        optimize=True,
    )
    imag = np.einsum(
        "bi,ij,bj->b", bits64, kernel_imag.astype(np.int64), bits64,
        optimize=True,
    )
    return real, imag


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("orbit", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--batch-size", type=int, default=65_536)
    parser.add_argument(
        "--epsilon-plus-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="global complements have the same Zpsi, so one sign half suffices",
    )
    args = parser.parse_args()

    if not torch.xpu.is_available():
        raise RuntimeError("Torch XPU is not available")
    import dpctl

    p = 13
    q = p * p
    packed = np.load(args.orbit, mmap_mode="r")
    if packed.ndim != 2 or packed.shape[1] != 3:
        raise RuntimeError(f"unexpected packed orbit shape {packed.shape}")
    kernel_real, kernel_imag = quartic_kernel(p)
    line_incidence = np.zeros((q, 7 * p), dtype=np.float32)
    for direction, (coordinate, _form) in enumerate(square_directions(p)):
        for value in range(p):
            line_incidence[coordinate == value, direction * p + value] = 1
    combined = np.concatenate(
        (kernel_real, kernel_imag, line_incidence), axis=1
    ).astype(np.float32)
    device = torch.device("xpu")
    kernel_xpu = torch.from_numpy(combined).to(device)

    histogram: Counter[int] = Counter()
    k_histogram: Counter[int] = Counter()
    k_abs_sq_sum: Counter[int] = Counter()
    k_abs_sq_histogram: Counter[tuple[int, int]] = Counter()
    total_abs_sq = 0
    total_rows = 0
    cpu_audit_rows = 0
    cpu_audit_ok = True
    started = time.monotonic()
    torch.set_num_threads(1)
    with torch.no_grad():
        for lo in range(0, len(packed), args.batch_size):
            hi = min(lo + args.batch_size, len(packed))
            packed_batch = np.asarray(packed[lo:hi])
            if args.epsilon_plus_only:
                packed_batch = packed_batch[(packed_batch[:, 0] & 1) == 0]
            if len(packed_batch) == 0:
                continue
            finite = unpack_finite_bits(packed_batch, q)
            bits_xpu = torch.from_numpy(finite).to(device)
            products = bits_xpu @ kernel_xpu
            real_xpu = torch.sum(products[:, :q] * bits_xpu, dim=1)
            imag_xpu = torch.sum(products[:, q : 2 * q] * bits_xpu, dim=1)
            real = np.rint(real_xpu.cpu().numpy()).astype(np.int64)
            imag = np.rint(imag_xpu.cpu().numpy()).astype(np.int64)
            line_negative_counts = np.rint(
                products[:, 2 * q :].cpu().numpy()
            ).astype(np.int16).reshape(-1, 7, p)
            infinity_negative = (packed_batch[:, 0] & 1).astype(np.int16)
            inactive_target = 6 + infinity_negative
            activity = np.any(
                line_negative_counts != inactive_target[:, None, None], axis=2
            )
            k_values = np.sum(activity, axis=1, dtype=np.int16)

            if cpu_audit_rows < 256:
                take = min(256 - cpu_audit_rows, len(finite))
                exact_real, exact_imag = direct_values(
                    finite[:take], kernel_real, kernel_imag
                )
                cpu_audit_ok = cpu_audit_ok and np.array_equal(
                    real[:take], exact_real
                ) and np.array_equal(imag[:take], exact_imag)
                cpu_audit_rows += take

            abs_sq = real * real + imag * imag
            values, counts = np.unique(abs_sq, return_counts=True)
            histogram.update(
                {int(value): int(count) for value, count in zip(values, counts)}
            )
            joint = np.stack((k_values.astype(np.int64), abs_sq), axis=1)
            joint_values, joint_counts = np.unique(
                joint, axis=0, return_counts=True
            )
            for (k_value, abs_value), count in zip(joint_values, joint_counts):
                k = int(k_value)
                value = int(abs_value)
                multiplicity = int(count)
                k_histogram[k] += multiplicity
                k_abs_sq_sum[k] += value * multiplicity
                k_abs_sq_histogram[(k, value)] += multiplicity
            total_abs_sq += int(np.sum(abs_sq, dtype=np.int64))
            total_rows += len(abs_sq)
            if lo == 0 or hi == len(packed) or hi % (10 * args.batch_size) == 0:
                print(
                    f"packed={hi:,}/{len(packed):,} evaluated={total_rows:,} "
                    f"elapsed={time.monotonic() - started:.1f}s",
                    flush=True,
                )
    torch.xpu.synchronize()
    elapsed = time.monotonic() - started
    if not cpu_audit_ok:
        raise RuntimeError("XPU quartic values failed the exact CPU audit")

    mean = Fraction(total_abs_sq, total_rows)
    threshold = Fraction(3 * q * (q - 1), 16)
    scalar = Fraction(32, q * (q - 1)) * mean
    by_activity = {}
    for k in sorted(k_histogram):
        k_mean = Fraction(k_abs_sq_sum[k], k_histogram[k])
        by_activity[str(k)] = {
            "count": k_histogram[k],
            "sum_abs_Zpsi_sq": k_abs_sq_sum[k],
            "E_abs_Zpsi_sq": str(k_mean),
            "clears_QVAR": k_mean >= threshold,
            "histogram": {
                str(value): k_abs_sq_histogram[(k, value)]
                for kk, value in sorted(k_abs_sq_histogram)
                if kk == k
            },
        }
    report = {
        "p": p,
        "packed_orbit": str(args.orbit),
        "packed_rows": len(packed),
        "epsilon_plus_only": args.epsilon_plus_only,
        "evaluated_rows": total_rows,
        "device": torch.xpu.get_device_name(0),
        "torch_version": torch.__version__,
        "dpctl_version": dpctl.__version__,
        "batch_size": args.batch_size,
        "float32_exact_integer_bound": True,
        "cpu_int64_audit_rows": cpu_audit_rows,
        "cpu_int64_audit_ok": cpu_audit_ok,
        "sum_abs_Zpsi_sq": total_abs_sq,
        "E_abs_Zpsi_sq": str(mean),
        "QVAR_threshold": str(threshold),
        "orbit_average_clears_QVAR": mean >= threshold,
        "induced_exceptional_scalar": str(scalar),
        "min_abs_Zpsi_sq": min(histogram),
        "max_abs_Zpsi_sq": max(histogram),
        "histogram": {
            str(value): histogram[value] for value in sorted(histogram)
        },
        "by_activity": by_activity,
        "elapsed_seconds": elapsed,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    main()
