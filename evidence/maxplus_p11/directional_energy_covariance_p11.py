#!/usr/bin/env python3
"""Directional profile-energy covariance on the full p=11 Max+ census.

This tests one possible route to the remaining exceptional QVAR inequality:
whether the quartic sign vector is forced to be a top covariance mode.  It is
diagnostic evidence only.  No finite-p covariance pattern is promoted to a
general theorem.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import cupy as cp
import numpy as np


P = 11
Q = P * P
N = Q + 1
COUNTS = {
    "k1": 2_772,
    "k3": 24_200,
    "k4": 58_080,
    "k5": 1_306_800,
    "k6": 36_065_260,
}


def _context():
    evidence = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(evidence))
    from quartic_profile_attack import (  # type: ignore
        field_context,
        projective_directions,
        quartic_character,
    )

    q, mul, trace, _add = field_context(P)
    psi = quartic_character(P, mul)
    directions = projective_directions(P, mul, trace)
    active = [(g, t_of) for g, t_of in directions if abs(psi[g].imag) < 0.5]
    weights = np.asarray([int(round(psi[g].real)) for g, _ in active], dtype=np.int64)
    if q != Q or len(active) != (P + 1) // 2 or weights.sum() != 0:
        raise RuntimeError("field/direction audit failed")
    return active, weights


def _empty_acc(m: int) -> dict:
    return {
        "count": 0,
        "sum": np.zeros(m, dtype=np.int64),
        "gram": np.zeros((m, m), dtype=np.float64),
        "imbalance_histogram": Counter(),
    }


def _finish(acc: dict, weights: np.ndarray) -> dict:
    count = acc["count"]
    mean = acc["sum"].astype(np.float64) / count
    second = acc["gram"] / count
    covariance = second - np.outer(mean, mean)
    covariance = (covariance + covariance.T) / 2
    eig = np.linalg.eigvalsh(covariance)
    quartic = float(weights @ covariance @ weights / (weights @ weights))
    return {
        "count": count,
        "mean_energy": mean.tolist(),
        "covariance": covariance.tolist(),
        "covariance_eigenvalues": eig.tolist(),
        "quartic_covariance_eigenvalue_if_mode": quartic,
        "quartic_is_top_mode": bool(abs(quartic - eig[-1]) < 1e-8),
        "quartic_to_top_ratio": quartic / eig[-1],
        "E_signed_imbalance_sq": quartic * float(weights @ weights),
        "imbalance_histogram": {
            str(k): int(v) for k, v in sorted(acc["imbalance_histogram"].items())
        },
    }


def main() -> None:
    root = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    rows = np.load(root / "maxplus_p11_eps1.npy", mmap_mode="r")
    if rows.shape != (sum(COUNTS.values()), N):
        raise RuntimeError(f"unexpected Max+ shape {rows.shape}")

    active, weights = _context()
    m = len(active)
    incidence = np.zeros((Q, m * P), dtype=np.float32)
    for j, (_g, t_of) in enumerate(active):
        incidence[np.arange(Q), j * P + t_of] = 1
    incidence_gpu = cp.asarray(incidence)

    chunk = int(os.environ.get("PROFILE_COV_CHUNK", "100000"))
    acc = {label: _empty_acc(m) for label in COUNTS}
    acc["total"] = _empty_acc(m)
    bounds = []
    start = 0
    for label, count in COUNTS.items():
        bounds.append((start, start + count, label))
        start += count

    started = time.time()
    for lo in range(0, len(rows), chunk):
        hi = min(len(rows), lo + chunk)
        y = cp.asarray(rows[lo:hi, 1:], dtype=cp.float32)
        sigma = (y @ incidence_gpu).reshape(hi - lo, m, P)
        energy_gpu = cp.rint(cp.sum((sigma - 1) ** 2, axis=2) / 4).astype(cp.int32)
        energy = cp.asnumpy(energy_gpu).astype(np.int64)
        del y, sigma, energy_gpu

        if not np.all(energy.sum(axis=1) == P * (Q - 1) // 4):
            raise RuntimeError("profile-energy conservation failed")
        activity = np.count_nonzero(energy, axis=1)
        imbalance = energy @ weights

        for blo, bhi, label in bounds:
            left, right = max(lo, blo), min(hi, bhi)
            if left >= right:
                continue
            sl = slice(left - lo, right - lo)
            if not np.all(activity[sl] == int(label[1:])):
                raise RuntimeError(f"stratum order failed in {label}")
            for key in (label, "total"):
                part = energy[sl]
                rec = acc[key]
                rec["count"] += len(part)
                rec["sum"] += part.sum(axis=0)
                part_gpu = cp.asarray(part, dtype=cp.float64)
                rec["gram"] += cp.asnumpy(part_gpu.T @ part_gpu)
                del part_gpu
                vals, freq = np.unique(imbalance[sl], return_counts=True)
                rec["imbalance_histogram"].update(
                    {int(v): int(c) for v, c in zip(vals, freq)}
                )
        if (lo // chunk + 1) % 25 == 0 or hi == len(rows):
            cp.get_default_memory_pool().free_all_blocks()
            print(f"rows {hi}/{len(rows)} elapsed={time.time()-started:.1f}s", flush=True)

    result = {
        "p": P,
        "weights": weights.tolist(),
        "strata": {label: _finish(acc[label], weights) for label in COUNTS},
        "total": _finish(acc["total"], weights),
        "interpretation_rule": (
            "A top-mode pattern at p=11 is diagnostic only; a general cyclic or "
            "group-covariance argument is required before it can support QVAR."
        ),
        "elapsed_seconds": time.time() - started,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    output = root / "directional_energy_covariance_p11.json"
    tracked_output = Path(__file__).with_suffix(".json")
    output.write_text(rendered)
    tracked_output.write_text(rendered)
    print(json.dumps(result, indent=2), flush=True)
    print(f"wrote {output}", flush=True)
    print(f"wrote {tracked_output}", flush=True)


if __name__ == "__main__":
    main()
