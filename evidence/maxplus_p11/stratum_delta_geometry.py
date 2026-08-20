#!/usr/bin/env python3
"""Decompose the p=11 four-point residual by profile stratum.

This is an exploratory calculation for the shared leftover-1/leftover-3
four-point tensor.  For every probability measure supported on Max+, the
master equation has the same particular solution ``m4_flat`` and

    m4 = m4_flat + delta,       delta in E_{4p},
    ||m4||^2 = ||m4_flat||^2 + ||delta||^2.

The p=11 enumeration is stored in profile order k=1,3,4,5,6.  We compute the
pair-moment Gram for k<=5 on the V100, obtain k=6 by exact subtraction from
the already verified total Gram, and report all pairwise inner products of
the residual vectors.  Large negative off-diagonal terms would identify the
cancellation responsible for the very small full-ensemble residual.

The script deliberately does not assert a general theorem or flip a flag.
It writes one JSON record next to the persistent p=11 data.
"""
from __future__ import annotations

import itertools
import json
import os
import time
from fractions import Fraction
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
ORDER = tuple(COUNTS)


def m4_flat_norm_sq(p: int) -> Fraction:
    """||(master particular four-point tensor)||^2."""
    return Fraction(
        (p - 1) * (p + 1) * (p * p + 1) * (3 * p * p + 17),
        24 * (p * p - 5),
    )


def activity_count_sample(rows: np.ndarray, finite_line_maps: np.ndarray) -> list[int]:
    """Return active-direction counts for a small eps=+1 sample."""
    out: list[int] = []
    for row in rows:
        yf = row[1:]
        active = 0
        for tmap in finite_line_maps:
            sums = np.bincount(tmap, weights=yf, minlength=P)
            active += int(np.any(sums != 1))
        out.append(active)
    return out


def accumulate_pair_gram(
    rows: np.ndarray,
    iu0_gpu: cp.ndarray,
    iu1_gpu: cp.ndarray,
    chunk: int,
    label: str,
) -> np.ndarray:
    """Return exact integer Q.T@Q for pair features Q_ab=y_a*y_b."""
    npair = int(iu0_gpu.size)
    acc = np.zeros((npair, npair), dtype=np.int64)
    started = time.time()
    nch = (len(rows) + chunk - 1) // chunk
    for ci, lo in enumerate(range(0, len(rows), chunk), start=1):
        yg = cp.asarray(rows[lo : lo + chunk], dtype=cp.int8)
        qg = (yg[:, iu0_gpu] * yg[:, iu1_gpu]).astype(cp.float32)
        # Every entry is an integer of magnitude <= chunk < 2^24, so V100
        # fp32 GEMM is exact before conversion to int64.
        gram = qg.T @ qg
        acc += cp.asnumpy(gram).astype(np.int64)
        del yg, qg, gram
        if ci % 10 == 0 or ci == nch:
            cp.get_default_memory_pool().free_all_blocks()
            print(
                f"{label}: chunk {ci}/{nch}, elapsed={time.time()-started:.1f}s",
                flush=True,
            )
    return acc


def main() -> None:
    data_root = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    output = data_root / "stratum_delta_geometry_p11.json"
    yall = np.load(data_root / "maxplus_p11_eps1.npy", mmap_mode="r")
    gtotal = np.load(data_root / "G_pairmoment_p11.npy", mmap_mode="r")

    if yall.shape != (sum(COUNTS.values()), N):
        raise RuntimeError(f"unexpected Max+ shape {yall.shape}")
    if not np.all(yall[:1000, 0] == 1):
        raise RuntimeError("eps=+1 gauge check failed")

    # The finite line maps are only needed to validate that the persisted
    # assembly still has the documented k-stratum order.
    import sys

    sys.path.insert(0, "/mnt/storage/e1work/scripts")
    from kgen import square_coords  # type: ignore

    _dirs, _forms, coords = square_coords(P)
    finite_line_maps = np.asarray(coords, dtype=np.int64)

    starts: dict[str, tuple[int, int]] = {}
    lo = 0
    for label in ORDER:
        hi = lo + COUNTS[label]
        starts[label] = (lo, hi)
        sample_idx = np.linspace(lo, hi - 1, min(12, hi - lo), dtype=np.int64)
        observed = activity_count_sample(np.asarray(yall[sample_idx]), finite_line_maps)
        expected = int(label[1:])
        if observed != [expected] * len(observed):
            raise RuntimeError(
                f"stratum-order check failed for {label}: observed {observed}"
            )
        lo = hi

    iu = np.triu_indices(N, 1)
    iu0_gpu = cp.asarray(iu[0].astype(np.int32))
    iu1_gpu = cp.asarray(iu[1].astype(np.int32))
    chunk = int(os.environ.get("PAIR_GRAM_CHUNK", "50000"))

    grams: dict[str, np.ndarray] = {}
    for label in ("k1", "k3", "k4", "k5"):
        slo, shi = starts[label]
        grams[label] = accumulate_pair_gram(
            yall[slo:shi], iu0_gpu, iu1_gpu, chunk, label
        )

    known_sum = sum(grams.values())
    grams["k6"] = np.asarray(gtotal, dtype=np.int64) - known_sum
    if not np.array_equal(sum(grams.values()), np.asarray(gtotal, dtype=np.int64)):
        raise RuntimeError("stratum Grams do not reconstruct the total Gram")

    pair_id = -np.ones((N, N), dtype=np.int32)
    pair_id[iu] = np.arange(len(iu[0]), dtype=np.int32)
    pair_id[(iu[1], iu[0])] = pair_id[iu]

    quads = np.asarray(list(itertools.combinations(range(N), 4)), dtype=np.int16)
    a, b, c, d = (quads[:, j].astype(np.int32) for j in range(4))
    ab_cd = (pair_id[a, b], pair_id[c, d])
    ac_bd = (pair_id[a, c], pair_id[b, d])
    ad_bc = (pair_id[a, d], pair_id[b, c])

    moments: dict[str, np.ndarray] = {}
    consistency: dict[str, int] = {}
    for label in ORDER:
        gram = grams[label]
        x = gram[ab_cd]
        y = gram[ac_bd]
        z = gram[ad_bc]
        consistency[label] = int(max(np.max(np.abs(x - y)), np.max(np.abs(x - z))))
        if consistency[label] != 0:
            raise RuntimeError(f"four-point pairing mismatch in {label}")
        moments[label] = x.astype(np.float64) / COUNTS[label]

    flat = float(m4_flat_norm_sq(P))
    residual_gram = np.empty((len(ORDER), len(ORDER)), dtype=np.float64)
    m4_gram = np.empty_like(residual_gram)
    for i, left in enumerate(ORDER):
        for j, right in enumerate(ORDER):
            ip = float(np.dot(moments[left], moments[right]))
            m4_gram[i, j] = ip
            residual_gram[i, j] = ip - flat

    weights = np.asarray([COUNTS[k] / len(yall) for k in ORDER])
    residual_norms = np.sqrt(np.maximum(np.diag(residual_gram), 0.0))
    residual_cosines = residual_gram / np.outer(residual_norms, residual_norms)
    global_delta_sq = float(weights @ residual_gram @ weights)
    direct_total = np.asarray(gtotal[ab_cd], dtype=np.float64) / len(yall)
    direct_total_delta_sq = float(np.dot(direct_total, direct_total) - flat)

    report = {
        "p": P,
        "counts": COUNTS,
        "weights": dict(zip(ORDER, weights.tolist())),
        "stratum_order_validated": True,
        "pairing_consistency_max_integer_error": consistency,
        "m4_flat_norm_sq": str(m4_flat_norm_sq(P)),
        "m4_norm_sq_by_stratum": {
            k: float(m4_gram[i, i]) for i, k in enumerate(ORDER)
        },
        "delta_norm_sq_by_stratum": {
            k: float(residual_gram[i, i]) for i, k in enumerate(ORDER)
        },
        "delta_inner_product_matrix_order": ORDER,
        "delta_inner_product_matrix": residual_gram.tolist(),
        "delta_cosine_matrix": residual_cosines.tolist(),
        "weighted_global_delta_norm_sq": global_delta_sq,
        "direct_global_delta_norm_sq": direct_total_delta_sq,
        "weighted_vs_direct_abs_error": abs(global_delta_sq - direct_total_delta_sq),
        "interpretation_rule": (
            "Only a general cancellation or norm mechanism counts as progress; "
            "these p=11 stratum numbers are diagnostic data, not a theorem."
        ),
    }
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
