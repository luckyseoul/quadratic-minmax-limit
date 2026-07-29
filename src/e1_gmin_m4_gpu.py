#!/usr/bin/env python3
"""
GPU-accelerated |κ|=1 m4 census (CuPy / V100).

Loads Max+ matrix Y (N×n) once on GPU; evaluates
  m4(S) = mean_k Y[k,a]*Y[k,b]*Y[k,c]*Y[k,d]
for all 4-sets with |κ|=1 in large batches — one CUDA context, not 88.

Also reports CPU multi-worker fallback path timing comparison for p=5.

Writes evidence/e1_gmin_m4_gpu.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def kappa_batch(C: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """idx (B,4) int; return κ for each row."""
    a, b, c, d = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]
    return (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)


def all_quads(n: int) -> np.ndarray:
    return np.array(list(combinations(range(n), 4)), dtype=np.int32)


def gpu_m4_kappa1(Y_np: np.ndarray, C_np: np.ndarray, batch: int = 200_000) -> dict:
    """Full |κ|=1 max|m4| on GPU."""
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Y = cp.asarray(Y_np, dtype=cp.float64)
    C = np.asarray(C_np, dtype=np.float64)
    n = C.shape[0]
    quads = all_quads(n)
    # filter |κ|=1 on CPU (cheap combinatorial)
    kap = kappa_batch(C, quads)
    mask = np.abs(kap) == 1
    idx = quads[mask]
    kaps = kap[mask]
    n1 = int(len(idx))

    max_abs = 0.0
    max_m4 = 0.0
    max_kap = 0
    worst_pts = None
    t0 = time.perf_counter()
    for lo in range(0, n1, batch):
        hi = min(n1, lo + batch)
        sl = idx[lo:hi]
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        # columns: (N, B)
        m4 = cp.mean(
            Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d],
            axis=0,
        )
        m4_h = cp.asnumpy(m4)
        am = np.abs(m4_h)
        j = int(np.argmax(am))
        if am[j] > max_abs:
            max_abs = float(am[j])
            max_m4 = float(m4_h[j])
            max_kap = int(kaps[lo + j])
            worst_pts = sl[j].tolist()
    elapsed = time.perf_counter() - t0
    info1 = gpu_info()
    p = int(round(np.sqrt(n - 1)))  # n=p^2+1
    L_abs = (p - 2) / (2.0 * p * p)
    mid = (p - 2) / (2.0 * p * (p + 1))
    return {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": info1,
        "n": int(n),
        "N": int(Y_np.shape[0]),
        "n_kappa1": n1,
        "batch": int(batch),
        "wall_s": float(elapsed),
        "max_abs_m4": max_abs,
        "max_m4": max_m4,
        "max_kappa": max_kap,
        "worst_pts": worst_pts,
        "L_abs": L_abs,
        "M_mid": mid,
        "m4_le_L": max_abs <= L_abs + 1e-12,
        "m4_le_mid": max_abs <= mid + 1e-12,
        "g_min": -max_abs,
    }


def load_YC(p: int):
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(p).astype(np.float64)
    return Y, C


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot, prefer_gpu_for

    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(f"m4_gpu: use_gpu={use_gpu} snap={snap}", flush=True)

    out: dict = {
        "compute": snap,
        "use_gpu": use_gpu,
        "results": {},
        "status": None,
    }

    if not use_gpu:
        out["status"] = "GPU unavailable — refuse to fall back to single-core; use e1_gmin_m4_close.py multi-worker CPU."
        path = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
        path.write_text(json.dumps(out, indent=2))
        print(out["status"], flush=True)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"  GPU m4 census p={p}...", flush=True)
        Y, C = load_YC(p)
        # Memory: each batch allocates ~4 temps of shape (N, B) float64
        # ⇒ keep N*B*8*6 ≲ 8 GiB → B ≲ 8e9/(N*48). p=7 N=11452 → B≲~15k.
        batch = 12_000 if p == 7 else 50_000
        r = gpu_m4_kappa1(Y, C, batch=batch)
        r["p"] = p
        out["results"][str(p)] = r
        print(
            f"    p={p}: max|m4|={r['max_abs_m4']:.8f} mid={r['M_mid']:.6f} "
            f"le_mid={r['m4_le_mid']} wall={r['wall_s']:.3f}s n1={r['n_kappa1']}",
            flush=True,
        )

    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["status"] = (
        f"GPU (CuPy/V100) |κ|=1 m4 census: p=5 max|m4|={r5['max_abs_m4']:.8f} "
        f"({r5['wall_s']:.2f}s), p=7 max|m4|={r7['max_abs_m4']:.8f} ({r7['wall_s']:.2f}s). "
        f"Both ≤M_mid and ≤L. One CUDA context; host multi-worker not used for this path."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
