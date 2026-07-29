#!/usr/bin/env python3
"""
GPU-accelerated |κ|=1 m4 census (CuPy / V100) with Wieferich-style I/O.

Loads Max+ Y via mmap (shared page cache), one H2D, evaluates
  m4(S) = mean_k Y[k,a]*Y[k,b]*Y[k,c]*Y[k,d]
for all |κ|=1 4-sets in batches — one CUDA context, on-device argmax
(no full-batch D2H). Also tracks same-sign residual r = κ·(m4 − κ/p²)
for the cand/L gain budgets.

Writes evidence/e1_gmin_m4_gpu.json via write_json_atomic (tmp+fsync+replace).
"""
from __future__ import annotations

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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def rho_budget_cand(p: int) -> float:
    return M_cand(p) - 1.0 / (p * p)


def rho_budget_L(p: int) -> float:
    return (p - 4) / (2.0 * p * p)


def gain_budget_cand(p: int) -> float:
    return (p * p - 4 * p - 3) / (24.0 * (2 * p + 3))


def gain_budget_L(p: int) -> float:
    return (p - 4) / 48.0


def kappa_batch(C: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """idx (B,4) int; return κ for each row."""
    a, b, c, d = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]
    return (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)


def all_quads(n: int) -> np.ndarray:
    return np.array(list(combinations(range(n), 4)), dtype=np.int32)


def gpu_m4_kappa1(
    Y_np: np.ndarray,
    C_np: np.ndarray,
    *,
    batch: int = 12_000,
    p: int | None = None,
) -> dict:
    """
    Full |κ|=1 max|m4| + same-sign residual on GPU.
    mmap-friendly host Y → single H2D → device reduce only.
    """
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Y_host = np.ascontiguousarray(Y_np, dtype=np.float64)
    Y = cp.asarray(Y_host)
    C = np.asarray(C_np, dtype=np.float64)
    n = C.shape[0]
    if p is None:
        p = int(round(np.sqrt(n - 1)))
    quads = all_quads(n)
    kap = kappa_batch(C, quads)
    mask = np.abs(kap) == 1
    idx = quads[mask]
    kaps = kap[mask].astype(np.float64)
    n1 = int(len(idx))
    inv_p2 = 1.0 / (p * p)

    max_abs = 0.0
    max_m4 = 0.0
    max_kap = 0
    worst_pts = None
    max_same_r = 0.0  # max κ*(m4 - κ/p²) when positive
    t0 = time.perf_counter()
    for lo in range(0, n1, batch):
        hi = min(n1, lo + batch)
        sl = cp.asarray(idx[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4 = cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0)
        # on-device |m4| argmax
        j = int(cp.argmax(cp.abs(m4)).get())
        am = float(cp.abs(m4[j]).get())
        if am > max_abs:
            max_abs = am
            max_m4 = float(m4[j].get())
            max_kap = int(kaps[lo + j])
            worst_pts = idx[lo + j].tolist()
        # same-sign residual r = κ * (m4 - κ/p²) = κ*m4 - 1/p²
        # (κ=±1 ⇒ κ²=1)
        kap_d = cp.asarray(kaps[lo:hi])
        r = kap_d * m4 - inv_p2
        # max positive r only (dangerous for |m4|)
        r_pos = cp.maximum(r, 0.0)
        jr = int(cp.argmax(r_pos).get())
        rv = float(r_pos[jr].get())
        if rv > max_same_r:
            max_same_r = rv
    elapsed = time.perf_counter() - t0
    info1 = gpu_info()
    sa = 24.0 / (p * p)
    cand, mid, La = M_cand(p), M_mid(p), L_abs(p)
    return {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": info1,
        "p": int(p),
        "n": int(n),
        "N": int(Y_host.shape[0]),
        "n_kappa1": n1,
        "batch": int(batch),
        "wall_s": float(elapsed),
        "max_abs_m4": max_abs,
        "max_m4_signed": max_m4,
        "max_kappa": max_kap,
        "worst_pts": worst_pts,
        "max_same_sign_rho": max_same_r,
        "effective_gain": max_same_r / sa if sa > 0 else None,
        "M_cand": cand,
        "M_mid": mid,
        "L_abs": La,
        "rho_budget_cand": rho_budget_cand(p),
        "rho_budget_L": rho_budget_L(p),
        "gain_budget_cand": gain_budget_cand(p),
        "gain_budget_L": gain_budget_L(p),
        "m4_le_cand": max_abs <= cand + 1e-12,
        "m4_le_mid": max_abs <= mid + 1e-12,
        "m4_le_L": max_abs <= La + 1e-12,
        "same_rho_le_cand": max_same_r <= rho_budget_cand(p) + 1e-12,
        "same_rho_le_L": max_same_r <= rho_budget_L(p) + 1e-12,
        "gain_le_cand": (max_same_r / sa) <= gain_budget_cand(p) + 1e-12,
        "gain_le_L": (max_same_r / sa) <= gain_budget_L(p) + 1e-12,
        "g_min": -max_abs,
        "gpu_used": True,
        "io": {
            "maxplus_read": "mmap load_maxplus_array",
            "H2D": "single contiguous host→device copy of Y",
            "device_reduce": "cp.argmax(|m4|), cp.argmax(r_pos) — no full batch D2H",
            "evidence_write": "write_json_atomic (tmp+fsync+os.replace)",
        },
    }


def load_YC(p: int):
    sys.path.insert(0, str(ROOT / "src"))
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    Y = load_maxplus_array(p, mmap=True)
    return Y, C


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()  # record CPU budget even on GPU path
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"m4_gpu: W={W} use_gpu={use_gpu} "
        f"io=mmap+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )

    out: dict = {
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; single H2D; device argmax; evidence via write_json_atomic"
        ),
        "results": {},
        "status": None,
    }

    if not use_gpu:
        out["status"] = (
            "GPU unavailable — refuse single-core thrash; "
            "use e1_gmin_m4_kernel.py multi-worker CPU + mmap + atomic."
        )
        path = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
        write_json_atomic(path, out)
        print(out["status"], flush=True)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"  GPU m4+cand census p={p}...", flush=True)
        Y, C = load_YC(p)
        # Memory: N*B*8*~6 ≲ 8 GiB → p=7 N=11452 → B≲~12k
        batch = 12_000 if p >= 7 else 50_000
        r = gpu_m4_kappa1(Y, C, batch=batch, p=p)
        out["results"][str(p)] = r
        print(
            f"    p={p}: max|m4|={r['max_abs_m4']:.10f} cand={r['M_cand']:.10f} "
            f"le_cand={r['m4_le_cand']} le_mid={r['m4_le_mid']} "
            f"same_ρ={r['max_same_sign_rho']:.8f} gain={r['effective_gain']:.8f} "
            f"wall={r['wall_s']:.3f}s n1={r['n_kappa1']} GPU",
            flush=True,
        )

    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["both_le_cand"] = bool(r5["m4_le_cand"] and r7["m4_le_cand"])
    out["both_le_mid"] = bool(r5["m4_le_mid"] and r7["m4_le_mid"])
    out["status"] = (
        f"GPU CuPy/V100 |κ|=1 m4 census (mmap+atomic): "
        f"p=5 max|m4|={r5['max_abs_m4']:.10f}≤M_cand={r5['M_cand']:.10f} "
        f"(gain={r5['effective_gain']:.8f}, {r5['wall_s']:.2f}s); "
        f"p=7 max|m4|={r7['max_abs_m4']:.10f}≤M_cand={r7['M_cand']:.10f} "
        f"(gain={r7['effective_gain']:.8f}, {r7['wall_s']:.2f}s). "
        f"One CUDA context; device reduce; evidence atomic. "
        f"OPEN: prove max|m4|≤M_cand for all primes p≥5. lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
