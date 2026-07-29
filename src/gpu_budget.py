"""GPU + CPU compute budget for this machine (V100 + multi-core).

Policy:
  - Always inventory free CPU *and* GPU before multi-minute work.
  - Prefer GPU for large dense batch matmuls / column products / sampling Φ.
  - Prefer ProcessPool (W≈nproc-2, OMP=1/worker) for independent irregular units
    (seeds, primes, combination shards that don't fit GPU well).
  - Do NOT spawn 88 CUDA contexts (one process owns the GPU).
  - Report workers + whether GPU was used in evidence JSON.

Fits GPU (CuPy on Tesla V100-SXM2-16GB):
  - Batch m4: mean(Y[:,a]*Y[:,b]*Y[:,c]*Y[:,d]) over many 4-sets
  - Batch Φ / Rayleigh on many A or many x
  - Dense eig / power on moderate matrices that fit in 16GB

Keep on CPU multi-worker:
  - Sparse irregular graph ops on binom(n,4) (unless reformulated dense)
  - pytest-xdist, Aut orbit enumerations with branching
"""
from __future__ import annotations

import os
import sys
from typing import Any


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def gpu_info() -> dict[str, Any]:
    """Return GPU status; never raises."""
    out: dict[str, Any] = {
        "available": False,
        "backend": None,
        "name": None,
        "mem_free_gib": None,
        "mem_total_gib": None,
        "util_pct": None,
    }
    try:
        import cupy as cp

        free, total = cp.cuda.runtime.memGetInfo()
        out["available"] = True
        out["backend"] = "cupy"
        out["mem_free_gib"] = free / 1e9
        out["mem_total_gib"] = total / 1e9
        try:
            props = cp.cuda.runtime.getDeviceProperties(0)
            name = props.get("name")
            if isinstance(name, bytes):
                name = name.decode()
            out["name"] = name
        except Exception:
            out["name"] = "cuda:0"
        return out
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def require_gpu() -> Any:
    """Import cupy or abort with a clear message."""
    try:
        import cupy as cp

        return cp
    except ImportError as e:
        raise SystemExit(
            "gpu_budget: CuPy required for this path but not importable. "
            f"({e})"
        ) from e


def compute_snapshot(headroom: int = 2) -> dict[str, Any]:
    """CPU workers + GPU in one dict for evidence logs."""
    _blas1()
    sys.path.insert(0, os.path.dirname(__file__))
    from workers import cpu_count_reliable, default_workers

    n = cpu_count_reliable()
    w = default_workers(headroom=headroom)
    g = gpu_info()
    snap = {
        "cpu_threads": n,
        "workers": w,
        "gpu": g,
        "policy": (
            "Use GPU for large dense batches if available; "
            f"else ProcessPool W={w}. Never single-core multi-minute thrash (F17)."
        ),
    }
    print(
        f"compute_snapshot: W={w} cpus={n} gpu="
        f"{'yes '+str(g.get('name')) if g.get('available') else 'no'} "
        f"mem_free={g.get('mem_free_gib')}",
        file=sys.stderr,
        flush=True,
    )
    return snap


def prefer_gpu_for(workload: str) -> bool:
    """Heuristic: should this workload try GPU first?"""
    g = gpu_info()
    if not g.get("available"):
        return False
    # free < 1 GiB → skip
    free = g.get("mem_free_gib") or 0
    if free < 1.0:
        return False
    dense = {
        "m4_batch",
        "phi_batch",
        "rayleigh_batch",
        "dense_matmul",
        "power_iter_dense",
        "y_column_products",
    }
    return workload in dense
