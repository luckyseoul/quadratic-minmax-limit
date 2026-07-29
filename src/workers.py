"""F17 multi-worker policy: never thrash a single core on multi-minute jobs.

Use require_workers() at the start of every heavy script/main.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor


def cpu_count_reliable() -> int:
    """Largest reliable CPU count (cgroup nproc can lie low)."""
    cands = [
        os.cpu_count() or 0,
    ]
    try:
        cands.append(len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    try:
        with open("/proc/cpuinfo") as f:
            cands.append(sum(1 for line in f if line.startswith("processor")))
    except OSError:
        pass
    return max(cands) if cands else 1


def default_workers(headroom: int = 2) -> int:
    n = cpu_count_reliable()
    return max(2, n - headroom)


def require_workers(min_workers: int = 4, headroom: int = 2) -> int:
    """
    Return worker count W = max(2, ncpu - headroom).
    Abort if the machine is too small to justify multi-worker policy (F17).
    Also force BLAS single-thread per process so ProcessPool does not oversubscribe.
    """
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[k] = "1"
    n = cpu_count_reliable()
    if n < min_workers:
        raise SystemExit(
            f"workers.py F17 FATAL: only {n} CPUs visible; refusing single-core heavy job "
            f"(need ≥{min_workers})"
        )
    w = max(2, n - headroom)
    # Allow explicit override only if still multi-worker
    env = os.environ.get("GROK_WORKERS") or os.environ.get("PYTEST_WORKERS")
    if env is not None:
        w = max(2, int(env))
    if w < 2:
        raise SystemExit("workers.py F17 FATAL: W < 2")
    return w


def pool(max_workers: int | None = None) -> ProcessPoolExecutor:
    w = max_workers if max_workers is not None else require_workers()
    return ProcessPoolExecutor(max_workers=w)
