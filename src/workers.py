"""F17 multi-worker + GPU policy: use the whole machine.

Hard rules (project + user-enforced):
  1. Call require_workers() / compute_snapshot() at the start of every heavy main.
  2. Fan out expensive *independent* units (quads, primes, seeds) with ProcessPool
     W≈nproc-2 and OMP=1 per worker — not toy tasks only.
  3. Prefer GPU (CuPy/V100) for large dense batches (m4 column products, Φ sampling,
     dense matmul). One CUDA context — do not spawn 88 GPU processes.
  4. NEVER serial leftover heavy work on the parent after a ProcessPool wave
     (classic bug: pool for p=3,5 then p=7 on main → 97% NLWP=1).
  5. After launch, self-check with assert_not_single_core_thrash() or ps;
     nvidia-smi if the path claimed GPU.

See evidence/E1_FAILURE_GRAPH.md F17, evidence/AGENT_BUG_F17_RECURRENCE.md,
and src/gpu_budget.py.
"""
from __future__ import annotations

import os
import subprocess
import sys
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
    Force BLAS single-thread per process so ProcessPool does not oversubscribe.
    """
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"
    n = cpu_count_reliable()
    if n < min_workers:
        raise SystemExit(
            f"workers.py F17 FATAL: only {n} CPUs visible; refusing single-core heavy job "
            f"(need ≥{min_workers})"
        )
    w = max(2, n - headroom)
    env = os.environ.get("GROK_WORKERS") or os.environ.get("PYTEST_WORKERS")
    if env is not None:
        w = max(2, int(env))
    if w < 2:
        raise SystemExit("workers.py F17 FATAL: W < 2")
    print(
        f"workers.py: F17 W={w} (cpus={n}, headroom={headroom}); "
        f"BLAS threads=1 per process",
        file=sys.stderr,
        flush=True,
    )
    return w


def pool(max_workers: int | None = None) -> ProcessPoolExecutor:
    w = max_workers if max_workers is not None else require_workers()
    if w < 2:
        raise SystemExit("workers.py F17 FATAL: refusing ProcessPool with W<2")
    return ProcessPoolExecutor(max_workers=w)


def assert_not_single_core_thrash(
    *,
    min_cpus: int = 8,
    pcpu_thresh: float = 80.0,
    our_pids: set[int] | None = None,
) -> None:
    """
    Best-effort live check: if this host has many CPUs and the hottest process is
    a single-threaded python at ~100% CPU, abort (F17).

    Call periodically from long jobs, or right after launching a background wave.
    """
    n = cpu_count_reliable()
    if n < min_cpus:
        return
    try:
        out = subprocess.check_output(
            ["ps", "-eo", "pid,pcpu,nlwp,comm,args", "--sort=-pcpu"],
            text=True,
            timeout=5,
        )
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return
    lines = out.strip().splitlines()
    if len(lines) < 2:
        return
    # skip header
    for line in lines[1:6]:
        parts = line.split(None, 4)
        if len(parts) < 5:
            continue
        pid_s, pcpu_s, nlwp_s, comm, args = parts
        try:
            pid = int(pid_s)
            pcpu = float(pcpu_s)
            nlwp = int(nlwp_s)
        except ValueError:
            continue
        if our_pids is not None and pid not in our_pids:
            # still flag any thrashing python owned by us via args path
            pass
        if pcpu < pcpu_thresh or nlwp != 1:
            continue
        if "python" not in comm.lower() and "python" not in args.lower():
            continue
        # ignore short-lived ps/self
        if "workers.py" in args or "assert_not_single_core" in args:
            continue
        raise SystemExit(
            f"workers.py F17 FATAL: single-core thrash detected: "
            f"pid={pid} pcpu={pcpu} nlwp={nlwp} cmd={args[:160]!r}\n"
            f"  Host has {n} CPUs. Rewrite to ProcessPool/xdist shards "
            f"(see evidence/E1_FAILURE_GRAPH.md F17)."
        )


def forbid_serial_heavy_on_main(stage: str) -> None:
    """
    Call before any parent-process heavy loop that is not sharded.
    Documents intent and aborts if GROK_ALLOW_SERIAL is not set.
    """
    if os.environ.get("GROK_ALLOW_SERIAL") == "1":
        print(
            f"workers.py: SERIAL ALLOWED by GROK_ALLOW_SERIAL=1 at stage={stage!r}",
            file=sys.stderr,
            flush=True,
        )
        return
    n = cpu_count_reliable()
    if n >= 8:
        raise SystemExit(
            f"workers.py F17 FATAL: refused serial heavy stage on main "
            f"({stage!r}) with {n} CPUs. Shard it or set GROK_ALLOW_SERIAL=1 "
            f"only if the algorithm is inherently serial (and log why)."
        )
