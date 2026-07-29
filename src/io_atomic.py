"""Atomic / high-throughput I/O helpers (Wieferich-style).

Patterns that buy orders of magnitude when jobs allow:

1. **mmap reads** — many ProcessPool workers open the same .npy with
   ``mmap_mode='r'``; page cache is shared (no 86× full copies of Max+).
2. **atomic writes** — write to ``path.tmp`` then ``os.replace`` (crash-safe
   checkpoints / evidence; concurrent readers never see partial JSON).
3. **GPU device reductions** — keep max/sum on-device (CuPy ``cp.max`` /
   custom atomic kernels) instead of shipping full batch arrays host-side.

Not every path needs these; use when:
  - multiple workers re-read the same large array, or
  - GPU batches are large enough that H2D/D2H dominates, or
  - evidence/checkpoints are written under parallel jobs.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


def write_bytes_atomic(path: str | Path, data: bytes, *, fsync: bool = True) -> None:
    """Write bytes atomically: tmp in same dir → fsync → os.replace."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            if fsync:
                os.fsync(f.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_text_atomic(path: str | Path, text: str, *, encoding: str = "utf-8") -> None:
    write_bytes_atomic(path, text.encode(encoding))


def _json_default(obj: Any) -> Any:
    """Numpy / Path → JSON-safe scalars."""
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def write_json_atomic(path: str | Path, obj: Any, *, indent: int = 2) -> None:
    write_text_atomic(
        path, json.dumps(obj, indent=indent, default=_json_default) + "\n"
    )


def write_npy_atomic(path: str | Path, arr: np.ndarray) -> None:
    """Atomic np.save via temp file in same directory."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Use a real .npy suffix so numpy does not append another .npy
    fd, tmp = tempfile.mkstemp(prefix=path.stem + ".", suffix=".npy", dir=path.parent)
    os.close(fd)
    tmp_path = Path(tmp)
    try:
        with open(tmp_path, "wb") as f:
            np.save(f, arr)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, path)
    except Exception:
        try:
            tmp_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise


def load_npy(
    path: str | Path,
    *,
    mmap_mode: str | None = "r",
    dtype: np.dtype | None = None,
) -> np.ndarray:
    """
    Load .npy with optional mmap (default read-only share-friendly).

    Use mmap_mode='r' for multi-worker ProcessPool consumers of the same Max+
    cache. Use mmap_mode=None when you need a writable private copy or will
    immediately push the whole array to GPU.
    """
    path = Path(path)
    arr = np.load(path, mmap_mode=mmap_mode)
    if dtype is not None and arr.dtype != dtype:
        # mmap cannot change dtype in place — force copy
        arr = np.asarray(arr, dtype=dtype)
    return arr


def load_maxplus_array(p: int, *, mmap: bool = True) -> np.ndarray:
    """Standard Max+ cache locations with mmap by default."""
    candidates = {
        5: [Path("/tmp/maxplus_p5.npy"), Path("/tmp/e1_p5/maxplus.npy")],
        7: [Path("/tmp/e1_p7/maxplus.npy"), Path("/tmp/maxplus_p7.npy")],
    }
    if p == 3:
        # often generated in-process; no default mmap path
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)

    for path in candidates.get(p, []):
        if path.is_file():
            return load_npy(path, mmap_mode="r" if mmap else None, dtype=None)
    raise FileNotFoundError(f"Max+ cache missing for p={p}")


def gpu_max_abs_m4_batches(
    Y_np: np.ndarray,
    idx_ab: np.ndarray,
    idx_cd: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None,
    *,
    batch: int = 12_000,
) -> tuple[float, float]:
    """
    On-device max |m4| over 4-set index rows without shipping every m4 to host.

    idx: (B,4) int32 of 4-set indices, processed in batches.
    Returns (max_abs_m4, signed_m4_at_max).
    """
    import cupy as cp

    if idx_cd is not None:
        raise ValueError("pass idx as (n1,4) array only")
    Y = cp.asarray(np.ascontiguousarray(Y_np, dtype=np.float64))
    idx = np.asarray(idx_ab, dtype=np.int32)
    n1 = idx.shape[0]
    max_abs = cp.float64(0.0)
    max_signed = cp.float64(0.0)
    for lo in range(0, n1, batch):
        hi = min(n1, lo + batch)
        sl = cp.asarray(idx[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4 = cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0)
        # device reduction — no full host transfer of the batch
        j = int(cp.argmax(cp.abs(m4)).get())
        am = float(cp.abs(m4[j]).get())
        if am > float(max_abs):
            max_abs = am
            max_signed = float(m4[j].get())
    return float(max_abs), float(max_signed)
