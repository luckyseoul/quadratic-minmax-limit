"""k=6 mesh stop flags and orbit claim.

Soft-soft: a worker that already holds `.lock_<tvidx>` finishes that orbit
and writes `orb<tvidx>.npy`. It does not mkdir a new lock once a stop flag
is up. Start of one worker never clears another worker's flag.
"""
from __future__ import annotations

import os
import time
from pathlib import Path

DEFAULT_ROOT = Path(os.environ.get("K6_ROOT", "/mnt/storage/e1work/maxplus_p13"))
WORKERS = ("v100", "nuka", "orin", "a380", "cpu44", "dash")
COMPUTE_WORKERS = ("v100", "nuka", "orin", "a380", "cpu44")
STOP_ALL = "ALL"
STARVE_MARK = "STARVE"
N_TASKS_DEFAULT = 17805


def stop_dir(root: Path | None = None) -> Path:
    r = Path(root) if root is not None else Path(
        os.environ.get("K6_STOP_DIR", str(DEFAULT_ROOT / "mesh_stop"))
    )
    return r


def stop_path(name: str, root: Path | None = None) -> Path:
    return stop_dir(root) / name


def request_stop(name: str, root: Path | None = None) -> Path:
    d = stop_dir(root)
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    p.write_text(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\nsoft\n")
    return p


def clear_stop(name: str, root: Path | None = None) -> None:
    """Remove one flag. ALL removes only the ALL file, not per-worker flags."""
    p = stop_path(name, root)
    try:
        p.unlink()
    except FileNotFoundError:
        pass


def prepare_start(name: str, root: Path | None = None) -> None:
    """Make *name* runnable without un-stopping the rest of the mesh.

    If ALL is set, convert it into per-worker flags on every other compute
    worker, then drop ALL and *name*.
    """
    d = stop_dir(root)
    d.mkdir(parents=True, exist_ok=True)
    all_p = d / STOP_ALL
    if all_p.is_file():
        try:
            all_p.unlink()
        except FileNotFoundError:
            pass
        for w in COMPUTE_WORKERS:
            if w != name:
                request_stop(w, root)
    clear_stop(name, root)


def should_stop(host: str | None = None, root: Path | None = None) -> bool:
    """True if this host must not start another orbit."""
    d = stop_dir(root)
    if (d / STOP_ALL).is_file():
        return True
    h = host if host is not None else os.environ.get("K6_HOST", "")
    if h and (d / h).is_file():
        return True
    return False


def stop_flag_state(root: Path | None = None) -> dict[str, bool]:
    d = stop_dir(root)
    out = {STOP_ALL: (d / STOP_ALL).is_file()}
    for w in WORKERS:
        if w == "dash":
            out[w] = (d / w).is_file()
        else:
            out[w] = (d / w).is_file() or out[STOP_ALL]
    return out


def try_claim_orbit(
    outdir: str | Path,
    tvidx: int,
    host: str | None = None,
    root: Path | None = None,
) -> int:
    """0 = claimed this lock, -1 = already done/locked, -2 = soft-stop skip.

    Must not mkdir when a stop flag is up. In-flight holders never reach
    this call for the orbit they already hold.
    """
    outdir = Path(outdir)
    if (outdir / f"orb{tvidx}.npy").is_file():
        return -1
    if should_stop(host, root):
        return -2
    try:
        (outdir / f".lock_{tvidx}").mkdir()
    except FileExistsError:
        return -1
    return 0


def _lock_tvidx(path: Path) -> int | None:
    name = path.name
    if not name.startswith(".lock_"):
        return None
    try:
        return int(name.split("_", 1)[1])
    except ValueError:
        return None


def starve_unclaimed(
    outdir: str | Path,
    n_tasks: int = N_TASKS_DEFAULT,
) -> dict:
    """Occupy remaining slots with marked dummy locks.

    In-flight `.lock_*` (no STARVE file) and existing `orb*.npy` are left
    alone so current orbits finish. Used so *old* run_kgauged parents that
    predate should_stop() still will not start a new outer.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    have: set[int] = set()
    for f in outdir.glob("orb*.npy"):
        try:
            have.add(int(f.name[3:-4]))
        except ValueError:
            pass
    inflight: list[int] = []
    for d in outdir.glob(".lock_*"):
        i = _lock_tvidx(d)
        if i is None:
            continue
        if not (d / STARVE_MARK).is_file():
            inflight.append(i)
    n = 0
    for i in range(n_tasks):
        if i in have or i in inflight:
            continue
        path = outdir / f".lock_{i}"
        try:
            path.mkdir()
        except FileExistsError:
            continue
        (path / STARVE_MARK).write_text("dummy\n")
        n += 1
    inflight.sort()
    return {"dummy_locks": n, "inflight": inflight, "npy": len(have)}


def drop_starve_locks(outdir: str | Path) -> int:
    """Remove only STARVE-marked dummy locks. Never rmdir an in-flight lock."""
    outdir = Path(outdir)
    n = 0
    for d in outdir.glob(".lock_*"):
        mark = d / STARVE_MARK
        if not mark.is_file():
            continue
        try:
            mark.unlink()
        except FileNotFoundError:
            pass
        try:
            d.rmdir()
            n += 1
        except OSError:
            pass
    return n
