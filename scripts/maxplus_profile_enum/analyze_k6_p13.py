#!/usr/bin/env python3
"""Post-enum: wait is external. Concat k=6 orbs, unique, Cy=py, call moments.py."""
from __future__ import annotations

import glob
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/mnt/storage/e1work/maxplus_p13")
LOG = ROOT / "enum_p13_k6.log"
OUT = ROOT / "k6_p13_full.npy"
TOOL = Path(__file__).resolve().parent


def log_complete(text: str) -> bool:
    return "[17805/17805]" in text and "k=6 gauged GPU TOTAL" in text


def concat_unique() -> tuple[int, int]:
    paths = sorted(glob.glob(str(ROOT / "k6_gpu_out" / "orb*.npy")))
    parts = []
    raw = 0
    for p in paths:
        a = np.load(p)
        raw += len(a)
        if len(a):
            parts.append(a)
    if not parts:
        raise SystemExit("no nonempty orbs")
    A = np.concatenate(parts, axis=0)
    view = A.reshape(A.shape[0], -1).view(
        np.dtype((np.void, A.dtype.itemsize * A.shape[1]))
    )
    _, idx = np.unique(view, return_index=True)
    B = A[np.sort(idx)].astype(np.int8)
    np.save(OUT, B)
    with open(OUT, "rb") as fh:
        os.fsync(fh.fileno())
    return raw, int(len(B))


def eigencheck(Y: np.ndarray, p: int = 13) -> int:
    sys.path.insert(0, str(TOOL.parents[1] / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.int64)
    # Y is finite q; prepend +1 at infinity (eps=+1)
    q = p * p
    ones = np.ones((len(Y), 1), dtype=np.int64)
    Yn = np.concatenate([ones, Y.astype(np.int64)], axis=1)
    fail = 0
    bs = 4096
    for i in range(0, len(Yn), bs):
        chunk = Yn[i : i + bs]
        got = chunk @ C.T
        if np.any(got != p * chunk):
            fail += int(np.sum(np.any(got != p * chunk, axis=1)))
    return fail


def main() -> None:
    text = LOG.read_text()
    if not log_complete(text):
        print("NOT FINISHED")
        sys.exit(2)
    raw, dist = concat_unique()
    Y = np.load(OUT)
    fail = eigencheck(Y)
    print(f"raw={raw} distinct={dist} eigen_fail={fail} wrote {OUT}")
    if fail:
        sys.exit(1)
    subprocess.check_call(
        [sys.executable, "-u", str(TOOL / "moments.py"), "13", str(OUT)]
    )


if __name__ == "__main__":
    main()
