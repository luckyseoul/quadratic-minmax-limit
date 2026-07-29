#!/usr/bin/env python3
"""
p=7 finer class probe + Max+ m4 census on |kappa|=1 (multi-worker + mmap).

Coarse (type6, ext_sum_hist) is NOT m4-constant at p=7 (max_std ~ 0.012).
This adds external sign-signature histogram and measures true max |m4|.

F17: ProcessPool W≈nproc-2; mmap Max+; atomic evidence JSON.
"""
from __future__ import annotations

import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_k] = "1"

_G: dict = {}


def _init(C: np.ndarray, mp_path: str) -> None:
    _G["C"] = C
    _G["Mp"] = np.load(mp_path, mmap_mode="r")
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def finer_key(C: np.ndarray, pts: tuple) -> tuple:
    from e1_gmin_moduli import type6, ext_sum_hist, kappa

    t6 = type6(C, pts)
    eh = ext_sum_hist(C, pts)
    kap = kappa(C, pts)
    pts_l = list(pts)
    S = set(pts_l)
    n = C.shape[0]
    ctr: Counter = Counter()
    for r in range(n):
        if r in S:
            continue
        sig = tuple(sorted(int(C[r, v]) for v in pts_l))
        ctr[sig] += 1
    ext_sig = tuple(sorted(ctr.items()))
    return (t6, eh, kap, ext_sig)


def _shard(quads: list) -> list:
    from e1_gmin_moduli import kappa

    C = _G["C"]
    Mp = _G["Mp"]
    out = []
    for q in quads:
        k = finer_key(C, q)
        a, b, c, d = q
        m4 = float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
        out.append((k, q, m4, kappa(C, q)))
    return out


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


def main() -> None:
    from e1_gmin_cbound_gen import trG2_gpu_allpairs
    from e1_gmin_moduli import L_p, T_p
    from io_atomic import load_maxplus_array, write_json_atomic
    from minmax_quadratic import paley_conference_prime_power
    from workers import require_workers

    W = require_workers()
    p = 7
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    mp_path = "/tmp/e1_p7/maxplus.npy"
    Mp = np.ascontiguousarray(load_maxplus_array(p, mmap=True), dtype=float)
    assert Path(mp_path).is_file()

    quads = list(combinations(range(n), 4))
    nq = len(quads)
    n_shards = min(W, max(1, nq // 500))
    ss = (nq + n_shards - 1) // n_shards
    shards = [quads[i : i + ss] for i in range(0, nq, ss)]
    print(
        f"p=7 n={n} N={Mp.shape[0]} quads={nq} shards={len(shards)} W={W}",
        flush=True,
    )

    t0 = time.perf_counter()
    classes: dict = defaultdict(list)
    m4s: dict = defaultdict(list)
    kaps: dict = {}
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init,
        initargs=(C, mp_path),
    ) as ex:
        futs = [ex.submit(_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            for k, q, m4, kap in fut.result():
                classes[k].append(q)
                m4s[k].append(m4)
                kaps[k] = kap
    wall = time.perf_counter() - t0

    n_cls = len(classes)
    stds = [float(np.std(v)) if len(v) > 1 else 0.0 for v in m4s.values()]
    n_const = sum(1 for s in stds if s < 1e-9)

    max_abs = 0.0
    n_kap1 = 0
    for k, vals in m4s.items():
        if abs(kaps[k]) != 1:
            continue
        n_kap1 += len(vals)
        max_abs = max(max_abs, max(abs(x) for x in vals))

    mid = M_mid(p)
    tr = trG2_gpu_allpairs(Mp)

    out = {
        "p": 7,
        "workers": W,
        "n_classes_finer": n_cls,
        "n_constant": n_const,
        "all_constant": n_const == n_cls,
        "max_std": max(stds) if stds else None,
        "wall_class_s": wall,
        "n_kappa1_quads": n_kap1,
        "max_abs_m4_kappa1": max_abs,
        "M_mid": mid,
        "L_abs": abs(L_p(p)),
        "T_abs": abs(T_p(p)),
        "max_abs_le_M_mid": max_abs <= mid + 1e-12,
        "max_abs_le_L": max_abs <= abs(L_p(p)) + 1e-12,
        "g_min_empirical": -max_abs,
        "Tr": tr,
        "io": "mmap Max+ in workers; write_json_atomic",
        "note": (
            "Finer key still may not be fully constant; empirical max|m4| on |kappa|=1 "
            "is the durable census for bi-tight at p=7."
        ),
    }
    print(
        f"classes={n_cls} const={n_const}/{n_cls} max_std={out['max_std']} "
        f"max|m4|={max_abs} mid={mid} le_mid={out['max_abs_le_M_mid']} "
        f"wall={wall:.1f}s Tr={tr['Tr_G2']:.4f}",
        flush=True,
    )
    path = ROOT / "evidence" / "e1_gmin_p7_finer_classes.json"
    write_json_atomic(path, out)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
