#!/usr/bin/env python3
"""
Prop 15.82: refine C-invariants so m4 is constant on classes for p≥7.

Compute policy (F17):
  - GPU all-m4 (mmap Max+, one CUDA context, H2D once)
  - ProcessPool W≈nproc-2 for class-key assignment + constancy
  - Atomic evidence JSON

Invariants tested (pure C / Paley):
  coarse = (type6, ext_sum_hist)   # Prop 15.53 — enough at p=5
  + CR (PGL-complete cross-ratio, e1_gmin_cr_classify.canon_cr_fn)
  + κ
  combinations thereof

Writes evidence/e1_gmin_m4_refine.json
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


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def gpu_all_m4(Y_np: np.ndarray, quads: np.ndarray, *, batch: int = 12_000):
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Y = cp.asarray(np.ascontiguousarray(Y_np, dtype=np.float64))
    m4 = np.empty(len(quads), dtype=np.float64)
    t0 = time.perf_counter()
    for lo in range(0, len(quads), batch):
        hi = min(len(quads), lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4[lo:hi] = cp.asnumpy(cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0))
    return m4, {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": gpu_info(),
        "wall_s": float(time.perf_counter() - t0),
        "n_quads": int(len(quads)),
        "batch": int(batch),
        "io": "mmap Max+; single H2D; D2H m4 only",
    }


# ── ProcessPool: assign class keys on shards ────────────────────────────────

_G: dict = {}


def _init_keys(C: np.ndarray, mode: str, p: int) -> None:
    _G["C"] = C
    _G["mode"] = mode
    _G["p"] = p
    _blas1()
    if "cr" in mode:
        from e1_gmin_cr_classify import canon_cr_fn

        _G["cr_fn"] = canon_cr_fn(p)
    from e1_gmin_moduli import type6, ext_sum_hist, kappa

    _G["type6"] = type6
    _G["ext_sum_hist"] = ext_sum_hist
    _G["kappa"] = kappa


def _key_of(pts: tuple) -> tuple:
    mode = _G["mode"]
    C = _G["C"]
    t6 = _G["type6"](C, pts)
    ext = _G["ext_sum_hist"](C, pts)
    kap = _G["kappa"](C, pts)
    if mode == "coarse":
        return (t6, ext)
    if mode == "coarse+kappa":
        return (t6, ext, kap)
    if mode == "cr":
        return (_G["cr_fn"](pts),)
    if mode == "cr+kappa":
        return (_G["cr_fn"](pts), kap)
    if mode == "type6+cr":
        return (t6, _G["cr_fn"](pts))
    if mode == "type6+cr+kappa":
        return (t6, _G["cr_fn"](pts), kap)
    if mode == "coarse+cr":
        return (t6, ext, _G["cr_fn"](pts))
    if mode == "coarse+cr+kappa":
        return (t6, ext, _G["cr_fn"](pts), kap)
    raise ValueError(mode)


def _shard_keys(args: tuple) -> list:
    """args = (lo, hi, quads_slice as list of tuples) -> list of (i, key)"""
    lo, quads_slice = args
    out = []
    for j, pts in enumerate(quads_slice):
        out.append((lo + j, _key_of(tuple(pts))))
    return out


def assign_keys_parallel(
    C: np.ndarray, quads: np.ndarray, mode: str, p: int, W: int
) -> list:
    """Return list key_of[i] for each quad index i."""
    nq = len(quads)
    # shard as lists of tuples for pickling
    nsh = min(W, max(1, nq // 500))
    ss = (nq + nsh - 1) // nsh
    shards = []
    for s in range(0, nq, ss):
        sl = quads[s : s + ss]
        shards.append((s, [tuple(map(int, row)) for row in sl]))
    keys: list = [None] * nq
    t0 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_keys,
        initargs=(C, mode, p),
    ) as ex:
        futs = [ex.submit(_shard_keys, sh) for sh in shards]
        for fut in as_completed(futs):
            for i, k in fut.result():
                keys[i] = k
    wall = time.perf_counter() - t0
    return keys, wall


def constancy_from_keys(keys: list, m4: np.ndarray) -> dict:
    """Group by key; compute m4 constancy (full arrays, host, fast)."""
    buckets: dict = defaultdict(list)
    for i, k in enumerate(keys):
        buckets[k].append(i)
    n_const = 0
    max_std = 0.0
    n_vals_max = 1
    nonconst_examples = []
    for k, idxs in buckets.items():
        vals = m4[np.array(idxs, dtype=np.int64)]
        std = float(np.std(vals)) if len(vals) > 1 else 0.0
        max_std = max(max_std, std)
        nuniq = len(set(np.round(vals, 10)))
        n_vals_max = max(n_vals_max, nuniq)
        if std < 1e-9:
            n_const += 1
        else:
            ctr = Counter(np.round(vals, 10))
            if len(nonconst_examples) < 6:
                nonconst_examples.append(
                    {
                        "n": len(idxs),
                        "std": std,
                        "n_distinct_m4": nuniq,
                        "m4_counts": {str(a): int(b) for a, b in ctr.most_common(4)},
                    }
                )
    return {
        "n_classes": len(buckets),
        "n_constant": n_const,
        "n_nonconstant": len(buckets) - n_const,
        "all_constant": n_const == len(buckets),
        "max_std": max_std,
        "max_distinct_m4_in_class": n_vals_max,
        "nonconst_examples": nonconst_examples,
    }


def max_abs_m4_kappa1(C, quads, m4, kap) -> float:
    mask = np.abs(kap) == 1
    return float(np.max(np.abs(m4[mask]))) if np.any(mask) else 0.0


def run_p(p: int, W: int) -> dict:
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power
    from gpu_budget import prefer_gpu_for

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = load_maxplus_array(p, mmap=True)
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    a, b, c, d = quads.T
    kap = (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)

    use_gpu = prefer_gpu_for("m4_batch")
    if not use_gpu:
        raise SystemExit("GPU required for m4 (F17)")
    batch = 12_000 if p >= 7 else 50_000
    print(f"  GPU all-m4 p={p} n_quads={len(quads)}...", flush=True)
    m4, gpu_meta = gpu_all_m4(Y, quads, batch=batch)
    print(f"    wall={gpu_meta['wall_s']:.3f}s", flush=True)

    modes = [
        "coarse",
        "coarse+kappa",
        "cr",
        "cr+kappa",
        "type6+cr",
        "type6+cr+kappa",
        "coarse+cr",
        "coarse+cr+kappa",
    ]
    strat = {}
    best = None
    for mode in modes:
        print(f"  stratify mode={mode} W={W}...", flush=True)
        keys, wall = assign_keys_parallel(C, quads, mode, p, W)
        cst = constancy_from_keys(keys, m4)
        cst["wall_key_s"] = float(wall)
        cst["mode"] = mode
        strat[mode] = cst
        print(
            f"    classes={cst['n_classes']} const={cst['n_constant']}/"
            f"{cst['n_classes']} max_std={cst['max_std']:.8f} "
            f"key_wall={wall:.2f}s",
            flush=True,
        )
        if cst["all_constant"] and (
            best is None or cst["n_classes"] < best["n_classes"]
        ):
            best = {"mode": mode, **cst}

    return {
        "p": int(p),
        "n": int(n),
        "N": int(np.asarray(Y).shape[0]),
        "n_quads": int(len(quads)),
        "gpu_m4": gpu_meta,
        "gpu_used": True,
        "max_abs_m4_kappa1": max_abs_m4_kappa1(C, quads, m4, kap),
        "M_cand": float(M_cand(p)),
        "m4_le_cand": max_abs_m4_kappa1(C, quads, m4, kap) <= M_cand(p) + 1e-12,
        "strategies": strat,
        "best_constant_strategy": best,
        "workers": int(W),
        "io": "mmap Max+; GPU all-m4; ProcessPool class keys; write_json_atomic",
    }


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"m4_refine Prop15.82: W={W} use_gpu={use_gpu} "
        f"io=mmap+GPU+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )
    out: dict = {
        "prop": "15.82",
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; GPU all-m4 (one CUDA context); "
            "ProcessPool class-key shards; write_json_atomic"
        ),
        "results": {},
        "status": None,
    }
    if not use_gpu:
        out["status"] = "GPU unavailable — refuse (F17)."
        write_json_atomic(ROOT / "evidence" / "e1_gmin_m4_refine.json", out)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"=== p={p} ===", flush=True)
        out["results"][str(p)] = run_p(p, W)

    r5, r7 = out["results"]["5"], out["results"]["7"]
    # p=5 coarse must stay constant; p=7 seek any all_constant strategy
    p5_ok = r5["strategies"]["coarse"]["all_constant"]
    p7_ok = any(s["all_constant"] for s in r7["strategies"].values())
    out["p5_coarse_constant"] = p5_ok
    out["p7_some_strategy_constant"] = p7_ok
    out["p7_best"] = r7.get("best_constant_strategy")
    out["status"] = (
        f"Prop 15.82: GPU+multi-worker class refinement. "
        f"p=5 coarse constant={p5_ok} ({r5['strategies']['coarse']['n_classes']} classes). "
        f"p=7 some strategy all-constant={p7_ok}"
        + (
            f" best={out['p7_best']['mode']} n_classes={out['p7_best']['n_classes']}"
            if out["p7_best"]
            else f" (best max_std among tries="
            f"{min(s['max_std'] for s in r7['strategies'].values()):.6f})"
        )
        + f". m4_le_cand p5={r5['m4_le_cand']} p7={r7['m4_le_cand']}. "
        f"OPEN: if p7 constant, run moduli c*≤c_GD; else find finer C/Aut invariant. "
        f"lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_refine.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
