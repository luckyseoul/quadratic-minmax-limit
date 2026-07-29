#!/usr/bin/env python3
"""
S1 sign cert at star_a=+1 on same-sign residual |κ|=1 sets.

Lead identity (Prop 15.76):
  p ρ − 2 star_a / p² = S1 + S3
If S1 ≤ 0 whenever star_a = +1 on r = κ·ρ > 0 sets, then
  p ρ ≤ 2/p² + S3  (star=+1),
and the cand attack reduces to bounding S3.

Compute policy (F17, non-negotiable):
  - Max+ via mmap (shared page cache)
  - All m4 on GPU (one CUDA context, batched column products, device reduce
    where useful; single H2D of Y; D2H only the compact m4 vector)
  - Structure walk (star / S1 / S3) on host ProcessPool over precomputed m4
  - Evidence via write_json_atomic (tmp + fsync + os.replace)

Writes evidence/e1_gmin_m4_S1_sign.json
"""
from __future__ import annotations

import os
import sys
import time
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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def rho_budget_cand(p: int) -> float:
    return M_cand(p) - 1.0 / (p * p)


def kappa_all(C: np.ndarray, quads: np.ndarray) -> np.ndarray:
    a, b, c, d = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]
    return (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)


def star_a(C: np.ndarray, S: tuple[int, ...], v: int) -> int:
    pr = 1
    for u in S:
        if u != v:
            pr *= int(C[v, u])
    return pr


def gpu_all_m4(
    Y_np: np.ndarray,
    quads: np.ndarray,
    *,
    batch: int = 12_000,
) -> tuple[np.ndarray, dict]:
    """
    Compute m4 for every 4-set on GPU.
    mmap-friendly host Y → single H2D → batched device products → D2H m4 only.
    """
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Y_host = np.ascontiguousarray(Y_np, dtype=np.float64)
    Y = cp.asarray(Y_host)
    n1 = int(quads.shape[0])
    m4_host = np.empty(n1, dtype=np.float64)
    t0 = time.perf_counter()
    for lo in range(0, n1, batch):
        hi = min(n1, lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4 = cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0)
        m4_host[lo:hi] = cp.asnumpy(m4)
    elapsed = time.perf_counter() - t0
    meta = {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": gpu_info(),
        "n_quads": n1,
        "batch": int(batch),
        "wall_s": float(elapsed),
        "N": int(Y_host.shape[0]),
        "n": int(Y_host.shape[1]),
        "io": {
            "maxplus_read": "mmap load_maxplus_array",
            "H2D": "single contiguous host→device copy of Y",
            "D2H": "m4 vector only (n_quads floats) — no full Y transfer back",
            "device_ops": "batched column products + mean; one CUDA context",
        },
    }
    return m4_host, meta


def build_lookup(
    quads: np.ndarray, m4: np.ndarray, kap: np.ndarray
) -> dict[tuple[int, int, int, int], tuple[float, int]]:
    """Map sorted 4-set → (m4, κ)."""
    out: dict[tuple[int, int, int, int], tuple[float, int]] = {}
    for i in range(quads.shape[0]):
        t = (int(quads[i, 0]), int(quads[i, 1]), int(quads[i, 2]), int(quads[i, 3]))
        out[t] = (float(m4[i]), int(kap[i]))
    return out


# ── ProcessPool structure walk over precomputed m4 ──────────────────────────
_G: dict = {}


def _init_walk(C: np.ndarray, lookup: dict, p: int) -> None:
    _G["C"] = C
    _G["lookup"] = lookup
    _G["p"] = p
    _blas1()


def _walk_shard(indices: list[tuple[int, int, int, int]]) -> tuple:
    """
    For each |κ|=1 4-set with same-sign residual r>0, and each a∈S with
    star_a=+1, accumulate S1 and S3. Return counts + maxes.
    """
    C = _G["C"]
    lookup = _G["lookup"]
    p = _G["p"]
    inv_p2 = 1.0 / (p * p)
    n = C.shape[0]
    n_pos = 0
    n_star_plus = 0
    n_S1_pos = 0
    max_S1 = -1e30
    max_S3 = -1e30
    max_r = 0.0
    max_rho = 0.0

    for S in indices:
        m4, k = lookup[S]
        if abs(k) != 1:
            continue
        rho = m4 - k * inv_p2
        r = rho * k
        if r <= 1e-15:
            continue
        n_pos += 1
        if r > max_r:
            max_r = r
        if abs(rho) > max_rho:
            max_rho = abs(rho)

        for v in S:
            if star_a(C, S, v) != 1:
                continue
            n_star_plus += 1
            others = tuple(x for x in S if x != v)
            S1 = 0.0
            S3 = 0.0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                m4p, kp = lookup[Sp]
                rp = m4p - kp * inv_p2
                w = float(C[v, rr]) * rp
                if abs(kp) == 1:
                    S1 += w
                else:
                    S3 += w
            if S1 > max_S1:
                max_S1 = S1
            if S3 > max_S3:
                max_S3 = S3
            if S1 > 1e-12:
                n_S1_pos += 1

    return n_pos, n_star_plus, n_S1_pos, max_S1, max_S3, max_r, max_rho


def cert_p(
    p: int,
    *,
    W: int,
    use_gpu: bool,
) -> dict:
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = load_maxplus_array(p, mmap=True)
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap = kappa_all(C, quads)

    t_gpu0 = time.perf_counter()
    if use_gpu:
        batch = 12_000 if p >= 7 else 50_000
        m4, gpu_meta = gpu_all_m4(Y, quads, batch=batch)
        gpu_used = True
    else:
        # multi-worker CPU fallback for m4 (still mmap Y)
        m4 = _cpu_all_m4(Y, quads, W=W)
        gpu_meta = {
            "backend": "cpu_processpool",
            "note": "GPU unavailable; m4 via ProcessPool + mmap",
        }
        gpu_used = False
    t_gpu = time.perf_counter() - t_gpu0

    lookup = build_lookup(quads, m4, kap)
    # only walk |κ|=1 candidates for sharding (still full for correctness)
    kappa1_sets = [
        (int(q[0]), int(q[1]), int(q[2]), int(q[3]))
        for q, k in zip(quads, kap)
        if abs(int(k)) == 1
    ]

    nsh = min(W, max(1, len(kappa1_sets) // 200))
    ss = (len(kappa1_sets) + nsh - 1) // nsh
    shards = [kappa1_sets[i : i + ss] for i in range(0, len(kappa1_sets), ss)]

    n_pos = n_sp = n_s1p = 0
    max_S1 = -1e30
    max_S3 = -1e30
    max_r = 0.0
    max_rho = 0.0
    t_walk0 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_walk,
        initargs=(C, lookup, p),
    ) as ex:
        futs = [ex.submit(_walk_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            a, b, c, d, e, f, g = fut.result()
            n_pos += a
            n_sp += b
            n_s1p += c
            max_S1 = max(max_S1, d)
            max_S3 = max(max_S3, e)
            max_r = max(max_r, f)
            max_rho = max(max_rho, g)
    t_walk = time.perf_counter() - t_walk0

    # budget checks
    rho_c = rho_budget_cand(p)
    # if S1≤0 and star=+1: pρ ≤ 2/p² + S3  ⇒  ρ ≤ 2/p³ + S3/p
    s3_implied_rho = (2.0 / (p * p * p)) + (max_S3 / p) if n_sp else None

    return {
        "p": int(p),
        "n": int(n),
        "N": int(np.asarray(Y).shape[0]),
        "n_quads": int(quads.shape[0]),
        "n_kappa1": int(len(kappa1_sets)),
        "n_pos_r_sets": int(n_pos),
        "n_star_plus_checks": int(n_sp),
        "n_S1_positive": int(n_s1p),
        "S1_always_nonpositive": n_s1p == 0,
        "max_S1_star_plus": float(max_S1) if n_sp else None,
        "max_S3_star_plus": float(max_S3) if n_sp else None,
        "max_same_sign_r": float(max_r),
        "max_abs_rho_on_pos": float(max_rho),
        "rho_budget_cand": float(rho_c),
        "M_cand": float(M_cand(p)),
        "L_abs": float(L_abs(p)),
        "s3_implied_rho_ub": s3_implied_rho,
        "s3_implied_le_cand": (
            s3_implied_rho is not None and s3_implied_rho <= rho_c + 1e-12
        ),
        "wall_m4_s": float(t_gpu),
        "wall_walk_s": float(t_walk),
        "wall_s": float(t_gpu + t_walk),
        "gpu_used": gpu_used,
        "gpu_m4": gpu_meta,
        "workers_walk": int(min(W, len(shards))),
        "io": (
            "mmap Max+; GPU all-m4 (single H2D, D2H m4 only); "
            "ProcessPool structure walk on precomputed lookup; "
            "write_json_atomic"
        ),
    }


def _cpu_all_m4(Y, quads: np.ndarray, *, W: int) -> np.ndarray:
    """Fallback: multi-worker m4 with mmap Y (no GPU)."""
    ypath = None
    # Y may be memmap — workers re-open via global path if possible
    # For safety, materialize contiguous if not path-backed
    Yc = np.ascontiguousarray(Y, dtype=np.float64)
    n1 = quads.shape[0]
    m4 = np.empty(n1, dtype=np.float64)
    # simple sequential mean in chunks (ProcessPool would need path; keep
    # vectorized numpy with BLAS=1 — still multi-core via numpy? No, OMP=1.
    # Use ProcessPool over shards with shared Yc via fork or init.
    global _G_m4

    def _init_m4(Yarr):
        global _G_m4
        _G_m4 = {"Y": Yarr}
        _blas1()

    def _shard_m4(lo_hi):
        lo, hi = lo_hi
        Yarr = _G_m4["Y"]
        sl = quads[lo:hi]
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        return lo, hi, np.mean(Yarr[:, a] * Yarr[:, b] * Yarr[:, c] * Yarr[:, d], axis=0)

    nsh = min(W, max(1, n1 // 2000))
    ss = (n1 + nsh - 1) // nsh
    shards = [(i, min(n1, i + ss)) for i in range(0, n1, ss)]
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_m4,
        initargs=(Yc,),
    ) as ex:
        for fut in as_completed([ex.submit(_shard_m4, sh) for sh in shards]):
            lo, hi, vals = fut.result()
            m4[lo:hi] = vals
    return m4


_G_m4: dict = {}


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"S1_sign: W={W} use_gpu={use_gpu} "
        f"io=mmap+GPU-m4+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )

    out: dict = {
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; GPU all-m4 (one CUDA context, H2D once, D2H m4 vector); "
            "ProcessPool S1/S3 walk on precomputed lookup; "
            "evidence via write_json_atomic (tmp+fsync+os.replace)"
        ),
        "results": {},
        "status": None,
    }

    if not use_gpu:
        print(
            "WARNING: GPU unavailable — falling back to ProcessPool m4 + mmap. "
            "Still atomic evidence; still multi-worker. Prefer free V100.",
            flush=True,
        )

    for p in (5, 7):
        print(f"  S1 cert p={p} (GPU m4 + multi-worker walk)...", flush=True)
        r = cert_p(p, W=W, use_gpu=use_gpu)
        out["results"][str(p)] = r
        print(
            f"    p={p}: pos_r={r['n_pos_r_sets']} star+={r['n_star_plus_checks']} "
            f"S1>0={r['n_S1_positive']} maxS1={r['max_S1_star_plus']:.8f} "
            f"maxS3={r['max_S3_star_plus']:.8f} "
            f"m4={r['wall_m4_s']:.3f}s walk={r['wall_walk_s']:.3f}s "
            f"gpu_m4={r['gpu_used']}",
            flush=True,
        )

    both = all(
        out["results"][str(p)]["S1_always_nonpositive"] for p in (5, 7)
    )
    out["S1_nonpos_both"] = both
    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["status"] = (
        f"S1≤0 at star_a=+1 on same-sign r>0 |κ|=1: p=5,7 => {both}. "
        f"GPU m4 p5={r5['wall_m4_s']:.2f}s p7={r7['wall_m4_s']:.2f}s; "
        f"walk p5={r5['wall_walk_s']:.2f}s p7={r7['wall_walk_s']:.2f}s; "
        f"mmap+atomic. If proved: pρ ≤ 2/p² + S3 (star=+1). "
        f"OPEN: prove S1≤0 general + bound S3. lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_S1_sign.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
