#!/usr/bin/env python3
"""Benchmark the exact Phi kernel (the scout-v5 scoring core) on any engine.

Kernel: Phi(K) = max over 2^(2n-1) states S of |sum((S @ K) * S, axis=1)/2|,
chunked over states; identical protocol on CUDA (cupy), ROCm (cupy), and numpy.

Usage:  python bench_phi.py --matrix K16_phi30.json --expected 30 --reps 3
"""
import argparse
import json
import time


def build_doubled(A, seed):
    """K = [[A, Z],[Z^T, -A]] with a deterministic float Z (the scout object)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    Z = rng.standard_normal((n, n))
    K = np.zeros((2 * n, 2 * n), dtype=np.float64)
    K[:n, :n] = A
    K[:n, n:] = Z
    K[n:, :n] = Z.T
    K[n:, n:] = -A
    return K


def build_chunk(u0, u1, n2, xp, int8, int64):
    u = xp.arange(u0, u1, dtype=int64)
    S = xp.ones((u1 - u0, n2), dtype=int8)
    for c in range(1, n2):
        S[:, c] = (1 - 2 * ((u >> (c - 1)) & 1)).astype(int8)
    return S


def phi_cupy(K, chunk_bits, cp):
    n2 = K.shape[0]
    m = 1 << (n2 - 1)
    step = 1 << chunk_bits
    Kt = cp.asarray(K, dtype=cp.float32)
    best = 0.0
    for lo in range(0, m, step):
        hi = min(lo + step, m)
        S = build_chunk(lo, hi, n2, cp, cp.int8, cp.int64)
        Sf = S.astype(cp.float32)
        q = cp.sum(cp.matmul(Sf, Kt) * Sf, axis=1) * cp.float32(0.5)
        v = float(cp.max(cp.abs(q)))
        if v > best:
            best = v
        del S, Sf, q
    return best


def phi_numpy(K, chunk_bits, np):
    n2 = K.shape[0]
    m = 1 << (n2 - 1)
    step = 1 << chunk_bits
    Kf = K.astype(np.float32)
    best = 0.0
    for lo in range(0, m, step):
        hi = min(lo + step, m)
        S = build_chunk(lo, hi, n2, np, np.int8, np.int64)
        Sf = S.astype(np.float32)
        q = np.einsum('mi,mi->m', Sf, Sf @ Kf) * 0.5
        v = float(np.abs(q).max())
        if v > best:
            best = v
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--expected", type=float, default=None)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--chunk-bits", type=int, default=21)
    ap.add_argument("--engine", choices=["auto", "cupy", "numpy"], default="auto")
    ap.add_argument("--doubled", action="store_true")
    ap.add_argument("--seed", type=int, default=999)
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    import numpy as np
    K = np.array(json.load(open(args.matrix))["A"], dtype=np.int64)
    if args.doubled:
        K = build_doubled(K, args.seed)
    n2 = K.shape[0]
    states = 1 << (n2 - 1)

    engine = args.engine
    cp = None
    dev = "cpu"
    ver = "numpy " + np.__version__
    if engine in ("auto", "cupy"):
        try:
            import cupy as cp  # noqa
            dev = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
            ver = "cupy " + cp.__version__
            engine = "cupy"
        except Exception as e:
            if engine == "cupy":
                raise
            engine = "numpy"
            cp = None

    times = []
    val = None
    for r in range(args.reps):
        t0 = time.time()
        if engine == "cupy":
            val = phi_cupy(K, args.chunk_bits, cp)
        else:
            val = phi_numpy(K, args.chunk_bits, np)
        times.append(time.time() - t0)
    ok = (args.expected is None) or (abs(val - args.expected) < args.tol)
    out = dict(n=n2, states=states, engine=engine, runtime=ver, device=dev,
               times=[round(t, 3) for t in times],
               best=round(min(times), 3), phi=val, expected=args.expected,
               tol=args.tol, ok=ok)
    print(json.dumps(out))
    if not ok:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
