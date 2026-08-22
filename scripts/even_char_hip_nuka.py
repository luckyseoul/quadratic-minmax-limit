#!/usr/bin/env python3
"""Even-character N-moments on nuka RX 9070 XT (HIP / hiprtc via CuPy-ROCm).

V100 NVRTC cannot compile sm_70 against CUDA 13 fp8 headers.  This path
uses HIP RawKernel (hiprtc) on gfx1201.

Streams Max+ from local npy (mmap) in batches; one HIP context.
Does not load the 4.3 GiB host array.  No leftover flag flip.

Env: ~/.venvs/rocm72  (cupy 13.5.1, ROCm 7.2)
Data: /home/nick/e1work/maxplus_p11/maxplus_p11_eps1.npy on nuka.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

HIP_SRC = r"""
extern "C" __global__
void n_from_shifts(const char* Z, const int* shift, float* N, int B, int q) {
    int i = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    int a1 = (int)blockIdx.y;
    if (i >= B || a1 >= q - 1) return;
    const char* row = Z + ((size_t)i * (size_t)q);
    const int* sh = shift + ((size_t)(a1 + 1) * (size_t)q);
    int s = 0;
    for (int x = 0; x < q; ++x) {
        char zx = row[x];
        char za = row[sh[x]];
        s += (int)((zx == (char)-1) & (za == (char)-1));
    }
    N[(size_t)i * (size_t)(q - 1) + (size_t)a1] = (float)s;
}
"""


def field_ops_15590(p: int):
    """e = p*a + b ↔ a + b t, t^2 = nonresidue.  MuLab / 15590."""
    r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)

    def fmul(e1, e2):
        a1, b1 = divmod(e1, p)
        a2, b2 = divmod(e2, p)
        return p * ((a1 * a2 + r * b1 * b2) % p) + ((a1 * b2 + a2 * b1) % p)

    def fadd(e1, e2):
        a1, b1 = divmod(e1, p)
        a2, b2 = divmod(e2, p)
        return p * ((a1 + a2) % p) + ((b1 + b2) % p)

    one = p
    return fmul, fadd, one


def field_ops_minmax(p: int):
    """e = a + b p ↔ a + b ω, ω^2 = ia ω + ib.  paley_conference_prime_power.

    First irreducible (ia, ib) in the same nested loop as minmax_quadratic.
    Column 1+e of maxplus_p11_eps1.npy is this encoding.
    """
    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def fmul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        e0 = (c0 * d0 + c1 * d1 * ib) % p
        e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        return e0 + e1 * p

    def fadd(u, v):
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    one = 1
    return fmul, fadd, one, (ia, ib)


def primitive_root(q, fmul, one):
    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    return next(e for e in range(2, q) if order_of(e) == q - 1)


def main():
    import cupy as cp

    p = 11
    q = p * p
    n = q + 1
    ypath = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "/home/nick/e1work/maxplus_p11/maxplus_p11_eps1.npy"
    )
    phipath = ypath.with_name("phiZ_p11.npy")
    batch = int(sys.argv[2]) if len(sys.argv) > 2 else 262144
    # p=11 Max+ is paley_conference_prime_power labeling (minmax).
    field = sys.argv[3] if len(sys.argv) > 3 else "minmax"

    free, tot = cp.cuda.runtime.memGetInfo()
    print(
        f"HIP cupy={cp.__version__}  device_mem={free/1e9:.2f}/{tot/1e9:.2f} GiB",
        flush=True,
    )
    t_compile = time.perf_counter()
    ker = cp.RawKernel(HIP_SRC, "n_from_shifts")
    # force compile
    z_dummy = cp.zeros((1, q), dtype=cp.int8)
    sh_dummy = cp.zeros((q, q), dtype=cp.int32)
    n_dummy = cp.zeros((1, q - 1), dtype=cp.float32)
    ker(((1 + 255) // 256, q - 1), (256,), (z_dummy, sh_dummy, n_dummy, 1, q))
    cp.cuda.runtime.deviceSynchronize()
    print(f"hiprtc RawKernel compiled in {time.perf_counter()-t_compile:.2f}s", flush=True)

    if field == "minmax":
        fmul, fadd, one, irr = field_ops_minmax(p)
        print(f"field=minmax irr={irr} one={one}", flush=True)
    elif field == "15590":
        fmul, fadd, one = field_ops_15590(p)
        print(f"field=15590 one={one}", flush=True)
    else:
        raise SystemExit(f"unknown field {field}")
    gen = primitive_root(q, fmul, one)
    dlog = np.full(q, -1, dtype=np.int32)
    x = one
    for k in range(q - 1):
        dlog[x] = k
        x = fmul(x, gen)
    shift = np.empty((q, q), dtype=np.int32)
    for a in range(q):
        for xx in range(q):
            shift[a, xx] = fadd(xx, a)
    half = (q - 1) // 2
    ks = [k for k in range(2, half, 2)]
    ang = np.empty((q - 1, len(ks)), dtype=np.complex64)
    for j, k in enumerate(ks):
        for i, a in enumerate(range(1, q)):
            ang[i, j] = np.exp(2j * np.pi * k * int(dlog[a]) / (q - 1))
    sh_g = cp.asarray(shift)
    ang_g = cp.asarray(ang)
    print(f"p={p} q={q} #k={len(ks)} gen={gen} batch={batch}", flush=True)

    Y = np.load(ypath, mmap_mode="r")
    M, nn = int(Y.shape[0]), int(Y.shape[1])
    assert nn == n
    print(f"mmap {ypath}  |Max+|={M}", flush=True)

    acc = cp.zeros(len(ks), dtype=cp.float64)
    nseen = 0
    t0 = time.perf_counter()
    t_h2d = t_gpu = 0.0
    peak_used = 0
    for lo in range(0, M, batch):
        hi = min(M, lo + batch)
        sl = np.asarray(Y[lo:hi, 1:], dtype=np.int8)
        B = sl.shape[0]
        t1 = time.perf_counter()
        Zg = cp.asarray(sl)
        Ng = cp.empty((B, q - 1), dtype=cp.float32)
        t2 = time.perf_counter()
        ker(((B + 255) // 256, q - 1), (256,), (Zg, sh_g, Ng, B, q))
        Sα = Ng @ ang_g
        acc += (cp.abs(Sα).astype(cp.float64) ** 2).sum(axis=0)
        cp.cuda.runtime.deviceSynchronize()
        t3 = time.perf_counter()
        t_h2d += t2 - t1
        t_gpu += t3 - t2
        nseen += B
        used = tot - cp.cuda.runtime.memGetInfo()[0]
        peak_used = max(peak_used, used)
        if (lo // batch) % 8 == 0 or hi == M:
            print(
                f"  {nseen}/{M}  wall={time.perf_counter()-t0:.1f}s  "
                f"h2d={t_h2d:.1f}s gpu={t_gpu:.1f}s  used={used/1e6:.0f}MB",
                flush=True,
            )
        del Zg, Ng, Sα, sl
        cp.get_default_memory_pool().free_all_blocks()

    e2 = cp.asnumpy(acc) / nseen
    c = 32 / (q * (q - 1))
    lams = c * e2
    thr = 3 * q * (q - 1) / 16
    print(f"DONE nseen={nseen}  gpu={t_gpu:.2f}s h2d={t_h2d:.2f}s  peak_used={peak_used/1e6:.0f}MB", flush=True)
    print(f"min E|Z|^2={e2.min():.6f}  thr={thr:.6f}  min λ={lams.min():.6f}", flush=True)
    if phipath.exists():
        w = np.sort(np.linalg.eigvalsh(np.load(phipath)))
        uniq = np.array(sorted(set(np.round(w, 6))))
        print(f"Φ min={w[0]:.6f} max={w[-1]:.6f} n_clusters={len(uniq)}", flush=True)
        print(f"{'k':>6} {'E|Z|^2':>14} {'λ':>12}  nearest Φ", flush=True)
        n_close = 0
        for k, e, lam in zip(ks, e2, lams):
            j = int(np.argmin(np.abs(uniq - lam)))
            d = abs(uniq[j] - lam)
            n_close += int(d < 0.05)
            print(f"{k:6d} {e:14.4f} {lam:12.6f}  Φ={uniq[j]:.6f} d={d:.3e}", flush=True)
        print(f"k within 0.05 of a Φ cluster: {n_close}/{len(ks)}", flush=True)
    out = Path.home() / "e1work/maxplus_p11/even_char_hip_nuka.npz"
    np.savez(out, ks=np.array(ks), e2=e2, lams=lams)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
