#!/usr/bin/env python3
"""Gate Step 5 at p=11 and record even-character Parseval at p=5,7.

For even α of F_q^* (α(-1)=1, α≠1, α²≠χ):
    λ(α) := 32 E|Z_α|² / [q(q-1)],  Z_α=∑_{a≠0} α(a) N(a)
should equal spec(Φ).  p=5,7 already exact on MuLab coordinates.  p=11 stored Max+ is
enumerated against paley_conference_prime_power (minmax_quadratic),
not 15590 field_ops — do not treat a mismatch as a kill of Step 5.

Also prints even-character Parseval:
    ∑_{even α} |Z_α|² = (q-1) ∑_{a≠0} N(a)²
and mean vs min vs QVAR threshold 3q(q-1)/16.

No flag flip.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import field_ops  # noqa: E402


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


def dlog_table(gen, q, fmul, one):
    tab = np.full(q, -1, dtype=np.int32)
    x = one
    for k in range(q - 1):
        tab[x] = k
        x = fmul(x, gen)
    return tab


def add_shift_index(q, fadd):
    """shift[a, x] = x+a in field encoding."""
    S = np.empty((q, q), dtype=np.int32)
    for a in range(q):
        for x in range(q):
            S[a, x] = fadd(x, a)
    return S


def moments_from_Z(Z, S, dlog, q, ks):
    """Z: (M,q) int8 ±1 on F_q. Returns E|ZN_k|² for each k in ks."""
    Dm = Z == -1
    M = Z.shape[0]
    # N[:, a] for a=1..q-1
    N = np.empty((M, q - 1), dtype=np.float64)
    for i, a in enumerate(range(1, q)):
        N[:, i] = (Dm & Dm[:, S[a]]).sum(axis=1)
    ang = np.empty((q - 1, len(ks)), dtype=np.complex128)
    for j, k in enumerate(ks):
        for i, a in enumerate(range(1, q)):
            ang[i, j] = np.exp(2j * np.pi * k * int(dlog[a]) / (q - 1))
    Sα = N @ ang  # (M, n_k)
    e2 = np.mean(np.abs(Sα) ** 2, axis=0)
    sumN2 = float(np.mean(np.sum(N * N, axis=1)))
    sumN = N.sum(axis=1)
    e_triv = float(np.mean(sumN.astype(np.float64) ** 2))
    return e2, sumN2, e_triv


# Fork-shared by parent after load (avoid 86× platter mmap).
_Y = None
_S = None
_ANG = None


def _p11_shard(lohi):
    import os

    os.environ["OMP_NUM_THREADS"] = "1"
    lo, hi = lohi
    Z = _Y[lo:hi, 1:]
    Dm = Z == -1
    q = _S.shape[0]
    N = np.empty((Z.shape[0], q - 1), dtype=np.float32)
    for i, a in enumerate(range(1, q)):
        N[:, i] = (Dm & Dm[:, _S[a]]).sum(axis=1)
    Sα = N @ _ANG
    return (np.abs(Sα) ** 2).sum(axis=0), Z.shape[0]


def p11_pool(Ypath, S, dlog, q, ks, n_workers=86):
    import multiprocessing as mp

    global _Y, _S, _ANG
    print("loading Max+ into RAM for fork-share...", flush=True)
    _Y = np.load(Ypath)
    M = int(_Y.shape[0])
    print(f"p=11 RAM |Max+|={M} n={_Y.shape[1]} workers={n_workers}", flush=True)
    _S = np.asarray(S, dtype=np.int32)
    _ANG = np.empty((q - 1, len(ks)), dtype=np.complex64)
    for j, k in enumerate(ks):
        for i, a in enumerate(range(1, q)):
            _ANG[i, j] = np.exp(2j * np.pi * k * int(dlog[a]) / (q - 1))
    step = (M + n_workers - 1) // n_workers
    shards = [(i * step, min(M, (i + 1) * step)) for i in range(n_workers)]
    shards = [s for s in shards if s[0] < s[1]]
    ctx = mp.get_context("fork")
    acc = np.zeros(len(ks), dtype=np.float64)
    nseen = 0
    with ctx.Pool(n_workers) as pool:
        for part, m in pool.imap_unordered(_p11_shard, shards, chunksize=1):
            acc += part
            nseen += m
            print(f"  rows {nseen}/{M}", flush=True)
    return acc / nseen


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    S = add_shift_index(q, fadd)
    half = (q - 1) // 2
    ks = [k for k in range(2, half, 2)]  # even in (0, half), includes ψ=(q-1)/4
    thr = 3 * q * (q - 1) / 16
    c = 32 / (q * (q - 1))
    print(f"p={p} q={q} n={n} #even-k={len(ks)} dimF={(q-5)//4} thr={thr:.6f}", flush=True)

    if p == 11:
        path = "/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy"
        e2 = p11_pool(path, S, dlog, q, ks)
        phi = np.load("/mnt/storage/e1work/maxplus_p11/phiZ_p11.npy")
        w = np.sort(np.linalg.eigvalsh(phi))
        from collections import Counter

        cl = sorted(Counter(np.round(w, 6)).items(), key=lambda kv: kv[0])
        print("Φ clusters:", cl, flush=True)
        lams = c * np.asarray(e2, dtype=np.float64)
        print(f"{'k':>6} {'E|Z|²':>14} {'32E/q(q-1)':>14}  nearest Φ", flush=True)
        uniq_w = np.array(sorted(set(np.round(w, 6))))
        for k, e, lam in zip(ks, e2, lams):
            j = int(np.argmin(np.abs(uniq_w - lam)))
            print(f"{k:6d} {float(e):14.6f} {float(lam):14.6f}  Φ={uniq_w[j]:.6f}  d={abs(uniq_w[j]-lam):.3e}", flush=True)
        print(f"min E|Z|²={float(np.min(e2)):.6f}  thr={thr:.6f}  gap={float(np.min(e2))-thr:.6f}", flush=True)
        print(f"min λ={float(np.min(lams)):.6f}  Φmin={w[0]:.6f}", flush=True)
        return

    from e1_gmin_m4_prop15590 import MuLab

    Y = MuLab(p, with_deg6=False).Yp.astype(np.int8)
    e2, sumN2, e_triv = moments_from_Z(Y[:, 1:], S, dlog, q, ks)
    lams = c * e2
    print(f"E∑_{{a≠0}} N(a)²={sumN2:.6f}  E|∑N|²={e_triv:.6f}", flush=True)
    parseval_rhs = (q - 1) * sumN2
    parseval_lhs = e_triv + float(np.sum(e2))  # missing quadratic k=half
    print(f"parseval even-without-χ: lhs~{parseval_lhs:.4f}  (q-1)∑N²={parseval_rhs:.4f}", flush=True)
    print(f"{'k':>6} {'ψ?':>3} {'E|Z|²':>12} {'λ':>12}", flush=True)
    psi = (q - 1) // 4
    for k, e, lam in zip(ks, e2, lams):
        print(f"{k:6d} {int(k==psi):3d} {e:12.4f} {lam:12.6f}", flush=True)
    print(f"min E|Z|²={float(np.min(e2)):.6f} thr={thr:.6f} gap={float(np.min(e2))-thr:.6f}", flush=True)
    print(f"mean leftover E|Z|²={float(np.mean(e2)):.6f}", flush=True)


if __name__ == "__main__":
    main()
