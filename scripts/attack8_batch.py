#!/usr/bin/env python3
"""Batch attack: harvest many distinct order-8 minimizers, run the active-set
completion optimizer on each, hunt for a completion <= 28.

m_16 <= 28 would settle the first open doubling instance (2*sqrt(2)*m_8 = 28.28).
Self-test first: the C14 completion must evaluate to 21.
"""
import json
import time

import numpy as np


def canon(X):
    k = X.shape[0]
    g = np.ones(k, dtype=np.int64)
    g[1:] = X[0, 1:]
    return X * (g[:, None] * g[None, :])


def main():
    import cupy as cp
    n = 8
    N2 = 1 << (2 * n - 1)

    # ---- self-test: C14 completion evaluates to 21
    j = json.load(open("/home/nick/quadratic-minmax-limit/evidence/"
                      "paley_two_block_c14_witness.json"))
    S2 = np.array(j["S2"], dtype=np.float32)
    K = cp.asarray(S2)
    uu = cp.arange(1 << 13, dtype=cp.int64)
    S = cp.ones((1 << 13, 14), dtype=cp.float32)
    for c in range(1, 14):
        S[:, c] = 1 - 2 * ((uu >> (c - 1)) & 1)
    vals = 0.5 * cp.sum((S @ K) * S, axis=1)
    selftest = float(cp.abs(vals).max())
    print(f"[selftest] C14 Phi = {selftest} (expect 21)", flush=True)
    assert abs(selftest - 21) < 1e-6

    # ---- harvest distinct order-8 minimizers by multi-start descent
    ux = cp.arange(1 << (n - 1), dtype=cp.int64)
    X8 = cp.ones((1 << (n - 1), n), dtype=cp.int8)
    for c in range(1, n):
        X8[:, c] = (1 - 2 * ((ux >> (c - 1)) & 1)).astype(cp.int8)
    pairs = [(i, k) for i in range(n) for k in range(i + 1, n)]
    P = len(pairs)
    G8 = cp.empty((P, 1 << (n - 1)), dtype=cp.float32)
    for a, (i, k) in enumerate(pairs):
        G8[a] = X8[:, i] * X8[:, k]
    del X8

    rng = np.random.default_rng(8)
    sources = {}
    t0 = time.time()
    starts = 0
    while len(sources) < 2000 and time.time() - t0 < 300:
        starts += 1
        Kv = cp.asarray(rng.choice(np.array([-1.0, 1.0], dtype=np.float32), P))
        Q = (G8 * Kv[:, None]).sum(axis=0)
        val = float(cp.abs(Q).max())
        # 1-flip descent
        while True:
            M = cp.abs(Q[None, :] - 2.0 * Kv[:, None] * G8)
            mv = M.max(axis=1)
            a = int(mv.argmin())
            if float(mv[a]) < val - 1e-9:
                Q = Q - 2.0 * Kv[a] * G8[a]
                Kv[a] = -Kv[a]
                val = float(mv[a])
            else:
                break
        if val == 10:
            Kmat = np.zeros((n, n), dtype=np.int64)
            Kvn = cp.asnumpy(Kv)
            for a, (i, k) in enumerate(pairs):
                Kmat[i, k] = Kmat[k, i] = int(Kvn[a])
            key = canon(Kmat).tobytes()
            if key not in sources:
                sources[key] = Kmat
    print(f"[harvest] {len(sources)} distinct order-8 minimizers from "
          f"{starts} starts [{time.time()-t0:.0f}s]", flush=True)

    # ---- batches over the completion family
    uxn = np.arange(1 << (n - 1), dtype=np.int64)
    Xn = np.ones((1 << (n - 1), n), dtype=np.float32)
    for c in range(1, n):
        Xn[:, c] = 1 - 2 * ((uxn >> (c - 1)) & 1)
    # y-states FULL (2^n)
    uyn = np.arange(1 << n, dtype=np.int64)
    Yn = np.ones((1 << n, n), dtype=np.float32)
    for c in range(1, n):
        Yn[:, c] = 1 - 2 * ((uyn >> (c - 1)) & 1)

    def active_set_attack(Af, iters=60, topk=400, rng_local=None):
        Qx = 0.5 * np.sum((Xn @ Af) * Xn, axis=1)
        Qy = 0.5 * np.sum((Yn @ Af) * Yn, axis=1)
        Dn = Qx[:, None] - Qy[None, :]
        B = (rng_local or np.random.default_rng(7)).choice(
            np.array([-1.0, 1.0], dtype=np.float32), (n, n))
        best = None
        for it in range(iters):
            if it > 0 and it % 40 == 0:
                B = (rng_local or np.random.default_rng(it)).choice(
                    np.array([-1.0, 1.0], dtype=np.float32), (n, n))
            V = Dn + Xn @ (B @ Yn.T)
            phi = float(np.abs(V).max())
            if best is None or phi < best:
                best = phi
            flat = np.abs(V).ravel()
            k = min(topk, flat.size)
            idx = np.argpartition(flat, flat.size - k)[flat.size - k:]
            i = idx // V.shape[1]
            j = idx % V.shape[1]
            Dk = Dn[i, j].astype(np.float64)
            Z = (Xn[i][:, :, None] * Yn[j][:, None, :]).reshape(len(i), -1)
            Bf = B.astype(np.float64)
            Bflat = Bf.reshape(-1)
            cross = Bflat @ Z.T
            improved = True
            sweeps = 0
            while improved and sweeps < 12:
                improved = False
                sweeps += 1
                for ab in range(n * n):
                    alt = cross - 2.0 * Bflat[ab] * Z[:, ab]
                    cur = np.abs(Dk + cross)
                    new = np.abs(Dk + alt)
                    if new.max() < cur.max() - 1e-9:
                        Bflat[ab] *= -1.0
                        cross = alt
                        improved = True
            B = Bf.astype(np.float32)
        return best, B

    t1 = time.time()
    results = []
    for idx, (key, Kmat) in enumerate(sources.items()):
        best, B = active_set_attack(Kmat.astype(np.float32),
                                    rng_local=np.random.default_rng(1000 + idx))
        results.append((best, idx, Kmat, B))
        if best <= 28:
            print(f"[JACKPOT] source #{idx}: completion Phi = {best}!!",
                  flush=True)
        if idx % 200 == 0:
            print(f"[attack] {idx}/{len(sources)} sources, best so far "
                  f"{min(r[0] for r in results)}, t={time.time()-t1:.0f}s",
                  flush=True)
    results.sort(key=lambda r: r[0])
    from collections import Counter
    dist = Counter(int(r[0]) for r in results)
    print(f"[attack] DONE over {len(results)} sources in "
          f"{time.time()-t1:.0f}s; distribution {dict(sorted(dist.items()))}",
          flush=True)
    print(f"[attack] best completion: {results[0][0]}", flush=True)
    # save all <= 30 for deep pushes
    winners = [{"phi": float(b), "A": K.tolist(), "B": B.tolist()}
               for b, i, K, B in results if b <= 30]
    json.dump(winners, open("/home/nick/scratch/completion_hits_8.json", "w"))
    print(f"[attack] saved {len(winners)} sources with completion <= 30",
          flush=True)


if __name__ == "__main__":
    main()
