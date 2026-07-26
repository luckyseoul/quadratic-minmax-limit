#!/usr/bin/env python3
"""Exact m_9 by sharding free-edge sign space across many processes."""
import os
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from numba import njit

N = 9
# free edges among vertices 1..8
FREE = [(i, j) for i in range(1, N) for j in range(i + 1, N)]
NF = len(FREE)  # 28
NPAT = 1 << (N - 1)  # 256


@njit(cache=True)
def prep_tables():
    n = 9
    nf = 28
    npat = 256
    free_i = np.empty(nf, np.int32)
    free_j = np.empty(nf, np.int32)
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            free_i[k] = i
            free_j[k] = j
            k += 1
    XX = np.empty((npat, nf), np.int32)
    Q0 = np.empty(npat, np.int32)
    for p in range(npat):
        s = 0
        for j in range(1, n):
            bit = (p >> (j - 1)) & 1
            xj = 1 if bit else -1
            s += xj
        Q0[p] = s
        for e in range(nf):
            i = free_i[e]
            j = free_j[e]
            xi = 1 if ((p >> (i - 1)) & 1) else -1
            xj = 1 if ((p >> (j - 1)) & 1) else -1
            XX[p, e] = xi * xj
    return XX, Q0


@njit(cache=True)
def shard_min_max(XX, Q0, start: int, count: int, nf: int) -> int:
    """Enumerate signs as integers in [start, start+count), return min over matrices of max|Q|."""
    npat = Q0.shape[0]
    best = 10**9
    signs = np.empty(nf, np.int32)
    Q = np.empty(npat, np.int32)
    end = start + count
    for bits in range(start, end):
        # decode signs
        for e in range(nf):
            signs[e] = 1 if ((bits >> e) & 1) == 0 else -1
        # Q = Q0 + XX @ signs
        for p in range(npat):
            s = Q0[p]
            for e in range(nf):
                s += XX[p, e] * signs[e]
            Q[p] = s
        mx = 0
        for p in range(npat):
            v = Q[p]
            if v < 0:
                v = -v
            if v > mx:
                mx = v
        if mx < best:
            best = mx
    return best


def worker(args):
    start, count = args
    XX, Q0 = prep_tables()
    return shard_min_max(XX, Q0, start, count, 28)


def main():
    ncpu = os.cpu_count() or 88
    workers = min(80, ncpu)
    total = 1 << 28
    chunk = (total + workers - 1) // workers
    shards = []
    for w in range(workers):
        start = w * chunk
        if start >= total:
            break
        count = min(chunk, total - start)
        shards.append((start, count))
    print(f"m_9: {workers} workers, {total} matrices, ~{chunk} each", flush=True)
    # warm numba in main
    XX, Q0 = prep_tables()
    _ = shard_min_max(XX, Q0, 0, 1000, 28)
    t0 = time.time()
    best = 10**9
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(worker, s) for s in shards]
        done = 0
        for fut in as_completed(futs):
            b = fut.result()
            best = min(best, b)
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"  progress {done}/{len(futs)} global_best={best} t={time.time()-t0:.1f}s", flush=True)
    print(f"RESULT m_9={best} alpha={best/9**1.5:.6f} time={time.time()-t0:.1f}s", flush=True)
    return best


if __name__ == "__main__":
    main()
