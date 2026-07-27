#!/usr/bin/env python3
"""Parallel SA for Max-covering perfect matchings at p=5; exact maxR + MITM Phi."""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import (  # noqa: E402
    maxR_matching_level,
    paley_conference_prime_power,
    phi_mitm,
)


def worker(seed: int) -> dict | None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    m = n // 2
    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.exists():
        # rebuild lightly via halfspace orbit not available — fail soft
        return None
    Max = np.load(cache).astype(np.float64)
    Max_full = np.vstack([Max, -Max])
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    pairs = [(int(perm[2 * i]), int(perm[2 * i + 1])) for i in range(m)]

    def minS() -> float:
        s = np.zeros(len(Max_full))
        for i, j in pairs:
            s += C[i, j] * Max_full[:, i] * Max_full[:, j]
        return float(s.min())

    cur = minS()
    best = cur
    bp = list(pairs)
    for t in range(40000):
        a, b = int(rng.integers(0, m)), int(rng.integers(0, m))
        if a == b:
            continue
        i, j = pairs[a]
        k, l = pairs[b]
        if rng.integers(0, 2) == 0:
            na, nb = (i, k), (j, l)
        else:
            na, nb = (i, l), (j, k)
        if len({na[0], na[1], nb[0], nb[1]}) < 4:
            continue
        old = pairs[a], pairs[b]
        pairs[a], pairs[b] = na, nb
        ns = minS()
        T = 2.0 * (1 - t / 40000) + 1e-6
        if ns >= cur or rng.random() < np.exp((ns - cur) / T):
            cur = ns
            if cur > best:
                best = cur
                bp = list(pairs)
                if best >= 1.0 - 1e-9:
                    break
        else:
            pairs[a], pairs[b] = old
    if best < 1.0 - 1e-9:
        return None
    mr = maxR_matching_level(C, bp, p)
    A = C.copy()
    for i, j in bp:
        A[i, j] *= -1
        A[j, i] *= -1
    ph = float(phi_mitm(A))
    return {
        "seed": seed,
        "minS": best,
        "maxR": mr,
        "phi": ph,
        "Phi": 65.0,
        "undercut": ph < 64.5,
        "matching": [[int(a), int(b)] for a, b in bp],
        "criterion": mr >= 60,
    }


def main():
    t0 = time.time()
    covers = []
    seeds = list(range(2000, 2400))
    with ProcessPoolExecutor(max_workers=40) as ex:
        futs = [ex.submit(worker, s) for s in seeds]
        for f in as_completed(futs):
            r = f.result()
            if r is not None:
                covers.append(r)
                print(
                    f"COVER seed={r['seed']} maxR={r['maxR']} phi={r['phi']} crit={r['criterion']}",
                    flush=True,
                )
    print(f"found {len(covers)} covers in {time.time() - t0:.1f}s", flush=True)
    out = {
        "n_seeds": len(seeds),
        "n_covers": len(covers),
        "covers": covers,
        "all_criterion": all(c["criterion"] for c in covers) if covers else None,
        "all_no_undercut": all(not c["undercut"] for c in covers) if covers else None,
        "maxR_values": sorted({c["maxR"] for c in covers}),
        "status": "OPEN — batch census of SA Max-covers; not forall-M proof",
    }
    path = ROOT / "evidence" / "e1_gamma_cover_batch.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path, flush=True)
    print(json.dumps({k: out[k] for k in out if k != "covers"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
