#!/usr/bin/env python3
"""Explicit Paley conference evaluation: Phi(C_{q+1}) by exact chunked enumeration.

Validated against C14 (expect 21) and C18 (expect 33):
    python scripts/paley_conference_eval.py --q 13
    python scripts/paley_conference_eval.py --q 37
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paley_structural_lift_family import paley_conference  # noqa: E402
from paley_two_block_decomposition import char_split  # noqa: E402


def phi_chunked(K, chunk_bits=21):
    import cupy as cp
    n2 = K.shape[0]
    m = 1 << (n2 - 1)
    step = 1 << chunk_bits
    Kt = cp.asarray(K, dtype=cp.float32)
    best = 0.0
    for lo in range(0, m, step):
        hi = min(lo + step, m)
        u = cp.arange(lo, hi, dtype=cp.int64)
        S = cp.ones((hi - lo, n2), dtype=cp.int8)
        for c in range(1, n2):
            S[:, c] = (1 - 2 * ((u >> (c - 1)) & 1)).astype(cp.int8)
        Sf = S.astype(cp.float32)
        q = cp.sum(cp.matmul(Sf, Kt) * Sf, axis=1) * cp.float32(0.5)
        v = float(cp.max(cp.abs(q)))
        if v > best:
            best = v
        del S, Sf, q
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--chunk-bits", type=int, default=21)
    args = ap.parse_args()

    t0 = time.time()
    C = paley_conference(args.q)
    n = C.shape[0]
    gram = C.astype(np.float64) @ C.astype(np.float64).T
    diag_err = float(np.abs(np.diag(gram) - args.q).max())
    off_err = float(np.abs(gram - np.diag(np.diag(gram))).max())
    print(f"q={args.q} order={n} Gram check: diag_err={diag_err} "
          f"offdiag_err={off_err}", flush=True)

    T, Tc, nu = char_split(args.q)
    S = np.ascontiguousarray(C[np.ix_(T + Tc, T + Tc)])
    m = n // 2
    block_ok = (bool((S[:m, :m] == -S[m:, m:]).all())
                and bool((S[:m, m:] == S[m:, :m].T).all()))
    print("two-block structure:", block_ok, flush=True)

    bphi = phi_chunked(np.ascontiguousarray(S[:m, :m]), args.chunk_bits)
    phi = phi_chunked(S, args.chunk_bits)
    print(f"Phi(block A{m}) = {bphi}", flush=True)
    print(f"Phi(C{n}) = {phi}  [{time.time()-t0:.0f}s]", flush=True)

    if args.out:
        json.dump({"q": args.q, "order": n, "phi": phi, "block_phi": bphi,
                   "gram_diag_err": diag_err, "gram_offdiag_err": off_err,
                   "two_block_ok": block_ok,
                   "seconds": round(time.time() - t0, 1)},
                  open(args.out, "w"), indent=1)
        print("->", args.out, flush=True)


if __name__ == "__main__":
    main()
