#!/usr/bin/env python3
"""Quartic exceptional-mode variance by p=11 profile stratum.

Prop 15.589 identifies the only sub-n PSL constituent with a quartic
character psi satisfying psi^2=chi and

    lambda_exc = 32 E|Z_psi|^2 / (q(q-1)),
    Z_psi(D) = sum_{a,b in D, a!=b} psi(b-a).

This script computes the exact Gaussian-integer Z_psi distribution on the
persisted p=11 Max+ enumeration, split by profile activity k=1,3,4,5,6.
It is diagnostic evidence only; it does not promote a p=11 pattern to a
general theorem.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import cupy as cp
import numpy as np


P = 11
Q = P * P
COUNTS = {
    "k1": 2_772,
    "k3": 24_200,
    "k4": 58_080,
    "k5": 1_306_800,
    "k6": 36_065_260,
}


def primitive_root(q: int, mul) -> int:
    for g in range(2, q):
        x = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(x)
            x = mul(x, g)
        if len(seen) == q - 1 and x == 1:
            return g
    raise RuntimeError("no primitive root")


def quartic_kernel(mul) -> tuple[np.ndarray, np.ndarray]:
    g = primitive_root(Q, mul)
    psi = np.zeros(Q, dtype=np.complex64)
    x = 1
    for j in range(Q - 1):
        psi[x] = (1j) ** (j % 4)
        x = mul(x, g)

    neg = np.zeros(Q, dtype=np.int32)
    for a in range(Q):
        neg[a] = ((-a) % P) + ((-(a // P)) % P) * P
    add = np.empty((Q, Q), dtype=np.int32)
    for a in range(Q):
        for b in range(Q):
            add[a, b] = (a % P + b % P) % P + ((a // P + b // P) % P) * P
    diff = add[np.arange(Q)[:, None], neg[np.arange(Q)][None, :]]
    K = psi[diff]
    np.fill_diagonal(K, 0)
    # psi^2=chi audit.
    for a in range(1, Q):
        chi = 1 if psi[a] ** 2 == 1 else -1
        if abs(psi[a] ** 2 - chi) > 1e-6:
            raise RuntimeError("quartic character audit failed")
    return K.real.astype(np.float32), K.imag.astype(np.float32)


def main() -> None:
    root = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    yall = np.load(root / "maxplus_p11_eps1.npy", mmap_mode="r")
    if yall.shape != (sum(COUNTS.values()), Q + 1):
        raise RuntimeError(f"unexpected Max+ shape {yall.shape}")

    sys.path.insert(0, "/mnt/storage/e1work/scripts")
    from kgen import field_ctx  # type: ignore

    q, mul, _chi, _tr = field_ctx(P)
    if q != Q:
        raise RuntimeError("field mismatch")
    Kr, Ki = quartic_kernel(mul)
    Krg, Kig = cp.asarray(Kr), cp.asarray(Ki)
    chunk = int(os.environ.get("QUARTIC_CHUNK", "100000"))

    report = {"p": P, "q": Q, "counts": COUNTS, "strata": {}}
    offset = 0
    total_sum2 = 0
    started = time.time()
    for label, count in COUNTS.items():
        stop = offset + count
        hist: Counter[tuple[int, int]] = Counter()
        sum2 = 0
        nonzero = 0
        for lo in range(offset, stop, chunk):
            hi = min(stop, lo + chunk)
            # eps=+1 gauge; D is the -1 set and has p(p-1)/2 points.
            Dn = ((1 - np.asarray(yall[lo:hi, 1:], dtype=np.int8)) // 2).astype(
                np.float32
            )
            if not np.all(Dn.sum(axis=1) == P * (P - 1) // 2):
                raise RuntimeError("D-size audit failed")
            D = cp.asarray(Dn, dtype=cp.float32)
            # Keep GPU work to cuBLAS GEMM.  CuPy's NVRTC reduction kernels
            # are incompatible with this node's newer CUDA headers on sm_70.
            # The two ~48 MB batch results are reduced exactly on the CPU.
            DKr = cp.asnumpy(D @ Krg)
            DKi = cp.asnumpy(D @ Kig)
            ar = np.rint(np.sum(DKr * Dn, axis=1)).astype(np.int32)
            ai = np.rint(np.sum(DKi * Dn, axis=1)).astype(np.int32)
            vals, freq = np.unique(np.stack([ar, ai], axis=1), axis=0, return_counts=True)
            for v, c in zip(vals, freq):
                hist[(int(v[0]), int(v[1]))] += int(c)
            sq = ar.astype(np.int64) ** 2 + ai.astype(np.int64) ** 2
            sum2 += int(sq.sum())
            nonzero += int(np.count_nonzero(sq))
            del D, Dn, DKr, DKi
        variance = sum2 / count
        lam = 32.0 * variance / (Q * (Q - 1))
        report["strata"][label] = {
            "count": count,
            "sum_abs_Zpsi_sq": sum2,
            "E_abs_Zpsi_sq": variance,
            "lambda_exc_if_stratum_measure": lam,
            "nonzero_fraction": nonzero / count,
            "n_values": len(hist),
            "histogram": {f"{a}{b:+d}i": c for (a, b), c in sorted(hist.items())},
        }
        total_sum2 += sum2
        offset = stop
        print(label, report["strata"][label], flush=True)

    total_n = len(yall)
    total_var = total_sum2 / total_n
    total_lam = 32.0 * total_var / (Q * (Q - 1))
    report["total"] = {
        "count": total_n,
        "sum_abs_Zpsi_sq": total_sum2,
        "E_abs_Zpsi_sq": total_var,
        "floor_threshold": 3 * Q * (Q - 1) / 16,
        "lambda_exc": total_lam,
        "matches_verified_p11_spectrum": abs(total_lam - 8.664378396) < 2e-7,
        "elapsed_seconds": time.time() - started,
    }
    output = root / "stratum_quartic_variance_p11.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["total"], indent=2), flush=True)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
