#!/usr/bin/env python3
"""CPU-only bench: serial odometer vs gpu_inner_fast (two-pass / CPU atomics).

Never imports gpu_inner / cupy / numba.cuda — those open a V100 context.
"""
from __future__ import annotations

import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "16")

sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from kgen import square_coords
from kgen3 import prep_subset
from kgen5 import _prep_tables
from gpu_inner_fast import (
    ATOMICS_AVAILABLE,
    _gen_candidates_serial,
    gen_candidates_parallel,
)

if "cupy" in sys.modules or "numba.cuda" in sys.modules:
    raise SystemExit("refusing to run: cupy/numba.cuda imported (would hit the live V100)")


def one_outer(p=13, k=6, nci=64):
    dirs, forms, coords = square_coords(p)
    sub = tuple(range(k))
    ctx = prep_subset(p, list(sub), forms, coords)
    upper = np.zeros((k, p), dtype=np.int64)
    deg = k - 2
    if deg in ctx["kern"] and len(ctx["kern"][deg]):
        upper = (3 * np.outer(ctx["kern"][deg][0], (np.arange(p) ** deg) % p)) % p
    bases, av, af, an, aull = _prep_tables(p, k, upper, 1)
    UU = ctx["UU"]
    q = p * p
    npr = min(10, q)
    probe_idx = np.linspace(0, q - 1, npr).astype(np.int64)
    # probes from bases/Tm like gpu_inner.load_outer host slice
    Tm = ctx["Tm"]
    cont = np.zeros((k, p * p, q), dtype=np.int16)
    cov = np.zeros((k, p * p, q), dtype=np.int8)
    for j in range(k):
        for u in range(p):
            W = bases[j, u]
            sig = 2 * ((W[None, :] + np.arange(p)[:, None]) % p) - p + 2
            cont[j, u * p : (u + 1) * p, :] = sig[:, Tm[j]].astype(np.int16)
            ind = ((W[None, :] + np.arange(p)[:, None]) % p) == (p - 1)
            cov[j, u * p : (u + 1) * p, :] = ind[:, Tm[j]].astype(np.int8)
    probes = np.ascontiguousarray(cont[:, :, probe_idx])
    cprobes = np.ascontiguousarray(cov[:, :, probe_idx].astype(np.int16))
    thi = (k - 1) + p
    tlo = (k - 1) - p
    u_lo, u_hi = 0, min(nci, UU.shape[0])
    cap = 40_000_000
    return dict(
        p=p, k=k, av=av, af=af, an=an, aull=aull, UU=UU, c0=ctx["c0"],
        u_lo=u_lo, u_hi=u_hi, probes=probes, cprobes=cprobes, thi=thi, tlo=tlo,
        cap=cap, nUU=int(UU.shape[0]),
    )


def run_serial(d):
    codes = np.zeros(d["cap"], np.int64)
    fsums = np.zeros(d["cap"], np.int64)
    count = np.zeros(1, np.int64)
    t0 = time.time()
    _gen_candidates_serial(
        d["p"], d["k"], d["av"], d["af"], d["an"], d["aull"], d["UU"], d["c0"],
        d["u_lo"], d["u_hi"], codes, fsums, count, d["probes"], d["cprobes"],
        d["thi"], d["tlo"],
    )
    return int(count[0]), time.time() - t0, codes, fsums


def run_par(d, mode="atomic"):
    codes = np.zeros(d["cap"], np.int64)
    fsums = np.zeros(d["cap"], np.int64)
    t0 = time.time()
    n = gen_candidates_parallel(
        d["p"], d["k"], d["av"], d["af"], d["an"], d["aull"], d["UU"], d["c0"],
        d["u_lo"], d["u_hi"], d["probes"], d["cprobes"], d["thi"], d["tlo"],
        codes, fsums, mode=mode,
    )
    return n, time.time() - t0, codes, fsums


def main():
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    nci = int(os.environ.get("NCI", "32"))
    threads = os.environ.get("NUMBA_NUM_THREADS", "16")
    os.environ["NUMBA_NUM_THREADS"] = threads
    print(
        f"CPU-ONLY NUMBA_NUM_THREADS={threads} NCI={nci} "
        f"ATOMICS_AVAILABLE={ATOMICS_AVAILABLE} cupy={('cupy' in sys.modules)}",
        flush=True,
    )
    d = one_outer(nci=nci)
    print(f"UU_rows={d['nUU']} chunk=[{d['u_lo']},{d['u_hi']})", flush=True)
    d_w = dict(d)
    d_w["u_hi"] = min(d["u_lo"] + 2, d["u_hi"])
    run_serial(d_w)
    run_par(d_w, mode="twopass")
    if ATOMICS_AVAILABLE:
        run_par(d_w, mode="atomic")
    ns, ts, cs, fs = run_serial(d)
    print(f"serial   n={ns} t={ts:.3f}s rate={ns/max(ts,1e-9):.3e}/s", flush=True)
    a = set(zip(cs[:ns].tolist(), fs[:ns].tolist())) if ns else set()
    for mode in ("twopass", "atomic"):
        if mode == "atomic" and not ATOMICS_AVAILABLE:
            print("atomic SKIP (not available)", flush=True)
            continue
        np_, tp, cp, fp = run_par(d, mode=mode)
        b = set(zip(cp[:np_].tolist(), fp[:np_].tolist())) if np_ else set()
        ok = a == b
        print(
            f"{mode:8s} n={np_} t={tp:.3f}s rate={np_/max(tp,1e-9):.3e}/s "
            f"speedup={ts/max(tp,1e-9):.2f}x codes_equal={ok}",
            flush=True,
        )
        if not ok:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
