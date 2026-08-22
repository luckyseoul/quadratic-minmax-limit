#!/usr/bin/env python3
"""CPU serial/twopass/atomic vs device atomic emit.

Run on nuka (RX 9070 XT / HIP). Do not run on soulkiller while k=6 owns the V100.
"""
from __future__ import annotations

import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "16")

sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from bench_gen_fast import one_outer, run_serial, run_par
from gpu_inner_fast import ATOMICS_AVAILABLE
from gpu_gen_device import gen_candidates_device


def run_gpu(d):
    codes = np.zeros(d["cap"], np.int64)
    fsums = np.zeros(d["cap"], np.int64)
    # warmup/compile not included: caller does a tiny pass first
    t0 = time.time()
    n = gen_candidates_device(
        d["p"], d["k"], d["av"], d["af"], d["an"], d["aull"], d["UU"], d["c0"],
        d["u_lo"], d["u_hi"], d["probes"], d["cprobes"], d["thi"], d["tlo"],
        codes, fsums,
    )
    return n, time.time() - t0, codes, fsums


def main():
    import cupy as cp
    nci = int(os.environ.get("NCI", "128"))
    threads = os.environ.get("NUMBA_NUM_THREADS", "16")
    free, tot = cp.cuda.runtime.memGetInfo()
    print(
        f"device mem_free={free/2**30:.2f}GiB tot={tot/2**30:.2f}GiB "
        f"cupy={cp.__version__} NUMBA_NUM_THREADS={threads} NCI={nci} "
        f"ATOMICS_AVAILABLE={ATOMICS_AVAILABLE}",
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
    run_gpu(d_w)
    ns, ts, cs, fs = run_serial(d)
    print(f"cpu_serial  n={ns} t={ts:.3f}s rate={ns/max(ts,1e-9):.3e}/s", flush=True)
    a = set(zip(cs[:ns].tolist(), fs[:ns].tolist())) if ns else set()
    for mode in ("twopass", "atomic"):
        if mode == "atomic" and not ATOMICS_AVAILABLE:
            print("cpu_atomic SKIP", flush=True)
            continue
        np_, tp, cp_, fp = run_par(d, mode=mode)
        b = set(zip(cp_[:np_].tolist(), fp[:np_].tolist())) if np_ else set()
        ok = a == b
        print(
            f"cpu_{mode:8s} n={np_} t={tp:.3f}s rate={np_/max(tp,1e-9):.3e}/s "
            f"speedup={ts/max(tp,1e-9):.2f}x codes_equal={ok}",
            flush=True,
        )
        if not ok:
            raise SystemExit(1)
    ng, tg, cg, fg = run_gpu(d)
    gset = set(zip(cg[:ng].tolist(), fg[:ng].tolist())) if ng else set()
    ok = a == gset
    print(
        f"gpu_atomic  n={ng} t={tg:.3f}s rate={ng/max(tg,1e-9):.3e}/s "
        f"speedup={ts/max(tg,1e-9):.2f}x vs_cpu80t_note=nuka_local "
        f"codes_equal={ok}",
        flush=True,
    )
    if not ok:
        print(f"|serial|={len(a)} |gpu|={len(gset)} only_cpu={len(a-gset)} only_gpu={len(gset-a)}", flush=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
