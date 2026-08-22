#!/usr/bin/env python3
"""CPU serial vs A380 SYCL emit. Run on jellyfin. Does not touch the V100."""
from __future__ import annotations

import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "8")
sys.path.insert(0, os.path.dirname(__file__))

from bench_gen_fast import one_outer, run_serial
from gpu_gen_sycl import gen_candidates_sycl, sycl_device_name


def run_sycl(d):
    codes = __import__("numpy").zeros(d["cap"], __import__("numpy").int64)
    fsums = __import__("numpy").zeros(d["cap"], __import__("numpy").int64)
    t0 = time.time()
    n = gen_candidates_sycl(
        d["p"], d["k"], d["av"], d["af"], d["an"], d["aull"], d["UU"], d["c0"],
        d["u_lo"], d["u_hi"], d["probes"], d["cprobes"], d["thi"], d["tlo"],
        codes, fsums,
    )
    return n, time.time() - t0, codes, fsums


def main():
    nci = int(os.environ.get("NCI", "128"))
    print(f"SYCL device={sycl_device_name()!r} NCI={nci}", flush=True)
    d = one_outer(nci=nci)
    print(f"UU_rows={d['nUU']} chunk=[{d['u_lo']},{d['u_hi']})", flush=True)
    d_w = dict(d)
    d_w["u_hi"] = min(d["u_lo"] + 2, d["u_hi"])
    run_serial(d_w)
    run_sycl(d_w)
    ns, ts, cs, fs = run_serial(d)
    print(f"cpu_serial n={ns} t={ts:.3f}s rate={ns/max(ts,1e-9):.3e}/s", flush=True)
    ng, tg, cg, fg = run_sycl(d)
    a = set(zip(cs[:ns].tolist(), fs[:ns].tolist())) if ns else set()
    b = set(zip(cg[:ng].tolist(), fg[:ng].tolist())) if ng else set()
    ok = a == b
    print(
        f"sycl_a380  n={ng} t={tg:.3f}s rate={ng/max(tg,1e-9):.3e}/s "
        f"vs_serial={ts/max(tg,1e-9):.2f}x codes_equal={ok}",
        flush=True,
    )
    if not ok:
        print(f"|serial|={len(a)} |sycl|={len(b)} only_cpu={len(a-b)} only_sycl={len(b-a)}", flush=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
