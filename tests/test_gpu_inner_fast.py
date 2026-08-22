"""Serial vs two-pass vs atomic-r/w candidate gen emit the same packed codes.

CPU only — does not construct GpuTester / steal V100 from the live k=6 job.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

TOOL = Path(__file__).resolve().parents[1] / "scripts" / "maxplus_profile_enum"
sys.path.insert(0, str(TOOL))

os.environ.setdefault("NUMBA_NUM_THREADS", "8")


def test_cpu_atomics_available_and_unique():
    from gpu_inner_fast import (
        ATOMICS_AVAILABLE,
        atomic_fetch_add_i64,
        atomic_load_i64,
        atomic_store_i64,
    )
    from numba import njit, prange

    assert ATOMICS_AVAILABLE, "CPU LLVM atomic load/store/fetch-add failed to compile"

    @njit(parallel=True)
    def k(n, count, out):
        for i in prange(n):
            slot = atomic_fetch_add_i64(count, 0, 1)
            atomic_store_i64(out, slot, i + 1)
        return atomic_load_i64(count, 0)

    n = 20000
    count = np.zeros(1, np.int64)
    out = np.zeros(n, np.int64)
    tot = int(k(n, count, out))
    assert tot == n == int(count[0])
    assert int(out.min()) == 1 and int(out.max()) == n
    assert len(set(out.tolist())) == n


@pytest.mark.parametrize("p,k,nci", [(7, 4, 16), (13, 6, 4)])
@pytest.mark.parametrize("mode", ["atomic", "twopass"])
def test_parallel_gen_matches_serial(p, k, nci, mode):
    from bench_gen_fast import one_outer, run_par, run_serial

    d = one_outer(p=p, k=k, nci=nci)
    w = dict(d)
    w["u_hi"] = w["u_lo"] + 1
    run_serial(w)
    run_par(w, mode=mode)
    ns, _, cs, fs = run_serial(d)
    np_, _, cp, fp = run_par(d, mode=mode)
    assert ns == np_
    a = set(zip(cs[:ns].tolist(), fs[:ns].tolist()))
    b = set(zip(cp[:np_].tolist(), fp[:np_].tolist()))
    assert a == b
