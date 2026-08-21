"""CPU / SYCL outer for k=6 mesh (no cupy). V100/nuka keep gpu_inner."""
from __future__ import annotations

import os
import time

import numpy as np
from kgen5 import _prep_tables


class NumpyTester:
    """Host copy of GpuTester.test_batch / load_outer. No CUDA."""

    def __init__(self, p, k, Tm, UU, eps=1):
        self.p, self.k, self.eps = p, k, eps
        q = Tm.shape[1]
        self.q = q
        self.thi = (k - 1) * eps + p
        self.tlo = (k - 1) * eps - p
        self.Tm = Tm
        self.UU = np.ascontiguousarray(UU.astype(np.int32))
        L = np.zeros((k, q, p), dtype=np.float32)
        for j in range(k):
            for x in range(q):
                L[j, x, Tm[j, x]] = 1.0
        self.L = L
        LIDX = np.zeros((k, p, p), dtype=np.int64)
        for j in range(k):
            for s in range(p):
                LIDX[j, s] = np.where(Tm[j] == s)[0]
        self.LIDX = LIDX

    def load_outer(self, av, af, bases):
        p, k, q = self.p, self.k, self.q
        cont = np.zeros((k, p * p, q), dtype=np.int16)
        cov = np.zeros((k, p * p, q), dtype=np.int8)
        for j in range(k):
            for u in range(p):
                W = bases[j, u]
                sig = 2 * ((W[None, :] + np.arange(p)[:, None]) % p) - p + 2
                cont[j, u * p:(u + 1) * p, :] = sig[:, self.Tm[j]].astype(np.int16)
                ind = (((W[None, :] + np.arange(p)[:, None]) % p) == p - 1)
                cov[j, u * p:(u + 1) * p, :] = ind[:, self.Tm[j]].astype(np.int8)
        self.cont = cont
        self.cov = cov
        self.av = np.ascontiguousarray(av.astype(np.int32))
        self.af = np.ascontiguousarray(af.astype(np.int32))
        return cont, cov

    def test_batch(self, codes, fsums):
        p, k, q = self.p, self.k, self.q
        B = len(codes)
        cg = np.asarray(codes)
        fg = np.asarray(fsums, dtype=np.int32)
        ci = (cg // (16 ** 6)).astype(np.int32)
        PS = np.zeros((B, q), dtype=np.int16)
        for j in range(k):
            vidx = ((cg // (16 ** j)) % 16).astype(np.int32)
            uj = self.UU[ci, j]
            vj = self.av[j, uj, vidx]
            rows = uj * p + vj
            PS += self.cont[j][rows]
        isf0 = fg == 0
        okpt = (PS == self.thi) | (PS == self.tlo)
        cnt = okpt.sum(axis=1)
        f0_pass = isf0 & (cnt == q)
        diff = PS.astype(np.int32) - self.thi
        g = diff // (2 * p)
        COV = np.zeros((B, q), dtype=np.int16)
        for j in range(k):
            vidx = ((cg // (16 ** j)) % 16).astype(np.int32)
            uj = self.UU[ci, j]
            vj = self.av[j, uj, vidx]
            rows = uj * p + vj
            COV += self.cov[j][rows]
        gok = (g >= -1) & (g <= fg[:, None]) & (g <= COV)
        gcnt = gok.sum(axis=1)
        flip_ok = (~isf0) & (gcnt == q)
        line_ok = np.ones(B, dtype=bool)
        F = fg.astype(np.float32)
        gf = g.astype(np.float32)
        for j in range(k):
            Gl = gf @ self.L[j]
            vidx = ((cg // (16 ** j)) % 16).astype(np.int32)
            uj = self.UU[ci, j]
            fj = self.af[j, uj, vidx].astype(np.float32)
            Fj = F - fj
            lo = Gl - Fj[:, None]
            hi = Gl + p - Fj[:, None]
            can0 = (lo <= 0) & (hi >= 0)
            can1 = (lo <= p) & (hi >= p)
            okl = can0 | can1
            line_ok &= okl.sum(axis=1) == p
        flip_pass = flip_ok & line_ok
        return np.where(f0_pass)[0], np.where(flip_pass)[0]


def process_outer_host(p, k, q, upper, UU, Tm, c0, eps, tester, sols,
                       gen_cap=None, emit="cpu"):
    """cpu (numba prange) or sycl emit, then NumpyTester."""
    from gpu_inner_fast import gen_candidates_parallel

    if gen_cap is None:
        gen_cap = int(os.environ.get("GEN_CAP", "8000000"))
    host = os.environ.get("K6_HOST", "")
    if emit == "sycl":
        from gpu_gen_sycl import sycl_device_name
        print(f"  {host} SYCL device={sycl_device_name()!r} gen_cap={gen_cap} "
              f"NUMBA={os.environ.get('NUMBA_NUM_THREADS','')} "
              f"OMP={os.environ.get('OMP_NUM_THREADS','')}",
              flush=True)
    else:
        print(f"  {host} CPU emit gen_cap={gen_cap} "
              f"NUMBA={os.environ.get('NUMBA_NUM_THREADS','')} "
              f"OMP={os.environ.get('OMP_NUM_THREADS','')}",
              flush=True)
    thi = (k - 1) * eps + p
    tlo = (k - 1) * eps - p
    bases, av, af, an, aull = _prep_tables(p, k, upper, eps)
    cont, cov = tester.load_outer(av, af, bases)
    npr = min(10, q)
    probe_idx = np.linspace(0, q - 1, npr).astype(np.int64)
    probes = np.ascontiguousarray(cont[:, :, probe_idx])
    cprobes = np.ascontiguousarray(cov[:, :, probe_idx].astype(np.int16))
    codes = np.zeros(gen_cap, np.int64)
    fsums = np.zeros(gen_cap, np.int64)
    CH = 250_000
    s_ar = np.arange(p, dtype=np.int64)

    def decode(cc_sub):
        B = len(cc_sub)
        ci = (cc_sub // 16 ** 6).astype(np.int64)
        SIG = np.zeros((B, k, p), np.int64)
        FV = np.zeros((B, k), np.int64)
        for j in range(k):
            vidx = ((cc_sub // 16 ** j) % 16).astype(np.int64)
            u = UU[ci, j]
            v = av[j, u, vidx]
            FV[:, j] = af[j, u, vidx]
            SIG[:, j, :] = 2 * ((upper[j][None, :] + u[:, None] * s_ar[None, :]
                                 + v[:, None]) % p) - p + 2
        return SIG, FV

    def resolve(cc, ff):
        f0_idx, fl_idx = tester.test_batch(cc, ff)
        if len(f0_idx):
            SIG, FV = decode(cc[f0_idx])
            PS = np.zeros((len(f0_idx), q), np.int64)
            for j in range(k):
                PS += SIG[:, j, :][:, Tm[j]]
            for r in range(len(f0_idx)):
                sols.append(np.where(PS[r] == thi, 1, -1).astype(np.int8))
        if len(fl_idx):
            _resolve_flips(cc[fl_idx], ff[fl_idx])

    def _resolve_flips(cc_sub, ff_sub):
        from flipnb import flip_batch
        SIG, FV = decode(cc_sub)
        cap = max(200_000, 8 * len(cc_sub) + 1000)
        ybuf = np.zeros((cap, q), np.int8)
        ycount = np.zeros(1, np.int64)
        flip_batch(p, k, q, SIG, FV, Tm, tester.LIDX, thi, tlo,
                   ybuf, ycount, 5_000_000)
        ny = int(ycount[0])
        if ny > cap:
            if len(cc_sub) <= 1:
                raise RuntimeError(f"flip ybuf overflow ny={ny}")
            mid = len(cc_sub) // 2
            _resolve_flips(cc_sub[:mid], ff_sub[:mid])
            _resolve_flips(cc_sub[mid:], ff_sub[mid:])
            return
        for r in range(ny):
            sols.append(ybuf[r].copy())

    def emit_chunk(ulo, uhi):
        if emit == "sycl":
            from gpu_gen_sycl import gen_candidates_sycl
            return gen_candidates_sycl(
                p, k, av, af, an, aull, UU, c0, ulo, uhi,
                probes, cprobes, thi, tlo, codes, fsums,
            )
        return gen_candidates_parallel(
            p, k, av, af, an, aull, UU, c0, ulo, uhi,
            probes, cprobes, thi, tlo, codes, fsums,
        )

    worst_per_uu = 1
    for _ in range(k - 1):
        worst_per_uu *= p
    UCH = max(1, min(2000, gen_cap // max(1, worst_per_uu)))
    ncand = 0
    t_em = time.time()
    nchunk = 0
    for ulo in range(0, UU.shape[0], UCH):
        nc = emit_chunk(ulo, min(ulo + UCH, UU.shape[0]))
        ncand += nc
        nchunk += 1
        for lo in range(0, nc, CH):
            resolve(codes[lo:lo + CH], fsums[lo:lo + CH])
        if nchunk == 1 or nchunk % 20 == 0:
            print(f"  {host} emit chunks={nchunk} ncand={ncand} "
                  f"{time.time()-t_em:.0f}s",
                  flush=True)
    return ncand


def make_tester(p, k, Tm, UU, backend=None):
    backend = backend or os.environ.get("K6_BACKEND", "cuda")
    if backend in ("cpu", "sycl"):
        return NumpyTester(p, k, Tm, UU)
    from gpu_inner import GpuTester
    return GpuTester(p, k, Tm, UU)


def process_outer(p, k, q, upper, UU, Tm, c0, eps, tester, sols, backend=None):
    backend = backend or os.environ.get("K6_BACKEND", "cuda")
    if backend == "sycl":
        return process_outer_host(
            p, k, q, upper, UU, Tm, c0, eps, tester, sols, emit="sycl")
    if backend == "cpu":
        return process_outer_host(
            p, k, q, upper, UU, Tm, c0, eps, tester, sols, emit="cpu")
    from gpu_inner import process_outer_gpu
    return process_outer_gpu(p, k, q, upper, UU, Tm, c0, eps, tester, sols)
