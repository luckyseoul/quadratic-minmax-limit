"""Device candidate emit: one thread per (UU-row, odometer combo).

Uses CuPy RawKernel + 64-bit atomicAdd (cursor) and atomicExch (payload).
HIP (nuka RX 9070 XT) and CUDA (V100 / H100) both lower these atomics.
CuPy RawKernel is compiled for the current device — no baked sm_70 cubin.
"""
from __future__ import annotations

import numpy as np

_KERNEL = r"""
extern "C" __global__
void emit_cands(
    const int p,
    const int k,
    const int u_lo,
    const int nci,
    const int c0,
    const int thi,
    const int tlo,
    const int npr,
    const int two_p,
    const unsigned long long ncombo,
    const int* __restrict__ av,
    const int* __restrict__ af,
    const int* __restrict__ an,
    const unsigned char* __restrict__ aull,
    const int* __restrict__ UU,
    const int* __restrict__ vmap,
    const short* __restrict__ probes,
    const short* __restrict__ cprobes,
    long long* __restrict__ codes,
    long long* __restrict__ fsums,
    unsigned long long* ctr,
    const long long cap
) {
    const int pp = p * p;
    const unsigned long long ntot = (unsigned long long)nci * ncombo;
    const unsigned long long gid =
        (unsigned long long)blockIdx.x * (unsigned long long)blockDim.x
        + (unsigned long long)threadIdx.x;
    const unsigned long long stride =
        (unsigned long long)gridDim.x * (unsigned long long)blockDim.x;
    int idx[8];
    for (unsigned long long t = gid; t < ntot; t += stride) {
        const int ci_off = (int)(t / ncombo);
        unsigned long long combo = t - (unsigned long long)ci_off * ncombo;
        const int ci = u_lo + ci_off;
        int ok = 1;
        for (int j = 0; j < k; ++j) {
            const int u = UU[ci * k + j];
            if (!aull[j * p + u]) { ok = 0; break; }
        }
        if (!ok) continue;
        int lens_ok = 1;
        for (int j = 0; j < k - 1; ++j) {
            idx[j] = (int)(combo % (unsigned long long)p);
            combo /= (unsigned long long)p;
            const int u = UU[ci * k + j];
            const int ln = an[j * p + u];
            if (idx[j] >= ln) { lens_ok = 0; break; }
        }
        if (!lens_ok) continue;
        int vs = 0;
        int fs = 0;
        for (int j = 0; j < k - 1; ++j) {
            const int u = UU[ci * k + j];
            const int base = (j * p + u) * p + idx[j];
            vs += av[base];
            fs += af[base];
        }
        int vlast = c0 - vs;
        vlast %= p;
        if (vlast < 0) vlast += p;
        const int ul = UU[ci * k + (k - 1)];
        const int pos = vmap[((k - 1) * p + ul) * p + vlast];
        if (pos < 0) continue;
        const int ftot = fs + af[((k - 1) * p + ul) * p + pos];
        int good = 1;
        for (int a = 0; a < npr; ++a) {
            int ps = 0;
            int cv = 0;
            for (int j = 0; j < k; ++j) {
                const int u = UU[ci * k + j];
                const int v = (j < k - 1) ? av[(j * p + u) * p + idx[j]] : vlast;
                const int uv = u * p + v;
                const int po = (j * pp + uv) * npr + a;
                ps += (int)probes[po];
                cv += (int)cprobes[po];
            }
            if (ftot == 0) {
                if (ps != thi && ps != tlo) { good = 0; break; }
            } else {
                const int d = ps - thi;
                int g = d / two_p;
                if (d < 0 && (d % two_p) != 0) g -= 1;
                if (g < -1 || g > ftot || g > cv) { good = 0; break; }
            }
        }
        if (!good) continue;
        long long code = (long long)ci * 16777216LL; /* 16**6 */
        long long place = 1;
        for (int j = 0; j < k - 1; ++j) {
            code += (long long)idx[j] * place;
            place *= 16;
        }
        code += (long long)pos * place;
        const unsigned long long slot = atomicAdd(ctr, 1ULL);
        if ((long long)slot < cap) {
            atomicExch((unsigned long long*)(codes + slot), (unsigned long long)code);
            atomicExch((unsigned long long*)(fsums + slot), (unsigned long long)ftot);
        }
    }
}
"""

_kern = None


def _kernel():
    global _kern
    if _kern is None:
        import cupy as cp
        _kern = cp.RawKernel(_KERNEL, "emit_cands")
    return _kern


def prepare_emit_tables(p, k, av, af, an, aull, UU, probes, cprobes):
    """Host→device once per outer. Returns a dict of cupy arrays + ncombo."""
    import cupy as cp
    from gpu_inner_fast import _vmap_from_av

    ncombo = 1
    for _ in range(k - 1):
        ncombo *= int(p)
    vmap = _vmap_from_av(p, k, av, an)
    return dict(
        ncombo=ncombo,
        av=cp.asarray(np.ascontiguousarray(av, dtype=np.int32)),
        af=cp.asarray(np.ascontiguousarray(af, dtype=np.int32)),
        an=cp.asarray(np.ascontiguousarray(an, dtype=np.int32)),
        aull=cp.asarray(np.ascontiguousarray(aull, dtype=np.uint8)),
        UU=cp.asarray(np.ascontiguousarray(UU, dtype=np.int32)),
        vmap=cp.asarray(np.ascontiguousarray(vmap, dtype=np.int32)),
        probes=cp.asarray(np.ascontiguousarray(probes, dtype=np.int16)),
        cprobes=cp.asarray(np.ascontiguousarray(cprobes, dtype=np.int16)),
    )


def emit_chunk_device(
    p, k, c0, thi, tlo, u_lo, u_hi, tab, codes_g, fsums_g, ctr,
):
    """Atomic emit into preallocated device buffers. Returns n_emitted."""
    import cupy as cp

    nci = int(u_hi - u_lo)
    if nci <= 0:
        return 0
    cap = int(codes_g.shape[0])
    ctr.fill(0)
    ncombo = tab["ncombo"]
    ntot = nci * ncombo
    tpb = 256
    nblocks = int(min((ntot + tpb - 1) // tpb, 262144))
    nblocks = max(nblocks, 1)
    npr = int(tab["probes"].shape[2])
    _kernel()(
        (nblocks,), (tpb,),
        (
            np.int32(p), np.int32(k), np.int32(u_lo), np.int32(nci),
            np.int32(c0), np.int32(thi), np.int32(tlo), np.int32(npr),
            np.int32(2 * p), np.uint64(ncombo),
            tab["av"], tab["af"], tab["an"], tab["aull"], tab["UU"], tab["vmap"],
            tab["probes"], tab["cprobes"],
            codes_g, fsums_g, ctr, np.int64(cap),
        ),
    )
    cp.cuda.runtime.deviceSynchronize()
    ntot_out = int(ctr.get()[0])
    if ntot_out > cap:
        raise RuntimeError(f"candidate overflow {ntot_out} > {cap}")
    return ntot_out


def gen_candidates_device(
    p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo,
    codes, fsums,
):
    """Fill host codes/fsums via device atomics. Returns n_emitted."""
    import cupy as cp

    tab = prepare_emit_tables(p, k, av, af, an, aull, UU, probes, cprobes)
    cap = int(codes.shape[0])
    codes_g = cp.zeros(cap, dtype=cp.int64)
    fsums_g = cp.zeros(cap, dtype=cp.int64)
    ctr = cp.zeros(1, dtype=cp.uint64)
    n = emit_chunk_device(
        p, k, c0, thi, tlo, u_lo, u_hi, tab, codes_g, fsums_g, ctr,
    )
    if n:
        codes[:n] = cp.asnumpy(codes_g[:n])
        fsums[:n] = cp.asnumpy(fsums_g[:n])
    return n
