"""Faster candidate generation for gpu_inner (do not swap into the live k=6 job).

Bottleneck at p=13 k=6: serial numba odometer emits ~3e8–4.5e8 codes/outer
(~1.3e6/s, ~280s) while V100 testing is a small fraction. This module:

- precomputes vlast→slot (no linear scan)
- single-pass `prange` emit via CPU LLVM **atomic read/write**
  (fetch-add the write cursor, atomic store of packed code + fsum)
- two-pass count/fill kept as fallback if atomics fail to compile

Same packed-code format as gpu_inner._gen_candidates so GpuTester.test_batch
is unchanged. Live PID 6838 keeps the old serial inner.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange, types
from numba.core import cgutils
from numba.extending import intrinsic


# ---------------------------------------------------------------------------
# CPU LLVM atomics (numba.cuda.atomic is GPU-only; V100 is owned by k=6).
# ---------------------------------------------------------------------------

def _ary1d_item_ptr(context, builder, signature, ary_val, idx_val, idx_pos=1):
    aryty = signature.args[0]
    lary = context.make_array(aryty)(context, builder, ary_val)
    idxp = context.cast(builder, idx_val, signature.args[idx_pos], types.intp)
    return cgutils.get_item_pointer(
        context, builder, aryty, lary, [idxp], wraparound=False
    )


@intrinsic
def atomic_fetch_add_i64(typingctx, ary, idx, val):
    """ary[idx] += val; return previous value. acq_rel."""
    if not isinstance(ary, types.Array) or ary.ndim != 1:
        return None
    sig = types.int64(ary, types.int64, types.int64)

    def codegen(context, builder, signature, args):
        ary_val, idx_val, val_val = args
        ptr = _ary1d_item_ptr(context, builder, signature, ary_val, idx_val)
        v = context.cast(builder, val_val, signature.args[2], signature.args[0].dtype)
        return builder.atomic_rmw("add", ptr, v, ordering="acq_rel")

    return sig, codegen


@intrinsic
def atomic_store_i64(typingctx, ary, idx, val):
    """Atomic *release* store of int64."""
    if not isinstance(ary, types.Array) or ary.ndim != 1:
        return None
    sig = types.none(ary, types.int64, types.int64)

    def codegen(context, builder, signature, args):
        ary_val, idx_val, val_val = args
        ptr = _ary1d_item_ptr(context, builder, signature, ary_val, idx_val)
        v = context.cast(builder, val_val, signature.args[2], signature.args[0].dtype)
        builder.store_atomic(v, ptr, "release", 8)
        return context.get_dummy_value()

    return sig, codegen


@intrinsic
def atomic_load_i64(typingctx, ary, idx):
    """Atomic *acquire* load of int64."""
    if not isinstance(ary, types.Array) or ary.ndim != 1:
        return None
    sig = types.int64(ary, types.int64)

    def codegen(context, builder, signature, args):
        ary_val, idx_val = args
        ptr = _ary1d_item_ptr(context, builder, signature, ary_val, idx_val)
        return builder.load_atomic(ptr, "acquire", 8)

    return sig, codegen


def _probe_cpu_atomics():
    """Compile+run a tiny prange kernel. False if this numba cannot lower atomics."""
    try:
        @njit(parallel=True)
        def _k(n, count, out):
            for i in prange(n):
                slot = atomic_fetch_add_i64(count, 0, 1)
                atomic_store_i64(out, slot, i + 1)
            return atomic_load_i64(count, 0)

        n = 4096
        count = np.zeros(1, np.int64)
        out = np.zeros(n, np.int64)
        tot = int(_k(n, count, out))
        if tot != n or int(count[0]) != n:
            return False
        if int(out.min()) != 1 or int(out.max()) != n:
            return False
        if len(set(out.tolist())) != n:
            return False
        return True
    except Exception:
        return False


ATOMICS_AVAILABLE = _probe_cpu_atomics()


# ---------------------------------------------------------------------------
# vmap + odometer walk
# ---------------------------------------------------------------------------

@njit(cache=True)
def _gen_candidates_serial(p, k, av, af, an, aull, UU, c0, u_lo, u_hi,
                           codes, fsums, count, probes, cprobes, thi, tlo):
    """Exact copy of gpu_inner._gen_candidates (CPU, no cupy)."""
    m = count[0]
    idx = np.zeros(k, np.int64)
    lens = np.zeros(k, np.int64)
    P6 = 16 ** 6
    npr = probes.shape[2]
    two_p = 2 * p
    for ci in range(u_lo, u_hi):
        ok = True
        for j in range(k):
            if not aull[j, UU[ci, j]]:
                ok = False
                break
        if not ok:
            continue
        for j in range(k):
            idx[j] = 0
            lens[j] = an[j, UU[ci, j]]
        while True:
            vs = 0
            fs = 0
            for j in range(k - 1):
                vs += av[j, UU[ci, j], idx[j]]
                fs += af[j, UU[ci, j], idx[j]]
            vlast = (c0 - vs) % p
            jl = k - 1
            ul = UU[ci, jl]
            pos = -1
            for t in range(lens[jl]):
                if av[jl, ul, t] == vlast:
                    pos = t
                    break
            if pos >= 0:
                ftot = fs + af[jl, ul, pos]
                good = True
                for a in range(npr):
                    ps = 0
                    cv = 0
                    for j in range(k):
                        u = UU[ci, j]
                        v = av[j, u, idx[j]] if j < k - 1 else vlast
                        ps += probes[j, u * p + v, a]
                        cv += cprobes[j, u * p + v, a]
                    if ftot == 0:
                        if ps != thi and ps != tlo:
                            good = False
                            break
                    else:
                        d = ps - thi
                        g = d // two_p
                        if g < -1 or g > ftot or g > cv:
                            good = False
                            break
                if good:
                    if m < codes.shape[0]:
                        code = ci * P6
                        for j in range(k - 1):
                            code += idx[j] * (16 ** j)
                        code += pos * (16 ** (k - 1))
                        codes[m] = code
                        fsums[m] = ftot
                    m += 1
            c2 = 0
            while c2 < k - 1:
                idx[c2] += 1
                if idx[c2] < lens[c2]:
                    break
                idx[c2] = 0
                c2 += 1
            if c2 >= k - 1:
                break
    count[0] = m


@njit(cache=True)
def _vmap_from_av(p, k, av, an):
    vmap = np.full((k, p, p), np.int64(-1))
    for j in range(k):
        for u in range(p):
            n = an[j, u]
            for t in range(n):
                vmap[j, u, av[j, u, t]] = t
    return vmap


@njit(cache=True)
def _walk_ci(
    p, k, ci, av, af, an, aull, UU, c0, vmap, probes, cprobes, thi, tlo, two_p, npr,
    write, codes, fsums, base,
):
    """Emit probe-pass candidates for one UU row (no atomics).

    write=0  count only; write=1  fill codes[base + m] (disjoint ranges).
    """
    for j in range(k):
        if not aull[j, UU[ci, j]]:
            return 0
    idx = np.zeros(k, np.int64)
    lens = np.zeros(k, np.int64)
    for j in range(k):
        lens[j] = an[j, UU[ci, j]]
        if lens[j] <= 0:
            return 0
    P6 = 16 ** 6
    m = 0
    while True:
        vs = 0
        fs = 0
        for j in range(k - 1):
            vs += av[j, UU[ci, j], idx[j]]
            fs += af[j, UU[ci, j], idx[j]]
        vlast = (c0 - vs) % p
        pos = vmap[k - 1, UU[ci, k - 1], vlast]
        if pos >= 0:
            ftot = fs + af[k - 1, UU[ci, k - 1], pos]
            good = True
            for a in range(npr):
                ps = 0
                cv = 0
                for j in range(k):
                    u = UU[ci, j]
                    v = av[j, u, idx[j]] if j < k - 1 else vlast
                    ps += probes[j, u * p + v, a]
                    cv += cprobes[j, u * p + v, a]
                if ftot == 0:
                    if ps != thi and ps != tlo:
                        good = False
                        break
                else:
                    d = ps - thi
                    g = d // two_p
                    if g < -1 or g > ftot or g > cv:
                        good = False
                        break
            if good:
                if write:
                    code = ci * P6
                    for j in range(k - 1):
                        code += idx[j] * (16 ** j)
                    code += pos * (16 ** (k - 1))
                    codes[base + m] = code
                    fsums[base + m] = ftot
                m += 1
        c2 = 0
        while c2 < k - 1:
            idx[c2] += 1
            if idx[c2] < lens[c2]:
                break
            idx[c2] = 0
            c2 += 1
        if c2 >= k - 1:
            break
    return m


@njit
def _walk_ci_atomic(
    p, k, ci, av, af, an, aull, UU, c0, vmap, probes, cprobes, thi, tlo, two_p, npr,
    codes, fsums, ctr,
):
    """Same odometer as _walk_ci; slot via atomic fetch-add, payload via atomic store."""
    for j in range(k):
        if not aull[j, UU[ci, j]]:
            return 0
    idx = np.zeros(k, np.int64)
    lens = np.zeros(k, np.int64)
    for j in range(k):
        lens[j] = an[j, UU[ci, j]]
        if lens[j] <= 0:
            return 0
    P6 = 16 ** 6
    cap = codes.shape[0]
    m = 0
    while True:
        vs = 0
        fs = 0
        for j in range(k - 1):
            vs += av[j, UU[ci, j], idx[j]]
            fs += af[j, UU[ci, j], idx[j]]
        vlast = (c0 - vs) % p
        pos = vmap[k - 1, UU[ci, k - 1], vlast]
        if pos >= 0:
            ftot = fs + af[k - 1, UU[ci, k - 1], pos]
            good = True
            for a in range(npr):
                ps = 0
                cv = 0
                for j in range(k):
                    u = UU[ci, j]
                    v = av[j, u, idx[j]] if j < k - 1 else vlast
                    ps += probes[j, u * p + v, a]
                    cv += cprobes[j, u * p + v, a]
                if ftot == 0:
                    if ps != thi and ps != tlo:
                        good = False
                        break
                else:
                    d = ps - thi
                    g = d // two_p
                    if g < -1 or g > ftot or g > cv:
                        good = False
                        break
            if good:
                code = ci * P6
                for j in range(k - 1):
                    code += idx[j] * (16 ** j)
                code += pos * (16 ** (k - 1))
                slot = atomic_fetch_add_i64(ctr, 0, 1)
                if slot < cap:
                    atomic_store_i64(codes, slot, code)
                    atomic_store_i64(fsums, slot, ftot)
                m += 1
        c2 = 0
        while c2 < k - 1:
            idx[c2] += 1
            if idx[c2] < lens[c2]:
                break
            idx[c2] = 0
            c2 += 1
        if c2 >= k - 1:
            break
    return m


@njit(parallel=True, cache=True)
def _count_all(
    p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes, thi, tlo, counts
):
    npr = probes.shape[2]
    two_p = 2 * p
    dummy = np.zeros(1, np.int64)
    for ci in prange(u_lo, u_hi):
        counts[ci - u_lo] = _walk_ci(
            p, k, ci, av, af, an, aull, UU, c0, vmap, probes, cprobes, thi, tlo,
            two_p, npr, False, dummy, dummy, 0,
        )


@njit(parallel=True, cache=True)
def _fill_all(
    p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes, thi, tlo,
    starts, codes, fsums,
):
    npr = probes.shape[2]
    two_p = 2 * p
    for ci in prange(u_lo, u_hi):
        base = starts[ci - u_lo]
        _walk_ci(
            p, k, ci, av, af, an, aull, UU, c0, vmap, probes, cprobes, thi, tlo,
            two_p, npr, True, codes, fsums, base,
        )


@njit(parallel=True)
def _emit_all_atomic(
    p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes, thi, tlo,
    ctr, codes, fsums,
):
    npr = probes.shape[2]
    two_p = 2 * p
    for ci in prange(u_lo, u_hi):
        _walk_ci_atomic(
            p, k, ci, av, af, an, aull, UU, c0, vmap, probes, cprobes, thi, tlo,
            two_p, npr, codes, fsums, ctr,
        )
    return atomic_load_i64(ctr, 0)


def gen_candidates_twopass(
    p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo, codes, fsums
):
    """prange count-then-fill (no atomics). Returns n_emitted."""
    nci = u_hi - u_lo
    if nci <= 0:
        return 0
    vmap = _vmap_from_av(p, k, av, an)
    counts = np.zeros(nci, dtype=np.int64)
    _count_all(
        p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes, thi, tlo, counts
    )
    starts = np.zeros(nci + 1, dtype=np.int64)
    for i in range(nci):
        starts[i + 1] = starts[i] + counts[i]
    ntot = int(starts[nci])
    if ntot > codes.shape[0]:
        raise RuntimeError(f"candidate overflow {ntot} > {codes.shape[0]}")
    if ntot:
        _fill_all(
            p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes, thi, tlo,
            starts[:-1], codes, fsums,
        )
    return ntot


def gen_candidates_atomic(
    p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo, codes, fsums
):
    """Single-pass prange emit with atomic fetch-add + atomic store."""
    nci = u_hi - u_lo
    if nci <= 0:
        return 0
    if not ATOMICS_AVAILABLE:
        raise RuntimeError("CPU LLVM atomics not available on this numba")
    vmap = _vmap_from_av(p, k, av, an)
    ctr = np.zeros(1, dtype=np.int64)
    ntot = int(
        _emit_all_atomic(
            p, k, av, af, an, aull, UU, c0, vmap, u_lo, u_hi, probes, cprobes,
            thi, tlo, ctr, codes, fsums,
        )
    )
    if ntot > codes.shape[0]:
        raise RuntimeError(f"candidate overflow {ntot} > {codes.shape[0]}")
    return ntot


def gen_candidates_parallel(
    p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo, codes, fsums,
    mode="atomic",
):
    """Drop-in replacement for gpu_inner._gen_candidates. Returns n_emitted.

    mode='atomic'   single-pass atomic r/w (default if ATOMICS_AVAILABLE)
    mode='twopass'  count then fill
    """
    if mode == "twopass" or not ATOMICS_AVAILABLE:
        return gen_candidates_twopass(
            p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo,
            codes, fsums,
        )
    if mode != "atomic":
        raise ValueError(f"unknown mode {mode!r}")
    return gen_candidates_atomic(
        p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo,
        codes, fsums,
    )


def process_outer_gpu_fast(p, k, q, upper, UU, Tm, c0, eps, tester, sols,
                           gen_cap=40_000_000):
    """Same as gpu_inner.process_outer_gpu with parallel/atomic candidate gen."""
    from kgen5 import _prep_tables

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

    worst_per_uu = 1
    for _ in range(k - 1):
        worst_per_uu *= p
    UCH = max(1, min(2000, gen_cap // max(1, worst_per_uu)))
    ncand = 0
    for ulo in range(0, UU.shape[0], UCH):
        nc = gen_candidates_parallel(
            p, k, av, af, an, aull, UU, c0, ulo,
            min(ulo + UCH, UU.shape[0]),
            probes, cprobes, thi, tlo, codes, fsums,
        )
        ncand += nc
        for lo in range(0, nc, CH):
            resolve(codes[lo:lo + CH], fsums[lo:lo + CH])
    return ncand
