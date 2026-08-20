"""GPU inner engine: numba generates packed candidate codes, cupy tests them.

Per outer:
  CPU (numba): odometer over (uu rows x v-lists), vlast lookup; emit packed
               codes (ci * 16^6 + v-index nibbles). No testing.
  GPU (cupy):  unpack codes, gather per-profile contribution vectors
               (int16, length q), sum -> PS.
               f0 candidates: pass iff all PS in {tlo, thi}.
               flip candidates: g = (PS - thi)/2p must satisfy
                 -1 <= g <= F pointwise (Xi in {g,g+1} subset [0,F])
                 plus per-direction line sums: p*xi_j(s) = sum_line Xi - F_j
                 with Xi in [g, g+1] => interval must meet {0, p}.
  CPU: surviving flips -> exact flipsolve.
V100/CUDA13 constraints respected: elementwise, gather, matmul only;
reductions via matmul with ones or float32 matmul.
"""
import numpy as np
from numba import njit
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen5 import _prep_tables
from kgen4 import _flip_solve
from flipnb import flip_batch

import cupy as cp
import os

_nw = int(os.environ.get("GPU_WORKERS", "3"))
cp.get_default_memory_pool().set_limit(size=int(15 * 1024**3 / _nw))


@njit(cache=True)
def _gen_candidates(p, k, av, af, an, aull, UU, c0, u_lo, u_hi,
                    codes, fsums, count, probes, cprobes, thi, tlo):
    """Emit packed candidate codes passing probe-point prefilters.
    probes: (k, p*p, npr) int16 contribution values at npr fixed points.
    f0 candidates must hit {thi, tlo} at every probe; flip candidates must
    satisfy the g-bounds -1 <= (PS - thi)/(2p) <= F at every probe."""
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
                # probe prefilter
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


class GpuTester:
    """Holds per-subset constants on GPU; tests candidate batches per outer."""

    def __init__(self, p, k, Tm, UU, eps=1):
        self.p, self.k, self.eps = p, k, eps
        q = Tm.shape[1]
        self.q = q
        self.thi = (k - 1) * eps + p
        self.tlo = (k - 1) * eps - p
        self.Tm = Tm
        self.UU_g = cp.asarray(UU.astype(np.int32))
        # line indicator matrices for the k directions: q x p float32
        L = np.zeros((k, q, p), dtype=np.float32)
        for j in range(k):
            for x in range(q):
                L[j, x, Tm[j, x]] = 1.0
        self.L_g = cp.asarray(L)
        self.ones_q = cp.ones((q, 1), dtype=cp.float32)
        p_ = p
        LIDX = np.zeros((k, p_, p_), dtype=np.int64)
        for j in range(k):
            for s in range(p_):
                LIDX[j, s] = np.where(Tm[j] == s)[0]
        self.LIDX = LIDX

    def load_outer(self, av, af, bases):
        """Build contribution dictionary: cont[j, u*p+v, x] = sigma value
        (std lift) of profile j with (u, v) at point x; also f table.
        Returns the host cont array (for probe slicing)."""
        p, k, q = self.p, self.k, self.q
        cont = np.zeros((k, p * p, q), dtype=np.int16)
        for j in range(k):
            for u in range(p):
                W = bases[j, u]              # length p
                sig = 2 * ((W[None, :] + np.arange(p)[:, None]) % p) - p + 2
                # sig[v, s]; map to points
                cont[j, u * p:(u + 1) * p, :] = sig[:, self.Tm[j]].astype(np.int16)
        cov = np.zeros((k, p * p, q), dtype=np.int8)
        for j in range(k):
            for u in range(p):
                W = bases[j, u]
                ind = (((W[None, :] + np.arange(p)[:, None]) % p) == p - 1)
                cov[j, u * p:(u + 1) * p, :] = ind[:, self.Tm[j]].astype(np.int8)
        self.cont_g = cp.asarray(cont)
        self.cov_g = cp.asarray(cov)
        self.av_g = cp.asarray(av.astype(np.int32))
        self.af_g = cp.asarray(af.astype(np.int32))
        return cont, cov

    def test_batch(self, codes, fsums):
        """codes int64 (B,), fsums int32 (B,). Returns (f0_pass_idx, flip_pass_idx)
        as host arrays of indices into the batch."""
        p, k, q = self.p, self.k, self.q
        B = len(codes)
        cg = cp.asarray(codes)
        fg = cp.asarray(fsums.astype(np.int32))
        ci = (cg // (16 ** 6)).astype(cp.int32)
        PS = cp.zeros((B, q), dtype=cp.int16)
        for j in range(k):
            vidx = ((cg // (16 ** j)) % 16).astype(cp.int32)
            uj = self.UU_g[ci, j]
            vj = self.av_g[j, uj, vidx]
            rows = uj * p + vj
            PS += self.cont_g[j][rows]
        isf0 = fg == 0
        okpt = (PS == self.thi) | (PS == self.tlo)
        # all along axis=1 without reductions: matmul trick
        cnt = okpt.astype(cp.float32) @ self.ones_q      # B x 1
        f0_pass = isf0 & (cnt[:, 0] == float(q))
        # flips: g bounds
        diff = PS - self.thi
        g = diff // (2 * p)
        COV = cp.zeros((B, q), dtype=cp.int8)
        for j in range(k):
            vidx = ((cg // (16 ** j)) % 16).astype(cp.int32)
            uj = self.UU_g[ci, j]
            vj = self.av_g[j, uj, vidx]
            rows = uj * p + vj
            COV += self.cov_g[j][rows]
        gok = (g >= -1) & (g <= fg[:, None]) & (g <= COV.astype(cp.int16))
        gcnt = gok.astype(cp.float32) @ self.ones_q
        flip_ok = (~isf0) & (gcnt[:, 0] == float(q))
        # line test on flip survivors only (keep full-batch math; cheap)
        gf = g.astype(cp.float32)
        line_ok = cp.ones(B, dtype=cp.bool_)
        F = fg.astype(cp.float32)
        for j in range(k):
            Gl = gf @ self.L_g[j]                         # B x p
            vidx = ((cg // (16 ** j)) % 16).astype(cp.int32)
            uj = self.UU_g[ci, j]
            fj = self.af_g[j, uj, vidx].astype(cp.float32)
            Fj = F - fj
            lo = Gl - Fj[:, None]
            hi = Gl + p - Fj[:, None]
            can0 = (lo <= 0) & (hi >= 0)
            can1 = (lo <= p) & (hi >= p)
            okl = can0 | can1
            okcnt = okl.astype(cp.float32) @ cp.ones((p, 1), dtype=cp.float32)
            line_ok &= (okcnt[:, 0] == float(p))
        flip_pass = flip_ok & line_ok
        f0_idx = cp.asnumpy(cp.where(f0_pass)[0])
        fl_idx = cp.asnumpy(cp.where(flip_pass)[0])
        del PS, COV, okpt, cnt, diff, g, gok, gcnt, gf, line_ok
        cp.get_default_memory_pool().free_all_blocks()
        return f0_idx, fl_idx


def process_outer_gpu(p, k, q, upper, UU, Tm, c0, eps, tester, sols,
                      gen_cap=40_000_000):
    """Full outer: streamed generation (uu-chunks) -> GPU test -> resolve."""
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
            assert (((PS == thi) | (PS == tlo)).all(axis=1)).all()
            for r in range(len(f0_idx)):
                sols.append(np.where(PS[r] == thi, 1, -1).astype(np.int8))
        if len(fl_idx):
            _resolve_flips(cc[fl_idx], ff[fl_idx])

    def _resolve_flips(cc_sub, ff_sub):
        # splits automatically if a chunk yields more solutions than the
        # output buffer holds, instead of crashing (larger p can produce
        # far more solutions per candidate batch than smaller p did)
        SIG, FV = decode(cc_sub)
        cap = max(200_000, 8 * len(cc_sub) + 1000)
        ybuf = np.zeros((cap, q), np.int8)
        ycount = np.zeros(1, np.int64)
        flip_batch(p, k, q, SIG, FV, Tm, tester.LIDX, thi, tlo,
                   ybuf, ycount, 5_000_000)
        ny = int(ycount[0])
        if ny > cap:
            if len(cc_sub) <= 1:
                raise RuntimeError(f"flip ybuf overflow on a single candidate "
                                    f"(ny={ny}); pathological, not a sizing issue")
            mid = len(cc_sub) // 2
            _resolve_flips(cc_sub[:mid], ff_sub[:mid])
            _resolve_flips(cc_sub[mid:], ff_sub[mid:])
            return
        for r in range(ny):
            sols.append(ybuf[r].copy())

    # stream over uu-chunks so a single huge outer cannot overflow buffers
    worst_per_uu = 1
    for _ in range(k - 1):
        worst_per_uu *= p
    UCH = max(1, min(2000, gen_cap // max(1, worst_per_uu)))
    ncand = 0
    for ulo in range(0, UU.shape[0], UCH):
        count = np.zeros(1, np.int64)
        _gen_candidates(p, k, av, af, an, aull, UU, c0, ulo,
                        min(ulo + UCH, UU.shape[0]),
                        codes, fsums, count, probes, cprobes, thi, tlo)
        nc = int(count[0])
        if nc > gen_cap:
            raise RuntimeError(f"candidate overflow {nc}")
        ncand += nc
        for lo in range(0, nc, CH):
            resolve(codes[lo:lo + CH], fsums[lo:lo + CH])
    return ncand
