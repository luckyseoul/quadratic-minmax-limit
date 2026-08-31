#!/usr/bin/env python3
"""Census: F2[D]-ideal of U-differences in W_0, and translation-stay generators.

Does not flip flags. Writes evidence/walsh_U_ideal_probe.json.
Backend: ProcessPool over p=3,5,7; F2 rref serial per prime. GPU unused.
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15606 import _W_basis  # noqa: E402
from e1_gmin_m4_prop15610 import _D_matrix, _dil_fn  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402


def _poly_gcd(a: list[int], b: list[int]) -> list[int]:
    def norm(p):
        p = list(p)
        while p and p[-1] == 0:
            p.pop()
        return p

    a, b = norm(a), norm(b)
    if len(a) < len(b):
        a, b = b, a
    while b:
        if len(a) < len(b):
            a, b = b, a
        shift = len(a) - len(b)
        for i, c in enumerate(b):
            a[i + shift] ^= c
        a = norm(a)
    return a or [0]


def _annihilator(x: np.ndarray, Dmat: np.ndarray) -> list[int]:
    """Minpoly of D on F2[D]x, as low-degree-first bits (monic)."""
    k = Dmat.shape[0]
    cols = [x.copy()]
    cur = x.copy()
    for d in range(1, k + 1):
        cur = (Dmat.astype(np.int32) @ cur.astype(np.int32) % 2).astype(np.uint8)
        M = np.stack(cols + [cur], axis=1)
        r = gf2_rref(M.copy())[2]
        if r <= len(cols):
            A = np.stack(cols, axis=1)
            # solve A c = cur
            Aug = np.concatenate(
                [A, cur.reshape(-1, 1)], axis=1
            )
            R, pivots, _ = gf2_rref(Aug.copy())
            c = np.zeros(len(cols), dtype=np.uint8)
            for i, pv in enumerate(pivots):
                if pv < len(cols):
                    c[pv] = R[i, -1]
            poly = [int(c[i]) for i in range(len(cols))] + [1]
            return poly
        cols.append(cur.copy())
    return [0] * k + [1]


def _im_power_member(x: np.ndarray, A: np.ndarray, t: int) -> bool:
    """True if x is in im(A^t)."""
    if t <= 0:
        return True
    At = np.eye(A.shape[0], dtype=np.uint8)
    for _ in range(t):
        At = (A.astype(np.int32) @ At.astype(np.int32) % 2).astype(np.uint8)
    M = np.concatenate([At, x.reshape(-1, 1)], axis=1)
    return gf2_rref(M.copy())[2] == gf2_rref(At.copy())[2]


def probe_one(p: int) -> dict:
    t0 = time.time()
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    n = B.shape[1]
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    Umask = fe < 0
    BU = B[Umask]
    nU = int(Umask.sum())
    direc = (BU ^ BU[0]) & 1
    dim_dirU = int(gf2_rref(direc.copy())[2])
    N = (p * p - 1) // 2
    target_W0 = N - 1
    target_V = N  # ⟨1⟩ ⊕ W_0

    WB, q, mul, add, chi = _W_basis(p)
    extra = (_qr_qnr(p)[4] if p % 4 == 1 else _qr_qnr(p)[5])[1:]
    omega = _primitive(mul, q)
    g = mul(omega, omega)
    Dmat = _D_matrix(WB, mul, g, q)
    k = WB.shape[1]
    Am = (Dmat.astype(np.int32) + np.eye(k, dtype=np.int32)) % 2

    def to_W0_coords(v_full: np.ndarray) -> np.ndarray:
        """v in F2^{P¹}; W_0-component in W-basis. v_∞ is coord 0."""
        wfn = np.zeros(q, dtype=np.uint8)
        # field e at column 1+e
        wfn[:] = v_full[1 : 1 + q]
        # strip ⟨1⟩: if v_∞=1, xor all-ones on P¹ including field
        if v_full[0]:
            wfn ^= 1
        # now w_∞=0 conceptually; W-basis already on F_q
        # express wfn in WB columns
        Aug = np.concatenate([WB, wfn.reshape(-1, 1)], axis=1)
        R, pivots, rank = gf2_rref(Aug.copy())
        x = np.zeros(k, dtype=np.uint8)
        if rank > k:
            return None  # not in W
        for i, pv in enumerate(pivots):
            if pv < k:
                x[pv] = R[i, k]
        return x

    # extra in dir(U)? extra as full vector: 0 at ∞, extra on F_q
    extra_full = np.zeros(n, dtype=np.uint8)
    extra_full[1 : 1 + q] = extra
    extra_in = gf2_rref(np.concatenate([direc.T, extra_full.reshape(-1, 1)], axis=1))[
        2
    ] == dim_dirU

    ones = np.ones(n, dtype=np.uint8)
    ones_in = gf2_rref(np.concatenate([direc.T, ones.reshape(-1, 1)], axis=1))[2] == dim_dirU

    # sample U-differences: first 64 (or all) xor base
    sample_n = min(64, nU)
    anns = []
    val_x1 = []
    in_im = []
    W0_cols = []
    for i in range(sample_n):
        d = (BU[i] ^ BU[0]) & 1
        x = to_W0_coords(d)
        if x is None:
            continue
        W0_cols.append(x)
        if x.max() == 0:
            val_x1.append(-1)  # pure ⟨1⟩
            continue
        anns.append(_annihilator(x, Dmat))
        v = 0
        while v < 8 and _im_power_member(x, Am, v + 1):
            v += 1
        val_x1.append(v)
        in_im.append(bool(_im_power_member(x, Am, 1)))

    gcd = [0, 1]
    if anns:
        gcd = anns[0]
        for a in anns[1:]:
            gcd = _poly_gcd(gcd, a)

    if W0_cols:
        dim_sample = int(gf2_rref(np.stack(W0_cols, axis=1).copy())[2])
    else:
        dim_sample = 0

    # translation-stay: for each a≠0, rows with B_0 == B_{-a}
    # (T_a y)_z = y_{z-a}; ℓ(T_a y)= y_∞ xor y_{-a}; in U iff y_0 = y_{-a}
    stay_dims = []
    stay_cols = []
    rng = np.random.default_rng(0)
    y0 = BU[0]
    A_y_count = 0
    for a in range(1, q):
        ma = 0
        # -a
        neg = 0
        for _ in range(p * p):  # add a, (p^2-1) times → -a
            # just compute add-inverse: find b: add(a,b)=0
            pass
        # inverse of a under add: (p-a%p) + (p-(a//p))%p * p with wrap
        # field add is componentwise mod p
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        col_neg = 1 + neg  # bits at field -a
        agree = BU[:, 1] == BU[:, col_neg]
        A_y_count += int(y0[1] == y0[col_neg])
        idx = np.flatnonzero(agree)
        if idx.size == 0:
            continue
        # sample up to 8 differences y xor T_a y
        take = idx if idx.size <= 8 else rng.choice(idx, 8, replace=False)
        # T_a permutation of coordinates
        perm = np.arange(n)
        perm[0] = 0
        for z in range(q):
            perm[1 + add(z, a)] = 1 + z  # (T_a y)_{z+a} = y_z so source col 1+z -> dest 1+add(z,a)
        # wait: (T_a y)_w = y_{w-a}, dest w = 1+add(z,a) has source 1+z
        # perm as source index for dest: dest_vec = src_vec[perm] with perm[dest]=src
        psrc = np.arange(n)
        psrc[0] = 0
        for z in range(q):
            psrc[1 + add(z, a)] = 1 + z
        for j in take:
            y = BU[j]
            ty = y[psrc]
            d = (y ^ ty) & 1
            x = to_W0_coords(d)
            if x is not None and x.max():
                stay_cols.append(x)
    if stay_cols:
        dim_stay = int(gf2_rref(np.stack(stay_cols, axis=1).copy())[2])
    else:
        dim_stay = 0

    # single-y translation span for BU[0]
    single_cols = []
    y = BU[0]
    for a in range(1, q):
        c0, c1 = a % p, a // p
        neg = ((p - c0) % p) + ((p - c1) % p) * p
        if y[1] != y[1 + neg]:
            continue
        psrc = np.arange(n)
        psrc[0] = 0
        for z in range(q):
            psrc[1 + add(z, a)] = 1 + z
        d = (y ^ y[psrc]) & 1
        x = to_W0_coords(d)
        if x is not None and x.max():
            single_cols.append(x)
    dim_single = (
        int(gf2_rref(np.stack(single_cols, axis=1).copy())[2]) if single_cols else 0
    )

    gcd_deg = len(gcd) - 1 if gcd and gcd[-1] else -1
    gcd_is_1 = gcd == [1]

    return {
        "p": p,
        "n": n,
        "nU": nU,
        "N": N,
        "dim_dirU": dim_dirU,
        "target_V": target_V,
        "ones_in_dirU": bool(ones_in),
        "extra_in_dirU": bool(extra_in),
        "sample_n": sample_n,
        "dim_sample_W0": dim_sample,
        "target_W0": target_W0,
        "val_x1_sample": val_x1[:20],
        "min_val_x1": min((v for v in val_x1 if v >= 0), default=None),
        "n_not_in_im_DI": int(sum(1 for v in val_x1 if v == 0)),
        "n_in_im_DI": int(sum(1 for v in val_x1 if v and v > 0)),
        "gcd_ann_is_1": bool(gcd_is_1),
        "gcd_ann_deg": gcd_deg,
        "gcd_ann": gcd[:16],
        "v2_N": _v2(N),
        "A_y0_count": A_y_count,
        "dim_stay_W0": dim_stay,
        "n_stay_cols": len(stay_cols),
        "dim_single_y_stay": dim_single,
        "n_single_cols": len(single_cols),
        "seconds": round(time.time() - t0, 3),
    }


def main():
    primes = (3, 5, 7)
    print("Walsh U-ideal probe", flush=True)
    rows = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        futs = {ex.submit(probe_one, p): p for p in primes}
        for fut in futs:
            rec = fut.result()
            rows[str(rec["p"])] = rec
            print(
                f"  p={rec['p']} dirU={rec['dim_dirU']}/{rec['target_V']} "
                f"sampleW0={rec['dim_sample_W0']}/{rec['target_W0']} "
                f"gcd1={rec['gcd_ann_is_1']} minval={rec['min_val_x1']} "
                f"stay={rec['dim_stay_W0']} single={rec['dim_single_y_stay']} "
                f"extra={rec['extra_in_dirU']} t={rec['seconds']}s",
                flush=True,
            )
    out = {
        "prop": "walsh_U_ideal_probe",
        "flags_not_flipped": True,
        "rows": rows,
        "backend": "ProcessPool W=3 over primes; F2 rref serial per prime; GPU unused",
    }
    dest = ROOT / "evidence" / "walsh_U_ideal_probe.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
