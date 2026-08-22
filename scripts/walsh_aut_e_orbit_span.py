#!/usr/bin/env python3
"""Aut_e-orbit affine span of actual Max- points in the xor-slice U.

Irreducibility on random H0-vectors is false. Special Max- points may
still have full-slice Aut_e orbits. Also: line-flip does not preserve
Max- (signs on S are determined by the exterior). Antipodal completion
restores dim H = n/2 at p=11.

No flag flip.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from walsh_linecode_rank import aut_e_generators, orbit_span_F2  # noqa: E402


def bits(Y):
    return ((1 - Y.astype(np.int8)) // 2).astype(np.uint8)


def affine_orbit_dim(gens, x):
    """Dim of affine span of {g x} = dim span{gx xor x0}."""
    d_lin, norb = orbit_span_F2(gens, x)
    # collect orbit, then dir
    seen = {x.tobytes()}
    dq = [x.copy()]
    orb = [x.copy()]
    while dq:
        u = dq.pop()
        for g in gens:
            v = u[g]
            k = v.tobytes()
            if k not in seen:
                seen.add(k)
                dq.append(v)
                orb.append(v)
    B = np.stack(orb, axis=0)
    direc = gf2_rref((B ^ B[0]) & 1)[2]
    return int(d_lin), int(direc), int(len(orb))


def line_flip_test(p: int) -> dict:
    """Negating y on S={∞}∪F_p: still Max-?"""
    from e1_gmin_m4_prop15590 import paley_conference

    # use same C as load_minus
    from minmax_quadratic import paley_conference_prime_power

    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int64)
    C = np.rint(C).astype(np.int64)
    S = np.array([0] + list(range(1, p + 1)))
    y = Y[0].copy()
    yp = y.copy()
    yp[S] *= -1
    pval = C.shape[0]
    # Cy' vs -p y'
    lhs = C @ yp
    rhs = -p * yp
    ok = bool(np.array_equal(lhs, rhs))
    # is yp in the ensemble?
    in_ens = bool(((Y == yp).all(axis=1)).any())
    return {"p": p, "flip_is_eigen": ok, "flip_in_ensemble": in_ens}


def analyse(p: int, n_pts: int = 20, seed: int = 0) -> dict:
    print(f"\n======== p={p} ========", flush=True)
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = bits(Y)
    n = B.shape[1]
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    U = fe < 0
    BU = B[U]
    gens, _, _ = aut_e_generators(p)
    print(f"  |Max-|={len(Y)} |U|={int(U.sum())} n={n}", flush=True)
    flip = line_flip_test(p)
    print(f"  line-flip eigen={flip['flip_is_eigen']} in_ens={flip['flip_in_ensemble']}", flush=True)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(BU), size=min(n_pts, len(BU)), replace=False)
    recs = []
    for i in idx:
        dlin, ddir, norb = affine_orbit_dim(gens, BU[i])
        recs.append({"lin": dlin, "dir": ddir, "norb": norb})
        print(f"    U-pt orbit lin={dlin} dir={ddir} |orb|={norb}", flush=True)
    dirs = [r["dir"] for r in recs]
    target = n // 2 - 1  # expected slice dir if H dim = n/2
    return {
        "p": p,
        "n": n,
        "nU": int(U.sum()),
        "line_flip": flip,
        "orbits": recs,
        "max_dir": max(dirs) if dirs else None,
        "min_dir": min(dirs) if dirs else None,
        "any_full_slice": any(d >= target for d in dirs),
        "all_full_slice": all(d >= target for d in dirs),
        "target_slice_dir": target,
    }


def p11_U_complements(nsamp=30000, seed=4) -> dict:
    print("\n======== p=11 U + antipodes ========", flush=True)
    from e1_gmin_m4_prop15598 import field_ctx
    from minmax_quadratic import paley_conference_prime_power
    from walsh_linecode_rank import square_line_matrix

    p, q, n = 11, 121, 122
    A = np.load("/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy", mmap_mode="r")
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(A.shape[0], size=nsamp, replace=False))
    C = paley_conference_prime_power(p)
    q_, mul, add, chi, frob, norm, ia, ib = field_ctx(p)

    def order_of(e):
        x, o = e, 1
        while x != 1:
            x = mul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    pi = np.zeros(n, dtype=np.int64)
    pi[0] = 0
    for e in range(q):
        pi[1 + e] = 1 + mul(e, gen)
    d = np.zeros(n, dtype=np.int64)
    d[0] = 1
    d[1:] = -np.rint(C[pi[0], pi[1:]]).astype(np.int64) * np.rint(C[0, 1:]).astype(
        np.int64
    )
    chunk = A[idx].astype(np.int64)
    Ym = d[None, :] * chunk[:, pi]
    B = ((1 - Ym) // 2).astype(np.uint8)
    ones = np.ones(n, dtype=np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Ym[:, 0] * Ym[:, 1]
    U = fe < 0
    BU = B[U]
    BUc = np.vstack([BU, (BU ^ ones) & 1])
    Bfull = np.vstack([B, (B ^ ones) & 1])

    def lin(M):
        return gf2_rref(M)[2] if len(M) < 80000 else None

    def direc(M):
        return gf2_rref((M ^ M[0]) & 1)[2]

    rec = {
        "nsamp": nsamp,
        "nU": int(U.sum()),
        "Max_lin_half": int(gf2_rref(B)[2]),
        "Max_dir_half": int(direc(B)),
        "Max_lin_anti": int(gf2_rref(Bfull)[2]),
        "Max_dir_anti": int(direc(Bfull)),
        "U_lin_half": int(gf2_rref(BU)[2]),
        "U_dir_half": int(direc(BU)),
        "U_lin_anti": int(gf2_rref(BUc)[2]),
        "U_dir_anti": int(direc(BUc)),
        "n_over_2": n // 2,
    }
    print(f"  {rec}", flush=True)
    return rec


def main():
    out = {"small": {}, "p11": None}
    for p in (5, 7):
        out["small"][str(p)] = analyse(p, n_pts=12)
    out["p11"] = p11_U_complements()
    dest = ROOT / "evidence" / "walsh_aut_e_orbit_span.json"
    dest.write_text(__import__("json").dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)


if __name__ == "__main__":
    main()
