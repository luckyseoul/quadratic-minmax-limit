#!/usr/bin/env python3
"""Aut-orbit values of L(r)=E ∑_{δ≠0} N(δ)N(rδ), and a cosine model of λ(k).

15.279 Q: L/Q is constant on ⟨Frob,inv⟩-orbits of squares/{±1};
leftover dofs = n_orb−2.  This script *names* those leftover values at
p=5,7 (MuLab) and tests whether a short cosine polynomial in
θ=2πk/(q−1) that fits p=5,7 predicts the p=11 even-character spectrum
(even_char_hip_minmax.npz).

No flag flip.  Not an identity file.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from e1_gmin_m4_prop15590 import MuLab, field_ops  # noqa: E402


def primitive_root(q, fmul, one):
    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    return next(e for e in range(2, q) if order_of(e) == q - 1)


def dlog_table(gen, q, fmul, one):
    tab = [-1] * q
    x = one
    for k in range(q - 1):
        tab[x] = k
        x = fmul(x, gen)
    return tab


def add_mul_tables(q, fadd, fmul):
    Add = np.empty((q, q), dtype=np.int32)
    Mul = np.empty((q, q), dtype=np.int32)
    for i in range(q):
        for j in range(q):
            Add[i, j] = fadd(i, j)
            Mul[i, j] = fmul(i, j)
    return Add, Mul


def inv_table(q, fmul, one):
    inv = np.zeros(q, dtype=np.int32)
    for a in range(1, q):
        x = one
        # a^{q-2}
        base, ee, r = a, q - 2, one
        while ee:
            if ee & 1:
                r = fmul(r, base)
            base = fmul(base, base)
            ee >>= 1
        inv[a] = r
        x = r
    return inv


def frob_pow(e, p, fmul, one):
    """e^p."""
    if e == 0:
        return 0
    r, base, ee = one, e, p
    while ee:
        if ee & 1:
            r = fmul(r, base)
        base = fmul(base, base)
        ee >>= 1
    return r


def n_aut_orbits_squares(p: int) -> int:
    if p % 4 == 1:
        return (p + 3) ** 2 // 16
    return (p + 1) * (p + 5) // 16


def orbit_invariants(p, r, fmul, fadd, fneg, one, chi, inv, dlog, gen):
    q = p * p
    # F_p embedding: e = p*a+b ↔ a+bt, F_p is b=0, e=p*a
    a, b = divmod(r, p)
    in_fp = b == 0
    # order in F_q^×
    ord_m = (q - 1) // np.gcd(dlog[r], q - 1) if r else 0
    rm1 = fadd(r, fneg(one))
    rp1 = fadd(r, one)
    chi_rm1 = 0 if rm1 == 0 else int(chi[rm1])
    chi_rp1 = 0 if rp1 == 0 else int(chi[rp1])
    # Norm_{q/p}(x)=x^{p+1}; 15590: a+bt, t^2=nonres
    # χ_p(Norm): r^{p+1} in F_p
    nr = r
    for _ in range(p):
        nr = fmul(nr, r)
    # r^{p+1}
    nr = frob_pow(r, p, fmul, one)
    nr = fmul(nr, r)
    na, nb = divmod(nr, p)
    # Tr = r + r^p
    rp = frob_pow(r, p, fmul, one)
    tr = fadd(r, rp)
    ta, tb = divmod(tr, p)
    return {
        "r": int(r),
        "in_Fp": bool(in_fp),
        "order": int(ord_m),
        "chi_rm1": chi_rm1,
        "chi_rp1": chi_rp1,
        "norm_a": int(na),
        "norm_b": int(nb),
        "tr_a": int(ta),
        "tr_b": int(tb),
        "dlog": int(dlog[r]),
    }


def ensemble_NL(p: int):
    q = p * p
    fmul, fadd, fneg, one = field_ops(p)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    Add, Mul = add_mul_tables(q, fadd, fmul)
    inv = inv_table(q, fmul, one)
    chi = np.zeros(q, dtype=np.int8)
    sq = {fmul(t, t) for t in range(1, q)}
    for e in range(1, q):
        chi[e] = 1 if e in sq else -1
    neg1 = fneg(one)

    lab = MuLab(p, with_deg6=False)
    Z = lab.Yp.astype(np.int8)[:, 1:]
    M = len(Z)
    Dm = Z == -1
    N = np.zeros((M, q), dtype=np.int32)
    for a in range(q):
        N[:, a] = (Dm & Dm[:, Add[:, a]]).sum(axis=1)
    nD_row = N[:, 0].astype(np.int64)
    nD_plus = p * (p - 1) // 2
    nD_minus = p * (p + 1) // 2
    sizes = set(int(x) for x in np.unique(nD_row))
    assert sizes <= {nD_plus, nD_minus}, sizes

    squares = [e for e in range(1, q) if chi[e] == 1]
    # T-rep min(r,-r)
    def t_rep(r):
        nr = fmul(r, neg1)
        return min(r, nr)

    # ⟨Frob, inv⟩ orbits on T
    unused = set(t_rep(r) for r in squares)
    orbits = []
    while unused:
        start = min(unused)
        stack = [start]
        seen = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.append(t_rep(frob_pow(x, p, fmul, one)))
            stack.append(t_rep(int(inv[x])))
        unused -= seen
        orbits.append(sorted(seen))

    L_of = {}
    for r in squares:
        idx = Mul[:, r]
        pair = (N.astype(np.int64) * N[:, idx]).sum(axis=1) - nD_row * nD_row
        L_of[r] = Fraction(int(pair.sum()), M)

    print(f"\n=== p={p} |Max+|={M} n_orbits={len(orbits)} formula={n_aut_orbits_squares(p)} ===")
    print(
        f"{'rep':>4} {'size':>4} {'inFp':>4} {'ord':>5} "
        f"{'χ(r-1)':>7} {'χ(r+1)':>7} {'L':>18} {'L/q^2':>10}"
    )
    rows = []
    for orb in sorted(orbits, key=lambda o: (len(o), o[0])):
        r0 = orb[0]
        Ls = {L_of[r] for r in orb}
        Ls |= {L_of[fmul(r, neg1)] for r in orb}
        invs = orbit_invariants(p, r0, fmul, fadd, fneg, one, chi, inv, dlog, gen)
        L = L_of[r0]
        const = len(Ls) == 1
        print(
            f"{r0:4d} {len(orb):4d} {int(invs['in_Fp']):4d} {invs['order']:5d} "
            f"{invs['chi_rm1']:7d} {invs['chi_rp1']:7d} "
            f"{str(L):>18} {float(L)/(q*q):10.4f}"
            f"{'' if const else '  NOT CONST '+str(Ls)}"
        )
        rows.append(
            {
                "rep": r0,
                "size": len(orb),
                "L": str(L),
                "L_float": float(L),
                **{k: invs[k] for k in invs if k != "r"},
                "const_on_orbit": const,
            }
        )

    # even-k λ
    half = (q - 1) // 2
    ang = np.zeros((half // 2, q), dtype=np.complex128)  # k=2,4,...,half-2 plus we'll do all even < half
    ks = list(range(2, half, 2))
    E2 = []
    for k in ks:
        a = np.zeros(q, dtype=np.complex128)
        for x in range(1, q):
            a[x] = np.exp(2j * np.pi * k * dlog[x] / (q - 1))
        S = N[:, 1:].astype(np.complex128) @ a[1:]
        E2.append(float(np.mean(np.abs(S) ** 2)))
    c = 32 / (q * (q - 1))
    lams = [c * e for e in E2]
    print(f"\n  even-k λ  (dim F={(q-5)//4})")
    print(f"  {'k':>4} {'λ':>12} {'λ-8':>12} {'θ/π':>8}")
    for k, lam in zip(ks, lams):
        print(f"  {k:4d} {lam:12.6f} {lam-8:12.6f} {2*k/(q-1):8.4f}")

    return {
        "p": p,
        "M": M,
        "n_orbits": len(orbits),
        "orbits": rows,
        "ks": ks,
        "E2": E2,
        "lams": lams,
    }


def cosine_design(ks, q, n_harm):
    """Columns: 1, cos(2π m k/(q-1)) for m=1..n_harm."""
    theta = 2 * np.pi * np.asarray(ks, dtype=np.float64) / (q - 1)
    cols = [np.ones(len(ks))]
    for m in range(1, n_harm + 1):
        cols.append(np.cos(m * theta))
    return np.column_stack(cols)


def fit_and_predict():
    recs = {}
    for p in (5, 7):
        recs[p] = ensemble_NL(p)

    # p=11 from HIP npz
    npz = np.load(ROOT / "evidence" / "even_char_hip_minmax.npz")
    recs[11] = {
        "p": 11,
        "ks": [int(x) for x in npz["ks"]],
        "lams": [float(x) for x in npz["lams"]],
        "E2": [float(x) for x in npz["e2"]],
    }

    print("\n=== cosine model of λ(k): fit p=5+7, predict p=11 ===")
    # use λ-8 (Wick baseline)
    for n_harm in range(1, 7):
        Xs, ys = [], []
        for p in (5, 7):
            q = p * p
            ks = recs[p]["ks"]
            y = np.asarray(recs[p]["lams"]) - 8
            X = cosine_design(ks, q, n_harm)
            Xs.append(X)
            ys.append(y)
        X = np.vstack(Xs)
        y = np.concatenate(ys)
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        # per-prime residual
        print(f"\n  n_harm={n_harm}  coef={np.round(coef, 5)}")
        for p in (5, 7, 11):
            q = p * p
            ks = recs[p]["ks"]
            ytrue = np.asarray(recs[p]["lams"]) - 8
            pred = cosine_design(ks, q, n_harm) @ coef
            err = np.max(np.abs(pred - ytrue))
            print(f"    p={p:2d} max|λ−pred|={err:.4f}  min_true={ytrue.min()+8:.4f} min_pred={pred.min()+8:.4f}")

    # Separate cosine fit per prime (sanity: how many harmonics to interpolate)
    print("\n=== per-prime interpolating cosine (how many harmonics?) ===")
    for p in (5, 7, 11):
        q = p * p
        ks = recs[p]["ks"]
        y = np.asarray(recs[p]["lams"]) - 8
        for n_harm in range(1, min(8, len(ks))):
            X = cosine_design(ks, q, n_harm)
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coef
            err = np.max(np.abs(pred - y))
            if err < 1e-6:
                print(f"  p={p} exact at n_harm={n_harm} coef={np.round(coef, 6)}")
                break
        else:
            print(f"  p={p} not exact by n_harm≤{n_harm} last maxerr={err:.3e}")

    # Chebyshev in x=cos(4πk/(q-1))  (period compatible with 4|k and even)
    print("\n=== polynomial in x=cos(4π k/(q-1)) ===")
    for deg in range(1, 6):
        print(f"  deg={deg}")
        for p in (5, 7, 11):
            q = p * p
            ks = np.asarray(recs[p]["ks"], dtype=np.float64)
            x = np.cos(4 * np.pi * ks / (q - 1))
            y = np.asarray(recs[p]["lams"]) - 8
            X = np.column_stack([x ** m for m in range(deg + 1)])
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coef
            err = np.max(np.abs(pred - y))
            print(f"    p={p} maxerr={err:.4f}")

    # dump
    out = {
        "p5_orbits": recs[5]["orbits"],
        "p7_orbits": recs[7]["orbits"],
        "p5_lams": list(zip(recs[5]["ks"], recs[5]["lams"])),
        "p7_lams": list(zip(recs[7]["ks"], recs[7]["lams"])),
    }
    path = ROOT / "evidence" / "aut_orbit_L_p5p7.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    fit_and_predict()
