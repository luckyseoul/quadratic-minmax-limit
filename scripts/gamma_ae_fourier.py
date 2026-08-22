#!/usr/bin/env python3
"""Identify A_e ⊂ PSL principal series and invert λ_α = ⟨Γ, χ_α⟩.

Uses the no-Frobenius (signed PSL) half of Aut(C), canonical lifts with
d[0]=+1.  χ_Z from Sym²(V_+) − (1+St) (15.589 A/B).  Principal series
ρ(α_k) of PSL(2,q) are even characters α_k of F_q^* modulo inverse.

Backend: p=5 CPU; p=7 Γ on V100 via CuPy.  No flags flipped.
"""
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15590 import MuLab, field_ops, paley_conference  # noqa: E402
from gamma_class_parameter import (  # noqa: E402
    chi_W_of,
    classify_perm,
    enum_aut,
    field_inv,
    frob_fn,
    gamma_cpu,
    gamma_gpu,
)


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


def sqrt_table(q, fmul):
    root = {0: 0}
    for x in range(1, q):
        sq = fmul(x, x)
        root.setdefault(sq, x)
    return root


def psl_ps_indices(q):
    """Even k in (0, (q-1)/2), excluding the quadratic k=(q-1)/2."""
    half = (q - 1) // 2
    return [k for k in range(2, half, 2)]


def compose(e1, e2):
    p1, d1 = np.asarray(e1[0]), np.asarray(e1[1])
    p2, d2 = np.asarray(e2[0]), np.asarray(e2[1])
    return (tuple(p2[p1]), tuple(d1 * d2[p1]))


def sl_t(cl, fmul, fadd, fneg, one, finv, sqrt, dlog):
    """Torus parameter t ∈ F_q^* for a split/involution PSL class, or None."""
    if cl["family"] not in ("split", "identity", "parabolic"):
        return None
    a, b, c, d = cl["A"]
    det = cl["det"]
    tr = cl["tr"]
    if det not in sqrt:
        return None
    u = sqrt[det]
    tr_sl = fmul(tr, finv[u])
    two = fadd(one, one)
    inv2 = finv[two]
    four = fadd(two, two)
    disc = fadd(fmul(tr_sl, tr_sl), fneg(four))  # tr² − 4
    if disc == 0:
        t = fmul(tr_sl, inv2)  # ±1
        return int(t)
    if disc not in sqrt:
        return None
    s = sqrt[disc]
    t = fmul(fadd(tr_sl, s), inv2)
    if t == 0:
        t = fmul(fadd(tr_sl, fneg(s)), inv2)
    return int(t) if t != 0 else None


def chi_ps(k, q, family, t, dlog):
    """ρ(α_k) on a PSL class.  t is the SL torus parameter or None."""
    if family == "identity":
        return q + 1.0
    if family == "parabolic":
        return 1.0
    if family == "elliptic":
        return 0.0
    if family == "split" and t is not None and t > 0 and dlog[t] >= 0:
        v = dlog[t]
        ang = 2 * math.pi * k * v / (q - 1)
        return 2.0 * math.cos(ang)
    return 0.0


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    q = p * p
    n = q + 1
    print(f"p={p} q={q} n={n}  |A_e|={(q-9)//8}  #PS={(q-5)//4}", flush=True)
    fmul, fadd, fneg, one = field_ops(p)
    finv = field_inv(p, q, fmul, one)
    frob = frob_fn(p, fmul)
    squares = {fmul(x, x) for x in range(1, q)}
    sqrt = sqrt_table(q, fmul)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    ks = psl_ps_indices(q)
    print(f"primitive gen={gen}  PS indices k={ks}", flush=True)

    C = paley_conference(p)
    lab = MuLab(p, with_deg6=False)
    Y = lab.Yp.astype(np.int8)
    print(f"|Max+|={len(Y)}", flush=True)
    elems = enum_aut(p, C)
    print(f"|G+|={len(elems)}", flush=True)
    gam = gamma_cpu(Y, elems, n) if p < 7 else gamma_gpu(Y, elems, n)
    print(f"Γ(e)={gam[0]:.4f}", flush=True)

    eidx = {e: i for i, e in enumerate(elems)}
    rows = []
    for i, (pi, d) in enumerate(elems):
        cl = classify_perm(pi, p, q, n, fmul, fadd, fneg, one, finv, frob, squares)
        if not cl.get("ok") or cl["use_frob"]:
            continue
        if d[0] != 1:
            continue
        chi = chi_W_of(pi, d, C, p)
        g2 = compose((pi, d), (pi, d))
        j2 = eidx.get(g2)
        chi2 = chi_W_of(*elems[j2], C, p) if j2 is not None else float("nan")
        t = sl_t(cl, fmul, fadd, fneg, one, finv, sqrt, dlog)
        rows.append((cl, float(gam[i]), chi, chi2, t, pi, d))
    print(f"canonical PSL lifts: {len(rows)} (expect |PSL|={q*(q*q-1)//2})", flush=True)

    # χ_Z = (χ_W² + χ_W(g²))/2 − nfix
    chiZ = []
    for cl, g, chi, chi2, t, pi, d in rows:
        chiZ.append(0.5 * (chi * chi + chi2) - cl["nfix"])
    chiZ = np.array(chiZ)
    Gsz = len(rows)

    print("\n=== ⟨χ_Z, ρ_k⟩  (A_e iff ≈1)  and  λ_k = ⟨Γ, ρ_k⟩ ===", flush=True)
    Ae = []
    lams = {}
    for k in ks:
        cz = 0.0
        cg = 0.0
        for (cl, g, chi, chi2, t, pi, d), z in zip(rows, chiZ):
            x = chi_ps(k, q, cl["family"], t, dlog)
            cz += z * x
            cg += g * x
        inn_z = cz / Gsz
        lam = cg / Gsz
        mark = "A_e" if abs(inn_z - 1) < 0.05 else ("0" if abs(inn_z) < 0.05 else "?")
        print(f"  k={k:4d}  ⟨χ_Z,ρ⟩={inn_z:8.4f}  λ=⟨Γ,ρ⟩={lam:12.6f}  {mark}", flush=True)
        if mark == "A_e":
            Ae.append(k)
            lams[k] = lam

    print(f"\nA_e k's: {Ae}   |A_e| computed={len(Ae)}  formula={(q-9)//8}", flush=True)
    print(f"  4|k hypothesis: all 4|k in PS = {[k for k in ks if k % 4 == 0]}", flush=True)
    print(f"  MATCH 4|k: {Ae == [k for k in ks if k % 4 == 0]}", flush=True)

    # W_e-in-Z character is χ_Z − ∑_{A_e} ρ, not tr(U|_{V_+})
    chiW = []
    for (cl, g, chi, chi2, t, pi, d), z in zip(rows, chiZ):
        w = z
        for k in Ae:
            w -= chi_ps(k, q, cl["family"], t, dlog)
        chiW.append(w)
    chiW = np.array(chiW)
    gamP = np.array([g for cl, g, chi, chi2, t, pi, d in rows])
    lexc = float(np.dot(gamP, chiW) / Gsz)
    nrm = float(np.dot(chiW, chiW) / Gsz)
    census = {5: 176 / 13, 7: 4320 / 409, 11: 8.664378396284}
    print(f"⟨χ_W_Z, χ_W_Z⟩={nrm:.4f} (expect 1)", flush=True)
    print(f"⟨Γ, χ_W_Z⟩={lexc:.6f}   census λ_exc={census.get(p)}", flush=True)

    print("\n=== reconstruction Γ vs λ_exc χ_W_Z + ∑_{A_e} λ_k χ_k ===", flush=True)
    by_fam = defaultdict(lambda: [0, 0])
    nchk = nbad = 0
    by_tau = defaultdict(list)
    for i, (cl, g, chi, chi2, t, pi, d) in enumerate(rows):
        rec = lexc * chiW[i]
        for k in Ae:
            rec += lams[k] * chi_ps(k, q, cl["family"], t, dlog)
        nchk += 1
        bad = abs(rec - g) > 0.05
        nbad += bad
        by_fam[cl["family"]][0] += 1
        by_fam[cl["family"]][1] += bad
        if cl["family"] == "split":
            by_tau[cl["tau"]].append((g, rec, chiW[i], t, dlog[t] if t else None))
    print(f"  all PSL checked={nchk}  mismatches={nbad}  by family {dict(by_fam)}", flush=True)
    for tau in sorted(by_tau):
        vs = by_tau[tau]
        g0, rec0, w0, t0, v0 = vs[0]
        print(
            f"  τ={tau:6d} t={t0} dlog={v0}  Γ={g0:.6f} rec={rec0:.6f}  "
            f"χ_WZ={w0:.2f}  n={len(vs)}"
        )

    # pattern in A_e
    print("\n=== A_e arithmetic ===", flush=True)
    print(f"  k mod (p-1)={p-1}: {[k % (p-1) for k in Ae]}")
    print(f"  k mod (p+1)={p+1}: {[k % (p+1) for k in Ae]}")
    print(f"  k / (p-1): {[k / (p-1) for k in Ae]}")
    print(f"  k even: {all(k % 2 == 0 for k in Ae)}")
    print(f"  4|k: {all(k % 4 == 0 for k in Ae)}  count 4|k in all PS: {sum(k%4==0 for k in ks)}")


if __name__ == "__main__":
    main()
