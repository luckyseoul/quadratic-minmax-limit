#!/usr/bin/env python3
"""|ẑ|² on Ω for Max+, p=5 (and p=7).  15.279: ẑ supported on {0}∪Ω,
E[|ẑ|²]=2q on Ω, λ=8+R̂_rest/q².  Few-valued |ẑ|² would name R̂_rest.

No flag flip.  CPU, MuLab ensembles.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from e1_gmin_m4_prop15590 import MuLab, field_ops, paley_conference  # noqa: E402


def tr_table(p, q):
    """Tr_{F_q/F_p}(a+bt)=2a with 15590 encoding e=p*a+b."""
    tr = np.empty(q, dtype=np.int32)
    for e in range(q):
        a, b = divmod(e, p)
        tr[e] = (2 * a) % p
    return tr


def add_table(q, fadd):
    A = np.empty((q, q), dtype=np.int32)
    for i in range(q):
        for j in range(q):
            A[i, j] = fadd(i, j)
    return A


def chi_table(q, fmul, one):
    chi = np.zeros(q, dtype=np.int8)
    for e in range(1, q):
        x, k = e, 1
        # e^{(q-1)/2}
        x = e
        exp = (q - 1) // 2
        r = one
        base = e
        ee = exp
        while ee:
            if ee & 1:
                r = fmul(r, base)
            base = fmul(base, base)
            ee >>= 1
        chi[e] = 1 if r == one else -1
    return chi


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    q = p * p
    fmul, fadd, fneg, one = field_ops(p)
    C = paley_conference(p)
    # Paley chi on F_q: C[1+x,1+y] = chi(x-y)
    chi = np.zeros(q, dtype=np.int8)
    for x in range(q):
        for y in range(q):
            if x != y:
                chi[fadd(x, fneg(y))] = int(C[1 + x, 1 + y])
                break
    # chi from squares is more reliable:
    chi = np.zeros(q, dtype=np.int8)
    sq = {fmul(t, t) for t in range(1, q)}
    for e in range(1, q):
        chi[e] = 1 if e in sq else -1
    tr = tr_table(p, q)
    # additive characters e_ξ(x)=exp(2πi Tr(ξx)/p)
    Add = add_table(q, fadd)
    # Tr(ξ x)
    Trxi = np.empty((q, q), dtype=np.int32)
    for xi in range(q):
        for x in range(q):
            Trxi[xi, x] = tr[fmul(xi, x)] if xi and x else (tr[0] if not xi else tr[fmul(xi, x)])
            Trxi[xi, x] = tr[fmul(xi, x)]
    # χ̂(ξ)=∑_x χ(x) e_ξ(x); Ω={ξ: χ̂=p}  (up to the usual p vs -p)
    w = np.exp(2j * np.pi * np.arange(p) / p)
    chihat = np.zeros(q, dtype=np.complex128)
    for xi in range(1, q):
        chihat[xi] = (chi[1:] * w[Trxi[xi, 1:]]).sum()  # χ(0)=0
        # include 0: χ(0)=0
    # G(χ)=∑ χ(x) e(x) ; Ω typically |chihat|=p
    mag = np.abs(chihat)
    print(f"p={p} chihat unique {sorted(set(np.round(chihat[1:].real, 4)))[:8]}", flush=True)
    Omega = np.where(np.abs(chihat.real - p) < 0.5)[0]
    Omega = Omega[Omega > 0]
    Omegam = np.where(np.abs(chihat.real + p) < 0.5)[0]
    Omegam = Omegam[Omegam > 0]
    print(f"|Ω+|={len(Omega)} |Ω-|={len(Omegam)} expect {(q-1)//2}", flush=True)

    lab = MuLab(p, with_deg6=False)
    Z = lab.Yp.astype(np.int8)[:, 1:]  # (M,q)
    M = len(Z)
    # ẑ(ξ)=∑_x z_x exp(2πi Tr(ξx)/p)
    E = w[Trxi[Omega][:, :]]  # (|Ω|, q)
    # Z @ E.T : (M, |Ω|)
    zhat = Z.astype(np.complex128) @ E.T
    u = np.abs(zhat) ** 2  # (M, |Ω|)
    # support check: energy off Ω+
    Em = w[Trxi[Omegam]]
    zhm = Z.astype(np.complex128) @ Em.T
    print(f"mean |ẑ|² on Ω- (should ~0 if supported on Ω+): {np.mean(np.abs(zhm)**2):.4e}", flush=True)
    print(f"E[|ẑ|²] mean on Ω+={u.mean():.4f} expect 2q={2*q}", flush=True)
    flat = u.ravel()
    rnd = np.round(flat, 4)
    cnt = Counter(rnd.tolist())
    print("unique |ẑ|² values (top 12):", cnt.most_common(12), flush=True)
    print(f"min={flat.min():.4f} max={flat.max():.4f}", flush=True)
    # per-vector variance on Ω
    vrow = u.var(axis=1)
    print(f"per-y Var_Ω(|ẑ|²): mean={vrow.mean():.2f} min={vrow.min():.2f} max={vrow.max():.2f}", flush=True)
    # R̂_rest from λ=8+R̂/q² using known mins
    q2 = q * q
    known_min = {5: 80 / 13, 7: 3072 / 409, 11: 8.054447}
    if p in known_min:
        lam = known_min[p]
        Rhat = (lam - 8) * q2
        print(f"λ_min={lam:.6f}  R̂_rest={Rhat:.4f}  budget -2q²={-2*q2}", flush=True)
        print(f"margin R̂-(-2q²)={Rhat+2*q2:.4f}", flush=True)


if __name__ == "__main__":
    main()
