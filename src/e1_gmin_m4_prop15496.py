#!/usr/bin/env python3
"""
Prop 15.496 — Fourier of the 3-linear A_r is affine in Re J(χ,ψ)
and the 15.302 Hasse–Davenport weight σ.  A4 / 4-distinct F_τ
are not in that span.  H-levels split A_r on ++ sub, not F_τ.

    A_r = E[1̂_D(ξ) 1̂_D(rξ) 1̂_D(−(1+r)ξ)]   (15.494)
    Ã(ψ) := ∑_{r∈T\\{±1}} ψ(r) A_r
    σ_k = (p−3)/2  if (p−1)|k,  else −1          (15.302)

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.495 flags.  Does **not** name A4 or Q_τ / F_τ.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import (
    H_legendre,
    _e_tr,
    _gauss,
    _kernel_field,
    _psi_vec,
)
from e1_gmin_m4_prop15290 import live_Q_on_T, omega_set, type_key
from e1_gmin_m4_prop15301 import jacobi_chi_psi, jacobi_inverted
from e1_gmin_m4_prop15494 import minus_one_plus_r

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15496.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
_AR_CACHE: dict = {}


def even_ks(q: int) -> list[int]:
    return [k for k in range(2, q - 1, 2) if k != (q - 1) // 2]


def sigma_hd(p: int, k: int) -> float:
    return (p - 3) / 2 if k % (p - 1) == 0 else -1.0


def coll_off(p: int) -> int:
    q = p * p
    return 3 * q * q - 10 * q


def a_r_table(p: int) -> dict:
    """Per-r A_r, H(r), type.  One Max+ fold."""
    if p in _AR_CACHE:
        return _AR_CACHE[p]
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    D = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0
    Omega = omega_set(F)
    xi = Omega[0]
    Toff = [r for r in F["squares"] if r not in (1, F["neg1"])]
    need = {xi}
    for r in Toff:
        need.add(F["mul"][r][xi])
        need.add(F["mul"][minus_one_plus_r(F, r)][xi])
    hats = {}
    mul = F["mul"]
    for a in need:
        e = np.array([_e_tr(F, mul[a][x]) for x in range(q)])
        hats[int(a)] = (D.astype(np.float64) * e).sum(axis=1)
    out = {}
    for r in Toff:
        rxi = mul[r][xi]
        sxi = mul[minus_one_plus_r(F, r)][xi]
        Ar = float(np.mean(hats[xi] * hats[rxi] * hats[sxi]).real)
        out[int(r)] = {
            "A_r": Ar,
            "H": int(H_legendre(F, r)),
            "type": type_key(F, r),
        }
    _AR_CACHE[p] = out
    return out


def Ahat_of(F, k: int, table: dict) -> complex:
    psi = _psi_vec(F, k)
    return sum(psi[r] * rec["A_r"] for r, rec in table.items())


def slack_of(F, k: int, Q: dict) -> float:
    q2 = float(F["q"] ** 2)
    Toff = [r for r in F["squares"] if r not in (1, F["neg1"])]
    psi = _psi_vec(F, k)
    sdelta = sum((4.0 - Q[r] / q2) * psi[r].real for r in Toff)
    return float(2.0 - sdelta)


def _lstsq_err(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    err = float(np.max(np.abs(X @ coef - y)))
    return coef, err


def prove_A() -> dict:
    """Ã is real and Ã = a + b Re J(χ,ψ) + c σ at p=5,7."""
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        tab = a_r_table(p)
        xs, ys, ims = [], [], []
        for k in even_ks(F["q"]):
            ah = Ahat_of(F, k, tab)
            J = jacobi_chi_psi(F, k)
            xs.append([1.0, float(J.real), sigma_hd(p, k)])
            ys.append(float(ah.real))
            ims.append(abs(ah.imag))
        X, y = np.asarray(xs), np.asarray(ys)
        coef, err = _lstsq_err(X, y)
        imax = float(max(ims))
        if err > 1e-8 or imax > 1e-8:
            ok = False
        rows[p] = {
            "coef": [float(c) for c in coef],
            "err": err,
            "im_max": imax,
            "n_k": int(len(ys)),
        }
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """Fail-when-wrong: drop σ; invert J; put slack in the same span."""
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        tab = a_r_table(p)
        Q = live_Q_on_T(p)
        xs, xdrop, xinv, yA, yS = [], [], [], [], []
        for k in even_ks(F["q"]):
            ah = Ahat_of(F, k, tab)
            J = jacobi_chi_psi(F, k)
            Jinv = jacobi_inverted(F, k)
            sig = sigma_hd(p, k)
            xs.append([1.0, float(J.real), sig])
            xdrop.append([1.0, float(J.real)])
            xinv.append([1.0, float(Jinv.real), sig])
            yA.append(float(ah.real))
            yS.append(slack_of(F, k, Q))
        X, Xd, Xi = map(np.asarray, (xs, xdrop, xinv))
        yA, yS = map(np.asarray, (yA, yS))
        _, e_full = _lstsq_err(X, yA)
        _, e_drop = _lstsq_err(Xd, yA)
        _, e_inv = _lstsq_err(Xi, yA)
        _, e_slack = _lstsq_err(X, yS)
        rec = {
            "Ahat_full": e_full,
            "Ahat_drop_sig": e_drop,
            "Ahat_Jinv": e_inv,
            "slack_full": e_slack,
        }
        rows[p] = rec
        # p=7: drop-σ and slack must fail; full Ahat stays tight
        if p == 7:
            if e_full > 1e-8 or e_drop < 0.5 or e_slack < 0.5:
                ok = False
        if p == 5:
            # p=5 is 3-class interpolant; slack also sits in the span (15.301)
            if e_full > 1e-8 or e_slack > 1e-8:
                ok = False
    # inverted J must miss at least one prime
    if rows[5]["Ahat_Jinv"] < 1e-6 and rows[7]["Ahat_Jinv"] < 1e-6:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """H-levels split A_r on ++ sub; F_τ is constant while H spreads."""
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        tab = a_r_table(p)
        Q = live_Q_on_T(p)
        q2 = float(F["q"] ** 2)
        Coll = coll_off(p)
        by = defaultdict(lambda: {"A": [], "H": [], "F": []})
        for r, rec in tab.items():
            by[str(rec["type"])]["A"].append(rec["A_r"])
            by[str(rec["type"])]["H"].append(rec["H"])
            by[str(rec["type"])]["F"].append(Q[r] - Coll)
        packed = {}
        for t, vs in by.items():
            A = np.asarray(vs["A"], float)
            H = np.asarray(vs["H"], float)
            Ff = np.asarray(vs["F"], float)
            packed[t] = {
                "n": int(len(A)),
                "A_spread": float(A.max() - A.min()),
                "H_spread": float(H.max() - H.min()),
                "F_spread": float(Ff.max() - Ff.min()),
                "H_vals": [int(h) for h in H],
                "A_by_H": {
                    str(int(h)): float(A[H == h].mean())
                    for h in sorted(set(H.astype(int)))
                },
            }
        # F type-constant
        if any(v["F_spread"] > 1e-6 for v in packed.values()):
            ok = False
        sub = packed.get("(1, 1, 'sub')")
        if sub is None:
            ok = False
        else:
            byH = defaultdict(list)
            for r, rec in tab.items():
                if rec["type"] == (1, 1, "sub"):
                    byH[rec["H"]].append(rec["A_r"])
            h_level_ok = all(max(vs) - min(vs) < 1e-8 for vs in byH.values())
            if not h_level_ok:
                ok = False
            packed["sub_H_levels"] = {
                str(h): {"n": len(vs), "spread": float(max(vs) - min(vs))}
                for h, vs in byH.items()
            }
            if p == 7:
                # H spreads, F does not: H is not a name of F_τ
                if sub["H_spread"] < 1.0 or sub["F_spread"] > 1e-6:
                    ok = False
                if len(byH) < 2:
                    ok = False
        rows[p] = packed
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_imported": False,
        "Q_tau_named": False,
        "A4_named": False,
        "F_tau_named": False,
        "note": (
            "Ã∈span{1, Re J(χ,ψ), σ} at p=5,7 names the Fourier of A_r, "
            "not A4=(slack) q G(ψ̄) and not the 4-distinct F=Q−Coll on "
            "15.290 types.  Coefs have den |H+|/(2p).  H-levels split "
            "A_r on ++ sub and do not split F.  ⟨δ,ψ⟩≤2 stays OPEN."
        ),
    }


def main() -> dict:
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
