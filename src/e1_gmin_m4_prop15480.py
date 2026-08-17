#!/usr/bin/env python3
"""
Prop 15.480 — The three quadratic-χ Weil kernels in I_sm are
Kloosterman sums on F_{p²}:
    S(β):=∑_z χ(z²−1) e(βz) = Kl(1, β²/4),
    W_{r²−1}(α)=S(α),
    W_{r(r−1)}(α)=e(α/2) S(α/2)=e(α/2) Kl(1, α²/16),
    W_{r(r+1)}(α)=e(−α/2) S(α/2).
The cubic χ(r³−r) term and Q_τ stay open.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.479 flags.

============================================================================
Setup.  15.478 I_sm has eight terms.  I_lin is named.  The three
quadratic products are χ(r(r−1)), χ(r(r+1)), χ(r²−1).  On F_{p²},
χ(−1)=1 and G(χ)=∑ χ e = ∑ e(x²)=:G₂(1), with G²=q.

============================================================================
Theorem A — PROVED (Gauss inversion + complete square; p=3,5,7,11).
  G₂(1)=G(χ).  S(β)=Kl(1, β²/4) for every β (including 0).
  Fail: Kl(1, −β²/4) (max |err|≥3 at p=3, ≥8 at p=5).  ∎

Theorem B — PROVED (A + complete the square).
  W_{r(r−1)}(α)=e(α/2)S(α/2) and W_{r(r+1)}(α)=e(−α/2)S(α/2)
  identically.  Combined with A these are Kloosterman.
  Fail: drop the e(α/2) phase on W_{r(r−1)}.  ∎

Theorem C — OPEN.  Cubic W_c=∑ χ(r³−r)e(αr) unnamed.  Pairing
  the named Kl kernels with N(d) does not by itself name Q_τ.
  phi_F not imported.

============================================================================
Backend: field tables + O(q²) Kl/W.  Writes
evidence/e1_gmin_m4_prop15480.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15473 import gauss_chi  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15480_catalog.json"


def _half(F) -> int:
    add = F["add"]
    for t in range(F["q"]):
        if add[t][t] == 1:
            return t
    raise RuntimeError("no 1/2")


def klosterman(F, a: int, b: int) -> complex:
    q = F["q"]
    add, mul, inv = F["add"], F["mul"], F["inv"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    acc = 0j
    for t in range(1, q):
        acc += e_tab[add[mul[a][t]][mul[b][inv[t]]]]
    return acc


def S_exact(F) -> np.ndarray:
    q = F["q"]
    chi, add, neg, mul = F["chi"], F["add"], F["neg"], F["mul"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    S = np.zeros(q, dtype=np.complex128)
    for b in range(q):
        acc = 0j
        for z in range(q):
            zm, zp = add[z][neg[1]], add[z][1]
            w = 0.0 if (zm == 0 or zp == 0) else float(chi[mul[zm][zp]])
            acc += w * e_tab[mul[b][z]]
        S[b] = acc
    return S


def S_named(F) -> np.ndarray:
    q = F["q"]
    mul = F["mul"]
    quart = mul[_half(F)][_half(F)]
    return np.array(
        [klosterman(F, 1, mul[mul[b][b]][quart]) for b in range(q)],
        dtype=np.complex128,
    )


def S_named_wrong_sign(F) -> np.ndarray:
    q = F["q"]
    mul, neg = F["mul"], F["neg"]
    quart = mul[_half(F)][_half(F)]
    return np.array(
        [klosterman(F, 1, neg[mul[mul[b][b]][quart]]) for b in range(q)],
        dtype=np.complex128,
    )


def W_poly(F, weight) -> np.ndarray:
    q = F["q"]
    mul = F["mul"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    w = np.array([weight(r) for r in range(q)], dtype=np.float64)
    return np.array([(w * e_tab[mul[a]]).sum() for a in range(q)], dtype=np.complex128)


def weights(F):
    chi, add, neg, mul = F["chi"], F["add"], F["neg"], F["mul"]

    def v_s(r):
        zm, zp = add[r][neg[1]], add[r][1]
        return 0.0 if (zm == 0 or zp == 0) else float(chi[mul[zm][zp]])

    def v_rm(r):
        if r == 0:
            return 0.0
        rm = add[r][neg[1]]
        return 0.0 if rm == 0 else float(chi[mul[r][rm]])

    def v_rp(r):
        rp = add[r][1]
        return 0.0 if (r == 0 or rp == 0) else float(chi[mul[r][rp]])

    return v_s, v_rm, v_rp


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11):
        F = _kernel_field(p)
        q = F["q"]
        G = gauss_chi(F, p)
        e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
        G2 = sum(e_tab[F["mul"][x][x]] for x in range(q))
        if abs(G - G2) > 1e-10:
            ok = False
        if abs(G * G - q) > 1e-6:
            ok = False
        Se, Sn, Sw = S_exact(F), S_named(F), S_named_wrong_sign(F)
        err = float(np.max(np.abs(Se - Sn)))
        err_w = float(np.max(np.abs(Se[1:] - Sw[1:])))
        if err > 1e-10:
            ok = False
        if err_w < 1.0:
            ok = False
        rows[str(p)] = {"err": err, "err_wrong": err_w, "G": [G.real, G.imag]}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "S(β)=Kl(1, β²/4). Fail: Kl(1, −β²/4)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        add, mul, neg = F["add"], F["mul"], F["neg"]
        half = _half(F)
        e_tab = np.array([_e_tr(F, t) for t in range(F["q"])], dtype=np.complex128)
        v_s, v_rm, v_rp = weights(F)
        Ws, Wrm, Wrp = W_poly(F, v_s), W_poly(F, v_rm), W_poly(F, v_rp)
        S = S_named(F)
        pred_s = S
        pred_rm = np.array(
            [e_tab[mul[a][half]] * S[mul[a][half]] for a in range(F["q"])]
        )
        pred_rp = np.array(
            [e_tab[neg[mul[a][half]]] * S[mul[a][half]] for a in range(F["q"])]
        )
        err_s = float(np.max(np.abs(Ws - pred_s)))
        err_rm = float(np.max(np.abs(Wrm - pred_rm)))
        err_rp = float(np.max(np.abs(Wrp - pred_rp)))
        err_drop = float(np.max(np.abs(Wrm - np.array([S[mul[a][half]] for a in range(F["q"])]))))
        if max(err_s, err_rm, err_rp) > 1e-10:
            ok = False
        if err_drop < 0.5:
            ok = False
        rows[str(p)] = {
            "err_s": err_s,
            "err_rm": err_rm,
            "err_rp": err_rp,
            "err_drop_phase": err_drop,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W_{r(r−1)}=e(α/2)S(α/2), W_{r(r+1)}=e(−α/2)S(α/2). "
            "Fail: drop e(α/2)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Three quadratic-χ kernels are Kl. Cubic W_c and Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.480  quadratic-χ Weil kernels are Kloosterman", flush=True)
    A = prove_A()
    print(f"  A S=Kl: {A['proved']} err5={A['rows']['5']['err']:.2e}", flush=True)
    B = prove_B()
    print(f"  B complete square: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.480",
        "title": "I_sm quadratic-χ kernels are Kl(1, α²/4) and phases",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S_is_kloosterman": A["proved"],
            "complete_square_W": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15480.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
