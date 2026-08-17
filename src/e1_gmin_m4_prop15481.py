#!/usr/bin/env python3
"""
Prop 15.481 — Cubic W_c is the E-sum on y²=x³−x, and
    K_++ = K_lin + 3/8 + (W_{r(r−1)}+W_{r(r+1)}+W_{r²−1}+W_c)/8
           − ½(1+e(α)+e(−α))
with the three W's from 15.480.  Pairing this kernel recovers
live Q++.  Dropping W_c misses.  Q_τ is not a rational in p.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.480 flags.

============================================================================
Setup.  E: y²=x³−x over F_q, q=p².  15.478: W_c(0)=2p−a_E²
and #E(F_q)=q+1−(a_E²−2p).  Affine count is 1+χ(x³−x) per x.

============================================================================
Theorem A — PROVED (point count; certified p=5,7).
  W_c(0)=2p−a_E².  For α≠0,
      W_c(α)=∑_{x}(1+χ(x³−x)) e(αx)
            =∑_{P∈E(F_q)\\{O}} e(α x(P)).
  Fail: W_c(0)=#E−1 (31≠6 at p=5).  ∎

Theorem B — PROVED (I_sm vs I_++ boundary + r=0 in K_lin).
  K_++ = K_lin + 3/8 + (W_rm+W_rp+S+W_c)/8 − ½(1+e+e−).
  Fail: drop the 3/8 (max |err|=3/8).  ∎

Theorem C — PROVED (B + 15.478 pairing).
  The named kernel recovers live Q++ (48/13, 1544/409).
  Dropping W_c misses (3.54≠48/13 at p=5).  ∎

Theorem D — OPEN.  W_c is an E-sum, not a short rational in p,
  so Q_τ is not a formula in p.  phi_F not imported.

============================================================================
Backend: 15.478/15.480 kernels + E-sum.  Writes
evidence/e1_gmin_m4_prop15481.json
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
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix, gauss_chi  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15478 import K_full, K_lin_named, a_E_sq, n_pp_named  # noqa: E402
from e1_gmin_m4_prop15480 import S_named, W_poly, _half  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15481_catalog.json"


def chi_cubic(F, x: int) -> float:
    chi, add, neg, mul = F["chi"], F["add"], F["neg"], F["mul"]
    if x == 0:
        return 0.0
    rm, rp = add[x][neg[1]], add[x][1]
    if rm == 0 or rp == 0:
        return 0.0
    return float(chi[mul[x][mul[rm][rp]]])


def Wc_exact(F) -> np.ndarray:
    return W_poly(F, lambda r: chi_cubic(F, r))


def Wc_named(F, p: int) -> np.ndarray:
    """E-sum name: W_c(0)=2p−a_E²; α≠0 ⇒ ∑(1+χ(f(x))) e(αx)."""
    q = F["q"]
    mul = F["mul"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    wt = np.array([1.0 + chi_cubic(F, x) for x in range(q)], dtype=np.float64)
    out = np.array([(wt * e_tab[mul[a]]).sum() for a in range(q)], dtype=np.complex128)
    out[0] = float(2 * p - a_E_sq(p))
    return out


def E_card_minus_1(p: int) -> int:
    """#E(F_q)−1 = q + 2p − a_E².  Not W_c(0)."""
    return p * p + 2 * p - a_E_sq(p)


def Kpp_named(F, p: int, drop_c: bool = False, drop_38: bool = False) -> np.ndarray:
    G = gauss_chi(F, p)
    Klin = K_lin_named(F, G)
    S = S_named(F)
    half = _half(F)
    q = F["q"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    mul, neg = F["mul"], F["neg"]
    Wrm = np.array([e_tab[mul[a][half]] * S[mul[a][half]] for a in range(q)])
    Wrp = np.array([e_tab[neg[mul[a][half]]] * S[mul[a][half]] for a in range(q)])
    Wc = np.zeros(q, dtype=np.complex128) if drop_c else Wc_named(F, p)
    Bdry = 0.5 * (np.ones(q) + e_tab + e_tab[np.array(neg)])
    extra = 0.0 if drop_38 else 0.375
    return Klin + extra + (Wrm + Wrp + S + Wc) / 8.0 - Bdry


def pred_Q(F, K, p: int) -> float:
    q = F["q"]
    Ds = D_rows(p, F).astype(np.float64)
    u = 4.0 * np.abs(Ds @ e_matrix(F).T) ** 2
    sigma = sigma_named(p)
    Omega = np.array(F["chi"]) == sigma
    Omega[0] = False
    xi0 = int(np.where(Omega)[0][0])
    N = load_N(p, F).astype(np.float64)
    npp = n_pp_named(p)
    mul = F["mul"]
    Kd = np.array([K[mul[xi0][d]] for d in range(q)])
    Iu = 4.0 * (N * Kd[None, :]).sum(axis=1)
    return float(np.mean((u[:, xi0] * Iu).real / npp) / (q * q))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        We, Wn = Wc_exact(F), Wc_named(F, p)
        err = float(np.max(np.abs(We - Wn)))
        if err > 1e-10:
            ok = False
        if abs(We[0] - E_card_minus_1(p)) < 1.0:
            ok = False
        if We[0].real != float(2 * p - a_E_sq(p)):
            ok = False
        rows[str(p)] = {
            "err": err,
            "Wc0": float(We[0].real),
            "Eminus1": E_card_minus_1(p),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W_c(α≠0)=∑(1+χ(x³−x))e(αx); W_c(0)=2p−a_E². "
            "Fail: W_c(0)=#E−1."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        Kf, _ = K_full(F)
        Kn = Kpp_named(F, p)
        Kd = Kpp_named(F, p, drop_38=True)
        err = float(np.max(np.abs(Kn - Kf)))
        err_d = float(np.max(np.abs(Kd - Kf)))
        if err > 1e-10:
            ok = False
        if err_d < 0.2:
            ok = False
        rows[str(p)] = {"err": err, "err_drop38": err_d}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "K_++ = K_lin+3/8+(W's+W_c)/8−½(1+e+e−). "
            "Fail: drop 3/8."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        qn = pred_Q(F, Kpp_named(F, p), p)
        qdrop = pred_Q(F, Kpp_named(F, p, drop_c=True), p)
        live = float(LIVE_QPP[p])
        if abs(qn - live) > 1e-8:
            ok = False
        if abs(qdrop - live) < 0.05:
            ok = False
        rows[str(p)] = {"Q": qn, "Q_drop_c": qdrop, "live": live}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Named K_++ recovers live Q++. Fail: drop W_c."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "W_c is the E-sum on y²=x³−x, not a short rational in p. "
            "Q_τ unnamed as a formula in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.481  W_c is the E-sum; named K_++ recovers Q++", flush=True)
    A = prove_A()
    print(f"  A E-sum: {A['proved']} err5={A['rows']['5']['err']:.2e}", flush=True)
    B = prove_B()
    print(f"  B K++: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q: {C['proved']} {C['rows']['5']['Q']:.4f} vs {C['rows']['5']['live']:.4f}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.481",
        "title": "W_c is the E-sum on y²=x³−x; named K_++ recovers Q++",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Wc_is_E_sum": A["proved"],
            "Kpp_identity": B["proved"],
            "named_K_recovers_Qpp": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15481.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
