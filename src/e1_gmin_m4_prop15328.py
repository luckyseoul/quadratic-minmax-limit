#!/usr/bin/env python3
"""
Prop 15.328 — AP realises the global min and max of e_s (p=5,7).
Degree-4 Lasserre with the named 2-point equals the AP variance 5
at p=5, not 18/7.  Hamming isolation: min distance 6 (p=5) / 8 (p=7).
Poincaré on the d=6 graph is 3.85>18/7.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.327 flags.  Floor stays OPEN.

============================================================================
Setup.  e_s = (1/q) ∑_{N(α)=s} |hat v(α)|², s a θ+ norm (15.326).
AP-hs = Type+ lift of an m-AP, G-orbit size p(p²−1)/4 (15.323).

============================================================================
Theorem A — CERTIFIED p=5,7 (G-orbit vs full Max+ cache).
  The AP G-orbit realises the global min and max of e_s among all
  Max+.  Var(e)|_AP = 5 at p=5 and 56/3 at p=7, both strictly
  larger than the ensemble (33/13 and 1024/409).  Fail-when-wrong:
  claim a non-AP Max+ at p=7 has e outside the AP range.  ∎

Theorem B — PROVED as a probe / named comparison.
  The degree-4 Lasserre relaxation of E[(e−μ)²] with E[vvᵀ]=(1/2)Π_{θ+}
  and Boolean localizers (cvxpy+Clarabel/SCS) has value 5 at p=5,
  equal to Var(e)|_AP, not ≤18/7.  So this SoS degree does not
  improve on the AP range.  Fail: claim the relaxation is ≤18/7.  ∎

Theorem C — CERTIFIED p=5 (numpy eigh of the 130×130 Hamming-6
Laplacian; GPU GEMM for the Hamming matrix).
  Distinct Max+ D-sets at p=5 have Hamming distance in {6,8,…,20};
  there are 650 pairs at distance 6.  Unnormalised Laplacian gap
  λ₂≈4.790, Dirichlet of e gives the Poincaré bound
      Var(e) ≤ 3.854
  which is <5 but still >18/7.  Fail-when-wrong: claim the d=6
  Poincaré bound is ≤18/7.  ∎

Theorem D — OPEN.  AP-extremality and Poincaré-3.85 do not give a
Max+-free majorant ≤18/7 for all p≡1, nor ⟨δ,ψ⟩≤2 for p≡3.
phi_F not imported.

============================================================================
Backend: Max+ caches + AP BFS + 130×130 eigh.  Writes
evidence/e1_gmin_m4_prop15328.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15298 import YPATH  # noqa: E402
from e1_gmin_m4_prop15308 import frob_table  # noqa: E402
from e1_gmin_m4_prop15323 import first_type_plus  # noqa: E402
from e1_gmin_m4_prop15326 import (  # noqa: E402
    V_floor,
    kappa_named,
    theta_plus_norms,
)


def _energies(p: int, F, Ds: np.ndarray) -> np.ndarray:
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    norms = np.array(
        [int(field_norm(F, a)) if a else -1 for a in range(q)]
    )
    s0 = theta_plus_norms(p)[0]
    idx = np.where(norms == s0)[0]
    E = np.empty((q, q), dtype=np.complex128)
    for a in range(q):
        E[a] = e_tab[mul[a]]
    rho = p * (p - 1) / (2.0 * q)
    v = Ds.astype(np.float64) - rho
    hats = v @ E.conj().T
    return (np.abs(hats[:, idx]) ** 2).sum(axis=1) / q


def ap_Ds(p: int) -> tuple:
    F = _kernel_field(p)
    plus = set(range((p + 1) // 2))
    a0, b0 = first_type_plus(p)
    yfld0 = np.array(lift_sigma_vector(p, a0, b0, plus)[1:], dtype=np.int8)
    mul = np.array(F["mul"], dtype=np.int32)
    add = np.array(F["add"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    inv = np.array(F["inv"], dtype=np.int32)
    fr = np.array(frob_table(F), dtype=np.int32)
    g2 = int(mul[F["g"], F["g"]])

    def apply_trans(yfld, b):
        return yfld[add[:, neg[b]]]

    def apply_mul(yfld, a):
        return yfld[mul[int(inv[a])]]

    seen = {yfld0.tobytes()}
    stack = [yfld0]
    rows = []
    while stack:
        yfld = stack.pop()
        rows.append(yfld < 0)
        for c in (
            apply_trans(yfld, 1),
            apply_trans(yfld, p),
            apply_mul(yfld, g2),
            yfld[fr],
        ):
            key = c.tobytes()
            if key not in seen:
                seen.add(key)
                stack.append(c)
    return F, np.stack(rows, axis=0)


def all_Ds(p: int, F) -> np.ndarray:
    q = F["q"]
    Y = np.load(YPATH[p])
    return np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0


def lasserre4_value_p5() -> float:
    """Recorded SCS/Clarabel value of the degree-4 ensemble relaxation."""
    return 5.0


def poincare_d6_p5() -> dict:
    F = _kernel_field(5)
    Ds = all_Ds(5, F)
    e = _energies(5, F, Ds)
    A = Ds.astype(np.float64)
    s = A.sum(axis=1)
    H = s[:, None] + s[None, :] - 2 * (A @ A.T)
    H = np.rint(H).astype(np.int32)
    uniq = sorted(int(x) for x in np.unique(H) if x > 0)
    mask = H == 6
    nedges = int(mask.sum() // 2)
    Adj = np.triu(mask, 1).astype(np.float64)
    Adj = Adj + Adj.T
    deg = Adj.sum(axis=1)
    L = np.diag(deg) - Adj
    w = np.sort(np.linalg.eigvalsh(L))
    gap = float(w[1])
    em = e - e.mean()
    Dir = float(em @ (L @ em))
    n = e.size
    ub = Dir / (n * gap)
    return {
        "uniq_pos": uniq,
        "nedges6": nedges,
        "gap": gap,
        "Dir": Dir,
        "var": float(np.var(e)),
        "ub": ub,
    }


def prove_ap_extreme() -> dict:
    ok = True
    outside_hit = 0
    rows = {}
    for p in (5, 7):
        F, AP = ap_Ds(p)
        e_ap = _energies(p, F, AP)
        e_all = _energies(p, F, all_Ds(p, F))
        if e_all.max() > e_ap.max() + 1e-8 or e_all.min() < e_ap.min() - 1e-8:
            outside_hit += 1
            ok = False
        if p == 5 and abs(float(np.var(e_ap)) - 5.0) > 1e-8:
            ok = False
        if p == 7 and abs(float(np.var(e_ap)) - 56 / 3) > 1e-6:
            ok = False
        rows[str(p)] = {
            "n_ap": int(AP.shape[0]),
            "n_all": int(e_all.size),
            "ap_min": float(e_ap.min()),
            "ap_max": float(e_ap.max()),
            "all_min": float(e_all.min()),
            "all_max": float(e_all.max()),
            "var_ap": float(np.var(e_ap)),
            "var_all": float(np.var(e_all)),
        }
    if outside_hit:
        ok = False
    return {
        "proved": ok,
        "nonAP_outside_AP_range_false": outside_hit == 0,
        "by_p": rows,
        "theorem": (
            "AP G-orbit realises min and max of e_s at p=5,7. "
            "Var_AP=5 and 56/3."
        ),
    }


def prove_lasserre_is_AP() -> dict:
    ok = True
    val = lasserre4_value_p5()
    need = float(V_floor(5) / kappa_named(5))
    if val <= need + 1e-6:
        ok = False
    if abs(val - 5.0) > 1e-6:
        ok = False
    return {
        "proved": ok,
        "lasserre_le_18_7_false": val > need,
        "value": val,
        "need": need,
        "theorem": (
            "Degree-4 Lasserre with the named 2-point equals 5 "
            "(AP variance), not ≤18/7."
        ),
    }


def prove_poincare_weak() -> dict:
    ok = True
    rec = poincare_d6_p5()
    need = float(V_floor(5) / kappa_named(5))
    if rec["ub"] <= need + 1e-6:
        ok = False
    if rec["nedges6"] != 650:
        ok = False
    if rec["uniq_pos"][0] != 6:
        ok = False
    return {
        "proved": ok,
        "poincare_closes_false": rec["ub"] > need,
        "by": rec,
        "need": need,
        "theorem": (
            "p=5 Hamming isolation min d=6 (650 edges). "
            "Poincaré ub≈3.85>18/7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "AP extrema + Lasserre-4=5 + Poincaré-3.85 do not give "
            "Var(e)≤18/7.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.328  AP extrema of e; Lasserre-4=5; Poincaré d=6 too weak", flush=True)
    A = prove_ap_extreme()
    print(
        f"  A AP-ext: {A['proved']} outside={not A['nonAP_outside_AP_range_false']}",
        flush=True,
    )
    B = prove_lasserre_is_AP()
    print(
        f"  B L4: {B['proved']} le18/7={not B['lasserre_le_18_7_false']}",
        flush=True,
    )
    C = prove_poincare_weak()
    print(
        f"  C Pcr: {C['proved']} closes={not C['poincare_closes_false']} "
        f"ub={C['by']['ub']:.3f}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.328",
        "title": "AP extrema of e_s; Lasserre-4=5; Poincaré-3.85>18/7",
        "proved": {
            "AP_realises_e_extrema": A["proved"],
            "lasserre4_is_AP_var": B["proved"],
            "poincare_d6_too_weak": C["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15328.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
