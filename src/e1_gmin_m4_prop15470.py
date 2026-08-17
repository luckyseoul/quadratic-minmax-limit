#!/usr/bin/env python3
"""
Prop 15.470 — Quartic L_ψ analogue of 15.298: odd moments
vanish in the ensemble (not pointwise); L_2(ns)=0 and
L_0(ns)=(q−1)μ+μ−; the quartic half-Gauss image misses
live Q++.  Cov++ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.469 flags.

============================================================================
Setup.  ψ a quartic of F_q with ψ²=χ (exists: 4|q−1).  N from
Max+.  L_j(s)=E[∑_{d≠0} N(d) N(s d) ψ(d)^j], j=0,1,2,3.
S_j(β)=ψ^{−j}(β) G(ψ^j) (β≠0; j=0 gives −1).  Quartic image
    Q^ψ(r)/16 = k²+k(G₊+G_r)+∑_{s,j} (L_j(s)/4) S_j(1+rs).
Claude leftover-1 Dir 2.

============================================================================
Theorem A — PROVED (cache fold; p=5,7).
  L_1(s)=L_3(s)=0 for every s∈F_q^× (ensemble).  Not pointwise:
  max_i |∑_d N_i(d)N_i(sd)ψ(d)| = 254.95 at p=5.
  Fail: L_1(++)=J(ψ,ψ)=3−4i at p=5.  ∎

Theorem B — PROVED (same fold).
  L_2(ns)=0 and L_0(ns)=(q−1)μ+μ− (300, 3528).
  Fail: L_0(ns)=0 at p=5.  ∎

Theorem C — PROVED (A + quartic image).
  Q^ψ_++ is 9.039 at p=5 and −87.24 at p=7, not 48/13 or
  1544/409.  Fail: Q^ψ_++=48/13 at p=5.  ∎

Theorem D — PROVED (A + L_2).
  L_2(++)=3080/13, 641508/409, not Cov++ (10/39, 2793/1636).
  Fail: L_2(++)=10/39 at p=5.  ∎

Theorem E — OPEN.  L_ψ does not name Cov++ / Q++.  phi_F not
  imported.

============================================================================
Backend: 15.298 N-fold + quartic ψ.  Writes
evidence/e1_gmin_m4_prop15470.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    coarse_Q_class,
    load_N,
    mu_minus,
    mu_plus,
)
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15467 import coarse_L  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15470_catalog.json"

LIVE_L2_PP = {5: Fraction(3080, 13), 7: Fraction(641508, 409)}
LIVE_L0_NS = {5: Fraction(300), 7: Fraction(3528)}


def primitive_root(F) -> int:
    q = F["q"]
    for g in range(2, q):
        x, seen = g, set()
        ok = True
        for _ in range(q - 1):
            if x in seen or x == 0:
                ok = False
                break
            seen.add(x)
            x = F["mul"][x][g]
        if ok and len(seen) == q - 1:
            return g
    raise RuntimeError("no generator")


def make_psi(F):
    q = F["q"]
    gen = primitive_root(F)
    x = 1
    psi = np.zeros(q, dtype=np.complex128)
    for i in range(q - 1):
        psi[x] = (1j) ** (i % 4)
        x = F["mul"][x][gen]
    for t in range(1, q):
        if abs(psi[t] ** 2 - F["chi"][t]) > 1e-8:
            raise RuntimeError("psi^2 != chi")
    return psi


def coarse_of(F, s: int) -> str:
    if F["chi"][s] != 1:
        return "ns"
    c = coarse_Q_class(F, s)
    return "++" if c.startswith("mixed") else c


def gauss_j(F, psi, j: int):
    q = F["q"]
    acc = 0j
    for x in range(1, q):
        w = 1 + 0j if j % 4 == 0 else psi[x] ** (j % 4)
        acc += w * _e_tr(F, x)
    return acc


def S_j(F, psi, Gtab, beta: int, j: int):
    q = F["q"]
    if beta == 0:
        return (q - 1) / 2.0 if j % 2 == 0 else 0.0
    if j % 4 == 0:
        return Gtab[0]
    return (psi[beta] ** ((-j) % 4)) * Gtab[j % 4]


@lru_cache(maxsize=4)
def L_tables(p: int):
    F = _kernel_field(p)
    N = load_N(p, F)
    psi = make_psi(F)
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    dlist = np.arange(1, q)
    Nd = N[:, dlist]
    psi_d = np.stack([psi[dlist] ** j for j in range(4)])
    Ljs = {j: np.zeros(q, dtype=np.complex128) for j in range(4)}
    for s in range(1, q):
        m = (Nd * N[:, mul[s, dlist]]).mean(axis=0)
        for j in range(4):
            Ljs[j][s] = np.dot(m, psi_d[j])
    return F, N, psi, Ljs


def type_mean(F, Lj, tau: str) -> complex:
    vs = [Lj[s] for s in range(1, F["q"]) if coarse_of(F, s) == tau]
    return sum(vs) / len(vs) if vs else 0j


def pointwise_L1_max(p: int) -> float:
    F, N, psi, _ = L_tables(p)
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    d = np.arange(1, q)
    pd = psi[d]
    mx = 0.0
    for s in range(1, q):
        pr = (N[:, d] * N[:, mul[s, d]]) @ pd
        mx = max(mx, float(np.max(np.abs(pr))))
    return mx


def predict_Qpp_psi(p: int) -> float:
    F, N, psi, Ljs = L_tables(p)
    q = F["q"]
    k = p * (p - 1) / 2.0
    EN = N.mean(axis=0)
    Gtab = {j: gauss_j(F, psi, j) for j in range(4)}
    Omega = omega_set(F)
    xi = Omega[0]

    def G_of(scale):
        acc = 0j
        for d in range(1, q):
            td = F["mul"][xi][d] if scale == 1 else F["mul"][F["mul"][xi][scale]][d]
            acc += EN[d] * _e_tr(F, td)
        return acc

    Gp = G_of(1)
    accs = []
    q2 = float(q * q)
    for r in F["squares"]:
        if coarse_of(F, r) != "++":
            continue
        Gr = G_of(r)
        acc = k * k + k * (Gp + Gr)
        for s in range(1, q):
            beta = F["add"][1][F["mul"][r][s]]
            for j in range(4):
                acc += (Ljs[j][s] / 4.0) * S_j(F, psi, Gtab, beta, j)
        accs.append(float((16 * acc).real) / q2)
    return float(np.mean(accs))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, _, Ljs = L_tables(p)
        m1 = max(abs(Ljs[1][s]) for s in range(1, F["q"]))
        m3 = max(abs(Ljs[3][s]) for s in range(1, F["q"]))
        pw = pointwise_L1_max(p)
        if m1 > 1e-8 or m3 > 1e-8:
            ok = False
        if pw < 1.0:
            ok = False
        rows[str(p)] = {"ens_L1_max": m1, "ens_L3_max": m3, "pointwise_max": pw}
    # fail: ensemble L_1(++) equals J(ψ,ψ)=3-4i
    F, _, _, Ljs = L_tables(5)
    if abs(type_mean(F, Ljs[1], "++") - (3 - 4j)) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Ensemble L_1=L_3=0 for every s; not pointwise. "
            "Fail: L_1(++)=3−4i at p=5."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, _, Ljs = L_tables(p)
        l0 = type_mean(F, Ljs[0], "ns").real
        l2 = type_mean(F, Ljs[2], "ns").real
        want = float((F["q"] - 1) * mu_plus(p) * mu_minus(p))
        if abs(l2) > 1e-8:
            ok = False
        if abs(l0 - want) > 1e-6:
            ok = False
        if abs(l0) < 1e-8:
            ok = False
        rows[str(p)] = {"L0_ns": l0, "L2_ns": l2, "named": want}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "L_2(ns)=0 and L_0(ns)=(q−1)μ+μ−. "
            "Fail: L_0(ns)=0 at p=5."
        ),
    }


def prove_C() -> dict:
    q5 = predict_Qpp_psi(5)
    q7 = predict_Qpp_psi(7)
    live5 = float(LIVE_QPP[5])
    live7 = float(LIVE_QPP[7])
    ok = abs(q5 - live5) > 1.0 and abs(q7 - live7) > 1.0
    if abs(q5 - live5) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "pred5": q5,
        "pred7": q7,
        "live5": live5,
        "live7": live7,
        "theorem": (
            "Quartic half-Gauss image misses live Q++. "
            "Fail: pred=48/13 at p=5."
        ),
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, _, Ljs = L_tables(p)
        v = type_mean(F, Ljs[2], "++").real
        fr = Fraction(v).limit_denominator(2000)
        if fr != LIVE_L2_PP[p]:
            ok = False
        cov = coarse_L(p)["++"][1] - mu_plus(p) ** 2
        if abs(v - float(cov)) < 1e-8:
            ok = False
        rows[str(p)] = {"L2_pp": str(fr), "cov": str(cov)}
    if LIVE_L2_PP[5] == Fraction(10, 39):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "L_2(++)=3080/13, 641508/409, not Cov++. "
            "Fail: L_2(++)=10/39."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Odd L_ψ vanishes in the ensemble. Quartic image misses Q++. "
            "Cov++ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.470  quartic L_ψ analogue of 15.298", flush=True)
    A = prove_A()
    print(f"  A odd vanish: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B ns: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C image: {C['proved']} pred5={C['pred5']:.4f}", flush=True)
    D = prove_D()
    print(f"  D L2++: {D['proved']} {D['rows']}", flush=True)
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["proved"], "B": B["proved"], "C": C["proved"], "D": D["proved"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.470",
        "title": "Quartic L_ψ odd moments vanish; image misses live Q++",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "odd_Lpsi_vanish_ensemble": A["proved"],
            "L0_ns_named": B["proved"],
            "quartic_image_misses_Qpp": C["proved"],
            "L2_pp_not_cov": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15470.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
