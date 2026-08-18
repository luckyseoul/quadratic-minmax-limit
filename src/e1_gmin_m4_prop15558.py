#!/usr/bin/env python3
"""
Prop 15.558 — J_all and J_T are named complete sums; Φ_free still
carries the unnamed Q3_T.
    J_all(v,w)=(1/8) Σ_{ε∈{0,1}³} χ_Ω^{|ε|} S_ε(v,w)
with 7 Gauss S_ε and S_{111}=G(χ,1) C(v,w),
    C=Σ_{t≠0,−1,−v/w} χ(t(t+1)(v+t w)).
    J_1line=Σ_{c∈F_p^×\\{-1}} S_Ω((1+c)v−c u),
    J_double=Σ_{r∈{1,−2,−1/2}} S_Ω((1+r)v−r u),
    J_n3=J_all−J_1line.
G3=G02+q^{-3} Σ_T Q3_T J_T.  Q3_T is not a p-law.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.557 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.557.

============================================================================
Setup.  1_Ω=(1+χ_Ω χ)/2.  S_Ω(δ)=|Ω| (δ=0) else −1/2+(χ_Ω/2)G(χ,δ).
G02=(N/p³)[1+2S_Ω(u)+2S_Ω(v)+2S_Ω(v−u)].  15.556 G3_1d.  15.554 Q3_T.

============================================================================
Theorem A — PROVED (Gauss inversion; p=5,7,11,13).
  The seven non-cubic S_ε are Gauss.  S_{111}=G(χ,1)C (w≠0;
  S_{111}(v,0)=−G(χ,v)).  J_all=(1/8)Σ χ_Ω^{|ε|} S_ε.
  Fail: drop S_{111}.  Fail: drop the 1/8.  ∎

Theorem B — PROVED (A + F_p-lines; p=5,7,11,13).
  η∈Ω ⇒ cη∈Ω for c∈F_p^×, so J_1line and J_double are S_Ω
  complete sums as above.  J_generic=J_1line−J_double,
  J_n3=J_all−J_1line.  Fail: J_1line=J_all.  Fail: J_double=J_1line.  ∎

Theorem C — PROVED (B + 15.554/555; p=5,7).
  G3(0,u,v)=G02+q^{-3} Σ_T Q3_T J_T on 3-distinct pairs.
  Hence Φ_free=G02−G3_1d+q^{-3} Σ_T Q3_T J_T.
  Fail: G02=Φ_free (p=5 maxerr 24).  Fail: one Q3 for all T.  ∎

Theorem D — OPEN.  Q3_T takes 2 (p=5) / 3 (p=7) values with
  Max+ den 13,409, not a Gauss/Jacobi integer in p.  Pointwise
  Φ_free / 16pA / Q_τ unnamed.  Do not interpolate (p−5)/15.

============================================================================
Backend: field Gauss O(q²); p=5,7 G3 fold.  Writes
evidence/e1_gmin_m4_prop15558.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15550 import _e_tab, omega_set_local
from e1_gmin_m4_prop15554 import LIVE_Q3, load_Y, triple_kind
from e1_gmin_m4_prop15556 import g3_1d_named, is_collinear

EV = ROOT / "evidence" / "e1_gmin_m4_prop15558.json"
LIVE_FOUR_G02 = {5: -5304.0, 7: -323928.0}
LIVE_N = {5: 130, 7: 5726}


def q302_named(p: int, N: int) -> float:
    """Q3_02 = −4N(2p²+1)/p from F=−2(3p²+2) (certified p=5..23)."""
    return -4.0 * N * (2 * p * p + 1) / p


def q302_drop_plus1(p: int, N: int) -> float:
    """Fail: drop the +1."""
    return -4.0 * N * (2 * p * p) / p


def gchi_table(F) -> np.ndarray:
    q = F["q"]
    G = np.zeros(q, dtype=np.complex128)
    for a in range(q):
        G[a] = sum(F["chi"][x] * _e_tr(F, F["mul"][a][x]) for x in range(1, q))
    return G


def S_Omega_named(F, chi_om: int, G, d: int) -> complex:
    if d == 0:
        return (F["q"] - 1) / 2.0
    return -0.5 + 0.5 * chi_om * G[d]


def S_eps_named(F, G, e1: int, e2: int, e3: int, v: int, w: int) -> complex:
    q = F["q"]
    add, neg, chi, mul, inv = F["add"], F["neg"], F["chi"], F["mul"], F["inv"]
    chi_m1 = chi[neg[1]]
    dv0 = 1.0 if v == 0 else 0.0
    dw0 = 1.0 if w == 0 else 0.0
    dvw = 1.0 if v == w else 0.0
    key = (e1, e2, e3)
    if key == (0, 0, 0):
        return (q * dv0 - 1) * (q * dw0 - 1) - (q * dvw - 1)
    if key == (1, 0, 0):
        return G[v] * (q * dw0 - 1) - G[add[v][neg[w]]]
    if key == (0, 1, 0):
        return G[w] * (q * dv0 - 1) - chi_m1 * G[add[v][neg[w]]]
    if key == (0, 0, 1):
        return G[w] * (q * dvw - 1) - G[v]
    if key == (1, 1, 0):
        return G[v] * G[w] - chi_m1 * (q * dvw - 1)
    if key == (1, 0, 1):
        return G[add[v][neg[w]]] * G[w] - (q * dv0 - 1)
    if key == (0, 1, 1):
        return G[v] * G[add[w][neg[v]]] - (q * dw0 - 1)
    if w == 0:
        return -G[v]
    if v == 0:
        return -G[w]
    C = 0j
    tw_skip = mul[neg[v]][inv[w]]
    for t in range(1, q):
        if t == neg[1] or t == tw_skip:
            continue
        C += chi[mul[mul[t][add[t][1]]][add[v][mul[t][w]]]]
    return G[1] * C


def S_eps_brute(F, e1, e2, e3, v, w) -> complex:
    q = F["q"]
    add, chi, mul = F["add"], F["chi"], F["mul"]
    acc = 0j
    for eta in range(1, q):
        for th in range(1, q):
            gam = add[eta][th]
            if gam == 0:
                continue
            c = 1
            if e1:
                c *= chi[eta]
            if e2:
                c *= chi[th]
            if e3:
                c *= chi[gam]
            acc += c * _e_tr(F, mul[eta][v]) * _e_tr(F, mul[th][w])
    return acc


def J_all_named(F, G, chi_om: int, v: int, w: int) -> complex:
    acc = 0j
    for e1 in (0, 1):
        for e2 in (0, 1):
            for e3 in (0, 1):
                acc += (chi_om ** (e1 + e2 + e3)) * S_eps_named(F, G, e1, e2, e3, v, w)
    return acc / 8.0


def J_all_named_drop111(F, G, chi_om: int, v: int, w: int) -> complex:
    acc = 0j
    for e1 in (0, 1):
        for e2 in (0, 1):
            for e3 in (0, 1):
                if (e1, e2, e3) == (1, 1, 1):
                    continue
                acc += (chi_om ** (e1 + e2 + e3)) * S_eps_named(F, G, e1, e2, e3, v, w)
    return acc / 8.0


def J_all_brute(F, Omega, v: int, w: int) -> complex:
    Om = set(Omega)
    add, mul = F["add"], F["mul"]
    acc = 0j
    for eta in Omega:
        for th in Omega:
            if add[eta][th] not in Om:
                continue
            acc += _e_tr(F, mul[eta][v]) * _e_tr(F, mul[th][w])
    return acc


def double_rs(p: int) -> tuple[int, ...]:
    half = pow(2, p - 2, p)
    return (1, (-2) % p, (-half) % p)


def J_1line_named(F, chi_om: int, G, u: int, v: int) -> complex:
    add, neg, mul = F["add"], F["neg"], F["mul"]
    acc = 0j
    for c in range(1, F["p"]):
        if (1 + c) % F["p"] == 0:
            continue
        arg = add[mul[neg[c]][u]][mul[(1 + c) % F["p"]][v]]
        acc += S_Omega_named(F, chi_om, G, arg)
    return acc


def J_double_named(F, chi_om: int, G, u: int, v: int) -> complex:
    add, neg, mul = F["add"], F["neg"], F["mul"]
    acc = 0j
    for c in double_rs(F["p"]):
        arg = add[mul[neg[c]][u]][mul[(1 + c) % F["p"]][v]]
        acc += S_Omega_named(F, chi_om, G, arg)
    return acc


def J_T_live(F, Omega, u: int, v: int) -> dict[str, complex]:
    add, mul, inv, neg = F["add"], F["mul"], F["inv"], F["neg"]
    e1 = _e_tab(F)
    Om = set(Omega)
    out = {"double": 0j, "generic": 0j, "n3": 0j, "all": 0j}
    for eta in Omega:
        for th in Omega:
            gam = add[eta][th]
            if gam not in Om:
                continue
            val = np.conj(e1[mul[th][u]]) * e1[mul[gam][v]]
            knd = triple_kind(F["p"], mul[th][inv[eta]])
            if knd != "degenerate":
                out[knd] += val
            out["all"] += val
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    samples = [(0, 0), (1, 0), (0, 1), (1, 1), (1, 2), (2, 3), (3, 1)]
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        G = gchi_table(F)
        Omega = omega_set_local(F)
        chi_om = int(F["chi"][Omega[0]])
        err_S = {f"{e1}{e2}{e3}": 0.0 for e1 in (0, 1) for e2 in (0, 1) for e3 in (0, 1)}
        err_J = err_drop = err_n8 = 0.0
        for v, w in samples:
            if v >= F["q"] or w >= F["q"]:
                continue
            for e1 in (0, 1):
                for e2 in (0, 1):
                    for e3 in (0, 1):
                        err_S[f"{e1}{e2}{e3}"] = max(
                            err_S[f"{e1}{e2}{e3}"],
                            abs(S_eps_named(F, G, e1, e2, e3, v, w) - S_eps_brute(F, e1, e2, e3, v, w)),
                        )
            jn = J_all_named(F, G, chi_om, v, w)
            jb = J_all_brute(F, Omega, v, w)
            err_J = max(err_J, abs(jn - jb))
            err_drop = max(err_drop, abs(J_all_named_drop111(F, G, chi_om, v, w) - jb))
            err_n8 = max(err_n8, abs(8 * jn - jb))
        row_ok = (
            all(err_S[k] < 1e-7 for k in err_S)
            and err_J < 1e-7
            and err_drop > 1.0
            and err_n8 > 1.0
        )
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "maxerr_S": {k: float(v) for k, v in err_S.items()},
            "maxerr_J": float(err_J),
            "drop111": float(err_drop),
            "drop_eighth": float(err_n8),
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "J_all=(1/8)Σ χ_Ω^{|ε|} S_ε. Fail drop S_{111}; fail drop 1/8.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    samples = [(1, 2), (1, 1), (3, 5), (2, 0), (4, 1)]
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        G = gchi_table(F)
        Omega = omega_set_local(F)
        chi_om = int(F["chi"][Omega[0]])
        e1 = e2 = e3 = eall = 0.0
        for u, v in samples:
            if u >= F["q"] or v >= F["q"]:
                continue
            live = J_T_live(F, Omega, u, v)
            j1 = J_1line_named(F, chi_om, G, u, v)
            jd = J_double_named(F, chi_om, G, u, v)
            jall = J_all_named(F, G, chi_om, u, v)
            e1 = max(e1, abs(j1 - live["double"] - live["generic"]))
            e2 = max(e2, abs(jd - live["double"]))
            e3 = max(e3, abs((jall - j1) - live["n3"]))
            eall = max(eall, abs(j1 - live["all"]))
        row_ok = e1 < 1e-7 and e2 < 1e-7 and e3 < 1e-7 and eall > 1.0
        # J_double vs J_1line: at p=5 generic is empty so they coincide
        if p != 5:
            F = _kernel_field(p)
            G = gchi_table(F)
            chi_om = int(F["chi"][omega_set_local(F)[0]])
            ddiff = abs(
                J_1line_named(F, chi_om, G, 1, 2) - J_double_named(F, chi_om, G, 1, 2)
            )
            if ddiff < 1.0:
                row_ok = False
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "err_1line": float(e1),
            "err_double": float(e2),
            "err_n3": float(e3),
            "J1_eq_Jall_err": float(eall),
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "J_1line=Σ_{c∈F_p^×\\{-1}} S_Ω((1+c)v−cu). "
            "Fail J_1line=J_all. Fail J_double=J_1line (p≥7)."
        ),
    }


@lru_cache(None)
def fold_G3(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    Y = load_Y(p, q)
    N = Y.shape[0]
    Omega = omega_set_local(F)
    chi_om = int(F["chi"][Omega[0]])
    G = gchi_table(F)
    add, neg, chi = F["add"], F["neg"], F["chi"]
    G3 = (Y * Y[:, 0:1]).T @ Y
    S = np.array([S_Omega_named(F, chi_om, G, d) for d in range(q)])
    coef = N / (p**3)
    err = err_g02 = err_oneQ = 0.0
    n = 0
    q3s = LIVE_Q3[p]
    for u in range(1, q):
        for v in range(1, q):
            if u == v:
                continue
            w = add[v][neg[u]]
            if w == 0:
                continue
            g02 = coef * (1 + 2 * S[u] + 2 * S[v] + 2 * S[w])
            live = J_T_live(F, Omega, u, v)
            rem = 0j
            rem1 = 0j
            q3d = q3s["double"]
            for knd, q3 in (
                ("double", q3s["double"]),
                ("generic", q3s["generic"]),
                ("n3", q3s["n3"]),
            ):
                if q3 is None:
                    continue
                rem += q3 * live[knd]
                rem1 += q3d * live[knd]
            pred = g02 + rem / (q**3)
            got = float(np.real(G3[u, v]))
            col = is_collinear(F, u, v)
            g31 = g3_1d_named(p, int(chi[u]), int(chi[v]), int(chi[w]), col)
            fr = got - g31
            err = max(err, abs(got - float(np.real(pred))))
            err_g02 = max(err_g02, abs(fr - float(np.real(g02 - g31))))
            err_oneQ = max(err_oneQ, abs(got - float(np.real(g02 + rem1 / (q**3)))))
            n += 1
    return {
        "N": int(N),
        "n": n,
        "err_identity": float(err),
        "err_G02_is_free": float(err_g02),
        "err_one_Q3": float(err_oneQ),
        "ok": bool(err < 1e-6 and err_g02 > 1.0 and err_oneQ > 1.0),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = fold_G3(p)
        if not rec["ok"]:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "G3=G02+q^{-3}Σ_T Q3_T J_T. Fail G02=Φ_free. Fail one Q3."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Phi_free_pointwise_named": False,
        "Q3_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "J_all and J_T are named complete sums. Q3_T splits "
            "(p=5: −5000/−25000; p=7: three values) with den 13,409, "
            "not a Gauss/Jacobi integer in p. Φ_free still carries Q3_T. "
            "Q3_02=−4N(2p²+1)/p follows from F=−2(3p²+2) certified "
            "field-only at p=5,7,11,13,17,19,23 (not a 2-point fit; "
            "not proved for all p). Do not interpolate (p−5)/15. "
            "16pA unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.558  J_all/J_T named; Q3_T not a p-law", flush=True)
    A = prove_A()
    print(f"  A J_all 8-term: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B J_T via S_Ω: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C G3 identity: {C['proved']}", flush=True)
    D = prove_open()
    out = {
        "prop": "15.558",
        "title": "J_all and J_T named; Φ_free still carries Q3_T",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "J_all_8_term": A["proved"],
            "J_T_S_Omega": B["proved"],
            "G3_Q3_identity": C["proved"],
            "Q3_named_in_p": False,
            "Phi_free_pointwise_named": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "J_all=(1/8)Σ χ_Ω^{|ε|} S_ε; "
            "J_1line=Σ_{c∈F_p^×\\{-1}} S_Ω((1+c)v−cu); "
            "G3=G02+q^{-3}Σ_T Q3_T J_T"
        ),
        "obstruction": "Q3_T not a p-law (den 13,409); Φ_free carries Q3_T",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
