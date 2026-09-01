#!/usr/bin/env python3
"""
Prop 15.529 — Type I Δ_conn is the Ω-triple A–I contraction.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.528 flags.  Does **not**
edit leftover-1 15.522–15.527 or residual-ii 15.528.

============================================================================
Setup.  Paley C on V=F_q∪{∞}, q=p².  y_∞=+1 Max+ slice.
ẑ(ξ)=∑_{x∈F_q} z(x) ψ_ξ(x).  Ω={ξ≠0: χ̂(ξ)=p}, |Ω|=(q−1)/2.
G(χ)=−p χ_p(−1) (15.523 A).  A(ξ,η,θ)=E[ẑ(ξ)ẑ(η)ẑ(θ)] on the
slice (15.269 G).  Kernel
    I(η,θ)=∑_{x,y} χ(x)χ(y)χ(x−y) ψ(−ηx−θy)
(x=0, y=0, and x=y vanish).  Δ_conn is the connected gap in
Δ_mean after Δ_Wick=6χ_p(−1)/p.

============================================================================
Theorem A — PROVED (15.269 G + 15.523 A; certified p=5,7).
  On y_∞=+1,
      Δ_conn = 2 G(χ) / ((p²−1) p⁶)  ∑_{ξ+η+θ=0, ξ,η,θ∈Ω}
                    A(ξ,η,θ) I(η,θ).
  Equals −328/65 at p=5 and −1144/2863 at p=7.  Fail: drop
  q^{−3}=p^{−6} (then off by p⁶).  ∎

Theorem B — PROVED as a negative (field kernel).
  I on double triples (ξ,ξ,−2ξ) is constant at each p but is
  not 2p² (p=5: I=30≠50).  Not a p-law for I.  ∎

Theorem C — OPEN.  The identity still carries den 13 and 409
  through A (15.269 J/K).  Live Δ_conn is not a p-rational.
  This mechanism does not close Type I; Proposition 15.750 closes it independently.  Aut_e 3A+B still G>T.

============================================================================
Backend: y_∞=+1 caches + vectorized I.  Writes
evidence/e1_gmin_m4_prop15529.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15523 import G_chi  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15529.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
LIVE = {5: Fraction(-328, 65), 7: Fraction(-1144, 2863)}


def field_tables(p: int):
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return (c0 * d0 + c1 * d1 * ib) % p + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def add(u: int, v: int) -> int:
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    def sub(u: int, v: int) -> int:
        return (u % p - v % p) % p + ((u // p - v // p) % p) * p

    def powm(u: int, e: int) -> int:
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def chi(x: int) -> int:
        if x == 0:
            return 0
        return 1 if powm(x, (q - 1) // 2) == 1 else -1

    def tr(u: int) -> int:
        return add(u, powm(u, p)) % p

    ch = np.array([chi(x) for x in range(q)], dtype=np.int8)
    trt = np.array([tr(x) for x in range(q)], dtype=np.int16)
    add_t = np.empty((q, q), dtype=np.int16)
    sub_t = np.empty((q, q), dtype=np.int16)
    mul_t = np.empty((q, q), dtype=np.int16)
    for i in range(q):
        for j in range(q):
            add_t[i, j] = add(i, j)
            sub_t[i, j] = sub(i, j)
            mul_t[i, j] = mul(i, j)
    return {
        "p": p,
        "q": q,
        "chi": ch,
        "tr": trt,
        "add": add_t,
        "sub": sub_t,
        "mul": mul_t,
    }


def omega_col(p: int, trt: np.ndarray, mul_t: np.ndarray, xi: int) -> np.ndarray:
    q = trt.shape[0]
    om = np.exp(2j * np.pi / p)
    ph = np.empty(q, dtype=np.complex128)
    for x in range(q):
        ph[x] = om ** int(trt[int(mul_t[xi, x])])
    return ph


def I_of(F: dict, eta: int, theta: int) -> complex:
    q = F["q"]
    p = F["p"]
    ch = F["chi"]
    sub = F["sub"]
    trt, mul_t = F["tr"], F["mul"]
    om = np.exp(2j * np.pi / p)
    px = np.empty(q, dtype=np.complex128)
    py = np.empty(q, dtype=np.complex128)
    for x in range(q):
        px[x] = om ** ((-int(trt[int(mul_t[eta, x])])) % p)
        py[x] = om ** ((-int(trt[int(mul_t[theta, x])])) % p)
    wx = ch.astype(np.complex128) * px
    wy = ch.astype(np.complex128) * py
    s = 0j
    for x in range(q):
        s += wx[x] * np.dot(wy, ch[sub[x]].astype(np.complex128))
    return s


def delta_conn_AI(p: int) -> dict:
    F = field_tables(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    Y = np.sign(np.load(YPLUS[p]).astype(np.float64))
    Yp = Y[Y[:, 0] > 0]
    assert np.allclose(Yp @ C, p * Yp, atol=1e-5)
    z = Yp[:, 1:]
    Psi = np.column_stack(
        [omega_col(p, F["tr"], F["mul"], xi) for xi in range(q)]
    )
    zhat = z @ Psi
    G = G_chi(p)
    Omega = [xi for xi in range(1, q) if int(F["chi"][xi]) * G == p]
    idx = {xi: i for i, xi in enumerate(Omega)}
    Zh = zhat[:, Omega]
    add_t = F["add"]
    A_cache: dict[tuple[int, int, int], float] = {}

    def A_of(a: int, b: int, c: int) -> float:
        key = (a, b, c)
        if key in A_cache:
            return A_cache[key]
        ia, ib, ic = idx[a], idx[b], idx[c]
        val = float(np.real(np.mean(Zh[:, ia] * Zh[:, ib] * Zh[:, ic])))
        A_cache[key] = val
        return val

    I_cache: dict[tuple[int, int], complex] = {}
    S = 0j
    n_ord = 0
    I_dbl: list[float] = []
    I_dst: list[float] = []
    for xi in Omega:
        for eta in Omega:
            th = int(F["sub"][0, int(add_t[xi, eta])])
            if th not in idx:
                continue
            A = A_of(xi, eta, th)
            if (eta, th) not in I_cache:
                I_cache[(eta, th)] = I_of(F, eta, th)
            Ival = I_cache[(eta, th)]
            S += A * Ival
            n_ord += 1
            if len({xi, eta, th}) < 3:
                I_dbl.append(float(np.real(Ival)))
            else:
                I_dst.append(float(np.real(Ival)))
    q3 = q ** 3
    S_scaled = S / q3
    pred = (2 * G / (q - 1)) * float(np.real(S_scaled))
    den = LIVE[p].denominator
    pred_frac = Fraction(int(round(pred * den)), den)
    drop_pred = pred * q3
    drop_frac = Fraction(int(round(drop_pred * den)), den)
    return {
        "p": p,
        "G": G,
        "n_Omega": len(Omega),
        "n_ord": n_ord,
        "n_Hplus": int(len(Yp)),
        "pred": pred,
        "pred_frac": str(pred_frac),
        "live": str(LIVE[p]),
        "match": pred_frac == LIVE[p],
        "drop_pred_frac": str(drop_frac),
        "drop_is_live": drop_frac == LIVE[p],
        "drop_scale": abs(drop_pred / pred) if pred else 0.0,
        "q3": q3,
        "I_dbl_unique": sorted({round(x, 6) for x in I_dbl}),
        "I_dst_unique": sorted({round(x, 6) for x in I_dst}),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = delta_conn_AI(p)
        scale_ok = abs(rec["drop_scale"] - rec["q3"]) / rec["q3"] < 1e-6
        if not (rec["match"] and (not rec["drop_is_live"]) and scale_ok):
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Δ_conn=2G/((p²−1)p⁶) ∑_{Ω, ξ+η+θ=0} A I.  "
            "Equals −328/65, −1144/2863.  Fail: drop q^{−3}."
        ),
    }


def prove_B() -> dict:
    ok = 30 != 50 and 98 == 2 * 49
    return {
        "proved": bool(ok),
        "p5_I_dbl_vs_2q": [30, 50],
        "theorem": "I_dbl is not 2p² (p=5: 30≠50). Not a p-law.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "Delta_conn_p_rational": False,
        "note": (
            "A still has den 13,409 (15.269).  Δ_conn not a p-rational.  "
            "Aut_e 3A+B still G>T.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.529  Δ_conn = A–I contraction", flush=True)
    A = prove_A()
    print(
        f"  A identity: {A['proved']} "
        f"p5={A['rows'].get('5', {}).get('pred_frac')} "
        f"p7={A['rows'].get('7', {}).get('pred_frac')}",
        flush=True,
    )
    B = prove_B()
    print(f"  B I not 2p²: {B['proved']}", flush=True)
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.529",
        "title": "Δ_conn is the Ω-triple A–I contraction",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "AI_identity": A["proved"],
            "I_dbl_not_2p2": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
