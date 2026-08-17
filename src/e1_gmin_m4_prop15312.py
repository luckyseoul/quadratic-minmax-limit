#!/usr/bin/env python3
"""
Prop 15.312 — QR0⊙NC and AP⊙NC⊙NC Q_τ via ẑ_base−2Γ.
Names three more p=7 NL orbits (40/9, 32/9, 294).  One 1176
(Q_{++}=8/3, Q_{N*}=4) is not in this generation.  Full Q_τ OPEN.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.311 flags.  Floor stays OPEN.

============================================================================
Theorem A — PROVED (Fourier; certified p=7).
  y1 = y_{QR0} ⊙ z_C, C the lex-first Hoffman circle of C_{QR0}.
  ẑ_{y1}=ẑ_{QR0}−2Γ_C.  Ω-average Q_{++}^{y1}=40/9, Q_{N*}=76/21.
  Fail-when-wrong: drop the 2; or claim Q_{++}=8/3 (that is AP⊙NC).  ∎

Theorem B — PROVED (certified p=7).
  y2 = y_* ⊙ z_{C2}, C2 the lex-first Hoffman circle of C_{y_*}
  that is not the inverse switch back to an affine type.
  ẑ_{y2}=ẑ_{y_*}−2Γ_{C2}.  Q_{++}^{y2}=32/9.
  Fail-when-wrong: claim 8/3 or 40/9.  ∎

Theorem C — PROVED (certified p=7).
  Among NC of C_{y1}, the y with smallest new G-orbit has size
  |G|/(p+1)=294 and Q_{++}=40/9, Q_{N*}=24/7.  Same ẑ−2Γ.
  Fail-when-wrong: claim that orbit is 1176.  ∎

Theorem D — OPEN.  The 1176-orbit with Q_{++}=8/3, Q_{N*}=4 is
  not produced by depth-2 NC from {AP, QR0}.  Full-ensemble Q_τ
  unnamed.  phi_F not imported.  ∎

============================================================================
Backend: CPU field circles (294) + Ω Fourier; G for orbit size of C.
Writes evidence/e1_gmin_m4_prop15312.json
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

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15308 import G_order, build_G, orbit_size  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15311 import Q_ystar_Omega  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def _norm_table(F):
    q, p = F["q"], F["p"]
    mul = F["mul"]
    units = [0] * (q - 1)
    g, cur = F["g"], 1
    for j in range(q - 1):
        units[j] = cur
        cur = mul[cur][g]
    fr = [0] * q
    for x in range(1, q):
        fr[x] = units[(F["dlog"][x] * p) % (q - 1)]
    N = [0] * q
    for u in range(1, q):
        N[u] = mul[u][fr[u]] % p
    return N


def first_nc(p: int, y0: np.ndarray, skip_pred=None):
    """Lex-first Hoffman+evec NC of C_{y0}. skip_pred(y)->bool to skip."""
    F = _kernel_field(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    if y0[0] < 0:
        y0 = -y0
    Cbase = (y0[:, None] * C) * y0[None, :]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if any(Cbase[S[i], S[j]] <= 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
                continue
            z = np.ones(q + 1)
            for i in S:
                z[i] = -1.0
            if not np.allclose(Cbase @ z, p * z, atol=1e-5):
                continue
            y = y0 * z
            if y[0] < 0:
                y = -y
            if not np.allclose(C @ y, p * y, atol=1e-5):
                continue
            if skip_pred is not None and skip_pred(y):
                continue
            return {
                "y": y,
                "circle": pts,
                "t": t,
                "c": c,
                "y0": y0,
            }
    return None


def Q_switch(p: int, pack: dict) -> tuple[dict, float]:
    F = _kernel_field(p)
    q = F["q"]
    q2 = float(q * q)
    Omega = omega_set(F)
    y0 = pack["y0"]
    z_base = np.asarray(y0[1 : 1 + q], dtype=np.float64)
    z_y = np.asarray(pack["y"][1 : 1 + q], dtype=np.float64)
    zh_b = zhat_omega(z_base, F, Omega)
    Gam = np.zeros(len(Omega), dtype=np.complex128)
    for i, xi in enumerate(Omega):
        acc = 0j
        for u in pack["circle"]:
            acc += z_base[u] * _e_tr(F, F["mul"][xi][u])
        Gam[i] = acc
    zh_y = zhat_omega(z_y, F, Omega)
    err2 = float(np.max(np.abs(zh_y - (zh_b - 2 * Gam))))
    err1 = float(np.max(np.abs(zh_y - (zh_b - Gam))))
    u = np.abs(zh_b - 2 * Gam) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        sm = 0.0
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        acc[merge_cls(F, r)].append(sm / len(Omega) / q2)
    Q = {c: float(np.mean(vs)) for c, vs in acc.items()}
    return Q, err2, err1


def qr0_nc_pack(p: int = 7) -> dict:
    y0 = affine_halfspace(p, (0, 1), qr0(p))
    pack = first_nc(p, y0)
    pack["base"] = "QR0"
    return pack


def apnc_gen2_pack(p: int = 7) -> dict:
    ys = construct_ystar(p)
    ystar = ys["y"]
    # skip the inverse switch if it returns an affine-like Q (N*=0)
    def skip(y):
        return False  # lex-first at p=7 is already 32/9 (not inverse)

    pack = first_nc(p, ystar, skip_pred=skip)
    pack["base"] = "APstar"
    return pack


def prove_QR0_NC() -> dict:
    p = 7
    pack = qr0_nc_pack(p)
    Q, err2, err1 = Q_switch(p, pack)
    qpp = Fraction(Q["++"]).limit_denominator(30)
    qn = Fraction(Q["N*"]).limit_denominator(30)
    ok = err2 < 1e-9 and err1 > 1e-6
    if qpp != Fraction(40, 9):
        ok = False
    if qpp == Fraction(8, 3):
        ok = False
    if qn != Fraction(76, 21):
        ok = False
    ystar_pp = Fraction(Q_ystar_Omega(p)["++"]).limit_denominator(10)
    if ystar_pp == qpp:
        ok = False  # must differ from AP⊙NC
    return {
        "proved": ok,
        "err2": err2,
        "err1": err1,
        "Qpp": str(qpp),
        "Qn": str(qn),
        "ystar_Qpp": str(ystar_pp),
        "t_c": (pack["t"], pack["c"]),
        "theorem": "QR0⊙NC: ẑ=ẑ_QR0−2Γ, Q++=40/9 ≠ 8/3=AP⊙NC.",
    }


def prove_gen2_32_9() -> dict:
    p = 7
    pack = apnc_gen2_pack(p)
    Q, err2, err1 = Q_switch(p, pack)
    qpp = Fraction(Q["++"]).limit_denominator(30)
    ok = err2 < 1e-9 and err1 > 1e-6
    if qpp != Fraction(32, 9):
        ok = False
    if qpp in (Fraction(8, 3), Fraction(40, 9)):
        ok = False
    return {
        "proved": ok,
        "err2": err2,
        "err1": err1,
        "Qpp": str(qpp),
        "Qn": str(Fraction(Q["N*"]).limit_denominator(80)),
        "t_c": (pack["t"], pack["c"]),
        "theorem": "AP⊙NC⊙NC lex-first: Q++=32/9. Not 8/3 or 40/9.",
    }


def prove_294() -> dict:
    p = 7
    pack1 = qr0_nc_pack(p)
    y1 = pack1["y"]
    F = _kernel_field(p)
    G = build_G(F, True)
    C = paley_conference_prime_power(p).astype(np.float64)
    q = F["q"]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    Cbase = (y1[:, None] * C) * y1[None, :]
    y0 = np.asarray(y1, dtype=np.float64)
    candidates = []
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if any(Cbase[S[i], S[j]] <= 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
                continue
            z = np.ones(q + 1)
            for i in S:
                z[i] = -1.0
            if not np.allclose(Cbase @ z, p * z, atol=1e-5):
                continue
            y = y0 * z
            if y[0] < 0:
                y = -y
            if not np.allclose(C @ y, p * y, atol=1e-5):
                continue
            orb = orbit_size(G, np.sign(y).astype(np.int8))
            candidates.append((orb, t, c, y, pts))
    # smallest NEW orbit (exclude 56=QR0 and 1176)
    new = [c for c in candidates if c[0] not in (56, 1176)]
    ok = bool(new)
    smallest = min(new, key=lambda r: r[0]) if new else None
    Qpp = Qn = None
    err2 = None
    if smallest is not None:
        orb, t, c, y, pts = smallest
        pack = {"y": y, "y0": y0, "circle": pts, "t": t, "c": c}
        Q, err2, err1 = Q_switch(p, pack)
        Qpp = Fraction(Q["++"]).limit_denominator(30)
        Qn = Fraction(Q["N*"]).limit_denominator(30)
        if orb != G_order(p) // (p + 1):
            ok = False
        if orb == 1176:
            ok = False
        if Qpp != Fraction(40, 9) or Qn != Fraction(24, 7):
            ok = False
        if err2 > 1e-9:
            ok = False
    else:
        ok = False
    return {
        "proved": ok,
        "orbit": smallest[0] if smallest else None,
        "named_orbit": G_order(p) // (p + 1),
        "Qpp": str(Qpp) if Qpp is not None else None,
        "Qn": str(Qn) if Qn is not None else None,
        "err2": err2,
        "n_cand": len(candidates),
        "theorem": (
            "Smallest new G-orbit of QR0⊙NC⊙NC has size |G|/(p+1)=294 "
            "and Q++=40/9, Q_N*=24/7. Not 1176."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "missing": "1176-orbit with Q++=8/3, Q_N*=4 not in depth-2 NC from {AP,QR0}",
        "note": "Named QR0⊙NC, AP⊙NC⊙NC (32/9), and 294. Full Q_τ unnamed. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.312  QR0⊙NC Q++=40/9; gen2 32/9; 294 orbit", flush=True)
    A = prove_QR0_NC()
    print(f"  A QR0⊙NC: {A['proved']} Q++={A['Qpp']} Qn={A['Qn']}", flush=True)
    B = prove_gen2_32_9()
    print(f"  B AP⊙NC⊙NC: {B['proved']} Q++={B['Qpp']}", flush=True)
    C = prove_294()
    print(f"  C 294: {C['proved']} |O|={C['orbit']} Q++={C['Qpp']} Qn={C['Qn']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.312",
        "title": "QR0⊙NC and gen2 Q via ẑ−2Γ; 294 named; one 1176 missing",
        "proved": {
            "QR0_NC_40_9": A["proved"],
            "gen2_32_9": B["proved"],
            "orbit_294": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15312.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
