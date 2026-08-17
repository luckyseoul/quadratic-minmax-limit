#!/usr/bin/env python3
"""
Prop 15.313 — Missing 1176-orbit is the 15.141 size-12 partner:
Q++=8/3, Q_N*=4 via ẑ_AP−2Γ_T.  p=7 eight-type mixture recovers
live Q_τ and ⟨δ,ψ⟩≤2.  Not a general-p name.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.312 flags.  Floor stays OPEN.

============================================================================
Theorem A — PROVED (Fourier; certified p=7).
  y_♯=y_AP ⊙ z_T, T=T_FIELD_SIZE12 (15.141).  ẑ_{y_♯}=ẑ_AP−2Γ_T.
  Ω-average Q_{++}=8/3, Q_{N*}=4, Q_{pm1}=8.  Fail-when-wrong: drop
  the 2; replace T by the fourth-power subgroup (not a C0-evec).  ∎

Theorem B — PROVED (certified p=7).
  Eight Max+-free types with named sizes |G|/|Stab| and Q_τ from
  Fejer/Gauss/ẑ−2Γ mix to the live coarse Q:
      Q++=1544/409, Q_N*=1496/409, Q_{--N1}=1376/409.
  Fail-when-wrong: drop the size-12 orbit.  ∎

Theorem C — CERTIFIED p=5,7; OPEN for all p≥5.
  Named Q_τ ⇒ min_k λ(ψ_k)>6 and ⟨δ,ψ⟩<2 at p=5 (15.311 mixture)
  and p=7 (Thm B).  No general-p type list.  phi_F not imported.  ∎

============================================================================
Backend: CPU ẑ on |Ω|=24; reuses 15.310–15.312 constructions.
Writes evidence/e1_gmin_m4_prop15313.json
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

from e1_gmin_m4_prop15141 import T_FIELD_SIZE12, construct_y_sharp  # noqa: E402
from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import (  # noqa: E402
    J_row,
    classes_of,
    dlog_T,
    primitive_root_T,
)
from e1_gmin_m4_prop15308 import G_order  # noqa: E402
from e1_gmin_m4_prop15310 import Q_AP_Fejer, Q_QR0_Gauss  # noqa: E402
from e1_gmin_m4_prop15311 import Q_ystar_Omega  # noqa: E402
from e1_gmin_m4_prop15312 import Q_switch, apnc_gen2_pack, qr0_nc_pack  # noqa: E402


def _fourths_is_evec() -> bool:
    from e1_gmin_m4_prop15139 import seidel_switch
    from minmax_quadratic import paley_conference_prime_power

    p = 7
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = affine_halfspace(p, (0, 1), {2, 3, 4, 5})
    n = p * p + 1
    z = np.ones(n)
    for u in fourths(p):
        z[1 + u] = -1.0
    C0 = seidel_switch(C, y0)
    return bool(np.allclose(C0 @ z, p * z, atol=1e-5))


def fourths(p: int) -> list[int]:
    F = _kernel_field(p)
    out = set()
    for x in range(1, F["q"]):
        y = x
        for _ in range(3):
            y = F["mul"][y][x]
        out.add(y)
    return sorted(out)


def Q_size12(p: int = 7) -> tuple[dict, float, float]:
    F = _kernel_field(p)
    q = F["q"]
    q2 = float(q * q)
    Omega = omega_set(F)
    rec = construct_y_sharp(p)
    y0 = affine_halfspace(p, (0, 1), {2, 3, 4, 5})
    z0 = np.asarray(y0[1 : 1 + q], dtype=np.float64)
    zy = np.asarray(rec["y"][1 : 1 + q], dtype=np.float64)
    T = list(T_FIELD_SIZE12)
    zh0 = zhat_omega(z0, F, Omega)
    Gam = np.zeros(len(Omega), dtype=np.complex128)
    for i, xi in enumerate(Omega):
        acc = 0j
        for u in T:
            acc += z0[u] * _e_tr(F, F["mul"][xi][u])
        Gam[i] = acc
    zhy = zhat_omega(zy, F, Omega)
    err2 = float(np.max(np.abs(zhy - (zh0 - 2 * Gam))))
    err1 = float(np.max(np.abs(zhy - (zh0 - Gam))))
    u = np.abs(zh0 - 2 * Gam) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    acc: dict[str, list[float]] = defaultdict(list)
    from e1_gmin_m4_prop15300 import merge_cls

    for r in F["squares"]:
        sm = 0.0
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        acc[merge_cls(F, r)].append(sm / len(Omega) / q2)
    Q = {c: float(np.mean(vs)) for c, vs in acc.items()}
    return Q, err2, err1


def prove_size12_Q() -> dict:
    Q, err2, err1 = Q_size12(7)
    qpp = Fraction(Q["++"]).limit_denominator(10)
    qn = Fraction(Q["N*"]).limit_denominator(10)
    ok = err2 < 1e-9 and err1 > 1e-6
    if qpp != Fraction(8, 3) or qn != 4:
        ok = False
    if abs(Q["pm1"] - 8.0) > 1e-6:
        ok = False
    if _fourths_is_evec():
        ok = False  # fail-when-wrong: fourths is not a C0-evec
    # fourths is not a C0 partner of this y0 (construct would fail or Q differs)
    return {
        "proved": ok,
        "err2": err2,
        "err1": err1,
        "Qpp": str(qpp),
        "Qn": str(qn),
        "Qpm1": Q["pm1"],
        "T": list(T_FIELD_SIZE12),
        "fourths": fourths(7),
        "T_is_fourths": set(T_FIELD_SIZE12) == set(fourths(7)),
        "fourths_is_C0_evec": _fourths_is_evec(),
        "theorem": (
            "Size-12 partner: ẑ=ẑ_AP−2Γ_T, Q++=8/3, Q_N*=4, Q_pm1=8. "
            "T is not the fourth-power subgroup."
        ),
    }


def _class_avg_from_pred(p: int, pred_fn) -> dict[str, float]:
    F = _kernel_field(p)
    from e1_gmin_m4_prop15290 import in_subfield
    from e1_gmin_m4_prop15300 import merge_cls

    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        val = pred_fn(r) if in_subfield(F, r) else 0.0
        acc[merge_cls(F, r)].append(float(val))
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def _nc_hits_of(p: int, y0: np.ndarray) -> list[dict]:
    from e1_gmin_m4_prop15139 import seidel_switch
    from e1_gmin_m4_prop15312 import _norm_table
    from minmax_quadratic import paley_conference_prime_power

    F = _kernel_field(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    if y0[0] < 0:
        y0 = -y0
    Cbase = (y0[:, None] * C) * y0[None, :]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    hits = []
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
            if np.allclose(C @ y, p * y, atol=1e-5):
                hits.append({"y": y, "y0": y0, "circle": pts, "t": t, "c": c})
    return hits


_TYPES_CACHE: list | None = None


def p7_named_types() -> list[dict]:
    """Eight types: size and Q dict (Q/q²) from constructions."""
    global _TYPES_CACHE
    if _TYPES_CACHE is not None:
        return _TYPES_CACHE
    p = 7
    g = G_order(p)
    from e1_gmin_m4_prop15308 import build_G, orbit_size

    Qap = _class_avg_from_pred(p, lambda r: Q_AP_Fejer(p, r % p))
    Qqr = _class_avg_from_pred(p, lambda r: Q_QR0_Gauss(p, r % p))
    Qys = Q_ystar_Omega(p)
    pack1 = qr0_nc_pack(p)
    Q1, _, _ = Q_switch(p, pack1)
    Q2, _, _ = Q_switch(p, apnc_gen2_pack(p))
    Qs, _, _ = Q_size12(p)
    F = _kernel_field(p)
    G = build_G(F, True)
    hits = _nc_hits_of(p, pack1["y"])
    Q294 = None
    Qcrrr = None
    for pack in hits:
        Qh, _, _ = Q_switch(p, pack)
        orb = orbit_size(G, np.sign(pack["y"]).astype(np.int8))
        qpm = Qh.get("pm1", 0.0)
        if orb == g // (p + 1):
            Q294 = Qh
        if abs(qpm - 8.0) < 1e-6 and orb == g // 2:
            Qcrrr = Qh
    if Q294 is None or Qcrrr is None:
        raise RuntimeError("missing 294 or CRRR among QR0⊙NC⊙NC hits")
    _TYPES_CACHE = [
        {"name": "AP", "n": g // (4 * p), "Q": Qap},
        {"name": "QR0", "n": g // (p * (p - 1)), "Q": Qqr},
        {"name": "APNC", "n": g // 4, "Q": Qys},
        {"name": "QR0NC", "n": g // 2, "Q": Q1},
        {"name": "size12", "n": g // 2, "Q": Qs},
        {"name": "APNCgen2", "n": g // 2, "Q": Q2},
        {"name": "CRRR", "n": g // 2, "Q": Qcrrr},
        {"name": "o294", "n": g // (p + 1), "Q": Q294},
    ]
    return _TYPES_CACHE


def mix_Q(types: list[dict]) -> dict[str, Fraction]:
    N = sum(t["n"] for t in types)
    out = {}
    for key in ("++", "N*", "--N1"):
        num = 0.0
        for t in types:
            num += t["n"] * float(t["Q"].get(key, 0.0))
        out[key] = Fraction(num / N).limit_denominator(2000)
    return out


def prove_p7_mixture() -> dict:
    types = p7_named_types()
    N = sum(t["n"] for t in types)
    M = mix_Q(types)
    dropped = mix_Q([t for t in types if t["name"] != "size12"])
    ok = N == 5726
    if M["++"] != Fraction(1544, 409):
        ok = False
    if M["N*"] != Fraction(1496, 409):
        ok = False
    if M["--N1"] != Fraction(1376, 409):
        ok = False
    if dropped["++"] == Fraction(1544, 409):
        ok = False
    return {
        "proved": ok,
        "N": N,
        "Qpp": str(M["++"]),
        "Qn": str(M["N*"]),
        "Qm": str(M["--N1"]),
        "drop12_Qpp": str(dropped["++"]),
        "sizes": {t["name"]: t["n"] for t in types},
        "theorem": (
            "p=7 eight-type mix recovers live Q_τ. "
            "Dropping the size-12 orbit misses 1544/409."
        ),
    }


def prove_floor_p57() -> dict:
    # p=7: λ from J @ Q_mix
    p = 7
    F = _kernel_field(p)
    C = classes_of(F)
    names = list(C)
    gen = primitive_root_T(F)
    nT = (F["q"] - 1) // 2
    dlog = dlog_T(F, gen)
    ks = (0, 2, 8, 12)
    Jr = np.real(np.array([J_row(F, C, dlog, nT, k) for k in ks], dtype=np.complex128))
    types = p7_named_types()
    M = mix_Q(types)
    Qv = []
    for name in names:
        if name == "pm1":
            Qv.append(8.0)
        elif name == "++":
            Qv.append(float(M["++"]))
        elif name == "N*":
            Qv.append(float(M["N*"]))
        else:
            Qv.append(float(M["--N1"]))
    Qv = np.array(Qv)
    S = Jr @ Qv  # S_k / q² scale: Qv is Q/q² so S = λ
    lam = {int(k): float(S[i]) for i, k in enumerate(ks)}
    # ⟨δ,ψ⟩ = ∑_{τ≠pm1} (4 - Q_τ) J_τ
    # S_k = 8 J_pm1 + ∑_{off} Q_τ J_τ
    # ∑_{off} J_τ = -J_pm1 for k≠0 (∑_T ψ=0, J_pm1=2)
    # ⟨δ,ψ⟩ = ∑_off (4-Q) J = 4(-J_pm1) - (S - 8 J_pm1) = -8 + 8*2 - S wait
    # J_pm1(k≠0)= ψ(1)+ψ(-1)=2
    # ∑_all J = 0 for k≠0 so ∑_off J = -2
    # ⟨δ,ψ⟩ = ∑_off 4 J - ∑_off Q J = 4(-2) - (S - 8*2) = -8 - S + 16 = 8 - S
    # λ = S (since Qv is Q/q² and S=∑ (Q/q²)J = S_□/q² = λ)
    # ⟨δ,ψ⟩ = 8 - λ
    # need ⟨δ,ψ⟩≤2 ⇔ λ≥6.  Yes that's the 15.290 restatement.
    lams = [lam[k] for k in (2, 8, 12)]
    dmin = min(lams)
    ok = dmin > 6 - 1e-9
    # p=5 from 15.311: 80/13 > 6
    return {
        "proved": ok,
        "lambda_p7": lam,
        "lambda_min_p7": dmin,
        "delta_max_p7": 8 - dmin,
        "p5_lambda_min": float(Fraction(80, 13)),
        "general": False,
        "theorem": (
            "Named Q_τ ⇒ λ_min=3072/409>6 at p=7 and 80/13>6 at p=5. "
            "Not a type list for all p≥5. phi_F not imported."
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
        "note": (
            "p=5 and p=7 Q_τ named by Fejer/Gauss/ẑ−2Γ mixtures. "
            "No general-p NL type list (size-12 T is p=7-only; NC empty at p=13). "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.313  size-12 Q++=8/3; p=7 mix 1544/409; floor p=5,7 only", flush=True)
    A = prove_size12_Q()
    print(f"  A size12: {A['proved']} Q++={A['Qpp']} Qn={A['Qn']} fourths={A['T_is_fourths']}", flush=True)
    B = prove_p7_mixture()
    print(f"  B mix: {B['proved']} Q++={B['Qpp']} drop={B['drop12_Qpp']}", flush=True)
    C = prove_floor_p57()
    print(f"  C floor p=7 λmin={C['lambda_min_p7']:.4f} δmax={C['delta_max_p7']:.4f}", flush=True)
    D = prove_open()
    print(f"  D general: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.313",
        "title": "Size-12 orbit Q++=8/3; p=7 mix names Q_τ; floor only p=5,7",
        "proved": {
            "size12_Q": A["proved"],
            "p7_mixture": B["proved"],
            "floor_p5_p7": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15313.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
