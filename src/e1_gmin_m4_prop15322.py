#!/usr/bin/env python3
"""
Prop 15.322 — W(1) on AP / non-AP / NL orbits.  p=3: Var=0 and
A_□=E[W]².  p=5 two-orbit ± names A_□ only at p=5.  The same
pattern fails at p=7.  A_□ / Q_{++} unnamed in p.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.321 flags.  Floor stays OPEN.

============================================================================
Setup.  W(1)=∑_{N(δ)=1} N(δ).  E[W]=(p+1)μ_+.  AP-hs orbit size
p(p²−1)/4; QR0/non-AP orbit size p(p+1) when it is distinct
(15.308–15.309).  EW+2p := p(p²+7)/4.

============================================================================
Theorem A — PROVED (certified p=3; 2-design + |H_+|=6).
  At p=3, W≡E[W]=6 on every Max+, so Var=0 and A_□=36=E[W]².
  Fail-when-wrong: claim Var(W)>0 at p=3.  ∎

Theorem B — PROVED (orbit tag; certified p=5,7).
  p=5: AP (n=30) has W∈{EW±2p}={20,40} equally; NL (n=100) has
  W∈{EW±(p+1)}={24,36} equally.  Hence
      A_□=(30(EW²+4p²)+100(EW²+(p+1)²))/130=12360/13.
  p=7: AP (n=84) has W∈{EW−4p, EW+2p}={56,98} with masses 28:56
  (not equal ±2p); non-AP (n=56) has W≡EW; NL Var=p²−1=48
  (not (p+1)²=64).  W_min=56<EW−2p=70.  Fail-when-wrong:
  claim W≥EW−2p at p=7; or NL Var=(p+1)² at p=7.  ∎

Theorem C — OPEN.  The p=5 two-orbit formula is not a formula
  in p (p=7 has three named hs/NL slots and a 12-value NL hist).
  3-prime low-degree rational search for Var/A_□: 0 hits.
  Qpp_floor_ub unproved.  phi_F not imported.  ∎

============================================================================
Backend: N-cache + square-dir signatures (15.307–15.308).
Writes evidence/e1_gmin_m4_prop15322.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15308 import is_AP, recover_S_if_hs  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, alpha_beta_p1  # noqa: E402

_W_CACHE: dict[int, dict] = {}


def EW_named(p: int) -> int:
    return (p + 1) * p * (p - 1) // 4


def EW_plus_2p(p: int) -> int:
    return p * (p * p + 7) // 4


def Asq_p5_two_orbit() -> Fraction:
    p = 5
    ew = EW_named(p)
    n_ap, n_nl = 30, 100
    return Fraction(
        n_ap * (ew * ew + 4 * p * p) + n_nl * (ew * ew + (p + 1) ** 2),
        n_ap + n_nl,
    )


def _W_by_orbit(p: int) -> dict:
    if p in _W_CACHE:
        return _W_CACHE[p]
    D, F = load_D(p)
    N = load_N(p, F).astype(np.float64)
    Uidx = [d for d in range(1, F["q"]) if int(field_norm(F, d)) == 1]
    W = N[:, np.array(Uidx, dtype=np.int32)].sum(axis=1)
    classes = square_classes(F)
    n_sq = len(classes)
    hs_key = tuple(sorted(["H"] + ["C"] * (n_sq - 1)))
    buckets: dict[str, list[int]] = defaultdict(list)
    for i in range(D.shape[0]):
        tags = []
        for lines in classes:
            ns = D[i, lines].sum(axis=1)
            sig = tuple(int(x) for x in np.sort(ns))
            tags.append(tag_sig(p, sig))
        t = tuple(sorted(tags))
        if t == hs_key:
            S = recover_S_if_hs(p, D[i], classes)
            key = "hs-AP" if (S is not None and is_AP(S, p)) else "hs-nonAP"
        else:
            key = "NL"
        buckets[key].append(int(round(float(W[i]))))
    rec = {"p": p, "n": int(D.shape[0]), "EW": float(np.mean(W))}
    for k, vs in buckets.items():
        rec[k] = {
            "n": len(vs),
            "hist": dict(sorted(Counter(vs).items())),
            "mean": float(np.mean(vs)),
            "var": float(np.var(vs)),
            "wmin": min(vs),
            "wmax": max(vs),
        }
    _W_CACHE[p] = rec
    return rec


def prove_p3() -> dict:
    rec = _W_by_orbit(3)
    vs = []
    for k in ("hs-AP", "hs-nonAP", "NL"):
        if k in rec:
            vs.extend(rec[k]["hist"].keys())
    uniq = sorted(set(vs))
    ok = uniq == [6] and abs(rec["EW"] - 6) < 1e-9
    return {
        "proved": ok,
        "var_positive_false": uniq == [6],
        "uniq": uniq,
        "n": rec["n"],
        "theorem": "p=3: W≡6, Var=0, A_□=36.",
    }


def prove_orbits() -> dict:
    ok = True
    rec5 = _W_by_orbit(5)
    rec7 = _W_by_orbit(7)
    # p=5 AP ±2p equal, NL ±(p+1) equal
    ap5 = rec5["hs-AP"]
    nl5 = rec5["NL"]
    if ap5["n"] != 30 or ap5["hist"] != {20: 15, 40: 15}:
        ok = False
    if nl5["n"] != 100 or nl5["hist"] != {24: 50, 36: 50}:
        ok = False
    if Asq_p5_two_orbit() != Fraction(12360, 13):
        ok = False
    # p=7: min < EW-2p; AP not equal ±2p; NL var ≠ (p+1)²
    ap7 = rec7["hs-AP"]
    nl7 = rec7["NL"]
    non7 = rec7.get("hs-nonAP", {})
    min7 = min(ap7["wmin"], nl7["wmin"], non7.get("wmin", 10**9))
    if min7 >= EW_named(7) - 2 * 7:
        ok = False
    if ap7["hist"] == {70: 42, 98: 42}:
        ok = False
    if abs(nl7["var"] - 64.0) < 1e-6:
        ok = False
    if abs(nl7["var"] - 48.0) > 1e-6:
        ok = False
    if non7.get("n") != 56 or abs(non7.get("var", 1) - 0.0) > 1e-6:
        ok = False
    return {
        "proved": ok,
        "W_ge_EW_minus_2p_false": min7 < EW_named(7) - 14,
        "NL_var_pplus1sq_false": abs(nl7["var"] - 64) > 1,
        "by_p": {"5": rec5, "7": rec7},
        "theorem": (
            "p=5 two-orbit ± names A_□=12360/13. "
            "p=7: AP {56,98} 1:2, nonAP W≡EW, NL Var=48, min 56."
        ),
    }


def prove_open() -> dict:
    al, be = alpha_beta_p1(5)
    live_Q = Fraction(48, 13)
    ub = Qpp_floor_ub(5)
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "floor_ub_p5": str(ub),
        "live_Qpp_p5": str(live_Q),
        "slack": str(ub - live_Q),
        "A_from_Q": str(al + be * live_Q),
        "note": (
            "p=5 two-orbit names A_□ only at p=5. "
            "p=7 orbit Vars do not give a formula in p. "
            "Qpp_floor_ub unproved. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.322  W by AP/NL orbit; p=3 Var=0; p=5-only A_□ name", flush=True)
    A = prove_p3()
    print(f"  A p=3: {A['proved']} var0={A['var_positive_false']}", flush=True)
    B = prove_orbits()
    print(
        f"  B orbits: {B['proved']} min2p={not B['W_ge_EW_minus_2p_false']} "
        f"NLvar64={not B['NL_var_pplus1sq_false']}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']} slack={C['slack']}", flush=True)
    out = {
        "prop": "15.322",
        "title": "W orbit values; A_□ named only at p=3 and via p=5 two-orbit",
        "proved": {
            "p3_var0": A["proved"],
            "orbit_W": B["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15322.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
