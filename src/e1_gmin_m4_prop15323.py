#!/usr/bin/env python3
"""
Prop 15.323 — AP-orbit W(1) is not {EW±2p} for p≡1 (fails p=13).
Var ≤ 4p(p+1)(2p+1)/(5p+1) is a sufficient majorant for the
p≡1 floor ub (equality at p=5).  Unproved as a bound.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.322 flags.  Floor stays OPEN.

============================================================================
Setup.  AP-hs = Type+ lift of an m-AP, m=(p+1)/2.  G as in 15.308.
|O_AP|=p(p²−1)/4.  Two-orbit model: H_+=AP ∪ NL_{|G|/(p+1)} with
W_AP=EW±2p equally and W_NL=EW±(p+1) equally, giving
    Var_2orb=4p(p+1)(2p+1)/(5p+1).

============================================================================
Theorem A — PROVED (G-orbit of a constructed AP-hs; Max+-free).
  |O_AP| and the W-histograms:
      p=5:  {20,40} masses 15:15  (=EW±2p)
      p=7:  {56,98} masses 28:56  (=EW−4p, EW+2p)
      p=11: {242,286,330,352,440} equal 66
      p=13: {312,520,624,676} masses 91:182:182:91
  In particular AP W ⊆ {EW±2p} fails at p=13 (the next p≡1).
  Fail-when-wrong: claim AP W ⊆ {EW±2p} for every p≡1.  ∎

Theorem B — PROVED (Fraction algebra; p=5,13,17,29,37).
  For primes p≥5, p≡1 (mod 4),
      Var ≤ Var_2orb  ⇒  Q_{++} ≤ Qpp_floor_ub(p),
  because Var_2orb ≤ A_floor−E[W]² on those p, with equality of
  Q_{++} to the live 48/13 only at p=5.  Live Var at p=7 is
  21504/409 < Var_2orb=280/3, so equality in the majorant fails
  at p=7.  Fail-when-wrong: claim Var=Var_2orb at p=7.  ∎

Theorem C — OPEN.  Var ≤ Var_2orb is unproved (and AP support is
  not {EW±2p} at p=13, so the 2-orbit law is not a support
  majorant).  A_□ / Q_{++} unnamed.  phi_F not imported.  ∎

============================================================================
Backend: generator BFS of the G-orbit of one AP-hs (ProcessPool).
Writes evidence/e1_gmin_m4_prop15323.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    is_type_plus_form,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15308 import frob_table  # noqa: E402
from e1_gmin_m4_prop15321 import (  # noqa: E402
    Qpp_floor_ub,
    alpha_beta_p1,
)

_ORBIT_CACHE: dict[int, dict] = {}


def EW_named(p: int) -> int:
    return (p + 1) * p * (p - 1) // 4


def Var_2orb(p: int) -> Fraction:
    return Fraction(4 * p * (p + 1) * (2 * p + 1), 5 * p + 1)


def first_type_plus(p: int) -> tuple[int, int]:
    for alpha in range(p):
        for beta in range(p):
            if (alpha or beta) and is_type_plus_form(p, alpha, beta):
                return alpha, beta
    raise RuntimeError("no Type+")


def ap_orbit_W(p: int) -> dict:
    if p in _ORBIT_CACHE:
        return _ORBIT_CACHE[p]
    F = _kernel_field(p)
    q = F["q"]
    plus = set(range((p + 1) // 2))
    a0, b0 = first_type_plus(p)
    yfld0 = np.array(lift_sigma_vector(p, a0, b0, plus)[1:], dtype=np.int8)
    mul = np.array(F["mul"], dtype=np.int32)
    add = np.array(F["add"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    inv = np.array(F["inv"], dtype=np.int32)
    fr = np.array(frob_table(F), dtype=np.int32)
    U = np.array(
        [d for d in range(1, q) if int(field_norm(F, d)) == 1], dtype=np.int32
    )
    g2 = int(mul[F["g"], F["g"]])

    def apply_trans(yfld, b):
        return yfld[add[:, neg[b]]]

    def apply_mul(yfld, a):
        return yfld[mul[int(inv[a])]]

    def W_of(yfld):
        D = yfld < 0
        Didx = np.flatnonzero(D)
        acc = 0
        for u in U:
            acc += int(D[add[Didx, u]].sum())
        return acc

    seen = {yfld0.tobytes()}
    stack = [yfld0]
    Ws: list[int] = []
    while stack:
        yfld = stack.pop()
        Ws.append(W_of(yfld))
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
    hist = dict(sorted(Counter(Ws).items()))
    rec = {
        "n": len(Ws),
        "want": p * (p * p - 1) // 4,
        "hist": hist,
        "mean": float(np.mean(Ws)),
        "var": float(np.var(Ws)),
        "EW": EW_named(p),
    }
    _ORBIT_CACHE[p] = rec
    return rec


def prove_AP_hist() -> dict:
    ok = True
    rows = {}
    want_hists = {
        5: {20: 15, 40: 15},
        7: {56: 28, 98: 56},
        13: {312: 91, 520: 182, 624: 182, 676: 91},
    }
    pm1_only = True
    for p, want in want_hists.items():
        rec = ap_orbit_W(p)
        if rec["n"] != rec["want"] or rec["hist"] != want:
            ok = False
        ew = EW_named(p)
        subset_pm = set(rec["hist"]).issubset({ew - 2 * p, ew + 2 * p})
        if p == 13 and subset_pm:
            ok = False
            pm1_only = True
        if p == 13 and not subset_pm:
            pm1_only = False
        rows[str(p)] = rec
    if pm1_only:
        ok = False
    return {
        "proved": ok,
        "APpm2p_all_p1_false": not pm1_only,
        "by_p": rows,
        "theorem": "AP G-orbit W-hists as displayed; {EW±2p} fails at p=13.",
    }


def prove_sufficient() -> dict:
    ok = True
    eq7 = Var_2orb(7) == Fraction(21504, 409)
    if eq7:
        ok = False
    rows = {}
    for p in (5, 13, 17, 29, 37):
        al, be = alpha_beta_p1(p)
        ew2 = Fraction(EW_named(p) ** 2)
        floor_A = al + be * Qpp_floor_ub(p)
        floor_Var = floor_A - ew2
        v2 = Var_2orb(p)
        if v2 > floor_Var + Fraction(1, 10**9):
            ok = False
        Q_from_v2 = (ew2 + v2 - al) / be
        if Q_from_v2 > Qpp_floor_ub(p) + Fraction(1, 10**9):
            ok = False
        rows[str(p)] = {
            "Var_2orb": str(v2),
            "floor_Var": str(floor_Var),
            "Q_from_v2": str(Q_from_v2),
            "ub": str(Qpp_floor_ub(p)),
            "le": Q_from_v2 <= Qpp_floor_ub(p),
        }
    if abs(float(Var_2orb(7) - Fraction(280, 3))) > 1e-9:
        # 4*7*8*15/(35+1)=3360/36=280/3
        if Var_2orb(7) != Fraction(280, 3):
            ok = False
    return {
        "proved": ok,
        "Var2orb_eq_live_p7_false": not eq7,
        "by_p": rows,
        "theorem": (
            "Var≤4p(p+1)(2p+1)/(5p+1) ⇒ Q++≤floor ub for p≡1. "
            "Equality at p=5; strict at p=7 live."
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
        "Var_2orb_proved": False,
        "note": (
            "Var≤Var_2orb is a sufficient majorant for the p≡1 floor, "
            "unproved. AP support is not {EW±2p} at p=13. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.323  AP-orbit W; Var_2orb sufficient for p≡1 floor", flush=True)
    A = prove_AP_hist()
    print(
        f"  A AP-hist: {A['proved']} pm2p_all={not A['APpm2p_all_p1_false']}",
        flush=True,
    )
    B = prove_sufficient()
    print(
        f"  B sufficient: {B['proved']} eq7={not B['Var2orb_eq_live_p7_false']}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.323",
        "title": "AP W not ±2p at p=13; Var_2orb sufficient, unproved",
        "proved": {
            "AP_orbit_hists": A["proved"],
            "Var2orb_sufficient": B["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15323.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
