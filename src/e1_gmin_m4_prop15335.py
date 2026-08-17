#!/usr/bin/env python3
"""
Prop 15.335 — QR0⊙NC coarse Q_τ via ẑ−2Γ at p=5,7,11,13.
Those values are not the full-ensemble Q_τ (dens 15,9,495,3003
≠ |H_+|/(2p) ∈ {13,409}).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.334 flags.  Floor stays OPEN.  Still 15.x.

============================================================================
Setup.  y = QR0⊙NC (15.334).  Q(r)=|Ω|⁻¹∑_ξ u(ξ)u(rξ), u=|ẑ|²,
reported /q² on coarse 15.290/15.300 classes.

============================================================================
Theorem A — PROVED (15.312 ẑ−2Γ; certified p=5,7,11,13).
  Q^{QR0⊙NC} equals the Ω-average of |ẑ_{QR0}−2Γ|².  Dropping the
  2 fails (err1 ≫ 0).  Class values:
      p=5:  Q++=64/15,  Q_N*=16/5     (coincides with ystar)
      p=7:  Q++=40/9,   Q_N*=76/21
      p=11: Q++=3976/495, Q_N*=144/55
      p=13: Q++=22688/3003, Q_N*=2992/1365
  Fail-when-wrong: drop the 2.  ∎

Theorem B — PROVED (denominator / live ensemble).
  These are not the H_+ ensemble Q_τ.  At p=5, 64/15 ≠ 48/13;
  at p=7, 40/9 ≠ 1544/409.  Dens {15,9,495,3003} ≠ {13,409}.
  Fail-when-wrong: claim QR0⊙NC Q++ equals the live ensemble.  ∎

Theorem C — OPEN.  Full-ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: first_nc + Ω Fourier.  Writes
evidence/e1_gmin_m4_prop15335.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15312 import Q_switch, qr0_nc_pack  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


def as_frac(x, cap=20000) -> Fraction:
    return Fraction(x).limit_denominator(cap)


@lru_cache(maxsize=8)
def Q_qr0nc(p: int) -> dict:
    pack = qr0_nc_pack(p)
    if pack is None or "y" not in pack:
        return {"found": False, "p": p}
    Q, err2, err1 = Q_switch(p, pack)
    return {
        "found": True,
        "p": p,
        "Q": {k: as_frac(v) for k, v in Q.items()},
        "err2": float(err2),
        "err1": float(err1),
        "t": pack.get("t"),
        "c": pack.get("c"),
    }


def prove_A() -> dict:
    expect = {
        5: {"++": Fraction(64, 15), "N*": Fraction(16, 5)},
        7: {"++": Fraction(40, 9), "N*": Fraction(76, 21)},
        11: {"++": Fraction(3976, 495), "N*": Fraction(144, 55)},
        13: {"++": Fraction(22688, 3003), "N*": Fraction(2992, 1365)},
    }
    ok = True
    rows = {}
    for p, exp in expect.items():
        rec = Q_qr0nc(p)
        if not rec["found"] or rec["err2"] >= 1e-9 or rec["err1"] <= 1e-3:
            ok = False
        hit = rec["found"] and rec["Q"]["++"] == exp["++"] and rec["Q"]["N*"] == exp["N*"]
        if not hit:
            ok = False
        rows[str(p)] = {
            "Qpp": str(rec["Q"]["++"]) if rec["found"] else None,
            "Qn": str(rec["Q"].get("N*")) if rec["found"] else None,
            "err2": rec.get("err2"),
            "err1": rec.get("err1"),
            "match": hit,
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "QR0⊙NC Q_τ = Ω-avg |ẑ_QR0−2Γ|².  Drop 2 fails.",
    }


def prove_B() -> dict:
    q5 = Q_qr0nc(5)["Q"]["++"]
    q7 = Q_qr0nc(7)["Q"]["++"]
    ok = q5 != LIVE_QPP[5] and q7 != LIVE_QPP[7]
    ok = ok and q5 == Fraction(64, 15) and q7 == Fraction(40, 9)
    ok = ok and q5.denominator != 13 and q7.denominator != 409
    return {
        "proved": ok,
        "qr0nc_p5": str(q5),
        "ensemble_p5": str(LIVE_QPP[5]),
        "qr0nc_p7": str(q7),
        "ensemble_p7": str(LIVE_QPP[7]),
        "dens_not_Hplus": True,
        "theorem": (
            "QR0⊙NC Q++ is 64/15 and 40/9, not the ensemble "
            "48/13 and 1544/409."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "Hplus_named_in_p": False,
        "note": (
            "QR0⊙NC type Q is certified, not the H_+ mixture.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.335  QR0⊙NC Q_τ ≠ ensemble Q_τ", flush=True)
    A = prove_A()
    print(f"  A ẑ−2Γ values: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B not ensemble: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.335",
        "title": "QR0⊙NC Q_τ named as rationals; not ensemble Q_τ",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "qr0nc_Q_via_zhat_minus_2Gamma": A["proved"],
            "not_ensemble_Q_tau": B["proved"],
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15335.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
